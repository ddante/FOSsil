MODULE petsc_driver

  USE Element_Class
  USE Geometry,       ONLY: N_dofs, N_seg, N_elements, elements

  IMPLICIT NONE

#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"

  !===================
  Mat  :: LHS
  KSP  :: ksp
  PC   :: pc
  Vec  :: b_rhs, x_sol
  !===================

  LOGICAL :: call_petsc = .FALSE.

  PRIVATE
  PUBLIC :: init_petsc, solve_sys, finalize_petsc
  PUBLIC :: LHS
  
CONTAINS

  !===========================
  SUBROUTINE init_petsc(N_eqn)
  !===========================

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: N_eqn
    !----------------------------

    INTEGER, DIMENSION(:), ALLOCATABLE :: nnz
    LOGICAL, DIMENSION(:), ALLOCATABLE :: seen_face
    INTEGER, DIMENSION(:), ALLOCATABLE :: idx

    INTEGER, DIMENSION(:), POINTER :: Nu

    INTEGER :: je, k, i
    !------------------------------------------

    PetscReal :: r_tol, a_tol, d_tol
    PetscInt :: max_res, max_ite

    PetscErrorCode :: ierr
    !------------------------------------------

    ! Count the #non-zero elements for each row
    ALLOCATE( nnz(N_dofs*N_eqn) ); nnz = 0
    ALLOCATE(  seen_face(N_seg) ); seen_face = .FALSE.

    DO je = 1, N_elements

       DO k = 1, elements(je)%p%N_faces

          Nu => elements(je)%p%faces(k)%f%NU

          IF( seen_face(elements(je)%p%faces(k)%f%g_seg) == .FALSE. ) THEN

             seen_face(elements(je)%p%faces(k)%f%g_seg) = .TRUE.

             ALLOCATE( idx(SIZE(Nu)) )
             
             DO i = 1, N_eqn

                idx = N_eqn*(Nu - 1) + i

                nnz(idx) = nnz(idx) + 1

             ENDDO

             DEALLOCATE( idx )
             NULLIFY( Nu )

          ENDIF

       ENDDO

    ENDDO

    ! diagonal term and system
    nnz = (nnz + 1)*N_eqn

    ! Init PETSC
    !---------------------------------------------
    CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)

    CALL MatCreateSeqAIJ(PETSC_COMM_WORLD,           &
                         N_eqn*N_dofs, N_eqn*N_dofs, &
                         PETSC_NULL_INTEGER, nnz,    &
                         LHS, ierr)
    
    CALL MatSetFromOptions(LHS, ierr)

    CALL VecCreate(PETSC_COMM_WORLD, &
                   b_rhs, ierr)

    CALL VecSetSizes(b_rhs, PETSC_DECIDE, &
                     N_eqn*N_dofs, ierr)

    CALL VecSetFromOptions(b_rhs, ierr)

    CALL VecDuplicate(b_rhs, x_sol, ierr)

    DEALLOCATE( nnz, seen_face )

    ! Solver option
    !---------------------------------------------
    r_tol = 1.0E-5; a_tol = 1.0E-17; d_tol = 1.0E5

    max_ite = 240; max_res = 60

    CALL KSPCreate(PETSC_COMM_WORLD, &
                   ksp, ierr)

    CALL KSPSetType(ksp, &
                    KSPGMRES, ierr)

    CALL KSPGMRESSetRestart(ksp, &
                            max_res, ierr)

    CALL KSPSetTolerances(ksp,                 &
                          r_tol, a_tol, d_tol, &
                          max_ite, ierr)

    CALL KSPSetOperators(ksp, LHS, &
                              LHS, DIFFERENT_NONZERO_PATTERN, ierr)

    
    CALL KSPSetFromOptions(ksp, ierr)

    ! Preconditioning
    !---------------------------------------------
    CALL KSPGetPC(ksp, pc, ierr)

    CALL PCSetType(pc, &
                       PCILU, ierr)

    CALL PCFactorSetReuseOrdering(pc, PETSC_TRUE)

    CALL PCFactorSetLevels(pc, 1)    

    call_petsc = .TRUE.

  END SUBROUTINE init_petsc
  !========================  

  !======================
  SUBROUTINE solve_sys(b)
  !======================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: b
    !------------------------------------------------

    REAL(KIND=8), DIMENSION(:), POINTER :: sol_v

    INTEGER,      DIMENSION(SIZE(b,1)) :: idx
    REAL(KIND=8), DIMENSION(SIZE(b,1)) :: vv
    
    INTEGER :: N_eqn, N_dof, i, j

    PetscErrorCode :: ierr
    !------------------------------------------------

    N_eqn = SIZE(b, 1)
    N_dof = SIZE(b, 2)
   
    ! Final assembly of the matrix
    CALL MatAssemblyBegin(LHS, MAT_FINAL_ASSEMBLY, ierr)
    CALL MatAssemblyEnd(  LHS, MAT_FINAL_ASSEMBLY, ierr)

    CALL MatSetOption(LHS,                               &
                           MAT_NEW_NONZERO_LOCATION_ERR, &
                      PETSC_TRUE, ierr)

    ! Assemby of rhs
    DO i = 1, N_dof

       idx = (/ (j, j = N_eqn*(i-1) + 1, i*N_eqn) /)
       idx = idx - 1

       vv = -b(:, i)

       CALL VecSetValues(b_rhs,                 &
                                N_eqn, idx, vv, &
                         INSERT_VALUES, ierr)

    ENDDO

    CALL VecAssemblyBegin(b_rhs, ierr)
    CALL VecAssemblyEnd(  b_rhs, ierr)

    ! Linear solver
    CALL KSPSetUp(ksp, ierr)

    CALL KSPSolve(ksp, b_rhs, x_sol, ierr)

    ! Copy the solution on the original array
    CALL VecGetArrayF90(x_sol, sol_v, ierr)
    DO i = 1, N_dof

       idx = (/ (j, j = N_eqn*(i-1) + 1, i*N_eqn) /)

       b(:, i) = sol_v(idx)

    ENDDO
    CALL VecRestoreArrayF90(x_sol, sol_v, ierr)

    ! Reset matrix
    CALL MatZeroEntries(LHS, ierr)
 
  END SUBROUTINE solve_sys
  !=======================

  !==========================
  SUBROUTINE finalize_petsc()
  !==========================

    IMPLICIT NONE

    PetscErrorCode :: ierr
    !--------------------------

    IF( call_petsc) THEN

       CALL MatDestroy(LHS, ierr)

       CALL VecDestroy(b_rhs, ierr)
       CALL VecDestroy(x_sol, ierr)

       CALL KSPDestroy(ksp, ierr)

       CALL PetscFinalize(ierr)

    ENDIF

  END SUBROUTINE finalize_petsc
  !============================
  
END MODULE petsc_driver
