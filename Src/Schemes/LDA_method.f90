MODULE LDA_method

  USE Element_class
  USE geometry,       ONLY: N_dim, elements
  USE init_problem,   ONLY: pb_type
  USE models,         ONLY: advection_speed
  USE init_problem,   ONLY: N_eqn, visc
  USE Lin_Algebra,    ONLY: inverse
  USE FOS_system

  USE petsc_driver,   ONLY: LHS

  IMPLICIT NONE

#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"

  PRIVATE

  PUBLIC :: LDA_scheme, LDA_scheme_imp
  
CONTAINS

  !====================================================
  SUBROUTINE LDA_scheme(ele, Phi_tot, u, Phi_i, inv_dt)
  !====================================================

    IMPLICIT NONE

    TYPE(element),                INTENT(IN)  :: ele
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: Phi_tot
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: u
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: Phi_i
    REAL(KIND=8),                 INTENT(OUT) :: inv_dt
    !---------------------------------------------------

    REAL(KIND=8), DIMENSION(N_eqn)        :: lambda
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: RR, LL
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: Beta
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: S_Kp, PP
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: inv_S_Kp

    REAL(KIND=8), DIMENSION(N_eqn, N_eqn, SIZE(u,2)) :: Kp
    !---------------------------------------------------

    REAL(KIND=8), DIMENSION(N_dim) :: v_n, a_m

    REAL(KIND=8) :: x_m, y_m, u_m, mod_n, r_Ns

    INTEGER :: Ns, i, j
    !---------------------------------------------------

    Ns = ele%N_points; r_Ns = REAL(Ns, 8)

    IF( Ns /= 3 ) THEN
       WRITE(*,*) 'ERROR: LDA with P1 triangles only'
       STOP
    ENDIF
    
    !-----------
    ! Mean state
    !-------------------------------------------------
    u_m = SUM(u(1,:)) / r_Ns

    x_m = SUM(ele%Coords(1, :)) / r_Ns
    y_m = SUM(ele%Coords(2, :)) / r_Ns

    a_m = advection_speed(pb_type, u_m, x_m, y_m)

    !------------
    ! LDA scheme
    !-------------------------------------------------
    S_Kp = 0.d0;  inv_dt = 0.d0

    DO i = 1, Ns

       mod_n = SQRT(SUM(ele%rd_n(:, i)**2))

       v_n = ele%rd_n(:, i)/mod_n
       
       lambda = FOS_eigenvalues(a_m, visc, v_n)

       RR = FOS_right_eigenvectors(a_m, visc, v_n)
       LL = FOS_left_eigenvectors(a_m, visc, v_n)

       PP = 0.d0
       DO j = 1, N_eqn
          PP(j, j) = MAX(0.d0, lambda(j))*mod_n
       ENDDO

       Kp(:,:, i) = 0.5d0 * MATMUL( RR, MATMUL(PP, LL) )

       S_Kp = S_Kp + Kp(:,:, i)

       ! Time step
       !inv_dt = MAX(inv_dt,  0.5d0*PP(2,2))
       inv_dt = MAX(inv_dt,  MAXVAL(ABS(PP)))
       
    ENDDO

    inv_S_Kp = inverse(S_Kp)
    
    DO i = 1, Ns

       Beta = MATMUL( Kp(:,:,i), inv_S_Kp )

       Phi_i(:,i) = MATMUL( Beta, Phi_tot )

    ENDDO
    
  END SUBROUTINE LDA_scheme
  !========================
  
  !===============================================================
  SUBROUTINE LDA_scheme_imp(ele, Phi_tot, u, J_tot, Phi_i, inv_dt)
  !===============================================================

    IMPLICIT NONE

    TYPE(element),                  INTENT(IN)  :: ele
    REAL(KIND=8), DIMENSION(:),     INTENT(IN)  :: Phi_tot
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: u
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: J_tot
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: Phi_i
    REAL(KIND=8),                   INTENT(OUT) :: inv_dt
    !----------------------------------------------------

    REAL(KIND=8), DIMENSION(N_eqn)        :: lambda
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: RR, LL
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: Beta
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: S_Kp, PP
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: inv_S_Kp

    REAL(KIND=8), DIMENSION(N_eqn, N_eqn, SIZE(u,2)) :: Kp

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: MM
    !------------------------------------------------------

    REAL(KIND=8), DIMENSION(N_dim) :: v_n, a_m

    REAL(KIND=8) :: x_m, y_m, u_m, mod_n, r_Ns
    !------------------------------------------------------

    INTEGER, DIMENSION(N_eqn) :: idr

    INTEGER, DIMENSION(:), ALLOCATABLE :: idc, idm

    INTEGER :: Ns, i, j, k

    PetscErrorCode :: ierr
    !------------------------------------------------------

    Ns = ele%N_points; r_Ns = REAL(Ns, 8)

    ALLOCATE( idc(Ns*N_eqn), &
              idm(N_eqn) )
   
    ALLOCATE( MM(N_eqn, N_eqn*Ns) )

    IF( Ns /= 3 ) THEN
       WRITE(*,*) 'ERROR: LDA with P1 triangles only'
       STOP
    ENDIF
    
    !-----------
    ! Mean state
    !-------------------------------------------------
    u_m = SUM(u(1,:)) / r_Ns

    x_m = SUM(ele%Coords(1, :)) / r_Ns
    y_m = SUM(ele%Coords(2, :)) / r_Ns

    a_m = advection_speed(pb_type, u_m, x_m, y_m)

    !------------
    ! LDA scheme
    !-------------------------------------------------
    S_Kp = 0.d0;  inv_dt = 0.d0

    DO i = 1, Ns

       mod_n = SQRT(SUM(ele%rd_n(:, i)**2))

       v_n = ele%rd_n(:, i)/mod_n
       
       lambda = FOS_eigenvalues(a_m, visc, v_n)

       RR = FOS_right_eigenvectors(a_m, visc, v_n)
       LL = FOS_left_eigenvectors(a_m, visc, v_n)

       PP = 0.d0
       DO j = 1, N_eqn
          PP(j, j) = MAX(0.d0, lambda(j))*mod_n
       ENDDO

       Kp(:,:, i) = 0.5d0 * MATMUL( RR, MATMUL(PP, LL) )

       S_Kp = S_Kp + Kp(:,:, i)

       ! Time step
       !inv_dt = MAX(inv_dt,  0.5d0*PP(2,2))
       inv_dt = MAX(inv_dt,  MAXVAL(ABS(PP)))
       
    ENDDO

    inv_S_Kp = inverse(S_Kp)
    
    DO i = 1, Ns

       Beta = MATMUL( Kp(:,:,i), inv_S_Kp )

       Phi_i(:,i) = MATMUL( Beta, Phi_tot )

       !----------------
       ! Jacobian matrix
       !-------------------------------------------------
       idr = N_eqn*(ele%Nu(i) - 1) + (/ (j, j = 1, N_eqn) /)
       idr = idr - 1
       
       DO k = 1, Ns

          idm = (/ (j, j = N_eqn*(k-1)+1, N_eqn*k) /)

          idc(idm) = N_eqn*(ele%Nu(k) - 1) + (/ (j, j = 1, N_eqn) /)
          idc(idm) = idc(idm) - 1

          MM(:, idm) = MATMUL( Beta, J_tot(:,:, k) )

       ENDDO

       DO j = 1, N_eqn

          CALL MatSetValues(LHS,     1,        idr(j), &
                                     N_eqn*Ns, idc,    &
                            MM(j,:), ADD_VALUES, ierr)

       ENDDO

    ENDDO

    DEALLOCATE(idc, idm, MM)
    
  END SUBROUTINE LDA_scheme_imp
  !============================


END MODULE LDA_method
