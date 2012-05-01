MODULE LLxFS_method

  USE Element_class
  USE geometry,       ONLY: N_dim, elements
  USE init_problem,   ONLY: pb_type
  USE models,         ONLY: advection_speed
  USE init_problem,   ONLY: N_eqn, pb_type, visc
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
  
  PUBLIC :: LLxFS_scheme, LLxFS_scheme_imp

CONTAINS

  !======================================================
  SUBROUTINE LLxFS_scheme(ele, Phi_tot, uu, Phi_i, alpha)
  !======================================================

    IMPLICIT NONE

    TYPE(element),                INTENT(IN)  :: ele
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: Phi_tot
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: uu
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: Phi_i
    REAL(KIND=8),                 INTENT(OUT) :: alpha
    !--------------------------------------------------

    REAL(KIND=8), DIMENSION(N_eqn) :: lambda
    !--------------------------------------------------

    REAL(KIND=8), DIMENSION(N_dim) :: v_n, a_m

    REAL(KIND=8), DIMENSION(N_eqn) :: uu_m

    REAL(KIND=8) :: x_m, y_m, r_N_dofs, mod_n, l_max

    INTEGER :: i, j, iq, id, N_dofs, N_verts
    !--------------------------------------------------

    N_dofs = ele%N_points; r_N_dofs = REAL(N_dofs, 8)
    N_verts = ele%N_verts

    !-----------
    ! LxF Scheme
    !----------------------------------------
    uu_m = 0.d0
    DO i = 1, N_dofs
       uu_m = uu_m + uu(:, i)
    ENDDO
    uu_m = uu_m / r_N_dofs

    x_m = SUM(ele%Coords(1, :)) / r_N_dofs
    y_m = SUM(ele%Coords(2, :)) / r_N_dofs

    a_m = advection_speed(pb_type, uu_m(1), x_m, y_m)

    !v_n = uu_m(2:3)/DSQRT(SUM(uu_m(2:3)**2))
    v_n = a_m/DSQRT(SUM(a_m**2))

    lambda = FOS_eigenvalues(a_m, visc, v_n)
    
    alpha = 0.d0
    DO i = 1, N_verts

       mod_n = DSQRT(SUM(ele%rd_n(:, i)**2))

       alpha = MAX( alpha, &
                    ABS(MAXVAL(lambda)*mod_n) )
       
    ENDDO
    alpha = 3.d0*alpha
    
    DO i = 1, N_dofs
       Phi_i(:, i) = Phi_tot/r_N_dofs + alpha*( uu(:, i) - uu_m )
    ENDDO

    !---------
    ! Limiting
    !-----------------------------------
    CALL Limitation(uu_m, a_m, Phi_i)

    !--------------
    ! Stabilization
    !-------------------------------------
    Phi_i = Phi_i + Stabilization(ele, uu)

  END SUBROUTINE LLxFS_scheme
  !==========================

  !=================================================================
  SUBROUTINE LLxFS_scheme_imp(ele, Phi_tot, uu, J_tot, Phi_i, alpha)
  !=================================================================

    IMPLICIT NONE

    TYPE(element),                  INTENT(IN)  :: ele
    REAL(KIND=8), DIMENSION(:),     INTENT(IN)  :: Phi_tot
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: uu
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: J_tot
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: Phi_i
    REAL(KIND=8),                   INTENT(OUT) :: alpha
    !---------------------------------------------------   

    REAL(KIND=8), DIMENSION(N_eqn) :: lambda
    !--------------------------------------------------

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Stab

    REAL(KIND=8), DIMENSION(N_dim) :: v_n, a_m

    REAL(KIND=8), DIMENSION(N_eqn) :: uu_m

    REAL(KIND=8) :: x_m, y_m, r_N_dofs, &
                    mod_n, l_max, const

    INTEGER :: i, j, k, iq, id, N_dofs, N_verts
    !--------------------------------------------------

    INTEGER, DIMENSION(N_eqn) :: idr

    INTEGER, DIMENSION(:), ALLOCATABLE :: idc, idm

    REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: J_LxF
    REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: MM

    PetscErrorCode :: ierr
    !----------------------------------------------------

    N_dofs = ele%N_points; r_N_dofs = REAL(N_dofs, 8)
    N_verts = ele%N_verts

    ALLOCATE( Stab(N_eqn, N_dofs) )

    ALLOCATE( J_LxF(N_eqn, N_eqn, N_dofs, N_dofs) )

    ALLOCATE( idc(N_dofs*N_eqn), &
              idm(N_eqn) )
   
    ALLOCATE( MM(N_eqn, N_eqn*N_dofs) )

    !-----------
    ! LxF Scheme
    !----------------------------------------
    uu_m = 0.d0
    DO i = 1, N_dofs
       uu_m = uu_m + uu(:, i)
    ENDDO
    uu_m = uu_m / r_N_dofs

    x_m = SUM(ele%Coords(1, :)) / r_N_dofs
    y_m = SUM(ele%Coords(2, :)) / r_N_dofs

    a_m = advection_speed(pb_type, uu_m(1), x_m, y_m)

    !v_n = uu_m(2:3)/DSQRT(SUM(uu_m(2:3)**2))
    v_n = a_m/DSQRT(SUM(a_m**2))

    lambda = FOS_eigenvalues(a_m, visc, v_n)
    
    alpha = 0.d0
    DO i = 1, N_verts

       mod_n = DSQRT(SUM(ele%rd_n(:, i)**2))

       alpha = MAX( alpha, &
                    ABS(MAXVAL(lambda)*mod_n) )
       
    ENDDO
    alpha = 3.d0*alpha
    
    DO i = 1, N_dofs

       Phi_i(:, i) = Phi_tot/r_N_dofs + alpha*( uu(:, i) - uu_m )

       !Jacobian
       DO k = 1, N_dofs
          
          IF( k == i ) THEN
             const = alpha*(r_N_dofs - 1.d0)/r_N_dofs
          ELSE
             const = -alpha/r_N_dofs
          ENDIF

          J_LxF(:,:, i, k) = J_tot(:,:, k)/r_N_dofs + const

       ENDDO 

    ENDDO

    !---------
    ! Limiting
    !-----------------------------------
    CALL Limitation(uu_m, a_m, Phi_i)

    !--------------
    ! Stabilization
    !-------------------------------------------
    CALL Stabilization_imp(ele, uu, J_LxF, Stab)

    Phi_i = Phi_i + Stab

    ! Assembly Jacobian
    DO i = 1, N_dofs

       idr = N_eqn*(ele%Nu(i) - 1) + (/ (j, j = 1, N_eqn) /)
       idr = idr - 1

       DO k = 1, N_dofs

          idm = (/ (j, j = N_eqn*(k-1)+1, N_eqn*k) /)

          idc(idm) = N_eqn*(ele%Nu(k) - 1) + (/ (j, j = 1, N_eqn) /)
          idc(idm) = idc(idm) - 1

          MM(:, idm) = J_LxF(:,:, i, k)

       ENDDO

       DO j = 1, N_eqn

          CALL MatSetValues(LHS,     1,            idr(j), &
                                     N_eqn*N_dofs, idc,    &
                            MM(j,:), ADD_VALUES, ierr)

       ENDDO

    ENDDO

    DEALLOCATE(idc, idm, MM, J_LxF, Stab)

  END SUBROUTINE LLxFS_scheme_imp
  !===============================

  !==================================
  SUBROUTINE Limitation(uu, a, Phi_i)
  !==================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: uu
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: a
    REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: Phi_i
    !----------------------------------------------------

    REAL(KIND=8) :: Phi_tot, Den

    REAL(KIND=8), DIMENSION(SIZE(Phi_i,2)) :: beta_p

    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: RR, LL

    REAL(KIND=8), DIMENSION(N_dim) :: v_n

    INTEGER :: i, j, N_dofs
    !----------------------------------------------------

    N_dofs = SIZE(Phi_i, 2)
    
    !v_n = uu(2:3)/DSQRT(SUM(uu(2:3)**2)
    v_n = a/DSQRT(SUM(a**2))

    RR = FOS_right_eigenvectors(a, visc, v_n)
    LL = FOS_left_eigenvectors(a, visc, v_n)

    ! Characteristic residual: 
    ! Projection on the left eigenvectors
    DO i = 1, N_dofs
       Phi_i(:, i) = MATMUL(LL, Phi_i(:, i))
    ENDDO

    ! Limitation on the characteristic space
    DO j = 1, N_eqn
       
       Phi_tot = SUM( Phi_i(j, :) )

       IF( ABS(Phi_tot) > 0.d0 ) THEN          

          beta_p = MAX( Phi_i(j,:)/Phi_tot, 0.d0 )

          Den = SUM(beta_p)

          Phi_i(j, :) = (beta_p/Den) * Phi_tot

       ELSE

          Phi_i(j, :) = 0.d0

       END IF         
       
    ENDDO
    
    ! Back to the physical residual
    DO i = 1, N_dofs
       Phi_i(:, i) = MATMUL(RR, Phi_i(:, i))
    ENDDO

  END SUBROUTINE Limitation
  !========================

  !===========================================
  FUNCTION Stabilization(ele, uu) RESULT(Stab)
  !===========================================

    IMPLICIT NONE

    TYPE(element),                INTENT(IN)  :: ele
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: uu

    REAL(KIND=8), DIMENSION(SIZE(uu,1), SIZE(uu,2)) :: Stab    
    !-------------------------------------------------------

    REAL(KIND=8), DIMENSION(N_eqn)        :: lambda
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: RR, LL
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: PP
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: Tau_, Tau
    !--------------------------------------------------

    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: D_phi_q
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: phi_q
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: xy
    REAL(KIND=8), DIMENSION(:),     POINTER :: w  
    !--------------------------------------------------
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn, N_dim) :: AA_q
    REAL(KIND=8), DIMENSION(N_dim, N_eqn) :: D_uu_q
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: Stab_left_i
    REAL(KIND=8), DIMENSION(N_eqn) :: Stab_right
    REAL(KIND=8), DIMENSION(N_eqn) :: uu_q
    !--------------------------------------------------

    REAL(KIND=8), DIMENSION(N_dim) :: v_n, a_m, a_q

    REAL(KIND=8) :: x_m, y_m, u_m, r_N_dofs, mod_n

    INTEGER :: i, j, iq, id, N_dofs, N_verts
    !--------------------------------------------------

    Stab = 0.d0

    N_dofs = ele%N_points; r_N_dofs = REAL(N_dofs, 8)
    N_verts = ele%N_verts
    
    !---------------
    ! Scaling matrix
    !---------------------------------------------------
    Tau_ = 0.d0

    u_m = SUM(uu(1,:)) / r_N_dofs

    x_m = SUM(ele%Coords(1, :)) / r_N_dofs
    y_m = SUM(ele%Coords(2, :)) / r_N_dofs

    a_m = advection_speed(pb_type, u_m, x_m, y_m)

    DO i = 1, N_dofs
       
       mod_n = DSQRT( SUM(ele%rd_n(:,i)**2) )
       v_n = ele%rd_n(:,i)/mod_n

       lambda = FOS_eigenvalues(a_m, visc, v_n)
       RR = FOS_right_eigenvectors(a_m, visc, v_n)
       LL = FOS_left_eigenvectors(a_m, visc, v_n)

       PP = 0.d0
       DO j = 1, N_eqn

         !PP(j, j) = MAX(0.d0, lambda(j))*mod_n
          PP(j, j) = ABS(lambda(j))*mod_n

       ENDDO

       Tau_ = Tau_ + 0.5d0 * MATMUL(RR, MATMUL(PP, LL))

    ENDDO

    Tau = 0.5d0 * ele%Volume * inverse(Tau_)

    !-------------------
    ! Stabilization term
    !---------------------------------------------------
    D_phi_q => ele%D_phi_q
      phi_q => ele%phi_q
          w => ele%w_q
         xy => ele%xx_q

    DO iq = 1, ele%N_quad

       uu_q = 0.d0; D_uu_q = 0.d0

       DO i = 1, N_dofs

          uu_q = uu_q + phi_q(i, iq)*uu(:,i)

          DO id = 1, N_dim
             D_uu_q(id, :) = D_uu_q(id, :) + &
                             D_phi_q(id, i, iq) * uu(:, i)
          ENDDO

       ENDDO !  # dofs -> i

       a_q = advection_speed(pb_type, uu_q(1), xy(1,iq), xy(2,iq))
       AA_q = FOS_Jacobian(a_q, visc)

       Stab_right = 0.d0
       DO id = 1, N_dim
          Stab_right = Stab_right + &
                       MATMUL( AA_q(:,:, id), D_uu_q(id, :) )
       ENDDO

       Stab_right = Stab_right - &
                    FOS_source(pb_type, uu_q, visc, xy(1,iq), xy(2,iq))
     
       DO i = 1, N_dofs

          Stab_left_i = 0.d0
          DO id = 1, N_dim
             Stab_left_i = Stab_left_i + &
                           AA_q(:,:, id)*D_phi_q(id, i, iq)
          ENDDO

          Stab(:,i) = Stab(:,i) + w(iq) * &
                      MATMUL( Stab_left_i, &
                              MATMUL(Tau, Stab_right) )

       ENDDO !  # dofs -> i

    ENDDO ! # quad -> iq

    NULLIFY( D_phi_q, phi_q, w, xy )

  END FUNCTION Stabilization
  !=========================

  !=================================================
  SUBROUTINE Stabilization_imp(ele, uu, J_LxF, Stab)
  !=================================================

    IMPLICIT NONE

    TYPE(element),                    INTENT(IN)    :: ele
    REAL(KIND=8), DIMENSION(:,:),     INTENT(IN)    :: uu
    REAL(KIND=8), DIMENSION(:,:,:,:), INTENT(INOUT) :: J_LxF
    REAL(KIND=8), DIMENSION(:,:),     INTENT(OUT)   :: Stab    
    !-------------------------------------------------------

    REAL(KIND=8), DIMENSION(N_eqn)        :: lambda
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: RR, LL
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: PP
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: Tau_, Tau
    !--------------------------------------------------

    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: D_phi_q
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: phi_q
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: xy
    REAL(KIND=8), DIMENSION(:),     POINTER :: w  
    !--------------------------------------------------
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn, N_dim) :: AA_q
    REAL(KIND=8), DIMENSION(N_dim, N_eqn) :: D_uu_q
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: Stab_left_i
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: Stab_right_k
    REAL(KIND=8), DIMENSION(N_eqn) :: Stab_right
    REAL(KIND=8), DIMENSION(N_eqn) :: uu_q
    !--------------------------------------------------

    REAL(KIND=8), DIMENSION(N_dim) :: v_n, a_m, a_q

    REAL(KIND=8) :: x_m, y_m, u_m, r_N_dofs, mod_n

    INTEGER :: i, j, k, iq, id, N_dofs, N_verts
    !--------------------------------------------------

    Stab = 0.d0

    N_dofs = ele%N_points; r_N_dofs = REAL(N_dofs, 8)
    N_verts = ele%N_verts
    
    !---------------
    ! Scaling matrix
    !---------------------------------------------------
    Tau_ = 0.d0

    u_m = SUM(uu(1,:)) / r_N_dofs

    x_m = SUM(ele%Coords(1, :)) / r_N_dofs
    y_m = SUM(ele%Coords(2, :)) / r_N_dofs

    a_m = advection_speed(pb_type, u_m, x_m, y_m)

    DO i = 1, N_dofs
       
       mod_n = DSQRT( SUM(ele%rd_n(:,i)**2) )
       v_n = ele%rd_n(:,i)/mod_n

       lambda = FOS_eigenvalues(a_m, visc, v_n)
       RR = FOS_right_eigenvectors(a_m, visc, v_n)
       LL = FOS_left_eigenvectors(a_m, visc, v_n)

       PP = 0.d0
       DO j = 1, N_eqn

         !PP(j, j) = MAX(0.d0, lambda(j))*mod_n
          PP(j, j) = ABS(lambda(j))*mod_n

       ENDDO

       Tau_ = Tau_ + 0.5d0 * MATMUL(RR, MATMUL(PP, LL))

    ENDDO

    Tau = 0.5d0 * ele%Volume * inverse(Tau_)

    !-------------------
    ! Stabilization term
    !---------------------------------------------------
    D_phi_q => ele%D_phi_q
      phi_q => ele%phi_q
          w => ele%w_q
         xy => ele%xx_q

    DO iq = 1, ele%N_quad

       uu_q = 0.d0; D_uu_q = 0.d0

       DO i = 1, N_dofs

          uu_q = uu_q + phi_q(i, iq)*uu(:,i)

          DO id = 1, N_dim
             D_uu_q(id, :) = D_uu_q(id, :) + &
                             D_phi_q(id, i, iq) * uu(:, i)
          ENDDO

       ENDDO !  # dofs -> i

       a_q = advection_speed(pb_type, uu_q(1), xy(1,iq), xy(2,iq))
       AA_q = FOS_Jacobian(a_q, visc)

       Stab_right = 0.d0
       DO id = 1, N_dim
          Stab_right = Stab_right + &
                       MATMUL( AA_q(:,:, id), D_uu_q(id, :) )
       ENDDO

       Stab_right = Stab_right - &
                    FOS_source(pb_type, uu_q, visc, xy(1,iq), xy(2,iq))
     
       DO i = 1, N_dofs

          Stab_left_i = 0.d0
          DO id = 1, N_dim
             Stab_left_i = Stab_left_i + &
                           AA_q(:,:, id)*D_phi_q(id, i, iq)
          ENDDO

          Stab(:,i) = Stab(:,i) + w(iq) * &
                      MATMUL( Stab_left_i, &
                              MATMUL(Tau, Stab_right) )

          !----------------
          ! Jacobian matrix
          !-------------------------------------------------
          DO k = 1, N_dofs

             Stab_right_k = 0.d0
             DO id = 1, N_dim
                 Stab_right_k = Stab_right_k + &
                                AA_q(:,:, id)*D_phi_q(id, k, iq)
             ENDDO

             Stab_right_k = Stab_right_k - &
                            phi_q(k, iq) * &
                            FOS_J_Source(pb_type, uu_q, visc, xy(1,iq), xy(2,iq))

             J_LxF(:,:, i,k) = J_LxF(:,:, i,k) + w(iq) * &
                               MATMUL( Stab_left_i, &
                                       MATMUL(Tau, Stab_right_k) )

          ENDDO  

       ENDDO !  # dofs -> i

    ENDDO ! # quad -> iq

    NULLIFY( D_phi_q, phi_q, w, xy )

  END SUBROUTINE Stabilization_imp  
  !===============================

END MODULE LLxFS_method
