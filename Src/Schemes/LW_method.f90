MODULE LW_method

  USE Element_class
  USE geometry,         ONLY: N_dim, elements
  USE init_problem,     ONLY: pb_type
  USE models,           ONLY: advection_speed
  USE init_problem,     ONLY: N_eqn, pb_type, visc
  USE Lin_Algebra,      ONLY: inverse
  USE FOS_system

  USE petsc_driver,   ONLY: LHS

  USE Quadrature_rules

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

  PUBLIC :: LW_scheme, LW_scheme_imp
  
CONTAINS

  !===================================================
  SUBROUTINE LW_scheme(ele, Phi_tot, u, Phi_i, inv_dt)
  !===================================================

    IMPLICIT NONE

    TYPE(element),                INTENT(IN)  :: ele
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: Phi_tot
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: u
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: Phi_i
    REAL(KIND=8),                 INTENT(OUT) :: inv_dt
    !--------------------------------------------------

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
    REAL(KIND=8), DIMENSION(N_dim, N_eqn) :: D_u_q
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: Stab_left_i
    REAL(KIND=8), DIMENSION(N_eqn) :: Stab_right
    REAL(KIND=8), DIMENSION(N_eqn) :: u_q
    !--------------------------------------------------

    REAL(KIND=8), DIMENSION(N_dim) :: v_n, a_m, a_q

    REAL(KIND=8) :: x_m, y_m, u_m, r_N_dofs, mod_n, l_max

    INTEGER :: i, j, iq, id, N_dofs, N_verts
    !--------------------------------------------------
    
    N_dofs = ele%N_points; r_N_dofs = REAL(N_dofs, 8)
    N_verts = ele%N_verts
    
    !---------------------
    ! Central distribution 
    !--------------------------------------------------
    DO j = 1, N_eqn

       Phi_i(j, :) = Phi_tot(j)/r_N_dofs

    ENDDO

    !---------------
    ! Scaling matrix
    !---------------------------------------------------
    Tau_ = 0.d0;  inv_dt = 0.d0

    u_m = SUM(u(1,:)) / r_N_dofs

    x_m = SUM(ele%Coords(1, :)) / r_N_dofs
    y_m = SUM(ele%Coords(2, :)) / r_N_dofs

    a_m = advection_speed(pb_type, u_m, x_m, y_m)

    DO i = 1, N_dofs
       
       mod_n = DSQRT( SUM(ele%rd_n(:,i)**2) )
       v_n = ele%rd_n(:,i)/mod_n

       lambda = FOS_eigenvalues(a_m, visc, v_n)
       RR = FOS_right_eigenvectors(a_m, visc, v_n)
       LL = FOS_left_eigenvectors(a_m, visc, v_n)

       PP = 0.d0; l_max = 0.d0
       DO j = 1, N_eqn

         !PP(j, j) = MAX(0.d0, lambda(j))*mod_n
          PP(j, j) = ABS(lambda(j))*mod_n

          l_max = MAX( l_max, &
                       MAX(0.d0, lambda(j))*mod_n )

       ENDDO

       Tau_ = Tau_ + 0.5d0 * MATMUL(RR, MATMUL(PP, LL))

       ! Time step
      !inv_dt = MAX( inv_dt,  MAXVAL(ABS(PP)) )
       inv_dt = MAX(inv_dt, l_max)

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

       u_q = 0.d0; D_u_q = 0.d0

       DO i = 1, N_dofs

          u_q = u_q + phi_q(i, iq)*u(:,i)

          DO id = 1, N_dim
             D_u_q(id, :) = D_u_q(id, :) + &
                            D_phi_q(id, i, iq) * u(:, i)
          ENDDO

       ENDDO !  # dofs -> i

       a_q = advection_speed(pb_type, u_q(1), xy(1,iq), xy(2,iq))
       AA_q = FOS_Jacobian(a_q, visc)

       Stab_right = 0.d0
       DO id = 1, N_dim
          Stab_right = Stab_right + &
                       MATMUL( AA_q(:,:, id), D_u_q(id, :) )
       ENDDO

       Stab_right = Stab_right - &
                    FOS_source(pb_type, u_q, visc, xy(1,iq), xy(2,iq))
     
       DO i = 1, N_dofs

          Stab_left_i = 0.d0
          DO id = 1, N_dim
             Stab_left_i = Stab_left_i + &
                           AA_q(:,:, id)*D_phi_q(id, i, iq)
          ENDDO

          Phi_i(:,i) = Phi_i(:,i) + w(iq) * &
                       MATMUL( Stab_left_i, &
                               MATMUL(Tau, Stab_right) )

       ENDDO !  # dofs -> i

    ENDDO ! # quad -> iq

    NULLIFY( D_phi_q, phi_q, w, xy )

  END SUBROUTINE LW_scheme
  !=======================

  !==============================================================
  SUBROUTINE LW_scheme_imp(ele, Phi_tot, u, J_tot, Phi_i, inv_dt)
  !==============================================================

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
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: PP
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: Tau_, Tau
    !--------------------------------------------------

    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: D_phi_q
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: phi_q
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: xy
    REAL(KIND=8), DIMENSION(:),     POINTER :: w  
    !--------------------------------------------------

    REAL(KIND=8), DIMENSION(N_eqn, N_eqn, N_dim) :: AA_q
    REAL(KIND=8), DIMENSION(N_dim, N_eqn) :: D_u_q
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: Stab_left_i
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: Stab_right_k
    REAL(KIND=8), DIMENSION(N_eqn) :: Stab_right
    REAL(KIND=8), DIMENSION(N_eqn) :: u_q
    !---------------------------------------------------

    INTEGER, DIMENSION(N_eqn) :: idr

    INTEGER, DIMENSION(:), ALLOCATABLE :: idc, idm

    REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: J_LW
    REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: MM
    !----------------------------------------------------

    REAL(KIND=8), DIMENSION(N_dim) :: v_n, a_m, a_q

    REAL(KIND=8) :: x_m, y_m, u_m, r_N_dofs, mod_n, l_max

    INTEGER :: i, j, k, iq, id, N_dofs, N_verts

    PetscErrorCode :: ierr
    !--------------------------------------------------
    
    N_dofs = ele%N_points; r_N_dofs = REAL(N_dofs, 8)
    N_verts = ele%N_verts

    ALLOCATE( J_LW(N_eqn, N_eqn, N_dofs, N_dofs) )

    ALLOCATE( idc(N_dofs*N_eqn), &
              idm(N_eqn) )
   
    ALLOCATE( MM(N_eqn, N_eqn*N_dofs) )

    !---------------------
    ! Central distribution 
    !--------------------------------------------------
    DO i = 1, N_dofs

       Phi_i(:, i) = Phi_tot/r_N_dofs

       !Jacobian
       DO k = 1, N_dofs
          J_LW(:,:, i, k) = J_tot(:,:, k)
       ENDDO       

    ENDDO

    !---------------
    ! Scaling matrix
    !---------------------------------------------------
    Tau_ = 0.d0;  inv_dt = 0.d0

    u_m = SUM(u(1,:)) / r_N_dofs

    x_m = SUM(ele%Coords(1, :)) / r_N_dofs
    y_m = SUM(ele%Coords(2, :)) / r_N_dofs

    a_m = advection_speed(pb_type, u_m, x_m, y_m)

    DO i = 1, N_dofs
       
       mod_n = DSQRT( SUM(ele%rd_n(:,i)**2) )
       v_n = ele%rd_n(:,i)/mod_n

       lambda = FOS_eigenvalues(a_m, visc, v_n)
       RR = FOS_right_eigenvectors(a_m, visc, v_n)
       LL = FOS_left_eigenvectors(a_m, visc, v_n)

       PP = 0.d0; l_max = 0.d0
       DO j = 1, N_eqn

         !PP(j, j) = MAX(0.d0, lambda(j))*mod_n
          PP(j, j) = ABS(lambda(j))*mod_n

          l_max = MAX( l_max, &
                       MAX(0.d0, lambda(j))*mod_n )

       ENDDO

       Tau_ = Tau_ + 0.5d0 * MATMUL(RR, MATMUL(PP, LL))

       ! Time step
      !inv_dt = MAX( inv_dt,  MAXVAL(ABS(PP)) )
       inv_dt = MAX(inv_dt, l_max)

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

       u_q = 0.d0; D_u_q = 0.d0

       DO i = 1, N_dofs

          u_q = u_q + phi_q(i, iq)*u(:,i)

          DO id = 1, N_dim
             D_u_q(id, :) = D_u_q(id, :) + &
                            D_phi_q(id, i, iq) * u(:, i)
          ENDDO

       ENDDO !  # dofs -> i

       a_q = advection_speed(pb_type, u_q(1), xy(1,iq), xy(2,iq))
       AA_q = FOS_Jacobian(a_q, visc)

       Stab_right = 0.d0
       DO id = 1, N_dim
          Stab_right = Stab_right + &
                       MATMUL( AA_q(:,:, id), D_u_q(id, :) )
       ENDDO

       Stab_right = Stab_right - &
                    FOS_source(pb_type, u_q, visc, xy(1,iq), xy(2,iq))
     
       DO i = 1, N_dofs

          Stab_left_i = 0.d0
          DO id = 1, N_dim
             Stab_left_i = Stab_left_i + &
                           AA_q(:,:, id)*D_phi_q(id, i, iq)
          ENDDO

          Phi_i(:,i) = Phi_i(:,i) + w(iq) * &
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
                            FOS_J_Source(pb_type, u_q, visc, xy(1,iq), xy(2,iq))

             J_LW(:,:, i,k) = J_LW(:,:, i,k) + w(iq) * &
                              MATMUL( Stab_left_i, &
                                      MATMUL(Tau, Stab_right_k) )

          ENDDO          

       ENDDO !  # dofs -> i

    ENDDO ! # quad -> iq

    NULLIFY( D_phi_q, phi_q, w, xy )

    ! Assembly Jacobian
    DO i = 1, N_dofs

       idr = N_eqn*(ele%Nu(i) - 1) + (/ (j, j = 1, N_eqn) /)
       idr = idr - 1

       DO k = 1, N_dofs

          idm = (/ (j, j = N_eqn*(k-1)+1, N_eqn*k) /)

          idc(idm) = N_eqn*(ele%Nu(k) - 1) + (/ (j, j = 1, N_eqn) /)
          idc(idm) = idc(idm) - 1

          MM(:, idm) = J_LW(:,:, i, k)

       ENDDO

       DO j = 1, N_eqn

          CALL MatSetValues(LHS,     1,            idr(j), &
                                     N_eqn*N_dofs, idc,    &
                            MM(j,:), ADD_VALUES, ierr)

       ENDDO

    ENDDO

    DEALLOCATE(idc, idm, MM, J_LW)

  END SUBROUTINE LW_scheme_imp
  !===========================

END MODULE LW_method
