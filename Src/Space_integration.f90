MODULE space_integration

  USE Element_Class

  USE geometry,       ONLY: N_dim, N_dofs, N_elements, elements

  USE init_problem,   ONLY: order, pb_type, bc_type, &
                            visc, CFL, N_eqn, CFL_l, R_0, R_1

  USE models,         ONLY: advection_flux, advection_speed, &
                            strong_bc, exact_solution, exact_grad

  USE Lin_Algebra,    ONLY: Identity_matrix

  USE petsc_driver,   ONLY: LHS

  USE FOS_system
  USE Num_scheme

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

  !==========================================
  INTEGER, PARAMETER :: STRONG_BC_TYPE = 0, &
                          WEAK_BC_TYPE = 1
  !==========================================


  PUBLIC :: compute_rhs, compute_lhs_rhs

CONTAINS

  !====================================
  SUBROUTINE compute_rhs(uu, rhs, Dt_V)
  !====================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: uu
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT)   :: rhs
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT)   :: Dt_V
    !-------------------------------------------------

    TYPE(element) :: loc_ele

    INTEGER,      DIMENSION(:),   ALLOCATABLE :: Nu
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: u1, Phi_i

    REAL(KIND=8), DIMENSION(N_eqn) :: Phi_tot
    REAL(KIND=8) :: inv_dt
    
    INTEGER :: Ns, je
    !----------------------------------------------

    rhs = 0.d0;  Dt_v = 0.d0

    DO je = 1, N_elements

       loc_ele = elements(je)%p

       Ns = loc_ele%N_points

       ALLOCATE( Nu(Ns) )
       ALLOCATE( u1(N_eqn, Ns), Phi_i(N_eqn, Ns) )                 

       Nu = loc_ele%NU

       u1 = uu(:, Nu)

       !-------------------------------
       ! Compute the total fluctuation
       !------------------------------------------------------
       Phi_tot = total_residual(loc_ele, u1)
    
       !---------------------------
       ! Distribute the fluctuation
       !------------------------------------------------------------
       CALL distribute_residual(loc_ele, Phi_tot, u1, Phi_i, inv_dt)

       ! Gather the nodal residual
       rhs(:, Nu) = rhs(:, Nu) + Phi_i
       
       Dt_V(Nu) = Dt_V(Nu) + inv_dt
       
       DEALLOCATE( Nu, u1, Phi_i )

    ENDDO

    Dt_V = CFL/DT_V
 
    !---------------------------
    ! Impose boundary conditions
    !-------------------------------------
    SELECT CASE(bc_type)

    CASE(STRONG_BC_TYPE)

       CALL strong_bc(pb_type, visc, uu, rhs)

    CASE(WEAK_BC_TYPE)

       rhs = rhs + Res_weak_bc(uu)

    CASE DEFAULT

       WRITE(*,*) 'ERROR: unknown boundary condition method'
       STOP

    END SELECT
    
  END SUBROUTINE compute_rhs
  !=========================

  !==================================
  SUBROUTINE compute_lhs_rhs(uu, rhs)
  !==================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: uu
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT)   :: rhs
    !-------------------------------------------------

    TYPE(element) :: loc_ele

    INTEGER,      DIMENSION(:),   ALLOCATABLE :: Nu
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: u1, Phi_i

    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: J_tot

    REAL(KIND=8), DIMENSION(N_eqn) :: Phi_tot

    REAL(KIND=8), DIMENSION(N_dofs) :: V_Dt

    REAL(KIND=8) :: inv_dt

    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: MM, II

    INTEGER, DIMENSION(N_eqn) :: idd
    
    INTEGER :: Ns, je, k, i

    PetscErrorCode :: ierr
    !----------------------------------------------

    rhs = 0.d0;  V_Dt = 0.d0

    DO je = 1, N_elements

       loc_ele = elements(je)%p

       Ns = loc_ele%N_points

       ALLOCATE( Nu(Ns) )
       ALLOCATE( u1(N_eqn, Ns), Phi_i(N_eqn, Ns) )
       ALLOCATE( J_tot(N_eqn, N_eqn, Ns) )

       Nu = loc_ele%NU

       u1 = uu(:, Nu)

       !-------------------------------
       ! Compute the total fluctuation
       !---------------------------------------------------
       CALL total_residual_imp(loc_ele, u1, Phi_tot, J_tot)

       !---------------------------
       ! Distribute the fluctuation
       !------------------------------------------------------------
        CALL distribute_residual_imp(loc_ele, Phi_tot, u1, J_tot, &
                                    Phi_i, inv_dt)

       ! Gather the nodal residual
       rhs(:, Nu) = rhs(:, Nu) + Phi_i
       
       V_Dt(Nu) = V_Dt(Nu) + inv_dt
       
       DEALLOCATE( Nu, u1, Phi_i, J_tot )

    ENDDO

    !---------------------------
    ! Impose boundary conditions
    !-------------------------------------
    SELECT CASE(bc_type)

    CASE(STRONG_BC_TYPE)

       !CALL strong_bc(pb_type, visc, uu, rhs)
       WRITE(*,*) 'ERROR: Strong BC not supported in implicit'
       STOP

    CASE(WEAK_BC_TYPE)

       rhs = rhs + Res_weak_bc_imp(uu)

    CASE DEFAULT

       WRITE(*,*) 'ERROR: unknown boundary condition method'
       STOP

    END SELECT

    !------------------------------
    ! Pseudo-transient continuation
    !-------------------------------------
    CFL_l = CFL_law(rhs, CFL_l)

    ! |C|/Dt
    V_Dt = V_Dt/CFL_l

    II = Identity_matrix(N_eqn)
    
    ! Jacobian diagonal contribution
    DO i = 1, N_dofs

       MM = II*V_Dt(i)
       idd = (/ (k, k = N_eqn*(i-1)+1, N_eqn*i) /)
       idd = idd - 1
       
       DO k = 1, N_eqn
       
          CALL MatSetValues(LHS,     1,     idd(k),  &
                                     N_eqn, idd,     &
                            MM(k,:), ADD_VALUES, ierr)

       ENDDO

    ENDDO
    
  END SUBROUTINE compute_lhs_rhs
  !==============================

  !-------------------------------------------------------
  !---------------------- EXPLICIT -----------------------
  !-------------------------------------------------------

  !===============================================
  FUNCTION total_residual(ele, uu) RESULT(Phi_tot)
  !===============================================

    IMPLICIT NONE

    TYPE(element),                INTENT(IN) :: ele
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: uu
    
    REAL(KIND=8), DIMENSION(N_eqn) :: Phi_tot
    !---------------------------------------------

    INTEGER,                      POINTER :: N_quad
    INTEGER,                      POINTER :: N_points
    INTEGER,      DIMENSION(:),   POINTER :: loc
    REAL(KIND=8), DIMENSION(:,:), POINTER :: p
    REAL(KIND=8), DIMENSION(:,:), POINTER :: n
    REAL(KIND=8), DIMENSION(:),   POINTER :: w
    REAL(KIND=8), DIMENSION(:,:), POINTER :: xy
    !---------------------------------------------

    REAL(KIND=8), DIMENSION(N_eqn) :: Phi_b, Phi_S
    REAL(KIND=8), DIMENSION(N_eqn) :: uu_q, S_q

    REAL(KIND=8), DIMENSION(N_eqn, N_dim) :: ff_q

    INTEGER :: i, j, iq, k, id
    !---------------------------------------------

    !------------------
    ! Boundary integral
    !----------------------------------------------
    phi_b = 0.d0

    DO j = 1, ele%N_faces

       N_quad   => ele%faces(j)%f%N_quad
       N_points => ele%faces(j)%f%N_points
       loc      => ele%faces(j)%f%l_nu
       p        => ele%faces(j)%f%phi_q
       w        => ele%faces(j)%f%w_q
       n        => ele%faces(j)%f%n_q
       xy       => ele%faces(j)%f%xx_q

       DO iq = 1, N_quad

          uu_q = 0.d0
          DO k = 1, N_points
             uu_q = uu_q + uu(:, loc(k)) * p(k, iq)
          ENDDO

          ff_q = FOS_advection_flux(pb_type, uu_q, visc, xy(1,iq), xy(2,iq))

          DO id = 1, N_dim
             phi_b = phi_b + w(iq)*ff_q(:,id)*n(id, iq)
          ENDDO

       ENDDO

       NULLIFY(N_quad, N_points, loc, p, w, n, xy)

    ENDDO
    !----------------------------------------------

    !----------------
    ! Domain integral
    !----------------------------------------------
    Phi_S = 0.d0

    p  => ele%phi_q
    w  => ele%w_q
    xy => ele%xx_q

    DO iq = 1, ele%N_quad

       uu_q = 0.d0
       DO k = 1, ele%N_points
          uu_q = uu_q + uu(:, k) * p(k, iq)
       ENDDO

       S_q = FOS_source(pb_type, uu_q, visc, xy(1,iq), xy(2,iq))

       Phi_S = Phi_S + w(iq)*S_q

    ENDDO
    
    NULLIFY(p, w, xy)
    !----------------------------------------------

    Phi_tot = Phi_b - Phi_S
                  !^^^

  END FUNCTION total_residual
  !==========================

  !=====================================
  FUNCTION Res_weak_bc(uu) RESULT(res_b)
  !=====================================
  !
  ! Impose the boundary conditions as 
  ! inflow/outflow type. The inflow/outflow
  ! state is the exact solution: (u, p, q)
  !
    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: uu

    REAL(KIND=8), DIMENSION(SIZE(uu,1),SIZE(uu,2)) :: res_b
    !------------------------------------------------------

    TYPE(element) :: ele

    INTEGER,                      POINTER :: N_quad
    INTEGER,                      POINTER :: N_points
    INTEGER,      DIMENSION(:),   POINTER :: loc
    INTEGER,      DIMENSION(:),   POINTER :: nu
    REAL(KIND=8), DIMENSION(:,:), POINTER :: p
    REAL(KIND=8), DIMENSION(:,:), POINTER :: n
    REAL(KIND=8), DIMENSION(:),   POINTER :: w
    REAL(KIND=8), DIMENSION(:,:), POINTER :: xy
    !---------------------------------------------------

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: u

    REAL(KIND=8), DIMENSION(N_eqn)        :: lambda
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: RR, LL
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: Neg

    REAL(KIND=8), DIMENSION(N_dim) :: a_q
    REAL(KIND=8), DIMENSION(N_eqn) :: u_ex, u_q, eta, alpha

    INTEGER :: je, jf, iq, k, j, i
    !---------------------------------------------------

    res_b = 0.d0

    DO je = 1, N_elements

       ele = elements(je)%p

       ALLOCATE( u(N_eqn, ele%N_points) )

       u = uu(:, ele%Nu)

       DO jf = 1, ele%N_faces

          ! Boundary face
          IF( ele%faces(jf)%f%c_ele == 0 ) THEN

             N_quad   => ele%faces(jf)%f%N_quad
             N_points => ele%faces(jf)%f%N_points
             nu       => ele%faces(jf)%f%nu
             loc      => ele%faces(jf)%f%l_nu
             p        => ele%faces(jf)%f%phi_q
             w        => ele%faces(jf)%f%w_q
             n        => ele%faces(jf)%f%n_q
             xy       => ele%faces(jf)%f%xx_q

             DO iq = 1, N_quad

                u_q = 0.d0
                DO k = 1, N_points
                   u_q = u_q + u(:, loc(k)) * p(k, iq)
                ENDDO

                a_q = advection_speed(pb_type, u_q(1), xy(1, iq), xy(2, iq))

                u_ex(1)   = exact_solution(pb_type, xy(:, iq), visc) ! u
                u_ex(2:3) = exact_grad(pb_type, xy(:, iq), visc)     ! p, q

                lambda = FOS_eigenvalues(a_q, visc, n(:,iq))
                RR = FOS_right_eigenvectors(a_q, visc, n(:,iq))
                LL = FOS_left_eigenvectors(a_q, visc, n(:,iq))

                Neg = 0.d0
                DO j = 1, N_eqn
                   Neg(j, j) = MIN(0.d0, lambda(j))
                ENDDO

                eta = MATMUL(LL, u_ex - u_q)
                
                alpha = MATMUL( RR, MATMUL(Neg, eta) )

                DO i = 1, N_points

                   res_b(:, Nu(i)) = res_b(:, Nu(i)) + &
                                     w(iq)*p(i, iq)*alpha

                ENDDO                

             ENDDO ! iq -> # quad

          ENDIF

       ENDDO ! jf -> # faces

       DEALLOCATE( u )

    ENDDO ! je -> # ele

  END FUNCTION Res_weak_bc 
  !=======================

  !-------------------------------------------------------
  !---------------------- IMPLICIT -----------------------
  !-------------------------------------------------------

  !=====================================================
  SUBROUTINE total_residual_imp(ele, uu, Phi_tot, J_tot)
  !=====================================================

    IMPLICIT NONE

    TYPE(element),                  INTENT(IN)  :: ele
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: uu
    REAL(KIND=8), DIMENSION(:),     INTENT(OUT) :: Phi_tot
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: J_tot
    !------------------------------------------------------

    INTEGER,                      POINTER :: N_quad
    INTEGER,                      POINTER :: N_points
    INTEGER,      DIMENSION(:),   POINTER :: loc
    REAL(KIND=8), DIMENSION(:,:), POINTER :: p
    REAL(KIND=8), DIMENSION(:,:), POINTER :: n
    REAL(KIND=8), DIMENSION(:),   POINTER :: w
    REAL(KIND=8), DIMENSION(:,:), POINTER :: xy
    !---------------------------------------------

    REAL(KIND=8), DIMENSION(N_eqn) :: Phi_b, Phi_S
    REAL(KIND=8), DIMENSION(N_eqn) :: uu_q, S_q

    REAL(KIND=8), DIMENSION(N_dim) :: a_q

    REAL(KIND=8), DIMENSION(N_eqn, N_dim) :: ff_q

    REAL(KIND=8), DIMENSION(N_eqn, N_eqn, N_dim) :: AA
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: AA_n

    INTEGER :: i, j, iq, k, l, id
    !---------------------------------------------

    !------------------
    ! Boundary integral
    !----------------------------------------------
    phi_b = 0.d0

    DO j = 1, ele%N_faces

       N_quad   => ele%faces(j)%f%N_quad
       N_points => ele%faces(j)%f%N_points
       loc      => ele%faces(j)%f%l_nu
       p        => ele%faces(j)%f%phi_q
       w        => ele%faces(j)%f%w_q
       n        => ele%faces(j)%f%n_q
       xy       => ele%faces(j)%f%xx_q

       DO iq = 1, N_quad

          uu_q = 0.d0
          DO k = 1, N_points
             uu_q = uu_q + uu(:, loc(k)) * p(k, iq)
          ENDDO

          ff_q = FOS_advection_flux(pb_type, uu_q, visc, xy(1,iq), xy(2,iq))

          DO id = 1, N_dim
             phi_b = phi_b + w(iq)*ff_q(:,id)*n(id, iq)
          ENDDO

       ENDDO

       NULLIFY(N_quad, N_points, loc, p, w, n, xy)

    ENDDO
    !----------------------------------------------

    !----------------
    ! Domain integral
    !----------------------------------------------
    Phi_S = 0.d0

    p  => ele%phi_q
    w  => ele%w_q
    xy => ele%xx_q

    DO iq = 1, ele%N_quad

       uu_q = 0.d0
       DO k = 1, ele%N_points
          uu_q = uu_q + uu(:, k) * p(k, iq)
       ENDDO

       S_q = FOS_source(pb_type, uu_q, visc, xy(1,iq), xy(2,iq))

       Phi_S = Phi_S + w(iq)*S_q

    ENDDO
    
    NULLIFY(p, w, xy)
    !----------------------------------------------
   
    Phi_tot = Phi_b - Phi_S
                  !^^^

    !-----------------
    ! Jacobian matrix
    !----------------------------------------------
    J_tot = 0.0

    DO l = 1, ele%N_points

       ! Boundary integral
       DO j = 1, ele%N_faces

          N_quad   => ele%faces(j)%f%N_quad
          N_points => ele%faces(j)%f%N_points
          loc      => ele%faces(j)%f%l_nu
          p        => ele%faces(j)%f%phi_q
          w        => ele%faces(j)%f%w_q
          n        => ele%faces(j)%f%n_q
          xy       => ele%faces(j)%f%xx_q

          DO i = 1, N_points

             IF( loc(i) == l ) THEN

                DO iq = 1, N_quad

                   uu_q = 0.d0
                   DO k = 1, N_points
                      uu_q = uu_q + uu(:, loc(k)) * p(k, iq)
                   ENDDO

                   a_q = advection_speed(pb_type, uu_q(1), xy(1,iq), xy(2,iq))
          
                   AA = FOS_Jacobian(a_q, visc)

                   AA_n = 0.d0
                   DO id = 1, ele%N_dim
                      AA_n = AA_n + AA(:,:, id)*n(id, iq)
                   ENDDO
                   
                   J_tot(:,:, l) = J_tot(:,:, l) + &
                                   w(iq) * AA_n * p(i, iq)

                ENDDO ! iq => face # N_quad

             ENDIF

          ENDDO ! i => ele # N_points

          NULLIFY(N_quad, N_points, loc, p, w, n, xy)

       ENDDO ! j => ele # N_faces

       ! Domain integral
       p  => ele%phi_q
       w  => ele%w_q
       xy => ele%xx_q

       DO iq = 1, ele%N_quad

          uu_q = 0.d0
          DO k = 1, ele%N_points
             uu_q = uu_q + uu(:, k) * p(k, iq)
          ENDDO

          J_tot(:,:, l) = J_tot(:,:, l) - &
                          w(iq) * p(l, iq) * &
                          FOS_J_Source(pb_type, uu_q, visc, xy(1,iq), xy(2,iq))

       ENDDO

       NULLIFY(p, w, xy)

    ENDDO ! l => ele # N_points

  END SUBROUTINE total_residual_imp
  !================================

  !=========================================
  FUNCTION Res_weak_bc_imp(uu) RESULT(res_b)
  !=========================================
  !
  ! Impose the boundary conditions as 
  ! inflow/outflow type. The inflow/outflow
  ! state is the exact solution: (u, p, q)
  !
    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: uu

    REAL(KIND=8), DIMENSION(SIZE(uu,1),SIZE(uu,2)) :: res_b
    !------------------------------------------------------

    TYPE(element) :: ele

    INTEGER,                      POINTER :: N_quad
    INTEGER,                      POINTER :: N_points
    INTEGER,      DIMENSION(:),   POINTER :: loc
    INTEGER,      DIMENSION(:),   POINTER :: nu
    REAL(KIND=8), DIMENSION(:,:), POINTER :: p
    REAL(KIND=8), DIMENSION(:,:), POINTER :: n
    REAL(KIND=8), DIMENSION(:),   POINTER :: w
    REAL(KIND=8), DIMENSION(:,:), POINTER :: xy
    !---------------------------------------------------

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: u

    REAL(KIND=8), DIMENSION(N_eqn)        :: lambda
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: RR, LL
    REAL(KIND=8), DIMENSION(N_eqn, N_eqn) :: Neg

    REAL(KIND=8), DIMENSION(N_dim) :: a_q
    REAL(KIND=8), DIMENSION(N_eqn) :: u_ex, u_q, eta, alpha
    !---------------------------------------------------

    REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: J_b
    REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: MM
    
    INTEGER, DIMENSION(N_eqn) :: idr

    INTEGER, DIMENSION(:), ALLOCATABLE :: idc, idm

    INTEGER :: je, jf, iq, k, j, i

    PetscErrorCode :: ierr
    !---------------------------------------------------

    res_b = 0.d0

    DO je = 1, N_elements

       ele = elements(je)%p

       ALLOCATE( u(N_eqn, ele%N_points) )

       u = uu(:, ele%Nu)

       DO jf = 1, ele%N_faces

          ! Boundary face
          IF( ele%faces(jf)%f%c_ele == 0 ) THEN

             N_quad   => ele%faces(jf)%f%N_quad
             N_points => ele%faces(jf)%f%N_points
             nu       => ele%faces(jf)%f%nu
             loc      => ele%faces(jf)%f%l_nu
             p        => ele%faces(jf)%f%phi_q
             w        => ele%faces(jf)%f%w_q
             n        => ele%faces(jf)%f%n_q
             xy       => ele%faces(jf)%f%xx_q

             ALLOCATE( J_b(N_eqn,    N_eqn, &
                           N_points, N_points) )
             
             J_b = 0.0
             
             ALLOCATE( idc(N_points*N_eqn), &
                       idm(N_eqn) )
             
             DO iq = 1, N_quad

                u_q = 0.d0
                DO k = 1, N_points
                   u_q = u_q + u(:, loc(k)) * p(k, iq)
                ENDDO

                a_q = advection_speed(pb_type, u_q(1), xy(1, iq), xy(2, iq))

                u_ex(1)   = exact_solution(pb_type, xy(:, iq), visc) ! u
                u_ex(2:3) = exact_grad(pb_type, xy(:, iq), visc)     ! p, q

                lambda = FOS_eigenvalues(a_q, visc, n(:,iq))
                RR = FOS_right_eigenvectors(a_q, visc, n(:,iq))
                LL = FOS_left_eigenvectors(a_q, visc, n(:,iq))

                Neg = 0.d0
                DO j = 1, N_eqn
                   Neg(j, j) = MIN(0.d0, lambda(j))
                ENDDO

                eta = MATMUL(LL, u_ex - u_q)
                
                alpha = MATMUL( RR, MATMUL(Neg, eta) )

                DO i = 1, N_points

                   res_b(:, Nu(i)) = res_b(:, Nu(i)) + &
                                     w(iq)*p(i, iq)*alpha

                   !----------------
                   ! Jacobian Matrix
                   !-----------------------------------------
                   DO k = 1, N_points

                      J_b(:,:, i, k) = J_b(:,:, i, k) -            &
                                       w(iq)*p(i, iq)*p(k, iq) *   &
                                       MATMUL( RR, MATMUL(Neg, LL) )

                   ENDDO

                ENDDO

             ENDDO ! iq -> # quad

             ! Assembly Jacobian
             ALLOCATE( MM(N_eqn, N_eqn*N_points) )

             DO i = 1, N_points

                idr = N_eqn*(ele%faces(jf)%f%Nu(i) - 1) + (/ (j, j = 1, N_eqn) /)
                idr = idr - 1
                
                DO k = 1, N_points
                   
                   idm = (/ (j, j = N_eqn*(k-1)+1, N_eqn*k) /)

                   idc(idm) = N_eqn*(ele%faces(jf)%f%Nu(k) - 1) + (/ (j, j = 1, N_eqn) /)
                   idc(idm) = idc(idm) - 1

                   MM(:, idm) = J_b(:,:, i, k)

                ENDDO
                
                DO j = 1, N_eqn

                   CALL MatSetValues(LHS, 1,              idr(j),  &
                                          N_eqn*N_points, idc,     &
                                          MM(j,:), ADD_VALUES, ierr)

                ENDDO

             ENDDO

             DEALLOCATE( J_b, MM, idc, idm )

          ENDIF

       ENDDO ! jf -> # faces

       DEALLOCATE( u )

    ENDDO ! je -> # ele

  END FUNCTION Res_weak_bc_imp
  !===========================

  !=========================================
  FUNCTION CFL_law(rhs, CFL_o) RESULT(CFL_n)
  !=========================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rhs
    REAL(KIND=8),                 INTENT(IN) :: CFL_o

    REAL(KIND=8) :: CFL_n
    !-----------------------------------------------
   
    REAL(KIND=8), PARAMETER :: CFL_max = 1.E+20
    !-----------------------------------------------

    ! L2 res u
    R_0 = SQRT(SUM(rhs(2, :)**2)) / REAL(N_dofs,8)
   
    IF( R_1 > R_0 ) THEN

       CFL_n = MIN(2.d0*CFL_o,  CFL_o * R_1/R_0)

    ELSE

       CFL_n = MAX(0.1d0*CFL_o,  CFL_o * R_1/R_0)

    ENDIF

    CFL_n = MIN(CFL_n, CFL_max)
       
    R_1 = R_0
    
  END FUNCTION CFL_law
  !===================
  
END MODULE space_integration
