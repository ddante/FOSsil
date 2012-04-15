MODULE space_integration

  USE Element_Class

  USE geometry,       ONLY: N_dim, N_dofs, N_elements, elements

  USE init_problem,   ONLY: order, pb_type, bc_type, &
                            visc, CFL, N_eqn

  USE models,         ONLY: advection_flux, advection_speed, &
                            strong_bc, exact_solution, exact_grad
  USE FOS_system
  USE Num_scheme

  IMPLICIT NONE

  PRIVATE  

  !==========================================
  INTEGER, PARAMETER :: STRONG_BC_TYPE = 0, &
                          WEAK_BC_TYPE = 1
  !==========================================

  PUBLIC :: compute_rhs

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


END MODULE space_integration
