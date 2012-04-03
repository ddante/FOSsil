MODULE space_integration

  USE Element_Class

  USE geometry,       ONLY: N_dim, N_dofs, N_elements, elements

  USE init_problem,   ONLY: order, pb_type, visc, CFL, N_eqn

  USE models,         ONLY: advection_flux, diffusion_flux, &
                            advection_speed, strong_bc
  USE FOS_system
  USE Num_scheme

  IMPLICIT NONE

  PRIVATE  
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

    !---------------------
    ! Impose strong bc 
    ! RHS = 0 on the inlet
    !-------------------------------------
    CALL strong_bc(pb_type, visc, uu, rhs)

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

END MODULE space_integration
