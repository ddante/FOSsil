MODULE space_integration

  USE Element_Class

  USE geometry,       ONLY: N_dim, N_dofs, N_elements, elements

  USE init_problem,   ONLY: order, pb_type, is_visc, visc, &
                            with_source, CFL

  USE models,         ONLY: advection_flux, diffusion_flux, &
                            advection_speed, strong_bc, source_term

  USE Num_scheme

  USE Gradient_Reconstruction

  USE Quadrature_rules, ONLY: oInt_n, Int_d

use test_gradient

  IMPLICIT NONE

  !================================================
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: D_uu
  !================================================

  PRIVATE  
  PUBLIC :: compute_rhs, Local_Pe, ClearGradients

CONTAINS

  !====================================
  SUBROUTINE compute_rhs(uu, rhs, Dt_V)
  !====================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: uu
    REAL(KIND=8), DIMENSION(:), INTENT(OUT)   :: rhs
    REAL(KIND=8), DIMENSION(:), INTENT(OUT)   :: Dt_V
    !-------------------------------------------------

    TYPE(element) :: loc_ele, adh_ele

    INTEGER,      DIMENSION(:),   ALLOCATABLE :: Nu
    REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: u1, Phi_i
    REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: u2
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: D_u1

    REAL(KIND=8) :: Phi_tot, inv_dt

    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: p_Dphi_1_q
    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: p_Dphi_2_q

    REAL(KIND=8), DIMENSION(N_dim) :: p_Du_1_q, p_Du_2_q
    
    INTEGER :: Ns, je, i_f, k, iq, n_ele, istat
    !----------------------------------------------

    rhs = 0.d0

    Dt_v = 0.d0

    IF (.NOT. ALLOCATED(D_uu) ) THEN
       ALLOCATE( D_uu(N_dim, SIZE(uu)) )
    ENDIF
    
    D_uu = Compute_gradient(uu)

    DO je = 1, N_elements

       loc_ele = elements(je)%p

       Ns = loc_ele%N_points

       ALLOCATE( Nu(Ns), u1(Ns), Phi_i(Ns), &
                 D_u1(N_dim, Ns)            )

       Nu = loc_ele%NU

       u1 = uu(Nu);  D_u1 = D_uu(:, Nu)

       !-------------------------------
       ! Compute the total fluctuation
       !------------------------------------------------------
       IF( is_visc ) THEN
          Phi_tot = total_residual(loc_ele, u1, D_u1)
       ELSE          
          Phi_tot = total_residual(loc_ele, u1)
       ENDIF
       
       !--------------------------
       ! Update the face gradients
       !------------------------------------------------------
       DO i_f = 1, loc_ele%N_faces

          p_Dphi_1_q => loc_ele%faces(i_f)%f%p_Dphi_1_q
          p_Dphi_2_q => loc_ele%faces(i_f)%f%p_Dphi_2_q
          
          n_ele = loc_ele%faces(i_f)%f%c_ele

          IF( n_ele /= 0 ) THEN ! * Internal face *

             adh_ele = elements(n_ele)%p

             ALLOCATE( u2(adh_ele%N_points) )
             
             u2 = uu(adh_ele%NU)

             DO iq = 1, loc_ele%faces(i_f)%f%N_quad

                p_Du_1_q = 0.d0
                DO k = 1, loc_ele%N_points
                   p_Du_1_q = p_Du_1_q + p_Dphi_1_q(:, k,iq)*u1(k)
                ENDDO

                p_Du_2_q = 0.d0
                DO k = 1, adh_ele%N_points
                   p_Du_2_q = p_Du_2_q + p_Dphi_2_q(:, k,iq)*u2(k)
                ENDDO

                loc_ele%faces(i_f)%f%p_Du_1_q(:,iq) = p_Du_1_q
                loc_ele%faces(i_f)%f%p_Du_2_q(:,iq) = p_Du_2_q

             ENDDO

             DEALLOCATE( u2 )

             NULLIFY( p_Dphi_1_q, p_Dphi_2_q )

          ELSE ! * Boundary Face *

             loc_ele%faces(i_f)%f%p_Du_1_q = 0.d0
             loc_ele%faces(i_f)%f%p_Du_2_q = 0.d0

!!$             DO iq = 1, loc_ele%faces(i_f)%f%N_quad                !*
!!$                                                                   !*
!!$                p_Du_1_q = 0.d0                                    !*
!!$                DO k = 1, loc_ele%N_points                         !*
!!$                   p_Du_1_q = p_Du_1_q + p_Dphi_1_q(:, k,iq)*u1(k) !*
!!$                ENDDO                                              !*
!!$                                                                   !*
!!$                loc_ele%faces(i_f)%f%p_Du_1_q(:,iq) = p_Du_1_q     !*
!!$                loc_ele%faces(i_f)%f%p_Du_2_q(:,iq) = 0.d0         !*
!!$                                                                   !*
!!$             ENDDO                                                 !*

          ENDIF

       ENDDO
       !------------------------------------------------------

       !---------------------------
       ! Distribute the fluctuation
       !------------------------------------------------------------
       CALL distribute_residual(loc_ele, Phi_tot, u1, D_u1, Phi_i, inv_dt)

       ! Gather the nodal residual
       rhs(Nu) = rhs(Nu) + Phi_i
       
       Dt_V(Nu) = Dt_V(Nu) + inv_dt
       
       DEALLOCATE( Nu, u1, D_u1, Phi_i )

    ENDDO

    Dt_V = CFL/DT_V

    !---------------------
    ! Impose strong bc 
    ! RHS = 0 on the inlet
    !-------------------------------------
    CALL strong_bc(pb_type, visc, uu, rhs)

  END SUBROUTINE compute_rhs
  !=========================

  !===================================================
  FUNCTION total_residual(ele, u, D_u) RESULT(Phi_tot)
  !===================================================

    IMPLICIT NONE

    TYPE(element),                INTENT(IN) :: ele
    REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: u
    REAL(KIND=8), DIMENSION(:,:), OPTIONAL, &
                                  INTENT(IN) :: D_u
    
    REAL(KIND=8) :: Phi_tot
    !---------------------------------------------

    REAL(KIND=8) :: x, y
    REAL(KIND=8) :: Phi_b, Source

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: ff_a
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: ff_v
    REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: SS

    INTEGER :: Ns, i, j
    !---------------------------------------------

    Phi_b  = 0.d0
    Source = 0.d0

    Ns = ele%N_points

    ALLOCATE( ff_a(N_dim, Ns), ff_v(N_dim, Ns), &
              SS(Ns)                            )

    ff_a = 0.d0; ff_v = 0.d0; SS = 0.d0

    DO i = 1, Ns

       x = ele%Coords(1, i)
       y = ele%Coords(2, i)

       ff_a(:, i) = advection_flux(pb_type, u(i), x, y)

       SS(i) = source_term(pb_type, u(i), visc, x, y)

       IF( PRESENT(D_u) ) THEN
          ff_v(:, i) = diffusion_flux(pb_type, visc, D_u(:,i), x, y)
       ENDIF
           
    ENDDO

    Phi_b = oInt_n(ele, ff_a - ff_v)

    IF( with_source ) THEN
       Source = Int_d(ele, SS)
    ENDIF
    
    DEALLOCATE( ff_a, ff_v, SS )

    Phi_tot = Phi_b - Source
                  !^^^ 

  END FUNCTION total_residual
  !==========================

!!$  !=================================================
!!$  FUNCTION Penalization_flux(ele, u) RESULT(Phi_pen)
!!$  !=================================================
!!$
!!$    IMPLICIT NONE
!!$    TYPE(element),                INTENT(IN) :: ele
!!$    REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: u
!!$    
!!$    REAL(KIND=8) :: Phi_pen
!!$    !-----------------------------------------------
!!$
!!$    INTEGER,                      POINTER :: N_quad
!!$    REAL(KIND=8), DIMENSION(:,:), POINTER :: n_1
!!$    REAL(KIND=8), DIMENSION(:),   POINTER :: w
!!$    REAL(KIND=8), DIMENSION(:,:), POINTER :: p_Du_1_q
!!$    REAL(KIND=8), DIMENSION(:,:), POINTER :: p_Du_2_q
!!$    INTEGER,      DIMENSION(:),   POINTER :: loc_con
!!$    !-----------------------------------------------
!!$
!!$    REAL(KIND=8) :: Jump_D_u, C_ip
!!$
!!$    INTEGER :: jf, iq
!!$    !-----------------------------------------------
!!$
!!$    C_ip = 0.05d0
!!$    
!!$    Phi_pen = 0.d0
!!$
!!$    DO jf = 1, ele%N_faces
!!$
!!$       N_quad     => ele%faces(jf)%f%N_quad      ! ok      
!!$       n_1        => ele%faces(jf)%f%n_q         ! ok
!!$       w          => ele%faces(jf)%f%w_q         ! ok      
!!$       p_Du_1_q   => ele%faces(jf)%f%p_Du_1_q    ! ok
!!$       p_Du_2_q   => ele%faces(jf)%f%p_Du_2_q    ! ok
!!$       loc_con    => ele%faces(jf)%f%loc_con     ! ok
!!$
!!$       DO iq = 1, N_quad
!!$
!!$          !---------------
!!$          ! Gradient Jump
!!$          !-----------------------------------------------------------
!!$          Jump_D_u = DOT_PRODUCT( p_Du_1_q(:, iq),  n_1(:, iq) ) + &
!!$                     DOT_PRODUCT( p_Du_2_q(:, iq), -n_1(:, iq) )
!!$
!!$          Phi_pen = Phi_pen + w(iq) * visc * Jump_D_u 
!!$
!!$       ENDDO
!!$
!!$     ENDDO
!!$
!!$    Phi_pen = C_ip * Phi_pen
!!$
!!$  END FUNCTION Penalization_flux
!!$  !=============================
  
  !=====================================
  FUNCTION Local_Pe(ele, u) RESULT(l_Pe)
  !=====================================

    IMPLICIT NONE

    TYPE(element),              INTENT(IN) :: ele
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u

    REAL(KIND=8) :: l_Pe
    !---------------------------------------------

    REAL(KIND=8) :: x_i, y_i, h, P

    REAL(KIND=8), DIMENSION(N_dim) :: a, gg, xx

    INTEGER :: N_v, i, j

    REAL(KIND=8), PARAMETER ::  Pi = DACOS(-1.d0)
    !-------------------------------------------------

    N_v = ele%N_verts

    h = 0.d0

    DO j = 1, N_dim
       gg(j) = SUM(ele%coords(j, 1:N_v)) / REAL(N_v)
    ENDDO    

    DO i = 1, N_v

       x_i = ele%coords(1, i)
       y_i = ele%coords(2, i)

       xx = (/ x_i, y_i /)

       h = h + SUM( (xx - gg)**2 )

       a = advection_speed(pb_type, u(i), x_i, y_i)
       
    ENDDO

!!$    h = 2.d0*ele%Volume / SQRT(N_v * h)

    P = 0.d0
    DO i = 1, ele%N_faces
       P = P + SUM(ele%faces(i)%f%w_q)
    ENDDO
    
    h = ele%Volume/P

    l_Pe = SQRT(SUM(a*a)) * h / visc

  END FUNCTION Local_Pe
  !====================
  
  !==========================
  SUBROUTINE ClearGradients()
  !==========================

    IMPLICIT NONE

    IF( ALLOCATED(D_uu) ) DEALLOCATE( D_uu )

    CALL DestroyPETSc()

  END SUBROUTINE ClearGradients
  !============================
  

END MODULE space_integration
