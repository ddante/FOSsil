MODULE LLxFS_method

  USE Element_class
  USE geometry,       ONLY: N_dim, elements
  USE init_problem,   ONLY: pb_type
  USE models,         ONLY: advection_speed, source_term
  USE init_problem,   ONLY: is_visc, visc

  IMPLICIT NONE
  PRIVATE
  
  PUBLIC :: LLxFS_scheme

CONTAINS

  !==========================================================
  SUBROUTINE LLxFS_scheme(ele, Phi_tot, u, D_u, Phi_i, alpha)
  !==========================================================

    IMPLICIT NONE

    TYPE(element),                INTENT(IN)  :: ele
    REAL(KIND=8),                 INTENT(IN)  :: Phi_tot
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: u
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: D_u    
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: Phi_i
    REAL(KIND=8),                 INTENT(OUT) :: alpha
    !--------------------------------------------------

    INTEGER :: Ns, k, j, Nu

    REAL(KIND=8) :: theta_e, xi, diff_u, k_s
    
    REAL(KIND=8)                   :: u_m, x_m, y_m
    REAL(KIND=8), DIMENSION(N_dim) :: a_m

    REAL(KIND=8), DIMENSION(SIZE(Phi_i)) :: Stab_i
    !--------------------------------------------------

    Ns = ele%N_points

    !-----------
    ! Mean state
    !----------------------------------------------
    u_m = SUM(u) / REAL(Ns)

    x_m = SUM(ele%Coords(1, :)) / REAL(Ns)
    y_m = SUM(ele%Coords(2, :)) / REAL(Ns)

    a_m = advection_speed(pb_type, u_m, x_m, y_m)

    !--------------------
    ! Diffusive parameter
    !-------------------------------------------------
    alpha = 0.d0

    DO k = 1, Ns
       alpha = 0.5d0 * MAX( alpha, ABS( DOT_PRODUCT(a_m, ele%rd_n(:, k)) ) )
    ENDDO

    alpha = 3.d0*alpha

    IF ( is_visc ) THEN
       alpha = alpha + 3.d0*visc
    ENDIF
    
    !-----------
    ! LxF Scheme
    !----------------------------------------
    Phi_i = Phi_tot/REAL(Ns) + alpha*(u - u_m)

    !---------
    ! Limiting
    !------------------------
    Phi_i = limitation(Phi_i)

    !--------------
    ! Stabilization
    !---------------------------------
    Stab_i = GLS_stabilization(ele, u, D_u)

    !Stab_i = CIP_stabilization2(ele, u)

    Phi_i = Phi_i + Stab_i

  END SUBROUTINE LLxFS_scheme
  !==========================

  !=============================================
  FUNCTION limitation(Phi_i__L) RESULT(Phi_i__H)
  !=============================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: Phi_i__L
 
    REAL(KIND=8), DIMENSION(SIZE(Phi_i__L)) :: Phi_i__H
    !---------------------------------------------------

    REAL(KIND=8), DIMENSION(SIZE(Phi_i__L)) :: beta

    REAL(KIND=8) :: Phi_tot

    REAL(KIND=8), PARAMETER :: eps = 1.0E-15
    !---------------------------------------------------

    Phi_tot = SUM(Phi_i__L)

    IF( ABS(Phi_tot) > eps ) THEN

       beta = Phi_i__L / Phi_tot 

       beta = MAX(0.d0, beta) + eps

       beta = beta / SUM(beta)

    ELSE

       beta = 0.d0

    ENDIF

    Phi_i__H = beta * Phi_tot
    
  END FUNCTION limitation
  !======================

  !===================================================
  FUNCTION GLS_stabilization(ele, u, D_u) RESULT(Stab)
  !===================================================

    IMPLICIT NONE

    TYPE(element),                INTENT(IN)  :: ele
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: u
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: D_u

    REAL(KIND=8), DIMENSION(SIZE(u)) :: Stab
    !---------------------------------------------

    REAL(KIND=8) :: x_i, y_i
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: a
    REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: SS
    
    REAL(KIND=8), DIMENSION(N_dim) :: D_u_q, &
                                      Dr_u_q, a_q

    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: D_muD_phi_q

    REAL(KIND=8) :: s_q

    REAL(KIND=8), DIMENSION(:),     POINTER :: w
    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: D_phi_q
    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: D_phi_k
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: phi_q

    REAL(KIND=8) :: D_muD_u_q,  xi_Re, Re_h
    REAL(KIND=8) :: tau, Stab_D_u, Stab_D_phi
    
    INTEGER :: iq, i, k, id, N_s, N_v 
    !---------------------------------------------

    Stab = 0.d0

    N_s = ele%N_points
    N_v = ele%N_verts

    ALLOCATE( a(N_dim, N_s), SS(N_s) )
    ALLOCATE( D_muD_phi_q(N_s) )
    
          w => ele%w_q
    D_phi_q => ele%D_phi_q
    D_phi_k => ele%D_phi_k
      phi_q => ele%phi_q

    SS = 0.d0

    DO i = 1, N_s
       
       x_i = ele%coords(1, i)
       y_i = ele%coords(2, i)

       a(:, i) = advection_speed(pb_type, u(i), x_i, y_i)
       SS(i)   = source_term( pb_type, u(i), visc, x_i, y_i) 
       
    ENDDO

    !------------------
    ! Scaling parameter
    !---------------------------------------------
    tau = 0.d0
    
    DO i = 1, N_v

       tau = tau + &
             0.5d0 * ABS( DOT_PRODUCT(a(:, i), ele%rd_n(:, i)) ) / ele%volume
    
    ENDDO

    tau = tau + ( REAL(N_v) * visc / ele%volume )
          
    tau = 1.d0 / tau

!    Re_h = local_Pe(ele, u)
!    xi_Re = MAX(0.d0, 1.d0 - 1.d0/Re_h)
!    tau = tau * xi_Re

    !-------------------
    ! Stabilization term
    !--------------------------------------------------
    DO iq = 1, ele%N_quad

       !------------------------
       ! Div( mu * Grad(Phi_i) )
       !-----------------------------------------------
       D_muD_phi_q = 0.d0       
       DO i = 1, N_s

          DO id = 1, N_dim

             DO k = 1, N_s
                D_muD_phi_q(i) = D_muD_phi_q(i) +     &
                                 visc*D_phi_k(id, i, k) * &
                                 D_phi_q(id, k, iq)
             ENDDO

          ENDDO

       ENDDO     
       !-------------------------------------------------
      
      Dr_u_q = 0.d0 ! Reconstructed gradient
       D_u_q = 0.d0 ! Continuous gradient
         a_q = 0.d0
         s_q = 0.d0

      D_muD_u_q = 0.d0
         
       DO i = 1, N_s

          Dr_u_q = Dr_u_q + phi_q(i, iq) * D_u(:, i)
           D_u_q = D_u_q  + D_phi_q(:, i, iq) * u(i)
             a_q = a_q    + phi_q(i, iq) * a(:, i)
             s_q = s_q    + phi_q(i, iq) * SS(i)

          D_muD_u_q = D_muD_u_q + D_muD_phi_q(i) * u(i)

       ENDDO

       Stab_D_u = DOT_PRODUCT(a_q, D_u_q) - s_q
  
       Stab_D_u = ( Stab_D_u - D_muD_u_q ) * tau

       DO i = 1, N_s
 
          Stab_D_phi = DOT_PRODUCT(a_q, D_phi_q(:, i, iq))

          Stab_D_phi = Stab_D_phi - D_muD_phi_q(i)
                   
          Stab(i) = Stab(i) + Stab_D_u * Stab_D_phi * w(iq)
          
         ! Lax-Wendroff like stabilization
         Stab(i) = Stab(i) + &
                   visc * DOT_PRODUCT(D_phi_q(:, i, iq), (D_u_q - Dr_u_q)) * w(iq) 

       ENDDO

   ENDDO

   NULLIFY( w, D_phi_q, phi_q, D_phi_k )
   DEALLOCATE( a, SS, D_muD_phi_q )

  END FUNCTION GLS_stabilization
  !=============================
  
  !==============================================
  FUNCTION CIP_stabilization(ele, u) RESULT(Stab)
  !==============================================

    IMPLICIT NONE

    TYPE(element),              INTENT(IN)  :: ele
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: u

    REAL(KIND=8), DIMENSION(SIZE(u)) :: Stab
    !---------------------------------------------   
    
    INTEGER,                        POINTER :: N_quad
    INTEGER,                        POINTER :: N_points_f
    INTEGER,      DIMENSION(:),     POINTER :: loc
    REAL(KIND=8),                   POINTER :: h_f
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: n_1
    REAL(KIND=8), DIMENSION(:),     POINTER :: w
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: p
    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: p_Dphi_1_q
    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: p_Dphi_2_q
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: p_Du_1_q
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: p_Du_2_q
    INTEGER,      DIMENSION(:),     POINTER :: loc_con

    REAL(KIND=8), DIMENSION(N_dim) :: n_2, tt
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: a_q

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: a    
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: p_Dphi_l_2_q 
    REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: Jump_grad_phi
    REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: Jump_grad_u

    REAL(KIND=8) :: g_1, b_n, b_t, x_i, y_i

    INTEGER :: i, if, iq, k
    !---------------------------------------------

    Stab = 0.d0

    g_1 = 0.01

    ALLOCATE ( a(N_dim, ele%N_points) )

    DO i = 1, ele%N_points

       x_i = ele%Coords(1, i)
       y_i = ele%Coords(2, i)

       a(:, i) = advection_speed(pb_type, u(i), x_i, y_i)

    ENDDO
    
    ! Construct the stabilization term for each node
    DO i = 1, SIZE(Stab)

       ! Split the boundary integral on the element
       ! as the sum of integral on each face
       DO if = 1, ele%N_faces

          N_quad     => ele%faces(if)%f%N_quad
          N_points_f => ele%faces(if)%f%N_points
          loc        => ele%faces(if)%f%l_nu
          w          => ele%faces(if)%f%w_q
          h_f        => ele%faces(if)%f%area
          n_1        => ele%faces(if)%f%n_q
          p          => ele%faces(if)%f%phi_q
          p_Dphi_1_q => ele%faces(if)%f%p_Dphi_1_q
          p_Dphi_2_q => ele%faces(if)%f%p_Dphi_2_q
          p_Du_1_q   => ele%faces(if)%f%p_Du_1_q
          p_Du_2_q   => ele%faces(if)%f%p_Du_2_q
          loc_con    => ele%faces(if)%f%loc_con

          ALLOCATE( Jump_grad_phi(N_quad), &
                    Jump_grad_u(N_quad)    )

          ALLOCATE( p_Dphi_l_2_q(N_dim, N_quad), &
                             a_q(N_dim, N_quad) )

          IF(loc_con(i) /= 0) THEN
             p_Dphi_l_2_q = p_Dphi_2_q(:, loc_con(i), :)             
          ELSE
             p_Dphi_l_2_q = 0.d0             
          ENDIF

          !------------------------
          ! Stabilization Parameter
          !----------------------------------------------
          b_n = 0.d0; b_t = 0.d0
          
          DO iq = 1, N_quad
             
             a_q(:,iq) = 0.0             
             DO k = 1, N_points_f                
                a_q(:, iq) = a_q(:, iq) + a(:, loc(k)) * p(k, iq)                
             ENDDO
             
             !b_n = MAX( b_n, ABS(DOT_PRODUCT(a_q(:,iq), n_1(:, iq))) )             
             !b_t = MAX( b_t, ABS(-a_q(1, iq)*n_1(2,iq) + a_q(2,iq)*n_1(1,iq)) )
                          
          ENDDO
          

          ! Perform the integration on the single face
          DO iq = 1, N_quad

            n_2 = -n_1(:, iq)
          
             ! Jump of the gradient of the i-th shape function
             ! at the quadrature points
             Jump_grad_phi(iq) = DOT_PRODUCT(p_Dphi_1_q(:, i,iq), n_1(:, iq) ) + &
                                 DOT_PRODUCT(p_Dphi_l_2_q(:, iq), n_2 )

             ! Jump of the gradient of the sulution
             ! at the quadrature points
             Jump_grad_u(iq) = DOT_PRODUCT(p_Du_1_q(:, iq), n_1(:, iq) ) + &
                               DOT_PRODUCT(p_Du_2_q(:, iq), n_2)

             Stab(i) = Stab(i) + &
                       g_1 *     &
                       ( DOT_PRODUCT(a_q(:,iq),n_1(:,iq))**2 ) * &
                       h_f*h_f * &
                       w(iq) * Jump_grad_phi(iq) * Jump_grad_u(iq)

          ENDDO

          NULLIFY( N_quad, N_points_f, loc,    &
                   w, h_f, n_1, p,             &
                   p_Dphi_1_q, p_Dphi_2_q,     &
                   p_Du_1_q, p_Du_2_q, loc_con )

          DEALLOCATE( Jump_grad_phi, Jump_grad_u, &
                      p_Dphi_l_2_q, a_q   )

       ENDDO

    ENDDO

    DEALLOCATE( a )

  END FUNCTION CIP_stabilization
  !=============================
    


  !==============================================
  FUNCTION CIP_stabilization2(ele, u) RESULT(Stab)
  !==============================================

    IMPLICIT NONE

    TYPE(element),              INTENT(IN)  :: ele
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: u

    REAL(KIND=8), DIMENSION(SIZE(u)) :: Stab
    !---------------------------------------------   
    
    INTEGER,                        POINTER :: N_quad
    INTEGER,                        POINTER :: N_points_f
    INTEGER,      DIMENSION(:),     POINTER :: loc
    REAL(KIND=8),                   POINTER :: h_f
    REAL(KIND=8), DIMENSION(:),     POINTER :: w
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: n_1
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: p
    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: p_Dphi_1_q
    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: p_Dphi_2_q
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: p_Du_1_q
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: p_Du_2_q
    INTEGER,      DIMENSION(:),     POINTER :: loc_con

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: a_q
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: a
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: p_Dphi_l_2_q 

    REAL(KIND=8) :: Jump_D_phi, Jump_D_u
    REAL(KIND=8) :: b_n, b_t, g_1, x_i, y_i, mod_a_q

    INTEGER :: i, if, iq, k
    !---------------------------------------------

    Stab = 0.d0

    g_1 = 0.008

    ALLOCATE ( a(N_dim, ele%N_points) )
    DO i = 1, ele%N_points

       x_i = ele%Coords(1, i)
       y_i = ele%Coords(2, i)

       a(:, i) = advection_speed(pb_type, u(i), x_i, y_i)

    ENDDO
   
    DO i = 1, SIZE(Stab)

       ! Split the boundary integral on the element
       ! as the sum of integral on each face
       DO if = 1, ele%N_faces

          N_quad     => ele%faces(if)%f%N_quad      ! ok
          N_points_f => ele%faces(if)%f%N_points    ! ok
          loc        => ele%faces(if)%f%l_nu        ! ok
          n_1        => ele%faces(if)%f%n_q         ! ok
          w          => ele%faces(if)%f%w_q         ! ok
          h_f        => ele%faces(if)%f%area        ! ok
          p          => ele%faces(if)%f%phi_q       ! ok
          p_Dphi_1_q => ele%faces(if)%f%p_Dphi_1_q  ! ok
          p_Dphi_2_q => ele%faces(if)%f%p_Dphi_2_q  ! ok
          p_Du_1_q   => ele%faces(if)%f%p_Du_1_q    ! ok
          p_Du_2_q   => ele%faces(if)%f%p_Du_2_q    ! ok
          loc_con    => ele%faces(if)%f%loc_con     ! ok

          ALLOCATE( p_Dphi_l_2_q(N_dim, N_quad), & 
                             a_q(N_dim, N_quad) )

          IF(loc_con(i) /= 0) THEN
             p_Dphi_l_2_q = p_Dphi_2_q(:, loc_con(i), :)
          ELSE
             p_Dphi_l_2_q = 0.d0
          ENDIF

          DO iq = 1, N_quad

             a_q(:, iq) = 0.d0
             DO k = 1, N_points_f
                a_q(:, iq) = a_q(:, iq) + a(:, loc(k)) * p(k, iq) 
             ENDDO              
             
          ENDDO
          

          DO iq = 1, N_quad

             !---------------
             ! Gradient Jump
             !-----------------------------------------------------------
             Jump_D_phi = DOT_PRODUCT( p_Dphi_1_q(:, i, iq), n_1(:, iq) ) + &
                          DOT_PRODUCT( p_Dphi_l_2_q(:, iq), -n_1(:, iq) )
             
             Jump_D_u   = DOT_PRODUCT( p_Du_1_q(:, iq),  n_1(:, iq) ) + &
                          DOT_PRODUCT( p_Du_2_q(:, iq), -n_1(:, iq) )

             b_n = ABS( DOT_PRODUCT(a_q(:, iq), n_1(:,iq)) )
             
             !--------------------
             ! Stabilization term
             !-----------------------------------------------------------
             Stab(i) = Stab(i) + g_1 * b_n * &
                       w(iq) * (h_f**2) * (Jump_D_phi * Jump_D_u)

          ENDDO

          !--------------------
          ! Destruc face stuff
          !----------------------------------------
          DEALLOCATE( p_Dphi_l_2_q, a_q )

          NULLIFY( N_quad, N_points_f, loc, w, n_1, &
                   h_f, p, p_Dphi_1_q, p_Dphi_2_q,  &
                   p_Du_1_q, p_Du_2_q, loc_con      )


       ENDDO

    ENDDO

    DEALLOCATE ( a )

  END FUNCTION CIP_stabilization2
  !=============================


  !=====================================
  FUNCTION Local_Pe(ele, u) RESULT(l_Pe)
  !=====================================

    IMPLICIT NONE

    TYPE(element),              INTENT(IN) :: ele
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u

    REAL(KIND=8) :: l_Pe
    !---------------------------------------------

    REAL(KIND=8) :: x_i, y_i, h, p

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

!!$    h =  ele%Volume / SQRT(N_v * h)

    P = 0.d0
    DO i = 1, ele%N_faces
       P = P + SUM(ele%faces(i)%f%w_q)
    ENDDO
    
    h = ele%Volume/P

    h = h*(2.d0*(2.d0 + SQRT(2.d0)))

    l_Pe = SQRT(SUM(a*a)) * h / visc

  END FUNCTION Local_Pe
  !====================

END MODULE LLxFS_method
