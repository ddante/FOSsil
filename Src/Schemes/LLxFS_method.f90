MODULE LLxFS_method

  USE Element_class
  USE geometry,       ONLY: N_dim, elements
  USE init_problem,   ONLY: pb_type
  USE models,         ONLY: advection_speed
  USE init_problem,   ONLY: N_eqn, pb_type, visc
  USE Lin_Algebra,    ONLY: inverse
  USE FOS_system

  IMPLICIT NONE
  PRIVATE
  
  PUBLIC :: LLxFS_scheme

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

END MODULE LLxFS_method
