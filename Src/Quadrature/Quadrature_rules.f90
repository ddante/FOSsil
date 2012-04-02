MODULE Quadrature_rules

  USE Element_Class
  
  IMPLICIT NONE

CONTAINS

  !====================================
  FUNCTION oInt_n(ele, ff) RESULT(i_ff)
  !====================================
  !
  ! \oint_{\partial E} {ff * n}
  !
    IMPLICIT NONE

    TYPE(element),                INTENT(IN) :: ele
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ff

    REAL(KIND=8) :: i_ff
    !----------------------------------------------

    INTEGER,                      POINTER :: N_quad
    INTEGER,                      POINTER :: N_points
    INTEGER,      DIMENSION(:),   POINTER :: loc
    REAL(KIND=8), DIMENSION(:,:), POINTER :: p
    REAL(KIND=8), DIMENSION(:,:), POINTER :: n
    REAL(KIND=8), DIMENSION(:),   POINTER :: w
    
    
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ff_q

    INTEGER :: j, iq, k
    !----------------------------------------------

    i_ff = 0.d0

    ALLOCATE( ff_q(SIZE(ff,1)) )

    DO j = 1, ele%N_faces

       N_quad   => ele%faces(j)%f%N_quad
       N_points => ele%faces(j)%f%N_points
       loc      => ele%faces(j)%f%l_nu
       p        => ele%faces(j)%f%phi_q
       w        => ele%faces(j)%f%w_q
       n        => ele%faces(j)%f%n_q
       
       DO iq = 1, N_quad

          ff_q = 0.d0
          DO k = 1, N_points 

             ff_q = ff_q + ff(:, loc(k)) * p(k, iq)

          ENDDO 

          i_ff = i_ff + w(iq) * DOT_PRODUCT(ff_q, n(:, iq) )

       ENDDO
       
       NULLIFY(N_quad, N_points, loc, p, w, n)
       
    ENDDO

    DEALLOCATE( ff_q )

  END FUNCTION oInt_n
  !==================  

  !===================================
  FUNCTION Int_d(ele, ff) RESULT(i_ff)
  !===================================
  !
  ! \int_{E} {ff}
  !
    IMPLICIT NONE

    TYPE(element),              INTENT(IN) :: ele
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: ff

    REAL(KIND=8) :: i_ff
    !----------------------------------------------

    INTEGER :: N_quad,  N_points
   
    REAL(KIND=8), DIMENSION(:,:), POINTER :: p
    REAL(KIND=8), DIMENSION(:),   POINTER :: w
        
    REAL(KIND=8) :: ff_q

    INTEGER :: iq
    !----------------------------------------------

    i_ff = 0.d0

    N_quad   =  ele%N_quad
    N_points =  ele%N_points
    p        => ele%phi_q
    w        => ele%w_q

    DO iq = 1, N_quad

       ff_q = SUM( ff * p(:, iq) )

       i_ff = i_ff + w(iq) * ff_q

    ENDDO
       
    NULLIFY( p, w )
       
  END FUNCTION Int_d
  !=================

  !==================================
  FUNCTION int_d_G(ele, k) RESULT(ig)
  !==================================
  !
  ! int_{E} {Grad phi_k}
  !
    IMPLICIT NONE

    TYPE(element), INTENT(IN) :: ele
    INTEGER,       INTENT(IN) :: k

    REAL, DIMENSION(ele%N_dim) :: ig
    !---------------------------------

    INTEGER :: N_quad

    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: d_phi
    REAL(KIND=8), DIMENSION(:),     POINTER :: w

    INTEGER :: iq
    !-----------------------------------------

    ig = 0.d0

    N_quad =  ele%N_quad
    w      => ele%w_q
    d_phi  => ele%D_phi_q

    DO iq = 1, N_quad

       ig = ig + d_phi(:, k, iq)*w(iq)

    ENDDO
    
    NULLIFY( w, d_phi )

  END FUNCTION int_d_G
  !===================  
  
  !=======================================
  FUNCTION int_d_Mij(ele, i, j) RESULT(ig)
  !=======================================
  !
  ! int_{E} {phi_i phi_j}
  !
    IMPLICIT NONE

    TYPE(element), INTENT(IN) :: ele
    INTEGER,       INTENT(IN) :: i, j

    REAL :: ig
    !---------------------------------

    INTEGER :: N_quad

    REAL(KIND=8), DIMENSION(:,:), POINTER :: phi
    REAL(KIND=8), DIMENSION(:),   POINTER :: w

    INTEGER :: iq
    !-----------------------------------------

    ig = 0.d0

    N_quad =  ele%N_quad
    w      => ele%w_q
    phi    => ele%phi_q

    DO iq = 1, N_quad

       ig = ig + &
            phi(i, iq)*phi(j, iq) * w(iq)

    ENDDO
    
    NULLIFY( w, phi )

  END FUNCTION int_d_Mij
  !=====================

  !========================================
  FUNCTION int_d_Gu_i(ele, u, i) RESULT(ig)
  !========================================
  !
  ! int_{E} {phi_i * Grad(u)}
  !
    IMPLICIT NONE

    TYPE(element),              INTENT(IN) :: ele
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u
    INTEGER,                    INTENT(IN) :: i

    REAL, DIMENSION(ele%N_dim) :: ig
    !---------------------------------

    INTEGER :: N_quad

    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: d_phi
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: phi
    REAL(KIND=8), DIMENSION(:),     POINTER :: w

    REAL, DIMENSION(ele%N_dim) :: Du_q

    INTEGER :: iq, id
    !-----------------------------------------

    ig = 0.d0

    N_quad =  ele%N_quad
    w      => ele%w_q
    phi    => ele%phi_q
    d_phi  => ele%D_phi_q

    DO iq = 1, N_quad

       Du_q = 0.d0
       DO id = 1, ele%N_dim
          Du_q(id) = Du_q(id) + SUM( d_phi(id, :, iq)*u )
       ENDDO

       ig = ig + phi(i, iq)*Du_q*w(iq)

    ENDDO
    
    NULLIFY( w, phi, d_phi )

  END FUNCTION int_d_Gu_i
  !======================


END MODULE Quadrature_rules
