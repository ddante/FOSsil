MODULE FOS_system
  
  USE geometry, ONLY: N_dim
  USE models,   ONLY: advection_speed, source_term
  USE Lin_Algebra

  IMPLICIT NONE

  PRIVATE

  REAL(KIND=8), PARAMETER :: PI = DACOS(-1.d0)

  PUBLIC :: FOS_advection_flux, FOS_source, &
            FOS_eigenvalues, FOS_right_eigenvectors, &
            FOS_left_eigenvectors, FOS_Jacobian

CONTAINS

  !============================================================
  FUNCTION FOS_advection_flux(type_pb, uu, nu, x, y) RESULT(ff)
  !============================================================

    IMPLICIT NONE
    
    INTEGER,                    INTENT(IN) :: type_pb
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: uu
    REAL(KIND=8),               INTENT(IN) :: nu
    REAL(KIND=8),               INTENT(IN) :: x, y
    
    REAL(KIND=8), DIMENSION(SIZE(uu), N_dim) :: ff
    !---------------------------------------------

    REAL(KIND=8), DIMENSION(N_dim) :: a

    REAL(KIND=8) :: u, p, q, TR, LR, mod_a
    REAL(KIND=8) :: b, c
    !---------------------------------------------

    u = uu(1); p = uu(2); q = uu(3)

    b = 0.d0; c = 0.d0

    a = advection_speed(type_pb, u, x, y)
    mod_a = DSQRT(SUM(a**2))
    
    Lr = FOS_Characteristic_length(a, nu)
    Tr = Lr/(mod_a + nu/Lr)
    
    ff(1, :) = (/  a(1)*u - nu*p,    a(2)*u - nu*q  /)
    ff(2, :) = (/ -u/Tr,            (b*p + c*q)/Tr  /)
    ff(3, :) = (/ -(b*p + c*q)/Tr,  -u/Tr           /)

  END FUNCTION FOS_advection_flux
  !==============================

  !====================================================
  FUNCTION fos_source(type_pb, uu, nu, x, y) RESULT(ss)
  !====================================================

    IMPLICIT NONE
    
    INTEGER,                    INTENT(IN) :: type_pb
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: uu
    REAL(KIND=8),               INTENT(IN) :: nu
    REAL(KIND=8),               INTENT(IN) :: x, y
    
    REAL(KIND=8), DIMENSION(SIZE(UU)) :: SS
    !---------------------------------------------

    REAL(KIND=8), DIMENSION(N_DIM) :: A

    REAL(KIND=8) :: u, p, q, Tr, Lr, mod_a
    !---------------------------------------------

    u = uu(1); p = uu(2); q = uu(3)

    a = advection_speed(type_pb, u, x, y)
    mod_a = DSQRT(SUM(a**2))

    Lr = fos_characteristic_length(a, nu)
    Tr = Lr/(mod_a + nu/Lr)

    ss(1) = source_term(type_pb, u, nu, x, y)    
    ss(2) = -p/Tr
    ss(3) = -q/Tr

  END FUNCTION FOS_SOURCE
  !======================

  !================================================
  FUNCTION FOS_eigenvalues(a, nu, n) RESULT(lambda)
  !================================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: a
    REAL(KIND=8),               INTENT(IN) :: nu
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: n

    REAL(KIND=8), DIMENSION(3) :: lambda
    !-------------------------------------------
    
    REAL(KINd=8) :: b, c

    REAL(KIND=8) :: An, An_m, An_p, Re_m, &
                    Re_p, Lr, at, Tr
    !---------------------------------------------
    
    b = 0.d0; c = 0.d0

    an = DOT_PRODUCT(a, n)

    an_m = MIN(0.d0, an)
    an_p = MAX(0.d0, an)

    Lr = FOS_characteristic_length(a, nu)

    Tr = Lr / (an_p - an_m + nu/Lr)

    !RE_M = AN_M * LR / NU
    !RE_P = AN_P * LR / NU

    lambda(1) = an_m - nu/Lr ! an_m*(1 - 1/re_m)
    lambda(2) = an_p + nu/Lr ! an_p*(1 + 1/re_p)
    lambda(3) = (b*n(2) - c*n(1))/Tr

  END FUNCTION FOS_eigenvalues
  !===========================  
  
  !===================================================
  FUNCTION FOS_right_eigenvectors(a, nu, n) RESULT(RR)
  !===================================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: a
    REAL(KIND=8),               INTENT(IN) :: nu
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: n

    REAL(KIND=8), DIMENSION(3,3) :: RR
    !-------------------------------------------

    REAL(KIND=8) :: b, c
    
    REAL(KIND=8), DIMENSION(3) :: ll
    
    REAL(KIND=8) :: An, An_m, An_p, &
                    Re_m, Re_p, Lr, Tr
    !---------------------------------------------
    
    b = 0.d0; c = 0.d0
    
    an = DOT_PRODUCT(a, n)

    an_m = MIN(0.d0, an)
    an_p = MAX(0.d0, an)

    Lr = fos_characteristic_length(a, nu)
    Tr = Lr/(an_p - an_m + nu/Lr)

    ll = FOS_eigenvalues(a, nu, n)

    RR(1, :) = (/  Tr*(ll(1) - ll(3)),   Tr*(ll(2) - ll(3)),  0.d0 /)
    RR(2, :) = (/  c/(ll(1)*Tr) + n(1),  c/(ll(2)*Tr) + n(1),-n(2) /)
    RR(3, :) = (/ -b/(ll(1)*Tr) + n(2), -b/(ll(2)*Tr) + n(2), n(1) /)

  END FUNCTION FOS_right_eigenvectors
  !==================================
  
  !==================================================
  FUNCTION FOS_left_eigenvectors(a, nu, n) RESULT(LL)
  !==================================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: a
    REAL(KIND=8), INTENT(IN) :: nu
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: n

    REAL(KIND=8), DIMENSION(3,3) :: LL
    !-------------------------------------------

    REAL(KIND=8), DIMENSION(3,3) :: RR

    REAL(KIND=8) :: an, an_m, an_p, Re_m, Re_p, Re_a, LR
    !---------------------------------------------

    RR = FOS_right_eigenvectors(a, nu, n)
    
    LL = inverse(RR)

  END FUNCTION FOS_left_eigenvectors
  !=================================

  !======================================
  FUNCTION FOS_Jacobian(a, nu) RESULT(AA)
  !======================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: a
    REAL(KIND=8),               INTENT(IN) :: nu

    REAL(KIND=8), DIMENSION(3,3, N_dim) :: AA
    !-------------------------------------------

    REAL(KIND=8) :: Tr, Lr, mod_a

    REAL(KIND=8) :: b, c
    !-------------------------------------------

    b = 0.d0; c = 0.d0

    mod_a = DSQRT(SUM(a**2))

    Lr = fos_characteristic_length(a, nu)
    Tr = Lr/(mod_a + nu/Lr)

    ! a_x
    !--------------------------------------
    AA(1,:, 1) = (/ a(1),    -nu,       0.d0  /)
    AA(2,:, 1) = (/-1.d0/Tr,  0.d0,     0.d0  /)
    AA(3,:, 1) = (/ 0.d0,    -b/Tr,     -c/Tr /)

    ! a_y
    !--------------------------------------
    AA(1,:, 2) = (/ a(2),    0.d0, -nu   /)
    AA(2,:, 2) = (/ 0.d0,    b/Tr,  c/Tr /)
    AA(3,:, 2) = (/-1.d0/Tr, 0.d0,  0.d0 /)

  END FUNCTION FOS_Jacobian
  !======================== 

  !==================================================
  FUNCTION FOS_Characteristic_length(a, nu) RESULT(L)
  !==================================================
  !
  ! H. Nishikawa: A first-order system approach for
  ! diffusion equation. II: Unification of advection 
  ! and diffusion. JCP 2010.
  !
    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: a
    REAL(KIND=8),               INTENT(IN) :: nu

    REAL(KIND=8) :: L
    !-------------------------------------------
    REAL(KIND=8) :: Re_pi
    !-------------------------------------------

    Re_pi = DSQRT(SUM(a**2)) / (nu * Pi)
    
    L = ( (RE_pi/( DSQRT(1.d0 + Re_pi**2) + 1.d0 )) + &
          (DSQRT(1.d0 + 2.d0/( DSQRT(1.d0 + Re_pi**2) +1.d0))) ) / (2.d0*PI)

  END FUNCTION FOS_Characteristic_length
  !=====================================


END MODULE FOS_system
