MODULE Face_Class

  IMPLICIT NONE

  !==================================
  INTEGER, PARAMETER :: SEG_P1 = 1
  INTEGER, PARAMETER :: SEG_P2 = 2
  !==================================

  !==================================
  TYPE, PUBLIC :: face
     
     INTEGER :: Type      ! Element type
     INTEGER :: N_dim = 2 ! # spatial dimension
     INTEGER :: N_verts   ! # vertices
     INTEGER :: N_points  ! # DoFs

     !--------------
     ! Face geometry
     !------------------------------------------------------------
     !
     ! Vertices coordinates
     ! Coords(id, j): id = 1, N_dim, j = 1, N_verts
     REAL(KIND=8), DIMENSION(:,:), POINTER :: Coords
     !------------------------------------------------------------
     !
     ! Area of the face
     REAL(KIND=8) :: Area

     !--------------
     ! Face Topology
     !------------------------------------------------------------
     !
     ! Nodes (global numeration of the mesh)
     ! NU(j), j = 1, N_points
     INTEGER,      DIMENSION(:),   POINTER :: NU
     !------------------------------------------------------------
     !
     ! Nodes (local numeration of the element)
     ! l_nu(j), j = 1, N_points
     INTEGER,      DIMENSION(:),   POINTER :: l_nu
     !------------------------------------------------------------
     ! Segment (global numeration)
     !
     INTEGER                               :: g_seg
     !------------------------------------------------------------
     ! Element wich share the same face
     !
     INTEGER                               :: c_ele
     
     !------------
     ! Quadrature
     !------------------------------------------------------------
     !
     ! # quadrature points
     INTEGER :: N_quad
     !------------------------------------------------------------
     !
     ! Quadrature weighs (multiplied by the jacobian)
     ! w_q(iq), iq = 1, N_quad
     REAL(KIND=8), DIMENSION(:),     POINTER :: w_q
     !------------------------------------------------------------
     !
     ! Reference Coordinates of the quadrature points
     ! xx_q(ik, iq), , ik = 1, N_dim_ele, iq = 1, N_quad
     REAL(kind=8), DIMENSION(:,:),   POINTER :: x_q
     !------------------------------------------------------------
     !
     ! Physical coordinates of the quadrature points
     ! xx_q(id, iq), , id = 1, N_dim, iq = 1, N_quad
     REAL(kind=8), DIMENSION(:,:),   POINTER :: xx_q
     !------------------------------------------------------------
     !
     ! Normal versor at quadrature points
     ! n_q(id, iq), id = 1, N_dim, iq = 1, N_quad
     REAL(KIND=8), DIMENSION(:,:),   POINTER :: n_q
     !------------------------------------------------------------
     !
     ! Value of basis functions at the quadrature point
     ! phi_q(k, iq), k = 1, N_points,  iq = 1, N_quad
     REAL(KIND=8), DIMENSION(:,:),   POINTER :: phi_q

     !--------------
     ! CIP structure
     !------------------------------------------------------------
     !
     ! 1 -> element to which the face belong to 
     ! 2 -> adiancen element which share the same face
     !
     ! Gradient of element basis functions at the quadrature point
     ! p_Dphi_x_q(id, k, iq)
     ! id = 1, N_dim, k = 1, N_points (ele), iq = 1, N_quad
     REAL(KIND=8), DIMENSION(:,:,:), POINTER :: p_Dphi_1_q
     REAL(KIND=8), DIMENSION(:,:,:), POINTER :: p_Dphi_2_q
     !------------------------------------------------------------
     !
     ! Gradient of the solution at the quadrature point
     ! p_u_x_q(id, iq),  id = 1, N_dim, iq = 1, N_quad
     REAL(KIND=8), DIMENSION(:,:), POINTER :: p_Du_1_q
     REAL(KIND=8), DIMENSION(:,:), POINTER :: p_Du_2_q
     !------------------------------------------------------------
     !
     ! Correspondence between the local numeration of two elements
     ! with share the same face
     ! loc_con(k), k = 1, N_points (ele)
     INTEGER, DIMENSION(:), POINTER :: loc_con

!   CONTAINS

  END TYPE face
  !=================================

CONTAINS

END MODULE Face_Class
