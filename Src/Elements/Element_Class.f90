MODULE Element_Class

  USE Face_Class
  
  IMPLICIT NONE

  !==================================  
  INTEGER, PARAMETER :: TRI_P1 = 10
  INTEGER, PARAMETER :: TRI_P2 = 11
  
  INTEGER, PARAMETER :: QUA_Q1 = 20
  INTEGER, PARAMETER :: QUA_Q2 = 21
  !==================================

  TYPE :: faces_ptr     
     CLASS(face), POINTER :: f     
  END TYPE faces_ptr
       
  !==================================
  TYPE, PUBLIC :: element

     INTEGER :: Type      ! Element type
     INTEGER :: N_dim = 2 ! # spatial dimension
     INTEGER :: N_verts   ! # vertices
     INTEGER :: N_points  ! # DoFs
     INTEGER :: N_Faces   ! # faces
     
     !------------------
     ! Element geometry
     !-------------------------------------------------
     ! Vertices coordinates
     ! Coords(id, j): id = 1, N_dim, j = 1, N_verts
     REAL(KIND=8), DIMENSION(:,:), POINTER :: Coords
     !-------------------------------------------------
     ! Volume of the face
     REAL(KIND=8) :: Volume

     !--------------
     ! Element Topology
     !-------------------------------------------------
     ! Nodes (global numeration of the mesh)
     ! NU(j), j = 1, N_points
     INTEGER,      DIMENSION(:),   POINTER :: NU
     !-------------------------------------------------
     ! Normals for the RD scheme
     ! rd_n(id, j), id = 1, N_dim, j = 1, N_points
     !-------------------------------------------------
     REAL(KIND=8), DIMENSION(:,:), POINTER :: rd_n

     !------------
     ! Quadrature
     !------------------------------------------------------------
     ! # quadrature points
     INTEGER :: N_quad
     !------------------------------------------------------------
     ! Quadrature weighs (multiplied by the jacobian)
     ! w_q(iq), iq = 1, N_quad
     !------------------------------------------------------------
     REAL(KIND=8), DIMENSION(:),     POINTER :: w_q
     ! Reference Coordinates of the quadrature points
     ! xx_q(ik, iq), , ik = 1, N_dim_ele, iq = 1, N_quad
     !------------------------------------------------------------
     REAL(kind=8), DIMENSION(:,:),   POINTER :: x_q    
     ! Physical coordinates of the quadrature points
     ! xx_q(id, iq), , id = 1, N_dim, iq = 1, N_quad
     REAL(kind=8), DIMENSION(:,:),   POINTER :: xx_q
     !------------------------------------------------------------
     ! Value of basis functions at the quadrature point
     ! phi_q(i, iq), i = 1, N_points,  iq = 1, N_quad
     REAL(KIND=8), DIMENSION(:,:),   POINTER :: phi_q
     !------------------------------------------------------------
     ! Gradient of basis functions at the quadrature point
     ! D_phi_q(id, i, iq), 
     ! id =1, N_dim, i = 1, N_points,  iq = 1, N_quad
     REAL(KIND=8), DIMENSION(:,:,:), POINTER :: D_phi_q

     ! Gradient of the basis functions at the DOFs
     !  D_phi_k(id, i, k)
     ! id =1, N_dim, i = 1, N_points,  k = 1, N_points
     REAL(KIND=8), DIMENSION(:,:,:), POINTER :: D_phi_k
     
     !---------------
     !Faces structure
     !---------------------------------------------------
     TYPE(faces_ptr), DIMENSION(:), ALLOCATABLE :: faces
     
!!$   CONTAINS

  END TYPE element  
  !==================================
    
CONTAINS
     
END MODULE Element_Class
