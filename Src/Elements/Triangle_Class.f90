MODULE Trianlge_Class

  USE Element_Class
  USE Segment_Class
  USE Lin_Algebra

  IMPLICIT NONE

  !=========================================
  TYPE, PUBLIC, EXTENDS(element) :: triangle

     ! < >

   CONTAINS

     PROCEDURE, PUBLIC :: initialize => initialize_sub

  END TYPE triangle
  !=========================================

  PRIVATE :: initialize_sub

  PRIVATE :: init_faces
  PRIVATE :: init_faces_TRI_P1
  PRIVATE :: init_faces_TRI_P2

  PRIVATE :: volume_quadrature

  PRIVATE :: basis_function
  PRIVATE :: basis_function_TRI_P1
  PRIVATE :: basis_function_TRI_P2
  
  PRIVATE :: gradient
  
  PRIVATE :: gradient_ref
  PRIVATE :: gradient_ref_TRI_P1
  PRIVATE :: gradient_ref_TRI_P2

  PRIVATE :: DOFs_ref
  PRIVATE :: DOFs_ref_TRI_P1
  PRIVATE :: DOFs_ref_TRI_P2
  
  PRIVATE :: gradient_trace
  PRIVATE :: face_trace
  PRIVATE :: face_trace_TRI_P1
  PRIVATE :: face_trace_TRI_P2  

  PRIVATE :: rd_normal
  PRIVATE :: rd_normal_TRI_P1
  PRIVATE :: rd_normal_TRI_P2

  PRIVATE :: Compute_Jacobian

CONTAINS

  !===============================================================
  SUBROUTINE initialize_sub(e, mode, Nodes, Coords, Nu_seg, n_ele)
  !===============================================================

    IMPLICIT NONE

    CLASS(triangle)                          :: e
    CHARACTER(*),                 INTENT(IN) :: mode
    INTEGER,      DIMENSION(:),   INTENT(IN) :: Nodes
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Coords
    INTEGER,      DIMENSION(:),   INTENT(IN) :: NU_seg
    INTEGER,      DIMENSION(:),   INTENT(IN) :: n_ele
    !-------------------------------------------------

    INTEGER :: i, id
    !-------------------------------------------------

    !---------------------------
    ! Attach data to the element
    !---------------------------------------------------
    ALLOCATE( e%Coords(SIZE(Coords,1), SIZE(Coords,2)) )
    ALLOCATE( e%NU(SIZE(Nodes)) )
     
    DO id = 1, SIZE(Coords,1)
       e%Coords(id, :) = Coords(id, :)
    ENDDO

    e%NU = Nodes

    SELECT CASE( SIZE(Nodes) )

    CASE(3)

       !---------
       ! ELELMENT
       !----------------------------------
       e%Type     = TRI_P1 ! Element type
       e%N_verts  = 3      ! # Vertices 
       e%N_points = 3      ! # DoFs       

    CASE(6)

       !---------
       ! ELELMENT
       !----------------------------------
       e%Type     = TRI_P2 ! Element type
       e%N_verts  = 3      ! # Vertices 
       e%N_points = 6      ! # DoFs

     CASE DEFAULT

        WRITE(*,*) 'ERROR: Unsupported Triangle type'
        STOP

     END SELECT

     IF( mode == "element") THEN

        ! Itialize the faces of the element
        CALL init_faces(e, Nodes, Coords, NU_seg, n_ele)

        ! Store the informations at the quadrature points
        CALL volume_quadrature(e)
        
        !--------------------------
        ! Normals for the RD scheme
        !--------------------------------------
        ALLOCATE( e%rd_n(e%N_dim, e%N_points) )

        DO i = 1, e%N_points
           e%rd_n(:, i) = rd_normal(e, i)
        ENDDO

        CALL nodal_gradients(e)

     ENDIF
         
  END SUBROUTINE initialize_sub
  !============================

  !=====================================================
  SUBROUTINE init_faces(e, Nodes, Coords, NU_seg, n_ele)
  !=====================================================

    IMPLICIT NONE
    
    CLASS(triangle)                          :: e
    INTEGER,      DIMENSION(:),   INTENT(IN) :: Nodes
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Coords
    INTEGER,      DIMENSION(:),   INTENT(IN) :: NU_seg
    INTEGER,      DIMENSION(:),   INTENT(IN) :: n_ele
    !-------------------------------------------------

    SELECT CASE(e%Type)
              
    CASE(TRI_P1)

       CALL init_faces_TRI_P1(e, Nodes, Coords, NU_seg, n_ele)
       
    CASE(TRI_P2)

       CALL init_faces_TRI_P2(e, Nodes, Coords, NU_seg, n_ele)
       
    CASE DEFAULT

       WRITE(*,*) 'Unknown Tiangle type for face initialization'
       WRITE(*,*) 'STOP'

    END SELECT 

  END SUBROUTINE init_faces
  !=========================

  !==============================
  SUBROUTINE volume_quadrature(e)
  !==============================
  !
  ! Store the following  quantities at the quadrature points
  !    - value of the basis function
  !    - normal versor
  !    - weight of the quadrature formula (multipplied by 
  !      the jacobian of the transformation)
  !    - value of the trace of the gradient of basis functions
  !    - physical coordiantes of the quadrature point
  !
  ! Compute the lenght of the segment
  !
    IMPLICIT NONE

    CLASS(triangle) :: e
    !-----------------------------------------------

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: JJ

    INTEGER :: iq, k
    !-----------------------------------------------

    SELECT CASE(e%Type)

    CASE(TRI_P1)

       CALL init_quadrature_TRI_P1(e)
       
    CASE(TRI_P2)

       CALL init_quadrature_TRI_P2(e)
       
    CASE DEFAULT

       WRITE(*,*) 'Unknown Triangle type for quadrature'
       WRITE(*,*) 'STOP'

    END SELECT

    !-------------------------------------    
    ! Attach data to the quadrature points
    !---------------------------------------------------
    ALLOCATE( JJ(e%N_dim, e%N_dim) )

    DO iq = 1, e%N_quad

       DO k = 1, e%N_points

          e%phi_q(k, iq) = basis_function( e, k, e%x_q(iq, :) )

          e%xx_q(:, iq) = e%xx_q(:, iq) + &               
                         basis_function( e, k, e%x_q(iq, :) ) * e%Coords(:, k)

          e%D_phi_q(:, k, iq) = gradient( e, k, e%x_q(iq, :) )

       ENDDO

       ! warning: not optimized, jacobian computed twice       
       JJ = Compute_Jacobian( e, e%x_q(iq, :) )

       e%w_q(iq) = e%w_q(iq) * determinant(JJ)
      
    ENDDO

    DEALLOCATE( JJ )

    ! Area of the element
    e%volume = SUM( e%w_q )    

  END SUBROUTINE volume_quadrature
  !===============================

  !============================
  SUBROUTINE nodal_gradients(e)
  !============================
  !
  ! Compute the gradient of the basis functions
  ! at the DOFs of the element
  !
    IMPLICIT NONE

    CLASS(triangle) :: e    
    !-----------------------------------------------

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: x_dof

    INTEGER :: i, k
    !-----------------------------------------------

    ALLOCATE( x_dof(e%N_points, 3) )
    
    x_dof = DOFs_ref(e)
    
    ALLOCATE( e%D_phi_k(e%N_dim, e%N_points, e%N_points) )

    DO i = 1, e%N_points

       DO k = 1, e%N_points

          e%D_phi_k(:, i, k) = gradient( e, i, x_dof(k, :) )

       ENDDO

    ENDDO
    
    DEALLOCATE(x_dof)

  END SUBROUTINE nodal_gradients
  !=============================

!*******************************************************************************
!*******************************************************************************
!                         COMPLEMENTARY FUNCTIONS                              !
!*******************************************************************************
!*******************************************************************************

  !==============================================
  FUNCTION basis_function(e, i, xi) RESULT(psi_i)
  !==============================================

    IMPLICIT NONE

    CLASS(triangle)                        :: e
    INTEGER,                    INTENT(IN) :: i
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8) :: psi_i
    !-------------------------------------------
      
    SELECT CASE(e%Type)
       
    CASE(TRI_P1)

       psi_i = basis_function_TRI_P1(e, i, xi)
       
    CASE(TRI_P2)

       psi_i = basis_function_TRI_P2(e, i, xi)
       
    CASE DEFAULT

       WRITE(*,*) 'Unknown Tiangle type for basis functions'
       WRITE(*,*) 'STOP'

    END SELECT

  END FUNCTION basis_function
  !==========================

  !==============================================
  FUNCTION gradient_ref(e, i, xi) RESULT(D_psi_i)
  !==============================================

    IMPLICIT NONE

    CLASS(triangle)                         :: e
    INTEGER,                     INTENT(IN) :: i    
    REAL(KIND=8),  DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8), DIMENSION(2) :: D_psi_i
    !-------------------------------------------
      
    SELECT CASE(e%Type)
       
    CASE(TRI_P1)

       D_psi_i = gradient_ref_TRI_P1(e, i, xi)
       
    CASE(TRI_P2)

       D_psi_i = gradient_ref_TRI_P2(e, i, xi)
     
    CASE DEFAULT

       WRITE(*,*) 'Unknown Triangle type for gradients'
       WRITE(*,*) 'STOP'

    END SELECT    

  END FUNCTION gradient_ref
  !========================

  !==========================================
  FUNCTION gradient(e, i, xi) RESULT(D_psi_i)
  !==========================================
  !
  ! Compute the gradient of the shape function i
  ! on the actual triangle at the point of baricentric
  ! coordinates xi
  !  
    IMPLICIT NONE

    CLASS(triangle)                         :: e
    INTEGER,                     INTENT(IN) :: i    
    REAL(KIND=8),  DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8), DIMENSION(e%N_dim) :: D_psi_i
    !-------------------------------------------

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: JJ, inv_J

    INTEGER :: k, l, m
    !---------------------------------------------

    ALLOCATE(    JJ(e%N_dim, e%N_dim) )
    ALLOCATE( inv_J(e%N_dim, e%N_dim) )
        
    JJ = compute_Jacobian( e, xi )
    
    inv_J = inverse(JJ)

    D_Psi_i = MATMUL( inv_J, gradient_ref(e, i, xi) )

    DEALLOCATE( JJ, inv_J )

  END FUNCTION gradient
  !====================
  
  !====================================
  FUNCTION rd_normal(e, i) RESULT(nn_i)
  !====================================

    IMPLICIT NONE

    CLASS(triangle)     :: e
    INTEGER, INTENT(IN) :: i

    REAL(KIND=8), DIMENSION(e%N_dim) :: nn_i
    !---------------------------------------    

    SELECT CASE(e%Type)

    CASE(TRI_P1)

       nn_i = rd_normal_TRI_P1(e, i)

    CASE(TRI_P2)

       nn_i = rd_normal_TRI_P2(e, i)

    CASE DEFAULT

       WRITE(*,*) 'ERROR: Unsupported Triangle type'
       STOP

    END SELECT
    
  END FUNCTION rd_normal
  !=====================

  !========================================
  SUBROUTINE gradient_trace(e, if, p_D_phi)
  !========================================
  !
  ! Compute the trace of the gradients of all shape
  ! functions on the face if
  !
    IMPLICIT NONE

    CLASS(triangle)                             :: e
    INTEGER,                        INTENT(IN)  :: if
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: p_D_phi
    !-----------------------------------------------------

    REAL(kind=8), DIMENSION(:,:), ALLOCATABLE :: xi_f
    INTEGER :: N_points_ele, N_points_face, k, i
    !-----------------------------------------------------

    N_points_ele  = e%N_points
    N_points_face = SIZE(p_D_phi, 2)

    ALLOCATE( xi_f(N_points_face, 3) )

    xi_f = face_trace(e, if)

    DO k = 1, N_points_face      

       DO i = 1, N_points_ele

          p_D_phi(:, k, i) = gradient( e, i, xi_f(k,:) )

       ENDDO

    ENDDO
    
    DEALLOCATE(xi_f)    

  END SUBROUTINE gradient_trace
  !============================

  !======================================
  FUNCTION face_trace(e, if) RESULT(xi_f)
  !======================================
  !
  ! Compute the baricentric coordinates on the
  ! face if
  !
    IMPLICIT NONE

    CLASS(triangle)          :: e
    INTEGER,      INTENT(IN) :: if
    
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: xi_f
    !---------------------------------------------
     
    SELECT CASE(e%Type)
       
    CASE(TRI_P1)

       ALLOCATE( xi_f(2, 3) )

       xi_f = face_trace_TRI_P1(e, if)
       
    CASE(TRI_P2)

       ALLOCATE( xi_f(3, 3) )

       xi_f = face_trace_TRI_P2(e, if)
       
    CASE DEFAULT

       WRITE(*,*) 'Unknown Tiangle type for face trace'
       WRITE(*,*) 'STOP'

    END SELECT 

  END FUNCTION face_trace
  !======================

  !=================================
  FUNCTION DOFs_ref(e) RESULT(x_dof)
  !=================================
  !
  ! Compute the coordinates of the DOFs
  ! on the reference element
  !
    IMPLICIT NONE

    CLASS(triangle) :: e
    
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: x_dof
    !-------------------------------------------------
     
    SELECT CASE(e%Type)
       
    CASE(TRI_P1)

       ALLOCATE( x_dof(3, 3) )

       x_dof = DOFs_ref_TRI_P1(e)
       
    CASE(TRI_P2)

       ALLOCATE( x_dof(6, 3) )

       x_dof = DOFs_ref_TRI_P2(e)
       
    CASE DEFAULT

       WRITE(*,*) 'Unknown Tiangle type for DOFs'
       WRITE(*,*) 'STOP'

    END SELECT 

  END FUNCTION DOFs_ref
  !====================

  !==========================================
  FUNCTION Compute_Jacobian(e, xi) RESULT(JJ)
  !==========================================

    IMPLICIT NONE

    CLASS(triangle)                         :: e
    REAL(KIND=8),  DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8),  DIMENSION(e%N_dim, e%N_dim) :: JJ
    !------------------------------------------------

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: d

    INTEGER :: k, l, m
    !------------------------------------------------

    ALLOCATE( d(e%N_dim, e%N_points) )

    DO m = 1, e%N_points
       d(:, m) = gradient_ref(e, m, xi)
    ENDDO

    ! Construction of the Jacobian matrix
    DO k = 1, e%N_dim

       DO l = 1, e%N_dim

          JJ(l, k) = SUM( e%Coords(k, :) * d(l, :) )

       ENDDO

    ENDDO

    DEALLOCATE( d )
    
  END FUNCTION Compute_Jacobian
  !============================  

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%% SPECIFIC FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%     TRIANGLE P1    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !============================================================
  SUBROUTINE init_faces_TRI_P1(e, Nodes, Coords, NU_seg, n_ele)
  !============================================================

    IMPLICIT NONE

    CLASS(triangle)                          :: e
    INTEGER,      DIMENSION(:),   INTENT(IN) :: Nodes
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Coords
    INTEGER,      DIMENSION(:),   INTENT(IN) :: NU_seg
    INTEGER,      DIMENSION(:),   INTENT(IN) :: n_ele
    !-------------------------------------------------

     TYPE(segment), POINTER :: seg

    INTEGER,      DIMENSION(:),   ALLOCATABLE :: loc
    INTEGER,      DIMENSION(:),   ALLOCATABLE :: VV
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: RR

    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: p_D_phi

    INTEGER :: id, if, i, j, k, N_ln, istat
    !-------------------------------------------------

    e%N_faces = 3 ! # Faces
       
    ALLOCATE( e%faces(e%N_faces) )

    ! # of local nodes    
    N_ln = 2
    
    ! DoFs on the faces
    ALLOCATE( VV(N_ln), loc(N_ln) )

    ! Coordinates on the faces
    ALLOCATE( RR(SIZE(Coords, 1), N_ln) )

    ALLOCATE( p_D_phi(e%N_dim, N_ln, e%N_points) )

    DO if = 1, e%N_faces

       i = MOD(if, MIN(if+1, e%N_faces)) + 1
       j = MOD(i,  MIN(if+2, e%N_faces)) + 1

       loc = (/ i, j /)

       VV = Nodes(loc)
       RR = Coords(:, loc)

       ALLOCATE(seg, STAT=istat)          
       IF(istat /=0) THEN
          WRITE(*,*) 'ERROR: failed trianlge faces allocation'
       ENDIF
       
       CALL seg%initialize( loc, VV, RR, NU_seg(if), n_ele(if) )

       CALL gradient_trace(e, if, p_D_phi)

       CALL seg%face_quadrature(p_D_phi)

       e%faces(if)%f => seg
          
    ENDDO

    DEALLOCATE( VV, RR, loc, p_D_phi )

  END SUBROUTINE init_faces_TRI_P1
  !===============================

  !=====================================================
  FUNCTION basis_function_TRI_P1(e, i, xi) RESULT(psi_i)
  !=====================================================

    IMPLICIT NONE

    CLASS(triangle)                        :: e
    INTEGER,                    INTENT(IN) :: i
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8) :: psi_i
    !-------------------------------------------

    SELECT CASE (i)
       
    CASE(1,2,3)
       
       psi_i = xi(i)

    CASE DEFAULT

       WRITE(*,*) 'ERROR: not supported Dof in Triangle  baisis function'
       STOP
       
    END SELECT

  END FUNCTION basis_function_TRI_P1
  !=================================

  !=====================================================
  FUNCTION gradient_ref_TRI_P1(e, i, xi) RESULT(D_psi_i)
  !=====================================================

    IMPLICIT NONE

    CLASS(triangle)                        :: e
    INTEGER,                    INTENT(IN) :: i
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8), DIMENSION(2) :: D_psi_i
    !-------------------------------------------

    SELECT CASE(i)

    CASE(1)

       D_psi_i = (/ -1.d0, -1.d0 /)
       
    CASE(2)

       D_psi_i = (/ 1.d0, 0.d0 /)
       
    CASE(3)

       D_psi_i = (/ 0.d0, 1.d0 /)

    CASE DEFAULT

       WRITE(*,*) 'ERROR: not supported Dof in Triangle gradient'
       STOP

    END SELECT
    
  END FUNCTION gradient_ref_TRI_P1
  !===============================

  !========================================
  FUNCTION DOFs_ref_TRI_P1(e) RESULT(x_dof)
  !========================================

    IMPLICIT NONE

    CLASS(triangle):: e

    REAL(KIND=8), DIMENSION(3,3) :: x_dof
    !------------------------------------

    x_dof(1, :) = (/ 1.d0,  0.d0,  0.d0  /)
    x_dof(2, :) = (/ 0.d0,  1.d0,  0.d0  /)
    x_dof(3, :) = (/ 0.d0,  0.d0,  1.d0  /)

  END FUNCTION DOFs_ref_TRI_P1
  !===========================
  
  !===========================================
  FUNCTION rd_normal_TRI_P1(e, k) RESULT(nn_k)
  !===========================================

    IMPLICIT NONE

    CLASS(triangle)     :: e
    INTEGER, INTENT(IN) :: k

    REAL(KIND=8), DIMENSION(e%N_dim) :: nn_k
    !---------------------------------------

    INTEGER :: i, j
    !---------------------------------------

    i = e%faces(k)%f%l_NU(1)
    j = e%faces(k)%f%l_NU(2)

    nn_k(1) = e%Coords(2,i) - e%Coords(2,j)
    nn_k(2) = e%Coords(1,j) - e%Coords(1,i)

  END FUNCTION rd_normal_TRI_P1
  !============================

  !=============================================
  FUNCTION face_trace_TRI_P1(e, if) RESULT(xi_f)
  !=============================================

    IMPLICIT NONE

    CLASS(triangle)          :: e
    INTEGER,      INTENT(IN) :: if

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: xi_f
    !----------------------------------------------

    ALLOCATE( xi_f(2, 3) )

    SELECT CASE(if)

    CASE(1)

       xi_f(1, :) = (/ 0.d0, 1.d0, 0.d0 /)
       xi_f(2, :) = (/ 0.d0, 0.d0, 1.d0 /)

    CASE(2)

       xi_f(1, :) = (/ 0.d0, 0.d0, 1.d0 /)
       xi_f(2, :) = (/ 1.d0, 0.d0, 0.d0 /)       
       
    CASE(3)
  
       xi_f(1, :) = (/ 1.d0, 0.d0, 0.d0 /)
       xi_f(2, :) = (/ 0.d0, 1.d0, 0.d0 /)       
            
    CASE DEFAULT

       WRITE(*,*) 'ERROR: wrong Triangle face in face trace P1'
       STOP

    END SELECT    
    
  END FUNCTION face_trace_TRI_P1
  !=============================

  !===================================
  SUBROUTINE init_quadrature_TRI_P1(e)
  !===================================

    IMPLICIT NONE

    CLASS(triangle) :: e
    !----------------------------------------------

    REAL(KIND=8) :: alpha_1, beta_1
    !-----------------------------------------------

    e%N_quad = 6

    ALLOCATE( e%phi_q(e%N_points, e%N_quad) )

    ALLOCATE( e%D_phi_q(e%N_dim, e%N_points, e%N_quad) )

    ALLOCATE( e%w_q(e%N_quad   ) )
    ALLOCATE( e%x_q(e%N_quad, 3) )

    ALLOCATE( e%xx_q(e%N_dim, e%N_quad) )
    
    !-------------------
    ! Quadrature formula
    !--------------------------------------
    alpha_1 = 0.816847572980459d0
    beta_1  = 0.091576213509771d0

    e%x_q(1,:) = (/ alpha_1, beta_1 , beta_1  /)
    e%x_q(2,:) = (/ beta_1 , alpha_1, beta_1  /)
    e%x_q(3,:) = (/ beta_1 , beta_1 , alpha_1 /)
    e%w_q(1:3) = 0.109951743655322d0
    
    alpha_1 = 0.108103018168070d0
    beta_1  = 0.445948490915965d0
    
    e%x_q(4,:) = (/ alpha_1, beta_1 , beta_1  /)
    e%x_q(5,:) = (/ beta_1 , alpha_1, beta_1  /)
    e%x_q(6,:) = (/ beta_1 , beta_1 , alpha_1 /)    
    e%w_q(4:6) = 0.223381589678011d0
    !--------------------------------------
    
    e%w_q = e%w_q*0.5d0 ! J -> 0.5*det|J|  
    
    e%xx_q = 0.d0
    
  END SUBROUTINE init_quadrature_TRI_P1
  !====================================

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%% SPECIFIC FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%     TRIANGLE P2    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  !============================================================
  SUBROUTINE init_faces_TRI_P2(e, Nodes, Coords, NU_seg, n_ele)
  !============================================================

    IMPLICIT NONE

    CLASS(triangle)                          :: e
    INTEGER,      DIMENSION(:),   INTENT(IN) :: Nodes
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Coords
    INTEGER,      DIMENSION(:),   INTENT(IN) :: NU_seg
    INTEGER,      DIMENSION(:),   INTENT(IN) :: n_ele
    !-------------------------------------------------

    TYPE(segment), POINTER :: seg

    INTEGER,      DIMENSION(:),   ALLOCATABLE :: loc
    INTEGER,      DIMENSION(:),   ALLOCATABLE :: VV
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: RR

    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: p_D_phi

    INTEGER :: id, if, i, j, k, N_ln, istat
    !-------------------------------------------------

    e%N_faces  = 3 ! # Faces

    ALLOCATE( e%faces(e%N_faces) )

    ! # of local nodes
    N_ln = 3

    ! Dofs on the faces
    ALLOCATE( VV(N_ln), loc(N_ln) )

    ! Coordinates on the faces
    ALLOCATE( RR(SIZE(Coords, 1), N_ln) )

    ALLOCATE( p_D_phi(e%N_dim, N_ln, e%N_points) )
    
    DO if = 1, e%N_faces

       i = MOD(if, MIN(if+1, e%N_faces)) + 1          
       j = MOD(i,  MIN(if+2, e%N_faces)) + 1

       k = i + 3

       loc = (/ i, j, k /)

       VV = Nodes(loc)
       RR = Coords(:, loc)

       ALLOCATE(seg, STAT=istat)          
       IF(istat /=0) THEN
          WRITE(*,*) 'ERROR: failed trianlge allocation'
       ENDIF

       CALL seg%initialize( loc, VV, RR, NU_seg(if), n_ele(if) )

       CALL gradient_trace(e, if, p_D_phi)

       CALL seg%face_quadrature(p_D_phi)
         
       e%faces(if)%f => seg
          
    ENDDO

    DEALLOCATE( VV, RR, loc, p_D_phi )

  END SUBROUTINE init_faces_TRI_P2
  !===============================

  !=====================================================
  FUNCTION basis_function_TRI_P2(e, i, xi) RESULT(psi_i)
  !=====================================================

    IMPLICIT NONE

    CLASS(triangle)                        :: e
    INTEGER,                    INTENT(IN) :: i
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8) :: psi_i
    !-------------------------------------------

    SELECT CASE (i)

    CASE(1,2,3)
       
       psi_i = (2.d0*xi(i) - 1.d0) * xi(i)
                     
    CASE(4)

       psi_i  = 4.d0*xi(1)*xi(2)
              
    CASE(5)

       psi_i  = 4.d0*xi(2)*xi(3)
              
    CASE(6)

       psi_i = 4.d0*xi(3)*xi(1)
       
    CASE DEFAULT

       WRITE(*,*) 'ERROR: not supported Dof in Triangle  baisis function'
       STOP
       
    END SELECT

  END FUNCTION basis_function_TRI_P2
  !=================================

  !=====================================================
  FUNCTION gradient_ref_TRI_P2(e, i, xi) RESULT(D_psi_i)
  !=====================================================

    IMPLICIT NONE

    CLASS(triangle)                        :: e
    INTEGER,                    INTENT(IN) :: i
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8), DIMENSION(2) :: D_psi_i
    !-------------------------------------------

    SELECT CASE(i)

    CASE(1)

       D_psi_i = -(/ 4.d0*xi(1) - 1.d0, 4.d0*xi(1) - 1.d0 /)
       
    CASE(2)

       D_psi_i = (/ 4.d0*xi(2) - 1.d0, 0.d0 /)
       
    CASE(3)

       D_psi_i = (/ 0.d0, 4.d0*xi(3) - 1.d0 /)

    CASE(4)

       D_psi_i = (/ 4.d0*(xi(1)-xi(2)), -4.d0*xi(2) /)

    CASE(5)

       D_psi_i = (/ 4.d0*xi(3), 4.d0*xi(2) /)

    CASE(6)

       D_psi_i = (/ -4.d0*xi(3), 4.d0*(xi(1) - xi(3)) /)

    CASE DEFAULT

       WRITE(*,*) 'ERROR: not supported Dof in Triangle gradient'
       STOP

    END SELECT
    
  END FUNCTION gradient_ref_TRI_P2
  !===============================

  !========================================
  FUNCTION DOFs_ref_TRI_P2(e) RESULT(x_dof)
  !========================================

    IMPLICIT NONE

    CLASS(triangle):: e

    REAL(KIND=8), DIMENSION(6,3) :: x_dof
    !-------------------------------------

    x_dof(1, :) = (/ 1.d0,  0.d0,  0.d0  /)
    x_dof(2, :) = (/ 0.d0,  1.d0,  0.d0  /)
    x_dof(3, :) = (/ 0.d0,  0.d0,  1.d0  /)
    x_dof(4, :) = (/ 0.5d0, 0.5d0, 0.d0  /)
    x_dof(5, :) = (/ 0.d0,  0.5d0, 0.5d0 /)
    x_dof(6, :) = (/ 0.5d0, 0.d0,  0.5d0 /)

  END FUNCTION DOFs_ref_TRI_P2
  !===========================  

  !===========================================
  FUNCTION rd_normal_TRI_P2(e, k) RESULT(nn_k)
  !===========================================

    IMPLICIT NONE

    CLASS(triangle)     :: e
    INTEGER, INTENT(IN) :: k

    REAL(KIND=8), DIMENSION(e%N_dim) :: nn_k
    !---------------------------------------

    INTEGER :: i, j, kk
    !---------------------------------------

    SELECT CASE(k)

    CASE(1, 2, 3)

       i = e%faces(k)%f%l_NU(1)
       j = e%faces(k)%f%l_NU(2)

       nn_k(1) = e%Coords(2,i) - e%Coords(2,j)
       nn_k(2) = e%Coords(1,j) - e%Coords(1,i)       

    CASE(4, 5, 6)

       kk = MOD(k+1, MIN(k+2, e%N_faces)) + 1

       i = e%faces(kk)%f%l_NU(1)      
       j = e%faces(kk)%f%l_NU(2)
       
       nn_k(1) = e%Coords(2,i) - e%Coords(2,j)       
       nn_k(2) = e%Coords(1,j) - e%Coords(1,i)

       nn_k = -nn_k       

    CASE DEFAULT

       WRITE(*,*) 'ERROR: wrong nodes for rd_normal'
       STOP

    END SELECT
       
  END FUNCTION rd_normal_TRI_P2
  !============================

  !=============================================
  FUNCTION face_trace_TRI_P2(e, if) RESULT(xi_f)
  !=============================================

    IMPLICIT NONE

    CLASS(triangle)          :: e
    INTEGER,      INTENT(IN) :: if

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: xi_f
    !----------------------------------------------

    ALLOCATE( xi_f(3, 3) )

    SELECT CASE(if)

    CASE(1)

       xi_f(1, :) = (/ 0.d0, 1.d0,  0.d0  /)
       xi_f(2, :) = (/ 0.d0, 0.d0,  1.d0  /)
       xi_f(3, :) = (/ 0.d0, 0.5d0, 0.5d0 /)

    CASE(2)

       xi_f(1, :) = (/ 0.d0,  0.d0, 1.d0  /)
       xi_f(2, :) = (/ 1.d0,  0.d0, 0.d0  /)       
       xi_f(3, :) = (/ 0.5d0, 0.d0, 0.5d0 /)

    CASE(3)
  
       xi_f(1, :) = (/ 1.d0,  0.d0,  0.d0 /)
       xi_f(2, :) = (/ 0.d0,  1.d0,  0.d0 /)
       xi_f(3, :) = (/ 0.5d0, 0.5d0, 0.d0 /)
            
    CASE DEFAULT

       WRITE(*,*) 'ERROR: wrong Triangle face in face trace P2'
       STOP

    END SELECT    
    
  END FUNCTION face_trace_TRI_P2
  !=============================

  !===================================
  SUBROUTINE init_quadrature_TRI_P2(e)
  !===================================

    IMPLICIT NONE

    CLASS(triangle) :: e
    !----------------------------------------------------

    REAL(KIND=8) :: alpha_1, beta_1, gamma_1
    !----------------------------------------------------

    e%N_quad = 12

    ALLOCATE( e%phi_q(e%N_points, e%N_quad) )

    ALLOCATE( e%D_phi_q(e%N_dim, e%N_points, e%N_quad) )
    
    ALLOCATE( e%w_q(e%N_quad   ) )
    ALLOCATE( e%x_q(e%N_quad, 3) )

    ALLOCATE( e%xx_q(e%N_dim, e%N_quad) )
    
    !-------------------
    ! Quadrature formula
    !--------------------------------------
    alpha_1 = 0.873821971016996d0
    beta_1  = 0.063089014491502d0

    e%x_q(1,:) = (/ alpha_1, beta_1,  beta_1  /)
    e%x_q(2,:) = (/ beta_1,  alpha_1, beta_1  /)
    e%x_q(3,:) = (/ beta_1,  beta_1,  alpha_1 /)
    e%w_q(1:3) = 0.050844906370207d0
    
    alpha_1 = 0.501426509658179d0
    beta_1  = 0.249286745170910d0

    e%x_q(4,:) = (/ beta_1,  beta_1,  alpha_1 /)
    e%x_q(5,:) = (/ beta_1,  alpha_1, beta_1  /)
    e%x_q(6,:) = (/ alpha_1, beta_1,  beta_1  /)
    e%w_q(4:6) = 0.116786275726379d0

    alpha_1 = 0.636502499121399d0 
    beta_1  = 0.310352451033785d0
    gamma_1 = 0.053145049844816d0 

    e%x_q(7,:)  = (/ alpha_1, beta_1,  gamma_1 /)
    e%x_q(8,:)  = (/ alpha_1, gamma_1, beta_1  /)
    e%x_q(9,:)  = (/ beta_1,  alpha_1, gamma_1 /)
    e%x_q(10,:) = (/ gamma_1, alpha_1, beta_1  /)
    e%x_q(11,:) = (/ beta_1,  gamma_1, alpha_1 /)
    e%x_q(12,:) = (/ gamma_1, beta_1,  alpha_1 /)    
    e%w_q(7:12) = 0.082851075618374d0
    !--------------------------------------
    e%w_q = e%w_q*0.5d0 ! J -> 0.5*det|J|
    
    e%xx_q = 0.d0
    
  END SUBROUTINE init_quadrature_TRI_P2
  !====================================

END MODULE Trianlge_Class
