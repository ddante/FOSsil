MODULE Quadrangle_Class

  USE Element_Class
  USE Segment_Class

  IMPLICIT NONE

  !=========================================
  TYPE, PUBLIC, EXTENDS(element) :: quadrangle

     ! < >

   CONTAINS

     PROCEDURE, PUBLIC :: initialize => initialize_sub
     PROCEDURE, PUBLIC :: basis_function => basis_function_fun
     PROCEDURE, PUBLIC :: gradient_ref   => gradient_ref_fun
     PROCEDURE, PUBLIC :: gradient       => gradient_fun
     PROCEDURE, PUBLIC :: face_trace     => face_trace_fun
     PROCEDURE, PUBLIC :: gradient_trace => gradient_trace_sub

  END TYPE quadrangle
  !=========================================

  PRIVATE :: initialize_sub
  PRIVATE :: init_faces_QUA_Q1
  PRIVATE :: init_faces_QUA_Q2

  PRIVATE :: basis_function_fun
  PRIVATE :: basis_function_QUA_Q1
  PRIVATE :: basis_function_QUA_Q2

  PRIVATE :: gradient_ref_fun
  PRIVATE :: gradient_ref_QUA_Q1
  PRIVATE :: gradient_ref_QUA_Q2

  PRIVATE :: gradient_fun

  PRIVATE :: face_trace_fun
  PRIVATE :: face_trace_QUA_Q1
  PRIVATE :: face_trace_QUA_Q2

  PRIVATE :: rd_normal_fun
  PRIVATE :: rd_normal_QUA_Q1
  PRIVATE :: rd_normal_QUA_Q2

CONTAINS

  !===============================================================
  SUBROUTINE initialize_sub(e, mode, Nodes, Coords, Nu_seg, n_ele)
  !===============================================================

    IMPLICIT NONE

    CLASS(quadrangle)                        :: e
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

    CASE(4)

       !---------
       ! ELELMENT
       !----------------------------------
       e%Type     = QUA_Q1 ! Element type
       e%N_verts  = 4      ! # Vertices
       e%N_points = 4      ! # DoFs

       IF( mode == "element") THEN
          ! FACES
          CALL init_faces_QUA_Q1(e, Nodes, Coords, NU_seg, n_ele)
       ENDIF

    CASE(9)

       !---------
       ! ELELMENT
       !----------------------------------       
       e%Type     = QUA_Q2 ! Element type
       e%N_verts  = 4      ! # Vertices
       e%N_points = 9      ! # DoFs

       IF( mode == "element") THEN
          ! FACES
          CALL init_faces_QUA_Q2(e, Nodes, Coords, NU_seg, n_ele)
       ENDIF

     CASE DEFAULT

        WRITE(*,*) 'ERROR: Unsupported Quadrangle type'
        STOP

     END SELECT

     !--------------------------
     ! Normals for the RD scheme
     !--------------------------------------
     ALLOCATE( e%rd_n(e%N_dim, e%N_points) )

     DO i = 1, e%N_points
        e%rd_n(:, i) = rd_normal_fun(e, i)
     ENDDO
     
  END SUBROUTINE initialize_sub
  !============================

  !==================================================
  FUNCTION basis_function_fun(e, i, xi) RESULT(psi_i)
  !==================================================

    IMPLICIT NONE

    CLASS(quadrangle)                      :: e
    INTEGER,                    INTENT(IN) :: i
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8) :: psi_i
    !-------------------------------------------
      
    SELECT CASE(e%Type)
       
    CASE(QUA_Q1)

       psi_i = basis_function_QUA_Q1(e, i, xi)
       
    CASE(QUA_Q2)

       psi_i = basis_function_QUA_Q2(e, i, xi)
       
    CASE DEFAULT

       WRITE(*,*) 'Unknown Quadrangle type for basis functions'
       WRITE(*,*) 'STOP'

    END SELECT

  END FUNCTION basis_function_fun
  !============================== 

  !==================================================
  FUNCTION gradient_ref_fun(e, i, xi) RESULT(D_psi_i)
  !==================================================

    IMPLICIT NONE

    CLASS(quadrangle)                       :: e
    INTEGER,                     INTENT(IN) :: i    
    REAL(KIND=8),  DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8), DIMENSION(2) :: D_psi_i
    !-------------------------------------------
      
    SELECT CASE(e%Type)
       
    CASE(QUA_Q1)

       D_psi_i = gradient_ref_QUA_Q1(e, i, xi)
       
    CASE(QUA_Q2)

       D_psi_i = gradient_ref_QUA_Q2(e, i, xi)
     
    CASE DEFAULT

       WRITE(*,*) 'Unknown Quadrangle type for gradients'
       WRITE(*,*) 'STOP'

    END SELECT    

  END FUNCTION gradient_ref_fun
  !============================

  !==============================================
  FUNCTION gradient_fun(e, i, xi) RESULT(D_psi_i)
  !==============================================
  !
  ! Compute the gradient of the shape function i
  ! on the actual quadrangle at the point of baricentric
  ! coordinates xi
  !  
    IMPLICIT NONE

    CLASS(quadrangle)                       :: e
    INTEGER,                     INTENT(IN) :: i    
    REAL(KIND=8),  DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8), DIMENSION(e%N_dim) :: D_psi_i
    !-------------------------------------------

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: d
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: JJ, inv_J

    INTEGER :: k, l, m
    !---------------------------------------------

    ALLOCATE( d(e%N_dim, e%N_points) )

    DO m = 1, e%N_points
       d(:, m) = e%gradient_ref(m, xi)
    ENDDO

    ! Construction of the Jacobian matrix
    ALLOCATE( JJ(e%N_dim, e%N_dim) )
    
    DO k = 1, e%N_dim

       DO l = 1, e%N_dim

          JJ(l, k) = SUM( e%Coords(k, :) * d(l, :) )

       ENDDO

    ENDDO

    ALLOCATE( inv_J(e%N_dim, e%N_dim) )
    
    inv_J = inverse(JJ)

    D_Psi_i = MATMUL( inv_J, d(:, i) )

    DEALLOCATE( d, JJ, inv_J )

  END FUNCTION gradient_fun
  !========================

  !========================================
  FUNCTION rd_normal_fun(e, i) RESULT(nn_i)
  !========================================

    IMPLICIT NONE

    CLASS(quadrangle)   :: e
    INTEGER, INTENT(IN) :: i

    REAL(KIND=8), DIMENSION(e%N_dim) :: nn_i
    !---------------------------------------    

    SELECT CASE(e%Type)

    CASE(QUA_Q1)

       nn_i = rd_normal_QUA_Q1(e, i)

    CASE(QUA_Q2)

       nn_i = rd_normal_QUA_Q2(e, i)

    CASE DEFAULT

       WRITE(*,*) 'ERROR: Unsupported Quadrangle type'
       STOP

    END SELECT
    
  END FUNCTION rd_normal_fun
  !=========================  

  !============================================
  SUBROUTINE gradient_trace_sub(e, if, p_D_phi)
  !============================================
  !
  ! Compute the trace of the gradients of all shape
  ! functions on the face if
  !
    IMPLICIT NONE

    CLASS(quadrangle)                           :: e
    INTEGER,                        INTENT(IN)  :: if
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: p_D_phi
    !-----------------------------------------------------

    REAL(kind=8), DIMENSION(:,:), ALLOCATABLE :: xi_f
    INTEGER :: N_points_ele, N_points_face, k, i
    !-----------------------------------------------------

    N_points_ele  = e%N_points
    N_points_face = SIZE(p_D_phi, 2)

    ALLOCATE( xi_f(N_points_face, 2) )

    xi_f = e%face_trace(if)    

    DO k = 1, N_points_face      

       DO i = 1, N_points_ele

          p_D_phi(:, k, i) = e%gradient( i, xi_f(k,:) )

       ENDDO

    ENDDO
    
    DEALLOCATE(xi_f)    

  END SUBROUTINE gradient_trace_sub
  !================================

  !==========================================
  FUNCTION face_trace_fun(e, if) RESULT(xi_f)
  !==========================================
  !
  ! Compute the coordinates on the face if
  !
    IMPLICIT NONE

    CLASS(quadrangle)        :: e
    INTEGER,      INTENT(IN) :: if
    
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: xi_f
    !---------------------------------------------
     
    SELECT CASE(e%Type)
       
    CASE(QUA_Q1)

       ALLOCATE( xi_f(2, 2) )

       xi_f = face_trace_QUA_Q1(e, if)
       
    CASE(QUA_Q2)

       ALLOCATE( xi_f(3, 2) )

       xi_f = face_trace_QUA_Q2(e, if)
       
    CASE DEFAULT

       WRITE(*,*) 'Unknown Quadrangle type for face trace'
       WRITE(*,*) 'STOP'

    END SELECT 

  END FUNCTION face_trace_fun
  !==========================

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%% SPECIFIC FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%    QUADRANGLE Q1   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !============================================================
  SUBROUTINE init_faces_QUA_Q1(e, Nodes, Coords, NU_seg, n_ele)
  !============================================================

    IMPLICIT NONE

    CLASS(quadrangle)                        :: e
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

    INTEGER :: id, if, i, j, k, istat
    !-------------------------------------------------

    e%N_faces  = 4 ! # Faces

    ALLOCATE( e%faces(e%N_faces) )

    ! Dofs on the faces
    ALLOCATE( VV(2), loc(2) )

    ! Coordinates on the faces
    ALLOCATE( RR(SIZE(Coords, 1), 2) )

    DO if = 1, e%N_faces

       i = MOD(if, MIN(if+1, e%N_faces)) + 1
       j = MOD(i,  MIN(if+2, e%N_faces)) + 1

       loc = (/ i, j /)

       VV = Nodes(loc)
       RR = Coords(:, loc)

       ALLOCATE(seg, STAT=istat)          
       IF(istat /=0) THEN
          WRITE(*,*) 'ERROR: failed quadrangle allocation'
       ENDIF

       CALL seg%initialize( loc, VV, RR, NU_seg(if), n_ele(if) )

       CALL e%gradient_trace(if, p_D_phi)
       
       CALL seg%face_quadrature(p_D_phi)
          
       e%faces(if)%f => seg
          
    ENDDO

    DEALLOCATE( VV, RR, loc, p_D_phi )

  END SUBROUTINE init_faces_QUA_Q1
  !===============================

  !=====================================================
  FUNCTION basis_function_QUA_Q1(e, i, xi) RESULT(psi_i)
  !=====================================================

    IMPLICIT NONE

    CLASS(quadrangle)                      :: e
    INTEGER,                    INTENT(IN) :: i
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8) :: psi_i
    !-------------------------------------------

    SELECT CASE (i)
       
    CASE(1)

       psi_i = (1.d0 - xi(1))*(1.d0 - xi(2))

    CASE(2)

       psi_i = xi(1)*(1.d0 - xi(2))

    CASE(3)
       psi_i = xi(1)*xi(2)
       
    CASE(4)       
       
       psi_i = (1.d0 - xi(1))*xi(2)

    CASE DEFAULT

       WRITE(*,*) 'ERROR: not supported Dof in Quadrangle  baisis function'
       STOP
       
    END SELECT

  END FUNCTION basis_function_QUA_Q1
  !=================================

  !=====================================================
  FUNCTION gradient_ref_QUA_Q1(e, i, xi) RESULT(D_psi_i)
  !=====================================================

    IMPLICIT NONE

    CLASS(quadrangle)                      :: e
    INTEGER,                    INTENT(IN) :: i
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8), DIMENSION(2) :: D_psi_i
    !-------------------------------------------

    SELECT CASE(i)

    CASE(1)

       D_psi_i = (/ -1.d0 + xi(2), &
                    -1.d0 + xi(1) /)
       
    CASE(2)

       D_psi_i = (/ -1.d0 + xi(2), &
                           -xi(1) /)
       
    CASE(3)

        D_psi_i = (/ xi(2), &
                     xi(1) /)

    CASE(4)

        D_psi_i = (/       -xi(2), &
                     1.d0 - xi(1) /)

    CASE DEFAULT

       WRITE(*,*) 'ERROR: not supported Dof in Quadrangle gradient'
       STOP

    END SELECT
    
  END FUNCTION gradient_ref_QUA_Q1
  !===============================  
  
  !===========================================
  FUNCTION rd_normal_QUA_Q1(e, k) RESULT(nn_k)
  !===========================================

    IMPLICIT NONE

    CLASS(quadrangle)   :: e
    INTEGER, INTENT(IN) :: k

    REAL(KIND=8), DIMENSION(e%N_dim) :: nn_k
    !---------------------------------------

    INTEGER :: i, j
    !---------------------------------------

    SELECT CASE(k)
    CASE(1)
       i = 2; j = 4
    CASE(2)
       i = 3; j = 1
    CASE(3)
       i = 4; j = 2
    CASE(4)
       i = 1; j = 3
    CASE DEFAULT
       WRITE(*,*) 'ERROR: wrong nodes for rd_normal'
       STOP       
    END SELECT    

    nn_k(1) = e%Coords(2,i) - e%Coords(2,j)
    nn_k(2) = e%Coords(1,j) - e%Coords(1,i)

  END FUNCTION rd_normal_QUA_Q1
  !============================

  !=============================================
  FUNCTION face_trace_QUA_Q1(e, if) RESULT(xi_f)
  !=============================================

    IMPLICIT NONE

    CLASS(quadrangle)        :: e
    INTEGER,      INTENT(IN) :: if

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: xi_f
    !------------------------------------------------

    ALLOCATE( xi_f(2, 2) )

    SELECT CASE(if)

    CASE(1)

       xi_f(1, :) = (/ 1.d0, 0.d0 /)
       xi_f(2, :) = (/ 1.d0, 1.d0 /)

    CASE(2)

       xi_f(1, :) = (/ 1.d0, 1.d0 /)
       xi_f(2, :) = (/ 0.d0, 1.d0 /)
       
    CASE(3)
  
       xi_f(1, :) = (/ 0.d0, 1.d0 /)
       xi_f(2, :) = (/ 0.d0, 0.d0 /)

    CASE(4)
  
       xi_f(1, :) = (/ 0.d0, 0.d0 /)
       xi_f(2, :) = (/ 1.d0, 0.d0 /)

    CASE DEFAULT

       WRITE(*,*) 'ERROR: wrong Quadrangle face in face trace Q1'
       STOP

    END SELECT    
    
  END FUNCTION face_trace_QUA_Q1
  !=============================

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%% SPECIFIC FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%    QUADRANGLE Q2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !============================================================
  SUBROUTINE init_faces_QUA_Q2(e, Nodes, Coords, NU_seg, n_ele)
  !============================================================

    IMPLICIT NONE

    CLASS(quadrangle)                        :: e
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
    
    INTEGER :: id, if, i, j, k, istat
    !-------------------------------------------------

    e%N_faces  = 4 ! # Faces

    ALLOCATE( e%faces(e%N_faces) )

    ! Dofs on the faces
    ALLOCATE( VV(3), loc(3) )

    ! Coordinates on the faces
    ALLOCATE( RR(SIZE(Coords, 1), 3) )

    DO if = 1, e%N_faces

       i = MOD(if, MIN(if+1, e%N_faces)) + 1
       j = MOD(i,  MIN(if+2, e%N_faces)) + 1

       k = i + 4

       loc = (/ i, j, k /)

       VV = Nodes(loc)
       RR = Coords(:, loc)

       ALLOCATE(seg, STAT=istat)          
       IF(istat /=0) THEN
          WRITE(*,*) 'ERROR: failed quadrangle allocation'
       ENDIF

       CALL seg%initialize( loc, VV, RR, NU_seg(if), n_ele(if) )

       CALL e%gradient_trace(if, p_D_phi)
       
       CALL seg%face_quadrature(p_D_phi)
          
       e%faces(if)%f => seg
          
    ENDDO

    DEALLOCATE( VV, RR, loc )

  END SUBROUTINE init_faces_QUA_Q2
  !===============================

  !=====================================================
  FUNCTION basis_function_QUA_Q2(e, i, xi) RESULT(psi_i)
  !=====================================================

    IMPLICIT NONE

    CLASS(quadrangle)                      :: e
    INTEGER,                    INTENT(IN) :: i
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8) :: psi_i
    !-------------------------------------------

    SELECT CASE (i)
       
    CASE(1)

       psi_i = (1.d0 - xi(1))*(1.d0 - 2.d0*xi(1)) * &
               (1.d0 - xi(2))*(1.d0 - 2.d0*xi(2))

    CASE(2)

       psi_i =      xi(1) *(2.d0*xi(1) - 1.d0) * &
               (1 - xi(2))*(1 - 2.d0*xi(2))

    CASE(3)
       
      psi_i = xi(1) *(2.d0*xi(1) - 1.d0) * &
              xi(2) *(2.d0*xi(2) - 1.d0)

    CASE(4)       
       
      psi_i = (1.d0 - xi(1))*(1.d0 - 2.d0*xi(1)) * &
                      xi(2) *(2.d0*xi(2) - 1.d0)

   CASE(5)

      psi_i = 4.d0*(1.d0 - xi(1))*xi(1) * &
              (1.d0 - xi(2))*(1.d0 - 2.d0*xi(2))

   CASE(6)

      psi_i = xi(1) *(2.d0*xi(1) - 1.d0) * &
              4.d0*(1.d0 - xi(2))*xi(2)

   CASE(7)

      psi_i = 4.d0*(1.d0 - xi(1))*xi(1) * &
              xi(2) *(2.d0*xi(2) - 1.d0)

   CASE(8)

      psi_i = (1.d0 - xi(1))*(1.d0 - 2.d0*xi(1)) * &
               4.d0*(1.d0 - xi(2))*xi(2)

   CASE(9)

      psi_i = 4.d0*(1.d0 - xi(1))*xi(1) * &
              4.d0*(1.d0 - xi(2))*xi(2)
            
   CASE DEFAULT
      
       WRITE(*,*) 'ERROR: not supported Dof in Quadrangle  baisis function'
       STOP
       
    END SELECT

  END FUNCTION basis_function_QUA_Q2
  !=================================

  !=====================================================
  FUNCTION gradient_ref_QUA_Q2(e, i, xi) RESULT(D_psi_i)
  !=====================================================

    IMPLICIT NONE

    CLASS(quadrangle)                      :: e
    INTEGER,                    INTENT(IN) :: i
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8), DIMENSION(2) :: D_psi_i
    !-------------------------------------------

    REAL(KIND=8) :: x, y
    !-------------------------------------------

    x = xi(1); y = xi(2)

    SELECT CASE(i)

    CASE(1)

       D_psi_i = (/ -3.d0 + 9.d0*y - 6.d0*y*y + 4.d0*x  - 12.d0*x*y + 8.d0*x*y*y, &
                    -3.d0 + 4.d0*y + 9.d0*x   -12.d0*x*y - 6.d0*x*x + 8.d0*x*x*y /)
       
    CASE(2)

       D_psi_i = (/ 4.d0*x - 12.d0*x*y + 8.d0*x*y*y - 1.d0 + 3.d0*y -2.d0*y*y, &
                   -6.d0*x*x + 8.d0*x*x*y + 3.d0*x -4.d0*x*y /)
       
    CASE(3)

        D_psi_i = (/ 8.d0*x*y*y - 4.d0*x*y - 2.d0*y*y + y, &
                     8.d0*x*x*y - 2.d0*x*x - 4.d0*x*y + x /)

    CASE(4)

        D_psi_i = (/ -6.d0*y*y + 3.d0*y + 8.d0*x*y*y - 4.d0*x*y, &
                      4.d0*y - 1.d0 - 12.d0*x*y + 3.d0*x + 8.d0*x*x*y - 2.d0*x*x /)

     CASE(5)

        D_psi_i = (/ -8.d0*x + 24.d0*x*y - 16.d0*x*y*y + 4.d0 - 12.d0*y + 8.d0*y*y, &
                    -12.d0*x + 16.d0*x*y + 12.d0*x*x - 16.d0*x*x*y /)
     CASE(6)

        D_psi_i = (/ 16.d0*x*y - 16.d0*x*y*y - 4.d0*y + 4.d0*y*y, &
                    -16.d0*x*x*y + 8.d0*x*y + 8.d0*x*x - 4.d0*x /)
        
     CASE(7)

        D_psi_i = (/ -16.d0*x*y*y + 8.d0*x*y + 8.d0*y*y - 4.d0*y, &
                      16.d0*x*y - 4.d0*x - 16.d0*x*x*y + 4.d0*x*x /)

     CASE(8)

        D_psi_i = (/ -12.d0*y + 12.d0*y*y + 16.d0*x*y - 16.d0*x*y*y, &
                     -8.d0*y + 24.d0*x*y - 16.d0*x*x*y + 4.d0 - 12.d0*x + 8.d0*x*x /)

     CASE(9)

        D_psi_i = (/ -32.d0*x*y + 32.d0*x*y*y + 16.d0*y - 16.d0*y*y, &
                     -32.d0*x*y + 32.d0*x*x*y + 16.d0*x - 16.d0*x*x /)
        
    CASE DEFAULT

       WRITE(*,*) 'ERROR: not supported Dof in Quadrangle gradient'
       STOP

    END SELECT
    
  END FUNCTION gradient_ref_QUA_Q2
  !===============================

  !===========================================
  FUNCTION rd_normal_QUA_Q2(e, k) RESULT(nn_k)
  !===========================================

    IMPLICIT NONE

    CLASS(quadrangle)   :: e
    INTEGER, INTENT(IN) :: k

    REAL(KIND=8), DIMENSION(e%N_dim) :: nn_k
    !---------------------------------------

    INTEGER :: i, j
    !---------------------------------------

    SELECT CASE(k)
    CASE(1)
       i = 2; j = 4
    CASE(2)
       i = 3; j = 1
    CASE(3)
       i = 4; j = 2
    CASE(4)
       i = 1; j = 3
    CASE(5)
       i = 2; j = 1
    CASE(6)
       i = 3; j = 2
    CASE(7)
       i = 4; j = 3
    CASE(8)
       i = 1; j = 4
    CASE(9)
       i = 1; j = 1        
    CASE DEFAULT
       WRITE(*,*) 'ERROR: wrong nodes for rd_normal'
       STOP       
    END SELECT    

    nn_k(1) = e%Coords(2,i) - e%Coords(2,j)
    nn_k(2) = e%Coords(1,j) - e%Coords(1,i)

  END FUNCTION rd_normal_QUA_Q2
  !============================

  !=============================================
  FUNCTION face_trace_QUA_Q2(e, if) RESULT(xi_f)
  !=============================================

    IMPLICIT NONE

    CLASS(quadrangle)        :: e
    INTEGER,      INTENT(IN) :: if

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: xi_f
    !----------------------------------------------

    ALLOCATE( xi_f(3, 2) )

    SELECT CASE(if)

    CASE(1)

       xi_f(1, :) = (/ 1.d0,  0.d0  /)
       xi_f(2, :) = (/ 1.d0,  1.d0  /)
       xi_f(3, :) = (/ 1.0d0, 0.5d0 /)

    CASE(2)

       xi_f(1, :) = (/ 1.d0,  1.d0 /)
       xi_f(2, :) = (/ 0.d0,  1.d0 /)
       xi_f(3, :) = (/ 0.5d0, 1.d0 /)

    CASE(3)
  
       xi_f(1, :) = (/ 0.d0,  1.d0  /)
       xi_f(2, :) = (/ 0.d0,  0.d0  /)
       xi_f(3, :) = (/ 0.0d0, 0.5d0 /)

     CASE(4)

       xi_f(1, :) = (/ 0.d0,  0.d0 /)
       xi_f(2, :) = (/ 1.d0,  0.d0 /)
       xi_f(3, :) = (/ 0.5d0, 0.d0 /)
           
    CASE DEFAULT

       WRITE(*,*) 'ERROR: wrong Quadrangle face in face trace Q2'
       STOP

    END SELECT    
    
  END FUNCTION face_trace_QUA_Q2
  !=============================

END MODULE Quadrangle_Class

