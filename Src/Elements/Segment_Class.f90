MODULE Segment_Class

  USE Face_Class
  USE Lin_Algebra

  IMPLICIT NONE
  
  !=========================================
  TYPE, PUBLIC, EXTENDS(face) :: segment

     ! < >

   CONTAINS

     PROCEDURE, PUBLIC :: initialize      => initialize_sub
     PROCEDURE, PUBLIC :: face_quadrature => face_quadrature_sub
     

  END TYPE segment
  !=========================================

  PRIVATE :: initialize_sub

  PRIVATE :: face_quadrature_sub
  PRIVATE :: init_quadrature_SEG_P1
  PRIVATE :: init_quadrature_SEG_P2

  PRIVATE :: normal
  
  PRIVATE :: basis_function
  PRIVATE :: basis_function_SEG_P1
  PRIVATE :: basis_function_SEG_P2

  PRIVATE :: gradient_ref
  PRIVATE :: gradient_ref_SEG_P1
  PRIVATE :: gradient_ref_SEG_P2

CONTAINS

  !==============================================================
  SUBROUTINE initialize_sub(e, loc, Nodes, Coords, NU_seg, n_ele)
  !==============================================================

    IMPLICIT NONE

    CLASS(segment)                           :: e
    INTEGER,      DIMENSION(:),   INTENT(IN) :: loc
    INTEGER,      DIMENSION(:),   INTENT(IN) :: Nodes
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Coords
    INTEGER,                      INTENT(IN) :: NU_seg
    INTEGER,                      INTENT(IN) :: n_ele
    !-------------------------------------------------

    INTEGER :: id
    !-------------------------------------------------

    SELECT CASE( SIZE(Nodes) )

    CASE(2)

       e%Type     = SEG_P1
       e%N_verts  = 2
       e%N_points = 2
              
    CASE(3)

       e%Type     = SEG_P2
       e%N_verts  = 2
       e%N_points = 3

     CASE DEFAULT

        WRITE(*,*) 'ERROR: Unsupported Segment type'
        STOP

     END SELECT

     !---------------------------
     ! Attach data to the element
     !---------------------------------------------------
     ALLOCATE( e%Coords(SIZE(Coords,1), SIZE(Coords,2)) )
     ALLOCATE( e%NU(SIZE(Nodes)) )
     ALLOCATE( e%l_NU(SIZE(loc)) )
          
     DO id = 1, SIZE(Coords,1)
        e%Coords(id, :) = Coords(id, :)
     ENDDO

     e%NU = Nodes

     e%l_NU = loc

     e%g_seg = NU_seg

     e%c_ele = n_ele
     
  END SUBROUTINE initialize_sub
  !============================

  !=========================================
  SUBROUTINE face_quadrature_sub(e, p_D_phi)
  !=========================================
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

    CLASS(segment)                             :: e
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN) :: p_D_phi
    !----------------------------------------------------

    REAL(KIND=8) :: Jac    
    INTEGER :: j, k, l
    !----------------------------------------------------

    SELECT CASE(e%Type)

    CASE(SEG_P1)

       CALL init_quadrature_SEG_P1(e, p_D_phi)
       
    CASE(SEG_P2)

       CALL init_quadrature_SEG_P2(e, p_D_phi)
       
    CASE DEFAULT

       WRITE(*,*) 'Unknown Segment type for quadrature'
       WRITE(*,*) 'STOP'

    END SELECT

    !-------------------------------------
    ! Attach data to the quadrature points
    !---------------------------------------------------
    DO j = 1, e%N_quad

       DO k = 1, e%N_points

          e%phi_q(k, j) = basis_function( e, k, e%x_q(j, 1) )

          e%xx_q(:, j) = e%xx_q(:, j) + &
                         basis_function( e, k, e%x_q(j, 1) ) * e%Coords(:, k)

       ENDDO

       e%n_q(:, j) = normal(e, e%x_q(j, 1) )

       Jac = SQRT( SUM(e%n_q(:, j)**2) )

       e%n_q(:, j) = e%n_q(:, j)/Jac

       e%w_q(j) = e%w_q(j) * Jac


       DO k = 1, SIZE(p_D_phi,3)

          DO l = 1, e%N_points
             
             e%p_Dphi_1_q(:, k, j) = e%p_Dphi_1_q(:, k, j) + &
                                     p_D_phi(:, l, k)*e%phi_q(l, j)
 
          ENDDO          
          
       ENDDO

    ENDDO

    ! Area of the element
    ! warning: straight segment only!!!
    e%area = SUM(e%w_q)
    
  END SUBROUTINE face_quadrature_sub
  !=================================


!*******************************************************************************
!*******************************************************************************
!                         COMPLEMENTARY FUNCTIONS                              !
!*******************************************************************************
!*******************************************************************************

  !==============================================
  FUNCTION basis_function(e, i, xi) RESULT(psi_i)
  !==============================================

    IMPLICIT NONE

    CLASS(segment)            :: e
    INTEGER,       INTENT(IN) :: i
    REAL(KIND=8),  INTENT(IN) :: xi

    REAL(KIND=8) :: psi_i
    !-------------------------------
      
    SELECT CASE(e%Type)
       
    CASE(SEG_P1)

       psi_i = basis_function_SEG_P1(e, i, xi)
       
    CASE(SEG_P2)

       psi_i = basis_function_SEG_P2(e, i, xi)
       
    CASE DEFAULT

       WRITE(*,*) 'Unknown Segment type for basis functions'
       WRITE(*,*) 'STOP'

    END SELECT

  END FUNCTION basis_function
  !==========================

  !==============================================
  FUNCTION gradient_ref(e, i, xi) RESULT(D_psi_i)
  !==============================================

    IMPLICIT NONE

    CLASS(segment)            :: e
    INTEGER,       INTENT(IN) :: i
    REAL(KIND=8),  INTENT(IN) :: xi

    REAL(KIND=8), DIMENSION(e%N_dim) :: D_psi_i
    !-------------------------------------------
      
    SELECT CASE(e%Type)
       
    CASE(SEG_P1)

       D_psi_i = gradient_ref_SEG_P1(e, i, xi)
       
    CASE(SEG_P2)

       D_psi_i = gradient_ref_SEG_P2(e, i, xi)
       
    CASE DEFAULT

       WRITE(*,*) 'Unknown Segment type for gradients'
       WRITE(*,*) 'STOP'

    END SELECT

  END FUNCTION gradient_ref
  !========================

  !=================================
  FUNCTION normal(e, xi) RESULT(n_i)
  !=================================
  !
  ! Compute the normal versor to the face 
  ! on the node with coordiantes xi
  !
    IMPLICIT NONE

    CLASS(segment)            :: e
    REAL(KIND=8),  INTENT(IN) :: xi

    REAL(KIND=8), DIMENSION(e%N_dim) :: n_i
    !---------------------------------------------

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: d
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Jb

    INTEGER :: i, k, l
    !---------------------------------------------

    ALLOCATE( d(e%N_dim, e%N_points) )

    DO i = 1, e%N_points
       d(:, i) = gradient_ref(e, i, xi)
    ENDDO    

    ALLOCATE( Jb(e%N_dim-1, e%N_dim) )

    DO k = 1, e%N_dim

       DO l = 1, e%N_dim - 1

          Jb(l, k) = SUM( e%Coords(k, :) * d(l, :) )

       ENDDO

    ENDDO

    n_i = Cross_product(Jb)

    DEALLOCATE( d, Jb )

  END FUNCTION normal
  !==================

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%% SPECIFIC FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%     SEGMENT P1     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  !=====================================================
  FUNCTION basis_function_SEG_P1(e, i, xi) RESULT(psi_i)
  !=====================================================

    IMPLICIT NONE

    CLASS(segment)            :: e
    INTEGER,       INTENT(IN) :: i
    REAL(KIND=8),  INTENT(IN) :: xi

    REAL(KIND=8) :: psi_i
    !------------------------------

    SELECT CASE(i)

    CASE(1)

       psi_i = 1.d0 - xi

    CASE(2)

       psi_i = xi

    CASE DEFAULT

       WRITE(*,*) 'ERROR: not supported Dof in Segment basis function'
       STOP

    END SELECT    

  END FUNCTION basis_function_SEG_P1
  !=================================

  !=====================================================
  FUNCTION gradient_ref_SEG_P1(e, i, xi) RESULT(D_psi_i)
  !=====================================================

    IMPLICIT NONE

    CLASS(segment)            :: e
    INTEGER,       INTENT(IN) :: i
    REAL(KIND=8),  INTENT(IN) :: xi

    REAL(KIND=8) :: D_psi_i
    !------------------------------

    SELECT CASE(i)

    CASE(1)

       D_psi_i = -1.d0

    CASE(2)

       D_psi_i = 1.d0

    CASE DEFAULT

       WRITE(*,*) 'ERROR: not supported Dof in Segment gradient'
       STOP

    END SELECT    

  END FUNCTION gradient_ref_SEG_P1
  !===============================
  
  !============================================
  SUBROUTINE init_quadrature_SEG_P1(e, p_D_phi)
  !============================================

    IMPLICIT NONE

    CLASS(segment)                             :: e
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN) :: p_D_phi
    !----------------------------------------------------

    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: x_q
    REAL(KIND=8) :: Jac
    INTEGER :: j, k, l
    !----------------------------------------------------

    e%N_quad = 2

    ALLOCATE(   e%n_q(e%N_dim,    e%N_quad) )
    ALLOCATE( e%phi_q(e%N_points, e%N_quad) )

    ALLOCATE( e%w_q(e%N_quad   ) )
    ALLOCATE( e%x_q(e%N_quad, 1) )

    ALLOCATE( e%p_Dphi_1_q(e%N_dim, SIZE(p_D_phi,3), e%N_quad) )
    ALLOCATE( e%p_Du_1_q(e%N_dim, e%N_quad) )

    ALLOCATE( e%xx_q(e%N_dim, e%N_quad) )

    !-------------------
    ! Quadrature formula
    !-----------------------------
    e%x_q(1, 1) = 0.5d0*( 1.d0 - SQRT(1.d0/3.d0) )
    e%x_q(2, 1) = 0.5d0*( 1.d0 + SQRT(1.d0/3.d0) )

    e%w_q(1) = 0.5d0
    e%w_q(2) = 0.5d0
    !-----------------------------

    e%p_Dphi_1_q = 0.d0
    e%p_Dphi_2_q = 0.d0

    e%xx_q = 0.d0
    
  END SUBROUTINE init_quadrature_SEG_P1
  !====================================

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%% SPECIFIC FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%     SEGMENT P2     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !=====================================================
  FUNCTION basis_function_SEG_P2(e, i, xi) RESULT(psi_i)
  !=====================================================

    IMPLICIT NONE

    CLASS(segment)            :: e
    INTEGER,       INTENT(IN) :: i
    REAL(KIND=8),  INTENT(IN) :: xi

    REAL(KIND=8) :: psi_i
    !------------------------------

    SELECT CASE(i)

    CASE(1)

       psi_i = (1.d0 - xi)*(1.d0 - 2.d0*xi)

    CASE(2)

       psi_i = xi*(2.d0*xi - 1.d0)

    CASE(3)       

       psi_i = 4.d0*(1.d0 - xi)*xi 

    CASE DEFAULT

       WRITE(*,*) 'ERROR: not supported Dof in Segment baisis function'
       STOP

    END SELECT    

  END FUNCTION basis_function_SEG_P2
  !=================================

  !=====================================================
  FUNCTION gradient_ref_SEG_P2(e, i, xi) RESULT(D_psi_i)
  !=====================================================

    IMPLICIT NONE

    CLASS(segment)            :: e
    INTEGER,       INTENT(IN) :: i
    REAL(KIND=8),  INTENT(IN) :: xi

    REAL(KIND=8) :: D_psi_i
    !------------------------------

    SELECT CASE(i)

    CASE(1)

       D_psi_i = 4.d0*xi - 3.d0

    CASE(2)

       D_psi_i = 4.d0*xi - 1.d0

    CASE(3)

       D_psi_i = 4.d0 - 8.d0*xi

    CASE DEFAULT

       WRITE(*,*) 'ERROR: not supported Dof in Segment gradient'
       STOP

    END SELECT    

  END FUNCTION gradient_ref_SEG_P2
  !================================  

  !============================================
  SUBROUTINE init_quadrature_SEG_P2(e, p_D_phi)
  !============================================

    IMPLICIT NONE

    CLASS(segment)                             :: e
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN) :: p_D_phi
    !----------------------------------------------------

    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: x_q
    REAL(KIND=8) :: Jac
    INTEGER :: j, k, l
    !----------------------------------------------------

    e%N_quad = 3

    ALLOCATE(   e%n_q(e%N_dim,    e%N_quad) )
    ALLOCATE( e%phi_q(e%N_points, e%N_quad) )

    ALLOCATE( e%w_q(e%N_quad   ) )
    ALLOCATE( e%x_q(e%N_quad, 1) )
    
    ALLOCATE( e%p_Dphi_1_q(e%N_dim, SIZE(p_D_phi,3), e%N_quad) )
    ALLOCATE( e%p_Du_1_q(e%N_dim, e%N_quad) )

    ALLOCATE( e%xx_q(e%N_dim, e%N_quad) )

    !-------------------
    ! Quadrature formula
    !----------------------------------------------
    e%x_q(1, 1) = 0.5d0 * (1.d0 - sqrt(3.d0/5.d0) )
    e%x_q(2, 1) = 0.5d0 * (1.d0 + sqrt(3.d0/5.d0) )
    e%x_q(3, 1) = 0.5d0

    e%w_q(1) = 5.d0/18.d0
    e%w_q(2) = 5.d0/18.d0
    e%w_q(3) = 8.d0/18.d0
    !-----------------------------------------------

    e%p_Dphi_1_q = 0.d0
    e%p_Dphi_2_q = 0.d0

    e%xx_q = 0.d0
    
  END SUBROUTINE init_quadrature_SEG_P2
  !====================================
 
END MODULE Segment_Class

