MODULE Lin_Algebra

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Cross_Product, inverse, determinant

CONTAINS

  !====================================
  FUNCTION Cross_Product(AA) RESULT(Vn)
  !====================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: AA

    REAL(KIND=8), DIMENSION(SIZE(AA,2)) :: Vn
    !----------------------------------------------

    REAL(KIND=8), DIMENSION(3)   :: a, b, c
    REAL(KIND=8), DIMENSION(3,3) :: a_cross
    !----------------------------------------------

    IF( (SIZE(AA,2) /= SIZE(AA,1)+1) .AND. &
        (SIZE(AA,2) /= 2 .OR. SIZE(AA,2) /= 3) ) THEN
        WRITE(*,*) 'ERROR: inconsistent sizes for the cross product'   
        STOP
    ENDIF

    a = 0.d0
    a(1:SIZE(AA,2)) = AA(1, :)

    b = 0.d0
    IF(SIZE(AA,1) == 1) THEN
       b(3) = 1.d0
    ELSE
       b(:) = AA(2, :)
    ENDIF
        
    a_cross(1, :) = (/  0.d0, -a(3),  a(2) /)
    a_cross(2, :) = (/  a(3),  0.d0, -a(1) /)
    a_cross(3, :) = (/ -a(2),  a(1),  0.d0 /)

    c = MATMUL(a_cross, b)   

    Vn = c(1:SIZE(AA,2))
    
  END FUNCTION Cross_Product
  !=========================

  !==============================
  FUNCTION inverse (A) RESULT (B)
  !==============================

      IMPLICIT NONE
     
      REAL(KIND=8), DIMENSION(:,:),     INTENT(IN)  ::  A
      REAL(KIND=8), DIMENSION(SIZE(A,1),SIZE(A,2))  ::  B

      REAL(KIND=8)  ::  det


      IF ( SIZE(A,1) .NE. SIZE(A,2) ) THEN
         WRITE(*,*) ' Matrix A is not square.' 
         WRITE(*,*) ' in FUNCTION inverse, MODULE lin_algebra' 
         WRITE(*,*) ' STOP. ' 
         STOP
      ENDIF

      det = determinant(A)

      IF ( det == 0 ) THEN
         WRITE(*,*) ' Matrix A is singular.' 
         WRITE(*,*) ' in FUNCTION inverse, MODULE lin_algebra' 
         WRITE(*,*) ' STOP. ' 
         STOP
      ENDIF

      SELECT CASE ( SIZE(A,1) )

      CASE (2) 
      ! ------
         B(1,1) =   A(2,2);  B(1,2) = - A(1,2) 
         B(2,1) = - A(2,1);  B(2,2) =   A(1,1) 
         B = B/det

      CASE (3)
      ! ------
         B(1,1) =   A(2,2)*A(3,3) - A(3,2)*A(2,3)
         B(1,2) = - A(1,2)*A(3,3) + A(3,2)*A(1,3)
         B(1,3) =   A(1,2)*A(2,3) - A(2,2)*A(1,3)

         B(2,1) = - A(2,1)*A(3,3) + A(3,1)*A(2,3)
         B(2,2) =   A(1,1)*A(3,3) - A(3,1)*A(1,3)
         B(2,3) = - A(1,1)*A(2,3) + A(2,1)*A(1,3)

         B(3,1) =   A(2,1)*A(3,2) - A(3,1)*A(2,2)
         B(3,2) = - A(1,1)*A(3,2) + A(3,1)*A(1,2)
         B(3,3) =   A(1,1)*A(2,2) - A(2,1)*A(1,2)

         B = B/det

      CASE DEFAULT
         WRITE(*,*) ' Matrix of size ', SIZE(A,1), ' not implemented' 
         WRITE(*,*) ' in FUNCTION inverse, MODULE lin_algebra' 
         WRITE(*,*) ' STOP. ' 
         STOP
         ! CALL ......

      END SELECT

   END FUNCTION inverse
   !===================

   !===================================
   FUNCTION determinant(A) RESULT (det)
   !=================================== 

      IMPLICIT NONE
     
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  ::  A
      REAL(KIND=8)  ::  det 


      IF ( SIZE(A,1) .NE. SIZE(A,2) ) THEN
         WRITE(*,*) ' Matrix A is not square.' 
         WRITE(*,*) ' in FUNCTION determinant, MODULE lin_algebra' 
         WRITE(*,*) ' STOP. ' 
         STOP
      ENDIF

      SELECT CASE ( SIZE(A,1) )

      CASE (2)
         det = A(1,1)*A(2,2) - A(1,2)*A(2,1)

      CASE (3)
         det = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3))  &
             - A(1,2)*(A(2,1)*A(3,3) - A(3,1)*A(2,3))  &
             + A(1,3)*(A(2,1)*A(3,2) - A(3,1)*A(2,2))

      CASE DEFAULT 
         WRITE(*,*) ' Matrix of size ', SIZE(A,1), ' not implemented' 
         WRITE(*,*) ' in FUNCTION determinant, MODULE lin_algebra' 
         WRITE(*,*) ' STOP. ' 
         STOP

      END SELECT

   END FUNCTION determinant 
   !=======================

END MODULE Lin_Algebra
