MODULE Post_pro

  USE Element_Class

  USE init_problem,  ONLY : pb_name, order, pb_type, visc, is_visc

  USE geometry,      ONLY : N_dim, N_elements, N_dofs, elements

  USE models,        ONLY : exact_solution, exact_grad

  USE Gradient_Reconstruction

  USE Quadrature_rules,  ONLY : Int_d

  IMPLICIT NONE
  PRIVATE

  !==========================================
  TYPE :: type_verts
     INTEGER, DIMENSION(:), POINTER :: verts
  END TYPE type_verts
  !==========================================

  PUBLIC :: plot_procedure, compute_error
  !==========================================

CONTAINS

  !============================
  SUBROUTINE plot_procedure(uu)
  !============================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: uu
    !------------------------------------------------------

    INTEGER, PARAMETER               :: TRI = 1, QUA = 2
    INTEGER, PARAMETER, DIMENSION(2) :: Nverts = (/ 3, 4 /)
    INTEGER,            DIMENSION(2) :: N_p_ele
    
    REAL(KIND=8), DIMENSION(N_dim, N_dofs) :: coord
    INTEGER,      DIMENSION(N_dofs)        :: Glo2Loc

    REAL :: uu_ex

    INTEGER :: N_tri, N_qua, ns
    INTEGER :: UNIT, ierror, type, jt, j, k
    !-------------------------------------------------------

    N_tri = 0
    N_qua = 0
    ns    = 0

    ! Number of elements (plus sub-element in Pk, k>=2)
    DO jt = 1, N_elements
       SELECT CASE(elements(jt)%p%Type)
       CASE(TRI_P1)
          N_tri = N_tri + 1
       CASE(TRI_P2)
          N_tri = N_tri + 4
       CASE(QUA_Q1)
          N_qua = N_qua + 1
       CASE(QUA_Q2)
          N_qua = N_qua + 4
       CASE DEFAULT
          WRITE(*,*) 'WRONG element. STOP!'
          STOP
       END SELECT
    ENDDO

    N_p_ele = (/ N_tri, N_qua /)

    UNIT = 1
    OPEN(UNIT, FILE = 'plot_'//TRIM(ADJUSTL(pb_name))//'.plt', & 
         ACTION = 'WRITE', IOSTAT = ierror)

    IF (ierror /= 0) THEN
       WRITE(*,*) 'ERROR: Impossible to write file .plt' 
       WRITE(*,*) 'STOP!'
       STOP
    ENDIF

    ! Expansion of nodes coordinates
    DO jt = 1, N_elements
       DO j = 1, elements(jt)%p%N_points
          coord(:, elements(jt)%p%Nu(j)) = elements(jt)%p%Coords(:, j)
       ENDDO
    ENDDO

   ! Tecplot's header
   WRITE(UNIT, 100)
   WRITE(UNIT, 110)

   Glo2Loc = 0 

   DO type = TRI, QUA

      IF (N_p_ele(type) > 0) THEN
         
         ! Nodes renumeration
         DO jt = 1, N_elements
              IF (elements(jt)%p%N_verts == Nverts(type)) THEN
                 DO j = 1, elements(jt)%p%N_points
                    IF ( Glo2Loc(elements(jt)%p%NU(j) ) == 0 ) THEN
                       ns = ns + 1
                       Glo2Loc(elements(jt)%p%NU(j)) = ns
                    ENDIF
                 ENDDO
              ENDIF
           ENDDO
           
           ! TECPLOT type headers
           SELECT CASE(type)
           CASE(TRI)
              WRITE(UNIT, 120) ns, N_p_ele(TRI)
           CASE(QUA)
              WRITE(UNIT, 121) ns, N_p_ele(QUA)
           END SELECT

           ns = 0

           !------------
           ! Coordiantes
           ! Solution
           !-------------------------------------           
           DO j = 1, N_dofs

              IF( Glo2Loc(j) > 0) THEN

                 ns = ns + 1
                 Glo2Loc(j) = ns

                 uu_ex = exact_solution(pb_type, coord(:, j), visc)

                 WRITE(UNIT, 200) coord(:, j), uu(j), uu_ex

              ENDIF

           ENDDO

           !-------------
           ! Connectivity
           !-----------------------------------------------
           DO jt = 1, N_elements

              IF(elements(jt)%p%N_verts == Nverts(type)) THEN
 
                 SELECT CASE(elements(jt)%p%Type)

                 CASE(TRI_P1)

                    WRITE(UNIT, *) ( Glo2Loc(elements(jt)%p%NU(k)), k = 1, elements(jt)%p%N_verts )

                 CASE(TRI_P2)

                     WRITE(UNIT, *)  Glo2Loc(elements(jt)%p%NU(1)), &
                                     Glo2Loc(elements(jt)%p%Nu(4)), &
                                     Glo2Loc(elements(jt)%p%Nu(6))

                     WRITE(UNIT, *)  Glo2Loc(elements(jt)%p%NU(2)), &
                                     Glo2Loc(elements(jt)%p%Nu(5)), &
                                     Glo2Loc(elements(jt)%p%Nu(4))

                     WRITE(UNIT, *)  Glo2Loc(elements(jt)%p%NU(3)), &
                                     Glo2Loc(elements(jt)%p%Nu(6)), &
                                     Glo2Loc(elements(jt)%p%Nu(5))

                     WRITE(UNIT, *)  Glo2Loc(elements(jt)%p%NU(6)), &
                                     Glo2Loc(elements(jt)%p%Nu(4)), &
                                     Glo2Loc(elements(jt)%p%Nu(5))

                  CASE(QUA_Q1)

                     WRITE(UNIT, *) ( Glo2Loc(elements(jt)%p%NU(k)), k = 1, elements(jt)%p%N_verts )

                  CASE(QUA_Q2)

                     WRITE(UNIT, *)  Glo2Loc(elements(jt)%p%NU(1)), &
                                     Glo2Loc(elements(jt)%p%Nu(5)), &
                                     Glo2Loc(elements(jt)%p%Nu(9)), &
                                     Glo2Loc(elements(jt)%p%Nu(8))


                     WRITE(UNIT, *)  Glo2Loc(elements(jt)%p%NU(5)), &
                                     Glo2Loc(elements(jt)%p%Nu(2)), &
                                     Glo2Loc(elements(jt)%p%Nu(6)), &
                                     Glo2Loc(elements(jt)%p%Nu(9))

                     WRITE(UNIT, *)  Glo2Loc(elements(jt)%p%NU(9)), &
                                     Glo2Loc(elements(jt)%p%Nu(6)), &
                                     Glo2Loc(elements(jt)%p%Nu(3)), &
                                     Glo2Loc(elements(jt)%p%Nu(7))

                     WRITE(UNIT, *)  Glo2Loc(elements(jt)%p%NU(8)), &
                                     Glo2Loc(elements(jt)%p%Nu(9)), &
                                     Glo2Loc(elements(jt)%p%Nu(7)), &
                                     Glo2Loc(elements(jt)%p%Nu(4))

                  CASE DEFAULT

                     WRITE(*,*) 'Wrong element type. STOP!'
                     STOP

                  END SELECT

               ENDIF

            ENDDO ! N_elements

         ENDIF
 
    ENDDO ! Type

   CLOSE(UNIT)


   
100 FORMAT('TITLE = "Scalar advection prblem 2d"')
110 FORMAT('VARIABLES = "X", "Y", "sol", "sol_ex"')
120 FORMAT('ZONE T="uu", F=FEPOINT, N=',I6, ', E=',I6 ', ET=TRIANGLE' )
121 FORMAT('ZONE T="uu", F=FEPOINT, N=',I6, ', E=',I6 ', ET=QUADRILATERAL' )
200 FORMAT(4(1x,e24.16))
300 FORMAT(3(1x,I6))   
      
    END SUBROUTINE plot_procedure
    !============================

    !===========================
    SUBROUTINE compute_error(uu)
    !===========================

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: uu      
      !---------------------------------------------

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: uu_ex
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: err_ele
      
      REAL(KIND=8) :: err_L2, err_Loo, int_uex
      REAL(KIND=8) :: err_Gx_L2, err_Gy_L2

      INTEGER :: je, i, ierror, UNIT
      !---------------------------------------------

      err_L2 = 0.d0; err_Loo = 0.d0; int_uex = 0.d0
      
      !---------------
      ! Exact solution
      !--------------------------------------------------
      DO je = 1, N_elements

         ALLOCATE(   uu_ex(elements(je)%p%N_points), &
                   err_ele(elements(je)%p%N_points)  )
     
         DO i = 1, elements(je)%p%N_points

            uu_ex(i) = exact_solution(pb_type, elements(je)%p%Coords(:, i), visc)

         ENDDO

         err_ele = uu(elements(je)%p%NU) - uu_ex
         
         !----------
         ! L_2 error
         !-------------------------------------
         err_L2  = err_L2  + Int_d(elements(je)%p, err_ele**2)
         int_uex = int_uex + Int_d(elements(je)%p, uu_ex**2)

         !-----------
         ! L_oo error
         !-------------------------------------
         err_Loo = MAX( err_Loo, &
                        MAXVAL(ABS(uu(elements(je)%p%NU) - uu_ex)) )

         DEALLOCATE( uu_ex, err_ele )

      ENDDO

      ! Normalize L2 error
      err_L2 = SQRT(err_L2) / SQRT(int_uex)

      ! Error of the gradients
      IF( is_visc ) THEN
         CALL compute_gradient_error(uu, err_Gx_L2, err_Gy_L2)
      ENDIF        

      WRITE(*,*) ' *** COMPUTE THE ERROR ***'
            
      WRITE(*,*) 'N. dof, err L2, err_Loo'
      WRITE(*,22) N_dofs, err_L2, err_Loo

      IF( is_visc ) THEN
         WRITE(*,*) 'L2 error of the gradients: E_G_x, E_G_y'
         WRITE(*,23) err_Gx_L2, err_Gy_L2
      ENDIF

      UNIT = 9
      
      OPEN(UNIT, FILE = 'error.'//TRIM(ADJUSTL(pb_name)), &
           ACTION= 'WRITE', IOSTAT = ierror)
           
      IF(ierror /= 0) THEN
         WRITE(*,*) 'ERROR: Impossible to open the file error'
         WRITE(*,*) 'STOP'
      
         STOP
      ENDIF

      WRITE(UNIT, *) 'N. dof, err L2, err_Loo'
      WRITE(UNIT,22)  N_dofs, err_L2, err_Loo

      IF( is_visc ) THEN
         WRITE(UNIT, *) 'L2 error of the gradients: E_G_x, E_G_y'
         WRITE(UNIT,23)  err_Gx_L2, err_Gy_L2
      ENDIF

      CLOSE(UNIT)
      
22 FORMAT(I6, 2(1x,e24.16))
23 FORMAT(2(1x,e24.16))
      
    END SUBROUTINE compute_error
    !===========================

    !==========================================================
    SUBROUTINE compute_gradient_error(uu, err_Gx_L2, err_Gy_L2)
    !==========================================================

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: uu
      REAL(KIND=8),               INTENT(OUT) :: err_Gx_L2
      REAL(KIND=8),               INTENT(OUT) :: err_Gy_L2
      !-------------------------------------------------

      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: D_uu

      REAL(KIND=8), DIMENSION(2) :: D_u_ex

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: err2_x, err2_y 
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: d2_x, d2_y

      INTEGER :: je, k, i

      REAL(KIND=8) :: int_d2_x, int_d2_y
      !---------------------------------------------

      err_Gx_L2 = 0.d0;  err_Gy_L2 = 0.d0
      int_d2_x = 0.d0;   int_d2_y = 0.d0

      ALLOCATE( D_uu(2, SIZE(uu)) )

      D_uu = Compute_gradient(uu)

      DO je = 1, N_elements

         ALLOCATE( err2_x(elements(je)%p%N_points), &
                   err2_y(elements(je)%p%N_points), &
                     d2_x(elements(je)%p%N_points), &
                     d2_y(elements(je)%p%N_points) )

         DO k = 1, elements(je)%p%N_points

            D_u_ex = exact_grad(pb_type, &
                                elements(je)%p%coords(:, k), &
                                visc )         

            err2_x(k) = ( D_uu(1, elements(je)%p%NU(k)) - D_u_ex(1) )**2
            err2_y(k) = ( D_uu(2, elements(je)%p%NU(k)) - D_u_ex(2) )**2

            d2_x(k) = D_u_ex(1)**2
            d2_y(k) = D_u_ex(2)**2

         ENDDO

         err_Gx_L2 = err_Gx_L2 + Int_d(elements(je)%p, err2_x)
         err_Gy_L2 = err_Gy_L2 + Int_d(elements(je)%p, err2_y)

         int_d2_x = int_d2_x + Int_d(elements(je)%p, d2_x)
         int_d2_y = int_d2_y + Int_d(elements(je)%p, d2_y)

         DEALLOCATE(  err2_x, err2_y, d2_x, d2_y )

      ENDDO
    
      err_Gx_L2 = DSQRT(err_Gx_L2)/DSQRT(int_d2_x)
      err_Gy_L2 = DSQRT(err_Gy_L2)/DSQRT(int_d2_y)

      DEALLOCATE( D_uu )

    END SUBROUTINE compute_gradient_error
    !====================================    
   

END MODULE Post_pro
