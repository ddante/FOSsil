MODULE models

  USE geometry,      ONLY: N_dim, N_elements, elements

  IMPLICIT NONE

  PRIVATE

   !===============================================
   INTEGER, PARAMETER :: LINEAR_ADVECTION   = 0, &
                         ROTATION           = 1, &
                         BURGER             = 2, &
                         LIN_VISC_ADVEC     = 3, &
                         SMITH_HUTTON       = 4, &
                         PERAIRE_ADVEC_DIFF = 5, &
                         ADVECTION_SOURCE   = 6, &
                         MANUFACTED_SOURCE  = 7, &
                         TRANSPORT_REACTION = 8, &
                         BOUNDARY_LAYER     = 9, &
                         PURE_DIFFUSION     = 10
   !===============================================
   
   REAL(KIND=8), PARAMETER :: PI = DACOS(-1.d0)
   REAL(KIND=8), PARAMETER :: theta = 100.d0
   REAL(KIND=8), PARAMETER :: theta_xi = 0.d0
   REAL(KIND=8), PARAMETER :: mu_br = 0.01d0
   REAL(KIND=8), PARAMETER :: c_bl  = 0.059d0
   !===============================================

   PUBLIC :: advection_speed, advection_flux, &
             diffusion_flux, exact_solution, exact_grad, &
             strong_bc, source_term, detect_source
   !===============================================

 CONTAINS
   
   !=====================================================
   FUNCTION advection_speed(type_pb, uu, x, y)  RESULT(a)
   !=====================================================
   !
   ! Advection speed of the bidimensional 
   ! problem written in quasi-linear form
   !
   IMPLICIT NONE
   
      INTEGER,      INTENT(IN) :: type_pb
      REAL(KIND=8), INTENT(IN) :: uu
      REAL(KIND=8), INTENT(IN) :: x, y
      
      REAL(KIND=8), DIMENSION(2) :: a
      !------------------------------------
      
      SELECT CASE (type_pb)
      
         CASE(LINEAR_ADVECTION)
         
            a(1) = 0.d0
            a(2) = 1.d0
      
         CASE(ROTATION)
         
            a(1) = -y
            a(2) =  x
      
         CASE(BURGER)
         
            a(1) = uu 
            a(2) = 1.d0
             
	 CASE(LIN_VISC_ADVEC)
	 
	    a(1) = 1.d0 * DSIN(theta_xi*PI/180.d0)
	    a(2) = 1.d0 * DCOS(theta_xi*PI/180.d0)

         CASE(SMITH_HUTTON)
         
            a(1) =  2.d0*(1 - x**2)*y
            a(2) = -2.d0*(1 - y**2)*x

         CASE(PERAIRE_ADVEC_DIFF)

            a(1) = 1.d0
            a(2) = 2.d0*x

         CASE(ADVECTION_SOURCE)

            a(1) = y
            a(2) = 1.d0 - x
            
         CASE(MANUFACTED_SOURCE)

            a(1) = 0.d0
            a(2) = 1.d0

         CASE(TRANSPORT_REACTION)

            a(1) =  y
            a(2) = -x

            a = a / SQRT(x*x + y+y)

         CASE(BOUNDARY_LAYER)

            a(1) = 1.d0
            a(2) = 0.d0

         CASE(PURE_DIFFUSION)

            a(1) = 0.d0
            a(2) = 0.d0

         CASE DEFAULT
         
            WRITE(*,*) 'Problem of unknow type.'
            WRITE(*,*) 'STOP!'
            
            STOP
      
      END SELECT
   
   END fUNCTION advection_speed
   !===========================

   !=======================================================
   FUNCTION advection_flux(type_pb, uu, x, y)  RESULT(flux)
   !=======================================================
   !
   ! Flux of the bidimensional problem 
   ! written in quasi-linear form
   !
   IMPLICIT NONE
   
      INTEGER,      INTENT(IN) :: type_pb
      REAL(KIND=8), INTENT(IN) :: uu
      REAL(KIND=8), INTENT(IN) :: x, y
      
      REAL(KIND=8), DIMENSION(2) :: flux
      !------------------------------------
      
      SELECT CASE (type_pb)
      
         CASE(LINEAR_ADVECTION)
         
            flux(1) = 0.d0
            flux(2) = uu
      
         CASE(ROTATION)
         
            flux(1) = -y*uu
            flux(2) =  x*uu
      
         CASE(BURGER)
         
            flux(1) = 0.5d0 * uu**2 
            flux(2) = uu

	 CASE(LIN_VISC_ADVEC)
	 
	    flux(1) = uu * DSIN(theta_xi*PI/180.d0)
            flux(2) = uu * DCOS(theta_xi*PI/180.d0)
                  
         CASE(SMITH_HUTTON)
         
            flux(1) =  2.d0*(1 - x**2)*y * uu
            flux(2) = -2.d0*(1 - y**2)*x * uu

         CASE(PERAIRE_ADVEC_DIFF)
            
            flux(1) = 1.d0*uu
            flux(2) = 2.d0*x*uu

         CASE(ADVECTION_SOURCE)

            flux(1) = y*uu
            flux(2) = (1.d0 - x)*uu

         CASE(MANUFACTED_SOURCE)            

            flux(1) = 0.d0
            flux(2) = uu

         CASE(TRANSPORT_REACTION)            

            flux(1) =  y * uu
            flux(2) = -x * uu

            flux = flux / SQRT(x*x + y*y)

         CASE(BOUNDARY_LAYER)

            flux(1) = uu
            flux(2) = 0.d0

         CASE(PURE_DIFFUSION)

            flux(1) = 0.d0
            flux(2) = 0.d0

         CASE DEFAULT
         
            WRITE(*,*) 'Problem of unknow type.'
            WRITE(*,*) 'STOP!'
            
            STOP
      
      END SELECT
   
   END FUNCTION advection_flux
   !==========================

   !=============================================================
   FUNCTION diffusion_flux(type_pb, mu, D_uu, x, y)  RESULT(flux)
   !=============================================================
   !
   ! Flux of the bidimensional problem 
   !
   IMPLICIT NONE
   
      INTEGER,                    INTENT(IN) :: type_pb
      REAL(KIND=8),               INTENT(IN) :: mu
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: D_uu
      REAL(KIND=8),               INTENT(IN) :: x, y
      
      REAL(KIND=8), DIMENSION(2) :: flux
      !------------------------------------------------
      
      SELECT CASE (type_pb)
         
         CASE(LIN_VISC_ADVEC)

            flux = mu * D_uu

         CASE(BOUNDARY_LAYER)

            flux = mu * D_uu

         CASE(PURE_DIFFUSION)

            flux = mu * D_uu

         CASE DEFAULT

!!$            WRITE(*,*) 'Problem of unknow type.'
!!$            WRITE(*,*) 'STOP!'
!!$            
!!$            STOP

            RETURN

      END SELECT         

   END FUNCTION diffusion_flux
   !==========================
       
   !====================================================
   FUNCTION source_term(type_pb, uu, mu, x, y) RESULT(S)
   !====================================================
   IMPLICIT NONE

      INTEGER,      INTENT(IN) :: type_pb
      REAL(KIND=8), INTENT(IN) :: uu
      REAL(KIND=8), INTENT(IN) :: mu
      REAL(KIND=8), INTENT(IN) :: x, y
      
      REAL(KIND=8) :: S
      !------------------------------------

      REAL(KIND=8) :: r, c, t1, t2, t3, t9, t13, &
                         t14, t15, t16, t24, t25, t27
      !------------------------------------

      S = 0.d0
      
      SELECT CASE (type_pb)

      CASE(ADVECTION_SOURCE)

         r = SQRT((x - 0.5d0)**2 + (y - 0.3)**2)

         IF ( r <= 0.25 ) THEN
            S = 10.d0
         ELSE
            S = 0.d0
         ENDIF
         
      CASE(MANUFACTED_SOURCE)
         
         s = DSIN(PI*x)*DCOS(PI*x)*DCOS(PI*y)*DCOS(PI*y)*PI - &
             DSIN(PI*x)*DCOS(PI*x)*DSIN(PI*y)*DSIN(PI*y)*PI

      CASE(TRANSPORT_REACTION)
         
         s = mu_br*uu

      CASE(BOUNDARY_LAYER)

         t1 = c_bl * mu
         t2 = t1 * x
         t3 = SQRT(t2)
         t9 = EXP(-y / t3)
         t13 = c_bl ** 2
         t14 = mu ** 2
         t15 = t13 * t14
         t16 = x ** 2
         t24 = y ** 2
         t25 = 0.1D1 / c_bl
         t27 = 0.1D1 / mu
         
         s = -y / t3 / t2 * t1 * t9 / 0.2D1 - mu * &
              (0.3D1 / 0.4D1 * y / t3 / t16 * t9 - &
              t24 * t25 * t27 / x / t16 * t9 / 0.4D1 - t25 * t27 / x * t9)


      CASE DEFAULT

!!$         WRITE(*,*) 'Problem of unknow type.'         
!!$         WRITE(*,*) 'STOP!'
!!$         
!!$         STOP

         RETURN

      END SELECT
      
   END FUNCTION source_term
   !=======================

   !============================================
   FUNCTION detect_source(type_pb) RESULT(logic)
   !============================================

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: type_pb

     LOGICAL :: logic
     !------------------------------

     SELECT CASE(type_pb)

        CASE( ADVECTION_SOURCE,   &
              MANUFACTED_SOURCE,  &
              TRANSPORT_REACTION, &
              BOUNDARY_LAYER )

           logic = .TRUE.

        CASE DEFAULT

           logic = .FALSE.

        END SELECT        
     
   END FUNCTION detect_source
   !=========================
 
   !===========================================
   SUBROUTINE strong_bc(type_pb, visc, uu, rhs)
   !===========================================
   !
   ! Impose the bondary conditions in the strong way,
   ! according to the chosen type of the problem
   !
   IMPLICIT NONE
   
      INTEGER,                    INTENT(IN)    :: type_pb
      REAL(KIND=8),               INTENT(IN)    :: visc
      REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: uu
      REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: rhs
      !======================================================
      
      INTEGER, DIMENSION(:), ALLOCATABLE :: is, is_loc
      
      REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: u_l, rhs_l
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: coord
      
      INTEGER :: n_loc, Nv
      
      LOGICAL :: b_flag
      
      INTEGER :: jt, i, k
      !=======================================================
      
      DO jt = 1, N_elements

         Nv = elements(jt)%p%N_points

         ALLOCATE( is(Nv), is_loc(Nv), u_l(Nv), rhs_l(Nv) )
         ALLOCATE( coord(N_dim, Nv) )
      
         is = elements(jt)%p%NU
         
         DO i = 1, N_dim
            coord(i, :) = elements(jt)%p%Coords(i, :)
         ENDDO

           u_l(:) =  uu(is)
         rhs_l(:) = rhs(is)

          n_loc = 0
         is_loc = 0
         b_flag = .FALSE.
         
         SELECT CASE(type_pb)
      
         CASE(LINEAR_ADVECTION)
      
            CALL bc_linadv()
               
         CASE(ROTATION)
            
            CALL bc_rotation()
               
         CASE(BURGER)
            
            CALL bc_burger()
               
	 CASE(LIN_VISC_ADVEC)
	    
	    CALL bc_linviscadv()

         CASE(SMITH_HUTTON)
            
            CALL bc_smithutton()

         CASE(PERAIRE_ADVEC_DIFF)

            CALL bc_peraire_advec_diff()

         CASE(ADVECTION_SOURCE)

            CALL bc_advection_source()

         CASE(MANUFACTED_SOURCE)

            CALL bc_manufacted_source()

         CASE(TRANSPORT_REACTION)

            CALL bc_transport_reaction()

         CASE(BOUNDARY_LAYER)

            CALL bc_boundary_layer()

         CASE(PURE_DIFFUSION)

            CALL bc_pure_diffusion()

         CASE DEFAULT
         
            WRITE(*,*) 'Problem of unknow type.'
            WRITE(*,*) 'STOP!'   
            STOP
            
         END SELECT
         
         IF (b_flag) THEN
               
            DO i = 1, n_loc
                uu( is(is_loc(i)) ) =   u_l( is_loc(i) )            
               rhs( is(is_loc(i)) ) = rhs_l( is_loc(i) )
            ENDDO
                     
         ENDIF

         DEALLOCATE( is, is_loc, coord, u_l, rhs_l )
      
      ENDDO

!=======
CONTAINS
!=======

      !......................
      SUBROUTINE bc_linadv()
      !
      ! Square [0,1]x[0,1]
      !
      IMPLICIT NONE

         DO k = 1, Nv
            IF ( ABS(coord(1, k)) <= 0.d0 .AND. & 
                 ABS(coord(2, k)) >  0.d0) THEN
               u_l(k) = initial_profile(coord(1, k))
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.
            ENDIF
         ENDDO
          
         DO k = 1, Nv
            IF ( (ABS(coord(1, k) - 1.d0)) <= 0.d0) THEN
               u_l(k) = initial_profile(coord(1, k))
               rhs_l(k) = 0.d0              
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.             
            ENDIF            
         ENDDO
          
         DO k = 1, Nv
            IF ( ABS(coord(2, k)) <= 0.d0 ) THEN
               u_l(k) = initial_profile(coord(1, k))
               rhs_l(k) = 0.d0              
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.             
            ENDIF            
         ENDDO        
      
      END SUBROUTINE bc_linadv
      !........................

      !........................
      SUBROUTINE bc_rotation()
      !
      ! Square [0,1]x[0,1]
      !
      IMPLICIT NONE

         DO k = 1, Nv
            IF ( (ABS(coord(1, k) - 1.d0)) <= 0.d0 ) THEN
               u_l(k) = 0.d0
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.
            ENDIF
         ENDDO
                    
         DO k = 1, Nv
            IF ( ABS(coord(2, k)) <= 0.d0 ) THEN
               IF (coord(1, k) >= 0.25d0 .AND. coord(1, k) <= 0.75d0) THEN
                  u_l(k) = DCOS(2.d0*PI*coord(1, k))**2
               ELSE
                  u_l(k) = 0
               ENDIF
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.           
            ENDIF            
         ENDDO        
      
      END SUBROUTINE bc_rotation
      !........................
      
      !......................
      SUBROUTINE bc_burger()
      !
      ! Square [0,1]x[0,1]
      !
      IMPLICIT NONE

         DO k = 1, Nv
            IF ( ABS(coord(2, k)) <= 0.0 ) THEN
               u_l(k) = 1.5d0 - 2.d0*coord(1, k)
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1              
               is_loc(n_loc) = k
               b_flag = .TRUE.  
            ENDIF
         ENDDO
                    
         DO k = 1, Nv
            IF ((ABS(coord(1, k)  - 1.d0)) <= 0.0 ) THEN
               u_l(k) = -0.5d0
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1              
               is_loc(n_loc) = k
               b_flag = .TRUE.             
            ENDIF            
         ENDDO 
                
         DO k = 1, Nv
            IF ( ABS(coord(1, k)) <= 0.d0 .AND. &
                 ABS(coord(2, k)) >= 0.d0) THEN
               u_l(k) = 1.5d0 - 2.d0*coord(1, k)
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1              
               is_loc(n_loc) = k
               b_flag = .TRUE.         
            ENDIF            
         ENDDO   
             
      END SUBROUTINE bc_burger
      !........................

      !.........................
      SUBROUTINE bc_linviscadv()
      !
      ! Square [0,1]x[0,1]
      !
      IMPLICIT NONE

         DO k = 1, Nv
            IF ( ABS(coord(1, k)) <= 0.d0 .AND. & 
                 ABS(coord(2, k)) >  0.d0 ) THEN
               u_l(k) = visc_advection(coord(1, k), coord(2, k), visc)
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.             
            ENDIF
         ENDDO

         DO k = 1, Nv           
            IF ( (ABS(coord(1, k)  - 1.d0)) <= 0.d0 ) THEN
               u_l(k) = visc_advection(coord(1, k), coord(2, k), visc)
               rhs_l(k) = 0.d0              
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.             
            ENDIF            
         ENDDO

         DO k = 1, Nv
            IF ( ABS(coord(2, k)) <= 0.d0 )THEN
               u_l(k) =  visc_advection(coord(1, k), coord(2, k), visc)
               rhs_l(k) = 0.d0              
               n_loc = n_loc + 1              
               is_loc(n_loc) = k
               b_flag = .TRUE.             
            ENDIF            
         ENDDO   
 
      END SUBROUTINE bc_linviscadv
      !...........................
      
      !.........................
      SUBROUTINE bc_smithutton()
      !
      ! Rectangle [-1,1]x[0,1]
      !
      IMPLICIT NONE

         DO k = 1, Nv
            IF ( ABS(coord(2, k)) <= 0.d0 .AND. &
                 coord(1, k) <= 0.d0 ) THEN
	       u_l(k) = 1 + TANH(theta*(2.d0*coord(1, k) + 1.d0))
	       rhs_l(k) = 0.d0
	       n_loc = n_loc + 1              
	       is_loc(n_loc) = k
	       b_flag = .TRUE.               
            ENDIF
         ENDDO
                    
         DO k = 1, Nv
            IF ((ABS(coord(2, k) - 1.d0)) <= 0.d0 ) THEN
               u_l(k) = 1.d0 - TANH(theta)
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1              
               is_loc(n_loc) = k
               b_flag = .TRUE.             
            ENDIF            
         ENDDO 
                
         DO k = 1, Nv
            IF ( (ABS(coord(1, k) - 1.d0)) <= 0.d0 .AND. &
                  ABS(coord(2, k))         >= 0.d0 ) THEN
               u_l(k) = 1.d0 - TANH(theta)
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1              
               is_loc(n_loc) = k
               b_flag = .TRUE.         
            ENDIF            
         ENDDO   
         
         DO k = 1, Nv
            IF ( (ABS(coord(1, k) + 1.d0)) <= 0.d0 .AND. &
                  ABS(coord(2, k) - 1.d0) >= 0.d0 ) THEN
               u_l(k) = 1.d0 - TANH(theta)
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1              
               is_loc(n_loc) = k
               b_flag = .TRUE.         
            ENDIF            
         ENDDO
         
      END SUBROUTINE bc_smithutton
      !............................

      !..................................
      SUBROUTINE  bc_peraire_advec_diff()
      !
      ! Square [0,1]x[0,1]
      !
      IMPLICIT NONE
     
        DO k = 1, Nv
            IF ( ABS(coord(1, k)) <= 0.d0 .AND. & 
                 ABS(coord(2, k)) >= 0.d0) THEN
               u_l(k) = 1.d0 - coord(2, k)
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.
            ENDIF
         ENDDO

        DO k = 1, Nv
            IF ( ABS(coord(1, k)) >=  0.d0 .AND. & 
                 ABS(coord(2, k)) <= 0.d0) THEN
               u_l(k) = coord(1, k) - 1.d0
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.
            ENDIF
         ENDDO         
           
      END SUBROUTINE bc_peraire_advec_diff
      !...................................

      !...............................
      SUBROUTINE bc_advection_source()
      !
      ! Square [0,1]x[0,1]
      !
      IMPLICIT NONE

        DO k = 1, Nv
            IF ( ABS(coord(1, k)) <= 0.d0 .AND. & 
                 ABS(coord(2, k)) >= 0.d0) THEN
               u_l(k) = 0.d0
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.
            ENDIF
         ENDDO       

        DO k = 1, Nv
            IF ( ABS(coord(1, k)) >= 0.d0 .AND. & 
                 ABS(coord(2, k)) <= 0.d0) THEN

               IF ( coord(1, k) >= 0.3d0 .AND. &
                    coord(1, k) <= 0.8d0 ) THEN
                  u_l(k) = -5.d0
               ELSE                  
                  u_l(k) = -0.d0                  
               ENDIF               

               rhs_l(k) = 0.d0
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.
            ENDIF
         ENDDO
      
      END SUBROUTINE bc_advection_source
      !.................................

      !................................
      SUBROUTINE bc_manufacted_source()
      !
      ! Square [0,1]x[0,1]
      !
      IMPLICIT NONE

         DO k = 1, Nv
            IF ( ABS(coord(1, k)) <= 0.d0 .AND. & 
                 ABS(coord(2, k)) >= 0.d0) THEN
               u_l(k) = 1.d0
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.
            ENDIF
         ENDDO
          
         DO k = 1, Nv
            IF ( (ABS(coord(1, k) - 1.d0)) <= 0.d0) THEN
               u_l(k) = 1.d0
               rhs_l(k) = 0.d0              
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.             
            ENDIF            
         ENDDO
          
         DO k = 1, Nv
            IF ( ABS(coord(2, k)) <= 0.d0 ) THEN
               u_l(k) = 1.d0
               rhs_l(k) = 0.d0              
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.             
            ENDIF            
         ENDDO        
      
      END SUBROUTINE bc_manufacted_source
      !..................................

      !.................................
      SUBROUTINE bc_transport_reaction()
      !
      ! arc circle r=[0.1, 1.0]
      !
      IMPLICIT NONE

         DO k = 1, Nv
            IF ( coord(1, k) <= 0.d0 ) THEN
               u_l(k) = u_ex_transport_reaction(coord(1, k), &
                                                coord(2, k) )
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.
            ENDIF
         ENDDO
      
      END SUBROUTINE bc_transport_reaction
      !...................................

      !.............................
      SUBROUTINE bc_boundary_layer()
      !
      ![0.05, 1.05]x[0, 0.001]
      ! 
      IMPLICIT NONE

         DO k = 1, Nv
            IF ( ABS(coord(2, k)) <= 0.d0 ) THEN
               u_l(k) = u_ex_bl(coord(1, k), coord(2, k), visc)
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.
            ENDIF
         ENDDO

         DO k = 1, Nv
            IF ( (ABS(coord(1, k) - 0.05d0)) < 0.d0 ) THEN
               u_l(k) = u_ex_bl(coord(1, k), coord(2, k), visc)
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.
            ENDIF
         ENDDO

      END SUBROUTINE bc_boundary_layer
      !...............................

      !................................
      SUBROUTINE  bc_pure_diffusion()
      !
      ! Square [0,1]x[0,1]
      !
      IMPLICIT NONE

         DO k = 1, Nv
            IF ( ABS(coord(1, k)) <= 0.d0 .AND. & 
                 ABS(coord(2, k)) >  0.d0) THEN
               u_l(k) = u_ex_pure_diff(coord(1, k), coord(2, k))
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.
            ENDIF
         ENDDO
          
         DO k = 1, Nv
            IF ( (ABS(coord(1, k) - 1.d0)) <= 0.d0 .AND. &
                  ABS(coord(2, k)) >  0.d0) THEN
               u_l(k) = u_ex_pure_diff(coord(1, k), coord(2, k))
               rhs_l(k) = 0.d0              
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.             
            ENDIF            
         ENDDO
          
         DO k = 1, Nv
            IF ( ABS(coord(2, k)) <= 0.d0 ) THEN
               u_l(k) = u_ex_pure_diff(coord(1, k), coord(2, k))
               rhs_l(k) = 0.d0              
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.             
            ENDIF            
         ENDDO
 
         DO k = 1, Nv
            IF ( ABS(coord(2, k) - 1.d0) <= 0.d0 ) THEN
               u_l(k) = u_ex_pure_diff(coord(1, k), coord(2, k))
               rhs_l(k) = 0.d0              
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.             
            ENDIF            
         ENDDO 

       END SUBROUTINE bc_pure_diffusion       
       !..................................

   END SUBROUTINE strong_bc
   !=======================
   
   !========================================================
   FUNCTION exact_solution(type_pb, coord, nu) RESULT(uu_ex)
   !========================================================
   !
   ! Compute the exact solution of the 2d problem, when possible
   !
   IMPLICIT NONE
   
      INTEGER,                    INTENT(IN) :: type_pb
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: coord
      REAL(KIND=8),               INTENT(IN) :: nu
      
      REAL(KIND=8) :: uu_ex
      !-------------------------------------------------
      
      REAL(KIND=8) :: x, y, r, phi    
      !-------------------------------------------------

      x = coord(1)
      y = coord(2)
      
      SELECT CASE (type_pb)
      
      CASE(LINEAR_ADVECTION)
                  
         uu_ex = initial_profile(x)
                                 
      CASE(ROTATION)
                  
         r = SQRT(x**2 + y**2)         
         IF( r >= 0.25 .AND. r <= 0.75 ) THEN
            uu_ex = DCOS(2*PI*r)**2
         ELSE            
            uu_ex = 0.d0            
         ENDIF
         
            
      CASE(BURGER)
                  
         IF (y >= 0.5) THEN
            IF (-2.d0*(x - 3.d0/4.d0) + (y - 0.5) <= 0) THEN
               uu_ex = -0.5
            ELSE
               uu_ex = 1.5   
            ENDIF               
         ELSE                        
            uu_ex = MAX( -0.5, MIN( 1.5, ((x - 3.d0/4.d0)/(y - 0.5)) ) )
         ENDIF


      CASE(LIN_VISC_ADVEC)
      
         IF (nu > 0.d0) THEN    
            uu_ex = visc_advection(x, y, nu)	          
         ELSE	    
            WRITE(*,*) 'ERROR: LIN_VISC_ADVEC'
            WRITE(*,*) 'STOP!'
            STOP
         ENDIF                 
            
      CASE(SMITH_HUTTON)
                       
         phi = -(1.d0 - x**2)*(1 - y**2)         
         uu_ex = 1.d0 + TANH(theta*(1.d0 - 2.d0*SQRT(1.d0 + phi)))

      CASE(MANUFACTED_SOURCE)

         uu_ex = 1.d0 + DSIN(PI*x)*DCOS(PI*x)*DSIN(PI*y)*DCOS(PI*y)
         
      CASE(TRANSPORT_REACTION)

         uu_ex = u_ex_transport_reaction(x, y)

      CASE(BOUNDARY_LAYER)

         IF (nu > 0.d0) THEN    
            uu_ex = u_ex_bl(x, y, nu)	          
         ELSE	    
            WRITE(*,*) 'ERROR: LIN_VISC_ADVEC'            
            WRITE(*,*) 'STOP!'
            STOP       
         ENDIF

      CASE(PURE_DIFFUSION)

         uu_ex = u_ex_pure_diff(x, y)

      CASE DEFAULT

         stop
         !uu_ex = 0.0
      
      END SELECT
   
   END FUNCTION exact_solution
   !==========================   

   !======================================================
   FUNCTION exact_grad(type_pb, coord, nu) RESULT(G_uu_ex)
   !======================================================
   !
   ! Compute the exact gradient solution
   !
   IMPLICIT NONE
   
      INTEGER,                    INTENT(IN) :: type_pb
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: coord
      REAL(KIND=8),               INTENT(IN) :: nu
      
      REAL(KIND=8), DIMENSION(2) :: G_uu_ex
      !-------------------------------------------------
      
      REAL(KIND=8) :: x, y, r, phi    
      !-------------------------------------------------

      x = coord(1)
      y = coord(2)
      
      SELECT CASE (type_pb)
      
      CASE(LIN_VISC_ADVEC)

         G_uu_ex = grad_visc_advection(x, y, nu)
         
      CASE(PURE_DIFFUSION)

         G_uu_ex = grad_pure_diff(x, y)

      CASE DEFAULT

         !stop
         !G_uu_ex = 0.0
         
      END SELECT

    END FUNCTION exact_grad
    !======================      

   !====================================
   FUNCTION initial_profile(x) RESULT(f)
   !====================================
   !
   ! Initial profile of the solution used as 
   ! inlet boundary condition in some problems
   !
   IMPLICIT NONE
   
      REAL(KIND=8), INTENT(IN) :: x
      
      REAL(KIND=8) :: f
      !--------------------------------------
      
      REAL(KIND=8), PARAMETER :: gamma = 100
      REAL(KIND=8) :: a, b      
      !--------------------------------------

      a = 1.d0
      b = a*gamma*(x - 0.5)**2
      
      f = EXP(-b)*SQRT(a)   
   
   END FUNCTION initial_profile
   !===========================
   
   !================================================
   FUNCTION grad_visc_advection(x, y, nu) RESULT(GG)
   !================================================
   !
   ! Gradient of the viscous linear advection problem
   !
   IMPLICIT NONE
   
      REAL(KIND=8), INTENT(IN) :: x
      REAL(KIND=8), INTENT(IN) :: y
      REAL(KIND=8), INTENT(IN) :: nu
      
      REAL(KIND=8), DIMENSION(2) :: GG
      !-----------------------------

      REAL(KIND=8) :: t5, t6, t11, t12, t16, t17, t19, t22, t26
      REAL(KIND=8) :: ax, ay
      !--------------------------------

      ax = DSIN(theta_xi*PI/180.d0)
      ay = DCOS(theta_xi*PI/180.d0)
      
      IF (nu > 0.d0) THEN

         t5 = 0.2D1 * PI * (ay*x - ax*y)
         t6 = DSIN(t5)
         t11 = PI ** 2
         t12 = nu ** 2
         t16 = DSQRT(0.1D1 + 0.16D2 * t11 * t12)
         t17 = 0.1D1 - t16
         t19 = 0.1D1 / nu
         t22 = EXP((ax*x + ay*y) * t17 * t19 / 0.2D1)
         t26 = DCOS(t5)
         GG(1) = 0.2D1 * t6 * PI * ay * t22 - t26 * ax * t17 * t19 * t22 / 0.2D1

         t5 = 0.2D1 * PI * (ay * x - ax * y)
         t6 = DSIN(t5)
         t11 = PI ** 2
         t12 = nu ** 2
         t16 = DSQRT(0.1D1 + 0.16D2 * t11 * t12)
         t17 = 0.1D1 - t16
         t19 = 0.1D1 / nu
         t22 = EXP((ax*x + ay*y) * t17 * t19 / 0.2D1)
         t26 = DCOS(t5)
         GG(2) = -0.2D1 * t6 * PI * ax * t22 - t26 * ay * t17 * t19 * t22 / 0.2D1

      ELSE

         WRITE(*,*) 'ERROR: Visc zero'
         WRITE(*,*) 'STOP!'
         STOP

      ENDIF
   
    END FUNCTION grad_visc_advection
    !===============================

   !==========================================
   FUNCTION visc_advection(x, y, nu) RESULT(f)
   !==========================================
   !
   ! Solution of the viscous linear advection problem
   !
   IMPLICIT NONE
   
      REAL(KIND=8), INTENT(IN) :: x
      REAL(KIND=8), INTENT(IN) :: y
      REAL(KIND=8), INTENT(IN) :: nu
      
      REAL(KIND=8) :: f
      REAL(KIND=8) :: eta, psi
      !-----------------------------
      
      eta = DCOS(theta_xi*PI/180.d0)*x - DSIN(theta_xi*PI/180.d0)*y
      psi = DSIN(theta_xi*PI/180.d0)*x + DCOS(theta_xi*PI/180.d0)*y

      IF (nu > 0.d0) THEN
	 f = -DCOS(2.d0*PI*eta)* & 
	     EXP(0.5d0*psi*( 1.d0 - SQRT(1.d0 + 16.d0*(PI**2)*nu**2) )/nu)
      ELSE
         WRITE(*,*) 'ERROR: Visc zero'
         WRITE(*,*) 'STOP!'
         STOP
      ENDIF
   
    END FUNCTION visc_advection
    !==========================

    !===============================================
    FUNCTION u_ex_transport_reaction(x, y) RESULT(f)
    !===============================================
    IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: x
      REAL(KIND=8), INTENT(IN) :: y
      
      REAL(KIND=8) :: f
      !--------------------------------------------

      REAL(KIND=8) :: r

      r = SQRT( x*x + y*y )
      
      f = EXP( mu_br*r * ASIN(y/r) ) * ATAN( (r - 0.5d0)/0.1 )
    
    END FUNCTION u_ex_transport_reaction
    !===================================

    !===================================
    FUNCTION u_ex_bl(x, y, nu) RESULT(f)
    !===================================

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: x
      REAL(KIND=8), INTENT(IN) :: y
      REAL(KIND=8), INTENT(IN) :: nu

      REAL(KIND=8) :: f
      !------------------------------

      f = 1.d0 - EXP( -y/SQRT(c_bl * nu * x) )
      
    END FUNCTION u_ex_bl
    !===================

    !=====================================
    FUNCTION u_ex_pure_diff(x,y) RESULT(f)
    !=====================================
      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: x
      REAL(KIND=8), INTENT(IN) :: y
      
      REAL(KIND=8) :: f
      !--------------------------------------------
      
      f = (SINH(Pi*x)*SIN(Pi*y) + SINH(Pi*y)*SIN(Pi*x)) / SINH(Pi)

     END FUNCTION u_ex_pure_diff
     !==========================

     !=======================================
     FUNCTION grad_pure_diff(x, y) RESULT(GG)
     !=======================================
   
       IMPLICIT NONE
   
       REAL(KIND=8), INTENT(IN) :: x
       REAL(KIND=8), INTENT(IN) :: y
            
       REAL(KIND=8), DIMENSION(2) :: GG
       !--------------------------------
      
       GG(1) = (COSH(Pi*x)*Pi*SIN(Pi*y) + SINH(Pi*y)*COS(Pi*x)*Pi) / SINH(Pi)
       GG(2) = (SINH(Pi*x)*COS(Pi*y)*Pi + COSH(Pi*y)*Pi*SIN(Pi*x)) / SINH(Pi)

     END FUNCTION grad_pure_diff
     !==========================

END MODULE models
