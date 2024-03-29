MODULE time_integration

  USE Element_Class
  
  USE geometry,          ONLY: N_dofs, N_elements, elements
  USE init_problem,      ONLY: pb_name, pb_type, visc, &
                               time_int, CFL, CFL_l, visc, N_eqn

  USE space_integration

  USE petsc_driver,      ONLY: solve_sys

  IMPLICIT NONE
  PRIVATE

  REAL(KIND=8), DIMENSION(3) :: res_0
  
  PUBLIC :: time_advance

CONTAINS
  
  !=========================================
  SUBROUTINE time_advance(ite, uu, rhs, res)
  !=========================================

    IMPLICIT NONE

    INTEGER,                        INTENT(INOUT) :: ite
    REAL(KIND=8), DIMENSION(:,:),   INTENT(INOUT) :: uu
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT)   :: rhs
    REAL(KIND=8), DIMENSION(N_eqn), INTENT(OUT)   :: res
    !----------------------------------------------------

    REAL(KIND=8), DIMENSION(N_dofs) :: Dt_V
    
    INTEGER :: i, UNIT, ierror
    !----------------------------------------------------
     
    IF ( time_int == 0 ) THEN

       CALL compute_rhs(uu, rhs, Dt_V)

       !---------------
       ! Esplicit Euler
       !----------------------------------
       DO i = 1, N_eqn
          uu(i,:) = uu(i,:) - Dt_V*rhs(i,:)
       ENDDO

    ELSEIF ( time_int == 1 ) THEN

       CALL compute_lhs_rhs(uu, rhs)

       !---------------
       ! Implicit Euler
       !---------------------------------
       CALL solve_sys(rhs)

       DO i = 1, N_eqn
          uu(i,:) = uu(i,:) + rhs(i,:)
       ENDDO
       
    ENDIF
   
    !-----------------------
    ! Normalized L2 Residual
    !-----------------------------------------------
    DO i = 1, N_eqn
       res(i) = SQRT(SUM(rhs(i,:)**2)) / REAL(N_dofs,8)
    ENDDO

    IF(ite == 1) res_0 = res
    res = res/res_0

    !------------------------
    ! Save convergence strory
    !------------------------------------------------
    UNIT = 3

    OPEN(UNIT, FILE = 'convergence.'//TRIM(ADJUSTL(pb_name)), &
         ACTION = 'WRITE', POSITION = 'APPEND', IOSTAT = ierror)
            
    IF(ierror /= 0) THEN
       
       WRITE(*,*) 'ERROR: Impossible to open the file converegence'
       WRITE(*,*) 'STOP'
       STOP
       
    ENDIF
    
    WRITE(UNIT, 500) ite, res
    
    CLOSE(UNIT)
     
    !--------------------
    ! Data on the screan
    !------------------------------------------------
    WRITE(*, 600) ite, MINVAL(uu(1,:)), MAXVAL(uu(1,:))
    WRITE(*, 601) res
    WRITE(*, 602) res_0
    IF ( time_int == 1 ) THEN
       WRITE(*, 603) CFL_l
    ENDIF
    WRITE(*, *)

500 FORMAT(I6, 3F24.16)      
600 FORMAT('Ite # ', I6, ';   uu min = ', F10.6, ', uu max =', F10.6)
601 FORMAT('Res   = ', 3(E12.5, ' | ') )
602 FORMAT('Res_0 = ', 3(E12.5, ' | ') )
603 FORMAT('CFL   = ', E12.5 )

  END SUBROUTINE time_advance
  !==========================

END MODULE time_integration
