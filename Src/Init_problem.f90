MODULE init_problem

  USE Geometry,       ONLY: N_dofs, elements,  N_elements
  USE Models,         ONLY: strong_bc

  IMPLICIT NONE

  PRIVATE

  !=======================================
  CHARACTER(len=64) :: pb_name
  INTEGER           :: order
  INTEGER           :: scheme_type
  INTEGER           :: bc_type
  INTEGER           :: mesh_format
  INTEGER           :: time_int
  INTEGER           :: pb_type
  INTEGER           :: ite_max
  REAL(KIND=8)      :: CFL
  REAL(KIND=8)      :: toll_res
  REAL(KIND=8)      :: visc

  INTEGER :: N_eqn
  !=======================================

  !=======================================
  
  PUBLIC :: read_param, initialization
  PUBLIC :: pb_name, order, time_int,    &
            scheme_type, mesh_format,    &
            pb_type, bc_type, visc, CFL, &
            ite_max, toll_res, N_eqn
  !=======================================

CONTAINS
  
   !=================================
   SUBROUTINE initialization(uu, rhs)
   !=================================
   !
   ! Initialize the solution at the first time step
   !
   IMPLICIT NONE
   
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: uu
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: rhs
      !-------------------------------------------------------------
      
      INTEGER :: ierror, UNIT
      !-------------------------------------------------------------

      N_eqn = 3
      
      ! Solution and RHS
      ALLOCATE( uu(N_eqn, N_dofs), rhs(N_eqn, N_dofs) )
      
      uu = 0.d0;  rhs = 0.d0

      CALL strong_bc(pb_type, visc, uu, rhs)

      ! Delete a possible previous convergence history...
      UNIT = 4
      
      OPEN(UNIT, FILE = 'convergence.'//TRIM(ADJUSTL(pb_name)), & 
          STATUS= 'REPLACE', IOSTAT = ierror)

      IF(ierror /= 0) THEN
         WRITE(*,*) 'ERROR: Impossible to open the file converegence'
         WRITE(*,*) 'STOP'
      
         STOP
      ENDIF     
      CLOSE (UNIT)
      
      ! ... and error file
      UNIT = 9
      
      OPEN(UNIT, FILE = 'error.'//TRIM(ADJUSTL(pb_name)), &
           STATUS= 'REPLACE', IOSTAT = ierror)

      IF(ierror /= 0) THEN
         WRITE(*,*) 'ERROR: Impossible to open the file error'
         WRITE(*,*) 'STOP'
      
         STOP
      ENDIF     
      CLOSE (UNIT)
   
   END SUBROUTINE initialization
   !============================

  !=======================================
  SUBROUTINE read_param (unit, param_file)
  !=======================================
  !
  ! Read the parameters of the simulation from the input file
  !
    IMPLICIT NONE

    INTEGER,           INTENT(IN) :: unit
    CHARACTER(len=64), INTENT(IN) :: param_file
    !===========================================

    INTEGER :: ierror
     
    OPEN(unit, FILE = param_file, ACTION = 'READ', IOSTAT = ierror)
    IF (ierror /= 0 ) THEN
       WRITE(*,*) 'ERROR: Impossible to open param file', TRIM(ADJUSTL(param_file))
       WRITE(*,*) 'STOP!'    
       STOP
    ENDIF
      
    READ(unit, *) pb_name
    READ(unit, *) order
    READ(unit, *) scheme_type
    READ(unit, *) mesh_format
    READ(unit, *) time_int
    READ(unit, *) bc_type
    READ(unit, *) pb_type
    READ(unit, *) ite_max
    READ(unit, *) toll_res
    READ(unit, *) CFL
    READ(unit, *) visc

    CLOSE (unit)
      
    WRITE(*,*)
    WRITE(*,*) 'Problem name: ',                   pb_name
    WRITE(*, '(" Order: " I2)')                    order
    WRITE(*, '(" Num scheme: " I2)')               scheme_type
    WRITE(*, '(" Mesh Format: " I2)')              mesh_format
    WRITE(*, '(" Time integration: " I2)')         time_int
    WRITE(*, '(" Boundary conditions type: " I2)') bc_type
    WRITE(*, '(" Problem type: " I2)')             pb_type
    WRITE(*, '(" Max Num. iterations: " I8)')      ite_max
    WRITE(*, '(" Residual tollerance: " E10.5)')   toll_res
    WRITE(*, '(" CFL Number: " F10.5)')            CFL
    WRITE(*, '(" Viscosity coefficient: " F10.5)') visc
   
  END SUBROUTINE read_param
  !========================

END MODULE init_problem
