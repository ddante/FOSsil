MODULE Num_scheme

  USE element_class
  USE init_problem,   ONLY: scheme_type
  USE LLxFS_method
  USE LW_method
  USE LDA_method

  IMPLICIT NONE
  PRIVATE

  !=====================================
  INTEGER, PARAMETER :: LN    = 1, &
                        LDA   = 2, &
                        LLXFS = 3, &
                        LW    = 4
  !=====================================

  PUBLIC :: distribute_residual
  !=====================================
  
CONTAINS

  !=============================================================
  SUBROUTINE distribute_residual(ele, Phi_tot, u, Phi_i, inv_dt)
  !=============================================================

    IMPLICIT NONE

    TYPE(element),                INTENT(IN)  :: ele
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: Phi_tot
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: u
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: Phi_i
    REAL(KIND=8),                 INTENT(OUT) :: inv_dt
    !--------------------------------------------------

    SELECT CASE(scheme_type)

     CASE(LDA)

       CALL LDA_scheme(ele, Phi_tot, u, Phi_i, inv_dt)
      
    CASE(LLXFS)
       
       CALL LLxFS_scheme(ele, Phi_tot, u, Phi_i, inv_dt)

    CASE(LW)
       
       CALL LW_scheme(ele, Phi_tot, u, Phi_i, inv_dt)

    CASE DEFAULT

       WRITE(*,*) 'ERROR: unknown scheme type'
       STOP

    END SELECT    
    
  END SUBROUTINE distribute_residual
  !=================================
       
END MODULE Num_scheme
