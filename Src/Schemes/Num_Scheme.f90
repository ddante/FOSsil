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

  PUBLIC :: distribute_residual, &
            distribute_residual_imp
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
       
  !=========================================================
  SUBROUTINE distribute_residual_imp(ele, Phi_tot, u, J_tot, &
                                     Phi_i, inv_dt)
  !=========================================================

    IMPLICIT NONE

    TYPE(element),                  INTENT(IN)  :: ele
    REAL(KIND=8), DIMENSION(:),     INTENT(IN)  :: Phi_tot
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: u
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: J_tot
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: Phi_i
    REAL(KIND=8),                   INTENT(OUT) :: inv_dt
    !-----------------------------------------------------

    SELECT CASE(scheme_type)

     CASE(LDA)

       CALL LDA_scheme_imp(ele, Phi_tot, u, J_tot, Phi_i, inv_dt)
      
    CASE(LLXFS)
stop       
       !CALL LLxFS_scheme_imp(ele, Phi_tot, u, Phi_i, inv_dt)

    CASE(LW)
stop       
       !CALL LW_scheme_imp(ele, Phi_tot, u, Phi_i, inv_dt)

    CASE DEFAULT

       WRITE(*,*) 'ERROR: unknown scheme type'
       STOP

    END SELECT    
    
  END SUBROUTINE distribute_residual_imp
  !=====================================

END MODULE Num_scheme
