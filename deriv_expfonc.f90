!---------------------------------------------------------------------------------------------------
! All the formulations  have been taken from the original
! article,the Span-Wagner EoS for CO2.
! The EoS has been published in J. Phys. Chem. Ref. Data, Vol.25,
! pp. 1509-1596, 1996.
!---------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
! SUBROUTINE   deriv_expfonc
! @brief   Compute expotential fuction Psi p.1545 Table 32
! @comments All exponential fuctions and their derivatives are involved
! only in termes 40-42 of the residual parts
! @authors Yu Fang
! @date    26-10-2016
!----------------------------------------------------------------------------------------------------
SUBROUTINE   deriv_expfonc( T, v, psi40, psi41, psi42,&
                            dpsi_ddelt40, dpsi_ddelt41, dpsi_ddelt42,&
                            dpsi_dtau40, dpsi_dtau41, dpsi_dtau42,&
                            ddpsi_dddelt40, ddpsi_dddelt41, ddpsi_dddelt42,&
                            ddpsi_ddtau40, ddpsi_ddtau41, ddpsi_ddtau42,&
                            ddpsi_ddeltdtau40, ddpsi_ddeltdtau41, ddpsi_ddeltdtau42 )





     USE def_constants

     IMPLICIT NONE

!IN/OUT

     REAL(pr) :: T,v
    
     REAL(pr) :: psi40, psi41, psi42

     REAL(pr) :: dpsi_ddelt40, dpsi_ddelt41, dpsi_ddelt42,&
                 dpsi_dtau40, dpsi_dtau41, dpsi_dtau42
  
     REAL(pr) :: ddpsi_dddelt40, ddpsi_dddelt41, ddpsi_dddelt42,&
                 ddpsi_ddtau40, ddpsi_ddtau41, ddpsi_ddtau42,&
                 ddpsi_ddeltdtau40, ddpsi_ddeltdtau41,ddpsi_ddeltdtau42

!LOCAL

    REAL(pr) :: tau, delt, rho
   
    REAL(pr) :: deltm1,deltm1_s, taum1,taum1_s

!Preparation of  variables

    rho        = 1_pr    / v
    delt       = 1_pr    / (rho_cr*v)
    tau        = T_cr    / T
!
!
    deltm1     = delt    - 1_pr
    deltm1_s   = deltm1  * deltm1
    taum1      = tau     - 1_pr
    taum1_s    = taum1   * taum1
!
!
    psi40 = exp(- CC40 * deltm1_s - DD40 * taum1_s)
    psi41 = exp(- CC41 * deltm1_s - DD40 * taum1_s)
    psi42 = exp(- CC42 * deltm1_s - DD40 * taum1_s)
!
!
!-----------------------------------------------------------------------------------------------------------
!Compute derivatives of psi 40-42
!-----------------------------------------------------------------------------------------------------------

    dpsi_ddelt40 = - 2_pr * CC40 * deltm1 * psi40

    dpsi_ddelt41 = - 2_pr * CC41 * deltm1 * psi41

    dpsi_ddelt42 = - 2_pr * CC42 * deltm1 * psi42
!
!
    dpsi_dtau40  = - 2_pr * DD40 * taum1 * psi40

    dpsi_dtau41  = - 2_pr * DD40 * taum1 * psi41

    dpsi_dtau42  = - 2_pr * DD40 * taum1 * psi42
!
!
    ddpsi_dddelt40 = (2_pr * CC40 * deltm1_s - 1_pr) * 2_pr * CC40 * psi40

    ddpsi_dddelt41 = (2_pr * CC41 * deltm1_s - 1_pr) * 2_pr * CC41 * psi41

    ddpsi_dddelt42 = (2_pr * CC42 * deltm1_s - 1_pr) * 2_pr * CC42 * psi42
!
!
    ddpsi_ddtau40 = (2_pr * DD40 * taum1_s - 1_pr) * 2_pr * DD40 * psi40

    ddpsi_ddtau41 = (2_pr * DD40 * taum1_s - 1_pr) * 2_pr * DD40 * psi41

    ddpsi_ddtau42 = (2_pr * DD40 * taum1_s - 1_pr) * 2_pr * DD40 * psi42
!
!
    ddpsi_ddeltdtau40 = 4_pr * CC40 * DD40 * deltm1 * taum1 * psi40

    ddpsi_ddeltdtau41 = 4_pr * CC41 * DD40 * deltm1 * taum1 * psi41

    ddpsi_ddeltdtau42 = 4_pr * CC42 * DD40 * deltm1 * taum1 * psi42
!
!
END SUBROUTINE deriv_expfonc
