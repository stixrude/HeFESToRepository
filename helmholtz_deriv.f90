!---------------------------------------------------------------------------------------------------
! All the formulations  have been taken from the original
! article,the Span-Wagner EoS for CO2.
! The EoS has been published in J. Phys. Chem. Ref. Data, Vol.25,
! pp. 1509-1596, 1996.
!---------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
! SUBROUTINE   helmholtz_deriv
! @brief   Compute derivatives of dimensionless helmotz energy p.1545 Table32:
!           phi_r_dtau,phi_ddelt,phi_r_ddtau,phi_r_dddelt,phi_r_ddeltdtau
! @routine deriv_disfonc.f90, deriv_expfonc.f90
! @authors Yu Fang
! @date    26-10-2016
!----------------------------------------------------------------------------------------------------
SUBROUTINE helmholtz_deriv( T, v, phi_r_ddelt, phi_r_dtau,&
                            phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau )
  
     USE def_constants

     IMPLICIT NONE

!IN/OUT

    REAL(pr) :: v, T
    REAL(pr) :: phi_r_ddelt, phi_r_dtau
    REAL(pr) :: phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau


!LOCAL
    REAL(pr) ::  tau,tau2,tau3,tau4,tau5,tau6,tau7,tau8,tau10,tau11,tau13,&
                 tau12,tau14,tau15,tau16,tau20,tau22,tau24,tau26,tau28,tau0_75,tau1_5,  &
                 tau2_5,taum1,taum1_s,tau21,tau23,tau27,tau0_75m1,tau1_5m1,&
                 rho,delt,delt2,delt3,delt4,delt5,delt6,delt7,         &
                 delt8,delt9,delt10,deltm1,deltm1_s,delt_m1
!
    REAL(pr) ::  tau_g2_35,  tau_g2_36, tau_g2_37, tau_g2_38,  tau_g2_39
!
!funciton distance trib
    REAL(pr) ::  trib40, trib41, trib42,&
                 dtrib_ddelt40, dtrib_ddelt41, dtrib_ddelt42,&
                 dtrib_dtau40, dtrib_dtau41, dtrib_dtau42,&
                 ddtrib_dddelt40, ddtrib_dddelt41, ddtrib_dddelt42,&
                 ddtrib_ddtau40,ddtrib_ddtau41, ddtrib_ddtau42,&
                 ddtrib_ddeltdtau40,ddtrib_ddeltdtau41,ddtrib_ddeltdtau42
!
!function exponential psi 
   REAL(pr) ::   psi40, psi41, psi42,&
                 dpsi_ddelt40, dpsi_ddelt41, dpsi_ddelt42,&
                 dpsi_dtau40, dpsi_dtau41, dpsi_dtau42,&
                 ddpsi_dddelt40, ddpsi_dddelt41, ddpsi_dddelt42,&
                 ddpsi_ddtau40, ddpsi_ddtau41, ddpsi_ddtau42,&
                 ddpsi_ddeltdtau40, ddpsi_ddeltdtau41,ddpsi_ddeltdtau42 
!
!
!
    REAL(pr) ::  sum1_phir_ddelt,sum2_phir_ddelt,sum3_phir_ddelt,sum4_phir_ddelt,&
                 sum1_phir_dtau,sum2_phir_dtau,sum3_phir_dtau,sum4_phir_dtau

    REAL(pr) ::  sum1_phir_dddelt,sum2_phir_dddelt,sum3_phir_dddelt,sum4_phir_dddelt,&
                 sum1_phir_ddtau,sum2_phir_ddtau,sum3_phir_ddtau,sum4_phir_ddtau,&
                 sum1_phir_ddeltdtau,sum2_phir_ddeltdtau,sum3_phir_ddeltdtau,sum4_phir_ddeltdtau
!--------------------------------------------------------------------------------------------------------
! Normalization of temperature and density and powers useful to speed up
! calculations
!-------------------------------------------------------------------------------------------------------
    tau     = T_cr / T
    tau2    = tau  * tau
    tau3    = tau2 * tau
    tau4    = tau2 * tau2
    tau5    = tau2 * tau3
    tau6    = tau3 * tau3
    tau7    = tau4 * tau3
    tau8    = tau4 * tau4
    tau10   = tau4 * tau6
    tau11   = tau8 * tau3
    tau12   = tau6 * tau6
    tau13   = tau6 * tau7
    tau14   = tau7 * tau7
    tau15   = tau7 * tau8
    tau16   = tau8 * tau8
    tau20   = tau10* tau10
    tau21   = tau13* tau8 
    tau22   = tau14* tau8
    tau23   = tau11* tau12
    tau24   = tau12* tau12
    tau26   = tau12* tau14
    tau27   = tau13* tau14
    tau28   = tau14* tau14
    tau0_75 = tau**(0.75_pr)
    tau1_5  = tau0_75 * tau0_75
    tau2_5  = tau1_5  * tau
    taum1   = tau  - 1.0_pr
    taum1_s = taum1 * taum1
    tau0_75m1 = tau0_75/tau
    tau1_5m1  = tau1_5/tau
!
!
!    rho      = 1.0_pr  / v
    delt     = 1.0_pr   / (rho_cr *v)
    delt_m1  = 1.0_pr  / delt
    delt2    = delt  * delt
    delt3    = delt2 * delt
    delt4    = delt2 * delt2
    delt5    = delt3 * delt2
    delt6    = delt3 * delt3
    delt7    = delt4 * delt3
    delt8    = delt4 * delt4
    delt9    = delt4 * delt5
    delt10   = delt5 * delt5
    deltm1   = delt  - 1.0_pr
    deltm1_s = deltm1 * deltm1
!
!
    tau_g2_35 = (tau-gam35)*(tau-gam35)
    tau_g2_36 = (tau-gam36)*(tau-gam36)
    tau_g2_37 = (tau-gam37)*(tau-gam37)
    tau_g2_38 = (tau-gam38)*(tau-gam38)
    tau_g2_39 = (tau-gam39)*(tau-gam39) 


    CALL deriv_disfonc ( T,v,trib40, trib41, trib42,&
                         dtrib_ddelt40, dtrib_ddelt41, dtrib_ddelt42,&
                         dtrib_dtau40, dtrib_dtau41, dtrib_dtau42,&
                         ddtrib_dddelt40, ddtrib_dddelt41, ddtrib_dddelt42,&
                         ddtrib_ddtau40,ddtrib_ddtau41, ddtrib_ddtau42,&
                         ddtrib_ddeltdtau40,ddtrib_ddeltdtau41,ddtrib_ddeltdtau42 )

    CALL deriv_expfonc( T, v, psi40, psi41, psi42,&
                        dpsi_ddelt40, dpsi_ddelt41, dpsi_ddelt42,&
                        dpsi_dtau40, dpsi_dtau41, dpsi_dtau42,&
                        ddpsi_dddelt40, ddpsi_dddelt41, ddpsi_dddelt42,&
                        ddpsi_ddtau40, ddpsi_ddtau41, ddpsi_ddtau42,&
                        ddpsi_ddeltdtau40, ddpsi_ddeltdtau41,ddpsi_ddeltdtau42 )

!-----------------------------------------------------------------------------------------------------
!Compute dphi_r/ddelt 
!-------------------------------------------------------------------------------------------------------
!sum1 terms 1-7

    sum1_phir_ddelt = n1 * d1                       +  n2 * d2         * tau0_75      + &
    &                 n3 * d3            * tau      +  n4 * d4         * tau2         + &
    &                 n5 * d5 * delt     * tau0_75  +  n6 * d6 * delt  * tau2         + &
    &                 n7 * d7 * delt2    * tau0_75

!sum2 terms 8-34
   
    sum2_phir_ddelt = n8         * tau1_5 *exp(-delt)  * (d8-delt)  + n9 * delt  * tau1_5 *exp(-delt) * (d9-delt)      +&
    &                 n10* delt3 * tau2_5 *exp(-delt)  * (d10-delt) + n11* delt4         *exp(-delt) * (d11-delt)      +&
    &                 n12* delt4 * tau1_5 *exp(-delt)  * (d12-delt) + n13* delt4 * tau2  *exp(-delt) * (d13-delt)      +&
    &                 n14* delt5          *exp(-delt)  * (d14-delt) + n15* delt5 * tau   *exp(-delt) * (d15-delt)      +&
    &                 n16* delt5 * tau2   *exp(-delt)  * (d16-delt) + n17        * tau3  *exp(-delt2)* (d17-2_pr*delt2)+&
    &                 n18        * tau6   *exp(-delt2) * (d18-2_pr*delt2) + &
    &                 n19* delt3 * tau3   *exp(-delt2) * (d19-2_pr*delt2) + &
    &                 n20* delt3 * tau6   *exp(-delt2) * (d20-2_pr*delt2) + &
    &                 n21* delt3 * tau8   *exp(-delt2) * (d21-2_pr*delt2) + &
    &                 n22* delt6 * tau6   *exp(-delt2) * (d22-2_pr*delt2) + & 
    &                 n23* delt7          *exp(-delt2) * (d23-2_pr*delt2) + &
    &                 n24* delt  * tau7   *exp(-delt3) * (d24-3_pr*delt3) + & 
    &                 n25* delt2 * tau12  *exp(-delt3) * (d25-3_pr*delt3) + &
    &                 n26* delt2 * tau16  *exp(-delt3) * (d26-3_pr*delt3) + &
    &                 n27* delt4 * tau22  *exp(-delt4) * (d27-4_pr*delt4) + &
    &                 n28* delt4 * tau24  *exp(-delt4) * (d28-4_pr*delt4) + &
    &                 n29* delt5 * tau16  *exp(-delt4) * (d29-4_pr*delt4) + &
    &                 n30* delt6 * tau24  *exp(-delt4) * (d30-4_pr*delt4) + &
    &                 n31* delt7 * tau8   *exp(-delt4) * (d31-4_pr*delt4) + &
    &                 n32* delt9 * tau2   *exp(-delt4) * (d32-4_pr*delt4) + &
    &                 n33* delt3 * tau28  *exp(-delt5) * (d33-5_pr*delt5) + &
    &                 n34* delt7 * tau14  *exp(-delt6) * (d34-6_pr*delt6)

!sum3 terms 35-39

    sum3_phir_ddelt =  n35 * delt2 * tau  *exp(-alp35*deltm1_s - bet35*tau_g2_35)  *(d35/delt - 2_pr*alp35*deltm1) +  &
    &                  n36 * delt2        *exp(-alp36*deltm1_s - bet36*tau_g2_36)  *(d36/delt - 2_pr*alp36*deltm1) +  &
    &                  n37 * delt2 * tau  *exp(-alp37*deltm1_s - bet37*tau_g2_37)  *(d37/delt - 2_pr*alp37*deltm1) +  &
    &                  n38 * delt3 * tau3 *exp(-alp38*deltm1_s - bet38*tau_g2_38)  *(d38/delt - 2_pr*alp38*deltm1) +  &
    &                  n39 * delt3 * tau3 *exp(-alp39*deltm1_s - bet39*tau_g2_39)  *(d39/delt - 2_pr*alp39*deltm1)

!sum4 terms 40-42 

    sum4_phir_ddelt = n40 * ( trib40 *(psi40 + delt*dpsi_ddelt40) + dtrib_ddelt40*delt*psi40  ) + &
    &                 n41 * ( trib41 *(psi41 + delt*dpsi_ddelt41) + dtrib_ddelt41*delt*psi41  ) + &
    &                 n42 * ( trib42 *(psi42 + delt*dpsi_ddelt42) + dtrib_ddelt42*delt*psi42  ) 
!
!   
    
!    print*, sum1_phir_ddelt, sum2_phir_ddelt,sum3_phir_ddelt,sum4_phir_ddelt
!    print*, trib40,psi40,dpsi_ddelt40,dtrib_ddelt40
    phi_r_ddelt = sum1_phir_ddelt + sum2_phir_ddelt + sum3_phir_ddelt + sum4_phir_ddelt
!
!
!-----------------------------------------------------------------------------------------------------------------------------
!Compute dphi_r/dtau
!-----------------------------------------------------------------------------------------------------------------------------
!sum1 terms 1-7

    sum1_phir_dtau =                                    n2 * delt * 0.75_pr * tau0_75m1      + &
    &                 n3 * delt                       + n4 * delt * 2.0_pr  * tau            + &
    &                 n5 * delt2* 0.75_pr * tau0_75m1 + n6 * delt2* 2.0_pr  * tau            + &
    &                 n7 * delt3* 0.75_pr * tau0_75m1
!
!
!sum2 terms 8-34
  
    sum2_phir_dtau = 0_pr
    sum2_phir_dtau =   n8 * delt  * 1.5_pr * tau1_5m1 *exp(-delt)  + n9 * delt2 * 1.5_pr * tau1_5m1 *exp(-delt) +&
    &                  n10* delt4 * 2.5_pr * tau1_5   *exp(-delt)                                               +&
    &                  n12* delt5 * 1.5_pr * tau1_5m1 *exp(-delt)  + n13* delt5 * 2.0_pr * tau      *exp(-delt) +&
    &                                                                n15* delt6                     *exp(-delt) +&
    &                  n16* delt6 * 2.0_pr * tau      *exp(-delt)  + n17* delt  * 3.0_pr * tau2     *exp(-delt2)+&
    &                  n18* delt  * 6.0_pr * tau5     *exp(-delt2) + &
    &                  n19* delt4 * 3.0_pr * tau2     *exp(-delt2)+& 
    &                  n20* delt4 * 6.0_pr * tau5     *exp(-delt2) + n21* delt4 * 8.0_pr * tau7     *exp(-delt2)+&
    &                  n22* delt7 * 6.0_pr * tau5     *exp(-delt2)                                              +&
    &                  n24* delt2 * 7.0_pr * tau6     *exp(-delt3) + n25* delt3 * 12.0_pr * tau11    *exp(-delt3)+&
    &                  n26* delt3 * 16.0_pr* tau15    *exp(-delt3) + n27* delt5 * 22.0_pr * tau21    *exp(-delt4)+&
    &                  n28* delt5 * 24.0_pr* tau23    *exp(-delt4) + n29* delt6 * 16.0_pr * tau15    *exp(-delt4)+&
    &                  n30* delt7 * 24.0_pr* tau23    *exp(-delt4) + n31* delt8 * 8.0_pr  * tau7     *exp(-delt4)+&
    &                  n32* delt10* 2.0_pr * tau      *exp(-delt4) + n33* delt4 * 28.0_pr *tau27     *exp(-delt5)+&
    &                  n34* delt8 * 14.0_pr* tau13    *exp(-delt6)
!
!        print*,' tau5  ', tau5  
!sum2 terms 35-39
  
    sum3_phir_dtau =  n35 * delt2 * tau  *exp(-alp35*deltm1_s - bet35*tau_g2_35) * ( 1.0_pr/tau - 2.0_pr * bet35 * (tau-gam35) ) + &
    &                 n36 * delt2        *exp(-alp36*deltm1_s - bet36*tau_g2_36) * (            - 2.0_pr * bet36 * (tau-gam36) ) + &
    &                 n37 * delt2 * tau  *exp(-alp37*deltm1_s - bet37*tau_g2_37) * ( 1.0_pr/tau - 2.0_pr * bet37 * (tau-gam37) ) + &
    &                 n38 * delt3 * tau3 *exp(-alp38*deltm1_s - bet38*tau_g2_38) * ( 3.0_pr/tau - 2.0_pr * bet38 * (tau-gam38) ) + &
    &                 n39 * delt3 * tau3 *exp(-alp39*deltm1_s - bet39*tau_g2_39) * ( 3.0_pr/tau - 2.0_pr * bet39 * (tau-gam39) )
!sum3 terms 40-42

    sum4_phir_dtau = n40 * delt * ( dtrib_dtau40*psi40 + trib40*dpsi_dtau40 ) + &
    &                 n41 * delt * ( dtrib_dtau41*psi41 + trib41*dpsi_dtau41 ) + &
    &                 n42 * delt * ( dtrib_dtau42*psi42 + trib42*dpsi_dtau42 ) 
!
!
    phi_r_dtau = sum1_phir_dtau +  sum2_phir_dtau + sum3_phir_dtau + sum4_phir_dtau
!     phi_r_dtau = sum1_phir_dtau + sum3_phir_dtau + sum4_phir_dtau
!print*,  'sum2 = ' , sum2_phir_dtau 
!
!
!-----------------------------------------------------------------------------------------------------------------------------------
!Compute ddphi_r / dddelt
!-------------------------------------------------------------------------------------------------------------------------------------
!sum1 terms 1-7 
 
    sum1_phir_dddelt =  n5 * d5                     * tau0_75  +  n6 * d6   * tau2  + &
    &                   n7 * d7 * (d7-1_pr) * delt  * tau0_75
!
!sum2 terms 8-34
    sum2_phir_dddelt =n8 * exp(-delt)  * ( tau1_5*delt**(d8-2_pr) * ((d8-delt) *(d8-1_pr-delt)  - delt ) )       + & 
    &                 n9 * exp(-delt)  * ( tau1_5*delt**(d9-2_pr) * ((d9-delt) *(d9-1_pr-delt)  - delt ) )       + &
    &                 n10* exp(-delt)  * ( tau2_5*delt**(d10-2_pr)* ((d10-delt)*(d10-1_pr-delt) - delt ) )       + &
    &                 n11* exp(-delt)  * (        delt**(d11-2_pr)* ((d11-delt)*(d11-1_pr-delt) - delt ) )       + &
    &                 n12* exp(-delt)  * ( tau1_5*delt**(d12-2_pr)* ((d12-delt)*(d12-1_pr-delt) - delt ) )       + &
    &                 n13* exp(-delt)  * ( tau2  *delt**(d13-2_pr)* ((d13-delt)*(d13-1_pr-delt) - delt ) )       + &
    &                 n14* exp(-delt)  * (        delt**(d14-2_pr)* ((d14-delt)*(d14-1_pr-delt) - delt ) )       + &
    &                 n15* exp(-delt)  * ( tau   *delt**(d15-2_pr)* ((d15-delt)*(d15-1_pr-delt) - delt ) )       + &
    &                 n16* exp(-delt)  * ( tau2  *delt**(d16-2_pr)* ((d16-delt)*(d16-1_pr-delt) - delt ) )       + &
                      n17* exp(-delt2) * ( tau3  *delt**(d17-2_pr)* ((d17-2_pr*delt2)*(d17-1_pr-2_pr*delt2) - 4_pr*delt2 )) +&
    &                 n18* exp(-delt2) * ( tau6  *delt**(d18-2_pr)* ((d18-2_pr*delt2)*(d18-1_pr-2_pr*delt2) - 4_pr*delt2 )) +&
    &                 n19* exp(-delt2) * ( tau3  *delt**(d19-2_pr)* ((d10-2_pr*delt2)*(d19-1_pr-2_pr*delt2) - 4_pr*delt2 )) +&
    &                 n20* exp(-delt2) * ( tau6  *delt**(d20-2_pr)* ((d20-2_pr*delt2)*(d20-1_pr-2_pr*delt2) - 4_pr*delt2 )) +&
    &                 n21* exp(-delt2) * ( tau8  *delt**(d21-2_pr)* ((d21-2_pr*delt2)*(d21-1_pr-2_pr*delt2) - 4_pr*delt2 )) +&
    &                 n22* exp(-delt2) * ( tau6  *delt**(d22-2_pr)* ((d22-2_pr*delt2)*(d22-1_pr-2_pr*delt2) - 4_pr*delt2 )) +&
    &                 n23* exp(-delt2) * (        delt**(d23-2_pr)* ((d23-2_pr*delt2)*(d23-1_pr-2_pr*delt2) - 4_pr*delt2 )) +&
    &                 n24* exp(-delt3) * ( tau7  *delt**(d24-2_pr)* ((d24-3_pr*delt3)*(d24-1_pr-3_pr*delt3) - 9_pr*delt3 )) +&
    &                 n25* exp(-delt3) * ( tau12 *delt**(d25-2_pr)* ((d25-3_pr*delt3)*(d25-1_pr-3_pr*delt3) - 9_pr*delt3 )) +&
    &                 n26* exp(-delt3) * ( tau16 *delt**(d26-2_pr)* ((d26-3_pr*delt3)*(d26-1_pr-3_pr*delt3) - 9_pr*delt3 )) +&
    &                 n27* exp(-delt4) * ( tau22 *delt**(d27-2_pr)* ((d27-4_pr*delt4)*(d27-1_pr-4_pr*delt4) - 16_pr*delt4 )) +&
    &                 n28* exp(-delt4) * ( tau24 *delt**(d28-2_pr)* ((d28-4_pr*delt4)*(d28-1_pr-4_pr*delt4) - 16_pr*delt4 )) +&
    &                 n29* exp(-delt4) * ( tau16 *delt**(d29-2_pr)* ((d29-4_pr*delt4)*(d29-1_pr-4_pr*delt4) - 16_pr*delt4 )) +&
    &                 n30* exp(-delt4) * ( tau24 *delt**(d30-2_pr)* ((d30-4_pr*delt4)*(d30-1_pr-4_pr*delt4) - 16_pr*delt4 )) +&
    &                 n31* exp(-delt4) * ( tau8  *delt**(d31-2_pr)* ((d31-4_pr*delt4)*(d31-1_pr-4_pr*delt4) - 16_pr*delt4 )) +&
    &                 n32* exp(-delt4) * ( tau2  *delt**(d32-2_pr)* ((d32-4_pr*delt4)*(d32-1_pr-4_pr*delt4) - 16_pr*delt4 )) +&
    &                 n33* exp(-delt5) * ( tau28 *delt**(d33-2_pr)* ((d33-5_pr*delt5)*(d33-1_pr-5_pr*delt5) - 25_pr*delt5 )) +&
    &                 n34* exp(-delt6) * ( tau14 *delt**(d34-2_pr)* ((d34-6_pr*delt6)*(d34-1_pr-6_pr*delt6) - 36_pr*delt6 ))
!
!sum3 terms 35-39
   
    sum3_phir_dddelt = n35  * tau  *exp(-alp35*deltm1_s - bet35*tau_g2_35)*  &
    &                  ( -2_pr*alp35*delt2 + 4_pr*alp35*alp35*delt2*deltm1_s - 4_pr*d35*alp35*delt*deltm1 + d35*(d35-1_pr) ) + &
    &                  n36         *exp(-alp36*deltm1_s - bet36*tau_g2_36)*  &
    &                  ( -2_pr*alp36*delt2 + 4_pr*alp36*alp36*delt2*deltm1_s - 4_pr*d36*alp36*delt*deltm1 + d36*(d36-1_pr) ) + &
    &                  n37  * tau  *exp(-alp37*deltm1_s - bet37*tau_g2_37)*  &
    &                  ( -2_pr*alp37*delt2 + 4_pr*alp37*alp37*delt2*deltm1_s - 4_pr*d37*alp37*delt*deltm1 + d37*(d37-1_pr) ) + &
    &                  n38  * tau3 *exp(-alp38*deltm1_s - bet38*tau_g2_38)*  &
    &                  ( -2_pr*alp38*delt3 + 4_pr*alp38*alp38*delt3*deltm1_s - 4_pr*d38*alp38*delt2*deltm1+ d38*(d38-1_pr)*delt )+&
    &                  n39  * tau3 *exp(-alp39*deltm1_s - bet39*tau_g2_39)*  &
                       ( -2_pr*alp39*delt3 + 4_pr*alp39*alp39*delt3*deltm1_s - 4_pr*d39*alp39*delt2*deltm1+ d39*(d39-1_pr)*delt )
!
!
!sum4 terms 40-42
   
    sum4_phir_dddelt = n40 * ( trib40*(2_pr*dpsi_ddelt40 + delt*ddpsi_dddelt40) + 2_pr*dtrib_ddelt40*(psi40 + delt*dpsi_ddelt40) + &
    &                           ddtrib_dddelt40*delt*psi40 ) +&
    &                  n41 * ( trib41*(2_pr*dpsi_ddelt41 + delt*ddpsi_dddelt41) + 2_pr*dtrib_ddelt41*(psi41 + delt*dpsi_ddelt41) + &
    &                           ddtrib_dddelt41*delt*psi41 ) +&
    &                  n42 * ( trib42*(2_pr*dpsi_ddelt42 + delt*ddpsi_dddelt42) + 2_pr*dtrib_ddelt42*(psi42 + delt*dpsi_ddelt42) + &
                                ddtrib_dddelt42*delt*psi42 )
!
!
    phi_r_dddelt = sum1_phir_dddelt + sum2_phir_dddelt + sum3_phir_dddelt + sum4_phir_dddelt
!
!
!-------------------------------------------------------------------------------------------------------------------------------------------------------------
!Compute dphi_r / ddtau
!-------------------------------------------------------------------------------------------------------------------------------------------------------------
!sum1 terms 1-7
   
     sum1_phir_ddtau = n2 * delt * 0.75_pr * (0.75_pr-1_pr)*tau**(0.75_pr-2_pr)  + &
    &                  n4 * delt * 2.0_pr                                        + &
    &                  n5 * delt2* 0.75_pr * (0.75_pr-1_pr)*tau**(0.75_pr-2_pr)  + &
    &                  n6 * delt2 * 2.0_pr                                       + &
    &                  n7 * delt3* 0.75_pr * (0.75_pr-1_pr)*tau**(0.75_pr-2_pr)
!
!
!sum2 terms 8-34
   
    sum2_phir_ddtau =  n8 * delt  * 1.5_pr *0.5_pr  * tau**(1.5_pr-2_pr) *exp(-delt)  +& 
    &                  n9 * delt2 * 1.5_pr *0.5_pr  * tau**(1.5_pr-2_pr) *exp(-delt)  +&
    &                  n10* delt4 * 2.5_pr *1.5_pr  * tau1_5m1           *exp(-delt)  +&
    &                  n12* delt5 * 1.5_pr *0.5_pr  * tau**(1.5_pr-2_pr) *exp(-delt)  +& 
    &                  n13* delt5 * 2.0_pr                               *exp(-delt)  +&
    &                  n16* delt6 * 2.0_pr                               *exp(-delt)  +& 
    &                  n17* delt  * 3.0_pr *2_pr    * tau                *exp(-delt2) +&
    &                  n18* delt  * 6.0_pr *5_pr    * tau4               *exp(-delt2) +& 
    &                  n19* delt4 * 3.0_pr *2_pr    * tau                *exp(-delt2 )+&
    &                  n20* delt4 * 6.0_pr *5_pr    * tau4               *exp(-delt2) +& 
    &                  n21* delt4 * 8.0_pr *7_pr    * tau6               *exp(-delt2) +&
    &                  n22* delt7 * 6.0_pr *5_pr    * tau4               *exp(-delt2) +&
    &                  n24* delt2 * 7.0_pr *6_pr    * tau5               *exp(-delt3) +& 
    &                  n25* delt3 * 12.0_pr*11_pr   * tau10              *exp(-delt3) +&
    &                  n26* delt3 * 16.0_pr*15_pr   * tau14              *exp(-delt3) +& 
    &                  n27* delt5 * 22.0_pr*21_pr   * tau20              *exp(-delt4) +&
    &                  n28* delt5 * 24.0_pr*23_pr   * tau22              *exp(-delt4) +& 
    &                  n29* delt6 * 16.0_pr*15_pr   * tau14              *exp(-delt4) +&
    &                  n30* delt7 * 24.0_pr*23_pr   * tau22              *exp(-delt4) +&
    &                  n31* delt8 * 8.0_pr *7_pr    * tau6               *exp(-delt4) +&
    &                  n32* delt10* 2.0_pr                               *exp(-delt4) +& 
    &                  n33* delt4 * 28.0_pr*27_pr   * tau26              *exp(-delt5) +&
    &                  n34* delt8 * 14.0_pr*13_pr   * tau12              *exp(-delt6)
!
!
!sum3 terms 35-39
 
    sum3_phir_ddtau = n35 * delt2 * tau  *exp(-alp35*deltm1_s - bet35*tau_g2_35) * &
    &                        ( (1.0_pr/tau - 2.0_pr * bet35 * (tau-gam35))**(2_pr) - 1_pr/tau2 - 2_pr*bet35 )   +& 
    &                 n36 * delt2        *exp(-alp36*deltm1_s - bet36*tau_g2_36) * &
    &                        ( (           - 2.0_pr * bet36 * (tau-gam36))**(2_pr)             - 2_pr*bet36)    +&
    &                 n37 * delt2 * tau  *exp(-alp37*deltm1_s - bet37*tau_g2_37) * &
    &                        ( (1.0_pr/tau - 2.0_pr * bet37 * (tau-gam37))**(2_pr) - 1_pr/tau2 - 2_pr*bet37  )  +&
    &                 n38 * delt3 * tau3 *exp(-alp38*deltm1_s - bet38*tau_g2_38) * &
    &                        ( (3.0_pr/tau - 2.0_pr * bet38 * (tau-gam38))**(2_pr) - 3_pr/tau2 - 2_pr*bet38   ) +&
    &                 n39 * delt3 * tau3 *exp(-alp39*deltm1_s - bet39*tau_g2_39) * &
    &                        ( (3.0_pr/tau - 2.0_pr * bet39 * (tau-gam39))**(2_pr) - 3_pr/tau2 - 2_pr*bet39  )
    
!    print*, sum1_phir_ddtau,sum2_phir_ddtau,sum3_phir_ddtau,sum4_phir_ddtau
!
!
!sum4 terms 40-42
    sum4_phir_ddtau = n40 * delt * ( ddtrib_ddtau40*psi40 + 2_pr*dtrib_dtau40*dpsi_dtau40 + trib40*ddpsi_ddtau40 ) + &
                      n41 * delt * ( ddtrib_ddtau41*psi41 + 2_pr*dtrib_dtau41*dpsi_dtau41 + trib41*ddpsi_ddtau41 ) + &
                      n42 * delt * ( ddtrib_ddtau42*psi42 + 2_pr*dtrib_dtau42*dpsi_dtau42 + trib42*ddpsi_ddtau42 ) 
!
!
    phi_r_ddtau = sum1_phir_ddtau + sum2_phir_ddtau + sum3_phir_ddtau + sum4_phir_ddtau
!
!
!------------------------------------------------------------------------------------------------------------------------------------------
!Compute ddphi_r / ddeltdtau
!------------------------------------------------------------------------------------------------------------------------------------------
!sum1 terms 1-7

    sum1_phir_ddeltdtau =  n2 * 0.75_pr  * tau0_75m1                                                + &
    &                      n3                                     +  n4 * 2_pr  * tau               + &
    &                      n5 * 2_pr * 0.75_pr * delt * tau0_75m1 +  n6 * 2_pr  *2_pr * delt * tau  + &
    &                      n7 * 3_pr * 0.75_pr * delt2* tau0_75m1

!
!sum2 terms 8-34

    sum2_phir_ddeltdtau = n8        * 1.5_pr * tau1_5m1  *exp(-delt)  * (d8-delt)  +&
    &                     n9 * delt * 1.5_pr * tau1_5m1  *exp(-delt)  * (d9-delt)  +&
    &                     n10* delt3* 2.5_pr * tau1_5    *exp(-delt)  * (d10-delt) +&
    &                     n12* delt4* 1.5_pr * tau1_5m1  *exp(-delt)  * (d12-delt) +&
    &                     n13* delt4* 2.0_pr * tau       *exp(-delt)  * (d13-delt) +&
    &                     n15* delt5                     *exp(-delt)  * (d15-delt) +&
    &                     n16* delt5* 2.0_pr * tau       *exp(-delt)  * (d16-delt) +&
    &                     n17       * 3.0_pr * tau2      *exp(-delt2) * (d17-2_pr*delt2)+&
    &                     n18       * 6.0_pr * tau5      *exp(-delt2) * (d18-2_pr*delt2)+&
    &                     n19* delt3* 3.0_pr * tau2      *exp(-delt2) * (d19-2_pr*delt2)+&
    &                     n20* delt3* 6.0_pr * tau5      *exp(-delt2) * (d20-2_pr*delt2)+&
    &                     n21* delt3* 8.0_pr * tau7      *exp(-delt2) * (d21-2_pr*delt2)+&
    &                     n22* delt6* 6.0_pr * tau5      *exp(-delt2) * (d22-2_pr*delt2)+&
    &                     n24* delt * 7.0_pr * tau6      *exp(-delt3) * (d24-3_pr*delt3)+&
    &                     n25* delt2* 12.0_pr* tau11     *exp(-delt3) * (d25-3_pr*delt3)+&
    &                     n26* delt2* 16.0_pr* tau15     *exp(-delt3) * (d26-3_pr*delt3)+&
    &                     n27* delt4* 22.0_pr* tau21     *exp(-delt4) * (d27-4_pr*delt4)+&
    &                     n28* delt4* 24.0_pr* tau23     *exp(-delt4) * (d28-4_pr*delt4)+&
    &                     n29* delt5* 16.0_pr* tau15     *exp(-delt4) * (d29-4_pr*delt4)+&
    &                     n30* delt6* 24.0_pr* tau23     *exp(-delt4) * (d30-4_pr*delt4)+&
    &                     n31* delt7* 8.0_pr * tau7      *exp(-delt4) * (d31-4_pr*delt4)+&
    &                     n32* delt9* 2.0_pr * tau       *exp(-delt4) * (d32-4_pr*delt4)+&
    &                     n33* delt3* 28.0_pr* tau27     *exp(-delt5) * (d33-5_pr*delt5)+&
    &                     n34* delt7* 14.0_pr* tau13     *exp(-delt6) * (d34-6_pr*delt6)
!
!sum3 terms 35-39

    sum3_phir_ddeltdtau = n35 * delt2 * tau  *exp(-alp35*deltm1_s - bet35*tau_g2_35)  *&
    &                         (d35/delt - 2_pr*alp35*deltm1) * (1_pr/tau - 2_pr*bet35*(tau-gam35)) +  &
    &                     n36 * delt2        *exp(-alp36*deltm1_s - bet36*tau_g2_36) *&
    &                         (d36/delt - 2_pr*alp36*deltm1) * (         - 2_pr*bet36*(tau-gam36)) +  &
    &                     n37 * delt2 * tau  *exp(-alp37*deltm1_s - bet37*tau_g2_37) *&
    &                         (d37/delt - 2_pr*alp37*deltm1) * (1_pr/tau - 2_pr*bet37*(tau-gam37)) +  &
    &                     n38 * delt3 * tau3 *exp(-alp38*deltm1_s - bet38*tau_g2_38) *&
    &                         (d38/delt - 2_pr*alp38*deltm1) * (3_pr/tau - 2_pr*bet38*(tau-gam38)) +  &
    &                     n39 * delt3 * tau3 *exp(-alp39*deltm1_s - bet39*tau_g2_39)*&
                              (d39/delt - 2_pr*alp39*deltm1) * (3_pr/tau - 2_pr*bet39*(tau-gam39))
!
!
!sum4 terms 40-42
    
    sum4_phir_ddeltdtau = n40 * ( trib40*(dpsi_dtau40+delt*ddpsi_ddeltdtau40) + delt*dtrib_ddelt40*dpsi_dtau40 +&
    &                             dtrib_dtau40*(psi40+delt*dpsi_ddelt40)      + ddtrib_ddeltdtau40*delt*psi40 )        +&
    &                     n41 * ( trib41*(dpsi_dtau41+delt*ddpsi_ddeltdtau41) + delt*dtrib_ddelt41*dpsi_dtau41 +&
    &                             dtrib_dtau41*(psi41+delt*dpsi_ddelt41)      + ddtrib_ddeltdtau41*delt*psi41 )        +&
    &                     n42 * ( trib42*(dpsi_dtau42+delt*ddpsi_ddeltdtau42) + delt*dtrib_ddelt42*dpsi_dtau42 +&
    &                             dtrib_dtau42*(psi42+delt*dpsi_ddelt42)      + ddtrib_ddeltdtau42*delt*psi42 )     
    
!
!
    phi_r_ddeltdtau =  sum1_phir_ddeltdtau + sum2_phir_ddeltdtau + sum3_phir_ddeltdtau + sum4_phir_ddeltdtau
!
!
END SUBROUTINE helmholtz_deriv
