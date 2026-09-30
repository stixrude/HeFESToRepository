!---------------------------------------------------------------------------------------------------
! All the formulations  have been taken from the original
! article,the Span-Wagner EoS for CO2.
! The EoS has been published in J. Phys. Chem. Ref. Data, Vol.25,
! pp. 1509-1596, 1996.
!---------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
! SUBROUTINE   helmholtz_dimless
! @brief   Compute dimensionless helmotz energy:  phi_0 (perfect gas part) and  phi_r (residul part)
!                                                 derivative of phi_0
! @authors Marco De Lorenzo, Yu Fang
! @date    25-10-2016
!----------------------------------------------------------------------------------------------------
SUBROUTINE helmholtz_dimless (T,v,phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                              phi0_dddelt,phi0_ddtau)

    USE def_constants
       
    IMPLICIT NONE

!IN/OUT
    REAL(pr) :: v, T
    REAL(pr) :: phi_0,phi_r
    REAL(pr) :: phi0_ddelt
    REAL(pr) :: phi0_dtau,phi0_dddelt,phi0_ddtau                    !,phi_0_ddeltdtau

!Local
    REAL(pr) ::  tau,   tau2,  tau3,  tau4,  tau6,  tau7,    tau8,           &
                tau12,  tau14, tau16, tau22, tau24, tau28, tau0_75, tau1_5,  &
                tau2_5, taum1, taum1_s,                                      &
                rho,    del,   del2,  del3,  del4,  del5,  del6,    del7,    &
                del8,   del10, delm1, delm1_s
!
    REAL(pr) :: phi_r_sum1, phi_r_sum2, phi_r_sum3, phi_r_sum4,              & 
                tau_g2_35,  tau_g2_36, tau_g2_37, tau_g2_38,  tau_g2_39,     &
                braces,     braces_s


!--------------------------------------------------------------------------------------------------------
! Normalization of temperature and density and powers useful to speed up calculations
! di in Table31 is computed in this part 
!-------------------------------------------------------------------------------------------------------
    tau     = T_cr / T
    tau2    = tau  * tau
    tau3    = tau2 * tau
    tau4    = tau2 * tau2
    tau6    = tau3 * tau3
    tau7    = tau4 * tau3
    tau8    = tau4 * tau4
    tau12   = tau6 * tau6
    tau14   = tau7 * tau7
    tau16   = tau8 * tau8
    tau22   = tau14* tau8
    tau24   = tau12* tau12
    tau28   = tau14* tau14
    tau0_75 = tau**0.75_pr
    tau1_5  = tau0_75 * tau0_75
    tau2_5  = tau1_5  * tau
    taum1   = tau  - 1_pr
    taum1_s = taum1 * taum1
!
!
!    rho     = 1_pr  / v
    del     = 1.0_pr / (rho_cr*v)
    del2    = del  * del
    del3    = del2 * del
    del4    = del2 * del2
    del5    = del3 * del2
    del6    = del3 * del3
    del7    = del4 * del3
    del8    = del4 * del4
    del10   = del5 * del5
    delm1   = del  - 1_pr
    delm1_s = delm1 * delm1
!
!-----------------------------------------------------------------------------------------------------
!
!Compute the perfect gas part of helhmoltz energy (p.1543 Eq.6.3)
!
!-----------------------------------------------------------------------------------------------------

    phi_0 = a1_0  +  a2_0 * tau + a3_0 * log(tau)  + log(del)+ &     !!!
    &       a4_0 * log(1_pr - exp(-tau*th4_0))                  +& 
    &       a5_0 * log(1_pr - exp(-tau*th5_0))                  +& 
    &       a6_0 * log(1_pr - exp(-tau*th6_0))                   + &
    &       a7_0 * log(1_pr - exp(-tau*th7_0))                   + &
    &       a8_0 * log(1_pr - exp(-tau*th8_0))

!-----------------------------------------------------------------------------------------------------
!
!Compute the  derivatives of the gas part of helhmoltz energy (p.1541 Table28)
!
!-----------------------------------------------------------------------------------------------------
    phi0_ddelt  = 1_pr/del
!
!
    phi0_dddelt = -1_pr/del2
!
!
!   phi_0_ddeltdtau 
!
     phi0_dtau  =  a2_0 + a3_0 / tau +&
     &             a4_0 * th4_0* ( (1.0_pr - exp(-tau*th4_0))**(-1.0_pr) - 1.0_pr )  + &
     &             a5_0 * th5_0* ( (1.0_pr - exp(-tau*th5_0))**(-1.0_pr) - 1.0_pr )  + &
     &             a6_0 * th6_0* ( (1.0_pr - exp(-tau*th6_0))**(-1.0_pr) - 1.0_pr )  + & 
     &             a7_0 * th7_0* ( (1.0_pr - exp(-tau*th7_0))**(-1.0_pr) - 1.0_pr )  + &
     &             a8_0 * th8_0* ( (1.0_pr - exp(-tau*th8_0))**(-1.0_pr) - 1.0_pr )
!
!
    phi0_ddtau = -a3_0/tau2 - (&
    &            a4_0 * th4_0*th4_0 * exp(-tau*th4_0) *  (1_pr - exp(-tau*th4_0))**(-2_pr)    + &
    &            a5_0 * th5_0*th5_0 * exp(-tau*th5_0) *  (1_pr - exp(-tau*th5_0))**(-2_pr)    + &
    &            a6_0 * th6_0*th6_0 * exp(-tau*th6_0) *  (1_pr - exp(-tau*th6_0))**(-2_pr)    + &
    &            a7_0 * th7_0*th7_0 * exp(-tau*th7_0) *  (1_pr - exp(-tau*th7_0))**(-2_pr)    + &
    &            a8_0 * th8_0*th8_0 * exp(-tau*th8_0) *  (1_pr - exp(-tau*th8_0))**(-2_pr))
!
!
!
!-----------------------------------------------------------------------------------------------------
!
!Compute the residul part of helhmoltz energy (p.1544 Eq.6.5)
!
!-----------------------------------------------------------------------------------------------------
! It is composed by 4 summations, here splitted in 4 different variables
! for easiness.
!
    phi_r_sum1 = n1 * del            +  n2 * del * tau0_75      + &
    &            n3 * del * tau      +  n4 * del * tau2         + &
    &            n5 * del2* tau0_75  +  n6 * del2* tau2         + &
    &            n7 * del3* tau0_75
!
!
    phi_r_sum2 =                                                      &
    & n8 * del  * tau1_5 *exp(-del)  + n9 * del2 * tau1_5 *exp(-del) +&
    & n10* del4 * tau2_5 *exp(-del)  + n11* del5          *exp(-del) +&
    & n12* del5 * tau1_5 *exp(-del)  + n13* del5 * tau2   *exp(-del) +&
    & n14* del6          *exp(-del)  + n15* del6 * tau    *exp(-del) +&
    & n16* del6 * tau2   *exp(-del)  + n17* del  * tau3   *exp(-del2)+&
    & n18* del  * tau6   *exp(-del2) + n19* del4 * tau3   *exp(-del2)+&
    & n20* del4 * tau6   *exp(-del2) + n21* del4 * tau8   *exp(-del2)+&
    & n22* del7 * tau6   *exp(-del2) + n23* del8          *exp(-del2)+&
    & n24* del2 * tau7   *exp(-del3) + n25* del3 * tau12  *exp(-del3)+&
    & n26* del3 * tau16  *exp(-del3) + n27* del5 * tau22  *exp(-del4)+&
    & n28* del5 * tau24  *exp(-del4) + n29* del6 * tau16  *exp(-del4)+&
    & n30* del7 * tau24  *exp(-del4) + n31* del8 * tau8   *exp(-del4)+&
    & n32* del10* tau2   *exp(-del4) + n33* del4 * tau28  *exp(-del5)+&
    & n34* del8 * tau14  *exp(-del6)
!
!
    tau_g2_35 = (tau-gam35)*(tau-gam35)
    tau_g2_36 = (tau-gam36)*(tau-gam36)
    tau_g2_37 = (tau-gam37)*(tau-gam37)
    tau_g2_38 = (tau-gam38)*(tau-gam38)
    tau_g2_39 = (tau-gam39)*(tau-gam39)
!
    phi_r_sum3 =                                                   &
    & n35 * del2 * tau  *exp(-alp35*delm1_s - bet35*tau_g2_35)  +  &
    & n36 * del2        *exp(-alp36*delm1_s - bet36*tau_g2_36)  +  &
    & n37 * del2 * tau  *exp(-alp37*delm1_s - bet37*tau_g2_37)  +  &
    & n38 * del3 * tau3 *exp(-alp38*delm1_s - bet38*tau_g2_38)  +  &
    & n39 * del3 * tau3 *exp(-alp39*delm1_s - bet39*tau_g2_39)
!
!
    braces   = - taum1 + AA40 * delm1_s**(0.5_pr/bet40)
    braces_s = braces * braces
!
!
    phi_r_sum4 =                                           &
    & n40 * del *exp(-CC40*delm1_s - DD40*taum1_s)         &
    &           *   (braces_s + BB40*delm1_s**a40)**b40  + &
    & n41 * del *exp(-CC41*delm1_s - DD40*taum1_s)         &
    &           *   (braces_s + BB41*delm1_s**a41)**b41  + &
    & n42 * del *exp(-CC42*delm1_s - DD40*taum1_s)         &
    &           *   (braces_s + BB42*delm1_s**a42)**b42
!
!
    phi_r =  phi_r_sum1 + phi_r_sum2 + phi_r_sum3 + phi_r_sum4
!
!
END SUBROUTINE helmholtz_dimless
