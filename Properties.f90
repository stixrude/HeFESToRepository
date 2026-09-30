!---------------------------------------------------------------------------------------------------
! All the formulations  have been taken from the original
! article,the Span-Wagner EoS for CO2.
! The EoS has been published in J. Phys. Chem. Ref. Data, Vol.25,
! pp. 1509-1596, 1996.
!---------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
! MODULE   properties
! @brief   Compute thermodynamic properties  in IS p.1517 Table3:
!          Pressure, specific internal energy , specific Cv, specific Cp, speed of sound,
!          specific entropy, specific helmholtz energy, specific gibbs energy
! @routine helmholtz_deriv.f90, helmholtz_dimless.f90
! @authors Yu Fang
! @date    10-02-2017
! NOTE: derivatives are added 27/11/2017
!--------------------------------------------------------------------------------------------------

MODULE properties
!
       USE def_constants 
!
       USE def_variables, ONLY: saturP, saturP_sat, vL_psat_spline, vV_psat_spline,&
&                             uL_psat_spline, uV_psat_spline, Tsat_psat_spline,  &
&                             saturP2, saturP_sat2, vL_psat_spline2, vV_psat_spline2,&
&                             uL_psat_spline2, uV_psat_spline2, Tsat_psat_spline2
       IMPLICIT NONE
! 
       PRIVATE
       PUBLIC :: pressure_sw, inter_energy, heat_cap_v, heat_cap_p,dpdv_T, dpdT_v,&
&                sound_speed, entropy, helmho, gibbs,dpdu_v, dpdr_u, satprop, satderiv, axlpress,&
&                dedr_T
       CONTAINS
!
!
!===============================================================================================
       SUBROUTINE pressure_sw(T,v,p)  !Pa
!===============================================================================================
        IMPLICIT NONE
!
!IN/OUT
        REAL(pr) :: T,v
        REAL(pr) :: p              
!LOCAL
        REAL(pr) :: rho,delt
!
!for helmholtz_deriv
!
        REAL(pr) :: phi_r_ddelt, phi_r_dtau,&
                    phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau
!Pre-compute
!
        rho      = 1_pr  / v
        delt     = 1_pr   / (rho_cr*v)
!
        CALL helmholtz_deriv( T, v, phi_r_ddelt, phi_r_dtau,&
                            phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau )
!       print*, phi_r_ddelt, phi_r_dtau,phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau
!
        p = (1_pr + delt * phi_r_ddelt) * R * T * rho
!
       END SUBROUTINE pressure_sw
!
!=================================================================================================
       SUBROUTINE inter_energy(T,v,e) !J/kg
!=================================================================================================
        IMPLICIT NONE
!
!IN/OUT
        REAL(pr) :: T,v,e
!LOCAL
        REAL(pr) :: tau
!for helmholtz_deriv

        REAL(pr) :: phi_r_ddelt, phi_r_dtau,&
                    phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau

!for helmholtz_dimless

        REAL(pr) :: phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                    phi0_dddelt,phi0_ddtau
!
!Pre-compute

        tau      = T_cr  / T
!initialisation
      phi_0 = 0_pr
      phi_r = 0_pr
      phi0_ddelt = 0_pr
      phi0_dtau = 0_pr
      phi0_dddelt = 0_pr
      phi0_ddtau = 0_pr
      phi_r_ddelt= 0_pr
      phi_r_dtau = 0_pr
      phi_r_dddelt = 0_pr
      phi_r_ddtau = 0_pr
      phi_r_ddeltdtau = 0_pr
!
        CALL helmholtz_dimless (T,v,phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                              phi0_dddelt,phi0_ddtau)
!       print*, phi_0,phi_r,phi0_ddelt,phi0_dtau,phi0_dddelt,phi0_ddtau
!
!        print*, T,v
        CALL helmholtz_deriv( T, v, phi_r_ddelt, phi_r_dtau,&
                            phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau )
!       print*, phi_r_ddelt, phi_r_dtau,phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau
!
!
        e = tau * ( phi0_dtau + phi_r_dtau) * R * T
!
!        print*, 'phi_r_dtau(properties) = ',phi_r_dtau
!
       END SUBROUTINE inter_energy
!
!============================================================================================================
       SUBROUTINE heat_cap_v(T,v,cv)  !J/kgK
!============================================================================================================
        IMPLICIT NONE
!
!IN/OUT
        REAL(pr) :: T,v,cv
!LOCAL
        REAL(pr) :: tau, tau2
!for helmholtz_deriv

        REAL(pr) :: phi_r_ddelt, phi_r_dtau,&
                    phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau

!for helmholtz_dimless

        REAL(pr) :: phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                    phi0_dddelt,phi0_ddtau
!
!Pre-compute

        tau      = T_cr  / T
        tau2     = tau   * tau
!initialisation
      phi_0 = 0_pr
      phi_r = 0_pr
      phi0_ddelt = 0_pr
      phi0_dtau = 0_pr
      phi0_dddelt = 0_pr
      phi0_ddtau = 0_pr
      phi_r_ddelt= 0_pr
      phi_r_dtau = 0_pr
      phi_r_dddelt = 0_pr
      phi_r_ddtau = 0_pr
      phi_r_ddeltdtau = 0_pr


        CALL helmholtz_dimless (T,v,phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                              phi0_dddelt,phi0_ddtau)
!       print*, phi_0,phi_r,phi0_ddelt,phi0_dtau,phi0_dddelt,phi0_ddtau

        CALL helmholtz_deriv( T, v, phi_r_ddelt, phi_r_dtau,&
                            phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau )
!       print*, phi_r_ddelt, phi_r_dtau,phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau
!
!
!
        cv = -tau2 * ( phi0_ddtau + phi_r_ddtau ) * R

!       print*, phi0_ddtau , phi_r_ddtau , tau2
!
!
       END SUBROUTINE heat_cap_v
!
!================================================================================================================        
       SUBROUTINE heat_cap_p(T,v,cp)  !J/kgK
!================================================================================================================
        IMPLICIT NONE
!IN/OUT
        REAL(pr) :: T,v,cp
!LOCAL
        REAL(pr) :: rho,delt,tau,tau2,delt2

!for helmholtz_deriv

        REAL(pr) :: phi_r_ddelt, phi_r_dtau,&
                    phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau

!for helmholtz_dimless

        REAL(pr) :: phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                    phi0_dddelt,phi0_ddtau


!Pre-compute

        rho      = 1_pr  / v
        delt     = rho   / rho_cr
        tau      = T_cr  / T
        tau2     = tau   * tau
        delt2    = delt  * delt
!initialisation
      phi_0 = 0_pr
      phi_r = 0_pr
      phi0_ddelt = 0_pr
      phi0_dtau = 0_pr
      phi0_dddelt = 0_pr
      phi0_ddtau = 0_pr
      phi_r_ddelt= 0_pr
      phi_r_dtau = 0_pr
      phi_r_dddelt = 0_pr
      phi_r_ddtau = 0_pr
      phi_r_ddeltdtau = 0_pr

        CALL helmholtz_dimless (T,v,phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                              phi0_dddelt,phi0_ddtau)
!       print*, phi_0,phi_r,phi0_ddelt,phi0_dtau,phi0_dddelt,phi0_ddtau

        CALL helmholtz_deriv( T, v, phi_r_ddelt, phi_r_dtau,&
                            phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau )
!       print*, phi_r_ddelt, phi_r_dtau,phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau
!
!
        cp = ( -tau2   * ( phi0_ddtau + phi_r_ddtau ) + &
             (1_pr+delt*phi_r_ddelt-delt*tau*phi_r_ddeltdtau)**(2_pr) / &
             (1_pr+2_pr*delt*phi_r_ddelt+delt2*phi_r_dddelt) ) * R
!
!
       END SUBROUTINE heat_cap_p
!
!=============================================================================================================
       SUBROUTINE sound_speed(T,v,c)   !m/s
!=============================================================================================================
        IMPLICIT NONE
!IN/OUT
        REAL(pr) :: T,v,c
!LOCAL
        REAL(pr) :: rho,delt,tau,tau2,delt2,c_2

!for helmholtz_deriv

        REAL(pr) :: phi_r_ddelt, phi_r_dtau,&
                    phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau

!for helmholtz_dimless

        REAL(pr) :: phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                    phi0_dddelt,phi0_ddtau


!Pre-compute

        rho      = 1_pr  / v
        delt     = 1_pr  / (rho_cr*v)
        tau      = T_cr  / T
        tau2     = tau   * tau
        delt2    = delt  * delt

!initialisation
      phi_0 = 0_pr
      phi_r = 0_pr
      phi0_ddelt = 0_pr
      phi0_dtau = 0_pr
      phi0_dddelt = 0_pr
      phi0_ddtau = 0_pr
      phi_r_ddelt= 0_pr
      phi_r_dtau = 0_pr
      phi_r_dddelt = 0_pr
      phi_r_ddtau = 0_pr
      phi_r_ddeltdtau = 0_pr

        CALL helmholtz_dimless (T,v,phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                              phi0_dddelt,phi0_ddtau)
!       print*, phi_0,phi_r,phi0_ddelt,phi0_dtau,phi0_dddelt,phi0_ddtau

        CALL helmholtz_deriv( T, v, phi_r_ddelt, phi_r_dtau,&
                            phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau )
!       print*, phi_r_ddelt, phi_r_dtau,phi_r_dddelt, phi_r_ddtau,
!       phi_r_ddeltdtau
!
!
        c_2 = ( 1_pr + 2_pr*delt*phi_r_ddelt + delt2*phi_r_dddelt - &
              (1_pr + delt*phi_r_ddelt  - delt*tau*phi_r_ddeltdtau)**(2_pr) / &
              ( tau2 * (phi0_ddtau+phi_r_ddtau) ) ) * R * T
!        print*, tau2 * (phi0_ddtau+phi_r_ddtau)
!         print*, delt*phi_r_dddelt,delt*tau* phi_r_ddeltdtau
!         print*, (1_pr + delt*phi_r_ddelt  - delt*tau*phi_r_ddeltdtau)
!         print*, delt*phi_r_ddelt
        c = abs( sqrt(c_2) )
!
!
       END SUBROUTINE sound_speed
!
!
!====================================================================================================================
       SUBROUTINE entropy(T,v,s)  !J/kgK
!====================================================================================================================
        IMPLICIT NONE
!IN/OUT
        REAL(pr) :: T,v,s
!LOCAL
        REAL(pr) :: tau
!for helmholtz_deriv

        REAL(pr) :: phi_r_ddelt, phi_r_dtau,&
                    phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau

!for helmholtz_dimless

        REAL(pr) :: phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                    phi0_dddelt,phi0_ddtau


!Pre-compute

        tau      = T_cr  / T
!
!initialisation
      phi_0 = 0_pr
      phi_r = 0_pr
      phi0_ddelt = 0_pr
      phi0_dtau = 0_pr
      phi0_dddelt = 0_pr
      phi0_ddtau = 0_pr
      phi_r_ddelt= 0_pr
      phi_r_dtau = 0_pr
      phi_r_dddelt = 0_pr
      phi_r_ddtau = 0_pr
      phi_r_ddeltdtau = 0_pr
 
        CALL helmholtz_dimless (T,v,phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                              phi0_dddelt,phi0_ddtau)
!       print*, phi_0,phi_r,phi0_ddelt,phi0_dtau,phi0_dddelt,phi0_ddtau

        CALL helmholtz_deriv( T, v, phi_r_ddelt, phi_r_dtau,&
                            phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau )
!       print*, phi_r_ddelt, phi_r_dtau,phi_r_dddelt, phi_r_ddtau,
!       phi_r_ddeltdtau
!
!
        s = (tau * ( phi0_dtau + phi_r_dtau) - phi_0 - phi_r) * R
!
!
       END SUBROUTINE entropy
!
!================================================================================================================
       SUBROUTINE helmho(T,v,h_helmho) !J/kg
!================================================================================================================
        IMPLICIT NONE
!IN/OUT
        REAL(pr) :: T,v,h_helmho
!LOCAL
        REAL(pr) :: phi
!for helmholtz_dimless

        REAL(pr) :: phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                    phi0_dddelt,phi0_ddtau
!
!initialisation
      phi_0 = 0_pr
      phi_r = 0_pr
      phi0_ddelt = 0_pr
      phi0_dtau = 0_pr
      phi0_dddelt = 0_pr
      phi0_ddtau = 0_pr

!
        CALL helmholtz_dimless (T,v,phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                              phi0_dddelt,phi0_ddtau)
!
!
        phi = phi_0 + phi_r
! print*, 'propoerties phi0, phir',phi_0, phi_r
! phi is a dimensionless magnitude, R is in [J/kg/K], T in [K]:
!
        h_helmho   = phi * R * T                      ! [J/kg]
!
!
       END SUBROUTINE helmho
!
!
!=================================================================================================================
       SUBROUTINE gibbs(T,v,g) !J/kg
!=================================================================================================================
        IMPLICIT NONE
!IN/OUT
        REAL(pr) :: T,v,g
!LOCAL
        REAL(pr) :: h_helmho, p
!
!
        CALL helmho(T,v,h_helmho)
        CALL pressure_sw(T,v,p)
!
!
        g = h_helmho + p * v
!
!
       END SUBROUTINE gibbs
!
!
!==================================================================================================================
      SUBROUTINE dpdr_u(T,v,deriv) !Pa.m3/kg
!=================================================================================================================
        IMPLICIT NONE
!IN/OUT
        REAL(pr) :: T,v,deriv
!LOCAL
        REAL(pr) :: rho, delt, phi_r_ddelt, phi_r_dtau,&
&                   phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau
        REAL(pr) :: phi_0,phi_r,phi0_ddelt,phi0_dtau,phi0_dddelt,phi0_ddtau,tau
        REAL(pr) :: dp_dr_T, dp_du_r, du_dr_T
!
        rho      = 1_pr  / v
        delt     = 1_pr  / (rho_cr*v)
        tau      = T_cr/T
        phi_r_ddelt = 0_pr
        phi_r_dtau  = 0_pr
        phi_r_dddelt= 0_pr
        phi_r_ddtau = 0_pr
        phi_r_ddeltdtau = 0_pr


!
        CALL helmholtz_dimless (T,v,phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                              phi0_dddelt,phi0_ddtau)

        CALL helmholtz_deriv( T, v, phi_r_ddelt, phi_r_dtau,&
&                             phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau )
!
        CALL dpdu_v(T,v,dp_du_r)
        dp_dr_T  = (1_pr + 2_pr*delt*phi_r_ddelt + delt*delt*phi_r_dddelt) *R *T
        du_dr_T  = R*T*tau/rho_cr * phi_r_ddeltdtau
!
        deriv = dp_dr_T - dp_du_r*du_dr_T

      END SUBROUTINE dpdr_u
!==================================================================================================================
      SUBROUTINE dpdv_T(T,v,deriv) !Pa.m3/kg
!=================================================================================================================
        IMPLICIT NONE
!IN/OUT
        REAL(pr) :: T,v,deriv
!LOCAL
        REAL(pr) :: rho, delt, phi_r_ddelt, phi_r_dtau,&
&                   phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau
        REAL(pr) :: phi_0,phi_r,phi0_ddelt,phi0_dtau,phi0_dddelt,phi0_ddtau,tau
        REAL(pr) :: dp_dr_T, dp_du_r, du_dr_T
!
        rho      = 1_pr  / v
        delt     = 1_pr  / (rho_cr*v)
        tau      = T_cr/T
        phi_r_ddelt = 0_pr
        phi_r_dtau  = 0_pr
        phi_r_dddelt= 0_pr
        phi_r_ddtau = 0_pr
        phi_r_ddeltdtau = 0_pr


!
        CALL helmholtz_dimless (T,v,phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                              phi0_dddelt,phi0_ddtau)

        CALL helmholtz_deriv( T, v, phi_r_ddelt, phi_r_dtau,&
&                             phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau )
!
        CALL dpdu_v(T,v,dp_du_r)
        dp_dr_T  = (1_pr + 2_pr*delt*phi_r_ddelt + delt*delt*phi_r_dddelt) *R *T
!
        deriv = dp_dr_T *(-rho*rho)

      END SUBROUTINE dpdv_T
!
!==================================================================================================================
      SUBROUTINE dpdT_v(T,v,deriv) 
!=================================================================================================================
        IMPLICIT NONE
!IN/OUT
        REAL(pr) :: T,v,deriv
!LOCAL
        REAL(pr) :: rho, delt, phi_r_ddelt, phi_r_dtau,&
&                   phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau
        REAL(pr) :: tau
!
        REAL(pr) :: dp_dT_v, du_dT_v
        rho      = 1_pr  / v
        delt     = 1_pr  / (rho_cr*v)
        tau      = T_cr/T
        phi_r_ddelt = 0_pr
        phi_r_dtau  = 0_pr
        phi_r_dddelt= 0_pr
        phi_r_ddtau = 0_pr
        phi_r_ddeltdtau = 0_pr
!
! 
        CALL helmholtz_deriv( T, v, phi_r_ddelt, phi_r_dtau,&
&                             phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau )
!        
        dp_dT_v = rho * R * (1_pr+delt*phi_r_ddelt-tau*delt*phi_r_ddeltdtau )


        deriv = dp_dT_v


      END SUBROUTINE dpdT_v
!

!
!==================================================================================================================
      SUBROUTINE dpdu_v(T,v,deriv) !Pa/J.kg/m3
!=================================================================================================================
        IMPLICIT NONE
!IN/OUT
        REAL(pr) :: T,v,deriv
!LOCAL
        REAL(pr) :: rho, delt, phi_r_ddelt, phi_r_dtau,&
&                   phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau
        REAL(pr) :: phi_0,phi_r,phi0_ddelt,phi0_dtau,phi0_dddelt,phi0_ddtau,tau
!
        REAL(pr) :: dp_dT_v, du_dT_v
        rho      = 1_pr  / v
        delt     = 1_pr  / (rho_cr*v)
        tau      = T_cr/T
        phi_r_ddelt = 0_pr
        phi_r_dtau  = 0_pr
        phi_r_dddelt= 0_pr
        phi_r_ddtau = 0_pr
        phi_r_ddeltdtau = 0_pr
!
!
!
        CALL helmholtz_dimless (T,v,phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                              phi0_dddelt,phi0_ddtau)
! 
        CALL helmholtz_deriv( T, v, phi_r_ddelt, phi_r_dtau,&
&                             phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau )
!        
        dp_dT_v = rho * R * (1_pr+delt*phi_r_ddelt-tau*delt*phi_r_ddeltdtau )
        du_dT_v = -R*tau*tau*( phi0_ddtau + phi_r_ddtau ) 


        deriv = dp_dT_v / du_dT_v


      END SUBROUTINE dpdu_v
!
!
!=================================================================================
!     
     SUBROUTINE satprop(mode, psat, Tsat, vvsat, vlsat, uvsat, ulsat)
!
!=================================================================================
      IMPLICIT NONE
!
      INTEGER  :: i, j, j_sat
!
      REAL(pr), INTENT(in)  :: psat
      INTEGER , INTENT(in)  :: mode
!
      REAL(pr), INTENT(out) :: Tsat, vvsat,vlsat,uvsat,ulsat
!
      REAL(pr) :: delta, pp, temp
      REAL(pr) :: vL, uL, vV, uV
!
!
!use 3e order spline
!
      IF (mode == 3) THEN
!
       delta = saturP(2) - saturP(1)
       j_sat = INT((psat  - saturP(1))/delta) + 1
!
!##computing saturation quantities
!
       vL = 0_pr; uL = 0_pr;  vV = 0_pr;  uV = 0_pr; temp = 0_pr
       DO i = 1, ord_spline+1
          pp = psat**(i-1)
          vL   = vL  + vL_psat_spline  (ord_spline+2-i, j_sat) * pp
          uL   = uL  + uL_psat_spline  (ord_spline+2-i, j_sat) * pp
          vV   = vV  + vV_psat_spline  (ord_spline+2-i, j_sat) * pp
          uV   = uV  + uV_psat_spline  (ord_spline+2-i, j_sat) * pp
          temp = temp+ Tsat_psat_spline(ord_spline+2-i, j_sat) * pp
       ENDDO
!
!##for output saturation quantities
!
       Tsat  = temp
       vvsat = vV
       vlsat = vL
       uvsat = uV
       ulsat = uL
!
!use 5e order spline
!
      ELSEIF (mode == 5) THEN
!
       delta = saturP2(2) - saturP2(1)
       j_sat = INT((psat  - saturP2(1))/delta) + 1
!
!##computing saturation quantities
!
       vL = 0_pr; uL = 0_pr;  vV = 0_pr;  uV = 0_pr
       DO i = 1, 6
          pp = psat**(i-1)
          vL   = vL  + vL_psat_spline2  (5+2-i, j_sat) * pp
          uL   = uL  + uL_psat_spline2  (5+2-i, j_sat) * pp
          vV   = vV  + vV_psat_spline2  (5+2-i, j_sat) * pp
          uV   = uV  + uV_psat_spline2  (5+2-i, j_sat) * pp
          temp = temp+ Tsat_psat_spline2(5+2-i, j_sat) * pp
       ENDDO
!
!##for output saturation quantities
!
       Tsat  = temp
       vvsat = vV
       vlsat = vL
       uvsat = uV
       ulsat = uL
!
      ENDIF
     END SUBROUTINE satprop
!
!=================================================================================
!     
     SUBROUTINE satderiv(mode, psat, duL_dp, duV_dp, dvL_dp, dvV_dp)
!
!=================================================================================
      IMPLICIT NONE
!
      INTEGER  :: i, j, j_sat
!
      REAL(pr), INTENT(in)  :: psat
      INTEGER , INTENT(in)  :: mode
!
      REAL(pr), INTENT(out) :: duL_dp, duV_dp, dvL_dp, dvV_dp
!
      REAL(pr) :: delta, pp
      REAL(pr) :: duLdp, duVdp, dvLdp, dvVdp
!
!
!use 3e order spline
!
      IF (mode == 3) THEN
!
       delta = saturP(2) - saturP(1)
       j_sat = INT((psat  - saturP(1))/delta) + 1
!
!##computing derivatives
!
       duLdp = 0_pr;  duVdp = 0_pr;  dvLdp = 0_pr;  dvVdp = 0_pr
       DO i = 1, ord_spline
          pp      = psat**(ord_spline - i)
          duLdp  = duLdp + (ord_spline+1-i) * uL_psat_spline(i,j_sat) *pp
          duVdp  = duVdp + (ord_spline+1-i) * uV_psat_spline(i,j_sat) *pp
          dvLdp  = dvLdp + (ord_spline+1-i) * vL_psat_spline(i,j_sat) *pp
          dvVdp  = dvVdp + (ord_spline+1-i) * vV_psat_spline(i,j_sat) *pp
       ENDDO
!
!##for output saturation quantities
!
       duL_dp = duLdp
       duV_dp = duVdp
       dvL_dp = dvLdp
       dvV_dp = dvVdp
!
!use 5e order spline
!
      ELSEIF (mode == 5) THEN
!
       delta = saturP2(2) - saturP2(1)
       j_sat = INT((psat  - saturP2(1))/delta) + 1
!
!##computing derivatives
!
       duLdp = 0_pr;  duVdp = 0_pr;  dvLdp = 0_pr;  dvVdp = 0_pr
       DO i = 1, 5
          pp      = psat**(5 - i)
          duLdp  = duLdp + (5+1-i) * uL_psat_spline2(i,j_sat) *pp
          duVdp  = duVdp + (5+1-i) * uV_psat_spline2(i,j_sat) *pp
          dvLdp  = dvLdp + (5+1-i) * vL_psat_spline2(i,j_sat) *pp
          dvVdp  = dvVdp + (5+1-i) * vV_psat_spline2(i,j_sat) *pp
       ENDDO
!
!##for output saturation quantities
!
       duL_dp = duLdp
       duV_dp = duVdp
       dvL_dp = dvLdp
       dvV_dp = dvVdp
!
      ENDIF
     END SUBROUTINE satderiv
!==================================================================================================================
      SUBROUTINE axlpress(T,vvap, vliq, p1, p2, p3) !Pa
!=================================================================================================================
        IMPLICIT NONE
!IN/OUT
        REAL(pr), INTENT(IN)  :: vvap,vliq,T
        REAL(pr), INTENT(OUT) :: p1,p2,p3
!LOCAL
        REAL(pr) :: deltp,deltpp,tau 
        REAL(pr) :: phi_r_ddelt1, phi_r_ddelt2, phi_r_dtau,&
&                   phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau
        REAL(pr) :: phi_0,phi_r1,phi_r2,phi0_ddelt,phi0_dtau,phi0_dddelt,phi0_ddtau
!
       
        deltp      = 1_pr  / (rho_cr*vliq)
        deltpp     = 1_pr  / (rho_cr*vvap)
        tau        = T_cr/T
!
        phi_r1 = 0_pr
        phi_r2 = 0_pr
        phi_r_ddelt1 = 0_pr
        phi_r_ddelt2 = 0_pr
        phi_r_dtau  = 0_pr
        phi_r_dddelt= 0_pr
        phi_r_ddtau = 0_pr
        phi_r_ddeltdtau = 0_pr
!
!
!
        CALL helmholtz_dimless (T,vliq,phi_0,phi_r1,phi0_ddelt,phi0_dtau,&
                              phi0_dddelt,phi0_ddtau)
! 
        CALL helmholtz_deriv( T, vliq, phi_r_ddelt1, phi_r_dtau,&
&                             phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau )
!        
        CALL helmholtz_dimless (T,vvap,phi_0,phi_r2,phi0_ddelt,phi0_dtau,&
                              phi0_dddelt,phi0_ddtau)
! 
        CALL helmholtz_deriv( T, vvap, phi_r_ddelt2, phi_r_dtau,&
&                             phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau )
!----------eq2.2a, 2.2b, 2.2C--------------
        p1 = R*T/vliq * (1.0+deltp*phi_r_ddelt1)
        p2 = R*T/vvap * (1.0+deltpp*phi_r_ddelt2)
        p3 = ( phi_r1 - phi_r2 + dlog(vvap/vliq) ) * R*T / (vvap-vliq)
!
!
      END SUBROUTINE axlpress
!
!==================================================================================================================
      SUBROUTINE dedr_T(T,v,deriv)
!=================================================================================================================
        IMPLICIT NONE
!IN/OUT
        REAL(pr) :: T,v,deriv
!LOCAL
        REAL(pr) :: rho, delt, phi_r_ddelt, phi_r_dtau,&
&                   phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau
!
        REAL(pr) :: de_dr_T, tau
        rho      = 1_pr  / v
        delt     = 1_pr  / (rho_cr*v)
        tau      = T_cr/T
        phi_r_ddelt = 0_pr
        phi_r_dtau  = 0_pr
        phi_r_dddelt= 0_pr
        phi_r_ddtau = 0_pr
        phi_r_ddeltdtau = 0_pr
!
! 
        CALL helmholtz_deriv( T, v, phi_r_ddelt, phi_r_dtau,&
&                             phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau )
!        
        de_dr_T = R * T_cr / rho_cr * phi_r_ddeltdtau


        deriv = de_dr_T


      END SUBROUTINE dedr_T
!
!
END MODULE properties
