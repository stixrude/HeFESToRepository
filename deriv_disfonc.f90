!---------------------------------------------------------------------------------------------------
! All the formulations  have been taken from the original
! article,the Span-Wagner EoS for CO2.
! The EoS has been published in J. Phys. Chem. Ref. Data, Vol.25,
! pp. 1509-1596, 1996.
!---------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
! SUBROUTINE   deriv_disfonc
! @brief   Compute distance fuction p.1545 Table 32
! @comments All distance fuctions and their derivatives are involved only in
!           termes 40-42 of the residual part
! @authors Yu Fang
! @date    26-10-2016
!----------------------------------------------------------------------------------------------------
SUBROUTINE   deriv_disfonc ( T,v,trib40, trib41, trib42,&
                             dtrib_ddelt40, dtrib_ddelt41, dtrib_ddelt42,&
                             dtrib_dtau40, dtrib_dtau41, dtrib_dtau42,&
                             ddtrib_dddelt40, ddtrib_dddelt41, ddtrib_dddelt42,&
                             ddtrib_ddtau40,ddtrib_ddtau41, ddtrib_ddtau42,&
                             ddtrib_ddeltdtau40, ddtrib_ddeltdtau41,ddtrib_ddeltdtau42 )
     
     USE def_constants

     IMPLICIT NONE
   
!IN/OUT
     REAL(pr) :: T,v

     REAL(pr) :: trib40, trib41, trib42

     REAL(pr) :: dtrib_ddelt40, dtrib_ddelt41, dtrib_ddelt42,&
                 dtrib_dtau40, dtrib_dtau41, dtrib_dtau42

     REAL(pr) :: ddtrib_dddelt40, ddtrib_dddelt41, ddtrib_dddelt42,&
                 ddtrib_ddtau40,ddtrib_ddtau41, ddtrib_ddtau42,&
                 ddtrib_ddeltdtau40,ddtrib_ddeltdtau41,ddtrib_ddeltdtau42

!LOCAL
    REAL(pr) :: tri40,tri41,tri42,th

    REAL(pr) :: tau, delt, rho
    
    REAL(pr) :: dtri_ddelt40, dtri_ddelt41, dtri_ddelt42,&
                 ddtri_dddelt40, ddtri_dddelt41, ddtri_dddelt42

    
    REAL(pr) :: deltm1,invdeltm1,deltm1_s,invbet,inv2betm1,inv2betm2


!Preparation of  variables
 
    rho        = 1_pr    / v
    delt       = 1_pr     / (rho_cr*v)   
    tau        = T_cr    / T
!
!   
    deltm1     = delt    - 1_pr
    invdeltm1  = 1_pr    / deltm1
    deltm1_s  = deltm1 * deltm1
    invbet     = 1_pr    / bet40 
    inv2betm1  = 1_pr    / (2_pr * bet40) - 1_pr
    inv2betm2  = 1_pr    / (2_pr * bet40) - 2_pr
    
    th       = 1_pr - tau + AA40 * deltm1_s**(invbet/2_pr)
    tri40    = th**2_pr + BB40 * deltm1_s**a40
    tri41    = th**2_pr + BB41 * deltm1_s**a41
    tri42    = th**2_pr + BB42 * deltm1_s**a42   
!
!
!------------------------------------------------------------------------------------------------------------------------
!Compute triangle (DELTA)^bi 40-42
!------------------------------------------------------------------------------------------------------------------------
    trib40 = tri40**(b40)
    trib41 = tri41**(b41)
    trib42 = tri42**(b42)
!
!
!--------------------------------------------------------------------------------------------------------------------------
!Compute derivatives of triangle (DELTA) 40-42
!--------------------------------------------------------------------------------------------------------------------------

    dtri_ddelt40 = deltm1 * (AA40 * th * 2_pr / bet40 * deltm1_s**inv2betm1 + 2_pr * BB40 * a40 * deltm1_s**(a40-1_pr))
!    print*, dtri_ddelt40
    dtri_ddelt41 = deltm1 * (AA40 * th * 2_pr / bet40 * deltm1_s**inv2betm1 + 2_pr * BB41 * a41 * deltm1_s**(a41-1_pr))
 
    dtri_ddelt42 = deltm1 * (AA40 * th * 2_pr / bet40 * deltm1_s**inv2betm1 + 2_pr * BB42 * a42 * deltm1_s**(a42-1_pr))

!
!
    ddtri_dddelt40 =  invdeltm1 * dtri_ddelt40 + deltm1_s * (4_pr * BB40 * a40 * (a40 - 1_pr) * deltm1_s**(a40-2_pr) + &
                                                              2_pr * AA40**2_pr * invbet**2_pr * (deltm1_s**inv2betm1)**2_pr +&
                                                              AA40 * th * 4_pr * invbet * inv2betm1 * deltm1_s**inv2betm2) 

    ddtri_dddelt41 =  invdeltm1 * dtri_ddelt41 + deltm1_s * (4_pr * BB41 * a41 * (a41 - 1_pr) * deltm1_s**(a41-2_pr) + &
                                                              2_pr * AA40**2_pr * invbet**2_pr * (deltm1_s**inv2betm1)**2_pr +&
                                                              AA40 * th * 4_pr* invbet * inv2betm1 * deltm1_s**inv2betm2)

    ddtri_dddelt42 =  invdeltm1 * dtri_ddelt42 + deltm1_s * (4_pr * BB42 * a42 * (a42 - 1_pr) * deltm1_s**(a42-2_pr) + &
                                                              2_pr * AA40**2_pr * invbet**2_pr * (deltm1_s**inv2betm1)**2_pr +&
                                                              AA40 * th *4_pr* invbet * inv2betm1 * deltm1_s**inv2betm2)
!
!
!--------------------------------------------------------------------------------------------------------------------------------
!Compute derivatives of triangle^bi  40-42
!--------------------------------------------------------------------------------------------------------------------------------

    dtrib_ddelt40 = b40 * tri40**(b40-1_pr) * dtri_ddelt40
!    print*, dtrib_ddelt40,b40,tri40,dtri_ddelt40

    dtrib_ddelt41 = b41 * tri41**(b41-1_pr) * dtri_ddelt41
  
    dtrib_ddelt42 = b42 * tri42**(b42-1_pr) * dtri_ddelt42
!
!
    dtrib_dtau40  =  - 2_pr * th * b40 * tri40**(b40-1_pr)

    dtrib_dtau41  =  - 2_pr * th * b41 * tri41**(b41-1_pr)

    dtrib_dtau42  =  - 2_pr * th * b42 * tri42**(b42-1_pr)
!
!
    ddtrib_dddelt40 = b40 * (tri40**(b40-1_pr) * ddtri_dddelt40 + (b40 - 1_pr) * tri40**(b40-2_pr) * dtri_ddelt40**2_pr)
 
    ddtrib_dddelt41 = b41 * (tri41**(b41-1_pr) * ddtri_dddelt41 + (b41 - 1_pr) * tri41**(b41-2_pr) * dtri_ddelt41**2_pr)
 
    ddtrib_dddelt42 = b42 * (tri42**(b42-1_pr) * ddtri_dddelt42 + (b42 - 1_pr) * tri42**(b42-2_pr) * dtri_ddelt42**2_pr)
!
!
    ddtrib_ddtau40 = 2_pr * b40 * tri40**(b40-1_pr) + 4_pr * th**2_pr * b40 * (b40 - 1_pr) * tri40**(b40-2_pr)
  
    ddtrib_ddtau41 = 2_pr * b41 * tri41**(b41-1_pr) + 4_pr * th**2_pr * b41 * (b41 - 1_pr) * tri41**(b41-2_pr)
  
    ddtrib_ddtau42 = 2_pr * b42 * tri42**(b42-1_pr) + 4_pr * th**2_pr * b42 * (b42 - 1_pr) * tri42**(b42-2_pr)
!
!
    ddtrib_ddeltdtau40 = - AA40 * b40 * 2_pr * invbet * tri40**(b40-1_pr) * deltm1 * deltm1_s**inv2betm1&
                         - 2_pr * th * b40 * (b40 - 1_pr) * tri40**(b40-2_pr) * dtri_ddelt40
  
    ddtrib_ddeltdtau41 = - AA40 * b41 * 2_pr * invbet * tri41**(b41-1_pr) * deltm1 * deltm1_s**inv2betm1&
                         - 2_pr * th * b41 * (b41 - 1_pr) * tri41**(b41-2_pr) * dtri_ddelt41

    ddtrib_ddeltdtau42 = - AA40 * b42 * 2_pr * invbet * tri42**(b42-1_pr) * deltm1 * deltm1_s**inv2betm1&
                         - 2_pr * th * b42 * (b42 - 1_pr) * tri42**(b42-2_pr) * dtri_ddelt42
  
!
!
END SUBROUTINE deriv_disfonc
