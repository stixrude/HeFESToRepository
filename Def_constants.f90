!---------------------------------------------------------------------------------------------------
! All the coefficients have been taken from the original
! article,the Span-Wagner EoS for CO2.
! The EoS has been published in J. Phys. Chem. Ref. Data, Vol.25,
! pp. 1509-1596, 1996.
!---------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
! MODULE   def_constants 
! @brief   Module for all coefficents in EoS
! @authors Marco De Lorenzo, Yu Fang 
! @date    25-10-2016
!----------------------------------------------------------------------------------------------------
MODULE def_constants

! USE mod_prec_defs, ONLY: pr

  IMPLICIT NONE

  INTEGER, PARAMETER :: pr = selected_real_kind(15,307)  !< working precision consistent with AVBP (here double precision)
!  INTEGER, PARAMETER :: pr = selected_real_kind(32)  !< working precision consistent with AVBP (here double precision)

!-----------------------------------------------------------------------------------------------------
!
! Coefficients for phi_0 (ideal gas part)    FROM p. 1540 of the paper, Table 27
!
!----------------------------------------------------------------------------------------------------

! Coefficient a_i_0
  REAL(pr), PARAMETER :: a1_0= 8.37304456_pr 
  REAL(pr), PARAMETER :: a2_0=-3.70454304_pr
  REAL(pr), PARAMETER :: a3_0= 2.50000000_pr
  REAL(pr), PARAMETER :: a4_0= 1.99427042_pr
  REAL(pr), PARAMETER :: a5_0= 0.62105248_pr
  REAL(pr), PARAMETER :: a6_0= 0.41195293_pr
  REAL(pr), PARAMETER :: a7_0= 1.04028922_pr
  REAL(pr), PARAMETER :: a8_0= 0.08327678_pr

! Coefficients thi_0
  REAL(pr), PARAMETER :: th4_0= 3.151630_pr
  REAL(pr), PARAMETER :: th5_0= 6.111900_pr
  REAL(pr), PARAMETER :: th6_0= 6.777080_pr
  REAL(pr), PARAMETER :: th7_0= 11.32384_pr
  REAL(pr), PARAMETER :: th8_0= 27.08792_pr

!-------------------------------------------------------------------------------------------------
! Physical parameters for CO2
!-------------------------------------------------------------------------------------------------
  REAL(pr), PARAMETER :: W      = 0.0440098_pr        ! Molar mass               [kg/mol]
  REAL(pr), PARAMETER :: R      = 188.9241_pr         ! specific gas constant    [J/kg/K]
! Critical Point 
  REAL(pr), PARAMETER :: T_cr   = 304.1282_pr         ! critical temperature     [K]
  REAL(pr), PARAMETER :: rho_cr = 467.6_pr            ! critical density         [kg/m3]
  REAL(pr), PARAMETER :: P_cr   = 7377300.0_pr        ! critical pressure        [Pa] 
  REAL(pr), PARAMETER :: e_cr   = -190311.39561_pr    ! critical internal energy [J/kg] 
  REAL(pr), PARAMETER :: c_cr   = 130_pr              ! singularity in EoS, average of the previous (T,v)
  REAL(pr), PARAMETER :: cp_cr   = 2000e3_pr              ! singularity in EoS, average of the previous (T,v)
! Triple Point 
  REAL(pr), PARAMETER :: T_tri  = 216.592_pr          ! triple point temperature [K]
  REAL(pr), PARAMETER :: P_tri  = 517950_pr       ! triple point pressure    [Pa]
! e = h-pv
!  REAL(pr), PARAMETER :: e_tri_L  = -427167.1751487_pr! liquid triple point internal energy      [J/kg]
  REAL(pr), PARAMETER :: e_tri_L  = -427183.310638_pr
  REAL(pr), PARAMETER :: e_tri_R  = -114003.706_pr    ! vapor triple point internal energy       [J/kg]
  REAL(pr), PARAMETER :: rho_tri_R  = 13.761_pr       ! vapor triple point density               [kg/m3]
  REAL(pr), PARAMETER :: rho_tri_L  = 1178.46_pr       ! liquid triple point density              [kg/m3]
  REAL(pr), PARAMETER :: v_tri_R  = 1_pr / 13.761_pr  ! vapor triple point specific volume      [m3/kg]
  REAL(pr), PARAMETER :: v_tri_L  = 1_pr / 1178.4_pr  ! liquid triple point specific volume      [m3/kg]
  REAL(pr), PARAMETER :: c_tri_L  = 975.85_pr         ! liquid triple point speed of sound       [m/s]
  REAL(pr), PARAMETER :: c_tri_R  = 222.78_pr         ! vapor triple point speed of sound        [m/s]
  REAL(pr), PARAMETER :: h_tri_L  = -426740_pr        ! liquid triple point specific enthalpy    [J/kg]
  REAL(pr), PARAMETER :: h_tri_R  = -76364_pr         ! vapor  point liquid specific enthalpy    [J/kg]
!  
!  The U maximum point along the saturation curve
  REAL(pr), PARAMETER :: T_umax   = 252.0757_pr       ! Umax temperature     [K]
  REAL(pr), PARAMETER :: rho_umax = 49.9239758_pr     ! Umax density         [kg/m3]
  REAL(pr), PARAMETER :: v_umax   = 2.003045597d-2    ! Umax specific volume [m3/kg]
  REAL(pr), PARAMETER :: P_umax   = 1905180.51_pr     ! Umax pressure        [Pa] 
  REAL(pr), PARAMETER :: e_umax   = -107978.7379_pr   ! Umax internal energy [J/kg] 
  REAL(pr), PARAMETER :: c_umax   = 220.697963_pr     ! Umax speed of sound  [m/s]
 
  REAL(pr), PARAMETER :: T_B      = 432_pr
  REAL(pr), PARAMETER :: v_B      = 1.502084e-3_pr
!  REAL(pr), PARAMETER :: T_B      = 336.94_pr         ! values computed by corner_B
!  REAL(pr), PARAMETER :: v_B      = 3.69e-3_pr
  
  REAL(pr), PARAMETER :: T_C      = 225.553_pr        ! values computed by corner_C
  REAL(pr), PARAMETER :: v_C      = 7.6597e-2_pr
!  
  REAL(pr), PARAMETER :: T_D      = 216.59199_pr        ! values computed by corner_D
  REAL(pr), PARAMETER :: v_D      = 5.51702987d-2
!
  REAL(pr), PARAMETER :: T_E      = 278.2886799_pr        ! values computed by corner_E
  REAL(pr), PARAMETER :: v_E      = 8.68527009d-3
!
  REAL(pr), PARAMETER :: e_F      = -111810.3809_pr   ! point F at saturation line before Umax
  REAL(pr), PARAMETER :: v_F      = 1.41627e-2_pr
  REAL(pr), PARAMETER :: P_F      = 1.905180d6
!  
  REAL(pr), PARAMETER :: e_G      = -170d3   ! point F at saturation line before Umax
  REAL(pr), PARAMETER :: v_G      = 2.126d-3
!------------------------------------------------------------------------------------------------
!
! Coefficients FOR phi_r (residual part)    FROM p. 1544 of the paper Table 31
!
! Some of the coefficients reported in Table 31 of p. 1544, have not
! been rewritten here as parameters but directly into the formulas.
!
! Namely, coefficients d_i, t_i and c_i are not needed since we will
! express differently the powers of 'tau' and 'delta'.
!
!------------------------------------------------------------------------------------------------
   REAL(pr), PARAMETER :: n1 =  0.38856823203161_pr
   REAL(pr), PARAMETER :: n2 =  0.29385475942740e+1_pr
   REAL(pr), PARAMETER :: n3 = -0.55867188534934e+1_pr
   REAL(pr), PARAMETER :: n4 = -0.76753199592477e+0_pr
   REAL(pr), PARAMETER :: n5 =  0.31729005580416e+0_pr
   REAL(pr), PARAMETER :: n6 =  0.54803315897767e+0_pr
   REAL(pr), PARAMETER :: n7 =  0.12279411220335e+0_pr
   REAL(pr), PARAMETER :: n8 =  0.21658961543220e+1_pr
   REAL(pr), PARAMETER :: n9 =  0.15841735109724e+1_pr
   REAL(pr), PARAMETER :: n10= -0.23132705405503e+0_pr
   REAL(pr), PARAMETER :: n11=  0.58116916431436e-1_pr
   REAL(pr), PARAMETER :: n12= -0.55369137205382e+0_pr
   REAL(pr), PARAMETER :: n13=  0.48946615909422e+0_pr
   REAL(pr), PARAMETER :: n14= -0.24275739843501e-1_pr
   REAL(pr), PARAMETER :: n15=  0.62494790501678e-1_pr
   REAL(pr), PARAMETER :: n16= -0.12175860225246e+0_pr
   REAL(pr), PARAMETER :: n17= -0.37055685270086e+0_pr
   REAL(pr), PARAMETER :: n18= -0.16775879700426e-1_pr
   REAL(pr), PARAMETER :: n19= -0.11960736637987e+0_pr
   REAL(pr), PARAMETER :: n20= -0.45619362508778e-1_pr
   REAL(pr), PARAMETER :: n21=  0.35612789270346e-1_pr
   REAL(pr), PARAMETER :: n22= -0.74427727132052e-2_pr
   REAL(pr), PARAMETER :: n23= -0.17395704902432e-2_pr
   REAL(pr), PARAMETER :: n24= -0.21810121289527e-1_pr
   REAL(pr), PARAMETER :: n25=  0.24332166559236e-1_pr
   REAL(pr), PARAMETER :: n26= -0.37440133423463e-1_pr
   REAL(pr), PARAMETER :: n27=  0.14338715756878e+0_pr
   REAL(pr), PARAMETER :: n28= -0.13491969083286e+0_pr
   REAL(pr), PARAMETER :: n29= -0.23151225053480e-1_pr
   REAL(pr), PARAMETER :: n30=  0.12363125492901e-1_pr
   REAL(pr), PARAMETER :: n31=  0.21058321972940e-2_pr
   REAL(pr), PARAMETER :: n32= -0.33958519026368e-3_pr
   REAL(pr), PARAMETER :: n33=  0.55993651771592e-2_pr
   REAL(pr), PARAMETER :: n34= -0.30335118055646e-3_pr
   REAL(pr), PARAMETER :: n35= -0.21365488688320e+3_pr
   REAL(pr), PARAMETER :: n36=  0.26641569149272e+5_pr
   REAL(pr), PARAMETER :: n37= -0.24027212204557e+5_pr
   REAL(pr), PARAMETER :: n38= -0.28341603423999e+3_pr
   REAL(pr), PARAMETER :: n39=  0.21247284400179e+3_pr
   REAL(pr), PARAMETER :: n40= -0.66642276540751e+0_pr
   REAL(pr), PARAMETER :: n41=  0.72608632349897e+0_pr
   REAL(pr), PARAMETER :: n42=  0.55068668612842e-1_pr
!
!
   REAL(pr), PARAMETER :: d1 =  1.0_pr
   REAL(pr), PARAMETER :: d2 =  1.0_pr
   REAL(pr), PARAMETER :: d3 =  1.0_pr
   REAL(pr), PARAMETER :: d4 =  1.0_pr
   REAL(pr), PARAMETER :: d5 =  2.0_pr
   REAL(pr), PARAMETER :: d6 =  2.0_pr
   REAL(pr), PARAMETER :: d7 =  3.0_pr
   REAL(pr), PARAMETER :: d8 =  1.0_pr
   REAL(pr), PARAMETER :: d9 =  2.0_pr
   REAL(pr), PARAMETER :: d10=  4.0_pr
   REAL(pr), PARAMETER :: d11=  5.0_pr
   REAL(pr), PARAMETER :: d12=  5.0_pr
   REAL(pr), PARAMETER :: d13=  5.0_pr
   REAL(pr), PARAMETER :: d14=  6.0_pr
   REAL(pr), PARAMETER :: d15=  6.0_pr
   REAL(pr), PARAMETER :: d16=  6.0_pr
   REAL(pr), PARAMETER :: d17=  1.0_pr
   REAL(pr), PARAMETER :: d18=  1.0_pr
   REAL(pr), PARAMETER :: d19=  4.0_pr
   REAL(pr), PARAMETER :: d20=  4.0_pr
   REAL(pr), PARAMETER :: d21=  4.0_pr
   REAL(pr), PARAMETER :: d22=  7.0_pr
   REAL(pr), PARAMETER :: d23=  8.0_pr
   REAL(pr), PARAMETER :: d24=  2.0_pr
   REAL(pr), PARAMETER :: d25=  3.0_pr
   REAL(pr), PARAMETER :: d26=  3.0_pr
   REAL(pr), PARAMETER :: d27=  5.0_pr
   REAL(pr), PARAMETER :: d28=  5.0_pr
   REAL(pr), PARAMETER :: d29=  6.0_pr
   REAL(pr), PARAMETER :: d30=  7.0_pr
   REAL(pr), PARAMETER :: d31=  8.0_pr
   REAL(pr), PARAMETER :: d32=  10.0_pr
   REAL(pr), PARAMETER :: d33=  4.0_pr
   REAL(pr), PARAMETER :: d34=  8.0_pr
   REAL(pr), PARAMETER :: d35=  2.0_pr
   REAL(pr), PARAMETER :: d36=  2.0_pr
   REAL(pr), PARAMETER :: d37=  2.0_pr
   REAL(pr), PARAMETER :: d38=  3.0_pr
   REAL(pr), PARAMETER :: d39=  3.0_pr
!   
!
!
! Coefficients alpha_i for summation i=35,39
   REAL(pr), PARAMETER :: alp35= 25_pr
   REAL(pr), PARAMETER :: alp36= 25_pr
   REAL(pr), PARAMETER :: alp37= 25_pr
   REAL(pr), PARAMETER :: alp38= 15_pr
   REAL(pr), PARAMETER :: alp39= 20_pr

! Coefficients beta_i for summation i=35,42
   REAL(pr), PARAMETER :: bet35= 325_pr
   REAL(pr), PARAMETER :: bet36= 300_pr
   REAL(pr), PARAMETER :: bet37= 300_pr
   REAL(pr), PARAMETER :: bet38= 275_pr
   REAL(pr), PARAMETER :: bet39= 275_pr
   REAL(pr), PARAMETER :: bet40= 0.3_pr
!   REAL(pr), PARAMETER :: bet41= 0.3d0)
!   REAL(pr), PARAMETER :: bet42= 0.3d0)

! Coefficients gamma_i for summation i=35,39
   REAL(pr), PARAMETER :: gam35= 1.16_pr
   REAL(pr), PARAMETER :: gam36= 1.19_pr
   REAL(pr), PARAMETER :: gam37= 1.19_pr
   REAL(pr), PARAMETER :: gam38= 1.25_pr
   REAL(pr), PARAMETER :: gam39= 1.22_pr

! Coefficients a_i for summation i=40,42
   REAL(pr), PARAMETER :: a40= 3.5_pr
   REAL(pr), PARAMETER :: a41= 3.5_pr
   REAL(pr), PARAMETER :: a42= 3.0_pr

! Coefficients b_i for summation i=40,42
   REAL(pr), PARAMETER :: b40= 0.875_pr
   REAL(pr), PARAMETER :: b41= 0.925_pr
   REAL(pr), PARAMETER :: b42= 0.875_Pr

! Coefficients A_i for summation i=40,42
   REAL(pr), PARAMETER :: AA40= 0.7_pr

! Coefficients B_i for summation i=40,42
   REAL(pr), PARAMETER :: BB40= 0.3_pr
   REAL(pr), PARAMETER :: BB41= 0.3_pr
   REAL(pr), PARAMETER :: BB42= 1.0_pr

! Coefficients C_i for summation i=40,42
   REAL(pr), PARAMETER :: CC40= 10.0_pr
   REAL(pr), PARAMETER :: CC41= 10.0_pr
   REAL(pr), PARAMETER :: CC42= 12.5_pr

! Coefficients D_i for summation i=40,42
   REAL(pr), PARAMETER :: DD40= 275_pr
!
!
! Definition of vertices of the cell in the (x,y) diagram for bicubic
! coefficient construction
!
!REAL(pr) :: x1 = 0d0, y1 = 0d0
!REAL(pr) :: x2 = 0d0, y2 = 1d0
!REAL(pr) :: x3 = 1d0, y3 = 1d0
!--------------------------------------------------------------------
!
! Parameters for the creation of the mesh in physical domains
!
!--------------------------------------------------------------------
! Order of polyfit for spline line
   INTEGER, PARAMETER :: ord_spline = 3
! Left boundary isobar  50MPa
   REAL(pr), PARAMETER :: p_max = 50.0d6
! Residu error
   REAL(pr), PARAMETER :: res_ref = 1.0d-7
! Maximum internal energy smaller than  1000K [J/kg]
   REAL(pr), PARAMETER :: u_end = 540.0e3_pr 
!   REAL(pr), PARAMETER :: e_ref = 500.0e3_pr 
   REAL(pr), PARAMETER :: e_ref = 506.78993e3_pr ! caliberated with NIST
! delta_R = u_min - u_triple in the Right region 
   REAL(pr), PARAMETER, PUBLIC  :: delta_R    = 5e0_pr
! X interval in transformed domain   
   REAL(pr), PARAMETER :: x_mesh_max = 101_pr
   REAL(pr), PARAMETER :: x_mesh_min = 1_pr
!   
   REAL(pr), PARAMETER :: x_mesh2_max = 101_pr
   REAL(pr), PARAMETER :: x_mesh2_min = 1_pr
! MESH GRID PARAMETERS
   INTEGER, PARAMETER :: NNN_LL = 100, MMM_LL = 100, NNN_sat_LL = 3*(NNN_LL-1)+1
   INTEGER, PARAMETER :: NNN_LH = 100, MMM_LH = 100, NNN_sat_LH = 3*(NNN_LH-1)+1
   INTEGER, PARAMETER :: NNN_R  = 100, MMM_R  = 100, NNN_sat_R  = 3*(NNN_R-1) +1
   INTEGER, PARAMETER :: NNN_HT = 100, MMM_HT = 100, NNN_sat_HT = 3*(NNN_HT-1)+1
!
! MESH GRID PARAMETERS TWO-PHASE REGION
   INTEGER, PARAMETER :: NNN_TPL = 100, MMM_TPL = 100, NNN_isop_TPL = 3*(NNN_TPL-1)+1
   INTEGER, PARAMETER :: NNN_TPM = 100, MMM_TPM = 100, NNN_sat_TPM = 3*(NNN_TPM-1)+1
   INTEGER, PARAMETER :: NNN_TPH = 100, MMM_TPH = 100, NNN_sat_TPH = 3*(NNN_TPH-1)+1
!
! NUMBER OF POINTS ON THE SATURATION CURVE FOR THE TWO-PHASE ANALYSIS
   INTEGER, PARAMETER :: NNN_TP = 600, NNN_sat_TP = 3*(NNN_TP-1)+1
   INTEGER, PARAMETER :: NNN_TP2 = 600, NNN_sat_TP2 = 5*(NNN_TP2-1)+1
!
!
!===========================================================================
!
! Parameters for transport properties
!
!============================================================================
!
!---------------------------------------------------------------
!
! Viscosity, project summary report: 
! Thermophysical Properties of CO2 and CO2-Rich Mixture
!
!--------------------------------------------------------------
REAL(pr), PARAMETER :: av_00= 1749.35489318835_pr
REAL(pr), PARAMETER :: av_01=-369.069300007128_pr
REAL(pr), PARAMETER :: av_02= 5423856.34887691_pr
REAL(pr), PARAMETER :: av_03=-2.21283852168356_pr
REAL(pr), PARAMETER :: av_04=-269503.247933569_pr
REAL(pr), PARAMETER :: av_05= 73145.021531826_pr
REAL(pr), PARAMETER :: av_06= 5.34368649509278_pr
!
REAL(pr), PARAMETER :: bv_00=-19.572881_pr
REAL(pr), PARAMETER :: bv_01= 219.73999_pr 
REAL(pr), PARAMETER :: bv_02=-1015.3226_pr
REAL(pr), PARAMETER :: bv_03= 2471.0125_pr
REAL(pr), PARAMETER :: bv_04=-3375.1717_pr
REAL(pr), PARAMETER :: bv_05= 2491.6597_pr
REAL(pr), PARAMETER :: bv_06=-787.26086_pr
REAL(pr), PARAMETER :: bv_07= 14.085455_pr
REAL(pr), PARAMETER :: bv_08=-0.34664158_pr
!
REAL(pr), PARAMETER :: tv_01= 0.25_pr
REAL(pr), PARAMETER :: tv_02= 0.5_pr
REAL(pr), PARAMETER :: tv_03= 0.75_pr
REAL(pr), PARAMETER :: tv_04= 1.0_pr 
REAL(pr), PARAMETER :: tv_05= 1.25_pr
REAL(pr), PARAMETER :: tv_06= 1.5_pr
REAL(pr), PARAMETER :: tv_07= 2.5_pr
REAL(pr), PARAMETER :: tv_08= 5.5_pr
!
REAL(pr), PARAMETER :: gamm_v= 8.06282737481277_pr
REAL(pr), PARAMETER :: cv_01 = 0.360603235428487_pr
REAL(pr), PARAMETER :: cv_02 = 0.121550806591497_pr
!
!---------------------------------------------------------------
!
! Thermal conductivity, project summary report: 
! Thermophysical Properties of CO2 and CO2-Rich Mixture
!
!--------------------------------------------------------------
!
!
REAL(pr), PARAMETER :: lc_00= 0.0151873407_pr
REAL(pr), PARAMETER :: lc_01= 0.0280674040_pr
REAL(pr), PARAMETER :: lc_02= 0.0228564190_pr
REAL(pr), PARAMETER :: lc_03=-0.0074162421_pr
!
REAL(pr), PARAMETER :: b1c_01= 0.0100128_pr
REAL(pr), PARAMETER :: b1c_02= 0.0560488_pr
REAL(pr), PARAMETER :: b1c_03=-0.0811620_pr
REAL(pr), PARAMETER :: b1c_04= 0.0624337_pr
REAL(pr), PARAMETER :: b1c_05=-0.0206336_pr
REAL(pr), PARAMETER :: b1c_06= 0.00253248_pr
!
REAL(pr), PARAMETER :: b2c_01= 0.00430829_pr
REAL(pr), PARAMETER :: b2c_02=-0.0358563_pr
REAL(pr), PARAMETER :: b2c_03= 0.0671480_pr
REAL(pr), PARAMETER :: b2c_04=-0.0522855_pr
REAL(pr), PARAMETER :: b2c_05= 0.0174571_pr
REAL(pr), PARAMETER :: b2c_06=-0.00196414_pr
!

END MODULE def_constants
