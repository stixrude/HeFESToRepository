!--------------------------------------------------------------------------------------------------
! MODULE   def_variables 
! @brief   Module for the definitions of all variables of tabulation 
! @authors Yu Fang, Marco De Lorenzo 
! @date    20-03-2017
!----------------------------------------------------------------------------------------------------
MODULE def_variables

  USE def_constants, ONLY: NNN_LL, MMM_LL, NNN_sat_LL,&
&                          NNN_LH, MMM_LH, NNN_sat_LH,&
&                          NNN_R , MMM_R , NNN_sat_R ,&
&                          NNN_HT, MMM_HT, NNN_sat_HT,&
&                          NNN_TPL,MMM_TPL,NNN_isop_TPL,&
&                          NNN_TPM,MMM_TPM,NNN_sat_TPM,&
&                          NNN_TPH,MMM_TPH,NNN_sat_TPH,&
&                          ord_spline, pr, NNN_sat_TP,NNN_TP,&
&                          NNN_sat_TP2,NNN_TP2 
  IMPLICIT NONE
!----------------------------------------------------------------------------------------
!
!Variables for the construction of the mesh
!
!----------------------------------------------------------------------------------------
!
! Derivatives for bicubic 
!
!REAL(pr), DIMENSION (NNN_LL,MMM_LL) :: dpdv_LL, dpdu_LL,   &
!&                  d2p_dv2_LL, d2p_du2_LL, d2p_dudv_LL
!REAL(pr), DIMENSION (NNN_LH,MMM_LH) :: dpdv_LH, dpdu_LH,   &
!&                  d2p_dv2_LH, d2p_du2_LH, d2p_dudv_LH
!REAL(pr), DIMENSION (NNN_R,MMM_R)   :: dpdv_R,  dpdu_R,    &
!&                   d2p_dv2_R,  d2p_du2_R, d2p_dudv_R
!REAL(pr), DIMENSION (NNN_HT,MMM_HT) :: dpdv_HT, dpdu_HT,   &
!&                  d2p_dv2_HT, d2p_du2_HT, d2p_dudv_HT
!
! Arrays for the definition of the grid (transformed space)
!
REAL(pr), DIMENSION (MMM_LL)  :: x_mesh_LL
REAL(pr), DIMENSION (MMM_LH)  :: x_mesh_LH
REAL(pr), DIMENSION (MMM_R)   :: x_mesh_R
REAL(pr), DIMENSION (MMM_HT)  :: x_mesh_HT
!
REAL(pr), DIMENSION (NNN_LL)  :: y_mesh_LL
REAL(pr), DIMENSION (NNN_LH)  :: y_mesh_LH
REAL(pr), DIMENSION (NNN_R)   :: y_mesh_R
REAL(pr), DIMENSION (NNN_HT)  :: y_mesh_HT
!
! Arrays for the grid in the Two-Phase Region
!
REAL(pr), DIMENSION (MMM_TPL)  :: x_mesh_TPL
REAL(pr), DIMENSION (MMM_TPM)  :: x_mesh_TPM
REAL(pr), DIMENSION (MMM_TPH)  :: x_mesh_TPH
!
REAL(pr), DIMENSION (NNN_TPL)  :: y_mesh_TPL
REAL(pr), DIMENSION (NNN_TPM)  :: y_mesh_TPM
REAL(pr), DIMENSION (NNN_TPH)  :: y_mesh_TPH
!
! Arrays for the left and right boundaries of LL, LH, R , HT, TPL, TPM, TPH
!
REAL(pr), DIMENSION (NNN_sat_LL)  :: y_mesh_sat_LL, v_Lsat_LL, v_Lpmax_LL,v_liq_meta
REAL(pr), DIMENSION (NNN_sat_LH)  :: y_mesh_sat_LH, v_Lsat_LH, v_Lpmax_LH
REAL(pr), DIMENSION (NNN_sat_R)   :: y_mesh_sat_R,  v_Vsat,    v_Vpmin
REAL(pr), DIMENSION (NNN_sat_HT)  :: y_mesh_sat_HT, v_left_HT, v_right_HT
REAL(pr), DIMENSION (NNN_isop_TPL):: y_mesh_Ptri_TPL, v_right_TPL
REAL(pr), DIMENSION (NNN_sat_TPM) :: y_mesh_sat_TPM, v_right_TPM,v_left_TPM
REAL(pr), DIMENSION (NNN_sat_TPH) :: y_mesh_sat_TPH, v_left_TPH
!
REAL(pr), DIMENSION (ord_spline+1,NNN_LL-1)  :: spline_Lsat_LL, spline_pmax_LL,spline_meta_LL
REAL(pr), DIMENSION (ord_spline+1,NNN_LH-1)  :: spline_Lsat_LH, spline_pmax_LH
REAL(pr), DIMENSION (ord_spline+1,NNN_R -1)  :: spline_pmin, spline_Vsat
REAL(pr), DIMENSION (ord_spline+1,NNN_HT-1)  :: spline_left_HT, spline_right_HT
REAL(pr), DIMENSION (ord_spline+1,NNN_TPL-1) :: spline_right_TPL
REAL(pr), DIMENSION (ord_spline+1,NNN_TPM-1) :: spline_right_TPM,spline_left_TPM
REAL(pr), DIMENSION (ord_spline+1,NNN_TPH-1) :: spline_left_TPH
!
! Values and spline coefficients for saturation curve and TP purposes
!
REAL(pr), DIMENSION (NNN_TP) :: saturP
REAL(pr), DIMENSION (NNN_sat_TP) :: saturP_sat
REAL(pr), DIMENSION (ord_spline+1,NNN_TP-1) :: vL_psat_spline, vV_psat_spline,&
&                                             uL_psat_spline, uV_psat_spline, Tsat_psat_spline
! 5 order spline for saturation curve
REAL(pr), DIMENSION (NNN_TP2) :: saturP2
REAL(pr), DIMENSION (NNN_sat_TP2) :: saturP_sat2
REAL(pr), DIMENSION (6,NNN_TP2-1) :: vL_psat_spline2, vV_psat_spline2,&
&                                             uL_psat_spline2, uV_psat_spline2, Tsat_psat_spline2
! Arrays of properties at the nodes
!
REAL(pr), DIMENSION (NNN_LL,MMM_LL) :: ppp_LL, TTT_LL, vvv_LL, ccc_LL,cpcp_LL
REAL(pr), DIMENSION (NNN_LH,MMM_LH) :: ppp_LH, TTT_LH, vvv_LH, ccc_LH,cpcp_LH
REAL(pr), DIMENSION (NNN_R,MMM_R)   :: ppp_R,  TTT_R,  vvv_R,  ccc_R ,cpcp_R
REAL(pr), DIMENSION (NNN_HT,MMM_HT) :: ppp_HT, TTT_HT, vvv_HT, ccc_HT,cpcp_HT
!
! Arrays of properties at the nodes in the Two-Phase region
!
REAL(pr), DIMENSION (NNN_TPL,MMM_TPL) :: ppp_TPL, vvv_TPL, xxx_TPL, ccc_TPL, TTT_TPL, cpcp_TPL
REAL(pr), DIMENSION (NNN_TPM,MMM_TPM) :: ppp_TPM, vvv_TPM, xxx_TPM, ccc_TPM, TTT_TPM, cpcp_TPM
REAL(pr), DIMENSION (NNN_TPH,MMM_TPH) :: ppp_TPH, vvv_TPH, xxx_TPH, ccc_TPH, TTT_TPH, cpcp_TPH
!
! Bicubic interpolation coefficients 
!
!REAL(pr), DIMENSION (NNN_LL-1,MMM_LL-2,16) :: coeffs_bicub_LL
!REAL(pr), DIMENSION (NNN_LH-1,MMM_LH-2,16) :: coeffs_bicub_LH
!REAL(pr), DIMENSION (NNN_R-1, MMM_R -2,16) :: coeffs_bicub_R
!REAL(pr), DIMENSION (NNN_HT-1,MMM_HT-2,16) :: coeffs_bicub_HT
!
REAL(pr) :: v_corner_A, T_corner_A, v_corner_B, T_corner_B, &
&          v_corner_C, T_corner_C, u_corner_C!,             &
!&          v_corner_D, T_corner_D, u_corner_D,             &
!&          v_corner_E, T_corner_E
!
!
END MODULE def_variables
