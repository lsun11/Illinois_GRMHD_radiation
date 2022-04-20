#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL eta;
  CCTK_REAL eta_falloff_dr;
  CCTK_REAL eta_falloff_radius;
  CCTK_REAL eta_final_value;
  CCTK_REAL hyper_Rscale;
  CCTK_REAL hyper_b1;
  CCTK_REAL hyper_b2;
  CCTK_REAL hyper_b3;
  CCTK_REAL hyper_phi_init;
  CCTK_REAL hyper_psi6_init;
  CCTK_REAL time_to_switch_to_eta_final_value;
  const char * shift_interp_operator;
  CCTK_INT eta_falloff_enable;
  CCTK_INT firstorder_shift_convert_Gammai_fisheye_to_physical;
  CCTK_INT hbpunc_advect_enable;
  CCTK_INT hbpuncture_shift_convert_Gammai_fisheye_to_physical;
  CCTK_INT hyper_Integration_Order;
  CCTK_INT punc_shift_bc_radial_falloff_power;
  CCTK_INT shift_bc_type;
} PRIVATE_SHIFT_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_SHIFT_STRUCT_PARAMS \
  CCTK_REAL const eta = PRIVATE_SHIFT_STRUCT.eta; \
  CCTK_REAL const eta_falloff_dr = PRIVATE_SHIFT_STRUCT.eta_falloff_dr; \
  CCTK_REAL const eta_falloff_radius = PRIVATE_SHIFT_STRUCT.eta_falloff_radius; \
  CCTK_REAL const eta_final_value = PRIVATE_SHIFT_STRUCT.eta_final_value; \
  CCTK_REAL const hyper_Rscale = PRIVATE_SHIFT_STRUCT.hyper_Rscale; \
  CCTK_REAL const hyper_b1 = PRIVATE_SHIFT_STRUCT.hyper_b1; \
  CCTK_REAL const hyper_b2 = PRIVATE_SHIFT_STRUCT.hyper_b2; \
  CCTK_REAL const hyper_b3 = PRIVATE_SHIFT_STRUCT.hyper_b3; \
  CCTK_REAL const hyper_phi_init = PRIVATE_SHIFT_STRUCT.hyper_phi_init; \
  CCTK_REAL const hyper_psi6_init = PRIVATE_SHIFT_STRUCT.hyper_psi6_init; \
  CCTK_REAL const time_to_switch_to_eta_final_value = PRIVATE_SHIFT_STRUCT.time_to_switch_to_eta_final_value; \
  const char * const shift_interp_operator = PRIVATE_SHIFT_STRUCT.shift_interp_operator; \
  CCTK_INT const eta_falloff_enable = PRIVATE_SHIFT_STRUCT.eta_falloff_enable; \
  CCTK_INT const firstorder_shift_convert_Gammai_fisheye_to_physical = PRIVATE_SHIFT_STRUCT.firstorder_shift_convert_Gammai_fisheye_to_physical; \
  CCTK_INT const hbpunc_advect_enable = PRIVATE_SHIFT_STRUCT.hbpunc_advect_enable; \
  CCTK_INT const hbpuncture_shift_convert_Gammai_fisheye_to_physical = PRIVATE_SHIFT_STRUCT.hbpuncture_shift_convert_Gammai_fisheye_to_physical; \
  CCTK_INT const hyper_Integration_Order = PRIVATE_SHIFT_STRUCT.hyper_Integration_Order; \
  CCTK_INT const punc_shift_bc_radial_falloff_power = PRIVATE_SHIFT_STRUCT.punc_shift_bc_radial_falloff_power; \
  CCTK_INT const shift_bc_type = PRIVATE_SHIFT_STRUCT.shift_bc_type; \
  enum { \
      dummy_PRIVATE_SHIFT_STRUCT_eta = sizeof( eta ) \
    , dummy_PRIVATE_SHIFT_STRUCT_eta_falloff_dr = sizeof( eta_falloff_dr ) \
    , dummy_PRIVATE_SHIFT_STRUCT_eta_falloff_radius = sizeof( eta_falloff_radius ) \
    , dummy_PRIVATE_SHIFT_STRUCT_eta_final_value = sizeof( eta_final_value ) \
    , dummy_PRIVATE_SHIFT_STRUCT_hyper_Rscale = sizeof( hyper_Rscale ) \
    , dummy_PRIVATE_SHIFT_STRUCT_hyper_b1 = sizeof( hyper_b1 ) \
    , dummy_PRIVATE_SHIFT_STRUCT_hyper_b2 = sizeof( hyper_b2 ) \
    , dummy_PRIVATE_SHIFT_STRUCT_hyper_b3 = sizeof( hyper_b3 ) \
    , dummy_PRIVATE_SHIFT_STRUCT_hyper_phi_init = sizeof( hyper_phi_init ) \
    , dummy_PRIVATE_SHIFT_STRUCT_hyper_psi6_init = sizeof( hyper_psi6_init ) \
    , dummy_PRIVATE_SHIFT_STRUCT_time_to_switch_to_eta_final_value = sizeof( time_to_switch_to_eta_final_value ) \
    , dummy_PRIVATE_SHIFT_STRUCT_shift_interp_operator = sizeof( shift_interp_operator ) \
    , dummy_PRIVATE_SHIFT_STRUCT_eta_falloff_enable = sizeof( eta_falloff_enable ) \
    , dummy_PRIVATE_SHIFT_STRUCT_firstorder_shift_convert_Gammai_fisheye_to_physical = sizeof( firstorder_shift_convert_Gammai_fisheye_to_physical ) \
    , dummy_PRIVATE_SHIFT_STRUCT_hbpunc_advect_enable = sizeof( hbpunc_advect_enable ) \
    , dummy_PRIVATE_SHIFT_STRUCT_hbpuncture_shift_convert_Gammai_fisheye_to_physical = sizeof( hbpuncture_shift_convert_Gammai_fisheye_to_physical ) \
    , dummy_PRIVATE_SHIFT_STRUCT_hyper_Integration_Order = sizeof( hyper_Integration_Order ) \
    , dummy_PRIVATE_SHIFT_STRUCT_punc_shift_bc_radial_falloff_power = sizeof( punc_shift_bc_radial_falloff_power ) \
    , dummy_PRIVATE_SHIFT_STRUCT_shift_bc_type = sizeof( shift_bc_type ) \
  }; \

