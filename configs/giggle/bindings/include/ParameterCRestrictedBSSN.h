#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL Xcenter;
  CCTK_REAL Ycenter;
  CCTK_REAL Zcenter;
  CCTK_REAL ah_excision_fac;
  CCTK_REAL binary_orb_freq;
  CCTK_REAL radmax_rot_interp_1;
  CCTK_REAL radmax_rot_interp_2;
  CCTK_REAL radmax_rot_interp_3;
  CCTK_REAL radmin_rot_interp_1;
  CCTK_REAL radmin_rot_interp_2;
  CCTK_REAL radmin_rot_interp_3;
  CCTK_REAL rah1_initial;
  CCTK_REAL rah2_initial;
  CCTK_REAL rotation_phase;
  CCTK_REAL xbh1_initial;
  CCTK_REAL xbh2_initial;
  CCTK_REAL y_offset;
  CCTK_INT BSSN_MaxNumConstrainedVars;
  CCTK_INT BSSN_MaxNumEvolvedVars;
  CCTK_INT Nlograd_rot_interp;
  CCTK_INT Nphi_rot_interp;
  CCTK_INT Ntheta_rot_interp;
  CCTK_INT bssn_enable_shift_upwind;
  CCTK_INT enable_lower_order_at_boundaries;
  CCTK_INT use_trans_fish_phys_new;
  CCTK_INT zero_out_matter_source_terms;
} RESTRICTED_BSSN_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_BSSN_STRUCT_PARAMS \
  CCTK_REAL const Xcenter = RESTRICTED_BSSN_STRUCT.Xcenter; \
  CCTK_REAL const Ycenter = RESTRICTED_BSSN_STRUCT.Ycenter; \
  CCTK_REAL const Zcenter = RESTRICTED_BSSN_STRUCT.Zcenter; \
  CCTK_REAL const ah_excision_fac = RESTRICTED_BSSN_STRUCT.ah_excision_fac; \
  CCTK_REAL const binary_orb_freq = RESTRICTED_BSSN_STRUCT.binary_orb_freq; \
  CCTK_REAL const radmax_rot_interp_1 = RESTRICTED_BSSN_STRUCT.radmax_rot_interp_1; \
  CCTK_REAL const radmax_rot_interp_2 = RESTRICTED_BSSN_STRUCT.radmax_rot_interp_2; \
  CCTK_REAL const radmax_rot_interp_3 = RESTRICTED_BSSN_STRUCT.radmax_rot_interp_3; \
  CCTK_REAL const radmin_rot_interp_1 = RESTRICTED_BSSN_STRUCT.radmin_rot_interp_1; \
  CCTK_REAL const radmin_rot_interp_2 = RESTRICTED_BSSN_STRUCT.radmin_rot_interp_2; \
  CCTK_REAL const radmin_rot_interp_3 = RESTRICTED_BSSN_STRUCT.radmin_rot_interp_3; \
  CCTK_REAL const rah1_initial = RESTRICTED_BSSN_STRUCT.rah1_initial; \
  CCTK_REAL const rah2_initial = RESTRICTED_BSSN_STRUCT.rah2_initial; \
  CCTK_REAL const rotation_phase = RESTRICTED_BSSN_STRUCT.rotation_phase; \
  CCTK_REAL const xbh1_initial = RESTRICTED_BSSN_STRUCT.xbh1_initial; \
  CCTK_REAL const xbh2_initial = RESTRICTED_BSSN_STRUCT.xbh2_initial; \
  CCTK_REAL const y_offset = RESTRICTED_BSSN_STRUCT.y_offset; \
  CCTK_INT const BSSN_MaxNumConstrainedVars = RESTRICTED_BSSN_STRUCT.BSSN_MaxNumConstrainedVars; \
  CCTK_INT const BSSN_MaxNumEvolvedVars = RESTRICTED_BSSN_STRUCT.BSSN_MaxNumEvolvedVars; \
  CCTK_INT const Nlograd_rot_interp = RESTRICTED_BSSN_STRUCT.Nlograd_rot_interp; \
  CCTK_INT const Nphi_rot_interp = RESTRICTED_BSSN_STRUCT.Nphi_rot_interp; \
  CCTK_INT const Ntheta_rot_interp = RESTRICTED_BSSN_STRUCT.Ntheta_rot_interp; \
  CCTK_INT const bssn_enable_shift_upwind = RESTRICTED_BSSN_STRUCT.bssn_enable_shift_upwind; \
  CCTK_INT const enable_lower_order_at_boundaries = RESTRICTED_BSSN_STRUCT.enable_lower_order_at_boundaries; \
  CCTK_INT const use_trans_fish_phys_new = RESTRICTED_BSSN_STRUCT.use_trans_fish_phys_new; \
  CCTK_INT const zero_out_matter_source_terms = RESTRICTED_BSSN_STRUCT.zero_out_matter_source_terms; \
  enum { \
      dummy_RESTRICTED_BSSN_STRUCT_Xcenter = sizeof( Xcenter ) \
    , dummy_RESTRICTED_BSSN_STRUCT_Ycenter = sizeof( Ycenter ) \
    , dummy_RESTRICTED_BSSN_STRUCT_Zcenter = sizeof( Zcenter ) \
    , dummy_RESTRICTED_BSSN_STRUCT_ah_excision_fac = sizeof( ah_excision_fac ) \
    , dummy_RESTRICTED_BSSN_STRUCT_binary_orb_freq = sizeof( binary_orb_freq ) \
    , dummy_RESTRICTED_BSSN_STRUCT_radmax_rot_interp_1 = sizeof( radmax_rot_interp_1 ) \
    , dummy_RESTRICTED_BSSN_STRUCT_radmax_rot_interp_2 = sizeof( radmax_rot_interp_2 ) \
    , dummy_RESTRICTED_BSSN_STRUCT_radmax_rot_interp_3 = sizeof( radmax_rot_interp_3 ) \
    , dummy_RESTRICTED_BSSN_STRUCT_radmin_rot_interp_1 = sizeof( radmin_rot_interp_1 ) \
    , dummy_RESTRICTED_BSSN_STRUCT_radmin_rot_interp_2 = sizeof( radmin_rot_interp_2 ) \
    , dummy_RESTRICTED_BSSN_STRUCT_radmin_rot_interp_3 = sizeof( radmin_rot_interp_3 ) \
    , dummy_RESTRICTED_BSSN_STRUCT_rah1_initial = sizeof( rah1_initial ) \
    , dummy_RESTRICTED_BSSN_STRUCT_rah2_initial = sizeof( rah2_initial ) \
    , dummy_RESTRICTED_BSSN_STRUCT_rotation_phase = sizeof( rotation_phase ) \
    , dummy_RESTRICTED_BSSN_STRUCT_xbh1_initial = sizeof( xbh1_initial ) \
    , dummy_RESTRICTED_BSSN_STRUCT_xbh2_initial = sizeof( xbh2_initial ) \
    , dummy_RESTRICTED_BSSN_STRUCT_y_offset = sizeof( y_offset ) \
    , dummy_RESTRICTED_BSSN_STRUCT_BSSN_MaxNumConstrainedVars = sizeof( BSSN_MaxNumConstrainedVars ) \
    , dummy_RESTRICTED_BSSN_STRUCT_BSSN_MaxNumEvolvedVars = sizeof( BSSN_MaxNumEvolvedVars ) \
    , dummy_RESTRICTED_BSSN_STRUCT_Nlograd_rot_interp = sizeof( Nlograd_rot_interp ) \
    , dummy_RESTRICTED_BSSN_STRUCT_Nphi_rot_interp = sizeof( Nphi_rot_interp ) \
    , dummy_RESTRICTED_BSSN_STRUCT_Ntheta_rot_interp = sizeof( Ntheta_rot_interp ) \
    , dummy_RESTRICTED_BSSN_STRUCT_bssn_enable_shift_upwind = sizeof( bssn_enable_shift_upwind ) \
    , dummy_RESTRICTED_BSSN_STRUCT_enable_lower_order_at_boundaries = sizeof( enable_lower_order_at_boundaries ) \
    , dummy_RESTRICTED_BSSN_STRUCT_use_trans_fish_phys_new = sizeof( use_trans_fish_phys_new ) \
    , dummy_RESTRICTED_BSSN_STRUCT_zero_out_matter_source_terms = sizeof( zero_out_matter_source_terms ) \
  }; \

