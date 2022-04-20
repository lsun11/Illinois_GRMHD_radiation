#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL hyper_a1;
  CCTK_REAL hyper_a2;
  CCTK_REAL hyper_a3;
  const char * slicing_type;
  CCTK_INT f_of_alpha;
  CCTK_INT hyper_Integration_Order;
  CCTK_INT lapse_MaxNumConstrainedVars;
  CCTK_INT lapse_MaxNumEvolvedVars;
  CCTK_INT opl_a0lap;
  CCTK_INT opl_advect_enable;
  CCTK_INT opl_alap;
} RESTRICTED_LAPSE_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_LAPSE_STRUCT_PARAMS \
  CCTK_REAL const hyper_a1 = RESTRICTED_LAPSE_STRUCT.hyper_a1; \
  CCTK_REAL const hyper_a2 = RESTRICTED_LAPSE_STRUCT.hyper_a2; \
  CCTK_REAL const hyper_a3 = RESTRICTED_LAPSE_STRUCT.hyper_a3; \
  const char * const slicing_type = RESTRICTED_LAPSE_STRUCT.slicing_type; \
  CCTK_INT const f_of_alpha = RESTRICTED_LAPSE_STRUCT.f_of_alpha; \
  CCTK_INT const hyper_Integration_Order = RESTRICTED_LAPSE_STRUCT.hyper_Integration_Order; \
  CCTK_INT const lapse_MaxNumConstrainedVars = RESTRICTED_LAPSE_STRUCT.lapse_MaxNumConstrainedVars; \
  CCTK_INT const lapse_MaxNumEvolvedVars = RESTRICTED_LAPSE_STRUCT.lapse_MaxNumEvolvedVars; \
  CCTK_INT const opl_a0lap = RESTRICTED_LAPSE_STRUCT.opl_a0lap; \
  CCTK_INT const opl_advect_enable = RESTRICTED_LAPSE_STRUCT.opl_advect_enable; \
  CCTK_INT const opl_alap = RESTRICTED_LAPSE_STRUCT.opl_alap; \
  enum { \
      dummy_RESTRICTED_LAPSE_STRUCT_hyper_a1 = sizeof( hyper_a1 ) \
    , dummy_RESTRICTED_LAPSE_STRUCT_hyper_a2 = sizeof( hyper_a2 ) \
    , dummy_RESTRICTED_LAPSE_STRUCT_hyper_a3 = sizeof( hyper_a3 ) \
    , dummy_RESTRICTED_LAPSE_STRUCT_slicing_type = sizeof( slicing_type ) \
    , dummy_RESTRICTED_LAPSE_STRUCT_f_of_alpha = sizeof( f_of_alpha ) \
    , dummy_RESTRICTED_LAPSE_STRUCT_hyper_Integration_Order = sizeof( hyper_Integration_Order ) \
    , dummy_RESTRICTED_LAPSE_STRUCT_lapse_MaxNumConstrainedVars = sizeof( lapse_MaxNumConstrainedVars ) \
    , dummy_RESTRICTED_LAPSE_STRUCT_lapse_MaxNumEvolvedVars = sizeof( lapse_MaxNumEvolvedVars ) \
    , dummy_RESTRICTED_LAPSE_STRUCT_opl_a0lap = sizeof( opl_a0lap ) \
    , dummy_RESTRICTED_LAPSE_STRUCT_opl_advect_enable = sizeof( opl_advect_enable ) \
    , dummy_RESTRICTED_LAPSE_STRUCT_opl_alap = sizeof( opl_alap ) \
  }; \

