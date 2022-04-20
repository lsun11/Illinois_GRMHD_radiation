#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_INT Spatial_Gauge;
  CCTK_INT shift_MaxNumConstrainedVars;
  CCTK_INT shift_MaxNumEvolvedVars;
} RESTRICTED_SHIFT_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_SHIFT_STRUCT_PARAMS \
  CCTK_INT const Spatial_Gauge = RESTRICTED_SHIFT_STRUCT.Spatial_Gauge; \
  CCTK_INT const shift_MaxNumConstrainedVars = RESTRICTED_SHIFT_STRUCT.shift_MaxNumConstrainedVars; \
  CCTK_INT const shift_MaxNumEvolvedVars = RESTRICTED_SHIFT_STRUCT.shift_MaxNumEvolvedVars; \
  enum { \
      dummy_RESTRICTED_SHIFT_STRUCT_Spatial_Gauge = sizeof( Spatial_Gauge ) \
    , dummy_RESTRICTED_SHIFT_STRUCT_shift_MaxNumConstrainedVars = sizeof( shift_MaxNumConstrainedVars ) \
    , dummy_RESTRICTED_SHIFT_STRUCT_shift_MaxNumEvolvedVars = sizeof( shift_MaxNumEvolvedVars ) \
  }; \

