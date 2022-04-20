#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_INT WaveMoL_MaxNumConstrainedVars;
  CCTK_INT WaveMoL_MaxNumEvolvedVars;
} RESTRICTED_SCALARWAVEMOL_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_SCALARWAVEMOL_STRUCT_PARAMS \
  CCTK_INT const WaveMoL_MaxNumConstrainedVars = RESTRICTED_SCALARWAVEMOL_STRUCT.WaveMoL_MaxNumConstrainedVars; \
  CCTK_INT const WaveMoL_MaxNumEvolvedVars = RESTRICTED_SCALARWAVEMOL_STRUCT.WaveMoL_MaxNumEvolvedVars; \
  enum { \
      dummy_RESTRICTED_SCALARWAVEMOL_STRUCT_WaveMoL_MaxNumConstrainedVars = sizeof( WaveMoL_MaxNumConstrainedVars ) \
    , dummy_RESTRICTED_SCALARWAVEMOL_STRUCT_WaveMoL_MaxNumEvolvedVars = sizeof( WaveMoL_MaxNumEvolvedVars ) \
  }; \

