#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_INT WaveToyMoL_MaxNumEvolvedVars;
} RESTRICTED_WAVETOYMOL_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_WAVETOYMOL_STRUCT_PARAMS \
  CCTK_INT const WaveToyMoL_MaxNumEvolvedVars = RESTRICTED_WAVETOYMOL_STRUCT.WaveToyMoL_MaxNumEvolvedVars; \
  enum { \
      dummy_RESTRICTED_WAVETOYMOL_STRUCT_WaveToyMoL_MaxNumEvolvedVars = sizeof( WaveToyMoL_MaxNumEvolvedVars ) \
  }; \

