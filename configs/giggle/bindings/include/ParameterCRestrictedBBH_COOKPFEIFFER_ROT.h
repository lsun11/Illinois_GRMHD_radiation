#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_INT genID_cmdline_output_enable;
} RESTRICTED_BBH_COOKPFEIFFER_ROT_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_BBH_COOKPFEIFFER_ROT_STRUCT_PARAMS \
  CCTK_INT const genID_cmdline_output_enable = RESTRICTED_BBH_COOKPFEIFFER_ROT_STRUCT.genID_cmdline_output_enable; \
  enum { \
      dummy_RESTRICTED_BBH_COOKPFEIFFER_ROT_STRUCT_genID_cmdline_output_enable = sizeof( genID_cmdline_output_enable ) \
  }; \

