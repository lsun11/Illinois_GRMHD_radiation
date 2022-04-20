#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  const char * cctk_run_title;
  const char * cctk_timer_output;
  const char * info_format;
  const char * recovery_mode;
  CCTK_INT allow_mixeddim_gfs;
  CCTK_INT cctk_brief_output;
  CCTK_INT cctk_full_warnings;
  CCTK_INT cctk_show_banners;
  CCTK_INT cctk_show_schedule;
  CCTK_INT cctk_strong_param_check;
  CCTK_INT highlight_warning_messages;
  CCTK_INT manual_cache_setup;
  CCTK_INT manual_cache_size;
  CCTK_INT manual_cacheline_bytes;
} PRIVATE_CACTUS_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_CACTUS_STRUCT_PARAMS \
  const char * const cctk_run_title = PRIVATE_CACTUS_STRUCT.cctk_run_title; \
  const char * const cctk_timer_output = PRIVATE_CACTUS_STRUCT.cctk_timer_output; \
  const char * const info_format = PRIVATE_CACTUS_STRUCT.info_format; \
  const char * const recovery_mode = PRIVATE_CACTUS_STRUCT.recovery_mode; \
  CCTK_INT const allow_mixeddim_gfs = PRIVATE_CACTUS_STRUCT.allow_mixeddim_gfs; \
  CCTK_INT const cctk_brief_output = PRIVATE_CACTUS_STRUCT.cctk_brief_output; \
  CCTK_INT const cctk_full_warnings = PRIVATE_CACTUS_STRUCT.cctk_full_warnings; \
  CCTK_INT const cctk_show_banners = PRIVATE_CACTUS_STRUCT.cctk_show_banners; \
  CCTK_INT const cctk_show_schedule = PRIVATE_CACTUS_STRUCT.cctk_show_schedule; \
  CCTK_INT const cctk_strong_param_check = PRIVATE_CACTUS_STRUCT.cctk_strong_param_check; \
  CCTK_INT const highlight_warning_messages = PRIVATE_CACTUS_STRUCT.highlight_warning_messages; \
  CCTK_INT const manual_cache_setup = PRIVATE_CACTUS_STRUCT.manual_cache_setup; \
  CCTK_INT const manual_cache_size = PRIVATE_CACTUS_STRUCT.manual_cache_size; \
  CCTK_INT const manual_cacheline_bytes = PRIVATE_CACTUS_STRUCT.manual_cacheline_bytes; \
  enum { \
      dummy_PRIVATE_CACTUS_STRUCT_cctk_run_title = sizeof( cctk_run_title ) \
    , dummy_PRIVATE_CACTUS_STRUCT_cctk_timer_output = sizeof( cctk_timer_output ) \
    , dummy_PRIVATE_CACTUS_STRUCT_info_format = sizeof( info_format ) \
    , dummy_PRIVATE_CACTUS_STRUCT_recovery_mode = sizeof( recovery_mode ) \
    , dummy_PRIVATE_CACTUS_STRUCT_allow_mixeddim_gfs = sizeof( allow_mixeddim_gfs ) \
    , dummy_PRIVATE_CACTUS_STRUCT_cctk_brief_output = sizeof( cctk_brief_output ) \
    , dummy_PRIVATE_CACTUS_STRUCT_cctk_full_warnings = sizeof( cctk_full_warnings ) \
    , dummy_PRIVATE_CACTUS_STRUCT_cctk_show_banners = sizeof( cctk_show_banners ) \
    , dummy_PRIVATE_CACTUS_STRUCT_cctk_show_schedule = sizeof( cctk_show_schedule ) \
    , dummy_PRIVATE_CACTUS_STRUCT_cctk_strong_param_check = sizeof( cctk_strong_param_check ) \
    , dummy_PRIVATE_CACTUS_STRUCT_highlight_warning_messages = sizeof( highlight_warning_messages ) \
    , dummy_PRIVATE_CACTUS_STRUCT_manual_cache_setup = sizeof( manual_cache_setup ) \
    , dummy_PRIVATE_CACTUS_STRUCT_manual_cache_size = sizeof( manual_cache_size ) \
    , dummy_PRIVATE_CACTUS_STRUCT_manual_cacheline_bytes = sizeof( manual_cacheline_bytes ) \
  }; \

