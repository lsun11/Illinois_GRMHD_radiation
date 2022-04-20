#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_INT truncate_files;
  CCTK_INT truncate_files_after_recovering;
} PRIVATE_IOUTIL_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_IOUTIL_STRUCT_PARAMS \
  CCTK_INT const truncate_files = PRIVATE_IOUTIL_STRUCT.truncate_files; \
  CCTK_INT const truncate_files_after_recovering = PRIVATE_IOUTIL_STRUCT.truncate_files_after_recovering; \
  enum { \
      dummy_PRIVATE_IOUTIL_STRUCT_truncate_files = sizeof( truncate_files ) \
    , dummy_PRIVATE_IOUTIL_STRUCT_truncate_files_after_recovering = sizeof( truncate_files_after_recovering ) \
  }; \

