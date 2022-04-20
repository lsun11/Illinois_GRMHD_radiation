#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_INT verbose;
} PRIVATE_CARPETREDUCE_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_CARPETREDUCE_STRUCT_PARAMS \
  CCTK_INT const verbose = PRIVATE_CARPETREDUCE_STRUCT.verbose; \
  enum { \
      dummy_PRIVATE_CARPETREDUCE_STRUCT_verbose = sizeof( verbose ) \
  }; \

