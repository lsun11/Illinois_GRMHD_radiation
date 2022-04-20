#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_INT surface[10];
  CCTK_INT verbose;
} PRIVATE_CARPETTRACKER_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_CARPETTRACKER_STRUCT_PARAMS \
  CCTK_INT const * const surface = PRIVATE_CARPETTRACKER_STRUCT.surface; \
  CCTK_INT const verbose = PRIVATE_CARPETTRACKER_STRUCT.verbose; \
  enum { \
      dummy_PRIVATE_CARPETTRACKER_STRUCT_surface = sizeof( surface ) \
    , dummy_PRIVATE_CARPETTRACKER_STRUCT_verbose = sizeof( verbose ) \
  }; \

