#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL courant_fac;
  CCTK_REAL dtfac;
  CCTK_REAL timestep;
  CCTK_INT timestep_outevery;
  CCTK_INT verbose;
} PRIVATE_TIME_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_TIME_STRUCT_PARAMS \
  CCTK_REAL const courant_fac = PRIVATE_TIME_STRUCT.courant_fac; \
  CCTK_REAL const dtfac = PRIVATE_TIME_STRUCT.dtfac; \
  CCTK_REAL const timestep = PRIVATE_TIME_STRUCT.timestep; \
  CCTK_INT const timestep_outevery = PRIVATE_TIME_STRUCT.timestep_outevery; \
  CCTK_INT const verbose = PRIVATE_TIME_STRUCT.verbose; \
  enum { \
      dummy_PRIVATE_TIME_STRUCT_courant_fac = sizeof( courant_fac ) \
    , dummy_PRIVATE_TIME_STRUCT_dtfac = sizeof( dtfac ) \
    , dummy_PRIVATE_TIME_STRUCT_timestep = sizeof( timestep ) \
    , dummy_PRIVATE_TIME_STRUCT_timestep_outevery = sizeof( timestep_outevery ) \
    , dummy_PRIVATE_TIME_STRUCT_verbose = sizeof( verbose ) \
  }; \

