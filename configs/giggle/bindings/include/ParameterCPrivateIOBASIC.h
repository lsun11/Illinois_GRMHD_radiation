#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL outInfo_dt;
  CCTK_REAL outScalar_dt;
  const char * outInfo_criterion;
  const char * outInfo_reductions;
  const char * outInfo_vars;
  const char * outScalar_criterion;
  const char * outScalar_reductions;
  const char * outScalar_style;
  const char * outScalar_vars;
  const char * out_dir;
  const char * out_format;
  CCTK_INT outInfo_every;
  CCTK_INT outScalar_every;
} PRIVATE_IOBASIC_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_IOBASIC_STRUCT_PARAMS \
  CCTK_REAL const outInfo_dt = PRIVATE_IOBASIC_STRUCT.outInfo_dt; \
  CCTK_REAL const outScalar_dt = PRIVATE_IOBASIC_STRUCT.outScalar_dt; \
  const char * const outInfo_criterion = PRIVATE_IOBASIC_STRUCT.outInfo_criterion; \
  const char * const outInfo_reductions = PRIVATE_IOBASIC_STRUCT.outInfo_reductions; \
  const char * const outInfo_vars = PRIVATE_IOBASIC_STRUCT.outInfo_vars; \
  const char * const outScalar_criterion = PRIVATE_IOBASIC_STRUCT.outScalar_criterion; \
  const char * const outScalar_reductions = PRIVATE_IOBASIC_STRUCT.outScalar_reductions; \
  const char * const outScalar_style = PRIVATE_IOBASIC_STRUCT.outScalar_style; \
  const char * const outScalar_vars = PRIVATE_IOBASIC_STRUCT.outScalar_vars; \
  const char * const out_dir = PRIVATE_IOBASIC_STRUCT.out_dir; \
  const char * const out_format = PRIVATE_IOBASIC_STRUCT.out_format; \
  CCTK_INT const outInfo_every = PRIVATE_IOBASIC_STRUCT.outInfo_every; \
  CCTK_INT const outScalar_every = PRIVATE_IOBASIC_STRUCT.outScalar_every; \
  enum { \
      dummy_PRIVATE_IOBASIC_STRUCT_outInfo_dt = sizeof( outInfo_dt ) \
    , dummy_PRIVATE_IOBASIC_STRUCT_outScalar_dt = sizeof( outScalar_dt ) \
    , dummy_PRIVATE_IOBASIC_STRUCT_outInfo_criterion = sizeof( outInfo_criterion ) \
    , dummy_PRIVATE_IOBASIC_STRUCT_outInfo_reductions = sizeof( outInfo_reductions ) \
    , dummy_PRIVATE_IOBASIC_STRUCT_outInfo_vars = sizeof( outInfo_vars ) \
    , dummy_PRIVATE_IOBASIC_STRUCT_outScalar_criterion = sizeof( outScalar_criterion ) \
    , dummy_PRIVATE_IOBASIC_STRUCT_outScalar_reductions = sizeof( outScalar_reductions ) \
    , dummy_PRIVATE_IOBASIC_STRUCT_outScalar_style = sizeof( outScalar_style ) \
    , dummy_PRIVATE_IOBASIC_STRUCT_outScalar_vars = sizeof( outScalar_vars ) \
    , dummy_PRIVATE_IOBASIC_STRUCT_out_dir = sizeof( out_dir ) \
    , dummy_PRIVATE_IOBASIC_STRUCT_out_format = sizeof( out_format ) \
    , dummy_PRIVATE_IOBASIC_STRUCT_outInfo_every = sizeof( outInfo_every ) \
    , dummy_PRIVATE_IOBASIC_STRUCT_outScalar_every = sizeof( outScalar_every ) \
  }; \

