#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL outInfo_dt;
  CCTK_REAL real_max;
  CCTK_REAL real_min;
  const char * outInfo_criterion;
  const char * outInfo_reductions;
  const char * outInfo_vars;
  CCTK_INT int_width;
  CCTK_INT iter_width;
  CCTK_INT outHeader_every;
  CCTK_INT outInfo_every;
  CCTK_INT real_prec;
  CCTK_INT real_prec_sci;
  CCTK_INT real_width;
  CCTK_INT time_prec;
  CCTK_INT time_width;
} PRIVATE_CARPETIOBASIC_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_CARPETIOBASIC_STRUCT_PARAMS \
  CCTK_REAL const outInfo_dt = PRIVATE_CARPETIOBASIC_STRUCT.outInfo_dt; \
  CCTK_REAL const real_max = PRIVATE_CARPETIOBASIC_STRUCT.real_max; \
  CCTK_REAL const real_min = PRIVATE_CARPETIOBASIC_STRUCT.real_min; \
  const char * const outInfo_criterion = PRIVATE_CARPETIOBASIC_STRUCT.outInfo_criterion; \
  const char * const outInfo_reductions = PRIVATE_CARPETIOBASIC_STRUCT.outInfo_reductions; \
  const char * const outInfo_vars = PRIVATE_CARPETIOBASIC_STRUCT.outInfo_vars; \
  CCTK_INT const int_width = PRIVATE_CARPETIOBASIC_STRUCT.int_width; \
  CCTK_INT const iter_width = PRIVATE_CARPETIOBASIC_STRUCT.iter_width; \
  CCTK_INT const outHeader_every = PRIVATE_CARPETIOBASIC_STRUCT.outHeader_every; \
  CCTK_INT const outInfo_every = PRIVATE_CARPETIOBASIC_STRUCT.outInfo_every; \
  CCTK_INT const real_prec = PRIVATE_CARPETIOBASIC_STRUCT.real_prec; \
  CCTK_INT const real_prec_sci = PRIVATE_CARPETIOBASIC_STRUCT.real_prec_sci; \
  CCTK_INT const real_width = PRIVATE_CARPETIOBASIC_STRUCT.real_width; \
  CCTK_INT const time_prec = PRIVATE_CARPETIOBASIC_STRUCT.time_prec; \
  CCTK_INT const time_width = PRIVATE_CARPETIOBASIC_STRUCT.time_width; \
  enum { \
      dummy_PRIVATE_CARPETIOBASIC_STRUCT_outInfo_dt = sizeof( outInfo_dt ) \
    , dummy_PRIVATE_CARPETIOBASIC_STRUCT_real_max = sizeof( real_max ) \
    , dummy_PRIVATE_CARPETIOBASIC_STRUCT_real_min = sizeof( real_min ) \
    , dummy_PRIVATE_CARPETIOBASIC_STRUCT_outInfo_criterion = sizeof( outInfo_criterion ) \
    , dummy_PRIVATE_CARPETIOBASIC_STRUCT_outInfo_reductions = sizeof( outInfo_reductions ) \
    , dummy_PRIVATE_CARPETIOBASIC_STRUCT_outInfo_vars = sizeof( outInfo_vars ) \
    , dummy_PRIVATE_CARPETIOBASIC_STRUCT_int_width = sizeof( int_width ) \
    , dummy_PRIVATE_CARPETIOBASIC_STRUCT_iter_width = sizeof( iter_width ) \
    , dummy_PRIVATE_CARPETIOBASIC_STRUCT_outHeader_every = sizeof( outHeader_every ) \
    , dummy_PRIVATE_CARPETIOBASIC_STRUCT_outInfo_every = sizeof( outInfo_every ) \
    , dummy_PRIVATE_CARPETIOBASIC_STRUCT_real_prec = sizeof( real_prec ) \
    , dummy_PRIVATE_CARPETIOBASIC_STRUCT_real_prec_sci = sizeof( real_prec_sci ) \
    , dummy_PRIVATE_CARPETIOBASIC_STRUCT_real_width = sizeof( real_width ) \
    , dummy_PRIVATE_CARPETIOBASIC_STRUCT_time_prec = sizeof( time_prec ) \
    , dummy_PRIVATE_CARPETIOBASIC_STRUCT_time_width = sizeof( time_width ) \
  }; \

