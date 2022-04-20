#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL out1D_xline_y;
  CCTK_REAL out1D_xline_z;
  CCTK_REAL out1D_yline_x;
  CCTK_REAL out1D_yline_z;
  CCTK_REAL out1D_zline_x;
  CCTK_REAL out1D_zline_y;
  CCTK_REAL out2D_xyplane_z;
  CCTK_REAL out2D_xzplane_y;
  CCTK_REAL out2D_yzplane_x;
  const char * out1D_dir;
  const char * out1D_style;
  const char * out1D_vars;
  const char * out2D_dir;
  const char * out2D_style;
  const char * out2D_vars;
  const char * out3D_dir;
  const char * out3D_style;
  const char * out3D_vars;
  const char * out_format;
  CCTK_INT out1D_d;
  CCTK_INT out1D_every;
  CCTK_INT out1D_x;
  CCTK_INT out1D_xline_yi;
  CCTK_INT out1D_xline_zi;
  CCTK_INT out1D_y;
  CCTK_INT out1D_yline_xi;
  CCTK_INT out1D_yline_zi;
  CCTK_INT out1D_z;
  CCTK_INT out1D_zline_xi;
  CCTK_INT out1D_zline_yi;
  CCTK_INT out2D_every;
  CCTK_INT out2D_xyplane_zi;
  CCTK_INT out2D_xzplane_yi;
  CCTK_INT out2D_yzplane_xi;
  CCTK_INT out3D_every;
} PRIVATE_IOASCII_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_IOASCII_STRUCT_PARAMS \
  CCTK_REAL const out1D_xline_y = PRIVATE_IOASCII_STRUCT.out1D_xline_y; \
  CCTK_REAL const out1D_xline_z = PRIVATE_IOASCII_STRUCT.out1D_xline_z; \
  CCTK_REAL const out1D_yline_x = PRIVATE_IOASCII_STRUCT.out1D_yline_x; \
  CCTK_REAL const out1D_yline_z = PRIVATE_IOASCII_STRUCT.out1D_yline_z; \
  CCTK_REAL const out1D_zline_x = PRIVATE_IOASCII_STRUCT.out1D_zline_x; \
  CCTK_REAL const out1D_zline_y = PRIVATE_IOASCII_STRUCT.out1D_zline_y; \
  CCTK_REAL const out2D_xyplane_z = PRIVATE_IOASCII_STRUCT.out2D_xyplane_z; \
  CCTK_REAL const out2D_xzplane_y = PRIVATE_IOASCII_STRUCT.out2D_xzplane_y; \
  CCTK_REAL const out2D_yzplane_x = PRIVATE_IOASCII_STRUCT.out2D_yzplane_x; \
  const char * const out1D_dir = PRIVATE_IOASCII_STRUCT.out1D_dir; \
  const char * const out1D_style = PRIVATE_IOASCII_STRUCT.out1D_style; \
  const char * const out1D_vars = PRIVATE_IOASCII_STRUCT.out1D_vars; \
  const char * const out2D_dir = PRIVATE_IOASCII_STRUCT.out2D_dir; \
  const char * const out2D_style = PRIVATE_IOASCII_STRUCT.out2D_style; \
  const char * const out2D_vars = PRIVATE_IOASCII_STRUCT.out2D_vars; \
  const char * const out3D_dir = PRIVATE_IOASCII_STRUCT.out3D_dir; \
  const char * const out3D_style = PRIVATE_IOASCII_STRUCT.out3D_style; \
  const char * const out3D_vars = PRIVATE_IOASCII_STRUCT.out3D_vars; \
  const char * const out_format = PRIVATE_IOASCII_STRUCT.out_format; \
  CCTK_INT const out1D_d = PRIVATE_IOASCII_STRUCT.out1D_d; \
  CCTK_INT const out1D_every = PRIVATE_IOASCII_STRUCT.out1D_every; \
  CCTK_INT const out1D_x = PRIVATE_IOASCII_STRUCT.out1D_x; \
  CCTK_INT const out1D_xline_yi = PRIVATE_IOASCII_STRUCT.out1D_xline_yi; \
  CCTK_INT const out1D_xline_zi = PRIVATE_IOASCII_STRUCT.out1D_xline_zi; \
  CCTK_INT const out1D_y = PRIVATE_IOASCII_STRUCT.out1D_y; \
  CCTK_INT const out1D_yline_xi = PRIVATE_IOASCII_STRUCT.out1D_yline_xi; \
  CCTK_INT const out1D_yline_zi = PRIVATE_IOASCII_STRUCT.out1D_yline_zi; \
  CCTK_INT const out1D_z = PRIVATE_IOASCII_STRUCT.out1D_z; \
  CCTK_INT const out1D_zline_xi = PRIVATE_IOASCII_STRUCT.out1D_zline_xi; \
  CCTK_INT const out1D_zline_yi = PRIVATE_IOASCII_STRUCT.out1D_zline_yi; \
  CCTK_INT const out2D_every = PRIVATE_IOASCII_STRUCT.out2D_every; \
  CCTK_INT const out2D_xyplane_zi = PRIVATE_IOASCII_STRUCT.out2D_xyplane_zi; \
  CCTK_INT const out2D_xzplane_yi = PRIVATE_IOASCII_STRUCT.out2D_xzplane_yi; \
  CCTK_INT const out2D_yzplane_xi = PRIVATE_IOASCII_STRUCT.out2D_yzplane_xi; \
  CCTK_INT const out3D_every = PRIVATE_IOASCII_STRUCT.out3D_every; \
  enum { \
      dummy_PRIVATE_IOASCII_STRUCT_out1D_xline_y = sizeof( out1D_xline_y ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out1D_xline_z = sizeof( out1D_xline_z ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out1D_yline_x = sizeof( out1D_yline_x ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out1D_yline_z = sizeof( out1D_yline_z ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out1D_zline_x = sizeof( out1D_zline_x ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out1D_zline_y = sizeof( out1D_zline_y ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out2D_xyplane_z = sizeof( out2D_xyplane_z ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out2D_xzplane_y = sizeof( out2D_xzplane_y ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out2D_yzplane_x = sizeof( out2D_yzplane_x ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out1D_dir = sizeof( out1D_dir ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out1D_style = sizeof( out1D_style ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out1D_vars = sizeof( out1D_vars ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out2D_dir = sizeof( out2D_dir ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out2D_style = sizeof( out2D_style ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out2D_vars = sizeof( out2D_vars ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out3D_dir = sizeof( out3D_dir ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out3D_style = sizeof( out3D_style ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out3D_vars = sizeof( out3D_vars ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out_format = sizeof( out_format ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out1D_d = sizeof( out1D_d ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out1D_every = sizeof( out1D_every ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out1D_x = sizeof( out1D_x ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out1D_xline_yi = sizeof( out1D_xline_yi ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out1D_xline_zi = sizeof( out1D_xline_zi ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out1D_y = sizeof( out1D_y ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out1D_yline_xi = sizeof( out1D_yline_xi ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out1D_yline_zi = sizeof( out1D_yline_zi ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out1D_z = sizeof( out1D_z ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out1D_zline_xi = sizeof( out1D_zline_xi ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out1D_zline_yi = sizeof( out1D_zline_yi ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out2D_every = sizeof( out2D_every ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out2D_xyplane_zi = sizeof( out2D_xyplane_zi ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out2D_xzplane_yi = sizeof( out2D_xzplane_yi ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out2D_yzplane_xi = sizeof( out2D_yzplane_xi ) \
    , dummy_PRIVATE_IOASCII_STRUCT_out3D_every = sizeof( out3D_every ) \
  }; \

