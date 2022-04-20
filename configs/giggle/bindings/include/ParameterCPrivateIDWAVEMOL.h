#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL amplitude;
  CCTK_REAL centrex;
  CCTK_REAL centrey;
  CCTK_REAL centrez;
  CCTK_REAL kx;
  CCTK_REAL ky;
  CCTK_REAL kz;
  CCTK_REAL offsett;
  CCTK_REAL offsetx;
  CCTK_REAL offsety;
  CCTK_REAL offsetz;
  CCTK_REAL radius;
  CCTK_REAL sigma;
  CCTK_REAL slopet;
  CCTK_REAL slopex;
  CCTK_REAL slopey;
  CCTK_REAL slopez;
} PRIVATE_IDWAVEMOL_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_IDWAVEMOL_STRUCT_PARAMS \
  CCTK_REAL const amplitude = PRIVATE_IDWAVEMOL_STRUCT.amplitude; \
  CCTK_REAL const centrex = PRIVATE_IDWAVEMOL_STRUCT.centrex; \
  CCTK_REAL const centrey = PRIVATE_IDWAVEMOL_STRUCT.centrey; \
  CCTK_REAL const centrez = PRIVATE_IDWAVEMOL_STRUCT.centrez; \
  CCTK_REAL const kx = PRIVATE_IDWAVEMOL_STRUCT.kx; \
  CCTK_REAL const ky = PRIVATE_IDWAVEMOL_STRUCT.ky; \
  CCTK_REAL const kz = PRIVATE_IDWAVEMOL_STRUCT.kz; \
  CCTK_REAL const offsett = PRIVATE_IDWAVEMOL_STRUCT.offsett; \
  CCTK_REAL const offsetx = PRIVATE_IDWAVEMOL_STRUCT.offsetx; \
  CCTK_REAL const offsety = PRIVATE_IDWAVEMOL_STRUCT.offsety; \
  CCTK_REAL const offsetz = PRIVATE_IDWAVEMOL_STRUCT.offsetz; \
  CCTK_REAL const radius = PRIVATE_IDWAVEMOL_STRUCT.radius; \
  CCTK_REAL const sigma = PRIVATE_IDWAVEMOL_STRUCT.sigma; \
  CCTK_REAL const slopet = PRIVATE_IDWAVEMOL_STRUCT.slopet; \
  CCTK_REAL const slopex = PRIVATE_IDWAVEMOL_STRUCT.slopex; \
  CCTK_REAL const slopey = PRIVATE_IDWAVEMOL_STRUCT.slopey; \
  CCTK_REAL const slopez = PRIVATE_IDWAVEMOL_STRUCT.slopez; \
  enum { \
      dummy_PRIVATE_IDWAVEMOL_STRUCT_amplitude = sizeof( amplitude ) \
    , dummy_PRIVATE_IDWAVEMOL_STRUCT_centrex = sizeof( centrex ) \
    , dummy_PRIVATE_IDWAVEMOL_STRUCT_centrey = sizeof( centrey ) \
    , dummy_PRIVATE_IDWAVEMOL_STRUCT_centrez = sizeof( centrez ) \
    , dummy_PRIVATE_IDWAVEMOL_STRUCT_kx = sizeof( kx ) \
    , dummy_PRIVATE_IDWAVEMOL_STRUCT_ky = sizeof( ky ) \
    , dummy_PRIVATE_IDWAVEMOL_STRUCT_kz = sizeof( kz ) \
    , dummy_PRIVATE_IDWAVEMOL_STRUCT_offsett = sizeof( offsett ) \
    , dummy_PRIVATE_IDWAVEMOL_STRUCT_offsetx = sizeof( offsetx ) \
    , dummy_PRIVATE_IDWAVEMOL_STRUCT_offsety = sizeof( offsety ) \
    , dummy_PRIVATE_IDWAVEMOL_STRUCT_offsetz = sizeof( offsetz ) \
    , dummy_PRIVATE_IDWAVEMOL_STRUCT_radius = sizeof( radius ) \
    , dummy_PRIVATE_IDWAVEMOL_STRUCT_sigma = sizeof( sigma ) \
    , dummy_PRIVATE_IDWAVEMOL_STRUCT_slopet = sizeof( slopet ) \
    , dummy_PRIVATE_IDWAVEMOL_STRUCT_slopex = sizeof( slopex ) \
    , dummy_PRIVATE_IDWAVEMOL_STRUCT_slopey = sizeof( slopey ) \
    , dummy_PRIVATE_IDWAVEMOL_STRUCT_slopez = sizeof( slopez ) \
  }; \

