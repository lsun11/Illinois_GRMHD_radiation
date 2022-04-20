#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL offsetx[10];
  CCTK_REAL offsety[10];
  CCTK_REAL offsetz[10];
  CCTK_INT num_offsets;
  CCTK_INT offset_firstlevel;
} RESTRICTED_CARPETREGRID_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_CARPETREGRID_STRUCT_PARAMS \
  CCTK_REAL const * const offsetx = RESTRICTED_CARPETREGRID_STRUCT.offsetx; \
  CCTK_REAL const * const offsety = RESTRICTED_CARPETREGRID_STRUCT.offsety; \
  CCTK_REAL const * const offsetz = RESTRICTED_CARPETREGRID_STRUCT.offsetz; \
  CCTK_INT const num_offsets = RESTRICTED_CARPETREGRID_STRUCT.num_offsets; \
  CCTK_INT const offset_firstlevel = RESTRICTED_CARPETREGRID_STRUCT.offset_firstlevel; \
  enum { \
      dummy_RESTRICTED_CARPETREGRID_STRUCT_offsetx = sizeof( offsetx ) \
    , dummy_RESTRICTED_CARPETREGRID_STRUCT_offsety = sizeof( offsety ) \
    , dummy_RESTRICTED_CARPETREGRID_STRUCT_offsetz = sizeof( offsetz ) \
    , dummy_RESTRICTED_CARPETREGRID_STRUCT_num_offsets = sizeof( num_offsets ) \
    , dummy_RESTRICTED_CARPETREGRID_STRUCT_offset_firstlevel = sizeof( offset_firstlevel ) \
  }; \

