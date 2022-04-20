#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL amplitude;
  CCTK_REAL width;
  const char * bound;
  CCTK_INT enable_moving_grid;
  CCTK_INT scalarwave_Symmetry;
} PRIVATE_SCALARWAVEMOL_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_SCALARWAVEMOL_STRUCT_PARAMS \
  CCTK_REAL const amplitude = PRIVATE_SCALARWAVEMOL_STRUCT.amplitude; \
  CCTK_REAL const width = PRIVATE_SCALARWAVEMOL_STRUCT.width; \
  const char * const bound = PRIVATE_SCALARWAVEMOL_STRUCT.bound; \
  CCTK_INT const enable_moving_grid = PRIVATE_SCALARWAVEMOL_STRUCT.enable_moving_grid; \
  CCTK_INT const scalarwave_Symmetry = PRIVATE_SCALARWAVEMOL_STRUCT.scalarwave_Symmetry; \
  enum { \
      dummy_PRIVATE_SCALARWAVEMOL_STRUCT_amplitude = sizeof( amplitude ) \
    , dummy_PRIVATE_SCALARWAVEMOL_STRUCT_width = sizeof( width ) \
    , dummy_PRIVATE_SCALARWAVEMOL_STRUCT_bound = sizeof( bound ) \
    , dummy_PRIVATE_SCALARWAVEMOL_STRUCT_enable_moving_grid = sizeof( enable_moving_grid ) \
    , dummy_PRIVATE_SCALARWAVEMOL_STRUCT_scalarwave_Symmetry = sizeof( scalarwave_Symmetry ) \
  }; \

