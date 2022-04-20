#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL amplitude;
  CCTK_REAL origin[3];
  CCTK_REAL phase_offset[3];
  CCTK_REAL pulse_direction[3];
  CCTK_REAL pulse_offset[3];
  CCTK_REAL radius;
  CCTK_REAL sigma;
  CCTK_REAL time_offset;
  CCTK_REAL wave_number[3];
  const char * initial_data;
} PRIVATE_IDSCALARWAVEMOL_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_IDSCALARWAVEMOL_STRUCT_PARAMS \
  CCTK_REAL const amplitude = PRIVATE_IDSCALARWAVEMOL_STRUCT.amplitude; \
  CCTK_REAL const * const origin = PRIVATE_IDSCALARWAVEMOL_STRUCT.origin; \
  CCTK_REAL const * const phase_offset = PRIVATE_IDSCALARWAVEMOL_STRUCT.phase_offset; \
  CCTK_REAL const * const pulse_direction = PRIVATE_IDSCALARWAVEMOL_STRUCT.pulse_direction; \
  CCTK_REAL const * const pulse_offset = PRIVATE_IDSCALARWAVEMOL_STRUCT.pulse_offset; \
  CCTK_REAL const radius = PRIVATE_IDSCALARWAVEMOL_STRUCT.radius; \
  CCTK_REAL const sigma = PRIVATE_IDSCALARWAVEMOL_STRUCT.sigma; \
  CCTK_REAL const time_offset = PRIVATE_IDSCALARWAVEMOL_STRUCT.time_offset; \
  CCTK_REAL const * const wave_number = PRIVATE_IDSCALARWAVEMOL_STRUCT.wave_number; \
  const char * const initial_data = PRIVATE_IDSCALARWAVEMOL_STRUCT.initial_data; \
  enum { \
      dummy_PRIVATE_IDSCALARWAVEMOL_STRUCT_amplitude = sizeof( amplitude ) \
    , dummy_PRIVATE_IDSCALARWAVEMOL_STRUCT_origin = sizeof( origin ) \
    , dummy_PRIVATE_IDSCALARWAVEMOL_STRUCT_phase_offset = sizeof( phase_offset ) \
    , dummy_PRIVATE_IDSCALARWAVEMOL_STRUCT_pulse_direction = sizeof( pulse_direction ) \
    , dummy_PRIVATE_IDSCALARWAVEMOL_STRUCT_pulse_offset = sizeof( pulse_offset ) \
    , dummy_PRIVATE_IDSCALARWAVEMOL_STRUCT_radius = sizeof( radius ) \
    , dummy_PRIVATE_IDSCALARWAVEMOL_STRUCT_sigma = sizeof( sigma ) \
    , dummy_PRIVATE_IDSCALARWAVEMOL_STRUCT_time_offset = sizeof( time_offset ) \
    , dummy_PRIVATE_IDSCALARWAVEMOL_STRUCT_wave_number = sizeof( wave_number ) \
    , dummy_PRIVATE_IDSCALARWAVEMOL_STRUCT_initial_data = sizeof( initial_data ) \
  }; \

