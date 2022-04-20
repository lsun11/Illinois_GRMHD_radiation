#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL Psi4imsumlw;
  CCTK_REAL Psi4resumlw;
  CCTK_REAL amplitude;
  CCTK_REAL corrector_iteration;
  CCTK_REAL time_shift;
  CCTK_REAL width;
  CCTK_INT mode;
} PRIVATE_LINEARIZED_WAVE_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_LINEARIZED_WAVE_STRUCT_PARAMS \
  CCTK_REAL const Psi4imsumlw = PRIVATE_LINEARIZED_WAVE_STRUCT.Psi4imsumlw; \
  CCTK_REAL const Psi4resumlw = PRIVATE_LINEARIZED_WAVE_STRUCT.Psi4resumlw; \
  CCTK_REAL const amplitude = PRIVATE_LINEARIZED_WAVE_STRUCT.amplitude; \
  CCTK_REAL const corrector_iteration = PRIVATE_LINEARIZED_WAVE_STRUCT.corrector_iteration; \
  CCTK_REAL const time_shift = PRIVATE_LINEARIZED_WAVE_STRUCT.time_shift; \
  CCTK_REAL const width = PRIVATE_LINEARIZED_WAVE_STRUCT.width; \
  CCTK_INT const mode = PRIVATE_LINEARIZED_WAVE_STRUCT.mode; \
  enum { \
      dummy_PRIVATE_LINEARIZED_WAVE_STRUCT_Psi4imsumlw = sizeof( Psi4imsumlw ) \
    , dummy_PRIVATE_LINEARIZED_WAVE_STRUCT_Psi4resumlw = sizeof( Psi4resumlw ) \
    , dummy_PRIVATE_LINEARIZED_WAVE_STRUCT_amplitude = sizeof( amplitude ) \
    , dummy_PRIVATE_LINEARIZED_WAVE_STRUCT_corrector_iteration = sizeof( corrector_iteration ) \
    , dummy_PRIVATE_LINEARIZED_WAVE_STRUCT_time_shift = sizeof( time_shift ) \
    , dummy_PRIVATE_LINEARIZED_WAVE_STRUCT_width = sizeof( width ) \
    , dummy_PRIVATE_LINEARIZED_WAVE_STRUCT_mode = sizeof( mode ) \
  }; \

