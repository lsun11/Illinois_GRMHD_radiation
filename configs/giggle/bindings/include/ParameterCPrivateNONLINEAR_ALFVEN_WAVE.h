#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL B1os4pi_l;
  CCTK_REAL B2os4pi_l;
  CCTK_REAL P0;
  CCTK_REAL npolywave;
  CCTK_REAL nx;
  CCTK_REAL ny;
  CCTK_REAL nz;
  CCTK_REAL rho_b0;
  CCTK_REAL rotations;
  CCTK_REAL u1_l;
  CCTK_REAL u2_l;
  CCTK_REAL u3_l;
  CCTK_REAL wtrans;
} PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT_PARAMS \
  CCTK_REAL const B1os4pi_l = PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT.B1os4pi_l; \
  CCTK_REAL const B2os4pi_l = PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT.B2os4pi_l; \
  CCTK_REAL const P0 = PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT.P0; \
  CCTK_REAL const npolywave = PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT.npolywave; \
  CCTK_REAL const nx = PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT.nx; \
  CCTK_REAL const ny = PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT.ny; \
  CCTK_REAL const nz = PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT.nz; \
  CCTK_REAL const rho_b0 = PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT.rho_b0; \
  CCTK_REAL const rotations = PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT.rotations; \
  CCTK_REAL const u1_l = PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT.u1_l; \
  CCTK_REAL const u2_l = PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT.u2_l; \
  CCTK_REAL const u3_l = PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT.u3_l; \
  CCTK_REAL const wtrans = PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT.wtrans; \
  enum { \
      dummy_PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT_B1os4pi_l = sizeof( B1os4pi_l ) \
    , dummy_PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT_B2os4pi_l = sizeof( B2os4pi_l ) \
    , dummy_PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT_P0 = sizeof( P0 ) \
    , dummy_PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT_npolywave = sizeof( npolywave ) \
    , dummy_PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT_nx = sizeof( nx ) \
    , dummy_PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT_ny = sizeof( ny ) \
    , dummy_PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT_nz = sizeof( nz ) \
    , dummy_PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT_rho_b0 = sizeof( rho_b0 ) \
    , dummy_PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT_rotations = sizeof( rotations ) \
    , dummy_PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT_u1_l = sizeof( u1_l ) \
    , dummy_PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT_u2_l = sizeof( u2_l ) \
    , dummy_PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT_u3_l = sizeof( u3_l ) \
    , dummy_PRIVATE_NONLINEAR_ALFVEN_WAVE_STRUCT_wtrans = sizeof( wtrans ) \
  }; \

