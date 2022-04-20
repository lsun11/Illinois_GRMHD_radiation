#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL B1os4pi_l;
  CCTK_REAL B1os4pi_r;
  CCTK_REAL B2os4pi_l;
  CCTK_REAL B2os4pi_r;
  CCTK_REAL B3os4pi_l;
  CCTK_REAL B3os4pi_r;
  CCTK_REAL Erad_l;
  CCTK_REAL Erad_r;
  CCTK_REAL Frad0_l;
  CCTK_REAL Frad0_r;
  CCTK_REAL Fradx_l;
  CCTK_REAL Fradx_r;
  CCTK_REAL M0_initial;
  CCTK_REAL P_l;
  CCTK_REAL P_r;
  CCTK_REAL npolyshock;
  CCTK_REAL nx;
  CCTK_REAL ny;
  CCTK_REAL nz;
  CCTK_REAL rho_l;
  CCTK_REAL rho_r;
  CCTK_REAL u1_l;
  CCTK_REAL u1_r;
  CCTK_REAL u2_l;
  CCTK_REAL u2_r;
  CCTK_REAL u3_l;
  CCTK_REAL u3_r;
  CCTK_REAL u_boost;
  CCTK_INT NS;
  CCTK_INT boost_shock;
  CCTK_INT integrate_flag;
  CCTK_INT read_radiation_shock_data;
} PRIVATE_SHOCK_MHD_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_SHOCK_MHD_STRUCT_PARAMS \
  CCTK_REAL const B1os4pi_l = PRIVATE_SHOCK_MHD_STRUCT.B1os4pi_l; \
  CCTK_REAL const B1os4pi_r = PRIVATE_SHOCK_MHD_STRUCT.B1os4pi_r; \
  CCTK_REAL const B2os4pi_l = PRIVATE_SHOCK_MHD_STRUCT.B2os4pi_l; \
  CCTK_REAL const B2os4pi_r = PRIVATE_SHOCK_MHD_STRUCT.B2os4pi_r; \
  CCTK_REAL const B3os4pi_l = PRIVATE_SHOCK_MHD_STRUCT.B3os4pi_l; \
  CCTK_REAL const B3os4pi_r = PRIVATE_SHOCK_MHD_STRUCT.B3os4pi_r; \
  CCTK_REAL const Erad_l = PRIVATE_SHOCK_MHD_STRUCT.Erad_l; \
  CCTK_REAL const Erad_r = PRIVATE_SHOCK_MHD_STRUCT.Erad_r; \
  CCTK_REAL const Frad0_l = PRIVATE_SHOCK_MHD_STRUCT.Frad0_l; \
  CCTK_REAL const Frad0_r = PRIVATE_SHOCK_MHD_STRUCT.Frad0_r; \
  CCTK_REAL const Fradx_l = PRIVATE_SHOCK_MHD_STRUCT.Fradx_l; \
  CCTK_REAL const Fradx_r = PRIVATE_SHOCK_MHD_STRUCT.Fradx_r; \
  CCTK_REAL const M0_initial = PRIVATE_SHOCK_MHD_STRUCT.M0_initial; \
  CCTK_REAL const P_l = PRIVATE_SHOCK_MHD_STRUCT.P_l; \
  CCTK_REAL const P_r = PRIVATE_SHOCK_MHD_STRUCT.P_r; \
  CCTK_REAL const npolyshock = PRIVATE_SHOCK_MHD_STRUCT.npolyshock; \
  CCTK_REAL const nx = PRIVATE_SHOCK_MHD_STRUCT.nx; \
  CCTK_REAL const ny = PRIVATE_SHOCK_MHD_STRUCT.ny; \
  CCTK_REAL const nz = PRIVATE_SHOCK_MHD_STRUCT.nz; \
  CCTK_REAL const rho_l = PRIVATE_SHOCK_MHD_STRUCT.rho_l; \
  CCTK_REAL const rho_r = PRIVATE_SHOCK_MHD_STRUCT.rho_r; \
  CCTK_REAL const u1_l = PRIVATE_SHOCK_MHD_STRUCT.u1_l; \
  CCTK_REAL const u1_r = PRIVATE_SHOCK_MHD_STRUCT.u1_r; \
  CCTK_REAL const u2_l = PRIVATE_SHOCK_MHD_STRUCT.u2_l; \
  CCTK_REAL const u2_r = PRIVATE_SHOCK_MHD_STRUCT.u2_r; \
  CCTK_REAL const u3_l = PRIVATE_SHOCK_MHD_STRUCT.u3_l; \
  CCTK_REAL const u3_r = PRIVATE_SHOCK_MHD_STRUCT.u3_r; \
  CCTK_REAL const u_boost = PRIVATE_SHOCK_MHD_STRUCT.u_boost; \
  CCTK_INT const NS = PRIVATE_SHOCK_MHD_STRUCT.NS; \
  CCTK_INT const boost_shock = PRIVATE_SHOCK_MHD_STRUCT.boost_shock; \
  CCTK_INT const integrate_flag = PRIVATE_SHOCK_MHD_STRUCT.integrate_flag; \
  CCTK_INT const read_radiation_shock_data = PRIVATE_SHOCK_MHD_STRUCT.read_radiation_shock_data; \
  enum { \
      dummy_PRIVATE_SHOCK_MHD_STRUCT_B1os4pi_l = sizeof( B1os4pi_l ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_B1os4pi_r = sizeof( B1os4pi_r ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_B2os4pi_l = sizeof( B2os4pi_l ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_B2os4pi_r = sizeof( B2os4pi_r ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_B3os4pi_l = sizeof( B3os4pi_l ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_B3os4pi_r = sizeof( B3os4pi_r ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_Erad_l = sizeof( Erad_l ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_Erad_r = sizeof( Erad_r ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_Frad0_l = sizeof( Frad0_l ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_Frad0_r = sizeof( Frad0_r ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_Fradx_l = sizeof( Fradx_l ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_Fradx_r = sizeof( Fradx_r ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_M0_initial = sizeof( M0_initial ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_P_l = sizeof( P_l ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_P_r = sizeof( P_r ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_npolyshock = sizeof( npolyshock ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_nx = sizeof( nx ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_ny = sizeof( ny ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_nz = sizeof( nz ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_rho_l = sizeof( rho_l ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_rho_r = sizeof( rho_r ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_u1_l = sizeof( u1_l ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_u1_r = sizeof( u1_r ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_u2_l = sizeof( u2_l ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_u2_r = sizeof( u2_r ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_u3_l = sizeof( u3_l ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_u3_r = sizeof( u3_r ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_u_boost = sizeof( u_boost ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_NS = sizeof( NS ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_boost_shock = sizeof( boost_shock ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_integrate_flag = sizeof( integrate_flag ) \
    , dummy_PRIVATE_SHOCK_MHD_STRUCT_read_radiation_shock_data = sizeof( read_radiation_shock_data ) \
  }; \

