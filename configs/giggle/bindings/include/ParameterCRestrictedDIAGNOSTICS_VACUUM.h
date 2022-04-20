#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL BH_Vol_Excise_Radius;
  CCTK_REAL dcostheta;
  CCTK_REAL ddrbddr;
  CCTK_REAL dphi;
  CCTK_REAL drbdr;
  CCTK_REAL inner_volInt_radius;
  CCTK_REAL rbr;
  CCTK_REAL rsurf2;
  CCTK_REAL surf_radius;
  CCTK_INT Compute_VolIntegrands_Every;
  CCTK_INT N_phi;
  CCTK_INT N_theta;
  CCTK_INT WhichIntegral;
  CCTK_INT arithsurf;
  CCTK_INT enable_Jz_constraint;
  CCTK_INT enable_M_constraint;
  CCTK_INT enable_P_constraint;
  CCTK_INT nsurf;
  CCTK_INT ntot;
  CCTK_INT num_BHs;
  CCTK_INT numphi;
  CCTK_INT numtheta;
  CCTK_INT scaledsurf;
  CCTK_INT sym_factor;
} RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_PARAMS \
  CCTK_REAL const BH_Vol_Excise_Radius = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.BH_Vol_Excise_Radius; \
  CCTK_REAL const dcostheta = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.dcostheta; \
  CCTK_REAL const ddrbddr = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.ddrbddr; \
  CCTK_REAL const dphi = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.dphi; \
  CCTK_REAL const drbdr = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.drbdr; \
  CCTK_REAL const inner_volInt_radius = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.inner_volInt_radius; \
  CCTK_REAL const rbr = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.rbr; \
  CCTK_REAL const rsurf2 = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.rsurf2; \
  CCTK_REAL const surf_radius = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.surf_radius; \
  CCTK_INT const Compute_VolIntegrands_Every = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.Compute_VolIntegrands_Every; \
  CCTK_INT const N_phi = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.N_phi; \
  CCTK_INT const N_theta = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.N_theta; \
  CCTK_INT const WhichIntegral = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.WhichIntegral; \
  CCTK_INT const arithsurf = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.arithsurf; \
  CCTK_INT const enable_Jz_constraint = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.enable_Jz_constraint; \
  CCTK_INT const enable_M_constraint = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.enable_M_constraint; \
  CCTK_INT const enable_P_constraint = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.enable_P_constraint; \
  CCTK_INT const nsurf = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.nsurf; \
  CCTK_INT const ntot = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.ntot; \
  CCTK_INT const num_BHs = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.num_BHs; \
  CCTK_INT const numphi = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.numphi; \
  CCTK_INT const numtheta = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.numtheta; \
  CCTK_INT const scaledsurf = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.scaledsurf; \
  CCTK_INT const sym_factor = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.sym_factor; \
  enum { \
      dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_BH_Vol_Excise_Radius = sizeof( BH_Vol_Excise_Radius ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_dcostheta = sizeof( dcostheta ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_ddrbddr = sizeof( ddrbddr ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_dphi = sizeof( dphi ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_drbdr = sizeof( drbdr ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_inner_volInt_radius = sizeof( inner_volInt_radius ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_rbr = sizeof( rbr ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_rsurf2 = sizeof( rsurf2 ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_surf_radius = sizeof( surf_radius ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_Compute_VolIntegrands_Every = sizeof( Compute_VolIntegrands_Every ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_N_phi = sizeof( N_phi ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_N_theta = sizeof( N_theta ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_WhichIntegral = sizeof( WhichIntegral ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_arithsurf = sizeof( arithsurf ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_enable_Jz_constraint = sizeof( enable_Jz_constraint ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_enable_M_constraint = sizeof( enable_M_constraint ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_enable_P_constraint = sizeof( enable_P_constraint ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_nsurf = sizeof( nsurf ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_ntot = sizeof( ntot ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_num_BHs = sizeof( num_BHs ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_numphi = sizeof( numphi ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_numtheta = sizeof( numtheta ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_scaledsurf = sizeof( scaledsurf ) \
    , dummy_RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT_sym_factor = sizeof( sym_factor ) \
  }; \

