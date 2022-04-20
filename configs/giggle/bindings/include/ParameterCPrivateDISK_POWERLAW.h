#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL A_const;
  CCTK_REAL BigMass;
  CCTK_REAL RoM_inner;
  CCTK_REAL bbh_r_const1;
  CCTK_REAL bbh_r_const2;
  CCTK_REAL bbh_r_const3;
  CCTK_REAL bbh_r_const4;
  CCTK_REAL bbh_rhosurf1;
  CCTK_REAL bbh_rhosurf2;
  CCTK_REAL bbh_rhosurf3;
  CCTK_REAL bbh_rhosurf4;
  CCTK_REAL betam1;
  CCTK_REAL disk_powerlaw_P_max;
  CCTK_REAL disk_powerlaw_rho_star_max;
  CCTK_REAL disk_powerlaw_tau_max;
  CCTK_REAL ell_inner;
  CCTK_REAL p_c;
  CCTK_REAL q;
  CCTK_REAL rho_cut;
  CCTK_REAL rho_cut_o_rho_max;
  CCTK_REAL sam_disk;
  CCTK_REAL surfxmax;
  CCTK_REAL surfxmin;
  CCTK_REAL surfzmax;
  CCTK_INT arithrhosurf;
  CCTK_INT em_field_type;
  CCTK_INT metric_type;
  CCTK_INT nsurfrho;
  CCTK_INT ntot_disk_powerlaw;
  CCTK_INT surfphinum;
  CCTK_INT surfxnum;
  CCTK_INT surfznum;
} PRIVATE_DISK_POWERLAW_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_DISK_POWERLAW_STRUCT_PARAMS \
  CCTK_REAL const A_const = PRIVATE_DISK_POWERLAW_STRUCT.A_const; \
  CCTK_REAL const BigMass = PRIVATE_DISK_POWERLAW_STRUCT.BigMass; \
  CCTK_REAL const RoM_inner = PRIVATE_DISK_POWERLAW_STRUCT.RoM_inner; \
  CCTK_REAL const bbh_r_const1 = PRIVATE_DISK_POWERLAW_STRUCT.bbh_r_const1; \
  CCTK_REAL const bbh_r_const2 = PRIVATE_DISK_POWERLAW_STRUCT.bbh_r_const2; \
  CCTK_REAL const bbh_r_const3 = PRIVATE_DISK_POWERLAW_STRUCT.bbh_r_const3; \
  CCTK_REAL const bbh_r_const4 = PRIVATE_DISK_POWERLAW_STRUCT.bbh_r_const4; \
  CCTK_REAL const bbh_rhosurf1 = PRIVATE_DISK_POWERLAW_STRUCT.bbh_rhosurf1; \
  CCTK_REAL const bbh_rhosurf2 = PRIVATE_DISK_POWERLAW_STRUCT.bbh_rhosurf2; \
  CCTK_REAL const bbh_rhosurf3 = PRIVATE_DISK_POWERLAW_STRUCT.bbh_rhosurf3; \
  CCTK_REAL const bbh_rhosurf4 = PRIVATE_DISK_POWERLAW_STRUCT.bbh_rhosurf4; \
  CCTK_REAL const betam1 = PRIVATE_DISK_POWERLAW_STRUCT.betam1; \
  CCTK_REAL const disk_powerlaw_P_max = PRIVATE_DISK_POWERLAW_STRUCT.disk_powerlaw_P_max; \
  CCTK_REAL const disk_powerlaw_rho_star_max = PRIVATE_DISK_POWERLAW_STRUCT.disk_powerlaw_rho_star_max; \
  CCTK_REAL const disk_powerlaw_tau_max = PRIVATE_DISK_POWERLAW_STRUCT.disk_powerlaw_tau_max; \
  CCTK_REAL const ell_inner = PRIVATE_DISK_POWERLAW_STRUCT.ell_inner; \
  CCTK_REAL const p_c = PRIVATE_DISK_POWERLAW_STRUCT.p_c; \
  CCTK_REAL const q = PRIVATE_DISK_POWERLAW_STRUCT.q; \
  CCTK_REAL const rho_cut = PRIVATE_DISK_POWERLAW_STRUCT.rho_cut; \
  CCTK_REAL const rho_cut_o_rho_max = PRIVATE_DISK_POWERLAW_STRUCT.rho_cut_o_rho_max; \
  CCTK_REAL const sam_disk = PRIVATE_DISK_POWERLAW_STRUCT.sam_disk; \
  CCTK_REAL const surfxmax = PRIVATE_DISK_POWERLAW_STRUCT.surfxmax; \
  CCTK_REAL const surfxmin = PRIVATE_DISK_POWERLAW_STRUCT.surfxmin; \
  CCTK_REAL const surfzmax = PRIVATE_DISK_POWERLAW_STRUCT.surfzmax; \
  CCTK_INT const arithrhosurf = PRIVATE_DISK_POWERLAW_STRUCT.arithrhosurf; \
  CCTK_INT const em_field_type = PRIVATE_DISK_POWERLAW_STRUCT.em_field_type; \
  CCTK_INT const metric_type = PRIVATE_DISK_POWERLAW_STRUCT.metric_type; \
  CCTK_INT const nsurfrho = PRIVATE_DISK_POWERLAW_STRUCT.nsurfrho; \
  CCTK_INT const ntot_disk_powerlaw = PRIVATE_DISK_POWERLAW_STRUCT.ntot_disk_powerlaw; \
  CCTK_INT const surfphinum = PRIVATE_DISK_POWERLAW_STRUCT.surfphinum; \
  CCTK_INT const surfxnum = PRIVATE_DISK_POWERLAW_STRUCT.surfxnum; \
  CCTK_INT const surfznum = PRIVATE_DISK_POWERLAW_STRUCT.surfznum; \
  enum { \
      dummy_PRIVATE_DISK_POWERLAW_STRUCT_A_const = sizeof( A_const ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_BigMass = sizeof( BigMass ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_RoM_inner = sizeof( RoM_inner ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_bbh_r_const1 = sizeof( bbh_r_const1 ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_bbh_r_const2 = sizeof( bbh_r_const2 ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_bbh_r_const3 = sizeof( bbh_r_const3 ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_bbh_r_const4 = sizeof( bbh_r_const4 ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_bbh_rhosurf1 = sizeof( bbh_rhosurf1 ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_bbh_rhosurf2 = sizeof( bbh_rhosurf2 ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_bbh_rhosurf3 = sizeof( bbh_rhosurf3 ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_bbh_rhosurf4 = sizeof( bbh_rhosurf4 ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_betam1 = sizeof( betam1 ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_disk_powerlaw_P_max = sizeof( disk_powerlaw_P_max ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_disk_powerlaw_rho_star_max = sizeof( disk_powerlaw_rho_star_max ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_disk_powerlaw_tau_max = sizeof( disk_powerlaw_tau_max ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_ell_inner = sizeof( ell_inner ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_p_c = sizeof( p_c ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_q = sizeof( q ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_rho_cut = sizeof( rho_cut ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_rho_cut_o_rho_max = sizeof( rho_cut_o_rho_max ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_sam_disk = sizeof( sam_disk ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_surfxmax = sizeof( surfxmax ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_surfxmin = sizeof( surfxmin ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_surfzmax = sizeof( surfzmax ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_arithrhosurf = sizeof( arithrhosurf ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_em_field_type = sizeof( em_field_type ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_metric_type = sizeof( metric_type ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_nsurfrho = sizeof( nsurfrho ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_ntot_disk_powerlaw = sizeof( ntot_disk_powerlaw ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_surfphinum = sizeof( surfphinum ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_surfxnum = sizeof( surfxnum ) \
    , dummy_PRIVATE_DISK_POWERLAW_STRUCT_surfznum = sizeof( surfznum ) \
  }; \

