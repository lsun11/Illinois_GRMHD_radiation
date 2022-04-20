#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL excision_radius;
  CCTK_REAL run_time;
  CCTK_INT Symmetry;
  CCTK_INT bssn_enable;
  CCTK_INT cowling_enable;
  CCTK_INT excision_enable;
  CCTK_INT fisheye_enable;
  CCTK_INT iter_count;
  CCTK_INT number_of_mol_ministeps;
  CCTK_INT rot_metric;
  CCTK_INT trA_detg_enforce;
} GLOBAL_PARAMETER_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_GLOBAL_PARAMETER_STRUCT_PARAMS \
  CCTK_REAL const excision_radius = GLOBAL_PARAMETER_STRUCT.excision_radius; \
  CCTK_REAL const run_time = GLOBAL_PARAMETER_STRUCT.run_time; \
  CCTK_INT const Symmetry = GLOBAL_PARAMETER_STRUCT.Symmetry; \
  CCTK_INT const bssn_enable = GLOBAL_PARAMETER_STRUCT.bssn_enable; \
  CCTK_INT const cowling_enable = GLOBAL_PARAMETER_STRUCT.cowling_enable; \
  CCTK_INT const excision_enable = GLOBAL_PARAMETER_STRUCT.excision_enable; \
  CCTK_INT const fisheye_enable = GLOBAL_PARAMETER_STRUCT.fisheye_enable; \
  CCTK_INT const iter_count = GLOBAL_PARAMETER_STRUCT.iter_count; \
  CCTK_INT const number_of_mol_ministeps = GLOBAL_PARAMETER_STRUCT.number_of_mol_ministeps; \
  CCTK_INT const rot_metric = GLOBAL_PARAMETER_STRUCT.rot_metric; \
  CCTK_INT const trA_detg_enforce = GLOBAL_PARAMETER_STRUCT.trA_detg_enforce; \
  enum { \
      dummy_GLOBAL_PARAMETER_STRUCT_excision_radius = sizeof( excision_radius ) \
    , dummy_GLOBAL_PARAMETER_STRUCT_run_time = sizeof( run_time ) \
    , dummy_GLOBAL_PARAMETER_STRUCT_Symmetry = sizeof( Symmetry ) \
    , dummy_GLOBAL_PARAMETER_STRUCT_bssn_enable = sizeof( bssn_enable ) \
    , dummy_GLOBAL_PARAMETER_STRUCT_cowling_enable = sizeof( cowling_enable ) \
    , dummy_GLOBAL_PARAMETER_STRUCT_excision_enable = sizeof( excision_enable ) \
    , dummy_GLOBAL_PARAMETER_STRUCT_fisheye_enable = sizeof( fisheye_enable ) \
    , dummy_GLOBAL_PARAMETER_STRUCT_iter_count = sizeof( iter_count ) \
    , dummy_GLOBAL_PARAMETER_STRUCT_number_of_mol_ministeps = sizeof( number_of_mol_ministeps ) \
    , dummy_GLOBAL_PARAMETER_STRUCT_rot_metric = sizeof( rot_metric ) \
    , dummy_GLOBAL_PARAMETER_STRUCT_trA_detg_enforce = sizeof( trA_detg_enforce ) \
  }; \

