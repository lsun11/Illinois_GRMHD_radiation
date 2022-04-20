#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL compute_Psi4_max_radius;
  CCTK_REAL compute_Psi4_min_radius;
  CCTK_REAL radius_GW_Psi4[101];
  const char * psif_vec;
  CCTK_INT enable_interp_onepoint;
  CCTK_INT nmodes_Psi4;
  CCTK_INT nmodes_ZM;
  CCTK_INT num_extraction_radii;
  CCTK_INT radius_power;
  CCTK_INT scale_with_radius;
  CCTK_INT use_Rij_from_compute_ricci;
} PRIVATE_GW_EXTRACTION_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_GW_EXTRACTION_STRUCT_PARAMS \
  CCTK_REAL const compute_Psi4_max_radius = PRIVATE_GW_EXTRACTION_STRUCT.compute_Psi4_max_radius; \
  CCTK_REAL const compute_Psi4_min_radius = PRIVATE_GW_EXTRACTION_STRUCT.compute_Psi4_min_radius; \
  CCTK_REAL const * const radius_GW_Psi4 = PRIVATE_GW_EXTRACTION_STRUCT.radius_GW_Psi4; \
  const char * const psif_vec = PRIVATE_GW_EXTRACTION_STRUCT.psif_vec; \
  CCTK_INT const enable_interp_onepoint = PRIVATE_GW_EXTRACTION_STRUCT.enable_interp_onepoint; \
  CCTK_INT const nmodes_Psi4 = PRIVATE_GW_EXTRACTION_STRUCT.nmodes_Psi4; \
  CCTK_INT const nmodes_ZM = PRIVATE_GW_EXTRACTION_STRUCT.nmodes_ZM; \
  CCTK_INT const num_extraction_radii = PRIVATE_GW_EXTRACTION_STRUCT.num_extraction_radii; \
  CCTK_INT const radius_power = PRIVATE_GW_EXTRACTION_STRUCT.radius_power; \
  CCTK_INT const scale_with_radius = PRIVATE_GW_EXTRACTION_STRUCT.scale_with_radius; \
  CCTK_INT const use_Rij_from_compute_ricci = PRIVATE_GW_EXTRACTION_STRUCT.use_Rij_from_compute_ricci; \
  enum { \
      dummy_PRIVATE_GW_EXTRACTION_STRUCT_compute_Psi4_max_radius = sizeof( compute_Psi4_max_radius ) \
    , dummy_PRIVATE_GW_EXTRACTION_STRUCT_compute_Psi4_min_radius = sizeof( compute_Psi4_min_radius ) \
    , dummy_PRIVATE_GW_EXTRACTION_STRUCT_radius_GW_Psi4 = sizeof( radius_GW_Psi4 ) \
    , dummy_PRIVATE_GW_EXTRACTION_STRUCT_psif_vec = sizeof( psif_vec ) \
    , dummy_PRIVATE_GW_EXTRACTION_STRUCT_enable_interp_onepoint = sizeof( enable_interp_onepoint ) \
    , dummy_PRIVATE_GW_EXTRACTION_STRUCT_nmodes_Psi4 = sizeof( nmodes_Psi4 ) \
    , dummy_PRIVATE_GW_EXTRACTION_STRUCT_nmodes_ZM = sizeof( nmodes_ZM ) \
    , dummy_PRIVATE_GW_EXTRACTION_STRUCT_num_extraction_radii = sizeof( num_extraction_radii ) \
    , dummy_PRIVATE_GW_EXTRACTION_STRUCT_radius_power = sizeof( radius_power ) \
    , dummy_PRIVATE_GW_EXTRACTION_STRUCT_scale_with_radius = sizeof( scale_with_radius ) \
    , dummy_PRIVATE_GW_EXTRACTION_STRUCT_use_Rij_from_compute_ricci = sizeof( use_Rij_from_compute_ricci ) \
  }; \

