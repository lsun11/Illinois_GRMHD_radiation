#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL compute_phi2_max_radius;
  CCTK_REAL compute_phi2_min_radius;
  CCTK_REAL radius_em_phi2[101];
  CCTK_INT num_extraction_radii_em;
  CCTK_INT radius_power_em;
  CCTK_INT scale_with_radius_em;
  CCTK_INT set_up_spherical_EM_wave_in_flat_spacetime;
} PRIVATE_EM_EXTRACTION_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_EM_EXTRACTION_STRUCT_PARAMS \
  CCTK_REAL const compute_phi2_max_radius = PRIVATE_EM_EXTRACTION_STRUCT.compute_phi2_max_radius; \
  CCTK_REAL const compute_phi2_min_radius = PRIVATE_EM_EXTRACTION_STRUCT.compute_phi2_min_radius; \
  CCTK_REAL const * const radius_em_phi2 = PRIVATE_EM_EXTRACTION_STRUCT.radius_em_phi2; \
  CCTK_INT const num_extraction_radii_em = PRIVATE_EM_EXTRACTION_STRUCT.num_extraction_radii_em; \
  CCTK_INT const radius_power_em = PRIVATE_EM_EXTRACTION_STRUCT.radius_power_em; \
  CCTK_INT const scale_with_radius_em = PRIVATE_EM_EXTRACTION_STRUCT.scale_with_radius_em; \
  CCTK_INT const set_up_spherical_EM_wave_in_flat_spacetime = PRIVATE_EM_EXTRACTION_STRUCT.set_up_spherical_EM_wave_in_flat_spacetime; \
  enum { \
      dummy_PRIVATE_EM_EXTRACTION_STRUCT_compute_phi2_max_radius = sizeof( compute_phi2_max_radius ) \
    , dummy_PRIVATE_EM_EXTRACTION_STRUCT_compute_phi2_min_radius = sizeof( compute_phi2_min_radius ) \
    , dummy_PRIVATE_EM_EXTRACTION_STRUCT_radius_em_phi2 = sizeof( radius_em_phi2 ) \
    , dummy_PRIVATE_EM_EXTRACTION_STRUCT_num_extraction_radii_em = sizeof( num_extraction_radii_em ) \
    , dummy_PRIVATE_EM_EXTRACTION_STRUCT_radius_power_em = sizeof( radius_power_em ) \
    , dummy_PRIVATE_EM_EXTRACTION_STRUCT_scale_with_radius_em = sizeof( scale_with_radius_em ) \
    , dummy_PRIVATE_EM_EXTRACTION_STRUCT_set_up_spherical_EM_wave_in_flat_spacetime = sizeof( set_up_spherical_EM_wave_in_flat_spacetime ) \
  }; \

