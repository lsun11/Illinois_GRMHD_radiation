#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL KO_strength_falloff_dr;
  CCTK_REAL KO_strength_falloff_radius;
  CCTK_REAL KO_strength_final_value;
  CCTK_REAL ah_max_epsdis;
  CCTK_REAL ah_radius_offset;
  CCTK_REAL ah_slope;
  CCTK_REAL epsdis;
  CCTK_REAL epsdis_for_level[32];
  CCTK_REAL ob_slope;
  CCTK_REAL outer_boundary_max_epsdis;
  const char * vars;
  CCTK_INT KO_strength_falloff_enable;
  CCTK_INT extra_dissipation_at_outerbound;
  CCTK_INT extra_dissipation_in_horizons;
  CCTK_INT horizon_number[100];
  CCTK_INT order;
  CCTK_INT outer_bound_npoints;
  CCTK_INT respect_emask;
  CCTK_INT surface_number[100];
  CCTK_INT update_ah_every;
  CCTK_INT use_dissipation_near_excision;
  CCTK_INT verbose;
} RESTRICTED_DISSIPATION_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_DISSIPATION_STRUCT_PARAMS \
  CCTK_REAL const KO_strength_falloff_dr = RESTRICTED_DISSIPATION_STRUCT.KO_strength_falloff_dr; \
  CCTK_REAL const KO_strength_falloff_radius = RESTRICTED_DISSIPATION_STRUCT.KO_strength_falloff_radius; \
  CCTK_REAL const KO_strength_final_value = RESTRICTED_DISSIPATION_STRUCT.KO_strength_final_value; \
  CCTK_REAL const ah_max_epsdis = RESTRICTED_DISSIPATION_STRUCT.ah_max_epsdis; \
  CCTK_REAL const ah_radius_offset = RESTRICTED_DISSIPATION_STRUCT.ah_radius_offset; \
  CCTK_REAL const ah_slope = RESTRICTED_DISSIPATION_STRUCT.ah_slope; \
  CCTK_REAL const epsdis = RESTRICTED_DISSIPATION_STRUCT.epsdis; \
  CCTK_REAL const * const epsdis_for_level = RESTRICTED_DISSIPATION_STRUCT.epsdis_for_level; \
  CCTK_REAL const ob_slope = RESTRICTED_DISSIPATION_STRUCT.ob_slope; \
  CCTK_REAL const outer_boundary_max_epsdis = RESTRICTED_DISSIPATION_STRUCT.outer_boundary_max_epsdis; \
  const char * const vars = RESTRICTED_DISSIPATION_STRUCT.vars; \
  CCTK_INT const KO_strength_falloff_enable = RESTRICTED_DISSIPATION_STRUCT.KO_strength_falloff_enable; \
  CCTK_INT const extra_dissipation_at_outerbound = RESTRICTED_DISSIPATION_STRUCT.extra_dissipation_at_outerbound; \
  CCTK_INT const extra_dissipation_in_horizons = RESTRICTED_DISSIPATION_STRUCT.extra_dissipation_in_horizons; \
  CCTK_INT const * const horizon_number = RESTRICTED_DISSIPATION_STRUCT.horizon_number; \
  CCTK_INT const order = RESTRICTED_DISSIPATION_STRUCT.order; \
  CCTK_INT const outer_bound_npoints = RESTRICTED_DISSIPATION_STRUCT.outer_bound_npoints; \
  CCTK_INT const respect_emask = RESTRICTED_DISSIPATION_STRUCT.respect_emask; \
  CCTK_INT const * const surface_number = RESTRICTED_DISSIPATION_STRUCT.surface_number; \
  CCTK_INT const update_ah_every = RESTRICTED_DISSIPATION_STRUCT.update_ah_every; \
  CCTK_INT const use_dissipation_near_excision = RESTRICTED_DISSIPATION_STRUCT.use_dissipation_near_excision; \
  CCTK_INT const verbose = RESTRICTED_DISSIPATION_STRUCT.verbose; \
  enum { \
      dummy_RESTRICTED_DISSIPATION_STRUCT_KO_strength_falloff_dr = sizeof( KO_strength_falloff_dr ) \
    , dummy_RESTRICTED_DISSIPATION_STRUCT_KO_strength_falloff_radius = sizeof( KO_strength_falloff_radius ) \
    , dummy_RESTRICTED_DISSIPATION_STRUCT_KO_strength_final_value = sizeof( KO_strength_final_value ) \
    , dummy_RESTRICTED_DISSIPATION_STRUCT_ah_max_epsdis = sizeof( ah_max_epsdis ) \
    , dummy_RESTRICTED_DISSIPATION_STRUCT_ah_radius_offset = sizeof( ah_radius_offset ) \
    , dummy_RESTRICTED_DISSIPATION_STRUCT_ah_slope = sizeof( ah_slope ) \
    , dummy_RESTRICTED_DISSIPATION_STRUCT_epsdis = sizeof( epsdis ) \
    , dummy_RESTRICTED_DISSIPATION_STRUCT_epsdis_for_level = sizeof( epsdis_for_level ) \
    , dummy_RESTRICTED_DISSIPATION_STRUCT_ob_slope = sizeof( ob_slope ) \
    , dummy_RESTRICTED_DISSIPATION_STRUCT_outer_boundary_max_epsdis = sizeof( outer_boundary_max_epsdis ) \
    , dummy_RESTRICTED_DISSIPATION_STRUCT_vars = sizeof( vars ) \
    , dummy_RESTRICTED_DISSIPATION_STRUCT_KO_strength_falloff_enable = sizeof( KO_strength_falloff_enable ) \
    , dummy_RESTRICTED_DISSIPATION_STRUCT_extra_dissipation_at_outerbound = sizeof( extra_dissipation_at_outerbound ) \
    , dummy_RESTRICTED_DISSIPATION_STRUCT_extra_dissipation_in_horizons = sizeof( extra_dissipation_in_horizons ) \
    , dummy_RESTRICTED_DISSIPATION_STRUCT_horizon_number = sizeof( horizon_number ) \
    , dummy_RESTRICTED_DISSIPATION_STRUCT_order = sizeof( order ) \
    , dummy_RESTRICTED_DISSIPATION_STRUCT_outer_bound_npoints = sizeof( outer_bound_npoints ) \
    , dummy_RESTRICTED_DISSIPATION_STRUCT_respect_emask = sizeof( respect_emask ) \
    , dummy_RESTRICTED_DISSIPATION_STRUCT_surface_number = sizeof( surface_number ) \
    , dummy_RESTRICTED_DISSIPATION_STRUCT_update_ah_every = sizeof( update_ah_every ) \
    , dummy_RESTRICTED_DISSIPATION_STRUCT_use_dissipation_near_excision = sizeof( use_dissipation_near_excision ) \
    , dummy_RESTRICTED_DISSIPATION_STRUCT_verbose = sizeof( verbose ) \
  }; \

