#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL radius_om[101];
  CCTK_INT num_interp_radii_om;
  CCTK_INT num_theta;
} PRIVATE_ANG_FREQ_OFBFIELD_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_ANG_FREQ_OFBFIELD_STRUCT_PARAMS \
  CCTK_REAL const * const radius_om = PRIVATE_ANG_FREQ_OFBFIELD_STRUCT.radius_om; \
  CCTK_INT const num_interp_radii_om = PRIVATE_ANG_FREQ_OFBFIELD_STRUCT.num_interp_radii_om; \
  CCTK_INT const num_theta = PRIVATE_ANG_FREQ_OFBFIELD_STRUCT.num_theta; \
  enum { \
      dummy_PRIVATE_ANG_FREQ_OFBFIELD_STRUCT_radius_om = sizeof( radius_om ) \
    , dummy_PRIVATE_ANG_FREQ_OFBFIELD_STRUCT_num_interp_radii_om = sizeof( num_interp_radii_om ) \
    , dummy_PRIVATE_ANG_FREQ_OFBFIELD_STRUCT_num_theta = sizeof( num_theta ) \
  }; \

