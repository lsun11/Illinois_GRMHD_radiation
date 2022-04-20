#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL KO_hydro_inside_BH;
  CCTK_REAL KO_inside_BH;
  CCTK_REAL KO_refbd;
  CCTK_REAL damp_lorenz;
  CCTK_INT number_of_points_apply_KO_outside_refboundaries;
  CCTK_INT really_want_gamma_th_of_1_point_3;
  CCTK_INT use_disk_bcs;
  CCTK_INT use_std_PPM_instead_of_efficient_3_gz_version;
} PRIVATE_MHD_EVOLVE_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_MHD_EVOLVE_STRUCT_PARAMS \
  CCTK_REAL const KO_hydro_inside_BH = PRIVATE_MHD_EVOLVE_STRUCT.KO_hydro_inside_BH; \
  CCTK_REAL const KO_inside_BH = PRIVATE_MHD_EVOLVE_STRUCT.KO_inside_BH; \
  CCTK_REAL const KO_refbd = PRIVATE_MHD_EVOLVE_STRUCT.KO_refbd; \
  CCTK_REAL const damp_lorenz = PRIVATE_MHD_EVOLVE_STRUCT.damp_lorenz; \
  CCTK_INT const number_of_points_apply_KO_outside_refboundaries = PRIVATE_MHD_EVOLVE_STRUCT.number_of_points_apply_KO_outside_refboundaries; \
  CCTK_INT const really_want_gamma_th_of_1_point_3 = PRIVATE_MHD_EVOLVE_STRUCT.really_want_gamma_th_of_1_point_3; \
  CCTK_INT const use_disk_bcs = PRIVATE_MHD_EVOLVE_STRUCT.use_disk_bcs; \
  CCTK_INT const use_std_PPM_instead_of_efficient_3_gz_version = PRIVATE_MHD_EVOLVE_STRUCT.use_std_PPM_instead_of_efficient_3_gz_version; \
  enum { \
      dummy_PRIVATE_MHD_EVOLVE_STRUCT_KO_hydro_inside_BH = sizeof( KO_hydro_inside_BH ) \
    , dummy_PRIVATE_MHD_EVOLVE_STRUCT_KO_inside_BH = sizeof( KO_inside_BH ) \
    , dummy_PRIVATE_MHD_EVOLVE_STRUCT_KO_refbd = sizeof( KO_refbd ) \
    , dummy_PRIVATE_MHD_EVOLVE_STRUCT_damp_lorenz = sizeof( damp_lorenz ) \
    , dummy_PRIVATE_MHD_EVOLVE_STRUCT_number_of_points_apply_KO_outside_refboundaries = sizeof( number_of_points_apply_KO_outside_refboundaries ) \
    , dummy_PRIVATE_MHD_EVOLVE_STRUCT_really_want_gamma_th_of_1_point_3 = sizeof( really_want_gamma_th_of_1_point_3 ) \
    , dummy_PRIVATE_MHD_EVOLVE_STRUCT_use_disk_bcs = sizeof( use_disk_bcs ) \
    , dummy_PRIVATE_MHD_EVOLVE_STRUCT_use_std_PPM_instead_of_efficient_3_gz_version = sizeof( use_std_PPM_instead_of_efficient_3_gz_version ) \
  }; \

