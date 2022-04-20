#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL radius_star1;
  CCTK_REAL radius_star2;
  CCTK_INT num_matter_refinement_box_centers;
  CCTK_INT track_bhns;
} RESTRICTED_MOVINGBOX_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_MOVINGBOX_STRUCT_PARAMS \
  CCTK_REAL const radius_star1 = RESTRICTED_MOVINGBOX_STRUCT.radius_star1; \
  CCTK_REAL const radius_star2 = RESTRICTED_MOVINGBOX_STRUCT.radius_star2; \
  CCTK_INT const num_matter_refinement_box_centers = RESTRICTED_MOVINGBOX_STRUCT.num_matter_refinement_box_centers; \
  CCTK_INT const track_bhns = RESTRICTED_MOVINGBOX_STRUCT.track_bhns; \
  enum { \
      dummy_RESTRICTED_MOVINGBOX_STRUCT_radius_star1 = sizeof( radius_star1 ) \
    , dummy_RESTRICTED_MOVINGBOX_STRUCT_radius_star2 = sizeof( radius_star2 ) \
    , dummy_RESTRICTED_MOVINGBOX_STRUCT_num_matter_refinement_box_centers = sizeof( num_matter_refinement_box_centers ) \
    , dummy_RESTRICTED_MOVINGBOX_STRUCT_track_bhns = sizeof( track_bhns ) \
  }; \

