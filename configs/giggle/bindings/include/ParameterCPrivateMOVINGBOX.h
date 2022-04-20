#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_INT BOXES_IN_FINAL_CONFIGURATION;
  CCTK_INT avoid_large_drift_in_ns_box;
  CCTK_INT moving_box_type;
} PRIVATE_MOVINGBOX_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_MOVINGBOX_STRUCT_PARAMS \
  CCTK_INT const BOXES_IN_FINAL_CONFIGURATION = PRIVATE_MOVINGBOX_STRUCT.BOXES_IN_FINAL_CONFIGURATION; \
  CCTK_INT const avoid_large_drift_in_ns_box = PRIVATE_MOVINGBOX_STRUCT.avoid_large_drift_in_ns_box; \
  CCTK_INT const moving_box_type = PRIVATE_MOVINGBOX_STRUCT.moving_box_type; \
  enum { \
      dummy_PRIVATE_MOVINGBOX_STRUCT_BOXES_IN_FINAL_CONFIGURATION = sizeof( BOXES_IN_FINAL_CONFIGURATION ) \
    , dummy_PRIVATE_MOVINGBOX_STRUCT_avoid_large_drift_in_ns_box = sizeof( avoid_large_drift_in_ns_box ) \
    , dummy_PRIVATE_MOVINGBOX_STRUCT_moving_box_type = sizeof( moving_box_type ) \
  }; \

