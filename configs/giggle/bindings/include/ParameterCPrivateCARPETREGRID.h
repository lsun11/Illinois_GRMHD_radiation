#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL l1xmax;
  CCTK_REAL l1xmin;
  CCTK_REAL l1ymax;
  CCTK_REAL l1ymin;
  CCTK_REAL l1zmax;
  CCTK_REAL l1zmin;
  CCTK_REAL l2xmax;
  CCTK_REAL l2xmin;
  CCTK_REAL l2ymax;
  CCTK_REAL l2ymin;
  CCTK_REAL l2zmax;
  CCTK_REAL l2zmin;
  CCTK_REAL l3xmax;
  CCTK_REAL l3xmin;
  CCTK_REAL l3ymax;
  CCTK_REAL l3ymin;
  CCTK_REAL l3zmax;
  CCTK_REAL l3zmin;
  CCTK_REAL maxerror;
  CCTK_REAL minfraction;
  CCTK_REAL moving_centre_x;
  CCTK_REAL moving_centre_y;
  CCTK_REAL moving_centre_z;
  CCTK_REAL moving_circle_frequency;
  CCTK_REAL moving_circle_radius;
  CCTK_REAL moving_region_radius;
  const char * activate_levels_on_regrid;
  const char * coordinates;
  const char * errorvar;
  const char * gridpoints;
  const char * moving_trajectory;
  const char * outerbounds;
  const char * refined_regions;
  CCTK_INT activate_next;
  CCTK_INT keep_same_grid_structure;
  CCTK_INT l1ixmax;
  CCTK_INT l1ixmin;
  CCTK_INT l1iymax;
  CCTK_INT l1iymin;
  CCTK_INT l1izmax;
  CCTK_INT l1izmin;
  CCTK_INT l2ixmax;
  CCTK_INT l2ixmin;
  CCTK_INT l2iymax;
  CCTK_INT l2iymin;
  CCTK_INT l2izmax;
  CCTK_INT l2izmin;
  CCTK_INT l3ixmax;
  CCTK_INT l3ixmin;
  CCTK_INT l3iymax;
  CCTK_INT l3iymin;
  CCTK_INT l3izmax;
  CCTK_INT l3izmin;
  CCTK_INT merge_overlapping_components;
  CCTK_INT minwidth;
  CCTK_INT num_new_levels;
  CCTK_INT refinement_levels;
  CCTK_INT regrid_every;
  CCTK_INT smart_outer_boundaries;
  CCTK_INT symmetry_x;
  CCTK_INT symmetry_y;
  CCTK_INT symmetry_z;
  CCTK_INT tracking;
  CCTK_INT verbose;
  CCTK_INT veryverbose;
} PRIVATE_CARPETREGRID_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_CARPETREGRID_STRUCT_PARAMS \
  CCTK_REAL const l1xmax = PRIVATE_CARPETREGRID_STRUCT.l1xmax; \
  CCTK_REAL const l1xmin = PRIVATE_CARPETREGRID_STRUCT.l1xmin; \
  CCTK_REAL const l1ymax = PRIVATE_CARPETREGRID_STRUCT.l1ymax; \
  CCTK_REAL const l1ymin = PRIVATE_CARPETREGRID_STRUCT.l1ymin; \
  CCTK_REAL const l1zmax = PRIVATE_CARPETREGRID_STRUCT.l1zmax; \
  CCTK_REAL const l1zmin = PRIVATE_CARPETREGRID_STRUCT.l1zmin; \
  CCTK_REAL const l2xmax = PRIVATE_CARPETREGRID_STRUCT.l2xmax; \
  CCTK_REAL const l2xmin = PRIVATE_CARPETREGRID_STRUCT.l2xmin; \
  CCTK_REAL const l2ymax = PRIVATE_CARPETREGRID_STRUCT.l2ymax; \
  CCTK_REAL const l2ymin = PRIVATE_CARPETREGRID_STRUCT.l2ymin; \
  CCTK_REAL const l2zmax = PRIVATE_CARPETREGRID_STRUCT.l2zmax; \
  CCTK_REAL const l2zmin = PRIVATE_CARPETREGRID_STRUCT.l2zmin; \
  CCTK_REAL const l3xmax = PRIVATE_CARPETREGRID_STRUCT.l3xmax; \
  CCTK_REAL const l3xmin = PRIVATE_CARPETREGRID_STRUCT.l3xmin; \
  CCTK_REAL const l3ymax = PRIVATE_CARPETREGRID_STRUCT.l3ymax; \
  CCTK_REAL const l3ymin = PRIVATE_CARPETREGRID_STRUCT.l3ymin; \
  CCTK_REAL const l3zmax = PRIVATE_CARPETREGRID_STRUCT.l3zmax; \
  CCTK_REAL const l3zmin = PRIVATE_CARPETREGRID_STRUCT.l3zmin; \
  CCTK_REAL const maxerror = PRIVATE_CARPETREGRID_STRUCT.maxerror; \
  CCTK_REAL const minfraction = PRIVATE_CARPETREGRID_STRUCT.minfraction; \
  CCTK_REAL const moving_centre_x = PRIVATE_CARPETREGRID_STRUCT.moving_centre_x; \
  CCTK_REAL const moving_centre_y = PRIVATE_CARPETREGRID_STRUCT.moving_centre_y; \
  CCTK_REAL const moving_centre_z = PRIVATE_CARPETREGRID_STRUCT.moving_centre_z; \
  CCTK_REAL const moving_circle_frequency = PRIVATE_CARPETREGRID_STRUCT.moving_circle_frequency; \
  CCTK_REAL const moving_circle_radius = PRIVATE_CARPETREGRID_STRUCT.moving_circle_radius; \
  CCTK_REAL const moving_region_radius = PRIVATE_CARPETREGRID_STRUCT.moving_region_radius; \
  const char * const activate_levels_on_regrid = PRIVATE_CARPETREGRID_STRUCT.activate_levels_on_regrid; \
  const char * const coordinates = PRIVATE_CARPETREGRID_STRUCT.coordinates; \
  const char * const errorvar = PRIVATE_CARPETREGRID_STRUCT.errorvar; \
  const char * const gridpoints = PRIVATE_CARPETREGRID_STRUCT.gridpoints; \
  const char * const moving_trajectory = PRIVATE_CARPETREGRID_STRUCT.moving_trajectory; \
  const char * const outerbounds = PRIVATE_CARPETREGRID_STRUCT.outerbounds; \
  const char * const refined_regions = PRIVATE_CARPETREGRID_STRUCT.refined_regions; \
  CCTK_INT const activate_next = PRIVATE_CARPETREGRID_STRUCT.activate_next; \
  CCTK_INT const keep_same_grid_structure = PRIVATE_CARPETREGRID_STRUCT.keep_same_grid_structure; \
  CCTK_INT const l1ixmax = PRIVATE_CARPETREGRID_STRUCT.l1ixmax; \
  CCTK_INT const l1ixmin = PRIVATE_CARPETREGRID_STRUCT.l1ixmin; \
  CCTK_INT const l1iymax = PRIVATE_CARPETREGRID_STRUCT.l1iymax; \
  CCTK_INT const l1iymin = PRIVATE_CARPETREGRID_STRUCT.l1iymin; \
  CCTK_INT const l1izmax = PRIVATE_CARPETREGRID_STRUCT.l1izmax; \
  CCTK_INT const l1izmin = PRIVATE_CARPETREGRID_STRUCT.l1izmin; \
  CCTK_INT const l2ixmax = PRIVATE_CARPETREGRID_STRUCT.l2ixmax; \
  CCTK_INT const l2ixmin = PRIVATE_CARPETREGRID_STRUCT.l2ixmin; \
  CCTK_INT const l2iymax = PRIVATE_CARPETREGRID_STRUCT.l2iymax; \
  CCTK_INT const l2iymin = PRIVATE_CARPETREGRID_STRUCT.l2iymin; \
  CCTK_INT const l2izmax = PRIVATE_CARPETREGRID_STRUCT.l2izmax; \
  CCTK_INT const l2izmin = PRIVATE_CARPETREGRID_STRUCT.l2izmin; \
  CCTK_INT const l3ixmax = PRIVATE_CARPETREGRID_STRUCT.l3ixmax; \
  CCTK_INT const l3ixmin = PRIVATE_CARPETREGRID_STRUCT.l3ixmin; \
  CCTK_INT const l3iymax = PRIVATE_CARPETREGRID_STRUCT.l3iymax; \
  CCTK_INT const l3iymin = PRIVATE_CARPETREGRID_STRUCT.l3iymin; \
  CCTK_INT const l3izmax = PRIVATE_CARPETREGRID_STRUCT.l3izmax; \
  CCTK_INT const l3izmin = PRIVATE_CARPETREGRID_STRUCT.l3izmin; \
  CCTK_INT const merge_overlapping_components = PRIVATE_CARPETREGRID_STRUCT.merge_overlapping_components; \
  CCTK_INT const minwidth = PRIVATE_CARPETREGRID_STRUCT.minwidth; \
  CCTK_INT const num_new_levels = PRIVATE_CARPETREGRID_STRUCT.num_new_levels; \
  CCTK_INT const refinement_levels = PRIVATE_CARPETREGRID_STRUCT.refinement_levels; \
  CCTK_INT const regrid_every = PRIVATE_CARPETREGRID_STRUCT.regrid_every; \
  CCTK_INT const smart_outer_boundaries = PRIVATE_CARPETREGRID_STRUCT.smart_outer_boundaries; \
  CCTK_INT const symmetry_x = PRIVATE_CARPETREGRID_STRUCT.symmetry_x; \
  CCTK_INT const symmetry_y = PRIVATE_CARPETREGRID_STRUCT.symmetry_y; \
  CCTK_INT const symmetry_z = PRIVATE_CARPETREGRID_STRUCT.symmetry_z; \
  CCTK_INT const tracking = PRIVATE_CARPETREGRID_STRUCT.tracking; \
  CCTK_INT const verbose = PRIVATE_CARPETREGRID_STRUCT.verbose; \
  CCTK_INT const veryverbose = PRIVATE_CARPETREGRID_STRUCT.veryverbose; \
  enum { \
      dummy_PRIVATE_CARPETREGRID_STRUCT_l1xmax = sizeof( l1xmax ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l1xmin = sizeof( l1xmin ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l1ymax = sizeof( l1ymax ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l1ymin = sizeof( l1ymin ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l1zmax = sizeof( l1zmax ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l1zmin = sizeof( l1zmin ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l2xmax = sizeof( l2xmax ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l2xmin = sizeof( l2xmin ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l2ymax = sizeof( l2ymax ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l2ymin = sizeof( l2ymin ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l2zmax = sizeof( l2zmax ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l2zmin = sizeof( l2zmin ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l3xmax = sizeof( l3xmax ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l3xmin = sizeof( l3xmin ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l3ymax = sizeof( l3ymax ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l3ymin = sizeof( l3ymin ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l3zmax = sizeof( l3zmax ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l3zmin = sizeof( l3zmin ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_maxerror = sizeof( maxerror ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_minfraction = sizeof( minfraction ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_moving_centre_x = sizeof( moving_centre_x ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_moving_centre_y = sizeof( moving_centre_y ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_moving_centre_z = sizeof( moving_centre_z ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_moving_circle_frequency = sizeof( moving_circle_frequency ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_moving_circle_radius = sizeof( moving_circle_radius ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_moving_region_radius = sizeof( moving_region_radius ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_activate_levels_on_regrid = sizeof( activate_levels_on_regrid ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_coordinates = sizeof( coordinates ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_errorvar = sizeof( errorvar ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_gridpoints = sizeof( gridpoints ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_moving_trajectory = sizeof( moving_trajectory ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_outerbounds = sizeof( outerbounds ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_refined_regions = sizeof( refined_regions ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_activate_next = sizeof( activate_next ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_keep_same_grid_structure = sizeof( keep_same_grid_structure ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l1ixmax = sizeof( l1ixmax ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l1ixmin = sizeof( l1ixmin ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l1iymax = sizeof( l1iymax ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l1iymin = sizeof( l1iymin ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l1izmax = sizeof( l1izmax ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l1izmin = sizeof( l1izmin ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l2ixmax = sizeof( l2ixmax ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l2ixmin = sizeof( l2ixmin ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l2iymax = sizeof( l2iymax ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l2iymin = sizeof( l2iymin ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l2izmax = sizeof( l2izmax ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l2izmin = sizeof( l2izmin ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l3ixmax = sizeof( l3ixmax ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l3ixmin = sizeof( l3ixmin ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l3iymax = sizeof( l3iymax ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l3iymin = sizeof( l3iymin ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l3izmax = sizeof( l3izmax ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_l3izmin = sizeof( l3izmin ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_merge_overlapping_components = sizeof( merge_overlapping_components ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_minwidth = sizeof( minwidth ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_num_new_levels = sizeof( num_new_levels ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_refinement_levels = sizeof( refinement_levels ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_regrid_every = sizeof( regrid_every ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_smart_outer_boundaries = sizeof( smart_outer_boundaries ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_symmetry_x = sizeof( symmetry_x ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_symmetry_y = sizeof( symmetry_y ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_symmetry_z = sizeof( symmetry_z ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_tracking = sizeof( tracking ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_verbose = sizeof( verbose ) \
    , dummy_PRIVATE_CARPETREGRID_STRUCT_veryverbose = sizeof( veryverbose ) \
  }; \

