/*@@
   @header  Dissipation_arguments.h
   @author  Automatically generated by GridFuncStuff.pl
   @desc
            Defines macros to declare/define/pass function arguments
            in calls from C to Fortran for thorn Dissipation
   @enddesc
 @@*/


#ifdef FCODE
#define DECLARE_DISSIPATION_PRIVATE_FARGUMENTS \
INTEGER X0epsdisA_group&&\
INTEGER X1epsdisA_group&&\
INTEGER X2epsdisA_group&&\
CCTK_REAL epsdisA(X0epsdisA_group,X1epsdisA_group,X2epsdisA_group)&&\


#define DISSIPATION_PRIVATE_FARGUMENTS \
X0epsdisA_group,X1epsdisA_group,X2epsdisA_group,epsdisA

#endif /* FCODE */

#ifdef CCODE
#define DECLARE_DISSIPATION_PRIVATE_CARGUMENTS \
CCTK_REAL * CCTK_RESTRICT epsdisA = (cctki_dummy_int = &epsdisA - &epsdisA, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "DISSIPATION", "epsdisA"));

#define DECLARE_DISSIPATION_PRIVATE_C2F \
static int CCTKARGNUM_epsdisA = -1; \
static int CCTKGROUPNUM_epsdisA_group = -1;

#define INITIALISE_DISSIPATION_PRIVATE_C2F \
if(CCTKARGNUM_epsdisA == -1) CCTKARGNUM_epsdisA = CCTK_VarIndex("Dissipation::epsdisA"); \
if(CCTKGROUPNUM_epsdisA_group == -1) CCTKGROUPNUM_epsdisA_group = CCTK_GroupIndex("Dissipation::epsdisA_group");

#define DISSIPATION_PRIVATE_C2F_PROTO \
const int *,const int *,const int *,CCTK_REAL *

#define PASS_DISSIPATION_PRIVATE_C2F(GH) \
PASS_GROUPSIZE(epsdisA_group, 0),\
PASS_GROUPSIZE(epsdisA_group, 1),\
PASS_GROUPSIZE(epsdisA_group, 2),\
(CCTK_REAL *)(PASS_REFERENCE(epsdisA, 0))

#endif /* CCODE */

#ifdef FCODE
#define DECLARE_DISSIPATION_PROTECTED_FARGUMENTS \


#define DISSIPATION_PROTECTED_FARGUMENTS \


#endif /* FCODE */

#ifdef CCODE
#define DECLARE_DISSIPATION_PROTECTED_CARGUMENTS \


#define DECLARE_DISSIPATION_PROTECTED_C2F \


#define INITIALISE_DISSIPATION_PROTECTED_C2F \


#define DISSIPATION_PROTECTED_C2F_PROTO \


#define PASS_DISSIPATION_PROTECTED_C2F(GH) \


#endif /* CCODE */

#ifdef FCODE
#define DECLARE_DISSIPATION_PUBLIC_FARGUMENTS \
INTEGER X0coordinates&&\
INTEGER X0mask&&\
INTEGER X0sf_radius&&\
INTEGER X0space_mask_group&&\
INTEGER X1coordinates&&\
INTEGER X1mask&&\
INTEGER X1sf_radius&&\
INTEGER X1space_mask_group&&\
INTEGER X2coordinates&&\
INTEGER X2mask&&\
INTEGER X2space_mask_group&&\
INTEGER sf_active_length&&\
INTEGER sf_coordinate_descriptors_length&&\
INTEGER sf_info_length&&\
INTEGER sf_maxreflevel_length&&\
INTEGER sf_minreflevel_length&&\
INTEGER sf_origin_length&&\
INTEGER sf_radius_length&&\
INTEGER sf_shape_descriptors_length&&\
INTEGER sf_valid_length&&\
CCTK_REAL coarse_dx&&\
CCTK_REAL coarse_dy&&\
CCTK_REAL coarse_dz&&\
CCTK_REAL emask(X0mask,X1mask,X2mask)&&\
CCTK_REAL r(X0coordinates,X1coordinates,X2coordinates)&&\
CCTK_INT sf_active(sf_active_length)&&\
CCTK_REAL sf_area(sf_info_length)&&\
CCTK_REAL sf_centroid_x(sf_info_length)&&\
CCTK_REAL sf_centroid_y(sf_info_length)&&\
CCTK_REAL sf_centroid_z(sf_info_length)&&\
CCTK_REAL sf_delta_phi(sf_coordinate_descriptors_length)&&\
CCTK_REAL sf_delta_theta(sf_coordinate_descriptors_length)&&\
CCTK_REAL sf_max_radius(sf_info_length)&&\
CCTK_REAL sf_max_x(sf_info_length)&&\
CCTK_REAL sf_max_y(sf_info_length)&&\
CCTK_REAL sf_max_z(sf_info_length)&&\
CCTK_INT sf_maxreflevel(sf_maxreflevel_length)&&\
CCTK_REAL sf_mean_radius(sf_info_length)&&\
CCTK_REAL sf_min_radius(sf_info_length)&&\
CCTK_REAL sf_min_x(sf_info_length)&&\
CCTK_REAL sf_min_y(sf_info_length)&&\
CCTK_REAL sf_min_z(sf_info_length)&&\
CCTK_INT sf_minreflevel(sf_minreflevel_length)&&\
CCTK_INT sf_nghostsphi(sf_shape_descriptors_length)&&\
CCTK_INT sf_nghoststheta(sf_shape_descriptors_length)&&\
CCTK_INT sf_nphi(sf_shape_descriptors_length)&&\
CCTK_INT sf_ntheta(sf_shape_descriptors_length)&&\
CCTK_REAL sf_origin_phi(sf_coordinate_descriptors_length)&&\
CCTK_REAL sf_origin_theta(sf_coordinate_descriptors_length)&&\
CCTK_REAL sf_origin_x(sf_origin_length)&&\
CCTK_REAL sf_origin_y(sf_origin_length)&&\
CCTK_REAL sf_origin_z(sf_origin_length)&&\
CCTK_REAL sf_quadrupole_xx(sf_info_length)&&\
CCTK_REAL sf_quadrupole_xy(sf_info_length)&&\
CCTK_REAL sf_quadrupole_xz(sf_info_length)&&\
CCTK_REAL sf_quadrupole_yy(sf_info_length)&&\
CCTK_REAL sf_quadrupole_yz(sf_info_length)&&\
CCTK_REAL sf_quadrupole_zz(sf_info_length)&&\
CCTK_REAL sf_radius(X0sf_radius,X1sf_radius,sf_radius_length)&&\
CCTK_INT sf_valid(sf_valid_length)&&\
CCTK_INT space_mask(X0space_mask_group,X1space_mask_group,X2space_mask_group)&&\
CCTK_REAL x(X0coordinates,X1coordinates,X2coordinates)&&\
CCTK_REAL y(X0coordinates,X1coordinates,X2coordinates)&&\
CCTK_REAL z(X0coordinates,X1coordinates,X2coordinates)&&\


#define DISSIPATION_PUBLIC_FARGUMENTS \
X0coordinates,X0mask,X0sf_radius,X0space_mask_group,X1coordinates,X1mask,X1sf_radius,X1space_mask_group,X2coordinates,X2mask,X2space_mask_group,sf_active_length,sf_coordinate_descriptors_length,sf_info_length,sf_maxreflevel_length,sf_minreflevel_length,sf_origin_length,sf_radius_length,sf_shape_descriptors_length,sf_valid_length,coarse_dx,coarse_dy,coarse_dz,emask,r,sf_active,sf_area,sf_centroid_x,sf_centroid_y,sf_centroid_z,sf_delta_phi,sf_delta_theta,sf_max_radius,sf_max_x,sf_max_y,sf_max_z,sf_maxreflevel,sf_mean_radius,sf_min_radius,sf_min_x,sf_min_y,sf_min_z,sf_minreflevel,sf_nghostsphi,sf_nghoststheta,sf_nphi,sf_ntheta,sf_origin_phi,sf_origin_theta,sf_origin_x,sf_origin_y,sf_origin_z,sf_quadrupole_xx,sf_quadrupole_xy,sf_quadrupole_xz,sf_quadrupole_yy,sf_quadrupole_yz,sf_quadrupole_zz,sf_radius,sf_valid,space_mask,x,y,z

#endif /* FCODE */

#ifdef CCODE
#define DECLARE_DISSIPATION_PUBLIC_CARGUMENTS \
CCTK_REAL * CCTK_RESTRICT coarse_dx = (cctki_dummy_int = &coarse_dx - &coarse_dx, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "GRID", "coarse_dx")); \
CCTK_REAL * CCTK_RESTRICT coarse_dy = (cctki_dummy_int = &coarse_dy - &coarse_dy, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "GRID", "coarse_dy")); \
CCTK_REAL * CCTK_RESTRICT coarse_dz = (cctki_dummy_int = &coarse_dz - &coarse_dz, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "GRID", "coarse_dz")); \
CCTK_REAL * CCTK_RESTRICT emask = (cctki_dummy_int = &emask - &emask, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPACEMASK", "emask")); \
CCTK_REAL * CCTK_RESTRICT r = (cctki_dummy_int = &r - &r, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "GRID", "r")); \
CCTK_INT * CCTK_RESTRICT sf_active = (cctki_dummy_int = &sf_active - &sf_active, (CCTK_INT *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_active[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_area = (cctki_dummy_int = &sf_area - &sf_area, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_area[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_centroid_x = (cctki_dummy_int = &sf_centroid_x - &sf_centroid_x, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_centroid_x[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_centroid_y = (cctki_dummy_int = &sf_centroid_y - &sf_centroid_y, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_centroid_y[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_centroid_z = (cctki_dummy_int = &sf_centroid_z - &sf_centroid_z, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_centroid_z[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_delta_phi = (cctki_dummy_int = &sf_delta_phi - &sf_delta_phi, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_delta_phi[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_delta_theta = (cctki_dummy_int = &sf_delta_theta - &sf_delta_theta, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_delta_theta[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_max_radius = (cctki_dummy_int = &sf_max_radius - &sf_max_radius, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_max_radius[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_max_x = (cctki_dummy_int = &sf_max_x - &sf_max_x, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_max_x[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_max_y = (cctki_dummy_int = &sf_max_y - &sf_max_y, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_max_y[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_max_z = (cctki_dummy_int = &sf_max_z - &sf_max_z, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_max_z[0]")); \
CCTK_INT * CCTK_RESTRICT sf_maxreflevel = (cctki_dummy_int = &sf_maxreflevel - &sf_maxreflevel, (CCTK_INT *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_maxreflevel[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_mean_radius = (cctki_dummy_int = &sf_mean_radius - &sf_mean_radius, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_mean_radius[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_min_radius = (cctki_dummy_int = &sf_min_radius - &sf_min_radius, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_min_radius[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_min_x = (cctki_dummy_int = &sf_min_x - &sf_min_x, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_min_x[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_min_y = (cctki_dummy_int = &sf_min_y - &sf_min_y, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_min_y[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_min_z = (cctki_dummy_int = &sf_min_z - &sf_min_z, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_min_z[0]")); \
CCTK_INT * CCTK_RESTRICT sf_minreflevel = (cctki_dummy_int = &sf_minreflevel - &sf_minreflevel, (CCTK_INT *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_minreflevel[0]")); \
CCTK_INT * CCTK_RESTRICT sf_nghostsphi = (cctki_dummy_int = &sf_nghostsphi - &sf_nghostsphi, (CCTK_INT *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_nghostsphi[0]")); \
CCTK_INT * CCTK_RESTRICT sf_nghoststheta = (cctki_dummy_int = &sf_nghoststheta - &sf_nghoststheta, (CCTK_INT *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_nghoststheta[0]")); \
CCTK_INT * CCTK_RESTRICT sf_nphi = (cctki_dummy_int = &sf_nphi - &sf_nphi, (CCTK_INT *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_nphi[0]")); \
CCTK_INT * CCTK_RESTRICT sf_ntheta = (cctki_dummy_int = &sf_ntheta - &sf_ntheta, (CCTK_INT *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_ntheta[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_origin_phi = (cctki_dummy_int = &sf_origin_phi - &sf_origin_phi, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_origin_phi[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_origin_theta = (cctki_dummy_int = &sf_origin_theta - &sf_origin_theta, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_origin_theta[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_origin_x = (cctki_dummy_int = &sf_origin_x - &sf_origin_x, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_origin_x[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_origin_y = (cctki_dummy_int = &sf_origin_y - &sf_origin_y, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_origin_y[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_origin_z = (cctki_dummy_int = &sf_origin_z - &sf_origin_z, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_origin_z[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_quadrupole_xx = (cctki_dummy_int = &sf_quadrupole_xx - &sf_quadrupole_xx, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_quadrupole_xx[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_quadrupole_xy = (cctki_dummy_int = &sf_quadrupole_xy - &sf_quadrupole_xy, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_quadrupole_xy[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_quadrupole_xz = (cctki_dummy_int = &sf_quadrupole_xz - &sf_quadrupole_xz, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_quadrupole_xz[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_quadrupole_yy = (cctki_dummy_int = &sf_quadrupole_yy - &sf_quadrupole_yy, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_quadrupole_yy[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_quadrupole_yz = (cctki_dummy_int = &sf_quadrupole_yz - &sf_quadrupole_yz, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_quadrupole_yz[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_quadrupole_zz = (cctki_dummy_int = &sf_quadrupole_zz - &sf_quadrupole_zz, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_quadrupole_zz[0]")); \
CCTK_REAL * CCTK_RESTRICT sf_radius = (cctki_dummy_int = &sf_radius - &sf_radius, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_radius[0]")); \
CCTK_INT * CCTK_RESTRICT sf_valid = (cctki_dummy_int = &sf_valid - &sf_valid, (CCTK_INT *) CCTKi_VarDataPtr(cctkGH, 0, "SPHERICALSURFACE", "sf_valid[0]")); \
CCTK_INT * CCTK_RESTRICT space_mask = (cctki_dummy_int = &space_mask - &space_mask, (CCTK_INT *) CCTKi_VarDataPtr(cctkGH, 0, "SPACEMASK", "space_mask")); \
CCTK_REAL * CCTK_RESTRICT x = (cctki_dummy_int = &x - &x, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "GRID", "x")); \
CCTK_REAL * CCTK_RESTRICT y = (cctki_dummy_int = &y - &y, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "GRID", "y")); \
CCTK_REAL * CCTK_RESTRICT z = (cctki_dummy_int = &z - &z, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "GRID", "z"));

#define DECLARE_DISSIPATION_PUBLIC_C2F \
static int CCTKARGNUM_coarse_dx = -1; \
static int CCTKGROUPNUM_gridspacings = -1; \
static int CCTKARGNUM_coarse_dy = -1; \
static int CCTKARGNUM_coarse_dz = -1; \
static int CCTKARGNUM_emask = -1; \
static int CCTKGROUPNUM_mask = -1; \
static int CCTKARGNUM_r = -1; \
static int CCTKGROUPNUM_coordinates = -1; \
static int CCTKARGNUM_sf_active = -1; \
static int CCTKGROUPNUM_sf_active = -1; \
static int CCTKARGNUM_sf_area = -1; \
static int CCTKGROUPNUM_sf_info = -1; \
static int CCTKARGNUM_sf_centroid_x = -1; \
static int CCTKARGNUM_sf_centroid_y = -1; \
static int CCTKARGNUM_sf_centroid_z = -1; \
static int CCTKARGNUM_sf_delta_phi = -1; \
static int CCTKGROUPNUM_sf_coordinate_descriptors = -1; \
static int CCTKARGNUM_sf_delta_theta = -1; \
static int CCTKARGNUM_sf_max_radius = -1; \
static int CCTKARGNUM_sf_max_x = -1; \
static int CCTKARGNUM_sf_max_y = -1; \
static int CCTKARGNUM_sf_max_z = -1; \
static int CCTKARGNUM_sf_maxreflevel = -1; \
static int CCTKGROUPNUM_sf_maxreflevel = -1; \
static int CCTKARGNUM_sf_mean_radius = -1; \
static int CCTKARGNUM_sf_min_radius = -1; \
static int CCTKARGNUM_sf_min_x = -1; \
static int CCTKARGNUM_sf_min_y = -1; \
static int CCTKARGNUM_sf_min_z = -1; \
static int CCTKARGNUM_sf_minreflevel = -1; \
static int CCTKGROUPNUM_sf_minreflevel = -1; \
static int CCTKARGNUM_sf_nghostsphi = -1; \
static int CCTKGROUPNUM_sf_shape_descriptors = -1; \
static int CCTKARGNUM_sf_nghoststheta = -1; \
static int CCTKARGNUM_sf_nphi = -1; \
static int CCTKARGNUM_sf_ntheta = -1; \
static int CCTKARGNUM_sf_origin_phi = -1; \
static int CCTKARGNUM_sf_origin_theta = -1; \
static int CCTKARGNUM_sf_origin_x = -1; \
static int CCTKGROUPNUM_sf_origin = -1; \
static int CCTKARGNUM_sf_origin_y = -1; \
static int CCTKARGNUM_sf_origin_z = -1; \
static int CCTKARGNUM_sf_quadrupole_xx = -1; \
static int CCTKARGNUM_sf_quadrupole_xy = -1; \
static int CCTKARGNUM_sf_quadrupole_xz = -1; \
static int CCTKARGNUM_sf_quadrupole_yy = -1; \
static int CCTKARGNUM_sf_quadrupole_yz = -1; \
static int CCTKARGNUM_sf_quadrupole_zz = -1; \
static int CCTKARGNUM_sf_radius = -1; \
static int CCTKGROUPNUM_sf_radius = -1; \
static int CCTKARGNUM_sf_valid = -1; \
static int CCTKGROUPNUM_sf_valid = -1; \
static int CCTKARGNUM_space_mask = -1; \
static int CCTKGROUPNUM_space_mask_group = -1; \
static int CCTKARGNUM_x = -1; \
static int CCTKARGNUM_y = -1; \
static int CCTKARGNUM_z = -1;

#define INITIALISE_DISSIPATION_PUBLIC_C2F \
if(CCTKARGNUM_coarse_dx == -1) CCTKARGNUM_coarse_dx = CCTK_VarIndex("GRID::coarse_dx"); \
if(CCTKGROUPNUM_gridspacings == -1) CCTKGROUPNUM_gridspacings = CCTK_GroupIndex("GRID::gridspacings"); \
if(CCTKARGNUM_coarse_dy == -1) CCTKARGNUM_coarse_dy = CCTK_VarIndex("GRID::coarse_dy"); \
if(CCTKARGNUM_coarse_dz == -1) CCTKARGNUM_coarse_dz = CCTK_VarIndex("GRID::coarse_dz"); \
if(CCTKARGNUM_emask == -1) CCTKARGNUM_emask = CCTK_VarIndex("SPACEMASK::emask"); \
if(CCTKGROUPNUM_mask == -1) CCTKGROUPNUM_mask = CCTK_GroupIndex("SPACEMASK::mask"); \
if(CCTKARGNUM_r == -1) CCTKARGNUM_r = CCTK_VarIndex("GRID::r"); \
if(CCTKGROUPNUM_coordinates == -1) CCTKGROUPNUM_coordinates = CCTK_GroupIndex("GRID::coordinates"); \
if(CCTKARGNUM_sf_active == -1) CCTKARGNUM_sf_active = CCTK_VarIndex("SPHERICALSURFACE::sf_active[0]"); \
if(CCTKGROUPNUM_sf_active == -1) CCTKGROUPNUM_sf_active = CCTK_GroupIndex("SPHERICALSURFACE::sf_active"); \
if(CCTKARGNUM_sf_area == -1) CCTKARGNUM_sf_area = CCTK_VarIndex("SPHERICALSURFACE::sf_area[0]"); \
if(CCTKGROUPNUM_sf_info == -1) CCTKGROUPNUM_sf_info = CCTK_GroupIndex("SPHERICALSURFACE::sf_info"); \
if(CCTKARGNUM_sf_centroid_x == -1) CCTKARGNUM_sf_centroid_x = CCTK_VarIndex("SPHERICALSURFACE::sf_centroid_x[0]"); \
if(CCTKARGNUM_sf_centroid_y == -1) CCTKARGNUM_sf_centroid_y = CCTK_VarIndex("SPHERICALSURFACE::sf_centroid_y[0]"); \
if(CCTKARGNUM_sf_centroid_z == -1) CCTKARGNUM_sf_centroid_z = CCTK_VarIndex("SPHERICALSURFACE::sf_centroid_z[0]"); \
if(CCTKARGNUM_sf_delta_phi == -1) CCTKARGNUM_sf_delta_phi = CCTK_VarIndex("SPHERICALSURFACE::sf_delta_phi[0]"); \
if(CCTKGROUPNUM_sf_coordinate_descriptors == -1) CCTKGROUPNUM_sf_coordinate_descriptors = CCTK_GroupIndex("SPHERICALSURFACE::sf_coordinate_descriptors"); \
if(CCTKARGNUM_sf_delta_theta == -1) CCTKARGNUM_sf_delta_theta = CCTK_VarIndex("SPHERICALSURFACE::sf_delta_theta[0]"); \
if(CCTKARGNUM_sf_max_radius == -1) CCTKARGNUM_sf_max_radius = CCTK_VarIndex("SPHERICALSURFACE::sf_max_radius[0]"); \
if(CCTKARGNUM_sf_max_x == -1) CCTKARGNUM_sf_max_x = CCTK_VarIndex("SPHERICALSURFACE::sf_max_x[0]"); \
if(CCTKARGNUM_sf_max_y == -1) CCTKARGNUM_sf_max_y = CCTK_VarIndex("SPHERICALSURFACE::sf_max_y[0]"); \
if(CCTKARGNUM_sf_max_z == -1) CCTKARGNUM_sf_max_z = CCTK_VarIndex("SPHERICALSURFACE::sf_max_z[0]"); \
if(CCTKARGNUM_sf_maxreflevel == -1) CCTKARGNUM_sf_maxreflevel = CCTK_VarIndex("SPHERICALSURFACE::sf_maxreflevel[0]"); \
if(CCTKGROUPNUM_sf_maxreflevel == -1) CCTKGROUPNUM_sf_maxreflevel = CCTK_GroupIndex("SPHERICALSURFACE::sf_maxreflevel"); \
if(CCTKARGNUM_sf_mean_radius == -1) CCTKARGNUM_sf_mean_radius = CCTK_VarIndex("SPHERICALSURFACE::sf_mean_radius[0]"); \
if(CCTKARGNUM_sf_min_radius == -1) CCTKARGNUM_sf_min_radius = CCTK_VarIndex("SPHERICALSURFACE::sf_min_radius[0]"); \
if(CCTKARGNUM_sf_min_x == -1) CCTKARGNUM_sf_min_x = CCTK_VarIndex("SPHERICALSURFACE::sf_min_x[0]"); \
if(CCTKARGNUM_sf_min_y == -1) CCTKARGNUM_sf_min_y = CCTK_VarIndex("SPHERICALSURFACE::sf_min_y[0]"); \
if(CCTKARGNUM_sf_min_z == -1) CCTKARGNUM_sf_min_z = CCTK_VarIndex("SPHERICALSURFACE::sf_min_z[0]"); \
if(CCTKARGNUM_sf_minreflevel == -1) CCTKARGNUM_sf_minreflevel = CCTK_VarIndex("SPHERICALSURFACE::sf_minreflevel[0]"); \
if(CCTKGROUPNUM_sf_minreflevel == -1) CCTKGROUPNUM_sf_minreflevel = CCTK_GroupIndex("SPHERICALSURFACE::sf_minreflevel"); \
if(CCTKARGNUM_sf_nghostsphi == -1) CCTKARGNUM_sf_nghostsphi = CCTK_VarIndex("SPHERICALSURFACE::sf_nghostsphi[0]"); \
if(CCTKGROUPNUM_sf_shape_descriptors == -1) CCTKGROUPNUM_sf_shape_descriptors = CCTK_GroupIndex("SPHERICALSURFACE::sf_shape_descriptors"); \
if(CCTKARGNUM_sf_nghoststheta == -1) CCTKARGNUM_sf_nghoststheta = CCTK_VarIndex("SPHERICALSURFACE::sf_nghoststheta[0]"); \
if(CCTKARGNUM_sf_nphi == -1) CCTKARGNUM_sf_nphi = CCTK_VarIndex("SPHERICALSURFACE::sf_nphi[0]"); \
if(CCTKARGNUM_sf_ntheta == -1) CCTKARGNUM_sf_ntheta = CCTK_VarIndex("SPHERICALSURFACE::sf_ntheta[0]"); \
if(CCTKARGNUM_sf_origin_phi == -1) CCTKARGNUM_sf_origin_phi = CCTK_VarIndex("SPHERICALSURFACE::sf_origin_phi[0]"); \
if(CCTKARGNUM_sf_origin_theta == -1) CCTKARGNUM_sf_origin_theta = CCTK_VarIndex("SPHERICALSURFACE::sf_origin_theta[0]"); \
if(CCTKARGNUM_sf_origin_x == -1) CCTKARGNUM_sf_origin_x = CCTK_VarIndex("SPHERICALSURFACE::sf_origin_x[0]"); \
if(CCTKGROUPNUM_sf_origin == -1) CCTKGROUPNUM_sf_origin = CCTK_GroupIndex("SPHERICALSURFACE::sf_origin"); \
if(CCTKARGNUM_sf_origin_y == -1) CCTKARGNUM_sf_origin_y = CCTK_VarIndex("SPHERICALSURFACE::sf_origin_y[0]"); \
if(CCTKARGNUM_sf_origin_z == -1) CCTKARGNUM_sf_origin_z = CCTK_VarIndex("SPHERICALSURFACE::sf_origin_z[0]"); \
if(CCTKARGNUM_sf_quadrupole_xx == -1) CCTKARGNUM_sf_quadrupole_xx = CCTK_VarIndex("SPHERICALSURFACE::sf_quadrupole_xx[0]"); \
if(CCTKARGNUM_sf_quadrupole_xy == -1) CCTKARGNUM_sf_quadrupole_xy = CCTK_VarIndex("SPHERICALSURFACE::sf_quadrupole_xy[0]"); \
if(CCTKARGNUM_sf_quadrupole_xz == -1) CCTKARGNUM_sf_quadrupole_xz = CCTK_VarIndex("SPHERICALSURFACE::sf_quadrupole_xz[0]"); \
if(CCTKARGNUM_sf_quadrupole_yy == -1) CCTKARGNUM_sf_quadrupole_yy = CCTK_VarIndex("SPHERICALSURFACE::sf_quadrupole_yy[0]"); \
if(CCTKARGNUM_sf_quadrupole_yz == -1) CCTKARGNUM_sf_quadrupole_yz = CCTK_VarIndex("SPHERICALSURFACE::sf_quadrupole_yz[0]"); \
if(CCTKARGNUM_sf_quadrupole_zz == -1) CCTKARGNUM_sf_quadrupole_zz = CCTK_VarIndex("SPHERICALSURFACE::sf_quadrupole_zz[0]"); \
if(CCTKARGNUM_sf_radius == -1) CCTKARGNUM_sf_radius = CCTK_VarIndex("SPHERICALSURFACE::sf_radius[0]"); \
if(CCTKGROUPNUM_sf_radius == -1) CCTKGROUPNUM_sf_radius = CCTK_GroupIndex("SPHERICALSURFACE::sf_radius"); \
if(CCTKARGNUM_sf_valid == -1) CCTKARGNUM_sf_valid = CCTK_VarIndex("SPHERICALSURFACE::sf_valid[0]"); \
if(CCTKGROUPNUM_sf_valid == -1) CCTKGROUPNUM_sf_valid = CCTK_GroupIndex("SPHERICALSURFACE::sf_valid"); \
if(CCTKARGNUM_space_mask == -1) CCTKARGNUM_space_mask = CCTK_VarIndex("SPACEMASK::space_mask"); \
if(CCTKGROUPNUM_space_mask_group == -1) CCTKGROUPNUM_space_mask_group = CCTK_GroupIndex("SPACEMASK::space_mask_group"); \
if(CCTKARGNUM_x == -1) CCTKARGNUM_x = CCTK_VarIndex("GRID::x"); \
if(CCTKARGNUM_y == -1) CCTKARGNUM_y = CCTK_VarIndex("GRID::y"); \
if(CCTKARGNUM_z == -1) CCTKARGNUM_z = CCTK_VarIndex("GRID::z");

#define DISSIPATION_PUBLIC_C2F_PROTO \
const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_INT *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_INT *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_INT *,CCTK_INT *,CCTK_INT *,CCTK_INT *,CCTK_INT *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_INT *,CCTK_INT *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *

#define PASS_DISSIPATION_PUBLIC_C2F(GH) \
PASS_GROUPSIZE(coordinates, 0),\
PASS_GROUPSIZE(mask, 0),\
PASS_GROUPSIZE(sf_radius, 0),\
PASS_GROUPSIZE(space_mask_group, 0),\
PASS_GROUPSIZE(coordinates, 1),\
PASS_GROUPSIZE(mask, 1),\
PASS_GROUPSIZE(sf_radius, 1),\
PASS_GROUPSIZE(space_mask_group, 1),\
PASS_GROUPSIZE(coordinates, 2),\
PASS_GROUPSIZE(mask, 2),\
PASS_GROUPSIZE(space_mask_group, 2),\
PASS_GROUPLEN(SPHERICALSURFACE, sf_active),\
PASS_GROUPLEN(SPHERICALSURFACE, sf_coordinate_descriptors),\
PASS_GROUPLEN(SPHERICALSURFACE, sf_info),\
PASS_GROUPLEN(SPHERICALSURFACE, sf_maxreflevel),\
PASS_GROUPLEN(SPHERICALSURFACE, sf_minreflevel),\
PASS_GROUPLEN(SPHERICALSURFACE, sf_origin),\
PASS_GROUPLEN(SPHERICALSURFACE, sf_radius),\
PASS_GROUPLEN(SPHERICALSURFACE, sf_shape_descriptors),\
PASS_GROUPLEN(SPHERICALSURFACE, sf_valid),\
(CCTK_REAL *)(PASS_REFERENCE(coarse_dx, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(coarse_dy, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(coarse_dz, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(emask, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(r, 0)),\
(CCTK_INT *)(PASS_REFERENCE(sf_active, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_area, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_centroid_x, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_centroid_y, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_centroid_z, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_delta_phi, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_delta_theta, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_max_radius, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_max_x, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_max_y, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_max_z, 0)),\
(CCTK_INT *)(PASS_REFERENCE(sf_maxreflevel, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_mean_radius, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_min_radius, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_min_x, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_min_y, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_min_z, 0)),\
(CCTK_INT *)(PASS_REFERENCE(sf_minreflevel, 0)),\
(CCTK_INT *)(PASS_REFERENCE(sf_nghostsphi, 0)),\
(CCTK_INT *)(PASS_REFERENCE(sf_nghoststheta, 0)),\
(CCTK_INT *)(PASS_REFERENCE(sf_nphi, 0)),\
(CCTK_INT *)(PASS_REFERENCE(sf_ntheta, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_origin_phi, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_origin_theta, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_origin_x, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_origin_y, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_origin_z, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_quadrupole_xx, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_quadrupole_xy, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_quadrupole_xz, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_quadrupole_yy, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_quadrupole_yz, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_quadrupole_zz, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(sf_radius, 0)),\
(CCTK_INT *)(PASS_REFERENCE(sf_valid, 0)),\
(CCTK_INT *)(PASS_REFERENCE(space_mask, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(x, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(y, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(z, 0))

#endif /* CCODE */

#ifdef FCODE
#define DISSIPATION_FARGUMENTS _CCTK_FARGUMENTS, DISSIPATION_PRIVATE_FARGUMENTS, DISSIPATION_PUBLIC_FARGUMENTS

#define DECLARE_DISSIPATION_FARGUMENTS _DECLARE_CCTK_FARGUMENTS DECLARE_DISSIPATION_PRIVATE_FARGUMENTS DECLARE_DISSIPATION_PUBLIC_FARGUMENTS

#endif /* FCODE */

#ifdef CCODE
#define DECLARE_DISSIPATION_CARGUMENTS _DECLARE_CCTK_CARGUMENTS DECLARE_DISSIPATION_PRIVATE_CARGUMENTS DECLARE_DISSIPATION_PUBLIC_CARGUMENTS

#define DISSIPATION_C2F_PROTO _CCTK_C2F_PROTO, DISSIPATION_PRIVATE_C2F_PROTO, DISSIPATION_PUBLIC_C2F_PROTO

#define PASS_DISSIPATION_C2F(GH) _PASS_CCTK_C2F(GH), PASS_DISSIPATION_PRIVATE_C2F(GH), PASS_DISSIPATION_PUBLIC_C2F(GH)

#define DECLARE_DISSIPATION_C2F _DECLARE_CCTK_C2F DECLARE_DISSIPATION_PRIVATE_C2F DECLARE_DISSIPATION_PUBLIC_C2F

#define INITIALISE_DISSIPATION_C2F _INITIALISE_CCTK_C2F INITIALISE_DISSIPATION_PRIVATE_C2F INITIALISE_DISSIPATION_PUBLIC_C2F

#define DISSIPATION_CARGUMENTS cGH *cctkGH

#endif /* CCODE */