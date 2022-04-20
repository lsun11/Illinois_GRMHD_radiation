#define DECLARE_CCTK_PARAMETERS \
CCTK_REAL  excision_radius&&\
CCTK_REAL  run_time&&\
CCTK_INT Symmetry&&\
CCTK_INT bssn_enable&&\
CCTK_INT cowling_enable&&\
CCTK_INT excision_enable&&\
CCTK_INT fisheye_enable&&\
CCTK_INT iter_count&&\
CCTK_INT number_of_mol_ministeps&&\
CCTK_INT rot_metric&&\
CCTK_INT trA_detg_enforce&&\
COMMON /cctk_params_global/excision_radius,run_time,Symmetry,bssn_enable,cowling_enable,excision_enable,fisheye_enable,iter_count,number_of_mol_ministeps,rot_metric,trA_detg_enforce&&\
CCTK_REAL  auto_res_ratio(42)&&\
CCTK_STRING  auto_res_grid(42)&&\
CCTK_INT auto_res(42)&&\
CCTK_INT maxnphi&&\
CCTK_INT maxntheta&&\
CCTK_INT nghostsphi(42)&&\
CCTK_INT nghoststheta(42)&&\
CCTK_INT nphi(42)&&\
CCTK_INT nsurfaces&&\
CCTK_INT ntheta(42)&&\
CCTK_INT symmetric_x(42)&&\
CCTK_INT symmetric_y(42)&&\
CCTK_INT symmetric_z(42)&&\
CCTK_INT verbose&&\
COMMON /SphericalSurfacerest/auto_res_ratio,auto_res_grid,auto_res,maxnphi,maxntheta,nghostsphi,nghoststheta,nphi,nsurfaces,ntheta,symmetric_x,symmetric_y,symmetric_z,verbose&&\
CCTK_REAL  origin_x(42)&&\
CCTK_REAL  origin_y(42)&&\
CCTK_REAL  origin_z(42)&&\
CCTK_REAL  radius(42)&&\
CCTK_REAL  radius_x(42)&&\
CCTK_REAL  radius_y(42)&&\
CCTK_REAL  radius_z(42)&&\
CCTK_INT set_elliptic(42)&&\
CCTK_INT set_spherical(42)&&\
COMMON /SphericalSurfacepriv/origin_x,origin_y,origin_z,radius,radius_x,radius_y,radius_z,set_elliptic,set_spherical&&\
