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
CCTK_REAL  ipoison&&\
CCTK_REAL  poison&&\
CCTK_INT barriers&&\
CCTK_INT check_tree_search&&\
CCTK_INT tree_search&&\
COMMON /CarpetInterppriv/ipoison,poison,barriers,check_tree_search,tree_search&&\
CCTK_REAL  CCTKH0&&\
CCTK_REAL  cctk_initial_time&&\
CCTK_REAL  CCTKH2&&\
CCTK_STRING  CCTKH3&&\
CCTK_INT CCTKH1&&\
CCTK_INT CCTKH4&&\
COMMON /CACTUSrest/CCTKH0,cctk_initial_time,CCTKH2,CCTKH3,CCTKH1,CCTKH4&&\

