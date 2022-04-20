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
CCTK_STRING  action_if_found&&\
CCTK_STRING  check_for&&\
CCTK_STRING  check_vars&&\
CCTK_STRING  verbose&&\
CCTK_INT check_after&&\
CCTK_INT check_every&&\
CCTK_INT out_NaNmask&&\
CCTK_INT report_max&&\
COMMON /NaNCheckerpriv/action_if_found,check_for,check_vars,verbose,check_after,check_every,out_NaNmask,report_max&&\

