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
CCTK_REAL  C_ko&&\
COMMON /excisionrest/C_ko&&\

