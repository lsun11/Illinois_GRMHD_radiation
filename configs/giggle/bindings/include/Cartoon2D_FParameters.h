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
CCTK_STRING  new_mask_excised_name&&\
CCTK_STRING  new_mask_field_name&&\
CCTK_STRING  new_style_excision_var&&\
CCTK_STRING  old_style_excision_var&&\
CCTK_INT allow_grid_resize&&\
CCTK_INT cartoon_active&&\
CCTK_INT eno_order&&\
CCTK_INT new_excision&&\
CCTK_INT old_excision&&\
CCTK_INT order&&\
CCTK_INT stencil&&\
CCTK_INT verbose&&\
COMMON /Cartoon2Dpriv/new_mask_excised_name,new_mask_field_name,new_style_excision_var,old_style_excision_var,allow_grid_resize,cartoon_active,eno_order,new_excision,old_excision,order,stencil,verbose&&\

