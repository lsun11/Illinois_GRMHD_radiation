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
CCTK_REAL  outInfo_dt&&\
CCTK_REAL  real_max&&\
CCTK_REAL  real_min&&\
CCTK_STRING  outInfo_criterion&&\
CCTK_STRING  outInfo_reductions&&\
CCTK_STRING  outInfo_vars&&\
CCTK_INT int_width&&\
CCTK_INT iter_width&&\
CCTK_INT outHeader_every&&\
CCTK_INT outInfo_every&&\
CCTK_INT real_prec&&\
CCTK_INT real_prec_sci&&\
CCTK_INT real_width&&\
CCTK_INT time_prec&&\
CCTK_INT time_width&&\
COMMON /CarpetIOBasicpriv/outInfo_dt,real_max,real_min,outInfo_criterion,outInfo_reductions,outInfo_vars,int_width,iter_width,outHeader_every,outInfo_every,real_prec,real_prec_sci,real_width,time_prec,time_width&&\
CCTK_REAL  CCTKH5&&\
CCTK_REAL  out_dt&&\
CCTK_REAL  CCTKH25&&\
CCTK_REAL  CCTKH27&&\
CCTK_REAL  CCTKH29&&\
CCTK_REAL  CCTKH31&&\
CCTK_REAL  CCTKH33&&\
CCTK_REAL  CCTKH35&&\
CCTK_REAL  CCTKH37&&\
CCTK_REAL  CCTKH39&&\
CCTK_REAL  CCTKH41&&\
CCTK_STRING  CCTKH2&&\
CCTK_STRING  CCTKH3&&\
CCTK_STRING  CCTKH6&&\
CCTK_STRING  CCTKH9&&\
CCTK_STRING  CCTKH10&&\
CCTK_STRING  CCTKH11&&\
CCTK_STRING  out_criterion&&\
CCTK_STRING  CCTKH14&&\
CCTK_STRING  CCTKH18&&\
CCTK_STRING  CCTKH19&&\
CCTK_STRING  CCTKH21&&\
CCTK_STRING  CCTKH43&&\
CCTK_STRING  CCTKH45&&\
CCTK_STRING  CCTKH47&&\
CCTK_STRING  CCTKH49&&\
CCTK_STRING  CCTKH50&&\
CCTK_STRING  CCTKH52&&\
CCTK_INT CCTKH0&&\
CCTK_INT CCTKH1&&\
CCTK_INT CCTKH4&&\
CCTK_INT CCTKH7&&\
CCTK_INT CCTKH8&&\
CCTK_INT CCTKH12&&\
CCTK_INT CCTKH13&&\
CCTK_INT CCTKH15&&\
CCTK_INT CCTKH16&&\
CCTK_INT CCTKH17&&\
CCTK_INT out_every&&\
CCTK_INT CCTKH20&&\
CCTK_INT CCTKH22&&\
CCTK_INT CCTKH23&&\
CCTK_INT CCTKH24&&\
CCTK_INT CCTKH26&&\
CCTK_INT CCTKH28&&\
CCTK_INT CCTKH30&&\
CCTK_INT CCTKH32&&\
CCTK_INT CCTKH34&&\
CCTK_INT CCTKH36&&\
CCTK_INT CCTKH38&&\
CCTK_INT CCTKH40&&\
CCTK_INT CCTKH42&&\
CCTK_INT CCTKH44&&\
CCTK_INT CCTKH46&&\
CCTK_INT CCTKH48&&\
CCTK_INT CCTKH51&&\
CCTK_INT strict_io_parameter_check&&\
COMMON /IOrest/CCTKH5,out_dt,CCTKH25,CCTKH27,CCTKH29,CCTKH31,CCTKH33,CCTKH35,CCTKH37,CCTKH39,CCTKH41,CCTKH2,CCTKH3,CCTKH6,CCTKH9,CCTKH10,CCTKH11,out_criterion,CCTKH14,CCTKH18,CCTKH19,CCTKH21,CCTKH43,CCTKH45,CCTKH47,CCTKH49,CCTKH50,CCTKH52,CCTKH0,CCTKH1,CCTKH4,CCTKH7,CCTKH8,CCTKH12,CCTKH13,CCTKH15,CCTKH16,CCTKH17,out_every,CCTKH20,CCTKH22,CCTKH23,CCTKH24,CCTKH26,CCTKH28,CCTKH30,CCTKH32,CCTKH34,CCTKH36,CCTKH38,CCTKH40,CCTKH42,CCTKH44,CCTKH46,CCTKH48,CCTKH51,strict_io_parameter_check&&\
