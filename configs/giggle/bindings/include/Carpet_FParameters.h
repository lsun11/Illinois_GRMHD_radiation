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
CCTK_INT periodic&&\
CCTK_INT periodic_x&&\
CCTK_INT periodic_y&&\
CCTK_INT periodic_z&&\
COMMON /Driverrest/periodic,periodic_x,periodic_y,periodic_z&&\
CCTK_REAL  aspect_ratio_x&&\
CCTK_REAL  aspect_ratio_y&&\
CCTK_REAL  aspect_ratio_z&&\
CCTK_STRING  base_extents&&\
CCTK_STRING  base_outerbounds&&\
CCTK_STRING  grid_coordinates_filename&&\
CCTK_STRING  grid_structure_filename&&\
CCTK_STRING  model&&\
CCTK_STRING  processor_topology&&\
CCTK_STRING  refinement_centering&&\
CCTK_STRING  space_refinement_factors&&\
CCTK_STRING  time_refinement_factors&&\
CCTK_STRING  timer_file&&\
CCTK_INT adaptive_stepsize&&\
CCTK_INT additional_buffer_zones&&\
CCTK_INT barriers&&\
CCTK_INT check_for_poison&&\
CCTK_INT checksum_timelevels&&\
CCTK_INT constant_load_per_processor&&\
CCTK_INT convergence_factor&&\
CCTK_INT convergence_level&&\
CCTK_INT deadbeef&&\
CCTK_INT domain_from_coordbase&&\
CCTK_INT domain_from_multipatch&&\
CCTK_INT enable_all_storage&&\
CCTK_INT ghost_size&&\
CCTK_INT ghost_size_x&&\
CCTK_INT ghost_size_y&&\
CCTK_INT ghost_size_z&&\
CCTK_INT global_nsize&&\
CCTK_INT global_nx&&\
CCTK_INT global_ny&&\
CCTK_INT global_nz&&\
CCTK_INT init_3_timelevels&&\
CCTK_INT init_each_timelevel&&\
CCTK_INT init_fill_timelevels&&\
CCTK_INT max_core_size_MB&&\
CCTK_INT max_memory_size_MB&&\
CCTK_INT max_poison_locations&&\
CCTK_INT max_refinement_levels&&\
CCTK_INT min_points_per_proc&&\
CCTK_INT num_convergence_levels&&\
CCTK_INT num_integrator_substeps&&\
CCTK_INT num_maps&&\
CCTK_INT num_threads&&\
CCTK_INT output_internal_data&&\
CCTK_INT output_timers_every&&\
CCTK_INT poison_new_timelevels&&\
CCTK_INT poison_value&&\
CCTK_INT print_timestats_every&&\
CCTK_INT processor_topology_3d_x&&\
CCTK_INT processor_topology_3d_y&&\
CCTK_INT processor_topology_3d_z&&\
CCTK_INT prolongate_initial_data&&\
CCTK_INT prolongation_order_space&&\
CCTK_INT prolongation_order_time&&\
CCTK_INT refine_timestep&&\
CCTK_INT refinement_factor&&\
CCTK_INT regrid_during_initialisation&&\
CCTK_INT regrid_during_recovery&&\
CCTK_INT regrid_in_level_mode&&\
CCTK_INT schedule_barriers&&\
CCTK_INT split_direction&&\
CCTK_INT storage_verbose&&\
CCTK_INT suppress_restriction&&\
CCTK_INT sync_during_time_integration&&\
CCTK_INT use_buffer_zones&&\
CCTK_INT use_tapered_grids&&\
CCTK_INT verbose&&\
CCTK_INT veryverbose&&\
COMMON /Carpetpriv/aspect_ratio_x,aspect_ratio_y,aspect_ratio_z,base_extents,base_outerbounds,grid_coordinates_filename,grid_structure_filename,model,processor_topology,refinement_centering,space_refinement_factors,time_refinement_factors,timer_file,adaptive_stepsize,additional_buffer_zones,barriers,check_for_poison,checksum_timelevels,constant_load_per_processor,convergence_factor,convergence_level,deadbeef,domain_from_coordbase,domain_from_multipatch,enable_all_storage,ghost_size,ghost_size_x,ghost_size_y,ghost_size_z,global_nsize,global_nx,global_ny,global_nz,init_3_timelevels,init_each_timelevel,init_fill_timelevels,max_core_size_MB,max_memory_size_MB,max_poison_locations,max_refinement_levels,min_points_per_proc,num_convergence_levels,num_integrator_substeps,num_maps,num_threads,output_internal_data,output_timers_every,poison_new_timelevels,poison_value,print_timestats_every,processor_topology_3d_x,processor_topology_3d_y,processor_topology_3d_z,prolongate_initial_data,prolongation_order_space,prolongation_order_time,refine_timestep,refinement_factor,regrid_during_initialisation,regrid_during_recovery,regrid_in_level_mode,schedule_barriers,split_direction,storage_verbose,suppress_restriction,sync_during_time_integration,use_buffer_zones,use_tapered_grids,verbose,veryverbose&&\
CCTK_REAL  cctk_final_time&&\
CCTK_REAL  cctk_initial_time&&\
CCTK_REAL  max_runtime&&\
CCTK_STRING  terminate&&\
CCTK_INT cctk_itlast&&\
CCTK_INT terminate_next&&\
COMMON /CACTUSrest/cctk_final_time,cctk_initial_time,max_runtime,terminate,cctk_itlast,terminate_next&&\
CCTK_STRING  initial_data_setup_method&&\
COMMON /INITBASErest/initial_data_setup_method&&\
CCTK_REAL  CCTKH5&&\
CCTK_REAL  CCTKH18&&\
CCTK_REAL  CCTKH27&&\
CCTK_REAL  CCTKH29&&\
CCTK_REAL  CCTKH31&&\
CCTK_REAL  CCTKH33&&\
CCTK_REAL  CCTKH35&&\
CCTK_REAL  CCTKH37&&\
CCTK_REAL  CCTKH39&&\
CCTK_REAL  CCTKH41&&\
CCTK_REAL  CCTKH43&&\
CCTK_STRING  CCTKH2&&\
CCTK_STRING  CCTKH3&&\
CCTK_STRING  CCTKH6&&\
CCTK_STRING  CCTKH9&&\
CCTK_STRING  CCTKH10&&\
CCTK_STRING  CCTKH11&&\
CCTK_STRING  CCTKH14&&\
CCTK_STRING  out_dir&&\
CCTK_STRING  CCTKH20&&\
CCTK_STRING  CCTKH21&&\
CCTK_STRING  CCTKH23&&\
CCTK_STRING  CCTKH45&&\
CCTK_STRING  CCTKH47&&\
CCTK_STRING  CCTKH49&&\
CCTK_STRING  CCTKH51&&\
CCTK_STRING  CCTKH52&&\
CCTK_STRING  CCTKH55&&\
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
CCTK_INT CCTKH19&&\
CCTK_INT CCTKH22&&\
CCTK_INT CCTKH24&&\
CCTK_INT CCTKH25&&\
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
CCTK_INT CCTKH50&&\
CCTK_INT CCTKH53&&\
CCTK_INT CCTKH54&&\
COMMON /IOrest/CCTKH5,CCTKH18,CCTKH27,CCTKH29,CCTKH31,CCTKH33,CCTKH35,CCTKH37,CCTKH39,CCTKH41,CCTKH43,CCTKH2,CCTKH3,CCTKH6,CCTKH9,CCTKH10,CCTKH11,CCTKH14,out_dir,CCTKH20,CCTKH21,CCTKH23,CCTKH45,CCTKH47,CCTKH49,CCTKH51,CCTKH52,CCTKH55,CCTKH0,CCTKH1,CCTKH4,CCTKH7,CCTKH8,CCTKH12,CCTKH13,CCTKH15,CCTKH16,CCTKH17,CCTKH19,CCTKH22,CCTKH24,CCTKH25,CCTKH26,CCTKH28,CCTKH30,CCTKH32,CCTKH34,CCTKH36,CCTKH38,CCTKH40,CCTKH42,CCTKH44,CCTKH46,CCTKH48,CCTKH50,CCTKH53,CCTKH54&&\
