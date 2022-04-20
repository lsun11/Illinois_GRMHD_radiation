#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL Aphi_power;
  CCTK_REAL BigM;
  CCTK_REAL Con_Pol;
  CCTK_REAL Con_Tor;
  CCTK_REAL Erad_over_rho;
  CCTK_REAL Iloop;
  CCTK_REAL Iloop2;
  CCTK_REAL MRI_wavelength_calculator_dxy;
  CCTK_REAL Omega_value;
  CCTK_REAL P_deplete_bhns;
  CCTK_REAL a2oa1_perturb;
  CCTK_REAL alpha_rho_cut_off;
  CCTK_REAL ampl_perturb;
  CCTK_REAL angle_to_tilt_magnetic_fields;
  CCTK_REAL angle_to_tilt_magnetic_fieldsII;
  CCTK_REAL beta_ext;
  CCTK_REAL betam1;
  CCTK_REAL betam1_2;
  CCTK_REAL bhns_P_max;
  CCTK_REAL bhns_P_max2;
  CCTK_REAL bhns_R_NS;
  CCTK_REAL bhns_R_NS2;
  CCTK_REAL bhns_avg_betam1;
  CCTK_REAL bhns_bh_filling_radius[100];
  CCTK_REAL bhns_fac_atm;
  CCTK_REAL bhns_max_b2;
  CCTK_REAL bhns_max_b2_posn_x;
  CCTK_REAL bhns_max_b2_posn_y;
  CCTK_REAL bhns_max_b2_posn_z;
  CCTK_REAL bhns_max_rho_b;
  CCTK_REAL bhns_rho_b_max_posn_x;
  CCTK_REAL bhns_rho_b_max_posn_y;
  CCTK_REAL bhns_rho_b_max_posn_z;
  CCTK_REAL bhns_rhob_fac;
  CCTK_REAL bhns_rhob_max;
  CCTK_REAL bhns_tracer_r;
  CCTK_REAL bhns_tracer_rin;
  CCTK_REAL bhns_tracer_x0;
  CCTK_REAL bhns_tracer_y0;
  CCTK_REAL bhns_tracer_z0;
  CCTK_REAL bhns_vx_CM;
  CCTK_REAL bhns_vx_CM2;
  CCTK_REAL bhns_vy_CM;
  CCTK_REAL bhns_vy_CM2;
  CCTK_REAL bhns_vz_CM;
  CCTK_REAL bhns_vz_CM2;
  CCTK_REAL conloop1;
  CCTK_REAL conloop2;
  CCTK_REAL horiz_radius;
  CCTK_REAL initial_ns2_coord_x;
  CCTK_REAL initial_ns2_coord_y;
  CCTK_REAL initial_ns_coord_x;
  CCTK_REAL initial_ns_coord_y;
  CCTK_REAL initial_wd_coord_x;
  CCTK_REAL initial_wd_coord_y;
  CCTK_REAL lambda_perturb;
  CCTK_REAL moncrief_radius_GW[11];
  CCTK_REAL mythbusters_boost_factor;
  CCTK_REAL p_c;
  CCTK_REAL particle_cone_angle;
  CCTK_REAL particle_cylinder_cone_zmax;
  CCTK_REAL particle_cylinder_cone_zmin;
  CCTK_REAL r0;
  CCTK_REAL r0_ah;
  CCTK_REAL rad_rhob_fac;
  CCTK_REAL radi_perturb;
  CCTK_REAL rhob_fac2;
  CCTK_REAL rhob_o_b2;
  CCTK_REAL rhobatm_falloff_power;
  CCTK_REAL rloop;
  CCTK_REAL rloop2;
  CCTK_REAL surfxmax;
  CCTK_REAL surfxmin;
  CCTK_REAL surfzmax;
  CCTK_REAL tracer_x_max;
  CCTK_REAL tracer_x_min;
  CCTK_REAL tracer_y_max;
  CCTK_REAL tracer_y_min;
  CCTK_REAL tracer_z_max;
  CCTK_REAL tracer_z_min;
  CCTK_REAL xh0;
  CCTK_REAL yh0;
  CCTK_REAL zh0;
  CCTK_INT ATM_TYPE;
  CCTK_INT INPUTARRAY_PHISIZE;
  CCTK_INT INPUTARRAY_THETASIZE;
  CCTK_INT ITERATION_TO_BOOST_MAGNETIC_FIELDS;
  CCTK_INT ITERATION_TO_INSERT_MAGNETIC_FIELDS;
  CCTK_INT ITERATION_TO_output_MRI_wavelength;
  CCTK_INT MRI_wavelength_calculator_Nxy;
  CCTK_INT NUM_ZERO_PTS;
  CCTK_INT N_particles_to_trace;
  CCTK_INT RADIAL_INTERP_ORDER;
  CCTK_INT RESET_RHO_B_ATM;
  CCTK_INT alpha_diagnostic;
  CCTK_INT angle_to_tilt_90;
  CCTK_INT bhns_B_v;
  CCTK_INT bhns_domain;
  CCTK_INT bhns_particle_tracer_start;
  CCTK_INT bhns_regrid_input_enable;
  CCTK_INT bhns_regrid_output_enable_iter;
  CCTK_INT center_Bfields_around_BH_instead;
  CCTK_INT chunklet_dump_every;
  CCTK_INT chunklet_procs_at_a_time;
  CCTK_INT em_field_type;
  CCTK_INT genID_cmdline_output_enable;
  CCTK_INT initial_particle_geometry;
  CCTK_INT moncrief_gw_num_radii;
  CCTK_INT nperturb;
  CCTK_INT ntot_bhns;
  CCTK_INT particle_center;
  CCTK_INT particle_tracer_substep_every;
  CCTK_INT piecewise;
  CCTK_INT refill_horizons_magfields_every;
  CCTK_INT reset_shift_lapse;
  CCTK_INT subtract_off_Omega_r;
  CCTK_INT superposition;
  CCTK_INT surfphinum;
  CCTK_INT surfxnum;
  CCTK_INT surfznum;
  CCTK_INT two_ns;
  CCTK_INT unequalmass;
  CCTK_INT use_new_bhns_initial_data;
} PRIVATE_BHNS_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_BHNS_STRUCT_PARAMS \
  CCTK_REAL const Aphi_power = PRIVATE_BHNS_STRUCT.Aphi_power; \
  CCTK_REAL const BigM = PRIVATE_BHNS_STRUCT.BigM; \
  CCTK_REAL const Con_Pol = PRIVATE_BHNS_STRUCT.Con_Pol; \
  CCTK_REAL const Con_Tor = PRIVATE_BHNS_STRUCT.Con_Tor; \
  CCTK_REAL const Erad_over_rho = PRIVATE_BHNS_STRUCT.Erad_over_rho; \
  CCTK_REAL const Iloop = PRIVATE_BHNS_STRUCT.Iloop; \
  CCTK_REAL const Iloop2 = PRIVATE_BHNS_STRUCT.Iloop2; \
  CCTK_REAL const MRI_wavelength_calculator_dxy = PRIVATE_BHNS_STRUCT.MRI_wavelength_calculator_dxy; \
  CCTK_REAL const Omega_value = PRIVATE_BHNS_STRUCT.Omega_value; \
  CCTK_REAL const P_deplete_bhns = PRIVATE_BHNS_STRUCT.P_deplete_bhns; \
  CCTK_REAL const a2oa1_perturb = PRIVATE_BHNS_STRUCT.a2oa1_perturb; \
  CCTK_REAL const alpha_rho_cut_off = PRIVATE_BHNS_STRUCT.alpha_rho_cut_off; \
  CCTK_REAL const ampl_perturb = PRIVATE_BHNS_STRUCT.ampl_perturb; \
  CCTK_REAL const angle_to_tilt_magnetic_fields = PRIVATE_BHNS_STRUCT.angle_to_tilt_magnetic_fields; \
  CCTK_REAL const angle_to_tilt_magnetic_fieldsII = PRIVATE_BHNS_STRUCT.angle_to_tilt_magnetic_fieldsII; \
  CCTK_REAL const beta_ext = PRIVATE_BHNS_STRUCT.beta_ext; \
  CCTK_REAL const betam1 = PRIVATE_BHNS_STRUCT.betam1; \
  CCTK_REAL const betam1_2 = PRIVATE_BHNS_STRUCT.betam1_2; \
  CCTK_REAL const bhns_P_max = PRIVATE_BHNS_STRUCT.bhns_P_max; \
  CCTK_REAL const bhns_P_max2 = PRIVATE_BHNS_STRUCT.bhns_P_max2; \
  CCTK_REAL const bhns_R_NS = PRIVATE_BHNS_STRUCT.bhns_R_NS; \
  CCTK_REAL const bhns_R_NS2 = PRIVATE_BHNS_STRUCT.bhns_R_NS2; \
  CCTK_REAL const bhns_avg_betam1 = PRIVATE_BHNS_STRUCT.bhns_avg_betam1; \
  CCTK_REAL const * const bhns_bh_filling_radius = PRIVATE_BHNS_STRUCT.bhns_bh_filling_radius; \
  CCTK_REAL const bhns_fac_atm = PRIVATE_BHNS_STRUCT.bhns_fac_atm; \
  CCTK_REAL const bhns_max_b2 = PRIVATE_BHNS_STRUCT.bhns_max_b2; \
  CCTK_REAL const bhns_max_b2_posn_x = PRIVATE_BHNS_STRUCT.bhns_max_b2_posn_x; \
  CCTK_REAL const bhns_max_b2_posn_y = PRIVATE_BHNS_STRUCT.bhns_max_b2_posn_y; \
  CCTK_REAL const bhns_max_b2_posn_z = PRIVATE_BHNS_STRUCT.bhns_max_b2_posn_z; \
  CCTK_REAL const bhns_max_rho_b = PRIVATE_BHNS_STRUCT.bhns_max_rho_b; \
  CCTK_REAL const bhns_rho_b_max_posn_x = PRIVATE_BHNS_STRUCT.bhns_rho_b_max_posn_x; \
  CCTK_REAL const bhns_rho_b_max_posn_y = PRIVATE_BHNS_STRUCT.bhns_rho_b_max_posn_y; \
  CCTK_REAL const bhns_rho_b_max_posn_z = PRIVATE_BHNS_STRUCT.bhns_rho_b_max_posn_z; \
  CCTK_REAL const bhns_rhob_fac = PRIVATE_BHNS_STRUCT.bhns_rhob_fac; \
  CCTK_REAL const bhns_rhob_max = PRIVATE_BHNS_STRUCT.bhns_rhob_max; \
  CCTK_REAL const bhns_tracer_r = PRIVATE_BHNS_STRUCT.bhns_tracer_r; \
  CCTK_REAL const bhns_tracer_rin = PRIVATE_BHNS_STRUCT.bhns_tracer_rin; \
  CCTK_REAL const bhns_tracer_x0 = PRIVATE_BHNS_STRUCT.bhns_tracer_x0; \
  CCTK_REAL const bhns_tracer_y0 = PRIVATE_BHNS_STRUCT.bhns_tracer_y0; \
  CCTK_REAL const bhns_tracer_z0 = PRIVATE_BHNS_STRUCT.bhns_tracer_z0; \
  CCTK_REAL const bhns_vx_CM = PRIVATE_BHNS_STRUCT.bhns_vx_CM; \
  CCTK_REAL const bhns_vx_CM2 = PRIVATE_BHNS_STRUCT.bhns_vx_CM2; \
  CCTK_REAL const bhns_vy_CM = PRIVATE_BHNS_STRUCT.bhns_vy_CM; \
  CCTK_REAL const bhns_vy_CM2 = PRIVATE_BHNS_STRUCT.bhns_vy_CM2; \
  CCTK_REAL const bhns_vz_CM = PRIVATE_BHNS_STRUCT.bhns_vz_CM; \
  CCTK_REAL const bhns_vz_CM2 = PRIVATE_BHNS_STRUCT.bhns_vz_CM2; \
  CCTK_REAL const conloop1 = PRIVATE_BHNS_STRUCT.conloop1; \
  CCTK_REAL const conloop2 = PRIVATE_BHNS_STRUCT.conloop2; \
  CCTK_REAL const horiz_radius = PRIVATE_BHNS_STRUCT.horiz_radius; \
  CCTK_REAL const initial_ns2_coord_x = PRIVATE_BHNS_STRUCT.initial_ns2_coord_x; \
  CCTK_REAL const initial_ns2_coord_y = PRIVATE_BHNS_STRUCT.initial_ns2_coord_y; \
  CCTK_REAL const initial_ns_coord_x = PRIVATE_BHNS_STRUCT.initial_ns_coord_x; \
  CCTK_REAL const initial_ns_coord_y = PRIVATE_BHNS_STRUCT.initial_ns_coord_y; \
  CCTK_REAL const initial_wd_coord_x = PRIVATE_BHNS_STRUCT.initial_wd_coord_x; \
  CCTK_REAL const initial_wd_coord_y = PRIVATE_BHNS_STRUCT.initial_wd_coord_y; \
  CCTK_REAL const lambda_perturb = PRIVATE_BHNS_STRUCT.lambda_perturb; \
  CCTK_REAL const * const moncrief_radius_GW = PRIVATE_BHNS_STRUCT.moncrief_radius_GW; \
  CCTK_REAL const mythbusters_boost_factor = PRIVATE_BHNS_STRUCT.mythbusters_boost_factor; \
  CCTK_REAL const p_c = PRIVATE_BHNS_STRUCT.p_c; \
  CCTK_REAL const particle_cone_angle = PRIVATE_BHNS_STRUCT.particle_cone_angle; \
  CCTK_REAL const particle_cylinder_cone_zmax = PRIVATE_BHNS_STRUCT.particle_cylinder_cone_zmax; \
  CCTK_REAL const particle_cylinder_cone_zmin = PRIVATE_BHNS_STRUCT.particle_cylinder_cone_zmin; \
  CCTK_REAL const r0 = PRIVATE_BHNS_STRUCT.r0; \
  CCTK_REAL const r0_ah = PRIVATE_BHNS_STRUCT.r0_ah; \
  CCTK_REAL const rad_rhob_fac = PRIVATE_BHNS_STRUCT.rad_rhob_fac; \
  CCTK_REAL const radi_perturb = PRIVATE_BHNS_STRUCT.radi_perturb; \
  CCTK_REAL const rhob_fac2 = PRIVATE_BHNS_STRUCT.rhob_fac2; \
  CCTK_REAL const rhob_o_b2 = PRIVATE_BHNS_STRUCT.rhob_o_b2; \
  CCTK_REAL const rhobatm_falloff_power = PRIVATE_BHNS_STRUCT.rhobatm_falloff_power; \
  CCTK_REAL const rloop = PRIVATE_BHNS_STRUCT.rloop; \
  CCTK_REAL const rloop2 = PRIVATE_BHNS_STRUCT.rloop2; \
  CCTK_REAL const surfxmax = PRIVATE_BHNS_STRUCT.surfxmax; \
  CCTK_REAL const surfxmin = PRIVATE_BHNS_STRUCT.surfxmin; \
  CCTK_REAL const surfzmax = PRIVATE_BHNS_STRUCT.surfzmax; \
  CCTK_REAL const tracer_x_max = PRIVATE_BHNS_STRUCT.tracer_x_max; \
  CCTK_REAL const tracer_x_min = PRIVATE_BHNS_STRUCT.tracer_x_min; \
  CCTK_REAL const tracer_y_max = PRIVATE_BHNS_STRUCT.tracer_y_max; \
  CCTK_REAL const tracer_y_min = PRIVATE_BHNS_STRUCT.tracer_y_min; \
  CCTK_REAL const tracer_z_max = PRIVATE_BHNS_STRUCT.tracer_z_max; \
  CCTK_REAL const tracer_z_min = PRIVATE_BHNS_STRUCT.tracer_z_min; \
  CCTK_REAL const xh0 = PRIVATE_BHNS_STRUCT.xh0; \
  CCTK_REAL const yh0 = PRIVATE_BHNS_STRUCT.yh0; \
  CCTK_REAL const zh0 = PRIVATE_BHNS_STRUCT.zh0; \
  CCTK_INT const ATM_TYPE = PRIVATE_BHNS_STRUCT.ATM_TYPE; \
  CCTK_INT const INPUTARRAY_PHISIZE = PRIVATE_BHNS_STRUCT.INPUTARRAY_PHISIZE; \
  CCTK_INT const INPUTARRAY_THETASIZE = PRIVATE_BHNS_STRUCT.INPUTARRAY_THETASIZE; \
  CCTK_INT const ITERATION_TO_BOOST_MAGNETIC_FIELDS = PRIVATE_BHNS_STRUCT.ITERATION_TO_BOOST_MAGNETIC_FIELDS; \
  CCTK_INT const ITERATION_TO_INSERT_MAGNETIC_FIELDS = PRIVATE_BHNS_STRUCT.ITERATION_TO_INSERT_MAGNETIC_FIELDS; \
  CCTK_INT const ITERATION_TO_output_MRI_wavelength = PRIVATE_BHNS_STRUCT.ITERATION_TO_output_MRI_wavelength; \
  CCTK_INT const MRI_wavelength_calculator_Nxy = PRIVATE_BHNS_STRUCT.MRI_wavelength_calculator_Nxy; \
  CCTK_INT const NUM_ZERO_PTS = PRIVATE_BHNS_STRUCT.NUM_ZERO_PTS; \
  CCTK_INT const N_particles_to_trace = PRIVATE_BHNS_STRUCT.N_particles_to_trace; \
  CCTK_INT const RADIAL_INTERP_ORDER = PRIVATE_BHNS_STRUCT.RADIAL_INTERP_ORDER; \
  CCTK_INT const RESET_RHO_B_ATM = PRIVATE_BHNS_STRUCT.RESET_RHO_B_ATM; \
  CCTK_INT const alpha_diagnostic = PRIVATE_BHNS_STRUCT.alpha_diagnostic; \
  CCTK_INT const angle_to_tilt_90 = PRIVATE_BHNS_STRUCT.angle_to_tilt_90; \
  CCTK_INT const bhns_B_v = PRIVATE_BHNS_STRUCT.bhns_B_v; \
  CCTK_INT const bhns_domain = PRIVATE_BHNS_STRUCT.bhns_domain; \
  CCTK_INT const bhns_particle_tracer_start = PRIVATE_BHNS_STRUCT.bhns_particle_tracer_start; \
  CCTK_INT const bhns_regrid_input_enable = PRIVATE_BHNS_STRUCT.bhns_regrid_input_enable; \
  CCTK_INT const bhns_regrid_output_enable_iter = PRIVATE_BHNS_STRUCT.bhns_regrid_output_enable_iter; \
  CCTK_INT const center_Bfields_around_BH_instead = PRIVATE_BHNS_STRUCT.center_Bfields_around_BH_instead; \
  CCTK_INT const chunklet_dump_every = PRIVATE_BHNS_STRUCT.chunklet_dump_every; \
  CCTK_INT const chunklet_procs_at_a_time = PRIVATE_BHNS_STRUCT.chunklet_procs_at_a_time; \
  CCTK_INT const em_field_type = PRIVATE_BHNS_STRUCT.em_field_type; \
  CCTK_INT const genID_cmdline_output_enable = PRIVATE_BHNS_STRUCT.genID_cmdline_output_enable; \
  CCTK_INT const initial_particle_geometry = PRIVATE_BHNS_STRUCT.initial_particle_geometry; \
  CCTK_INT const moncrief_gw_num_radii = PRIVATE_BHNS_STRUCT.moncrief_gw_num_radii; \
  CCTK_INT const nperturb = PRIVATE_BHNS_STRUCT.nperturb; \
  CCTK_INT const ntot_bhns = PRIVATE_BHNS_STRUCT.ntot_bhns; \
  CCTK_INT const particle_center = PRIVATE_BHNS_STRUCT.particle_center; \
  CCTK_INT const particle_tracer_substep_every = PRIVATE_BHNS_STRUCT.particle_tracer_substep_every; \
  CCTK_INT const piecewise = PRIVATE_BHNS_STRUCT.piecewise; \
  CCTK_INT const refill_horizons_magfields_every = PRIVATE_BHNS_STRUCT.refill_horizons_magfields_every; \
  CCTK_INT const reset_shift_lapse = PRIVATE_BHNS_STRUCT.reset_shift_lapse; \
  CCTK_INT const subtract_off_Omega_r = PRIVATE_BHNS_STRUCT.subtract_off_Omega_r; \
  CCTK_INT const superposition = PRIVATE_BHNS_STRUCT.superposition; \
  CCTK_INT const surfphinum = PRIVATE_BHNS_STRUCT.surfphinum; \
  CCTK_INT const surfxnum = PRIVATE_BHNS_STRUCT.surfxnum; \
  CCTK_INT const surfznum = PRIVATE_BHNS_STRUCT.surfznum; \
  CCTK_INT const two_ns = PRIVATE_BHNS_STRUCT.two_ns; \
  CCTK_INT const unequalmass = PRIVATE_BHNS_STRUCT.unequalmass; \
  CCTK_INT const use_new_bhns_initial_data = PRIVATE_BHNS_STRUCT.use_new_bhns_initial_data; \
  enum { \
      dummy_PRIVATE_BHNS_STRUCT_Aphi_power = sizeof( Aphi_power ) \
    , dummy_PRIVATE_BHNS_STRUCT_BigM = sizeof( BigM ) \
    , dummy_PRIVATE_BHNS_STRUCT_Con_Pol = sizeof( Con_Pol ) \
    , dummy_PRIVATE_BHNS_STRUCT_Con_Tor = sizeof( Con_Tor ) \
    , dummy_PRIVATE_BHNS_STRUCT_Erad_over_rho = sizeof( Erad_over_rho ) \
    , dummy_PRIVATE_BHNS_STRUCT_Iloop = sizeof( Iloop ) \
    , dummy_PRIVATE_BHNS_STRUCT_Iloop2 = sizeof( Iloop2 ) \
    , dummy_PRIVATE_BHNS_STRUCT_MRI_wavelength_calculator_dxy = sizeof( MRI_wavelength_calculator_dxy ) \
    , dummy_PRIVATE_BHNS_STRUCT_Omega_value = sizeof( Omega_value ) \
    , dummy_PRIVATE_BHNS_STRUCT_P_deplete_bhns = sizeof( P_deplete_bhns ) \
    , dummy_PRIVATE_BHNS_STRUCT_a2oa1_perturb = sizeof( a2oa1_perturb ) \
    , dummy_PRIVATE_BHNS_STRUCT_alpha_rho_cut_off = sizeof( alpha_rho_cut_off ) \
    , dummy_PRIVATE_BHNS_STRUCT_ampl_perturb = sizeof( ampl_perturb ) \
    , dummy_PRIVATE_BHNS_STRUCT_angle_to_tilt_magnetic_fields = sizeof( angle_to_tilt_magnetic_fields ) \
    , dummy_PRIVATE_BHNS_STRUCT_angle_to_tilt_magnetic_fieldsII = sizeof( angle_to_tilt_magnetic_fieldsII ) \
    , dummy_PRIVATE_BHNS_STRUCT_beta_ext = sizeof( beta_ext ) \
    , dummy_PRIVATE_BHNS_STRUCT_betam1 = sizeof( betam1 ) \
    , dummy_PRIVATE_BHNS_STRUCT_betam1_2 = sizeof( betam1_2 ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_P_max = sizeof( bhns_P_max ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_P_max2 = sizeof( bhns_P_max2 ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_R_NS = sizeof( bhns_R_NS ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_R_NS2 = sizeof( bhns_R_NS2 ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_avg_betam1 = sizeof( bhns_avg_betam1 ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_bh_filling_radius = sizeof( bhns_bh_filling_radius ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_fac_atm = sizeof( bhns_fac_atm ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_max_b2 = sizeof( bhns_max_b2 ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_max_b2_posn_x = sizeof( bhns_max_b2_posn_x ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_max_b2_posn_y = sizeof( bhns_max_b2_posn_y ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_max_b2_posn_z = sizeof( bhns_max_b2_posn_z ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_max_rho_b = sizeof( bhns_max_rho_b ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_rho_b_max_posn_x = sizeof( bhns_rho_b_max_posn_x ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_rho_b_max_posn_y = sizeof( bhns_rho_b_max_posn_y ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_rho_b_max_posn_z = sizeof( bhns_rho_b_max_posn_z ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_rhob_fac = sizeof( bhns_rhob_fac ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_rhob_max = sizeof( bhns_rhob_max ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_tracer_r = sizeof( bhns_tracer_r ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_tracer_rin = sizeof( bhns_tracer_rin ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_tracer_x0 = sizeof( bhns_tracer_x0 ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_tracer_y0 = sizeof( bhns_tracer_y0 ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_tracer_z0 = sizeof( bhns_tracer_z0 ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_vx_CM = sizeof( bhns_vx_CM ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_vx_CM2 = sizeof( bhns_vx_CM2 ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_vy_CM = sizeof( bhns_vy_CM ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_vy_CM2 = sizeof( bhns_vy_CM2 ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_vz_CM = sizeof( bhns_vz_CM ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_vz_CM2 = sizeof( bhns_vz_CM2 ) \
    , dummy_PRIVATE_BHNS_STRUCT_conloop1 = sizeof( conloop1 ) \
    , dummy_PRIVATE_BHNS_STRUCT_conloop2 = sizeof( conloop2 ) \
    , dummy_PRIVATE_BHNS_STRUCT_horiz_radius = sizeof( horiz_radius ) \
    , dummy_PRIVATE_BHNS_STRUCT_initial_ns2_coord_x = sizeof( initial_ns2_coord_x ) \
    , dummy_PRIVATE_BHNS_STRUCT_initial_ns2_coord_y = sizeof( initial_ns2_coord_y ) \
    , dummy_PRIVATE_BHNS_STRUCT_initial_ns_coord_x = sizeof( initial_ns_coord_x ) \
    , dummy_PRIVATE_BHNS_STRUCT_initial_ns_coord_y = sizeof( initial_ns_coord_y ) \
    , dummy_PRIVATE_BHNS_STRUCT_initial_wd_coord_x = sizeof( initial_wd_coord_x ) \
    , dummy_PRIVATE_BHNS_STRUCT_initial_wd_coord_y = sizeof( initial_wd_coord_y ) \
    , dummy_PRIVATE_BHNS_STRUCT_lambda_perturb = sizeof( lambda_perturb ) \
    , dummy_PRIVATE_BHNS_STRUCT_moncrief_radius_GW = sizeof( moncrief_radius_GW ) \
    , dummy_PRIVATE_BHNS_STRUCT_mythbusters_boost_factor = sizeof( mythbusters_boost_factor ) \
    , dummy_PRIVATE_BHNS_STRUCT_p_c = sizeof( p_c ) \
    , dummy_PRIVATE_BHNS_STRUCT_particle_cone_angle = sizeof( particle_cone_angle ) \
    , dummy_PRIVATE_BHNS_STRUCT_particle_cylinder_cone_zmax = sizeof( particle_cylinder_cone_zmax ) \
    , dummy_PRIVATE_BHNS_STRUCT_particle_cylinder_cone_zmin = sizeof( particle_cylinder_cone_zmin ) \
    , dummy_PRIVATE_BHNS_STRUCT_r0 = sizeof( r0 ) \
    , dummy_PRIVATE_BHNS_STRUCT_r0_ah = sizeof( r0_ah ) \
    , dummy_PRIVATE_BHNS_STRUCT_rad_rhob_fac = sizeof( rad_rhob_fac ) \
    , dummy_PRIVATE_BHNS_STRUCT_radi_perturb = sizeof( radi_perturb ) \
    , dummy_PRIVATE_BHNS_STRUCT_rhob_fac2 = sizeof( rhob_fac2 ) \
    , dummy_PRIVATE_BHNS_STRUCT_rhob_o_b2 = sizeof( rhob_o_b2 ) \
    , dummy_PRIVATE_BHNS_STRUCT_rhobatm_falloff_power = sizeof( rhobatm_falloff_power ) \
    , dummy_PRIVATE_BHNS_STRUCT_rloop = sizeof( rloop ) \
    , dummy_PRIVATE_BHNS_STRUCT_rloop2 = sizeof( rloop2 ) \
    , dummy_PRIVATE_BHNS_STRUCT_surfxmax = sizeof( surfxmax ) \
    , dummy_PRIVATE_BHNS_STRUCT_surfxmin = sizeof( surfxmin ) \
    , dummy_PRIVATE_BHNS_STRUCT_surfzmax = sizeof( surfzmax ) \
    , dummy_PRIVATE_BHNS_STRUCT_tracer_x_max = sizeof( tracer_x_max ) \
    , dummy_PRIVATE_BHNS_STRUCT_tracer_x_min = sizeof( tracer_x_min ) \
    , dummy_PRIVATE_BHNS_STRUCT_tracer_y_max = sizeof( tracer_y_max ) \
    , dummy_PRIVATE_BHNS_STRUCT_tracer_y_min = sizeof( tracer_y_min ) \
    , dummy_PRIVATE_BHNS_STRUCT_tracer_z_max = sizeof( tracer_z_max ) \
    , dummy_PRIVATE_BHNS_STRUCT_tracer_z_min = sizeof( tracer_z_min ) \
    , dummy_PRIVATE_BHNS_STRUCT_xh0 = sizeof( xh0 ) \
    , dummy_PRIVATE_BHNS_STRUCT_yh0 = sizeof( yh0 ) \
    , dummy_PRIVATE_BHNS_STRUCT_zh0 = sizeof( zh0 ) \
    , dummy_PRIVATE_BHNS_STRUCT_ATM_TYPE = sizeof( ATM_TYPE ) \
    , dummy_PRIVATE_BHNS_STRUCT_INPUTARRAY_PHISIZE = sizeof( INPUTARRAY_PHISIZE ) \
    , dummy_PRIVATE_BHNS_STRUCT_INPUTARRAY_THETASIZE = sizeof( INPUTARRAY_THETASIZE ) \
    , dummy_PRIVATE_BHNS_STRUCT_ITERATION_TO_BOOST_MAGNETIC_FIELDS = sizeof( ITERATION_TO_BOOST_MAGNETIC_FIELDS ) \
    , dummy_PRIVATE_BHNS_STRUCT_ITERATION_TO_INSERT_MAGNETIC_FIELDS = sizeof( ITERATION_TO_INSERT_MAGNETIC_FIELDS ) \
    , dummy_PRIVATE_BHNS_STRUCT_ITERATION_TO_output_MRI_wavelength = sizeof( ITERATION_TO_output_MRI_wavelength ) \
    , dummy_PRIVATE_BHNS_STRUCT_MRI_wavelength_calculator_Nxy = sizeof( MRI_wavelength_calculator_Nxy ) \
    , dummy_PRIVATE_BHNS_STRUCT_NUM_ZERO_PTS = sizeof( NUM_ZERO_PTS ) \
    , dummy_PRIVATE_BHNS_STRUCT_N_particles_to_trace = sizeof( N_particles_to_trace ) \
    , dummy_PRIVATE_BHNS_STRUCT_RADIAL_INTERP_ORDER = sizeof( RADIAL_INTERP_ORDER ) \
    , dummy_PRIVATE_BHNS_STRUCT_RESET_RHO_B_ATM = sizeof( RESET_RHO_B_ATM ) \
    , dummy_PRIVATE_BHNS_STRUCT_alpha_diagnostic = sizeof( alpha_diagnostic ) \
    , dummy_PRIVATE_BHNS_STRUCT_angle_to_tilt_90 = sizeof( angle_to_tilt_90 ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_B_v = sizeof( bhns_B_v ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_domain = sizeof( bhns_domain ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_particle_tracer_start = sizeof( bhns_particle_tracer_start ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_regrid_input_enable = sizeof( bhns_regrid_input_enable ) \
    , dummy_PRIVATE_BHNS_STRUCT_bhns_regrid_output_enable_iter = sizeof( bhns_regrid_output_enable_iter ) \
    , dummy_PRIVATE_BHNS_STRUCT_center_Bfields_around_BH_instead = sizeof( center_Bfields_around_BH_instead ) \
    , dummy_PRIVATE_BHNS_STRUCT_chunklet_dump_every = sizeof( chunklet_dump_every ) \
    , dummy_PRIVATE_BHNS_STRUCT_chunklet_procs_at_a_time = sizeof( chunklet_procs_at_a_time ) \
    , dummy_PRIVATE_BHNS_STRUCT_em_field_type = sizeof( em_field_type ) \
    , dummy_PRIVATE_BHNS_STRUCT_genID_cmdline_output_enable = sizeof( genID_cmdline_output_enable ) \
    , dummy_PRIVATE_BHNS_STRUCT_initial_particle_geometry = sizeof( initial_particle_geometry ) \
    , dummy_PRIVATE_BHNS_STRUCT_moncrief_gw_num_radii = sizeof( moncrief_gw_num_radii ) \
    , dummy_PRIVATE_BHNS_STRUCT_nperturb = sizeof( nperturb ) \
    , dummy_PRIVATE_BHNS_STRUCT_ntot_bhns = sizeof( ntot_bhns ) \
    , dummy_PRIVATE_BHNS_STRUCT_particle_center = sizeof( particle_center ) \
    , dummy_PRIVATE_BHNS_STRUCT_particle_tracer_substep_every = sizeof( particle_tracer_substep_every ) \
    , dummy_PRIVATE_BHNS_STRUCT_piecewise = sizeof( piecewise ) \
    , dummy_PRIVATE_BHNS_STRUCT_refill_horizons_magfields_every = sizeof( refill_horizons_magfields_every ) \
    , dummy_PRIVATE_BHNS_STRUCT_reset_shift_lapse = sizeof( reset_shift_lapse ) \
    , dummy_PRIVATE_BHNS_STRUCT_subtract_off_Omega_r = sizeof( subtract_off_Omega_r ) \
    , dummy_PRIVATE_BHNS_STRUCT_superposition = sizeof( superposition ) \
    , dummy_PRIVATE_BHNS_STRUCT_surfphinum = sizeof( surfphinum ) \
    , dummy_PRIVATE_BHNS_STRUCT_surfxnum = sizeof( surfxnum ) \
    , dummy_PRIVATE_BHNS_STRUCT_surfznum = sizeof( surfznum ) \
    , dummy_PRIVATE_BHNS_STRUCT_two_ns = sizeof( two_ns ) \
    , dummy_PRIVATE_BHNS_STRUCT_unequalmass = sizeof( unequalmass ) \
    , dummy_PRIVATE_BHNS_STRUCT_use_new_bhns_initial_data = sizeof( use_new_bhns_initial_data ) \
  }; \

