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
CCTK_REAL  A_const&&\
CCTK_REAL  BigMass&&\
CCTK_REAL  RoM_inner&&\
CCTK_REAL  bbh_r_const1&&\
CCTK_REAL  bbh_r_const2&&\
CCTK_REAL  bbh_r_const3&&\
CCTK_REAL  bbh_r_const4&&\
CCTK_REAL  bbh_rhosurf1&&\
CCTK_REAL  bbh_rhosurf2&&\
CCTK_REAL  bbh_rhosurf3&&\
CCTK_REAL  bbh_rhosurf4&&\
CCTK_REAL  betam1&&\
CCTK_REAL  disk_powerlaw_P_max&&\
CCTK_REAL  disk_powerlaw_rho_star_max&&\
CCTK_REAL  disk_powerlaw_tau_max&&\
CCTK_REAL  ell_inner&&\
CCTK_REAL  p_c&&\
CCTK_REAL  q&&\
CCTK_REAL  rho_cut&&\
CCTK_REAL  rho_cut_o_rho_max&&\
CCTK_REAL  sam_disk&&\
CCTK_REAL  surfxmax&&\
CCTK_REAL  surfxmin&&\
CCTK_REAL  surfzmax&&\
CCTK_INT arithrhosurf&&\
CCTK_INT em_field_type&&\
CCTK_INT metric_type&&\
CCTK_INT nsurfrho&&\
CCTK_INT ntot_disk_powerlaw&&\
CCTK_INT surfphinum&&\
CCTK_INT surfxnum&&\
CCTK_INT surfznum&&\
COMMON /disk_powerlawpriv/A_const,BigMass,RoM_inner,bbh_r_const1,bbh_r_const2,bbh_r_const3,bbh_r_const4,bbh_rhosurf1,bbh_rhosurf2,bbh_rhosurf3,bbh_rhosurf4,betam1,disk_powerlaw_P_max,disk_powerlaw_rho_star_max,disk_powerlaw_tau_max,ell_inner,p_c,q,rho_cut,rho_cut_o_rho_max,sam_disk,surfxmax,surfxmin,surfzmax,arithrhosurf,em_field_type,metric_type,nsurfrho,ntot_disk_powerlaw,surfphinum,surfxnum,surfznum&&\
CCTK_REAL  CCTKH5&&\
CCTK_REAL  CCTKH6&&\
CCTK_REAL  CCTKH7&&\
CCTK_REAL  CCTKH8&&\
CCTK_REAL  binary_orb_freq&&\
CCTK_REAL  CCTKH10&&\
CCTK_REAL  CCTKH11&&\
CCTK_REAL  CCTKH12&&\
CCTK_REAL  CCTKH13&&\
CCTK_REAL  CCTKH14&&\
CCTK_REAL  CCTKH15&&\
CCTK_REAL  rah1_initial&&\
CCTK_REAL  rah2_initial&&\
CCTK_REAL  CCTKH16&&\
CCTK_REAL  xbh1_initial&&\
CCTK_REAL  xbh2_initial&&\
CCTK_REAL  CCTKH18&&\
CCTK_INT CCTKH0&&\
CCTK_INT CCTKH1&&\
CCTK_INT CCTKH2&&\
CCTK_INT CCTKH3&&\
CCTK_INT CCTKH4&&\
CCTK_INT CCTKH9&&\
CCTK_INT enable_lower_order_at_boundaries&&\
CCTK_INT CCTKH17&&\
CCTK_INT CCTKH19&&\
COMMON /BSSNrest/CCTKH5,CCTKH6,CCTKH7,CCTKH8,binary_orb_freq,CCTKH10,CCTKH11,CCTKH12,CCTKH13,CCTKH14,CCTKH15,rah1_initial,rah2_initial,CCTKH16,xbh1_initial,xbh2_initial,CCTKH18,CCTKH0,CCTKH1,CCTKH2,CCTKH3,CCTKH4,CCTKH9,enable_lower_order_at_boundaries,CCTKH17,CCTKH19&&\
CCTK_REAL  CCTKH20&&\
CCTK_REAL  CCTKH21&&\
CCTK_REAL  CCTKH27&&\
CCTK_REAL  inner_lum_rad_ratio&&\
CCTK_REAL  lum_outer_rad&&\
CCTK_REAL  CCTKH29&&\
CCTK_REAL  CCTKH30&&\
CCTK_REAL  CCTKH31&&\
CCTK_REAL  CCTKH32&&\
CCTK_REAL  CCTKH33&&\
CCTK_REAL  CCTKH34&&\
CCTK_REAL  rhosurf_rmax&&\
CCTK_REAL  rhosurf_rmin&&\
CCTK_REAL  CCTKH35&&\
CCTK_INT N_rad_ray&&\
CCTK_INT CCTKH22&&\
CCTK_INT CCTKH23&&\
CCTK_INT CCTKH24&&\
CCTK_INT CCTKH25&&\
CCTK_INT CCTKH26&&\
CCTK_INT CCTKH28&&\
COMMON /DIAGNOSTICS_MHDrest/CCTKH20,CCTKH21,CCTKH27,inner_lum_rad_ratio,lum_outer_rad,CCTKH29,CCTKH30,CCTKH31,CCTKH32,CCTKH33,CCTKH34,rhosurf_rmax,rhosurf_rmin,CCTKH35,N_rad_ray,CCTKH22,CCTKH23,CCTKH24,CCTKH25,CCTKH26,CCTKH28&&\
CCTK_REAL  CCTKH36&&\
CCTK_REAL  CCTKH39&&\
CCTK_REAL  ddrbddr&&\
CCTK_REAL  CCTKH40&&\
CCTK_REAL  drbdr&&\
CCTK_REAL  CCTKH44&&\
CCTK_REAL  rbr&&\
CCTK_REAL  rsurf2&&\
CCTK_REAL  surf_radius&&\
CCTK_INT CCTKH37&&\
CCTK_INT N_phi&&\
CCTK_INT N_theta&&\
CCTK_INT CCTKH38&&\
CCTK_INT arithsurf&&\
CCTK_INT CCTKH41&&\
CCTK_INT CCTKH42&&\
CCTK_INT CCTKH43&&\
CCTK_INT nsurf&&\
CCTK_INT CCTKH45&&\
CCTK_INT CCTKH46&&\
CCTK_INT CCTKH47&&\
CCTK_INT CCTKH48&&\
CCTK_INT scaledsurf&&\
CCTK_INT sym_factor&&\
COMMON /DIAGNOSTICS_VACUUMrest/CCTKH36,CCTKH39,ddrbddr,CCTKH40,drbdr,CCTKH44,rbr,rsurf2,surf_radius,CCTKH37,N_phi,N_theta,CCTKH38,arithsurf,CCTKH41,CCTKH42,CCTKH43,nsurf,CCTKH45,CCTKH46,CCTKH47,CCTKH48,scaledsurf,sym_factor&&\
CCTK_REAL  CCTKH49&&\
CCTK_REAL  CCTKH50&&\
CCTK_REAL  CCTKH51&&\
CCTK_REAL  CCTKH52&&\
CCTK_REAL  CCTKH62&&\
CCTK_REAL  CCTKH63&&\
CCTK_REAL  CCTKH64&&\
CCTK_REAL  CCTKH65&&\
CCTK_REAL  CCTKH66&&\
CCTK_REAL  CCTKH67&&\
CCTK_REAL  CCTKH68&&\
CCTK_REAL  CCTKH69&&\
CCTK_STRING  bitant_plane&&\
CCTK_STRING  domain&&\
CCTK_STRING  CCTKH53&&\
CCTK_STRING  CCTKH54&&\
CCTK_STRING  CCTKH61&&\
CCTK_INT CCTKH55&&\
CCTK_INT CCTKH56&&\
CCTK_INT CCTKH57&&\
CCTK_INT CCTKH58&&\
CCTK_INT CCTKH59&&\
CCTK_INT CCTKH60&&\
COMMON /GRIDrest/CCTKH49,CCTKH50,CCTKH51,CCTKH52,CCTKH62,CCTKH63,CCTKH64,CCTKH65,CCTKH66,CCTKH67,CCTKH68,CCTKH69,bitant_plane,domain,CCTKH53,CCTKH54,CCTKH61,CCTKH55,CCTKH56,CCTKH57,CCTKH58,CCTKH59,CCTKH60&&\
CCTK_REAL  CCTKH70&&\
CCTK_REAL  CCTKH71&&\
CCTK_REAL  phi_GW&&\
CCTK_REAL  radius_GW&&\
CCTK_REAL  CCTKH72&&\
CCTK_REAL  theta_GW&&\
COMMON /GW_EXTRACTIONrest/CCTKH70,CCTKH71,phi_GW,radius_GW,CCTKH72,theta_GW&&\
CCTK_REAL  CCTKH78&&\
CCTK_REAL  CCTKH92&&\
CCTK_REAL  CCTKH100&&\
CCTK_REAL  CCTKH102&&\
CCTK_REAL  CCTKH104&&\
CCTK_REAL  CCTKH106&&\
CCTK_REAL  CCTKH108&&\
CCTK_REAL  CCTKH110&&\
CCTK_REAL  CCTKH112&&\
CCTK_REAL  CCTKH114&&\
CCTK_REAL  CCTKH116&&\
CCTK_STRING  CCTKH75&&\
CCTK_STRING  CCTKH76&&\
CCTK_STRING  CCTKH79&&\
CCTK_STRING  CCTKH82&&\
CCTK_STRING  CCTKH83&&\
CCTK_STRING  CCTKH84&&\
CCTK_STRING  CCTKH87&&\
CCTK_STRING  CCTKH88&&\
CCTK_STRING  CCTKH93&&\
CCTK_STRING  CCTKH94&&\
CCTK_STRING  CCTKH96&&\
CCTK_STRING  CCTKH118&&\
CCTK_STRING  CCTKH120&&\
CCTK_STRING  CCTKH122&&\
CCTK_STRING  CCTKH124&&\
CCTK_STRING  CCTKH125&&\
CCTK_STRING  CCTKH128&&\
CCTK_INT CCTKH73&&\
CCTK_INT CCTKH74&&\
CCTK_INT CCTKH77&&\
CCTK_INT CCTKH80&&\
CCTK_INT CCTKH81&&\
CCTK_INT CCTKH85&&\
CCTK_INT CCTKH86&&\
CCTK_INT CCTKH89&&\
CCTK_INT CCTKH90&&\
CCTK_INT CCTKH91&&\
CCTK_INT out_every&&\
CCTK_INT CCTKH95&&\
CCTK_INT CCTKH97&&\
CCTK_INT CCTKH98&&\
CCTK_INT CCTKH99&&\
CCTK_INT CCTKH101&&\
CCTK_INT CCTKH103&&\
CCTK_INT CCTKH105&&\
CCTK_INT CCTKH107&&\
CCTK_INT CCTKH109&&\
CCTK_INT CCTKH111&&\
CCTK_INT CCTKH113&&\
CCTK_INT CCTKH115&&\
CCTK_INT CCTKH117&&\
CCTK_INT CCTKH119&&\
CCTK_INT CCTKH121&&\
CCTK_INT CCTKH123&&\
CCTK_INT CCTKH126&&\
CCTK_INT CCTKH127&&\
COMMON /IOrest/CCTKH78,CCTKH92,CCTKH100,CCTKH102,CCTKH104,CCTKH106,CCTKH108,CCTKH110,CCTKH112,CCTKH114,CCTKH116,CCTKH75,CCTKH76,CCTKH79,CCTKH82,CCTKH83,CCTKH84,CCTKH87,CCTKH88,CCTKH93,CCTKH94,CCTKH96,CCTKH118,CCTKH120,CCTKH122,CCTKH124,CCTKH125,CCTKH128,CCTKH73,CCTKH74,CCTKH77,CCTKH80,CCTKH81,CCTKH85,CCTKH86,CCTKH89,CCTKH90,CCTKH91,out_every,CCTKH95,CCTKH97,CCTKH98,CCTKH99,CCTKH101,CCTKH103,CCTKH105,CCTKH107,CCTKH109,CCTKH111,CCTKH113,CCTKH115,CCTKH117,CCTKH119,CCTKH121,CCTKH123,CCTKH126,CCTKH127&&\
CCTK_REAL  CCTKH131&&\
CCTK_REAL  CCTKH132&&\
CCTK_REAL  K_poly&&\
CCTK_REAL  CCTKH135&&\
CCTK_REAL  CCTKH137&&\
CCTK_REAL  P_fact&&\
CCTK_REAL  CCTKH138&&\
CCTK_REAL  CCTKH139&&\
CCTK_REAL  Sym_Bz&&\
CCTK_REAL  CCTKH141&&\
CCTK_REAL  CCTKH145&&\
CCTK_REAL  CCTKH146&&\
CCTK_REAL  CCTKH149&&\
CCTK_REAL  CCTKH154&&\
CCTK_REAL  CCTKH155&&\
CCTK_REAL  CCTKH157&&\
CCTK_REAL  gamma_th&&\
CCTK_REAL  CCTKH161&&\
CCTK_REAL  CCTKH163&&\
CCTK_REAL  CCTKH164&&\
CCTK_REAL  CCTKH165&&\
CCTK_REAL  CCTKH166&&\
CCTK_REAL  CCTKH168&&\
CCTK_REAL  pfloor&&\
CCTK_REAL  CCTKH171&&\
CCTK_REAL  CCTKH172&&\
CCTK_REAL  CCTKH173&&\
CCTK_REAL  CCTKH174&&\
CCTK_REAL  CCTKH175&&\
CCTK_REAL  CCTKH177&&\
CCTK_REAL  CCTKH178&&\
CCTK_REAL  CCTKH182&&\
CCTK_REAL  CCTKH183&&\
CCTK_REAL  rho_b_atm&&\
CCTK_REAL  rho_b_max&&\
CCTK_REAL  rho_fact&&\
CCTK_REAL  CCTKH186&&\
CCTK_REAL  sdots_o_rhot&&\
CCTK_REAL  CCTKH187&&\
CCTK_REAL  tau_atm&&\
CCTK_REAL  tau_fact&&\
CCTK_INT CCTKH129&&\
CCTK_INT CCTKH130&&\
CCTK_INT CCTKH133&&\
CCTK_INT CCTKH134&&\
CCTK_INT CCTKH136&&\
CCTK_INT CCTKH140&&\
CCTK_INT CCTKH142&&\
CCTK_INT CCTKH143&&\
CCTK_INT CCTKH144&&\
CCTK_INT CCTKH147&&\
CCTK_INT constrained_transport_scheme&&\
CCTK_INT CCTKH148&&\
CCTK_INT em_evolve_enable&&\
CCTK_INT CCTKH150&&\
CCTK_INT enable_HARM_energyvariable&&\
CCTK_INT CCTKH151&&\
CCTK_INT CCTKH152&&\
CCTK_INT enable_primitives_disk&&\
CCTK_INT enable_shocktest_primitive_mode&&\
CCTK_INT CCTKH153&&\
CCTK_INT eps_flag&&\
CCTK_INT CCTKH156&&\
CCTK_INT force_font_fix_fail&&\
CCTK_INT CCTKH158&&\
CCTK_INT CCTKH159&&\
CCTK_INT CCTKH160&&\
CCTK_INT CCTKH162&&\
CCTK_INT CCTKH167&&\
CCTK_INT neos&&\
CCTK_INT CCTKH169&&\
CCTK_INT CCTKH170&&\
CCTK_INT primitives_solver&&\
CCTK_INT CCTKH176&&\
CCTK_INT CCTKH179&&\
CCTK_INT CCTKH180&&\
CCTK_INT CCTKH181&&\
CCTK_INT CCTKH184&&\
CCTK_INT CCTKH185&&\
CCTK_INT tau_stildefix_enable&&\
CCTK_INT use_HARM_primitives&&\
CCTK_INT CCTKH188&&\
CCTK_INT CCTKH189&&\
COMMON /MHD_EVOLVErest/CCTKH131,CCTKH132,K_poly,CCTKH135,CCTKH137,P_fact,CCTKH138,CCTKH139,Sym_Bz,CCTKH141,CCTKH145,CCTKH146,CCTKH149,CCTKH154,CCTKH155,CCTKH157,gamma_th,CCTKH161,CCTKH163,CCTKH164,CCTKH165,CCTKH166,CCTKH168,pfloor,CCTKH171,CCTKH172,CCTKH173,CCTKH174,CCTKH175,CCTKH177,CCTKH178,CCTKH182,CCTKH183,rho_b_atm,rho_b_max,rho_fact,CCTKH186,sdots_o_rhot,CCTKH187,tau_atm,tau_fact,CCTKH129,CCTKH130,CCTKH133,CCTKH134,CCTKH136,CCTKH140,CCTKH142,CCTKH143,CCTKH144,CCTKH147,constrained_transport_scheme,CCTKH148,em_evolve_enable,CCTKH150,enable_HARM_energyvariable,CCTKH151,CCTKH152,enable_primitives_disk,enable_shocktest_primitive_mode,CCTKH153,eps_flag,CCTKH156,force_font_fix_fail,CCTKH158,CCTKH159,CCTKH160,CCTKH162,CCTKH167,neos,CCTKH169,CCTKH170,primitives_solver,CCTKH176,CCTKH179,CCTKH180,CCTKH181,CCTKH184,CCTKH185,tau_stildefix_enable,use_HARM_primitives,CCTKH188,CCTKH189&&\

