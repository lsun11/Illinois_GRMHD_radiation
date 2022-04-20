subroutine BHNS_setup_emfield_part2_local_setBi(cctk_dim,cctk_gsh,cctk_lsh,cctk_lbnd,cctk_ubnd,cctk_lssh,cctk_from,cctk_to,cctk_bbo&
  &x,cctk_delta_time,cctk_time,cctk_delta_space,cctk_origin_space,cctk_levfac,cctk_levoff,cctk_levoffdenom,cctk_timefac,cctk_convle&
  &vel,cctk_convfac,cctk_nghostzones,cctk_iteration,cctkGH, X0Axyzintrgroup,X0Tracer_rhs_group,X0Tracergroup,X0last_ah_data,X1Axyzi&
  &ntrgroup,X1Tracer_rhs_group,X2Axyzintrgroup,Axintr,Ayintr,Azintr,M0dot_last,Pintr,ah_radii_last,bhnsEthermalNOW,bhns_rk4_particl&
  &e_tracer_step,bhnsinsertBNOW,bhnsinsertRADNOW,int_M0dot,local_density,m1_minus_u_lower_zero,rhobintr,t_last,taul_rhs,tracer_taul&
  &,tracer_x,tracer_xc,tracer_y,tracer_yc,tracer_z,u_x_local,u_y_local,x_rhs,xc_rhs,xh_last,y_rhs,yc_rhs,yh_last,z_rhs,zh_last, X0B&
  &SSN_AH,X0BSSN_aux_private,X0BSSN_aux_restrict2,X0BSSN_diag_restrict,X0BSSN_gupij,X0BSSN_matter,X0BSSN_refbd,X0BSSN_rhs,X0BSSN_va&
  &rs,X0NPphi2_em,X0OmegF,X0Psi4_gw,X0Stagger_Bs,X0ah_radius,X0bh_posns,X0diag_vacuum_gf,X0disk_atmosphere,X0em_Ax,X0em_Ay,X0em_Az,&
  &X0em_Blagrangemultiplier,X0em_Blagrangemultiplier_rhs,X0em_Phi,X0em_Phi_rhs,X0em_conservativex,X0em_conservativey,X0em_conservat&
  &ivez,X0em_rhsx,X0em_rhsy,X0em_rhsz,X0eos_params1,X0eos_params2,X0excision_int_gfs,X0field_line_variables,X0field_line_variables_&
  &rhs,X0fisheye_vars,X0gw_moment_arrays,X0have_global_bdry,X0lapse_derivatives,X0lapse_rhs,X0lapse_vars,X0lapse_vars_aux,X0metric_&
  &spher_pol_1,X0metric_spher_pol_2,X0metric_spher_pol_3,X0mhd_conservatives,X0mhd_nosync,X0mhd_primitives,X0mhd_rhs,X0mhd_sync_lr_&
  &B_quantities,X0mhd_sync_lr_hydro_quantities,X0mhd_sync_metric_facevals,X0mhd_sync_nablas,X0mhd_sync_nablas_diag,X0mhd_sync_rho_b&
  &r_rho_bl,X0mhd_temps,X0mhd_vs,X0micphys_conservatives,X0micphys_conservatives_rhs,X0micphys_sync_lr_hydro_quantities,X0micphys_s&
  &ync_nablas,X0microphys_primitives,X0phi_derivs,X0rad_conservatives,X0rad_conservatives_rhs,X0rad_pressure,X0rad_primitives,X0rad&
  &_sync_lr,X0rhovecs,X0shift_rhs,X0shift_vars,X0shift_vars_temp,X0surf_params,X1BSSN_AH,X1BSSN_aux_private,X1BSSN_aux_restrict2,X1&
  &BSSN_diag_restrict,X1BSSN_gupij,X1BSSN_matter,X1BSSN_refbd,X1BSSN_rhs,X1BSSN_vars,X1NPphi2_em,X1OmegF,X1Psi4_gw,X1Stagger_Bs,X1a&
  &h_radius,X1diag_vacuum_gf,X1disk_atmosphere,X1em_Ax,X1em_Ay,X1em_Az,X1em_Blagrangemultiplier,X1em_Blagrangemultiplier_rhs,X1em_P&
  &hi,X1em_Phi_rhs,X1em_conservativex,X1em_conservativey,X1em_conservativez,X1em_rhsx,X1em_rhsy,X1em_rhsz,X1excision_int_gfs,X1fiel&
  &d_line_variables,X1field_line_variables_rhs,X1fisheye_vars,X1have_global_bdry,X1lapse_derivatives,X1lapse_rhs,X1lapse_vars,X1lap&
  &se_vars_aux,X1mhd_conservatives,X1mhd_nosync,X1mhd_primitives,X1mhd_rhs,X1mhd_sync_lr_B_quantities,X1mhd_sync_lr_hydro_quantitie&
  &s,X1mhd_sync_metric_facevals,X1mhd_sync_nablas,X1mhd_sync_nablas_diag,X1mhd_sync_rho_br_rho_bl,X1mhd_temps,X1mhd_vs,X1micphys_co&
  &nservatives,X1micphys_conservatives_rhs,X1micphys_sync_lr_hydro_quantities,X1micphys_sync_nablas,X1microphys_primitives,X1phi_de&
  &rivs,X1rad_conservatives,X1rad_conservatives_rhs,X1rad_pressure,X1rad_primitives,X1rad_sync_lr,X1shift_rhs,X1shift_vars,X1shift_&
  &vars_temp,X2BSSN_AH,X2BSSN_aux_private,X2BSSN_aux_restrict2,X2BSSN_diag_restrict,X2BSSN_gupij,X2BSSN_matter,X2BSSN_refbd,X2BSSN_&
  &rhs,X2BSSN_vars,X2NPphi2_em,X2OmegF,X2Psi4_gw,X2Stagger_Bs,X2ah_radius,X2diag_vacuum_gf,X2disk_atmosphere,X2em_Ax,X2em_Ay,X2em_A&
  &z,X2em_Blagrangemultiplier,X2em_Blagrangemultiplier_rhs,X2em_Phi,X2em_Phi_rhs,X2em_conservativex,X2em_conservativey,X2em_conserv&
  &ativez,X2em_rhsx,X2em_rhsy,X2em_rhsz,X2excision_int_gfs,X2field_line_variables,X2field_line_variables_rhs,X2fisheye_vars,X2lapse&
  &_derivatives,X2lapse_rhs,X2lapse_vars,X2lapse_vars_aux,X2mhd_conservatives,X2mhd_nosync,X2mhd_primitives,X2mhd_rhs,X2mhd_sync_lr&
  &_B_quantities,X2mhd_sync_lr_hydro_quantities,X2mhd_sync_metric_facevals,X2mhd_sync_nablas,X2mhd_sync_nablas_diag,X2mhd_sync_rho_&
  &br_rho_bl,X2mhd_temps,X2mhd_vs,X2micphys_conservatives,X2micphys_conservatives_rhs,X2micphys_sync_lr_hydro_quantities,X2micphys_&
  &sync_nablas,X2microphys_primitives,X2phi_derivs,X2rad_conservatives,X2rad_conservatives_rhs,X2rad_pressure,X2rad_primitives,X2ra&
  &d_sync_lr,X2shift_rhs,X2shift_vars,X2shift_vars_temp,ah_centroid_length,ah_flags_length,ah_origin_length,ah_radius_length,Aupxx,&
  &Aupxy,Aupxz,Aupyy,Aupyz,Aupzz,Ax,Ax_p,Ax_p_p,Axx,Axx_p,Axx_p_p,Axx_rhs,Axy,Axy_p,Axy_p_p,Axy_rhs,Axz,Axz_p,Axz_p_p,Axz_rhs,Ay,Ay&
  &_p,Ay_p_p,Ayy,Ayy_p,Ayy_p_p,Ayy_rhs,Ayz,Ayz_p,Ayz_p_p,Ayz_rhs,Az,Az_p,Az_p_p,Azz,Azz_p,Azz_p_p,Azz_rhs,Bfreq1,Bfreq1_p,Bfreq1_p_&
  &p,Bfreq2,Bfreq2_p,Bfreq2_p_p,Blagrangemultiplier,Blagrangemultiplier_p,Blagrangemultiplier_p_p,Blagrangemultiplier_rhs,Box1X_Vol&
  &Int,Box1X_VolInt1,Box1X_VolInt2,Box1Y_VolInt,Box1Y_VolInt1,Box1Y_VolInt2,Box1Z_VolInt,Box1Z_VolInt1,Box1Z_VolInt2,Box1denom_VolI&
  &nt,Box1denom_VolInt1,Box1denom_VolInt2,Bx,Bx_stagger,Bxl,Bxr,Bxtilde,Bxtilde_p,Bxtilde_p_p,Bxtilde_or_Ax_rhs,By,By_stagger,Byl,B&
  &yr,Bytilde,Bytilde_p,Bytilde_p_p,Bytilde_or_Ay_rhs,Bz,Bz_stagger,Bzl,Bzr,Bztilde,Bztilde_p,Bztilde_p_p,Bztilde_or_Az_rhs,CoM_Vol&
  &Int_denominator,CoMx_VolInt,CoMy_VolInt,CoMz_VolInt,E_rad,E_rad_nue,E_rad_nuel,E_rad_nuer,E_rad_nux,E_rad_nuxl,E_rad_nuxr,E_radl&
  &,E_radr,Ex,Ey,Ez,F_rad0,F_rad0_nue,F_rad0_nux,F_rad_scalar,F_rad_scalar_nue,F_rad_scalar_nux,F_radx,F_radx_nue,F_radx_nuele,F_ra&
  &dx_nuer,F_radx_nux,F_radx_nuxle,F_radx_nuxr,F_radxle,F_radxr,F_rady,F_rady_nue,F_rady_nuele,F_rady_nuer,F_rady_nux,F_rady_nuxle,&
  &F_rady_nuxr,F_radyle,F_radyr,F_radz,F_radz_nue,F_radz_nuele,F_radz_nuer,F_radz_nux,F_radz_nuxle,F_radz_nuxr,F_radzle,F_radzr,FaF&
  &al,FaFal_nue,FaFal_nux,FaFar,FaFar_nue,FaFar_nux,Fontfix_tracker_gf,Gammax,Gammax_p,Gammax_p_p,Gammax_drive,Gammax_rhs,Gammaxxx,&
  &Gammaxxy,Gammaxxz,Gammaxyy,Gammaxyz,Gammaxzz,Gammay,Gammay_p,Gammay_p_p,Gammay_drive,Gammay_rhs,Gammayxx,Gammayxy,Gammayxz,Gamma&
  &yyy,Gammayyz,Gammayzz,Gammaz,Gammaz_p,Gammaz_p_p,Gammaz_drive,Gammaz_rhs,Gammazxx,Gammazxy,Gammazxz,Gammazyy,Gammazyz,Gammazzz,G&
  &amx_const_VolInt,Gamy_const_VolInt,Gamz_const_VolInt,Ham_const_VolIntD,Ham_const_VolIntN,Ham_const_excised_VolIntD,Ham_const_exc&
  &ised_VolIntN,Ham_const_excised_innerregion_VolIntD,Ham_const_excised_innerregion_VolIntN,J_ADM_VolInt,J_ADM_VolInt_inner,J_ADM_V&
  &olInt_inner2,Jz_constraint,Jz_constraint_inner,KK,K_pp_rot1,K_pp_rot2,K_pp_rot3,K_rp_rot1,K_rp_rot2,K_rp_rot3,K_rr_rot1,K_rr_rot&
  &2,K_rr_rot3,K_rth_rot1,K_rth_rot2,K_rth_rot3,K_thp_rot1,K_thp_rot2,K_thp_rot3,K_thth_rot1,K_thth_rot2,K_thth_rot3,M0_AH_VolInt,M&
  &0_VolInt,M0_escape100M,M0_escape30M,M0_escape50M,M0_escape70M,M0_horiz_VolInt,M0_r1_VolInt,M0_r2_VolInt,M0_r3_VolInt,MNorm,MONOP&
  &OLE,MRsx,MRsy,MRsz,M_ADM_VolInt,M_constraint,M_constraint_inner,Minternal_VolInt,Minternal_cold_VolInt,NPphi0im,NPphi0im_p,NPphi&
  &0im_p_p,NPphi0re,NPphi0re_p,NPphi0re_p_p,NPphi2im,NPphi2im_p,NPphi2im_p_p,NPphi2re,NPphi2re_p,NPphi2re_p_p,OS_surf_rad,P,P_cld,P&
  &_constraint,P_constraint_inner,P_floor,P_radxx,P_radxx_nue,P_radxx_nux,P_radxy,P_radxy_nue,P_radxy_nux,P_radxz,P_radxz_nue,P_rad&
  &xz_nux,P_radyy,P_radyy_nue,P_radyy_nux,P_radyz,P_radyz_nue,P_radyz_nux,P_radzz,P_radzz_nue,P_radzz_nux,P_tab,P_thermal,PhysicalR&
  &adius,Pl,Pr,PsiNorm,PsiRes,PsiTau,Pvec,RadiusDerivative,RadiusDerivative2,Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,S,SPx,SPx_p,SPx_p_p,SPy,SPy_p,&
  &SPy_p_p,SPz,SPz_p,SPz_p_p,S_rad_x,S_rad_x_p,S_rad_x_p_p,S_rad_x_nue,S_rad_x_nue_p,S_rad_x_nue_p_p,S_rad_x_nue_rhs,S_rad_x_nux,S_&
  &rad_x_nux_p,S_rad_x_nux_p_p,S_rad_x_nux_rhs,S_rad_x_rhs,S_rad_y,S_rad_y_p,S_rad_y_p_p,S_rad_y_nue,S_rad_y_nue_p,S_rad_y_nue_p_p,&
  &S_rad_y_nue_rhs,S_rad_y_nux,S_rad_y_nux_p,S_rad_y_nux_p_p,S_rad_y_nux_rhs,S_rad_y_rhs,S_rad_z,S_rad_z_p,S_rad_z_p_p,S_rad_z_nue,&
  &S_rad_z_nue_p,S_rad_z_nue_p_p,S_rad_z_nue_rhs,S_rad_z_nux,S_rad_z_nux_p,S_rad_z_nux_p_p,S_rad_z_nux_rhs,S_rad_z_rhs,Sx,Sxx,Sxy,S&
  &xz,Sy,Syy,Syz,Sz,Szz,T_0x,T_0x_p,T_0x_p_p,T_0y,T_0y_p,T_0y_p_p,T_0z,T_0z_p,T_0z_p_p,T_VolInt,T_fluid,T_fluidle,T_fluidr,Tem0_0_V&
  &olInt,Tem0_0_outsideBH_VolInt,Tfluid0_0_VolInt,Tfluid0_0_outsideBH_VolInt,VolIntegrand,VolIntegrand2,VolIntegrand3,VolIntegrand4&
  &,Y_e,Y_el,Y_er,ah_centroid_iteration,ah_centroid_iteration_p,ah_centroid_t,ah_centroid_t_p,ah_centroid_valid,ah_centroid_valid_p&
  &,ah_centroid_x,ah_centroid_x_p,ah_centroid_y,ah_centroid_y_p,ah_centroid_z,ah_centroid_z_p,ah_found_flag,ah_initial_find_flag,ah&
  &_origin_x,ah_origin_y,ah_origin_z,ah_radius,ah_really_initial_find_flag,ah_search_flag,b_phi_VolInt,bh_posn_x,bh_posn_y,bh_posn_&
  &z,bh_radius_z,brem_qee_VolInt,brem_qei_VolInt,chi,chi_p,chi_p_p,chi_rad,chi_rad_nue,chi_rad_nux,chi_rhs,cmax,cmax_rad,cmax_rad_n&
  &ue,cmax_rad_nux,cmin,cmin_rad,cmin_rad_nue,cmin_rad_nux,dBx_m,dBy_m,dBz_m,dP_m,dT_fluid_m,ddrbrvec,density_modes_i1,density_mode&
  &s_i2,density_modes_i3,density_modes_i4,density_modes_i5,density_modes_i6,density_modes_r0,density_modes_r1,density_modes_r2,dens&
  &ity_modes_r3,density_modes_r4,density_modes_r5,density_modes_r6,div_beta,drbrvec,drhoYe_m,drhoYe_m_x,drhoYe_m_xp1,drho_b_m,drho_&
  &b_m_x,drho_b_m_xp1,dvx_m,dvx_m_x,dvx_m_xp1,dvy_m,dvz_m,em_J_between_VolInt,em_energy2_VolInt,em_energy_VolInt,em_energy_between_&
  &VolInt,em_energy_outsideBH_VolInt,em_energy_outsideradius1_VolInt,em_energy_outsideradius2_VolInt,emission_gf,emission_gf_nue,em&
  &ission_gf_nux,eps_cld,eps_tab,eps_thermal,eps_tot,eta_nue,excision_zone_gf,fluid_J_horiz_VolInt,fluid_J_r1_VolInt,fluid_J_r2_Vol&
  &Int,fluid_J_r3_VolInt,fluid_energy_horiz_VolInt,fluid_energy_r1_VolInt,fluid_energy_r2_VolInt,fluid_energy_r3_VolInt,gamma_tab,g&
  &conx,gcony,gconz,gupxx,gupxx_f,gupxy,gupxy_f,gupxz,gupxz_f,gupyy,gupyy_f,gupyz,gupyz_f,gupzz,gupzz_f,gxx,gxx_p,gxx_p_p,gxx_f,gxx&
  &_rhs,gxxx,gxxy,gxxz,gxy,gxy_p,gxy_p_p,gxy_f,gxy_rhs,gxyx,gxyy,gxyz,gxz,gxz_p,gxz_p_p,gxz_f,gxz_rhs,gxzx,gxzy,gxzz,gyy,gyy_p,gyy_&
  &p_p,gyy_f,gyy_rhs,gyyx,gyyy,gyyz,gyz,gyz_p,gyz_p_p,gyz_f,gyz_rhs,gyzx,gyzy,gyzz,gzz,gzz_p,gzz_p_p,gzz_f,gzz_rhs,gzzx,gzzy,gzzz,h&
  &,h_p,half_b2_u0_VolInt,half_b2_u0_outsideBH_VolInt,have_global_bdry_max,have_global_bdry_min,int_momentsi,int_momentsr,k_tab,ka_&
  &gf,ka_gf_nue,ka_gf_nux,ks_gf,ks_gf_nue,ks_gf_nux,kxx,kxy,kxz,kyy,kyz,kzz,lapm1,lapm1_p,lapm1_p_p,lapm1_f,lapm1_rhs,lapm1_rot1,la&
  &pm1_rot2,lapm1_rot3,lapset,lapset_p,lapset_p_p,lapset_rhs,lapset_timederiv,lapsex,lapsey,lapsez,mass_star1,mass_star2,mhd_chi_li&
  &ne,mhd_chi_line_p,mhd_chi_line_p_p,mhd_chi_line_rhs,mhd_psi_line,mhd_psi_line_p,mhd_psi_line_p_p,mhd_psi_line_rhs,mhd_st_x,mhd_s&
  &t_x_p,mhd_st_x_p_p,mhd_st_x_rhs,mhd_st_y,mhd_st_y_p,mhd_st_y_p_p,mhd_st_y_rhs,mhd_st_z,mhd_st_z_p,mhd_st_z_p_p,mhd_st_z_rhs,mhd_&
  &u_chi,mhd_u_chi_p,mhd_u_chi_p_p,mhd_u_chi_rhs,mhd_u_psi,mhd_u_psi_p,mhd_u_psi_p_p,mhd_u_psi_rhs,minternal_cold_horiz_VolInt,mint&
  &ernal_cold_r1_VolInt,minternal_cold_r2_VolInt,minternal_cold_r3_VolInt,minternal_horiz_VolInt,minternal_r1_VolInt,minternal_r2_V&
  &olInt,minternal_r3_VolInt,mom_const_VolIntD,mom_const_excised_VolIntD,momentsi,momentsi_old,momentsr,momentsr_old,momx_const_Vol&
  &IntN,momx_const_excised_VolIntN,momy_const_VolIntN,momy_const_excised_VolIntN,momz_const_VolIntN,momz_const_excised_VolIntN,mono&
  &pole_VolInt,monopole_outsideBH_VolInt,mu_nu,n_poly,odd_momentsi,odd_momentsi_old,odd_momentsr,odd_momentsr_old,optd,optd_x,optd_&
  &y,optd_z,pfloor_gf,phi,phi_p,phi_p_p,phi_f,phi_rhs,phi_rot1,phi_rot2,phi_rot3,phix,phiy,phiz,psi,psi0im,psi0im_p,psi0im_p_p,psi0&
  &re,psi0re_p,psi0re_p_p,psi6phi,psi6phi_p,psi6phi_p_p,psi6phi_rhs,rad_const,rad_energy_VolInt,rad_energy_nue_VolInt,rad_energy_nu&
  &x_VolInt,rbrvec,refbd,refbd_p,refbd_p_p,rho,rhoYe,rhoYe_p,rhoYe_p_p,rhoYe_rhs,rho_b,rho_b_atm_gf,rho_bl,rho_br,rho_star,rho_star&
  &_p,rho_star_p_p,rho_star_rhs,rho_tab,rhob_floor,rhovec,sbt,sbx,sby,sbz,shiftp_rot1,shiftp_rot2,shiftp_rot3,shiftr_rot1,shiftr_ro&
  &t2,shiftr_rot3,shiftth_rot1,shiftth_rot2,shiftth_rot3,shiftx,shiftx_p,shiftx_p_p,shiftx_f,shiftx_rhs,shiftxt,shiftxt_p,shiftxt_p&
  &_p,shiftxt_rhs,shiftxt_timederiv,shifty,shifty_p,shifty_p_p,shifty_f,shifty_rhs,shiftyt,shiftyt_p,shiftyt_p_p,shiftyt_rhs,shifty&
  &t_timederiv,shiftz,shiftz_p,shiftz_p_p,shiftz_f,shiftz_rhs,shiftzt,shiftzt_p,shiftzt_p_p,shiftzt_rhs,shiftzt_timederiv,smallb2,s&
  &t_x,st_y,st_z,surfvec,synch_VolInt,tau,tau_p,tau_p_p,tau_rad,tau_rad_p,tau_rad_p_p,tau_rad_nue,tau_rad_nue_p,tau_rad_nue_p_p,tau&
  &_rad_nue_rhs,tau_rad_nux,tau_rad_nux_p,tau_rad_nux_p_p,tau_rad_nux_rhs,tau_rad_rhs,tau_rhs,temp1,temp10,temp11,temp12,temp13,tem&
  &p14,temp15,temp16,temp17,temp18,temp19,temp2,temp20,temp21,temp22,temp23,temp24,temp25,temp26,temp3,temp4,temp5,temp6,temp7,temp&
  &8,temp9,temp_g00,trK,trK_p,trK_p_p,trK_rhs,trR,trRtilde,u0,v02_rad_nuel,v02_rad_nuer,v02_rad_nuxl,v02_rad_nuxr,v02_radl,v02_radr&
  &,v02l,v02r,vvec,vx,vxl,vxr,vy,vyl,vyr,vz,vzl,vzr,w,zeta_rad,zeta_rad_nue,zeta_rad_nux, X0coordinates,X0mask,X0radii,X0sf_radius,&
  &X0space_mask_group,X1coordinates,X1mask,X1sf_radius,X1space_mask_group,X2coordinates,X2mask,X2space_mask_group,active_length,num&
  &_levels_length,positions_length,radii_length,sf_active_length,sf_coordinate_descriptors_length,sf_info_length,sf_maxreflevel_len&
  &gth,sf_minreflevel_length,sf_origin_length,sf_radius_length,sf_shape_descriptors_length,sf_valid_length,active,coarse_dx,coarse_&
  &dy,coarse_dz,emask,num_levels,position_x,position_y,position_z,r,radius,sf_active,sf_area,sf_centroid_x,sf_centroid_y,sf_centroi&
  &d_z,sf_delta_phi,sf_delta_theta,sf_max_radius,sf_max_x,sf_max_y,sf_max_z,sf_maxreflevel,sf_mean_radius,sf_min_radius,sf_min_x,sf&
  &_min_y,sf_min_z,sf_minreflevel,sf_nghostsphi,sf_nghoststheta,sf_nphi,sf_ntheta,sf_origin_phi,sf_origin_theta,sf_origin_x,sf_orig&
  &in_y,sf_origin_z,sf_quadrupole_xx,sf_quadrupole_xy,sf_quadrupole_xz,sf_quadrupole_yy,sf_quadrupole_yz,sf_quadrupole_zz,sf_radius&
  &,sf_valid,space_mask,x,y,z)
    implicit none
    INTEGER cctk_dim
  INTEGER cctk_gsh(cctk_dim),cctk_lsh(cctk_dim)
  INTEGER cctk_lbnd(cctk_dim),cctk_ubnd(cctk_dim)
  INTEGER cctk_lssh(3*cctk_dim)
  INTEGER cctk_from(cctk_dim),cctk_to(cctk_dim)
  INTEGER cctk_bbox(2*cctk_dim)
  REAL*8 cctk_delta_time, cctk_time
  REAL*8 cctk_delta_space(cctk_dim)
  REAL*8 cctk_origin_space(cctk_dim)
  INTEGER cctk_levfac(cctk_dim)
  INTEGER cctk_levoff(cctk_dim)
  INTEGER cctk_levoffdenom(cctk_dim)
  INTEGER cctk_timefac
  INTEGER cctk_convlevel
  INTEGER cctk_convfac
  INTEGER cctk_nghostzones(cctk_dim)
  INTEGER cctk_iteration
  integer*8 cctkGH
  INTEGER X0Axyzintrgroup
  INTEGER X0Tracer_rhs_group
  INTEGER X0Tracergroup
  INTEGER X0last_ah_data
  INTEGER X1Axyzintrgroup
  INTEGER X1Tracer_rhs_group
  INTEGER X2Axyzintrgroup
  REAL*8 Axintr(X0Axyzintrgroup,X1Axyzintrgroup,X2Axyzintrgroup)
  REAL*8 Ayintr(X0Axyzintrgroup,X1Axyzintrgroup,X2Axyzintrgroup)
  REAL*8 Azintr(X0Axyzintrgroup,X1Axyzintrgroup,X2Axyzintrgroup)
  REAL*8 M0dot_last
  REAL*8 Pintr(X0Axyzintrgroup,X1Axyzintrgroup,X2Axyzintrgroup)
  REAL*8 ah_radii_last(X0last_ah_data)
  INTEGER*4 bhnsEthermalNOW
  INTEGER*4 bhns_rk4_particle_tracer_step
  INTEGER*4 bhnsinsertBNOW
  INTEGER*4 bhnsinsertRADNOW
  REAL*8 int_M0dot
  REAL*8 local_density(X0Tracergroup)
  REAL*8 m1_minus_u_lower_zero(X0Tracergroup)
  REAL*8 rhobintr(X0Axyzintrgroup,X1Axyzintrgroup,X2Axyzintrgroup)
  REAL*8 t_last
  REAL*8 taul_rhs(X0Tracer_rhs_group,X1Tracer_rhs_group)
  REAL*8 tracer_taul(X0Tracergroup)
  REAL*8 tracer_x(X0Tracergroup)
  REAL*8 tracer_xc(X0Tracergroup)
  REAL*8 tracer_y(X0Tracergroup)
  REAL*8 tracer_yc(X0Tracergroup)
  REAL*8 tracer_z(X0Tracergroup)
  REAL*8 u_x_local(X0Tracergroup)
  REAL*8 u_y_local(X0Tracergroup)
  REAL*8 x_rhs(X0Tracer_rhs_group,X1Tracer_rhs_group)
  REAL*8 xc_rhs(X0Tracer_rhs_group,X1Tracer_rhs_group)
  REAL*8 xh_last
  REAL*8 y_rhs(X0Tracer_rhs_group,X1Tracer_rhs_group)
  REAL*8 yc_rhs(X0Tracer_rhs_group,X1Tracer_rhs_group)
  REAL*8 yh_last
  REAL*8 z_rhs(X0Tracer_rhs_group,X1Tracer_rhs_group)
  REAL*8 zh_last
  INTEGER X0BSSN_AH
  INTEGER X0BSSN_aux_private
  INTEGER X0BSSN_aux_restrict2
  INTEGER X0BSSN_diag_restrict
  INTEGER X0BSSN_gupij
  INTEGER X0BSSN_matter
  INTEGER X0BSSN_refbd
  INTEGER X0BSSN_rhs
  INTEGER X0BSSN_vars
  INTEGER X0NPphi2_em
  INTEGER X0OmegF
  INTEGER X0Psi4_gw
  INTEGER X0Stagger_Bs
  INTEGER X0ah_radius
  INTEGER X0bh_posns
  INTEGER X0diag_vacuum_gf
  INTEGER X0disk_atmosphere
  INTEGER X0em_Ax
  INTEGER X0em_Ay
  INTEGER X0em_Az
  INTEGER X0em_Blagrangemultiplier
  INTEGER X0em_Blagrangemultiplier_rhs
  INTEGER X0em_Phi
  INTEGER X0em_Phi_rhs
  INTEGER X0em_conservativex
  INTEGER X0em_conservativey
  INTEGER X0em_conservativez
  INTEGER X0em_rhsx
  INTEGER X0em_rhsy
  INTEGER X0em_rhsz
  INTEGER X0eos_params1
  INTEGER X0eos_params2
  INTEGER X0excision_int_gfs
  INTEGER X0field_line_variables
  INTEGER X0field_line_variables_rhs
  INTEGER X0fisheye_vars
  INTEGER X0gw_moment_arrays
  INTEGER X0have_global_bdry
  INTEGER X0lapse_derivatives
  INTEGER X0lapse_rhs
  INTEGER X0lapse_vars
  INTEGER X0lapse_vars_aux
  INTEGER X0metric_spher_pol_1
  INTEGER X0metric_spher_pol_2
  INTEGER X0metric_spher_pol_3
  INTEGER X0mhd_conservatives
  INTEGER X0mhd_nosync
  INTEGER X0mhd_primitives
  INTEGER X0mhd_rhs
  INTEGER X0mhd_sync_lr_B_quantities
  INTEGER X0mhd_sync_lr_hydro_quantities
  INTEGER X0mhd_sync_metric_facevals
  INTEGER X0mhd_sync_nablas
  INTEGER X0mhd_sync_nablas_diag
  INTEGER X0mhd_sync_rho_br_rho_bl
  INTEGER X0mhd_temps
  INTEGER X0mhd_vs
  INTEGER X0micphys_conservatives
  INTEGER X0micphys_conservatives_rhs
  INTEGER X0micphys_sync_lr_hydro_quantities
  INTEGER X0micphys_sync_nablas
  INTEGER X0microphys_primitives
  INTEGER X0phi_derivs
  INTEGER X0rad_conservatives
  INTEGER X0rad_conservatives_rhs
  INTEGER X0rad_pressure
  INTEGER X0rad_primitives
  INTEGER X0rad_sync_lr
  INTEGER X0rhovecs
  INTEGER X0shift_rhs
  INTEGER X0shift_vars
  INTEGER X0shift_vars_temp
  INTEGER X0surf_params
  INTEGER X1BSSN_AH
  INTEGER X1BSSN_aux_private
  INTEGER X1BSSN_aux_restrict2
  INTEGER X1BSSN_diag_restrict
  INTEGER X1BSSN_gupij
  INTEGER X1BSSN_matter
  INTEGER X1BSSN_refbd
  INTEGER X1BSSN_rhs
  INTEGER X1BSSN_vars
  INTEGER X1NPphi2_em
  INTEGER X1OmegF
  INTEGER X1Psi4_gw
  INTEGER X1Stagger_Bs
  INTEGER X1ah_radius
  INTEGER X1diag_vacuum_gf
  INTEGER X1disk_atmosphere
  INTEGER X1em_Ax
  INTEGER X1em_Ay
  INTEGER X1em_Az
  INTEGER X1em_Blagrangemultiplier
  INTEGER X1em_Blagrangemultiplier_rhs
  INTEGER X1em_Phi
  INTEGER X1em_Phi_rhs
  INTEGER X1em_conservativex
  INTEGER X1em_conservativey
  INTEGER X1em_conservativez
  INTEGER X1em_rhsx
  INTEGER X1em_rhsy
  INTEGER X1em_rhsz
  INTEGER X1excision_int_gfs
  INTEGER X1field_line_variables
  INTEGER X1field_line_variables_rhs
  INTEGER X1fisheye_vars
  INTEGER X1have_global_bdry
  INTEGER X1lapse_derivatives
  INTEGER X1lapse_rhs
  INTEGER X1lapse_vars
  INTEGER X1lapse_vars_aux
  INTEGER X1mhd_conservatives
  INTEGER X1mhd_nosync
  INTEGER X1mhd_primitives
  INTEGER X1mhd_rhs
  INTEGER X1mhd_sync_lr_B_quantities
  INTEGER X1mhd_sync_lr_hydro_quantities
  INTEGER X1mhd_sync_metric_facevals
  INTEGER X1mhd_sync_nablas
  INTEGER X1mhd_sync_nablas_diag
  INTEGER X1mhd_sync_rho_br_rho_bl
  INTEGER X1mhd_temps
  INTEGER X1mhd_vs
  INTEGER X1micphys_conservatives
  INTEGER X1micphys_conservatives_rhs
  INTEGER X1micphys_sync_lr_hydro_quantities
  INTEGER X1micphys_sync_nablas
  INTEGER X1microphys_primitives
  INTEGER X1phi_derivs
  INTEGER X1rad_conservatives
  INTEGER X1rad_conservatives_rhs
  INTEGER X1rad_pressure
  INTEGER X1rad_primitives
  INTEGER X1rad_sync_lr
  INTEGER X1shift_rhs
  INTEGER X1shift_vars
  INTEGER X1shift_vars_temp
  INTEGER X2BSSN_AH
  INTEGER X2BSSN_aux_private
  INTEGER X2BSSN_aux_restrict2
  INTEGER X2BSSN_diag_restrict
  INTEGER X2BSSN_gupij
  INTEGER X2BSSN_matter
  INTEGER X2BSSN_refbd
  INTEGER X2BSSN_rhs
  INTEGER X2BSSN_vars
  INTEGER X2NPphi2_em
  INTEGER X2OmegF
  INTEGER X2Psi4_gw
  INTEGER X2Stagger_Bs
  INTEGER X2ah_radius
  INTEGER X2diag_vacuum_gf
  INTEGER X2disk_atmosphere
  INTEGER X2em_Ax
  INTEGER X2em_Ay
  INTEGER X2em_Az
  INTEGER X2em_Blagrangemultiplier
  INTEGER X2em_Blagrangemultiplier_rhs
  INTEGER X2em_Phi
  INTEGER X2em_Phi_rhs
  INTEGER X2em_conservativex
  INTEGER X2em_conservativey
  INTEGER X2em_conservativez
  INTEGER X2em_rhsx
  INTEGER X2em_rhsy
  INTEGER X2em_rhsz
  INTEGER X2excision_int_gfs
  INTEGER X2field_line_variables
  INTEGER X2field_line_variables_rhs
  INTEGER X2fisheye_vars
  INTEGER X2lapse_derivatives
  INTEGER X2lapse_rhs
  INTEGER X2lapse_vars
  INTEGER X2lapse_vars_aux
  INTEGER X2mhd_conservatives
  INTEGER X2mhd_nosync
  INTEGER X2mhd_primitives
  INTEGER X2mhd_rhs
  INTEGER X2mhd_sync_lr_B_quantities
  INTEGER X2mhd_sync_lr_hydro_quantities
  INTEGER X2mhd_sync_metric_facevals
  INTEGER X2mhd_sync_nablas
  INTEGER X2mhd_sync_nablas_diag
  INTEGER X2mhd_sync_rho_br_rho_bl
  INTEGER X2mhd_temps
  INTEGER X2mhd_vs
  INTEGER X2micphys_conservatives
  INTEGER X2micphys_conservatives_rhs
  INTEGER X2micphys_sync_lr_hydro_quantities
  INTEGER X2micphys_sync_nablas
  INTEGER X2microphys_primitives
  INTEGER X2phi_derivs
  INTEGER X2rad_conservatives
  INTEGER X2rad_conservatives_rhs
  INTEGER X2rad_pressure
  INTEGER X2rad_primitives
  INTEGER X2rad_sync_lr
  INTEGER X2shift_rhs
  INTEGER X2shift_vars
  INTEGER X2shift_vars_temp
  INTEGER ah_centroid_length
  INTEGER ah_flags_length
  INTEGER ah_origin_length
  INTEGER ah_radius_length
  REAL*8 Aupxx(X0BSSN_aux_private,X1BSSN_aux_private,X2BSSN_aux_private)
  REAL*8 Aupxy(X0BSSN_aux_private,X1BSSN_aux_private,X2BSSN_aux_private)
  REAL*8 Aupxz(X0BSSN_aux_private,X1BSSN_aux_private,X2BSSN_aux_private)
  REAL*8 Aupyy(X0BSSN_aux_private,X1BSSN_aux_private,X2BSSN_aux_private)
  REAL*8 Aupyz(X0BSSN_aux_private,X1BSSN_aux_private,X2BSSN_aux_private)
  REAL*8 Aupzz(X0BSSN_aux_private,X1BSSN_aux_private,X2BSSN_aux_private)
  REAL*8 Ax(X0em_Ax,X1em_Ax,X2em_Ax)
  REAL*8 Ax_p(X0em_Ax,X1em_Ax,X2em_Ax)
  REAL*8 Ax_p_p(X0em_Ax,X1em_Ax,X2em_Ax)
  REAL*8 Axx(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Axx_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Axx_p_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Axx_rhs(X0BSSN_rhs,X1BSSN_rhs,X2BSSN_rhs)
  REAL*8 Axy(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Axy_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Axy_p_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Axy_rhs(X0BSSN_rhs,X1BSSN_rhs,X2BSSN_rhs)
  REAL*8 Axz(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Axz_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Axz_p_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Axz_rhs(X0BSSN_rhs,X1BSSN_rhs,X2BSSN_rhs)
  REAL*8 Ay(X0em_Ay,X1em_Ay,X2em_Ay)
  REAL*8 Ay_p(X0em_Ay,X1em_Ay,X2em_Ay)
  REAL*8 Ay_p_p(X0em_Ay,X1em_Ay,X2em_Ay)
  REAL*8 Ayy(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Ayy_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Ayy_p_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Ayy_rhs(X0BSSN_rhs,X1BSSN_rhs,X2BSSN_rhs)
  REAL*8 Ayz(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Ayz_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Ayz_p_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Ayz_rhs(X0BSSN_rhs,X1BSSN_rhs,X2BSSN_rhs)
  REAL*8 Az(X0em_Az,X1em_Az,X2em_Az)
  REAL*8 Az_p(X0em_Az,X1em_Az,X2em_Az)
  REAL*8 Az_p_p(X0em_Az,X1em_Az,X2em_Az)
  REAL*8 Azz(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Azz_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Azz_p_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Azz_rhs(X0BSSN_rhs,X1BSSN_rhs,X2BSSN_rhs)
  REAL*8 Bfreq1(X0OmegF,X1OmegF,X2OmegF)
  REAL*8 Bfreq1_p(X0OmegF,X1OmegF,X2OmegF)
  REAL*8 Bfreq1_p_p(X0OmegF,X1OmegF,X2OmegF)
  REAL*8 Bfreq2(X0OmegF,X1OmegF,X2OmegF)
  REAL*8 Bfreq2_p(X0OmegF,X1OmegF,X2OmegF)
  REAL*8 Bfreq2_p_p(X0OmegF,X1OmegF,X2OmegF)
  REAL*8 Blagrangemultiplier(X0em_Blagrangemultiplier,X1em_Blagrangemultiplier,X2em_Blagrangemultiplier)
  REAL*8 Blagrangemultiplier_p(X0em_Blagrangemultiplier,X1em_Blagrangemultiplier,X2em_Blagrangemultiplier)
  REAL*8 Blagrangemultiplier_p_p(X0em_Blagrangemultiplier,X1em_Blagrangemultiplier,X2em_Blagrangemultiplier)
  REAL*8 Blagrangemultiplier_rhs(X0em_Blagrangemultiplier_rhs,X1em_Blagrangemultiplier_rhs,X2em_Blagrangemultiplier_rhs)
  REAL*8 Box1X_VolInt
  REAL*8 Box1X_VolInt1
  REAL*8 Box1X_VolInt2
  REAL*8 Box1Y_VolInt
  REAL*8 Box1Y_VolInt1
  REAL*8 Box1Y_VolInt2
  REAL*8 Box1Z_VolInt
  REAL*8 Box1Z_VolInt1
  REAL*8 Box1Z_VolInt2
  REAL*8 Box1denom_VolInt
  REAL*8 Box1denom_VolInt1
  REAL*8 Box1denom_VolInt2
  REAL*8 Bx(X0mhd_primitives,X1mhd_primitives,X2mhd_primitives)
  REAL*8 Bx_stagger(X0Stagger_Bs,X1Stagger_Bs,X2Stagger_Bs)
  REAL*8 Bxl(X0mhd_sync_lr_B_quantities,X1mhd_sync_lr_B_quantities,X2mhd_sync_lr_B_quantities)
  REAL*8 Bxr(X0mhd_sync_lr_B_quantities,X1mhd_sync_lr_B_quantities,X2mhd_sync_lr_B_quantities)
  REAL*8 Bxtilde(X0em_conservativex,X1em_conservativex,X2em_conservativex)
  REAL*8 Bxtilde_p(X0em_conservativex,X1em_conservativex,X2em_conservativex)
  REAL*8 Bxtilde_p_p(X0em_conservativex,X1em_conservativex,X2em_conservativex)
  REAL*8 Bxtilde_or_Ax_rhs(X0em_rhsx,X1em_rhsx,X2em_rhsx)
  REAL*8 By(X0mhd_primitives,X1mhd_primitives,X2mhd_primitives)
  REAL*8 By_stagger(X0Stagger_Bs,X1Stagger_Bs,X2Stagger_Bs)
  REAL*8 Byl(X0mhd_sync_lr_B_quantities,X1mhd_sync_lr_B_quantities,X2mhd_sync_lr_B_quantities)
  REAL*8 Byr(X0mhd_sync_lr_B_quantities,X1mhd_sync_lr_B_quantities,X2mhd_sync_lr_B_quantities)
  REAL*8 Bytilde(X0em_conservativey,X1em_conservativey,X2em_conservativey)
  REAL*8 Bytilde_p(X0em_conservativey,X1em_conservativey,X2em_conservativey)
  REAL*8 Bytilde_p_p(X0em_conservativey,X1em_conservativey,X2em_conservativey)
  REAL*8 Bytilde_or_Ay_rhs(X0em_rhsy,X1em_rhsy,X2em_rhsy)
  REAL*8 Bz(X0mhd_primitives,X1mhd_primitives,X2mhd_primitives)
  REAL*8 Bz_stagger(X0Stagger_Bs,X1Stagger_Bs,X2Stagger_Bs)
  REAL*8 Bzl(X0mhd_sync_lr_B_quantities,X1mhd_sync_lr_B_quantities,X2mhd_sync_lr_B_quantities)
  REAL*8 Bzr(X0mhd_sync_lr_B_quantities,X1mhd_sync_lr_B_quantities,X2mhd_sync_lr_B_quantities)
  REAL*8 Bztilde(X0em_conservativez,X1em_conservativez,X2em_conservativez)
  REAL*8 Bztilde_p(X0em_conservativez,X1em_conservativez,X2em_conservativez)
  REAL*8 Bztilde_p_p(X0em_conservativez,X1em_conservativez,X2em_conservativez)
  REAL*8 Bztilde_or_Az_rhs(X0em_rhsz,X1em_rhsz,X2em_rhsz)
  REAL*8 CoM_VolInt_denominator
  REAL*8 CoMx_VolInt
  REAL*8 CoMy_VolInt
  REAL*8 CoMz_VolInt
  REAL*8 E_rad(X0rad_primitives,X1rad_primitives,X2rad_primitives)
  REAL*8 E_rad_nue(X0rad_primitives,X1rad_primitives,X2rad_primitives)
  REAL*8 E_rad_nuel(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 E_rad_nuer(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 E_rad_nux(X0rad_primitives,X1rad_primitives,X2rad_primitives)
  REAL*8 E_rad_nuxl(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 E_rad_nuxr(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 E_radl(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 E_radr(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 Ex(X0mhd_primitives,X1mhd_primitives,X2mhd_primitives)
  REAL*8 Ey(X0mhd_primitives,X1mhd_primitives,X2mhd_primitives)
  REAL*8 Ez(X0mhd_primitives,X1mhd_primitives,X2mhd_primitives)
  REAL*8 F_rad0(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 F_rad0_nue(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 F_rad0_nux(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 F_rad_scalar(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 F_rad_scalar_nue(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 F_rad_scalar_nux(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 F_radx(X0rad_primitives,X1rad_primitives,X2rad_primitives)
  REAL*8 F_radx_nue(X0rad_primitives,X1rad_primitives,X2rad_primitives)
  REAL*8 F_radx_nuele(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 F_radx_nuer(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 F_radx_nux(X0rad_primitives,X1rad_primitives,X2rad_primitives)
  REAL*8 F_radx_nuxle(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 F_radx_nuxr(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 F_radxle(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 F_radxr(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 F_rady(X0rad_primitives,X1rad_primitives,X2rad_primitives)
  REAL*8 F_rady_nue(X0rad_primitives,X1rad_primitives,X2rad_primitives)
  REAL*8 F_rady_nuele(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 F_rady_nuer(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 F_rady_nux(X0rad_primitives,X1rad_primitives,X2rad_primitives)
  REAL*8 F_rady_nuxle(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 F_rady_nuxr(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 F_radyle(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 F_radyr(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 F_radz(X0rad_primitives,X1rad_primitives,X2rad_primitives)
  REAL*8 F_radz_nue(X0rad_primitives,X1rad_primitives,X2rad_primitives)
  REAL*8 F_radz_nuele(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 F_radz_nuer(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 F_radz_nux(X0rad_primitives,X1rad_primitives,X2rad_primitives)
  REAL*8 F_radz_nuxle(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 F_radz_nuxr(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 F_radzle(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 F_radzr(X0rad_sync_lr,X1rad_sync_lr,X2rad_sync_lr)
  REAL*8 FaFal(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 FaFal_nue(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 FaFal_nux(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 FaFar(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 FaFar_nue(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 FaFar_nux(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 Fontfix_tracker_gf(X0disk_atmosphere,X1disk_atmosphere,X2disk_atmosphere)
  REAL*8 Gammax(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Gammax_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Gammax_p_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Gammax_drive(X0shift_vars_temp,X1shift_vars_temp,X2shift_vars_temp)
  REAL*8 Gammax_rhs(X0BSSN_rhs,X1BSSN_rhs,X2BSSN_rhs)
  REAL*8 Gammaxxx(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 Gammaxxy(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 Gammaxxz(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 Gammaxyy(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 Gammaxyz(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 Gammaxzz(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 Gammay(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Gammay_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Gammay_p_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Gammay_drive(X0shift_vars_temp,X1shift_vars_temp,X2shift_vars_temp)
  REAL*8 Gammay_rhs(X0BSSN_rhs,X1BSSN_rhs,X2BSSN_rhs)
  REAL*8 Gammayxx(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 Gammayxy(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 Gammayxz(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 Gammayyy(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 Gammayyz(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 Gammayzz(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 Gammaz(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Gammaz_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Gammaz_p_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 Gammaz_drive(X0shift_vars_temp,X1shift_vars_temp,X2shift_vars_temp)
  REAL*8 Gammaz_rhs(X0BSSN_rhs,X1BSSN_rhs,X2BSSN_rhs)
  REAL*8 Gammazxx(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 Gammazxy(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 Gammazxz(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 Gammazyy(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 Gammazyz(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 Gammazzz(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 Gamx_const_VolInt
  REAL*8 Gamy_const_VolInt
  REAL*8 Gamz_const_VolInt
  REAL*8 Ham_const_VolIntD
  REAL*8 Ham_const_VolIntN
  REAL*8 Ham_const_excised_VolIntD
  REAL*8 Ham_const_excised_VolIntN
  REAL*8 Ham_const_excised_innerregion_VolIntD
  REAL*8 Ham_const_excised_innerregion_VolIntN
  REAL*8 J_ADM_VolInt
  REAL*8 J_ADM_VolInt_inner
  REAL*8 J_ADM_VolInt_inner2
  REAL*8 Jz_constraint
  REAL*8 Jz_constraint_inner
  REAL*8 KK(X0BSSN_diag_restrict,X1BSSN_diag_restrict,X2BSSN_diag_restrict)
  REAL*8 K_pp_rot1(X0metric_spher_pol_1)
  REAL*8 K_pp_rot2(X0metric_spher_pol_2)
  REAL*8 K_pp_rot3(X0metric_spher_pol_3)
  REAL*8 K_rp_rot1(X0metric_spher_pol_1)
  REAL*8 K_rp_rot2(X0metric_spher_pol_2)
  REAL*8 K_rp_rot3(X0metric_spher_pol_3)
  REAL*8 K_rr_rot1(X0metric_spher_pol_1)
  REAL*8 K_rr_rot2(X0metric_spher_pol_2)
  REAL*8 K_rr_rot3(X0metric_spher_pol_3)
  REAL*8 K_rth_rot1(X0metric_spher_pol_1)
  REAL*8 K_rth_rot2(X0metric_spher_pol_2)
  REAL*8 K_rth_rot3(X0metric_spher_pol_3)
  REAL*8 K_thp_rot1(X0metric_spher_pol_1)
  REAL*8 K_thp_rot2(X0metric_spher_pol_2)
  REAL*8 K_thp_rot3(X0metric_spher_pol_3)
  REAL*8 K_thth_rot1(X0metric_spher_pol_1)
  REAL*8 K_thth_rot2(X0metric_spher_pol_2)
  REAL*8 K_thth_rot3(X0metric_spher_pol_3)
  REAL*8 M0_AH_VolInt
  REAL*8 M0_VolInt
  REAL*8 M0_escape100M
  REAL*8 M0_escape30M
  REAL*8 M0_escape50M
  REAL*8 M0_escape70M
  REAL*8 M0_horiz_VolInt
  REAL*8 M0_r1_VolInt
  REAL*8 M0_r2_VolInt
  REAL*8 M0_r3_VolInt
  REAL*8 MNorm(X0BSSN_diag_restrict,X1BSSN_diag_restrict,X2BSSN_diag_restrict)
  REAL*8 MONOPOLE(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 MRsx(X0BSSN_diag_restrict,X1BSSN_diag_restrict,X2BSSN_diag_restrict)
  REAL*8 MRsy(X0BSSN_diag_restrict,X1BSSN_diag_restrict,X2BSSN_diag_restrict)
  REAL*8 MRsz(X0BSSN_diag_restrict,X1BSSN_diag_restrict,X2BSSN_diag_restrict)
  REAL*8 M_ADM_VolInt
  REAL*8 M_constraint
  REAL*8 M_constraint_inner
  REAL*8 Minternal_VolInt
  REAL*8 Minternal_cold_VolInt
  REAL*8 NPphi0im(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 NPphi0im_p(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 NPphi0im_p_p(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 NPphi0re(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 NPphi0re_p(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 NPphi0re_p_p(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 NPphi2im(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 NPphi2im_p(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 NPphi2im_p_p(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 NPphi2re(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 NPphi2re_p(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 NPphi2re_p_p(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 OS_surf_rad
  REAL*8 P(X0mhd_primitives,X1mhd_primitives,X2mhd_primitives)
  REAL*8 P_cld(X0microphys_primitives,X1microphys_primitives,X2microphys_primitives)
  REAL*8 P_constraint
  REAL*8 P_constraint_inner
  REAL*8 P_floor(X0mhd_nosync,X1mhd_nosync,X2mhd_nosync)
  REAL*8 P_radxx(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 P_radxx_nue(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 P_radxx_nux(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 P_radxy(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 P_radxy_nue(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 P_radxy_nux(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 P_radxz(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 P_radxz_nue(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 P_radxz_nux(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 P_radyy(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 P_radyy_nue(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 P_radyy_nux(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 P_radyz(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 P_radyz_nue(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 P_radyz_nux(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 P_radzz(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 P_radzz_nue(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 P_radzz_nux(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 P_tab(X0eos_params1)
  REAL*8 P_thermal(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 PhysicalRadius(X0fisheye_vars,X1fisheye_vars,X2fisheye_vars)
  REAL*8 Pl(X0mhd_sync_lr_hydro_quantities,X1mhd_sync_lr_hydro_quantities,X2mhd_sync_lr_hydro_quantities)
  REAL*8 Pr(X0mhd_sync_lr_hydro_quantities,X1mhd_sync_lr_hydro_quantities,X2mhd_sync_lr_hydro_quantities)
  REAL*8 PsiNorm(X0BSSN_diag_restrict,X1BSSN_diag_restrict,X2BSSN_diag_restrict)
  REAL*8 PsiRes(X0BSSN_diag_restrict,X1BSSN_diag_restrict,X2BSSN_diag_restrict)
  REAL*8 PsiTau(X0BSSN_diag_restrict,X1BSSN_diag_restrict,X2BSSN_diag_restrict)
  REAL*8 Pvec(X0rhovecs)
  REAL*8 RadiusDerivative(X0fisheye_vars,X1fisheye_vars,X2fisheye_vars)
  REAL*8 RadiusDerivative2(X0fisheye_vars,X1fisheye_vars,X2fisheye_vars)
  REAL*8 Rxx(X0BSSN_diag_restrict,X1BSSN_diag_restrict,X2BSSN_diag_restrict)
  REAL*8 Rxy(X0BSSN_diag_restrict,X1BSSN_diag_restrict,X2BSSN_diag_restrict)
  REAL*8 Rxz(X0BSSN_diag_restrict,X1BSSN_diag_restrict,X2BSSN_diag_restrict)
  REAL*8 Ryy(X0BSSN_diag_restrict,X1BSSN_diag_restrict,X2BSSN_diag_restrict)
  REAL*8 Ryz(X0BSSN_diag_restrict,X1BSSN_diag_restrict,X2BSSN_diag_restrict)
  REAL*8 Rzz(X0BSSN_diag_restrict,X1BSSN_diag_restrict,X2BSSN_diag_restrict)
  REAL*8 S(X0BSSN_matter,X1BSSN_matter,X2BSSN_matter)
  REAL*8 SPx(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 SPx_p(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 SPx_p_p(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 SPy(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 SPy_p(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 SPy_p_p(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 SPz(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 SPz_p(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 SPz_p_p(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 S_rad_x(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_x_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_x_p_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_x_nue(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_x_nue_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_x_nue_p_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_x_nue_rhs(X0rad_conservatives_rhs,X1rad_conservatives_rhs,X2rad_conservatives_rhs)
  REAL*8 S_rad_x_nux(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_x_nux_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_x_nux_p_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_x_nux_rhs(X0rad_conservatives_rhs,X1rad_conservatives_rhs,X2rad_conservatives_rhs)
  REAL*8 S_rad_x_rhs(X0rad_conservatives_rhs,X1rad_conservatives_rhs,X2rad_conservatives_rhs)
  REAL*8 S_rad_y(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_y_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_y_p_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_y_nue(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_y_nue_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_y_nue_p_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_y_nue_rhs(X0rad_conservatives_rhs,X1rad_conservatives_rhs,X2rad_conservatives_rhs)
  REAL*8 S_rad_y_nux(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_y_nux_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_y_nux_p_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_y_nux_rhs(X0rad_conservatives_rhs,X1rad_conservatives_rhs,X2rad_conservatives_rhs)
  REAL*8 S_rad_y_rhs(X0rad_conservatives_rhs,X1rad_conservatives_rhs,X2rad_conservatives_rhs)
  REAL*8 S_rad_z(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_z_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_z_p_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_z_nue(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_z_nue_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_z_nue_p_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_z_nue_rhs(X0rad_conservatives_rhs,X1rad_conservatives_rhs,X2rad_conservatives_rhs)
  REAL*8 S_rad_z_nux(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_z_nux_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_z_nux_p_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 S_rad_z_nux_rhs(X0rad_conservatives_rhs,X1rad_conservatives_rhs,X2rad_conservatives_rhs)
  REAL*8 S_rad_z_rhs(X0rad_conservatives_rhs,X1rad_conservatives_rhs,X2rad_conservatives_rhs)
  REAL*8 Sx(X0BSSN_matter,X1BSSN_matter,X2BSSN_matter)
  REAL*8 Sxx(X0BSSN_matter,X1BSSN_matter,X2BSSN_matter)
  REAL*8 Sxy(X0BSSN_matter,X1BSSN_matter,X2BSSN_matter)
  REAL*8 Sxz(X0BSSN_matter,X1BSSN_matter,X2BSSN_matter)
  REAL*8 Sy(X0BSSN_matter,X1BSSN_matter,X2BSSN_matter)
  REAL*8 Syy(X0BSSN_matter,X1BSSN_matter,X2BSSN_matter)
  REAL*8 Syz(X0BSSN_matter,X1BSSN_matter,X2BSSN_matter)
  REAL*8 Sz(X0BSSN_matter,X1BSSN_matter,X2BSSN_matter)
  REAL*8 Szz(X0BSSN_matter,X1BSSN_matter,X2BSSN_matter)
  REAL*8 T_0x(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 T_0x_p(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 T_0x_p_p(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 T_0y(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 T_0y_p(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 T_0y_p_p(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 T_0z(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 T_0z_p(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 T_0z_p_p(X0NPphi2_em,X1NPphi2_em,X2NPphi2_em)
  REAL*8 T_VolInt
  REAL*8 T_fluid(X0microphys_primitives,X1microphys_primitives,X2microphys_primitives)
  REAL*8 T_fluidle(X0micphys_sync_lr_hydro_quantities,X1micphys_sync_lr_hydro_quantities,X2micphys_sync_lr_hydro_quantities)
  REAL*8 T_fluidr(X0micphys_sync_lr_hydro_quantities,X1micphys_sync_lr_hydro_quantities,X2micphys_sync_lr_hydro_quantities)
  REAL*8 Tem0_0_VolInt
  REAL*8 Tem0_0_outsideBH_VolInt
  REAL*8 Tfluid0_0_VolInt
  REAL*8 Tfluid0_0_outsideBH_VolInt
  REAL*8 VolIntegrand(X0diag_vacuum_gf,X1diag_vacuum_gf,X2diag_vacuum_gf)
  REAL*8 VolIntegrand2(X0diag_vacuum_gf,X1diag_vacuum_gf,X2diag_vacuum_gf)
  REAL*8 VolIntegrand3(X0diag_vacuum_gf,X1diag_vacuum_gf,X2diag_vacuum_gf)
  REAL*8 VolIntegrand4(X0diag_vacuum_gf,X1diag_vacuum_gf,X2diag_vacuum_gf)
  REAL*8 Y_e(X0microphys_primitives,X1microphys_primitives,X2microphys_primitives)
  REAL*8 Y_el(X0micphys_sync_lr_hydro_quantities,X1micphys_sync_lr_hydro_quantities,X2micphys_sync_lr_hydro_quantities)
  REAL*8 Y_er(X0micphys_sync_lr_hydro_quantities,X1micphys_sync_lr_hydro_quantities,X2micphys_sync_lr_hydro_quantities)
  INTEGER*4 ah_centroid_iteration(ah_flags_length)
  INTEGER*4 ah_centroid_iteration_p(ah_flags_length)
  REAL*8 ah_centroid_t(ah_centroid_length)
  REAL*8 ah_centroid_t_p(ah_centroid_length)
  INTEGER*4 ah_centroid_valid(ah_flags_length)
  INTEGER*4 ah_centroid_valid_p(ah_flags_length)
  REAL*8 ah_centroid_x(ah_centroid_length)
  REAL*8 ah_centroid_x_p(ah_centroid_length)
  REAL*8 ah_centroid_y(ah_centroid_length)
  REAL*8 ah_centroid_y_p(ah_centroid_length)
  REAL*8 ah_centroid_z(ah_centroid_length)
  REAL*8 ah_centroid_z_p(ah_centroid_length)
  INTEGER*4 ah_found_flag(ah_flags_length)
  INTEGER*4 ah_initial_find_flag(ah_flags_length)
  REAL*8 ah_origin_x(ah_origin_length)
  REAL*8 ah_origin_y(ah_origin_length)
  REAL*8 ah_origin_z(ah_origin_length)
  REAL*8 ah_radius(X0ah_radius,X1ah_radius,X2ah_radius,ah_radius_length)
  INTEGER*4 ah_really_initial_find_flag(ah_flags_length)
  INTEGER*4 ah_search_flag(ah_flags_length)
  REAL*8 b_phi_VolInt
  REAL*8 bh_posn_x(X0bh_posns)
  REAL*8 bh_posn_y(X0bh_posns)
  REAL*8 bh_posn_z(X0bh_posns)
  REAL*8 bh_radius_z(X0bh_posns)
  REAL*8 brem_qee_VolInt
  REAL*8 brem_qei_VolInt
  REAL*8 chi(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 chi_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 chi_p_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 chi_rad(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 chi_rad_nue(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 chi_rad_nux(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 chi_rhs(X0BSSN_rhs,X1BSSN_rhs,X2BSSN_rhs)
  REAL*8 cmax(X0mhd_nosync,X1mhd_nosync,X2mhd_nosync)
  REAL*8 cmax_rad(X0mhd_nosync,X1mhd_nosync,X2mhd_nosync)
  REAL*8 cmax_rad_nue(X0mhd_nosync,X1mhd_nosync,X2mhd_nosync)
  REAL*8 cmax_rad_nux(X0mhd_nosync,X1mhd_nosync,X2mhd_nosync)
  REAL*8 cmin(X0mhd_nosync,X1mhd_nosync,X2mhd_nosync)
  REAL*8 cmin_rad(X0mhd_nosync,X1mhd_nosync,X2mhd_nosync)
  REAL*8 cmin_rad_nue(X0mhd_nosync,X1mhd_nosync,X2mhd_nosync)
  REAL*8 cmin_rad_nux(X0mhd_nosync,X1mhd_nosync,X2mhd_nosync)
  REAL*8 dBx_m(X0mhd_sync_nablas,X1mhd_sync_nablas,X2mhd_sync_nablas)
  REAL*8 dBy_m(X0mhd_sync_nablas,X1mhd_sync_nablas,X2mhd_sync_nablas)
  REAL*8 dBz_m(X0mhd_sync_nablas,X1mhd_sync_nablas,X2mhd_sync_nablas)
  REAL*8 dP_m(X0mhd_sync_nablas,X1mhd_sync_nablas,X2mhd_sync_nablas)
  REAL*8 dT_fluid_m(X0micphys_sync_nablas,X1micphys_sync_nablas,X2micphys_sync_nablas)
  REAL*8 ddrbrvec(X0surf_params)
  REAL*8 density_modes_i1
  REAL*8 density_modes_i2
  REAL*8 density_modes_i3
  REAL*8 density_modes_i4
  REAL*8 density_modes_i5
  REAL*8 density_modes_i6
  REAL*8 density_modes_r0
  REAL*8 density_modes_r1
  REAL*8 density_modes_r2
  REAL*8 density_modes_r3
  REAL*8 density_modes_r4
  REAL*8 density_modes_r5
  REAL*8 density_modes_r6
  REAL*8 div_beta(X0BSSN_aux_private,X1BSSN_aux_private,X2BSSN_aux_private)
  REAL*8 drbrvec(X0surf_params)
  REAL*8 drhoYe_m(X0micphys_sync_nablas,X1micphys_sync_nablas,X2micphys_sync_nablas)
  REAL*8 drhoYe_m_x(X0mhd_sync_nablas_diag,X1mhd_sync_nablas_diag,X2mhd_sync_nablas_diag)
  REAL*8 drhoYe_m_xp1(X0mhd_sync_nablas_diag,X1mhd_sync_nablas_diag,X2mhd_sync_nablas_diag)
  REAL*8 drho_b_m(X0mhd_sync_nablas,X1mhd_sync_nablas,X2mhd_sync_nablas)
  REAL*8 drho_b_m_x(X0mhd_sync_nablas_diag,X1mhd_sync_nablas_diag,X2mhd_sync_nablas_diag)
  REAL*8 drho_b_m_xp1(X0mhd_sync_nablas_diag,X1mhd_sync_nablas_diag,X2mhd_sync_nablas_diag)
  REAL*8 dvx_m(X0mhd_sync_nablas,X1mhd_sync_nablas,X2mhd_sync_nablas)
  REAL*8 dvx_m_x(X0mhd_sync_nablas_diag,X1mhd_sync_nablas_diag,X2mhd_sync_nablas_diag)
  REAL*8 dvx_m_xp1(X0mhd_sync_nablas_diag,X1mhd_sync_nablas_diag,X2mhd_sync_nablas_diag)
  REAL*8 dvy_m(X0mhd_sync_nablas,X1mhd_sync_nablas,X2mhd_sync_nablas)
  REAL*8 dvz_m(X0mhd_sync_nablas,X1mhd_sync_nablas,X2mhd_sync_nablas)
  REAL*8 em_J_between_VolInt
  REAL*8 em_energy2_VolInt
  REAL*8 em_energy_VolInt
  REAL*8 em_energy_between_VolInt
  REAL*8 em_energy_outsideBH_VolInt
  REAL*8 em_energy_outsideradius1_VolInt
  REAL*8 em_energy_outsideradius2_VolInt
  REAL*8 emission_gf(X0microphys_primitives,X1microphys_primitives,X2microphys_primitives)
  REAL*8 emission_gf_nue(X0microphys_primitives,X1microphys_primitives,X2microphys_primitives)
  REAL*8 emission_gf_nux(X0microphys_primitives,X1microphys_primitives,X2microphys_primitives)
  REAL*8 eps_cld(X0microphys_primitives,X1microphys_primitives,X2microphys_primitives)
  REAL*8 eps_tab(X0eos_params1)
  REAL*8 eps_thermal(X0microphys_primitives,X1microphys_primitives,X2microphys_primitives)
  REAL*8 eps_tot(X0microphys_primitives,X1microphys_primitives,X2microphys_primitives)
  REAL*8 eta_nue(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  INTEGER*4 excision_zone_gf(X0excision_int_gfs,X1excision_int_gfs,X2excision_int_gfs)
  REAL*8 fluid_J_horiz_VolInt
  REAL*8 fluid_J_r1_VolInt
  REAL*8 fluid_J_r2_VolInt
  REAL*8 fluid_J_r3_VolInt
  REAL*8 fluid_energy_horiz_VolInt
  REAL*8 fluid_energy_r1_VolInt
  REAL*8 fluid_energy_r2_VolInt
  REAL*8 fluid_energy_r3_VolInt
  REAL*8 gamma_tab(X0eos_params2)
  REAL*8 gconx(X0BSSN_diag_restrict,X1BSSN_diag_restrict,X2BSSN_diag_restrict)
  REAL*8 gcony(X0BSSN_diag_restrict,X1BSSN_diag_restrict,X2BSSN_diag_restrict)
  REAL*8 gconz(X0BSSN_diag_restrict,X1BSSN_diag_restrict,X2BSSN_diag_restrict)
  REAL*8 gupxx(X0BSSN_gupij,X1BSSN_gupij,X2BSSN_gupij)
  REAL*8 gupxx_f(X0mhd_sync_metric_facevals,X1mhd_sync_metric_facevals,X2mhd_sync_metric_facevals)
  REAL*8 gupxy(X0BSSN_gupij,X1BSSN_gupij,X2BSSN_gupij)
  REAL*8 gupxy_f(X0mhd_sync_metric_facevals,X1mhd_sync_metric_facevals,X2mhd_sync_metric_facevals)
  REAL*8 gupxz(X0BSSN_gupij,X1BSSN_gupij,X2BSSN_gupij)
  REAL*8 gupxz_f(X0mhd_sync_metric_facevals,X1mhd_sync_metric_facevals,X2mhd_sync_metric_facevals)
  REAL*8 gupyy(X0BSSN_gupij,X1BSSN_gupij,X2BSSN_gupij)
  REAL*8 gupyy_f(X0mhd_sync_metric_facevals,X1mhd_sync_metric_facevals,X2mhd_sync_metric_facevals)
  REAL*8 gupyz(X0BSSN_gupij,X1BSSN_gupij,X2BSSN_gupij)
  REAL*8 gupyz_f(X0mhd_sync_metric_facevals,X1mhd_sync_metric_facevals,X2mhd_sync_metric_facevals)
  REAL*8 gupzz(X0BSSN_gupij,X1BSSN_gupij,X2BSSN_gupij)
  REAL*8 gupzz_f(X0mhd_sync_metric_facevals,X1mhd_sync_metric_facevals,X2mhd_sync_metric_facevals)
  REAL*8 gxx(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 gxx_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 gxx_p_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 gxx_f(X0mhd_sync_metric_facevals,X1mhd_sync_metric_facevals,X2mhd_sync_metric_facevals)
  REAL*8 gxx_rhs(X0BSSN_rhs,X1BSSN_rhs,X2BSSN_rhs)
  REAL*8 gxxx(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 gxxy(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 gxxz(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 gxy(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 gxy_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 gxy_p_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 gxy_f(X0mhd_sync_metric_facevals,X1mhd_sync_metric_facevals,X2mhd_sync_metric_facevals)
  REAL*8 gxy_rhs(X0BSSN_rhs,X1BSSN_rhs,X2BSSN_rhs)
  REAL*8 gxyx(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 gxyy(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 gxyz(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 gxz(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 gxz_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 gxz_p_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 gxz_f(X0mhd_sync_metric_facevals,X1mhd_sync_metric_facevals,X2mhd_sync_metric_facevals)
  REAL*8 gxz_rhs(X0BSSN_rhs,X1BSSN_rhs,X2BSSN_rhs)
  REAL*8 gxzx(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 gxzy(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 gxzz(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 gyy(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 gyy_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 gyy_p_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 gyy_f(X0mhd_sync_metric_facevals,X1mhd_sync_metric_facevals,X2mhd_sync_metric_facevals)
  REAL*8 gyy_rhs(X0BSSN_rhs,X1BSSN_rhs,X2BSSN_rhs)
  REAL*8 gyyx(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 gyyy(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 gyyz(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 gyz(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 gyz_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 gyz_p_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 gyz_f(X0mhd_sync_metric_facevals,X1mhd_sync_metric_facevals,X2mhd_sync_metric_facevals)
  REAL*8 gyz_rhs(X0BSSN_rhs,X1BSSN_rhs,X2BSSN_rhs)
  REAL*8 gyzx(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 gyzy(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 gyzz(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 gzz(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 gzz_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 gzz_p_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 gzz_f(X0mhd_sync_metric_facevals,X1mhd_sync_metric_facevals,X2mhd_sync_metric_facevals)
  REAL*8 gzz_rhs(X0BSSN_rhs,X1BSSN_rhs,X2BSSN_rhs)
  REAL*8 gzzx(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 gzzy(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 gzzz(X0BSSN_aux_restrict2,X1BSSN_aux_restrict2,X2BSSN_aux_restrict2)
  REAL*8 h(X0mhd_primitives,X1mhd_primitives,X2mhd_primitives)
  REAL*8 h_p(X0mhd_primitives,X1mhd_primitives,X2mhd_primitives)
  REAL*8 half_b2_u0_VolInt
  REAL*8 half_b2_u0_outsideBH_VolInt
  INTEGER*4 have_global_bdry_max(X0have_global_bdry,X1have_global_bdry)
  INTEGER*4 have_global_bdry_min(X0have_global_bdry,X1have_global_bdry)
  REAL*8 int_momentsi(X0gw_moment_arrays)
  REAL*8 int_momentsr(X0gw_moment_arrays)
  REAL*8 k_tab(X0eos_params2)
  REAL*8 ka_gf(X0microphys_primitives,X1microphys_primitives,X2microphys_primitives)
  REAL*8 ka_gf_nue(X0microphys_primitives,X1microphys_primitives,X2microphys_primitives)
  REAL*8 ka_gf_nux(X0microphys_primitives,X1microphys_primitives,X2microphys_primitives)
  REAL*8 ks_gf(X0microphys_primitives,X1microphys_primitives,X2microphys_primitives)
  REAL*8 ks_gf_nue(X0microphys_primitives,X1microphys_primitives,X2microphys_primitives)
  REAL*8 ks_gf_nux(X0microphys_primitives,X1microphys_primitives,X2microphys_primitives)
  REAL*8 kxx(X0BSSN_AH,X1BSSN_AH,X2BSSN_AH)
  REAL*8 kxy(X0BSSN_AH,X1BSSN_AH,X2BSSN_AH)
  REAL*8 kxz(X0BSSN_AH,X1BSSN_AH,X2BSSN_AH)
  REAL*8 kyy(X0BSSN_AH,X1BSSN_AH,X2BSSN_AH)
  REAL*8 kyz(X0BSSN_AH,X1BSSN_AH,X2BSSN_AH)
  REAL*8 kzz(X0BSSN_AH,X1BSSN_AH,X2BSSN_AH)
  REAL*8 lapm1(X0lapse_vars,X1lapse_vars,X2lapse_vars)
  REAL*8 lapm1_p(X0lapse_vars,X1lapse_vars,X2lapse_vars)
  REAL*8 lapm1_p_p(X0lapse_vars,X1lapse_vars,X2lapse_vars)
  REAL*8 lapm1_f(X0mhd_sync_metric_facevals,X1mhd_sync_metric_facevals,X2mhd_sync_metric_facevals)
  REAL*8 lapm1_rhs(X0lapse_rhs,X1lapse_rhs,X2lapse_rhs)
  REAL*8 lapm1_rot1(X0metric_spher_pol_1)
  REAL*8 lapm1_rot2(X0metric_spher_pol_2)
  REAL*8 lapm1_rot3(X0metric_spher_pol_3)
  REAL*8 lapset(X0lapse_vars,X1lapse_vars,X2lapse_vars)
  REAL*8 lapset_p(X0lapse_vars,X1lapse_vars,X2lapse_vars)
  REAL*8 lapset_p_p(X0lapse_vars,X1lapse_vars,X2lapse_vars)
  REAL*8 lapset_rhs(X0lapse_rhs,X1lapse_rhs,X2lapse_rhs)
  REAL*8 lapset_timederiv(X0lapse_vars_aux,X1lapse_vars_aux,X2lapse_vars_aux)
  REAL*8 lapsex(X0lapse_derivatives,X1lapse_derivatives,X2lapse_derivatives)
  REAL*8 lapsey(X0lapse_derivatives,X1lapse_derivatives,X2lapse_derivatives)
  REAL*8 lapsez(X0lapse_derivatives,X1lapse_derivatives,X2lapse_derivatives)
  REAL*8 mass_star1
  REAL*8 mass_star2
  REAL*8 mhd_chi_line(X0field_line_variables,X1field_line_variables,X2field_line_variables)
  REAL*8 mhd_chi_line_p(X0field_line_variables,X1field_line_variables,X2field_line_variables)
  REAL*8 mhd_chi_line_p_p(X0field_line_variables,X1field_line_variables,X2field_line_variables)
  REAL*8 mhd_chi_line_rhs(X0field_line_variables_rhs,X1field_line_variables_rhs,X2field_line_variables_rhs)
  REAL*8 mhd_psi_line(X0field_line_variables,X1field_line_variables,X2field_line_variables)
  REAL*8 mhd_psi_line_p(X0field_line_variables,X1field_line_variables,X2field_line_variables)
  REAL*8 mhd_psi_line_p_p(X0field_line_variables,X1field_line_variables,X2field_line_variables)
  REAL*8 mhd_psi_line_rhs(X0field_line_variables_rhs,X1field_line_variables_rhs,X2field_line_variables_rhs)
  REAL*8 mhd_st_x(X0mhd_conservatives,X1mhd_conservatives,X2mhd_conservatives)
  REAL*8 mhd_st_x_p(X0mhd_conservatives,X1mhd_conservatives,X2mhd_conservatives)
  REAL*8 mhd_st_x_p_p(X0mhd_conservatives,X1mhd_conservatives,X2mhd_conservatives)
  REAL*8 mhd_st_x_rhs(X0mhd_rhs,X1mhd_rhs,X2mhd_rhs)
  REAL*8 mhd_st_y(X0mhd_conservatives,X1mhd_conservatives,X2mhd_conservatives)
  REAL*8 mhd_st_y_p(X0mhd_conservatives,X1mhd_conservatives,X2mhd_conservatives)
  REAL*8 mhd_st_y_p_p(X0mhd_conservatives,X1mhd_conservatives,X2mhd_conservatives)
  REAL*8 mhd_st_y_rhs(X0mhd_rhs,X1mhd_rhs,X2mhd_rhs)
  REAL*8 mhd_st_z(X0mhd_conservatives,X1mhd_conservatives,X2mhd_conservatives)
  REAL*8 mhd_st_z_p(X0mhd_conservatives,X1mhd_conservatives,X2mhd_conservatives)
  REAL*8 mhd_st_z_p_p(X0mhd_conservatives,X1mhd_conservatives,X2mhd_conservatives)
  REAL*8 mhd_st_z_rhs(X0mhd_rhs,X1mhd_rhs,X2mhd_rhs)
  REAL*8 mhd_u_chi(X0field_line_variables,X1field_line_variables,X2field_line_variables)
  REAL*8 mhd_u_chi_p(X0field_line_variables,X1field_line_variables,X2field_line_variables)
  REAL*8 mhd_u_chi_p_p(X0field_line_variables,X1field_line_variables,X2field_line_variables)
  REAL*8 mhd_u_chi_rhs(X0field_line_variables_rhs,X1field_line_variables_rhs,X2field_line_variables_rhs)
  REAL*8 mhd_u_psi(X0field_line_variables,X1field_line_variables,X2field_line_variables)
  REAL*8 mhd_u_psi_p(X0field_line_variables,X1field_line_variables,X2field_line_variables)
  REAL*8 mhd_u_psi_p_p(X0field_line_variables,X1field_line_variables,X2field_line_variables)
  REAL*8 mhd_u_psi_rhs(X0field_line_variables_rhs,X1field_line_variables_rhs,X2field_line_variables_rhs)
  REAL*8 minternal_cold_horiz_VolInt
  REAL*8 minternal_cold_r1_VolInt
  REAL*8 minternal_cold_r2_VolInt
  REAL*8 minternal_cold_r3_VolInt
  REAL*8 minternal_horiz_VolInt
  REAL*8 minternal_r1_VolInt
  REAL*8 minternal_r2_VolInt
  REAL*8 minternal_r3_VolInt
  REAL*8 mom_const_VolIntD
  REAL*8 mom_const_excised_VolIntD
  REAL*8 momentsi(X0gw_moment_arrays)
  REAL*8 momentsi_old(X0gw_moment_arrays)
  REAL*8 momentsr(X0gw_moment_arrays)
  REAL*8 momentsr_old(X0gw_moment_arrays)
  REAL*8 momx_const_VolIntN
  REAL*8 momx_const_excised_VolIntN
  REAL*8 momy_const_VolIntN
  REAL*8 momy_const_excised_VolIntN
  REAL*8 momz_const_VolIntN
  REAL*8 momz_const_excised_VolIntN
  REAL*8 monopole_VolInt
  REAL*8 monopole_outsideBH_VolInt
  REAL*8 mu_nu(X0microphys_primitives,X1microphys_primitives,X2microphys_primitives)
  REAL*8 n_poly
  REAL*8 odd_momentsi(X0gw_moment_arrays)
  REAL*8 odd_momentsi_old(X0gw_moment_arrays)
  REAL*8 odd_momentsr(X0gw_moment_arrays)
  REAL*8 odd_momentsr_old(X0gw_moment_arrays)
  REAL*8 optd(X0microphys_primitives,X1microphys_primitives,X2microphys_primitives)
  REAL*8 optd_x(X0microphys_primitives,X1microphys_primitives,X2microphys_primitives)
  REAL*8 optd_y(X0microphys_primitives,X1microphys_primitives,X2microphys_primitives)
  REAL*8 optd_z(X0microphys_primitives,X1microphys_primitives,X2microphys_primitives)
  REAL*8 pfloor_gf(X0disk_atmosphere,X1disk_atmosphere,X2disk_atmosphere)
  REAL*8 phi(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 phi_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 phi_p_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 phi_f(X0mhd_sync_metric_facevals,X1mhd_sync_metric_facevals,X2mhd_sync_metric_facevals)
  REAL*8 phi_rhs(X0BSSN_rhs,X1BSSN_rhs,X2BSSN_rhs)
  REAL*8 phi_rot1(X0metric_spher_pol_1)
  REAL*8 phi_rot2(X0metric_spher_pol_2)
  REAL*8 phi_rot3(X0metric_spher_pol_3)
  REAL*8 phix(X0phi_derivs,X1phi_derivs,X2phi_derivs)
  REAL*8 phiy(X0phi_derivs,X1phi_derivs,X2phi_derivs)
  REAL*8 phiz(X0phi_derivs,X1phi_derivs,X2phi_derivs)
  REAL*8 psi(X0BSSN_AH,X1BSSN_AH,X2BSSN_AH)
  REAL*8 psi0im(X0Psi4_gw,X1Psi4_gw,X2Psi4_gw)
  REAL*8 psi0im_p(X0Psi4_gw,X1Psi4_gw,X2Psi4_gw)
  REAL*8 psi0im_p_p(X0Psi4_gw,X1Psi4_gw,X2Psi4_gw)
  REAL*8 psi0re(X0Psi4_gw,X1Psi4_gw,X2Psi4_gw)
  REAL*8 psi0re_p(X0Psi4_gw,X1Psi4_gw,X2Psi4_gw)
  REAL*8 psi0re_p_p(X0Psi4_gw,X1Psi4_gw,X2Psi4_gw)
  REAL*8 psi6phi(X0em_Phi,X1em_Phi,X2em_Phi)
  REAL*8 psi6phi_p(X0em_Phi,X1em_Phi,X2em_Phi)
  REAL*8 psi6phi_p_p(X0em_Phi,X1em_Phi,X2em_Phi)
  REAL*8 psi6phi_rhs(X0em_Phi_rhs,X1em_Phi_rhs,X2em_Phi_rhs)
  REAL*8 rad_const
  REAL*8 rad_energy_VolInt
  REAL*8 rad_energy_nue_VolInt
  REAL*8 rad_energy_nux_VolInt
  REAL*8 rbrvec(X0surf_params)
  REAL*8 refbd(X0BSSN_refbd,X1BSSN_refbd,X2BSSN_refbd)
  REAL*8 refbd_p(X0BSSN_refbd,X1BSSN_refbd,X2BSSN_refbd)
  REAL*8 refbd_p_p(X0BSSN_refbd,X1BSSN_refbd,X2BSSN_refbd)
  REAL*8 rho(X0BSSN_matter,X1BSSN_matter,X2BSSN_matter)
  REAL*8 rhoYe(X0micphys_conservatives,X1micphys_conservatives,X2micphys_conservatives)
  REAL*8 rhoYe_p(X0micphys_conservatives,X1micphys_conservatives,X2micphys_conservatives)
  REAL*8 rhoYe_p_p(X0micphys_conservatives,X1micphys_conservatives,X2micphys_conservatives)
  REAL*8 rhoYe_rhs(X0micphys_conservatives_rhs,X1micphys_conservatives_rhs,X2micphys_conservatives_rhs)
  REAL*8 rho_b(X0mhd_primitives,X1mhd_primitives,X2mhd_primitives)
  REAL*8 rho_b_atm_gf(X0disk_atmosphere,X1disk_atmosphere,X2disk_atmosphere)
  REAL*8 rho_bl(X0mhd_sync_rho_br_rho_bl,X1mhd_sync_rho_br_rho_bl,X2mhd_sync_rho_br_rho_bl)
  REAL*8 rho_br(X0mhd_sync_rho_br_rho_bl,X1mhd_sync_rho_br_rho_bl,X2mhd_sync_rho_br_rho_bl)
  REAL*8 rho_star(X0mhd_conservatives,X1mhd_conservatives,X2mhd_conservatives)
  REAL*8 rho_star_p(X0mhd_conservatives,X1mhd_conservatives,X2mhd_conservatives)
  REAL*8 rho_star_p_p(X0mhd_conservatives,X1mhd_conservatives,X2mhd_conservatives)
  REAL*8 rho_star_rhs(X0mhd_rhs,X1mhd_rhs,X2mhd_rhs)
  REAL*8 rho_tab(X0eos_params1)
  REAL*8 rhob_floor(X0mhd_nosync,X1mhd_nosync,X2mhd_nosync)
  REAL*8 rhovec(X0rhovecs)
  REAL*8 sbt(X0mhd_primitives,X1mhd_primitives,X2mhd_primitives)
  REAL*8 sbx(X0mhd_primitives,X1mhd_primitives,X2mhd_primitives)
  REAL*8 sby(X0mhd_primitives,X1mhd_primitives,X2mhd_primitives)
  REAL*8 sbz(X0mhd_primitives,X1mhd_primitives,X2mhd_primitives)
  REAL*8 shiftp_rot1(X0metric_spher_pol_1)
  REAL*8 shiftp_rot2(X0metric_spher_pol_2)
  REAL*8 shiftp_rot3(X0metric_spher_pol_3)
  REAL*8 shiftr_rot1(X0metric_spher_pol_1)
  REAL*8 shiftr_rot2(X0metric_spher_pol_2)
  REAL*8 shiftr_rot3(X0metric_spher_pol_3)
  REAL*8 shiftth_rot1(X0metric_spher_pol_1)
  REAL*8 shiftth_rot2(X0metric_spher_pol_2)
  REAL*8 shiftth_rot3(X0metric_spher_pol_3)
  REAL*8 shiftx(X0shift_vars,X1shift_vars,X2shift_vars)
  REAL*8 shiftx_p(X0shift_vars,X1shift_vars,X2shift_vars)
  REAL*8 shiftx_p_p(X0shift_vars,X1shift_vars,X2shift_vars)
  REAL*8 shiftx_f(X0mhd_sync_metric_facevals,X1mhd_sync_metric_facevals,X2mhd_sync_metric_facevals)
  REAL*8 shiftx_rhs(X0shift_rhs,X1shift_rhs,X2shift_rhs)
  REAL*8 shiftxt(X0shift_vars,X1shift_vars,X2shift_vars)
  REAL*8 shiftxt_p(X0shift_vars,X1shift_vars,X2shift_vars)
  REAL*8 shiftxt_p_p(X0shift_vars,X1shift_vars,X2shift_vars)
  REAL*8 shiftxt_rhs(X0shift_rhs,X1shift_rhs,X2shift_rhs)
  REAL*8 shiftxt_timederiv(X0shift_vars_temp,X1shift_vars_temp,X2shift_vars_temp)
  REAL*8 shifty(X0shift_vars,X1shift_vars,X2shift_vars)
  REAL*8 shifty_p(X0shift_vars,X1shift_vars,X2shift_vars)
  REAL*8 shifty_p_p(X0shift_vars,X1shift_vars,X2shift_vars)
  REAL*8 shifty_f(X0mhd_sync_metric_facevals,X1mhd_sync_metric_facevals,X2mhd_sync_metric_facevals)
  REAL*8 shifty_rhs(X0shift_rhs,X1shift_rhs,X2shift_rhs)
  REAL*8 shiftyt(X0shift_vars,X1shift_vars,X2shift_vars)
  REAL*8 shiftyt_p(X0shift_vars,X1shift_vars,X2shift_vars)
  REAL*8 shiftyt_p_p(X0shift_vars,X1shift_vars,X2shift_vars)
  REAL*8 shiftyt_rhs(X0shift_rhs,X1shift_rhs,X2shift_rhs)
  REAL*8 shiftyt_timederiv(X0shift_vars_temp,X1shift_vars_temp,X2shift_vars_temp)
  REAL*8 shiftz(X0shift_vars,X1shift_vars,X2shift_vars)
  REAL*8 shiftz_p(X0shift_vars,X1shift_vars,X2shift_vars)
  REAL*8 shiftz_p_p(X0shift_vars,X1shift_vars,X2shift_vars)
  REAL*8 shiftz_f(X0mhd_sync_metric_facevals,X1mhd_sync_metric_facevals,X2mhd_sync_metric_facevals)
  REAL*8 shiftz_rhs(X0shift_rhs,X1shift_rhs,X2shift_rhs)
  REAL*8 shiftzt(X0shift_vars,X1shift_vars,X2shift_vars)
  REAL*8 shiftzt_p(X0shift_vars,X1shift_vars,X2shift_vars)
  REAL*8 shiftzt_p_p(X0shift_vars,X1shift_vars,X2shift_vars)
  REAL*8 shiftzt_rhs(X0shift_rhs,X1shift_rhs,X2shift_rhs)
  REAL*8 shiftzt_timederiv(X0shift_vars_temp,X1shift_vars_temp,X2shift_vars_temp)
  REAL*8 smallb2(X0mhd_primitives,X1mhd_primitives,X2mhd_primitives)
  REAL*8 st_x(X0mhd_primitives,X1mhd_primitives,X2mhd_primitives)
  REAL*8 st_y(X0mhd_primitives,X1mhd_primitives,X2mhd_primitives)
  REAL*8 st_z(X0mhd_primitives,X1mhd_primitives,X2mhd_primitives)
  REAL*8 surfvec(X0surf_params)
  REAL*8 synch_VolInt
  REAL*8 tau(X0mhd_conservatives,X1mhd_conservatives,X2mhd_conservatives)
  REAL*8 tau_p(X0mhd_conservatives,X1mhd_conservatives,X2mhd_conservatives)
  REAL*8 tau_p_p(X0mhd_conservatives,X1mhd_conservatives,X2mhd_conservatives)
  REAL*8 tau_rad(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 tau_rad_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 tau_rad_p_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 tau_rad_nue(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 tau_rad_nue_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 tau_rad_nue_p_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 tau_rad_nue_rhs(X0rad_conservatives_rhs,X1rad_conservatives_rhs,X2rad_conservatives_rhs)
  REAL*8 tau_rad_nux(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 tau_rad_nux_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 tau_rad_nux_p_p(X0rad_conservatives,X1rad_conservatives,X2rad_conservatives)
  REAL*8 tau_rad_nux_rhs(X0rad_conservatives_rhs,X1rad_conservatives_rhs,X2rad_conservatives_rhs)
  REAL*8 tau_rad_rhs(X0rad_conservatives_rhs,X1rad_conservatives_rhs,X2rad_conservatives_rhs)
  REAL*8 tau_rhs(X0mhd_rhs,X1mhd_rhs,X2mhd_rhs)
  REAL*8 temp1(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp10(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp11(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp12(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp13(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp14(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp15(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp16(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp17(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp18(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp19(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp2(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp20(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp21(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp22(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp23(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp24(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp25(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp26(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp3(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp4(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp5(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp6(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp7(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp8(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp9(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 temp_g00(X0mhd_temps,X1mhd_temps,X2mhd_temps)
  REAL*8 trK(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 trK_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 trK_p_p(X0BSSN_vars,X1BSSN_vars,X2BSSN_vars)
  REAL*8 trK_rhs(X0BSSN_rhs,X1BSSN_rhs,X2BSSN_rhs)
  REAL*8 trR(X0BSSN_aux_private,X1BSSN_aux_private,X2BSSN_aux_private)
  REAL*8 trRtilde(X0BSSN_diag_restrict,X1BSSN_diag_restrict,X2BSSN_diag_restrict)
  REAL*8 u0(X0mhd_nosync,X1mhd_nosync,X2mhd_nosync)
  REAL*8 v02_rad_nuel(X0mhd_nosync,X1mhd_nosync,X2mhd_nosync)
  REAL*8 v02_rad_nuer(X0mhd_nosync,X1mhd_nosync,X2mhd_nosync)
  REAL*8 v02_rad_nuxl(X0mhd_nosync,X1mhd_nosync,X2mhd_nosync)
  REAL*8 v02_rad_nuxr(X0mhd_nosync,X1mhd_nosync,X2mhd_nosync)
  REAL*8 v02_radl(X0mhd_nosync,X1mhd_nosync,X2mhd_nosync)
  REAL*8 v02_radr(X0mhd_nosync,X1mhd_nosync,X2mhd_nosync)
  REAL*8 v02l(X0mhd_nosync,X1mhd_nosync,X2mhd_nosync)
  REAL*8 v02r(X0mhd_nosync,X1mhd_nosync,X2mhd_nosync)
  REAL*8 vvec(X0rhovecs)
  REAL*8 vx(X0mhd_vs,X1mhd_vs,X2mhd_vs)
  REAL*8 vxl(X0mhd_sync_lr_hydro_quantities,X1mhd_sync_lr_hydro_quantities,X2mhd_sync_lr_hydro_quantities)
  REAL*8 vxr(X0mhd_sync_lr_hydro_quantities,X1mhd_sync_lr_hydro_quantities,X2mhd_sync_lr_hydro_quantities)
  REAL*8 vy(X0mhd_vs,X1mhd_vs,X2mhd_vs)
  REAL*8 vyl(X0mhd_sync_lr_hydro_quantities,X1mhd_sync_lr_hydro_quantities,X2mhd_sync_lr_hydro_quantities)
  REAL*8 vyr(X0mhd_sync_lr_hydro_quantities,X1mhd_sync_lr_hydro_quantities,X2mhd_sync_lr_hydro_quantities)
  REAL*8 vz(X0mhd_vs,X1mhd_vs,X2mhd_vs)
  REAL*8 vzl(X0mhd_sync_lr_hydro_quantities,X1mhd_sync_lr_hydro_quantities,X2mhd_sync_lr_hydro_quantities)
  REAL*8 vzr(X0mhd_sync_lr_hydro_quantities,X1mhd_sync_lr_hydro_quantities,X2mhd_sync_lr_hydro_quantities)
  REAL*8 w(X0mhd_primitives,X1mhd_primitives,X2mhd_primitives)
  REAL*8 zeta_rad(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 zeta_rad_nue(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  REAL*8 zeta_rad_nux(X0rad_pressure,X1rad_pressure,X2rad_pressure)
  INTEGER X0coordinates
  INTEGER X0mask
  INTEGER X0radii
  INTEGER X0sf_radius
  INTEGER X0space_mask_group
  INTEGER X1coordinates
  INTEGER X1mask
  INTEGER X1sf_radius
  INTEGER X1space_mask_group
  INTEGER X2coordinates
  INTEGER X2mask
  INTEGER X2space_mask_group
  INTEGER active_length
  INTEGER num_levels_length
  INTEGER positions_length
  INTEGER radii_length
  INTEGER sf_active_length
  INTEGER sf_coordinate_descriptors_length
  INTEGER sf_info_length
  INTEGER sf_maxreflevel_length
  INTEGER sf_minreflevel_length
  INTEGER sf_origin_length
  INTEGER sf_radius_length
  INTEGER sf_shape_descriptors_length
  INTEGER sf_valid_length
  INTEGER*4 active(active_length)
  REAL*8 coarse_dx
  REAL*8 coarse_dy
  REAL*8 coarse_dz
  REAL*8 emask(X0mask,X1mask,X2mask)
  INTEGER*4 num_levels(num_levels_length)
  REAL*8 position_x(positions_length)
  REAL*8 position_y(positions_length)
  REAL*8 position_z(positions_length)
  REAL*8 r(X0coordinates,X1coordinates,X2coordinates)
  REAL*8 radius(X0radii,radii_length)
  INTEGER*4 sf_active(sf_active_length)
  REAL*8 sf_area(sf_info_length)
  REAL*8 sf_centroid_x(sf_info_length)
  REAL*8 sf_centroid_y(sf_info_length)
  REAL*8 sf_centroid_z(sf_info_length)
  REAL*8 sf_delta_phi(sf_coordinate_descriptors_length)
  REAL*8 sf_delta_theta(sf_coordinate_descriptors_length)
  REAL*8 sf_max_radius(sf_info_length)
  REAL*8 sf_max_x(sf_info_length)
  REAL*8 sf_max_y(sf_info_length)
  REAL*8 sf_max_z(sf_info_length)
  INTEGER*4 sf_maxreflevel(sf_maxreflevel_length)
  REAL*8 sf_mean_radius(sf_info_length)
  REAL*8 sf_min_radius(sf_info_length)
  REAL*8 sf_min_x(sf_info_length)
  REAL*8 sf_min_y(sf_info_length)
  REAL*8 sf_min_z(sf_info_length)
  INTEGER*4 sf_minreflevel(sf_minreflevel_length)
  INTEGER*4 sf_nghostsphi(sf_shape_descriptors_length)
  INTEGER*4 sf_nghoststheta(sf_shape_descriptors_length)
  INTEGER*4 sf_nphi(sf_shape_descriptors_length)
  INTEGER*4 sf_ntheta(sf_shape_descriptors_length)
  REAL*8 sf_origin_phi(sf_coordinate_descriptors_length)
  REAL*8 sf_origin_theta(sf_coordinate_descriptors_length)
  REAL*8 sf_origin_x(sf_origin_length)
  REAL*8 sf_origin_y(sf_origin_length)
  REAL*8 sf_origin_z(sf_origin_length)
  REAL*8 sf_quadrupole_xx(sf_info_length)
  REAL*8 sf_quadrupole_xy(sf_info_length)
  REAL*8 sf_quadrupole_xz(sf_info_length)
  REAL*8 sf_quadrupole_yy(sf_info_length)
  REAL*8 sf_quadrupole_yz(sf_info_length)
  REAL*8 sf_quadrupole_zz(sf_info_length)
  REAL*8 sf_radius(X0sf_radius,X1sf_radius,sf_radius_length)
  INTEGER*4 sf_valid(sf_valid_length)
  INTEGER*4 space_mask(X0space_mask_group,X1space_mask_group,X2space_mask_group)
  REAL*8 x(X0coordinates,X1coordinates,X2coordinates)
  REAL*8 y(X0coordinates,X1coordinates,X2coordinates)
  REAL*8 z(X0coordinates,X1coordinates,X2coordinates)
  
    REAL*8  excision_radius
  REAL*8  run_time
  INTEGER*4 Symmetry
  INTEGER*4 bssn_enable
  INTEGER*4 cowling_enable
  INTEGER*4 excision_enable
  INTEGER*4 fisheye_enable
  INTEGER*4 iter_count
  INTEGER*4 number_of_mol_ministeps
  INTEGER*4 rot_metric
  INTEGER*4 trA_detg_enforce
  COMMON /cctk_params_global/excision_radius,run_time,Symmetry,bssn_enable,cowling_enable,excision_enable,fisheye_enable,iter_count&
  &,number_of_mol_ministeps,rot_metric,trA_detg_enforce
  INTEGER*4 iteration_to_compute_temp
  INTEGER*4 iteration_to_insert_rad
  INTEGER*4 num_CO
  COMMON /bhnsrest/iteration_to_compute_temp,iteration_to_insert_rad,num_CO
  REAL*8  Aphi_power
  REAL*8  BigM
  REAL*8  Con_Pol
  REAL*8  Con_Tor
  REAL*8  Erad_over_rho
  REAL*8  Iloop
  REAL*8  Iloop2
  REAL*8  MRI_wavelength_calculator_dxy
  REAL*8  Omega_value
  REAL*8  P_deplete_bhns
  REAL*8  a2oa1_perturb
  REAL*8  alpha_rho_cut_off
  REAL*8  ampl_perturb
  REAL*8  angle_to_tilt_magnetic_fields
  REAL*8  angle_to_tilt_magnetic_fieldsII
  REAL*8  beta_ext
  REAL*8  betam1
  REAL*8  betam1_2
  REAL*8  bhns_P_max
  REAL*8  bhns_P_max2
  REAL*8  bhns_R_NS
  REAL*8  bhns_R_NS2
  REAL*8  bhns_avg_betam1
  REAL*8  bhns_bh_filling_radius(100)
  REAL*8  bhns_fac_atm
  REAL*8  bhns_max_b2
  REAL*8  bhns_max_b2_posn_x
  REAL*8  bhns_max_b2_posn_y
  REAL*8  bhns_max_b2_posn_z
  REAL*8  bhns_max_rho_b
  REAL*8  bhns_rho_b_max_posn_x
  REAL*8  bhns_rho_b_max_posn_y
  REAL*8  bhns_rho_b_max_posn_z
  REAL*8  bhns_rhob_fac
  REAL*8  bhns_rhob_max
  REAL*8  bhns_tracer_r
  REAL*8  bhns_tracer_rin
  REAL*8  bhns_tracer_x0
  REAL*8  bhns_tracer_y0
  REAL*8  bhns_tracer_z0
  REAL*8  bhns_vx_CM
  REAL*8  bhns_vx_CM2
  REAL*8  bhns_vy_CM
  REAL*8  bhns_vy_CM2
  REAL*8  bhns_vz_CM
  REAL*8  bhns_vz_CM2
  REAL*8  conloop1
  REAL*8  conloop2
  REAL*8  horiz_radius
  REAL*8  initial_ns2_coord_x
  REAL*8  initial_ns2_coord_y
  REAL*8  initial_ns_coord_x
  REAL*8  initial_ns_coord_y
  REAL*8  initial_wd_coord_x
  REAL*8  initial_wd_coord_y
  REAL*8  lambda_perturb
  REAL*8  moncrief_radius_GW(11)
  REAL*8  mythbusters_boost_factor
  REAL*8  p_c
  REAL*8  particle_cone_angle
  REAL*8  particle_cylinder_cone_zmax
  REAL*8  particle_cylinder_cone_zmin
  REAL*8  r0
  REAL*8  r0_ah
  REAL*8  rad_rhob_fac
  REAL*8  radi_perturb
  REAL*8  rhob_fac2
  REAL*8  rhob_o_b2
  REAL*8  rhobatm_falloff_power
  REAL*8  rloop
  REAL*8  rloop2
  REAL*8  surfxmax
  REAL*8  surfxmin
  REAL*8  surfzmax
  REAL*8  tracer_x_max
  REAL*8  tracer_x_min
  REAL*8  tracer_y_max
  REAL*8  tracer_y_min
  REAL*8  tracer_z_max
  REAL*8  tracer_z_min
  REAL*8  xh0
  REAL*8  yh0
  REAL*8  zh0
  INTEGER*4 ATM_TYPE
  INTEGER*4 INPUTARRAY_PHISIZE
  INTEGER*4 INPUTARRAY_THETASIZE
  INTEGER*4 ITERATION_TO_BOOST_MAGNETIC_FIELDS
  INTEGER*4 ITERATION_TO_INSERT_MAGNETIC_FIELDS
  INTEGER*4 ITERATION_TO_output_MRI_wavelength
  INTEGER*4 MRI_wavelength_calculator_Nxy
  INTEGER*4 NUM_ZERO_PTS
  INTEGER*4 N_particles_to_trace
  INTEGER*4 RADIAL_INTERP_ORDER
  INTEGER*4 RESET_RHO_B_ATM
  INTEGER*4 alpha_diagnostic
  INTEGER*4 angle_to_tilt_90
  INTEGER*4 bhns_B_v
  INTEGER*4 bhns_domain
  INTEGER*4 bhns_particle_tracer_start
  INTEGER*4 bhns_regrid_input_enable
  INTEGER*4 bhns_regrid_output_enable_iter
  INTEGER*4 center_Bfields_around_BH_instead
  INTEGER*4 chunklet_dump_every
  INTEGER*4 chunklet_procs_at_a_time
  INTEGER*4 em_field_type
  INTEGER*4 genID_cmdline_output_enable
  INTEGER*4 initial_particle_geometry
  INTEGER*4 moncrief_gw_num_radii
  INTEGER*4 nperturb
  INTEGER*4 ntot_bhns
  INTEGER*4 particle_center
  INTEGER*4 particle_tracer_substep_every
  INTEGER*4 piecewise
  INTEGER*4 refill_horizons_magfields_every
  INTEGER*4 reset_shift_lapse
  INTEGER*4 subtract_off_Omega_r
  INTEGER*4 superposition
  INTEGER*4 surfphinum
  INTEGER*4 surfxnum
  INTEGER*4 surfznum
  INTEGER*4 two_ns
  INTEGER*4 unequalmass
  INTEGER*4 use_new_bhns_initial_data
  COMMON /bhnspriv/Aphi_power,BigM,Con_Pol,Con_Tor,Erad_over_rho,Iloop,Iloop2,MRI_wavelength_calculator_dxy,Omega_value,P_deplete_b&
  &hns,a2oa1_perturb,alpha_rho_cut_off,ampl_perturb,angle_to_tilt_magnetic_fields,angle_to_tilt_magnetic_fieldsII,beta_ext,betam1,b&
  &etam1_2,bhns_P_max,bhns_P_max2,bhns_R_NS,bhns_R_NS2,bhns_avg_betam1,bhns_bh_filling_radius,bhns_fac_atm,bhns_max_b2,bhns_max_b2_&
  &posn_x,bhns_max_b2_posn_y,bhns_max_b2_posn_z,bhns_max_rho_b,bhns_rho_b_max_posn_x,bhns_rho_b_max_posn_y,bhns_rho_b_max_posn_z,bh&
  &ns_rhob_fac,bhns_rhob_max,bhns_tracer_r,bhns_tracer_rin,bhns_tracer_x0,bhns_tracer_y0,bhns_tracer_z0,bhns_vx_CM,bhns_vx_CM2,bhns&
  &_vy_CM,bhns_vy_CM2,bhns_vz_CM,bhns_vz_CM2,conloop1,conloop2,horiz_radius,initial_ns2_coord_x,initial_ns2_coord_y,initial_ns_coor&
  &d_x,initial_ns_coord_y,initial_wd_coord_x,initial_wd_coord_y,lambda_perturb,moncrief_radius_GW,mythbusters_boost_factor,p_c,part&
  &icle_cone_angle,particle_cylinder_cone_zmax,particle_cylinder_cone_zmin,r0,r0_ah,rad_rhob_fac,radi_perturb,rhob_fac2,rhob_o_b2,r&
  &hobatm_falloff_power,rloop,rloop2,surfxmax,surfxmin,surfzmax,tracer_x_max,tracer_x_min,tracer_y_max,tracer_y_min,tracer_z_max,tr&
  &acer_z_min,xh0,yh0,zh0,ATM_TYPE,INPUTARRAY_PHISIZE,INPUTARRAY_THETASIZE,ITERATION_TO_BOOST_MAGNETIC_FIELDS,ITERATION_TO_INSERT_M&
  &AGNETIC_FIELDS,ITERATION_TO_output_MRI_wavelength,MRI_wavelength_calculator_Nxy,NUM_ZERO_PTS,N_particles_to_trace,RADIAL_INTERP_&
  &ORDER,RESET_RHO_B_ATM,alpha_diagnostic,angle_to_tilt_90,bhns_B_v,bhns_domain,bhns_particle_tracer_start,bhns_regrid_input_enable&
  &,bhns_regrid_output_enable_iter,center_Bfields_around_BH_instead,chunklet_dump_every,chunklet_procs_at_a_time,em_field_type,genI&
  &D_cmdline_output_enable,initial_particle_geometry,moncrief_gw_num_radii,nperturb,ntot_bhns,particle_center,particle_tracer_subst&
  &ep_every,piecewise,refill_horizons_magfields_every,reset_shift_lapse,subtract_off_Omega_r,superposition,surfphinum,surfxnum,surf&
  &znum,two_ns,unequalmass,use_new_bhns_initial_data
  REAL*8  CCTKH5
  REAL*8  CCTKH6
  REAL*8  CCTKH7
  REAL*8  CCTKH8
  REAL*8  binary_orb_freq
  REAL*8  CCTKH11
  REAL*8  CCTKH12
  REAL*8  CCTKH13
  REAL*8  CCTKH14
  REAL*8  CCTKH15
  REAL*8  CCTKH16
  REAL*8  CCTKH17
  REAL*8  CCTKH18
  REAL*8  rotation_phase
  REAL*8  CCTKH20
  REAL*8  CCTKH21
  REAL*8  CCTKH22
  INTEGER*4 CCTKH0
  INTEGER*4 CCTKH1
  INTEGER*4 CCTKH2
  INTEGER*4 CCTKH3
  INTEGER*4 CCTKH4
  INTEGER*4 CCTKH9
  INTEGER*4 CCTKH10
  INTEGER*4 CCTKH19
  INTEGER*4 CCTKH23
  COMMON /BSSNrest/CCTKH5,CCTKH6,CCTKH7,CCTKH8,binary_orb_freq,CCTKH11,CCTKH12,CCTKH13,CCTKH14,CCTKH15,CCTKH16,CCTKH17,CCTKH18,rota&
  &tion_phase,CCTKH20,CCTKH21,CCTKH22,CCTKH0,CCTKH1,CCTKH2,CCTKH3,CCTKH4,CCTKH9,CCTKH10,CCTKH19,CCTKH23
  REAL*8  CCTKH24
  REAL*8  M_ADM
  REAL*8  initial_monopole_value
  REAL*8  CCTKH30
  REAL*8  CCTKH31
  REAL*8  CCTKH33
  REAL*8  CCTKH34
  REAL*8  CCTKH35
  REAL*8  CCTKH36
  REAL*8  rhob_cutoff
  REAL*8  CCTKH37
  REAL*8  CCTKH38
  REAL*8  CCTKH39
  REAL*8  CCTKH40
  INTEGER*4 CCTKH25
  INTEGER*4 Nphi_points
  INTEGER*4 CCTKH26
  INTEGER*4 CCTKH27
  INTEGER*4 CCTKH28
  INTEGER*4 CCTKH29
  INTEGER*4 CCTKH32
  COMMON /DIAGNOSTICS_MHDrest/CCTKH24,M_ADM,initial_monopole_value,CCTKH30,CCTKH31,CCTKH33,CCTKH34,CCTKH35,CCTKH36,rhob_cutoff,CCTK&
  &H37,CCTKH38,CCTKH39,CCTKH40,CCTKH25,Nphi_points,CCTKH26,CCTKH27,CCTKH28,CCTKH29,CCTKH32
  REAL*8  BH_Vol_Excise_Radius
  REAL*8  CCTKH43
  REAL*8  ddrbddr
  REAL*8  CCTKH44
  REAL*8  drbdr
  REAL*8  CCTKH48
  REAL*8  rbr
  REAL*8  rsurf2
  REAL*8  surf_radius
  INTEGER*4 CCTKH41
  INTEGER*4 N_phi
  INTEGER*4 N_theta
  INTEGER*4 CCTKH42
  INTEGER*4 arithsurf
  INTEGER*4 CCTKH45
  INTEGER*4 CCTKH46
  INTEGER*4 CCTKH47
  INTEGER*4 nsurf
  INTEGER*4 CCTKH49
  INTEGER*4 CCTKH50
  INTEGER*4 CCTKH51
  INTEGER*4 CCTKH52
  INTEGER*4 scaledsurf
  INTEGER*4 sym_factor
  COMMON /DIAGNOSTICS_VACUUMrest/BH_Vol_Excise_Radius,CCTKH43,ddrbddr,CCTKH44,drbdr,CCTKH48,rbr,rsurf2,surf_radius,CCTKH41,N_phi,N_&
  &theta,CCTKH42,arithsurf,CCTKH45,CCTKH46,CCTKH47,nsurf,CCTKH49,CCTKH50,CCTKH51,CCTKH52,scaledsurf,sym_factor
  INTEGER*4 CCTKH53
  INTEGER*4 CCTKH54
  INTEGER*4 CCTKH55
  INTEGER*4 CCTKH56
  COMMON /DRIVERrest/CCTKH53,CCTKH54,CCTKH55,CCTKH56
  REAL*8  CCTKH57
  REAL*8  CCTKH58
  REAL*8  phi_GW
  REAL*8  CCTKH59
  REAL*8  CCTKH60
  REAL*8  theta_GW
  COMMON /GW_EXTRACTIONrest/CCTKH57,CCTKH58,phi_GW,CCTKH59,CCTKH60,theta_GW
  REAL*8  CCTKH66
  REAL*8  CCTKH80
  REAL*8  CCTKH88
  REAL*8  CCTKH90
  REAL*8  out_xyplane_z
  REAL*8  out_xzplane_y
  REAL*8  CCTKH94
  REAL*8  CCTKH96
  REAL*8  out_yzplane_x
  REAL*8  CCTKH99
  REAL*8  CCTKH101
  integer*8  CCTKH63
  integer*8  CCTKH64
  integer*8  CCTKH67
  integer*8  CCTKH70
  integer*8  CCTKH71
  integer*8  CCTKH72
  integer*8  CCTKH75
  integer*8  CCTKH76
  integer*8  CCTKH81
  integer*8  CCTKH82
  integer*8  CCTKH84
  integer*8  CCTKH103
  integer*8  CCTKH105
  integer*8  CCTKH107
  integer*8  CCTKH109
  integer*8  CCTKH110
  integer*8  CCTKH113
  INTEGER*4 CCTKH61
  INTEGER*4 CCTKH62
  INTEGER*4 CCTKH65
  INTEGER*4 CCTKH68
  INTEGER*4 CCTKH69
  INTEGER*4 CCTKH73
  INTEGER*4 CCTKH74
  INTEGER*4 CCTKH77
  INTEGER*4 CCTKH78
  INTEGER*4 CCTKH79
  INTEGER*4 out_every
  INTEGER*4 CCTKH83
  INTEGER*4 CCTKH85
  INTEGER*4 CCTKH86
  INTEGER*4 CCTKH87
  INTEGER*4 CCTKH89
  INTEGER*4 CCTKH91
  INTEGER*4 CCTKH92
  INTEGER*4 CCTKH93
  INTEGER*4 CCTKH95
  INTEGER*4 CCTKH97
  INTEGER*4 CCTKH98
  INTEGER*4 CCTKH100
  INTEGER*4 CCTKH102
  INTEGER*4 CCTKH104
  INTEGER*4 CCTKH106
  INTEGER*4 CCTKH108
  INTEGER*4 CCTKH111
  INTEGER*4 CCTKH112
  COMMON /IOrest/CCTKH66,CCTKH80,CCTKH88,CCTKH90,out_xyplane_z,out_xzplane_y,CCTKH94,CCTKH96,out_yzplane_x,CCTKH99,CCTKH101,CCTKH63&
  &,CCTKH64,CCTKH67,CCTKH70,CCTKH71,CCTKH72,CCTKH75,CCTKH76,CCTKH81,CCTKH82,CCTKH84,CCTKH103,CCTKH105,CCTKH107,CCTKH109,CCTKH110,CC&
  &TKH113,CCTKH61,CCTKH62,CCTKH65,CCTKH68,CCTKH69,CCTKH73,CCTKH74,CCTKH77,CCTKH78,CCTKH79,out_every,CCTKH83,CCTKH85,CCTKH86,CCTKH87&
  &,CCTKH89,CCTKH91,CCTKH92,CCTKH93,CCTKH95,CCTKH97,CCTKH98,CCTKH100,CCTKH102,CCTKH104,CCTKH106,CCTKH108,CCTKH111,CCTKH112
  REAL*8  CCTKH114
  REAL*8  CCTKH115
  REAL*8  CCTKH116
  REAL*8  CCTKH120
  REAL*8  CCTKH122
  REAL*8  CCTKH123
  REAL*8  CCTKH124
  REAL*8  CCTKH125
  integer*8  CCTKH121
  INTEGER*4 CCTKH117
  INTEGER*4 CCTKH118
  INTEGER*4 CCTKH119
  COMMON /MAGNETARrest/CCTKH114,CCTKH115,CCTKH116,CCTKH120,CCTKH122,CCTKH123,CCTKH124,CCTKH125,CCTKH121,CCTKH117,CCTKH118,CCTKH119
  REAL*8  Erad_atm_cut
  REAL*8  CCTKH127
  REAL*8  K_poly
  REAL*8  M_B
  REAL*8  CCTKH131
  REAL*8  P_fact
  REAL*8  CCTKH132
  REAL*8  Psi6threshold
  REAL*8  Sym_Bz
  REAL*8  T_fluid_cgs_atm
  REAL*8  CCTKH137
  REAL*8  CCTKH138
  REAL*8  CCTKH140
  REAL*8  eps_thermal_bhns
  REAL*8  ergo_sigma
  REAL*8  CCTKH145
  REAL*8  gamma_th
  REAL*8  CCTKH148
  REAL*8  CCTKH150
  REAL*8  CCTKH151
  REAL*8  CCTKH152
  REAL*8  CCTKH153
  REAL*8  CCTKH154
  REAL*8  pfloor
  REAL*8  CCTKH156
  REAL*8  rad_T_cutoff
  REAL*8  rad_T_fac
  REAL*8  rad_T_floor
  REAL*8  rad_T_pow
  REAL*8  CCTKH157
  REAL*8  CCTKH158
  REAL*8  CCTKH160
  REAL*8  CCTKH161
  REAL*8  rho_b_atm
  REAL*8  rho_b_max
  REAL*8  rho_fact
  REAL*8  rhos_max
  REAL*8  CCTKH164
  REAL*8  CCTKH165
  REAL*8  tau_atm
  REAL*8  tau_fact
  INTEGER*4 CCTKH126
  INTEGER*4 EM_BC
  INTEGER*4 CCTKH128
  INTEGER*4 CCTKH129
  INTEGER*4 CCTKH130
  INTEGER*4 CCTKH133
  INTEGER*4 CCTKH134
  INTEGER*4 CCTKH135
  INTEGER*4 CCTKH136
  INTEGER*4 compute_microphysics
  INTEGER*4 constrained_transport_scheme
  INTEGER*4 CCTKH139
  INTEGER*4 em_evolve_enable
  INTEGER*4 em_gauge
  INTEGER*4 CCTKH141
  INTEGER*4 enable_OS_collapse
  INTEGER*4 CCTKH142
  INTEGER*4 CCTKH143
  INTEGER*4 enable_shocktest_primitive_mode
  INTEGER*4 enable_trace_field_line
  INTEGER*4 eps_flag
  INTEGER*4 ergo_star
  INTEGER*4 CCTKH144
  INTEGER*4 horizon_enforce_rho_profile
  INTEGER*4 CCTKH146
  INTEGER*4 CCTKH147
  INTEGER*4 CCTKH149
  INTEGER*4 microphysics_scheme
  INTEGER*4 neos
  INTEGER*4 CCTKH155
  INTEGER*4 primitives_debug
  INTEGER*4 primitives_solver
  INTEGER*4 rad_closure_scheme
  INTEGER*4 rad_evolve_enable
  INTEGER*4 rad_fix
  INTEGER*4 CCTKH159
  INTEGER*4 CCTKH162
  INTEGER*4 CCTKH163
  INTEGER*4 tau_stildefix_enable
  INTEGER*4 CCTKH166
  INTEGER*4 CCTKH167
  INTEGER*4 CCTKH168
  COMMON /MHD_EVOLVErest/Erad_atm_cut,CCTKH127,K_poly,M_B,CCTKH131,P_fact,CCTKH132,Psi6threshold,Sym_Bz,T_fluid_cgs_atm,CCTKH137,CC&
  &TKH138,CCTKH140,eps_thermal_bhns,ergo_sigma,CCTKH145,gamma_th,CCTKH148,CCTKH150,CCTKH151,CCTKH152,CCTKH153,CCTKH154,pfloor,CCTKH&
  &156,rad_T_cutoff,rad_T_fac,rad_T_floor,rad_T_pow,CCTKH157,CCTKH158,CCTKH160,CCTKH161,rho_b_atm,rho_b_max,rho_fact,rhos_max,CCTKH&
  &164,CCTKH165,tau_atm,tau_fact,CCTKH126,EM_BC,CCTKH128,CCTKH129,CCTKH130,CCTKH133,CCTKH134,CCTKH135,CCTKH136,compute_microphysics&
  &,constrained_transport_scheme,CCTKH139,em_evolve_enable,em_gauge,CCTKH141,enable_OS_collapse,CCTKH142,CCTKH143,enable_shocktest_&
  &primitive_mode,enable_trace_field_line,eps_flag,ergo_star,CCTKH144,horizon_enforce_rho_profile,CCTKH146,CCTKH147,CCTKH149,microp&
  &hysics_scheme,neos,CCTKH155,primitives_debug,primitives_solver,rad_closure_scheme,rad_evolve_enable,rad_fix,CCTKH159,CCTKH162,CC&
  &TKH163,tau_stildefix_enable,CCTKH166,CCTKH167,CCTKH168
  REAL*8  CCTKH170
  REAL*8  CCTKH171
  INTEGER*4 CCTKH169
  INTEGER*4 track_bhns
  COMMON /MOVINGBOXrest/CCTKH170,CCTKH171,CCTKH169,track_bhns
  
    external     CCTK_PointerTo
  integer*8 CCTK_PointerTo
  interface
  integer function CCTK_Equals (arg1, arg2)
  implicit none
  integer*8 arg1
  character(*) arg2
  end function CCTK_Equals
  integer function CCTK_MyProc (cctkGH)
  implicit none
  integer*8 cctkGH
  end function CCTK_MyProc
  integer function CCTK_nProcs (cctkGH)
  implicit none
  integer*8 cctkGH
  end function CCTK_nProcs
  integer function CCTK_IsThornActive (name)
  implicit none
  character(*) name
  end function CCTK_IsThornActive
  integer*8 function CCTK_NullPointer ()
  implicit none
  end function CCTK_NullPointer
  end interface
  interface
  INTEGER*4 function HorizonRadiusInDirection (horizon_number, N_points, x, y, z, radius)
  implicit none
  INTEGER*4 horizon_number
  INTEGER*4 N_points
  REAL*8 x(*)
  REAL*8 y(*)
  REAL*8 z(*)
  REAL*8 radius(*)
  end function HorizonRadiusInDirection
  end interface
  
    interface
       subroutine gderivs_oct(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
         implicit none
         integer, dimension(3)                    :: ex
         real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
         real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
       end subroutine gderivs_oct
       subroutine gderivs_eq(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
         implicit none
         integer, dimension(3)                    :: ex
         real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
         real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
       end subroutine gderivs_eq
       subroutine gderivs_axi(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
         implicit none
         integer, dimension(3)                    :: ex
         real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
         real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
       end subroutine gderivs_axi
    end interface
    integer, dimension(3)                    :: ext
    real*8                                   :: dX,dY,dZ,Ab
    integer                                      :: handle,index,ierr
    integer                                  :: imin,imax,jmin,jmax,kmin,kmax
    integer                                  :: i,j,k,dummy
    integer                                  :: im1,jm1,km1,ip1,jp1,kp1
    real*8,allocatable,dimension(:,:,:)      :: A_phi,A_phix,A_phiy,A_phiz
    real*8                                   :: psim6,psim6_s,pomega2,al,sqrtg,sqrtg4,B2s
    real*8                                   :: fs4pi,B_xs,B_ys,B_zs,sb0,sb2,xn,yn
    real*8                                   :: sb_x,sb_y,sb_z,psi4,u_x,u_y,u_z
    real*8                                   :: psin,x_NS_CoM_coord,y_NS_CoM_coord
    real*8                                   :: Yijk,Yijkp1,Yijp1k,Yijp1kp1,Yip1jk,Yip1jkp1,Yip1jp1k,Yip1jp1kp1
    real*8                                   :: Xijk,Xijkp1,Xijp1k,Xijp1kp1,Xip1jk,Xip1jkp1,Xip1jp1k,Xip1jp1kp1
    real*8                                     :: Ax000,Ax001,Ax010,Ax011,Ax100,Ax101,Ax110,Ax111
    real*8                                   :: Ay000,Ay001,Ay010,Ay011,Ay100,Ay101,Ay110,Ay111
    integer                                  :: AXISYM,EQUATORIAL
    integer                                  :: OCTANT
    real*8                                   :: x_NS_CoM_coord2, y_NS_CoM_coord2
    integer, dimension(3) :: U_syms
    real*8, parameter                        :: SYM = 1.d0, ANTI = -1.d0
    REAL*8 reduction_value
    parameter(EQUATORIAL = 1, OCTANT = 2, AXISYM = 4)
    !
    print *, "Start BHNS_setup_emfield_part2_local_setBi"
    if(CCTK_ITERATION.eq.ITERATION_TO_INSERT_MAGNETIC_FIELDS .and. ITERATION_TO_INSERT_MAGNETIC_FIELDS.eq.0) then
       fs4pi = sqrt(4.d0*acos(-1.d0))
       ext = cctk_lsh
       dX = X(2,1,1) - X(1,1,1)
       dY = Y(1,2,1) - Y(1,1,1)
       dZ = Z(1,1,2) - Z(1,1,1)
       allocate(A_phi(ext(1),ext(2),ext(3)))
       allocate(A_phix(ext(1),ext(2),ext(3)))
       allocate(A_phiy(ext(1),ext(2),ext(3)))
       allocate(A_phiz(ext(1),ext(2),ext(3)))
       imin = 1
       jmin = 1
       kmin = 1
       imax = ext(1)
       jmax = ext(2)
       kmax = ext(3)
       ! Compute magnetic vector potential A_{\phi} with Ab=1
       ! Note that A_{\phi} is calculated from the center of each NS
       !     x_NS_CoM_coord = CoMx_VolInt/CoM_VolInt_denominator
       !     y_NS_CoM_coord = CoMy_VolInt/CoM_VolInt_denominator
       do k=1,cctk_lsh(3)
          do j=1,cctk_lsh(2)
             do i=1,cctk_lsh(1)
                if(isnan(Ax(i,j,k)).eqv. .true.) then
                   write(*,*) "1. Inside driver_setup_emfields_part2-local-setBi.F90, Ax is nan!!!"
                end if
             end do
          end do
       end do
       x_NS_CoM_coord = initial_ns_coord_x
       y_NS_CoM_coord = initial_ns_coord_y
       write(*,*) "INSIDE EMFIELDS SETUP: x,y of ns:",x_NS_CoM_coord,y_NS_CoM_coord
       if((use_new_bhns_initial_data.eq.3).or.(use_new_bhns_initial_data.eq.4)) then
          x_NS_CoM_coord2 = initial_ns2_coord_x
          y_NS_CoM_coord2 = initial_ns2_coord_y
          write(*,*) "INSIDE EMFIELDS SETUP: x,y of second ns:",x_NS_CoM_coord2,y_NS_CoM_coord2
       else
          x_NS_CoM_coord2 = x_NS_CoM_coord
          y_NS_CoM_coord2 = y_NS_CoM_coord
       end if
       write(*,*) "INSIDE EMFIELDS SETUP: Pmax:",bhns_P_max,p_c
       if (em_field_type==1 .and. Sym_Bz .gt. 0.d0) then
          write(*,*) "Sorry, em_field_type==1 (toroidal fields) not yet supported with Sym_Bz==1.  You can set Sym_Bz=-1 !"
          stop
       end if
       if (enable_trace_field_line==1 .and. em_field_type==1) then
          write(*,*) "Sorry, the initial data for the field line tracer variables currently only supports the poloidal initial data&
  &. You have to write your own initial data for the toroidal configuration."
       end if
       !
       ! Compute the vector potential A_phi 
       !
       if (em_field_type .ne. 1) then
          if(Symmetry==1 .and. angle_to_tilt_magnetic_fields.ne.0) then
             write(*,*) "CANNOT TILT THE MAGNETIC FIELDS WITH Symmetry==1."
             stop
          end if
          call BHNS_compute_Aphi(ext,X,Y,Z,PhysicalRadius,P,A_phi,Ax,Ay,Az,        &
               mhd_psi_line, mhd_u_psi, mhd_chi_line, mhd_u_chi,                   &
               betam1,p_c,x_NS_CoM_coord,y_NS_CoM_coord,bhns_P_max,Sym_Bz,         &
               enable_trace_field_line,constrained_transport_scheme,em_field_type, &
               angle_to_tilt_magnetic_fields,angle_to_tilt_magnetic_fieldsII,      &
               Aphi_power,rloop,Iloop,                                             &
               use_new_bhns_initial_data,x_NS_CoM_coord2,y_NS_CoM_coord2,          &
               bhns_P_max2,unequalmass,betam1_2,rloop2,Iloop2,conloop1,conloop2)
               !phi,Con_Pol,Con_Tor,superposition,r0)
       write(*,*) "after BHNS_compute_Aphi: x,y of ns:",x_NS_CoM_coord,y_NS_CoM_coord
       end if
       !
       ! Compute the vector potential Ax and Ay for near toroidal configuration
       !
       if (em_field_type==1) then
          if(constrained_transport_scheme==3) then
             call BHNS_compute_A_toroidal_constrained_transport_eq_3(ext,X,Y,Z,P,Ax,Ay,Az,phi,betam1, &
                  p_c,x_NS_CoM_coord,y_NS_CoM_coord,r0,bhns_P_max)
          else if(constrained_transport_scheme==1) then
             call BHNS_compute_A_toroidal(ext,X,Y,Z,P,Ax,Ay,Az,phi,betam1, &
                  p_c,x_NS_CoM_coord,y_NS_CoM_coord,r0,bhns_P_max)
          else
             write(*,*) "Sorry, em_field_type==1 (toroidal fields) only supported for constrained_transport_scheme == 1 or 3."
             stop
          end if
          write(*,*) "after BHNS_compute_A_toroidal: x,y of ns:",x_NS_CoM_coord,y_NS_CoM_coord
       end if
       !
       ! Compute the derivatives of A_phi, if we are not using constrained_transport_scheme==3
       !
       if (constrained_transport_scheme .ne. 3) then
          if (Symmetry==EQUATORIAL) then
             call gderivs_eq(ext,A_phi,A_phix,A_phiy,A_phiz,dX,dY,dZ,SYM,SYM ,Sym_Bz )
          else
             write(*,*) 'Symmetry type not supported in BHNS_setup_emfield'
             stop
          end if
       end if
       call BHNS_compute_Bi(ext,X,Y,Z,PhysicalRadius, &
            phi, &
            A_phi,A_phix,A_phiy,A_phiz, &
            Ax,Ay,Az, &
            Bx,By,Bz, &
            x_NS_CoM_coord,y_NS_CoM_coord,rho_b_atm, &
            constrained_transport_scheme,em_field_type)
       write(*,*) "after BHNS_compute_Bi: x,y of ns:",x_NS_CoM_coord,y_NS_CoM_coord
       do k=1,cctk_lsh(3)
          do j=1,cctk_lsh(2)
             do i=1,cctk_lsh(1)
                if(isnan(Ax(i,j,k)).eqv. .true.) then
                   write(*,*) "2. Inside driver_setup_emfields_part2-local-setBi.F90, Ax is nan!!!"
                end if
             end do
          end do
       end do
       ! Reset the atmospheric density to be variable and not constant if desired,
       ! falling off according to some prescription specified by parameter atm_type
       if (RESET_RHO_B_ATM.eq.1) then
          !        if (atm_type.eq.1) then
          call bhns_compute_b2_cpp(cctkGH,cctk_lsh, phi, lapm1, &
               shiftx,shifty,shiftz,vx,vy,vz,Bx,By,Bz, &
               gxx, gxy, gxz, gyy, gyz, gzz, temp2)
          !        end if
          call reset_rhob_atm_f90(cctk_lsh,X,Y,Z,temp2,rho_b,P,vx,vy,vz,shiftx,shifty,bhns_B_v,atm_type,&
               x_NS_CoM_coord,y_NS_CoM_coord,bhns_vx_CM,bhns_vy_CM,bhns_vz_CM,bhns_rhob_fac,rad_rhob_fac, &
               bhns_rhob_max,rhobatm_falloff_power,rhob_o_b2,bhns_R_NS,rho_b_atm,bhns_fac_atm,beta_ext,k_tab,gamma_tab,neos,x_NS_Co&
  &M_coord2,y_NS_CoM_coord2,use_new_bhns_initial_data)
       write(*,*) "after BHNS_reset_compute_b2_cpp: x,y of ns:",x_NS_CoM_coord,y_NS_CoM_coord
       end if
       temp2 = 0.D0
       temp3 = 0.D0
       if(1==0) then
          ! Compute B^i on staggered grid and temporarily store them in A_phix,A_phiy,A_phiz
          !
          if (constrained_transport_scheme==3) then
             do k=1,ext(3)
                do j=1,ext(2)
                   do i=1,ext(1)
                      im1 = max(i-1,1)
                      jm1 = max(j-1,1)
                      km1 = max(k-1,1)
                      ip1 = min(i+1,ext(1))
                      jp1 = min(j+1,ext(2))
                      kp1 = min(k+1,ext(3))
                      psim6_s = exp(-3.d0 * (phi(i,j,k) + phi(ip1,j,k)) )
                      A_phix(i,j,k) = ( (Az(i,j,k)-Az(i,jm1,k))/dY   &
                           - (Ay(i,j,k)-Ay(i,j,km1))/dZ ) * psim6_s
                      psim6_s = exp(-3.d0 * (phi(i,j,k) + phi(i,jp1,k)) )
                      A_phiy(i,j,k) = ( (Ax(i,j,k)-Ax(i,j,km1))/dZ &
                           - (Az(i,j,k)-Az(im1,j,k))/dX ) * psim6_s
                      psim6_s = exp(-3.d0 * (phi(i,j,k) + phi(i,j,kp1)) )
                      A_phiz(i,j,k) = ( (Ay(i,j,k)-Ay(im1,j,k))/dX &
                           - (Ax(i,j,k)-Ax(i,jm1,k))/dY ) * psim6_s
                   end do
                end do
             end do
             ! Now compute B^i on unstaggered grid by simple averge
             do k=1,ext(3)
                do j=1,ext(2)
                   do i=1,ext(1)
                      im1 = max(i-1,1)
                      jm1 = max(j-1,1)
                      km1 = max(k-1,1)
                      Bx(i,j,k) = 0.5d0* (A_phix(i,j,k) + A_phix(im1,j,k))
                      By(i,j,k) = 0.5d0* (A_phiy(i,j,k) + A_phiy(i,jm1,k))
                      Bz(i,j,k) = 0.5d0* (A_phiz(i,j,k) + A_phiz(i,j,km1))
                   end do
                end do
             end do
          end if
          do k = kmin,kmax
             do j = jmin,jmax
                do i = imin,imax
                   xn = X(i,j,k) - x_NS_CoM_coord
                   yn = Y(i,j,k) - y_NS_CoM_coord
                   psim6 = exp(-6.d0*phi(i,j,k))
                   pomega2 = xn**2 + yn**2
                   if (constrained_transport_scheme==3) then
                      ! do nothing since B^i has been computed 
                   elseif (em_field_type==0 .or. i==imax .or. j==jmax .or. k==kmax) then
                      Bx(i,j,k) = -xn/pomega2 * psim6 * A_phiz(i,j,k)
                      By(i,j,k) = -yn/pomega2 * psim6 * A_phiz(i,j,k)
                      Bz(i,j,k) = psim6/pomega2 * (xn*A_phix(i,j,k) +  &
                           yn*A_phiy(i,j,k))
                   else
                      if(i==imin .and. j==jmin .and. k==kmin) write(*,*) "CONSTRAINED TRANSPORT SCHEME != 3. = ",constrained_transp&
  &ort_scheme
                      ! Compute Bx, By, Bz that satisfy div(B)=0 to machine precision.
                      ! Here is the recipe, which has been verified by Mathematica: 
                      !
                      ! \tilde{B}^x(i,j,k)=-(Ay(i,j,k+1)-Ay(i,j,k) + Ay(i+1,j,k+1)-Ay(i+1,j,k) + Ay(i,j+1,k+1)-Ay(i,j+1,k) + Ay(i+1,
                      ! \tilde{B}^y(i,j,k)=(Ax(i,j,k+1)-Ax(i,j,k) + Ax(i+1,j,k+1)-Ax(i+1,j,k) + Ax(i,j+1,k+1)-Ax(i,j+1,k) + Ax(i+1,j
                      ! \tilde{B}^z(i,j,k):=(Ay(i+1,j,k)-Ay(i,j,k) + Ay(i+1,j+1,k)-Ay(i,j+1,k) + Ay(i+1,j,k+1)-Ay(i,j,k+1) + Ay(i+1,
                      ! Here \tilde{B}^i = psi^6 B^i, and 
                      ! Ax(i,j,k), Ay(i,j,k), and Az(i,j,k) are the 
                      ! 3 (covariant) components of the vector potential at the corner 
                      ! point (x_i-dx/2, y_j-dy/2, z_k-dz/2).
                      ! In the present case, we set Ax = -A_phi * y/(x^2+y^2), 
                      !  Ay = A_phi * x/(x^2+y^2), and Az = 0.
                      ! Note that B^i cant be set this way at the boundary points
                      ! becuase of the sentcil structure. Presumably, this 
                      ! is not a problem when we have enough ghostzones.
                      !
                      ! In the following, Axabc denotes Ax(i+a,j+b,k+c) and 
                      ! similarly for Ayabc.
                      ! 
                      ! Here we set the coordinates, relative to the origin, which in this case is the center of the NS.
                      Yijk      = yn
                      Yijkp1    = yn
                      Yijp1k    = yn + dY
                      Yijp1kp1  = yn + dY
                      Yip1jk    = yn
                      Yip1jkp1  = yn
                      Yip1jp1k  = yn + dY
                      Yip1jp1kp1= yn + dY
                      Xijk      = xn
                      Xijkp1    = xn
                      Xijp1k    = xn
                      Xijp1kp1  = xn
                      Xip1jk    = xn + dX
                      Xip1jkp1  = xn + dX
                      Xip1jp1k  = xn + dX
                      Xip1jp1kp1= xn + dX
                      Ax000 = -A_phi(i,j,k)* (Yijk-0.5d0*dY)/ &
                           ((Xijk-0.5d0*dX)**2 + (Yijk-0.5d0*dY)**2)
                      Ax001 = -A_phi(i,j,k+1)*(Yijkp1-0.5d0*dY)/ &
                           ((Xijkp1-0.5d0*dX)**2 + (Yijkp1-0.5d0*dY)**2)
                      Ax010 = -A_phi(i,j+1,k)* (Yijp1k-0.5d0*dY)/ &
                           ((Xijp1k-0.5d0*dX)**2 + (Yijp1k-0.5d0*dY)**2)
                      Ax011 = -A_phi(i,j+1,k+1)*(Yijp1kp1-0.5d0*dY)/ &
                           ((Xijp1kp1-0.5d0*dX)**2 + (Yijp1kp1-0.5d0*dY)**2)
                      Ax100 = -A_phi(i+1,j,k)*(Yip1jk-0.5d0*dY)/ &
                           ((Xip1jk-0.5d0*dX)**2 + (Yip1jk-0.5d0*dY)**2)
                      Ax101 = -A_phi(i+1,j,k+1)*(Yip1jkp1-0.5d0*dY)/ &
                           ((Xip1jkp1-0.5d0*dX)**2 + (Yip1jkp1-0.5d0*dY)**2)
                      Ax110 = -A_phi(i+1,j+1,k)*(Yip1jp1k-0.5d0*dY)/ &
                           ((Xip1jp1k-0.5d0*dX)**2 + (Yip1jp1k-0.5d0*dY)**2)
                      Ax111 = -A_phi(i+1,j+1,k+1)*(Yip1jp1kp1-0.5d0*dY)/ &
                           ((Xip1jp1kp1-0.5d0*dX)**2 + (Yip1jp1kp1-0.5d0*dY)**2)
                      Ay000 = A_phi(i,j,k)* (Xijk-0.5d0*dX)/ &
                           ((Xijk-0.5d0*dX)**2 + (Yijk-0.5d0*dY)**2)
                      Ay001 = A_phi(i,j,k+1)*(Xijkp1-0.5d0*dX)/ &
                           ((Xijkp1-0.5d0*dX)**2 + (Yijkp1-0.5d0*dY)**2)
                      Ay010 = A_phi(i,j+1,k)* (Xijp1k-0.5d0*dX)/ &
                           ((Xijp1k-0.5d0*dX)**2 + (Yijp1k-0.5d0*dY)**2)
                      Ay011 = A_phi(i,j+1,k+1)*(Xijp1kp1-0.5d0*dX)/ &
                           ((Xijp1kp1-0.5d0*dX)**2 + (Yijp1kp1-0.5d0*dY)**2)
                      Ay100 = A_phi(i+1,j,k)*(Xip1jk-0.5d0*dX)/ &
                           ((Xip1jk-0.5d0*dX)**2 + (Yip1jk-0.5d0*dY)**2)
                      Ay101 = A_phi(i+1,j,k+1)*(Xip1jkp1-0.5d0*dX)/ &
                           ((Xip1jkp1-0.5d0*dX)**2 + (Yip1jkp1-0.5d0*dY)**2)
                      Ay110 = A_phi(i+1,j+1,k)*(Xip1jp1k-0.5d0*dX)/ &
                           ((Xip1jp1k-0.5d0*dX)**2 + (Yip1jp1k-0.5d0*dY)**2)
                      Ay111 = A_phi(i+1,j+1,k+1)*(Xip1jp1kp1-0.5d0*dX)/ &
                           ((Xip1jp1kp1-0.5d0*dX)**2 + (Yip1jp1kp1-0.5d0*dY)**2)
                      Bx(i,j,k) = -( (Ay001-Ay000) + (Ay101-Ay100) + (Ay011-Ay010) + (Ay111-Ay110) ) * 0.25d0/dZ * psim6
                      By(i,j,k) = ( (Ax001-Ax000) + (Ax101-Ax100) + (Ax011-Ax010) + (Ax111-Ax110) ) * 0.25d0/dZ * psim6
                      Bz(i,j,k) = ( (Ay100-Ay000) + (Ay110-Ay010) + (Ay101-Ay001) + (Ay111-Ay011) ) * 0.25d0/dX * psim6 - ( (Ax010-&
  &Ax000) + (Ax110-Ax100) + (Ax011-Ax001) + (Ax111-Ax101) ) * 0.25d0/dY * psim6
                   end if
                   if(1==0) then
                      ! Compute b^0 and b_i
                      al = 1.d0 + lapm1(i,j,k)
                      sqrtg = 1.d0/psim6
                      sqrtg4 = al * sqrtg
                      B2s = exp(4.d0*phi(i,j,k))*(gxx(i,j,k)*Bx(i,j,k)**2 + &
                           2.d0*gxy(i,j,k)*Bx(i,j,k)*By(i,j,k) + &
                           2.d0*gxz(i,j,k)*Bx(i,j,k)*Bz(i,j,k) + &
                           gyy(i,j,k)*By(i,j,k)**2 + 2.d0*gyz(i,j,k)*By(i,j,k)*Bz(i,j,k) + &
                           gzz(i,j,k)*Bz(i,j,k)**2)/(fs4pi*al)**2
                      psin = exp(4.d0*phi(i,j,k))/al/fs4pi
                      B_xs  = psin * (gxx(i,j,k) * Bx(i,j,k) + gxy(i,j,k) * By(i,j,k) + &
                           gxz(i,j,k) * Bz(i,j,k))
                      B_ys  = psin * (gxy(i,j,k) * Bx(i,j,k) + gyy(i,j,k) * By(i,j,k) + &
                           gyz(i,j,k) * Bz(i,j,k))
                      B_zs  = psin * (gxz(i,j,k) * Bx(i,j,k) + gyz(i,j,k) * By(i,j,k) + &
                           gzz(i,j,k) * Bz(i,j,k))
                      psin = psi4*u0(i,j,k)
                      u_x = ( gxx(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) +  &
                           gxy(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) +  &
                           gxz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) )*psin
                      u_y = ( gxy(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) +  &
                           gyy(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) +  &
                           gyz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) )*psin
                      u_z = ( gxz(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) +  &
                           gyz(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) +  &
                           gzz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) )*psin
                      sb0 = (u_x*Bx(i,j,k) + u_y*By(i,j,k) + &
                           u_z*Bz(i,j,k))/fs4pi/al
                      sb2 = (B2s + sb0**2)/u0(i,j,k)**2
                      !Bransons way of ensuring b^2/P is what we want: 
                      !  Call this routine once to calibrate sb2/P, then average it.
                      !   Scale sb2 accordingly, then call this function again to obtain the desired sb2/P
                      ! Heres the necessary line of code.  Youll need to add the input parameters as well.
                      !bsq(i,j,k) = sb2/P(i,j,k)
                      sb_x = (B_xs + u_x*sb0)/u0(i,j,k)
                      sb_y = (B_ys + u_y*sb0)/u0(i,j,k)
                      sb_z = (B_zs + u_z*sb0)/u0(i,j,k)
                      ! Now compute mhd_st_i and tau
                      mhd_st_x(i,j,k) = st_x(i,j,k)+sqrtg4*(u0(i,j,k)*sb2*u_x-sb0*sb_x)
                      mhd_st_y(i,j,k) = st_y(i,j,k)+sqrtg4*(u0(i,j,k)*sb2*u_y-sb0*sb_y)
                      mhd_st_z(i,j,k) = st_z(i,j,k)+sqrtg4*(u0(i,j,k)*sb2*u_z-sb0*sb_z)
                      if(i==43 .and. j==40 .and. k==5) then
                         write(*,*) "hey! setting tau:",tau(i,j,k)
                      end if
                      tau(i,j,k) = tau(i,j,k) + sqrtg*( sb2*(al*u0(i,j,k))**2 &
                           - sb2*0.5d0 - (al*sb0)**2 )
                      if(i==43 .and. j==40 .and. k==5) then
                         write(*,*) "hey! setting tau:",tau(i,j,k)
                      end if
                      ! CHECK NEXT FUNCTION IN SCHEDULE.CCL for reason why temp1 and temp2 are defined.
                      ! Temporarily store e^(6 phi)*b^2/2 (b^2/2 = P_mag) to temp1, and e^(6 phi)*P(i,j,k) to temp2, e^(6 phi) to te
                      temp1(i,j,k) = sb2*0.5d0*sqrtg
                      !if(sb2.gt.0.D0) write(*,*) sb2
                      temp2(i,j,k) = P(i,j,k)*sqrtg
                      if (rho_b(i,j,k) .gt. rho_b_atm*1.d5) then
                         temp3(i,j,k) = 1.d0
                      else
                         temp3(i,j,k) = 0.d0
                      end if
                   end if
                end do
             end do
          end do
       end if
       deallocate(A_phi, A_phix, A_phiy, A_phiz)
    end if
 end subroutine BHNS_setup_emfield_part2_local_setBi
  !------------------------------------------------------------------------
  ! Compute A_phi = pm^2 * sqrt(8 pi bhns_P_max * betam1) * max( P/bhns_P_max - p_c, 0) for Sym_Bz = 1  and 
  !         A_phi = pm^2 * z/r * sqrt(8 pi bhns_P_max * betam1) * max( P/bhns_P_max - p_c, 0) for Sym_Bz = -1 
  !   A_phi = (x^2+y^2) * sqrt(8 pi P_max * betam1) * 
  !                max( P/P_max - p_c, 0) 
  ! where pm^2 = (x-x_ns_com)^2 + (y-y_ns_com)^2 
  ! x_ns_com is the position of the center of the neutron star.
  !------------------------------------------------------------------------
  subroutine BHNS_compute_Aphi(ext,X,Y,Z,PhysR,P,A_phi,Ax,Ay,Az,           &
       mhd_psi_line, mhd_u_psi, mhd_chi_line, mhd_u_chi,                   &
       betam1,p_c,x_NS_CoM_coord,y_NS_CoM_coord,bhns_P_max,Sym_Bz,         &
       enable_trace_field_line,constrained_transport_scheme,em_field_type, &
       angle_to_tilt_magnetic_fields,angle_to_tilt_magnetic_fieldsII,      &
       Aphi_power,a0,I0,two_ns,                                            &
       x_NS_CoM_coord2,y_NS_CoM_coord2,                                    &
       bhns_P_max2,unequalmass,betam1_2,a02,I02,con1,con2)
    implicit none
    integer, dimension(3)                          :: ext
    real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z,PhysR
    real*8, dimension(ext(1),ext(2),ext(3))        :: A_phi, P, Ax,Ay,Az
    real*8, dimension(ext(1),ext(2),ext(3))        :: mhd_psi_line, mhd_u_psi, mhd_chi_line, mhd_u_chi
    real*8                                         :: x_NS_CoM_coord,y_NS_CoM_coord
    real*8                                         :: x_NS_CoM_coord2, y_NS_CoM_coord2, pomega2b
    real*8                                           :: betam1,betam1_2,p_c,fac,fac2,pomega2,a0,I0,a02,I02,con1,con2
    real*8                                         :: bhns_P_max,bhns_P_max2,xp,yp,zp,r
    real*8                                            :: xp_s,yp_s,zp_s,Aphi_s,P_s,phi
    integer                                           :: i,j,k,em_field_type,two_ns,unequalmass
    integer                                        :: fisheye_enable,ip1,jp1,kp1,im1,jm1,km1
    integer                                        :: constrained_transport_scheme
    integer                                        :: enable_trace_field_line,tilt_star1,tilt_star2
    real*8                                         :: Sym_Bz,P_corner,hdX,hdY,hdZ, sintheta,costheta, Ax_temp, Ay_temp, Az_temp, x_&
  &before_tilt
    real*8                                            :: angle_to_tilt_magnetic_fields,angle_to_tilt_magnetic_fieldsII,Aphi_power,a&
  &0tmp,I0tmp,bhns_P_max_tmp,fac_tmp
!    real*8                                            :: Ax_tmpI, Ay_tmpI, Az_tmpI, Ax_tmpII, Ay_tmpII, Az_tmpII
    !  real*8, parameter :: Aphi_power = 1.d0
    !
    ! fac = sqrt(8.d0*acos(-1.d0)*bhns_P_max*betam1)/Riso**2
    fac = sqrt(8.d0*acos(-1.d0)*bhns_P_max*betam1)
    if(unequalmass.eq.1) then
       fac2 = sqrt(8.d0*acos(-1.d0)*bhns_P_max2*betam1_2)
    end if
    print *, "bhns_P_max2,bhns_P_max",bhns_P_max2,bhns_P_max
    print *, "fac2,a02,I02",fac2,a02,I02
    print *, "unequalmass,betam1_2",unequalmass,betam1_2
    print *, "Angle tilt",angle_to_tilt_magnetic_fields,angle_to_tilt_magnetic_fieldsII
    !n_b = 2.d0/3.d0
    hdX = 0.5d0*(X(2,1,1) - X(1,1,1))
    hdY = 0.5d0*(Y(1,2,1) - Y(1,1,1))
    hdZ = 0.5d0*(Z(1,1,2) - Z(1,1,1))
    if ( abs(sqrt(X(1,1,1)**2 + Y(1,1,1)**2 + Z(1,1,1)**2)-PhysR(1,1,1)) .gt. hdX*1.d-3) then
       fisheye_enable = 1
    else
       fisheye_enable = 0
    end if
    if (fisheye_enable==1 .and. (em_field_type==2 .or. constrained_transport_scheme==3)) then
       write(*,*) 'em_field_type=2 and constrained_transport_scheme=3 does not support fisheye (yet).'
       stop
    end if
    do k = 1,ext(3)
       do j = 1,ext(2)
          do i = 1,ext(1)
             ip1 = min(i+1,ext(1))
             jp1 = min(j+1,ext(2))
             kp1 = min(k+1,ext(3))
             im1 = max(i-1,1)
             jm1 = max(j-1,1)
             km1 = max(k-1,1)
             r = sqrt(X(i,1,1)**2 + Y(1,j,1)**2 + Z(1,1,k)**2)
             xp = X(i,1,1)/r * PhysR(i,j,k)
             yp = Y(1,j,1)/r * PhysR(i,j,k)
             !!pomega2 = xp**2 + yp**2 + hdX*1.d-13
             pomega2 = (xp-x_NS_CoM_coord)**2 + (yp-y_NS_CoM_coord)**2 + hdX*1.d-13
             if (constrained_transport_scheme==3) then
                ! setup A_i on staggered grid
                xp_s = X(i,1,1)+hdX
                yp_s = Y(1,j,1)+hdY
                zp_s = Z(1,1,k)+hdZ
                !x-component of A_i :
                xp = X(i,1,1)
                yp = yp_s
                zp = zp_s
                !!pomega2 = xp**2 + yp**2 + hdX*1.d-13        
                pomega2 = (xp-x_NS_CoM_coord)**2 + (yp-y_NS_CoM_coord)**2 + hdX*1.d-13
                if((two_ns.eq.3).or.(two_ns.eq.4)) then
                   pomega2b = (xp-x_NS_CoM_coord2)**2 + (yp-y_NS_CoM_coord2)**2 + hdX*1.d-13
                end if
                P_s = 0.25d0*( P(i,j,k) + P(i,jp1,k) + P(i,j,kp1) + P(i,jp1,kp1) )
                if (em_field_type.eq.0) then
                   if((two_ns.eq.3).or.(two_ns.eq.4)) then
                      ! We need to choose which star we are inside (If were not inside, Aphi_s should be zero anyway.)
                      if (pomega2 < pomega2b) then
                         if(unequalmass.eq.1) then
                            fac_tmp=fac2
                            bhns_P_max_tmp=bhns_P_max2
                         else
                            fac_tmp=fac
                            bhns_P_max_tmp=bhns_P_max
                         end if
                         Aphi_s = pomega2 * fac_tmp * max(P_s/bhns_P_max_tmp - p_c, 0.d0)**Aphi_power
                         if (Sym_Bz .lt. 0.d0) Aphi_s = Aphi_s*zp/sqrt(pomega2)
                         Ax(i,j,k) = -(yp-y_NS_CoM_coord)/pomega2*Aphi_s
                      else
                         Aphi_s = pomega2b * fac * max(P_s/bhns_P_max - p_c, 0.d0)**Aphi_power
                         if (Sym_Bz .lt. 0.d0) Aphi_s = Aphi_s*zp/sqrt(pomega2b)
                         Ax(i,j,k) = -(yp-y_NS_CoM_coord2)/pomega2b*Aphi_s
                      end if
                   else
                      Aphi_s = pomega2 * fac * max(P_s/bhns_P_max - p_c, 0.d0)**Aphi_power
                      if (Sym_Bz .lt. 0.d0) Aphi_s = Aphi_s*zp/sqrt(pomega2)
                      Ax(i,j,k) = -(yp-y_NS_CoM_coord)/pomega2*Aphi_s
                   end if
                else if (em_field_type.eq.3) then
                   if(unequalmass.eq.1) then
                      a0tmp=a02
                      I0tmp=I02
                   else
                      a0tmp=a0
                      I0tmp=I0
                   end if
                   call compute_Aphi_bhns_current_loop_app(xp,yp,zp,x_NS_CoM_coord,y_NS_CoM_coord,a0tmp,I0tmp,Aphi_s)
                   if (Sym_Bz .lt. 0.d0) Aphi_s = Aphi_s*zp/sqrt(pomega2)
                   Ax(i,j,k) = -(yp-y_NS_CoM_coord)/pomega2*Aphi_s
                   if((two_ns.eq.3).or.(two_ns.eq.4)) then
                      call compute_Aphi_bhns_current_loop_app(xp,yp,zp,x_NS_CoM_coord2,y_NS_CoM_coord2,a0,I0,Aphi_s)
                      if (Sym_Bz .lt. 0.d0) Aphi_s = Aphi_s*zp/sqrt(pomega2b)
                      Ax(i,j,k) = Ax(i,j,k) - (yp-y_NS_CoM_coord2)/pomega2b*Aphi_s
                   end if
                else if (em_field_type.eq.4) then
                   call compute_Aphi_bhns_current_loop_app(xp,yp,zp,x_NS_CoM_coord,y_NS_CoM_coord,a0,I0,Aphi_s)
                   Aphi_s = Aphi_s*(1.d0-exp(-(r*r/(con2*con2))**con1))&
                        +(pomega2 * fac * max(P_s/bhns_P_max - p_c, 0.d0)**Aphi_power)*exp(-(r*r/(con2*con2))**con1)
                   if (Sym_Bz .lt. 0.d0) Aphi_s = Aphi_s*zp/sqrt(pomega2)
                   Ax(i,j,k) = -(yp-y_NS_CoM_coord)/pomega2*Aphi_s
                   if(two_ns == 3) then
                      write(*,*) "Stopping: bhns_initial_data=3 and em_field_type=4 is not supported at the moment"
                      exit
                   end if
                end if
                !y-component of A_i :
                xp = xp_s
                yp = Y(1,j,1)
                zp = zp_s
                !!pomega2 = xp**2 + yp**2 + hdX*1.d-13
                pomega2 = (xp-x_NS_CoM_coord)**2 + (yp-y_NS_CoM_coord)**2 + hdX*1.d-13
                if((two_ns.eq.3).or.(two_ns.eq.4)) then
                   pomega2b = (xp-x_NS_CoM_coord2)**2 + (yp-y_NS_CoM_coord2)**2 + hdX*1.d-13
                end if
                P_s = 0.25d0*( P(i,j,k) + P(ip1,j,k) + P(i,j,kp1) + P(ip1,j,kp1) )
                if (em_field_type.eq.0) then
                   if((two_ns.eq.3).or.(two_ns.eq.4)) then
                      ! We need to choose which star we are inside (If were not inside, Aphi_s should be zero anyway.)
                      if (pomega2 < pomega2b) then
                         if(unequalmass.eq.1) then
                            fac_tmp=fac2
                            bhns_P_max_tmp=bhns_P_max2
                         else
                            fac_tmp=fac
                            bhns_P_max_tmp=bhns_P_max
                         end if
                         Aphi_s = pomega2 * fac_tmp * max(P_s/bhns_P_max_tmp - p_c, 0.d0)**Aphi_power
                         if (Sym_Bz .lt. 0.d0) Aphi_s = Aphi_s*zp/sqrt(pomega2)
                         Ay(i,j,k) = (xp-x_NS_CoM_coord)/pomega2*Aphi_s
                      else
                         Aphi_s = pomega2b * fac * max(P_s/bhns_P_max - p_c, 0.d0)**Aphi_power
                         if (Sym_Bz .lt. 0.d0) Aphi_s = Aphi_s*zp/sqrt(pomega2b)
                         Ay(i,j,k) = (xp-x_NS_CoM_coord2)/pomega2b*Aphi_s
                      end if
                   else
                      Aphi_s = pomega2 * fac * max(P_s/bhns_P_max - p_c, 0.d0)**Aphi_power
                      if (Sym_Bz .lt. 0.d0) Aphi_s = Aphi_s*zp/sqrt(pomega2)
                      Ay(i,j,k) = (xp-x_NS_CoM_coord)/pomega2*Aphi_s
                   end if
                else if (em_field_type.eq.3) then
                   if(unequalmass.eq.1) then
                      a0tmp=a02
                      I0tmp=I02
                   else
                      a0tmp=a0
                      I0tmp=I0
                   end if
                   call compute_Aphi_bhns_current_loop_app(xp,yp,zp,x_NS_CoM_coord,y_NS_CoM_coord,a0tmp,I0tmp,Aphi_s)
                   if (Sym_Bz .lt. 0.d0) Aphi_s = Aphi_s*zp/sqrt(pomega2)
                   Ay(i,j,k) = (xp-x_NS_CoM_coord)/pomega2*Aphi_s
                   if((two_ns.eq.3).or.(two_ns.eq.4)) then
                      call compute_Aphi_bhns_current_loop_app(xp,yp,zp,x_NS_CoM_coord2,y_NS_CoM_coord2,a0,I0,Aphi_s)
                      if (Sym_Bz .lt. 0.d0) Aphi_s = Aphi_s*zp/sqrt(pomega2b)
                      Ay(i,j,k) = Ay(i,j,k) + (xp-x_NS_CoM_coord2)/pomega2b*Aphi_s
                   end if
                else if (em_field_type.eq.4) then
                   call compute_Aphi_bhns_current_loop_app(xp,yp,zp,x_NS_CoM_coord,y_NS_CoM_coord,a0,I0,Aphi_s)
                   Aphi_s = Aphi_s*(1.d0-exp(-(r*r/(con2*con2))**con1))&
                         +(pomega2 * fac * max(P_s/bhns_P_max - p_c, 0.d0)**Aphi_power)*exp(-(r*r/(con2*con2))**con1)
                   if (Sym_Bz .lt. 0.d0) Aphi_s = Aphi_s*zp/sqrt(pomega2)
                   Ay(i,j,k) = (xp-x_NS_CoM_coord)/pomega2*Aphi_s
                   if(two_ns == 3) then
                      write(*,*) "Stopping: bhns_initial_data=3 and em_field_type=4 is not supported at the moment"
                      exit
                   end if
                end if
                !z-component of A_i :
                Az(i,j,k) = 0.d0
                !****************************************
                !*** do you want to tilt the Bfield?  ***
                !****************************************
                if((angle_to_tilt_magnetic_fields.gt.0.D0).or. &
                      (angle_to_tilt_magnetic_fieldsII.gt.0.D0)) then
                   ! reset values
                   Ax(i,j,k) = 0.0d0
                   Ay(i,j,k) = 0.0d0
                   Az(i,j,k) = 0.0d0
                   r = sqrt(X(i,1,1)**2 + Y(1,j,1)**2 + Z(1,1,k)**2)
                   xp = X(i,1,1)/r * PhysR(i,j,k)
                   yp = Y(1,j,1)/r * PhysR(i,j,k)
                   xp_s = X(i,1,1)+hdX
                   yp_s = Y(1,j,1)+hdY
                   zp_s = Z(1,1,k)+hdZ
                   pomega2  = (zp_s)**2   + (yp-y_NS_CoM_coord )**2 + hdX*1.d-13
                   pomega2b = (zp_s)**2   + (yp-y_NS_CoM_coord2)**2 + hdX*1.d-13
                   !x-component of A_i :
                   xp = xp_s
                   yp = yp_s
                   zp = Z(1,1,k)
                   call compute_Aphi_bhns_current_loop_app(zp,yp,xp-x_NS_CoM_coord, 0.0,y_NS_CoM_coord, a0,I0,Aphi_s)
                   Az_temp =         - (yp - y_NS_CoM_coord )*Aphi_s/pomega2
                   call compute_Aphi_bhns_current_loop_app(zp,yp,xp-x_NS_CoM_coord2,0.0, y_NS_CoM_coord2,a0,I0,Aphi_s)
                   Az_temp = Az_temp - (yp - y_NS_CoM_coord2)*Aphi_s/pomega2b
                   !y-component of A_i :
                   xp = xp_s
                   yp = Y(1,j,1)
                   zp = zp_s
                   call compute_Aphi_bhns_current_loop_app(zp,yp,xp-x_NS_CoM_coord, 0.0,y_NS_CoM_coord, a0,I0,Aphi_s)
                   Ay_temp =          zp*Aphi_s/pomega2
                   call compute_Aphi_bhns_current_loop_app(zp,yp,xp-x_NS_CoM_coord2, 0.0,y_NS_CoM_coord2,a0,I0,Aphi_s)
                   Ay_temp= Ay_temp + zp*Aphi_s/pomega2b
                   !z-component of A_i :
                   Ax_temp = 0.d0
                   Ax(i,j,k) =  Ax_temp
                   Ay(i,j,k) =  Ay_temp
                   Az(i,j,k) =  Az_temp
                end if
                A_phi(i,j,k) = 0.d0
                if (enable_trace_field_line==1) then
                   pomega2 = (X(i,j,k)-x_NS_CoM_coord)**2 + (Y(i,j,k)-y_NS_CoM_coord)**2 + hdX*1.d-13
                   if (em_field_type.eq.0) then
                      !Aphi_s = pomega2 * fac * max(P_s/bhns_P_max - p_c, 0.d0)**1
                      Aphi_s = pomega2 * fac * max((P(i,j,k)/bhns_P_max - p_c), 0.d0)**Aphi_power
                   else if (em_field_type.eq.3) then
                      call compute_Aphi_bhns_current_loop_app(xp,yp,zp,x_NS_CoM_coord,y_NS_CoM_coord,a0,I0,Aphi_s)
                   end if
                   phi = atan2(Y(i,j,k)-y_NS_CoM_coord,X(i,j,k)-x_NS_CoM_coord)
                   mhd_psi_line(i,j,k) = Aphi_s*cos(phi)
                   mhd_chi_line(i,j,k) = Aphi_s*sin(phi)
                   mhd_u_psi(i,j,k) = 0.d0
                   mhd_u_chi(i,j,k) = 0.d0
                end if
             elseif (em_field_type==0 .or. i==1 .or. j==1 .or. k==1) then
                !A_phi(i,j,k) = pomega2 * fac * max(P(i,j,k)/bhns_P_max - p_c, 0.d0)**1
                A_phi(i,j,k) = pomega2 * fac * max(P(i,j,k)/bhns_P_max - p_c, 0.d0)**Aphi_power
                if (enable_trace_field_line==1) then
                   phi = atan2(Y(i,j,k)-y_NS_CoM_coord,X(i,j,k)-x_NS_CoM_coord)
                   mhd_psi_line(i,j,k) = A_phi(i,j,k)*cos(phi)
                   mhd_chi_line(i,j,k) = A_phi(i,j,k)*sin(phi)
                   mhd_u_psi(i,j,k) = 0.d0
                   mhd_u_chi(i,j,k) = 0.d0
                end if
                ! Bransons alternate ways of setting Aphi.  n_b sets the power-law dependence on input parameter P.
                ! Alt way 1: A_phi(i,j,k) = pomega2 * fac * max(P(i,j,k)**(1.d0/n_b) - (bhns_P_max*p_c)**(1.d0/n_b), 0.d0)
                ! Alt way 2: A_phi(i,j,k) = pomega2 * fac * max((P(i,j,k)/(p_c*bhns_P_max))**(1.d0/n_b) - 1.d0, 0.d0)
                if (Sym_Bz .lt. 0.d0) A_phi(i,j,k) = A_phi(i,j,k)*Z(1,1,k)/r*PhysR(i,j,k)/sqrt(pomega2)
             else
                ! Note that for this initial data thorn, A_phi(i,j,k) 
                ! stores A_phi at the corner point
                ! x_i - dx/2, y_j-dy/2, z_j-dz/2
                ! We can only do that for i>1, j>1 and k>1. I assume 
                ! this will be fine as long as we have enough ghostzones.
                P_corner = 0.125d0*(P(i,j,k) + P(im1,j,k) +  &
                     P(i,jm1,k) + P(im1,jm1,k) + &
                     P(i,j,km1) + P(im1,j,km1) + &
                     P(i,jm1,km1) + P(im1,jm1,km1) )
                xp = X(i,j,k) - hdX
                yp = Y(i,j,k) - hdY
                zp = Z(i,j,k) - hdZ
                !!pomega2 = xp**2 + yp**2 + hdX*1.d-13
                pomega2 = (xp-x_NS_CoM_coord)**2 + (yp-y_NS_CoM_coord)**2 + hdX*1.d-13
                !A_phi(i,j,k) = pomega2 * fac * max(P_corner/bhns_P_max - p_c, 0.d0)**1
                A_phi(i,j,k) = pomega2 * fac * max(P_corner/bhns_P_max - p_c, 0.d0)**Aphi_power
                if (Sym_Bz .lt. 0.d0) A_phi(i,j,k) = A_phi(i,j,k)*zp/sqrt(pomega2)
                if (enable_trace_field_line==1) then
                   pomega2 = (X(i,j,k)-x_NS_CoM_coord)**2 + (Y(i,j,k)-y_NS_CoM_coord)**2 + hdX*1.d-13
                   !Aphi_s = pomega2 * fac * max((P(i,j,k)/bhns_P_max - p_c), 0.d0)**1
                   Aphi_s = pomega2 * fac * max((P(i,j,k)/bhns_P_max - p_c), 0.d0)**Aphi_power
                   phi = atan2(Y(i,j,k)-y_NS_CoM_coord,X(i,j,k)-x_NS_CoM_coord)
                   mhd_psi_line(i,j,k) = Aphi_s*cos(phi)
                   mhd_chi_line(i,j,k) = Aphi_s*sin(phi)
                   mhd_u_psi(i,j,k) = 0.d0
                   mhd_u_chi(i,j,k) = 0.d0
                end if
             end if
          end do
       end do
    end do
 end subroutine BHNS_compute_Aphi
  ! Set up vector potential for toroidal B field according to 
  ! Ax = ((x-x_ns)z/r0) sqrt[ 8 pi betam1 psi^4 max(P - p_c Pmax , 0) ]
  ! Ay = ((y-y_ns)z/r0) sqrt[ 8 pi betam1 psi^4 max(P - p_c Pmax , 0) ]
  ! Az = 0
  subroutine BHNS_compute_A_toroidal(ext,X,Y,Z,P,Ax,Ay,Az,phi,betam1, &
       p_c,x_NS_CoM_coord,y_NS_CoM_coord,r0,bhns_P_max)
    implicit none
    integer, dimension(3) :: ext
    real*8, dimension(ext(1),ext(2),ext(3)) :: X,Y,Z,P,Ax,Ay,Az,phi
    real*8 :: betam1,p_c,x_NS_CoM_coord,y_NS_CoM_coord,bhns_P_max,r0
    real*8 :: zs,xi,yj,hdx,hdy,hdz,psi4s,Ps,Ps_Pcut
    integer :: i,j,k,ip1,jp1,kp1
    real*8, parameter :: f8pi = 25.13274122871834590768d0
    !
    hdX = 0.5d0*(X(2,1,1) - X(1,1,1))
    hdY = 0.5d0*(Y(1,2,1) - Y(1,1,1))
    hdZ = 0.5d0*(Z(1,1,2) - Z(1,1,1))
    Az = 0.d0
    do k=1,ext(3)
       do j=1,ext(2)
          do i=1,ext(1)
             xi = X(i,j,k) - x_NS_CoM_coord
             yj = Y(i,j,k) - y_NS_CoM_coord
             zs = Z(i,j,k) + hdz
             ip1 = min(i+1,ext(1))
             jp1 = min(j+1,ext(2))
             kp1 = min(k+1,ext(3))
             Ps = 0.25d0*(P(i,j,k) + P(i,jp1,k) + P(i,j,kp1) + P(i,jp1,kp1) )
             Ps_Pcut = max(Ps - p_c*bhns_P_max, 0.d0)
             psi4s = exp( phi(i,j,k) + phi(i,jp1,k) + phi(i,j,kp1) + phi(i,jp1,kp1) )
             Ax(i,j,k) = xi*zs/r0*sqrt(f8pi*betam1*psi4s*Ps_Pcut)
             Ps = 0.25d0*(P(i,j,k) + P(ip1,j,k) + P(i,j,kp1) + P(ip1,j,kp1) )
             Ps_Pcut = max(Ps - p_c*bhns_P_max, 0.d0)
             psi4s = exp( phi(i,j,k) + phi(ip1,j,k) + phi(i,j,kp1) + phi(ip1,j,kp1) )
             Ay(i,j,k) = yj*zs/r0*sqrt(f8pi*betam1*psi4s*Ps_Pcut)
          end do
       end do
    end do
  end subroutine BHNS_compute_A_toroidal
  ! Set up vector potential for toroidal B field according to 
  ! Ax = ((x-x_ns)z/r0) sqrt[ 8 pi betam1 psi^4 max(P - p_c Pmax , 0) ]
  ! Ay = ((y-y_ns)z/r0) sqrt[ 8 pi betam1 psi^4 max(P - p_c Pmax , 0) ]
  ! Az = 0
  subroutine BHNS_compute_A_toroidal_constrained_transport_eq_3(ext,X,Y,Z,P,Ax,Ay,Az,phi,betam1, &
       p_c,x_NS_CoM_coord,y_NS_CoM_coord,r0,bhns_P_max)
    implicit none
    integer, dimension(3) :: ext
    real*8, dimension(ext(1),ext(2),ext(3)) :: X,Y,Z,P,Ax,Ay,Az,phi
    real*8 :: betam1,p_c,x_NS_CoM_coord,y_NS_CoM_coord,bhns_P_max,r0
    real*8 :: xp_s,yp_s,zp_s,xp,yp,zp,hdx,hdy,hdz,psi4s,P_s,P_s_Pcut
    integer :: i,j,k,ip1,jp1,kp1
    real*8, parameter :: f8pi = 25.13274122871834590768d0
    hdX = 0.5d0*(X(2,1,1) - X(1,1,1))
    hdY = 0.5d0*(Y(1,2,1) - Y(1,1,1))
    hdZ = 0.5d0*(Z(1,1,2) - Z(1,1,1))
    do k=1,ext(3)
       do j=1,ext(2)
          do i=1,ext(1)
             xp_s = X(i,1,1)+hdX
             yp_s = Y(1,j,1)+hdY
             zp_s = Z(1,1,k)+hdZ
             ip1 = min(i+1,ext(1))
             jp1 = min(j+1,ext(2))
             kp1 = min(k+1,ext(3))
             !x-component of A_i :
             !           Ps = 0.25d0*(P(i,j,k) + P(i,jp1,k) + P(i,j,kp1) + P(i,jp1,kp1) )
             !           Ps_Pcut = max(Ps - p_c*bhns_P_max, 0.d0)
             !           psi4s = exp( phi(i,j,k) + phi(i,jp1,k) + phi(i,j,kp1) + phi(i,jp1,kp1) )
             !           Ax(i,j,k) = xi*zs/r0*sqrt(f8pi*betam1*psi4s*Ps_Pcut)
             xp = X(i,1,1)
             yp = yp_s
             zp = zp_s
             P_s = 0.25d0*( P(i,j,k) + P(i,jp1,k) + P(i,j,kp1) + P(i,jp1,kp1) )
             P_s_Pcut = max(P_s - p_c*bhns_P_max, 0.d0)
             psi4s = exp( phi(i,j,k) + phi(i,jp1,k) + phi(i,j,kp1) + phi(i,jp1,kp1) )
             Ax(i,j,k) = (xp-x_NS_CoM_coord)*zp/r0*sqrt(f8pi*betam1*psi4s*P_s_Pcut)
             !y-component of A_i :
             !           Ps = 0.25d0*(P(i,j,k) + P(ip1,j,k) + P(i,j,kp1) + P(ip1,j,kp1) )
             !           Ps_Pcut = max(Ps - p_c*bhns_P_max, 0.d0)
             !           psi4s = exp( phi(i,j,k) + phi(ip1,j,k) + phi(i,j,kp1) + phi(ip1,j,kp1) )
             !           Ay(i,j,k) = yj*zs/r0*sqrt(f8pi*betam1*psi4s*Ps_Pcut)
             xp = xp_s
             yp = Y(1,j,1)
             zp = zp_s
             P_s = 0.25d0*(P(i,j,k) + P(ip1,j,k) + P(i,j,kp1) + P(ip1,j,kp1) )
             P_s_Pcut = max(P_s - p_c*bhns_P_max, 0.d0)
             psi4s = exp( phi(i,j,k) + phi(ip1,j,k) + phi(i,j,kp1) + phi(ip1,j,kp1) )
             Ay(i,j,k) = (yp-y_NS_CoM_coord)*zp/r0*sqrt(f8pi*betam1*psi4s*P_s_Pcut)
             !z-component of A_i :
             Az(i,j,k) = 0.d0
          end do
       end do
    end do
  end subroutine BHNS_compute_A_toroidal_constrained_transport_eq_3
  subroutine BHNS_compute_Bi(ext,X,Y,Z,PhysicalRadius, &
       phi, &
       A_phi,A_phix,A_phiy,A_phiz, &
       Ax,Ay,Az, &
       Bx,By,Bz, &
       x_NS_CoM_coord,y_NS_CoM_coord,rho_b_atm, &
       constrained_transport_scheme,em_field_type)
    ! INPUT/OUTPUT PARAMETERS
    integer, dimension(3)                    :: ext
    real*8, dimension(ext(1),ext(2),ext(3))  :: X,Y,Z,PhysicalRadius
    real*8, dimension(ext(1),ext(2),ext(3))  :: phi
    real*8, dimension(ext(1),ext(2),ext(3))  :: A_phi,A_phix,A_phiy,A_phiz
    real*8, dimension(ext(1),ext(2),ext(3))  :: Ax,Ay,Az
    real*8, dimension(ext(1),ext(2),ext(3))  :: Bx,By,Bz
    real*8                                   :: x_NS_CoM_coord,y_NS_CoM_coord,rho_b_atm
    integer                                    :: constrained_transport_scheme,em_field_type
    ! INTERNAL PARAMETERS
    integer                                  :: i,j,k
    integer                                  :: imin,imax,jmin,jmax,kmin,kmax
    real*8                                   :: dX,dY,dZ,psim6,psim6_s,pomega2,al,sqrtg,sqrtg4,B2s
    real*8                                   :: fs4pi,B_xs,B_ys,B_zs,sb0,sb2,xn,yn
    real*8                                   :: sb_x,sb_y,sb_z,psi4,u_x,u_y,u_z
    real*8                                   :: psin
    real*8                                   :: Yijk,Yijkp1,Yijp1k,Yijp1kp1,Yip1jk,Yip1jkp1,Yip1jp1k,Yip1jp1kp1
    real*8                                   :: Xijk,Xijkp1,Xijp1k,Xijp1kp1,Xip1jk,Xip1jkp1,Xip1jp1k,Xip1jp1kp1
    real*8                                   :: Ax000,Ax001,Ax010,Ax011,Ax100,Ax101,Ax110,Ax111
    real*8                                   :: Ay000,Ay001,Ay010,Ay011,Ay100,Ay101,Ay110,Ay111
    integer                                  :: AXISYM,EQUATORIAL,OCTANT
    parameter(EQUATORIAL = 1, OCTANT = 2, AXISYM = 4)
    real*8, parameter                        :: SYM = 1.d0, ANTI = -1.d0
    dX = X(2,1,1) - X(1,1,1)
    dY = Y(1,2,1) - Y(1,1,1)
    dZ = Z(1,1,2) - Z(1,1,1)
    imin = 1
    jmin = 1
    kmin = 1
    imax = ext(1)
    jmax = ext(2)
    kmax = ext(3)
    fs4pi = sqrt(4.d0*acos(-1.d0))
    ! Compute the derivatives of A_phi, if we are not using constrained_transport_scheme==3
    !
    if (constrained_transport_scheme .ne. 3) then
       write(*,*) 'SORRY, CONSTRAINED_TRANSPORT_SCHEME != 3 NOT SUPPORTED WHEN SETTING B FIELDS AFTER CHECKPOINT RESTART.'
       write(*,*) 'You will need to compute Bitildes at the end of the ...setBi_after_checkpoint.F90 routine!'
       stop
    end if
    ! Compute B^i on staggered grid and temporarily store them in A_phix,A_phiy,A_phiz
    !
    if (constrained_transport_scheme==3) then
       do k=1,ext(3)
          do j=1,ext(2)
             do i=1,ext(1)
                im1 = max(i-1,1)
                jm1 = max(j-1,1)
                km1 = max(k-1,1)
                ip1 = min(i+1,ext(1))
                jp1 = min(j+1,ext(2))
                kp1 = min(k+1,ext(3))
                psim6_s = exp(-3.d0 * (phi(i,j,k) + phi(ip1,j,k)) )
                A_phix(i,j,k) = ( (Az(i,j,k)-Az(i,jm1,k))/dY   &
                     - (Ay(i,j,k)-Ay(i,j,km1))/dZ ) * psim6_s
                psim6_s = exp(-3.d0 * (phi(i,j,k) + phi(i,jp1,k)) )
                A_phiy(i,j,k) = ( (Ax(i,j,k)-Ax(i,j,km1))/dZ &
                     - (Az(i,j,k)-Az(im1,j,k))/dX ) * psim6_s
                psim6_s = exp(-3.d0 * (phi(i,j,k) + phi(i,j,kp1)) )
                A_phiz(i,j,k) = ( (Ay(i,j,k)-Ay(im1,j,k))/dX &
                     - (Ax(i,j,k)-Ax(i,jm1,k))/dY ) * psim6_s
             end do
          end do
       end do
       ! Now compute B^i on unstaggered grid by simple averge
       do k=1,ext(3)
          do j=1,ext(2)
             do i=1,ext(1)
                im1 = max(i-1,1)
                jm1 = max(j-1,1)
                km1 = max(k-1,1)
                Bx(i,j,k) = 0.5d0* (A_phix(i,j,k) + A_phix(im1,j,k))
                By(i,j,k) = 0.5d0* (A_phiy(i,j,k) + A_phiy(i,jm1,k))
                Bz(i,j,k) = 0.5d0* (A_phiz(i,j,k) + A_phiz(i,j,km1))
             end do
          end do
       end do
    end if
    !
    ! Now compute B^i according to (exercise for the readers)
    !  B^x = (-x/pomega^2) e^(-6 phi) * A_{phi,z}; 
    !  B^y = (-y/pomega^2) e^(-6 phi) * A_{phi,z};
    !  B^z = e^(-6 phi) * (x A_{phi,x} + y A_{phi,y})/pomega^2; 
    !  pomega^2 = x^2 + y^2
    !
    ! and then calculate mhd_st_i and tau
    !
    do k = kmin,kmax
       do j = jmin,jmax
          do i = imin,imax
             ip1 = min(i+1,ext(1))
             jp1 = min(j+1,ext(2))
             kp1 = min(k+1,ext(3))
             xn = X(i,j,k) - x_NS_CoM_coord
             yn = Y(i,j,k) - y_NS_CoM_coord
             psim6 = exp(-6.d0*phi(i,j,k))
             pomega2 = xn**2 + yn**2
             if (constrained_transport_scheme==3) then
                ! do nothing since B^i has been computed 
             elseif (em_field_type==0 .or. em_field_type==3 .or. i==imax .or. j==jmax .or. k==kmax) then
                Bx(i,j,k) = -xn/pomega2 * psim6 * A_phiz(i,j,k)
                By(i,j,k) = -yn/pomega2 * psim6 * A_phiz(i,j,k)
                Bz(i,j,k) = psim6/pomega2 * (xn*A_phix(i,j,k) +  &
                     yn*A_phiy(i,j,k))
             else
                ! Compute Bx, By, Bz that satisfy div(B)=0 to machine precision.
                ! Here is the recipe, which has been verified by Mathematica: 
                !
                ! \tilde{B}^x(i,j,k)=-(Ay(i,j,kp1)-Ay(i,j,k) + Ay(ip1,j,kp1)-Ay(ip1,j,k) + Ay(i,jp1,kp1)-Ay(i,jp1,k) + Ay(ip1,jp1,kp
                ! \tilde{B}^y(i,j,k)=(Ax(i,j,kp1)-Ax(i,j,k) + Ax(ip1,j,kp1)-Ax(ip1,j,k) + Ax(i,jp1,kp1)-Ax(i,jp1,k) + Ax(ip1,jp1,kp1
                ! \tilde{B}^z(i,j,k):=(Ay(ip1,j,k)-Ay(i,j,k) + Ay(ip1,jp1,k)-Ay(i,jp1,k) + Ay(ip1,j,kp1)-Ay(i,j,kp1) + Ay(ip1,jp1,kp
                ! Here \tilde{B}^i = psi^6 B^i, and 
                ! Ax(i,j,k), Ay(i,j,k), and Az(i,j,k) are the 
                ! 3 (covariant) components of the vector potential at the corner 
                ! point (x_i-dx/2, y_j-dy/2, z_k-dz/2).
                ! In the present case, we set Ax = -A_phi * y/(x^2+y^2), 
                !  Ay = A_phi * x/(x^2+y^2), and Az = 0.
                ! Note that B^i cant be set this way at the boundary points
                ! becuase of the sentcil structure. Presumably, this 
                ! is not a problem when we have enough ghostzones.
                !
                ! In the following, Axabc denotes Ax(i+a,j+b,k+c) and 
                ! similarly for Ayabc.
                ! 
                ! Here we set the coordinates, relative to the origin, which in this case is the center of the NS.
                Yijk      = yn
                Yijkp1    = yn
                Yijp1k    = yn + dY
                Yijp1kp1  = yn + dY
                Yip1jk    = yn
                Yip1jkp1  = yn
                Yip1jp1k  = yn + dY
                Yip1jp1kp1= yn + dY
                Xijk      = xn
                Xijkp1    = xn
                Xijp1k    = xn
                Xijp1kp1  = xn
                Xip1jk    = xn + dX
                Xip1jkp1  = xn + dX
                Xip1jp1k  = xn + dX
                Xip1jp1kp1= xn + dX
                Ax000 = -A_phi(i,j,k)* (Yijk-0.5d0*dY)/ &
                     ((Xijk-0.5d0*dX)**2 + (Yijk-0.5d0*dY)**2)
                Ax001 = -A_phi(i,j,kp1)*(Yijkp1-0.5d0*dY)/ &
                     ((Xijkp1-0.5d0*dX)**2 + (Yijkp1-0.5d0*dY)**2)
                Ax010 = -A_phi(i,jp1,k)* (Yijp1k-0.5d0*dY)/ &
                     ((Xijp1k-0.5d0*dX)**2 + (Yijp1k-0.5d0*dY)**2)
                Ax011 = -A_phi(i,jp1,kp1)*(Yijp1kp1-0.5d0*dY)/ &
                     ((Xijp1kp1-0.5d0*dX)**2 + (Yijp1kp1-0.5d0*dY)**2)
                Ax100 = -A_phi(ip1,j,k)*(Yip1jk-0.5d0*dY)/ &
                     ((Xip1jk-0.5d0*dX)**2 + (Yip1jk-0.5d0*dY)**2)
                Ax101 = -A_phi(ip1,j,kp1)*(Yip1jkp1-0.5d0*dY)/ &
                     ((Xip1jkp1-0.5d0*dX)**2 + (Yip1jkp1-0.5d0*dY)**2)
                Ax110 = -A_phi(ip1,jp1,k)*(Yip1jp1k-0.5d0*dY)/ &
                     ((Xip1jp1k-0.5d0*dX)**2 + (Yip1jp1k-0.5d0*dY)**2)
                Ax111 = -A_phi(ip1,jp1,kp1)*(Yip1jp1kp1-0.5d0*dY)/ &
                     ((Xip1jp1kp1-0.5d0*dX)**2 + (Yip1jp1kp1-0.5d0*dY)**2)
                Ay000 = A_phi(i,j,k)* (Xijk-0.5d0*dX)/ &
                     ((Xijk-0.5d0*dX)**2 + (Yijk-0.5d0*dY)**2)
                Ay001 = A_phi(i,j,kp1)*(Xijkp1-0.5d0*dX)/ &
                     ((Xijkp1-0.5d0*dX)**2 + (Yijkp1-0.5d0*dY)**2)
                Ay010 = A_phi(i,jp1,k)* (Xijp1k-0.5d0*dX)/ &
                     ((Xijp1k-0.5d0*dX)**2 + (Yijp1k-0.5d0*dY)**2)
                Ay011 = A_phi(i,jp1,kp1)*(Xijp1kp1-0.5d0*dX)/ &
                     ((Xijp1kp1-0.5d0*dX)**2 + (Yijp1kp1-0.5d0*dY)**2)
                Ay100 = A_phi(ip1,j,k)*(Xip1jk-0.5d0*dX)/ &
                     ((Xip1jk-0.5d0*dX)**2 + (Yip1jk-0.5d0*dY)**2)
                Ay101 = A_phi(ip1,j,kp1)*(Xip1jkp1-0.5d0*dX)/ &
                     ((Xip1jkp1-0.5d0*dX)**2 + (Yip1jkp1-0.5d0*dY)**2)
                Ay110 = A_phi(ip1,jp1,k)*(Xip1jp1k-0.5d0*dX)/ &
                     ((Xip1jp1k-0.5d0*dX)**2 + (Yip1jp1k-0.5d0*dY)**2)
                Ay111 = A_phi(ip1,jp1,kp1)*(Xip1jp1kp1-0.5d0*dX)/ &
                     ((Xip1jp1kp1-0.5d0*dX)**2 + (Yip1jp1kp1-0.5d0*dY)**2)
                Bx(i,j,k) = -( (Ay001-Ay000) + (Ay101-Ay100) + (Ay011-Ay010) + (Ay111-Ay110) ) * 0.25d0/dZ * psim6
                By(i,j,k) = ( (Ax001-Ax000) + (Ax101-Ax100) + (Ax011-Ax010) + (Ax111-Ax110) ) * 0.25d0/dZ * psim6
                Bz(i,j,k) = ( (Ay100-Ay000) + (Ay110-Ay010) + (Ay101-Ay001) + (Ay111-Ay011) ) * 0.25d0/dX * psim6 - ( (Ax010-Ax000)&
  & + (Ax110-Ax100) + (Ax011-Ax001) + (Ax111-Ax101) ) * 0.25d0/dY * psim6
             end if
          end do
       end do
    end do
  end subroutine BHNS_compute_Bi
  subroutine compute_Aphi_bhns_current_loop_app(x,y,z,x_ns_init, y_ns_init,a,I0,Aphi)
    implicit none
    real*8, intent(in)  :: x,y,z,x_ns_init, y_ns_init, a,I0
    real*8, intent(out) :: Aphi
    real*8 :: xn,yn,r,pm2, w2,w
    real*8, parameter :: pi = 3.14159265358979323846d0
    real*8, parameter :: fiftno23 = 0.6521739130434782608695652d0
    xn = x - x_ns_init
    yn = y - y_ns_init
    pm2 = xn*xn + yn*yn
    w2 = a*a+z*z+pm2
    w = sqrt(w2)
    Aphi = fiftno23*pi*a*a*pm2*I0*(1.d0 + 1.875d0*(pm2+a*a)*a*a/(w2*w2))/(w2*w)
  end subroutine compute_Aphi_bhns_current_loop_app
