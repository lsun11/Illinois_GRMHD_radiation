subroutine particle_tracer(cctk_dim,cctk_gsh,cctk_lsh,cctk_lbnd,cctk_ubnd,cctk_lssh,cctk_from,cctk_to,cctk_bbox,cctk_delta_time,cct&
  &k_time,cctk_delta_space,cctk_origin_space,cctk_levfac,cctk_levoff,cctk_levoffdenom,cctk_timefac,cctk_convlevel,cctk_convfac,cctk&
  &_nghostzones,cctk_iteration,cctkGH, X0analytic,X0h17,X0more_pcle_stuff,X0particle_tracer,X0pcle_stuff,X0v_previous,X1analytic,X1&
  &h17,X1more_pcle_stuff,X1particle_tracer,X1v_previous,X2analytic,X2h17,X2v_previous,E_anal,E_rad_gradient,E_rad_pcle,Ec_anal,Ec_p&
  &cle,F_anal,F_coeff,F_pcle,Fc_anal,Fc_pcle,P_pcle,Q,areal_radius,coord,coordt,eta,gxx_pcle,gxz_pcle,gzz_pcle,h17_rhs,phi_pcle,pos&
  &,rho_b_anal,rho_b_analytic,rho_b_pcle,slope,slopem,slopet,tau_center_OS,u0_p,u0_pcle,vx_p,vy_p,vz_p,z_anal, X0BSSN_AH,X0BSSN_aux&
  &_private,X0BSSN_aux_restrict2,X0BSSN_diag_restrict,X0BSSN_gupij,X0BSSN_matter,X0BSSN_refbd,X0BSSN_rhs,X0BSSN_vars,X0NPphi2_em,X0&
  &OmegF,X0Psi4_gw,X0Stagger_Bs,X0ah_radius,X0bh_posns,X0diag_vacuum_gf,X0disk_atmosphere,X0em_Ax,X0em_Ay,X0em_Az,X0em_Blagrangemul&
  &tiplier,X0em_Blagrangemultiplier_rhs,X0em_Phi,X0em_Phi_rhs,X0em_conservativex,X0em_conservativey,X0em_conservativez,X0em_rhsx,X0&
  &em_rhsy,X0em_rhsz,X0eos_params1,X0eos_params2,X0excision_int_gfs,X0field_line_variables,X0field_line_variables_rhs,X0fisheye_var&
  &s,X0gw_moment_arrays,X0have_global_bdry,X0lapse_derivatives,X0lapse_rhs,X0lapse_vars,X0lapse_vars_aux,X0metric_spher_pol_1,X0met&
  &ric_spher_pol_2,X0metric_spher_pol_3,X0mhd_conservatives,X0mhd_nosync,X0mhd_primitives,X0mhd_rhs,X0mhd_sync_lr_B_quantities,X0mh&
  &d_sync_lr_hydro_quantities,X0mhd_sync_metric_facevals,X0mhd_sync_nablas,X0mhd_sync_nablas_diag,X0mhd_sync_rho_br_rho_bl,X0mhd_te&
  &mps,X0mhd_vs,X0micphys_conservatives,X0micphys_conservatives_rhs,X0micphys_sync_lr_hydro_quantities,X0micphys_sync_nablas,X0micr&
  &ophys_primitives,X0phi_derivs,X0rad_conservatives,X0rad_conservatives_rhs,X0rad_pressure,X0rad_primitives,X0rad_sync_lr,X0rhovec&
  &s,X0shift_rhs,X0shift_vars,X0shift_vars_temp,X0surf_params,X1BSSN_AH,X1BSSN_aux_private,X1BSSN_aux_restrict2,X1BSSN_diag_restric&
  &t,X1BSSN_gupij,X1BSSN_matter,X1BSSN_refbd,X1BSSN_rhs,X1BSSN_vars,X1NPphi2_em,X1OmegF,X1Psi4_gw,X1Stagger_Bs,X1ah_radius,X1diag_v&
  &acuum_gf,X1disk_atmosphere,X1em_Ax,X1em_Ay,X1em_Az,X1em_Blagrangemultiplier,X1em_Blagrangemultiplier_rhs,X1em_Phi,X1em_Phi_rhs,X&
  &1em_conservativex,X1em_conservativey,X1em_conservativez,X1em_rhsx,X1em_rhsy,X1em_rhsz,X1excision_int_gfs,X1field_line_variables,&
  &X1field_line_variables_rhs,X1fisheye_vars,X1have_global_bdry,X1lapse_derivatives,X1lapse_rhs,X1lapse_vars,X1lapse_vars_aux,X1mhd&
  &_conservatives,X1mhd_nosync,X1mhd_primitives,X1mhd_rhs,X1mhd_sync_lr_B_quantities,X1mhd_sync_lr_hydro_quantities,X1mhd_sync_metr&
  &ic_facevals,X1mhd_sync_nablas,X1mhd_sync_nablas_diag,X1mhd_sync_rho_br_rho_bl,X1mhd_temps,X1mhd_vs,X1micphys_conservatives,X1mic&
  &phys_conservatives_rhs,X1micphys_sync_lr_hydro_quantities,X1micphys_sync_nablas,X1microphys_primitives,X1phi_derivs,X1rad_conser&
  &vatives,X1rad_conservatives_rhs,X1rad_pressure,X1rad_primitives,X1rad_sync_lr,X1shift_rhs,X1shift_vars,X1shift_vars_temp,X2BSSN_&
  &AH,X2BSSN_aux_private,X2BSSN_aux_restrict2,X2BSSN_diag_restrict,X2BSSN_gupij,X2BSSN_matter,X2BSSN_refbd,X2BSSN_rhs,X2BSSN_vars,X&
  &2NPphi2_em,X2OmegF,X2Psi4_gw,X2Stagger_Bs,X2ah_radius,X2diag_vacuum_gf,X2disk_atmosphere,X2em_Ax,X2em_Ay,X2em_Az,X2em_Blagrangem&
  &ultiplier,X2em_Blagrangemultiplier_rhs,X2em_Phi,X2em_Phi_rhs,X2em_conservativex,X2em_conservativey,X2em_conservativez,X2em_rhsx,&
  &X2em_rhsy,X2em_rhsz,X2excision_int_gfs,X2field_line_variables,X2field_line_variables_rhs,X2fisheye_vars,X2lapse_derivatives,X2la&
  &pse_rhs,X2lapse_vars,X2lapse_vars_aux,X2mhd_conservatives,X2mhd_nosync,X2mhd_primitives,X2mhd_rhs,X2mhd_sync_lr_B_quantities,X2m&
  &hd_sync_lr_hydro_quantities,X2mhd_sync_metric_facevals,X2mhd_sync_nablas,X2mhd_sync_nablas_diag,X2mhd_sync_rho_br_rho_bl,X2mhd_t&
  &emps,X2mhd_vs,X2micphys_conservatives,X2micphys_conservatives_rhs,X2micphys_sync_lr_hydro_quantities,X2micphys_sync_nablas,X2mic&
  &rophys_primitives,X2phi_derivs,X2rad_conservatives,X2rad_conservatives_rhs,X2rad_pressure,X2rad_primitives,X2rad_sync_lr,X2shift&
  &_rhs,X2shift_vars,X2shift_vars_temp,ah_centroid_length,ah_flags_length,ah_origin_length,ah_radius_length,Aupxx,Aupxy,Aupxz,Aupyy&
  &,Aupyz,Aupzz,Ax,Ax_p,Ax_p_p,Axx,Axx_p,Axx_p_p,Axx_rhs,Axy,Axy_p,Axy_p_p,Axy_rhs,Axz,Axz_p,Axz_p_p,Axz_rhs,Ay,Ay_p,Ay_p_p,Ayy,Ayy&
  &_p,Ayy_p_p,Ayy_rhs,Ayz,Ayz_p,Ayz_p_p,Ayz_rhs,Az,Az_p,Az_p_p,Azz,Azz_p,Azz_p_p,Azz_rhs,Bfreq1,Bfreq1_p,Bfreq1_p_p,Bfreq2,Bfreq2_p&
  &,Bfreq2_p_p,Blagrangemultiplier,Blagrangemultiplier_p,Blagrangemultiplier_p_p,Blagrangemultiplier_rhs,Box1X_VolInt,Box1X_VolInt1&
  &,Box1X_VolInt2,Box1Y_VolInt,Box1Y_VolInt1,Box1Y_VolInt2,Box1Z_VolInt,Box1Z_VolInt1,Box1Z_VolInt2,Box1denom_VolInt,Box1denom_VolI&
  &nt1,Box1denom_VolInt2,Bx,Bx_stagger,Bxl,Bxr,Bxtilde,Bxtilde_p,Bxtilde_p_p,Bxtilde_or_Ax_rhs,By,By_stagger,Byl,Byr,Bytilde,Bytild&
  &e_p,Bytilde_p_p,Bytilde_or_Ay_rhs,Bz,Bz_stagger,Bzl,Bzr,Bztilde,Bztilde_p,Bztilde_p_p,Bztilde_or_Az_rhs,CoM_VolInt_denominator,C&
  &oMx_VolInt,CoMy_VolInt,CoMz_VolInt,E_rad,E_rad_nue,E_rad_nuel,E_rad_nuer,E_rad_nux,E_rad_nuxl,E_rad_nuxr,E_radl,E_radr,Ex,Ey,Ez,&
  &F_rad0,F_rad0_nue,F_rad0_nux,F_rad_scalar,F_rad_scalar_nue,F_rad_scalar_nux,F_radx,F_radx_nue,F_radx_nuele,F_radx_nuer,F_radx_nu&
  &x,F_radx_nuxle,F_radx_nuxr,F_radxle,F_radxr,F_rady,F_rady_nue,F_rady_nuele,F_rady_nuer,F_rady_nux,F_rady_nuxle,F_rady_nuxr,F_rad&
  &yle,F_radyr,F_radz,F_radz_nue,F_radz_nuele,F_radz_nuer,F_radz_nux,F_radz_nuxle,F_radz_nuxr,F_radzle,F_radzr,FaFal,FaFal_nue,FaFa&
  &l_nux,FaFar,FaFar_nue,FaFar_nux,Fontfix_tracker_gf,Gammax,Gammax_p,Gammax_p_p,Gammax_drive,Gammax_rhs,Gammaxxx,Gammaxxy,Gammaxxz&
  &,Gammaxyy,Gammaxyz,Gammaxzz,Gammay,Gammay_p,Gammay_p_p,Gammay_drive,Gammay_rhs,Gammayxx,Gammayxy,Gammayxz,Gammayyy,Gammayyz,Gamm&
  &ayzz,Gammaz,Gammaz_p,Gammaz_p_p,Gammaz_drive,Gammaz_rhs,Gammazxx,Gammazxy,Gammazxz,Gammazyy,Gammazyz,Gammazzz,Gamx_const_VolInt,&
  &Gamy_const_VolInt,Gamz_const_VolInt,Ham_const_VolIntD,Ham_const_VolIntN,Ham_const_excised_VolIntD,Ham_const_excised_VolIntN,Ham_&
  &const_excised_innerregion_VolIntD,Ham_const_excised_innerregion_VolIntN,J_ADM_VolInt,J_ADM_VolInt_inner,J_ADM_VolInt_inner2,Jz_c&
  &onstraint,Jz_constraint_inner,KK,K_pp_rot1,K_pp_rot2,K_pp_rot3,K_rp_rot1,K_rp_rot2,K_rp_rot3,K_rr_rot1,K_rr_rot2,K_rr_rot3,K_rth&
  &_rot1,K_rth_rot2,K_rth_rot3,K_thp_rot1,K_thp_rot2,K_thp_rot3,K_thth_rot1,K_thth_rot2,K_thth_rot3,M0_AH_VolInt,M0_VolInt,M0_escap&
  &e100M,M0_escape30M,M0_escape50M,M0_escape70M,M0_horiz_VolInt,M0_r1_VolInt,M0_r2_VolInt,M0_r3_VolInt,MNorm,MONOPOLE,MRsx,MRsy,MRs&
  &z,M_ADM_VolInt,M_constraint,M_constraint_inner,Minternal_VolInt,Minternal_cold_VolInt,NPphi0im,NPphi0im_p,NPphi0im_p_p,NPphi0re,&
  &NPphi0re_p,NPphi0re_p_p,NPphi2im,NPphi2im_p,NPphi2im_p_p,NPphi2re,NPphi2re_p,NPphi2re_p_p,OS_surf_rad,P,P_cld,P_constraint,P_con&
  &straint_inner,P_floor,P_radxx,P_radxx_nue,P_radxx_nux,P_radxy,P_radxy_nue,P_radxy_nux,P_radxz,P_radxz_nue,P_radxz_nux,P_radyy,P_&
  &radyy_nue,P_radyy_nux,P_radyz,P_radyz_nue,P_radyz_nux,P_radzz,P_radzz_nue,P_radzz_nux,P_tab,P_thermal,PhysicalRadius,Pl,Pr,PsiNo&
  &rm,PsiRes,PsiTau,Pvec,RadiusDerivative,RadiusDerivative2,Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,S,SPx,SPx_p,SPx_p_p,SPy,SPy_p,SPy_p_p,SPz,SPz_p&
  &,SPz_p_p,S_rad_x,S_rad_x_p,S_rad_x_p_p,S_rad_x_nue,S_rad_x_nue_p,S_rad_x_nue_p_p,S_rad_x_nue_rhs,S_rad_x_nux,S_rad_x_nux_p,S_rad&
  &_x_nux_p_p,S_rad_x_nux_rhs,S_rad_x_rhs,S_rad_y,S_rad_y_p,S_rad_y_p_p,S_rad_y_nue,S_rad_y_nue_p,S_rad_y_nue_p_p,S_rad_y_nue_rhs,S&
  &_rad_y_nux,S_rad_y_nux_p,S_rad_y_nux_p_p,S_rad_y_nux_rhs,S_rad_y_rhs,S_rad_z,S_rad_z_p,S_rad_z_p_p,S_rad_z_nue,S_rad_z_nue_p,S_r&
  &ad_z_nue_p_p,S_rad_z_nue_rhs,S_rad_z_nux,S_rad_z_nux_p,S_rad_z_nux_p_p,S_rad_z_nux_rhs,S_rad_z_rhs,Sx,Sxx,Sxy,Sxz,Sy,Syy,Syz,Sz,&
  &Szz,T_0x,T_0x_p,T_0x_p_p,T_0y,T_0y_p,T_0y_p_p,T_0z,T_0z_p,T_0z_p_p,T_VolInt,T_fluid,T_fluidle,T_fluidr,Tem0_0_VolInt,Tem0_0_outs&
  &ideBH_VolInt,Tfluid0_0_VolInt,Tfluid0_0_outsideBH_VolInt,VolIntegrand,VolIntegrand2,VolIntegrand3,VolIntegrand4,Y_e,Y_el,Y_er,ah&
  &_centroid_iteration,ah_centroid_iteration_p,ah_centroid_t,ah_centroid_t_p,ah_centroid_valid,ah_centroid_valid_p,ah_centroid_x,ah&
  &_centroid_x_p,ah_centroid_y,ah_centroid_y_p,ah_centroid_z,ah_centroid_z_p,ah_found_flag,ah_initial_find_flag,ah_origin_x,ah_orig&
  &in_y,ah_origin_z,ah_radius,ah_really_initial_find_flag,ah_search_flag,b_phi_VolInt,bh_posn_x,bh_posn_y,bh_posn_z,bh_radius_z,bre&
  &m_qee_VolInt,brem_qei_VolInt,chi,chi_p,chi_p_p,chi_rad,chi_rad_nue,chi_rad_nux,chi_rhs,cmax,cmax_rad,cmax_rad_nue,cmax_rad_nux,c&
  &min,cmin_rad,cmin_rad_nue,cmin_rad_nux,dBx_m,dBy_m,dBz_m,dP_m,dT_fluid_m,ddrbrvec,density_modes_i1,density_modes_i2,density_mode&
  &s_i3,density_modes_i4,density_modes_i5,density_modes_i6,density_modes_r0,density_modes_r1,density_modes_r2,density_modes_r3,dens&
  &ity_modes_r4,density_modes_r5,density_modes_r6,div_beta,drbrvec,drhoYe_m,drhoYe_m_x,drhoYe_m_xp1,drho_b_m,drho_b_m_x,drho_b_m_xp&
  &1,dvx_m,dvx_m_x,dvx_m_xp1,dvy_m,dvz_m,em_J_between_VolInt,em_energy2_VolInt,em_energy_VolInt,em_energy_between_VolInt,em_energy_&
  &outsideBH_VolInt,em_energy_outsideradius1_VolInt,em_energy_outsideradius2_VolInt,emission_gf,emission_gf_nue,emission_gf_nux,eps&
  &_cld,eps_tab,eps_thermal,eps_tot,eta_nue,excision_zone_gf,fluid_J_horiz_VolInt,fluid_J_r1_VolInt,fluid_J_r2_VolInt,fluid_J_r3_Vo&
  &lInt,fluid_energy_horiz_VolInt,fluid_energy_r1_VolInt,fluid_energy_r2_VolInt,fluid_energy_r3_VolInt,gamma_tab,gconx,gcony,gconz,&
  &gupxx,gupxx_f,gupxy,gupxy_f,gupxz,gupxz_f,gupyy,gupyy_f,gupyz,gupyz_f,gupzz,gupzz_f,gxx,gxx_p,gxx_p_p,gxx_f,gxx_rhs,gxxx,gxxy,gx&
  &xz,gxy,gxy_p,gxy_p_p,gxy_f,gxy_rhs,gxyx,gxyy,gxyz,gxz,gxz_p,gxz_p_p,gxz_f,gxz_rhs,gxzx,gxzy,gxzz,gyy,gyy_p,gyy_p_p,gyy_f,gyy_rhs&
  &,gyyx,gyyy,gyyz,gyz,gyz_p,gyz_p_p,gyz_f,gyz_rhs,gyzx,gyzy,gyzz,gzz,gzz_p,gzz_p_p,gzz_f,gzz_rhs,gzzx,gzzy,gzzz,h,h_p,half_b2_u0_V&
  &olInt,half_b2_u0_outsideBH_VolInt,have_global_bdry_max,have_global_bdry_min,int_momentsi,int_momentsr,k_tab,ka_gf,ka_gf_nue,ka_g&
  &f_nux,ks_gf,ks_gf_nue,ks_gf_nux,kxx,kxy,kxz,kyy,kyz,kzz,lapm1,lapm1_p,lapm1_p_p,lapm1_f,lapm1_rhs,lapm1_rot1,lapm1_rot2,lapm1_ro&
  &t3,lapset,lapset_p,lapset_p_p,lapset_rhs,lapset_timederiv,lapsex,lapsey,lapsez,mass_star1,mass_star2,mhd_chi_line,mhd_chi_line_p&
  &,mhd_chi_line_p_p,mhd_chi_line_rhs,mhd_psi_line,mhd_psi_line_p,mhd_psi_line_p_p,mhd_psi_line_rhs,mhd_st_x,mhd_st_x_p,mhd_st_x_p_&
  &p,mhd_st_x_rhs,mhd_st_y,mhd_st_y_p,mhd_st_y_p_p,mhd_st_y_rhs,mhd_st_z,mhd_st_z_p,mhd_st_z_p_p,mhd_st_z_rhs,mhd_u_chi,mhd_u_chi_p&
  &,mhd_u_chi_p_p,mhd_u_chi_rhs,mhd_u_psi,mhd_u_psi_p,mhd_u_psi_p_p,mhd_u_psi_rhs,minternal_cold_horiz_VolInt,minternal_cold_r1_Vol&
  &Int,minternal_cold_r2_VolInt,minternal_cold_r3_VolInt,minternal_horiz_VolInt,minternal_r1_VolInt,minternal_r2_VolInt,minternal_r&
  &3_VolInt,mom_const_VolIntD,mom_const_excised_VolIntD,momentsi,momentsi_old,momentsr,momentsr_old,momx_const_VolIntN,momx_const_e&
  &xcised_VolIntN,momy_const_VolIntN,momy_const_excised_VolIntN,momz_const_VolIntN,momz_const_excised_VolIntN,monopole_VolInt,monop&
  &ole_outsideBH_VolInt,mu_nu,n_poly,odd_momentsi,odd_momentsi_old,odd_momentsr,odd_momentsr_old,optd,optd_x,optd_y,optd_z,pfloor_g&
  &f,phi,phi_p,phi_p_p,phi_f,phi_rhs,phi_rot1,phi_rot2,phi_rot3,phix,phiy,phiz,psi,psi0im,psi0im_p,psi0im_p_p,psi0re,psi0re_p,psi0r&
  &e_p_p,psi6phi,psi6phi_p,psi6phi_p_p,psi6phi_rhs,rad_const,rad_energy_VolInt,rad_energy_nue_VolInt,rad_energy_nux_VolInt,rbrvec,r&
  &efbd,refbd_p,refbd_p_p,rho,rhoYe,rhoYe_p,rhoYe_p_p,rhoYe_rhs,rho_b,rho_b_atm_gf,rho_bl,rho_br,rho_star,rho_star_p,rho_star_p_p,r&
  &ho_star_rhs,rho_tab,rhob_floor,rhovec,sbt,sbx,sby,sbz,shiftp_rot1,shiftp_rot2,shiftp_rot3,shiftr_rot1,shiftr_rot2,shiftr_rot3,sh&
  &iftth_rot1,shiftth_rot2,shiftth_rot3,shiftx,shiftx_p,shiftx_p_p,shiftx_f,shiftx_rhs,shiftxt,shiftxt_p,shiftxt_p_p,shiftxt_rhs,sh&
  &iftxt_timederiv,shifty,shifty_p,shifty_p_p,shifty_f,shifty_rhs,shiftyt,shiftyt_p,shiftyt_p_p,shiftyt_rhs,shiftyt_timederiv,shift&
  &z,shiftz_p,shiftz_p_p,shiftz_f,shiftz_rhs,shiftzt,shiftzt_p,shiftzt_p_p,shiftzt_rhs,shiftzt_timederiv,smallb2,st_x,st_y,st_z,sur&
  &fvec,synch_VolInt,tau,tau_p,tau_p_p,tau_rad,tau_rad_p,tau_rad_p_p,tau_rad_nue,tau_rad_nue_p,tau_rad_nue_p_p,tau_rad_nue_rhs,tau_&
  &rad_nux,tau_rad_nux_p,tau_rad_nux_p_p,tau_rad_nux_rhs,tau_rad_rhs,tau_rhs,temp1,temp10,temp11,temp12,temp13,temp14,temp15,temp16&
  &,temp17,temp18,temp19,temp2,temp20,temp21,temp22,temp23,temp24,temp25,temp26,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp_g00,&
  &trK,trK_p,trK_p_p,trK_rhs,trR,trRtilde,u0,v02_rad_nuel,v02_rad_nuer,v02_rad_nuxl,v02_rad_nuxr,v02_radl,v02_radr,v02l,v02r,vvec,v&
  &x,vxl,vxr,vy,vyl,vyr,vz,vzl,vzr,w,zeta_rad,zeta_rad_nue,zeta_rad_nux, X0coordinates,X1coordinates,X2coordinates,coarse_dx,coarse&
  &_dy,coarse_dz,r,x,y,z)
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
  INTEGER X0analytic
  INTEGER X0h17
  INTEGER X0more_pcle_stuff
  INTEGER X0particle_tracer
  INTEGER X0pcle_stuff
  INTEGER X0v_previous
  INTEGER X1analytic
  INTEGER X1h17
  INTEGER X1more_pcle_stuff
  INTEGER X1particle_tracer
  INTEGER X1v_previous
  INTEGER X2analytic
  INTEGER X2h17
  INTEGER X2v_previous
  REAL*8 E_anal(X0pcle_stuff)
  REAL*8 E_rad_gradient(X0pcle_stuff)
  REAL*8 E_rad_pcle(X0pcle_stuff)
  REAL*8 Ec_anal(X0pcle_stuff)
  REAL*8 Ec_pcle(X0pcle_stuff)
  REAL*8 F_anal(X0pcle_stuff)
  REAL*8 F_coeff(X0pcle_stuff)
  REAL*8 F_pcle(X0pcle_stuff)
  REAL*8 Fc_anal(X0pcle_stuff)
  REAL*8 Fc_pcle(X0pcle_stuff)
  REAL*8 P_pcle(X0pcle_stuff)
  REAL*8 Q(X0pcle_stuff)
  REAL*8 areal_radius(X0pcle_stuff)
  REAL*8 coord(X0particle_tracer,X1particle_tracer)
  REAL*8 coordt(X0particle_tracer,X1particle_tracer)
  REAL*8 eta(X0pcle_stuff)
  REAL*8 gxx_pcle(X0pcle_stuff)
  REAL*8 gxz_pcle(X0pcle_stuff)
  REAL*8 gzz_pcle(X0pcle_stuff)
  REAL*8 h17_rhs(X0h17,X1h17,X2h17)
  REAL*8 phi_pcle(X0pcle_stuff)
  REAL*8 pos(X0more_pcle_stuff,X1more_pcle_stuff)
  REAL*8 rho_b_anal(X0pcle_stuff)
  REAL*8 rho_b_analytic(X0analytic,X1analytic,X2analytic)
  REAL*8 rho_b_pcle(X0pcle_stuff)
  REAL*8 slope(X0particle_tracer,X1particle_tracer)
  REAL*8 slopem(X0particle_tracer,X1particle_tracer)
  REAL*8 slopet(X0particle_tracer,X1particle_tracer)
  REAL*8 tau_center_OS
  REAL*8 u0_p(X0v_previous,X1v_previous,X2v_previous)
  REAL*8 u0_pcle(X0pcle_stuff)
  REAL*8 vx_p(X0v_previous,X1v_previous,X2v_previous)
  REAL*8 vy_p(X0v_previous,X1v_previous,X2v_previous)
  REAL*8 vz_p(X0v_previous,X1v_previous,X2v_previous)
  REAL*8 z_anal(X0pcle_stuff)
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
  INTEGER X1coordinates
  INTEGER X2coordinates
  REAL*8 coarse_dx
  REAL*8 coarse_dy
  REAL*8 coarse_dz
  REAL*8 r(X0coordinates,X1coordinates,X2coordinates)
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
  REAL*8  R_edge
  REAL*8  betam1
  REAL*8  moncrief_radius_GW(11)
  integer*8  OS_filename
  INTEGER*4 em_field_type
  INTEGER*4 moncrief_gw_num_radii
  INTEGER*4 narr
  INTEGER*4 set_lapse_one
  COMMON /OS_collapsepriv/R_edge,betam1,moncrief_radius_GW,OS_filename,em_field_type,moncrief_gw_num_radii,narr,set_lapse_one
  REAL*8  CCTKH0
  REAL*8  CCTKH1
  REAL*8  CCTKH7
  REAL*8  CCTKH8
  REAL*8  CCTKH9
  REAL*8  CCTKH11
  REAL*8  CCTKH12
  REAL*8  CCTKH13
  REAL*8  CCTKH14
  REAL*8  CCTKH15
  REAL*8  CCTKH16
  REAL*8  rhosurf_rmax
  REAL*8  rhosurf_rmin
  REAL*8  CCTKH17
  INTEGER*4 N_rad_ray
  INTEGER*4 CCTKH2
  INTEGER*4 CCTKH3
  INTEGER*4 CCTKH4
  INTEGER*4 CCTKH5
  INTEGER*4 CCTKH6
  INTEGER*4 CCTKH10
  COMMON /DIAGNOSTICS_MHDrest/CCTKH0,CCTKH1,CCTKH7,CCTKH8,CCTKH9,CCTKH11,CCTKH12,CCTKH13,CCTKH14,CCTKH15,CCTKH16,rhosurf_rmax,rhosu&
  &rf_rmin,CCTKH17,N_rad_ray,CCTKH2,CCTKH3,CCTKH4,CCTKH5,CCTKH6,CCTKH10
  REAL*8  CCTKH18
  REAL*8  CCTKH21
  REAL*8  ddrbddr
  REAL*8  CCTKH22
  REAL*8  drbdr
  REAL*8  CCTKH26
  REAL*8  rbr
  REAL*8  rsurf2
  REAL*8  surf_radius
  INTEGER*4 CCTKH19
  INTEGER*4 N_phi
  INTEGER*4 N_theta
  INTEGER*4 CCTKH20
  INTEGER*4 arithsurf
  INTEGER*4 CCTKH23
  INTEGER*4 CCTKH24
  INTEGER*4 CCTKH25
  INTEGER*4 nsurf
  INTEGER*4 CCTKH27
  INTEGER*4 CCTKH28
  INTEGER*4 CCTKH29
  INTEGER*4 CCTKH30
  INTEGER*4 scaledsurf
  INTEGER*4 sym_factor
  COMMON /DIAGNOSTICS_VACUUMrest/CCTKH18,CCTKH21,ddrbddr,CCTKH22,drbdr,CCTKH26,rbr,rsurf2,surf_radius,CCTKH19,N_phi,N_theta,CCTKH20&
  &,arithsurf,CCTKH23,CCTKH24,CCTKH25,nsurf,CCTKH27,CCTKH28,CCTKH29,CCTKH30,scaledsurf,sym_factor
  INTEGER*4 CCTKH31
  INTEGER*4 CCTKH32
  INTEGER*4 CCTKH33
  INTEGER*4 CCTKH34
  COMMON /DRIVERrest/CCTKH31,CCTKH32,CCTKH33,CCTKH34
  REAL*8  CCTKH35
  REAL*8  CCTKH36
  REAL*8  phi_GW
  REAL*8  CCTKH37
  REAL*8  CCTKH38
  REAL*8  theta_GW
  COMMON /GW_EXTRACTIONrest/CCTKH35,CCTKH36,phi_GW,CCTKH37,CCTKH38,theta_GW
  REAL*8  CCTKH44
  REAL*8  CCTKH58
  REAL*8  CCTKH66
  REAL*8  CCTKH68
  REAL*8  CCTKH70
  REAL*8  CCTKH72
  REAL*8  CCTKH74
  REAL*8  CCTKH76
  REAL*8  CCTKH78
  REAL*8  CCTKH80
  REAL*8  CCTKH82
  integer*8  CCTKH41
  integer*8  CCTKH42
  integer*8  CCTKH45
  integer*8  CCTKH48
  integer*8  CCTKH49
  integer*8  CCTKH50
  integer*8  CCTKH53
  integer*8  CCTKH54
  integer*8  CCTKH59
  integer*8  CCTKH60
  integer*8  CCTKH62
  integer*8  CCTKH84
  integer*8  CCTKH86
  integer*8  CCTKH88
  integer*8  CCTKH90
  integer*8  CCTKH91
  integer*8  CCTKH94
  INTEGER*4 CCTKH39
  INTEGER*4 CCTKH40
  INTEGER*4 CCTKH43
  INTEGER*4 CCTKH46
  INTEGER*4 CCTKH47
  INTEGER*4 CCTKH51
  INTEGER*4 CCTKH52
  INTEGER*4 CCTKH55
  INTEGER*4 CCTKH56
  INTEGER*4 CCTKH57
  INTEGER*4 out_every
  INTEGER*4 CCTKH61
  INTEGER*4 CCTKH63
  INTEGER*4 CCTKH64
  INTEGER*4 CCTKH65
  INTEGER*4 CCTKH67
  INTEGER*4 CCTKH69
  INTEGER*4 CCTKH71
  INTEGER*4 CCTKH73
  INTEGER*4 CCTKH75
  INTEGER*4 CCTKH77
  INTEGER*4 CCTKH79
  INTEGER*4 CCTKH81
  INTEGER*4 CCTKH83
  INTEGER*4 CCTKH85
  INTEGER*4 CCTKH87
  INTEGER*4 CCTKH89
  INTEGER*4 CCTKH92
  INTEGER*4 CCTKH93
  COMMON /IOrest/CCTKH44,CCTKH58,CCTKH66,CCTKH68,CCTKH70,CCTKH72,CCTKH74,CCTKH76,CCTKH78,CCTKH80,CCTKH82,CCTKH41,CCTKH42,CCTKH45,CC&
  &TKH48,CCTKH49,CCTKH50,CCTKH53,CCTKH54,CCTKH59,CCTKH60,CCTKH62,CCTKH84,CCTKH86,CCTKH88,CCTKH90,CCTKH91,CCTKH94,CCTKH39,CCTKH40,CC&
  &TKH43,CCTKH46,CCTKH47,CCTKH51,CCTKH52,CCTKH55,CCTKH56,CCTKH57,out_every,CCTKH61,CCTKH63,CCTKH64,CCTKH65,CCTKH67,CCTKH69,CCTKH71,&
  &CCTKH73,CCTKH75,CCTKH77,CCTKH79,CCTKH81,CCTKH83,CCTKH85,CCTKH87,CCTKH89,CCTKH92,CCTKH93
  REAL*8  CCTKH97
  REAL*8  CCTKH98
  REAL*8  CCTKH99
  integer*8  slicing_type
  INTEGER*4 CCTKH95
  INTEGER*4 CCTKH96
  INTEGER*4 CCTKH100
  INTEGER*4 CCTKH101
  INTEGER*4 opl_a0lap
  INTEGER*4 opl_advect_enable
  INTEGER*4 opl_alap
  COMMON /LAPSErest/CCTKH97,CCTKH98,CCTKH99,slicing_type,CCTKH95,CCTKH96,CCTKH100,CCTKH101,opl_a0lap,opl_advect_enable,opl_alap
  REAL*8  CCTKH104
  REAL*8  CCTKH105
  REAL*8  K_poly
  REAL*8  CCTKH108
  REAL*8  CCTKH110
  REAL*8  P_fact
  REAL*8  CCTKH111
  REAL*8  CCTKH112
  REAL*8  Sym_Bz
  REAL*8  CCTKH114
  REAL*8  CCTKH118
  REAL*8  CCTKH119
  REAL*8  CCTKH123
  REAL*8  CCTKH132
  REAL*8  CCTKH133
  REAL*8  CCTKH136
  REAL*8  gamma_th
  REAL*8  CCTKH140
  REAL*8  CCTKH142
  REAL*8  CCTKH143
  REAL*8  CCTKH144
  REAL*8  CCTKH145
  REAL*8  CCTKH147
  REAL*8  pfloor
  REAL*8  CCTKH151
  REAL*8  CCTKH152
  REAL*8  CCTKH153
  REAL*8  CCTKH154
  REAL*8  CCTKH155
  REAL*8  CCTKH157
  REAL*8  CCTKH158
  REAL*8  CCTKH162
  REAL*8  CCTKH163
  REAL*8  rho_b_atm
  REAL*8  CCTKH166
  REAL*8  rho_fact
  REAL*8  CCTKH167
  REAL*8  CCTKH168
  REAL*8  CCTKH169
  REAL*8  tau_atm
  REAL*8  tau_fact
  INTEGER*4 CCTKH102
  INTEGER*4 CCTKH103
  INTEGER*4 CCTKH106
  INTEGER*4 CCTKH107
  INTEGER*4 CCTKH109
  INTEGER*4 CCTKH113
  INTEGER*4 CCTKH115
  INTEGER*4 CCTKH116
  INTEGER*4 CCTKH117
  INTEGER*4 CCTKH120
  INTEGER*4 CCTKH121
  INTEGER*4 CCTKH122
  INTEGER*4 CCTKH124
  INTEGER*4 CCTKH125
  INTEGER*4 CCTKH126
  INTEGER*4 CCTKH127
  INTEGER*4 CCTKH128
  INTEGER*4 CCTKH129
  INTEGER*4 CCTKH130
  INTEGER*4 CCTKH131
  INTEGER*4 eps_flag
  INTEGER*4 CCTKH134
  INTEGER*4 CCTKH135
  INTEGER*4 CCTKH137
  INTEGER*4 CCTKH138
  INTEGER*4 CCTKH139
  INTEGER*4 CCTKH141
  INTEGER*4 CCTKH146
  INTEGER*4 neos
  INTEGER*4 CCTKH148
  INTEGER*4 CCTKH149
  INTEGER*4 CCTKH150
  INTEGER*4 CCTKH156
  INTEGER*4 CCTKH159
  INTEGER*4 CCTKH160
  INTEGER*4 CCTKH161
  INTEGER*4 CCTKH164
  INTEGER*4 CCTKH165
  INTEGER*4 CCTKH170
  INTEGER*4 CCTKH171
  INTEGER*4 CCTKH172
  INTEGER*4 CCTKH173
  COMMON /MHD_EVOLVErest/CCTKH104,CCTKH105,K_poly,CCTKH108,CCTKH110,P_fact,CCTKH111,CCTKH112,Sym_Bz,CCTKH114,CCTKH118,CCTKH119,CCTK&
  &H123,CCTKH132,CCTKH133,CCTKH136,gamma_th,CCTKH140,CCTKH142,CCTKH143,CCTKH144,CCTKH145,CCTKH147,pfloor,CCTKH151,CCTKH152,CCTKH153&
  &,CCTKH154,CCTKH155,CCTKH157,CCTKH158,CCTKH162,CCTKH163,rho_b_atm,CCTKH166,rho_fact,CCTKH167,CCTKH168,CCTKH169,tau_atm,tau_fact,C&
  &CTKH102,CCTKH103,CCTKH106,CCTKH107,CCTKH109,CCTKH113,CCTKH115,CCTKH116,CCTKH117,CCTKH120,CCTKH121,CCTKH122,CCTKH124,CCTKH125,CCT&
  &KH126,CCTKH127,CCTKH128,CCTKH129,CCTKH130,CCTKH131,eps_flag,CCTKH134,CCTKH135,CCTKH137,CCTKH138,CCTKH139,CCTKH141,CCTKH146,neos,&
  &CCTKH148,CCTKH149,CCTKH150,CCTKH156,CCTKH159,CCTKH160,CCTKH161,CCTKH164,CCTKH165,CCTKH170,CCTKH171,CCTKH172,CCTKH173
  INTEGER*4 Spatial_Gauge
  INTEGER*4 CCTKH174
  INTEGER*4 CCTKH175
  COMMON /SHIFTrest/Spatial_Gauge,CCTKH174,CCTKH175
  
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
  INTEGER*4 function HorizonLocalCoordinateOrigin (horizon_number, origin_x, origin_y, origin_z)
  implicit none
  INTEGER*4 horizon_number
  REAL*8 origin_x
  REAL*8 origin_y
  REAL*8 origin_z
  end function HorizonLocalCoordinateOrigin
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
  INTEGER*4 function HorizonWasFound (horizon_number)
  implicit none
  INTEGER*4 horizon_number
  end function HorizonWasFound
  end interface
  
  integer, dimension(3)                       :: ext,global_ext
  real*8                                      :: dT,dTo2,dTo6,dX,dY,dZ,xmax
  integer                                     :: ii,jj
  real*8, allocatable,dimension(:)            :: u0_p_int,vx_p_int,vy_p_int,vz_p_int
  real*8, allocatable,dimension(:)            :: u0_int,vx_int,vy_int,vz_int
  integer                                     :: vindex
  ext = cctk_lsh
  global_ext = cctk_gsh
  dT = (cctk_delta_time/cctk_timefac)
  dX = (cctk_delta_space(1)/cctk_levfac(1))
  dY = (cctk_delta_space(2)/cctk_levfac(2))
  dZ = (cctk_delta_space(3)/cctk_levfac(3))
  allocate(u0_p_int(narr))
  allocate(vx_p_int(narr))
  allocate(vy_p_int(narr))
  allocate(vz_p_int(narr))
  allocate(u0_int(narr))
  allocate(vx_int(narr))
  allocate(vy_int(narr))
  allocate(vz_int(narr))
  if (CCTK_TIME .ne. 0) then
     call CCTK_VarIndex(vindex,"OS_collapse::u0_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,u0_p_int)
     call CCTK_VarIndex(vindex,"OS_collapse::vx_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vx_p_int)
     call CCTK_VarIndex(vindex,"OS_collapse::vy_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vy_p_int)
     call CCTK_VarIndex(vindex,"OS_collapse::vz_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vz_p_int)
     slope(:,1)=1.d0/u0_p_int
     slope(:,2)=vx_p_int
     slope(:,3)=vy_p_int
     slope(:,4)=vz_p_int
     ! call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !      ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !      Xglobal,Yglobal,Zglobal,pos,1.d0/u0_p,slope(:,1))
     ! call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !      ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !      Xglobal,Yglobal,Zglobal,pos,vx_p,slope(:,2))
     ! call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !      ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !      Xglobal,Yglobal,Zglobal,pos,vy_p,slope(:,3)) 
     ! call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !      ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !      Xglobal,Yglobal,Zglobal,pos,vz_p,slope(:,4)) 
     dTo2=dT*0.5d0
     dTo6=dT/6.0d0
     coordt=coord+dTo2*slope
     write(*,*) "dTo2: ",dTo2
     do ii=1,3
        do jj=1,narr
           pos(jj,ii) = coordt(jj,ii+1)
        end do
     end do
     call CCTK_VarIndex(vindex,"OS_collapse::u0_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,u0_p_int)
     call CCTK_VarIndex(vindex,"OS_collapse::vx_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vx_p_int)
     call CCTK_VarIndex(vindex,"OS_collapse::vy_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vy_p_int)
     call CCTK_VarIndex(vindex,"OS_collapse::vz_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vz_p_int)
     call CCTK_VarIndex(vindex,"mhd_evolve::u0")
     call interp_driver_carp(cctkGH,narr,pos,vindex,u0_int)
     call CCTK_VarIndex(vindex,"mhd_evolve::vx")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vx_int)
     call CCTK_VarIndex(vindex,"mhd_evolve::vy")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vy_int)
     call CCTK_VarIndex(vindex,"mhd_evolve::vz")
     slopet(:,1)=0.5d0*(1.d0/u0_p_int+1.d0/u0_int)
     slopet(:,2)=0.5d0*(vx_p_int+vx_int)
     slopet(:,3)=0.5d0*(vy_p_int+vy_int)
     slopet(:,4)=0.5d0*(vz_p_int+vz_int)
     ! call interp_driver_carp(cctkGH,narr,pos,vindex,vz_p_int)    
     ! call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !      ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !      Xglobal,Yglobal,Zglobal,pos,0.5d0*(1.d0/u0_p+1.d0/u0),slopet(:,1))
     ! call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !      ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !      Xglobal,Yglobal,Zglobal,pos,0.5d0*(vx_p+vx),slopet(:,2))
     ! call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !      ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !      Xglobal,Yglobal,Zglobal,pos,0.5d0*(vy_p+vy),slopet(:,3)) 
     ! call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !      ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !      Xglobal,Yglobal,Zglobal,pos,0.5d0*(vz_p+vz),slopet(:,4)) 
     coordt=coord+dTo2*slopet
     do ii=1,3
        do jj=1,narr
           pos(jj,ii) = coordt(jj,ii+1)
        end do
     end do
     call CCTK_VarIndex(vindex,"OS_collapse::u0_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,u0_p_int)
     call CCTK_VarIndex(vindex,"OS_collapse::vx_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vx_p_int)
     call CCTK_VarIndex(vindex,"OS_collapse::vy_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vy_p_int)
     call CCTK_VarIndex(vindex,"OS_collapse::vz_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vz_p_int)
     call CCTK_VarIndex(vindex,"mhd_evolve::u0")
     call interp_driver_carp(cctkGH,narr,pos,vindex,u0_int)
     call CCTK_VarIndex(vindex,"mhd_evolve::vx")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vx_int)
     call CCTK_VarIndex(vindex,"mhd_evolve::vy")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vy_int)
     call CCTK_VarIndex(vindex,"mhd_evolve::vz")
     slopem(:,1)=0.5d0*(1.d0/u0_p_int+1.d0/u0_int)
     slopem(:,2)=0.5d0*(vx_p_int+vx_int)
     slopem(:,3)=0.5d0*(vy_p_int+vy_int)
     slopem(:,4)=0.5d0*(vz_p_int+vz_int)
     !call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !     ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !     Xglobal,Yglobal,Zglobal,pos,0.5d0*(1.d0/u0_p+1.d0/u0),slopem(:,1))
     !call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !     ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !     Xglobal,Yglobal,Zglobal,pos,0.5d0*(vx_p+vx),slopem(:,2))
     !call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !     ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !     Xglobal,Yglobal,Zglobal,pos,0.5d0*(vy_p+vy),slopem(:,3)) 
     !call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !     ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !     Xglobal,Yglobal,Zglobal,pos,0.5d0*(vz_p+vz),slopem(:,4)) 
     coordt=coord+dT*slopem
     slopem=slopet+slopem
     do ii=1,3
        do jj=1,narr
           pos(jj,ii) = coordt(jj,ii+1)
        end do
     end do
     call CCTK_VarIndex(vindex,"mhd_evolve::u0")
     call interp_driver_carp(cctkGH,narr,pos,vindex,u0_int)
     call CCTK_VarIndex(vindex,"mhd_evolve::vx")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vx_int)
     call CCTK_VarIndex(vindex,"mhd_evolve::vy")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vy_int)
     call CCTK_VarIndex(vindex,"mhd_evolve::vz")
     slopet(:,1)=1.d0/u0_int
     slopet(:,2)=vx_int
     slopet(:,3)=vy_int
     slopet(:,4)=vz_int
     !call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !  ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !  Xglobal,Yglobal,Zglobal,pos,1.d0/u0,slopet(:,1))
     !call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !     ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !     Xglobal,Yglobal,Zglobal,pos,vx,slopet(:,2))
     !call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !     ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !     Xglobal,Yglobal,Zglobal,pos,vy,slopet(:,3)) 
     !call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !     ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !     Xglobal,Yglobal,Zglobal,pos,vz,slopet(:,4)) 
     coord=coord+dTo6*(slope+slopet+2.0d0*slopem)
  endif
  vx_p=vx
  vy_p=vy
  vz_p=vz
  u0_p=u0
end subroutine particle_tracer
