!------------------------------------------------------------
! Update all auxiliary MHD variables after an AMR regridding
!------------------------------------------------------------
!
subroutine mhd_post_regrid(cctk_dim,cctk_gsh,cctk_lsh,cctk_lbnd,cctk_ubnd,cctk_lssh,cctk_from,cctk_to,cctk_bbox,cctk_delta_time,cct&
  &k_time,cctk_delta_space,cctk_origin_space,cctk_levfac,cctk_levoff,cctk_levoffdenom,cctk_timefac,cctk_convlevel,cctk_convfac,cctk&
  &_nghostzones,cctk_iteration,cctkGH, X0em_fijs,X0em_ftijs,X0rad_fijs,X0rad_fijs_diag,X1em_fijs,X1em_ftijs,X1rad_fijs,X1rad_fijs_d&
  &iag,X2em_fijs,X2em_ftijs,X2rad_fijs,X2rad_fijs_diag,S_radx_flux,S_radx_flux_x,S_radx_flux_xp1,S_radx_nue_flux,S_radx_nux_flux,S_&
  &rady_flux,S_rady_nue_flux,S_rady_nux_flux,S_radz_flux,S_radz_nue_flux,S_radz_nux_flux,ftxy,ftxz,ftyx,ftyz,ftzx,ftzy,fxy,fxz,fyx,&
  &fyz,fzx,fzy,tau_rad_advect_diag,tau_rad_flux,tau_rad_flux_diag,tau_rad_flux_x,tau_rad_nue_flux,tau_rad_nux_flux,tau_rad_scalar_d&
  &iag, X0BSSN_AH,X0BSSN_aux_private,X0BSSN_aux_restrict2,X0BSSN_diag_restrict,X0BSSN_gupij,X0BSSN_matter,X0BSSN_refbd,X0BSSN_rhs,X&
  &0BSSN_vars,X0NPphi2_em,X0OmegF,X0Psi4_gw,X0Stagger_Bs,X0ah_radius,X0bh_posns,X0diag_vacuum_gf,X0disk_atmosphere,X0em_Ax,X0em_Ay,&
  &X0em_Az,X0em_Blagrangemultiplier,X0em_Blagrangemultiplier_rhs,X0em_Phi,X0em_Phi_rhs,X0em_conservativex,X0em_conservativey,X0em_c&
  &onservativez,X0em_rhsx,X0em_rhsy,X0em_rhsz,X0eos_params1,X0eos_params2,X0excision_int_gfs,X0field_line_variables,X0field_line_va&
  &riables_rhs,X0fisheye_vars,X0gw_moment_arrays,X0have_global_bdry,X0lapse_derivatives,X0lapse_rhs,X0lapse_vars,X0lapse_vars_aux,X&
  &0metric_spher_pol_1,X0metric_spher_pol_2,X0metric_spher_pol_3,X0mhd_conservatives,X0mhd_nosync,X0mhd_primitives,X0mhd_rhs,X0mhd_&
  &sync_lr_B_quantities,X0mhd_sync_lr_hydro_quantities,X0mhd_sync_metric_facevals,X0mhd_sync_nablas,X0mhd_sync_nablas_diag,X0mhd_sy&
  &nc_rho_br_rho_bl,X0mhd_temps,X0mhd_vs,X0micphys_conservatives,X0micphys_conservatives_rhs,X0micphys_sync_lr_hydro_quantities,X0m&
  &icphys_sync_nablas,X0microphys_primitives,X0phi_derivs,X0rad_conservatives,X0rad_conservatives_rhs,X0rad_pressure,X0rad_primitiv&
  &es,X0rad_sync_lr,X0rhovecs,X0shift_rhs,X0shift_vars,X0shift_vars_temp,X0surf_params,X1BSSN_AH,X1BSSN_aux_private,X1BSSN_aux_rest&
  &rict2,X1BSSN_diag_restrict,X1BSSN_gupij,X1BSSN_matter,X1BSSN_refbd,X1BSSN_rhs,X1BSSN_vars,X1NPphi2_em,X1OmegF,X1Psi4_gw,X1Stagge&
  &r_Bs,X1ah_radius,X1diag_vacuum_gf,X1disk_atmosphere,X1em_Ax,X1em_Ay,X1em_Az,X1em_Blagrangemultiplier,X1em_Blagrangemultiplier_rh&
  &s,X1em_Phi,X1em_Phi_rhs,X1em_conservativex,X1em_conservativey,X1em_conservativez,X1em_rhsx,X1em_rhsy,X1em_rhsz,X1excision_int_gf&
  &s,X1field_line_variables,X1field_line_variables_rhs,X1fisheye_vars,X1have_global_bdry,X1lapse_derivatives,X1lapse_rhs,X1lapse_va&
  &rs,X1lapse_vars_aux,X1mhd_conservatives,X1mhd_nosync,X1mhd_primitives,X1mhd_rhs,X1mhd_sync_lr_B_quantities,X1mhd_sync_lr_hydro_q&
  &uantities,X1mhd_sync_metric_facevals,X1mhd_sync_nablas,X1mhd_sync_nablas_diag,X1mhd_sync_rho_br_rho_bl,X1mhd_temps,X1mhd_vs,X1mi&
  &cphys_conservatives,X1micphys_conservatives_rhs,X1micphys_sync_lr_hydro_quantities,X1micphys_sync_nablas,X1microphys_primitives,&
  &X1phi_derivs,X1rad_conservatives,X1rad_conservatives_rhs,X1rad_pressure,X1rad_primitives,X1rad_sync_lr,X1shift_rhs,X1shift_vars,&
  &X1shift_vars_temp,X2BSSN_AH,X2BSSN_aux_private,X2BSSN_aux_restrict2,X2BSSN_diag_restrict,X2BSSN_gupij,X2BSSN_matter,X2BSSN_refbd&
  &,X2BSSN_rhs,X2BSSN_vars,X2NPphi2_em,X2OmegF,X2Psi4_gw,X2Stagger_Bs,X2ah_radius,X2diag_vacuum_gf,X2disk_atmosphere,X2em_Ax,X2em_A&
  &y,X2em_Az,X2em_Blagrangemultiplier,X2em_Blagrangemultiplier_rhs,X2em_Phi,X2em_Phi_rhs,X2em_conservativex,X2em_conservativey,X2em&
  &_conservativez,X2em_rhsx,X2em_rhsy,X2em_rhsz,X2excision_int_gfs,X2field_line_variables,X2field_line_variables_rhs,X2fisheye_vars&
  &,X2lapse_derivatives,X2lapse_rhs,X2lapse_vars,X2lapse_vars_aux,X2mhd_conservatives,X2mhd_nosync,X2mhd_primitives,X2mhd_rhs,X2mhd&
  &_sync_lr_B_quantities,X2mhd_sync_lr_hydro_quantities,X2mhd_sync_metric_facevals,X2mhd_sync_nablas,X2mhd_sync_nablas_diag,X2mhd_s&
  &ync_rho_br_rho_bl,X2mhd_temps,X2mhd_vs,X2micphys_conservatives,X2micphys_conservatives_rhs,X2micphys_sync_lr_hydro_quantities,X2&
  &micphys_sync_nablas,X2microphys_primitives,X2phi_derivs,X2rad_conservatives,X2rad_conservatives_rhs,X2rad_pressure,X2rad_primiti&
  &ves,X2rad_sync_lr,X2shift_rhs,X2shift_vars,X2shift_vars_temp,ah_centroid_length,ah_flags_length,ah_origin_length,ah_radius_lengt&
  &h,Aupxx,Aupxy,Aupxz,Aupyy,Aupyz,Aupzz,Ax,Ax_p,Ax_p_p,Axx,Axx_p,Axx_p_p,Axx_rhs,Axy,Axy_p,Axy_p_p,Axy_rhs,Axz,Axz_p,Axz_p_p,Axz_r&
  &hs,Ay,Ay_p,Ay_p_p,Ayy,Ayy_p,Ayy_p_p,Ayy_rhs,Ayz,Ayz_p,Ayz_p_p,Ayz_rhs,Az,Az_p,Az_p_p,Azz,Azz_p,Azz_p_p,Azz_rhs,Bfreq1,Bfreq1_p,B&
  &freq1_p_p,Bfreq2,Bfreq2_p,Bfreq2_p_p,Blagrangemultiplier,Blagrangemultiplier_p,Blagrangemultiplier_p_p,Blagrangemultiplier_rhs,B&
  &ox1X_VolInt,Box1X_VolInt1,Box1X_VolInt2,Box1Y_VolInt,Box1Y_VolInt1,Box1Y_VolInt2,Box1Z_VolInt,Box1Z_VolInt1,Box1Z_VolInt2,Box1de&
  &nom_VolInt,Box1denom_VolInt1,Box1denom_VolInt2,Bx,Bx_stagger,Bxl,Bxr,Bxtilde,Bxtilde_p,Bxtilde_p_p,Bxtilde_or_Ax_rhs,By,By_stagg&
  &er,Byl,Byr,Bytilde,Bytilde_p,Bytilde_p_p,Bytilde_or_Ay_rhs,Bz,Bz_stagger,Bzl,Bzr,Bztilde,Bztilde_p,Bztilde_p_p,Bztilde_or_Az_rhs&
  &,CoM_VolInt_denominator,CoMx_VolInt,CoMy_VolInt,CoMz_VolInt,E_rad,E_rad_nue,E_rad_nuel,E_rad_nuer,E_rad_nux,E_rad_nuxl,E_rad_nux&
  &r,E_radl,E_radr,Ex,Ey,Ez,F_rad0,F_rad0_nue,F_rad0_nux,F_rad_scalar,F_rad_scalar_nue,F_rad_scalar_nux,F_radx,F_radx_nue,F_radx_nu&
  &ele,F_radx_nuer,F_radx_nux,F_radx_nuxle,F_radx_nuxr,F_radxle,F_radxr,F_rady,F_rady_nue,F_rady_nuele,F_rady_nuer,F_rady_nux,F_rad&
  &y_nuxle,F_rady_nuxr,F_radyle,F_radyr,F_radz,F_radz_nue,F_radz_nuele,F_radz_nuer,F_radz_nux,F_radz_nuxle,F_radz_nuxr,F_radzle,F_r&
  &adzr,FaFal,FaFal_nue,FaFal_nux,FaFar,FaFar_nue,FaFar_nux,Fontfix_tracker_gf,Gammax,Gammax_p,Gammax_p_p,Gammax_drive,Gammax_rhs,G&
  &ammaxxx,Gammaxxy,Gammaxxz,Gammaxyy,Gammaxyz,Gammaxzz,Gammay,Gammay_p,Gammay_p_p,Gammay_drive,Gammay_rhs,Gammayxx,Gammayxy,Gammay&
  &xz,Gammayyy,Gammayyz,Gammayzz,Gammaz,Gammaz_p,Gammaz_p_p,Gammaz_drive,Gammaz_rhs,Gammazxx,Gammazxy,Gammazxz,Gammazyy,Gammazyz,Ga&
  &mmazzz,Gamx_const_VolInt,Gamy_const_VolInt,Gamz_const_VolInt,Ham_const_VolIntD,Ham_const_VolIntN,Ham_const_excised_VolIntD,Ham_c&
  &onst_excised_VolIntN,Ham_const_excised_innerregion_VolIntD,Ham_const_excised_innerregion_VolIntN,J_ADM_VolInt,J_ADM_VolInt_inner&
  &,J_ADM_VolInt_inner2,Jz_constraint,Jz_constraint_inner,KK,K_pp_rot1,K_pp_rot2,K_pp_rot3,K_rp_rot1,K_rp_rot2,K_rp_rot3,K_rr_rot1,&
  &K_rr_rot2,K_rr_rot3,K_rth_rot1,K_rth_rot2,K_rth_rot3,K_thp_rot1,K_thp_rot2,K_thp_rot3,K_thth_rot1,K_thth_rot2,K_thth_rot3,M0_AH_&
  &VolInt,M0_VolInt,M0_escape100M,M0_escape30M,M0_escape50M,M0_escape70M,M0_horiz_VolInt,M0_r1_VolInt,M0_r2_VolInt,M0_r3_VolInt,MNo&
  &rm,MONOPOLE,MRsx,MRsy,MRsz,M_ADM_VolInt,M_constraint,M_constraint_inner,Minternal_VolInt,Minternal_cold_VolInt,NPphi0im,NPphi0im&
  &_p,NPphi0im_p_p,NPphi0re,NPphi0re_p,NPphi0re_p_p,NPphi2im,NPphi2im_p,NPphi2im_p_p,NPphi2re,NPphi2re_p,NPphi2re_p_p,OS_surf_rad,P&
  &,P_cld,P_constraint,P_constraint_inner,P_floor,P_radxx,P_radxx_nue,P_radxx_nux,P_radxy,P_radxy_nue,P_radxy_nux,P_radxz,P_radxz_n&
  &ue,P_radxz_nux,P_radyy,P_radyy_nue,P_radyy_nux,P_radyz,P_radyz_nue,P_radyz_nux,P_radzz,P_radzz_nue,P_radzz_nux,P_tab,P_thermal,P&
  &hysicalRadius,Pl,Pr,PsiNorm,PsiRes,PsiTau,Pvec,RadiusDerivative,RadiusDerivative2,Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,S,SPx,SPx_p,SPx_p_p,SP&
  &y,SPy_p,SPy_p_p,SPz,SPz_p,SPz_p_p,S_rad_x,S_rad_x_p,S_rad_x_p_p,S_rad_x_nue,S_rad_x_nue_p,S_rad_x_nue_p_p,S_rad_x_nue_rhs,S_rad_&
  &x_nux,S_rad_x_nux_p,S_rad_x_nux_p_p,S_rad_x_nux_rhs,S_rad_x_rhs,S_rad_y,S_rad_y_p,S_rad_y_p_p,S_rad_y_nue,S_rad_y_nue_p,S_rad_y_&
  &nue_p_p,S_rad_y_nue_rhs,S_rad_y_nux,S_rad_y_nux_p,S_rad_y_nux_p_p,S_rad_y_nux_rhs,S_rad_y_rhs,S_rad_z,S_rad_z_p,S_rad_z_p_p,S_ra&
  &d_z_nue,S_rad_z_nue_p,S_rad_z_nue_p_p,S_rad_z_nue_rhs,S_rad_z_nux,S_rad_z_nux_p,S_rad_z_nux_p_p,S_rad_z_nux_rhs,S_rad_z_rhs,Sx,S&
  &xx,Sxy,Sxz,Sy,Syy,Syz,Sz,Szz,T_0x,T_0x_p,T_0x_p_p,T_0y,T_0y_p,T_0y_p_p,T_0z,T_0z_p,T_0z_p_p,T_VolInt,T_fluid,T_fluidle,T_fluidr,&
  &Tem0_0_VolInt,Tem0_0_outsideBH_VolInt,Tfluid0_0_VolInt,Tfluid0_0_outsideBH_VolInt,VolIntegrand,VolIntegrand2,VolIntegrand3,VolIn&
  &tegrand4,Y_e,Y_el,Y_er,ah_centroid_iteration,ah_centroid_iteration_p,ah_centroid_t,ah_centroid_t_p,ah_centroid_valid,ah_centroid&
  &_valid_p,ah_centroid_x,ah_centroid_x_p,ah_centroid_y,ah_centroid_y_p,ah_centroid_z,ah_centroid_z_p,ah_found_flag,ah_initial_find&
  &_flag,ah_origin_x,ah_origin_y,ah_origin_z,ah_radius,ah_really_initial_find_flag,ah_search_flag,b_phi_VolInt,bh_posn_x,bh_posn_y,&
  &bh_posn_z,bh_radius_z,brem_qee_VolInt,brem_qei_VolInt,chi,chi_p,chi_p_p,chi_rad,chi_rad_nue,chi_rad_nux,chi_rhs,cmax,cmax_rad,cm&
  &ax_rad_nue,cmax_rad_nux,cmin,cmin_rad,cmin_rad_nue,cmin_rad_nux,dBx_m,dBy_m,dBz_m,dP_m,dT_fluid_m,ddrbrvec,density_modes_i1,dens&
  &ity_modes_i2,density_modes_i3,density_modes_i4,density_modes_i5,density_modes_i6,density_modes_r0,density_modes_r1,density_modes&
  &_r2,density_modes_r3,density_modes_r4,density_modes_r5,density_modes_r6,div_beta,drbrvec,drhoYe_m,drhoYe_m_x,drhoYe_m_xp1,drho_b&
  &_m,drho_b_m_x,drho_b_m_xp1,dvx_m,dvx_m_x,dvx_m_xp1,dvy_m,dvz_m,em_J_between_VolInt,em_energy2_VolInt,em_energy_VolInt,em_energy_&
  &between_VolInt,em_energy_outsideBH_VolInt,em_energy_outsideradius1_VolInt,em_energy_outsideradius2_VolInt,emission_gf,emission_g&
  &f_nue,emission_gf_nux,eps_cld,eps_tab,eps_thermal,eps_tot,eta_nue,excision_zone_gf,fluid_J_horiz_VolInt,fluid_J_r1_VolInt,fluid_&
  &J_r2_VolInt,fluid_J_r3_VolInt,fluid_energy_horiz_VolInt,fluid_energy_r1_VolInt,fluid_energy_r2_VolInt,fluid_energy_r3_VolInt,gam&
  &ma_tab,gconx,gcony,gconz,gupxx,gupxx_f,gupxy,gupxy_f,gupxz,gupxz_f,gupyy,gupyy_f,gupyz,gupyz_f,gupzz,gupzz_f,gxx,gxx_p,gxx_p_p,g&
  &xx_f,gxx_rhs,gxxx,gxxy,gxxz,gxy,gxy_p,gxy_p_p,gxy_f,gxy_rhs,gxyx,gxyy,gxyz,gxz,gxz_p,gxz_p_p,gxz_f,gxz_rhs,gxzx,gxzy,gxzz,gyy,gy&
  &y_p,gyy_p_p,gyy_f,gyy_rhs,gyyx,gyyy,gyyz,gyz,gyz_p,gyz_p_p,gyz_f,gyz_rhs,gyzx,gyzy,gyzz,gzz,gzz_p,gzz_p_p,gzz_f,gzz_rhs,gzzx,gzz&
  &y,gzzz,h,h_p,half_b2_u0_VolInt,half_b2_u0_outsideBH_VolInt,have_global_bdry_max,have_global_bdry_min,int_momentsi,int_momentsr,k&
  &_tab,ka_gf,ka_gf_nue,ka_gf_nux,ks_gf,ks_gf_nue,ks_gf_nux,kxx,kxy,kxz,kyy,kyz,kzz,lapm1,lapm1_p,lapm1_p_p,lapm1_f,lapm1_rhs,lapm1&
  &_rot1,lapm1_rot2,lapm1_rot3,lapset,lapset_p,lapset_p_p,lapset_rhs,lapset_timederiv,lapsex,lapsey,lapsez,mass_star1,mass_star2,mh&
  &d_chi_line,mhd_chi_line_p,mhd_chi_line_p_p,mhd_chi_line_rhs,mhd_psi_line,mhd_psi_line_p,mhd_psi_line_p_p,mhd_psi_line_rhs,mhd_st&
  &_x,mhd_st_x_p,mhd_st_x_p_p,mhd_st_x_rhs,mhd_st_y,mhd_st_y_p,mhd_st_y_p_p,mhd_st_y_rhs,mhd_st_z,mhd_st_z_p,mhd_st_z_p_p,mhd_st_z_&
  &rhs,mhd_u_chi,mhd_u_chi_p,mhd_u_chi_p_p,mhd_u_chi_rhs,mhd_u_psi,mhd_u_psi_p,mhd_u_psi_p_p,mhd_u_psi_rhs,minternal_cold_horiz_Vol&
  &Int,minternal_cold_r1_VolInt,minternal_cold_r2_VolInt,minternal_cold_r3_VolInt,minternal_horiz_VolInt,minternal_r1_VolInt,minter&
  &nal_r2_VolInt,minternal_r3_VolInt,mom_const_VolIntD,mom_const_excised_VolIntD,momentsi,momentsi_old,momentsr,momentsr_old,momx_c&
  &onst_VolIntN,momx_const_excised_VolIntN,momy_const_VolIntN,momy_const_excised_VolIntN,momz_const_VolIntN,momz_const_excised_VolI&
  &ntN,monopole_VolInt,monopole_outsideBH_VolInt,mu_nu,n_poly,odd_momentsi,odd_momentsi_old,odd_momentsr,odd_momentsr_old,optd,optd&
  &_x,optd_y,optd_z,pfloor_gf,phi,phi_p,phi_p_p,phi_f,phi_rhs,phi_rot1,phi_rot2,phi_rot3,phix,phiy,phiz,psi,psi0im,psi0im_p,psi0im_&
  &p_p,psi0re,psi0re_p,psi0re_p_p,psi6phi,psi6phi_p,psi6phi_p_p,psi6phi_rhs,rad_const,rad_energy_VolInt,rad_energy_nue_VolInt,rad_e&
  &nergy_nux_VolInt,rbrvec,refbd,refbd_p,refbd_p_p,rho,rhoYe,rhoYe_p,rhoYe_p_p,rhoYe_rhs,rho_b,rho_b_atm_gf,rho_bl,rho_br,rho_star,&
  &rho_star_p,rho_star_p_p,rho_star_rhs,rho_tab,rhob_floor,rhovec,sbt,sbx,sby,sbz,shiftp_rot1,shiftp_rot2,shiftp_rot3,shiftr_rot1,s&
  &hiftr_rot2,shiftr_rot3,shiftth_rot1,shiftth_rot2,shiftth_rot3,shiftx,shiftx_p,shiftx_p_p,shiftx_f,shiftx_rhs,shiftxt,shiftxt_p,s&
  &hiftxt_p_p,shiftxt_rhs,shiftxt_timederiv,shifty,shifty_p,shifty_p_p,shifty_f,shifty_rhs,shiftyt,shiftyt_p,shiftyt_p_p,shiftyt_rh&
  &s,shiftyt_timederiv,shiftz,shiftz_p,shiftz_p_p,shiftz_f,shiftz_rhs,shiftzt,shiftzt_p,shiftzt_p_p,shiftzt_rhs,shiftzt_timederiv,s&
  &mallb2,st_x,st_y,st_z,surfvec,synch_VolInt,tau,tau_p,tau_p_p,tau_rad,tau_rad_p,tau_rad_p_p,tau_rad_nue,tau_rad_nue_p,tau_rad_nue&
  &_p_p,tau_rad_nue_rhs,tau_rad_nux,tau_rad_nux_p,tau_rad_nux_p_p,tau_rad_nux_rhs,tau_rad_rhs,tau_rhs,temp1,temp10,temp11,temp12,te&
  &mp13,temp14,temp15,temp16,temp17,temp18,temp19,temp2,temp20,temp21,temp22,temp23,temp24,temp25,temp26,temp3,temp4,temp5,temp6,te&
  &mp7,temp8,temp9,temp_g00,trK,trK_p,trK_p_p,trK_rhs,trR,trRtilde,u0,v02_rad_nuel,v02_rad_nuer,v02_rad_nuxl,v02_rad_nuxr,v02_radl,&
  &v02_radr,v02l,v02r,vvec,vx,vxl,vxr,vy,vyl,vyr,vz,vzl,vzr,w,zeta_rad,zeta_rad_nue,zeta_rad_nux, X0coordinates,X1coordinates,X2coo&
  &rdinates,coarse_dx,coarse_dy,coarse_dz,r,x,y,z)
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
  INTEGER X0em_fijs
  INTEGER X0em_ftijs
  INTEGER X0rad_fijs
  INTEGER X0rad_fijs_diag
  INTEGER X1em_fijs
  INTEGER X1em_ftijs
  INTEGER X1rad_fijs
  INTEGER X1rad_fijs_diag
  INTEGER X2em_fijs
  INTEGER X2em_ftijs
  INTEGER X2rad_fijs
  INTEGER X2rad_fijs_diag
  REAL*8 S_radx_flux(X0rad_fijs,X1rad_fijs,X2rad_fijs)
  REAL*8 S_radx_flux_x(X0rad_fijs_diag,X1rad_fijs_diag,X2rad_fijs_diag)
  REAL*8 S_radx_flux_xp1(X0rad_fijs_diag,X1rad_fijs_diag,X2rad_fijs_diag)
  REAL*8 S_radx_nue_flux(X0rad_fijs,X1rad_fijs,X2rad_fijs)
  REAL*8 S_radx_nux_flux(X0rad_fijs,X1rad_fijs,X2rad_fijs)
  REAL*8 S_rady_flux(X0rad_fijs,X1rad_fijs,X2rad_fijs)
  REAL*8 S_rady_nue_flux(X0rad_fijs,X1rad_fijs,X2rad_fijs)
  REAL*8 S_rady_nux_flux(X0rad_fijs,X1rad_fijs,X2rad_fijs)
  REAL*8 S_radz_flux(X0rad_fijs,X1rad_fijs,X2rad_fijs)
  REAL*8 S_radz_nue_flux(X0rad_fijs,X1rad_fijs,X2rad_fijs)
  REAL*8 S_radz_nux_flux(X0rad_fijs,X1rad_fijs,X2rad_fijs)
  REAL*8 ftxy(X0em_ftijs,X1em_ftijs,X2em_ftijs)
  REAL*8 ftxz(X0em_ftijs,X1em_ftijs,X2em_ftijs)
  REAL*8 ftyx(X0em_ftijs,X1em_ftijs,X2em_ftijs)
  REAL*8 ftyz(X0em_ftijs,X1em_ftijs,X2em_ftijs)
  REAL*8 ftzx(X0em_ftijs,X1em_ftijs,X2em_ftijs)
  REAL*8 ftzy(X0em_ftijs,X1em_ftijs,X2em_ftijs)
  REAL*8 fxy(X0em_fijs,X1em_fijs,X2em_fijs)
  REAL*8 fxz(X0em_fijs,X1em_fijs,X2em_fijs)
  REAL*8 fyx(X0em_fijs,X1em_fijs,X2em_fijs)
  REAL*8 fyz(X0em_fijs,X1em_fijs,X2em_fijs)
  REAL*8 fzx(X0em_fijs,X1em_fijs,X2em_fijs)
  REAL*8 fzy(X0em_fijs,X1em_fijs,X2em_fijs)
  REAL*8 tau_rad_advect_diag(X0rad_fijs_diag,X1rad_fijs_diag,X2rad_fijs_diag)
  REAL*8 tau_rad_flux(X0rad_fijs,X1rad_fijs,X2rad_fijs)
  REAL*8 tau_rad_flux_diag(X0rad_fijs_diag,X1rad_fijs_diag,X2rad_fijs_diag)
  REAL*8 tau_rad_flux_x(X0rad_fijs_diag,X1rad_fijs_diag,X2rad_fijs_diag)
  REAL*8 tau_rad_nue_flux(X0rad_fijs,X1rad_fijs,X2rad_fijs)
  REAL*8 tau_rad_nux_flux(X0rad_fijs,X1rad_fijs,X2rad_fijs)
  REAL*8 tau_rad_scalar_diag(X0rad_fijs_diag,X1rad_fijs_diag,X2rad_fijs_diag)
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
  REAL*8  Erad_atm_cut
  REAL*8  Erad_cut
  REAL*8  K_poly
  REAL*8  M_B
  REAL*8  OS_T_Rmax
  REAL*8  P_fact
  REAL*8  P_max
  REAL*8  Psi6threshold
  REAL*8  Sym_Bz
  REAL*8  T_fluid_cgs_atm
  REAL*8  c_h_default
  REAL*8  c_p_default
  REAL*8  dxvec
  REAL*8  eps_thermal_bhns
  REAL*8  ergo_sigma
  REAL*8  gamma_OS
  REAL*8  gamma_th
  REAL*8  lambda_line
  REAL*8  mhd_bigp
  REAL*8  mhd_mbh
  REAL*8  mhd_mdot
  REAL*8  mhd_r_crit
  REAL*8  min_BH_radius
  REAL*8  pfloor
  REAL*8  r_ex_line
  REAL*8  rad_T_cutoff
  REAL*8  rad_T_fac
  REAL*8  rad_T_floor
  REAL*8  rad_T_pow
  REAL*8  rad_emissivity_abs
  REAL*8  rad_emissivity_sct
  REAL*8  rad_opacity_abs
  REAL*8  rad_opacity_sct
  REAL*8  rho_b_atm
  REAL*8  rho_b_max
  REAL*8  rho_fact
  REAL*8  rhos_max
  REAL*8  sdots_o_rhot
  REAL*8  t_cool
  REAL*8  tau_atm
  REAL*8  tau_fact
  INTEGER*4 A_BC_rad_mag_bondi_enable
  INTEGER*4 EM_BC
  INTEGER*4 MHD_MaxNumConstrainedVars
  INTEGER*4 MHD_MaxNumEvolvedVars
  INTEGER*4 Matter_BC
  INTEGER*4 Reconstruction
  INTEGER*4 adm_ppm_b
  INTEGER*4 allow_negative_eps_th
  INTEGER*4 artificial_cooling_enable
  INTEGER*4 compute_microphysics
  INTEGER*4 constrained_transport_scheme
  INTEGER*4 cooling_in_St_eq
  INTEGER*4 em_evolve_enable
  INTEGER*4 em_gauge
  INTEGER*4 enable_HARM_energyvariable
  INTEGER*4 enable_OS_collapse
  INTEGER*4 enable_disk_em_flux_induction
  INTEGER*4 enable_primitives_disk
  INTEGER*4 enable_shocktest_primitive_mode
  INTEGER*4 enable_trace_field_line
  INTEGER*4 eps_flag
  INTEGER*4 ergo_star
  INTEGER*4 force_font_fix_fail
  INTEGER*4 horizon_enforce_rho_profile
  INTEGER*4 hyperbolic_divergence_cleaning_centered_differencing
  INTEGER*4 hyperbolic_divergence_cleaning_enable
  INTEGER*4 m
  INTEGER*4 microphysics_scheme
  INTEGER*4 neos
  INTEGER*4 nrhovec
  INTEGER*4 primitives_debug
  INTEGER*4 primitives_solver
  INTEGER*4 rad_closure_scheme
  INTEGER*4 rad_evolve_enable
  INTEGER*4 rad_fix
  INTEGER*4 rad_fourforce_enable
  INTEGER*4 reconstruct_Bitildes_instead_of_Bis
  INTEGER*4 reconstruct_Pthermal_instead_of_P
  INTEGER*4 tau_stildefix_enable
  INTEGER*4 use_HARM_primitives
  INTEGER*4 use_central_scheme_instead_of_hll
  INTEGER*4 use_new_code
  COMMON /mhd_evolverest/Erad_atm_cut,Erad_cut,K_poly,M_B,OS_T_Rmax,P_fact,P_max,Psi6threshold,Sym_Bz,T_fluid_cgs_atm,c_h_default,c&
  &_p_default,dxvec,eps_thermal_bhns,ergo_sigma,gamma_OS,gamma_th,lambda_line,mhd_bigp,mhd_mbh,mhd_mdot,mhd_r_crit,min_BH_radius,pf&
  &loor,r_ex_line,rad_T_cutoff,rad_T_fac,rad_T_floor,rad_T_pow,rad_emissivity_abs,rad_emissivity_sct,rad_opacity_abs,rad_opacity_sc&
  &t,rho_b_atm,rho_b_max,rho_fact,rhos_max,sdots_o_rhot,t_cool,tau_atm,tau_fact,A_BC_rad_mag_bondi_enable,EM_BC,MHD_MaxNumConstrain&
  &edVars,MHD_MaxNumEvolvedVars,Matter_BC,Reconstruction,adm_ppm_b,allow_negative_eps_th,artificial_cooling_enable,compute_microphy&
  &sics,constrained_transport_scheme,cooling_in_St_eq,em_evolve_enable,em_gauge,enable_HARM_energyvariable,enable_OS_collapse,enabl&
  &e_disk_em_flux_induction,enable_primitives_disk,enable_shocktest_primitive_mode,enable_trace_field_line,eps_flag,ergo_star,force&
  &_font_fix_fail,horizon_enforce_rho_profile,hyperbolic_divergence_cleaning_centered_differencing,hyperbolic_divergence_cleaning_e&
  &nable,m,microphysics_scheme,neos,nrhovec,primitives_debug,primitives_solver,rad_closure_scheme,rad_evolve_enable,rad_fix,rad_fou&
  &rforce_enable,reconstruct_Bitildes_instead_of_Bis,reconstruct_Pthermal_instead_of_P,tau_stildefix_enable,use_HARM_primitives,use&
  &_central_scheme_instead_of_hll,use_new_code
  REAL*8  KO_hydro_inside_BH
  REAL*8  KO_inside_BH
  REAL*8  KO_refbd
  REAL*8  damp_lorenz
  INTEGER*4 number_of_points_apply_KO_outside_refboundaries
  INTEGER*4 really_want_gamma_th_of_1_point_3
  INTEGER*4 use_disk_bcs
  INTEGER*4 use_std_PPM_instead_of_efficient_3_gz_version
  COMMON /mhd_evolvepriv/KO_hydro_inside_BH,KO_inside_BH,KO_refbd,damp_lorenz,number_of_points_apply_KO_outside_refboundaries,reall&
  &y_want_gamma_th_of_1_point_3,use_disk_bcs,use_std_PPM_instead_of_efficient_3_gz_version
  INTEGER*4 CCTKH0
  INTEGER*4 iteration_to_insert_rad
  INTEGER*4 CCTKH1
  COMMON /BHNSrest/CCTKH0,iteration_to_insert_rad,CCTKH1
  REAL*8  CCTKH7
  REAL*8  CCTKH8
  REAL*8  CCTKH9
  REAL*8  CCTKH10
  REAL*8  CCTKH11
  REAL*8  CCTKH14
  REAL*8  CCTKH15
  REAL*8  CCTKH16
  REAL*8  CCTKH17
  REAL*8  CCTKH18
  REAL*8  CCTKH19
  REAL*8  CCTKH20
  REAL*8  CCTKH21
  REAL*8  CCTKH22
  REAL*8  CCTKH24
  REAL*8  CCTKH25
  REAL*8  CCTKH26
  INTEGER*4 CCTKH2
  INTEGER*4 CCTKH3
  INTEGER*4 CCTKH4
  INTEGER*4 CCTKH5
  INTEGER*4 CCTKH6
  INTEGER*4 CCTKH12
  INTEGER*4 CCTKH13
  INTEGER*4 CCTKH23
  INTEGER*4 zero_out_matter_source_terms
  COMMON /BSSNrest/CCTKH7,CCTKH8,CCTKH9,CCTKH10,CCTKH11,CCTKH14,CCTKH15,CCTKH16,CCTKH17,CCTKH18,CCTKH19,CCTKH20,CCTKH21,CCTKH22,CCT&
  &KH24,CCTKH25,CCTKH26,CCTKH2,CCTKH3,CCTKH4,CCTKH5,CCTKH6,CCTKH12,CCTKH13,CCTKH23,zero_out_matter_source_terms
  REAL*8  CCTKH27
  REAL*8  CCTKH33
  REAL*8  CCTKH34
  REAL*8  CCTKH35
  REAL*8  CCTKH36
  REAL*8  CCTKH40
  REAL*8  CCTKH45
  REAL*8  CCTKH46
  REAL*8  CCTKH48
  INTEGER*4 CCTKH28
  INTEGER*4 CCTKH29
  INTEGER*4 CCTKH30
  INTEGER*4 CCTKH31
  INTEGER*4 CCTKH32
  INTEGER*4 CCTKH37
  INTEGER*4 CCTKH38
  INTEGER*4 CCTKH39
  INTEGER*4 CCTKH41
  INTEGER*4 CCTKH42
  INTEGER*4 num_BHs
  INTEGER*4 CCTKH43
  INTEGER*4 CCTKH44
  INTEGER*4 CCTKH47
  INTEGER*4 CCTKH49
  COMMON /DIAGNOSTICS_VACUUMrest/CCTKH27,CCTKH33,CCTKH34,CCTKH35,CCTKH36,CCTKH40,CCTKH45,CCTKH46,CCTKH48,CCTKH28,CCTKH29,CCTKH30,CC&
  &TKH31,CCTKH32,CCTKH37,CCTKH38,CCTKH39,CCTKH41,CCTKH42,num_BHs,CCTKH43,CCTKH44,CCTKH47,CCTKH49
  REAL*8  CCTKH50
  REAL*8  CCTKH51
  REAL*8  CCTKH52
  REAL*8  CCTKH53
  REAL*8  CCTKH63
  REAL*8  CCTKH64
  REAL*8  CCTKH65
  REAL*8  CCTKH66
  REAL*8  CCTKH67
  REAL*8  CCTKH68
  REAL*8  CCTKH69
  REAL*8  CCTKH70
  integer*8  bitant_plane
  integer*8  domain
  integer*8  CCTKH54
  integer*8  CCTKH55
  integer*8  CCTKH62
  INTEGER*4 CCTKH56
  INTEGER*4 CCTKH57
  INTEGER*4 CCTKH58
  INTEGER*4 CCTKH59
  INTEGER*4 CCTKH60
  INTEGER*4 CCTKH61
  COMMON /GRIDrest/CCTKH50,CCTKH51,CCTKH52,CCTKH53,CCTKH63,CCTKH64,CCTKH65,CCTKH66,CCTKH67,CCTKH68,CCTKH69,CCTKH70,bitant_plane,dom&
  &ain,CCTKH54,CCTKH55,CCTKH62,CCTKH56,CCTKH57,CCTKH58,CCTKH59,CCTKH60,CCTKH61
  REAL*8  CCTKH73
  REAL*8  CCTKH74
  REAL*8  CCTKH75
  integer*8  slicing_type
  INTEGER*4 CCTKH71
  INTEGER*4 CCTKH72
  INTEGER*4 CCTKH76
  INTEGER*4 CCTKH77
  INTEGER*4 CCTKH78
  INTEGER*4 CCTKH79
  INTEGER*4 CCTKH80
  COMMON /LAPSErest/CCTKH73,CCTKH74,CCTKH75,slicing_type,CCTKH71,CCTKH72,CCTKH76,CCTKH77,CCTKH78,CCTKH79,CCTKH80
  INTEGER*4 CCTKH81
  INTEGER*4 CCTKH82
  INTEGER*4 CCTKH83
  INTEGER*4 CCTKH84
  INTEGER*4 CCTKH85
  INTEGER*4 CCTKH86
  INTEGER*4 CCTKH87
  INTEGER*4 CCTKH88
  INTEGER*4 CCTKH89
  INTEGER*4 CCTKH90
  INTEGER*4 CCTKH91
  INTEGER*4 MoL_Num_Constrained_Vars
  INTEGER*4 MoL_Num_Evolved_Vars
  INTEGER*4 MoL_Num_SaveAndRestore_Vars
  INTEGER*4 CCTKH92
  COMMON /METHODOFLINESrest/CCTKH81,CCTKH82,CCTKH83,CCTKH84,CCTKH85,CCTKH86,CCTKH87,CCTKH88,CCTKH89,CCTKH90,CCTKH91,MoL_Num_Constra&
  &ined_Vars,MoL_Num_Evolved_Vars,MoL_Num_SaveAndRestore_Vars,CCTKH92
  INTEGER*4 Spatial_Gauge
  INTEGER*4 CCTKH93
  INTEGER*4 CCTKH94
  COMMON /SHIFTrest/Spatial_Gauge,CCTKH93,CCTKH94
  
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
  INTEGER*4 function Boundary_SelectGroupForBC (GH, faces, boundary_width, table_handle, var_name, bc_name)
  implicit none
  integer*8 GH
  INTEGER*4 faces
  INTEGER*4 boundary_width
  INTEGER*4 table_handle
  character(*) var_name
  character(*) bc_name
  end function Boundary_SelectGroupForBC
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
  INTEGER*4 function MoLRegisterConstrained (ConstrainedIndex)
  implicit none
  INTEGER*4 ConstrainedIndex
  end function MoLRegisterConstrained
  end interface
  interface
  INTEGER*4 function MoLRegisterEvolvedGroup (EvolvedIndex, RHSIndex)
  implicit none
  INTEGER*4 EvolvedIndex
  INTEGER*4 RHSIndex
  end function MoLRegisterEvolvedGroup
  end interface
  
  integer, dimension(3) :: ext,fake_ext,thisproc_have_global_bdry_min,thisproc_have_global_bdry_max
  real*8                :: dX,dY,dZ,levelnumber
  real*8                :: Xglobmin,Yglobmin,Zglobmin,Xglobmax,Yglobmax,Zglobmax
  real*8                :: rho_fail_max_step,M_fail_step,b2bt
  real*8                :: xNS1,yNS1,xNS2,yNS2
  real*8                :: tau_radl,S_rad_xl,S_rad_yl,S_rad_zl,Psi4,Psi6,alpn1,gupxx_phys,gupyy_phys,gupzz_phys,gupxy_phys,gupxz_ph&
  &ys,gupyz_phys
  real*8                :: au0,F_rad0_denom,F_rad0_num,F_rad0l,P_radl
  integer               :: proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax
  integer               :: glob_imax,glob_jmax,glob_kmax
  integer               :: i,j,k,vindex,NFont,NFont_L
  integer               :: handle,dummy,ierr,repairs_needed,repairs_rad_needed
  REAL*8             :: reduction_value
  integer :: STATUS = 0
  integer               :: AXISYM
  parameter(AXISYM = 4)
  ext = cctk_lsh
  levelnumber = cctk_levfac(1)
  levelnumber = log(levelnumber)/log(2.D0)+1.D0
  write(*,*) "At the beginning of driver_post_regrid. levelnumber is ", levelnumber
!  write(*,*) cctk_levfac(1) is , cctk_levfac(1)
  !write(*,*) cctk_levfac is , cctk_levfac 
  ! IMPORTANT NOTE:
  ! After an AMR regridding (i.e., AMR grid movement), Carpet sets only the evolved variables (e.g., gijs, Aijs) 
  !   at all points on all grids.
  ! However, auxiliary variables (such as, e.g., matter source terms, Rij, etc.) are set only in regions where the grid
  !   _was_, and not in new regions.  
  !
  ! To properly explain this, its best to draw a 1D picture of a gridfunction defined along a line parallel to the x-axis.
  ! Let + denote the gridpoints where the gridfunction is defined, X denote the gridpoints where the gridfunction is
  !   not defined, and | the grid boundary
  ! Before a regrid, heres what our grid looks like:
  ! |++++++++++++++++++++++++++|
  !  <---      x-axis      --->  
  ! After a regrid, the grid has moved a little to the right:
  !      |+++++++++++++++++++++XXXXX|
  !  <---      x-axis      --->
  !  Attempting to access a gridfunction at the X points will result in memory errors, and ultimately the simulation may crash.
  ! This function is designed to fill in all the gridpoints after a regrid, including the X points.
  if(cctk_iteration.ge.0) then
     dX = (cctk_delta_space(1)/cctk_levfac(1))
     dY = (cctk_delta_space(2)/cctk_levfac(2))
     dZ = (cctk_delta_space(3)/cctk_levfac(3))
 !    write(*,*) Start driver_post_regird.F90!!!
     ! Should check if there are NaNs first...
     !call fill_bssn_symmetry_gz_bssn_vars(ext,X,Y,Z,Symmetry,phi,chi,trK,Gammax,Gammay,Gammaz,gxx,gxy,gxz,gyy,gyz,gzz,Axx,Axy,Axz,A
     !     call fill_shift_symmetry_gz(ext,X,Y,Z,Symmetry,shiftx,shifty,shiftz,shiftxt,shiftyt,shiftzt)
     ! FOLLOWING CARTSYMS ARE NEEDED.  Otherwise, well end up with memory errors inside ghostzones (even BSSN quantities)!
     call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')
     call CartSymGN(dummy,cctkGH,'bssn::BSSN_vars')
     call CartSymGN(dummy,cctkGH,'shift::shift_vars')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_conservatives')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativex')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativey')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativez')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::em_Blagrangemultiplier')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_primitives')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vs')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_conservatives')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_primitives')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::microphys_primitives')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::micphys_conservatives')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_pressure')
!     call CartSymGN(dummy,cctkGH,mhd_evolve::temperatures)
     !write(*,*) checkpoint (2). levelnumber is , levelnumber
     ! gupij is not evolved, and so is not defined anywhere that the grid has moved.
     !  Here we recompute gupij from the newly updated gijs:
     call BSSN_compute_gupij(cctkGH,cctk_lsh, &
          gxx,gxy,gxz,gyy,gyz,gzz, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)
     !write(*,*) checkpoint (3). levelnumber is , levelnumber
     if (constrained_transport_scheme .ne. 3) then
        !$omp parallel do
        do k=1,cctk_lsh(3)
           do j=1,cctk_lsh(2)
              do i=1,cctk_lsh(1)
                 !Note that the primitives solver requires Bx,By,Bz as input.
                 Bx(i,j,k) = Bxtilde(i,j,k)
                 By(i,j,k) = Bytilde(i,j,k)
                 Bz(i,j,k) = Bztilde(i,j,k)
              end do
           end do
        end do
        !$omp end parallel do
        ! Notice that we set B^i = tildeB^i above.
        b2bt = -1.D0
        call convert_b(ext,Bx,By,Bz,phi,b2bt)
     end if
     !write(*,*) checkpoint (4). levelnumber is , levelnumber
     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              !Set face values of field values to zero, so that valgrind is happy.
              !Should probably comment these out if youre debugging the code
              gupxx_f(i,j,k) = 0.D0
              gupyy_f(i,j,k) = 0.D0
              gupzz_f(i,j,k) = 0.D0
              lapm1_f(i,j,k) = 0.D0
              shiftx_f(i,j,k) = 0.D0
              shifty_f(i,j,k) = 0.D0
              shiftz_f(i,j,k) = 0.D0
              gxx_f(i,j,k) = 0.D0
              gxy_f(i,j,k) = 0.D0
              gxz_f(i,j,k) = 0.D0
              gyy_f(i,j,k) = 0.D0
              gyz_f(i,j,k) = 0.D0
              gzz_f(i,j,k) = 0.D0
              phi_f(i,j,k) = 0.D0
              if (isnan(Ax(i,j,k)) .or. isnan(Ay(i,j,k)) .or. isnan(Az(i,j,k))) then
                 write (*,*) "Inside mhd_postinitialdata.F90, Ax(i,j,k) is nan!!!!"
              end if
           end do
        end do
     end do
     !$omp end parallel do
     !-----------------------------------------------------------------------------------------
     ! PRIMITIVE RECONSTRUCTION STEP:
     !write(*,*) checkpoint (5). levelnumber is , levelnumber
     glob_imax = ext(1)
     glob_jmax = ext(2)
     glob_kmax = ext(3)
     proc_imin = -100
     proc_jmin = -100
     proc_kmin = -100
     proc_imax = -1
     proc_jmax = -1
     proc_kmax = -1
     thisproc_have_global_bdry_min = have_global_bdry_min(int(levelnumber),:)
     thisproc_have_global_bdry_max = have_global_bdry_max(int(levelnumber),:)
     !write(*,*) checkpoint (6). levelnumber is , levelnumber
     if(thisproc_have_global_bdry_min(1)==1) proc_imin = 0
     if(thisproc_have_global_bdry_min(2)==1) proc_jmin = 0
     if(thisproc_have_global_bdry_min(3)==1) proc_kmin = 0
     if(thisproc_have_global_bdry_max(1)==1) proc_imax = glob_imax
     if(thisproc_have_global_bdry_max(2)==1) proc_jmax = glob_jmax
     if(thisproc_have_global_bdry_max(3)==1) proc_kmax = glob_kmax
! calcualte CoM for tempearture computation.
        if (Box1denom_VolInt == 0.0) then
        Box1denom_VolInt = 1.0e-10
        end if
        if (Box1denom_VolInt1 == 0.0) then
        Box1denom_VolInt2 = 1.0e-10
        end if
        if (Box1denom_VolInt1 == 0.0) then
        Box1denom_VolInt2 = 1.0e-10
        end if
     write(*,*) "HELLO in MHD POST_REGRID",levelnumber,thisproc_have_global_bdry_min,thisproc_have_global_bdry_max
     if (enable_primitives_disk == 0) then
        if(primitives_solver==11) then
           !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
                    ! Here we provide approximations to rho_b and P, as first guesses for the primitives solver
                    rho_b(i,j,k) = max(rho_star(i,j,k) * exp(-6.D0*phi(i,j,k)), rho_b_atm)
                    call compute_pcold_epscold(rho_b(i,j,k),P(i,j,k),temp1(i,j,k), &
                         neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab,enable_OS_collapse)
                    vx(i,j,k) = -shiftx(i,j,k)
                    vy(i,j,k) = -shifty(i,j,k)
                    vz(i,j,k) = -shiftz(i,j,k)
                    u0(i,j,k) = 1.D0/(lapm1(i,j,k)+1.d0)
                    !! Add h(i,j,k) before prim solver to aviod zero values at boundary.                                            
                    !h(i,j,k) = 1.d0 + (P(i,j,k) + (gamma_th-1.d0) * rho_b(i,j,k) * temp1(i,j,k)*(eps_thermal_bhns -1.d0))/rho_b(i,j
                 end do
              end do
           end do
           !$omp end parallel do
           !Here we use the HARM primitives solver, with a new prescription that minimizes changes to
           !  conservatives without applying the Font fix.
           !We intend to make the below function truly generic, but currently it only supports Gamma=2.
           !  It can be trivially extended for arbitrary Gamma-law EOSs, but will require work for more
           !  generic EOSs.
           !We hope the below function will eventually be used in place of other primitives solver options,
           !  since it can be easily extended to use your desired technique for solving primitives.
!           write(6,*) (1)Inside driver_post_regrid, before primitives_generic, gamma_tab =, gamma_tab
           !if (levelnumber == max_refinement_level) then
           !write(*,*) (1)levelnumber is at max refinement level!!!
           !else if (levelnumber.le.max_refinement_level-1) then
           call primitives_generic(cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z, &
                phi,gxx,gxy,gxz,gyy,gyz,gzz, &
                gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                lapm1,shiftx,shifty,shiftz, &
                Bx,By,Bz, &
                eps_tot, eps_thermal, eps_cld, P_cld,&
                ka_gf, ks_gf, emission_gf, chi_rad, chi_rad_nue, Y_e, optd, eta_nue,&
                ka_gf_nue, ks_gf_nue, emission_gf_nue, ka_gf_nux,ks_gf_nux, emission_gf_nux,&
                mhd_st_x,mhd_st_y,mhd_st_z,tau,rho_star, rhoYe, &
                vx,vy,vz,P,rho_b,h,u0, &
                rho_b_atm,tau_atm, rho_b_max,&
                xNS1,yNS1, xNS2, yNS2, M_B, rad_T_fac, rad_T_cutoff, rad_T_pow, rad_T_floor, T_fluid,&
                neos,ergo_star, ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
                temp1,temp2,temp3,temp4,temp5, &
                primitives_debug,Psi6threshold,rad_const,horizon_enforce_rho_profile,enable_OS_collapse,&
                rad_evolve_enable, compute_microphysics, &
                microphysics_scheme, T_fluid_cgs_atm, rad_fix)
           ! write(*,*) levelnumber is less than max_refinement level, use prim solver 
           !end if                
           !print *,after primitives_generic2 ***************
           !The goal of the below function is to just update the metric source terms.  We should
           !  really get rid of the old h,w,st_i,Ei gridfunctions...
           !Currently this only works for Gamma-law EOSs but could be extended for arbitary EOSs
           !  with some work.
           !if (levelnumber == max_refinement_level) then
           !write(*,*) max-level reached!
           !else if (levelnumber.le.max_refinement_level-1) then
           call metric_source_terms_and_misc_vars(cctkGH,cctk_lsh, &
                rho,Sx,Sy,Sz, &
                Sxx,Sxy,Sxz,Syy,Syz,Szz, &
                tau_rad, S_rad_x,S_rad_y, S_rad_z,&
                rhoYe,Y_e,&
                E_rad, F_radx, F_rady, F_radz, F_rad0, F_rad_scalar,&
                h,w,st_x,st_y,st_z, &
                Ex,Ey,Ez, &
                sbt,sbx,sby,sbz, &
                lapm1,shiftx,shifty,shiftz, &
                phi,gxx,gxy,gxz,gyy,gyz,gzz, &
                gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                P,rho_b,u0,vx,vy,vz, &
                Bx,By,Bz, &
                rho_star,mhd_st_x,mhd_st_y,mhd_st_z,&
                neos,ergo_star, ergo_sigma,&
                rho_tab, P_tab,eps_tab,k_tab,gamma_tab,&
                rho_b_atm, enable_OS_collapse, rad_evolve_enable, rad_closure_scheme, Psi6threshold, Erad_atm_cut)
!                write(*,*) IN driver_post_regrid.F90, after metric_source_terms_and_misc_vars, level number is,levelnumber
           !end if
!          write(*,*) HI tau=,tau(1,1,1),levelnumber
!           call compute_optical_depth(ext, X,Y,Z, dX,dY,dZ, optd, optd_x, optd_y, optd_z, ka_gf, ks_gf)  
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
           if (T_fluid(i,j,k) .gt. 1.0e-30 .or. eps_tot(i,j,k) .eq. 0.0 .or. rho_star(i,j,k) .le. 0.0) then
              write(*,*) " In driver_post_regrid (2), T_fluid too large or eps_tot = 0.0 or rho_star < 0.0!!!!"
              write(*,*) "i,j,k X,Y,Z, T_fluid(i,j,k) = ", i,j,k, X(i,j,k), Y(i,j,k), Z(i,j,k), T_fluid(i,j,k)
              write(*,*) "rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k)=",rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), v&
  &z(i,j,k)
              write(*,*) "eps_tot(i,j,k), eps_cld(i,j,k), eps_thermal(i,j,k), P_thermal(i,j,k)=",eps_tot(i,j,k), eps_cld(i,j,k), ep&
  &s_thermal(i,j,k), P_thermal(i,j,k)
              write(*,*) "rho_star(i,j,k), tau(i,j,k), mhd_st_x(i,j,k), mhd_st_y(i,j,k), mhd_st_z(i,j,k), u0(i,j,k)=", rho_star(i,j&
  &,k), tau(i,j,k), mhd_st_x(i,j,k), mhd_st_y(i,j,k), mhd_st_z(i,j,k), u0(i,j,k)
              write(*,*) "Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)=", Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)
!              call EXIT(STATUS)                                                                                                    
           end if
        end do
     end do
  end do
  if(rad_evolve_enable==1 .and. rad_closure_scheme == 1) then
!     call primitive_rad_guess(ext, tau_rad, S_rad_x, S_rad_y, S_rad_z, phi, lapm1, gupxx, gupyy, gupzz, gupxy, gupxz, gupyz, shiftx
 !    if(1==0) then
     do k=1,ext(3)
        do j=1,ext(2)
           do i=1,ext(1)
              ! Here we provide approximations to Erad and Frad using Eddington factor = 2/3 (as the midpoint between thick and thin
              tau_radl = tau_rad(i,j,k)
              S_rad_xl = S_rad_x(i,j,k)
              S_rad_yl = S_rad_y(i,j,k)
              S_rad_zl = S_rad_z(i,j,k)
              Psi4 = exp(4.d0*phi(i,j,k))
              Psi6 = exp(6.d0*phi(i,j,k))
              alpn1 = 1.d0 + lapm1(i,j,k)
              gupxx_phys = gupxx(i,j,k)/Psi4
              gupyy_phys = gupyy(i,j,k)/Psi4
              gupzz_phys = gupzz(i,j,k)/Psi4
              gupxy_phys = gupxy(i,j,k)/Psi4
              gupxz_phys = gupxz(i,j,k)/Psi4
              gupyz_phys = gupyz(i,j,k)/Psi4
              if (tau_radl < 0.0) then
                 E_rad(i,j,k) = Erad_atm_cut
                 F_radx(i,j,k) = 0.0
                 F_rady(i,j,k) = 0.0
                 F_radz(i,j,k) = 0.0
              else
                 au0 = alpn1*u0(i,j,k)
                 F_rad0_denom = alpn1 * Psi6 * (2.0*au0*au0 + 1.0)/(4.0*au0*au0 - 1.0)
                 F_rad0_num = 4.0*au0*(1.0 - au0*au0)/(4.0*au0*au0 - 1.0)*tau_radl +( (shiftx(i,j,k)*u0(i,j,k)+u0(i,j,k)*vx(i,j,k))&
  &*S_rad_xl +(shifty(i,j,k)*u0(i,j,k)+u0(i,j,k)*vy(i,j,k))*S_rad_yl+(shiftz(i,j,k)*u0(i,j,k)+u0(i,j,k)*vz(i,j,k))*S_rad_zl)
                 F_rad0l = F_rad0_num/F_rad0_denom
                 E_rad(i,j,k) = 3.0*(tau_radl/Psi6 - 2.0*alpn1*au0*F_rad0l)/(4.0*au0*au0 - 1.0)
                 P_radl = 2.0*E_rad(i,j,k)/3.0
                 F_radx(i,j,k) = ((gupxx_phys*S_rad_xl + gupxy_phys*S_rad_yl + gupxz_phys*S_rad_zl))/(alpn1*Psi6*u0(i,j,k))-(E_rad(&
  &i,j,k)+P_radl)*u0(i,j,k)*(vx(i,j,k)+shiftx(i,j,k)) - 2.0*(F_rad0l * shiftx(i,j,k)) - F_rad0l*vx(i,j,k)
                 F_rady(i,j,k) = ((gupxy_phys*S_rad_xl + gupyy_phys*S_rad_yl + gupyz_phys*S_rad_zl))/(alpn1*Psi6*u0(i,j,k))-(E_rad(&
  &i,j,k)+P_radl)*u0(i,j,k)*(vy(i,j,k)+shifty(i,j,k)) - 2.0*(F_rad0l * shifty(i,j,k)) - F_rad0l*vy(i,j,k)
                 F_radz(i,j,k) = ((gupxz_phys*S_rad_xl + gupyz_phys*S_rad_yl + gupzz_phys*S_rad_zl))/(alpn1*Psi6*u0(i,j,k))-(E_rad(&
  &i,j,k)+P_radl)*u0(i,j,k)*(vz(i,j,k)+shiftz(i,j,k)) - 2.0*(F_rad0l * shiftz(i,j,k)) - F_rad0l*vz(i,j,k)
              end if
              if (abs(F_radx(i,j,k)) .gt. 1.0e2 .or. abs(F_rady(i,j,k)) .gt. 1.0e2 .or. abs(F_radz(i,j,k)) .gt. 1.0e2) then
                 E_rad(i,j,k) = Erad_atm_cut
                 F_radx(i,j,k) = 0.0
                 F_rady(i,j,k) = 0.0
                 F_radz(i,j,k) = 0.0
              end if
           end do
        end do
     end do
!     end if
     do k=1, ext(3)
        do j=1, ext(2)
           do i=1, ext(1)
              if (abs(F_radx(i,j,k)) .gt. 1.0e2 .or. abs(F_rady(i,j,k)) .gt. 1.0e2 .or. abs(F_radz(i,j,k)) .gt. 1.0e2) then
                 write(*,*) "1. In driver_post_regrid.F90, BEFORE primitive_vars_rad_pr, F_radx too large!!!!"
                 write(*,*) "i,j,k X,Y,Z = ", i,j,k, X(i,j,k), Y(i,j,k), Z(i,j,k)
                 write(*,*) "E_rad(i,j,k), F_radx(i,j,k), F_rady(i,j,k), F_radz(i,j,k), F_rad0(i,j,k)=",E_rad(i,j,k), F_radx(i,j,k)&
  &, F_rady(i,j,k), F_radz(i,j,k), F_rad0(i,j,k)
                 write(*,*) "tau_rad(i,j,k), S_rad_x(i,j,k), S_rad_y(i,j,k), S_rad_z(i,j,k)=",tau_rad(i,j,k), S_rad_x(i,j,k), S_rad&
  &_y(i,j,k), S_rad_z(i,j,k)
                 write(*,*) "rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k), u0(i,j,k)=", rho_b(i,j,k), P(i,j,k), vx(i,j,k&
  &), vy(i,j,k), vz(i,j,k), u0(i,j,k)
                 write(*,*) "Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)=", Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)
                 write(*,*) "Resettting!!!!!"
                 E_rad(i,j,k) = Erad_atm_cut
                 F_radx(i,j,k) = 0.0
                 F_rady(i,j,k) = 0.0
                 F_radz(i,j,k) = 0.0
              end if
           end do
        end do
     end do
     print *, "^^^^^^^^^^^^^^^^^^^^^^ Inside driver_post_regrid, START primitive_vars_rad_cpp (anit-electron)^^^^^^^^^^^^^^^^^^^^^^&
  &^^"
     call primitive_vars_rad_cpp(ext,cctk_nghostzones, X, Y, Z, &
                   rho_star, &
                   tau_rad,S_rad_x,S_rad_y,S_rad_z, &
                   Sx, Sy, Sz, rho, &
                   Sxx, Sxy, Sxz, Syy, Syz, Szz,&
                   E_rad,F_radx,F_rady,F_radz, F_rad0, F_rad_scalar,&
                   P_radxx, P_radyy, P_radzz, P_radxy, P_radxz, P_radyz, &
                   phi, lapm1, shiftx, shifty, shiftz, gxx, gxy, gxz, gyy, gyz, gzz, &
                   gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
                   vx, vy, vz, u0, chi_rad, zeta_rad,&
                   cctkGH,temp6,0, repairs_rad_needed, Psi6threshold, Erad_atm_cut)
              do k=1, ext(3)
                 do j=1, ext(2)
                    do i=1, ext(1)
                       if (abs(F_radx(i,j,k)) .gt. 1.0e2 .or. abs(F_rady(i,j,k)) .gt. 1.0e2 .or. abs(F_radz(i,j,k)) .gt. 1.0e2) then
                          write(*,*) " In driver_post_regrid.F90, AFTER primitive_vars_rad_pr, F_radx too large!!!!"
                          write(*,*) "i,j,k X,Y,Z = ", i,j,k, X(i,j,k), Y(i,j,k), Z(i,j,k)
                          write(*,*) "E_rad(i,j,k), F_radx(i,j,k), F_rady(i,j,k), F_radz(i,j,k), F_rad0(i,j,k)=",E_rad(i,j,k), F_ra&
  &dx(i,j,k), F_rady(i,j,k), F_radz(i,j,k), F_rad0(i,j,k)
                          write(*,*) "tau_rad(i,j,k), S_rad_x(i,j,k), S_rad_y(i,j,k), S_rad_z(i,j,k)=",tau_rad(i,j,k), S_rad_x(i,j,&
  &k), S_rad_y(i,j,k), S_rad_z(i,j,k)
                          write(*,*) "rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k), u0(i,j,k)=", rho_b(i,j,k), P(i,j,k),&
  & vx(i,j,k), vy(i,j,k), vz(i,j,k), u0(i,j,k)
                          write(*,*) "Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)=", Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i&
  &,j,k)
                       end if
                    end do
                 end do
              end do
              print *, "^^^^^^^^^^^^^^^^^^^^^^ Inside driver_post_regrid, END primitive_vars_rad_cpp^^^^^^^^^^^^^^^^^^^^^^^^"
           ! Here we only use this function to repair radiation!!!                                                                  
           if (repairs_rad_needed==1) then
              call repair_failures_mhd_hybrid(ext,Z,gamma_th, &
                   temp10, temp6, rho_b, P, &
                   vx, vy, vz, u0, w, h, rho_star, tau, &
                   st_x, st_y, st_z, mhd_st_x, mhd_st_y, mhd_st_z, &
                   rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
                   E_rad,F_radx,F_rady,F_radz,F_rad0, F_rad_scalar,&
                   P_radxx, P_radyy, P_radzz, P_radxy, P_radxz, P_radyz,&
                   tau_rad, S_rad_x, S_rad_y, S_rad_z, &
                   Erad_atm_cut, rad_closure_scheme, enable_OS_collapse, &
                   lapm1, shiftx, shifty, shiftz, phi, &
                   gxx, gxy, gxz, gyy, gyz, gzz, &
                   gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Bx, By, Bz, &
                   sbt, sbx, sby, sbz, rho_b_atm, &
                   neos,ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
                   proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
                   glob_imax,glob_jmax,glob_kmax,Symmetry, primitives_solver)
           end if
           if (microphysics_scheme== 1) then
              !!! Now let evolve other neutrino specices!!!
!              call primitive_rad_guess(ext, tau_rad_nue, S_rad_x_nue, S_rad_y_nue, S_rad_z_nue, phi, lapm1, gupxx, gupyy, gupzz, gu
     do k=1,ext(3)
        do j=1,ext(2)
           do i=1,ext(1)
              ! Here we provide approximations to Erad and Frad using Eddington factor = 2/3 (as the midpoint between thick and thin
              tau_radl = tau_rad_nue(i,j,k)
              S_rad_xl = S_rad_x_nue(i,j,k)
              S_rad_yl = S_rad_y_nue(i,j,k)
              S_rad_zl = S_rad_z_nue(i,j,k)
              Psi4 = exp(4.d0*phi(i,j,k))
              Psi6 = exp(6.d0*phi(i,j,k))
              alpn1 = 1.d0 + lapm1(i,j,k)
              gupxx_phys = gupxx(i,j,k)/Psi4
              gupyy_phys = gupyy(i,j,k)/Psi4
              gupzz_phys = gupzz(i,j,k)/Psi4
              gupxy_phys = gupxy(i,j,k)/Psi4
              gupxz_phys = gupxz(i,j,k)/Psi4
              gupyz_phys = gupyz(i,j,k)/Psi4
              if (tau_radl < 0.0) then
                 E_rad_nue(i,j,k) = Erad_atm_cut
                 F_radx_nue(i,j,k) = 0.0
                 F_rady_nue(i,j,k) = 0.0
                 F_radz_nue(i,j,k) = 0.0
              else
                 au0 = alpn1*u0(i,j,k)
                 F_rad0_denom = alpn1 * Psi6 * (2.0*au0*au0 + 1.0)/(4.0*au0*au0 - 1.0)
                 F_rad0_num = 4.0*au0*(1.0 - au0*au0)/(4.0*au0*au0 - 1.0)*tau_radl +( (shiftx(i,j,k)*u0(i,j,k)+u0(i,j,k)*vx(i,j,k))&
  &*S_rad_xl +(shifty(i,j,k)*u0(i,j,k)+u0(i,j,k)*vy(i,j,k))*S_rad_yl+(shiftz(i,j,k)*u0(i,j,k)+u0(i,j,k)*vz(i,j,k))*S_rad_zl)
                 F_rad0l = F_rad0_num/F_rad0_denom
                 E_rad_nue(i,j,k) = 3.0*(tau_radl/Psi6 - 2.0*alpn1*au0*F_rad0l)/(4.0*au0*au0 - 1.0)
                 P_radl = 2.0*E_rad_nue(i,j,k)/3.0
                 F_radx_nue(i,j,k) = ((gupxx_phys*S_rad_xl + gupxy_phys*S_rad_yl + gupxz_phys*S_rad_zl))/(alpn1*Psi6*u0(i,j,k))-(E_&
  &rad_nue(i,j,k)+P_radl)*u0(i,j,k)*(vx(i,j,k)+shiftx(i,j,k)) - 2.0*(F_rad0l * shiftx(i,j,k)) - F_rad0l*vx(i,j,k)
                 F_rady_nue(i,j,k) = ((gupxy_phys*S_rad_xl + gupyy_phys*S_rad_yl + gupyz_phys*S_rad_zl))/(alpn1*Psi6*u0(i,j,k))-(E_&
  &rad_nue(i,j,k)+P_radl)*u0(i,j,k)*(vy(i,j,k)+shifty(i,j,k)) - 2.0*(F_rad0l * shifty(i,j,k)) - F_rad0l*vy(i,j,k)
                 F_radz_nue(i,j,k) = ((gupxz_phys*S_rad_xl + gupyz_phys*S_rad_yl + gupzz_phys*S_rad_zl))/(alpn1*Psi6*u0(i,j,k))-(E_&
  &rad_nue(i,j,k)+P_radl)*u0(i,j,k)*(vz(i,j,k)+shiftz(i,j,k)) - 2.0*(F_rad0l * shiftz(i,j,k)) - F_rad0l*vz(i,j,k)
              end if
              if (abs(F_radx(i,j,k)) .gt. 1.0e2 .or. abs(F_rady(i,j,k)) .gt. 1.0e2 .or. abs(F_radz(i,j,k)) .gt. 1.0e2) then
                 E_rad_nue(i,j,k) = Erad_atm_cut
                 F_radx_nue(i,j,k) = 0.0
                 F_rady_nue(i,j,k) = 0.0
                 F_radz_nue(i,j,k) = 0.0
              end if
           end do
        end do
     end do
              do k=1, ext(3)
                 do j=1, ext(2)
                    do i=1, ext(1)
                       if (abs(F_radx_nue(i,j,k)) .gt. 1.0e2 .or. abs(F_rady_nue(i,j,k)) .gt. 1.0e2 .or. abs(F_radz_nue(i,j,k)) .gt&
  &. 1.0e2) then
                          write(*,*) "2. In driver_post_regrid.F90, BEFORE primitive_vars_rad_pr, F_radx too large!!!!"
                          write(*,*) "i,j,k X,Y,Z = ", i,j,k, X(i,j,k), Y(i,j,k), Z(i,j,k)
                          write(*,*) "E_rad_nue(i,j,k), F_radx_nue(i,j,k), F_rady_nue(i,j,k), F_radz_nue(i,j,k), F_rad0_nue(i,j,k)=&
  &",E_rad_nue(i,j,k), F_radx_nue(i,j,k), F_rady_nue(i,j,k), F_radz_nue(i,j,k), F_rad0_nue(i,j,k)
                          write(*,*) "tau_rad_nue(i,j,k), S_rad_x_nue(i,j,k), S_rad_y_nue(i,j,k), S_rad_z_nue(i,j,k)=",tau_rad_nue(&
  &i,j,k), S_rad_x_nue(i,j,k), S_rad_y_nue(i,j,k), S_rad_z_nue(i,j,k)
                          write(*,*) "rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k), u0(i,j,k)=", rho_b(i,j,k), P(i,j,k),&
  & vx(i,j,k), vy(i,j,k), vz(i,j,k), u0(i,j,k)
                          write(*,*) "Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)=", Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i&
  &,j,k)
                          write(*,*) "Resettting!!!!!"
                          E_rad_nue(i,j,k) = Erad_atm_cut
                          F_radx_nue(i,j,k) = 0.0
                          F_rady_nue(i,j,k) = 0.0
                          F_radz_nue(i,j,k) = 0.0
                       end if
                    end do
                 end do
              end do
              print *, "^^^^^^^^^^^^^^^^^^^^^^ Inside driver_post_regrid, START primitive_vars_rad_cpp for NUE^^^^^^^^^^^^^^^^^^^^^&
  &^^^"
              call  primitive_vars_rad_cpp(ext,cctk_nghostzones, X, Y, Z, &
                   rho_star, &
                   tau_rad_nue,S_rad_x_nue,S_rad_y_nue,S_rad_z_nue, &
                   Sx, Sy, Sz, rho, &
                   Sxx, Sxy, Sxz, Syy, Syz, Szz,&
                   E_rad_nue,F_radx_nue,F_rady_nue,F_radz_nue, F_rad0_nue, F_rad_scalar_nue,&
                   P_radxx_nue, P_radyy_nue, P_radzz_nue, P_radxy_nue, P_radxz_nue, P_radyz_nue,&
                   phi, lapm1, shiftx, shifty, shiftz, gxx, gxy, gxz, gyy, gyz, gzz, &
                   gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
                   vx, vy, vz, u0, chi_rad_nue, zeta_rad_nue,&
                   cctkGH, temp6, 0, repairs_rad_needed, Psi6threshold, Erad_atm_cut)
              ! Here we only use this function to repair radiation!!!                                                               
              if (repairs_rad_needed==1) then
                 call repair_failures_mhd_hybrid(ext,Z,gamma_th, &
                      temp10, temp6, rho_b, P, &
                      vx, vy, vz, u0, w, h, rho_star, tau, &
                      st_x, st_y, st_z, mhd_st_x, mhd_st_y, mhd_st_z, &
                      rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
                      E_rad_nue,F_radx_nue,F_rady_nue,F_radz_nue,F_rad0_nue, F_rad_scalar_nue,&
                      P_radxx_nue, P_radyy_nue, P_radzz_nue, P_radxy_nue, P_radxz_nue, P_radyz_nue,&
                      tau_rad_nue, S_rad_x_nue, S_rad_y_nue, S_rad_z_nue, &
                      Erad_atm_cut, rad_closure_scheme, enable_OS_collapse, &
                      lapm1, shiftx, shifty, shiftz, phi, &
                      gxx, gxy, gxz, gyy, gyz, gzz, &
                      gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Bx, By, Bz, &
                      sbt, sbx, sby, sbz, rho_b_atm, &
                      neos,ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
                      proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
                      glob_imax,glob_jmax,glob_kmax,Symmetry, primitives_solver)
              end if
!              call primitive_rad_guess(ext, tau_rad_nux, S_rad_x_nux, S_rad_y_nux, S_rad_z_nux, phi, lapm1, gupxx, gupyy, gupzz, gu
     do k=1,ext(3)
        do j=1,ext(2)
           do i=1,ext(1)
              ! Here we provide approximations to Erad and Frad using Eddington factor = 2/3 (as the midpoint between thick and thin
              tau_radl = tau_rad_nux(i,j,k)
              S_rad_xl = S_rad_x_nux(i,j,k)
              S_rad_yl = S_rad_y_nux(i,j,k)
              S_rad_zl = S_rad_z_nux(i,j,k)
              Psi4 = exp(4.d0*phi(i,j,k))
              Psi6 = exp(6.d0*phi(i,j,k))
              alpn1 = 1.d0 + lapm1(i,j,k)
              gupxx_phys = gupxx(i,j,k)/Psi4
              gupyy_phys = gupyy(i,j,k)/Psi4
              gupzz_phys = gupzz(i,j,k)/Psi4
              gupxy_phys = gupxy(i,j,k)/Psi4
              gupxz_phys = gupxz(i,j,k)/Psi4
              gupyz_phys = gupyz(i,j,k)/Psi4
              if (tau_radl < 0.0) then
                 E_rad_nux(i,j,k) = Erad_atm_cut
                 F_radx_nux(i,j,k) = 0.0
                 F_rady_nux(i,j,k) = 0.0
                 F_radz_nux(i,j,k) = 0.0
              else
                 au0 = alpn1*u0(i,j,k)
                 F_rad0_denom = alpn1 * Psi6 * (2.0*au0*au0 + 1.0)/(4.0*au0*au0 - 1.0)
                 F_rad0_num = 4.0*au0*(1.0 - au0*au0)/(4.0*au0*au0 - 1.0)*tau_radl +( (shiftx(i,j,k)*u0(i,j,k)+u0(i,j,k)*vx(i,j,k))&
  &*S_rad_xl +(shifty(i,j,k)*u0(i,j,k)+u0(i,j,k)*vy(i,j,k))*S_rad_yl+(shiftz(i,j,k)*u0(i,j,k)+u0(i,j,k)*vz(i,j,k))*S_rad_zl)
                 F_rad0l = F_rad0_num/F_rad0_denom
                 E_rad_nux(i,j,k) = 3.0*(tau_radl/Psi6 - 2.0*alpn1*au0*F_rad0l)/(4.0*au0*au0 - 1.0)
                 P_radl = 2.0*E_rad_nux(i,j,k)/3.0
                 F_radx_nux(i,j,k) = ((gupxx_phys*S_rad_xl + gupxy_phys*S_rad_yl + gupxz_phys*S_rad_zl))/(alpn1*Psi6*u0(i,j,k))-(E_&
  &rad_nux(i,j,k)+P_radl)*u0(i,j,k)*(vx(i,j,k)+shiftx(i,j,k)) - 2.0*(F_rad0l * shiftx(i,j,k)) - F_rad0l*vx(i,j,k)
                 F_rady_nux(i,j,k) = ((gupxy_phys*S_rad_xl + gupyy_phys*S_rad_yl + gupyz_phys*S_rad_zl))/(alpn1*Psi6*u0(i,j,k))-(E_&
  &rad_nux(i,j,k)+P_radl)*u0(i,j,k)*(vy(i,j,k)+shifty(i,j,k)) - 2.0*(F_rad0l * shifty(i,j,k)) - F_rad0l*vy(i,j,k)
                 F_radz_nux(i,j,k) = ((gupxz_phys*S_rad_xl + gupyz_phys*S_rad_yl + gupzz_phys*S_rad_zl))/(alpn1*Psi6*u0(i,j,k))-(E_&
  &rad_nux(i,j,k)+P_radl)*u0(i,j,k)*(vz(i,j,k)+shiftz(i,j,k)) - 2.0*(F_rad0l * shiftz(i,j,k)) - F_rad0l*vz(i,j,k)
              end if
              if (abs(F_radx(i,j,k)) .gt. 1.0e2 .or. abs(F_rady(i,j,k)) .gt. 1.0e2 .or. abs(F_radz(i,j,k)) .gt. 1.0e2) then
                 E_rad_nux(i,j,k) = Erad_atm_cut
                 F_radx_nux(i,j,k) = 0.0
                 F_rady_nux(i,j,k) = 0.0
                 F_radz_nux(i,j,k) = 0.0
              end if
           end do
        end do
     end do
              do k=1, ext(3)
                 do j=1, ext(2)
                    do i=1, ext(1)
                       if (abs(F_radx_nux(i,j,k)) .gt. 1.0e2 .or. abs(F_rady_nux(i,j,k)) .gt. 1.0e2 .or. abs(F_radz_nux(i,j,k)) .gt&
  &. 1.0e2) then
                          write(*,*) "3. In driver_post_regrid.F90, BEFORE primitive_vars_rad_pr, F_radx too large!!!!"
                          write(*,*) "i,j,k X,Y,Z = ", i,j,k, X(i,j,k), Y(i,j,k), Z(i,j,k)
                          write(*,*) "E_rad_nux(i,j,k), F_radx_nux(i,j,k), F_rady_nux(i,j,k), F_radz_nux(i,j,k), F_rad0_nux(i,j,k)=&
  &",E_rad_nux(i,j,k), F_radx_nux(i,j,k), F_rady_nux(i,j,k), F_radz_nux(i,j,k), F_rad0_nux(i,j,k)
                          write(*,*) "tau_rad_nux(i,j,k), S_rad_x_nux(i,j,k), S_rad_y_nux(i,j,k), S_rad_z_nux(i,j,k)=",tau_rad_nux(&
  &i,j,k), S_rad_x_nux(i,j,k), S_rad_y_nux(i,j,k), S_rad_z_nux(i,j,k)
                          write(*,*) "rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k), u0(i,j,k)=", rho_b(i,j,k), P(i,j,k),&
  & vx(i,j,k), vy(i,j,k), vz(i,j,k), u0(i,j,k)
                          write(*,*) "Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)=", Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i&
  &,j,k)
                          write(*,*) "Resettting!!!!!"
                          E_rad_nux(i,j,k) = Erad_atm_cut
                          F_radx_nux(i,j,k) = 0.0
                          F_rady_nux(i,j,k) = 0.0
                          F_radz_nux(i,j,k) = 0.0
                       end if
                    end do
                 end do
              end do
              print *, "^^^^^^^^^^^^^^^^^^^^^^ Inside driver_post_regrid, START primitive_vars_rad_cpp for NUX^^^^^^^^^^^^^^^^^^^^^&
  &^^^"
              call  primitive_vars_rad_cpp(ext,cctk_nghostzones, X, Y, Z, &
                   rho_star, &
                   tau_rad_nux,S_rad_x_nux,S_rad_y_nux,S_rad_z_nux, &
                   Sx, Sy, Sz, rho, &
                   Sxx, Sxy, Sxz, Syy, Syz, Szz,&
                   E_rad_nux,F_radx_nux,F_rady_nux,F_radz_nux, F_rad0_nux, F_rad_scalar_nux,&
                   P_radxx_nux, P_radyy_nux, P_radzz_nux, P_radxy_nux, P_radxz_nux, P_radyz_nux,&
                   phi, lapm1, shiftx, shifty, shiftz, gxx, gxy, gxz, gyy, gyz, gzz, &
                   gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
                   vx, vy, vz, u0, chi_rad_nux, zeta_rad_nux,&
                   cctkGH, temp6, 0, repairs_rad_needed, Psi6threshold, Erad_atm_cut)
              ! Here we only use this function to repair radiation!!!                                                               
              if (repairs_rad_needed==1) then
                 call repair_failures_mhd_hybrid(ext,Z,gamma_th, &
                      temp10, temp6, rho_b, P, &
                      vx, vy, vz, u0, w, h, rho_star, tau, &
                      st_x, st_y, st_z, mhd_st_x, mhd_st_y, mhd_st_z, &
                      rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
                      E_rad_nux,F_radx_nux,F_rady_nux,F_radz_nux,F_rad0_nux, F_rad_scalar_nux,&
                      P_radxx_nux, P_radyy_nux, P_radzz_nux, P_radxy_nux, P_radxz_nux, P_radyz_nux,&
                      tau_rad_nux, S_rad_x_nux, S_rad_y_nux, S_rad_z_nux, &
                      Erad_atm_cut, rad_closure_scheme, enable_OS_collapse, &
                      lapm1, shiftx, shifty, shiftz, phi, &
                      gxx, gxy, gxz, gyy, gyz, gzz, &
                      gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Bx, By, Bz, &
                      sbt, sbx, sby, sbz, rho_b_atm, &
                      neos,ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
                      proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
                      glob_imax,glob_jmax,glob_kmax,Symmetry, primitives_solver)
              end if
           end if
        end if
        else if(primitives_solver==1) then
           !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
                    !Note that the correct eqn for w is: 
                    ! w=rho_star*(lapm1+1.0)*u0
                    ! We must fill in w everywhere here since, for primitive_vars_hybrid2(), we need the old value of w
                    !    everywhere on the grid.  This is a problem since w is not prolongated and hence is not set
                    !    everywhere after regridding.  The below line just provides an initial guess for primitive_vars_hybrid2()
                    !Here, temp2 is supposed to be w_p
                    temp2(i,j,k) = rho_star(i,j,k)
                 end do
              end do
           end do
           !$omp end parallel do
           !write(*,*) PHI=,phi(79,80,3),ext
           !  Zach says: Warning: if ignore_ghostzones=1, we get undefined values post-regridding
           call hydro_primitives(cctkGH,ext,cctk_nghostzones, X, Y, Z, rho_star, tau,&
                mhd_st_x, mhd_st_y, mhd_st_z,&
                u0,vx,vy,vz,&
                w, temp2, rho_b,rho, P, h, &
                Sx, Sy, Sz, &
                Sxx, Sxy, Sxz, Syy, Syz, Szz, &
                phi, lapm1, shiftx,shifty,shiftz, &
                gxx, gxy, gxz, gyy, gyz, gzz, &
                gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,&
                tau_stildefix_enable,tau_atm,enable_shocktest_primitive_mode,&
                rho_b_max, rho_fail_max_step, M_fail_step, rho_b_atm,&
                gamma_th,K_poly,sdots_o_rhot,Symmetry,0)
           !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
                    ! Following lines are necessary
                    st_x(i,j,k) = mhd_st_x(i,j,k)
                    st_y(i,j,k) = mhd_st_y(i,j,k)
                    st_z(i,j,k) = mhd_st_z(i,j,k)
                    !MAYBE FIXME: Put in a warning if B!=0.  However, if B!=0, you shouldnt be using this primitives solver anyway!
                    sbt(i,j,k) = 0.D0
                    sbx(i,j,k) = 0.D0
                    sby(i,j,k) = 0.D0
                    sbz(i,j,k) = 0.D0
                    Bx(i,j,k) = 0.D0
                    By(i,j,k) = 0.D0
                    Bz(i,j,k) = 0.D0
                    Ex(i,j,k) = 0.D0
                    Ey(i,j,k) = 0.D0
                    Ez(i,j,k) = 0.D0
                 end do
              end do
           end do
           !$omp end parallel do
        else if(primitives_solver==2) then
           !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
                    ! Here we provide approximations to rho_b and P, as first guesses for the primitives solver
                    rho_b(i,j,k) = max(rho_star(i,j,k) * exp(-6.D0*phi(i,j,k)), rho_b_atm)
                    call compute_pcold_epscold(rho_b(i,j,k),P(i,j,k),temp1(i,j,k), &
                         neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab,enable_OS_collapse)
                    ! Finally, an approximation of h_p, so that eps = 1.D0:
                    !   (Recall that eps = h_old[index]-1.0-P_old[index]/rho_b_old)
                    h_p(i,j,k) = 1.d0 + P(i,j,k)/rho_b(i,j,k) + temp1(i,j,k)
                    vx(i,j,k) = -shiftx(i,j,k)
                    vy(i,j,k) = -shifty(i,j,k)
                    vz(i,j,k) = -shiftz(i,j,k)
                    u0(i,j,k) = 1.D0/(lapm1(i,j,k)+1.d0)
                    if(isnan(tau_rad(i,j,k)).or.isnan(S_rad_x(i,j,k))) then
                       write(*,*) "In driver_primitives_postmolstep, before primitive_vars_hybrid2_cpp, tau_rad and S_rad_x=", tau_&
  &rad, S_rad_x
                       write(*,*) "X,Y,Z=", X(i,j,k), Y(i,j,k), Z(i,j,k)
                    end if
                 end do
              end do
           end do
           !$omp end parallel do
           !TODO: remove temp1 (unused) parameter.
           !  Zach says: Warning: if ignore_ghostzones=1, we get undefined values post-regridding
!           print *, ^^^^^^^^^^^^^^^^^^^^^^ Inside driver_post_regrid.F90 start primitive_vars_hybrid2_cpp ^^^^^^^^^^^^^^^^^^^^^^^^
           call primitive_vars_hybrid2_cpp(ext,cctk_nghostzones,X,Y,Z, &
                rho_star,tau,st_x,st_y,st_z, &
                mhd_st_x,mhd_st_y,mhd_st_z, rhoYe, &
                tau_rad,S_rad_x,S_rad_y,S_rad_z,&
                neos, ergo_star, ergo_sigma,&
                rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, &
                w,temp1,rho_b,rho,P,h,Sx,Sy,Sz, &
                Sxx,Sxy,Sxz,Syy,Syz,Szz, &
                E_rad,F_radx,F_rady,F_radz, F_rad0, F_rad_scalar, Y_e, eps_thermal,&
                phi,lapm1,shiftx,shifty,shiftz, &
                gxx,gxy,gxz,gyy,gyz,gzz, &
                gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                h_p,u0,rho_b_max,rho_b_atm, &
                rho_fail_max_step,M_fail_step,rhos_max, &
                Bx,By,Bz,Ex,Ey,Ez, &
                vx,vy,vz, &
                sbt,sbx,sby,sbz, &
                xNS1,yNS1, xNS2, yNS2, M_B, rad_T_fac,rad_T_cutoff, rad_T_pow, rad_T_floor, T_fluid,&
                proc_imin,proc_jmin,proc_kmin, &
                proc_imax,proc_jmax,proc_kmax, &
                glob_imax,glob_jmax,glob_kmax, &
                Symmetry,pfloor,excision_enable, &
                excision_zone_gf, tau_stildefix_enable,tau_atm,temp4,cctkGH,0, &
                enable_shocktest_primitive_mode,repairs_needed,rad_closure_scheme, &
                rad_const, enable_OS_collapse, compute_microphysics, Psi6threshold, Erad_atm_cut)
           if(rad_evolve_enable==1 .and. rad_closure_scheme == 1) then
!              call CartSymGN(dummy,cctkGH,mhd_evolve::rad_pressure)
              call primitive_vars_rad_cpp(ext,cctk_nghostzones, X, Y, Z, &
                   rho_star, &
                   tau_rad,S_rad_x,S_rad_y,S_rad_z, &
                   Sx, Sy, Sz, rho, &
                   Sxx, Sxy, Sxz, Syy, Syz, Szz,&
                   E_rad,F_radx,F_rady,F_radz, F_rad0, F_rad_scalar,&
                   P_radxx, P_radyy, P_radzz, P_radxy, P_radxz, P_radyz, &
                   phi, lapm1, shiftx, shifty, shiftz, gxx, gxy, gxz, gyy, gyz, gzz, &
                   gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
                   vx, vy, vz, u0, chi_rad, zeta_rad,&
                   cctkGH, temp5, 0, repairs_rad_needed, Psi6threshold, Erad_atm_cut)
              end if
           if (repairs_needed==1) then
              call repair_failures_mhd_hybrid(ext,Z,gamma_th, &
                  temp4, temp5, rho_b, P, &
                  vx, vy, vz, u0, w, h, rho_star, tau, &
                  st_x, st_y, st_z, mhd_st_x, mhd_st_y, mhd_st_z, &
                  rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
                  E_rad,F_radx,F_rady,F_radz,F_rad0, F_rad_scalar,&
                  P_radxx, P_radyy, P_radzz, P_radxy, P_radxz, P_radyz,&
                  tau_rad, S_rad_x, S_rad_y, S_rad_z, &
                  Erad_atm_cut, rad_closure_scheme, enable_OS_collapse, &
                  lapm1, shiftx, shifty, shiftz, phi, &
                  gxx, gxy, gxz, gyy, gyz, gzz, &
                  gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Bx, By, Bz, &
                  sbt, sbx, sby, sbz, rho_b_atm, &
                  neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
                  proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
                  glob_imax,glob_jmax,glob_kmax,Symmetry, primitives_solver)
           end if
           if(excision_enable==1) then
             call hydro_ezbc_hybrid(ext,X,Y,Z,rho_star,tau, &
                  mhd_st_x,mhd_st_y,mhd_st_z,st_x,st_y,st_z, &
                  rho_b,P,h,vx,vy,vz,w,&
                  sbt,sbx,sby,sbz,Bx,By,Bz, &
                  lapm1,shiftx,shifty,shiftz,phi,&
                  gxx,gxy,gxz,gyy,gyz,gzz, &
                  gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,&
                  Symmetry,excision_zone_gf,gamma_th,neos,ergo_star, ergo_sigma, rho_tab,&
                  P_tab,eps_tab,k_tab,gamma_tab)
             if (hyperbolic_divergence_cleaning_enable==1) then
                call scalar_excision_bc(ext,X,Y,Z,Blagrangemultiplier,Symmetry,excision_zone_gf)
             end if
             call remove_interior2(ext,X,Y,Z,sbt,excision_zone_gf,Symmetry)
             call remove_interior2(ext,X,Y,Z,sbx,excision_zone_gf,Symmetry)
             call remove_interior2(ext,X,Y,Z,sby,excision_zone_gf,Symmetry)
             call remove_interior2(ext,X,Y,Z,sbz,excision_zone_gf,Symmetry)
           end if
        else if(primitives_solver==0) then
           !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
                    !if(i.le.cctk_nghostzones(1) .or. j.le.cctk_nghostzones(2) .or. k.le.cctk_nghostzones(3)) then
                    ! Here we provide approximations to rho_b and P, as first guesses for the primitives solver
                    rho_b(i,j,k) = max(rho_star(i,j,k) * exp(-6.D0*phi(i,j,k)), rho_b_atm)
                    call compute_pcold_epscold(rho_b(i,j,k),P(i,j,k),temp1(i,j,k), &
                                                   neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,enable_OS_collap&
  &se)
                    ! Finally, an approximation of h_p, so that eps = 1.D0:
                    !   (Recall that eps = h_old[index]-1.0-P_old[index]/rho_b_old)
                    h_p(i,j,k) = 1.d0 + P(i,j,k)/rho_b(i,j,k) + temp1(i,j,k)
                    vx(i,j,k) = 0.D0
                    vy(i,j,k) = 0.D0
                    vz(i,j,k) = 0.D0
                 end do
              end do
           end do
            !$omp end parallel do
           !TODO: remove temp1 (unused) parameter.
           !  Zach says: Warning: if ignore_ghostzones=1, we get undefined values post-regridding
           call primitive_vars_hybrid2(ext,cctk_nghostzones,X,Y,Z, &
                rho_star,tau,st_x,st_y,st_z, &
                mhd_st_x,mhd_st_y,mhd_st_z,neos, ergo_star, ergo_sigma,&
                rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, &
                w,temp1,rho_b,rho,P,h,Sx,Sy,Sz, &
                Sxx,Sxy,Sxz,Syy,Syz,Szz, &
                phi,lapm1,shiftx,shifty,shiftz, &
                gxx,gxy,gxz,gyy,gyz,gzz, &
                gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                h_p,u0,rho_b_max,rho_b_atm, &
                rho_fail_max_step,M_fail_step,rhos_max, &
                Bx,By,Bz,Ex,Ey,Ez, &
                vx,vy,vz, &
                sbt,sbx,sby,sbz, temp4, &
                proc_imin,proc_jmin,proc_kmin, &
                proc_imax,proc_jmax,proc_kmax, &
                glob_imax,glob_jmax,glob_kmax, &
                Symmetry,pfloor,excision_enable, &
                excision_zone_gf, tau_stildefix_enable,tau_atm,0)
        end if
     end if
     write(*,*) "HELLO postregrid:tau,dX=",tau(1,1,1),dX
  end if
end subroutine mhd_post_regrid
subroutine primitive_rad_guess(ext, tau_rad, S_rad_x, S_rad_y, S_rad_z, phi, lapm1, gupxx, gupyy, gupzz, gupxy, gupxz, gupyz, shift&
  &x, shifty, shiftz, vx, vy, vz, E_rad, F_radx, F_rady, F_radz, Erad_atm_cut)
implicit none
integer, dimension(3)                          :: ext
real*8, dimension(ext(1),ext(2),ext(3))        :: tau_rad, S_rad_x, S_rad_y, S_rad_z
real*8, dimension(ext(1),ext(2),ext(3))        :: phi, lapm1, u0, gupxx, gupyy, gupzz, gupxy, gupxz, gupyz
real*8, dimension(ext(1),ext(2),ext(3))        :: shiftx, shifty, shiftz
real*8, dimension(ext(1),ext(2),ext(3))        :: vx, vy, vz
real*8, dimension(ext(1),ext(2),ext(3))        :: E_rad, F_radx, F_rady, F_radz
real*8                                         :: tau_radl, S_rad_xl, S_rad_yl, S_rad_zl
real*8                                         :: Psi4, Psi6, alpn1
real*8                                         :: gupxx_phys, gupyy_phys, gupzz_phys, gupxy_phys, gupxz_phys, gupyz_phys
real*8                                         :: au0, F_rad0_denom, F_rad0_num, F_rad0l, P_radl
real*8                                         :: Erad_atm_cut
integer                                        :: i,j,k
do k=1,ext(3)
   do j=1,ext(2)
      do i=1,ext(1)
         ! Here we provide approximations to Erad and Frad using Eddington factor = 2/3 (as the midpoint between thick and thin limi
         tau_radl = tau_rad(i,j,k)
         S_rad_xl = S_rad_x(i,j,k)
         S_rad_yl = S_rad_y(i,j,k)
         S_rad_zl = S_rad_z(i,j,k)
         Psi4 = exp(4.d0*phi(i,j,k))
         Psi6 = exp(6.d0*phi(i,j,k))
         alpn1 = 1.d0 + lapm1(i,j,k)
         gupxx_phys = gupxx(i,j,k)/Psi4
         gupyy_phys = gupyy(i,j,k)/Psi4
         gupzz_phys = gupzz(i,j,k)/Psi4
         gupxy_phys = gupxy(i,j,k)/Psi4
         gupxz_phys = gupxz(i,j,k)/Psi4
         gupyz_phys = gupyz(i,j,k)/Psi4
         if (tau_radl < 0.0) then
            E_rad(i,j,k) = Erad_atm_cut
            F_radx(i,j,k) = 0.0
            F_rady(i,j,k) = 0.0
            F_radz(i,j,k) = 0.0
         else
            au0 = alpn1*u0(i,j,k)
            F_rad0_denom = alpn1 * Psi6 * (2.0*au0*au0 + 1.0)/(4.0*au0*au0 - 1.0)
            F_rad0_num = 4.0*au0*(1.0 - au0*au0)/(4.0*au0*au0 - 1.0)*tau_radl +( (shiftx(i,j,k)*u0(i,j,k)+u0(i,j,k)*vx(i,j,k))*S_ra&
  &d_xl +(shifty(i,j,k)*u0(i,j,k)+u0(i,j,k)*vy(i,j,k))*S_rad_yl+(shiftz(i,j,k)*u0(i,j,k)+u0(i,j,k)*vz(i,j,k))*S_rad_zl)
            F_rad0l = F_rad0_num/F_rad0_denom
            E_rad(i,j,k) = 3.0*(tau_radl/Psi6 - 2.0*alpn1*au0*F_rad0l)/(4.0*au0*au0 - 1.0)
            P_radl = 2.0*E_rad(i,j,k)/3.0
            F_radx(i,j,k) = ((gupxx_phys*S_rad_xl + gupxy_phys*S_rad_yl + gupxz_phys*S_rad_zl))/(alpn1*Psi6*u0(i,j,k))-(E_rad(i,j,k&
  &)+P_radl)*u0(i,j,k)*(vx(i,j,k)+shiftx(i,j,k)) - 2.0*(F_rad0l * shiftx(i,j,k)) - F_rad0l*vx(i,j,k)
            F_rady(i,j,k) = ((gupxy_phys*S_rad_xl + gupyy_phys*S_rad_yl + gupyz_phys*S_rad_zl))/(alpn1*Psi6*u0(i,j,k))-(E_rad(i,j,k&
  &)+P_radl)*u0(i,j,k)*(vy(i,j,k)+shifty(i,j,k)) - 2.0*(F_rad0l * shifty(i,j,k)) - F_rad0l*vy(i,j,k)
            F_radz(i,j,k) = ((gupxz_phys*S_rad_xl + gupyz_phys*S_rad_yl + gupzz_phys*S_rad_zl))/(alpn1*Psi6*u0(i,j,k))-(E_rad(i,j,k&
  &)+P_radl)*u0(i,j,k)*(vz(i,j,k)+shiftz(i,j,k)) - 2.0*(F_rad0l * shiftz(i,j,k)) - F_rad0l*vz(i,j,k)
         end if
         if (abs(F_radx(i,j,k)) .gt. 1.0e2 .or. abs(F_rady(i,j,k)) .gt. 1.0e2 .or. abs(F_radz(i,j,k)) .gt. 1.0e2) then
!            write(*,*)  In driver_post_regrid.F90, BEFORE primitive_vars_rad, F_radx too large!!!!
!            write(*,*) i,j,k X,Y,Z = , i,j,k, X(i,j,k), Y(i,j,k), Z(i,j,k)
!            write(*,*) E_rad(i,j,k), F_radx(i,j,k), F_rady(i,j,k), F_radz(i,j,k), F_rad0(i,j,k)=,E_rad(i,j,k), F_radx(i,j,k), F_rad
!            write(*,*) tau_rad(i,j,k), S_rad_x(i,j,k), S_rad_y(i,j,k), S_rad_z(i,j,k)=,tau_rad(i,j,k), S_rad_x(i,j,k), S_rad_y(i,j,
!            write(*,*) rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k), u0(i,j,k)=, rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i
!            write(*,*) Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)=, Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)
!            write(*,*) Resettting!!!!!
            E_rad(i,j,k) = Erad_atm_cut
            F_radx(i,j,k) = 0.0
            F_rady(i,j,k) = 0.0
            F_radz(i,j,k) = 0.0
         end if
      end do
   end do
end do
end subroutine primitive_rad_guess
