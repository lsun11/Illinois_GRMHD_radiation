!-------------------------------------------------------------------
! Here we set up all required variables for shock tests, including:
!  emfields, BSSN variables, primitives, radiations, etc.
!-------------------------------------------------------------------
subroutine mhd_shock_initialdata_local(cctk_dim,cctk_gsh,cctk_lsh,cctk_lbnd,cctk_ubnd,cctk_lssh,cctk_from,cctk_to,cctk_bbox,cctk_de&
  &lta_time,cctk_time,cctk_delta_space,cctk_origin_space,cctk_levfac,cctk_levoff,cctk_levoffdenom,cctk_timefac,cctk_convlevel,cctk_&
  &convfac,cctk_nghostzones,cctk_iteration,cctkGH, X0mhd_shock_private,X1mhd_shock_private,X2mhd_shock_private,shock_rest_mass_inte&
  &grand,shock_restmass_VolInt, X0BSSN_AH,X0BSSN_aux_private,X0BSSN_aux_restrict2,X0BSSN_diag_restrict,X0BSSN_gupij,X0BSSN_matter,X&
  &0BSSN_refbd,X0BSSN_rhs,X0BSSN_vars,X0NPphi2_em,X0OmegF,X0Psi4_gw,X0Stagger_Bs,X0ah_radius,X0bh_posns,X0diag_vacuum_gf,X0disk_atm&
  &osphere,X0em_Ax,X0em_Ay,X0em_Az,X0em_Blagrangemultiplier,X0em_Blagrangemultiplier_rhs,X0em_Phi,X0em_Phi_rhs,X0em_conservativex,X&
  &0em_conservativey,X0em_conservativez,X0em_rhsx,X0em_rhsy,X0em_rhsz,X0eos_params1,X0eos_params2,X0excision_int_gfs,X0field_line_v&
  &ariables,X0field_line_variables_rhs,X0fisheye_vars,X0gw_moment_arrays,X0have_global_bdry,X0lapse_derivatives,X0lapse_rhs,X0lapse&
  &_vars,X0lapse_vars_aux,X0metric_spher_pol_1,X0metric_spher_pol_2,X0metric_spher_pol_3,X0mhd_conservatives,X0mhd_nosync,X0mhd_pri&
  &mitives,X0mhd_rhs,X0mhd_sync_lr_B_quantities,X0mhd_sync_lr_hydro_quantities,X0mhd_sync_metric_facevals,X0mhd_sync_nablas,X0mhd_s&
  &ync_nablas_diag,X0mhd_sync_rho_br_rho_bl,X0mhd_temps,X0mhd_vs,X0micphys_conservatives,X0micphys_conservatives_rhs,X0micphys_sync&
  &_lr_hydro_quantities,X0micphys_sync_nablas,X0microphys_primitives,X0phi_derivs,X0rad_conservatives,X0rad_conservatives_rhs,X0rad&
  &_pressure,X0rad_primitives,X0rad_sync_lr,X0rhovecs,X0shift_rhs,X0shift_vars,X0shift_vars_temp,X0surf_params,X1BSSN_AH,X1BSSN_aux&
  &_private,X1BSSN_aux_restrict2,X1BSSN_diag_restrict,X1BSSN_gupij,X1BSSN_matter,X1BSSN_refbd,X1BSSN_rhs,X1BSSN_vars,X1NPphi2_em,X1&
  &OmegF,X1Psi4_gw,X1Stagger_Bs,X1ah_radius,X1diag_vacuum_gf,X1disk_atmosphere,X1em_Ax,X1em_Ay,X1em_Az,X1em_Blagrangemultiplier,X1e&
  &m_Blagrangemultiplier_rhs,X1em_Phi,X1em_Phi_rhs,X1em_conservativex,X1em_conservativey,X1em_conservativez,X1em_rhsx,X1em_rhsy,X1e&
  &m_rhsz,X1excision_int_gfs,X1field_line_variables,X1field_line_variables_rhs,X1fisheye_vars,X1have_global_bdry,X1lapse_derivative&
  &s,X1lapse_rhs,X1lapse_vars,X1lapse_vars_aux,X1mhd_conservatives,X1mhd_nosync,X1mhd_primitives,X1mhd_rhs,X1mhd_sync_lr_B_quantiti&
  &es,X1mhd_sync_lr_hydro_quantities,X1mhd_sync_metric_facevals,X1mhd_sync_nablas,X1mhd_sync_nablas_diag,X1mhd_sync_rho_br_rho_bl,X&
  &1mhd_temps,X1mhd_vs,X1micphys_conservatives,X1micphys_conservatives_rhs,X1micphys_sync_lr_hydro_quantities,X1micphys_sync_nablas&
  &,X1microphys_primitives,X1phi_derivs,X1rad_conservatives,X1rad_conservatives_rhs,X1rad_pressure,X1rad_primitives,X1rad_sync_lr,X&
  &1shift_rhs,X1shift_vars,X1shift_vars_temp,X2BSSN_AH,X2BSSN_aux_private,X2BSSN_aux_restrict2,X2BSSN_diag_restrict,X2BSSN_gupij,X2&
  &BSSN_matter,X2BSSN_refbd,X2BSSN_rhs,X2BSSN_vars,X2NPphi2_em,X2OmegF,X2Psi4_gw,X2Stagger_Bs,X2ah_radius,X2diag_vacuum_gf,X2disk_a&
  &tmosphere,X2em_Ax,X2em_Ay,X2em_Az,X2em_Blagrangemultiplier,X2em_Blagrangemultiplier_rhs,X2em_Phi,X2em_Phi_rhs,X2em_conservativex&
  &,X2em_conservativey,X2em_conservativez,X2em_rhsx,X2em_rhsy,X2em_rhsz,X2excision_int_gfs,X2field_line_variables,X2field_line_vari&
  &ables_rhs,X2fisheye_vars,X2lapse_derivatives,X2lapse_rhs,X2lapse_vars,X2lapse_vars_aux,X2mhd_conservatives,X2mhd_nosync,X2mhd_pr&
  &imitives,X2mhd_rhs,X2mhd_sync_lr_B_quantities,X2mhd_sync_lr_hydro_quantities,X2mhd_sync_metric_facevals,X2mhd_sync_nablas,X2mhd_&
  &sync_nablas_diag,X2mhd_sync_rho_br_rho_bl,X2mhd_temps,X2mhd_vs,X2micphys_conservatives,X2micphys_conservatives_rhs,X2micphys_syn&
  &c_lr_hydro_quantities,X2micphys_sync_nablas,X2microphys_primitives,X2phi_derivs,X2rad_conservatives,X2rad_conservatives_rhs,X2ra&
  &d_pressure,X2rad_primitives,X2rad_sync_lr,X2shift_rhs,X2shift_vars,X2shift_vars_temp,ah_centroid_length,ah_flags_length,ah_origi&
  &n_length,ah_radius_length,Aupxx,Aupxy,Aupxz,Aupyy,Aupyz,Aupzz,Ax,Ax_p,Ax_p_p,Axx,Axx_p,Axx_p_p,Axx_rhs,Axy,Axy_p,Axy_p_p,Axy_rhs&
  &,Axz,Axz_p,Axz_p_p,Axz_rhs,Ay,Ay_p,Ay_p_p,Ayy,Ayy_p,Ayy_p_p,Ayy_rhs,Ayz,Ayz_p,Ayz_p_p,Ayz_rhs,Az,Az_p,Az_p_p,Azz,Azz_p,Azz_p_p,A&
  &zz_rhs,Bfreq1,Bfreq1_p,Bfreq1_p_p,Bfreq2,Bfreq2_p,Bfreq2_p_p,Blagrangemultiplier,Blagrangemultiplier_p,Blagrangemultiplier_p_p,B&
  &lagrangemultiplier_rhs,Box1X_VolInt,Box1X_VolInt1,Box1X_VolInt2,Box1Y_VolInt,Box1Y_VolInt1,Box1Y_VolInt2,Box1Z_VolInt,Box1Z_VolI&
  &nt1,Box1Z_VolInt2,Box1denom_VolInt,Box1denom_VolInt1,Box1denom_VolInt2,Bx,Bx_stagger,Bxl,Bxr,Bxtilde,Bxtilde_p,Bxtilde_p_p,Bxtil&
  &de_or_Ax_rhs,By,By_stagger,Byl,Byr,Bytilde,Bytilde_p,Bytilde_p_p,Bytilde_or_Ay_rhs,Bz,Bz_stagger,Bzl,Bzr,Bztilde,Bztilde_p,Bztil&
  &de_p_p,Bztilde_or_Az_rhs,CoM_VolInt_denominator,CoMx_VolInt,CoMy_VolInt,CoMz_VolInt,E_rad,E_rad_nue,E_rad_nuel,E_rad_nuer,E_rad_&
  &nux,E_rad_nuxl,E_rad_nuxr,E_radl,E_radr,Ex,Ey,Ez,F_rad0,F_rad0_nue,F_rad0_nux,F_rad_scalar,F_rad_scalar_nue,F_rad_scalar_nux,F_r&
  &adx,F_radx_nue,F_radx_nuele,F_radx_nuer,F_radx_nux,F_radx_nuxle,F_radx_nuxr,F_radxle,F_radxr,F_rady,F_rady_nue,F_rady_nuele,F_ra&
  &dy_nuer,F_rady_nux,F_rady_nuxle,F_rady_nuxr,F_radyle,F_radyr,F_radz,F_radz_nue,F_radz_nuele,F_radz_nuer,F_radz_nux,F_radz_nuxle,&
  &F_radz_nuxr,F_radzle,F_radzr,FaFal,FaFal_nue,FaFal_nux,FaFar,FaFar_nue,FaFar_nux,Fontfix_tracker_gf,Gammax,Gammax_p,Gammax_p_p,G&
  &ammax_drive,Gammax_rhs,Gammaxxx,Gammaxxy,Gammaxxz,Gammaxyy,Gammaxyz,Gammaxzz,Gammay,Gammay_p,Gammay_p_p,Gammay_drive,Gammay_rhs,&
  &Gammayxx,Gammayxy,Gammayxz,Gammayyy,Gammayyz,Gammayzz,Gammaz,Gammaz_p,Gammaz_p_p,Gammaz_drive,Gammaz_rhs,Gammazxx,Gammazxy,Gamma&
  &zxz,Gammazyy,Gammazyz,Gammazzz,Gamx_const_VolInt,Gamy_const_VolInt,Gamz_const_VolInt,Ham_const_VolIntD,Ham_const_VolIntN,Ham_con&
  &st_excised_VolIntD,Ham_const_excised_VolIntN,Ham_const_excised_innerregion_VolIntD,Ham_const_excised_innerregion_VolIntN,J_ADM_V&
  &olInt,J_ADM_VolInt_inner,J_ADM_VolInt_inner2,Jz_constraint,Jz_constraint_inner,KK,K_pp_rot1,K_pp_rot2,K_pp_rot3,K_rp_rot1,K_rp_r&
  &ot2,K_rp_rot3,K_rr_rot1,K_rr_rot2,K_rr_rot3,K_rth_rot1,K_rth_rot2,K_rth_rot3,K_thp_rot1,K_thp_rot2,K_thp_rot3,K_thth_rot1,K_thth&
  &_rot2,K_thth_rot3,M0_AH_VolInt,M0_VolInt,M0_escape100M,M0_escape30M,M0_escape50M,M0_escape70M,M0_horiz_VolInt,M0_r1_VolInt,M0_r2&
  &_VolInt,M0_r3_VolInt,MNorm,MONOPOLE,MRsx,MRsy,MRsz,M_ADM_VolInt,M_constraint,M_constraint_inner,Minternal_VolInt,Minternal_cold_&
  &VolInt,NPphi0im,NPphi0im_p,NPphi0im_p_p,NPphi0re,NPphi0re_p,NPphi0re_p_p,NPphi2im,NPphi2im_p,NPphi2im_p_p,NPphi2re,NPphi2re_p,NP&
  &phi2re_p_p,OS_surf_rad,P,P_cld,P_constraint,P_constraint_inner,P_floor,P_radxx,P_radxx_nue,P_radxx_nux,P_radxy,P_radxy_nue,P_rad&
  &xy_nux,P_radxz,P_radxz_nue,P_radxz_nux,P_radyy,P_radyy_nue,P_radyy_nux,P_radyz,P_radyz_nue,P_radyz_nux,P_radzz,P_radzz_nue,P_rad&
  &zz_nux,P_tab,P_thermal,PhysicalRadius,Pl,Pr,PsiNorm,PsiRes,PsiTau,Pvec,RadiusDerivative,RadiusDerivative2,Rxx,Rxy,Rxz,Ryy,Ryz,Rz&
  &z,S,SPx,SPx_p,SPx_p_p,SPy,SPy_p,SPy_p_p,SPz,SPz_p,SPz_p_p,S_rad_x,S_rad_x_p,S_rad_x_p_p,S_rad_x_nue,S_rad_x_nue_p,S_rad_x_nue_p_&
  &p,S_rad_x_nue_rhs,S_rad_x_nux,S_rad_x_nux_p,S_rad_x_nux_p_p,S_rad_x_nux_rhs,S_rad_x_rhs,S_rad_y,S_rad_y_p,S_rad_y_p_p,S_rad_y_nu&
  &e,S_rad_y_nue_p,S_rad_y_nue_p_p,S_rad_y_nue_rhs,S_rad_y_nux,S_rad_y_nux_p,S_rad_y_nux_p_p,S_rad_y_nux_rhs,S_rad_y_rhs,S_rad_z,S_&
  &rad_z_p,S_rad_z_p_p,S_rad_z_nue,S_rad_z_nue_p,S_rad_z_nue_p_p,S_rad_z_nue_rhs,S_rad_z_nux,S_rad_z_nux_p,S_rad_z_nux_p_p,S_rad_z_&
  &nux_rhs,S_rad_z_rhs,Sx,Sxx,Sxy,Sxz,Sy,Syy,Syz,Sz,Szz,T_0x,T_0x_p,T_0x_p_p,T_0y,T_0y_p,T_0y_p_p,T_0z,T_0z_p,T_0z_p_p,T_VolInt,T_f&
  &luid,T_fluidle,T_fluidr,Tem0_0_VolInt,Tem0_0_outsideBH_VolInt,Tfluid0_0_VolInt,Tfluid0_0_outsideBH_VolInt,VolIntegrand,VolIntegr&
  &and2,VolIntegrand3,VolIntegrand4,Y_e,Y_el,Y_er,ah_centroid_iteration,ah_centroid_iteration_p,ah_centroid_t,ah_centroid_t_p,ah_ce&
  &ntroid_valid,ah_centroid_valid_p,ah_centroid_x,ah_centroid_x_p,ah_centroid_y,ah_centroid_y_p,ah_centroid_z,ah_centroid_z_p,ah_fo&
  &und_flag,ah_initial_find_flag,ah_origin_x,ah_origin_y,ah_origin_z,ah_radius,ah_really_initial_find_flag,ah_search_flag,b_phi_Vol&
  &Int,bh_posn_x,bh_posn_y,bh_posn_z,bh_radius_z,brem_qee_VolInt,brem_qei_VolInt,chi,chi_p,chi_p_p,chi_rad,chi_rad_nue,chi_rad_nux,&
  &chi_rhs,cmax,cmax_rad,cmax_rad_nue,cmax_rad_nux,cmin,cmin_rad,cmin_rad_nue,cmin_rad_nux,dBx_m,dBy_m,dBz_m,dP_m,dT_fluid_m,ddrbrv&
  &ec,density_modes_i1,density_modes_i2,density_modes_i3,density_modes_i4,density_modes_i5,density_modes_i6,density_modes_r0,densit&
  &y_modes_r1,density_modes_r2,density_modes_r3,density_modes_r4,density_modes_r5,density_modes_r6,div_beta,drbrvec,drhoYe_m,drhoYe&
  &_m_x,drhoYe_m_xp1,drho_b_m,drho_b_m_x,drho_b_m_xp1,dvx_m,dvx_m_x,dvx_m_xp1,dvy_m,dvz_m,em_J_between_VolInt,em_energy2_VolInt,em_&
  &energy_VolInt,em_energy_between_VolInt,em_energy_outsideBH_VolInt,em_energy_outsideradius1_VolInt,em_energy_outsideradius2_VolIn&
  &t,emission_gf,emission_gf_nue,emission_gf_nux,eps_cld,eps_tab,eps_thermal,eps_tot,eta_nue,excision_zone_gf,fluid_J_horiz_VolInt,&
  &fluid_J_r1_VolInt,fluid_J_r2_VolInt,fluid_J_r3_VolInt,fluid_energy_horiz_VolInt,fluid_energy_r1_VolInt,fluid_energy_r2_VolInt,fl&
  &uid_energy_r3_VolInt,gamma_tab,gconx,gcony,gconz,gupxx,gupxx_f,gupxy,gupxy_f,gupxz,gupxz_f,gupyy,gupyy_f,gupyz,gupyz_f,gupzz,gup&
  &zz_f,gxx,gxx_p,gxx_p_p,gxx_f,gxx_rhs,gxxx,gxxy,gxxz,gxy,gxy_p,gxy_p_p,gxy_f,gxy_rhs,gxyx,gxyy,gxyz,gxz,gxz_p,gxz_p_p,gxz_f,gxz_r&
  &hs,gxzx,gxzy,gxzz,gyy,gyy_p,gyy_p_p,gyy_f,gyy_rhs,gyyx,gyyy,gyyz,gyz,gyz_p,gyz_p_p,gyz_f,gyz_rhs,gyzx,gyzy,gyzz,gzz,gzz_p,gzz_p_&
  &p,gzz_f,gzz_rhs,gzzx,gzzy,gzzz,h,h_p,half_b2_u0_VolInt,half_b2_u0_outsideBH_VolInt,have_global_bdry_max,have_global_bdry_min,int&
  &_momentsi,int_momentsr,k_tab,ka_gf,ka_gf_nue,ka_gf_nux,ks_gf,ks_gf_nue,ks_gf_nux,kxx,kxy,kxz,kyy,kyz,kzz,lapm1,lapm1_p,lapm1_p_p&
  &,lapm1_f,lapm1_rhs,lapm1_rot1,lapm1_rot2,lapm1_rot3,lapset,lapset_p,lapset_p_p,lapset_rhs,lapset_timederiv,lapsex,lapsey,lapsez,&
  &mass_star1,mass_star2,mhd_chi_line,mhd_chi_line_p,mhd_chi_line_p_p,mhd_chi_line_rhs,mhd_psi_line,mhd_psi_line_p,mhd_psi_line_p_p&
  &,mhd_psi_line_rhs,mhd_st_x,mhd_st_x_p,mhd_st_x_p_p,mhd_st_x_rhs,mhd_st_y,mhd_st_y_p,mhd_st_y_p_p,mhd_st_y_rhs,mhd_st_z,mhd_st_z_&
  &p,mhd_st_z_p_p,mhd_st_z_rhs,mhd_u_chi,mhd_u_chi_p,mhd_u_chi_p_p,mhd_u_chi_rhs,mhd_u_psi,mhd_u_psi_p,mhd_u_psi_p_p,mhd_u_psi_rhs,&
  &minternal_cold_horiz_VolInt,minternal_cold_r1_VolInt,minternal_cold_r2_VolInt,minternal_cold_r3_VolInt,minternal_horiz_VolInt,mi&
  &nternal_r1_VolInt,minternal_r2_VolInt,minternal_r3_VolInt,mom_const_VolIntD,mom_const_excised_VolIntD,momentsi,momentsi_old,mome&
  &ntsr,momentsr_old,momx_const_VolIntN,momx_const_excised_VolIntN,momy_const_VolIntN,momy_const_excised_VolIntN,momz_const_VolIntN&
  &,momz_const_excised_VolIntN,monopole_VolInt,monopole_outsideBH_VolInt,mu_nu,n_poly,odd_momentsi,odd_momentsi_old,odd_momentsr,od&
  &d_momentsr_old,optd,optd_x,optd_y,optd_z,pfloor_gf,phi,phi_p,phi_p_p,phi_f,phi_rhs,phi_rot1,phi_rot2,phi_rot3,phix,phiy,phiz,psi&
  &,psi0im,psi0im_p,psi0im_p_p,psi0re,psi0re_p,psi0re_p_p,psi6phi,psi6phi_p,psi6phi_p_p,psi6phi_rhs,rad_const,rad_energy_VolInt,rad&
  &_energy_nue_VolInt,rad_energy_nux_VolInt,rbrvec,refbd,refbd_p,refbd_p_p,rho,rhoYe,rhoYe_p,rhoYe_p_p,rhoYe_rhs,rho_b,rho_b_atm_gf&
  &,rho_bl,rho_br,rho_star,rho_star_p,rho_star_p_p,rho_star_rhs,rho_tab,rhob_floor,rhovec,sbt,sbx,sby,sbz,shiftp_rot1,shiftp_rot2,s&
  &hiftp_rot3,shiftr_rot1,shiftr_rot2,shiftr_rot3,shiftth_rot1,shiftth_rot2,shiftth_rot3,shiftx,shiftx_p,shiftx_p_p,shiftx_f,shiftx&
  &_rhs,shiftxt,shiftxt_p,shiftxt_p_p,shiftxt_rhs,shiftxt_timederiv,shifty,shifty_p,shifty_p_p,shifty_f,shifty_rhs,shiftyt,shiftyt_&
  &p,shiftyt_p_p,shiftyt_rhs,shiftyt_timederiv,shiftz,shiftz_p,shiftz_p_p,shiftz_f,shiftz_rhs,shiftzt,shiftzt_p,shiftzt_p_p,shiftzt&
  &_rhs,shiftzt_timederiv,smallb2,st_x,st_y,st_z,surfvec,synch_VolInt,tau,tau_p,tau_p_p,tau_rad,tau_rad_p,tau_rad_p_p,tau_rad_nue,t&
  &au_rad_nue_p,tau_rad_nue_p_p,tau_rad_nue_rhs,tau_rad_nux,tau_rad_nux_p,tau_rad_nux_p_p,tau_rad_nux_rhs,tau_rad_rhs,tau_rhs,temp1&
  &,temp10,temp11,temp12,temp13,temp14,temp15,temp16,temp17,temp18,temp19,temp2,temp20,temp21,temp22,temp23,temp24,temp25,temp26,te&
  &mp3,temp4,temp5,temp6,temp7,temp8,temp9,temp_g00,trK,trK_p,trK_p_p,trK_rhs,trR,trRtilde,u0,v02_rad_nuel,v02_rad_nuer,v02_rad_nux&
  &l,v02_rad_nuxr,v02_radl,v02_radr,v02l,v02r,vvec,vx,vxl,vxr,vy,vyl,vyr,vz,vzl,vzr,w,zeta_rad,zeta_rad_nue,zeta_rad_nux, X0coordin&
  &ates,X1coordinates,X2coordinates,coarse_dx,coarse_dy,coarse_dz,r,x,y,z)
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
  INTEGER X0mhd_shock_private
  INTEGER X1mhd_shock_private
  INTEGER X2mhd_shock_private
  REAL*8 shock_rest_mass_integrand(X0mhd_shock_private,X1mhd_shock_private,X2mhd_shock_private)
  REAL*8 shock_restmass_VolInt
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
  REAL*8  B1os4pi_l
  REAL*8  B1os4pi_r
  REAL*8  B2os4pi_l
  REAL*8  B2os4pi_r
  REAL*8  B3os4pi_l
  REAL*8  B3os4pi_r
  REAL*8  Erad_l
  REAL*8  Erad_r
  REAL*8  Frad0_l
  REAL*8  Frad0_r
  REAL*8  Fradx_l
  REAL*8  Fradx_r
  REAL*8  M0_initial
  REAL*8  P_l
  REAL*8  P_r
  REAL*8  npolyshock
  REAL*8  nx
  REAL*8  ny
  REAL*8  nz
  REAL*8  rho_l
  REAL*8  rho_r
  REAL*8  u1_l
  REAL*8  u1_r
  REAL*8  u2_l
  REAL*8  u2_r
  REAL*8  u3_l
  REAL*8  u3_r
  REAL*8  u_boost
  INTEGER*4 NS
  INTEGER*4 boost_shock
  INTEGER*4 integrate_flag
  INTEGER*4 read_radiation_shock_data
  COMMON /shock_mhdpriv/B1os4pi_l,B1os4pi_r,B2os4pi_l,B2os4pi_r,B3os4pi_l,B3os4pi_r,Erad_l,Erad_r,Frad0_l,Frad0_r,Fradx_l,Fradx_r,M&
  &0_initial,P_l,P_r,npolyshock,nx,ny,nz,rho_l,rho_r,u1_l,u1_r,u2_l,u2_r,u3_l,u3_r,u_boost,NS,boost_shock,integrate_flag,read_radia&
  &tion_shock_data
  REAL*8  CCTKH0
  REAL*8  CCTKH1
  REAL*8  CCTKH2
  integer*8  CCTKH3
  INTEGER*4 cctk_itlast
  INTEGER*4 CCTKH4
  COMMON /CACTUSrest/CCTKH0,CCTKH1,CCTKH2,CCTKH3,cctk_itlast,CCTKH4
  REAL*8  CCTKH7
  REAL*8  CCTKH8
  REAL*8  CCTKH9
  REAL*8  CCTKH10
  REAL*8  CCTKH20
  REAL*8  CCTKH21
  REAL*8  CCTKH22
  REAL*8  CCTKH23
  REAL*8  CCTKH24
  REAL*8  CCTKH25
  REAL*8  CCTKH26
  REAL*8  CCTKH27
  integer*8  CCTKH5
  integer*8  CCTKH6
  integer*8  CCTKH11
  integer*8  CCTKH12
  integer*8  CCTKH19
  INTEGER*4 CCTKH13
  INTEGER*4 CCTKH14
  INTEGER*4 CCTKH15
  INTEGER*4 CCTKH16
  INTEGER*4 CCTKH17
  INTEGER*4 CCTKH18
  COMMON /GRIDrest/CCTKH7,CCTKH8,CCTKH9,CCTKH10,CCTKH20,CCTKH21,CCTKH22,CCTKH23,CCTKH24,CCTKH25,CCTKH26,CCTKH27,CCTKH5,CCTKH6,CCTKH&
  &11,CCTKH12,CCTKH19,CCTKH13,CCTKH14,CCTKH15,CCTKH16,CCTKH17,CCTKH18
  REAL*8  CCTKH33
  REAL*8  CCTKH47
  REAL*8  CCTKH55
  REAL*8  CCTKH57
  REAL*8  CCTKH59
  REAL*8  CCTKH61
  REAL*8  CCTKH63
  REAL*8  CCTKH65
  REAL*8  CCTKH67
  REAL*8  CCTKH69
  REAL*8  CCTKH71
  integer*8  CCTKH30
  integer*8  CCTKH31
  integer*8  CCTKH34
  integer*8  CCTKH37
  integer*8  CCTKH38
  integer*8  CCTKH39
  integer*8  CCTKH42
  integer*8  CCTKH43
  integer*8  CCTKH48
  integer*8  CCTKH49
  integer*8  CCTKH51
  integer*8  CCTKH73
  integer*8  CCTKH75
  integer*8  CCTKH77
  integer*8  CCTKH79
  integer*8  CCTKH80
  integer*8  CCTKH83
  INTEGER*4 CCTKH28
  INTEGER*4 CCTKH29
  INTEGER*4 CCTKH32
  INTEGER*4 CCTKH35
  INTEGER*4 CCTKH36
  INTEGER*4 CCTKH40
  INTEGER*4 CCTKH41
  INTEGER*4 CCTKH44
  INTEGER*4 CCTKH45
  INTEGER*4 CCTKH46
  INTEGER*4 out_every
  INTEGER*4 CCTKH50
  INTEGER*4 CCTKH52
  INTEGER*4 CCTKH53
  INTEGER*4 CCTKH54
  INTEGER*4 CCTKH56
  INTEGER*4 CCTKH58
  INTEGER*4 CCTKH60
  INTEGER*4 CCTKH62
  INTEGER*4 CCTKH64
  INTEGER*4 CCTKH66
  INTEGER*4 CCTKH68
  INTEGER*4 CCTKH70
  INTEGER*4 CCTKH72
  INTEGER*4 CCTKH74
  INTEGER*4 CCTKH76
  INTEGER*4 CCTKH78
  INTEGER*4 CCTKH81
  INTEGER*4 CCTKH82
  COMMON /IOrest/CCTKH33,CCTKH47,CCTKH55,CCTKH57,CCTKH59,CCTKH61,CCTKH63,CCTKH65,CCTKH67,CCTKH69,CCTKH71,CCTKH30,CCTKH31,CCTKH34,CC&
  &TKH37,CCTKH38,CCTKH39,CCTKH42,CCTKH43,CCTKH48,CCTKH49,CCTKH51,CCTKH73,CCTKH75,CCTKH77,CCTKH79,CCTKH80,CCTKH83,CCTKH28,CCTKH29,CC&
  &TKH32,CCTKH35,CCTKH36,CCTKH40,CCTKH41,CCTKH44,CCTKH45,CCTKH46,out_every,CCTKH50,CCTKH52,CCTKH53,CCTKH54,CCTKH56,CCTKH58,CCTKH60,&
  &CCTKH62,CCTKH64,CCTKH66,CCTKH68,CCTKH70,CCTKH72,CCTKH74,CCTKH76,CCTKH78,CCTKH81,CCTKH82
  REAL*8  CCTKH86
  REAL*8  CCTKH87
  REAL*8  K_poly
  REAL*8  M_B
  REAL*8  CCTKH91
  REAL*8  P_fact
  REAL*8  CCTKH92
  REAL*8  CCTKH93
  REAL*8  Sym_Bz
  REAL*8  CCTKH95
  REAL*8  CCTKH99
  REAL*8  CCTKH100
  REAL*8  CCTKH103
  REAL*8  CCTKH109
  REAL*8  CCTKH110
  REAL*8  CCTKH112
  REAL*8  gamma_th
  REAL*8  CCTKH116
  REAL*8  CCTKH118
  REAL*8  CCTKH119
  REAL*8  CCTKH120
  REAL*8  CCTKH121
  REAL*8  CCTKH123
  REAL*8  pfloor
  REAL*8  CCTKH126
  REAL*8  CCTKH127
  REAL*8  CCTKH128
  REAL*8  CCTKH129
  REAL*8  CCTKH130
  REAL*8  CCTKH131
  REAL*8  CCTKH132
  REAL*8  rad_opacity_abs
  REAL*8  rad_opacity_sct
  REAL*8  rho_b_atm
  REAL*8  CCTKH137
  REAL*8  rho_fact
  REAL*8  CCTKH138
  REAL*8  CCTKH139
  REAL*8  CCTKH140
  REAL*8  tau_atm
  REAL*8  tau_fact
  INTEGER*4 CCTKH84
  INTEGER*4 CCTKH85
  INTEGER*4 CCTKH88
  INTEGER*4 CCTKH89
  INTEGER*4 CCTKH90
  INTEGER*4 CCTKH94
  INTEGER*4 CCTKH96
  INTEGER*4 CCTKH97
  INTEGER*4 CCTKH98
  INTEGER*4 CCTKH101
  INTEGER*4 constrained_transport_scheme
  INTEGER*4 CCTKH102
  INTEGER*4 CCTKH104
  INTEGER*4 CCTKH105
  INTEGER*4 enable_HARM_energyvariable
  INTEGER*4 CCTKH106
  INTEGER*4 CCTKH107
  INTEGER*4 enable_primitives_disk
  INTEGER*4 enable_shocktest_primitive_mode
  INTEGER*4 CCTKH108
  INTEGER*4 eps_flag
  INTEGER*4 CCTKH111
  INTEGER*4 force_font_fix_fail
  INTEGER*4 CCTKH113
  INTEGER*4 CCTKH114
  INTEGER*4 CCTKH115
  INTEGER*4 CCTKH117
  INTEGER*4 CCTKH122
  INTEGER*4 neos
  INTEGER*4 CCTKH124
  INTEGER*4 CCTKH125
  INTEGER*4 primitives_solver
  INTEGER*4 rad_closure_scheme
  INTEGER*4 rad_evolve_enable
  INTEGER*4 CCTKH133
  INTEGER*4 CCTKH134
  INTEGER*4 CCTKH135
  INTEGER*4 CCTKH136
  INTEGER*4 tau_stildefix_enable
  INTEGER*4 CCTKH141
  INTEGER*4 CCTKH142
  INTEGER*4 CCTKH143
  COMMON /MHD_EVOLVErest/CCTKH86,CCTKH87,K_poly,M_B,CCTKH91,P_fact,CCTKH92,CCTKH93,Sym_Bz,CCTKH95,CCTKH99,CCTKH100,CCTKH103,CCTKH10&
  &9,CCTKH110,CCTKH112,gamma_th,CCTKH116,CCTKH118,CCTKH119,CCTKH120,CCTKH121,CCTKH123,pfloor,CCTKH126,CCTKH127,CCTKH128,CCTKH129,CC&
  &TKH130,CCTKH131,CCTKH132,rad_opacity_abs,rad_opacity_sct,rho_b_atm,CCTKH137,rho_fact,CCTKH138,CCTKH139,CCTKH140,tau_atm,tau_fact&
  &,CCTKH84,CCTKH85,CCTKH88,CCTKH89,CCTKH90,CCTKH94,CCTKH96,CCTKH97,CCTKH98,CCTKH101,constrained_transport_scheme,CCTKH102,CCTKH104&
  &,CCTKH105,enable_HARM_energyvariable,CCTKH106,CCTKH107,enable_primitives_disk,enable_shocktest_primitive_mode,CCTKH108,eps_flag,&
  &CCTKH111,force_font_fix_fail,CCTKH113,CCTKH114,CCTKH115,CCTKH117,CCTKH122,neos,CCTKH124,CCTKH125,primitives_solver,rad_closure_s&
  &cheme,rad_evolve_enable,CCTKH133,CCTKH134,CCTKH135,CCTKH136,tau_stildefix_enable,CCTKH141,CCTKH142,CCTKH143
  
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
  
  integer, dimension(3) :: ext
  real*8                :: dT,dX,dY,dZ,P_max,rho_max,tau_max,rhos_max,rho_fail_max_step,M_fail_step
  real*8                :: xmin,ymin,zmin,xmax,ymax,zmax
  integer               :: proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax
  integer               :: ierr,index,handle,dummy,glob_imax,glob_jmax,glob_kmax
  REAL*8             :: reduction_value
  integer               :: i,j,k
  real*8                :: norm, w_l,w_r,h_l,h_r,ut_l,ut_r,hm1_l,hm1_r,tau_l,tau_r
  real*8                :: tau_rad_l, tau_rad_r
  real*8                :: ux_l,ux_r,uy_l,uy_r,uz_l,uz_r, F_rad0l, P_radl
  real*8                :: xs,ys,yx,yy,yz,zx,zy,zz
  real*8                :: sbx_l,sbx_r,sby_l,sby_r,sbz_l,sbz_r
  real*8                :: Bx_l,Bx_r,By_l,By_r,Bz_l,Bz_r
  real*8                :: B1_l,B2_l,B3_l,B1_r,B2_r,B3_r
  real*8                 :: A1,A2,A3
  real*8, parameter :: pi = 3.1415926535897932d0
  real*8 :: s4pi,sb0_l,sb0_r,sb1_l,sb1_r,sb2_l,sb2_r,sb3_r,sb3_l,sbs_l,sbs_r
  real*8                :: w_boost, P_save, P_savep1, rho_save, rho_savep1, vx_save, vx_savep1, u0_save, u0_savep1, X_p
  real*8                :: zeta, chi_shock, Fksq, C_A, C_B, P_rad0xl, P_rad0yl, P_rad0zl, P_rad00l
  real*8, dimension(1)  :: xsave, tausave, rhosave,Psave,uxsave,Eradsave,Fradsave
  real*8, dimension(1) :: xread,rho0read,Pread,uxread,Eread,Fxread,tauread
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  integer :: ONE,ZERO
  parameter(ONE = 1.D0, ZERO = 0.D0)
  if(primitives_solver.eq.0 .or. enable_shocktest_primitive_mode.ne.1) then
     write(*,*) "ERROR.  YOU MUST SET primitives_solver=1 or 2 AND enable_shocktest_primitive_mode=1 FOR SHOCK TESTS"
     stop
  end if
  ! sqrt(4 pi) 
  s4pi = 2.d0*sqrt(pi)
  ext = cctk_lsh
  write(*,*) ext(1), ext(2), ext(3)
  if(cctk_lsh(1).le.2*cctk_nghostzones(1)) then
     write(*,*) "ERROR: MUST SET NX>NGHOSTZONESX*2!"
     stop
  else if(cctk_lsh(2).le.2*cctk_nghostzones(2)) then
     write(*,*) "ERROR: MUST SET NY>NGHOSTZONESY*2!"
     stop
  else if(cctk_lsh(3).le.2*cctk_nghostzones(3)) then
     write(*,*) "ERROR: MUST SET NZ>NGHOSTZONESZ*2!"
     stop
  end if
  !!write(*,*) Okay... Looks like you set up the grid correctly EXTENTS:,cctk_lsh
  !!write(*,*) Okay... Looks like you set up the grid correctly GZS:,cctk_nghostzones
  dT = (cctk_delta_time/cctk_timefac)
  dX = (cctk_delta_space(1)/cctk_levfac(1))
  dY = (cctk_delta_space(2)/cctk_levfac(2))
  dZ = (cctk_delta_space(3)/cctk_levfac(3))
  ! u^0 on the left (w_l) and right (w_r) side of the shock
  w_l = sqrt(1.0 + u1_l**2 + u2_l**2 + u3_l**2)
  w_r = sqrt(1.0 + u1_r**2 + u2_r**2 + u3_r**2)
  B1_l = B1os4pi_l * s4pi
  B2_l = B2os4pi_l * s4pi
  B3_l = B3os4pi_l * s4pi
  B1_r = B1os4pi_r * s4pi
  B2_r = B2os4pi_r * s4pi
  B3_r = B3os4pi_r * s4pi
  ! Compute b^mu
  sb0_l = u1_l*B1os4pi_l + u2_l*B2os4pi_l + u3_l*B3os4pi_l
  sb0_r = u1_r*B1os4pi_r + u2_r*B2os4pi_r + u3_r*B3os4pi_r
  sb1_l = (B1os4pi_l + sb0_l*u1_l)/w_l
  sb1_r = (B1os4pi_r + sb0_r*u1_r)/w_r
  sb2_l = (B2os4pi_l + sb0_l*u2_l)/w_l
  sb2_r = (B2os4pi_r + sb0_r*u2_r)/w_r
  sb3_l = (B3os4pi_l + sb0_l*u3_l)/w_l
  sb3_r = (B3os4pi_r + sb0_r*u3_r)/w_r
  sbs_l = sb1_l**2 + sb2_l**2 + sb3_l**2 - sb0_l**2
  sbs_r = sb1_r**2 + sb2_r**2 + sb3_r**2 - sb0_r**2
  hm1_l = ( npolyshock + 1.d0 ) * P_l / rho_l
  hm1_r = ( npolyshock + 1.d0 ) * P_r / rho_r
  h_l = 1.d0 + hm1_l
  h_r = 1.d0 + hm1_r
  tau_l = rho_l*w_l*(h_l*w_l - 1.d0) - P_l + sbs_l*(w_l*w_l-0.5d0) - sb0_l**2
  tau_r = rho_r*w_r*(h_r*w_r - 1.d0) - P_r + sbs_r*(w_r*w_r-0.5d0) - sb0_r**2
  ! just to make sure, normalize n:
  !
  norm = sqrt(nx*nx + ny*ny + nz*nz)
  nx = nx/norm
  ny = ny/norm
  nz = nz/norm
  if (nx==0.d0) then
     ! In this case, the relationship between the coordinates in the shock frame 
     ! and the lab frame is \hat{xs}=\hat{n}, \hat{ys}=-\hat{x}, 
     ! \hat{zs} = \hat{xs} cross \hat{ys}
     yx = -1.d0
     yy = 0.d0
     yz = 0.d0
     zx = 0.d0
     zy = -nz
     zz = ny
     ux_l = -u2_l
     ux_r = -u2_r
     uy_l = ny*u1_l - nz*u3_l
     uy_r = ny*u1_r - nz*u3_r
     uz_l = nz*u1_l + ny*u3_l
     uz_r = nz*u1_r + ny*u3_r
     sbx_l = -sb2_l
     sbx_r = -sb2_r
     sby_l = ny*sb1_l - nz*sb3_l
     sby_r = ny*sb1_r - nz*sb3_r
     sbz_l = nz*sb1_l + ny*sb3_l
     sbz_r = nz*sb1_r + ny*sb3_r
     Bx_l = -B2os4pi_l*s4pi
     Bx_r = -B2os4pi_r*s4pi
     By_l = (ny*B1os4pi_l - nz*B3os4pi_l)*s4pi
     By_r = (ny*B1os4pi_r - nz*B3os4pi_r)*s4pi
     Bz_l = (nz*B1os4pi_l + ny*B3os4pi_l)*s4pi
     Bz_r = (nz*B1os4pi_r + ny*B3os4pi_r)*s4pi
  else if (ny==0.d0) then
     ! In this case, the relationship between the coordinates in the shock frame
     ! and the lab frame is \hat{xs}=\hat{n}, \hat{ys}=\hat{y},
     ! \hat{zs} = \hat{xs} cross \hat{ys}
     yx = 0.d0
     yy = 1.d0
     yz = 0.d0
     zx = -nz
     zy = 0.d0
     zz = nx
     ux_l = nx*u1_l - nz*u3_l
     ux_r = nx*u1_r - nz*u3_r
     uy_l = u2_l
     uy_r = u2_r
     uz_l = nx*u3_l + nz*u1_l
     uz_r = nx*u3_r + nz*u1_r
     sbx_l = nx*sb1_l - nz*sb3_l
     sbx_r = nx*sb1_r - nz*sb3_r
     sby_l = sb2_l
     sby_r = sb2_r
     sbz_l = nx*sb3_l + nz*sb1_l
     sbz_r = nx*sb3_r + nz*sb1_r
     Bx_l = (nx*B1os4pi_l - nz*B3os4pi_l)*s4pi
     Bx_r = (nx*B1os4pi_r - nz*B3os4pi_r)*s4pi
     By_l = B2os4pi_l*s4pi
     By_r = B2os4pi_r*s4pi
     Bz_l = (nx*B3os4pi_l + nz*B1os4pi_l)*s4pi
     Bz_r = (nx*B3os4pi_r + nz*B1os4pi_r)*s4pi
  else if (nz==0.d0) then
     ! In this case, the relationship between the coordinates in the shock frame
     ! and the lab frame is \hat{xs}=\hat{n}, \hat{zs}=\hat{z},
     ! \hat{ys} = \hat{z} cross \hat{n}
     yx = -ny
     yy = nx
     yz = 0.d0
     zx = 1.d0
     zy = 0.d0
     zz = 0.d0
     ux_l = nx*u1_l - ny*u2_l
     ux_r = nx*u1_r - ny*u2_r
     uy_l = nx*u2_l + ny*u1_l
     uy_r = nx*u2_r + ny*u1_r
     uz_l = u3_l
     uz_r = u3_r
     sbx_l = nx*sb1_l - ny*sb2_l
     sbx_r = nx*sb1_r - ny*sb2_r
     sby_l = nx*sb2_l + ny*sb1_l
     sby_r = nx*sb2_r + ny*sb1_r
     sbz_l = sb3_l
     sbz_r = sb3_r
     Bx_l = (nx*B1os4pi_l - ny*B2os4pi_l)*s4pi
     Bx_r = (nx*B1os4pi_r - ny*B2os4pi_r)*s4pi
     By_l = (nx*B2os4pi_l + ny*B1os4pi_l)*s4pi
     By_r = (nx*B2os4pi_r + ny*B1os4pi_r)*s4pi
     Bz_l = B3os4pi_l*s4pi
     Bz_r = B3os4pi_r*s4pi
  else
     write(*,*) 'At least one of the nx, ny or ny must be set to 0'
     stop
  end if
if (read_radiation_shock_data .eq. 1) then
 write(*,*) 'DIRECTLY READ SHOCK DATA FROM INPUT FILE!!'
 write(*,*) 'rad_const =', rad_const
! here we call radiationshock_readfiles to read in rho0, P, E_rad, F_rad, vx, and tau.
! so at this point the 6 GFs should already be filled.
! we just need to fill in others.
        do k=1,ext(3)
          do j=1,ext(2)
            do i=1,ext(1)
                if (rho_b(i,j,k) == 0.d0) then
                    write(*,*) " At the begining of initialdata_local.F90, zero rho is found!!!!", i,j,k
                end if
                ! here we need to distinguish wether the shock is boosted!!!!
                w(i,j,k) = 1.d0/sqrt(1.d0 - vx(i,j,k)**2)
                if (boost_shock .eq. 1) then
                ! u0(i,j,k) = w(i,j,k)/sqrt(1.d0 - u_boost**2)*(1.d0 - u_boost*vx(i,j,k))
                ! The following expression is correct because vx(i,j,k) is already boosted in readinputfile_read.C                  
                u0(i,j,k) = 1.d0/sqrt(1.d0 - vx(i,j,k)**2)
                 ! F_radx(i,j,k) here is already boosted in readinputfile_read.C so only a factor of -u_boost is needed.
                !F_rad0(i,j,k) = u_boost* F_radx(i,j,k)
                F_rad0(i,j,k) = F_radx(i,j,k)*vx(i,j,k)
                else
                u0(i,j,k) = w(i,j,k)
                ! F_rad0 = - F_radx*u_x/u_0 = - F_radx*v_x = F_radx*vx (u_0 = - u0 / v_x = -vx) 
                F_rad0(i,j,k) = F_radx(i,j,k)*vx(i,j,k)
                end if
if (i==4.and.j==4.and.k==4) then
    write(*,*) " Inside shock_mhd initial_data_local  u0(i,j,k) is ", u0(i,j,k)
    write(*,*) " w(i,j,k), u_boost, vx(i,j,k) are",  w(i,j,k), u_boost, vx(i,j,k)
    write(*,*) "E_rad and F_radx are", E_rad(i,j,k), F_radx(i,j,k)
end if
                !vx(i,j,k) = ux(i,j,k)/w(i,j,k)
                vy(i,j,k) = uy_l/u0(i,j,k)
                vz(i,j,k) = uz_l/u0(i,j,k)
                F_rady(i,j,k) = ZERO
                F_radz(i,j,k) = ZERO
                P_radl = E_rad(i,j,k)/3.0
                rho_star(i,j,k) = rho_b(i,j,k) * u0(i,j,k)
                ! Lunan : add this if statement to aviod dividing by zero.
                if (rho_b(i,j,k) == 0.d0 .and. P(i,j,k) == 0.d0) then
                h(i,j,k) = 1.d0 + ( npolyshock + 1.d0 )
                else
                h(i,j,k) = 1.d0 + ( npolyshock + 1.d0 ) * P(i,j,k) / rho_b(i,j,k)
                end if
                tau(i,j,k) = rho_star(i,j,k)*(h(i,j,k)*u0(i,j,k)-1.d0)-P(i,j,k)+sbs_l*(u0(i,j,k)*u0(i,j,k)-0.5d0)-sb0_l**2
                st_x(i,j,k)  = rho_star(i,j,k)*h(i,j,k)*vx(i,j,k)/sqrt(1.d0-vx(i,j,k)**2)
                !st_x(i,j,k)  = rho_star(i,j,k)*h(i,j,k)*vx(i,j,k)/(1.d0-vx(i,j,k)**2)
                st_y(i,j,k)  = rho_star(i,j,k)*h(i,j,k)*uy_l
                st_z(i,j,k)  = rho_star(i,j,k)*h(i,j,k)*uz_l
                mhd_st_x(i,j,k) = st_x(i,j,k) + u0(i,j,k)*sbs_l*vx(i,j,k)/sqrt(1-vx(i,j,k)**2) - sb0_l*sbx_l
                mhd_st_y(i,j,k) = st_y(i,j,k) + u0(i,j,k)*sbs_l*uy_l - sb0_l*sby_l
                mhd_st_z(i,j,k) = st_z(i,j,k) + u0(i,j,k)*sbs_l*uz_l - sb0_l*sbz_l
                Bx(i,j,k)   = Bx_l
                By(i,j,k)   = By_l
                Bz(i,j,k)   = Bz_l
                Bxtilde(i,j,k) = Bx_l
                Bytilde(i,j,k) = By_l
                Bztilde(i,j,k) = Bz_l
                sbt(i,j,k) = sb0_l
                sbx(i,j,k) = sbx_l
                sby(i,j,k) = sby_l
                sbz(i,j,k) = sbz_l
                !T_fluid(i,j,k) = (E_rad(i,j,k)/rad_const)**(0.25)
                T_fluid(i,j,k) = P(i,j,k)/rho(i,j,k)
                if (rad_evolve_enable .eq. 1) then
                   if(rad_closure_scheme .eq. 0) then
                      tau_rad(i,j,k) = (E_rad(i,j,k) + P_radl)*u0(i,j,k)*u0(i,j,k) + 2.0*u0(i,j,k)*F_rad0(i,j,k)-P_radl
                      ! In flat space, F_rad_alpha = g_{alpha beta} F_rad^alpha
                      ! Notice that: F_rad_0 = -F_rad0, F_rad_i = F_radi, u_0 = -u0, v_i = -vi, u_i = ui
                      S_rad_x(i,j,k) = (E_rad(i,j,k) + P_radl)*u0(i,j,k)*u0(i,j,k)*vx(i,j,k) + F_rad0(i,j,k)*u0(i,j,k)*vx(i,j,k) + &
  &F_radx(i,j,k)*u0(i,j,k)
                      S_rad_y(i,j,k) = (E_rad(i,j,k) + P_radl)*u0(i,j,k)*uy_l + F_rad0(i,j,k)*uy_l + F_rady(i,j,k)*u0(i,j,k)
                      S_rad_z(i,j,k) = (E_rad(i,j,k) + P_radl)*u0(i,j,k)*uz_l + F_rad0(i,j,k)*uz_l + F_radz(i,j,k)*u0(i,j,k)
                   else
                      zeta = sqrt((-F_rad0(i,j,k)**2 + F_radx(i,j,k)**2 + F_rady(i,j,k)**2 + F_radz(i,j,k)**2)/E_rad(i,j,k)**2)
                      chi_shock  = 1.0/3.0 + zeta**2*(6.0-2.0*zeta+6.0*zeta**2)/15.0
                      chi_rad(i,j,k)=chi_shock
                      Fksq = F_radx(i,j,k)**2 + F_rady(i,j,k)**2 + F_radz(i,j,k)**2
                      C_A = (3.0*chi_shock -1.0)/2.0
                      C_B = 1.5*(1.0-chi_shock)
                      P_radxx(i,j,k) = F_radx(i,j,k)**2*E_rad(i,j,k)/Fksq*C_A + E_rad(i,j,k)*(1.0+(u0(i,j,k)*vx(i,j,k))**2)/3.0*C_B
                      P_radyy(i,j,k) = F_rady(i,j,k)**2*E_rad(i,j,k)/Fksq*C_A + E_rad(i,j,k)*(1.0+uy_l**2)/3.0*C_B
                      P_radzz(i,j,k) = F_radz(i,j,k)**2*E_rad(i,j,k)/Fksq*C_A + E_rad(i,j,k)*(1.0+uz_l**2)/3.0*C_B
                      P_radxy(i,j,k) = F_radx(i,j,k)*F_rady(i,j,k)/Fksq*C_A + E_rad(i,j,k)* u0(i,j,k)*vx(i,j,k)*uy_l/3.0*C_B
                      P_radxz(i,j,k) = F_radx(i,j,k)*F_radz(i,j,k)/Fksq*C_A + E_rad(i,j,k)* u0(i,j,k)*vz(i,j,k)*uz_l/3.0*C_B
                      P_radyz(i,j,k) = F_rady(i,j,k)*F_radz(i,j,k)/Fksq*C_A + E_rad(i,j,k)* uy_l*uz_l/3.0*C_B
                    ! In flat space, P_rad^i_j = P_rad^{ij} 
                      P_rad0xl = (P_radxx(i,j,k)*u0(i,j,k)*vx(i,j,k) + P_radxy(i,j,k) * uy_l + P_radxz(i,j,k) * uz_l)/u0(i,j,k)
                      P_rad0yl = (P_radxy(i,j,k)*u0(i,j,k)*vx(i,j,k) + P_radyy(i,j,k) * uy_l + P_radyz(i,j,k) * uz_l)/u0(i,j,k)
                      P_rad0yl = (P_radxz(i,j,k)*u0(i,j,k)*vx(i,j,k) + P_radyz(i,j,k) * uy_l + P_radzz(i,j,k) * uz_l)/u0(i,j,k)
                      P_rad00l = (P_rad0xl*u0(i,j,k)*vx(i,j,k) + P_rad0yl * uy_l + P_rad0zl * uz_l)/u0(i,j,k)
                      tau_rad(i,j,k) = E_rad(i,j,k)*u0(i,j,k)*u0(i,j,k) + 2.0*u0(i,j,k)*F_rad0(i,j,k) + P_rad00l
                      S_rad_x(i,j,k) = E_rad(i,j,k)*u0(i,j,k)*u0(i,j,k)*vx(i,j,k) + F_rad0(i,j,k)*u0(i,j,k)*vx(i,j,k) &
                                       + F_radx(i,j,k)*u0(i,j,k) + P_rad0xl
                      S_rad_y(i,j,k) = E_rad(i,j,k)*u0(i,j,k)*uy_l + F_rad0(i,j,k)*uy_l + F_rady(i,j,k)*u0(i,j,k) + P_rad0yl
                      S_rad_z(i,j,k) = E_rad(i,j,k)*u0(i,j,k)*uz_l + F_rad0(i,j,k)*uz_l + F_radz(i,j,k)*u0(i,j,k) + P_rad0zl
                   end if
                else
                   tau_rad(i,j,k) = 0.d0
                   S_rad_x(i,j,k) = 0.d0
                   S_rad_y(i,j,k) = 0.d0
                   S_rad_z(i,j,k) = 0.d0
                end if
           ! D.J. Griffths Introduction to Electrodynamics]
           ! Note that Ai are staggered
                A1 = 0.d0
                xs = X(i,j,k)*nx + (Y(i,j,k)+0.5d0*dY)*ny + (Z(i,j,k)+0.5d0*dZ)*nz
                ys = X(i,j,k)*yx + (Y(i,j,k)+0.5d0*dY)*yy + (Z(i,j,k)+0.5d0*dZ)*yz
                if (xs .le. 0.d0) then
                   A2 = B3_l*xs
                         A3 = B1_l*ys - B2_l*xs
                else
                   A2 = B3_r*xs
                         A3 = B1_r*ys - B2_r*xs
                end if
                   Ax(i,j,k) = nx*A1 + yx*A2 + zx*A3
                   xs = (X(i,j,k)+0.5d0*dX)*nx + Y(i,j,k)*ny + (Z(i,j,k)+0.5d0*dZ)*nz
                   ys = (X(i,j,k)+0.5d0*dX)*yx + Y(i,j,k)*yy + (Z(i,j,k)+0.5d0*dZ)*yz
                   if (xs .le. 0.d0) then
                         A2 = B3_l*xs
                         A3 = B1_l*ys - B2_l*xs
                   else
                   A2 = B3_r*xs
                         A3 = B1_r*ys - B2_r*xs
                   end if
                      Ay(i,j,k) = ny*A1 + yy*A2 + zy*A3
                      xs = (X(i,j,k)+0.5d0*dX)*nx + (Y(i,j,k)+0.5d0*dY)*ny + Z(i,j,k)*nz
                     ys = (X(i,j,k)+0.5d0*dX)*yx + (Y(i,j,k)+0.5d0*dY)*yy + Z(i,j,k)*yz
                   if (xs .le. 0.d0) then
                   A2 = B3_l*xs
                         A3 = B1_l*ys - B2_l*xs
                else
                   A2 = B3_r*xs
                         A3 = B1_r*ys - B2_r*xs
                end if
                Az(i,j,k) = nz*A1 + yz*A3 + zz*A3
        end do
     end do
  end do
else
 write(*,*) 'CALCULATE RAD SHOCK DATA USING PARFILE!'
 if (rad_evolve_enable .eq. 1) then
       write(*,*) 'integrate_flag is', integrate_flag
  call  initialdata_radiation (ext, x, npolyshock, &
        rad_opacity_abs, rad_opacity_sct, &
        rho_l, P_l, u1_l, Erad_l, Fradx_l, Fradx_r, &
        cctkGH, cctk_lsh, integrate_flag, xsave, tausave, &
        rhosave, Psave, uxsave, Eradsave, Fradsave)
          do k=1,ext(3)
               do j=1,ext(2)
            do i=1,ext(1)
              xs = X(i,j,k)*nx + Y(i,j,k)*ny + Z(i,j,k)*nz
              if (xs .le. 0.D0) then
             !
             ! left side of partition
                   w(i,j,k) = sqrt(1.d0 + uxsave(i)**2)
                rho_b(i,j,k) = rhosave(i)
                   P(i,j,k) = Psave(i)
                   vx(i,j,k) = uxsave(i)/w(i,j,k)
                   vy(i,j,k) = uy_l/w(i,j,k)
                   vz(i,j,k) = uz_l/w(i,j,k)
                   E_rad(i,j,k) = Eradsave(i)
                   F_radx(i,j,k) = Fradsave(i)
                F_rad0(i,j,k) = ZERO
                F_rady(i,j,k) = ZERO
                F_radz(i,j,k) = ZERO
                P_radl = E_rad(i,j,k)/3.0
                   tau (i,j,k) = tausave(i)
                   rho_star(i,j,k) = rho_b(i,j,k) * w(i,j,k)
                   h(i,j,k) = 1.d0 + ( npolyshock + 1.d0 ) * P(i,j,k) / rho_b(i,j,k)
                   st_x(i,j,k)  = rho_star(i,j,k)*h(i,j,k)*uxsave(i)
                   st_y(i,j,k)  = rho_star(i,j,k)*h(i,j,k)*uy_l
                  st_z(i,j,k)  = rho_star(i,j,k)*h(i,j,k)*uz_l
                   mhd_st_x(i,j,k) = st_x(i,j,k) + w(i,j,k)*sbs_l*uxsave(i) - sb0_l*sbx_l
                   mhd_st_y(i,j,k) = st_y(i,j,k) + w(i,j,k)*sbs_l*uy_l - sb0_l*sby_l
                   mhd_st_z(i,j,k) = st_z(i,j,k) + w(i,j,k)*sbs_l*uz_l - sb0_l*sbz_l
                   u0(i,j,k) = w (i,j,k)
                Bx(i,j,k)   = Bx_l
                By(i,j,k)   = By_l
                   Bz(i,j,k)   = Bz_l
                   Bxtilde(i,j,k) = Bx_l
                   Bytilde(i,j,k) = By_l
                   Bztilde(i,j,k) = Bz_l
                   sbt(i,j,k) = sb0_l
                   sbx(i,j,k) = sbx_l
                   sby(i,j,k) = sby_l
                   sbz(i,j,k) = sbz_l
                !T_fluid(i,j,k) = (E_rad(i,j,k)/rad_const)**(0.25)
                T_fluid(i,j,k) = P(i,j,k)/rho(i,j,k)
                tau_rad(i,j,k) = w(i,j,k)*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*u0(i,j,k) + 2.0*u0(i,j,k)*F_rad0(i,j,k)-P_radl)
                S_rad_x(i,j,k) = w(i,j,k)*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*uxsave(i) + F_rad0(i,j,k)*uxsave(i) + F_radx(i,j,k)*u0&
  &(i,j,k))
                S_rad_y(i,j,k) = w(i,j,k)*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*uy_l + F_rad0(i,j,k)*uy_l + F_rady(i,j,k)*u0(i,j,k))
                S_rad_z(i,j,k) = w(i,j,k)*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*uz_l + F_rad0(i,j,k)*uz_l + F_radz(i,j,k)*u0(i,j,k))
           else
              !
              ! right side of partition
              !
              P(i,j,k)    = P_r
              rho_b(i,j,k) = rho_r
              rho_star(i,j,k) = rho_r * w_r
              tau(i,j,k)   = tau_r
              h(i,j,k) = h_r
              w(i,j,k)    = rho_star(i,j,k) *w_r
              st_x(i,j,k)  = rho_star(i,j,k) *h_r*ux_r
              st_y(i,j,k)  = rho_star(i,j,k) *h_r*uy_r
              st_z(i,j,k)  = rho_star(i,j,k) *h_r*uz_r
              mhd_st_x(i,j,k) = st_x(i,j,k) + w_r*sbs_r*ux_r - sb0_r*sbx_r
              mhd_st_y(i,j,k) = st_y(i,j,k) + w_r*sbs_r*uy_r - sb0_r*sby_r
              mhd_st_z(i,j,k) = st_z(i,j,k) + w_r*sbs_r*uz_r - sb0_r*sbz_r
              vx(i,j,k)   = ux_r/w_r
              vy(i,j,k)   = uy_r/w_r
              vz(i,j,k)   = uz_r/w_r
              u0(i,j,k)   = w_r
              Bx(i,j,k)   = Bx_r
              By(i,j,k)   = By_r
              Bz(i,j,k)   = Bz_r
              Bxtilde(i,j,k) = Bx_r
              Bytilde(i,j,k) = By_r
              Bztilde(i,j,k) = Bz_r
              sbt(i,j,k) = sb0_r
              sbx(i,j,k) = sbx_r
              sby(i,j,k) = sby_r
              sbz(i,j,k) = sbz_r
              E_rad(i,j,k) = Erad_r
              F_rad0(i,j,k) = ZERO
              F_radx(i,j,k) = Fradx_r
              F_rady(i,j,k) = ZERO
              F_radz(i,j,k) = ZERO
              P_radl = Erad_r/3.0
              !T_fluid(i,j,k) = P(i,j,k)*M_B/rho_star(i,j,k)
              T_fluid(i,j,k) = (E_rad(i,j,k)/rad_const)**(0.25)
              tau_rad(i,j,k) = w_r*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*u0(i,j,k) + 2.0*u0(i,j,k)*F_rad0(i,j,k)-P_radl)
              S_rad_x(i,j,k) = w_r*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*ux_r + F_rad0(i,j,k)*ux_r + F_radx(i,j,k)*u0(i,j,k))
              S_rad_y(i,j,k) = w_r*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*uy_r + F_rad0(i,j,k)*uy_r + F_rady(i,j,k)*u0(i,j,k))
              S_rad_z(i,j,k) = w_r*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*uz_r + F_rad0(i,j,k)*uz_r + F_radz(i,j,k)*u0(i,j,k))
           end if
           ! D.J. Griffths Introduction to Electrodynamics]
           ! Note that Ai are staggered
           A1 = 0.d0
           xs = X(i,j,k)*nx + (Y(i,j,k)+0.5d0*dY)*ny + (Z(i,j,k)+0.5d0*dZ)*nz
           ys = X(i,j,k)*yx + (Y(i,j,k)+0.5d0*dY)*yy + (Z(i,j,k)+0.5d0*dZ)*yz
           if (xs .le. 0.d0) then
              A2 = B3_l*xs
              A3 = B1_l*ys - B2_l*xs
           else
              A2 = B3_r*xs
              A3 = B1_r*ys - B2_r*xs
           end if
           Ax(i,j,k) = nx*A1 + yx*A2 + zx*A3
           xs = (X(i,j,k)+0.5d0*dX)*nx + Y(i,j,k)*ny + (Z(i,j,k)+0.5d0*dZ)*nz
           ys = (X(i,j,k)+0.5d0*dX)*yx + Y(i,j,k)*yy + (Z(i,j,k)+0.5d0*dZ)*yz
           if (xs .le. 0.d0) then
              A2 = B3_l*xs
              A3 = B1_l*ys - B2_l*xs
           else
              A2 = B3_r*xs
              A3 = B1_r*ys - B2_r*xs
           end if
           Ay(i,j,k) = ny*A1 + yy*A2 + zy*A3
           xs = (X(i,j,k)+0.5d0*dX)*nx + (Y(i,j,k)+0.5d0*dY)*ny + Z(i,j,k)*nz
           ys = (X(i,j,k)+0.5d0*dX)*yx + (Y(i,j,k)+0.5d0*dY)*yy + Z(i,j,k)*yz
           if (xs .le. 0.d0) then
              A2 = B3_l*xs
              A3 = B1_l*ys - B2_l*xs
           else
              A2 = B3_r*xs
              A3 = B1_r*ys - B2_r*xs
           end if
           Az(i,j,k) = nz*A1 + yz*A3 + zz*A3
           end do
     end do
  end do
  else
  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
              xs = X(i,j,k)*nx + Y(i,j,k)*ny + Z(i,j,k)*nz
              if (xs .le. 0.D0) then
             !
             ! left side of partition
             ! If radiation is involved, we calcuate the 1D profile (in x-direction) using initialdata_radiation and fill
             !the 3D arrays with this 1D array.
              P(i,j,k)    = P_l
              rho_b(i,j,k) = rho_l
              rho_star(i,j,k) = rho_l * w_l
              tau(i,j,k)   = tau_l
              h(i,j,k) = h_l
              w(i,j,k)    = rho_star(i,j,k) *w_l
              st_x(i,j,k)  = rho_star(i,j,k) *h_l*ux_l
              st_y(i,j,k)  = rho_star(i,j,k) *h_l*uy_l
              st_z(i,j,k)  = rho_star(i,j,k) *h_l*uz_l
              mhd_st_x(i,j,k) = st_x(i,j,k) + w_l*sbs_l*ux_l - sb0_l*sbx_l
              mhd_st_y(i,j,k) = st_y(i,j,k) + w_l*sbs_l*uy_l - sb0_l*sby_l
              mhd_st_z(i,j,k) = st_z(i,j,k) + w_l*sbs_l*uz_l - sb0_l*sbz_l
              vx(i,j,k)   = ux_l/w_l
              vy(i,j,k)   = uy_l/w_l
              vz(i,j,k)   = uz_l/w_l
              u0(i,j,k)   = w_l
                 Bx(i,j,k)   = Bx_l
              By(i,j,k)   = By_l
              Bz(i,j,k)   = Bz_l
              Bxtilde(i,j,k) = Bx_l
              Bytilde(i,j,k) = By_l
              Bztilde(i,j,k) = Bz_l
              sbt(i,j,k) = sb0_l
              sbx(i,j,k) = sbx_l
              sby(i,j,k) = sby_l
              sbz(i,j,k) = sbz_l
              E_rad(i,j,k) = Erad_l
              F_rad0(i,j,k) = ZERO
              F_radx(i,j,k) = Fradx_l
              F_rady(i,j,k) = ZERO
              F_radz(i,j,k) = ZERO
              P_radl = Erad_l/3.0
              !T_fluid(i,j,k) = (E_rad(i,j,k)/rad_const)**(0.25)
              T_fluid(i,j,k) = P(i,j,k)/rho(i,j,k)
              tau_rad(i,j,k) = w_l*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*u0(i,j,k) + 2.0*u0(i,j,k)*F_rad0(i,j,k)-P_radl)
              S_rad_x(i,j,k) = w_l*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*ux_l + F_rad0(i,j,k)*ux_l + F_radx(i,j,k)*u0(i,j,k))
              S_rad_y(i,j,k) = w_l*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*uy_l + F_rad0(i,j,k)*uy_l + F_rady(i,j,k)*u0(i,j,k))
              S_rad_z(i,j,k) = w_l*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*uz_l + F_rad0(i,j,k)*uz_l + F_radz(i,j,k)*u0(i,j,k))
           else
              !
              ! right side of partition
              !
              P(i,j,k)    = P_r
              rho_b(i,j,k) = rho_r
              rho_star(i,j,k) = rho_r * w_r
              tau(i,j,k)   = tau_r
              h(i,j,k) = h_r
              w(i,j,k)    = rho_star(i,j,k) *w_r
              st_x(i,j,k)  = rho_star(i,j,k) *h_r*ux_r
              st_y(i,j,k)  = rho_star(i,j,k) *h_r*uy_r
              st_z(i,j,k)  = rho_star(i,j,k) *h_r*uz_r
              mhd_st_x(i,j,k) = st_x(i,j,k) + w_r*sbs_r*ux_r - sb0_r*sbx_r
              mhd_st_y(i,j,k) = st_y(i,j,k) + w_r*sbs_r*uy_r - sb0_r*sby_r
              mhd_st_z(i,j,k) = st_z(i,j,k) + w_r*sbs_r*uz_r - sb0_r*sbz_r
              vx(i,j,k)   = ux_r/w_r
              vy(i,j,k)   = uy_r/w_r
              vz(i,j,k)   = uz_r/w_r
              u0(i,j,k)   = w_r
              Bx(i,j,k)   = Bx_r
              By(i,j,k)   = By_r
              Bz(i,j,k)   = Bz_r
              Bxtilde(i,j,k) = Bx_r
              Bytilde(i,j,k) = By_r
              Bztilde(i,j,k) = Bz_r
              sbt(i,j,k) = sb0_r
              sbx(i,j,k) = sbx_r
              sby(i,j,k) = sby_r
              sbz(i,j,k) = sbz_r
              E_rad(i,j,k) = Erad_r
              F_rad0(i,j,k) = ZERO
              F_radx(i,j,k) = Fradx_r
                F_rady(i,j,k) = ZERO
              F_radz(i,j,k) = ZERO
              P_radl = Erad_r/3.0
              !T_fluid(i,j,k) = (E_rad(i,j,k)/rad_const)**(0.25)
              T_fluid(i,j,k) = P(i,j,k)/rho(i,j,k)
              tau_rad(i,j,k) = w_r*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*u0(i,j,k) + 2.0*u0(i,j,k)*F_rad0(i,j,k)-P_radl)
              S_rad_x(i,j,k) = w_r*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*ux_r + F_rad0(i,j,k)*ux_r + F_radx(i,j,k)*u0(i,j,k))
              S_rad_y(i,j,k) = w_r*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*uy_r + F_rad0(i,j,k)*uy_r + F_rady(i,j,k)*u0(i,j,k))
              S_rad_z(i,j,k) = w_r*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*uz_r + F_rad0(i,j,k)*uz_r + F_radz(i,j,k)*u0(i,j,k))
        end if
           ! Setup vector potential [use the expression of problem 5.33 in
           ! D.J. Griffths Introduction to Electrodynamics]
           ! Note that Ai are staggered
           A1 = 0.d0
           xs = X(i,j,k)*nx + (Y(i,j,k)+0.5d0*dY)*ny + (Z(i,j,k)+0.5d0*dZ)*nz
           ys = X(i,j,k)*yx + (Y(i,j,k)+0.5d0*dY)*yy + (Z(i,j,k)+0.5d0*dZ)*yz
           if (xs .le. 0.d0) then
              A2 = B3_l*xs
              A3 = B1_l*ys - B2_l*xs
           else
              A2 = B3_r*xs
              A3 = B1_r*ys - B2_r*xs
           end if
           Ax(i,j,k) = nx*A1 + yx*A2 + zx*A3
           xs = (X(i,j,k)+0.5d0*dX)*nx + Y(i,j,k)*ny + (Z(i,j,k)+0.5d0*dZ)*nz
           ys = (X(i,j,k)+0.5d0*dX)*yx + Y(i,j,k)*yy + (Z(i,j,k)+0.5d0*dZ)*yz
           if (xs .le. 0.d0) then
              A2 = B3_l*xs
              A3 = B1_l*ys - B2_l*xs
           else
              A2 = B3_r*xs
              A3 = B1_r*ys - B2_r*xs
           end if
           Ay(i,j,k) = ny*A1 + yy*A2 + zy*A3
           xs = (X(i,j,k)+0.5d0*dX)*nx + (Y(i,j,k)+0.5d0*dY)*ny + Z(i,j,k)*nz
           ys = (X(i,j,k)+0.5d0*dX)*yx + (Y(i,j,k)+0.5d0*dY)*yy + Z(i,j,k)*yz
           if (xs .le. 0.d0) then
              A2 = B3_l*xs
              A3 = B1_l*ys - B2_l*xs
           else
              A2 = B3_r*xs
              A3 = B1_r*ys - B2_r*xs
           end if
           Az(i,j,k) = nz*A1 + yz*A3 + zz*A3
           end do
             end do
      end do
  end if
end if
  !tau_p = tau
  !h_p = h
  neos = 1
  write(*,*) "Assuming polytropic EOS..."
  write(*,*) "Polytropic index = ",npolyshock
  write(*,*) "Polytropic constant K = ",K_poly
  rho_tab(1) = 0.03d0
  P_tab(1) = K_poly * rho_tab(1)**(1.0+1.0/npolyshock)
  eps_tab(1) = npolyshock*P_tab(1)/rho_tab(1)
  gamma_tab(1) = 1.d0 + 1.d0/npolyshock
  k_tab(1) = K_poly
  gamma_tab(2) = 1.d0 + 1.d0/npolyshock
  k_tab(2) = K_poly
  gamma_th = 1.d0 + 1.d0/npolyshock
  !======================================
  ! Set metric data to Minkowski
  !======================================
  phi = ZERO
  phix = ZERO
  phiy = ZERO
  phiz = ZERO
  psi = exp(phi)
  lapm1 = ZERO
  gxx = ONE
  gyy = ONE
  gzz = ONE
  gxy = ZERO
  gxz = ZERO
  gyz = ZERO
  gupxx = ONE
  gupyy = ONE
  gupzz = ONE
  gupxy = ZERO
  gupxz = ZERO
  gupyz = ZERO
  Gammax = ZERO
  Gammay = ZERO
  Gammaz = ZERO
  !======================================
  ! Set K_ij to Zero!
  !======================================
  trK = ZERO
  Axx = ZERO
  Axy = ZERO
  Axz = ZERO
  Ayy = ZERO
  Ayz = ZERO
  Azz = ZERO
  !======================================
  ! Set everything else to Zero!
  !======================================
  shiftx = ZERO
  shifty = ZERO
  shiftz = ZERO
  !======================================
  ! Set initial EM fields to zero:
  !======================================
  Ex = By*vz - Bz*vy
  Ey = Bz*vx - Bx*vz
  Ez = Bx*vy - By*vx
  !!Ex = ZERO
  !!Ey = ZERO
  !!Ez = ZERO
  !FIXME: read note:
  ! NOTE: The following lines are not strictly correct, but this should NOT matter, since Cowling is turned on!
  Sxx = ZERO
  Sxy = ZERO
  Sxz = ZERO
  Syy = ZERO
  Syz = ZERO
  Szz = ZERO
  write(*,*) "checkpoint5"
  Sx = ZERO
  Sy = ZERO
  Sz = ZERO
  rho = ZERO
  write(*,*) "LOCAL INITIALDATA FINISHED!",dx,mhd_st_x(4,4,4), rho_star(4,4,4), h(4,4,4), vx(4,4,4), npolyshock, P(4,4,4),rho_b(4,4&
  &,4)
  write(*,*) "LOCAL INITIALDATA FINISHED! y=",y(1,:,1)
end subroutine mhd_shock_initialdata_local
