           call	primitive_vars_hybrid2_cpp(ext,cctk_nghostzones,X,Y,Z, &
		rho_star,tau,st_x,st_y,st_z, &
		mhd_st_x,mhd_st_y,mhd_st_z,neos, ergo_star,ergo_sigma,&
		rho_tab, P_tab, eps_tab, k_tab,	gamma_tab, gamma_th, &
		w,temp1,rho_b,rho,P,h,Sx,Sy,Sz, &
		Sxx,Sxy,Sxz,Syy,Syz,Szz, &
		phi,lapm1,shiftx,shifty,shiftz, &
		gxx,gxy,gxz,gyy,gyz,gzz, &
		gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
		h_p,u0,rho_b_max,rho_b_atm, &
		rho_fail_max_step,M_fail_step,rhos_max, &
		Bx,By,Bz,Ex,Ey,Ez, &
		vx,vy,vz, &
                sbt,sbx,sby,sbz, &
		xNS1,yNS1, xNS2, yNS2, num_CO, M_B, rad_T_fac, rad_T_cutoff, rad_T_pow, rad_T_floor, T_fluid,&
		proc_imin,proc_jmin,proc_kmin, &
	       	proc_imax,proc_jmax,proc_kmax, &
		glob_imax,glob_jmax,glob_kmax, &
	 	Symmetry,pfloor,excision_enable, &
                excision_zone_gf, tau_stildefix_enable,tau_atm,temp4,cctkGH,ignore_ghostzones, &
		enable_shocktest_primitive_mode,repairs_needed, rad_closure_scheme, rad_const)

