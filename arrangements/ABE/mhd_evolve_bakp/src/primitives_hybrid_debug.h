#define PRINT_SHOCKTUBE_DEBUG printf("Primitive solver fails at (i,j,k): %d %d %d\n",i,j,k); \
  printf("(x,y,z): %e %e %e\n",X[index],Y[index],Z[index]);		\
  printf("rho_star: %e\n",rho_star[index]);				\
  printf("tau: %e\n",tau[index]);					\
  printf("mhd_st_x: %e\n",mhd_st_x[index]);				\
  printf("mhd_st_y: %e\n",mhd_st_y[index]);				\
  printf("mhd_st_z: %e\n",mhd_st_z[index]);				\
  printf("Bx: %e\n",Bx[index]);						\
  printf("By: %e\n",By[index]);						\
  printf("Bz: %e\n",Bz[index]);						\
  printf("B2: %e\n",AUX.B2);						\
  printf("gamma^xx: %e\n",AUX.gupxx_phys);				\
  printf("gamma^xy: %e\n",AUX.gupxy_phys);				\
  printf("gamma^xz: %e\n",AUX.gupxz_phys);				\
  printf("gamma^yy: %e\n",AUX.gupyy_phys);				\
  printf("gamma^yz: %e\n",AUX.gupyz_phys);				\
  printf("gamma^zz: %e\n",AUX.gupzz_phys);				\
  printf("exp(6 phi): %e\n",AUX.Psi6);

#define PRINT_TAU_STILDE_DEBUG \
  printf("ERROR. TAU STILDE FIX GUARANTEES A SOLUTION (analytically, at least).  THERE WAS NO SOLUTION FOUND.\n"); \
  printf("To read more about the tau stilde fix, check out: arXiv:0708.2436v3, or Phys.Rev.D76:104021,2007 (published version)\n"); \
  printf("Remember that you can use the tau stilde fix only if B=0.\n"); \
  printf("SOLUTION: USE THE QUARTIC PRIMITIVES SOLVER INSTEAD! (set mhd_evolve::primitives_solver=1 in your .par file)!\n\n"); \
  printf("coords of bad point: %d %d %d\n",i,j,k);			\
  printf("AFTER: UU = %e %e %e %e\n",UU[1],UU[2],UU[3],UU[4]);		\
									\
  if(check==1) printf("\n\n\nNote: the reason the primitives solver failed is because it hit the max count on iterations\n\n\n"); \
									\
  printf("OUTPUT UU[1] %.16e, check %d\n",UU[1],check);			\
  printf("OUTPUT UU[2] %.16e, check %d\n",UU[2],check);			\
  printf("OUTPUT UU[3] %.16e, check %d\n",UU[3],check);			\
  printf("OUTPUT UU[4] %.16e, check %d\n",UU[4],check);			\
									\
  printf("UUguess[1]: %.16e\n",UUguess[1]);				\
  printf("UUguess[2]: %.16e\n",UUguess[2]);				\
  printf("UUguess[3]: %.16e\n",UUguess[3]);				\
  printf("UUguess[4]: %.16e\n",UUguess[4]);				\
									\
  printf("rhobatm: %.16e\n",rhobatm);					\
  printf("tau_atm: %.16e\n",tau_atm);					\
									\
  printf("lapse: %.16e\n",alpn1);					\
  printf("shiftx: %.16e\n",shiftxL);					\
  printf("shifty: %.16e\n",shiftyL);					\
  printf("shiftz: %.16e\n",shiftzL);					\
  printf("tau: %.16e\n",tau[index]);					\
  printf("rho_star: %.16e\n",rho_star[index]);				\
  printf("mhd_st_x: %.16e\n",mhd_st_x[index]);				\
  printf("mhd_st_y: %.16e\n",mhd_st_y[index]);				\
  printf("mhd_st_z: %.16e\n",mhd_st_z[index]);				\
  printf("gamma^xx_phys: %.16e\n",AUX.gupxx_phys);			\
  printf("gamma^xy_phys: %.16e\n",AUX.gupxy_phys);			\
  printf("gamma^xz_phys: %.16e\n",AUX.gupxz_phys);			\
  printf("gamma^yy_phys: %.16e\n",AUX.gupyy_phys);			\
  printf("gamma^yz_phys: %.16e\n",AUX.gupyz_phys);			\
  printf("gamma^zz_phys: %.16e\n",AUX.gupzz_phys);			\
									\
  printf("gamma^xx: %.16e\n",gupxx[index]);				\
  printf("gamma^xy: %.16e\n",gupxy[index]);				\
  printf("gamma^xz: %.16e\n",gupxz[index]);				\
  printf("gamma^yy: %.16e\n",gupyy[index]);				\
  printf("gamma^yz: %.16e\n",gupyz[index]);				\
  printf("gamma^zz: %.16e\n",gupzz[index]);				\
									\
  printf("exp(6 phi): %.16e\n",AUX.Psi6);				\
  printf("phi: %.16e\n",phi[index]);					\
									\
  printf("u_scal: %.16e\n",u_scal);					\
  printf("eps_scal: %.16e\n",eps_scal);					\
  printf("sti_scal: %.16e\n",sti_scal);					\
  printf("tau_scal: %.16e\n",tau_scal);					\
									\
  printf("gxx: %.16e\n",gxxL);						\
  printf("gxy: %.16e\n",gxyL);						\
  printf("gxz: %.16e\n",gxzL);						\
  printf("gyy: %.16e\n",gyyL);						\
  printf("gyz: %.16e\n",gyzL);						\
  printf("gzz: %.16e\n",gzzL);						\
  printf("P_l: %.16e\n",h_l);						\
  printf("h_l: %.16e\n",P_l);						\
									\
  printf("\n\n VARIABLES THAT SHOULD BE ZERO:\n");			\
  printf("Bx: %e\n",Bx[index]);						\
  printf("By: %e\n",By[index]);						\
  printf("Bz: %e\n",Bz[index]);						\
  printf("B2: %e\n",AUX.B2);

#define PRINT_FONT_FIX_DEBUG printf("ERROR: FONT FIX (secondary solver) JUST FAILED\n"); \
	         printf("Problem at (x,y,z) = %e %e %e, %d %d %d\n", X[index],Y[index],Z[index],i,j,k); \
	         printf("rho_star = %e\n",rho_star[index]); \
	         printf("tau = %e\n",tau[index]); \
	         printf("mhd_st_x = %e\n",mhd_st_x[index]); \
	         printf("mhd_st_y = %e\n",mhd_st_y[index]); \
	         printf("mhd_st_z = %e\n",mhd_st_z[index]); \
	         printf("Bx = %e\n",BxL); \
	         printf("By = %e\n",ByL); \
	         printf("Bz = %e\n",BzL); \
	         printf("B2 = %e\n",AUX.B2); \
	         printf("gamma^xx = %e\n",AUX.gupxx_phys); \
	         printf("gamma^xy = %e\n",AUX.gupxy_phys); \
	         printf("gamma^xz = %e\n",AUX.gupxz_phys); \
	         printf("gamma^yy = %e\n",AUX.gupyy_phys); \
	         printf("gamma^yz = %e\n",AUX.gupyz_phys); \
	         printf("gamma^zz = %e\n",AUX.gupzz_phys); \
	         printf("exp(6 phi) = %e\n",AUX.Psi6);
                       
