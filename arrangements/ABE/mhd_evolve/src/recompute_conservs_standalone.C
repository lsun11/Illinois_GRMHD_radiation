#include "cctk.h"
#include "harm_primitives_headers.h"


inline void recompute_conservatives_fast_standalone_singlept(double &rho_b,double &P,double &vx,double &vy,double &vz,
							     double &phi,double &lapm1,
							     double &shiftx,double &shifty,double &shiftz,
							     double &gxx,double &gxy,double &gxz,double &gyy,double &gyz,double &gzz,
							     double &gupxx,double &gupxy,double &gupxz,double &gupyy,double &gupyz,double &gupzz,
							     double &Bx,double &By,double &Bz, double &T_fluid,double &eps_tot, double &eps_cld, double &P_cld,
							     double &rho_star,double &tau,double &mhd_st_x,double &mhd_st_y,double &mhd_st_z,
							     int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,
							     double *k_tab,double *gamma_tab,double &rho_max, double &rho_b_atm, int &enable_OS_collapse,int &compute_microphysics);



extern "C" void CCTK_FCALL CCTK_FNAME(recompute_conservatives_fast_standalone_gf)
  (const cGH **cctkGH,int *ext,
   double *rho_b,double *P,double *vx,double *vy,double *vz,
   double *phi,double *lapm1,
   double *shiftx,double *shifty,double *shiftz,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
   double *Bx,double *By,double *Bz, double *T_fluid,double *eps_tot, double *eps_cld, double *P_cld,
   double *rho_star,double *tau,double *mhd_st_x,double *mhd_st_y,double *mhd_st_z,
   int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,
   double *k_tab,double *gamma_tab,double &rho_max, double &rho_b_atm, int &enable_OS_collapse, int &compute_microphysics);

extern "C" void recompute_conservatives_fast_standalone_gf(const cGH *cctkGH,int *ext,
							   double *rho_b,double *P,double *vx,double *vy,double *vz,
							   double *phi,double *lapm1,
							   double *shiftx,double *shifty,double *shiftz,
							   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
							   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
							   double *Bx,double *By,double *Bz, double *T_fluid,double *eps_tot, double *eps_cld, double *P_cld,
							   double *rho_star,double *tau,double *mhd_st_x,double *mhd_st_y,double *mhd_st_z,
							   int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,
							   double *k_tab,double *gamma_tab,double &rho_max, double &rho_b_atm, int &enable_OS_collapse, int &compute_microphysics) {
  
#pragma omp parallel for
  for(int k=0;k<ext[2];k++) for(int j=0;j<ext[1];j++) for(int i=0;i<ext[0];i++) {
	int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
	recompute_conservatives_fast_standalone_singlept(rho_b[index],P[index],vx[index],vy[index],vz[index],
							 phi[index],lapm1[index],
							 shiftx[index],shifty[index],shiftz[index],
							 gxx[index],gxy[index],gxz[index],gyy[index],gyz[index],gzz[index],
							 gupxx[index],gupxy[index],gupxz[index],gupyy[index],gupyz[index],gupzz[index],
							 Bx[index],By[index],Bz[index], T_fluid[index], eps_tot[index], eps_cld[index], P_cld[index],
							 rho_star[index],tau[index],mhd_st_x[index],mhd_st_y[index],mhd_st_z[index],
							 neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab,
							 k_tab,gamma_tab,rho_max,rho_b_atm, enable_OS_collapse, compute_microphysics);
      }
}

extern "C" void CCTK_FCALL CCTK_FNAME(recompute_conservatives_fast_standalone_gf)
  (const cGH **cctkGH,int *ext,
   double *rho_b,double *P,double *vx,double *vy,double *vz,
   double *phi,double *lapm1,
   double *shiftx,double *shifty,double *shiftz,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
   double *Bx,double *By,double *Bz, double *T_fluid,double *eps_tot, double *eps_cld, double *P_cld,
   double *rho_star,double *tau,double *mhd_st_x,double *mhd_st_y,double *mhd_st_z,
   int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,
   double *k_tab,double *gamma_tab,double &rho_max, double &rho_b_atm, int &enable_OS_collapse, int &compute_microphysics) 
{
  recompute_conservatives_fast_standalone_gf
    (*cctkGH,ext,
     rho_b,P,vx,vy,vz,
     phi,lapm1,
     shiftx,shifty,shiftz,
     gxx,gxy,gxz,gyy,gyz,gzz,
     gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,
     Bx,By,Bz,T_fluid,eps_tot,eps_cld,P_cld,
     rho_star,tau,mhd_st_x,mhd_st_y,mhd_st_z,
     neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab,
     k_tab,gamma_tab,rho_max,rho_b_atm, enable_OS_collapse, compute_microphysics);
}



inline void recompute_conservatives_fast_standalone_singlept(double &rho_b,double &P,double &vx,double &vy,double &vz,
						      double &phi,double &lapm1,
						      double &shiftx,double &shifty,double &shiftz,
						      double &gxx,double &gxy,double &gxz,double &gyy,double &gyz,double &gzz,
						      double &gupxx,double &gupxy,double &gupxz,double &gupyy,double &gupyz,double &gupzz,
						      double &Bx,double &By,double &Bz, double &T_fluid,double &eps_tot, double &eps_cld, double &P_cld,
        					      double &rho_star,double &tau,double &mhd_st_x,double &mhd_st_y,double &mhd_st_z,
	       					      int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,
						      double *k_tab,double *gamma_tab,double &rho_max, double &rho_b_atm, int &enable_OS_collapse, int &compute_microphysics) {

  void recompute_conservatives_fast(struct output_primitives &new_primitives,
				    int &enable_OS_collapse, int &compute_microphysics,
				    double &Psi6,double &alpha,
				    double &gxx_phys,double &gxy_phys,double &gxz_phys,double &gyy_phys,double &gyz_phys,double &gzz_phys,
				    double &gupxx_phys,double &gupxy_phys,double &gupxz_phys,double &gupyy_phys,double &gupyz_phys,double &gupzz_phys,
				    double &Bx,double &By,double &Bz, double &T_fluid,double &eps_tot, double &eps_cld, double &P_cld,
				    double &rho_star,double &tau,double &mhd_st_x,double &mhd_st_y,double &mhd_st_z,
				    double &u_xl,double &u_yl, double &u_zl,double &au0m1,
                                    int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,
  			            double *k_tab,double *gamma_tab,double &rho_max, double &rho_b_atm);

  struct output_primitives input_primitives;


  double alpha = lapm1 + 1.0;

  input_primitives.rho_b_new = rho_b;
  input_primitives.P_new = P;

  input_primitives.vx_new = vx;
  input_primitives.vy_new = vy;
  input_primitives.vz_new = vz;

  double Psi4 = exp(4.0*phi);
  double Psim4 = 1.0/Psi4;

  double gxx_physL = gxx*Psi4;
  double gxy_physL = gxy*Psi4;
  double gxz_physL = gxz*Psi4;
  double gyy_physL = gyy*Psi4;
  double gyz_physL = gyz*Psi4;
  double gzz_physL = gzz*Psi4;
  
  double gupxx_physL = gupxx*Psim4;
  double gupxy_physL = gupxy*Psim4;
  double gupxz_physL = gupxz*Psim4;
  double gupyy_physL = gupyy*Psim4;
  double gupyz_physL = gupyz*Psim4;
  double gupzz_physL = gupzz*Psim4;

  double alpha_inv = 1.0/(lapm1+1.0);

  //  printf("recompute conserv vatiable gxx_physL= %e gxy_physL= %e gxz_physL= %e  gyy_physL= %e gyz_physL= %e gzz_physL= %e\n ",
  //	 gxx_physL,gxy_physL,gxz_physL,gyy_physL,gyz_physL,gzz_physL);

  double v2 = (gxx_physL*(vx+shiftx)*(vx+shiftx) 
	       + 2.0*gxy_physL*(vx+shiftx)*(vy+shifty) 
	       + 2.0*gxz_physL*(vx+shiftx)*(vz+shiftz)
	       + gyy_physL*(vy+shifty)*(vy+shifty) 
	       + 2.0*gyz_physL*(vy+shifty)*(vz+shiftz) 
	       + gzz_physL*(vz+shiftz)*(vz+shiftz) )*alpha_inv*alpha_inv;

  // velocity limiter

  if(v2 >= 1.0){
    double max = 100;
    double A   =  (1.0 -  1.0/(max*max))/v2;
    vx         =  sqrt(A)*vx;
    vy         =  sqrt(A)*vy;
    vz         =  sqrt(A)*vz;
    v2         =  A*v2;
      }

  double alp_u02 = 1.0/(1.0-v2);
  double u02 = alp_u02*alpha_inv*alpha_inv;
  double u0 = sqrt(u02);

  input_primitives.u0_new = u0;

  //u0 = (au0m1+1.0)*alpha_inv implies:
  double au0m1 = u0*alpha - 1.0;

  double u_xl = u0*(gxx_physL*(vx+shiftx) + gxy_physL*(vy+shifty) + gxz_physL*(vz+shiftz) );
  double u_yl = u0*(gxy_physL*(vx+shiftx) + gyy_physL*(vy+shifty) + gyz_physL*(vz+shiftz) );
  double u_zl = u0*(gxz_physL*(vx+shiftx) + gyz_physL*(vy+shifty) + gzz_physL*(vz+shiftz) );

  double Psi6 = sqrt(Psi4)*Psi4;


  recompute_conservatives_fast(input_primitives,
			       enable_OS_collapse, compute_microphysics,
			       Psi6,alpha,
			       gxx_physL,gxy_physL,gxz_physL,gyy_physL,gyz_physL,gzz_physL,
			       gupxx_physL,gupxy_physL,gupxz_physL,gupyy_physL,gupyz_physL,gupzz_physL,
			       Bx,By,Bz, T_fluid,eps_tot,eps_cld,P_cld,
			       rho_star,tau,mhd_st_x,mhd_st_y,mhd_st_z,
			       u_xl,u_yl, u_zl,au0m1,
			       neos,ergo_star,ergo_sigma, rho_tab, P_tab, eps_tab,
			       k_tab,gamma_tab,rho_max,rho_b_atm);
}

//#include "harm_primitives_recompute_conservs_rho_star_tau_floor.C"
