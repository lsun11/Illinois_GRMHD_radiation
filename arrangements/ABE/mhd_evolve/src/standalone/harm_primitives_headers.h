#include "harm_u2p_defs.h"
#include "harm_u2p_util.h"



struct output_primitives {
  double rho_b_new,P_new,vx_new,vy_new,vz_new,u0_new;
};

struct output_stats {
  int font_fixed;
};



/* these variables need to be shared between the functions
   Utoprim_1D, residual, and utsq */
//FTYPE Bsq,QdotBsq,Qtsq,Qdotn,D ;


/********************************************************************************************/
// Function prototype declarations:



int Utoprim_2d(FTYPE U[NPR], FTYPE gcov[NDIM][NDIM], FTYPE gcon[NDIM][NDIM], 
	       FTYPE gdet, FTYPE prim[NPR], 
	       int &neos, int &ergo_star, FTYPE &ergo_sigma, FTYPE &T_fluid, FTYPE *rho_tab, FTYPE *P_tab,  FTYPE *eps_tab,  FTYPE *k_tab,  FTYPE *gamma_tab, int &enable_OS_collapse, int &compute_microphysics);



/*
void recompute_conservatives_fast(struct output_primitives &new_primitives,
				  double &phi,double &alpha,double &shiftx,double &shifty,double &shiftz,
				  double &gxx,double &gxy,double &gxz,double &gyy,double &gyz,double &gzz,
				  double &gupxx,double &gupxy,double &gupxz,double &gupyy,double &gupyz,double &gupzz,
				  double &Bx,double &By,double &Bz,
				  double &rho_star,double &tau,double &mhd_st_x,double &mhd_st_y,double &mhd_st_z);
*/

void recompute_conservatives_fast(struct output_primitives &new_primitives,
				  int &enable_OS_collapse, int &compute_microphysics,
				  double &Psi6,double &alpha,
				  double &gxx_phys,double &gxy_phys,double &gxz_phys,double &gyy_phys,double &gyz_phys,double &gzz_phys,
				  double &gupxx_phys,double &gupxy_phys,double &gupxz_phys,double &gupyy_phys,double &gupyz_phys,double &gupzz_phys,
				  double &Bx,double &By,double &Bz, double &T_fluid, double &eps_tot, double &eps_cld, double &P_cld,
				  double &rho_star,double &tau,double &mhd_st_x,double &mhd_st_y,double &mhd_st_z,
				  double &u_xl,double &u_yl, double &u_zl,double &au0m1,
				  int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,
				  double *k_tab,double *gamma_tab,double &rho_max, double &rho_b_atm);



int apply_tau_floor(int &index,double &tau_atm,double &rhobatm,
		    double *Bx,double *By,double *Bz,
		    double *tau,double *rho_star,double *mhd_st_x,double *mhd_st_y,double *mhd_st_z,
		    double *phi,double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
		    double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
		    struct output_primitives &new_primitives,
		    double *lapm1,double *shiftx,double *shifty,double *shiftz,
		    int &neos,int &ergo_star, double &ergo_sigma,
		    double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab,
		    double &Psi6threshold);


int harm_primitives_gammalaw_lowlevel(int &index,double *X,double *Y,double *Z,
				      int &enable_OS_collapse,
				      double *phi,double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
				      double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
				      double *lapm1,double *shiftx,double *shifty,double *shiftz,
				      double *Bx,double *By,double *Bz,
				      double *T_fluid,				      
				      double *mhd_st_x,double *mhd_st_y,double *mhd_st_z,double *tau,double *rho_star,
				      double *vx,double *vy,double *vz,double *P,double *rho_b,double *h,double *u0,
				      double *ka_gf, double *ks_gf, double *emission_gf, double *chi_rad, double *chi_rad_nue, double *Y_e, double *optd,
				      double *ka_gf_nue, double *ks_gf_nue, double *emission_gf_nue, double *ka_gf_nux, double *ks_gf_nux, double *emission_gf_nux,
				      double *eps_tot, double *eps_thermal, double *eps_cld, double *P_cld,
				      double &rhobatm,double &tau_atm,
				      int &neos,int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab,
				      double &rho_max, 
				      struct output_primitives &new_primitives, 
				      struct output_stats &stats,double &Psi6threshold,
				      int &horizon_enforce_rho_profile, int &compute_mircophysics, int &microphysics_scheme, double &T_fluid_cgs_atm, int &rad_fix);




int font_fix(double *UU_font_fix,
             double &rho_starL,double &tauL,double &mhd_st_xL,double &mhd_st_yL,double &mhd_st_zL,
             double &gamma_th,double &BxL,double &ByL,double &BzL,
             double &alphaL,double &gxx_physL,double &gxy_physL,double &gxz_physL,double &gyy_physL,double &gyz_physL,double &gzz_physL,
             double &gupxx_physL,double &gupxy_physL,double &gupxz_physL,double &gupyy_physL,double &gupyz_physL,double &gupzz_physL,
             double &Psi2,double &Psi6,
             int &neos,double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab);




int font_fix_gamma_equals2(double &u_x, double &u_y, double &u_z, double &rhob,
             double &rho_starL,double &mhd_st_xL,double &mhd_st_yL,double &mhd_st_zL,
             double &BxL,double &ByL,double &BzL,
             double &gxx_physL,double &gxy_physL,double &gxz_physL,double &gyy_physL,double &gyz_physL,double &gzz_physL,
             double &gupxx_physL,double &gupxy_physL,double &gupxz_physL,double &gupyy_physL,double &gupyz_physL,double &gupzz_physL,
			   double &Psi6, double &kpoly,int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab, int &enable_OS_collapse);



int font_fix_general_gamma(double &u_x, double &u_y, double &u_z, double &rhob,
             double &rho_starL,double &mhd_st_xL,double &mhd_st_yL,double &mhd_st_zL,
             double &BxL,double &ByL,double &BzL,
             double &gxx_physL,double &gxy_physL,double &gxz_physL,double &gyy_physL,double &gyz_physL,double &gzz_physL,
             double &gupxx_physL,double &gupxy_physL,double &gupxz_physL,double &gupyy_physL,double &gupyz_physL,double &gupzz_physL,
			   double &Psi6, double &gamma, double &kpoly,int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab, int &enable_OS_collapse);



void eigenvalues_3by3_real_sym_matrix(double & lam1, double & lam2, double & lam3,
         double M11, double M12, double M13, double M22, double M23, double M33);











/********************************************************************************************/
