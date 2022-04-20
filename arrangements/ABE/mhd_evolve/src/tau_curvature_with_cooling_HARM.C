//--------------------------------------------------------------------

// Add the extrinsic curvature terms to tau_rhs, as well as radiation source terms to tau_rhs and mhd_st_i_rhs, i=x,y,z
//--------------------------------------------------------------------

#define THIRD 0.333333333333333333333333333
#define SQR(x) ((x) * (x))

#include "math.h"
#include "cctk.h"


void compute_pcold_epscold_cpp(double &rhob, double &P_cold, double &eps_cold, 
			       int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab, int &enable_OS_collapse);

extern "C" void CCTK_FCALL CCTK_FNAME(mhd_tau_curvature_cpp_with_cooling_HARM)
  (const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry,
   double *tau_rhs,double *st_x_rhs,double *st_y_rhs,double *st_z_rhs,
   double *rho_star,double *h,double *P, double *rho_b,
   int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab, double *P_tab, double *eps_tab, double *k_tab, double *gamma_tab, double &gamma_th,
   double *sbt,double *sbx,double *sby,double *sbz,
   double *u0,double *vx,double *vy,double *vz,
   double *alpha,double *betax,double *betay,double *betaz,
   double *phi,double *trK,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
   double *Axx,double *Axy,double *Axz,double *Ayy,double *Ayz,double *Azz,
   double *dX,double *dY,double *dZ, double &t_cool, int &cooling_in_St_eq, int &allow_negative_eps_th,
   double *dt_al,double *dt_bx,double *dt_by,double *dt_bz,double *dt_phi,
   double *dt_gxx,double *dt_gxy,double *dt_gxz,double *dt_gyy,double *dt_gyz,double *dt_gzz, int *primitives_solver, int &enable_OS_collapse);

extern "C" void mhd_tau_curvature_cpp_with_cooling_HARM(const cGH *cctkGH,int *cctk_lsh, int *nghostzones, int Symmetry,
						   double *tau_rhs,double *st_x_rhs,double *st_y_rhs,double *st_z_rhs,
						   double *rho_star,double *h,double *P, double *rho_b,
						   int &neos, int &ergo_star, double &ergo_sigma, 
                                                   double *rho_tab, double *P_tab, double *eps_tab, double *k_tab, double *gamma_tab, double &gamma_th,
						   double *sbt,double *sbx,double *sby,double *sbz,
						   double *u0,double *vx,double *vy,double *vz,
						   double *alpha,double *betax,double *betay,double *betaz,
						   double *phi,double *trK,
						   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
						   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
						   double *Axx,double *Axy,double *Axz,double *Ayy,double *Ayz,double *Azz,
						   double dX,double dY,double dZ,double &t_cool, int &cooling_in_St_eq, int &allow_negative_eps_th,
						   double *dt_al,double *dt_bx,double *dt_by,double *dt_bz,double *dt_phi,
							double *dt_gxx,double *dt_gxy,double *dt_gxz,double *dt_gyy,double *dt_gyz,double *dt_gzz,int primitives_solver, int &enable_OS_collapse) {

  double sqrtmfourpi = 1.0/sqrt(4.0*M_PI);

  int AXISYM = 4;

  /* Set up variables used in the grid loop for the physical grid points */
  int istart = 0;
  int jstart = 0;
  int kstart = 0;
  int iend = cctk_lsh[0];
  int jend = cctk_lsh[1];
  int kend = cctk_lsh[2];
  /*
    int istart = nghostzones[0];
    int jstart = nghostzones[1];
    int kstart = nghostzones[2];
    int iend = cctk_lsh[0] - nghostzones[0];
    int jend = cctk_lsh[1] - nghostzones[1];
    int kend = cctk_lsh[2] - nghostzones[2];
  */

  /*
  //Following lines needed since nghostzones[0] = ORDER, and 
  //   not ORDER-1 in axisymmetry 
  //   (so that rotation can be done on multiprocessor runs)
  if(Symmetry==4) {
    istart--;
    iend++;
  }
  */

  if(Symmetry==4) {
    jstart = 1; jend = 2;
  }


  double eps, eps_cold, P_cold, eps_th;

#pragma omp parallel for
  for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    double alphaL = 1.0+alpha[index];
    double alphaLm1 = 1.0/alphaL;
    double Psi4 = exp(4.0*phi[index]);
    //    double Psi6 = exp(6.0*phi[index]);
    double Psim4 = 1.0/Psi4;     // Psi^{-4}
    double sqrtg = alphaL *exp(6.0*phi[index]); // Psi4 * sqrt(Psi4); // == alphaL*exp(6.0*phi)
    
    // Divide b^a by alpha sqrt(4 pi)
    double sbtL = sbt[index]*alphaLm1*sqrtmfourpi;
    double sbxL = sbx[index]*alphaLm1*sqrtmfourpi;
    double sbyL = sby[index]*alphaLm1*sqrtmfourpi;
    double sbzL = sbz[index]*alphaLm1*sqrtmfourpi;

    double betaxL = betax[index];
    double betayL = betay[index];
    double betazL = betaz[index];

    double gxxL = gxx[index];
    double gxyL = gxy[index];
    double gxzL = gxz[index];
    double gyyL = gyy[index];
    double gyzL = gyz[index];
    double gzzL = gzz[index];

    double rho_starL = rho_star[index];
    double u0L = u0[index];
    double hL = h[index];
    double PL = P[index];

    double rho_bL = rho_b[index];

    double vxL = vx[index];
    double vyL = vy[index];
    double vzL = vz[index];


    // Compute eps_cold and P_cold
    
    compute_pcold_epscold_cpp(rho_bL, P_cold, eps_cold,neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab, enable_OS_collapse);
    
    // Use total pressure P_l and cold pressure P_cold to compute the thermal specific energy
    if (enable_OS_collapse == 1){
    eps_th = 0.0;
    }else{
    eps_th = (PL - P_cold)/(gamma_th-1.0)/rho_bL;
    }

    if(allow_negative_eps_th==0 && eps_th < 0.) eps_th = 0.;
    

    // Compute b^2 
    double b2 = -SQR(alphaL*sbtL) + Psi4*( gxxL*SQR(sbxL+betaxL*sbtL) + 
					   2.0*gxyL*(sbxL+betaxL*sbtL)*(sbyL+betayL*sbtL) + 
					   2.0*gxzL*(sbxL+betaxL*sbtL)*(sbzL+betazL*sbtL) + 
					   gyyL*SQR(sbyL+betayL*sbtL) + 
					   2.0*gyzL*(sbyL+betayL*sbtL)*(sbzL+betazL*sbtL) + 
					   gzzL*SQR(sbzL+betazL*sbtL) );


    // Calculate \beta_i and beta^2 = beta_i*beta^i
    double beta_xL = Psi4*(gxxL*betaxL + gxyL*betayL + gxzL*betazL );    
    double beta_yL = Psi4*(gxyL*betaxL + gyyL*betayL + gyzL*betazL );
    double beta_zL = Psi4*(gxzL*betaxL + gyzL*betayL + gzzL*betazL );
    double beta2 = beta_xL*betaxL + beta_yL*betayL + beta_zL*betazL;

    double gupxxL = gupxx[index];
    double gupxyL = gupxy[index];
    double gupxzL = gupxz[index];
    double gupyyL = gupyy[index];
    double gupyzL = gupyz[index];
    double gupzzL = gupzz[index];


    double dt_alL = dt_al[index];
    double dt_bxL = dt_bx[index];
    double dt_byL = dt_by[index];
    double dt_bzL = dt_bz[index];
    double dt_phiL = dt_phi[index];
    double dt_gxxL = dt_gxx[index];
    double dt_gxyL = dt_gxy[index];
    double dt_gxzL = dt_gxz[index];
    double dt_gyyL = dt_gyy[index];
    double dt_gyzL = dt_gyz[index];
    double dt_gzzL = dt_gzz[index];

    // Calculate time derivatives of the 4-metric

    double dt_g4tt = -2.0*alphaL*dt_alL + 2.0*(beta_xL*dt_bxL + beta_yL*dt_byL + beta_zL*dt_bzL) + 
                     4.0*beta2*dt_phiL + Psi4*(
					       betaxL*betaxL*dt_gxxL + 2.0*betaxL*betayL*dt_gxyL +
					       2.0*betaxL*betazL*dt_gxzL + betayL*betayL*dt_gyyL +
					       2.0*betayL*betazL*dt_gyzL + betazL*betazL*dt_gzzL
					       );
    
    double dt_g4tx = 4.0*beta_xL*dt_phiL + Psi4*( betaxL*dt_gxxL + betayL*dt_gxyL + betazL*dt_gxzL) + 
                     Psi4*(gxxL*dt_bxL + gxyL*dt_byL + gxzL*dt_bzL);

    double dt_g4ty = 4.0*beta_yL*dt_phiL + Psi4*( betaxL*dt_gxyL + betayL*dt_gyyL + betazL*dt_gyzL) + 
                     Psi4*(gxyL*dt_bxL + gyyL*dt_byL + gyzL*dt_bzL);

    double dt_g4tz = 4.0*beta_zL*dt_phiL + Psi4*( betaxL*dt_gxzL + betayL*dt_gyzL + betazL*dt_gzzL) + 
                     Psi4*(gxzL*dt_bxL + gyzL*dt_byL + gzzL*dt_bzL);

    double dt_g4xx = Psi4*(4.0*gxxL*dt_phiL + dt_gxxL);
    double dt_g4xy = Psi4*(4.0*gxyL*dt_phiL + dt_gxyL);
    double dt_g4xz = Psi4*(4.0*gxzL*dt_phiL + dt_gxzL);
    double dt_g4yy = Psi4*(4.0*gyyL*dt_phiL + dt_gyyL);
    double dt_g4yz = Psi4*(4.0*gyzL*dt_phiL + dt_gyzL);
    double dt_g4zz = Psi4*(4.0*gzzL*dt_phiL + dt_gzzL);



    //    calculate sqrt{-g}T^\mu\nu
    double Tuptt = rho_starL*hL*u0L + sqrtg * ( b2*SQR(u0L) - (PL + 0.5*b2)/(SQR(alphaL))-SQR(sbtL) );
    double Tuptx = rho_starL*hL*u0L*vxL + sqrtg * ( b2*SQR(u0L)*vxL + (PL + 0.5*b2)*betaxL/(SQR(alphaL))-sbtL*sbxL );
    double Tupty = rho_starL*hL*u0L*vyL + sqrtg * ( b2*SQR(u0L)*vyL + (PL + 0.5*b2)*betayL/(SQR(alphaL))-sbtL*sbyL );
    double Tuptz = rho_starL*hL*u0L*vzL + sqrtg * ( b2*SQR(u0L)*vzL + (PL + 0.5*b2)*betazL/(SQR(alphaL))-sbtL*sbzL );
    double Tupxx = rho_starL*hL*u0L*vxL*vxL + 
      sqrtg * ( b2*SQR(u0L)*vxL*vxL + (PL + 0.5*b2)*(Psim4*gupxxL - betaxL*betaxL/(SQR(alphaL)))-sbxL*sbxL );
    double Tupxy = rho_starL*hL*u0L*vxL*vyL + 
      sqrtg * ( b2*SQR(u0L)*vxL*vyL + (PL + 0.5*b2)*(Psim4*gupxyL - betaxL*betayL/(SQR(alphaL)))-sbxL*sbyL );
    double Tupxz = rho_starL*hL*u0L*vxL*vzL + 
      sqrtg * ( b2*SQR(u0L)*vxL*vzL + (PL + 0.5*b2)*(Psim4*gupxzL - betaxL*betazL/(SQR(alphaL)))-sbxL*sbzL );
    double Tupyy = rho_starL*hL*u0L*vyL*vyL + 
      sqrtg * ( b2*SQR(u0L)*vyL*vyL + (PL + 0.5*b2)*(Psim4*gupyyL - betayL*betayL/(SQR(alphaL)))-sbyL*sbyL );
    double Tupyz = rho_starL*hL*u0L*vyL*vzL + 
      sqrtg * ( b2*SQR(u0L)*vyL*vzL + (PL + 0.5*b2)*(Psim4*gupyzL - betayL*betazL/(SQR(alphaL)))-sbyL*sbzL );
    double Tupzz = rho_starL*hL*u0L*vzL*vzL + 
      sqrtg * ( b2*SQR(u0L)*vzL*vzL + (PL + 0.5*b2)*(Psim4*gupzzL - betazL*betazL/(SQR(alphaL)))-sbzL*sbzL );
      
      
      
    // Now add time derivative source terms -0.5*sqrt{-g}T^{\mu\nu}\partial_t g_{\mu\nu} to tau_rhs: 
    double tau_rhsL = -(Tuptt*dt_g4tt + 2.0*Tuptx*dt_g4tx + 2.0*Tupty*dt_g4ty + 2.0*Tuptz*dt_g4tz + 
                        Tupxx*dt_g4xx + 2.0*Tupxy*dt_g4xy + 2.0*Tupxz*dt_g4xz + 
                        Tupyy*dt_g4yy + 2.0*Tupyz*dt_g4yz + Tupzz*dt_g4zz); 

    tau_rhs[index] += tau_rhsL*0.5;



    ///////////////////////////////////////////////////////////////////////////
    // Cooling Source terms for exponential cooling of eps_th
    ///////////////////////////////////////////////////////////////////////////

    // Compute cooling (i.e. integrated emissivity)
    double Gamma = rho_bL*eps_th/t_cool;
        
    // Add Cooling source to tau_rhs
    tau_rhs[index] += (-alphaL*sqrtg*u0L*Gamma);



    //////////////////////////////////////////////////////////////////////////////////
    // If the primitives_solver is 3 (i.e., isothermal the energy variable 
    // does not matter. In case a random evolution crashes the code we set tau_rhs =0
    //////////////////////////////////////////////////////////////////////////////////

    if (primitives_solver==3) {
      tau_rhs[index] = 0.0;
    }

    //printf("tau_rhs=%e, Gamma=%e", (-alphaL*alphaL*Psi6*u0L*Gamma), Gamma);

    
   // Add Cooling source to st_x_rhs, st_y_rhs,st_z_rhs, if cooling_in_St_eq = 1

    // Compute u_x, u_y, u_z low indices
    double u_xll, u_yll, u_zll;

      if (cooling_in_St_eq==1){
	u_xll = (gxxL*(betaxL+vx[index]) + 
		 gxyL*(betayL+vy[index]) + 
		 gxzL*(betazL+vz[index]))*Psi4*u0L;
	u_yll = (gxyL*(betaxL+vx[index]) + 
		 gyyL*(betayL+vy[index]) + 
		 gyzL*(betazL+vz[index]))*Psi4*u0L;
	u_zll = (gxzL*(betaxL+vx[index]) + 
		 gyzL*(betayL+vy[index]) + 
		 gzzL*(betazL+vz[index]))*Psi4*u0L;

	st_x_rhs[index] += (-sqrtg*u_xll*Gamma);
	st_y_rhs[index] += (-sqrtg*u_yll*Gamma);
	st_z_rhs[index] += (-sqrtg*u_zll*Gamma);
      }
  }
}



extern "C" void CCTK_FCALL CCTK_FNAME(mhd_tau_curvature_cpp_with_cooling_HARM)
  (const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry,
   double *tau_rhs,double *st_x_rhs,double *st_y_rhs,double *st_z_rhs,
   double *rho_star,double *h,double *P, double *rho_b,
   int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab, double *P_tab, double *eps_tab, double *k_tab, double *gamma_tab, double &gamma_th,
   double *sbt,double *sbx,double *sby,double *sbz,
   double *u0,double *vx,double *vy,double *vz,
   double *alpha,double *betax,double *betay,double *betaz,
   double *phi,double *trK,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
   double *Axx,double *Axy,double *Axz,double *Ayy,double *Ayz,double *Azz,
   double *dX,double *dY,double *dZ, double &t_cool, int &cooling_in_St_eq, int &allow_negative_eps_th,
   double *dt_al,double *dt_bx,double *dt_by,double *dt_bz,double *dt_phi,
   double *dt_gxx,double *dt_gxy,double *dt_gxz,double *dt_gyy,double *dt_gyz,double *dt_gzz,int *primitives_solver, int &enable_OS_collapse)
{
  mhd_tau_curvature_cpp_with_cooling_HARM(*cctkGH,cctk_lsh, nghostzones, *Symmetry,
                        tau_rhs,st_x_rhs,st_y_rhs,st_z_rhs,rho_star,h,P,rho_b,
			neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th,
			sbt,sbx,sby,sbz,
			u0,vx,vy,vz,
			alpha,betax,betay,betaz,
			phi,trK,
			gxx,gxy,gxz,gyy,gyz,gzz,
			gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,
			Axx,Axy,Axz,Ayy,Ayz,Azz,
			*dX,*dY,*dZ,t_cool,cooling_in_St_eq,allow_negative_eps_th,
			dt_al,dt_bx,dt_by,dt_bz,dt_phi,
			dt_gxx,dt_gxy,dt_gxz,dt_gyy,dt_gyz,dt_gzz, *primitives_solver, enable_OS_collapse);
}
