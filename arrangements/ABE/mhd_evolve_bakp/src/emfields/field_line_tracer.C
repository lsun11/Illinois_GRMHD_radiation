//--------------------------------------------------------------------
// Add the extrinsic curvature terms to tau_rhs 
//--------------------------------------------------------------------

#define SQR(x) ((x) * (x))

#include "math.h"
#include "stdio.h"
#include "cctk.h"

extern "C" void CCTK_FCALL CCTK_FNAME(field_line_tracer_cpp)
  (const cGH **cctkGH,int *cctk_lsh, int *nghostzones, 
   double *alpm1, double *shiftx, double *shifty, double *shiftz, 
   double *phi,double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *vx,double *vy,double *vz, double *Bxtilde,double *Bytilde,double *Bztilde, 
   double *mhd_psi_line, double *mhd_u_psi, double *mhd_chi_line, double *mhd_u_chi, 
   double *mhd_psi_line_rhs, double *mhd_u_psi_rhs, double *mhd_chi_line_rhs, double *mhd_u_chi_rhs, 
   double *lambda_line, double *dX, double *dY, double *dZ, 
   double *X, double *Y, double *Z, double *r_ah, 
   double *x_ah, double *y_ah, double *z_ah);

extern "C" void field_line_tracer_cpp(const cGH *cctkGH,int *cctk_lsh, int *nghostzones,
				      double *alpm1, double *shiftx, double *shifty, double *shiftz,
                                      double *phi,double *gxx,double *gxy,double *gxz,double *gyy,
				      double *gyz,double *gzz, double *vx,double *vy,double *vz,
				      double *Bxtilde,double *Bytilde,double *Bztilde, double *mhd_psi_line, 
				      double *mhd_u_psi, double *mhd_chi_line, double *mhd_u_chi,
				      double *mhd_psi_line_rhs, double *mhd_u_psi_rhs,
				      double *mhd_chi_line_rhs, double *mhd_u_chi_rhs, 
				      double lambda_line, double dX,double dY,double dZ, 
				      double *X, double *Y, double *Z, double r_ah,
				      double x_ah, double y_ah, double z_ah) {

  printf("IF YOU WANT TO USE THE FIELD LINE TRACER, YOU MUST RE-ENABLE THE SYNC: field_line... stuff in mhd_evolve/schedule.ccl!\n");

  /* Set up variables used in the grid loop for the physical grid points */
  int istart = nghostzones[0];
  int jstart = nghostzones[1];
  int kstart = nghostzones[2];
  int iend = cctk_lsh[0]-nghostzones[0];
  int jend = cctk_lsh[1]-nghostzones[1];
  int kend = cctk_lsh[2]-nghostzones[2];     
  
  double hdXm1 = 0.5/dX;
  double hdYm1 = 0.5/dY;
  double hdZm1 = 0.5/dZ;

#pragma omp parallel for
  for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    double alphaL = 1.0+alpm1[index];
    double alphaLm1 = 1.0/alphaL;
    double Psi4 = exp(4.0*phi[index]);
    double betaxL = shiftx[index];
    double betayL = shifty[index];
    double betazL = shiftz[index];
    double gxxL = Psi4*gxx[index];
    double gxyL = Psi4*gxy[index];
    double gxzL = Psi4*gxz[index];
    double gyyL = Psi4*gyy[index];
    double gyzL = Psi4*gyz[index];
    double gzzL = Psi4*gzz[index];
    double BxtL = Bxtilde[index];
    double BytL = Bytilde[index];
    double BztL = Bztilde[index];
    double vxL = vx[index];
    double vyL = vy[index];
    double vzL = vz[index];
    double cVxL = (vxL + betaxL)*alphaLm1;
    double cVyL = (vyL + betayL)*alphaLm1;
    double cVzL = (vzL + betazL)*alphaLm1;
    double V2 = gxxL*SQR(cVxL) + 2.0*gxyL*cVxL*cVyL + 2.0*gxzL*cVxL*cVzL
		+ gyyL*SQR(cVyL) + 2.0*gyzL*cVyL*cVzL + gzzL*SQR(cVzL);
    // Check for superluminal velocity
    double V2_reset = 0.99;
    if (V2>1.0) {
       double fac = sqrt(V2_reset/V2);
       V2 = V2_reset;
       cVxL = fac*cVxL;
       cVyL = fac*cVyL;
       cVzL = fac*cVzL;
       vxL = alphaL*cVxL - betaxL;
       vyL = alphaL*cVyL - betayL;
       vzL = alphaL*cVzL - betazL;
    }
    double V_plus_beta = alphaL*sqrt(V2);

    double beta = sqrt( fabs( gxxL*SQR(betaxL) + 2.0*gxyL*betaxL*betayL + 2.0*gxzL*betaxL*betazL
                + gyyL*SQR(betayL) + 2.0*gyzL*betayL*betazL + gzzL*SQR(betazL) ) );
    
    double Bt2 = gxxL*SQR(BxtL) + 2.0*gxyL*BxtL*BytL + 2.0*gxzL*BxtL*BztL 
		+ gyyL*SQR(BytL) + 2.0*gyzL*BytL*BztL + gzzL*SQR(BztL);

    double mu_Btx, mu_Bty, mu_Btz;
    double hatBxt, hatByt, hatBzt;
    if (Bt2 > 0.0) {
       double Bt = sqrt(Bt2);
       hatBxt = BxtL/Bt; 
       hatByt = BytL/Bt; 
       hatBzt = BztL/Bt;
       double tmp1 = (alphaL - beta)*(alphaL-V_plus_beta); 
       double common_fact;
       if (tmp1 > 0) {
	  common_fact = sqrt(tmp1);
       } else {
	  common_fact = 0.0;
       }
       mu_Btx = common_fact*hatBxt;
       mu_Bty = common_fact*hatByt;
       mu_Btz = common_fact*hatBzt;
    } else {
      hatBxt=0.0; hatByt=0.0; hatBzt=0.0;
      mu_Btx=0.0; mu_Bty=0.0; mu_Btz=0.0;
    }

    int upwind_enable = 0; 
    if (SQR(mu_Btx) + SQR(mu_Bty) + SQR(mu_Btz) == 0 ) upwind_enable = 1;

    // Compute derivatives of mhd_psi_line, mhd_u_psi, mhd_chi_line, and mhd_u_chi
    int index1,indexm1;
    double dx_mhd_psi_line, dx_mhd_u_psi, dx_mhd_chi_line, dx_mhd_u_chi;
    if (upwind_enable==0) {
       indexm1 = CCTK_GFINDEX3D(cctkGH,i-1,j,k);
       index1  = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
       dx_mhd_psi_line = hdXm1*(mhd_psi_line[index1] - mhd_psi_line[indexm1]);
       dx_mhd_chi_line = hdXm1*(mhd_chi_line[index1] - mhd_chi_line[indexm1]);
       dx_mhd_u_psi    = hdXm1*(mhd_u_psi[index1] - mhd_u_psi[indexm1]);
       dx_mhd_u_chi    = hdXm1*(mhd_u_chi[index1] - mhd_u_chi[indexm1]);
    }


    int index2, indexm2;
    if (upwind_enable==1) {
       if (vxL > 0) {
          indexm1 = CCTK_GFINDEX3D(cctkGH,i-1,j,k);
          indexm2 = CCTK_GFINDEX3D(cctkGH,i-2,j,k);
          double cm2 = 0.5/dX;
          double cm1 = -4.0*cm2;
          double c0 = 3.0*cm2;
          dx_mhd_psi_line = c0*mhd_psi_line[index] + cm1*mhd_psi_line[indexm1] + cm2*mhd_psi_line[indexm2];
          dx_mhd_u_psi = c0*mhd_u_psi[index] + cm1*mhd_u_psi[indexm1] + cm2*mhd_u_psi[indexm2];
          dx_mhd_chi_line = c0*mhd_chi_line[index] + cm1*mhd_chi_line[indexm1] + cm2*mhd_chi_line[indexm2];
          dx_mhd_u_chi = c0*mhd_u_chi[index] + cm1*mhd_u_chi[indexm1] + cm2*mhd_u_chi[indexm2];
       } else {
          index1 = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
          index2 = CCTK_GFINDEX3D(cctkGH,i+2,j,k);
          double c2 = -0.5/dX;
          double c1 = -4.0*c2;
          double c0 = 3.0*c2;
          dx_mhd_psi_line = c0*mhd_psi_line[index] + c1*mhd_psi_line[index1] + c2*mhd_psi_line[index2];
          dx_mhd_u_psi = c0*mhd_u_psi[index] + c1*mhd_u_psi[index1] + c2*mhd_u_psi[index2];
          dx_mhd_chi_line = c0*mhd_chi_line[index] + c1*mhd_chi_line[index1] + c2*mhd_chi_line[index2];
          dx_mhd_u_chi = c0*mhd_u_chi[index] + c1*mhd_u_chi[index1] + c2*mhd_u_chi[index2];
       }
    }

    double dy_mhd_psi_line, dy_mhd_u_psi, dy_mhd_chi_line, dy_mhd_u_chi;
    if (upwind_enable==0) {
       indexm1 = CCTK_GFINDEX3D(cctkGH,i,j-1,k);
       index1  = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
       dy_mhd_psi_line = hdYm1*(mhd_psi_line[index1] - mhd_psi_line[indexm1]);
       dy_mhd_chi_line = hdYm1*(mhd_chi_line[index1] - mhd_chi_line[indexm1]);
       dy_mhd_u_psi    = hdYm1*(mhd_u_psi[index1] - mhd_u_psi[indexm1]);
       dy_mhd_u_chi    = hdYm1*(mhd_u_chi[index1] - mhd_u_chi[indexm1]);
    }

    if (upwind_enable==1) {
      if (vyL > 0) {
         indexm1 = CCTK_GFINDEX3D(cctkGH,i,j-1,k);
         indexm2 = CCTK_GFINDEX3D(cctkGH,i,j-2,k);
         double cm2 = 0.5/dY;
         double cm1 = -4.0*cm2;
         double c0 = 3.0*cm2;
         dy_mhd_psi_line = c0*mhd_psi_line[index] + cm1*mhd_psi_line[indexm1] + cm2*mhd_psi_line[indexm2];
         dy_mhd_u_psi = c0*mhd_u_psi[index] + cm1*mhd_u_psi[indexm1] + cm2*mhd_u_psi[indexm2];
         dy_mhd_chi_line = c0*mhd_chi_line[index] + cm1*mhd_chi_line[indexm1] + cm2*mhd_chi_line[indexm2];
         dy_mhd_u_chi = c0*mhd_u_chi[index] + cm1*mhd_u_chi[indexm1] + cm2*mhd_u_chi[indexm2];
      } else {
         index1 = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
         index2 = CCTK_GFINDEX3D(cctkGH,i,j+2,k);
         double c2 = -0.5/dY;
         double c1 = -4.0*c2;
         double c0 = 3.0*c2;
         dy_mhd_psi_line = c0*mhd_psi_line[index] + c1*mhd_psi_line[index1] + c2*mhd_psi_line[index2];
         dy_mhd_u_psi = c0*mhd_u_psi[index] + c1*mhd_u_psi[index1] + c2*mhd_u_psi[index2];
         dy_mhd_chi_line = c0*mhd_chi_line[index] + c1*mhd_chi_line[index1] + c2*mhd_chi_line[index2];
         dy_mhd_u_chi = c0*mhd_u_chi[index] + c1*mhd_u_chi[index1] + c2*mhd_u_chi[index2];
      }
    }

    double dz_mhd_psi_line, dz_mhd_u_psi, dz_mhd_chi_line, dz_mhd_u_chi;
    if (upwind_enable==0) {
       indexm1 = CCTK_GFINDEX3D(cctkGH,i,j,k-1);
       index1  = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
       dz_mhd_psi_line = hdZm1*(mhd_psi_line[index1] - mhd_psi_line[indexm1]);
       dz_mhd_chi_line = hdZm1*(mhd_chi_line[index1] - mhd_chi_line[indexm1]);
       dz_mhd_u_psi    = hdZm1*(mhd_u_psi[index1] - mhd_u_psi[indexm1]);
       dz_mhd_u_chi    = hdZm1*(mhd_u_chi[index1] - mhd_u_chi[indexm1]);
    }

    if (upwind_enable==1) { 
      if (vzL > 0) {
         indexm1 = CCTK_GFINDEX3D(cctkGH,i,j,k-1);
         indexm2 = CCTK_GFINDEX3D(cctkGH,i,j,k-2);
         double cm2 = 0.5/dZ;
         double cm1 = -4.0*cm2;
         double c0 = 3.0*cm2;
         dz_mhd_psi_line = c0*mhd_psi_line[index] + cm1*mhd_psi_line[indexm1] + cm2*mhd_psi_line[indexm2];
         dz_mhd_u_psi = c0*mhd_u_psi[index] + cm1*mhd_u_psi[indexm1] + cm2*mhd_u_psi[indexm2];
         dz_mhd_chi_line = c0*mhd_chi_line[index] + cm1*mhd_chi_line[indexm1] + cm2*mhd_chi_line[indexm2];
         dz_mhd_u_chi = c0*mhd_u_chi[index] + cm1*mhd_u_chi[indexm1] + cm2*mhd_u_chi[indexm2];
      } else {
         index1 = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
         index2 = CCTK_GFINDEX3D(cctkGH,i,j,k+2);
         double c2 = -0.5/dZ;
         double c1 = -4.0*c2;
         double c0 = 3.0*c2;
         dz_mhd_psi_line = c0*mhd_psi_line[index] + c1*mhd_psi_line[index1] + c2*mhd_psi_line[index2];
         dz_mhd_u_psi = c0*mhd_u_psi[index] + c1*mhd_u_psi[index1] + c2*mhd_u_psi[index2];
         dz_mhd_chi_line = c0*mhd_chi_line[index] + c1*mhd_chi_line[index1] + c2*mhd_chi_line[index2];
         dz_mhd_u_chi = c0*mhd_u_chi[index] + c1*mhd_u_chi[index1] + c2*mhd_u_chi[index2];
      }
    }


    // Final, compute mhd_psi_line_rhs, mhd_u_psi_rhs, mhd_chi_line_rhs, and mhd_u_chi_rhs
    mhd_psi_line_rhs[index] = -vxL*dx_mhd_psi_line - vyL*dy_mhd_psi_line - vzL*dz_mhd_psi_line 
        			+ mu_Btx*dx_mhd_u_psi + mu_Bty*dy_mhd_u_psi + mu_Btz*dz_mhd_u_psi;

    mhd_chi_line_rhs[index] = -vxL*dx_mhd_chi_line - vyL*dy_mhd_chi_line - vzL*dz_mhd_chi_line
                                + mu_Btx*dx_mhd_u_chi + mu_Bty*dy_mhd_u_chi + mu_Btz*dz_mhd_u_chi;

    mhd_u_psi_rhs[index] = (mu_Btx*dx_mhd_psi_line + mu_Bty*dy_mhd_psi_line + mu_Btz*dz_mhd_psi_line )
    				- lambda_line*mhd_u_psi[index];

    mhd_u_chi_rhs[index] = (mu_Btx*dx_mhd_chi_line + mu_Bty*dy_mhd_chi_line + mu_Btz*dz_mhd_chi_line)
                                - lambda_line*mhd_u_chi[index];


  } // end of for-loop

}

extern "C" void CCTK_FCALL CCTK_FNAME(field_line_tracer_cpp)
  (const cGH **cctkGH,int *cctk_lsh, int *nghostzones,
   double *alpm1, double *shiftx, double *shifty, double *shiftz,
   double *phi,double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *vx,double *vy,double *vz, double *Bxtilde,double *Bytilde,double *Bztilde,
   double *mhd_psi_line, double *mhd_u_psi, double *mhd_chi_line, double *mhd_u_chi,
   double *mhd_psi_line_rhs, double *mhd_u_psi_rhs, double *mhd_chi_line_rhs, double *mhd_u_chi_rhs,
   double *lambda_line, double *dX, double *dY, double *dZ, 
   double *X, double *Y, double *Z, double *r_ah,
   double *x_ah, double *y_ah, double *z_ah)
{
  field_line_tracer_cpp(*cctkGH,cctk_lsh, nghostzones, 
 		        alpm1, shiftx,shifty,shiftz, 
		        phi,gxx,gxy,gxz,gyy,gyz,gzz,
			vx,vy,vz,Bxtilde,Bytilde,Bztilde,
			mhd_psi_line, mhd_u_psi, mhd_chi_line, mhd_u_chi,
			mhd_psi_line_rhs, mhd_u_psi_rhs, mhd_chi_line_rhs, mhd_u_chi_rhs,
			*lambda_line, *dX,*dY,*dZ, X,Y,Z, *r_ah, 
                        *x_ah, *y_ah, *z_ah);
}
