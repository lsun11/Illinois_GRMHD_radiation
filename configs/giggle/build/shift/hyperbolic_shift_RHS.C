//--------------------------------------------------------
// Right-hand side of hyperbolic shift evolution equation
//--------------------------------------------------------

#include "stdio.h"
#include "math.h"
#include "cctk.h"

#define SQR(x) ((x) * (x))

extern "C" void CCTK_FCALL hyperbolic_shift_rhs_
  (const cGH **cctkGH,int *nghostzones,int *cctk_lsh,
   double *k1,double *k2,double *k3, 
   double *RadiusDerivativeG,double *lapseG,
   double *shiftx_rhsG,double *shifty_rhsG,double *shiftz_rhsG,
   double *dtshiftx_oldG,double *dtshifty_oldG,double *dtshiftz_oldG,
   double *dtshiftx_rhsG,double *dtshifty_rhsG,double *dtshiftz_rhsG,
   double *Gammax_rhsG,double *Gammay_rhsG,double *Gammaz_rhsG,
   double *GammaxG,double *GammayG,double *GammazG,
   double *Gammax_iG,double *Gammay_iG,double *Gammaz_iG);


extern "C" void hyperbolic_shift_rhs(const cGH *cctkGH,int *nghostzones,int *cctk_lsh,
				     double k1,double k2,double k3, 
				     double *RadiusDerivativeG,double *lapseG,
				     double *shiftx_rhsG,double *shifty_rhsG,double *shiftz_rhsG,
				     double *dtshiftx_oldG,double *dtshifty_oldG,double *dtshiftz_oldG,
				     double *dtshiftx_rhsG,double *dtshifty_rhsG,double *dtshiftz_rhsG,
				     double *Gammax_rhsG,double *Gammay_rhsG,double *Gammaz_rhsG,
				     double *GammaxG,double *GammayG,double *GammazG,
				     double *Gammax_iG,double *Gammay_iG,double *Gammaz_iG) {
//

  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

    double lapse = lapseG[index];
    double RadiusDerivative = RadiusDerivativeG[index];
    double Gammax_rhs = Gammax_rhsG[index];
    double Gammay_rhs = Gammay_rhsG[index];
    double Gammaz_rhs = Gammaz_rhsG[index];
    
    double Gammax = GammaxG[index];
    double Gammay = GammayG[index];
    double Gammaz = GammazG[index];
    
    double Gammax_i = 0.0;
    double Gammay_i = 0.0;
    double Gammaz_i = 0.0;

    double dtshiftx_old = dtshiftx_oldG[index];
    double dtshifty_old = dtshifty_oldG[index];
    double dtshiftz_old = dtshiftz_oldG[index];

    shiftx_rhsG[index] = dtshiftx_old;
    shifty_rhsG[index] = dtshifty_old;
    shiftz_rhsG[index] = dtshiftz_old;

    // Note: the equations for the dtshifti_rhs have been changed to take
    // into account the fisheye transformation. When fisheye coordinates
    // are not used, RadiusDerivative = 1, and the equations reduce to their
    // original form.

    dtshiftx_rhsG[index] = k1 * (lapse + 1.0) / SQR(RadiusDerivative) * Gammax_rhs 
      - k2 * dtshiftx_old
      + k3 * (lapse + 1.0) * (Gammax - Gammax_i);

    dtshifty_rhsG[index] = k1 * (lapse + 1.0) / SQR(RadiusDerivative) * Gammay_rhs 
      - k2 * dtshifty_old
      + k3 * (lapse + 1.0) * (Gammay - Gammay_i);

    dtshiftz_rhsG[index] = k1 * (lapse + 1.0) / SQR(RadiusDerivative) * Gammaz_rhs 
      - k2 * dtshiftz_old
      + k3 * (lapse + 1.0) * (Gammaz - Gammaz_i);

    //if(i==1 && j==1 && k==1) printf("   %.15e\t%.15e\t%.15e\n",Gammax_rhs,Gammax,Gammax_i);
    //if(i==1 && j==1 && k==1) printf("   %.15e\t%.15e\t%.15e\n",dtshiftx_rhsG[index],dtshifty_rhsG[index],dtshiftz_rhsG[index]);

  }

}

extern "C" void CCTK_FCALL hyperbolic_shift_rhs_
  (const cGH **cctkGH,int *nghostzones,int *cctk_lsh,
   double *k1,double *k2,double *k3, 
   double *RadiusDerivativeG,double *lapseG,
   double *shiftx_rhsG,double *shifty_rhsG,double *shiftz_rhsG,
   double *dtshiftx_oldG,double *dtshifty_oldG,double *dtshiftz_oldG,
   double *dtshiftx_rhsG,double *dtshifty_rhsG,double *dtshiftz_rhsG,
   double *Gammax_rhsG,double *Gammay_rhsG,double *Gammaz_rhsG,
   double *GammaxG,double *GammayG,double *GammazG,
   double *Gammax_iG,double *Gammay_iG,double *Gammaz_iG)
{
  hyperbolic_shift_rhs(*cctkGH,nghostzones,cctk_lsh,
		       *k1,*k2,*k3, 
		       RadiusDerivativeG,lapseG,
		       shiftx_rhsG,shifty_rhsG,shiftz_rhsG,
		       dtshiftx_oldG,dtshifty_oldG,dtshiftz_oldG,
		       dtshiftx_rhsG,dtshifty_rhsG,dtshiftz_rhsG,
		       Gammax_rhsG,Gammay_rhsG,Gammaz_rhsG,
		       GammaxG,GammayG,GammazG,
		       Gammax_iG,Gammay_iG,Gammaz_iG);
}
