//--------------------------------------------------------
// Right-hand side of hyperbolic lapse evolution equation
//--------------------------------------------------------

#include "stdio.h"
#include "math.h"
#include "cctk.h"

#define SQR(x) ((x) * (x))

extern "C" void CCTK_FCALL hyperbolic_lapse_rhs_
  (const cGH **cctkGH,int *nghostzones,int *cctk_lsh,
   double *a1, double *a2, double *a3, 
   double *rG,double *PhysicalRadiusG,double *RadiusDerivativeG,
   double *phiG,double *lapse_oldG,double *lapse_rhsG,
   double *dtlapse_oldG,double *dtlapse_rhsG,
   double *K_rhsG,double *trKG,double *K_initG);

extern "C" void hyperbolic_lapse_rhs(const cGH *cctkGH,int *nghostzones,int *cctk_lsh,
				     double a1, double a2, double a3, 
				     double *rG,double *PhysicalRadiusG,double *RadiusDerivativeG,
				     double *phiG,double *lapse_oldG,double *lapse_rhsG,
				     double *dtlapse_oldG,double *dtlapse_rhsG,
				     double *K_rhsG,double *trKG,double *K_initG) {

//

  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

    double lapse_old = lapse_oldG[index];
    double RadiusDerivative = RadiusDerivativeG[index];
    double PhysicalRadius = PhysicalRadiusG[index];

    double K_init = 0.0;
    double K_rhs = K_rhsG[index];
    double trK = trKG[index];

    double phi = phiG[index];

    double dtlapse_rhs = dtlapse_rhsG[index];
    double dtlapse_old = dtlapse_oldG[index];

    lapse_rhsG[index] = dtlapse_old * (lapse_old + 1.0);

    dtlapse_rhs = rG[index];
    dtlapse_rhs = pow(RadiusDerivative,(2.0/3.0)) * pow(PhysicalRadius/dtlapse_rhs,(4.0/3.0));

    dtlapse_rhs = - a1 * ( (lapse_old + 1.0) * K_rhs + a2 * dtlapse_old + 
				   a3 * dtlapse_rhs * (trK - K_init) * exp(-4.0*phi) * 
				   (lapse_old + 1.0) );

    dtlapse_rhsG[index] = dtlapse_rhs;

  }

}

extern "C" void CCTK_FCALL hyperbolic_lapse_rhs_
  (const cGH **cctkGH,int *nghostzones,int *cctk_lsh,
   double *a1, double *a2, double *a3, 
   double *rG,double *PhysicalRadiusG,double *RadiusDerivativeG,
   double *phiG,double *lapse_oldG,double *lapse_rhsG,
   double *dtlapse_oldG,double *dtlapse_rhsG,
   double *K_rhsG,double *trKG,double *K_initG)
{
  hyperbolic_lapse_rhs(*cctkGH,nghostzones,cctk_lsh,
		       *a1, *a2, *a3, 
		       rG,PhysicalRadiusG,RadiusDerivativeG,
		       phiG,lapse_oldG,lapse_rhsG,
		       dtlapse_oldG,dtlapse_rhsG,
		       K_rhsG,trKG,K_initG);
}
