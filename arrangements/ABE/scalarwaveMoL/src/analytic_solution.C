//------------------------------------
// Analytic Solution Subroutines
//------------------------------------

#include <stdio.h>
#include <math.h>
#include <stddef.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "Symmetry.h"

#define val(gridfunc,i,j,k)  gridfunc[CCTK_GFINDEX3D(cctkGH,i,j,k)]

//======================================================
// Compute analytic solution at time t, radius r:
//======================================================
extern "C" double Analytic(double t, double r, double Amp, double Width)
{
  r = r+1e-12;
  double u = r+t;
  double v = r-t;
  return 0.5*Amp*( u/r*exp(-Width*u*u) + v/r*exp(-Width*v*v));
}

//======================================================================
// Find (analytical solution - numerical solution), set this to phi_anal
//======================================================================

extern "C" void Compute_Difference_Numerical_Analytic(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  int istart = 0;
  int jstart = 0;
  int kstart = 0;

  int iend = cctk_lsh[0];
  int jend = cctk_lsh[1];
  int kend = cctk_lsh[2];

  double dx = CCTK_DELTA_SPACE(0);
  double dy = CCTK_DELTA_SPACE(1);
  double dz = CCTK_DELTA_SPACE(2);

  //Uncomment following lines to set phi_anal.  Note that this is an expensive operation, requiring N exp() operations, so
  //this WILL slow down the code (by a factor of ~1.5)

  //  phi[CCTK_GFINDEX3D(cctkGH,20,20,14)] = phi_p[CCTK_GFINDEX3D(cctkGH,20,20,14)];

  for (int k=kstart; k<kend; k++) for (int j=jstart; j<jend; j++) for (int i=istart; i<iend; i++) {
    int vindex = CCTK_GFINDEX3D(cctkGH,i,j,k);
    double radius=r[vindex];
    //double r=sqrt(x[vindex]*x[vindex]+y[vindex]*y[vindex]+z[vindex]*z[vindex]);
    double time = cctk_time;

    if(cctk_iteration==0) {
      //phi[vindex] = Analytic(time,radius,amplitude,width);
    }

    phi_analytic[vindex] = Analytic(time,radius,amplitude,width);
    phi_analytic_minus_numeric[vindex] = phi_analytic[vindex] - phi[vindex];

    double xl = x[vindex] + 0.5*dx; 
    double yl = y[vindex] + 0.5*dy;
    double zl = z[vindex];
    //radius = sqrt(xl*xl+yl*yl+zl*zl);
    //double phi_stagger_analytic = Analytic(time,radius,amplitude,width);
    double phi_stagger_analytic = 1.0 + xl + 0.2*xl*xl + 0.5*yl + 0.7*yl*yl 
      + 1.3*zl*zl + time*(-0.1 + 0.02*xl + 0.01*xl*xl 
			  - 0.03*yl + 0.002*yl*yl + 0.008*zl*zl);
    //- 0.03*yl + 0.002*yl*yl + 0.008*zl*zl + 3.0*xl*xl*xl - 1.2*xl*zl*zl + 4*yl*yl*xl);
    //			  - 0.03*yl + 0.002*yl*yl + 0.008*zl*zl + 3.0*xl*xl*xl - 1.2*xl*zl*zl + 4*yl*yl*xl - xl*xl*xl*xl);
    phi_stagger_analytic_minus_numeric[vindex] = (phi_stagger_analytic - phi_stagger[vindex])/(fabs(phi_stagger_analytic)+1.e-10);

    //if(i==20 && j==20 && k==14) printf("guh %d %e: %e %e %e, %e - %e = %e\n",k,time,radius,amplitude,width,Analytic(time,radius,amplitude,width),phi_p[vindex],phi_anal[vindex]);

    /*
    if(cctk_iteration!=0) {
      phi_anal[vindex] = Analytic(time,radius,amplitude,width) - phi_p[vindex];
    } else {
      phi_anal[vindex] = Analytic(time,radius,amplitude,width) - phi_p[vindex];
    }
    */
    //if(fabs(phi_anal[vindex])>1e-1) printf("hi %d %d %d , %e %e %e %e, %e - %e -- %e\n",i,j,k,time,radius,amplitude,width,Analytic(time,radius,amplitude,width),phi[vindex],phi_anal[vindex]);
    //if(i==30 && j==30 && k==4) printf("hi %e %e %e %e, %e - %e\n",time,radius,amplitude,width,Analytic(time,radius,amplitude,width),phi_anal[vindex]);

  }

}
