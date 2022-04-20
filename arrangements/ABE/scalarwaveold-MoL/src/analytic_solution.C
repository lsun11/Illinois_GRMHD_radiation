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

  //Uncomment following lines to set phi_anal.  Note that this is an expensive operation, requiring N exp() operations, so
  //this WILL slow down the code (by a factor of ~1.5)

  for (int k=kstart; k<kend; k++) for (int j=jstart; j<jend; j++) for (int i=istart; i<iend; i++) {
    int vindex = CCTK_GFINDEX3D(cctkGH,i,j,k);
    double r=sqrt(x[vindex]*x[vindex]+y[vindex]*y[vindex]+z[vindex]*z[vindex]);
    double time = cctk_time;

    phi_anal[vindex] = Analytic(time,r,amplitude,width) - phi[vindex];

  }

}
