//-----------------------------------------------------------------------------
//
// Initial Data, Analytical Solution Subroutines
//
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <math.h>
#include <stddef.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "Symmetry.h"

#define val(gridfunc,i,j,k)  gridfunc[CCTK_GFINDEX3D(cctkGH,i,j,k)]

double Analytic(double t, double r, double Amp, double Width);

//======================================================
//
// Setup Initial Data
//
//======================================================
static char *rcsid = "$Meow...$";

CCTK_FILEVERSION(ABE_ScalarWave_Setup_InitialData)

extern "C" void ABE_ScalarWave_Setup_InitialData(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  int istart = 0;
  int jstart = 0;
  int kstart = 0;
  int iend = cctk_lsh[0];
  int jend = cctk_lsh[1];
  int kend = cctk_lsh[2];

  for (int k=kstart; k<kend; k++) for (int j=jstart; j<jend; j++) for (int i=istart; i<iend; i++) {
    int vindex = CCTK_GFINDEX3D(cctkGH,i,j,k);
    double r=sqrt(x[vindex]*x[vindex]+y[vindex]*y[vindex]+z[vindex]*z[vindex]);

    //Need this if timelevels=3 is set.
    phi_p[vindex] = Analytic(0.0,r,amplitude,width);
    phidot_p[vindex] = 0;

    //Following must be set so that t=0 data are correct
    phi[vindex] = Analytic(0.0,r,amplitude,width);
    phidot[vindex] = 0;
  }
}

//======================================================
//
// Find (analytical solution - numerical solution)
//
//======================================================

double Analytic(double t, double r, double Amp, double Width)
{
  double u = r+t;
  double v = r-t;
  return 0.5*Amp*( u/r*exp(-Width*u*u) + v/r*exp(-Width*v*v));
}

extern "C" void Compute_Anal(CCTK_ARGUMENTS)
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
    
    //if(i==2 && j==1 && k==0)
    //    if(j<2 && k==0 && i==2 && time>0) printf("%e\t%e\t%e\t%e\t%d\t%d\n",x[vindex] ,z[vindex],Analytic(time,r,amplitude,width),phi[vindex],i,k);
    //    if(j==1) printf("%e\t%e\t%e\t%e\t%d\t%d\n",x[vindex] ,z[vindex],Analytic(time,r,amplitude,width),phi[vindex],i,k);

    phi_anal[vindex] = Analytic(time,r,amplitude,width) - phi[vindex];
  }
  //  printf("TIEEME! %e\n",cctk_time);
}
