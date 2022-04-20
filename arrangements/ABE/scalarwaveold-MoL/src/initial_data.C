//----------------------------
// Initial Data Subroutine
//----------------------------

#include <stdio.h>
#include <math.h>
#include <stddef.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "Symmetry.h"

extern "C" double Analytic(double t, double r, double Amp, double Width);

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

  for(int k=kstart; k<kend; k++) for(int j=jstart; j<jend; j++) for(int i=istart; i<iend; i++) {
    int vindex = CCTK_GFINDEX3D(cctkGH,i,j,k);
    double radius=r[vindex];
    double time = 0.0;

    //Need this if timelevels=3 is set.
    /*
    phi_p_p[vindex] = Analytic(time,radius,amplitude,width);
    phidot_p_p[vindex] = 0.0;

    phi_p[vindex] = Analytic(time,radius,amplitude,width);
    phidot_p[vindex] = 0.0;
    */
    //Following must be set so that t=0 data are correct
    phi[vindex] = Analytic(time,radius,amplitude,width);
    phidot[vindex] = 0.0;

    phi_p[vindex] = phi[vindex];
    phidot_p[vindex] = phidot[vindex];

    phi_p_p[vindex] = phi[vindex];
    phidot_p_p[vindex] = phidot[vindex];
  }
}
