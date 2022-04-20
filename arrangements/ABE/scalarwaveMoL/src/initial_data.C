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
  double dx = CCTK_DELTA_SPACE(0);
  double dy = CCTK_DELTA_SPACE(1);
  double dz = CCTK_DELTA_SPACE(2);
  double dt = CCTK_DELTA_TIME;

  printf("computing ID: %d %d %d\n",iend,jend,kend);
  //printf("ijkminmax %d %d %d %d %d %d %e %e\n",istart,jstart,kstart,iend,jend,kend,amplitude,width);

  for(int k=kstart; k<kend; k++) for(int j=jstart; j<jend; j++) for(int i=istart; i<iend; i++) {
    int vindex = CCTK_GFINDEX3D(cctkGH,i,j,k);
    double radius=r[vindex];
    double time=0.0;

    //Following must be set so that t=0 data are correct
    //if(radius==0) radius=1e-10;
    phi[vindex] = Analytic(time,radius,amplitude,width);
    //if(i==28 && j==27 && k==0) printf("bb %d %d %d , %e %e %e %e, %e\n",i,j,k,time,radius,amplitude,width,Analytic(time,radius,amplitude,width));
    phidot[vindex] = 0.0;

    // *** TEST ***
    double xl = x[vindex] + 0.5*dx;
    double yl = y[vindex] + 0.5*dy;
    double zl = z[vindex];
    phi_stagger[vindex] = 1.0 + xl + 0.2*xl*xl 
        		  + 0.5*yl + 0.7*yl*yl 
        		  + 1.3*zl*zl;
    phidot_stagger[vindex] = -0.1 + 0.02*xl + 0.01*xl*xl 
        		     -0.03*yl + 0.002*yl*yl + 0.008*zl*zl;
    phi_stagger_p[vindex] = phi_stagger[vindex];
    phidot_stagger_p[vindex] = phidot_stagger[vindex];
    phi_stagger_p_p[vindex] = phi_stagger[vindex];
    phidot_stagger_p_p[vindex] = phidot_stagger[vindex];
    // ************
    
    //double xl = x[vindex];
    //double yl = y[vindex] + 0.5*dy;
    //double zl = z[vindex] + 0.5*dz;
    //radius = sqrt(xl*xl + yl*yl + zl*zl);
    //phi_stagger[vindex] = Analytic(time,radius,amplitude,width);
    //phidot_stagger[vindex] = 0.0;
    //phi_stagger_p[vindex] = phi_stagger[vindex];
    //phidot_stagger_p[vindex] = phidot_stagger[vindex];

    //if(j==1 && k==1) printf("okay: %d %e\n",i,x[vindex]);
    //if(i==20 && j==20 && k==14) printf("hi %d %e %e %e %e, %e - %e\n",k,z[vindex],x[vindex],y[vindex],radius,Analytic(time,radius,amplitude,width),phi[vindex]);
    phi_p_p[vindex] = phi[vindex];
    phidot_p_p[vindex] = phidot[vindex];

  }

  //printf("computing ID at time: %e, %e %e %d %d %d %e\n",cctk_time,amplitude,width,iend,jend,kend,phi[CCTK_GFINDEX3D(cctkGH,28,27,0)]);
  
}
