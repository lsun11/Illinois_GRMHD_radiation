//-----------------------------------------------------------------------
// Read WDNS bin. files and do fancy things with them...
//-----------------------------------------------------------------------
#include <iostream>
#include <sstream> 
#include <string>  

#include <unistd.h>

#include <cmath>
#include <iomanip> 
#include <fstream> 
#include <stdlib.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

using namespace std;

extern "C" void CCTK_FCALL CCTK_FNAME(calc_u_0_PoPcold)
  (const cGH **cctkGH, double *xmin, double *ymin, double *zmin, 
   double *dx, double *dy, double *dz, int *ext, 
   double *lapm1, double *shiftx, double *shifty, double *shiftz, 
   double *psi, double *gmxx, double *gmxy, double *gmxz, double *gmyy, double *gmyz, double *gmzz, 
   double *rho_b, double *u_up0, double *Pr, double *vx, double *vy, double *vz, double *K, double *u_0, 
   double &gamma1, double &gamma2, double &gamma3, double &kappa1, double &kappa2, double &kappa3, double &rhoo1, double &rhoo2);


extern "C" void calc_u_0_PoPcold(const cGH *cctkGH, double xmin, double ymin, double zmin, 
				 double dx, double dy, double dz, int *ext, 
				 double *lapm1, double *shiftx, double *shifty, double *shiftz, 
				 double *psi, double *gmxx, double *gmxy, double *gmxz, double *gmyy, double *gmyz, double *gmzz, 
				 double *rho_b, double *u_up0, double *Pr, double *vx, double *vy, double *vz, double *K, double *u_0, 
				 double gamma1, double gamma2, double gamma3, double kappa1, double kappa2, double kappa3, double rhoo1, double rhoo2)
{


    int nx = ext[0];
    int ny = ext[1];
    int nz = ext[2];



    double Pcold;

    for(int k=0;k<nz;k++){
      //      double z = zmin + k*dz;
      for(int j=0;j<ny;j++){
	//	double y = ymin + j*dy;
	for(int i=0;i<nx;i++){
	  //	  double x = xmin + i*dx;

	  int ind = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  double rhoo = rho_b[ind];
	  double pres  = Pr[ind];
	  double u0 = u_up0[ind];
	  double conf = psi[ind];
	  double lap = lapm1[ind]+1.0;
	  double betx = shiftx[ind];
	  double bety = shifty[ind];
	  double betz = shiftz[ind];    
	  double gxx = gmxx[ind];
	  double gxy = gmxy[ind];
	  double gxz = gmxz[ind];
	  double gyy = gmyy[ind];
	  double gyz = gmyz[ind];
	  double gzz = gmzz[ind];
	  double v3x = vx[ind];
	  double v3y = vy[ind];
	  double v3z = vz[ind];

	  u_0[ind] = u0*(-lap*lap + conf*conf*conf*conf*( gxx*betx*(betx+v3x) + 2.*gxy*betx*(bety+v3y) + 2.*gxz*betx*(betz+v3z) + gyy*bety*(bety+v3y) + 2.*gyz*bety*(betz+v3z) + gzz*betz*(betz+v3z)   )  );
  

	  if ( (rhoo - rhoo1)/rhoo1 <= 1.e-15){
	    Pcold = kappa1*pow(rhoo,gamma1);
	  }
	  else if ( ((rhoo1 - rhoo)/rhoo < 1.e-15) && ((rhoo - rhoo2)/rhoo2 <= 1.e-15) ){
	    Pcold = kappa2*pow(rhoo,gamma2);
	  }
	  else{
	    Pcold = kappa3*pow(rhoo,gamma3);
	  }


	  K[ind] = pres/Pcold;
	  //K[ind] = Pcold;          


	}
      }
    }
    
  
}



extern "C" void CCTK_FCALL CCTK_FNAME(calc_u_0_PoPcold)
  (const cGH **cctkGH, double *xmin, double *ymin, double *zmin, 
   double *dx, double *dy, double *dz, int *ext, 
   double *lapm1, double *shiftx, double *shifty, double *shiftz, 
   double *psi, double *gmxx, double *gmxy, double *gmxz, double *gmyy, double *gmyz, double *gmzz, 
   double *rho_b, double *u_up0, double *Pr, double *vx, double *vy, double *vz, double *K, double *u_0, 
   double &gamma1, double &gamma2, double &gamma3, double &kappa1, double &kappa2, double &kappa3, double &rhoo1, double &rhoo2)
{
  calc_u_0_PoPcold(*cctkGH, *xmin, *ymin, *zmin, *dx, *dy, *dz, ext, 
   lapm1, shiftx, shifty, shiftz, 
   psi, gmxx, gmxy, gmxz, gmyy, gmyz, gmzz, 
   rho_b, u_up0, Pr, vx, vy, vz, K, u_0, 
   gamma1, gamma2, gamma3, kappa1, kappa2, kappa3, rhoo1, rhoo2);
}
