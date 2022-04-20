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
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

using namespace std;

extern "C" void CCTK_FCALL bhns_checksum_
  (const cGH **cctkGH,
   double *xmin,double *ymin,double *zmin,
   double *dx,double *dy,double *dz,int *ext, int &checksum);

extern "C" void bhns_checksum(const cGH *cctkGH,
			      double xmin,double ymin,double zmin,
			      double dx,double dy,double dz,int *ext, int &checksum){



  static int nStatic;
  nStatic += 1;
  double xmax=xmin+ext[0]*dx;
  double ymax=ymin+ext[1]*dy;
  double zmax=zmin+ext[2]*dz;

checksum = ext[0]*301+ext[1]*10+ext[2]*1000000  + 1241*CCTK_MyProc(cctkGH) + (int)(999999.*dx) + (int)(100100100.*dy) + (int)(800000.*dz) + 
    (int)(fabs(xmin)*10.)     + (int)(fabs(ymin)*101.)      + (int)(fabs(zmin)*1001.) +
    (int)(fabs(xmax)*999999.) + (int)(fabs(ymax)*30000.)    + (int)(fabs(zmax)*10100.)+
    (int)(xmax*xmax*191919. ) + (int)(ymax*ymax*30220. )    + (int)(zmax*zmax*20202.) +nStatic;    

}

extern "C" void CCTK_FCALL bhns_checksum_
  (const cGH **cctkGH,
   double *xmin,double *ymin,double *zmin,
   double *dx,double *dy,double *dz,int *ext, int &checksum)
{  
  bhns_checksum(*cctkGH,
		*xmin,*ymin,*zmin,*dx,*dy,*dz,ext,checksum);
}
