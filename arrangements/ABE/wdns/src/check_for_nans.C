// Here is a nan checker routine
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <sys/time.h>
#include <math.h>
#include "cctk.h"

extern "C" void CCTK_FCALL CCTK_FNAME(check_for_nans)
  (const cGH **cctkGH,int *ext,
   double *Bx,double *By,double *Bz);

void check_for_nans(const cGH *cctkGH,int *ext,
				       double *Bx,double *By,double *Bz
				       ) {
#pragma omp parallel for
  for(int k=0;k<ext[2];k++)
    for(int j=0;j<ext[1];j++)
      for(int i=0;i<ext[0];i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);


	if(isnan(Bx[index]) || isnan(By[index]) || isnan(Bz[index])) { 
	  printf("Nans Detected: %e. %e %e %\n",Bx[index],By[index],Bz[index]); 
	  exit(1); 
	}
      }
}


extern "C" void CCTK_FCALL CCTK_FNAME(check_for_nans)
  (const cGH **cctkGH,int *ext,
   double *Bx,double *By,double *Bz) {
  check_for_nans (*cctkGH,ext,
				     Bx,By,Bz);

}
