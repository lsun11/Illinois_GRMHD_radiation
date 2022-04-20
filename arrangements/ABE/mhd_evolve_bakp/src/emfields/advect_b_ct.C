//------------------------------------------------------//
// Advect B^i					       //
//------------------------------------------------------//

#include "cctk.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#define NO_SYMM 0
#define OCTANT  2
#define AXISYM  4

extern "C" void CCTK_FCALL CCTK_FNAME(advect_b_ct_cpp)
  (const cGH **cctkGH,int *ext,double *X,double *Y,double *Z,double *Bx_rhs,double *By_rhs,double *Bz_rhs,
   double *fxy,double *fxz,double *fyx,double *fyz,double *fzx,double *fzy, int &Symmetry);


//
extern "C" void advect_b_ct_cpp(const cGH *cctkGH,int *ext,double *X,double *Y,double *Z,double *Bx_rhs,double *By_rhs,double *Bz_rhs,
				double *fxy,double *fxz,double *fyx,double *fyz,double *fzx,double *fzy, int &Symmetry) {


  //
  double ddx=1.0/(X[CCTK_GFINDEX3D(cctkGH,1,0,0)]-X[CCTK_GFINDEX3D(cctkGH,0,0,0)]);
  double ddy=1.0/(Y[CCTK_GFINDEX3D(cctkGH,0,1,0)]-Y[CCTK_GFINDEX3D(cctkGH,0,0,0)]);
  double ddz=1.0/(Z[CCTK_GFINDEX3D(cctkGH,0,0,1)]-Z[CCTK_GFINDEX3D(cctkGH,0,0,0)]);
  //

  /* Set up variables used in the grid loop for the physical grid points */
  int imin = 1;
  int jmin = 1;
  int kmin = 1;
  int imax = ext[0] - 1;
  int jmax = ext[1] - 1;
  int kmax = ext[2] - 1;

  //
  //

  // Testing:  Initialize everything to zero first:
#pragma omp parallel for
  for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    int indexip1 = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
    int indexjp1 = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
    int indexkp1 = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
    if (Symmetry!=AXISYM) {
      Bx_rhs[index]=ddy*(fyx[index]-fyx[indexjp1]) + 
	ddz*(fzx[index]-fzx[indexkp1]);
      By_rhs[index]=ddx*(fxy[index]-fxy[indexip1]) + 
	ddz*(fzy[index]-fzy[indexkp1]);
      Bz_rhs[index]=ddx*(fxz[index]-fxz[indexip1]) + 
	ddy*(fyz[index]-fyz[indexjp1]);
    } else {
      Bx_rhs[index]=ddz*(fzx[index]-fzx[indexkp1])/X[index];
      By_rhs[index]=ddx*(fxy[index]-fxy[indexip1]) + 
	ddz*(fzy[index]-fzy[indexkp1]);
      Bz_rhs[index]=ddx*(fxz[index]-fxz[indexip1])/X[index];
    }
  }
  //write(*,*) "advect_b_ct:",Bx_rhs(2,2,2),fyx(2,2,2),fyx(2,2+1,2),fzx(2,2,2),fzx(2,2,2+1)
}

extern "C" void CCTK_FCALL CCTK_FNAME(advect_b_ct_cpp)
  (const cGH **cctkGH,int *ext,double *X,double *Y,double *Z,double *Bx_rhs,double *By_rhs,double *Bz_rhs,
   double *fxy,double *fxz,double *fyx,double *fyz,double *fzx,double *fzy, int &Symmetry)
{
  advect_b_ct_cpp(*cctkGH,ext,X,Y,Z,Bx_rhs,By_rhs,Bz_rhs,
   fxy,fxz,fyx,fyz,fzx,fzy,Symmetry);
}

