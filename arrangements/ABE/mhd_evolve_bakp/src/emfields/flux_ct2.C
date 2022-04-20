//----------------------------------------------------------//
// Constraint Transport: f^{ij} -> \tilde{f}^{ij}	    //
//----------------------------------------------------------//

#include "cctk.h"

#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#define NO_SYMM 0
#define OCTANT  2
#define AXISYM  4

extern "C" void CCTK_FCALL CCTK_FNAME(flux_ct_cpp)
  (const cGH **cctkGH,int *ext,double *X,double *Y,double *Z,
   double *fxy,double *fxz,double *fyx,double *fyz,double *fzx,double *fzy,
   double *ftxy,double *ftxz,double *ftyx,double *ftyz,double *ftzx,double *ftzy,
   double &Sym_Bz,int &Symmetry);

void flux_ct_cpp(const cGH *cctkGH,int *ext,double *X,double *Y,double *Z,
	     double *fxy,double *fxz,double *fyx,double *fyz,double *fzx,double *fzy,
	     double *ftxy,double *ftxz,double *ftyx,double *ftyz,double *ftzx,double *ftzy,
	     double &Sym_Bz,int &Symmetry) {

  /*
  if (Symmetry==NO_SYMM || Z(1,1,kmin) < 0.d0) {
     Symz = 1.d0
  } else {
     Symz = Sym_Bz
     //     Symz = -Sym_Bz
  }
  */

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

    ftxy[index] = 0.0;
    ftxz[index] = 0.0;
    ftyx[index] = 0.0;
    ftyz[index] = 0.0;
    ftzx[index] = 0.0;
    ftzy[index] = 0.0;
    
    ftxy[index] = fxy[index];
    ftxz[index] = fxz[index];
    ftyx[index] = fyx[index];
    ftyz[index] = fyz[index];
    ftzx[index] = fzx[index];
    ftzy[index] = fzy[index];
  }


  /* 
     IS THIS NECESSARY??????????????????????????????????????????????
     //Enforce fzx, fzy = 0 on equatorial plane//
     if (Sym_Bz < 0.0) {
     fzx(:,:,kmin) = 0.0;
     if (Symmetry != AXISYM) fzy(:,:,kmin) = 0.0;
     }
  */

  if (Symmetry != AXISYM) {
#pragma omp parallel for
    for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) {
      int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
      int indexip1 = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
      int indexim1 = CCTK_GFINDEX3D(cctkGH,i-1,j,k);
      int indexjp1 = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
      int indexjm1 = CCTK_GFINDEX3D(cctkGH,i,j-1,k);
      int indexkp1 = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
      int indexkm1 = CCTK_GFINDEX3D(cctkGH,i,j,k-1);
      int indexip1jm1 = CCTK_GFINDEX3D(cctkGH,i+1,j-1,k);
      int indexim1jp1 = CCTK_GFINDEX3D(cctkGH,i-1,j+1,k);
      int indexjp1km1 = CCTK_GFINDEX3D(cctkGH,i,j+1,k-1);
      int indexjm1kp1 = CCTK_GFINDEX3D(cctkGH,i,j-1,k+1);
      int indexip1km1 = CCTK_GFINDEX3D(cctkGH,i+1,j,k-1);
      int indexim1kp1 = CCTK_GFINDEX3D(cctkGH,i-1,j,k+1);
      // Calculate \tilde{f}^{xy}(i-1/2,j,k) and store in ftxy[index],
      //	   \tilde{f}^{yx}(i,j-1/2,k) and store in ftyx[index]
      ftxy[index]=2.0*fxy[index]-fyx[index];
      ftxy[index]=ftxy[index]+fxy[indexjp1]-fyx[indexjp1];
      ftxy[index]=ftxy[index]+fxy[indexjm1];
      ftxy[index]=ftxy[index]-fyx[indexim1];
      ftxy[index]=ftxy[index]-fyx[indexim1jp1];

      ftyx[index]=2.0*fyx[index]-fxy[index];
      ftyx[index]=ftyx[index]+fyx[indexip1]-fxy[indexip1];
      ftyx[index]=ftyx[index]+fyx[indexim1];
      ftyx[index]=ftyx[index]-fxy[indexjm1];
      ftyx[index]=ftyx[index]-fxy[indexip1jm1];
      
      ftxy[index]=0.125*ftxy[index];
      ftyx[index]=0.125*ftyx[index];

      // Calculate \tilde{f}^{yz}(i,j-1/2,k) and store in ftyz[index],
      //           \tilde{f}^{zy}(i,j,k-1/2) and store in ftzy[index]
      ftyz[index]=2.0*fyz[index]-fzy[index];
      ftyz[index]=ftyz[index]+fyz[indexkp1]-fzy[indexkp1];
      ftyz[index]=ftyz[index]+fyz[indexkm1];
      ftyz[index]=ftyz[index]-fzy[indexjm1];
      ftyz[index]=ftyz[index]-fzy[indexjm1kp1];

      ftzy[index]=2.0*fzy[index]-fyz[index];
      ftzy[index]=ftzy[index]+fzy[indexjp1]-fyz[indexjp1];
      ftzy[index]=ftzy[index]+fzy[indexjm1];
      ftzy[index]=ftzy[index]-fyz[indexkm1];
      ftzy[index]=ftzy[index]-fyz[indexjp1km1];

      ftyz[index]=0.125*ftyz[index];
      ftzy[index]=0.125*ftzy[index];

      // Calculate \tilde{f}^{xz}(i-1/2,j,k) and store in ftxz[index],
      //           \tilde{f}^{zx}(i,j,k-1/2) and store in ftzx[index]
      ftxz[index]=2.0*fxz[index]-fzx[index];
      ftxz[index]=ftxz[index]+fxz[indexkp1]-fzx[indexkp1];
      ftxz[index]=ftxz[index]+fxz[indexkm1];
      ftxz[index]=ftxz[index]-fzx[indexim1];
      ftxz[index]=ftxz[index]-fzx[indexim1kp1];

      ftzx[index]=2.0*fzx[index]-fxz[index];
      ftzx[index]=ftzx[index]+fzx[indexip1]-fxz[indexip1];
      ftzx[index]=ftzx[index]+fzx[indexim1];
      ftzx[index]=ftzx[index]-fxz[indexkm1];
      ftzx[index]=ftzx[index]-fxz[indexip1km1];

      ftxz[index]=0.125*ftxz[index];
      ftzx[index]=0.125*ftzx[index];
    }
  } else {
    imin=jmin=kmin=0;
    imax++; jmax++; kmax++;
    for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) {
      int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
      ftxy[index]=fxy[index];
      ftzy[index]=fzy[index];
    }
  }
}

extern "C" void CCTK_FCALL CCTK_FNAME(flux_ct_cpp)
  (const cGH **cctkGH,int *ext,double *X,double *Y,double *Z,
   double *fxy,double *fxz,double *fyx,double *fyz,double *fzx,double *fzy,
   double *ftxy,double *ftxz,double *ftyx,double *ftyz,double *ftzx,double *ftzy,
   double &Sym_Bz,int &Symmetry)
{
  flux_ct_cpp(*cctkGH,ext,X,Y,Z,
	      fxy,fxz,fyx,fyz,fzx,fzy,
	      ftxy,ftxz,ftyx,ftyz,ftzx,ftzy,
	      Sym_Bz,Symmetry);
}
