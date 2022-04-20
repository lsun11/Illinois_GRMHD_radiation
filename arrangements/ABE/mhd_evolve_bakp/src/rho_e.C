//-----------------------------------------------------------------------------
// Advect rho_star and e
//-----------------------------------------------------------------------------

#include "math.h"
#include "cctk.h"

extern "C" void CCTK_FCALL CCTK_FNAME(advect_rho_e_cpp)
  (int *flux_direction, const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry,
   double *rho_star_rhs,double *e_rhs,double *frm,double *fem,
   double *dX,double *dY,double *dZ,double *X,double *Y,double *Z);

extern "C" void advect_rho_e_cpp(int flux_direction,const cGH *cctkGH,int *cctk_lsh,int *nghostzones,int Symmetry,
				 double *rho_star_rhs,double *e_rhs,double *frm,double *fem,
				 double dX,double dY,double dZ,double *X,double *Y,double *Z) {
  int AXISYM = 4;

  /* Set up variables used in the grid loop for the physical grid points */
  int istart = 1;
  int jstart = 1;
  int kstart = 1;
  int iend = cctk_lsh[0]-1;
  int jend = cctk_lsh[1]-1;
  int kend = cctk_lsh[2]-1;

  /*
    int istart = nghostzones[0];
    int jstart = nghostzones[1];
    int kstart = nghostzones[2];
    int iend = cctk_lsh[0] - nghostzones[0];
    int jend = cctk_lsh[1] - nghostzones[1];
    int kend = cctk_lsh[2] - nghostzones[2];
  */

  double dxi = 1.0/dX;
  double dyi = 1.0/dY;
  double dzi = 1.0/dZ;

  /*
  //Following lines needed since nghostzones[0] = ORDER, and 
  //   not ORDER-1 in axisymmetry 
  //   (so that rotation can be done on multiprocessor runs)
  if(Symmetry==4) {
    istart--;
    iend++;
  }
  */

  if(Symmetry==4) {
    jstart = 0;
    jend = cctk_lsh[1];
    jstart++;
    jend--;
  }

  if (flux_direction==1) {
    if (Symmetry==AXISYM) {
      for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
	int index   = CCTK_GFINDEX3D(cctkGH,i,  j,k);
	int indexp1 = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
	
	rho_star_rhs[index] += (frm[index]-frm[indexp1]) * dxi / X[index];
	e_rhs[index] += (fem[index] - fem[indexp1]) * dxi / X[index];
      }
    } else {
#pragma omp parallel for
      for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
	int index   = CCTK_GFINDEX3D(cctkGH,i,  j,k);
	int indexp1 = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
	
	rho_star_rhs[index] += (frm[index]-frm[indexp1]) * dxi;
	e_rhs[index] += (fem[index] - fem[indexp1]) * dxi;
      }
    }
  } else if (flux_direction==2) {
    if (Symmetry != AXISYM) {
#pragma omp parallel for
      for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
 	int index   = CCTK_GFINDEX3D(cctkGH,i,j,  k);
	int indexp1 = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
	rho_star_rhs[index] += (frm[index]-frm[indexp1]) * dyi;
	e_rhs[index] += (fem[index] - fem[indexp1]) * dyi;
      }
    }
  } else {
#pragma omp parallel for
    for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
      int index   = CCTK_GFINDEX3D(cctkGH,i,j,k  );
      int indexp1 = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
      rho_star_rhs[index] += (frm[index]-frm[indexp1]) * dzi;
      e_rhs[index] += (fem[index] - fem[indexp1]) * dzi;
    }
  }
}

extern "C" void CCTK_FCALL CCTK_FNAME(advect_rho_e_cpp)
  (int *flux_direction, const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry,
   double *rho_star_rhs,double *e_rhs,double *frm,double *fem,
   double *dX,double *dY,double *dZ,double *X,double *Y,double *Z) 
{
  advect_rho_e_cpp(*flux_direction,*cctkGH,cctk_lsh, nghostzones, *Symmetry,
		   rho_star_rhs,e_rhs,frm,fem,
		   *dX,*dY,*dZ,X,Y,Z);
}
