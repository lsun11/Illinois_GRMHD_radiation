//-----------------------------------------------------------------------------
// Advect rho_star and e
//-----------------------------------------------------------------------------

#include "math.h"
#include "cctk.h"
#include <stdio.h>

extern "C" void CCTK_FCALL CCTK_FNAME(advect_rhoYe_cpp)
  (int *flux_direction, const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry,
   double *rhoYe_rhs,double *frYem, 
   double *dX,double *dY,double *dZ,double *X,double *Y,double *Z, double *drhoYe_m_x, double *drhoYe_m_xp1);

extern "C" void advect_rhoYe_cpp(int flux_direction,const cGH *cctkGH,int *cctk_lsh,int *nghostzones,int Symmetry,
				 double *rhoYe_rhs,double *frYem, 
				 double dX,double dY,double dZ,double *X,double *Y,double *Z, double *drhoYe_m_x, double *drhoYe_m_xp1) {
  int AXISYM = 4;

  /* Set up variables used in the grid loop for the physical grid points */
  int istart = 1;
  int jstart = 1;
  int kstart = 1;
  int iend = cctk_lsh[0]-1;
  int jend = cctk_lsh[1]-1;
  int kend = cctk_lsh[2]-1;

  double dxi = 1.0/dX;
  double dyi = 1.0/dY;
  double dzi = 1.0/dZ;

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
	
	rhoYe_rhs[index] += (frYem[index]-frYem[indexp1]) * dxi / X[index];

      }
    } else {
#pragma omp parallel for
      for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
	    int index   = CCTK_GFINDEX3D(cctkGH,i,  j,k);
	    int indexp1 = CCTK_GFINDEX3D(cctkGH,i+1,j,k);	
	    rhoYe_rhs[index] += (frYem[index]-frYem[indexp1]) * dxi;
	    
	    drhoYe_m_x[index] = frYem[index];
	    drhoYe_m_xp1[index] = frYem[indexp1];
      }
    }
  } else if (flux_direction==2) {
    if (Symmetry != AXISYM) {
#pragma omp parallel for
      for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
	    int index   = CCTK_GFINDEX3D(cctkGH,i,j,  k);
	    int indexp1 = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
	    rhoYe_rhs[index] += (frYem[index]-frYem[indexp1]) * dyi;

	  }
    } else if (flux_direction==3){
#pragma omp parallel for
      for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
	    int index   = CCTK_GFINDEX3D(cctkGH,i,j,k  );
	    int indexp1 = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
	    rhoYe_rhs[index] += (frYem[index]-frYem[indexp1]) * dzi;
    
      }
    }
  }
}

extern "C" void CCTK_FCALL CCTK_FNAME(advect_rhoYe_cpp)
  (int *flux_direction, const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry,
   double *rhoYe_rhs, double *frYem, 
   double *dX,double *dY,double *dZ,double *X,double *Y,double *Z, double *drhoYe_m_x, double *drhoYe_m_xp1) 
{
  advect_rhoYe_cpp(*flux_direction,*cctkGH,cctk_lsh, nghostzones, *Symmetry,
		   rhoYe_rhs, frYem, 
		   *dX,*dY,*dZ,X,Y,Z,drhoYe_m_x, drhoYe_m_xp1);
}
