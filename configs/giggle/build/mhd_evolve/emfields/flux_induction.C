//---------------------------
// fij = v^i B^j - v^j B^i
//---------------------------
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "cctk.h"

#define SQR(x) ((x) * (x))
#define AXISYM 4

extern "C" void CCTK_FCALL flux_induction_cpp_
  (const cGH **cctkGH,int *ext,double *X,double *fij,
   double *vir,double *vjr,double *Bir,double *Bjr,double *vil,double *vjl,double *Bil,double *Bjl,
   double *cmax,double *cmin,int &pow_axi,int &Symmetry);

void flux_hll_cpp_finduction(int *ext,double &Ur,double &Ul,double &Fr,double &Fl,double &F,double &cmax,double &cmin);

extern "C" void flux_induction_cpp(const cGH *cctkGH,int *ext,double *X,double *fij,
				   double *vir,double *vjr,double *Bir,double *Bjr,double *vil,double *vjl,double *Bil,double *Bjl,
				   double *cmax,double *cmin,int &pow_axi,int &Symmetry) {

// Auxiliary arrays
// real*8, dimension(ext(1),ext(2),ext(3))        :: Fr,Fl, X_f
//

  int istart = 0;
  int jstart = 0;
  int kstart = 0;
  int iend = ext[0];
  int jend = ext[1];
  int kend = ext[2];

//Next apply the HLL approximate Riemann solver to the induction equation:
#pragma omp parallel for
    for(int k=kstart;k<kend;k++)
      for(int j=jstart;j<jend;j++)
	for(int i=istart;i<iend;i++) {
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  //printf("hi %d %d %d %e %e %e %e\n",i,j,k,vir[index],Bjr[index],vjr[index],Bir[index]);
	  double Fr = vir[index]*Bjr[index] - vjr[index]*Bir[index];
	  double Fl = vil[index]*Bjl[index] - vjl[index]*Bil[index];

	  /*
	  if(isnan(Fr) || isnan(Fl)) {
	    printf("BadFr: %d %d %d %e %e %e %e\n",i,j,k,vir[index],Bjr[index],vjr[index],Bir[index]);
	    exit(1);
	  }
	  */

	  flux_hll_cpp_finduction(ext,Bjr[index],Bjl[index],Fr,Fl,fij[index],cmax[index],cmin[index]);
	}	  

  if (Symmetry==AXISYM && pow_axi != 0) {
    double dX2 = 0.5*(X[CCTK_GFINDEX3D(cctkGH,1,0,0)]-X[CCTK_GFINDEX3D(cctkGH,0,0,0)]);
    
#pragma omp parallel for
    for(int k=kstart;k<kend;k++)
      for(int j=jstart;j<jend;j++)
	for(int i=istart;i<iend;i++) {
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  double X_f;
	  if (pow_axi == 1) {
	    X_f = X[index] - dX2;
	  } else if (pow_axi == 2) {
	    //pow_axi == 2 for fzx//
	    X_f = X[index];
	  } else {
	    printf("BAD POW_AXI!\n"); exit(1);
	  }
	  fij[index] *= X_f;
	}
  }
}

extern "C" void CCTK_FCALL flux_induction_cpp_
  (const cGH **cctkGH,int *ext,double *X,double *fij,
   double *vir,double *vjr,double *Bir,double *Bjr,double *vil,double *vjl,double *Bil,double *Bjl,
   double *cmax,double *cmin,int &pow_axi,int &Symmetry)
{
  flux_induction_cpp
    (*cctkGH,ext,X,fij,
     vir,vjr,Bir,Bjr,vil,vjl,Bil,Bjl,
     cmax,cmin,pow_axi,Symmetry);
}

void flux_hll_cpp_finduction(int *ext,double &Ur,double &Ul,double &Fr,double &Fl,double &F,double &cmax,double &cmin) {
  F = (cmin*Fr + cmax*Fl - cmin*cmax*(Ur-Ul) )/(cmax + cmin);
  //We find that cmax==cmin==0 when we are in the ghostzone sometimes when B fields are being evolved.
  //   In the above line of code, this means that F=NaN, which leads to instability.
  //if(cmax+cmin==0) F=0;
}
