//---------------------------
// Implement 2D HLL flux 
// [see Del Zanna, Bucciantini & Londrillo A&A 400, 397 (2003), Eq. (44)]
//---------------------------
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "cctk.h"

extern "C" void CCTK_FCALL flux_hll2d_cpp_
  (const cGH **cctkGH,int *ext,double *cmax1, 
   double *cmin1, double *cmax2, double *cmin2,
   double *frr, double *frl, double *flr,
   double *fll, double *Ur_Ul1, double *Ur_Ul2, double *fhll, 
   int *use_central_scheme_instead_of_hll);

extern "C" void flux_hll2d_cpp(const cGH *cctkGH,int *ext,double *cmax1,
			    	double *cmin1, double *cmax2, double *cmin2, 
				double *frr, double *frl, double *flr, 
				double *fll, double *Ur_Ul1, double *Ur_Ul2, 
				double *fhll, int use_central_scheme_instead_of_hll)
{

  int istart = 0;
  int jstart = 0;
  int kstart = 0;
  int iend = ext[0];
  int jend = ext[1];
  int kend = ext[2];

#pragma omp parallel for
    for(int k=kstart;k<kend;k++)
      for(int j=jstart;j<jend;j++)
	for(int i=istart;i<iend;i++) {
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

          double cmax_1 = cmax1[index];
	  double cmax_2 = cmax2[index];
	  double cmin_1 = cmin1[index];
	  double cmin_2 = cmin2[index];
          double frr_ijk = frr[index];
	  double frl_ijk = frl[index];
	  double flr_ijk = flr[index];
	  double fll_ijk = fll[index];
	  double Ur_Ul1_ijk = Ur_Ul1[index];
	  double Ur_Ul2_ijk = Ur_Ul2[index];

          if (use_central_scheme_instead_of_hll==1) {
             double clax1, clax2;
             clax1 = cmax_1; 
             if (cmin_1 > clax1) clax1 = cmin_1; 
             cmax_1 = clax1; cmin_1 = clax1;
             clax2 = cmax_2; 
             if (cmin_2 > clax2) clax2 = cmin_2; 
             cmax_2 = clax2; cmin_2 = clax2;
          }
	  
	  fhll[index] = (cmax_1*cmax_2*fll_ijk + cmax_1*cmin_2*flr_ijk 
			    + cmin_1*cmax_2*frl_ijk + cmin_1*cmin_2*frr_ijk) 
			    /( (cmax_1+cmin_1)*(cmax_2+cmin_2) ) 
			    - cmax_1*cmin_1*(Ur_Ul2_ijk)/(cmax_1+cmin_1) 
			    + cmax_2*cmin_2*(Ur_Ul1_ijk)/(cmax_2+cmin_2);
	  if(isnan(fhll[index]))
	    {
	      printf("Inside flux_hll2d.C, fhll is nan, cmax_1=%e,cmax_2=%e,cmin_1=%e,cmin_2=%e", cmax_1, cmax_2, cmin_1, cmin_2);
	    }
	}	  
}

extern "C" void CCTK_FCALL flux_hll2d_cpp_
  (const cGH **cctkGH,int *ext,double *cmax1,
   double *cmin1, double *cmax2, double *cmin2,
   double *frr, double *frl, double *flr,
   double *fll, double *Ur_Ul1, double *Ur_Ul2, double *fhll, 
   int *use_central_scheme_instead_of_hll)
{
  flux_hll2d_cpp(*cctkGH,ext,cmax1,cmin1,cmax2,cmin2,frr,frl,flr,fll,
		 Ur_Ul1,Ur_Ul2, fhll, *use_central_scheme_instead_of_hll);
}
