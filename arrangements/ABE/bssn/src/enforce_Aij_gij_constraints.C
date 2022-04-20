#define Fm1o3 -0.3333333333333333333333333333333
#define TWO 2.0


#include <stdio.h>
#include <stdlib.h>
#include "cctk.h"
#include <math.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

extern "C" void CCTK_FCALL CCTK_FNAME(enforce_Aij_gij_constraints)
  (const cGH **cctkGH,int *cctk_lsh,
   double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,
   double *Axx, double *Axy, double *Axz, double *Ayy, double *Ayz, double *Azz,
   double *x, double *y, double *z);

void enforce_Aij_gij_constraints(const cGH *cctkGH, int *cctk_lsh,
				 double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,
				 double *Axx, double *Axy, double *Axz, double *Ayy, double *Ayz, double *Azz,
				 double *x, double *y, double *z) {

  double dX = x[CCTK_GFINDEX3D(cctkGH,1,0,0)] - x[CCTK_GFINDEX3D(cctkGH,0,0,0)];

#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++)
    for(int j=0;j<cctk_lsh[1];j++)
      for(int i=0;i<cctk_lsh[0];i++) {

        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        double gxxL = gxx[index];
        double gxyL = gxy[index];
        double gxzL = gxz[index];
        double gyyL = gyy[index];
        double gyzL = gyz[index];
        double gzzL = gzz[index];

        double AxxL = Axx[index];
        double AxyL = Axy[index];
        double AxzL = Axz[index];
        double AyyL = Ayy[index];
        double AyzL = Ayz[index];
        double AzzL = Azz[index];

	double det =  gxxL * gyyL * gzzL + gxyL * gyzL * gxzL + gxzL * gxyL * gyzL 
	  - gxzL * gyyL * gxzL - gxyL * gxyL * gzzL - gxxL * gyzL * gyzL;
	double sign=1.0;
	if(det<=0.0) {
	  if(gxx[index]==0.0 && gxy[index]==0.0 && gxz[index]==0.0 && gyy[index]==0.0 && gyz[index]==0.0 && gzz[index]==0.0) {
	  } else {
	    printf("enforce_Aij_gij_constraints: WARNING: NEGATIVE DETERMINANT (det=%.15e) OF g_{ij} AT %d %d %d.\n dx = %e \txyz = %e %e %e \t gij_old=%e %e %e %e %e %e\n",
		   det,i,j,k,dX,x[index],y[index],z[index],gxx[index],gxy[index],gxz[index],gyy[index],gyz[index],gzz[index]);
	    sign=-1.0;
	  }
	  /*
	    // ZACH SAYS: This should happen RARELY.
	    double averagegxx=0;
	    double averagegxy=0;
	    double averagegxz=0;
	    double averagegyy=0;
	    double averagegyz=0;
	    double averagegzz=0;
	    int numpts=0;
	    if(i>0) { 
	      int idx=CCTK_GFINDEX3D(cctkGH,i-1,j,k); gxxL=gxx[idx]; gxyL=gxy[idx]; gxzL=gxz[idx]; gyyL=gyy[idx]; gyzL=gyz[idx]; gzzL=gzz[idx]; 
	      double detL=gxxL * gyyL * gzzL + gxyL * gyzL * gxzL + gxzL * gxyL * gyzL - gxzL * gyyL * gxzL - gxyL * gxyL * gzzL - gxxL * gyzL * gyzL;
	      if(detL>0) { averagegxx+=gxxL; averagegxy+=gxyL; averagegxz+=gxzL; averagegyy+=gyyL; averagegyz+=gyzL; averagegzz+=gzzL; numpts++; } }
	    if(j>0) { 
	      int idx=CCTK_GFINDEX3D(cctkGH,i,j-1,k); gxxL=gxx[idx]; gxyL=gxy[idx]; gxzL=gxz[idx]; gyyL=gyy[idx]; gyzL=gyz[idx]; gzzL=gzz[idx]; 
	      double detL=gxxL * gyyL * gzzL + gxyL * gyzL * gxzL + gxzL * gxyL * gyzL - gxzL * gyyL * gxzL - gxyL * gxyL * gzzL - gxxL * gyzL * gyzL;
	      if(detL>0) { averagegxx+=gxxL; averagegxy+=gxyL; averagegxz+=gxzL; averagegyy+=gyyL; averagegyz+=gyzL; averagegzz+=gzzL; numpts++; } }
	    if(k>0) { 
	      int idx=CCTK_GFINDEX3D(cctkGH,i,j,k-1); gxxL=gxx[idx]; gxyL=gxy[idx]; gxzL=gxz[idx]; gyyL=gyy[idx]; gyzL=gyz[idx]; gzzL=gzz[idx]; 
	      double detL=gxxL * gyyL * gzzL + gxyL * gyzL * gxzL + gxzL * gxyL * gyzL - gxzL * gyyL * gxzL - gxyL * gxyL * gzzL - gxxL * gyzL * gyzL;
	      if(detL>0) { averagegxx+=gxxL; averagegxy+=gxyL; averagegxz+=gxzL; averagegyy+=gyyL; averagegyz+=gyzL; averagegzz+=gzzL; numpts++; } }
	    if(i<cctk_lsh[0]-1) { 
	      int idx=CCTK_GFINDEX3D(cctkGH,i+1,j,k); gxxL=gxx[idx]; gxyL=gxy[idx]; gxzL=gxz[idx]; gyyL=gyy[idx]; gyzL=gyz[idx]; gzzL=gzz[idx]; 
	      double detL=gxxL * gyyL * gzzL + gxyL * gyzL * gxzL + gxzL * gxyL * gyzL - gxzL * gyyL * gxzL - gxyL * gxyL * gzzL - gxxL * gyzL * gyzL;
	      if(detL>0) { averagegxx+=gxxL; averagegxy+=gxyL; averagegxz+=gxzL; averagegyy+=gyyL; averagegyz+=gyzL; averagegzz+=gzzL; numpts++; } }
	    if(j<cctk_lsh[1]-1) { 
	      int idx=CCTK_GFINDEX3D(cctkGH,i,j+1,k); gxxL=gxx[idx]; gxyL=gxy[idx]; gxzL=gxz[idx]; gyyL=gyy[idx]; gyzL=gyz[idx]; gzzL=gzz[idx]; 
	      double detL=gxxL * gyyL * gzzL + gxyL * gyzL * gxzL + gxzL * gxyL * gyzL - gxzL * gyyL * gxzL - gxyL * gxyL * gzzL - gxxL * gyzL * gyzL;
	      if(detL>0) { averagegxx+=gxxL; averagegxy+=gxyL; averagegxz+=gxzL; averagegyy+=gyyL; averagegyz+=gyzL; averagegzz+=gzzL; numpts++; } }
	    if(k<cctk_lsh[2]-1) { 
	      int idx=CCTK_GFINDEX3D(cctkGH,i,j,k+1); gxxL=gxx[idx]; gxyL=gxy[idx]; gxzL=gxz[idx]; gyyL=gyy[idx]; gyzL=gyz[idx]; gzzL=gzz[idx]; 
	      double detL=gxxL * gyyL * gzzL + gxyL * gyzL * gxzL + gxzL * gxyL * gyzL - gxzL * gyyL * gxzL - gxyL * gxyL * gzzL - gxxL * gyzL * gyzL;
	      if(detL>0) { averagegxx+=gxxL; averagegxy+=gxyL; averagegxz+=gxzL; averagegyy+=gyyL; averagegyz+=gyzL; averagegzz+=gzzL; numpts++; } }
	    if(averagegxx==0) { printf("BAD POINTS EVERYWHERE!\n"); }
	    gxx[index] = averagegxx/((double)numpts);
	    gxy[index] = averagegxy/((double)numpts);
	    gxz[index] = averagegxz/((double)numpts);
	    gyy[index] = averagegyy/((double)numpts);
	    gyz[index] = averagegyz/((double)numpts);
	    gzz[index] = averagegzz/((double)numpts);

	    gxxL = gxx[index];
	    gxyL = gxy[index];
	    gxzL = gxz[index];
	    gyyL = gyy[index];
	    gyzL = gyz[index];
	    gzzL = gzz[index];

	    printf("FIXED A DET! old = %e\t",det);
	    det =  gxxL * gyyL * gzzL + gxyL * gyzL * gxzL + gxzL * gxyL * gyzL 
	      - gxzL * gyyL * gxzL - gxyL * gxyL * gzzL - gxxL * gyzL * gyzL;
	    printf("new = %e\n",det);
	    if(det<=0) { printf("STILL THE DET IS NEGATIVE! ARGH!\n"); sign=-1; }
	  */
	}


	double det_Fm1o3 = sign*pow(fabs(det),Fm1o3);

	gxxL = gxxL * det_Fm1o3;
	gxyL = gxyL * det_Fm1o3; 
	gxzL = gxzL * det_Fm1o3; 
	gyyL = gyyL * det_Fm1o3; 
	gyzL = gyzL * det_Fm1o3; 
	gzzL = gzzL * det_Fm1o3;

	double gupxx =   ( gyyL * gzzL - gyzL * gyzL ); 
	double gupxy = - ( gxyL * gzzL - gyzL * gxzL );
	double gupxz =   ( gxyL * gyzL - gyyL * gxzL ); 
	double gupyy =   ( gxxL * gzzL - gxzL * gxzL ); 
	double gupyz = - ( gxxL * gyzL - gxyL * gxzL ); 
	double gupzz =   ( gxxL * gyyL - gxyL * gxyL ); 
	
	double trA3 = Fm1o3*(gupxx * AxxL + gupyy * AyyL + gupzz * AzzL + 
	  TWO * ( gupxy * AxyL + gupxz * AxzL + gupyz * AyzL ));

	AxxL = AxxL + trA3 * gxxL;
	AxyL = AxyL + trA3 * gxyL;
	AxzL = AxzL + trA3 * gxzL;
	AyyL = AyyL + trA3 * gyyL;
	AyzL = AyzL + trA3 * gyzL;
	AzzL = AzzL + trA3 * gzzL;

	//printf("HI %d %d %d %e\n",i,j,k,AxxL);

	//FIXME: Read the printf statement below for description of problem.  Carpet may not obey checkpointing rules as rigidly as the old PUGH, it seems...
	if(gxx[index]==0.0 || gyy[index]==0.0 || gzz[index]==0.0) {
	  //if(i==0 && j==0 && k==0) {
	    //printf("enforce_Aij_gij_constraints: LOOKS LIKE YOU'RE RECOVERING FROM CHECKPOINT IN CARPET.  I'll ignore the gijs here.\n");
	    printf("ERROR.  One of the diagonal gij's equals zero!  ijk=(%d, %d, %d)  gxx=%e, gyy=%e, gzz=%e\n",i,j,k,gxx[index],gyy[index],gzz[index]);
	    //exit(0);
	    //}
	} else {
	  if(isnan(gxxL)) { printf("gxxL: enforce_Aij_gij_constraints: FOUND NAN AT %d/%d %d/%d %d/%d.\n dx = %e \txyz = %e %e %e \t gij_old=%e %e %e %e %e %e\n",i,cctk_lsh[0],j,cctk_lsh[1],k,cctk_lsh[2],dX,x[index],y[index],z[index],gxx[index],gxy[index],gxz[index],gyy[index],gyz[index],gzz[index]); } //exit(1);}
	  if(isnan(gxxL)) { printf("gxyL: enforce_Aij_gij_constraints: FOUND NAN AT %d/%d %d/%d %d/%d.\n dx = %e \txyz = %e %e %e \t gij_old=%e %e %e %e %e %e\n",i,cctk_lsh[0],j,cctk_lsh[1],k,cctk_lsh[2],dX,x[index],y[index],z[index],gxx[index],gxy[index],gxz[index],gyy[index],gyz[index],gzz[index]); } //exit(1);}
	  if(isnan(gxxL)) { printf("gxzL: enforce_Aij_gij_constraints: FOUND NAN AT %d/%d %d/%d %d/%d.\n dx = %e \txyz = %e %e %e \t gij_old=%e %e %e %e %e %e\n",i,cctk_lsh[0],j,cctk_lsh[1],k,cctk_lsh[2],dX,x[index],y[index],z[index],gxx[index],gxy[index],gxz[index],gyy[index],gyz[index],gzz[index]); } //exit(1);}
	  if(isnan(gxxL)) { printf("gyyL: enforce_Aij_gij_constraints: FOUND NAN AT %d/%d %d/%d %d/%d.\n dx = %e \txyz = %e %e %e \t gij_old=%e %e %e %e %e %e\n",i,cctk_lsh[0],j,cctk_lsh[1],k,cctk_lsh[2],dX,x[index],y[index],z[index],gxx[index],gxy[index],gxz[index],gyy[index],gyz[index],gzz[index]); } //exit(1);}
	  if(isnan(gxxL)) { printf("gyzL: enforce_Aij_gij_constraints: FOUND NAN AT %d/%d %d/%d %d/%d.\n dx = %e \txyz = %e %e %e \t gij_old=%e %e %e %e %e %e\n",i,cctk_lsh[0],j,cctk_lsh[1],k,cctk_lsh[2],dX,x[index],y[index],z[index],gxx[index],gxy[index],gxz[index],gyy[index],gyz[index],gzz[index]); } //exit(1);}
	  if(isnan(gxxL)) { printf("gzzL: enforce_Aij_gij_constraints: FOUND NAN AT %d/%d %d/%d %d/%d.\n dx = %e \txyz = %e %e %e \t gij_old=%e %e %e %e %e %e\n",i,cctk_lsh[0],j,cctk_lsh[1],k,cctk_lsh[2],dX,x[index],y[index],z[index],gxx[index],gxy[index],gxz[index],gyy[index],gyz[index],gzz[index]); } //exit(1);}
	  if(isnan(gxxL) || isnan(gxyL) || isnan(gxzL) || isnan(gyyL) || isnan(gyzL) || isnan(gzzL)) {
	  } else {
	  gxx[index] = gxxL;
	  gxy[index] = gxyL;
	  gxz[index] = gxzL;
	  gyy[index] = gyyL;
	  gyz[index] = gyzL;
	  gzz[index] = gzzL;

	  Axx[index] = AxxL;
	  Axy[index] = AxyL;
	  Axz[index] = AxzL;
	  Ayy[index] = AyyL;
	  Ayz[index] = AyzL;
	  Azz[index] = AzzL;
	  }
	}

      }
}

extern "C" void CCTK_FCALL CCTK_FNAME(enforce_Aij_gij_constraints)
  (const cGH **cctkGH,int *cctk_lsh,
   double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,
   double *Axx, double *Axy, double *Axz, double *Ayy, double *Ayz, double *Azz,
   double *x, double *y, double *z)
{
  enforce_Aij_gij_constraints(*cctkGH, cctk_lsh, gxx, gxy, gxz, gyy, gyz, gzz,
			      Axx, Axy, Axz, Ayy, Ayz, Azz,x,y,z);
}
