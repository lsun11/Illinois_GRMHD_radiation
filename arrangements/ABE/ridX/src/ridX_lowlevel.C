
#define INDEXVAL(i, j, k)    ((i) + cctk_lsh[0] *		\
			      ((j) + cctk_lsh[1] * (k)))

#include <stdio.h>
#include "cctk.h"
#include <math.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"
#include "nrheaders.h"

extern "C" void CCTK_FCALL CCTK_FNAME(ridX_lowlevel)
  (const cGH **cctkGH, int *cctk_lsh, double &dx, double &dy, double &dz,
   double *Bx, double *By,double *Bz,double *monopole);

/*
  ridX Monopole Cleaner
  
  Solves underdetermined linear system of equations A x = b, for one possible x, using
  Numerical Recipes' SVD (Singular Value Decomposition) routine.

  Here, x_i denotes the perturbation of a particular B-field component (x,y, or z)
  at a particular location in space (i,j,k).  After all the perturbations are applied,
  we should have no monopoles on the grid, except at the boundary points.

*/

extern "C" void ridX_lowlevel(const cGH *cctkGH, int *cctk_lsh, double &dx, double &dy, double &dz,
			      double *Bx, double *By,double *Bz,double *monopole) {
  
  DECLARE_CCTK_PARAMETERS;

  /* Set up variables used in the grid loop for the physical grid points */
  int istart = 1;
  int jstart = 1;
  int kstart = 1;
  int iend = cctk_lsh[0];
  int jend = cctk_lsh[1];
  int kend = cctk_lsh[2];

  printf("hi dxyz = %e %e %e\n",dx,dy,dz);

  int *hasmpole = (int *)malloc(sizeof(int)*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]);
  int *mustfixB = (int *)malloc(sizeof(int)*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]);

  for(int k=0;k<cctk_lsh[2];k++)
    for(int j=0;j<cctk_lsh[1];j++)
      for(int i=0;i<cctk_lsh[0];i++) {
	hasmpole[INDEXVAL(i,j,k)] = 1;
	mustfixB[INDEXVAL(i,j,k)] = 1;
      }

  for(int k=kstart;k<kend;k++)
    for(int j=jstart;j<jend;j++)
      for(int i=istart;i<iend;i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	if(fabs(monopole[index])<1e-15) {
	  //If there is NO monopole at i,j,k, we 
	  // a) Zero out that element of our monopole tracking array:
	  hasmpole[INDEXVAL(i,j,k)] = 0;

	  // b) To preserve divB=0 at this non-monopole point, we mark 
	  //    the magnetic fields in its stencil un-repairable:
	  mustfixB[INDEXVAL(i,j,k)]=0;
	  mustfixB[INDEXVAL(i-1,j,k)]=0;
	  mustfixB[INDEXVAL(i-1,j-1,k)]=0;
	  mustfixB[INDEXVAL(i-1,j,k-1)]=0;
	  
	  mustfixB[INDEXVAL(i,j-1,k)]=0;
	  mustfixB[INDEXVAL(i,j-1,k-1)]=0;
	  mustfixB[INDEXVAL(i,j,k-1)]=0;
	  mustfixB[INDEXVAL(i-1,j-1,k-1)]=0;
	}
      }

  // Next, we count the total number of points where a monopole exists (nummonopoles),
  // and the total number of magnetic field components that are modifyable.
  //
  // In addition, we set an array monopoleindx[] that numbers each monopole, which is
  // a necessary step for filling the full monopole array:
  int *monopoleindx = (int *)malloc(sizeof(int)*(cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]+1));
  int nummonopoles=0,numBfixes=0;
  for(int k=kstart;k<kend;k++)
    for(int j=jstart;j<jend;j++)
      for(int i=istart;i<iend;i++) {
	if(mustfixB[INDEXVAL(i,j,k)]==1) numBfixes+=3;
	if(hasmpole[INDEXVAL(i,j,k)]==1) { 
	  nummonopoles++;
	  monopoleindx[INDEXVAL(i,j,k)] = nummonopoles;
	}
      }

  printf("numB = %d\t nummonopoles = %d\n",numBfixes,nummonopoles);

  // bigmatrix is the big MxN matrix -- where:
  // M=number of equations (i.e., LHSs of divB equations)
  // N=number of unknowns (i.e., number of fixable B-field components over entire grid)
  double **bigmatrix=dmatrix_ridX(1,nummonopoles,1,numBfixes);

  //First initialize values of bigmatrix to zero:
  for(int j=1;j<=numBfixes;j++)
    for(int i=1;i<=nummonopoles;i++) {
      bigmatrix[i][j] = 0.0;
    }

  // Next, we go to each point where the Bi's are fixable.  Each Bi that needs to be fixed
  // will cause 8 monopoles to be created (since 8 monopoles have a stencil that depends on 
  // Bi at this point), except maybe when Bi is located at a boundary.
  int whichB=1;
  for(int k=kstart;k<kend;k++)
    for(int j=jstart;j<jend;j++)
      for(int i=istart;i<iend;i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	if(mustfixB[INDEXVAL(i,j,k)]==1) {
	  for(int w=1;w<=3;w++) {
	    //Remember convention: bigmatrix[row][col]
	    //The idea is that each row corresponds to a monopole equation,
	    //  and each column a location,component where the i'th component of B^i field must be fixed

	    //monopole at i,j,k:
	    bigmatrix[monopoleindx[INDEXVAL(i,j,k)]][whichB] = 1;

	    //monopole at i+1,j,k:
	    if(i<iend-1) {
	      bigmatrix[monopoleindx[INDEXVAL(i+1,j,k)]][whichB] = 1;
	      if(w==1) bigmatrix[monopoleindx[INDEXVAL(i+1,j,k)]][whichB]*=-1;
	      //printf("%d %d %d %d %d\n",i,j,k,whichB,monopoleindx[INDEXVAL(i+1,j,k)]);
	    }

	    //monopole at i,j+1,k:
	    if(j<jend-1) {
	      bigmatrix[monopoleindx[INDEXVAL(i,j+1,k)]][whichB] = 1;
	      if(w==2) bigmatrix[monopoleindx[INDEXVAL(i,j+1,k)]][whichB]*=-1;
	    }

	    //monopole at i+1,j+1,k:
	    if(i<iend-1 && j<jend-1) {
	      bigmatrix[monopoleindx[INDEXVAL(i+1,j+1,k)]][whichB] = 1;
	      if(w==1 || w==2) bigmatrix[monopoleindx[INDEXVAL(i+1,j+1,k)]][whichB]*=-1;
	    }

	    //monopole at i,j,k+1:
	    if(k<kend-1) {
	      bigmatrix[monopoleindx[INDEXVAL(i,j,k+1)]][whichB] = 1;
	      if(w==3) bigmatrix[monopoleindx[INDEXVAL(i,j,k+1)]][whichB]*=-1;
	    }

	    //monopole at i+1,j,k+1:
	    if(i<iend-1 && k<kend-1) {
	      bigmatrix[monopoleindx[INDEXVAL(i+1,j,k+1)]][whichB] = 1;
	      if(w==1 || w==3) bigmatrix[monopoleindx[INDEXVAL(i+1,j,k+1)]][whichB]*=-1;
	    }

	    //monopole at i,j+1,k+1:
	    if(j<jend-1 && k<kend-1) {
	      bigmatrix[monopoleindx[INDEXVAL(i,j+1,k+1)]][whichB] = 1;
	      if(w==2 || w==3) bigmatrix[monopoleindx[INDEXVAL(i,j+1,k+1)]][whichB]*=-1;
	    }

	    //monopole at i+1,j+1,k+1:
	    if(i<iend-1 && j<jend-1 && k<kend-1) {
	      bigmatrix[monopoleindx[INDEXVAL(i+1,j+1,k+1)]][whichB] = -1;
	    }

	    whichB++;
	  }
	}
      }

  //Next, we set up variables for the SVD algorithms (see Numerical Recipes for details)
  double *wvec,**vvv;
  wvec=dvector_ridX(1,numBfixes);
  vvv=dmatrix_ridX(1,numBfixes,1,numBfixes);
  int tmp2=NR_END;

  printf("BEFORE DSVDCMP NR_END=%d\n",tmp2);

  // Perform the SVD... This is an N^3 operation, so it will take a LONG time (~5 minutes) for
  // a ~1500x3000 matrix.
  dsvdcmp_ridX(bigmatrix,nummonopoles,numBfixes,wvec,vvv);
 
  // Next, we zero out very small diagonal elements given by the wvec vector.  This is necessary
  // for getting a reasonable solution vector, according to Numerical Recipes.
  for(int i=1;i<=numBfixes;i++) {
    if(fabs(wvec[i])<1e-10) wvec[i]=0;
    printf("%d\t%e\n",i,wvec[i]);
  }

  // Finally, we use the SVD'ed bigmatrix and obtain a solution vector, x.  But first, we set
  // b=bbb, as in A x = b.
  double *bbb=dvector_ridX(1,numBfixes);
  double *x=dvector_ridX(1,numBfixes);

  nummonopoles=0;
  for(int k=kstart;k<kend;k++)
    for(int j=jstart;j<jend;j++)
      for(int i=istart;i<iend;i++) {
	if(hasmpole[INDEXVAL(i,j,k)]==1) { 
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  nummonopoles++;
	  // Yuk Tung: Please double-check the coefficient here (-4.0*dx).
	  bbb[nummonopoles] = -monopole[index]*4.0*dx;
	}
      }

  dsvbksb_ridX(bigmatrix,wvec,vvv,nummonopoles,numBfixes,bbb,x);

  // Now that we have the solution vector x, we can check for monopoles!
  // First we calculate how many monopoles before applying the solution vector x.
  int newnumms=0;
  for(int k=kstart;k<kend-1;k++)
    for(int j=jstart;j<jend-1;j++)
      for(int i=istart;i<iend-1;i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	double monopooo=0.5*(0.5*(1.0/dx*( 
					  (Bx[index] - Bx[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]) + 
					  (Bx[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] - Bx[CCTK_GFINDEX3D(cctkGH,i-1,j-1,k)]) + 
					  (Bx[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] - Bx[CCTK_GFINDEX3D(cctkGH,i-1,j,k-1)]) + 
					  (Bx[CCTK_GFINDEX3D(cctkGH,i,j-1,k-1)] - Bx[CCTK_GFINDEX3D(cctkGH,i-1,j-1,k-1)]) 
					  ))) + 
	  0.5*(0.5*(1.0/dy*( 
			    (By[index] - By[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]) + 
			    (By[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] - By[CCTK_GFINDEX3D(cctkGH,i-1,j-1,k)]) + 
			    (By[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] - By[CCTK_GFINDEX3D(cctkGH,i,j-1,k-1)]) + 
			    (By[CCTK_GFINDEX3D(cctkGH,i-1,j,k-1)] - By[CCTK_GFINDEX3D(cctkGH,i-1,j-1,k-1)]) 
			    ))) + 
	  0.5*(0.5*(1.0/dz*( 
			    (Bz[index] - Bz[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]) + 
			    (Bz[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] - Bz[CCTK_GFINDEX3D(cctkGH,i,j-1,k-1)]) + 
			    (Bz[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] - Bz[CCTK_GFINDEX3D(cctkGH,i-1,j,k-1)]) + 
			    (Bz[CCTK_GFINDEX3D(cctkGH,i-1,j-1,k)] - Bz[CCTK_GFINDEX3D(cctkGH,i-1,j-1,k-1)]) 
			    )));
	
	
	if(fabs(monopooo)>1e-14) {
	  newnumms++;
	}
      }

  printf("BEFORE FIX: total num of monopoles remaining: %d\n",newnumms);

  // We apply the solution vector x.
  whichB=1;
  for(int k=kstart;k<kend;k++)
    for(int j=jstart;j<jend;j++)
      for(int i=istart;i<iend;i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	if(mustfixB[INDEXVAL(i,j,k)]==1) {
	  printf("%d\t%d\t%d\t%e %e %e\n",i,j,k,x[whichB],x[whichB+1],x[whichB+2]);
	  
	  Bx[index] -= x[whichB]*0.25/dx; whichB++;
	  By[index] -= x[whichB]*0.25/dy; whichB++;
	  Bz[index] -= x[whichB]*0.25/dz; whichB++;
	}
      }  

  // Finally, we count the total number of monopoles remaining.
  newnumms=0;
  for(int k=kstart;k<kend-1;k++)
    for(int j=jstart;j<jend-1;j++)
      for(int i=istart;i<iend-1;i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	double monopooo=0.5*(0.5*(1.0/dx*( 
					  (Bx[index] - Bx[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]) + 
					  (Bx[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] - Bx[CCTK_GFINDEX3D(cctkGH,i-1,j-1,k)]) + 
					  (Bx[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] - Bx[CCTK_GFINDEX3D(cctkGH,i-1,j,k-1)]) + 
					  (Bx[CCTK_GFINDEX3D(cctkGH,i,j-1,k-1)] - Bx[CCTK_GFINDEX3D(cctkGH,i-1,j-1,k-1)]) 
					  ))) + 
	  0.5*(0.5*(1.0/dy*( 
			    (By[index] - By[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]) + 
			    (By[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] - By[CCTK_GFINDEX3D(cctkGH,i-1,j-1,k)]) + 
			    (By[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] - By[CCTK_GFINDEX3D(cctkGH,i,j-1,k-1)]) + 
			    (By[CCTK_GFINDEX3D(cctkGH,i-1,j,k-1)] - By[CCTK_GFINDEX3D(cctkGH,i-1,j-1,k-1)]) 
			    ))) + 
	  0.5*(0.5*(1.0/dz*( 
			    (Bz[index] - Bz[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]) + 
			    (Bz[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] - Bz[CCTK_GFINDEX3D(cctkGH,i,j-1,k-1)]) + 
			    (Bz[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] - Bz[CCTK_GFINDEX3D(cctkGH,i-1,j,k-1)]) + 
			    (Bz[CCTK_GFINDEX3D(cctkGH,i-1,j-1,k)] - Bz[CCTK_GFINDEX3D(cctkGH,i-1,j-1,k-1)]) 
			    )));
	
	
	if(fabs(monopooo)>1e-14) {
	  //printf("%d\t%d\t%d\t%e\n",i,j,k,monopooo);
	  newnumms++;
	}
      }

  printf("AFTER FIX: total num of monopoles remaining: %d\n",newnumms);
  exit(1);

  free_dmatrix_ridX(bigmatrix,1,nummonopoles,1,nummonopoles);
  free(hasmpole);
  free(mustfixB);
  free(monopoleindx);
  free(hasmpole);
  free(mustfixB);

}

extern "C" void CCTK_FCALL CCTK_FNAME(ridX_lowlevel)
  (const cGH **cctkGH, int *cctk_lsh, double &dx, double &dy, double &dz,
   double *Bx, double *By,double *Bz,double *monopole)
{
  ridX_lowlevel(*cctkGH,cctk_lsh, dx, dy, dz,
                Bx, By,Bz,monopole);
}


#include "nr_lowlevel_functs.h"
