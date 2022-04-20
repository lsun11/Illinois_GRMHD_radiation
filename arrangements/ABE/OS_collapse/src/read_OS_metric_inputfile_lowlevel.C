//-----------------------------------------------------------------------
// $Id$
//-----------------------------------------------------------------------
// Read binary files and do fancy things with them...
//-----------------------------------------------------------------------
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
#include "cctk_Parameters.h"
#include "nr.h"
#include "nrutil.h"

using namespace std;

extern "C" void CCTK_FCALL CCTK_FNAME(read_OS_metric_inputfile)
  (const cGH **cctkGH,
   double *psi,double *lapm1,double *X,double *Y,double *Z,
   int *nx, int *ny, int *nz, double *coord, int *narr,double *R_edge);

extern "C" void read_OS_metric_inputfile(const cGH *cctkGH,
					 double *psi,double *lapm1,double *X,
					 double *Y,double *Z,
					 int nx, int ny, int nz,
					 double *coord, int narr, double R_edge) {
  
  char filename[100];
  double time;
  int mark,i;
  int rpts;
  //read in from file
  double *rad,*A,*r_areal,*rho_aux,*Sr,*Srr,*lapse,*shift;
  //derivs obtained from spline
  double *ddA,*ddlapse, *ddrad;
  //radius
  double riso;
  //interpolated values
  double A_interp,lapse_interp;
  double HALF=0.5,ONE=1.0,ZERO=0.0,TEN=10.0,f4o3=4.0/3.0,f3o4=3.0/4.0;
  ifstream data_file;
  data_file.open("data_r5");
  if(!data_file) {
    cerr << "\a Can't open " << filename << " for input." << endl;
    exit(1);
  }
  data_file.read((char*)&rpts,sizeof(int));
  //Set up vectors
  rad = dvector(1,rpts);
  A = dvector(1,rpts);
  rho_aux = dvector(1,rpts);
  Sr = dvector(1,rpts);
  Srr = dvector(1,rpts);
  lapse = dvector(1,rpts);
  shift = dvector(1,rpts);
  r_areal = dvector(1,rpts);
  ddA = dvector(1,rpts);
  ddlapse = dvector(1,rpts);
  ddrad = dvector(1,rpts);
 //read in initial data
  data_file.read((char*)&time,sizeof(double));
  for(i=1;i<=rpts;i++) {
    data_file.read((char*)&rad[i],sizeof(double));
    data_file.read((char*)&lapse[i],sizeof(double));
    data_file.read((char*)&shift[i],sizeof(double));
    data_file.read((char*)&rho_aux[i],sizeof(double));
    data_file.read((char*)&Sr[i],sizeof(double));
    data_file.read((char*)&Srr[i],sizeof(double));
    data_file.read((char*)&A[i],sizeof(double));  
  }
  data_file.read((char*)&mark,sizeof(int));

  for(i=1;i<=rpts;i++) {
    r_areal[i] = A[i]*rad[i];
  }
  //spline
  spline(rad,A,rpts,0.0,1.e31,ddA);
  spline(rad,lapse,rpts,0.0,1.e31,ddlapse);
  spline(rad,rad,rpts,0.0,1.e31,ddrad);
  //
  // Now that we've read the matter in, we want to ditch the exterior:
  // it's zero and splining through the sudden drop-off causes severe
  // problems.
  //          
  
  //int nx = ext[0];
  //int ny = ext[1];
  //int nz = ext[2];

  //particle tracer stuff
  for (int index=0;index<=narr-1;index++) {
    coord[index + 0*narr] = ZERO;
    splint(r_areal,rad,ddrad,rpts,(index+1)*(0.97*R_edge/narr),&coord[index + 1*narr]); 
    coord[index + 2*narr] = ZERO;
    coord[index + 3*narr] = ZERO;
  } 
  
  for (int k = 0; k < nz; k++) for (int j = 0; j < ny; j++) for (int i = 0; i < nx; i++) {
    int vindex = CCTK_GFINDEX3D(cctkGH,i,j,k);
    //physical isotropic radius
    riso=sqrt(X[vindex]*X[vindex]+Y[vindex]*Y[vindex]+Z[vindex]*Z[vindex]);
    
    //do interpolations
    splint(rad,A,ddA,rpts,riso,&A_interp);
    splint(rad,lapse,ddlapse,rpts,riso,&lapse_interp);
    
    lapm1[vindex] = lapse_interp-ONE;
    psi[vindex] = sqrt(A_interp);
    //   printf("lapm1[vindex]: %10.6e\n",lapm1[vindex]);
  }
  data_file.close();  
}

extern "C" void CCTK_FCALL CCTK_FNAME(read_OS_metric_inputfile)
  (const cGH **cctkGH,
   double *psi,double *lapm1,double *X,double *Y,double *Z,
   int *nx,int *ny, int *nz, double *coord, int *narr,double *R_edge)
{  
  read_OS_metric_inputfile(*cctkGH,psi,lapm1,X,Y,Z,*nx,*ny,*nz,
			   coord,*narr,*R_edge);
}

