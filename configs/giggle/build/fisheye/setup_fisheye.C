//----------------------------------------------------------------------------
//
// $Id: $
//
//----------------------------------------------------------------------------
//
// Set up fisheye coords
//
//----------------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <fstream>

#include "Symmetry.h"

void MultFishEyePhysicalRadius(int n,double *ai, double *r0i, double *si,double r, double& R);
void MultFishEyedRdr(int n,double *ai, double *r0i, double *si,double r, double& dR);
void MultFishEyed2Rdr2(int n,double *ai, double *r0i, double *si,double r, double& ddR);

void LogCoordPhysicalRadius(double r_min, double r, double& R);
void LogCoorddRdr(double r_min, double r, double& dR);
void LogCoordd2Rdr2(double r_min, double r, double& ddR);

static char *rcsid = "$Meow...$";
CCTK_FILEVERSION(Setup_Fisheye_Coords)

//
// Fisheye transformation
//
extern "C" void Setup_Fisheye_Coords(CCTK_ARGUMENTS) 
{
  using namespace std;

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int ext[3];

  double aFE_n;

  int nFE = 0;
  double *aFE,*r0FE,*sFE;

  ifstream file2;
  int log_transform_flag = 0;
  double r_min;

  ext[0] = cctk_lsh[0];
  ext[1] = cctk_lsh[1];
  ext[2] = cctk_lsh[2];

  if (fisheye_enable==0) {
    cout << "    No Fisheye " << endl;
    aFE = new double[nFE+1];
    r0FE = new double[nFE+1];
    sFE = new double[nFE+1];
    aFE[0] = 1.0;
  }
  else {
    ifstream file;
    file.open("RadialCoordinate_Input");
    if (!file) {
       cerr << "Can't open RadialCoordinate_Input" << endl;
       exit(1);
    }
    char buf[100],c;
    file.get(buf,100,'='); file.get(c); file >> nFE;
    aFE = new double[nFE+1];
    r0FE = new double[nFE+1];
    sFE = new double[nFE+1];

    aFE[0]=0;r0FE[0]=0;sFE[0]=0;

    file.get(buf,100,'='); file.get(c); file >> aFE[0];
    if(nFE>0) {
      for(int i=1; i<=nFE; i++) file >> aFE[i];
      file.get(buf,100,'='); file.get(c); file >> r0FE[1];
      for(int i=2; i<=nFE; i++) file >> r0FE[i];
      file.get(buf,100,'='); file.get(c); file >> sFE[1];
      for(int i=2; i<=nFE; i++) file >> sFE[i];
    }
    file.close();
    cout << "     Using fisheye transformation, n = "
	 << nFE << endl;
    cout << "     a = ";
    for(int i=0; i<=nFE; i++) cout << aFE[i] << "  ";
    cout << endl;
    cout << "     r0 = ";
    for(int i=1; i<=nFE; i++) cout << r0FE[i] << "  ";
    cout << endl;
    cout << "     s = ";
    for(int i=1; i<=nFE; i++) cout << sFE[i] << "  ";
    cout << endl;
  }
  aFE_n = aFE[nFE];

  //
  // Log radius transformation.  Do not use with fisheye.
  //

  file2.open( "LogRadialCoordinate_Input");
  if (!file2) {
    cout << "    No Log_radius file found" << endl;
    r_min = 0.0;
    if(fisheye_enable==0) {
#pragma omp parallel for
      for(int k=0; k<ext[2]; k++)
	for(int j=0; j<ext[1]; j++)
	  for(int i=0; i<ext[0]; i++)
	    {	  
	      int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	      PhysicalRadius[index] = r[index];
	      RadiusDerivative[index] = 1.0;  
	      RadiusDerivative2[index] = 0.0;
	    }
      delete[] aFE; delete[] r0FE; delete[] sFE;
      cout << "Finished with Fisheye" << endl;
      return;
    }
  } else {
    if(fisheye_enable==0) {
      printf("ERROR.  I found a LogRadialCoordinate_Input file in this directory, but you have fisheye_enable==0\n");
      exit(1);
    }
    log_transform_flag = 1;
    char buf[100],c;
    file2.get(buf,100,'='); file2.get(c); file2 >> r_min;
  }

#pragma omp parallel for
  for(int k=0; k<ext[2]; k++)
    for(int j=0; j<ext[1]; j++)
      for(int i=0; i<ext[0]; i++)
	{	  
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  if (log_transform_flag==1) {
	    LogCoordPhysicalRadius(r_min,r[index],PhysicalRadius[index]);
	    LogCoorddRdr(r_min,r[index],RadiusDerivative[index]);
	    LogCoordd2Rdr2(r_min,r[index],RadiusDerivative2[index]);
	  } else {
	    if(i==0 && j==0 && k==0) printf("hi2. %d %d %d\n",i,j,k);
	    MultFishEyePhysicalRadius(nFE,aFE,r0FE,sFE,r[index],PhysicalRadius[index]);
	    MultFishEyedRdr(nFE,aFE,r0FE,sFE,r[index],RadiusDerivative[index]);
	    MultFishEyed2Rdr2(nFE,aFE,r0FE,sFE,r[index],RadiusDerivative2[index]);
	  }
	}

  printf("Setting up fisheye in region where dX=%e\n",CCTK_DELTA_SPACE(0));
  delete[] aFE; delete[] r0FE; delete[] sFE;
  
  return;
}

/************************************************/
/*  MULTIPLE FISHEYE TRANSFORMATION COORDINATES */
/************************************************/
void MultFishEyePhysicalRadius(int n,double *ai, double *r0i, double *si,double r, double& R) {
  double *a = new double[n+1];
  double *r0 = new double[n+1];
  double *s = new double[n+1];
  double *kappa = new double[n+1];
  a[0] = ai[0];
  for(int i=1; i<=n; i++) {
    a[i] = ai[i];
    r0[i] = r0i[i];
    s[i] = si[i];
    kappa[i] = (a[i-1]-a[i])*s[i]/(2.0*tanh(r0[i]/s[i]));
  }
  R = a[n]*r;
  for (int i=1; i<=n; i++) {
    if (r > r0[i]) {
      R += kappa[i] * ( 2.0*r0[i]/s[i] + 
			log( (1.0+exp(-2.0*(r+r0[i])/s[i])) /
			     (1.0+exp(-2.0*(r-r0[i])/s[i])) ) );
    } else {
      R += kappa[i] * ( 2.0*r/s[i] + 
			log( (1.0+exp(-2.0*(r+r0[i])/s[i])) /
			     (1.0+exp(-2.0*(r0[i]-r)/s[i])) ) );
    }
  }
  //printf("hifisha: %e\n",R);
  delete[] a; delete[] r0;delete[] s;delete[] kappa;
}

void MultFishEyedRdr(int n,double *ai, double *r0i, double *si,double r, double& dR) {
  double *a = new double[n+1];
  double *r0 = new double[n+1];
  double *s = new double[n+1];
  double *kappa = new double[n+1];
  a[0] = ai[0];
  for(int i=1; i<=n; i++) {
    a[i] = ai[i];
    r0[i] = r0i[i];
    s[i] = si[i];
    kappa[i] = (a[i-1]-a[i])*s[i]/(2.0*tanh(r0[i]/s[i]));
  }
  dR = a[n];
  for(int i=1; i<=n; i++) {
    if (r < r0[i]) {
      dR += (a[i-1]-a[i]) * ( 1.0 + 2.0*exp(-2.0*r0[i]/s[i]) + 
			      exp(-4.0*r0[i]/s[i]) ) /
	( 1.0 + exp(-4.0*r0[i]/s[i]) + 
	  exp(-2.0*(r0[i]-r)/s[i]) + 
	  exp(-2.0*(r+r0[i])/s[i]) );
    } else {
      dR += (a[i-1]-a[i]) * exp(-2.0*(r-r0[i])/s[i]) * 
	( 1.0 + exp(-4.0*r0[i]/s[i]) + 
	  2.0*exp(-2.0*r0[i]/s[i]) ) /
	( 1.0 + exp(-4.0*r/s[i]) + 
	  exp(-2.0*(r-r0[i])/s[i]) + 
	  exp(-2.0*(r+r0[i])/s[i]) );	
    }
  }
  delete[] a; delete[] r0;delete[] s;delete[] kappa;
}

void MultFishEyed2Rdr2(int n,double *ai, double *r0i, double *si,double r, double& ddR) {
  double *a = new double[n+1];
  double *r0 = new double[n+1];
  double * s = new double[n+1];
  double *kappa = new double[n+1];
  a[0] = ai[0];
  for(int i=1; i<=n; i++) {
    a[i] = ai[i];
    r0[i] = r0i[i];
    s[i] = si[i];
    kappa[i] = (a[i-1]-a[i])*s[i]/(2.0*tanh(r0[i]/s[i]));
  }
  ddR = 0.0;
  for(int i=1; i<=n; i++) {
    if (r < r0[i]) {
      ddR += 2.0*(a[i-1]-a[i])/s[i] * ( 1.0 + 2.0*exp(-2.0*r0[i]/s[i]) +
					exp(-4.0*r0[i]/s[i]) ) /
	( 1.0 + exp(-4.0*r0[i]/s[i]) +
	  exp(-2.0*(r0[i]-r)/s[i]) +
	  exp(-2.0*(r+r0[i])/s[i]) ) / 
	( 1.0 + exp(-4.0*r0[i]/s[i]) +
	  exp(-2.0*(r0[i]-r)/s[i]) +
	  exp(-2.0*(r+r0[i])/s[i]) ) * 
	( -exp(-2.0*(r0[i]-r)/s[i]) 
	  + exp(-2.0*(r+r0[i])/s[i]) );
    } else {
      ddR += 2.0*(a[i-1]-a[i])/s[i] * exp(-2.0*(r-r0[i])/s[i]) *
	( 1.0 + exp(-4.0*r0[i]/s[i]) +
	  2.0*exp(-2.0*r0[i]/s[i]) ) /
	( 1.0 + exp(-4.0*r/s[i]) +
	  exp(-2.0*(r-r0[i])/s[i]) +
	  exp(-2.0*(r+r0[i])/s[i]) ) * 
	( ( exp(-2.0*(r-r0[i])/s[i]) + 
	    exp(-2.0*(r+r0[i])/s[i]) + 
	    2.0*exp(-4.0*r/s[i]) ) / 
	  ( 1.0 + exp(-4.0*r/s[i]) +
	    exp(-2.0*(r-r0[i])/s[i]) +
	    exp(-2.0*(r+r0[i])/s[i]) ) - 1.0 );
    }
  }
  delete[] a; delete[] r0;delete[] s;delete[] kappa;
}


/********************/
/*  LOG COORDINATES */
/********************/

void LogCoordPhysicalRadius(double r_min, double r, double& R) {
  R = exp(r-r_min);
}

void LogCoorddRdr(double r_min, double r, double& dR) {
  dR = exp(r-r_min);
}

void LogCoordd2Rdr2(double r_min, double r, double& ddR) {
  ddR = exp(r-r_min);
}
