#include <cmath>            
#include <iostream>         
#include <iomanip>          
#include <sstream>
#include <string>           
#include <fstream>          
#include <stdlib.h>         


#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"
#include "cctk_FortranString.h"
using namespace std;
static char *rcsid="$mew. $";

// Achtung, you need to restructure the function calls starting with rtbis.  This is confusing.

void MultFishEyePhysicalRadius(int n,double *ai, double *r0i, double *si,double r, double& R);

void LogCoordPhysicalRadius(double r_min, double r, double& R);

double physical_radius_function(int n, double *ai, double *r0i, double *si, double R_input, 
				double x_in) {
  double R;
  
  MultFishEyePhysicalRadius(n,ai,r0i,si,x_in,R);
  return (R-R_input);
}

double physical_radius_function_log(double r_min, double R_input, double x_in) {
  double R;
  LogCoordPhysicalRadius(r_min,x_in,R);
  return (R-R_input);
}

double rtbis(double (*func)(int, double *, double *, double *, double, double),
	     int n, double *ai, double *r0i, double *si,
	     double x1, 
	     double x2, double xacc, double R_input);

double rtbis_log(double (*func)(double, double, double), double r_min, double x1, double x2, 
		 double xacc, double R_input);


// Note that grid_rmin and grid_rmax are meant to be the smallest and largest coordinate
// radius values on the grid.  Hopefully, they will bracket the root we're looking for.

CCTK_FILEVERSION(Convert_Physical_Radius)

extern "C" CCTK_INT Convert_Physical_Radius(double R_input, double grid_rmin, double grid_rmax, int fisheye_enable, double &r_out) {

    // first you need to read in all the fisheye parameters and log parameters.
    if (fisheye_enable == 1) {

	double aFE_n;
	
	int nFE = 0;
	double *aFE,*r0FE,*sFE;
	
	ifstream file2;
	int log_transform_flag = 0;
	double r_min;
	
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
	
	aFE_n = aFE[nFE];
	
	//
	// Log radius transformation.  Do not use with fisheye.
	//
	
	file2.open( "LogRadialCoordinate_Input");
	if (!file2) {
	    cout << "    No Log_radius file found" << endl;
	    r_min = 0.0;
	} else {
	    log_transform_flag = 1;
	    char buf[100],c;
	    file2.get(buf,100,'='); file2.get(c); file2 >> r_min;
	}
	
	double xacc;
	xacc = 1.e-10;
	
	//	double r_out;
	if (log_transform_flag) {
	    r_out = rtbis_log(physical_radius_function_log,r_min,
			      grid_rmin,grid_rmax,xacc,R_input);
	} else {

	    r_out = rtbis(physical_radius_function, nFE, aFE, r0FE, sFE,
			  grid_rmin,grid_rmax,xacc,R_input);
	}
	//	return r_out;
	
    } else {  // no fisheye or log coordinate transform being used.
      //	return R_input;
      r_out = 0.0;
    }
    return 0;
}

extern "C" void CCTK_FCALL CCTK_FNAME(Convert_Physical_Radius) 
  (CCTK_INT *retval, double *R_input, double *grid_rmin, double *grid_rmax, int *fisheye_enable, double *R_output)
{
  *retval = Convert_Physical_Radius(*R_input,*grid_rmin,*grid_rmax,*fisheye_enable, *R_output);
}

#define JMAX 100
double rtbis(double (*func)(int, double *, double *, double *, double, double), 
	     int n,double *ai, double *r0i, double *si,
	     double x1, double x2, double xacc, double R_input)
{
        void nrerror(char error_text[]);
        int j;
        double dx,f,fmid,xmid,rtb;

        f=(*func)(n,ai,r0i,si,R_input,x1);
        fmid=(*func)(n,ai,r0i,si,R_input,x2);
        if (f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis");
        rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
        for (j=1;j<=JMAX;j++) {
                fmid=(*func)(n,ai,r0i,si,R_input,xmid=rtb+(dx *= 0.5));
                if (fmid <= 0.0) rtb=xmid;
                if (fabs(dx) < xacc || fmid == 0.0) return rtb;
        }
        nrerror("Too many bisections in rtbis");
        return 0.0;
}

double rtbis_log(double (*func)(double, double, double), double r_min,
		 double x1, double x2, double xacc, double R_input)
{
        void nrerror(char error_text[]);
        int j;
        double dx,f,fmid,xmid,rtb;

        f=(*func)(r_min,R_input,x1);
        fmid=(*func)(r_min,R_input,x2);
        if (f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis");
        rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
        for (j=1;j<=JMAX;j++) {
                fmid=(*func)(r_min,R_input,xmid=rtb+(dx *= 0.5));
                if (fmid <= 0.0) rtb=xmid;
                if (fabs(dx) < xacc || fmid == 0.0) return rtb;
        }
        nrerror("Too many bisections in rtbis");
        return 0.0;
}
#undef JMAX


