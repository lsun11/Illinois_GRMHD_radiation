// Read data code for CTS WDNS initial data

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "utilities.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <time.h>
#include <ctype.h>
using namespace std;

void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
void polin2(double x1a[], double x2a[], double **ya, int m, int n, double x1, double x2, double *y, double *dy);
void polin3(double x1a[], double x2a[], double x3a[], double ***ya, int m, int n, int l, double x1, double x2, double x3, double *y, double *dy);
void interpolate(double *xco, double *yco, double *zco, double *var, int nx, int ny, int nz, double x1, double x2, double x3, double &val, int order);
void locate(double *xx, unsigned long n, double x, long& j);

int main() 
{
  using namespace std;
  
  double xmin = 0.0;
  double ymin = 0.0;
  double zmin = 0.0;
  
  double dx = 1.0;
  double dy = 1.0;
  double dz = 1.0;
  
  int nx = 100;
  int ny = 100;
  int nz = 100;
  int ntotal = nx*ny*nz;

  double x_arr[nx];
  double y_arr[ny];
  double z_arr[nz];
  double function_arr[nx*ny*nz];
    for(int k=0; k<nz; k++){
      z_arr[k]=zmin+dz*k;
      for(int j=0; j<ny; j++){
	y_arr[j]=ymin+dy*j;
	for(int i=0; i<nx; i++){
	  x_arr[i]=xmin+dx*i;
	  function_arr[i+nx*(j+ny*k)]=x_arr[i]*y_arr[j]*z_arr[k];
	  //function_arr[i+nx*(j+ny*k)]=5.0;
	}
      }
    }
  int order = 2;  
  
  double x, y, z;
  x=25.5;
  y=25.5;
  z=25.5;
    
  double function_val;
 
  int count=0;
	
  interpolate(x_arr,y_arr,z_arr,function_arr,nx,ny,nz,x,y,z,function_val,order);

  cout <<"function_val: " << function_val << endl;
}

// -------------------------------------------------------------------------------

void locate(double *xx, unsigned long n, double x, long& j)
{

// Adapted from numerical recipes for a zero-offset array which is assumed monotonic
// n here is the size of the array xx[0,1,...n-1]
// It returns j so that x[j]< x < x[j+1] and returns j=-1 or j=n-1 when x is out of range

	long ju,jm,jl;
	int ascnd;

	jl=-1;
	ju=n;
	ascnd=(xx[n-1] >= xx[0]);
	while (ju-jl > 1) {
		jm=(ju+jl) >> 1;
		if (x >= xx[jm] == ascnd)
			jl=jm;
		else
			ju=jm;
	}
	if (x == xx[0]) j=0;
	else if(x == xx[n-1]) j=n-2;
	else j=jl;
}

// -------------------------------------------------------------------------------


//===========================================================
// find corresponding i,j and k given the index indx
//===========================================================

void indx_to_ijk(int &i, int &j, int &k, int indx, int nx, int ny)
{
  int II=indx;
  int nxny = nx*ny;  
  k = II/nxny;
  int II2 = II - k*nxny;
  j = II2/nx;
  i = II2 - j*nx;
}

// -------------------------------------------------------------------------------





// -------------------------------------------------------------------------------


void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[1]);
	c=dvector(1,n);
	d=dvector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_dvector(d,1,n);
	free_dvector(c,1,n);
}
//#undef NRANSI

// -------------------------------------------------------------------------------




// -------------------------------------------------------------------------------
void polin2(double x1a[], double x2a[], double **ya, int m, int n, double x1,
	double x2, double *y, double *dy)
{
  double *ymtmp;

	ymtmp=dvector(1,m);
	for (int j=1;j<=m;j++) {
		polint(x2a,ya[j],n,x2,&ymtmp[j],dy);
	}
	polint(x1a,ymtmp,m,x1,y,dy);
	free_dvector(ymtmp,1,m);
}
//#undef NRANSI
// -------------------------------------------------------------------------------



// -------------------------------------------------------------------------------
void polin3(double x1a[], double x2a[], double x3a[], double ***ya, int m, int n, int l, double x1,
	double x2, double x3, double *y, double *dy)
{

  /* Provide a m x n x l block of coordinates containing the given point (x1,x2,x3)
     and also a pointer (ya) to this block and the routine will provide the 
     intepolated value y(x1,x2,x3) via (m-1)th-order polylomial interp. in the x1
     direction, (n-1)th-order polylomial interp. in the x2 and 
     (l-1)th-order polylomial interp. in the x3 */


	// Call polin2 to do 2D interpolations along the line x1a=x1, x2a =y1 for all x3a in the block m x n x l
	// and store the result in the vector ytmp1, then use ytmp1 to do 1D  interpolation

  int i,j,k;
  double *ymtmp1,**ymtmp2;


	ymtmp1=dvector(1,l);
	ymtmp2=dmatrix(1,m,1,n);


	for (k=1;k<l+1;k++) {
	  for (i=1;i<m+1;i++){
	    for (j=1;j<n+1;j++){
	      ymtmp2[i][j] = ya[i][j][k];
	    }
	  }
	  polin2(x1a,x2a, ymtmp2 , m, n, x1, x2, &ymtmp1[k],dy);
	}
	polint(x3a, ymtmp1, l, x3 , y , dy);

	free_dvector(ymtmp1,1,m);
	free_dmatrix(ymtmp2,1,m,1,n);
}
// -------------------------------------------------------------------------------




//------------------------------------------------------------------------------------

void interpolate(double *xco, double *yco, double *zco, double *var, int nx, int ny, int nz, double x1, double x2, double x3, double &val, int order)
{

  int  ierr, i,j,k, nxny=nx*ny;
  long i5,j5,k5;

  double *x1a, *x2a, *x3a, y_val, dely, ***ya;
  int m1, n1, l1; // 2nd-order polynomial interpolation in all directions
  int ind;
  int *i1, *j1, *k1;


  if (order == 2){
    m1 = 3; n1 = 3; l1 = 3;
  }
  else if(order == 3){
    m1 = 4; n1 = 4; l1 = 4;
  }
  else if(order == 4){
    m1 = 5; n1 = 5; l1 = 5;
  }
  else if(order == 5){
    m1 = 6; n1 = 6; l1 = 6;
  }
  else{
    printf("Stopping: Selected interpolation order= %d not supported. Choose order = 2 or 3 or 4", order);
    exit(1);
  }

  x1a = dvector(1,m1);
  x2a = dvector(1,n1);
  x3a = dvector(1,l1);
  ya = d3tensor(1,m1,1,n1,1,l1);


  i1 = ivector(0,m1-1);
  j1 = ivector(0,n1-1);
  k1 = ivector(0,l1-1);

  locate(xco, nx, x1, i5);
  locate(yco, ny, x2, j5);
  locate(zco, nz, x3, k5);


  // Second order interpolation
  if (order == 2){    
    if (i5 ==0){
      i1[0] = i5; i1[1] = i5 + 1; i1[2] = i5 + 2;
    }
    else{
      i1[0] = i5 - 1; i1[1] = i5; i1[2] = i5 + 1;
    }
    
    
    if (j5 ==0){
      j1[0] = j5; j1[1] = j5 + 1; j1[2] = j5 + 2;
    }
    else{
      j1[0] = j5 -1 ; j1[1] = j5; j1[2] = j5 + 1;
    }
    
    
    if (k5 ==0){
      k1[0] = k5; k1[1] = k5 + 1; k1[2] = k5 + 2;
    }
    else{
      k1[0] = k5-1; k1[1] = k5; k1[2] = k5 + 1;
    }
    
    x1a[1] = xco[i1[0]]; x1a[2] = xco[i1[1]]; x1a[3] = xco[i1[2]];
    x2a[1] = yco[j1[0]]; x2a[2] = yco[j1[1]]; x2a[3] = yco[j1[2]];
    x3a[1] = zco[k1[0]]; x3a[2] = zco[k1[1]]; x3a[3] = zco[k1[2]];
  }

  // Third order interpolation
  if(order == 3){
    if (i5 ==0){
      i1[0] = i5; i1[1] = i5 + 1; i1[2] = i5 + 2; i1[3] = i5 + 3;
    }
    else if (i5 == nx -2){
      i1[0] = i5-2; i1[1] = i5 - 1; i1[2] = i5; i1[3] = i5 + 1;    
    }
    else{
      i1[0] = i5 - 1; i1[1] = i5; i1[2] = i5 + 1; i1[3] = i5 + 2;
    }
    
    
    if (j5 ==0){
      j1[0] = j5; j1[1] = j5 + 1; j1[2] = j5 + 2; j1[3] = j5 + 3;
    }
    else if (j5 == ny-2){
      j1[0] = j5-2; j1[1] = j5 - 1; j1[2] = j5; j1[3] = j5 + 1;
    }
    else{
      j1[0] = j5 -1 ; j1[1] = j5; j1[2] = j5 + 1; j1[3] = j5 + 2;
    }
    
    
    if (k5 ==0){
      k1[0] = k5; k1[1] = k5 + 1; k1[2] = k5 + 2; k1[3] = k5 + 3;
    }
    else if (k5 == nz - 2){
      k1[0] = k5-2; k1[1] = k5 - 1; k1[2] = k5; k1[3] = k5 + 1;
    }
    else{
      k1[0] = k5-1; k1[1] = k5; k1[2] = k5 + 1; k1[3] = k5 + 2;
    }
    
    x1a[1] = xco[i1[0]]; x1a[2] = xco[i1[1]]; x1a[3] = xco[i1[2]]; x1a[4] = xco[i1[3]];
    x2a[1] = yco[j1[0]]; x2a[2] = yco[j1[1]]; x2a[3] = yco[j1[2]]; x2a[4] = yco[j1[3]];
    x3a[1] = zco[k1[0]]; x3a[2] = zco[k1[1]]; x3a[3] = zco[k1[2]]; x3a[4] = zco[k1[3]];
  }


  // Fourth order interpolation
  if(order == 4){
    if (i5 ==0){
      i1[0] = i5; i1[1] = i5 + 1; i1[2] = i5 + 2; i1[3] = i5 + 3; i1[4] = i5 + 4;
    }
    else if (i5 == nx-2){
      i1[0] = i5-3; i1[1] = i5 -2 ; i1[2] = i5 -1 ; i1[3] = i5; i1[4] = i5 + 1;
    }
    else if (i5 == nx-3){
      i1[0] = i5-2; i1[1] = i5 -1 ; i1[2] = i5 ; i1[3] = i5+1; i1[4] = i5 + 2;
    }
    else{
      i1[0] = i5 - 1; i1[1] = i5; i1[2] = i5 + 1; i1[3] = i5 + 2; i1[4] = i5 + 3;
    }

    
    
    if (j5 ==0){
      j1[0] = j5; j1[1] = j5 + 1; j1[2] = j5 + 2; j1[3] = j5 + 3; j1[4] = j5 + 4;
    }
    else if (j5 == ny -2 ){
      j1[0] = j5-3; j1[1] = j5 - 2; j1[2] = j5 -1; j1[3] = j5 ; j1[4] = j5 + 1;
    }
    else if (j5 == ny -3 ){
      j1[0] = j5-2; j1[1] = j5 - 1; j1[2] = j5; j1[3] = j5+1 ; j1[4] = j5 + 2;
    }
    else{
      j1[0] = j5 -1 ; j1[1] = j5; j1[2] = j5 + 1; j1[3] = j5 + 2; j1[4] = j5 + 3;
    }
    
    
    if (k5 ==0){
      k1[0] = k5; k1[1] = k5 + 1; k1[2] = k5 + 2; k1[3] = k5 + 3; k1[4] = k5 + 4;
    }
    else if (k5 == nz-2){
      k1[0] = k5-3; k1[1] = k5 -2; k1[2] = k5 -1; k1[3] = k5; k1[4] = k5 + 1;
    }
    else if (k5 == nz-3){
      k1[0] = k5-2; k1[1] = k5 -1; k1[2] = k5; k1[3] = k5+1; k1[4] = k5 + 2;
    }
    else{
      k1[0] = k5-1; k1[1] = k5; k1[2] = k5 + 1; k1[3] = k5 + 2; k1[4] = k5 + 3;
    }
    
    x1a[1] = xco[i1[0]]; x1a[2] = xco[i1[1]]; x1a[3] = xco[i1[2]]; x1a[4] = xco[i1[3]]; x1a[5] = xco[i1[4]];
    x2a[1] = yco[j1[0]]; x2a[2] = yco[j1[1]]; x2a[3] = yco[j1[2]]; x2a[4] = yco[j1[3]]; x2a[5] = yco[j1[4]];
    x3a[1] = zco[k1[0]]; x3a[2] = zco[k1[1]]; x3a[3] = zco[k1[2]]; x3a[4] = zco[k1[3]]; x3a[5] = zco[k1[4]];
  }


  // Fifth order interpolation
  if(order == 5){
    if (i5 ==0){
      i1[0] = i5; i1[1] = i5 + 1; i1[2] = i5 + 2; i1[3] = i5 + 3; i1[4] = i5 + 4; i1[5] = i5 + 5;
    }
    else if (i5 == 1){
      i1[0] = i5 - 1; i1[1] = i5; i1[2] = i5 + 1; i1[3] = i5 + 2; i1[4] = i5 + 3; i1[5] = i5 + 4;
    }
    else if (i5 == nx-2){
      i1[0] = i5 - 4; i1[1] = i5-3; i1[2] = i5 -2; i1[3] = i5 -1; i1[4] = i5; i1[5] = i5 + 1;
    }
    else if (i5 == nx-3){
      i1[0] = i5 - 3; i1[1] = i5-2; i1[2] = i5 -1; i1[3] = i5; i1[4] = i5+1; i1[5] = i5 + 2;
    }
    else{
      i1[0] = i5 - 2; i1[1] = i5-1; i1[2] = i5; i1[3] = i5 + 1; i1[4] = i5 + 2; i1[5] = i5 + 3;
    }

    
    
    if (j5 ==0){
      j1[0] = j5; j1[1] = j5 + 1; j1[2] = j5 + 2; j1[3] = j5 + 3; j1[4] = j5 + 4; j1[5] = j5 + 5;
    }
    else if (j5 ==1){
      j1[0] = j5 -1 ; j1[1] = j5; j1[2] = j5 + 1; j1[3] = j5 + 2; j1[4] = j5 + 3; j1[5] = j5 + 4;
    }
    else if (j5 == ny-2){
      j1[0] = j5 -4 ; j1[1] = j5-3; j1[2] = j5 -2; j1[3] = j5 -1; j1[4] = j5; j1[5] = j5 + 1;
    }
    else if (j5 == ny-3){
      j1[0] = j5 -3 ; j1[1] = j5-2; j1[2] = j5 -1; j1[3] = j5; j1[4] = j5+1; j1[5] = j5 + 2;
    }
    else{
      j1[0] = j5 -2 ; j1[1] = j5-1; j1[2] = j5; j1[3] = j5 + 1; j1[4] = j5 + 2; j1[5] = j5 + 3;
    }
    
    
    if (k5 ==0){
      k1[0] = k5; k1[1] = k5 + 1; k1[2] = k5 + 2; k1[3] = k5 + 3; k1[4] = k5 + 4; k1[5] = k5 + 5;
    }
    else if (k5 ==1){
      k1[0] = k5-1; k1[1] = k5; k1[2] = k5 + 1; k1[3] = k5 + 2; k1[4] = k5 + 3;  k1[5] = k5 + 4;
    }
    else if (k5 == nz-2){
      k1[0] = k5-4; k1[1] = k5-3; k1[2] = k5 -2; k1[3] = k5 -1; k1[4] = k5;  k1[5] = k5 + 1;
    }
    else if (k5 == nz-3){
      k1[0] = k5-3; k1[1] = k5-2; k1[2] = k5 -1; k1[3] = k5; k1[4] = k5+1;  k1[5] = k5 + 2;
    }
    else{
      k1[0] = k5-2; k1[1] = k5-1; k1[2] = k5; k1[3] = k5 + 1; k1[4] = k5 + 2; k1[5] = k5 + 3;
    }
    
    x1a[1] = xco[i1[0]]; x1a[2] = xco[i1[1]]; x1a[3] = xco[i1[2]]; x1a[4] = xco[i1[3]]; x1a[5] = xco[i1[4]]; x1a[6] = xco[i1[5]];
    x2a[1] = yco[j1[0]]; x2a[2] = yco[j1[1]]; x2a[3] = yco[j1[2]]; x2a[4] = yco[j1[3]]; x2a[5] = yco[j1[4]]; x2a[6] = yco[j1[5]];
    x3a[1] = zco[k1[0]]; x3a[2] = zco[k1[1]]; x3a[3] = zco[k1[2]]; x3a[4] = zco[k1[3]]; x3a[5] = zco[k1[4]]; x3a[6] = zco[k1[5]];
  }

  

  double value;
  for (i = 0; i < m1; i++){
    for (j = 0; j < n1; j++){
      for (k = 0; k < l1; k++){
	ind = i1[i] + nx*j1[j] + nxny*k1[k];
	value = var[ind];
	ya[i+1][j+1][k+1] = value;
	//	if (fabs(x2-39.0625) < 1.e-10 && fabs(x1+176.5625) < 1.e-10 && fabs(x3-26.5625) < 1.e-10){
	//	  cout << "value[" << i << "," << j << "," << k << "]=" << value << ", ind=" << ind << endl; 
	//	}
      }
    }
  }

  polin3(x1a, x2a, x3a, ya,  m1, n1, l1, x1,  x2, x3, &y_val, &dely);

  val = y_val;


  free_dvector(x1a, 1,m1);
  free_dvector(x2a, 1,n1);
  free_dvector(x3a, 1,l1);
  free_d3tensor(ya, 1,m1,1,n1,1,l1);

  free_ivector(i1, 0,m1-1);
  free_ivector(j1, 0,n1-1);
  free_ivector(k1, 0,l1-1);
  
}
//------------------------------------------------------------------------------------

