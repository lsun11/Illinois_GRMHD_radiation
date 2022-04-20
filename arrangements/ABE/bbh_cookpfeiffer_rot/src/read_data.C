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

typedef struct LEVEL { 
  double *phi, *lapse, *betax, *betay, *betaz, *A11, *A12, *A13, *A22, *A23, *A33, 
    *rhoo, *ux, *uy, *uz, *x, *y, *z, xmin, xmax, ymin, ymax, zmin, zmax;
} level;


void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
void polin2(double x1a[], double x2a[], double **ya, int m, int n, double x1, double x2, double *y, double *dy);
void polin3(double x1a[], double x2a[], double x3a[], double ***ya, int m, int n, int l, double x1, double x2, double x3, double *y, double *dy);
void interpolate(double *xco, double *yco, double *zco, double *var, int nx, int ny, int nz, double x1, double x2, double x3, double &val, int order);
void locate(double *xx, unsigned long n, double x, long& j);

void read_inputfile_wdns(level &lv, int &nx, int &ny, int &nz,			 
			 double &an1, double &an2, double &an, double &kappa1, double &kappa2, double &kappa, double &rhoo1, double &rhoo2, 
			 double &M_ADM_surf, double &M_K_surf, double &J_ADM_surf, double &MNS_b, double &MWD_b, double &Omega, double &xrot,
			 int levelnumber, int print_params);



void output_asci_file(level lv, int nx_lv, int ny_lv, int nz_lv, int levelnumber);


int main(int argc, char** argv) 
{
  using namespace std;
  
  if (argv[1]==NULL || argv[12]==NULL) { 
    printf("Sorry.  I was expecting 12 command line arguments:\n");
    printf("Usage: ./CTS_output <xmin> <ymin> <zmin> <dx> <dy> <dz> <nx> <ny> <nz> <processor number> <output_asci> <interpolation order>\n");
    exit(1);
  }
  
  double xmin = strtod(argv[1],NULL) ;
  double ymin = strtod(argv[2],NULL) ;
  double zmin = strtod(argv[3],NULL) ;
  
  double dx = strtod(argv[4],NULL) ;
  double dy = strtod(argv[5],NULL) ;
  double dz = strtod(argv[6],NULL) ;
  
  int nx = (int)strtod(argv[7],NULL) ;
  int ny = (int)strtod(argv[8],NULL) ;
  int nz = (int)strtod(argv[9],NULL) ;
  int ntotal = nx*ny*nz;

  int output_asci = (int)strtod(argv[11],NULL) ;
  int order = (int)strtod(argv[12],NULL) ;  
  
  printf("\n");
  printf("Given Input parameters: %e %e %e %e %e %e %d %d %d %d %d\n",xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,output_asci,order);
  printf("xmax = %f, ymax =%f, zmax = %f\n", xmin+(nx-1)*dx, ymin+(ny-1)*dy, zmin+(nz-1)*dz);

  int nx_lv, ny_lv, nz_lv, i,j,k;
  double an1, an2, an, kappa1, kappa2, kappa, rhoo1, rhoo2; // equation of state parameters
  double M_ADM_surf, M_K_surf, J_ADM_surf, MNS_b, MWD_b, Omega, xrot;
  double x, y, z;
  level lv1, lv2, lv3, lv4, lv5;
  int procnumber = (int)strtod(argv[10],NULL) ;


  // Read WDNS data reading each level at a time
  read_inputfile_wdns(lv1, nx_lv, ny_lv, nz_lv,			 
		      an1, an2, an, kappa1, kappa2, kappa, rhoo1, rhoo2, 
		      M_ADM_surf, M_K_surf, J_ADM_surf, MNS_b, MWD_b, Omega, xrot, 1, 0); 

  read_inputfile_wdns(lv2, nx_lv, ny_lv, nz_lv,			 
		      an1, an2, an, kappa1, kappa2, kappa, rhoo1, rhoo2, 
		      M_ADM_surf, M_K_surf, J_ADM_surf, MNS_b, MWD_b, Omega, xrot, 2, 0); 

  read_inputfile_wdns(lv3, nx_lv, ny_lv, nz_lv,			 
		      an1, an2, an, kappa1, kappa2, kappa, rhoo1, rhoo2, 
		      M_ADM_surf, M_K_surf, J_ADM_surf, MNS_b, MWD_b, Omega, xrot, 3, 0); 

  read_inputfile_wdns(lv4, nx_lv, ny_lv, nz_lv,			 
		      an1, an2, an, kappa1, kappa2, kappa, rhoo1, rhoo2, 
		      M_ADM_surf, M_K_surf, J_ADM_surf, MNS_b, MWD_b, Omega, xrot, 4, 0); 

  read_inputfile_wdns(lv5, nx_lv, ny_lv, nz_lv,			 
		      an1, an2, an, kappa1, kappa2, kappa, rhoo1, rhoo2, 
		      M_ADM_surf, M_K_surf, J_ADM_surf, MNS_b, MWD_b, Omega, xrot, 5, 1); 


  // Set coordinates of boundaries of each level 
  lv1.xmin = lv1.x[0];  lv1.xmax = lv1.x[nx_lv-1];  lv1.ymin = lv1.y[0];  lv1.ymax = lv1.y[ny_lv-1];  lv1.zmin = lv1.z[0];  lv1.zmax = lv1.z[nz_lv-1];
  lv2.xmin = lv2.x[0];  lv2.xmax = lv2.x[nx_lv-1];  lv2.ymin = lv2.y[0];  lv2.ymax = lv2.y[ny_lv-1];  lv2.zmin = lv2.z[0];  lv2.zmax = lv2.z[nz_lv-1];
  lv3.xmin = lv3.x[0];  lv3.xmax = lv3.x[nx_lv-1];  lv3.ymin = lv3.y[0];  lv3.ymax = lv3.y[ny_lv-1];  lv3.zmin = lv3.z[0];  lv3.zmax = lv3.z[nz_lv-1];
  lv4.xmin = lv4.x[0];  lv4.xmax = lv4.x[nx_lv-1];  lv4.ymin = lv4.y[0];  lv4.ymax = lv4.y[ny_lv-1];  lv4.zmin = lv4.z[0];  lv4.zmax = lv4.z[nz_lv-1];
  lv5.xmin = lv5.x[0];  lv5.xmax = lv5.x[nx_lv-1];  lv5.ymin = lv5.y[0];  lv5.ymax = lv5.y[ny_lv-1];  lv5.zmin = lv5.z[0];  lv5.zmax = lv5.z[nz_lv-1];

  printf("\n");
  printf("Boundaries of grid levels from CTS data read\n");
  printf("Level 1\n");
  printf("lv1.xmin, lv1.xmax = %f %f\n", lv1.xmin, lv1.xmax);
  printf("lv1.ymin, lv1.ymax = %f %f\n", lv1.ymin, lv1.ymax);
  printf("lv1.zmin, lv1.zmax = %f %f\n", lv1.zmin, lv1.zmax);
  printf("Level 2\n");
  printf("lv2.xmin, lv2.xmax = %f %f\n", lv2.xmin, lv2.xmax);
  printf("lv2.ymin, lv2.ymax = %f %f\n", lv2.ymin, lv2.ymax);
  printf("lv2.zmin, lv2.zmax = %f %f\n", lv2.zmin, lv2.zmax);
  printf("Level 3\n");
  printf("lv3.xmin, lv3.xmax = %f %f\n", lv3.xmin, lv3.xmax);
  printf("lv3.ymin, lv3.ymax = %f %f\n", lv3.ymin, lv3.ymax);
  printf("lv3.zmin, lv3.zmax = %f %f\n", lv3.zmin, lv3.zmax);
  printf("Level 4\n");
  printf("lv4.xmin, lv4.xmax = %f %f\n", lv4.xmin, lv4.xmax);
  printf("lv4.ymin, lv4.ymax = %f %f\n", lv4.ymin, lv4.ymax);
  printf("lv4.zmin, lv4.zmax = %f %f\n", lv4.zmin, lv4.zmax);
  printf("Level 5\n");
  printf("lv5.xmin, lv5.xmax = %f %f\n", lv5.xmin, lv5.xmax);
  printf("lv5.ymin, lv5.ymax = %f %f\n", lv5.ymin, lv5.ymax);
  printf("lv5.zmin, lv5.zmax = %f %f\n", lv5.zmin, lv5.zmax);


  // Export variables to a file 
  ofstream outfilepsi; 
  std::ostringstream filenamepsi; 
  char fname_psi[100]; 


  // Set filename, using procnumber 
  /* 
     Note that procnumber below is an integer that uniquely specifies the chunk 
     of grid we want to output.  For an N processor run, with M levels, there 
     are N*M initial data files and unique procnumber's.  We have found that 
     this works best when setting up initial data on grids.  All the initial 
     data binary file output routine should need is: 
     xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,procnumber 
  */ 
  sprintf(fname_psi,"CTS_bin-proc%d.d",procnumber); 
  outfilepsi.open(fname_psi,ios::out | ios::ate); 
  
  /* 
     Output full dataset to file first.  After this we will write all the 
     parameters.  We do this so that, in case the file is not fully written, 
     our initial data thorn will give an error message when it does not find 
     parameters.  Otherwise, there is a chance that the initial data thorn will 
     not fill the appropriate arrays, and our evolution code will mysteriously 
     crash. 
     
  */

  double lapse_val, phi_val, betax_val, betay_val, betaz_val, A11_val, A12_val, A13_val, A22_val, A23_val, A33_val;
  double rhoo_val, ux_val, uy_val, uz_val;

  ofstream solfile;
  std::ostringstream asci_filename; 
  char asci_fname[100]; 
  sprintf(asci_fname,"inter_sol%d.dat",procnumber); 
  solfile.open(asci_fname); 


  outfilepsi.write((char *) &(nx),sizeof(int)); 
  outfilepsi.write((char *) &(ny),sizeof(int)); 
  outfilepsi.write((char *) &(nz),sizeof(int)); 
  outfilepsi.write((char *) &(ntotal),sizeof(int)); 


 
  //  solfile.open("inter_sol.dat");
  int count=0;

  for(k=0; k<nz; k++){
    z = zmin + k*dz;          
    //z = 0.;
    for(j=0; j<ny; j++){
      y = ymin + j*dy;
      //y = 0.;
      for(i=0; i<nx; i++){      
	x = xmin + i*dx;	
	if ( (x > lv5.xmin && x < lv5.xmax) && (y > lv5.ymin && y < lv5.ymax) && (z > lv5.zmin && z < lv5.zmax) ){ // Inside level 5
	  //	 count +=1;
	  interpolate(lv5.x, lv5.y, lv5.z, lv5.lapse, nx_lv, ny_lv, nz_lv, x, y, z, lapse_val, order);
	  interpolate(lv5.x, lv5.y, lv5.z, lv5.phi, nx_lv, ny_lv, nz_lv, x, y, z, phi_val, order);
	  interpolate(lv5.x, lv5.y, lv5.z, lv5.betax, nx_lv, ny_lv, nz_lv, x, y, z, betax_val, order);
	  interpolate(lv5.x, lv5.y, lv5.z, lv5.betay, nx_lv, ny_lv, nz_lv, x, y, z, betay_val, order);
	  interpolate(lv5.x, lv5.y, lv5.z, lv5.betaz, nx_lv, ny_lv, nz_lv, x, y, z, betaz_val, order);
	  interpolate(lv5.x, lv5.y, lv5.z, lv5.A11, nx_lv, ny_lv, nz_lv, x, y, z, A11_val, order);
	  interpolate(lv5.x, lv5.y, lv5.z, lv5.A12, nx_lv, ny_lv, nz_lv, x, y, z, A12_val, order);
	  interpolate(lv5.x, lv5.y, lv5.z, lv5.A13, nx_lv, ny_lv, nz_lv, x, y, z, A13_val, order);
	  interpolate(lv5.x, lv5.y, lv5.z, lv5.A22, nx_lv, ny_lv, nz_lv, x, y, z, A22_val, order);
	  interpolate(lv5.x, lv5.y, lv5.z, lv5.A23, nx_lv, ny_lv, nz_lv, x, y, z, A23_val, order);
	  interpolate(lv5.x, lv5.y, lv5.z, lv5.A33, nx_lv, ny_lv, nz_lv, x, y, z, A33_val, order);
	  interpolate(lv5.x, lv5.y, lv5.z, lv5.rhoo, nx_lv, ny_lv, nz_lv, x, y, z, rhoo_val, order);
	  interpolate(lv5.x, lv5.y, lv5.z, lv5.ux, nx_lv, ny_lv, nz_lv, x, y, z, ux_val, order);
	  interpolate(lv5.x, lv5.y, lv5.z, lv5.uy, nx_lv, ny_lv, nz_lv, x, y, z, uy_val, order);
	  interpolate(lv5.x, lv5.y, lv5.z, lv5.uz, nx_lv, ny_lv, nz_lv, x, y, z, uz_val, order);
	}       
	else  if ( (x > lv4.xmin && x < lv4.xmax) && (y > lv4.ymin && y < lv4.ymax) && (z > lv4.zmin && z < lv4.zmax) ){ // Inside level 5
	  //	 count +=1;
	  interpolate(lv4.x, lv4.y, lv4.z, lv4.lapse, nx_lv, ny_lv, nz_lv, x, y, z, lapse_val, order);
	  interpolate(lv4.x, lv4.y, lv4.z, lv4.phi, nx_lv, ny_lv, nz_lv, x, y, z, phi_val, order);
	  interpolate(lv4.x, lv4.y, lv4.z, lv4.betax, nx_lv, ny_lv, nz_lv, x, y, z, betax_val, order);
	  interpolate(lv4.x, lv4.y, lv4.z, lv4.betay, nx_lv, ny_lv, nz_lv, x, y, z, betay_val, order);
	  interpolate(lv4.x, lv4.y, lv4.z, lv4.betaz, nx_lv, ny_lv, nz_lv, x, y, z, betaz_val, order);
	  interpolate(lv4.x, lv4.y, lv4.z, lv4.A11, nx_lv, ny_lv, nz_lv, x, y, z, A11_val, order);
	  interpolate(lv4.x, lv4.y, lv4.z, lv4.A12, nx_lv, ny_lv, nz_lv, x, y, z, A12_val, order);
	  interpolate(lv4.x, lv4.y, lv4.z, lv4.A13, nx_lv, ny_lv, nz_lv, x, y, z, A13_val, order);
	  interpolate(lv4.x, lv4.y, lv4.z, lv4.A22, nx_lv, ny_lv, nz_lv, x, y, z, A22_val, order);
	  interpolate(lv4.x, lv4.y, lv4.z, lv4.A23, nx_lv, ny_lv, nz_lv, x, y, z, A23_val, order);
	  interpolate(lv4.x, lv4.y, lv4.z, lv4.A33, nx_lv, ny_lv, nz_lv, x, y, z, A33_val, order);
	  interpolate(lv4.x, lv4.y, lv4.z, lv4.rhoo, nx_lv, ny_lv, nz_lv, x, y, z, rhoo_val, order);
	  interpolate(lv4.x, lv4.y, lv4.z, lv4.ux, nx_lv, ny_lv, nz_lv, x, y, z, ux_val, order);
	  interpolate(lv4.x, lv4.y, lv4.z, lv4.uy, nx_lv, ny_lv, nz_lv, x, y, z, uy_val, order);
	  interpolate(lv4.x, lv4.y, lv4.z, lv4.uz, nx_lv, ny_lv, nz_lv, x, y, z, uz_val, order);
	}       
	else  if ( (x > lv3.xmin && x < lv3.xmax) && (y > lv3.ymin && y < lv3.ymax) && (z > lv3.zmin && z < lv3.zmax) ){ // Inside level 5
	  //	 count +=1;
	  interpolate(lv3.x, lv3.y, lv3.z, lv3.lapse, nx_lv, ny_lv, nz_lv, x, y, z, lapse_val, order);
	  interpolate(lv3.x, lv3.y, lv3.z, lv3.phi, nx_lv, ny_lv, nz_lv, x, y, z, phi_val, order);
	  interpolate(lv3.x, lv3.y, lv3.z, lv3.betax, nx_lv, ny_lv, nz_lv, x, y, z, betax_val, order);
	  interpolate(lv3.x, lv3.y, lv3.z, lv3.betay, nx_lv, ny_lv, nz_lv, x, y, z, betay_val, order);
	  interpolate(lv3.x, lv3.y, lv3.z, lv3.betaz, nx_lv, ny_lv, nz_lv, x, y, z, betaz_val, order);
	  interpolate(lv3.x, lv3.y, lv3.z, lv3.A11, nx_lv, ny_lv, nz_lv, x, y, z, A11_val, order);
	  interpolate(lv3.x, lv3.y, lv3.z, lv3.A12, nx_lv, ny_lv, nz_lv, x, y, z, A12_val, order);
	  interpolate(lv3.x, lv3.y, lv3.z, lv3.A13, nx_lv, ny_lv, nz_lv, x, y, z, A13_val, order);
	  interpolate(lv3.x, lv3.y, lv3.z, lv3.A22, nx_lv, ny_lv, nz_lv, x, y, z, A22_val, order);
	  interpolate(lv3.x, lv3.y, lv3.z, lv3.A23, nx_lv, ny_lv, nz_lv, x, y, z, A23_val, order);
	  interpolate(lv3.x, lv3.y, lv3.z, lv3.A33, nx_lv, ny_lv, nz_lv, x, y, z, A33_val, order);
	  interpolate(lv3.x, lv3.y, lv3.z, lv3.rhoo, nx_lv, ny_lv, nz_lv, x, y, z, rhoo_val, order);
	  interpolate(lv3.x, lv3.y, lv3.z, lv3.ux, nx_lv, ny_lv, nz_lv, x, y, z, ux_val, order);
	  interpolate(lv3.x, lv3.y, lv3.z, lv3.uy, nx_lv, ny_lv, nz_lv, x, y, z, uy_val, order);
	  interpolate(lv3.x, lv3.y, lv3.z, lv3.uz, nx_lv, ny_lv, nz_lv, x, y, z, uz_val, order);
	}       
	else  if ( (x > lv2.xmin && x < lv2.xmax) && (y > lv2.ymin && y < lv2.ymax) && (z > lv2.zmin && z < lv2.zmax) ){ // Inside level 5
	  //	 count +=1;
	  interpolate(lv2.x, lv2.y, lv2.z, lv2.lapse, nx_lv, ny_lv, nz_lv, x, y, z, lapse_val, order);
	  interpolate(lv2.x, lv2.y, lv2.z, lv2.phi, nx_lv, ny_lv, nz_lv, x, y, z, phi_val, order);
	  interpolate(lv2.x, lv2.y, lv2.z, lv2.betax, nx_lv, ny_lv, nz_lv, x, y, z, betax_val, order);
	  interpolate(lv2.x, lv2.y, lv2.z, lv2.betay, nx_lv, ny_lv, nz_lv, x, y, z, betay_val, order);
	  interpolate(lv2.x, lv2.y, lv2.z, lv2.betaz, nx_lv, ny_lv, nz_lv, x, y, z, betaz_val, order);
	  interpolate(lv2.x, lv2.y, lv2.z, lv2.A11, nx_lv, ny_lv, nz_lv, x, y, z, A11_val, order);
	  interpolate(lv2.x, lv2.y, lv2.z, lv2.A12, nx_lv, ny_lv, nz_lv, x, y, z, A12_val, order);
	  interpolate(lv2.x, lv2.y, lv2.z, lv2.A13, nx_lv, ny_lv, nz_lv, x, y, z, A13_val, order);
	  interpolate(lv2.x, lv2.y, lv2.z, lv2.A22, nx_lv, ny_lv, nz_lv, x, y, z, A22_val, order);
	  interpolate(lv2.x, lv2.y, lv2.z, lv2.A23, nx_lv, ny_lv, nz_lv, x, y, z, A23_val, order);
	  interpolate(lv2.x, lv2.y, lv2.z, lv2.A33, nx_lv, ny_lv, nz_lv, x, y, z, A33_val, order);
	  interpolate(lv2.x, lv2.y, lv2.z, lv2.rhoo, nx_lv, ny_lv, nz_lv, x, y, z, rhoo_val, order);
	  interpolate(lv2.x, lv2.y, lv2.z, lv2.ux, nx_lv, ny_lv, nz_lv, x, y, z, ux_val, order);
	  interpolate(lv2.x, lv2.y, lv2.z, lv2.uy, nx_lv, ny_lv, nz_lv, x, y, z, uy_val, order);
	  interpolate(lv2.x, lv2.y, lv2.z, lv2.uz, nx_lv, ny_lv, nz_lv, x, y, z, uz_val, order);
	}       
	else  if ( (x > lv1.xmin && x < lv1.xmax) && (y > lv1.ymin && y < lv1.ymax) && (z > lv1.zmin && z < lv1.zmax) ){ // Inside level 5
	  //	 count +=1;
	  interpolate(lv1.x, lv1.y, lv1.z, lv1.lapse, nx_lv, ny_lv, nz_lv, x, y, z, lapse_val, order);
	  interpolate(lv1.x, lv1.y, lv1.z, lv1.phi, nx_lv, ny_lv, nz_lv, x, y, z, phi_val, order);
	  interpolate(lv1.x, lv1.y, lv1.z, lv1.betax, nx_lv, ny_lv, nz_lv, x, y, z, betax_val, order);
	  interpolate(lv1.x, lv1.y, lv1.z, lv1.betay, nx_lv, ny_lv, nz_lv, x, y, z, betay_val, order);
	  interpolate(lv1.x, lv1.y, lv1.z, lv1.betaz, nx_lv, ny_lv, nz_lv, x, y, z, betaz_val, order);
	  interpolate(lv1.x, lv1.y, lv1.z, lv1.A11, nx_lv, ny_lv, nz_lv, x, y, z, A11_val, order);
	  interpolate(lv1.x, lv1.y, lv1.z, lv1.A12, nx_lv, ny_lv, nz_lv, x, y, z, A12_val, order);
	  interpolate(lv1.x, lv1.y, lv1.z, lv1.A13, nx_lv, ny_lv, nz_lv, x, y, z, A13_val, order);
	  interpolate(lv1.x, lv1.y, lv1.z, lv1.A22, nx_lv, ny_lv, nz_lv, x, y, z, A22_val, order);
	  interpolate(lv1.x, lv1.y, lv1.z, lv1.A23, nx_lv, ny_lv, nz_lv, x, y, z, A23_val, order);
	  interpolate(lv1.x, lv1.y, lv1.z, lv1.A33, nx_lv, ny_lv, nz_lv, x, y, z, A33_val, order);
	  interpolate(lv1.x, lv1.y, lv1.z, lv1.rhoo, nx_lv, ny_lv, nz_lv, x, y, z, rhoo_val, order);
	  interpolate(lv1.x, lv1.y, lv1.z, lv1.ux, nx_lv, ny_lv, nz_lv, x, y, z, ux_val, order);
	  interpolate(lv1.x, lv1.y, lv1.z, lv1.uy, nx_lv, ny_lv, nz_lv, x, y, z, uy_val, order);
	  interpolate(lv1.x, lv1.y, lv1.z, lv1.uz, nx_lv, ny_lv, nz_lv, x, y, z, uz_val, order);
	}       
	else{
	  printf("Stopping: Grid point (x,y,z) = (%f,%f,%f) is outside of CTS grid boundary!\n",x,y,z);
	  exit(1);
	}
	
	if (rhoo_val < 1.e-16){
	  rhoo_val = 0.; 
	  ux_val = 0.; 
	  uy_val = 0.; 
	  uz_val = 0.;
	}

	outfilepsi.write((char *) &phi_val,sizeof(double)); 
	outfilepsi.write((char *) &lapse_val,sizeof(double)); 
	outfilepsi.write((char *) &betax_val,sizeof(double)); 
	outfilepsi.write((char *) &betay_val,sizeof(double)); 
	outfilepsi.write((char *) &betaz_val,sizeof(double)); 
	outfilepsi.write((char *) &A11_val,sizeof(double)); 
	outfilepsi.write((char *) &A12_val,sizeof(double)); 
	outfilepsi.write((char *) &A13_val,sizeof(double)); 
	outfilepsi.write((char *) &A22_val,sizeof(double)); 
	outfilepsi.write((char *) &A23_val,sizeof(double)); 
	outfilepsi.write((char *) &A33_val,sizeof(double)); 
	outfilepsi.write((char *) &rhoo_val,sizeof(double)); 
	outfilepsi.write((char *) &ux_val,sizeof(double)); 
	outfilepsi.write((char *) &uy_val,sizeof(double)); 
	outfilepsi.write((char *) &uz_val,sizeof(double)); 
	

	if (phi_val-1. < 1.e-16 || lapse_val-1. > 1.e-16){
	  printf("WARNING!!!: BAD value for either psi or lapse = %f %f at x,y,z = %f %f %f\n", phi_val, lapse_val,x,y,z);
	  //	  exit(1);
	}
	
	//	if( fabs(z) <=2000.){
	solfile << x <<" " << y  << " " << z << " " <<  rhoo_val <<  " " << phi_val << " " << lapse_val << " " << betax_val << " " << betay_val << " " << betaz_val << " " << A11_val  << " " << A12_val << " " << A13_val << " " << A22_val << " " << A23_val << " " << A33_val << endl;
	//	}
      }
    }
  }
  // cout << count << endl;
  // cout << nx*ny*nz << endl;
  solfile.close();
  
  // Finally, output the binary parameters to a file, and close the file. 
  outfilepsi.write((char *) &(Omega),sizeof(double)); 
  outfilepsi.write((char *) &(xrot),sizeof(double)); 
  outfilepsi.write((char *) &(MNS_b),sizeof(double)); 
  outfilepsi.write((char *) &(MWD_b),sizeof(double)); 
  outfilepsi.write((char *) &(M_ADM_surf),sizeof(double)); 
  outfilepsi.write((char *) &(M_K_surf),sizeof(double)); 
  outfilepsi.write((char *) &(J_ADM_surf),sizeof(double)); 
  outfilepsi.write((char *) &(an1),sizeof(double)); 
  outfilepsi.write((char *) &(an2),sizeof(double)); 
  outfilepsi.write((char *) &(an),sizeof(double)); 
  outfilepsi.write((char *) &(kappa1),sizeof(double)); 
  outfilepsi.write((char *) &(kappa2),sizeof(double)); 
  outfilepsi.write((char *) &(kappa),sizeof(double)); 
  outfilepsi.write((char *) &(rhoo1),sizeof(double)); 
  outfilepsi.write((char *) &(rhoo2),sizeof(double)); 

      cout << "\n";
      cout<<"Energy and Angular momentum" <<endl;
      cout<<"MNS_baryon="<< MNS_b<<endl;
      cout<<"MWD_baryon="<< MWD_b<<endl;
      cout<<"ADM mass="<<M_ADM_surf<<endl;
      cout<<"Komar mass="<< M_K_surf <<endl;
      cout<<"Angular momentum ="<< J_ADM_surf <<endl;
      printf("\n");
      cout<<"EOS parameters" <<endl;
      cout<<"an1 ="<< an1 <<endl;
      cout<<"an2 ="<< an2 <<endl;
      cout<<"an3 ="<< an <<endl;
      cout<<"Gamma1 ="<< 1. + 1./an1 <<endl;
      cout<<"Gamma2 ="<< 1. + 1./an2 <<endl;
      cout<<"Gamma3 ="<< 1. + 1./an <<endl;
      cout<<"kappa1 ="<< kappa1 <<endl;
      cout<<"kappa2 ="<< kappa2 <<endl;
      cout<<"kappa3 ="<< kappa <<endl;
      cout<<"rhoo1 ="<< rhoo1 <<endl;
      cout<<"rhoo2 ="<< rhoo2 <<endl;


  
  outfilepsi.close(); 
  
  if (output_asci == 1){
    output_asci_file(lv1, nx_lv, ny_lv, nz_lv, 1);
    output_asci_file(lv2, nx_lv, ny_lv, nz_lv, 2);
    output_asci_file(lv3, nx_lv, ny_lv, nz_lv, 3);
    output_asci_file(lv4, nx_lv, ny_lv, nz_lv, 4);
    output_asci_file(lv5, nx_lv, ny_lv, nz_lv, 5);
  }
  
  delete [] lv1.phi; delete [] lv1.lapse; delete [] lv1.betax; delete [] lv1.betay; delete [] lv1.betaz;
  delete [] lv1.A11;  delete [] lv1.A12;  delete [] lv1.A13;  delete [] lv1.A22;  delete [] lv1.A23;  delete [] lv1.A33;
  delete [] lv1.rhoo;  delete [] lv1.ux;  delete [] lv1.uy;  delete [] lv1.uz;  delete [] lv1.x;  delete [] lv1.y;  delete [] lv1.z;
  
  delete [] lv2.phi; delete [] lv2.lapse; delete [] lv2.betax; delete [] lv2.betay; delete [] lv2.betaz;
  delete [] lv2.A11;  delete [] lv2.A12;  delete [] lv2.A13;  delete [] lv2.A22;  delete [] lv2.A23;  delete [] lv2.A33;
  delete [] lv2.rhoo;  delete [] lv2.ux;  delete [] lv2.uy;  delete [] lv2.uz;  delete [] lv2.x;  delete [] lv2.y;  delete [] lv2.z;
  
  delete [] lv3.phi; delete [] lv3.lapse; delete [] lv3.betax; delete [] lv3.betay; delete [] lv3.betaz;
  delete [] lv3.A11;  delete [] lv3.A12;  delete [] lv3.A13;  delete [] lv3.A22;  delete [] lv3.A23;  delete [] lv3.A33;
  delete [] lv3.rhoo;  delete [] lv3.ux;  delete [] lv3.uy;  delete [] lv3.uz;  delete [] lv3.x;  delete [] lv3.y;  delete [] lv3.z;
  
  delete [] lv4.phi; delete [] lv4.lapse; delete [] lv4.betax; delete [] lv4.betay; delete [] lv4.betaz;
  delete [] lv4.A11;  delete [] lv4.A12;  delete [] lv4.A13;  delete [] lv4.A22;  delete [] lv4.A23;  delete [] lv4.A33;
  delete [] lv4.rhoo;  delete [] lv4.ux;  delete [] lv4.uy;  delete [] lv4.uz;  delete [] lv4.x;  delete [] lv4.y;  delete [] lv4.z;
  
  delete [] lv5.phi; delete [] lv5.lapse; delete [] lv5.betax; delete [] lv5.betay; delete [] lv5.betaz;
  delete [] lv5.A11;  delete [] lv5.A12;  delete [] lv5.A13;  delete [] lv5.A22;  delete [] lv5.A23;  delete [] lv5.A33;
  delete [] lv5.rhoo;  delete [] lv5.ux;  delete [] lv5.uy;  delete [] lv5.uz;  delete [] lv5.x;  delete [] lv5.y;  delete [] lv5.z;
  
  
}



void output_asci_file(level lv, int nx_lv, int ny_lv, int nz_lv, int levelnumber)
{
  int i,j,k, ind;
  double lapse_val, phi_val, betax_val, betay_val, betaz_val, A11_val, A12_val, A13_val, A22_val, A23_val, A33_val;
  double rhoo_val, ux_val, uy_val, uz_val, x, y, z, u0;
  
  ofstream solfile;
  std::ostringstream asci_filename; 
  char asci_fname[100]; 
  sprintf(asci_fname,"asci_wdns%d.dat",levelnumber); 
  solfile.open(asci_fname); 

  for(i=0; i<nx_lv; i++){      
    for(j=0; j<ny_lv; j++){
      for(k=0; k<nz_lv; k++){

	ind = i + nx_lv*j + nx_lv*ny_lv*k;

	lapse_val = lv.lapse[ind];
	phi_val = lv.phi[ind];
	betax_val = lv.betax[ind];
	betay_val = lv.betay[ind];
	betaz_val = lv.betaz[ind];
	A11_val = lv.A11[ind];
	A12_val = lv.A12[ind];
	A13_val = lv.A13[ind];
	A22_val = lv.A22[ind];
	A23_val = lv.A23[ind];
	rhoo_val = lv.rhoo[ind];
	ux_val = lv.ux[ind];
	uy_val = lv.uy[ind];
	uz_val = lv.uz[ind];
	x = lv.x[i];
	y = lv.y[j];
	z = lv.z[k];
       
	u0=sqrt(1.0+pow(phi_val,-4)*(ux_val*ux_val+uy_val*uy_val+uz_val*uz_val))/(lapse_val);
	ux_val=pow(phi_val,-4)*ux_val/u0-betax_val;
	uy_val=pow(phi_val,-4)*uy_val/u0-betay_val;
	uz_val=pow(phi_val,-4)*uz_val/u0-betaz_val;

	//	if(fabs(y) < 10. && fabs(z) < 10.){
       solfile << x <<" " << y  << " " << z << " " <<  rhoo_val <<  " " << phi_val << " " << lapse_val << " " << betax_val << " " << betay_val << " " << betaz_val << " " << A11_val  << " " << A12_val << " " << A13_val << " " << A22_val << " " << A23_val << " " << A33_val << " " << ux_val << " " << uy_val << " " << uz_val << 
endl;
       //	}
      }
    }
  }

}




void read_inputfile_wdns(level &lv, int &nx, int &ny, int &nz,			 
			 double &an1, double &an2, double &an, double &kappa1, double &kappa2, double &kappa, double &rhoo1, double &rhoo2, 
			 double &M_ADM_surf, double &M_K_surf, double &J_ADM_surf, double &MNS_b, double &MWD_b, double &Omega, double &xrot,
			 int levelnumber, int print_params)
{

  char filename[100];
  sprintf(filename,"CTS_bin-proc%d.dat",levelnumber);


  ifstream infile1;
  infile1.open(filename);
  if(!infile1) {
    cerr << "\a Can't open " << filename << " for input." << endl;
    exit(1);
  }

  double doubledum, xmin, ymin, zmin, dx, dy, dz;
  int ntot,intdum;
  
    // Read  grid parameters
    infile1.read((char *) &nx, sizeof(int));
    infile1.read((char *) &ny, sizeof(int));
    infile1.read((char *) &nz, sizeof(int));
    infile1.read((char *) &xmin, sizeof(double));
    infile1.read((char *) &ymin, sizeof(double));
    infile1.read((char *) &zmin, sizeof(double));
    infile1.read((char *) &dx, sizeof(double));
    infile1.read((char *) &dy, sizeof(double));
    infile1.read((char *) &dz, sizeof(double));
    
    // Allocate and save grid

    lv.x = new (nothrow) double [nx]; 
    lv.y = new (nothrow) double [ny]; 
    lv.z = new (nothrow) double [nz]; 

    for(int i=0;i<nx;i++)  lv.x[i] = xmin + i*dx;
    for(int j=0;j<ny;j++)  lv.y[j] = ymin + j*dy;
    for(int k=0;k<nz;k++)  lv.z[k] = zmin + k*dz;
    

    ntot = nx*ny*nz;

    // Allocate grid variables
    lv.lapse = new (nothrow) double [ntot]; 
    lv.betax = new (nothrow) double [ntot]; 
    lv.betay = new (nothrow) double [ntot]; 
    lv.betaz = new (nothrow) double [ntot]; 
    lv.phi = new (nothrow) double [ntot]; 
    lv.A11 = new (nothrow) double [ntot]; 
    lv.A12 = new (nothrow) double [ntot]; 
    lv.A13 = new (nothrow) double [ntot]; 
    lv.A22 = new (nothrow) double [ntot]; 
    lv.A23 = new (nothrow) double [ntot]; 
    lv.A33 = new (nothrow) double [ntot]; 
    lv.rhoo = new (nothrow) double [ntot]; 
    lv.ux = new (nothrow) double [ntot]; 
    lv.uy = new (nothrow) double [ntot]; 
    lv.uz = new (nothrow) double [ntot]; 

    for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){
	for(int k=0;k<nz;k++){

	  int ind = i + j*nx + k*nx*ny;

	  infile1.read((char *) &lv.rhoo[ind],sizeof(double));
	  infile1.read((char *) &lv.ux[ind],sizeof(double));
	  infile1.read((char *) &lv.uy[ind],sizeof(double));
	  infile1.read((char *) &lv.uz[ind],sizeof(double));
	  infile1.read((char *) &lv.phi[ind],sizeof(double));
	  infile1.read((char *) &lv.lapse[ind],sizeof(double));
	  infile1.read((char *) &lv.betax[ind],sizeof(double));
	  infile1.read((char *) &lv.betay[ind],sizeof(double));
	  infile1.read((char *) &lv.betaz[ind],sizeof(double));
	  infile1.read((char *) &lv.A11[ind],sizeof(double));
	  infile1.read((char *) &lv.A12[ind],sizeof(double));
	  infile1.read((char *) &lv.A13[ind],sizeof(double));	  
	  infile1.read((char *) &lv.A22[ind],sizeof(double));
	  infile1.read((char *) &lv.A23[ind],sizeof(double));
	  infile1.read((char *) &lv.A33[ind],sizeof(double));       

	}
      }
    }

    infile1.read((char *) &Omega, sizeof(double));
    infile1.read((char *) &xrot, sizeof(double));
    infile1.read((char *) &MNS_b, sizeof(double));
    infile1.read((char *) &MWD_b, sizeof(double));
    infile1.read((char *) &M_ADM_surf, sizeof(double));
    infile1.read((char *) &M_K_surf, sizeof(double));
    infile1.read((char *) &J_ADM_surf, sizeof(double));
    infile1.read((char *) &doubledum, sizeof(double));
    infile1.read((char *) &doubledum, sizeof(double));
    infile1.read((char *) &doubledum, sizeof(double));
    infile1.read((char *) &an1, sizeof(double));
    infile1.read((char *) &an2, sizeof(double));
    infile1.read((char *) &an, sizeof(double));
    infile1.read((char *) &kappa1, sizeof(double));
    infile1.read((char *) &kappa2, sizeof(double));
    infile1.read((char *) &kappa, sizeof(double));
    infile1.read((char *) &rhoo1, sizeof(double));
    infile1.read((char *) &rhoo2, sizeof(double));
    infile1.close();  

    if (print_params == 1){
      printf("\n");
      cout<<"Grid points" <<endl;
      cout<<"nx="<<nx<<endl;
      cout<<"ny="<<nx<<endl;
      cout<<"nz="<<nx<<endl;
      printf("\n");
      cout<<"Kinematic parameters" <<endl;
      cout<<"x coord of rotation axis="<<xrot<<endl;
      cout<<"Omega="<<Omega<<endl;
      printf("\n");
      cout<<"Energy and Angular momentum" <<endl;
      cout<<"MNS_baryon="<< MNS_b<<endl;
      cout<<"MWD_baryon="<< MWD_b<<endl;
      cout<<"ADM mass="<<M_ADM_surf<<endl;
      cout<<"Komar mass="<< M_K_surf <<endl;
      cout<<"Angular momentum ="<< J_ADM_surf <<endl;
      printf("\n");
      cout<<"EOS parameters" <<endl;
      cout<<"an1 ="<< an1 <<endl;
      cout<<"an2 ="<< an2 <<endl;
      cout<<"an3 ="<< an <<endl;
      cout<<"Gamma1 ="<< 1. + 1./an1 <<endl;
      cout<<"Gamma2 ="<< 1. + 1./an2 <<endl;
      cout<<"Gamma3 ="<< 1. + 1./an <<endl;
      cout<<"kappa1 ="<< kappa1 <<endl;
      cout<<"kappa2 ="<< kappa2 <<endl;
      cout<<"kappa3 ="<< kappa <<endl;
      cout<<"rhoo1 ="<< rhoo1 <<endl;
      cout<<"rhoo2 ="<< rhoo2 <<endl;
    }
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

