//-----------------------------------------------------------------------
// Read BHNS bin. files and do fancy things with them...
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

using namespace std;

extern "C" void CCTK_FCALL CCTK_FNAME(read_inputfile_bhns)
  (const cGH **cctkGH,int *genID_cmdline_output_enable,int *fisheye_enable,
   int *reset_shift_lapse,double *xmin,double *ymin,double *zmin,
   double *dx,double *dy,double *dz,int *ext,
   double *lapse,double *shiftx,double *shifty,double *shiftz,
   double *psi,double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,
   double *rho_nbar,double *ux,double *uy, double *uz,
   double &gamma_th,double &xbh,double &mirr);

extern "C" void read_inputfile_bhns(const cGH *cctkGH,int genID_cmdline_output_enable,int fisheye_enable,
				    int reset_shift_lapse,double xmin,double ymin,double zmin,
				    double dx,double dy,double dz,int *ext,
				    double *lapse,double *shiftx,double *shifty,double *shiftz,
				    double *psi,double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,
				    double *rho_nbar,double *ux,double *uy, double *uz,
				    double &gamma_th,double &xbh,double &mirr) {

  int filename_checksum = ext[1]*301+ext[2]*1+ext[0]*1000000 + 1241*CCTK_MyProc(cctkGH)+(int)(999999*dx)+(int)(fabs(xmin)*10.)+(int)(fabs(zmin)*100.);
  char filename[100];
  sprintf(filename,"ini_bhns-proc%d.d",filename_checksum);

  if(genID_cmdline_output_enable==1) {
    printf("./bhns_initialdata_isotropic %.16e %.16e %.16e %.16e %.16e %.16e %d %d %d %d\n",xmin,ymin,zmin,dx,dy,dz,ext[0],ext[1],ext[2],filename_checksum);
    printf("blah %d %d %d %d %d %d %d\n",ext[1]*301,ext[2]*1,ext[0]*1000000,1241*CCTK_MyProc(cctkGH),(int)(999999*dx),(int)(fabs(xmin)*10.),(int)(fabs(zmin)*100.));
  } else {

    printf("Attempting to read %s now...  \nNote that you'll need to store the initial data files so that ALL processors can see them!\n",filename);

    ifstream infile1;
    infile1.open(filename);
    if(!infile1) {
      cerr << "\a Can't open " << filename << " for input." << endl;
      exit(1);
    }

    double doubledum;
    int ntot,intdum;
  
    // Read BHNS parameters
    infile1.read((char *) &gamma_th, sizeof(double));
    infile1.read((char *) &doubledum, sizeof(double));
    infile1.read((char *) &doubledum, sizeof(double));
    infile1.read((char *) &doubledum, sizeof(double));
    infile1.read((char *) &doubledum, sizeof(double));
    cout<<"Omega="<<doubledum<<endl;
    infile1.read((char *) &doubledum, sizeof(double));
    cout<<"ADM mass="<<doubledum<<endl;
    infile1.read((char *) &doubledum, sizeof(double));
    cout<<"Ang mom="<<doubledum<<endl;
    infile1.read((char *) &mirr, sizeof(double));
    infile1.read((char *) &doubledum, sizeof(double));
    cout<<"r_AH="<<doubledum<<endl;
    infile1.read((char *) &doubledum, sizeof(double));
    cout<<"M_baryon="<<doubledum<<endl;
    infile1.read((char *) &doubledum, sizeof(double));
    infile1.read((char *) &xbh, sizeof(double));
    infile1.read((char *) &doubledum, sizeof(double));
    infile1.read((char *) &doubledum, sizeof(double));
    infile1.read((char *) &doubledum, sizeof(double));
    infile1.read((char *) &doubledum, sizeof(double));
    infile1.read((char *) &doubledum, sizeof(double));
    infile1.read((char *) &doubledum, sizeof(double));
    infile1.read((char *) &doubledum, sizeof(double));
    infile1.read((char *) &doubledum, sizeof(double));
    infile1.read((char *) &intdum, sizeof(int));
    infile1.read((char *) &ntot, sizeof(int));

    int nx = ext[0];
    int ny = ext[1];
    int nz = ext[2];

    if(ntot!=nx*ny*nz) {
      printf("MAYBE PROBLEM: This data file has %d elements, but current grid has %d!\n",ntot,nx*ny*nz);
      //exit(1);
    } else {
      printf("This data file has %d elements, like we expect!\n",ntot,nx*ny*nz);
    }

    //    if(reset_shift_lapse==0) {
      for(int k=0;k<nz;k++) for(int j=0;j<ny;j++) for(int i=0;i<nx;i++) infile1.read((char *) &lapse[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
      for(int k=0;k<nz;k++) for(int j=0;j<ny;j++) for(int i=0;i<nx;i++) infile1.read((char *) &shiftx[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
      for(int k=0;k<nz;k++) for(int j=0;j<ny;j++) for(int i=0;i<nx;i++) infile1.read((char *) &shifty[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
      for(int k=0;k<nz;k++) for(int j=0;j<ny;j++) for(int i=0;i<nx;i++) infile1.read((char *) &shiftz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
      //} else {
      for(int k=0;k<nz;k++) for(int j=0;j<ny;j++) for(int i=0;i<nx;i++) infile1.read((char *) &psi[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
      for(int k=0;k<nz;k++) for(int j=0;j<ny;j++) for(int i=0;i<nx;i++) infile1.read((char *) &kxx[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
      for(int k=0;k<nz;k++) for(int j=0;j<ny;j++) for(int i=0;i<nx;i++) infile1.read((char *) &kxy[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
      for(int k=0;k<nz;k++) for(int j=0;j<ny;j++) for(int i=0;i<nx;i++) infile1.read((char *) &kxz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
      for(int k=0;k<nz;k++) for(int j=0;j<ny;j++) for(int i=0;i<nx;i++) infile1.read((char *) &kyy[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
      for(int k=0;k<nz;k++) for(int j=0;j<ny;j++) for(int i=0;i<nx;i++) infile1.read((char *) &kyz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));

      //NEXT GF IS A DUMMY.  We know Kzz from the other Kij's and the fact that trK==0!
      for(int k=0;k<nz;k++) for(int j=0;j<ny;j++) for(int i=0;i<nx;i++) infile1.read((char *) &rho_nbar[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));

      for(int k=0;k<nz;k++) for(int j=0;j<ny;j++) for(int i=0;i<nx;i++) infile1.read((char *) &rho_nbar[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));

      for(int k=0;k<nz;k++) for(int j=0;j<ny;j++) for(int i=0;i<nx;i++) infile1.read((char *) &ux[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
      for(int k=0;k<nz;k++) for(int j=0;j<ny;j++) for(int i=0;i<nx;i++) infile1.read((char *) &uy[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
      for(int k=0;k<nz;k++) for(int j=0;j<ny;j++) for(int i=0;i<nx;i++) infile1.read((char *) &uz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
      //    }

    infile1.close();  
  }

}

extern "C" void CCTK_FCALL CCTK_FNAME(read_inputfile_bhns)
  (const cGH **cctkGH,int *genID_cmdline_output_enable,int *fisheye_enable,
   int *reset_shift_lapse,double *xmin,double *ymin,double *zmin,
   double *dx,double *dy,double *dz,int *ext,
   double *lapse,double *shiftx,double *shifty,double *shiftz,
   double *psi,double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,
   double *rho_nbar,double *ux,double *uy, double *uz,
   double &gamma_th,double &xbh,double &mirr)
{  
  read_inputfile_bhns(*cctkGH,*genID_cmdline_output_enable,*fisheye_enable,
		      *reset_shift_lapse,*xmin,*ymin,*zmin,*dx,*dy,*dz,ext,
		      lapse,shiftx,shifty,shiftz,
		      psi,kxx,kxy,kxz,kyy,kyz,
		      rho_nbar,ux,uy, uz,
		      gamma_th,xbh,mirr);
}

