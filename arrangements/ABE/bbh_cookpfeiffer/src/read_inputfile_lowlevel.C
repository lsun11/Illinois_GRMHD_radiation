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

using namespace std;

extern "C" void CCTK_FCALL CCTK_FNAME(read_inputfile_bbhcp)
  (const cGH **cctkGH,int *genID_cmdline_output_enable,int *fisheye_enable,
   double *phi,double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,double *shiftx,double *shifty,double *shiftz,double *lapse,
   double *xmin,double *ymin,double *zmin,
   double *dx,double *dy,double *dz,int *ext);

extern "C" void read_inputfile_bbhcp(const cGH *cctkGH,int genID_cmdline_output_enable,int fisheye_enable,
				     double *phi,double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,double *shiftx,double *shifty,double *shiftz,double *lapse,
				     double xmin,double ymin,double zmin,
				     double dx,double dy,double dz,int *ext) {

  char filename[100];
  //sprintf(filename,"outfields_bbhcp_proc-%d.dat",ext[1]+ext[2]+ext[0]*CCTK_MyProc(cctkGH)+(int)(119831*dx));
  sprintf(filename,"outfields_bbhcp_proc-%d.dat",ext[1]*3015+ext[2]*9121+ext[0]*CCTK_MyProc(cctkGH)+(int)(1198311*dx)-(int)(ymin));

  if(genID_cmdline_output_enable==1) {
    printf("./getbbhID-fish-parallel %e %e %e %e %e %e %d %d %d %d\n",xmin,ymin,zmin,dx,dy,dz,ext[0],ext[1],ext[2],ext[1]*3015+ext[2]*9121+ext[0]*CCTK_MyProc(cctkGH)+(int)(1198311*dx)-(int)(ymin));
  } else {
    /*
      printf("Now generating %s datafile for local processor...\n",filename);
      char cmdline[100];
      sprintf(cmdline,"./bbhnew %.15e %.15e %.15e %.15e %.15e %.15e %d %d %d %d",xmin,ymin,zmin,dx,dy,dz,ext[0],ext[1],ext[2],CCTK_MyProc(cctkGH));

      system(cmdline);

      sleep(100);
      exit(0);
    */

    printf("Attempting to read %s now...  \nNote that you'll need to store the initial data files so that ALL processors can see them!\n",filename);

    ifstream infile1;
    infile1.open(filename);
    if(!infile1) {
      cerr << "\a Can't open " << filename << " for input." << endl;
      exit(1);
    }

    int nx = ext[0];
    int ny = ext[1];
    int nz = ext[2];
  
    for (int k = 0; k < nz; k++) for (int j = 0; j < ny; j++) for (int i = 0; i < nx; i++) {
      infile1.read((char *) &phi[CCTK_GFINDEX3D(cctkGH,i,j,k)], sizeof(double));
      infile1.read((char *) &kxx[CCTK_GFINDEX3D(cctkGH,i,j,k)], sizeof(double));
      infile1.read((char *) &kxy[CCTK_GFINDEX3D(cctkGH,i,j,k)], sizeof(double));
      infile1.read((char *) &kxz[CCTK_GFINDEX3D(cctkGH,i,j,k)], sizeof(double));
      infile1.read((char *) &kyy[CCTK_GFINDEX3D(cctkGH,i,j,k)], sizeof(double));
      infile1.read((char *) &kyz[CCTK_GFINDEX3D(cctkGH,i,j,k)], sizeof(double));
      infile1.read((char *) &shiftx[CCTK_GFINDEX3D(cctkGH,i,j,k)], sizeof(double));
      infile1.read((char *) &shifty[CCTK_GFINDEX3D(cctkGH,i,j,k)], sizeof(double));
      infile1.read((char *) &shiftz[CCTK_GFINDEX3D(cctkGH,i,j,k)], sizeof(double));
      infile1.read((char *) &lapse[CCTK_GFINDEX3D(cctkGH,i,j,k)], sizeof(double));
      lapse[CCTK_GFINDEX3D(cctkGH,i,j,k)] = lapse[CCTK_GFINDEX3D(cctkGH,i,j,k)] - 1.0;
    }

    printf("==========================================\n");
    printf("Okay, the initial data has been read in...\n");
    printf("==========================================\n");

    infile1.close();  
  }

}

extern "C" void CCTK_FCALL CCTK_FNAME(read_inputfile_bbhcp)
  (const cGH **cctkGH,int *genID_cmdline_output_enable,int *fisheye_enable,
   double *phi,double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,double *shiftx,double *shifty,double *shiftz,double *lapse,
   double *xmin,double *ymin,double *zmin,
   double *dx,double *dy,double *dz,int *ext)
{  
  read_inputfile_bbhcp(*cctkGH,*genID_cmdline_output_enable,*fisheye_enable,
			   phi,kxx,kxy,kxz,kyy,kyz,shiftx,shifty,shiftz,lapse,*xmin,*ymin,*zmin,*dx,*dy,*dz,ext);
}

