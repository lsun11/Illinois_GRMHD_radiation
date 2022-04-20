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

extern "C" void CCTK_FCALL CCTK_FNAME(read_inputfile_NPunctures)
  (const cGH **cctkGH,int &genID_cmdline_output_enable,int &fisheye_enable,
   double *f,double &xmin,double &ymin,double &zmin,
   double &dx,double &dy,double &dz,int *ext,
   double & half_binary_separation, double & x_offset, int &npunctures_numbhs,
   double *bhmass,
   double *bhxpos, double *bhypos,
   double *bh_px, double *bh_py, double *bh_pz,
   double *bh_sx, double *bh_sy, double *bh_sz,
   double & BigM);

extern "C" void read_inputfile_NPunctures(const cGH *cctkGH,int &genID_cmdline_output_enable,int &fisheye_enable,
					  double *f,double &xmin,double &ymin,double &zmin,
					  double &dx,double &dy,double &dz,int *ext,
					  double & half_binary_separation, double & x_offset, int &npunctures_numbhs,
					  double *bhmass,
					  double *bhxpos, double *bhypos,
					  double *bh_px, double *bh_py, double *bh_pz,
					  double *bh_sx, double *bh_sy, double *bh_sz,
					  double & BigM)
{
  char filename[100];
  sprintf(filename,"psi_bbh-proc%d.dat",ext[1]*3015+ext[2]*91211+ext[0] + 1241*CCTK_MyProc(cctkGH)+(int)(1198311*dx)+(int)(xmin*1198)+(int)(ymin*2919));

  if(genID_cmdline_output_enable==1) {
    printf("./bbh_NPunctures %e %e %e %e %e %e %d %d %d %d\n",xmin,ymin,zmin,dx,dy,dz,ext[0],ext[1],ext[2],ext[1]*3015+ext[2]*91211+ext[0] + 1241*CCTK_MyProc(cctkGH)+(int)(1198311*dx)+(int)(xmin*1198)+(int)(ymin*2919));
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
      infile1.read((char *) &f[CCTK_GFINDEX3D(cctkGH,i,j,k)], sizeof(double));
    }

    printf("=========================================================================================\n");
    printf("Okay, the initial data has been read in...  Now I'm going to check additional parameters!\n");
    printf("=========================================================================================\n");

    double b, mplus,pplusx,pplusy,pplusz,splusx,splusy,splusz;
    double mminus,pminusx,pminusy,pminusz,sminusx,sminusy,sminusz;

    infile1.read((char *) &half_binary_separation, sizeof(double));
    infile1.read((char *) &x_offset, sizeof(double));

    infile1.read((char *) &npunctures_numbhs, sizeof(int));
    
    for(int i=0;i<npunctures_numbhs;i++) {
      infile1.read((char *) &bhxpos[i], sizeof(double));

      infile1.read((char *) &bhxpos[i], sizeof(double));
      infile1.read((char *) &bhypos[i], sizeof(double));

      infile1.read((char *) &bh_px[i], sizeof(double));
      infile1.read((char *) &bh_py[i], sizeof(double));
      infile1.read((char *) &bh_pz[i], sizeof(double));
      
      infile1.read((char *) &bh_sx[i], sizeof(double));
      infile1.read((char *) &bh_sy[i], sizeof(double));
      infile1.read((char *) &bh_sz[i], sizeof(double));
      
      if(bh_sz[i]!=0) {
	printf("ERROR. BH number %d has bh_sz!=0!\n",bh_sz[i]);
	exit(1);
      }
    }

    infile1.read((char *) &BigM, sizeof(double));
    
    printf("ADM mass of the binary = %f\n",BigM);

    if(fisheye_enable==1) {
	printf("Fisheye not supported in NPunctures!");
	exit(1);
    }
  
    infile1.close();  
  }

}

extern "C" void CCTK_FCALL CCTK_FNAME(read_inputfile_NPunctures)
  (const cGH **cctkGH,int &genID_cmdline_output_enable,int &fisheye_enable,
   double *f,double &xmin,double &ymin,double &zmin,
   double &dx,double &dy,double &dz,int *ext,
   double & half_binary_separation, double & x_offset, int &npunctures_numbhs, 
   double *bhmass,
   double *bhxpos, double *bhypos,
   double *bh_px, double *bh_py, double *bh_pz,
   double *bh_sx, double *bh_sy, double *bh_sz,
   double & BigM)
{  
  read_inputfile_NPunctures(*cctkGH,genID_cmdline_output_enable,fisheye_enable,
			    f,xmin,ymin,zmin,
			    dx,dy,dz,ext,
			    half_binary_separation,  x_offset, npunctures_numbhs,
			    bhmass,
			    bhxpos, bhypos,
			    bh_px, bh_py, bh_pz,
			    bh_sx, bh_sy, bh_sz,
			    BigM);
}

