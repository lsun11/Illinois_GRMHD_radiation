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

extern "C" void CCTK_FCALL CCTK_FNAME(read_inputfile_TwoPunctures)
  (const cGH **cctkGH,int *genID_cmdline_output_enable,int *fisheye_enable,
   double *f,double *xmin,double *ymin,double *zmin,
   double *dx,double *dy,double *dz,int *ext,
   double & bh_mass_plus, double & bh_mass_minus, double & bh_px_plus, 
   double & bh_px_minus, double & bh_py_plus, double & bh_py_minus, 
   double & bh_spin_plus, double & bh_spin_minus, 
   double & half_binary_separation, double & x_offset, double & BigM);

extern "C" void read_inputfile_TwoPunctures(const cGH *cctkGH,int genID_cmdline_output_enable,int fisheye_enable,
					 double *f,double xmin,double ymin,double zmin,
					 double dx,double dy,double dz,int *ext,
					 double & bh_mass_plus, double & bh_mass_minus, double & bh_px_plus,
   					 double & bh_px_minus, double & bh_py_plus, double & bh_py_minus,
   					 double & bh_spin_plus, double & bh_spin_minus,
    				 	 double & half_binary_separation, double & x_offset,
					 double & BigM)
{
  char filename[100];

  sprintf(filename,"psi_bbh-proc%d.dat",ext[1]*3015+ext[2]*9121+ext[0]*CCTK_MyProc(cctkGH)+(int)(1198311*dx)-(int)(xmin*5.1295812));
  //sprintf(filename,"psi_bbh-proc%d.dat",ext[1]*3015+ext[2]*91211+ext[0] + 1241*CCTK_MyProc(cctkGH)+(int)(1198311*dx));

  if(genID_cmdline_output_enable==1) {
    printf("./bbh_TwoPunctures %.17e %.17e %.17e %.17e %.17e %.17e %d %d %d %d\n",xmin,ymin,zmin,dx,dy,dz,ext[0],ext[1],ext[2],ext[1]*3015+ext[2]*9121+ext[0]*CCTK_MyProc(cctkGH)+(int)(1198311*dx)-(int)(xmin*5.1295812));
    //ext[1]*3015+ext[2]*91211+ext[0] + 1241*CCTK_MyProc(cctkGH)+(int)(1198311*dx));
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
/*
    infile1.read((char *) &b, sizeof(double));
    infile1.read((char *) &mplus, sizeof(double));
    infile1.read((char *) &pplusx, sizeof(double));
    infile1.read((char *) &pplusy, sizeof(double));
    infile1.read((char *) &pplusz, sizeof(double));
    infile1.read((char *) &splusx, sizeof(double));
    infile1.read((char *) &splusy, sizeof(double));
    infile1.read((char *) &splusz, sizeof(double));
    infile1.read((char *) &mminus, sizeof(double));
    infile1.read((char *) &pminusx, sizeof(double));
    infile1.read((char *) &pminusy, sizeof(double));
    infile1.read((char *) &pminusz, sizeof(double));
    infile1.read((char *) &sminusx, sizeof(double));
    infile1.read((char *) &sminusy, sizeof(double));
    infile1.read((char *) &sminusz, sizeof(double));
    infile1.read((char *) &x_offset, sizeof(double));
    infile1.read((char *) &BigM, sizeof(double));
*/

    infile1.read((char *) &half_binary_separation, sizeof(double));
    infile1.read((char *) &bh_mass_plus, sizeof(double));
    infile1.read((char *) &bh_px_plus, sizeof(double));
    infile1.read((char *) &bh_py_plus, sizeof(double));
    infile1.read((char *) &pplusz, sizeof(double));
    infile1.read((char *) &splusx, sizeof(double));
    infile1.read((char *) &splusy, sizeof(double));
    infile1.read((char *) &bh_spin_plus, sizeof(double));
    infile1.read((char *) &bh_mass_minus, sizeof(double));
    infile1.read((char *) &bh_px_minus, sizeof(double));
    infile1.read((char *) &bh_py_minus, sizeof(double));
    infile1.read((char *) &pminusz, sizeof(double));
    infile1.read((char *) &sminusx, sizeof(double));
    infile1.read((char *) &sminusy, sizeof(double));
    infile1.read((char *) &bh_spin_minus, sizeof(double));
    infile1.read((char *) &x_offset, sizeof(double));
    infile1.read((char *) &BigM, sizeof(double));
    
    printf("ADM mass of the binary = %f\n",BigM);

    if(fisheye_enable==1) {
	printf("Fisheye not supported in TwoPunctures!");
	exit(1);
    }

/*
    if (fabs(b-half_binary_separation) > b*1.e-5) {
       printf("half_binary_separation from data file = %f\n",b);
       printf("half_binary_separation from par file = %f\n",half_binary_separation);
       printf("These two values don't agree. Fix it you silly human being!");
       exit(1);
    }

    if (fabs(mplus - bh_mass_plus) > bh_mass_plus*1.e-5) {
       printf("bh_mass_plus from data file = %f\n",mplus);
       printf("bh_mass_plus from par file = %f\n",bh_mass_plus);
       printf("These two values don't agree. Fix it you silly human being!");
       exit(1);
    }

    if (fabs(pplusx - bh_px_plus) > fabs(bh_px_plus)*1.e-5) { 
       printf("bh_px_plus from data file = %e\n",pplusx);
       printf("bh_px_plus from par file = %e\n",bh_px_plus);
       printf("These two values don't agree. Fix it you silly human being!");
       exit(1);
    }

    if (fabs(pplusy - bh_py_plus) > fabs(bh_py_plus)*1.e-5) {
       printf("bh_py_plus from data file = %e\n",pplusy);
       printf("bh_py_plus from par file = %e\n",bh_py_plus);
       printf("These two values don't agree. Fix it you silly human being!");
       exit(1);
    }
*/

    if (fabs(pplusz) > 1.e-5*BigM) {
       printf("bh_pz_plus from data file = %e\n",pplusz);
       printf("This is inconsistent with equatorial symmetry. Fix it you silly human being!");
       exit(1);
    }

    if (fabs(splusx) > 1.e-5*BigM*BigM) { 
       printf("bh_sx_plus from data file = %e\n",splusx);
       printf("This is inconsistent with equatorial symmetry. Fix it you silly human being!");
       exit(1);
    }

    if (fabs(splusy) > 1.e-5*BigM*BigM) {
       printf("bh_sy_plus from data file = %e\n",splusy);
       printf("This is inconsistent with equatorial symmetry. Fix it you silly human being!");
       exit(1);
    }

/*
    if (fabs(splusz - bh_spin_plus) > fabs(bh_spin_plus)*1.e-5) { 
       printf("bh_spin_plus from data file = %e\n",splusz);
       printf("bh_spin_plus from par file = %e\n",bh_spin_plus);
       printf("These two values don't agree. Fix it you silly human being!");
       exit(1);
    }

    if (fabs(mminus - bh_mass_minus) > bh_mass_minus*1.e-5) {
       printf("bh_mass_minus from data file = %f\n",mminus);
       printf("bh_mass_minus from par file = %f\n",bh_mass_minus);
       printf("These two values don't agree. Fix it you silly human being!");
       exit(1);
    }

    if (fabs(pminusx - bh_px_minus) > fabs(bh_px_minus)*1.e-5) {
       printf("bh_px_minus from data file = %e\n",pminusx);
       printf("bh_px_minus from par file = %e\n",bh_px_minus);
       printf("These two values don't agree. Fix it you silly human being!");
       exit(1);
    }

    if (fabs(pminusy - bh_py_minus) > fabs(bh_py_minus)*1.e-5) {
       printf("bh_py_minus from data file = %e\n",pminusy);
       printf("bh_py_minus from par file = %e\n",bh_py_minus);
       printf("These two values don't agree. Fix it you silly human being!");
       exit(1);
    }
*/

    if (fabs(pminusz) > 1.e-5*BigM) {
       printf("bh_pz_minus from data file = %e\n",pminusz);
       printf("This is inconsistent with equatorial symmetry. Fix it you silly human being!");
       exit(1);
    }

    if (fabs(sminusx) > 1.e-5*BigM*BigM) {
       printf("bh_sx_minus from data file = %e\n",sminusx);
       printf("This is inconsistent with equatorial symmetry. Fix it you silly human being!");
       exit(1);
    }

    if (fabs(sminusy) > 1.e-5*BigM*BigM) {
       printf("bh_sy_minus from data file = %e\n",sminusy);
       printf("This is inconsistent with equatorial symmetry. Fix it you silly human being!");
       exit(1);
    }

/*
    if (fabs(sminusz - bh_spin_minus) > fabs(bh_spin_minus)*1.e-5) {
       printf("bh_spin_minus from data file = %e\n",sminusz);
       printf("bh_spin_minus from par file = %e\n",bh_spin_minus);
       printf("These two values don't agree. Fix it you silly human being!");
       exit(1);
    }
*/
  
    infile1.close();  
  }

}

extern "C" void CCTK_FCALL CCTK_FNAME(read_inputfile_TwoPunctures)
  (const cGH **cctkGH,int *genID_cmdline_output_enable,int *fisheye_enable,
   double *f,double *xmin,double *ymin,double *zmin,
   double *dx,double *dy,double *dz,int *ext,
   double & bh_mass_plus, double & bh_mass_minus, double & bh_px_plus,
   double & bh_px_minus, double & bh_py_plus, double & bh_py_minus,
   double & bh_spin_plus, double & bh_spin_minus,
   double & half_binary_separation, double & x_offset, double & BigM)
{  
  read_inputfile_TwoPunctures(*cctkGH,*genID_cmdline_output_enable,*fisheye_enable,
			   f,*xmin,*ymin,*zmin,*dx,*dy,*dz,ext,
			   bh_mass_plus, bh_mass_minus, bh_px_plus, bh_px_minus,
			   bh_py_plus, bh_py_minus, bh_spin_plus, bh_spin_minus,
		 	   half_binary_separation, x_offset, BigM);
}

