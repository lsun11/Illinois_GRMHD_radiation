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

extern "C" void CCTK_FCALL read_inputfile_bbhlorene_
  (const cGH **cctkGH,int *genID_cmdline_output_enable,int *fisheye_enable,
   double *f,double *xmin,double *ymin,double *zmin,
   double *dx,double *dy,double *dz,int *ext,
   double &mbhin, double &pbhin, double &sep0in, double &sbhin);

extern "C" void read_inputfile_bbhlorene(const cGH *cctkGH,int genID_cmdline_output_enable,int fisheye_enable,
					 double *f,double xmin,double ymin,double zmin,
					 double dx,double dy,double dz,int *ext,
					 double &mbhin, double &pbhin, double &sep0in, double &sbhin) {

  char filename[100];
  sprintf(filename,"phi_bbh-proc%d.dat",ext[1]*3015+ext[2]*9121+ext[0]*CCTK_MyProc(cctkGH)+(int)(1198311*dx)-(int)(ymin));

  if(genID_cmdline_output_enable==1) {
    printf("./bbhnew %e %e %e %e %e %e %d %d %d %d\n",xmin,ymin,zmin,dx,dy,dz,ext[0],ext[1],ext[2],ext[1]*3015+ext[2]*9121+ext[0]*CCTK_MyProc(cctkGH)+(int)(1198311*dx)-(int)(ymin));
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

    infile1.read((char *) &mbhin, sizeof(double));
    // Make sure to enter the following two values as P/M_tot and r/M_tot!!
    infile1.read((char *) &pbhin, sizeof(double));
    infile1.read((char *) &sep0in, sizeof(double));
    // Make sure to enter the following value as S/M_tot^2!!
    infile1.read((char *) &sbhin, sizeof(double));
    cout << "Mbh/M_tot="<<mbhin <<" P/M_tot="<<pbhin<<" r/M_tot="<<sep0in<<" S/M_tot^2="<<sbhin<< endl;

    if(fisheye_enable==1) {
      int nFEbin = 0;
      double *aFEbin,*r0FEbin,*sFEbin;
      infile1.read((char *) &nFEbin, sizeof(int));
      aFEbin = new double[nFEbin+1];
      r0FEbin = new double[nFEbin+1];
      sFEbin = new double[nFEbin+1];
      infile1.read((char *) &aFEbin[0], sizeof(double));
      for(int i=1; i<=nFEbin; i++)infile1.read((char *) &aFEbin[i], sizeof(double));
      for(int i=1; i<=nFEbin; i++)infile1.read((char *) &r0FEbin[i], sizeof(double));
      for(int i=1; i<=nFEbin; i++)infile1.read((char *) &sFEbin[i], sizeof(double));
  

      int nFEasc = 0;
      double *aFEasc,*r0FEasc,*sFEasc;
      ifstream file;
      char buf[100],c;
 
      file.open("RadialCoordinate_Input");
      if (!file) {
	cerr << "Can't open RadialCoordinate_Input" << endl;
	exit(1);
      }
      file.get(buf,100,'='); file.get(c); file >> nFEasc;
      if(nFEasc != nFEbin) {
	printf("Either I've hit the end of the file prematurely, or \n");
	printf("there is an error with n in RadialCoordinate_Input: File value n = %d, Local RadialCoordinate_Input value n = %d\n",nFEasc,nFEbin);
	exit(1);
      }
      aFEasc = new double[nFEasc+1];
      r0FEasc = new double[nFEasc+1];
      sFEasc = new double[nFEasc+1];
  
      aFEasc[0]=0;r0FEasc[0]=0;sFEasc[0]=0;
  
      file.get(buf,100,'='); file.get(c); file >> aFEasc[0];
      if(aFEasc[0] != aFEbin[0]) {
	printf("Error with aFE[0] in Radial: %.15e %.15e\n",aFEasc[0],aFEbin[0]);
	exit(1);
      }
      if(nFEasc>0) {
	for(int i=1; i<=nFEasc; i++){
	  file >> aFEasc[i];
	  if(aFEasc[i] != aFEbin[i]) {
	    printf("Error with aFE[i] in Radial: %.15e %.15e %.15e\n",i,aFEasc[i],aFEbin[i]);
	    exit(1);
	  }
	}
	file.get(buf,100,'='); file.get(c); file >> r0FEasc[1];
	if(fabs(r0FEasc[1] - r0FEbin[1]) > r0FEbin[1]*1e-12) {
	  printf("Error with r0FE[1] in Radial: %.15e %.15e\n",r0FEasc[1],r0FEbin[1]);
	  exit(1);
	}
	for(int i=2; i<=nFEasc; i++){
	  file >> r0FEasc[i];
	  if(fabs(r0FEasc[i] - r0FEbin[i]) > r0FEbin[i]*1e-12) {
	    printf("Error with r0FE[i] in Radial: %.15e %.15e %.15e\n",i,r0FEasc[i],r0FEbin[i]);
	    exit(1);
	  }
	}
	file.get(buf,100,'='); file.get(c); file >> sFEasc[1];
	if(fabs(sFEasc[1] - sFEbin[1]) > sFEbin[1]*1e-12) { 
	  printf("Error with sFE[1] in Radial: %.15e %.15e\n",sFEasc[1],sFEbin[1]);
	  exit(1);
	}
	for(int i=2; i<=nFEasc; i++){
	  file >> sFEasc[i];
	  if(fabs(sFEasc[i] - sFEbin[i]) > sFEbin[i]*1e-12) { 
	    printf("Error with sFE[i] in Radial: %.15e %.15e %.15e\n",i,sFEasc[i],sFEbin[i]);
	    exit(1);
	  }
	}
      }
      file.close();
    }

    infile1.close();  
  }

}

extern "C" void CCTK_FCALL read_inputfile_bbhlorene_
  (const cGH **cctkGH,int *genID_cmdline_output_enable,int *fisheye_enable,
   double *f,double *xmin,double *ymin,double *zmin,
   double *dx,double *dy,double *dz,int *ext,
   double &mbhin, double &pbhin, double &sep0in, double &sbhin)
{  
  read_inputfile_bbhlorene(*cctkGH,*genID_cmdline_output_enable,*fisheye_enable,
			   f,*xmin,*ymin,*zmin,*dx,*dy,*dz,ext,
			   mbhin,pbhin,sep0in,sbhin);
}
