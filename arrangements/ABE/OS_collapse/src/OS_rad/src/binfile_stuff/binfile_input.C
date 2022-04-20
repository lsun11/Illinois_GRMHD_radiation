//-----------------------------------------------------------------------
// $Id$
//-----------------------------------------------------------------------
// Read binary files and do fancy things with them...
//-----------------------------------------------------------------------
#include <iostream>
#include <sstream> 
#include <string>  

#include <cmath>
#include <iomanip> 
#include <fstream> 
#include <stdlib.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

using namespace std;

extern "C" void CCTK_FCALL CCTK_FNAME(read_binfile)
  (char *gridfuncname,int *binfile_checkpoint_iteration,double *f,const cGH **cctkGH,double *time,
   int *symm_iadj_glob,int *symm_jadj_glob,int *symm_kadj_glob,
   int *local_origin,int *loc_ext);

extern "C" void read_binfile(char *gridfuncname,int binfile_checkpoint_iteration, double *f,
			     const cGH *cctkGH,double *time,
			     int symm_iadj_glob,int symm_jadj_glob,int symm_kadj_glob, int *local_origin,int *loc_ext) {

  int glob_imin,glob_imax,glob_jmin;
  int glob_jmax,glob_kmin,glob_kmax;

  glob_imin = local_origin[0];
  glob_jmin = local_origin[1];
  glob_kmin = local_origin[2];
  
  glob_imax = local_origin[0]+loc_ext[0];
  glob_jmax = local_origin[1]+loc_ext[1];
  glob_kmax = local_origin[2]+loc_ext[2];

  //Clever trick to remove the uglies inside a FORTRAN-created string:
  int index=0;
  while(gridfuncname[index]!='\0') {
    if(gridfuncname[index]==' ') gridfuncname[index]='\0';
    index++;
  }

  char filename[30];

  sprintf(filename,"%s.it",gridfuncname);
  //  printf("hi.... %s\n",filename);

  if     (binfile_checkpoint_iteration<10) sprintf(filename,"%s00000%d.bin",filename,binfile_checkpoint_iteration);
  else if(binfile_checkpoint_iteration<100) sprintf(filename,"%s0000%d.bin",filename,binfile_checkpoint_iteration);
  else if(binfile_checkpoint_iteration<1000) sprintf(filename,"%s000%d.bin",filename,binfile_checkpoint_iteration);
  else if(binfile_checkpoint_iteration<10000) sprintf(filename,"%s00%d.bin",filename,binfile_checkpoint_iteration);
  else if(binfile_checkpoint_iteration<100000) sprintf(filename,"%s0%d.bin",filename,binfile_checkpoint_iteration);
  else if(binfile_checkpoint_iteration<1000000) sprintf(filename,"%s%d.bin",filename,binfile_checkpoint_iteration);

  printf("Reading %s now...\n",filename);

  ifstream infile1;
  infile1.open(filename);
  if(!infile1) {
    cerr << "\a Can't open " << filename << " for input." << endl;
    exit(1);
  }
  double data;
  // read properties of the binary file
  int nx, ny, nz;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  infile1.seekg(0, ios::beg);
  infile1.read((char *) time, sizeof(double));
  infile1.read((char *) &nx, sizeof(int));
  infile1.read((char *) &ny, sizeof(int));
  infile1.read((char *) &nz, sizeof(int));
  infile1.read((char *) &xmin, sizeof(double));
  infile1.read((char *) &xmax, sizeof(double));
  infile1.read((char *) &ymin, sizeof(double));
  infile1.read((char *) &ymax, sizeof(double));
  infile1.read((char *) &zmin, sizeof(double));
  infile1.read((char *) &zmax, sizeof(double));
  int i = 0, j = 0, k = 0; 

  char g[100];

  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {
	infile1.read((char *) &data, sizeof(double));
	if(i>=glob_imin && i<=glob_imax  &&
	   j>=glob_jmin && j<=glob_jmax  &&
	   k>=glob_kmin && k<=glob_kmax) {
	  f[CCTK_GFINDEX3D(cctkGH,i+symm_iadj_glob-glob_imin,j+symm_jadj_glob-glob_jmin,k+symm_kadj_glob-glob_kmin)] = data;
	}
      }
    }
  }
  infile1.close();
}
extern "C" void CCTK_FCALL CCTK_FNAME(read_binfile)
  (char *gridfuncname,int *binfile_checkpoint_iteration,double *f,const cGH **cctkGH,double *time,
   int *symm_iadj_glob,int *symm_jadj_glob,int *symm_kadj_glob, int *local_origin,int *loc_ext) {
  
  read_binfile(gridfuncname,*binfile_checkpoint_iteration,f,*cctkGH,time,*symm_iadj_glob,*symm_jadj_glob,*symm_kadj_glob,local_origin,loc_ext);

}
