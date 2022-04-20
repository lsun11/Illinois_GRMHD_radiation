/*
  DAGH-based binfile bns_output thorn
*/

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
CCTK_FILEVERSION(bns_outputbinfile)

extern "C" CCTK_INT bns_outputbinfile(const cGH *cctkGH,int *glob_ext,int *loc_ext,int *loc_origin,double cctk_time,double *xyz_glob_min,double *xyz_glob_max,
				  int cctk_iteration,double *f,char *gridfuncname,int Symmetry);

extern "C" void CCTK_FCALL CCTK_FNAME(bns_outputbinfile)
  (CCTK_INT *retval,const cGH **cctkGH, int *glob_ext,int *loc_ext,int *loc_origin,double *cctk_time,double *xyz_glob_min,double *xyz_glob_max,
   int *cctk_iteration,double *f,char *gridfuncname,int *Symmetry);

extern "C" CCTK_INT bns_outputbinfile(const cGH *cctkGH,int *glob_ext,int *loc_ext,int *loc_origin,double cctk_time,double *xyz_glob_min,double *xyz_glob_max,
				  int cctk_iteration,double *f,char *gridfuncname,int Symmetry)
{

  int upper[3],lower[3];
  double coordulbound[2];
  //strstream filename;
  std::ostringstream filename;

  //Clever trick to remove the uglies inside a FORTRAN-created string:
  if(CCTK_MyProc(cctkGH)==0) {
    int i=0;
    while(gridfuncname[i]!='\0') {
      if(gridfuncname[i]==' ') gridfuncname[i]='\0';
      i++;
    }
  }

  filename << gridfuncname << ".it" << setfill('0') << setw(6) 
  	   << cctk_iteration << ".bin" << ends;

  double xmin,ymin,zmin;
  double xmax,ymax,zmax;
  int handle,index;

  int NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4;

  int ext[3];

  int sizex_glob,sizey_glob,sizez_glob;
  int sizex_loc,sizey_loc,sizez_loc;
  int originx_loc,originy_loc,originz_loc;

  sizex_glob = glob_ext[0];
  sizey_glob = glob_ext[1];
  sizez_glob = glob_ext[2];

  sizex_loc = loc_ext[0];
  sizey_loc = loc_ext[1];
  sizez_loc = loc_ext[2];

  originx_loc = loc_origin[0];
  originy_loc = loc_origin[1];
  originz_loc = loc_origin[2];

  if(Symmetry == AXISYM) {
    ext[0] = sizex_glob;
    ext[1] = sizey_glob;
    ext[2] = sizez_glob;
    //Following accounts for ghostzones.  If following lines excluded, the CCTK_Reduce function call would double count across proc boundaries.
    if(originx_loc+sizex_loc < sizex_glob) sizex_loc-=4;
    if(originz_loc+sizez_loc < sizez_glob) sizez_loc-=2;

  } 
  else if(Symmetry == EQUATORIAL) { 
    ext[0] = sizex_glob;
    ext[1] = sizey_glob;
    ext[2] = sizez_glob-1;
    //Following accounts for ghostzones.  If following lines excluded, the CCTK_Reduce function call would double count across proc boundaries.
    if(originx_loc+sizex_loc < sizex_glob) sizex_loc-=2;
    if(originy_loc+sizey_loc < sizey_glob) sizey_loc-=2;
    if(originz_loc+sizez_loc < sizez_glob) sizez_loc-=2;

  }

  upper[0] = ext[0];
  upper[1] = ext[1];
  upper[2] = ext[2];
  lower[0] = 0;
  lower[1] = 0;
  lower[2] = 0;

  if(glob_ext[1] == 3) {
    upper[1]++;
  }

  //  xmin = x[CCTK_GFINDEX3D(cctkGH,0,0,0)];
  xmin = xyz_glob_min[0];
  ymin = xyz_glob_min[1];
  zmin = xyz_glob_min[2];

  xmax = xyz_glob_max[0];
  ymax = xyz_glob_max[1];
  zmax = xyz_glob_max[2];

  //First combine data from multiple processors into one big data array
  double *onebigdataarray;
  //onebigdataarray = malloc(glob_ext[0]*glob_ext[1]*glob_ext[2]*sizeof(double));
  onebigdataarray = new double[glob_ext[0]*glob_ext[1]*glob_ext[2]];

  //Initialize to zero
  for(int i=0;i<glob_ext[0]*glob_ext[1]*glob_ext[2];i++) onebigdataarray[i]=0.0;
  
  for(int i=0;i<sizex_loc;i++) 
    for(int j=0;j<sizey_loc;j++) 
      for(int k=0;k<sizez_loc;k++) {
	int iglob = i+originx_loc;
	int jglob = j+originy_loc;
	int kglob = k+originz_loc;
	
	int vindex_loc = CCTK_GFINDEX3D(cctkGH,i,j,k);
	long vindex_glob = iglob + jglob*sizex_glob+ kglob*sizex_glob*sizey_glob;
	onebigdataarray[vindex_glob] = f[vindex_loc];
      }
  /*
  double *onebigdataarray_o;
  onebigdataarray_o = malloc(glob_ext[0]*glob_ext[1]*glob_ext[2]*sizeof(double));
  //Initialize to zero
  for(int i=0;i<glob_ext[0]*glob_ext[1]*glob_ext[2];i++) onebigdataarray_o[i]=0.0;
  */

  CCTK_ReduceLocArrayToArray1D(cctkGH,-1,CCTK_ReductionHandle("sum"),onebigdataarray,onebigdataarray,glob_ext[0]*glob_ext[1]*glob_ext[2],CCTK_VARIABLE_REAL);

  //---------------------------------
  // Write to file ON PROCESSOR ZERO  >
  //---------------------------------

  if(CCTK_MyProc(cctkGH) == 0) {

    ofstream file;
    printf("outputting binfile %s at time %e\n",filename.str().c_str(),cctk_time);

    file.open(filename.str().c_str(),ios::out | ios::ate);
    //------------\
    // Write time  >
    //------------/
    file.write((char *) &cctk_time, sizeof(double));
    //----------------------\
    // Write integer extents >
    //----------------------/
    int i;
    for(i=0;i<3;i++) {
      if(lower[i] != upper[i]) {
	int extent = upper[i];
	file.write((char *) &extent,sizeof(int));
      }
    }
    //-----------------------------
    // Write world coord boundaries >
    //-----------------------------/
    if(Symmetry == AXISYM) { // || Symmetry == EQUATORIAL) {
      coordulbound[0] = fabs(xmin); 
    } 
    else coordulbound[0] = xmin; 
    coordulbound[1] = xmax; 
    file.write(((char *) (coordulbound)),2*sizeof(double)); 

    coordulbound[1] = ymax; 
    if(Symmetry == AXISYM || Symmetry == EQUATORIAL) {
      coordulbound[0] = -ymax;
    } 
    else coordulbound[0] = ymin; 
    file.write(((char *) (coordulbound)),2*sizeof(double)); 

    if(Symmetry == AXISYM || Symmetry==EQUATORIAL) {
      coordulbound[0] = fabs(zmin); 
    } 
    else coordulbound[0] = zmin;
    coordulbound[1] = zmax; 
    file.write(((char *) (coordulbound)),2*sizeof(double));

    //--------------------------
    // Write data and close file >
    //--------------------------/
    
    int imin=0,jmin=0,kmin=0;
    if(Symmetry == AXISYM) {
      imin = imin + 1;
      kmin = kmin + 1;
    }
    else if(Symmetry == OCTANT) {
      imin = imin + 1;
      jmin = jmin + 1;
      kmin = kmin + 1;
    }
    else if(Symmetry == EQUATORIAL) {
      kmin = kmin + 1;
    }
    for(int k = kmin;k<ext[2];k++) {
      for(int j = jmin;j<ext[1];j++) {
	for(int i = imin;i<ext[0];i++) {
	  //int vindex = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  int vindex = i + j*sizex_glob+ k*sizex_glob*sizey_glob;
	  //	  printf("Hello. %d %d %d %e\n",i,j,k,onebigdataarray[vindex]);
	  file.write((char *) &onebigdataarray[vindex],sizeof(double));
	}
      }
    }
    file.close();
  }

  delete [] onebigdataarray;
  return(0);
}


extern "C" void CCTK_FCALL CCTK_FNAME(bns_outputbinfile) 
  (CCTK_INT *retval,const cGH **cctkGH,int *glob_ext,int *loc_ext,int *loc_origin,double *cctk_time,double *xyz_glob_min,double *xyz_glob_max,
   int *cctk_iteration,double *f,char *gridfuncname,int *Symmetry)
{
  *retval = bns_outputbinfile(*cctkGH, glob_ext, loc_ext, loc_origin,*cctk_time,xyz_glob_min,xyz_glob_max,
			  *cctk_iteration, f, gridfuncname, *Symmetry);
}
