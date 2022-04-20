//-----------------------------------------------------------------------
// Read BHNS_BHSPIN bin. files and do fancy things with them...
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

extern "C" void CCTK_FCALL get_ah_radius_in_dirn_
  (const cGH **cctkGH,const double &horizdirn_x,const double &horizdirn_y,const double &horizdirn_z,double &output_radius);

extern "C" void get_ah_radius_in_dirn(const cGH *cctkGH,const double &horizdirn_x,const double &horizdirn_y,const double &horizdirn_z,double &output_radius) {

  //###########################################################################
  //# HorizonRadiusInDirection() computes the horizon radius in the direction
  //# of each (x,y,z) point, or -1.0 if this horizon wasn't found the most
  //# recent time AHFinderDirect searched for it.  More precisely, for each
  //# (x,y,z), consider the ray from the local coordinate origin through
  //# (x,y,z).  This function computes the Euclidean distance between the
  //# local coordinate origin and this ray's intersection with the horizon,
  //# or -1.0 if this horizon wasn't found the most recent time AHFinderDirect
  //# searched for it.  
  //#

  int horizon_number = 1;
  int N_points = 1;
  HorizonRadiusInDirection(horizon_number,N_points,&horizdirn_x, &horizdirn_y, &horizdirn_z,&output_radius);
}

extern "C" void CCTK_FCALL get_ah_radius_in_dirn_
  (const cGH **cctkGH,const double &horizdirn_x,const double &horizdirn_y,const double &horizdirn_z,double &output_radius)
{  
  get_ah_radius_in_dirn(*cctkGH,horizdirn_x,horizdirn_y,horizdirn_z,output_radius);
}
