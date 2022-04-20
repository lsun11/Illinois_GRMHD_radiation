//-----------------------------------------------------------------------
// Regridder: Read one line of BHNS grid input file.
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

extern "C" void CCTK_FCALL CCTK_FNAME(bhns_regridder_out_1gf)
  (int &numpoints,double *output,int &checksum);

extern "C" void CCTK_FCALL CCTK_FNAME(bhns_regridder_in_1gf)
  (int &numpoints,double *input,int &checksum);

extern "C" void bhns_regridder_out_1gf(int &numpoints,double *output,int &checksum) {
  char filename[100];
  sprintf(filename,"ini_bhns-regrid_%d.dat",checksum);

  printf("Attempting to write %s now...  \n",filename);

  ofstream outfile;
  outfile.open (filename, ios::out | ios::app | ios::binary);

  outfile.write ((char *)&numpoints, sizeof(int));
  outfile.write ((char *)&output, sizeof(double)*numpoints);

  outfile.close();
}

extern "C" void CCTK_FCALL CCTK_FNAME(bhns_regridder_out_1gf) 
  (int &numpoints,double *output,int &checksum)
{  
  bhns_regridder_out_1gf(numpoints,output,checksum);
}

extern "C" void bhns_regridder_in_1gf(int &numpoints,double *input,int &checksum) {
  char filename[100];
  sprintf(filename,"ini_bhns-regrid_%d.dat",checksum);

  printf("Attempting to write %s now...  \n",filename);

  ifstream infile;
  infile.open (filename, ios::out | ios::app | ios::binary);

  infile.read ((char *)&numpoints, sizeof(int));
  infile.read ((char *)&input, sizeof(double)*numpoints);

  infile.close();
}

extern "C" void CCTK_FCALL CCTK_FNAME(bhns_regridder_in_1gf)
  (int &numpoints,double *input,int &checksum)
{
  bhns_regridder_in_1gf(numpoints,input,checksum);
}
