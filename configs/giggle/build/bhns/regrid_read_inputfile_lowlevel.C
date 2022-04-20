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

extern "C" void CCTK_FCALL bhns_regrid_read_one_line_of_inputfile_
  (int &which_line,
   double &xmin_interp,double &ymin_interp,double &zmin_interp,
   double &dx_interp,double &dy_interp,double &dz_interp,
   int &Nx_interp,int &Ny_interp,int &Nz_interp,
   int &checksum);

extern "C" void bhns_regrid_read_one_line_of_inputfile(int &which_line,
						       double &xmin_interp,double &ymin_interp,double &zmin_interp,
						       double &dx_interp,double &dy_interp,double &dz_interp,
						       int &Nx_interp,int &Ny_interp,int &Nz_interp,
						       int &checksum) {
  

  char filename[100];
  sprintf(filename,"ini_bhns-grids.d",checksum);

  printf("Attempting to read grid input file %s now...  \n",filename);

  ifstream infile1;
  infile1.open(filename);
  if(!infile1) {
    cerr << "\a Can't open " << filename << " for input." << endl;
    exit(1);
  }

  int current_line=1;
  while (!infile1.eof()) {
    /* Each line of file has this format:
     * ./initdata_bhns -2.1250000000000000e+01 -2.0400000000000006e+01 -1.6999999999999993e+00 4.2499999999999999e-01 4.2499999999999999e-01 4.2499999999999999e-01 57 53 53 4825353
     */

    char dummy[100];
    infile1 >> dummy;
    // In the case that there's an empty line at the end of the file.
    if(infile1.eof()) { checksum=-1; return; }

    infile1 >> xmin_interp;
    infile1 >> ymin_interp;
    infile1 >> zmin_interp;

    infile1 >> dx_interp;
    infile1 >> dy_interp;
    infile1 >> dz_interp;

    infile1 >> Nx_interp;
    infile1 >> Ny_interp;
    infile1 >> Nz_interp;

    infile1 >> checksum;

    printf("which line=%d ; current line = %d\n",current_line,which_line);

    if(current_line==which_line) {     
      printf("Read line: ./initdata_bhns %.15e %.15e %.15e %.15e %.15e %.15e %d %d %d %d\n",xmin_interp,ymin_interp,zmin_interp,dx_interp,dy_interp,dz_interp,Nx_interp,Ny_interp,Nz_interp,checksum);
      infile1.close();  return; }
    current_line++;
  }

  infile1.close();
  
  checksum=-1;
  return;
}

extern "C" void CCTK_FCALL bhns_regrid_read_one_line_of_inputfile_
  (int &which_line,
   double &xmin_interp,double &ymin_interp,double &zmin_interp,
   double &dx_interp,double &dy_interp,double &dz_interp,
   int &Nx_interp,int &Ny_interp,int &Nz_interp,
   int &checksum)
{  
  bhns_regrid_read_one_line_of_inputfile(which_line,
					 xmin_interp,ymin_interp,zmin_interp,
					 dx_interp,dy_interp,dz_interp,
					 Nx_interp,Ny_interp,Nz_interp,
					 checksum);
}
