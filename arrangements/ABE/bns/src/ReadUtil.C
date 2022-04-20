//-----------------------------------------------------------------------------
//
// $Id: ReadUtil.C,v 1.1.1.1 2006/02/23 17:48:41 zetienne Exp $
//
//-----------------------------------------------------------------------------
//
// Contains several routines to read data from files
//
//-----------------------------------------------------------------------------
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <sstream>

void convert_double(double *d)
{
  char temp[8];
  int i;


  for(i=0;i<8;i++)
    temp[i]=((char*)d)[7-i];
  for(i=0;i<8;i++)
    ((char*)d)[i]=temp[i];
}

void convert_int(int *d)
{
  char temp[4];
  int i;

  for(i=0;i<4;i++)
    temp[i]=((char*)d)[3-i];
  for(i=0;i<4;i++)
    ((char*)d)[i]=temp[i];
}

//====================================================================
//
// Function reads data from file into array data.
//
// This routine reads the size of the array nx, ny, nz from "filename", 
// allocates the array "data" accordingly, and reads the data from the 
// file "filename".
//
//====================================================================
/*
		 double & xmin, double & xmax, double & ymin, 
		 double & ymax, double & zmin, double & zmax, 
*/
void ReadtoArray(const char * filename, double ** data,
		 double & xmin, double & xmax, double & ymin, 
		 double & ymax, double & zmin, double & zmax, 
		 int & nx, int & ny, int & nz)
{
  //==================================================================
  // open file
  //==================================================================
  fstream file;
  file.open(filename, ios::in);
  if (!file) { printf("Exiting: Cannot open file %s \n",filename); exit(0); }
  //==================================================================
  // Get time or re and grid parameters from file
  //==================================================================
  file.seekg(0, ios::beg);
  double re = 0.0;
  file.read((char *) &re, sizeof(double));
  file.read((char *) &nx, sizeof(int));
  file.read((char *) &ny, sizeof(int));
  file.read((char *) &nz, sizeof(int));
  file.read((char *) &xmin, sizeof(double));
  file.read((char *) &xmax, sizeof(double));
  file.read((char *) &ymin, sizeof(double));
  file.read((char *) &ymax, sizeof(double));
  file.read((char *) &zmin, sizeof(double));
  file.read((char *) &zmax, sizeof(double));
//
// CONVERSIONS 
// 
  convert_double(&re);
  convert_int(&nx);
  convert_int(&ny);
  convert_int(&nz);
  convert_double(&xmin);
  convert_double(&xmax);
  convert_double(&ymin);
  convert_double(&ymax);
  convert_double(&zmin);
  convert_double(&zmax);

  //==================================================================
  // Rescale limits
  //==================================================================
  xmin *= re;
  ymin *= re;
  zmin *= re;
  xmax *= re;
  ymax *= re;
  zmax *= re;
  //==================================================================
  // Allocate Array and read data
  //==================================================================
  *data = new double[nx*ny*nz];
  file.read((char *) *data, nx*ny*nz*sizeof(double));
// CONVERSIONS
  for(int i = 0; i < nx*ny*nz; i++) convert_double(&((*data)[i]));
//
  file.close();
} 

//====================================================================
//
// Function reads data from file into gridfunction.
//
// It compares the extends of the gridfunction grid and the data in
// file and then decides if data should be reflected across coordinate
// planes.  If so, it does so using the symmetry properties specified
// by xsym, ysym and zsym.  It also returns the rescaling paramter re,
// which is needed to rescale the coordinates.
//
//====================================================================

// function prototype:
void index(const int i, const int j, const int k,
	  const int Nx, const int Ny, const int Nz,
	  const int xflag, const int yflag, const int zflag,
	  const double xsym, const double ysym, const double zsym, 
	  double sign, int iind, int jind, int kind);

//====================================================================
//
// Function determines grid size from file.
//
// It opens one file to find nx, ny, nz from file, compares these
// with the specified values, decides if data should be reflected 
// across coordinate planes, and accordingly sets xmin, xmax...
//
// To be called from Read_Input, runs only on processor 0 (me == 0)
//
//====================================================================

void FindGridSize(const char * filename, const int me, 
		  const int Nx, const int Ny, const int Nz, 
		  double Xmin, double Xmax, double Ymin,
		  double Ymax, double Zmin, double Zmax)
{
  if (me == 0) {
    //==================================================================
    // open file
    //==================================================================
    fstream file;
    file.open(filename, ios::in);
    if (!file) { printf("Exiting: Cannot open file %s \n",filename); exit(0); }
    //==================================================================
    // Get time or re and grid parameters from file
    //==================================================================
    int nx, ny, nz;
    double re;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    file.seekg(0, ios::beg);
    file.read((char *) &re, sizeof(double));
    file.read((char *) &nx, sizeof(int));
    file.read((char *) &ny, sizeof(int));
    file.read((char *) &nz, sizeof(int));
    file.read((char *) &xmin, sizeof(double));
    file.read((char *) &xmax, sizeof(double));
    file.read((char *) &ymin, sizeof(double));
    file.read((char *) &ymax, sizeof(double));
    file.read((char *) &zmin, sizeof(double));
    file.read((char *) &zmax, sizeof(double));
    //==================================================================
    // Now look at what we got in x, y and z direction...
    //==================================================================
    if (xmin != 0.0) {
      cout << " Ooops -- xmin = " << xmin << endl;
      printf("xmin not what I thought it would be in file  %s \n",filename); exit(0);
    }
    if (ymin != 0.0) {
      cout << " Ooops -- ymin = " << ymin << endl;
      printf("ymin not what I thought it would be in file  %s \n",filename); exit(0);
    }
    if (zmin != 0.0) {
      cout << " Ooops -- zmin = " << zmin << endl;
      printf("zmin not what I thought it would be in file  %s \n",filename); exit(0);
    }
    //==================================================================
    // Set maxs...
    //==================================================================
    Xmax = xmax;
    Ymax = ymax;
    Zmax = zmax;
    //==================================================================
    // Set mins...
    //==================================================================
    if (Nx == nx)        { Xmin = xmin; }
    else if (Nx == 2*nx) { Xmin = - xmax; }
    else { printf("Nx and nx do not match!\n"); exit(0); }
    if (Ny == ny)        { Ymin = ymin; }
    else if (Ny == 2*ny) { Ymin = - ymax; }
    else { printf("Ny and ny do not match!\n"); exit(0); }
    if (Nz == nz)        { Zmin = zmin; }
    else if (Nz == 2*nz) { Zmin = - zmax; }
    else { printf("Nz and nz do not match!\n"); exit(0); }
    //==================================================================
    // ... and rescale everything
    //==================================================================
    Xmin *= re;
    Xmax *= re;
    Ymin *= re;
    Ymax *= re;
    Zmin *= re;
    Zmax *= re;
  }
}
  
//====================================================================
//
// Finds the correct index, taking into account reflections and such...
//
// NOTE: THIS ROUTINE ASSUMES A CELL-CENTERED GRID!!!
//
//====================================================================

void index(const int i, const int j, const int k,
	  const int Nx, const int Ny, const int Nz,
	  const int xflag, const int yflag, const int zflag,
	  const double xsym, const double ysym, const double zsym, 
	  double sign, int iind, int jind, int kind)
{
  sign = 1.0;
  //=======================================================
  // x:
  //=======================================================
  if (i >= Nx) { printf("i >= Nx in ReadUtil index!\n"); exit(0); } 
  if (xflag == 1) {           // fine
    iind = i; 
  } else if (xflag == 2) {    // refection possible
    int Nxhalf = Nx/2;
    if (i < Nxhalf) {         // need to reflect
      sign *= xsym;
      iind = Nx - 1 - i;
    } else {                  // don't need to reflect
      iind = i - Nxhalf;
    }
  }
  //=======================================================
  // y:
  //=======================================================
  if (j >= Ny) { printf("j >= Ny in ReadUtil index!\n"); exit(0); }
  if (yflag == 1) {           // fine
    jind = j; 
  } else if (yflag == 2) {    // refection possible
    int Nyhalf = Ny/2;
    if (j < Nyhalf) {         // need to reflect
      sign *= ysym;
      jind = Ny - 1 - j;
    } else {                  // don't need to reflect
      jind = j - Nyhalf;
    }
  }
  //=======================================================
  // z:
  //=======================================================
  if (k >= Nz) { printf("k >= Nz in ReadUtil index!\n"); exit(0); } 
  if (zflag == 1) {           // fine
    kind = k; 
  } else if (zflag == 2) {    // refection possible
    int Nzhalf = Nz/2;
    if (k < Nzhalf) {         // need to reflect
      sign *= zsym;
      kind = Nz - 1 - k;
    } else {                  // don't need to reflect
      kind = k - Nzhalf;
    }
  }
}
