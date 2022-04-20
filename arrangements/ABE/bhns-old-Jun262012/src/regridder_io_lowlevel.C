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
#include "cctk_Functions.h"
#include "util_Table.h"

using namespace std;



extern "C" void CCTK_FCALL CCTK_FNAME(bhns_regridder_out_allgfs)
  (const cGH **cctkGH,int &checksum,
   double &xmin_interp,double &ymin_interp,double &zmin_interp,
   double &dx_interp,double &dy_interp,double &dz_interp,
   int &Nx_interp,int &Ny_interp,int &Nz_interp,
   double &rho_b_atm,double &tau_atm);

extern "C" void bhns_regridder_out_allgfs(const cGH *cctkGH,int &checksum,
					  double &xmin_interp,double &ymin_interp,double &zmin_interp,
					  double &dx_interp,double &dy_interp,double &dz_interp,
					  int &Nx_interp,int &Ny_interp,int &Nz_interp,
					  double &rho_b_atm,double &tau_atm) {

  void bhns_regridder_out_1gf(int &numpoints,double *output,int &checksum);
  void bhns_interp_driver_carp(const cGH *cctkGH,int &num_points,double *pointcoordsx,double *pointcoordsy,double *pointcoordsz,int &gridfunc_varindex,double *output);


  int numpoints = Nx_interp*Ny_interp*Nz_interp;

  printf("hello! numpoints = %d\n",numpoints);


  double *output = (double *)malloc(sizeof(double)*numpoints);
  double *x_interp = (double *)malloc(sizeof(double)*numpoints);
  double *y_interp = (double *)malloc(sizeof(double)*numpoints);
  double *z_interp = (double *)malloc(sizeof(double)*numpoints);

  int whichpoint,vi;

  /* 
     Restagger A^{\mu}. Before this function was called, we scheduled regridder_unstagger_driver() and
     saved the unstaggered A^{\mu} = Ax,Ay,Az,psi6phi to Bxtilde,Bytilde,Bztilde, and 
     Blagrangemultiplier, respectively. Here we restagger these gridfunctions.
     
     Note that restaggering is necessary because with the staggered A^{\mu}'s, the carpet interpolator assumes
     that coarse gridpoints overlap with fine gridpoints, which is WRONG. Thus without restaggering, we get
     WRONG values for A^{\mu}
   */
  double xoff,yoff,zoff;
  double xmax=xmin_interp+dx_interp*Nx_interp;
  double ymax=ymin_interp+dy_interp*Ny_interp;
  double zmax=zmin_interp+dz_interp*Nz_interp;
  // Restagger Ax:
  xoff=0.0*dx_interp;yoff=0.5*dy_interp;zoff=0.5*dz_interp;
  whichpoint=0; for(int k=0;k<Nz_interp;k++) for(int j=0;j<Ny_interp;j++) for(int i=0;i<Nx_interp;i++) {
  	x_interp[whichpoint] = xmin_interp + i*dx_interp+xoff; y_interp[whichpoint] = ymin_interp + j*dy_interp+yoff; z_interp[whichpoint] = zmin_interp + k*dz_interp+zoff; whichpoint++; }
  vi=CCTK_VarIndex("mhd_evolve::Bxtilde");bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);

  // Restagger Ay:
  xoff=0.5*dx_interp;yoff=0.0*dy_interp;zoff=0.5*dz_interp;
  whichpoint=0; for(int k=0;k<Nz_interp;k++) for(int j=0;j<Ny_interp;j++) for(int i=0;i<Nx_interp;i++) {
  	x_interp[whichpoint] = xmin_interp + i*dx_interp+xoff; y_interp[whichpoint] = ymin_interp + j*dy_interp+yoff; z_interp[whichpoint] = zmin_interp + k*dz_interp+zoff; whichpoint++; }
  vi=CCTK_VarIndex("mhd_evolve::Bytilde"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  //for(int i=0;i<numpoints;i++) { if(fabs(output[i])>1e-4) printf("C:: %d\t%e\n",i,output[i]); }

  // Restagger Az:
  xoff=0.5*dx_interp;yoff=0.5*dy_interp;zoff=0.0*dz_interp;
  whichpoint=0; for(int k=0;k<Nz_interp;k++) for(int j=0;j<Ny_interp;j++) for(int i=0;i<Nx_interp;i++) {
  	x_interp[whichpoint] = xmin_interp + i*dx_interp+xoff; y_interp[whichpoint] = ymin_interp + j*dy_interp+yoff; z_interp[whichpoint] = zmin_interp + k*dz_interp+zoff; whichpoint++; }
  vi=CCTK_VarIndex("mhd_evolve::Bztilde"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);

  // Restagger psi6phi:
  xoff=0.5*dx_interp;yoff=0.5*dy_interp;zoff=0.5*dz_interp;
  whichpoint=0; for(int k=0;k<Nz_interp;k++) for(int j=0;j<Ny_interp;j++) for(int i=0;i<Nx_interp;i++) {
  	x_interp[whichpoint] = xmin_interp + i*dx_interp+xoff; y_interp[whichpoint] = ymin_interp + j*dy_interp+yoff; z_interp[whichpoint] = zmin_interp + k*dz_interp+zoff; whichpoint++; }
  vi=CCTK_VarIndex("mhd_evolve::Blagrangemultiplier"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  printf("CHECK: %e\t%e\t%e! %e\t%e\t%e\n",x_interp[1],y_interp[1],z_interp[1],xmin_interp,ymin_interp,zmin_interp);
  printf("CHECK2: %e\t%e\t%e!\n",xoff,yoff,zoff);


  whichpoint=0; for(int k=0;k<Nz_interp;k++) for(int j=0;j<Ny_interp;j++) for(int i=0;i<Nx_interp;i++) {
	x_interp[whichpoint] = xmin_interp + i*dx_interp; y_interp[whichpoint] = ymin_interp + j*dy_interp; z_interp[whichpoint] = zmin_interp + k*dz_interp; whichpoint++; }
  vi=CCTK_VarIndex("mhd_evolve::rho_star");bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  whichpoint=0; for(int k=0;k<Nz_interp;k++) for(int j=0;j<Ny_interp;j++) for(int i=0;i<Nx_interp;i++) {
	x_interp[whichpoint] = xmin_interp + i*dx_interp; y_interp[whichpoint] = ymin_interp + j*dy_interp; z_interp[whichpoint] = zmin_interp + k*dz_interp; 
	if(i==42 && j==39 && k==4) printf("interped rho star  tau: %.15e\n",output[whichpoint]);
	whichpoint++; }
  vi=CCTK_VarIndex("mhd_evolve::tau"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  whichpoint=0; for(int k=0;k<Nz_interp;k++) for(int j=0;j<Ny_interp;j++) for(int i=0;i<Nx_interp;i++) {
	x_interp[whichpoint] = xmin_interp + i*dx_interp; y_interp[whichpoint] = ymin_interp + j*dy_interp; z_interp[whichpoint] = zmin_interp + k*dz_interp; 
	if(i==42 && j==39 && k==4) printf("interped tau: %.15e\n",output[whichpoint]);
	whichpoint++; }
  vi=CCTK_VarIndex("mhd_evolve::mhd_st_x"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  vi=CCTK_VarIndex("mhd_evolve::mhd_st_y"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  vi=CCTK_VarIndex("mhd_evolve::mhd_st_z"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);

  vi=CCTK_VarIndex("shift::shiftx"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  vi=CCTK_VarIndex("shift::shifty"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  vi=CCTK_VarIndex("shift::shiftz"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  vi=CCTK_VarIndex("lapse::lapm1"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);

  vi=CCTK_VarIndex("bssn::phi"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output);  bhns_regridder_out_1gf(numpoints,output,checksum);
  vi=CCTK_VarIndex("bssn::chi"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  vi=CCTK_VarIndex("bssn::trK"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);

  vi=CCTK_VarIndex("bssn::gxx"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  vi=CCTK_VarIndex("bssn::gxy"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  vi=CCTK_VarIndex("bssn::gxz"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  vi=CCTK_VarIndex("bssn::gyy"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  vi=CCTK_VarIndex("bssn::gyz"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  vi=CCTK_VarIndex("bssn::gzz"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);

  vi=CCTK_VarIndex("bssn::Axx"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  vi=CCTK_VarIndex("bssn::Axy"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  vi=CCTK_VarIndex("bssn::Axz"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  vi=CCTK_VarIndex("bssn::Ayy"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  vi=CCTK_VarIndex("bssn::Ayz"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  vi=CCTK_VarIndex("bssn::Azz"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);

  vi=CCTK_VarIndex("bssn::Gammax"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  vi=CCTK_VarIndex("bssn::Gammay"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);
  vi=CCTK_VarIndex("bssn::Gammaz"); bhns_interp_driver_carp(cctkGH,numpoints,x_interp,y_interp,z_interp,vi,output); bhns_regridder_out_1gf(numpoints,output,checksum);

  numpoints=1;
  bhns_regridder_out_1gf(numpoints,&rho_b_atm,checksum);
  bhns_regridder_out_1gf(numpoints,&tau_atm,checksum);


  free(output);
  free(x_interp);
  free(y_interp);
  free(z_interp);
}
extern "C" void CCTK_FCALL CCTK_FNAME(bhns_regridder_out_allgfs)
  (const cGH **cctkGH,int &checksum,
   double &xmin_interp,double &ymin_interp,double &zmin_interp,
   double &dx_interp,double &dy_interp,double &dz_interp,
   int &Nx_interp,int &Ny_interp,int &Nz_interp,
   double &rho_b_atm,double &tau_atm) {

  bhns_regridder_out_allgfs(*cctkGH,checksum,
			    xmin_interp,ymin_interp,zmin_interp,
			    dx_interp,dy_interp,dz_interp,
			    Nx_interp,Ny_interp,Nz_interp,
			    rho_b_atm,tau_atm);
}

void bhns_regridder_out_1gf(int &numpoints,double *output,int &checksum) {

  char filename[100];
  sprintf(filename,"ini_bhns-regrid_%d.dat",checksum);

  printf("Attempting to write %d points to %s now. You'll see this message N times for N gridfunctions.\n",numpoints,filename);

  ofstream outfile;
  outfile.open (filename, ios::out | ios::app | ios::binary);

  outfile.write ((char *)&numpoints, sizeof(int));
  outfile.write ((char *)output, sizeof(double)*numpoints);

  outfile.close();
}

extern "C" void CCTK_FCALL CCTK_FNAME(bhns_regridder_in_allgfs)
  (int &numpoints,int &checksum,double *Ax,double *Ay,double *Az,double *psi6phi,
   double *rho_star,double *tau,double *mhd_st_x,double *mhd_st_y,double *mhd_st_z,
   double *shiftx,double *shifty,double *shiftz,double *lapm1,
   double *phi,double *chi,double *trK,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *Axx,double *Axy,double *Axz,double *Ayy,double *Ayz,double *Azz,
   double *Gammax,double *Gammay,double *Gammaz,
   double &rho_b_atm,double &tau_atm);

extern "C" void bhns_regridder_in_allgfs(int &numpoints,int &checksum,double *Ax,double *Ay,double *Az,double *psi6phi,
					 double *rho_star,double *tau,double *mhd_st_x,double *mhd_st_y,double *mhd_st_z,
					 double *shiftx,double *shifty,double *shiftz,double *lapm1,
					 double *phi,double *chi,double *trK,
					 double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
					 double *Axx,double *Axy,double *Axz,double *Ayy,double *Ayz,double *Azz,
					 double *Gammax,double *Gammay,double *Gammaz,
					 double &rho_b_atm,double &tau_atm) {
  char filename[100];
  sprintf(filename,"ini_bhns-regrid_%d.dat",checksum);

  printf("Attempting to read %s now...  \n",filename);

  ifstream infile;
  infile.open (filename, ios::in | ios::binary);

  int numpoints_file;

  printf("Attempting to read magnetic field data...\n");
  // First do a sanity check: do the files contain the same # of datapoints as this grid?
  infile.read ((char *)&numpoints_file, sizeof(int)); 
  if(numpoints_file!=numpoints) { 
    printf("ERROR. This file contains %d points, was expecting %d!\n",numpoints_file,numpoints); exit(1); 
  }
  printf("Good. The files contain the same # of datapoints as this grid %d!\n",numpoints);
                                                 infile.read ((char *)Ax, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)Ay, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)Az, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)psi6phi, sizeof(double)*numpoints);

  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)rho_star, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)tau, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)mhd_st_x, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)mhd_st_y, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)mhd_st_z, sizeof(double)*numpoints);

  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)shiftx, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)shifty, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)shiftz, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)lapm1, sizeof(double)*numpoints);

  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)phi, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)chi, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)trK, sizeof(double)*numpoints);

  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)gxx, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)gxy, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)gxz, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)gyy, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)gyz, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)gzz, sizeof(double)*numpoints);

  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)Axx, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)Axy, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)Axz, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)Ayy, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)Ayz, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)Azz, sizeof(double)*numpoints);

  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)Gammax, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)Gammay, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)Gammaz, sizeof(double)*numpoints);

  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)&rho_b_atm, sizeof(double)*numpoints);
  infile.read ((char *)&numpoints, sizeof(int)); infile.read ((char *)&tau_atm, sizeof(double)*numpoints);

  printf("hi. rho_b_atm = %e, tau_atm = %e\n",rho_b_atm,tau_atm);

  printf("Finished reading %s ! \n",filename);

  infile.close();
  printf("Finished reading %s ! \n",filename);

  return;

}

extern "C" void CCTK_FCALL CCTK_FNAME(bhns_regridder_in_allgfs)
  (int &numpoints,int &checksum,double *Ax,double *Ay,double *Az,double *psi6phi,
   double *rho_star,double *tau,double *mhd_st_x,double *mhd_st_y,double *mhd_st_z,
   double *shiftx,double *shifty,double *shiftz,double *lapm1,
   double *phi,double *chi,double *trK,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *Axx,double *Axy,double *Axz,double *Ayy,double *Ayz,double *Azz,
   double *Gammax,double *Gammay,double *Gammaz,
   double &rho_b_atm,double &tau_atm)
{
  bhns_regridder_in_allgfs(numpoints,checksum,Ax,Ay,Az,psi6phi,
			   rho_star,tau,mhd_st_x,mhd_st_y,mhd_st_z,
			   shiftx,shifty,shiftz,lapm1,
			   phi,chi,trK,
			   gxx,gxy,gxz,gyy,gyz,gzz,
			   Axx,Axy,Axz,Ayy,Ayz,Azz,
			   Gammax,Gammay,Gammaz,
			   rho_b_atm,tau_atm);
}


void bhns_interp_driver_carp(const cGH *cctkGH,int &num_points,double *pointcoordsx,double *pointcoordsy,double *pointcoordsz,int &gridfunc_varindex,double *output) {
  int N_dims = 3;
  int N_interp_points = num_points;
  int N_input_arrays = 1;
  int N_output_arrays = 1;

  CCTK_INT input_array_indices[1];
  input_array_indices[0] = gridfunc_varindex;

  int interp_handle = CCTK_InterpHandle ("Lagrange polynomial interpolation");
  if (interp_handle < 0) {
    CCTK_WARN(0,"Cannot get handle for interpolation // Forgot to activate an implementation providing interpolation operators ??");
  }
  
  int coord_system_handle = CCTK_CoordSystemHandle("cart3d");
  if (coord_system_handle < 0) {
    CCTK_WARN(0,"Cannot get handle for cart3d coordinate system // Forgot to activate an implementation providing coordinates ??");
  }

  static const CCTK_INT output_array_type_codes[1] = { CCTK_VARIABLE_REAL };  

  void *output_arrays[1];
  output_arrays[0] = (void *) output;

  const void *interp_coords[3];
  interp_coords[0] = (const void *) pointcoordsx;
  interp_coords[1] = (const void *) pointcoordsy;
  interp_coords[2] = (const void *) pointcoordsz;

  int status =
    CCTK_InterpGridArrays(cctkGH,
                          N_dims,
                          interp_handle, Util_TableCreateFromString("order=4"),
                          coord_system_handle,
                          N_interp_points,
                          CCTK_VARIABLE_REAL,
                          interp_coords,
                          N_input_arrays,
                          input_array_indices,
                          N_output_arrays,
                          output_array_type_codes,
                          output_arrays);
}

