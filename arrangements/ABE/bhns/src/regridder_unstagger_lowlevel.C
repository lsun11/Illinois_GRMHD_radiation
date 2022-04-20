//----------------------
// Regridder: Unstagger A's
//----------------------
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

void bhns_local_interp(double *input,
		       const cGH *cctkGH,
		       int *cctk_lsh,
		       double &dX,double &dY,double &dZ,
		       double * X,double * Y,double * Z,
		       double *x_interp, double *y_interp, double *z_interp,
		       double *output);


extern "C" void CCTK_FCALL CCTK_FNAME(bhns_regridder_unstagger)
  (const cGH **cctkGH,int *cctk_lsh,int *cctk_nghostzones,
   double *X,double *Y,double *Z,double &dX,double &dY,double &dZ,
   double *Ax_s,double *Ay_s,double *Az_s,double *psi6phi_s,
   double *Ax_u,double *Ay_u,double *Az_u,double *psi6phi_u);

extern "C" void bhns_regridder_unstagger(const cGH *cctkGH,int *cctk_lsh,int *cctk_nghostzones,
					 double *X,double *Y,double *Z,double &dX,double &dY,double &dZ,
					 double *Ax_s,double *Ay_s,double *Az_s,double *psi6phi_s,
					 double *Ax_u,double *Ay_u,double *Az_u,double *psi6phi_u) {

  printf("INSIDE REGRIDDER_UNSTAGGER!\n");
  
  int numpoints = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];

  double *x_interp = (double *)malloc(sizeof(double)*numpoints);
  double *y_interp = (double *)malloc(sizeof(double)*numpoints);
  double *z_interp = (double *)malloc(sizeof(double)*numpoints);

  int whichpoint;
  double xoff,yoff,zoff;

  xoff = 0.0*dX,yoff = 0.5*dY,zoff = 0.5*dZ;
  whichpoint=0;  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
	double x_origin = X[CCTK_GFINDEX3D(cctkGH,0,0,0)] - xoff, y_origin = Y[CCTK_GFINDEX3D(cctkGH,0,0,0)] - yoff, z_origin = Z[CCTK_GFINDEX3D(cctkGH,0,0,0)] - zoff;
	x_interp[whichpoint] = x_origin + dX*i;
	y_interp[whichpoint] = y_origin + dY*j;
	z_interp[whichpoint] = z_origin + dZ*k;
	if(x_interp[whichpoint] < X[CCTK_GFINDEX3D(cctkGH,0,0,0)]) x_interp[whichpoint] = X[CCTK_GFINDEX3D(cctkGH,0,0,0)];
	if(y_interp[whichpoint] < Y[CCTK_GFINDEX3D(cctkGH,0,0,0)]) y_interp[whichpoint] = Y[CCTK_GFINDEX3D(cctkGH,0,0,0)];
	if(z_interp[whichpoint] < Z[CCTK_GFINDEX3D(cctkGH,0,0,0)]) z_interp[whichpoint] = Z[CCTK_GFINDEX3D(cctkGH,0,0,0)];
	whichpoint++;
      }
  bhns_local_interp(Ax_s,cctkGH,cctk_lsh,dX,dY,dZ,X,Y,Z,x_interp,y_interp,z_interp,Ax_u);

  xoff = 0.5*dX,yoff = 0.0*dY,zoff = 0.5*dZ;
  whichpoint=0;  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
	double x_origin = X[CCTK_GFINDEX3D(cctkGH,0,0,0)] - xoff, y_origin = Y[CCTK_GFINDEX3D(cctkGH,0,0,0)] - yoff, z_origin = Z[CCTK_GFINDEX3D(cctkGH,0,0,0)] - zoff;
	x_interp[whichpoint] = x_origin + dX*i;
	y_interp[whichpoint] = y_origin + dY*j;
	z_interp[whichpoint] = z_origin + dZ*k;
	if(x_interp[whichpoint] < X[CCTK_GFINDEX3D(cctkGH,0,0,0)]) x_interp[whichpoint] = X[CCTK_GFINDEX3D(cctkGH,0,0,0)];
	if(y_interp[whichpoint] < Y[CCTK_GFINDEX3D(cctkGH,0,0,0)]) y_interp[whichpoint] = Y[CCTK_GFINDEX3D(cctkGH,0,0,0)];
	if(z_interp[whichpoint] < Z[CCTK_GFINDEX3D(cctkGH,0,0,0)]) z_interp[whichpoint] = Z[CCTK_GFINDEX3D(cctkGH,0,0,0)];
	whichpoint++;
      }
  bhns_local_interp(Ay_s,cctkGH,cctk_lsh,dX,dY,dZ,X,Y,Z,x_interp,y_interp,z_interp,Ay_u);

  xoff = 0.5*dX,yoff = 0.5*dY,zoff = 0.0*dZ;
  whichpoint=0;  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
	double x_origin = X[CCTK_GFINDEX3D(cctkGH,0,0,0)] - xoff, y_origin = Y[CCTK_GFINDEX3D(cctkGH,0,0,0)] - yoff, z_origin = Z[CCTK_GFINDEX3D(cctkGH,0,0,0)] - zoff;
	x_interp[whichpoint] = x_origin + dX*i;
	y_interp[whichpoint] = y_origin + dY*j;
	z_interp[whichpoint] = z_origin + dZ*k;
	if(x_interp[whichpoint] < X[CCTK_GFINDEX3D(cctkGH,0,0,0)]) x_interp[whichpoint] = X[CCTK_GFINDEX3D(cctkGH,0,0,0)];
	if(y_interp[whichpoint] < Y[CCTK_GFINDEX3D(cctkGH,0,0,0)]) y_interp[whichpoint] = Y[CCTK_GFINDEX3D(cctkGH,0,0,0)];
	if(z_interp[whichpoint] < Z[CCTK_GFINDEX3D(cctkGH,0,0,0)]) z_interp[whichpoint] = Z[CCTK_GFINDEX3D(cctkGH,0,0,0)];
	whichpoint++;
      }
  bhns_local_interp(Az_s,cctkGH,cctk_lsh,dX,dY,dZ,X,Y,Z,x_interp,y_interp,z_interp,Az_u);

  xoff = 0.5*dX,yoff = 0.5*dY,zoff = 0.5*dZ;
  whichpoint=0;  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
	double x_origin = X[CCTK_GFINDEX3D(cctkGH,0,0,0)] - xoff, y_origin = Y[CCTK_GFINDEX3D(cctkGH,0,0,0)] - yoff, z_origin = Z[CCTK_GFINDEX3D(cctkGH,0,0,0)] - zoff;
	x_interp[whichpoint] = x_origin + dX*i;
	y_interp[whichpoint] = y_origin + dY*j;
	z_interp[whichpoint] = z_origin + dZ*k;
	if(x_interp[whichpoint] < X[CCTK_GFINDEX3D(cctkGH,0,0,0)]) x_interp[whichpoint] = X[CCTK_GFINDEX3D(cctkGH,0,0,0)];
	if(y_interp[whichpoint] < Y[CCTK_GFINDEX3D(cctkGH,0,0,0)]) y_interp[whichpoint] = Y[CCTK_GFINDEX3D(cctkGH,0,0,0)];
	if(z_interp[whichpoint] < Z[CCTK_GFINDEX3D(cctkGH,0,0,0)]) z_interp[whichpoint] = Z[CCTK_GFINDEX3D(cctkGH,0,0,0)];
	whichpoint++;
      }
  bhns_local_interp(psi6phi_s,cctkGH,cctk_lsh,dX,dY,dZ,X,Y,Z,x_interp,y_interp,z_interp,psi6phi_u);

  printf("check: %e\t%e\t%e! %e\t%e\t%e\n",x_interp[1],y_interp[1],z_interp[1],
	 X[CCTK_GFINDEX3D(cctkGH,0,0,0)],
	 Y[CCTK_GFINDEX3D(cctkGH,0,0,0)],
	 Z[CCTK_GFINDEX3D(cctkGH,0,0,0)]);
  printf("check2: %e\t%e\t%e!\n",xoff,yoff,zoff);  

  free(x_interp);
  free(y_interp);
  free(z_interp);
  printf("howdy2!\n");
}

extern "C" void CCTK_FCALL CCTK_FNAME(bhns_regridder_unstagger)
  (const cGH **cctkGH,int *cctk_lsh,int *cctk_nghostzones,
   double *X,double *Y,double *Z,double &dX,double &dY,double &dZ,
   double *Ax_s,double *Ay_s,double *Az_s,double *psi6phi_s,
   double *Ax_u,double *Ay_u,double *Az_u,double *psi6phi_u)
{
  bhns_regridder_unstagger
    (*cctkGH,cctk_lsh,cctk_nghostzones,
     X,Y,Z,dX,dY,dZ,
     Ax_s,Ay_s,Az_s,psi6phi_s,
     Ax_u,Ay_u,Az_u,psi6phi_u);
}


void bhns_local_interp(double *input,
		       const cGH *cctkGH,
		       int *cctk_lsh,
		       double &dX,double &dY,double &dZ,
		       double * X,double * Y,double * Z,
		       double *x_interp, double *y_interp, double *z_interp,
		       double *output) {
  int N_dims = 3;

  int operator_handle = CCTK_InterpHandle ("Lagrange polynomial interpolation");
  if (operator_handle < 0) {
    CCTK_WARN(0,"Cannot get handle for interpolation // Forgot to activate an implementation providing interpolation operators ??");
  }

  const CCTK_REAL origin[3] = { X[CCTK_GFINDEX3D(cctkGH,0,0,0)], Y[CCTK_GFINDEX3D(cctkGH,0,0,0)], Z[CCTK_GFINDEX3D(cctkGH,0,0,0)] };
  const CCTK_REAL delta [3] = { dX, dY, dZ };

  int N_interp_points = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];

  const void *const interp_coords[3] = { (const void *)x_interp, (const void *)y_interp, (const void *)z_interp };

  int N_input_arrays = 1;

  const CCTK_INT input_array_dims[3] = { cctk_lsh[0], cctk_lsh[1], cctk_lsh[2] };

  const CCTK_INT input_array_type_codes[1] = { CCTK_VARIABLE_REAL };

  const void *input_arrays[1] = { (const void *) input } ;       /* see note above */

  int N_output_arrays = 1;

  const CCTK_INT output_array_type_codes[1] = { CCTK_VARIABLE_REAL };

  void *const output_arrays[1] = { (void *)output };

  /*CCTK_InterpLocalUniform(int N_dims,
                             int operator_handle,
                             int param_table_handle,
                             const CCTK_REAL coord_origin[],
                             const CCTK_REAL coord_delta[],
                             int N_interp_points,
                               int interp_coords_type_code,
                               const void *const interp_coords[],
                             int N_input_arrays,
                               const CCTK_INT input_array_dims[],
                               const CCTK_INT input_array_type_codes[],
                               const void *const input_arrays[],
                             int N_output_arrays,
                               const CCTK_INT output_array_type_codes[],
                               void *const output_arrays[]);
  */

  /* do the actual interpolation, and check for error returns */
  CCTK_InterpLocalUniform(N_dims,
			      operator_handle, Util_TableCreateFromString("order=4"),
			      origin, delta,
			      N_interp_points,
			      CCTK_VARIABLE_REAL,
			      interp_coords,
			      N_input_arrays,
			      input_array_dims,
			      input_array_type_codes,
			      input_arrays,
			      N_output_arrays,
			      output_array_type_codes,
			      output_arrays);
  /*
  if (CCTK_InterpLocalUniform(N_dims,
			      operator_handle, Util_TableCreateFromString("order=4"),
			      origin, delta,
			      N_interp_points,
			      CCTK_VARIABLE_REAL,
			      interp_coords,
			      N_input_arrays,
			      input_array_dims,
			      input_array_type_codes,
			      input_arrays,
			      N_output_arrays,
			      output_array_type_codes,
			      output_arrays) < 0)
    CCTK_WARN(CCTK_WARN_ABORT, "error return from interpolator!");
*/
}
