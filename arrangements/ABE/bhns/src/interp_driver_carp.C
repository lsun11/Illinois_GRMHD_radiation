#include "stdio.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

//------------------------------------
// Interpolation 
//------------------------------------

extern "C" void CCTK_FCALL CCTK_FNAME(bhns_interp_driver_carp)
  (const cGH **cctkGH,int &num_points,double *pointcoordsx,double *pointcoordsy,double *pointcoordsz,int &gridfunc_varindex,double *output);


void bhns_interp_driver_carp(const cGH *cctkGH,int &num_points,double *pointcoordsx,double *pointcoordsy,double *pointcoordsz,int &gridfunc_varindex,double *output) {
  int N_dims = 3;
  int N_interp_points = num_points;
  int N_input_arrays = 1;
  int N_output_arrays = 1;

  const CCTK_INT input_array_indices[1] = { gridfunc_varindex };

  int interp_handle = CCTK_InterpHandle ("Lagrange polynomial interpolation");
  if (interp_handle < 0) {
    CCTK_WARN(0,"Cannot get handle for interpolation // Forgot to activate an implementation providing interpolation operators ??");
  }
  
  int coord_system_handle = CCTK_CoordSystemHandle("cart3d");
  if (coord_system_handle < 0) {
    CCTK_WARN(0,"Cannot get handle for cart3d coordinate system // Forgot to activate an implementation providing coordinates ??");
  }

  const int input_array_type_codes = CCTK_VARIABLE_REAL;
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
			  input_array_type_codes,
			  interp_coords,
			  N_input_arrays,
			  input_array_indices,
			  N_output_arrays,
			  output_array_type_codes,
			  output_arrays);
}

extern "C" void CCTK_FCALL CCTK_FNAME(bhns_interp_driver_carp)
  (const cGH **cctkGH,int &num_points,double *pointcoordsx,double *pointcoordsy,double *pointcoordsz,int &gridfunc_varindex,double *output) {
  bhns_interp_driver_carp(*cctkGH,num_points,pointcoordsx,pointcoordsy,pointcoordsz,gridfunc_varindex,output);
}
