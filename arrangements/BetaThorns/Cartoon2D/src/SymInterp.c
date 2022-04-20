/* $Header: /cactusdevcvs/BetaThorns/Cartoon2D/src/SymInterp.c,v 1.4 2005/12/07 09:56:11 hawke Exp $ */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "Cartoon2D.h"



/*
 * NOTE:
 *
 * The contents of this file have been largely copied from the thorn
 * AEIDevelopment/RotatingSymmetry90.  Please do not invest much
 * effort in changing this file here -- changes to its ancestor should
 * be easy to port over.  The only difference between this file and
 * its ancestor are a few lines in the coordinate rotation routines.
 * I hope that a generic tensor infrastructure will come up, and then
 * this duplication will go away.
 *
 * Erik Schnetter, 2004-05-22
 */




/* This is pretty hard coded into all the tensor types and cannot be
   changed easily.  */
#define DIM 3 /* spatial dimension */



/* These can be increased if necessary, but not all necessary tensor
   types may be defined .  */
#define MAX_RANK                9 /* maximum tensor rank */
#define MAX_TIME_LEVELS         3 /* maximum number of time levels */
#define MAX_SPATIAL_DERIV_ORDER 2 /* maximum spatial derivative order */
#define MAX_TIME_DERIV_ORDER    2 /* maximum time derivative order */



/* A tensor description */
struct tensor
{
  const char * name;            /* description */
  int dim;
  int rank;
  int ncomps;                   /* dim^rank */
  int * restrict vars;          /* map component to variable */
  int * restrict parity;        /* parity for the above */
  int nvars;                    /* depends on symmetries */
  int * restrict comps;         /* map variable to component */
};



/* A scalar */
static int scalar_vars[] = {
  0
};
static int scalar_parity[] = {
  +1
};
static int scalar_comps[] = {
  0
};
static struct tensor scalar = {
  "scalar",
  3, 0, 1, scalar_vars, scalar_parity, 1, scalar_comps
};



/* A vector */
static int vector_vars[] = {
  0,1,2
};
static int vector_parity[] = {
  +1,+1,+1
};
static int vector_comps[] = {
  0,1,2
};
static struct tensor vector = {
  "vector",
  3, 1, 3, vector_vars, vector_parity, 3, vector_comps
};



/* A second rank tensor without symmetries */
static int tensor_vars[] = {
  0,1,2,   3,4,5,   6,7,8
};
static int tensor_parity[] = {
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1
};
static int tensor_comps[] = {
  0,1,2,   3,4,5,   6,7,8
};
static struct tensor tensor = {
  "tensor",
  3, 2, 9, tensor_vars, tensor_parity, 9, tensor_comps
};



/* A symmetric second rank tensor */
static int symmtensor_vars[] = {
  0,1,2,   1,3,4,   2,4,5
};
static int symmtensor_parity[] = {
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1
};
static int symmtensor_comps[] = {
  0,1,2,   4,5,   8
};
static struct tensor symmtensor = {
  "symmetric tensor T_(ij)",
  3, 2, 9, symmtensor_vars, symmtensor_parity, 6, symmtensor_comps
};



/* A third rank tensor with symmetry in the first two indices */
static int symmtensor3a_vars[] = {
   0, 1, 2,
   3, 4, 5,
   6, 7, 8,
  
   3, 4, 5,
   9,10,11,
  12,13,14,
  
   6, 7, 8,
  12,13,14,
  15,16,17
};
static int symmtensor3a_parity[] = {
  +1,+1,+1,
  +1,+1,+1,
  +1,+1,+1,
  
  +1,+1,+1,
  +1,+1,+1,
  +1,+1,+1,
  
  +1,+1,+1,
  +1,+1,+1,
  +1,+1,+1,
};
static int symmtensor3a_comps[] = {
   0, 1, 2,
   3, 4, 5,
   6, 7, 8,
  
  12,13,14,
  15,16,17,
  
  24,25,26
};
static struct tensor symmtensor3a = {
  "symmetric tensor T_(ij)k",
  3, 3, 27, symmtensor3a_vars, symmtensor3a_parity, 18, symmtensor3a_comps
};



/* A third rank tensor with symmetry in the last two indices */
static int symmtensor3b_vars[] = {
   0, 1, 2,    1, 3, 4,    2, 4, 5,
   6, 7, 8,    7, 9,10,    8,10,11,
  12,13,14,   13,15,16,   14,16,17
};
static int symmtensor3b_parity[] = {
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1
};
static int symmtensor3b_comps[] = {
   0, 1, 2,    4, 5,    8,
   9,10,11,   13,14,   17,
  18,19,20,   22,23,   26
};
static struct tensor symmtensor3b = {
  "symmetric tensor T_i(jk)",
  3, 3, 27, symmtensor3b_vars, symmtensor3b_parity, 18, symmtensor3b_comps
};



/* A fourth rank tensor with symmetries both in its first and last two
   indices */
static int symmtensor4_vars[] = {
   0, 1, 2,    1, 3, 4,    2, 4, 5,
   6, 7, 8,    7, 9,10,    8,10,11,
  12,13,14,   13,15,16,   14,16,17,
  
   6, 7, 8,    7, 9,10,    8,10,11,
  18,19,20,   19,21,22,   20,22,23,
  24,25,26,   25,27,28,   26,28,29,
  
  12,13,14,   13,15,16,   14,16,17,
  24,25,26,   25,27,28,   26,28,29,
  30,31,32,   31,33,34,   32,34,35
};
static int symmtensor4_parity[] = {
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1
};
static int symmtensor4_comps[] = {
   0, 1, 2,    4, 5,    8,
   9,10,11,   13,14,   17,
  18,19,20,   22,23,   26,
  
  36,37,38,   40,41,   44,
  45,46,47,   49,50,   53,
  
  72,73,74,   76,77,   80
};
static struct tensor symmtensor4 = {
  "symmetric tensor T_(ij)(kl)",
  3, 4, 81, symmtensor4_vars, symmtensor4_parity, 36, symmtensor4_comps
};


static int
ipow (int base, int expo)
{
  int res;
  assert (expo >= 0);
  res = expo & 1 ? base : 1;
  while (expo >>= 1)
  {
    base *= base;
    if (expo & 1) res *= base;
  }
  return res;
}



static int
ilog (int res, int const base)
{
  int expo;
  assert (base > 0);
  assert (res >= 0);
  for (expo = 0; res > 0; ++ expo) {
    res /= base;
  }
  return expo;
}



/* Ensure that all tensor declarations are internally consistent */
static void
CheckTensorType (struct tensor const * restrict const atensor)
{
  int i, n;
  assert (atensor->name);
  assert (atensor->dim>=0);
  assert (atensor->rank>=0);
  assert (atensor->ncomps>=0);
  assert (atensor->ncomps == floor(pow(atensor->dim, atensor->rank) + 0.5));
  assert (atensor->nvars>=0 && atensor->nvars<=atensor->ncomps);
  assert (atensor->vars);
  for (i=0; i<atensor->ncomps; ++i) {
    assert (atensor->vars[i]>=0 && atensor->vars[i]<atensor->nvars);
    assert (abs(atensor->parity[i]) <= 1);
  }
  assert (atensor->comps);
  for (n=0; n<atensor->nvars; ++n) {
    assert (atensor->comps[n]>=0 && atensor->comps[n]<atensor->ncomps);
    assert (atensor->vars[atensor->comps[n]] == n);
  }
}

void
Cartoon2D_CheckTensorTypes (CCTK_ARGUMENTS)
{
  int gi;
  
  /* Check all internal tensor type definitions */
  CheckTensorType (&scalar);
  CheckTensorType (&vector);
  CheckTensorType (&tensor);
  CheckTensorType (&symmtensor);
  CheckTensorType (&symmtensor3a);
  CheckTensorType (&symmtensor3b);
  CheckTensorType (&symmtensor4);
  
  /* Check tensor types of all groups */
  for (gi=0; gi<CCTK_NumGroups(); ++gi) {
    
    char tensortypealias[100];
    int numvars, firstvar;
    int table;
    int ierr;
    
    numvars = CCTK_NumVarsInGroupI(gi);
    if (numvars == 0) continue;
    assert (numvars>0);
    firstvar = CCTK_FirstVarIndexI(gi);
    assert (firstvar>=0);
    table = CCTK_GroupTagsTableI(gi);
    assert (table>=0);
    
    ierr = Util_TableGetString
      (table, sizeof tensortypealias, tensortypealias, "tensortypealias");
    if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
      char * groupname = CCTK_GroupName(gi);
      assert (groupname);
      CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Tensor type alias not declared for group \"%s\" -- assuming a scalar",
                  groupname);
      free (groupname);
      strcpy (tensortypealias, "");
    } else if (ierr<0) {
      char * groupname = CCTK_GroupName(gi);
      assert (groupname);
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Error in tensor type alias declaration for group \"%s\"",
                  groupname);
      free (groupname);
    }
    
    if (CCTK_EQUALS (tensortypealias, "")) {
      /* do nothing */
    } else if (CCTK_EQUALS (tensortypealias, "scalar")) {
      /* scalar */
      if (numvars != 1) {
        char * groupname = CCTK_GroupName(gi);
        assert (groupname);
        CCTK_VWarn (3, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Group \"%s\" has the tensor type alias \"scalar\", but contains more than 1 element",
                    groupname);
        free (groupname);
      }
    } else if (CCTK_EQUALS (tensortypealias, "4scalar")) {
      /* 4-scalar */
      if (numvars != 1) {
        char * groupname = CCTK_GroupName(gi);
        assert (groupname);
        CCTK_VWarn (3, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Group \"%s\" has the tensor type alias \"4scalar\", but contains more than 1 element",
                    groupname);
        free (groupname);
      }
    } else if (CCTK_EQUALS (tensortypealias, "u")
               || CCTK_EQUALS (tensortypealias, "d"))
    {
      /* vector */
      assert (numvars == DIM);
    } else if (CCTK_EQUALS (tensortypealias, "4u")
               || CCTK_EQUALS (tensortypealias, "4d"))
    {
      /* 4-vector */
      assert (numvars == DIM+1);
    } else if (CCTK_EQUALS (tensortypealias, "uu_sym")
               || CCTK_EQUALS (tensortypealias, "dd_sym")) {
      /* symmetric tensor */
      assert (numvars == DIM*(DIM+1)/2);
    } else if (CCTK_EQUALS (tensortypealias, "4uu_sym")
               || CCTK_EQUALS (tensortypealias, "4dd_sym")) {
      /* symmetric 4-tensor */
      assert (numvars == (DIM+1)*(DIM+2)/2);
    } else {
      char * groupname = CCTK_GroupName(gi);
      assert (groupname);
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Illegal tensor type alias for group \"%s\"",
                  groupname);
      free (groupname);
    }
    
  }
}
  


/* Symmetry interpolation */
CCTK_INT
Cartoon2D_SymmetryInterpolate (CCTK_POINTER_TO_CONST const cctkGH_,
                               CCTK_INT const N_dims,
                               CCTK_INT const local_interp_handle,
                               CCTK_INT const param_table_handle,
                               CCTK_INT const coord_system_handle,
                               CCTK_INT const N_interp_points,
                               CCTK_INT const interp_coords_type,
                               CCTK_POINTER_TO_CONST const interp_coords[],
                               CCTK_INT const N_input_arrays,
                               CCTK_INT const input_array_indices[],
                               CCTK_INT const N_output_arrays,
                               CCTK_INT const output_array_types[],
                               CCTK_POINTER const output_arrays[],
                               CCTK_INT const faces)
{
  cGH const * restrict const cctkGH = cctkGH_;
  
  CCTK_POINTER new_interp_coords[DIM];
  CCTK_POINTER * restrict new_output_arrays;
  CCTK_INT new_faces;
  CCTK_INT * restrict input_array_time_levels;
  CCTK_INT * restrict operand_indices;
  CCTK_INT * restrict operation_codes;
  CCTK_INT * restrict time_deriv_order;
  CCTK_INT * restrict output_array_indices;
  
  struct tensor const * restrict * restrict thetensor
    [MAX_TIME_DERIV_ORDER+1][MAX_SPATIAL_DERIV_ORDER+1][MAX_TIME_LEVELS+1];
  int * restrict * restrict thevars
    [MAX_TIME_DERIV_ORDER+1][MAX_SPATIAL_DERIV_ORDER+1][MAX_TIME_LEVELS+1];
  int * restrict thebase;
  int * restrict thevar;
  
  int m;                        /* output array */
  int n;                        /* point */
  int i;                        /* var */
  int f;                        /* face */
  int d;                        /* dimension */
  int p;                        /* number of time derivatives */
  int q;                        /* number of spatial derivatives */
  int tl;                       /* time level */
  int r;                        /* rank */
  
  int iret;                     /* interpolator return value */
  int ierr;
  
  /* Check arguments */
  assert (N_dims == DIM);
  assert (N_interp_points >= 0);
  assert (interp_coords_type >= 0);
  for (d=0; d<DIM; ++d) {
    assert (N_interp_points == 0 || interp_coords[d]);
  }
  assert (N_output_arrays >= 0);
  
  /* Coordinates must be CCTK_REAL */
  assert (interp_coords_type == CCTK_VARIABLE_REAL);
  for (m=0; m<N_output_arrays; ++m) {
    assert (output_array_types[m] == CCTK_VARIABLE_REAL);
  }
  
  
  
  /* Claim faces */
  assert (faces & (1 << 0));
  assert (faces & (1 << 2));
  assert (faces & (1 << 3));
  new_faces = faces;
  new_faces &= ~ (1 << 0);
  new_faces &= ~ (1 << 2);
  new_faces &= ~ (1 << 3);
  
  
  
  /* Copy coordinates */
  for (d=0; d<DIM; ++d) {
    new_interp_coords[d] = malloc (N_interp_points * sizeof(CCTK_REAL));
    assert (new_interp_coords[d]);
  }
  
  /* Fold coordinates */
  for (n=0; n<N_interp_points; ++n) {
    /* Rotate the point onto the positive x axis */
    CCTK_REAL x = ((CCTK_REAL const *)interp_coords[0])[n];
    CCTK_REAL y = ((CCTK_REAL const *)interp_coords[1])[n];
    CCTK_REAL z = ((CCTK_REAL const *)interp_coords[2])[n];
    CCTK_REAL const r = sqrt(x*x + y*y);
    /* CCTK_REAL const w = atan2(y, x); */
    ((CCTK_REAL *)new_interp_coords[0])[n] = r;
    ((CCTK_REAL *)new_interp_coords[1])[n] = 0;
    ((CCTK_REAL *)new_interp_coords[2])[n] = z;
  }
  
  
  
  /* Allocate new output arrays */
  new_output_arrays = malloc (N_output_arrays * sizeof *new_output_arrays);
  assert (new_output_arrays);
  for (m=0; m<N_output_arrays; ++m) {
    new_output_arrays[m] = malloc (N_interp_points * sizeof(CCTK_REAL));
    assert (new_output_arrays[m]);
  }
  
  
  
  /* Recursive call */
  iret = SymmetryInterpolateFaces
    (cctkGH_, N_dims,
     local_interp_handle, param_table_handle, coord_system_handle,
     N_interp_points, interp_coords_type, new_interp_coords,
     N_input_arrays, input_array_indices,
     N_output_arrays, output_array_types, new_output_arrays,
     new_faces);
  
  
  
  /* Free coordinates */
  for (d=0; d<DIM; ++d) {
    free (new_interp_coords[d]);
  }
  
  
  
  /* Find output variable indices */
  input_array_time_levels
    = malloc (N_input_arrays * sizeof *input_array_time_levels);
  assert (input_array_time_levels);
  ierr = Util_TableGetIntArray
    (param_table_handle, N_input_arrays,
     input_array_time_levels, "input_array_time_levels");
  if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    for (m=0; m<N_input_arrays; ++m) {
      input_array_time_levels[m] = 0; /* time level is 0 */
    }
  } else {
    assert (ierr == N_input_arrays);
  }
  
  operand_indices = malloc (N_output_arrays * sizeof *operand_indices);
  assert (operand_indices);
  ierr = Util_TableGetIntArray
    (param_table_handle, N_output_arrays, operand_indices, "operand_indices");
  if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    assert (N_output_arrays == N_input_arrays);
    for (m=0; m<N_output_arrays; ++m) {
      operand_indices[m] = m;   /* output index equals input index */
    }
  } else {
    assert (ierr == N_output_arrays);
  }
  
  operation_codes = malloc (N_output_arrays * sizeof *operation_codes);
  assert (operation_codes);
  ierr = Util_TableGetIntArray
    (param_table_handle, N_output_arrays, operation_codes, "operation_codes");
  if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    assert (N_output_arrays == N_input_arrays); /* why? */
    for (m=0; m<N_output_arrays; ++m) {
      operation_codes[m] = 0;     /* do not take derivatives */
    }
  } else {
    assert (ierr == N_output_arrays);
  }
    for (m=0; m<N_output_arrays; ++m) {
    assert (operation_codes[m] >= 0
            && operation_codes[m] < ipow (10, MAX_SPATIAL_DERIV_ORDER));
  }
  
  time_deriv_order = malloc (N_output_arrays * sizeof *time_deriv_order);
  assert (time_deriv_order);
  ierr = Util_TableGetIntArray
    (param_table_handle, N_output_arrays, time_deriv_order, "time_deriv_order");
  if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    for (m=0; m<N_output_arrays; ++m) {
      time_deriv_order[m] = 0;  /* do not take time derivatives */
    }
  } else {
    assert (ierr == N_output_arrays);
  }
  for (m=0; m<N_output_arrays; ++m) {
    assert (time_deriv_order[m] >= 0
            && time_deriv_order[m] <= MAX_TIME_DERIV_ORDER);
  }

  output_array_indices
    = malloc (N_output_arrays * sizeof *output_array_indices);
  assert (output_array_indices);
  for (m=0; m<N_output_arrays; ++m) {
    assert (operand_indices[m]>=0 && operand_indices[m]<N_input_arrays);
    output_array_indices[m] = input_array_indices[operand_indices[m]];
    assert (output_array_indices[m]>=0
            && output_array_indices[m]<CCTK_NumVars());
  }
  
  
  
  /* Map Cactus variables to tensor objects */
  for (p=0; p<=MAX_TIME_DERIV_ORDER; ++p) {
    for (q=0; q<=MAX_SPATIAL_DERIV_ORDER; ++q) {
      for (tl=0; tl<=MAX_TIME_LEVELS; ++tl) {
        thetensor[p][q][tl]
          = malloc (CCTK_NumVars() * sizeof *thetensor[p][q][tl]);
        assert (thetensor[p][q][tl]);
        thevars[p][q][tl] = malloc (CCTK_NumVars() * sizeof *thevars[p][q][tl]);
        assert (thevars[p][q][tl]);
        for (n=0; n<CCTK_NumVars(); ++n) {
          thetensor[p][q][tl][n] = NULL;
          thevars[p][q][tl][n] = NULL;
        }
      }
    }
  }
  
  /* Map output arrays to the base Cactus variable (i.e. the first
     component in the tensor objects) */
  thebase = malloc (N_output_arrays * sizeof *thebase);
  assert (thebase);
  thevar = malloc (N_output_arrays * sizeof *thevar);
  assert (thevar);
  
  for (m=0; m<N_output_arrays; ++m) {
    
    int vi, gi;
    int numvars, firstvar;
    int table;
    
    char tensortypealias[100];
    
    struct tensor const * tensortype;
    int basevar;
    int var;
    int indices[MAX_RANK+1];
    int oldrank;
    int num_time_derivs;
    int num_derivs;
    int time_level;
    
    /* Get some variable information */
    vi = output_array_indices[m];
    assert (vi>=0 && vi<CCTK_NumVars());
    gi = CCTK_GroupIndexFromVarI (vi);
    assert (gi>=0 && gi<CCTK_NumGroups());
    numvars = CCTK_NumVarsInGroupI(gi);
    assert (numvars>0);
    firstvar = CCTK_FirstVarIndexI(gi);
    assert (firstvar>=0);
    table = CCTK_GroupTagsTableI(gi);
    assert (table>=0);
    
    /* Get the tensor type alias */
    ierr = Util_TableGetString
      (table, sizeof tensortypealias, tensortypealias, "tensortypealias");
    if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
      char * groupname = CCTK_GroupName(gi);
      assert (groupname);
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Tensor type alias not declared for group \"%s\" -- assuming a scalar",
                  groupname);
      free (groupname);
      strcpy (tensortypealias, "scalar");
    } else if (ierr<0) {
      char * groupname = CCTK_GroupName(gi);
      assert (groupname);
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Error in tensor type alias declaration for group \"%s\"",
                  groupname);
      free (groupname);
    }
    
    /* Find the tensor type */
    tensortype = NULL;
    basevar = -1;
    var = -1;
    if (CCTK_EQUALS (tensortypealias, "scalar")) {
      /* scalar */
      if (numvars != 1) {
        char * groupname = CCTK_GroupName(gi);
        assert (groupname);
        CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Group \"%s\" has the tensor type alias \"scalar\", but contains more than 1 element",
                    groupname);
        free (groupname);
      }
      tensortype = &scalar;
      basevar = vi;
      var = 0;
    } else if ((CCTK_EQUALS (tensortypealias, "u")) ||(CCTK_EQUALS (tensortypealias, "d"))) {
      /* vector */
      assert (numvars == DIM);
      tensortype = &vector;
      basevar = firstvar;
      var = vi - firstvar;
    } else if (CCTK_EQUALS (tensortypealias, "4u")
               || CCTK_EQUALS (tensortypealias, "4d"))
    {
      /* 4-vector */
      assert (numvars == DIM+1);
      if (vi == firstvar) {
        /* temporal component */
        int const off = 0;
        tensortype = &scalar;
        basevar = firstvar + off;
        var = vi - basevar;
      } else {
        /* spatial components */
        int const off = 1;
        tensortype = &vector;
        basevar = firstvar + off;
        var = vi - basevar;
      }
    } else if (CCTK_EQUALS (tensortypealias, "dd_sym")) {
      /* symmetric tensor */
      assert (numvars == DIM*(DIM+1)/2);
      tensortype = &symmtensor;
      basevar = firstvar;
      var = vi - firstvar;
    } else if (CCTK_EQUALS (tensortypealias, "4uu_sym")
               || CCTK_EQUALS (tensortypealias, "4dd_sym")) {
      /* symmetric 4-tensor */
      assert (numvars == (DIM+1)*(DIM+2)/2);
      if (vi == firstvar) {
        /* temporal-temporal component */
        int const off = 0;
        tensortype = &scalar;
        basevar = firstvar + off;
        var = vi - basevar;
      } else if (vi < firstvar+DIM+1) {
        /* temporal-spatial components */
        int const off = 1;
        tensortype = &vector;
        basevar = firstvar + off;
        var = vi - basevar;
      } else {
        /* spatial-spatial components */
        int const off = DIM+1;
        tensortype = &symmtensor;
        basevar = firstvar + off;
        var = vi - basevar;
      }
    } else {
      char * groupname = CCTK_GroupName(gi);
      assert (groupname);
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Illegal tensor type alias for group \"%s\"",
                  groupname);
      free (groupname);
    }
    thebase[m] = basevar;
    
    /* Find my component */
    {
      int component;
      assert (var>=0 && var<tensortype->nvars);
      component = tensortype->comps[var];
      assert (tensortype->rank>=0 && tensortype->rank<=MAX_RANK);
      for (r=tensortype->rank-1; r>=0; --r) {
        indices[r] = component % tensortype->dim;
        assert (indices[r]>=0 && indices[r]<tensortype->dim);
        component /= tensortype->dim;
      }
      assert (component == 0);
      oldrank = tensortype->rank;
    }
        
    /* Take time derivative order (i.e., time derivatives) into account */
    num_time_derivs = time_deriv_order[m];
    assert (num_time_derivs>=0 && num_time_derivs<=MAX_TIME_DERIV_ORDER);
    
    /* Take operation code (i.e., spatial derivatives) into account */
    num_derivs = ilog (operation_codes[m], 10);
    assert (num_derivs>=0 && num_derivs<=MAX_SPATIAL_DERIV_ORDER);
    
    /* Get the time level */
    time_level = input_array_time_levels[operand_indices[m]];
    assert (time_level>=0 && time_level<=MAX_TIME_LEVELS);
    
    if (tensortype == &scalar) {
      switch (num_derivs) {
      case 0: break;
      case 1: tensortype = &vector; break;
      case 2: tensortype = &symmtensor; break;
      default: assert (0);
      }
    } else if (tensortype == &vector) {
      switch (num_derivs) {
      case 0: break;
      case 1: tensortype = &tensor; break;
      case 2: tensortype = &symmtensor3b; break;
      default: assert (0);
      }
    } else if (tensortype == &symmtensor) {
      switch (num_derivs) {
      case 0: break;
      case 1: tensortype = &symmtensor3a; break;
      case 2: tensortype = &symmtensor4; break;
      default: assert (0);
      }
    } else {
      assert (0);
    }
    
    /* Find the additional indices */
    {
      int code;
      assert (tensortype->rank>=0 && tensortype->rank<=MAX_RANK);
      assert (tensortype->rank >= oldrank);
      code = operation_codes[m];
      for (r=oldrank; r<tensortype->rank; ++r) {
        const int thedir = code % 10 - 1;
        code /= 10;
        assert (thedir>=0 && thedir<tensortype->dim);
        indices[r] = thedir;
      }
    }
    
    /* Re-calculate component */
    {
      int component = 0;
      for (r=0; r<tensortype->rank; ++r) {
        component = component * tensortype->dim + indices[r];
      }
      assert (component>=0 && component<tensortype->ncomps);
      var = tensortype->vars[component];
      assert (var>=0 && var<tensortype->nvars);
      thevar[m] = var;
    }
    
    /* Create or cross-check the tensor object */
    if (! thetensor[num_time_derivs][num_derivs][time_level][basevar]) {
      thetensor[num_time_derivs][num_derivs][time_level][basevar] = tensortype;
      assert (! thevars[num_time_derivs][num_derivs][time_level][basevar]);
      thevars[num_time_derivs][num_derivs][time_level][basevar]
        = malloc (tensortype->nvars
                  * sizeof *thevars[num_time_derivs][num_derivs][time_level][basevar]);
      assert (thevars[num_time_derivs][num_derivs][time_level][basevar]);
      for (i=0; i<tensortype->nvars; ++i) {
        thevars[num_time_derivs][num_derivs][time_level][basevar][i] = -1;
      }
    }
    assert (thetensor[num_time_derivs][num_derivs][time_level][basevar]
            == tensortype);
    assert (thevars[num_time_derivs][num_derivs][time_level][basevar]);
    /* This does not hold if the caller requests the same
       interpolation to be done into different output arrays.  This
       may happen e.g. when CarpetInterp needs to differentiate in
       time.  This is arguably a performance bug in CarpetInterp.  */
    /* See whether this goes away now.  */
    assert (thevars[num_time_derivs][num_derivs][time_level][basevar][var]
            == -1);
    thevars[num_time_derivs][num_derivs][time_level][basevar][var] = m;
    
  } /* for m */
  
  
  
  /* Loop over all output arrays */
  for (m=0; m<N_output_arrays; ++m) {
    
    int num_time_derivs;
    int num_derivs;
    int time_level;
    int basevar;
    struct tensor const * restrict tensortype;
    int const * restrict vars;
    int var;
    int indices[MAX_RANK+1];
    
    /* Take time derivative order (i.e., time derivatives) into account */
    num_time_derivs = time_deriv_order[m];
    assert (num_time_derivs>=0 && num_time_derivs<=MAX_TIME_DERIV_ORDER);
    
    /* Take operation code (i.e., spatial derivatives) into account */
    num_derivs = ilog (operation_codes[m], 10);
    assert (num_derivs>=0 && num_derivs<=MAX_SPATIAL_DERIV_ORDER);
    
    /* Get the time level */
    time_level = input_array_time_levels[operand_indices[m]];
    assert (time_level>=0 && time_level<=MAX_TIME_LEVELS);
        
    /* Get the tensor type */
    basevar = thebase[m];
    assert (basevar>=0 && basevar<=CCTK_NumVars());
    tensortype = thetensor[num_time_derivs][num_derivs][time_level][basevar];
    assert (tensortype);
    vars = thevars[num_time_derivs][num_derivs][time_level][basevar];
    assert (vars);
    var = thevar[m];
    assert (var>=0 && var<tensortype->nvars);
    
    /* Transform into indices */
    {
      int component = tensortype->comps[var];
      assert (tensortype->rank>=0 && tensortype->rank<=MAX_RANK);
      for (r=tensortype->rank-1; r>=0; --r) {
        indices[r] = component % tensortype->dim;
        assert (indices[r]>=0 && indices[r]<tensortype->dim);
        component /= tensortype->dim;
      }
    }
    
    /* Loop over all grid points */
    for (n=0; n<N_interp_points; ++n) {
      
      CCTK_REAL const x = ((CCTK_REAL const *)interp_coords[0])[n];
      CCTK_REAL const y = ((CCTK_REAL const *)interp_coords[1])[n];
/*       CCTK_REAL const z = ((CCTK_REAL const *)interp_coords[2])[n]; */
      
      CCTK_REAL const r = sqrt(x*x + y*y);
/*       CCTK_REAL const w = atan2(y, x); */
      
/*       CCTK_REAL const cw = cos(w); */
/*       CCTK_REAL const sw = sin(w); */
      CCTK_REAL const cw = (x+1.0e-10)/(r+1.0e-10);
      CCTK_REAL const sw = y/(r+1.0e-10);
      
      CCTK_REAL rot[3][3];
      rot[0][0] =  cw;
      rot[0][1] = -sw;
      rot[0][2] =  0;
      rot[1][0] =  sw;
      rot[1][1] =  cw;
      rot[1][2] =  0;
      rot[2][0] =  0;
      rot[2][1] =  0;
      rot[2][2] =  1;
      
      CCTK_REAL myval = 0;
      
      int thecomponent;
      for (thecomponent=0; thecomponent<tensortype->ncomps; ++thecomponent) {
        int theindices[MAX_RANK+1];
        int thevar;
        int theparity;
        int mm;
        assert (tensortype->rank<MAX_RANK+1);
        {
          int r;
          int component = thecomponent;
          for (r=tensortype->rank-1; r>=0; --r) {
            theindices[r] = component % tensortype->dim;
            assert (theindices[r]>=0 && theindices[r]<tensortype->dim);
            component /= tensortype->dim;
          }
          assert (component == 0);
        }
        theparity = tensortype->parity[thecomponent];
        assert (abs(theparity) <= 1);
        thevar = tensortype->vars[thecomponent];
        assert (thevar>=0 && thevar<tensortype->nvars);
        mm = vars[thevar];
        assert (mm>=0 && mm<N_output_arrays);
        {
          int r;
          CCTK_REAL val
            = theparity * ((CCTK_REAL const *)new_output_arrays[mm])[n];
          for (r=0; r<tensortype->rank; ++r) {
            val *= rot[indices[r]][theindices[r]];
          }
          myval += val;
        }
      }
      
      ((CCTK_REAL *)output_arrays[m])[n] = myval;
      
    } /* for n */
    
  } /* for m */
  
  
  
  /* Free output variable descriptions */
  for (m=0; m<N_output_arrays; ++m) {
    free (new_output_arrays[m]);
  }
  free (new_output_arrays);
  
  free (input_array_time_levels);
  free (operand_indices);
  free (operation_codes);
  free (output_array_indices);
  
  
  for (p=0; p<=MAX_TIME_DERIV_ORDER; ++p) {
    for (q=0; q<=MAX_SPATIAL_DERIV_ORDER; ++q) {
      for (tl=0; tl<=MAX_TIME_LEVELS; ++tl) {
        free ((void *) (thetensor[p][q][tl]));
        free ((void *) (thevars[p][q][tl]));
      }
    }
  }
  free (thebase);
  free (thevar);
  
  

  return iret;
}
