/* $Header: /numrelcvs/AEIThorns/AEILocalInterp/src/Hermite/1d.cube.order4.c,v 1.1.1.1 2003/07/06 11:16:11 jthorn Exp $ */

#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "util_ErrorCodes.h"
#include "cctk.h"
#include "../InterpLocalUniform.h"
#include "../common/structs.h"
#include "../common/load.h"
#include "../common/evaluate.h"
#include "../common/store.h"

/* function prototype */
#define FUNCTION_NAME			AEILocalInterp_U_Herm_1cube_4
#include "../template.h"

#define N_DIMS				1
#define MOLECULE_MIN_M			-2
#define MOLECULE_MAX_M			3
#define MOLECULE_SIZE			6

/* which derivative ops do we support? */
#define HAVE_OP_I
#define HAVE_OP_DX
#define HAVE_OP_DXX

#define XYZ				x
#define FP_XYZ				fp x
#define STRIDE_IJK			stride_i
#define JACOBIAN_MIJK_STRIDE		Jacobian_mi_stride

#define DATA_STRUCT			data_struct_1d_cube_size6
#define COEFFS_STRUCT			coeffs_struct_1d_cube_size6

#define LOAD_DATA_REAL			AEILocalInterp_load_1dcube6_r
#define LOAD_DATA_REAL4			AEILocalInterp_load_1dcube6_r4
#define LOAD_DATA_REAL8			AEILocalInterp_load_1dcube6_r8
#define LOAD_DATA_REAL16		AEILocalInterp_load_1dcube6_r16
#define LOAD_DATA_COMPLEX		AEILocalInterp_load_1dcube6_c
#define LOAD_DATA_COMPLEX8		AEILocalInterp_load_1dcube6_c8
#define LOAD_DATA_COMPLEX16		AEILocalInterp_load_1dcube6_c16
#define LOAD_DATA_COMPLEX32		AEILocalInterp_load_1dcube6_c32

#define EVALUATE_MOLECULE		AEILocalInterp_eval_1dcube6

#define STORE_COEFFS			AEILocalInterp_store_1dcube6

/* note pathnames are all relative to "../template.c" */
#define COEFFS_I_COMPUTE_FILE_NAME	"Hermite/1d.coeffs/1d.cube.order4/coeffs-I.compute.c"
#define COEFFS_DX_COMPUTE_FILE_NAME	"Hermite/1d.coeffs/1d.cube.order4/coeffs-dx.compute.c"
#define COEFFS_DXX_COMPUTE_FILE_NAME	"Hermite/1d.coeffs/1d.cube.order4/coeffs-dxx.compute.c"

/* actual code */
#include "../template.c"
