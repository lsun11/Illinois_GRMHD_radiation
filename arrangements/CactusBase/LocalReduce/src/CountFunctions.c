 /*@@
   @file      CountFunctions.c
   @date      
   @author    Tom Goodale, Yaakoub Y El Khamra
   @desc
              The functions responsible for performing the actual iteration.
              Having cascaded switch statements broke some compilers.
   @enddesc
   @version   $Id: CountFunctions.c,v 1.5 2005/10/25 20:30:15 yye00 Exp $
 @@*/

#include "cctk.h"
#include "local_reductions.h"
#include "Count_Functions.h"

#ifdef __cplusplus
extern "C" {
#endif

int LocalReduce_Count_BYTE(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle)
{
  int iter = 0;
  int sum_indices = 0;
  int flag, product, j, k, ierr;

  /* Weight variables */
  CCTK_REAL weight_sum = 0.0;
  CCTK_REAL weight_value = 1.0;

#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

  switch (output_number_type_codes[i])
  {
    /* out values type switches*/
    case CCTK_VARIABLE_BYTE:
      ITERATE_ON_ARRAY(i,CCTK_BYTE, input_arrays[i], CCTK_BYTE, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    case CCTK_VARIABLE_INT:
      ITERATE_ON_ARRAY(i,CCTK_BYTE, input_arrays[i], CCTK_INT, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_INT1
    case CCTK_VARIABLE_INT1:
      ITERATE_ON_ARRAY(i,CCTK_BYTE, input_arrays[i], CCTK_INT1, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif                                                              
    #ifdef HAVE_CCTK_INT2
    case CCTK_VARIABLE_INT2:
      ITERATE_ON_ARRAY(i,CCTK_BYTE, input_arrays[i], CCTK_INT2, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT4
    case CCTK_VARIABLE_INT4:
      ITERATE_ON_ARRAY(i,CCTK_BYTE, input_arrays[i], CCTK_INT4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT8
    case CCTK_VARIABLE_INT8:
      ITERATE_ON_ARRAY(i,CCTK_BYTE, input_arrays[i], CCTK_INT8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;                                                                
    #endif
    case CCTK_VARIABLE_REAL:
      ITERATE_ON_ARRAY(i,CCTK_BYTE, input_arrays[i], CCTK_REAL, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
      ITERATE_ON_ARRAY(i,CCTK_BYTE, input_arrays[i], CCTK_REAL4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      ITERATE_ON_ARRAY(i,CCTK_BYTE, input_arrays[i], CCTK_REAL8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
      ITERATE_ON_ARRAY(i,CCTK_BYTE, input_arrays[i], CCTK_REAL16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    case CCTK_VARIABLE_COMPLEX:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif

#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_BYTE, input_arrays[i], CCTK_COMPLEX, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #ifdef HAVE_CCTK_COMPLEX8
    case CCTK_VARIABLE_COMPLEX8:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_BYTE, input_arrays[i], CCTK_COMPLEX8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX16
    case CCTK_VARIABLE_COMPLEX16:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_BYTE, input_arrays[i], CCTK_COMPLEX16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX32
    case CCTK_VARIABLE_COMPLEX32:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_BYTE, input_arrays[i], CCTK_COMPLEX32, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)
    break;
    #endif
  }
  return num_points;
}

int LocalReduce_Count_INT(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle)
{
  int iter = 0;
  int sum_indices = 0;
  int flag, product, j, k, ierr;

  /* Weight variables */
  CCTK_REAL weight_sum = 0.0;
  CCTK_REAL weight_value = 1.0;

#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

  switch (output_number_type_codes[i])
  {
    /* out values type switches*/
    case CCTK_VARIABLE_BYTE:
      ITERATE_ON_ARRAY(i,CCTK_INT, input_arrays[i], CCTK_BYTE, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    case CCTK_VARIABLE_INT:
      ITERATE_ON_ARRAY(i,CCTK_INT, input_arrays[i], CCTK_INT, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_INT1
    case CCTK_VARIABLE_INT1:
      ITERATE_ON_ARRAY(i,CCTK_INT, input_arrays[i], CCTK_INT1, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif                                                              
    #ifdef HAVE_CCTK_INT2
    case CCTK_VARIABLE_INT2:
      ITERATE_ON_ARRAY(i,CCTK_INT, input_arrays[i], CCTK_INT2, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT4
    case CCTK_VARIABLE_INT4:
      ITERATE_ON_ARRAY(i,CCTK_INT, input_arrays[i], CCTK_INT4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT8
    case CCTK_VARIABLE_INT8:
      ITERATE_ON_ARRAY(i,CCTK_INT, input_arrays[i], CCTK_INT8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;                                                                
    #endif
    case CCTK_VARIABLE_REAL:
      ITERATE_ON_ARRAY(i,CCTK_INT, input_arrays[i], CCTK_REAL, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
      ITERATE_ON_ARRAY(i,CCTK_INT, input_arrays[i], CCTK_REAL4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      ITERATE_ON_ARRAY(i,CCTK_INT, input_arrays[i], CCTK_REAL8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
      ITERATE_ON_ARRAY(i,CCTK_INT, input_arrays[i], CCTK_REAL16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    case CCTK_VARIABLE_COMPLEX:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_INT, input_arrays[i], CCTK_COMPLEX, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #ifdef HAVE_CCTK_COMPLEX8
    case CCTK_VARIABLE_COMPLEX8:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_INT, input_arrays[i], CCTK_COMPLEX8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX16
    case CCTK_VARIABLE_COMPLEX16:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_INT, input_arrays[i], CCTK_COMPLEX16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX32
    case CCTK_VARIABLE_COMPLEX32:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_INT, input_arrays[i], CCTK_COMPLEX32, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)
    break;
    #endif
  }
  return num_points;
}

#ifdef HAVE_CCTK_INT1
int LocalReduce_Count_INT1(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle)
{
  int iter = 0;
  int sum_indices = 0;
  int flag, product, j, k, ierr;

  /* Weight variables */
  CCTK_REAL weight_sum = 0.0;
  CCTK_REAL weight_value = 1.0;

#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

  switch (output_number_type_codes[i])
  {
    /* out values type switches*/
    case CCTK_VARIABLE_BYTE:
      ITERATE_ON_ARRAY(i,CCTK_INT1, input_arrays[i], CCTK_BYTE, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    case CCTK_VARIABLE_INT:
      ITERATE_ON_ARRAY(i,CCTK_INT1, input_arrays[i], CCTK_INT, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_INT11
    case CCTK_VARIABLE_INT1:
      ITERATE_ON_ARRAY(i,CCTK_INT1, input_arrays[i], CCTK_INT1, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif                                                              
    #ifdef HAVE_CCTK_INT12
    case CCTK_VARIABLE_INT2:
      ITERATE_ON_ARRAY(i,CCTK_INT1, input_arrays[i], CCTK_INT2, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT14
    case CCTK_VARIABLE_INT4:
      ITERATE_ON_ARRAY(i,CCTK_INT1, input_arrays[i], CCTK_INT4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT18
    case CCTK_VARIABLE_INT8:
      ITERATE_ON_ARRAY(i,CCTK_INT1, input_arrays[i], CCTK_INT8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;                                                                
    #endif
    case CCTK_VARIABLE_REAL:
      ITERATE_ON_ARRAY(i,CCTK_INT1, input_arrays[i], CCTK_REAL, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
      ITERATE_ON_ARRAY(i,CCTK_INT1, input_arrays[i], CCTK_REAL4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      ITERATE_ON_ARRAY(i,CCTK_INT1, input_arrays[i], CCTK_REAL8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
      ITERATE_ON_ARRAY(i,CCTK_INT1, input_arrays[i], CCTK_REAL16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    case CCTK_VARIABLE_COMPLEX:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_INT1, input_arrays[i], CCTK_COMPLEX, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #ifdef HAVE_CCTK_COMPLEX8
    case CCTK_VARIABLE_COMPLEX8:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_INT1, input_arrays[i], CCTK_COMPLEX8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX16
    case CCTK_VARIABLE_COMPLEX16:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_INT1, input_arrays[i], CCTK_COMPLEX16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX32
    case CCTK_VARIABLE_COMPLEX32:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_INT1, input_arrays[i], CCTK_COMPLEX32, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)
    break;
    #endif
  }
  return num_points;
}
#endif                                                              

#ifdef HAVE_CCTK_INT2
int LocalReduce_Count_INT2(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle)
{
  int iter = 0;
  int sum_indices = 0;
  int flag, product, j, k, ierr;

  /* Weight variables */
  CCTK_REAL weight_sum = 0.0;
  CCTK_REAL weight_value = 1.0;

#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

  switch (output_number_type_codes[i])
  {
    /* out values type switches*/
    case CCTK_VARIABLE_BYTE:
      ITERATE_ON_ARRAY(i,CCTK_INT2, input_arrays[i], CCTK_BYTE, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    case CCTK_VARIABLE_INT:
      ITERATE_ON_ARRAY(i,CCTK_INT2, input_arrays[i], CCTK_INT, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_INT21
    case CCTK_VARIABLE_INT1:
      ITERATE_ON_ARRAY(i,CCTK_INT2, input_arrays[i], CCTK_INT1, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif                                                              
    #ifdef HAVE_CCTK_INT22
    case CCTK_VARIABLE_INT2:
      ITERATE_ON_ARRAY(i,CCTK_INT2, input_arrays[i], CCTK_INT2, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT24
    case CCTK_VARIABLE_INT4:
      ITERATE_ON_ARRAY(i,CCTK_INT2, input_arrays[i], CCTK_INT4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT28
    case CCTK_VARIABLE_INT8:
      ITERATE_ON_ARRAY(i,CCTK_INT2, input_arrays[i], CCTK_INT8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;                                                                
    #endif
    case CCTK_VARIABLE_REAL:
      ITERATE_ON_ARRAY(i,CCTK_INT2, input_arrays[i], CCTK_REAL, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
      ITERATE_ON_ARRAY(i,CCTK_INT2, input_arrays[i], CCTK_REAL4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      ITERATE_ON_ARRAY(i,CCTK_INT2, input_arrays[i], CCTK_REAL8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
      ITERATE_ON_ARRAY(i,CCTK_INT2, input_arrays[i], CCTK_REAL16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    case CCTK_VARIABLE_COMPLEX:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_INT2, input_arrays[i], CCTK_COMPLEX, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #ifdef HAVE_CCTK_COMPLEX8
    case CCTK_VARIABLE_COMPLEX8:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_INT2, input_arrays[i], CCTK_COMPLEX8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX16
    case CCTK_VARIABLE_COMPLEX16:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_INT2, input_arrays[i], CCTK_COMPLEX16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX32
    case CCTK_VARIABLE_COMPLEX32:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_INT2, input_arrays[i], CCTK_COMPLEX32, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)
    break;
    #endif
  }
  return num_points;
}
#endif

#ifdef HAVE_CCTK_INT4
int LocalReduce_Count_INT4(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle)
{
  int iter = 0;
  int sum_indices = 0;
  int flag, product, j, k, ierr;

  /* Weight variables */
  CCTK_REAL weight_sum = 0.0;
  CCTK_REAL weight_value = 1.0;

#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

  switch (output_number_type_codes[i])
  {
    /* out values type switches*/
    case CCTK_VARIABLE_BYTE:
      ITERATE_ON_ARRAY(i,CCTK_INT4, input_arrays[i], CCTK_BYTE, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    case CCTK_VARIABLE_INT:
      ITERATE_ON_ARRAY(i,CCTK_INT4, input_arrays[i], CCTK_INT, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_INT1
    case CCTK_VARIABLE_INT1:
      ITERATE_ON_ARRAY(i,CCTK_INT4, input_arrays[i], CCTK_INT1, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif                                                              
    #ifdef HAVE_CCTK_INT2
    case CCTK_VARIABLE_INT2:
      ITERATE_ON_ARRAY(i,CCTK_INT4, input_arrays[i], CCTK_INT2, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT4
    case CCTK_VARIABLE_INT4:
      ITERATE_ON_ARRAY(i,CCTK_INT4, input_arrays[i], CCTK_INT4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT8
    case CCTK_VARIABLE_INT8:
      ITERATE_ON_ARRAY(i,CCTK_INT4, input_arrays[i], CCTK_INT8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;                                                                
    #endif
    case CCTK_VARIABLE_REAL:
      ITERATE_ON_ARRAY(i,CCTK_INT4, input_arrays[i], CCTK_REAL, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
      ITERATE_ON_ARRAY(i,CCTK_INT4, input_arrays[i], CCTK_REAL4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      ITERATE_ON_ARRAY(i,CCTK_INT4, input_arrays[i], CCTK_REAL8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
      ITERATE_ON_ARRAY(i,CCTK_INT4, input_arrays[i], CCTK_REAL16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    case CCTK_VARIABLE_COMPLEX:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_INT4, input_arrays[i], CCTK_COMPLEX, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #ifdef HAVE_CCTK_COMPLEX8
    case CCTK_VARIABLE_COMPLEX8:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_INT4, input_arrays[i], CCTK_COMPLEX8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX16
    case CCTK_VARIABLE_COMPLEX16:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_INT4, input_arrays[i], CCTK_COMPLEX16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX32
    case CCTK_VARIABLE_COMPLEX32:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_INT4, input_arrays[i], CCTK_COMPLEX32, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)
    break;
    #endif
  }
  return num_points;
}
#endif

#ifdef HAVE_CCTK_INT8
int LocalReduce_Count_INT8(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle)
{
  int iter = 0;
  int sum_indices = 0;
  int flag, product, j, k, ierr;

  /* Weight variables */
  CCTK_REAL weight_sum = 0.0;
  CCTK_REAL weight_value = 1.0;

#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

  switch (output_number_type_codes[i])
  {
    /* out values type switches*/
    case CCTK_VARIABLE_BYTE:
      ITERATE_ON_ARRAY(i,CCTK_INT8, input_arrays[i], CCTK_BYTE, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    case CCTK_VARIABLE_INT:
      ITERATE_ON_ARRAY(i,CCTK_INT8, input_arrays[i], CCTK_INT, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_INT1
    case CCTK_VARIABLE_INT1:
      ITERATE_ON_ARRAY(i,CCTK_INT8, input_arrays[i], CCTK_INT1, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif                                                              
    #ifdef HAVE_CCTK_INT2
    case CCTK_VARIABLE_INT2:
      ITERATE_ON_ARRAY(i,CCTK_INT8, input_arrays[i], CCTK_INT2, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT4
    case CCTK_VARIABLE_INT4:
      ITERATE_ON_ARRAY(i,CCTK_INT8, input_arrays[i], CCTK_INT4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT8
    case CCTK_VARIABLE_INT8:
      ITERATE_ON_ARRAY(i,CCTK_INT8, input_arrays[i], CCTK_INT8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;                                                                
    #endif
    case CCTK_VARIABLE_REAL:
      ITERATE_ON_ARRAY(i,CCTK_INT8, input_arrays[i], CCTK_REAL, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
      ITERATE_ON_ARRAY(i,CCTK_INT8, input_arrays[i], CCTK_REAL4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      ITERATE_ON_ARRAY(i,CCTK_INT8, input_arrays[i], CCTK_REAL8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
      ITERATE_ON_ARRAY(i,CCTK_INT8, input_arrays[i], CCTK_REAL16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    case CCTK_VARIABLE_COMPLEX:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_INT8, input_arrays[i], CCTK_COMPLEX, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #ifdef HAVE_CCTK_COMPLEX8
    case CCTK_VARIABLE_COMPLEX8:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_INT8, input_arrays[i], CCTK_COMPLEX8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX16
    case CCTK_VARIABLE_COMPLEX16:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_INT8, input_arrays[i], CCTK_COMPLEX16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX32
    case CCTK_VARIABLE_COMPLEX32:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_INT8, input_arrays[i], CCTK_COMPLEX32, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)
    break;
    #endif
  }
  return num_points;
}
#endif

int LocalReduce_Count_REAL(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle)
{
  int iter = 0;
  int sum_indices = 0;
  int flag, product, j, k, ierr;

  /* Weight variables */
  CCTK_REAL weight_sum = 0.0;
  CCTK_REAL weight_value = 1.0;

#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

  switch (output_number_type_codes[i])
  {
    /* out values type switches*/
    case CCTK_VARIABLE_BYTE:
      ITERATE_ON_ARRAY(i,CCTK_REAL, input_arrays[i], CCTK_BYTE, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    case CCTK_VARIABLE_INT:
      ITERATE_ON_ARRAY(i,CCTK_REAL, input_arrays[i], CCTK_INT, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_INT1
    case CCTK_VARIABLE_INT1:
      ITERATE_ON_ARRAY(i,CCTK_REAL, input_arrays[i], CCTK_INT1, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif                                                              
    #ifdef HAVE_CCTK_INT2
    case CCTK_VARIABLE_INT2:
      ITERATE_ON_ARRAY(i,CCTK_REAL, input_arrays[i], CCTK_INT2, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT4
    case CCTK_VARIABLE_INT4:
      ITERATE_ON_ARRAY(i,CCTK_REAL, input_arrays[i], CCTK_INT4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT8
    case CCTK_VARIABLE_INT8:
      ITERATE_ON_ARRAY(i,CCTK_REAL, input_arrays[i], CCTK_INT8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;                                                                
    #endif
    case CCTK_VARIABLE_REAL:
      ITERATE_ON_ARRAY(i,CCTK_REAL, input_arrays[i], CCTK_REAL, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
      ITERATE_ON_ARRAY(i,CCTK_REAL, input_arrays[i], CCTK_REAL4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      ITERATE_ON_ARRAY(i,CCTK_REAL, input_arrays[i], CCTK_REAL8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
      ITERATE_ON_ARRAY(i,CCTK_REAL, input_arrays[i], CCTK_REAL16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    case CCTK_VARIABLE_COMPLEX:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_REAL, input_arrays[i], CCTK_COMPLEX, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #ifdef HAVE_CCTK_COMPLEX8
    case CCTK_VARIABLE_COMPLEX8:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_REAL, input_arrays[i], CCTK_COMPLEX8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX16
    case CCTK_VARIABLE_COMPLEX16:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_REAL, input_arrays[i], CCTK_COMPLEX16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX32
    case CCTK_VARIABLE_COMPLEX32:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_REAL, input_arrays[i], CCTK_COMPLEX32, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)
    break;
    #endif
  }
  return num_points;
}

#ifdef HAVE_CCTK_REAL4
int LocalReduce_Count_REAL4(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle)
{
  int iter = 0;
  int sum_indices = 0;
  int flag, product, j, k, ierr;

  /* Weight variables */
  CCTK_REAL weight_sum = 0.0;
  CCTK_REAL weight_value = 1.0;

#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

  switch (output_number_type_codes[i])
  {
    /* out values type switches*/
    case CCTK_VARIABLE_BYTE:
      ITERATE_ON_ARRAY(i,CCTK_REAL4, input_arrays[i], CCTK_BYTE, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    case CCTK_VARIABLE_INT:
      ITERATE_ON_ARRAY(i,CCTK_REAL4, input_arrays[i], CCTK_INT, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_INT1
    case CCTK_VARIABLE_INT1:
      ITERATE_ON_ARRAY(i,CCTK_REAL4, input_arrays[i], CCTK_INT1, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif                                                              
    #ifdef HAVE_CCTK_INT2
    case CCTK_VARIABLE_INT2:
      ITERATE_ON_ARRAY(i,CCTK_REAL4, input_arrays[i], CCTK_INT2, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT4
    case CCTK_VARIABLE_INT4:
      ITERATE_ON_ARRAY(i,CCTK_REAL4, input_arrays[i], CCTK_INT4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT8
    case CCTK_VARIABLE_INT8:
      ITERATE_ON_ARRAY(i,CCTK_REAL4, input_arrays[i], CCTK_INT8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;                                                                
    #endif
    case CCTK_VARIABLE_REAL:
      ITERATE_ON_ARRAY(i,CCTK_REAL4, input_arrays[i], CCTK_REAL, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
      ITERATE_ON_ARRAY(i,CCTK_REAL4, input_arrays[i], CCTK_REAL4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      ITERATE_ON_ARRAY(i,CCTK_REAL4, input_arrays[i], CCTK_REAL8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
      ITERATE_ON_ARRAY(i,CCTK_REAL4, input_arrays[i], CCTK_REAL16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    case CCTK_VARIABLE_COMPLEX:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_REAL4, input_arrays[i], CCTK_COMPLEX, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #ifdef HAVE_CCTK_COMPLEX8
    case CCTK_VARIABLE_COMPLEX8:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_REAL4, input_arrays[i], CCTK_COMPLEX8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX16
    case CCTK_VARIABLE_COMPLEX16:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_REAL4, input_arrays[i], CCTK_COMPLEX16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX32
    case CCTK_VARIABLE_COMPLEX32:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_REAL4, input_arrays[i], CCTK_COMPLEX32, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)
    break;
    #endif
  }
  return num_points;
}
#endif

#ifdef HAVE_CCTK_REAL8
int LocalReduce_Count_REAL8(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle)
{
  int iter = 0;
  int sum_indices = 0;
  int flag, product, j, k, ierr;

  /* Weight variables */
  CCTK_REAL weight_sum = 0.0;
  CCTK_REAL weight_value = 1.0;

#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

  switch (output_number_type_codes[i])
  {
    /* out values type switches*/
    case CCTK_VARIABLE_BYTE:
      ITERATE_ON_ARRAY(i,CCTK_REAL8, input_arrays[i], CCTK_BYTE, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    case CCTK_VARIABLE_INT:
      ITERATE_ON_ARRAY(i,CCTK_REAL8, input_arrays[i], CCTK_INT, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_INT1
    case CCTK_VARIABLE_INT1:
      ITERATE_ON_ARRAY(i,CCTK_REAL8, input_arrays[i], CCTK_INT1, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif                                                              
    #ifdef HAVE_CCTK_INT2
    case CCTK_VARIABLE_INT2:
      ITERATE_ON_ARRAY(i,CCTK_REAL8, input_arrays[i], CCTK_INT2, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT4
    case CCTK_VARIABLE_INT4:
      ITERATE_ON_ARRAY(i,CCTK_REAL8, input_arrays[i], CCTK_INT4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT8
    case CCTK_VARIABLE_INT8:
      ITERATE_ON_ARRAY(i,CCTK_REAL8, input_arrays[i], CCTK_INT8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;                                                                
    #endif
    case CCTK_VARIABLE_REAL:
      ITERATE_ON_ARRAY(i,CCTK_REAL8, input_arrays[i], CCTK_REAL, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
      ITERATE_ON_ARRAY(i,CCTK_REAL8, input_arrays[i], CCTK_REAL4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      ITERATE_ON_ARRAY(i,CCTK_REAL8, input_arrays[i], CCTK_REAL8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
      ITERATE_ON_ARRAY(i,CCTK_REAL8, input_arrays[i], CCTK_REAL16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    case CCTK_VARIABLE_COMPLEX:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_REAL8, input_arrays[i], CCTK_COMPLEX, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #ifdef HAVE_CCTK_COMPLEX8
    case CCTK_VARIABLE_COMPLEX8:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_REAL8, input_arrays[i], CCTK_COMPLEX8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX16
    case CCTK_VARIABLE_COMPLEX16:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_REAL8, input_arrays[i], CCTK_COMPLEX16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX32
    case CCTK_VARIABLE_COMPLEX32:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_REAL8, input_arrays[i], CCTK_COMPLEX32, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];

#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)
    break;
    #endif
  }
  return num_points;
}
#endif

#ifdef HAVE_CCTK_REAL16
int LocalReduce_Count_REAL16(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle)
{
  int iter = 0;
  int sum_indices = 0;
  int flag, product, j, k, ierr;

  /* Weight variables */
  CCTK_REAL weight_sum = 0.0;
  CCTK_REAL weight_value = 1.0;

#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

  switch (output_number_type_codes[i])
  {
    /* out values type switches*/
    case CCTK_VARIABLE_BYTE:
      ITERATE_ON_ARRAY(i,CCTK_REAL16, input_arrays[i], CCTK_BYTE, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    case CCTK_VARIABLE_INT:
      ITERATE_ON_ARRAY(i,CCTK_REAL16, input_arrays[i], CCTK_INT, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_INT1
    case CCTK_VARIABLE_INT1:
      ITERATE_ON_ARRAY(i,CCTK_REAL16, input_arrays[i], CCTK_INT1, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif                                                              
    #ifdef HAVE_CCTK_INT2
    case CCTK_VARIABLE_INT2:
      ITERATE_ON_ARRAY(i,CCTK_REAL16, input_arrays[i], CCTK_INT2, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT4
    case CCTK_VARIABLE_INT4:
      ITERATE_ON_ARRAY(i,CCTK_REAL16, input_arrays[i], CCTK_INT4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT8
    case CCTK_VARIABLE_INT8:
      ITERATE_ON_ARRAY(i,CCTK_REAL16, input_arrays[i], CCTK_INT8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;                                                                
    #endif
    case CCTK_VARIABLE_REAL:
      ITERATE_ON_ARRAY(i,CCTK_REAL16, input_arrays[i], CCTK_REAL, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
      ITERATE_ON_ARRAY(i,CCTK_REAL16, input_arrays[i], CCTK_REAL4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      ITERATE_ON_ARRAY(i,CCTK_REAL16, input_arrays[i], CCTK_REAL8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
      ITERATE_ON_ARRAY(i,CCTK_REAL16, input_arrays[i], CCTK_REAL16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    case CCTK_VARIABLE_COMPLEX:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_REAL16, input_arrays[i], CCTK_COMPLEX, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #ifdef HAVE_CCTK_COMPLEX8
    case CCTK_VARIABLE_COMPLEX8:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_REAL16, input_arrays[i], CCTK_COMPLEX8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX16
    case CCTK_VARIABLE_COMPLEX16:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_REAL16, input_arrays[i], CCTK_COMPLEX16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX32
    case CCTK_VARIABLE_COMPLEX32:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL4) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL8) typed_vdata[sum_indices];
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval.Re = (CCTK_REAL16) typed_vdata[sum_indices];
#endif
#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_REAL16, input_arrays[i], CCTK_COMPLEX32, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP
#undef REDUCTION_PREOP_CAST

#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (out_type) typed_vdata[sum_indices];
#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)
    break;
    #endif
  }
  return num_points;
}
#endif

int LocalReduce_Count_COMPLEX(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle)
{
  int iter = 0;
  int sum_indices = 0;
  int flag, product, j, k, ierr;

  /* Weight variables */
  CCTK_REAL weight_sum = 0.0;
  CCTK_REAL weight_value = 1.0;

#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL4) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL8) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL16) typed_vdata[sum_indices].Re;
#endif

#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

  switch (output_number_type_codes[i])
  {
    /* out values type switches*/
    case CCTK_VARIABLE_BYTE:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX, input_arrays[i], CCTK_BYTE, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    case CCTK_VARIABLE_INT:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX, input_arrays[i], CCTK_INT, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_INT1
    case CCTK_VARIABLE_INT1:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX, input_arrays[i], CCTK_INT1, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif                                                              
    #ifdef HAVE_CCTK_INT2
    case CCTK_VARIABLE_INT2:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX, input_arrays[i], CCTK_INT2, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT4
    case CCTK_VARIABLE_INT4:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX, input_arrays[i], CCTK_INT4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT8
    case CCTK_VARIABLE_INT8:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX, input_arrays[i], CCTK_INT8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;                                                                
    #endif
    case CCTK_VARIABLE_REAL:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX, input_arrays[i], CCTK_REAL, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX, input_arrays[i], CCTK_REAL4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX, input_arrays[i], CCTK_REAL8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX, input_arrays[i], CCTK_REAL16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    case CCTK_VARIABLE_COMPLEX:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        (inval).Re = (CCTK_REAL) (typed_vdata[sum_indices]).Re; (inval).Im = (CCTK_REAL) (typed_vdata[sum_indices]).Im;

#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX, input_arrays[i], CCTK_COMPLEX, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL4) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL8) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL16) typed_vdata[sum_indices].Re;
#endif

#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #ifdef HAVE_CCTK_COMPLEX8
    case CCTK_VARIABLE_COMPLEX8:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        (inval).Re = (CCTK_REAL4) (typed_vdata[sum_indices]).Re; (inval).Im = (CCTK_REAL4) (typed_vdata[sum_indices]).Im;

#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX, input_arrays[i], CCTK_COMPLEX8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL4) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL8) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL16) typed_vdata[sum_indices].Re;
#endif

#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX16
    case CCTK_VARIABLE_COMPLEX16:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        (inval).Re = (CCTK_REAL8) (typed_vdata[sum_indices]).Re; (inval).Im = (CCTK_REAL8) (typed_vdata[sum_indices]).Im;

#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_COMPLEX, input_arrays[i], CCTK_COMPLEX16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL4) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL8) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL16) typed_vdata[sum_indices].Re;
#endif

#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX32
    case CCTK_VARIABLE_COMPLEX32:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        (inval).Re = (CCTK_REAL16) (typed_vdata[sum_indices]).Re; (inval).Im = (CCTK_REAL16) (typed_vdata[sum_indices]).Im;

#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_COMPLEX, input_arrays[i], CCTK_COMPLEX32, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL4) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL8) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL16) typed_vdata[sum_indices].Re;
#endif

#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)
    break;
    #endif
  }
  return num_points;
}


#ifdef HAVE_CCTK_COMPLEX8
int LocalReduce_Count_COMPLEX8(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle)
{
  int iter = 0;
  int sum_indices = 0;
  int flag, product, j, k, ierr;

  /* Weight variables */
  CCTK_REAL weight_sum = 0.0;
  CCTK_REAL weight_value = 1.0;

#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL4) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL8) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL16) typed_vdata[sum_indices].Re;
#endif

#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

  switch (output_number_type_codes[i])
  {
    /* out values type switches*/
    case CCTK_VARIABLE_BYTE:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX8, input_arrays[i], CCTK_BYTE, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    case CCTK_VARIABLE_INT:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX8, input_arrays[i], CCTK_INT, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_INT1
    case CCTK_VARIABLE_INT1:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX8, input_arrays[i], CCTK_INT1, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif                                                              
    #ifdef HAVE_CCTK_INT2
    case CCTK_VARIABLE_INT2:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX8, input_arrays[i], CCTK_INT2, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT4
    case CCTK_VARIABLE_INT4:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX8, input_arrays[i], CCTK_INT4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT8
    case CCTK_VARIABLE_INT8:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX8, input_arrays[i], CCTK_INT8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;                                                                
    #endif
    case CCTK_VARIABLE_REAL:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX8, input_arrays[i], CCTK_REAL, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX8, input_arrays[i], CCTK_REAL4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX8, input_arrays[i], CCTK_REAL8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX8, input_arrays[i], CCTK_REAL16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    case CCTK_VARIABLE_COMPLEX:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        (inval).Re = (CCTK_REAL) (typed_vdata[sum_indices]).Re; (inval).Im = (CCTK_REAL) (typed_vdata[sum_indices]).Im;

#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX8, input_arrays[i], CCTK_COMPLEX, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL4) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL8) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL16) typed_vdata[sum_indices].Re;
#endif

#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #ifdef HAVE_CCTK_COMPLEX8
    case CCTK_VARIABLE_COMPLEX8:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        (inval).Re = (CCTK_REAL4) (typed_vdata[sum_indices]).Re; (inval).Im = (CCTK_REAL4) (typed_vdata[sum_indices]).Im;

#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX8, input_arrays[i], CCTK_COMPLEX8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL4) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL8) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL16) typed_vdata[sum_indices].Re;
#endif

#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX16
    case CCTK_VARIABLE_COMPLEX16:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        (inval).Re = (CCTK_REAL8) (typed_vdata[sum_indices]).Re; (inval).Im = (CCTK_REAL8) (typed_vdata[sum_indices]).Im;

#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_COMPLEX8, input_arrays[i], CCTK_COMPLEX16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL4) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL8) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL16) typed_vdata[sum_indices].Re;
#endif

#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX32
    case CCTK_VARIABLE_COMPLEX32:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        (inval).Re = (CCTK_REAL16) (typed_vdata[sum_indices]).Re; (inval).Im = (CCTK_REAL16) (typed_vdata[sum_indices]).Im;

#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_COMPLEX8, input_arrays[i], CCTK_COMPLEX32, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL4) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL8) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL16) typed_vdata[sum_indices].Re;
#endif

#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)
    break;
    #endif
  }
  return num_points;
}
#endif

#ifdef HAVE_CCTK_COMPLEX16
int LocalReduce_Count_COMPLEX16(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle)
{
  int iter = 0;
  int sum_indices = 0;
  int flag, product, j, k, ierr;

  /* Weight variables */
  CCTK_REAL weight_sum = 0.0;
  CCTK_REAL weight_value = 1.0;

#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL4) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL8) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL16) typed_vdata[sum_indices].Re;
#endif

#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

  switch (output_number_type_codes[i])
  {
    /* out values type switches*/
    case CCTK_VARIABLE_BYTE:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX16, input_arrays[i], CCTK_BYTE, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    case CCTK_VARIABLE_INT:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX16, input_arrays[i], CCTK_INT, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_INT1
    case CCTK_VARIABLE_INT1:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX16, input_arrays[i], CCTK_INT1, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif                                                              
    #ifdef HAVE_CCTK_INT2
    case CCTK_VARIABLE_INT2:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX16, input_arrays[i], CCTK_INT2, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT4
    case CCTK_VARIABLE_INT4:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX16, input_arrays[i], CCTK_INT4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT8
    case CCTK_VARIABLE_INT8:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX16, input_arrays[i], CCTK_INT8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;                                                                
    #endif
    case CCTK_VARIABLE_REAL:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX16, input_arrays[i], CCTK_REAL, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX16, input_arrays[i], CCTK_REAL4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX16, input_arrays[i], CCTK_REAL8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX16, input_arrays[i], CCTK_REAL16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    case CCTK_VARIABLE_COMPLEX:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        (inval).Re = (CCTK_REAL) (typed_vdata[sum_indices]).Re; (inval).Im = (CCTK_REAL) (typed_vdata[sum_indices]).Im;

#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX16, input_arrays[i], CCTK_COMPLEX, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL4) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL8) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL16) typed_vdata[sum_indices].Re;
#endif

#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #ifdef HAVE_CCTK_COMPLEX8
    case CCTK_VARIABLE_COMPLEX8:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        (inval).Re = (CCTK_REAL4) (typed_vdata[sum_indices]).Re; (inval).Im = (CCTK_REAL4) (typed_vdata[sum_indices]).Im;

#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX16, input_arrays[i], CCTK_COMPLEX8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL4) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL8) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL16) typed_vdata[sum_indices].Re;
#endif

#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX16
    case CCTK_VARIABLE_COMPLEX16:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        (inval).Re = (CCTK_REAL8) (typed_vdata[sum_indices]).Re; (inval).Im = (CCTK_REAL8) (typed_vdata[sum_indices]).Im;

#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_COMPLEX16, input_arrays[i], CCTK_COMPLEX16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL4) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL8) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL16) typed_vdata[sum_indices].Re;
#endif

#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX32
    case CCTK_VARIABLE_COMPLEX32:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        (inval).Re = (CCTK_REAL16) (typed_vdata[sum_indices]).Re; (inval).Im = (CCTK_REAL16) (typed_vdata[sum_indices]).Im;

#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_COMPLEX16, input_arrays[i], CCTK_COMPLEX32, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL4) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL8) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL16) typed_vdata[sum_indices].Re;
#endif

#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)
    break;
    #endif
  }
  return num_points;
}
#endif

#ifdef HAVE_CCTK_COMPLEX32
int LocalReduce_Count_COMPLEX32(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle)
{
  int iter = 0;
  int sum_indices = 0;
  int flag, product, j, k, ierr;

  /* Weight variables */
  CCTK_REAL weight_sum = 0.0;
  CCTK_REAL weight_value = 1.0;

#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL4) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL8) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL16) typed_vdata[sum_indices].Re;
#endif

#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

  switch (output_number_type_codes[i])
  {
    /* out values type switches*/
    case CCTK_VARIABLE_BYTE:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX32, input_arrays[i], CCTK_BYTE, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    case CCTK_VARIABLE_INT:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX32, input_arrays[i], CCTK_INT, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_INT1
    case CCTK_VARIABLE_INT1:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX32, input_arrays[i], CCTK_INT1, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif                                                              
    #ifdef HAVE_CCTK_INT2
    case CCTK_VARIABLE_INT2:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX32, input_arrays[i], CCTK_INT2, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT4
    case CCTK_VARIABLE_INT4:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX32, input_arrays[i], CCTK_INT4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_INT8
    case CCTK_VARIABLE_INT8:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX32, input_arrays[i], CCTK_INT8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;                                                                
    #endif
    case CCTK_VARIABLE_REAL:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX32, input_arrays[i], CCTK_REAL, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX32, input_arrays[i], CCTK_REAL4, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX32, input_arrays[i], CCTK_REAL8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    #ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX32, input_arrays[i], CCTK_REAL16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
    break;
    #endif
    case CCTK_VARIABLE_COMPLEX:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        (inval).Re = (CCTK_REAL) (typed_vdata[sum_indices]).Re; (inval).Im = (CCTK_REAL) (typed_vdata[sum_indices]).Im;

#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX32, input_arrays[i], CCTK_COMPLEX, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL4) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL8) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL16) typed_vdata[sum_indices].Re;
#endif

#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #ifdef HAVE_CCTK_COMPLEX8
    case CCTK_VARIABLE_COMPLEX8:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        (inval).Re = (CCTK_REAL4) (typed_vdata[sum_indices]).Re; (inval).Im = (CCTK_REAL4) (typed_vdata[sum_indices]).Im;

#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)
      ITERATE_ON_ARRAY(i,CCTK_COMPLEX32, input_arrays[i], CCTK_COMPLEX8, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL4) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL8) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL16) typed_vdata[sum_indices].Re;
#endif

#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX16
    case CCTK_VARIABLE_COMPLEX16:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        (inval).Re = (CCTK_REAL8) (typed_vdata[sum_indices]).Re; (inval).Im = (CCTK_REAL8) (typed_vdata[sum_indices]).Im;

#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_COMPLEX32, input_arrays[i], CCTK_COMPLEX16, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL4) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL8) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL16) typed_vdata[sum_indices].Re;
#endif

#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)

    break;
    #endif
    #ifdef HAVE_CCTK_COMPLEX32
    case CCTK_VARIABLE_COMPLEX32:
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        (inval).Re = (CCTK_REAL16) (typed_vdata[sum_indices]).Re; (inval).Im = (CCTK_REAL16) (typed_vdata[sum_indices]).Im;

#define REDUCTION_OPERATION(Count, scalar) (Count).Re = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) (Count).Re = num_points;
#define REDUCTION_INITIAL(num) (num).Re = 0.0; (num).Im = 0.0;
#define EXTRA_STEP(a, b)

      ITERATE_ON_ARRAY(i,CCTK_COMPLEX32, input_arrays[i], CCTK_COMPLEX32, output_numbers[i], weight_on, weight, input_array_offsets[i], indices, sum_indices, max_iter, iter, flag, actual_indices,input_array_strides, input_array_min_subscripts,input_array_dims,product)
#undef REDUCTION_OPERATION
#undef WEIGHTED_REDUCTION_OPERATION
#undef REDUCTION_INITIAL
#undef EXTRA_STEP

#ifdef  CCTK_REAL_PRECISION_4
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL4) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_8
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL8) typed_vdata[sum_indices].Re;
#elif   CCTK_REAL_PRECISION_16
#undef  REDUCTION_PREOP_CAST
#define REDUCTION_PREOP_CAST(inval, typed_vdata,sum_indices, out_type) \
        inval = (CCTK_REAL16) typed_vdata[sum_indices].Re;
#endif

#define REDUCTION_OPERATION(Count, scalar) Count = num_points;
#define WEIGHTED_REDUCTION_OPERATION(Count, scalar, weight) Count = num_points;
#define REDUCTION_INITIAL(num) num = 0;
#define EXTRA_STEP(a, b)
    break;
    #endif
  }
  return num_points;
}
#endif


#ifdef __cplusplus
}
#endif
