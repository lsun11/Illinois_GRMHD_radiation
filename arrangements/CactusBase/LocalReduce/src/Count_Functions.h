 /*@@
   @header    Count_Functions.h
   @date      
   @author    Tom Goodale, Yaakoub Y El Khamra
   @desc
              Prototypes for Count reduction operators
   @enddesc
   @version   $Header: /cactusdevcvs/CactusBase/LocalReduce/src/Count_Functions.h,v 1.5 2005/10/25 20:30:15 yye00 Exp $
 @@*/

#ifndef _Count_FUNCTIONS_H_
#define _Count_FUNCTIONS_H_

#include "cctk.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Count reduction functions */
int LocalReduce_Count_BYTE(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle);

int LocalReduce_Count_INT(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle);

#ifdef HAVE_CCTK_INT1
int LocalReduce_Count_INT1(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle);
#endif                                                              

#ifdef HAVE_CCTK_INT2
int LocalReduce_Count_INT2(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle);
#endif

#ifdef HAVE_CCTK_INT4
int LocalReduce_Count_INT4(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle);
#endif

#ifdef HAVE_CCTK_INT8
int LocalReduce_Count_INT8(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle);
#endif

int LocalReduce_Count_REAL(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle);

#ifdef HAVE_CCTK_REAL4
int LocalReduce_Count_REAL4(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle);
#endif

#ifdef HAVE_CCTK_REAL8
int LocalReduce_Count_REAL8(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle);
#endif

#ifdef HAVE_CCTK_REAL16
int LocalReduce_Count_REAL16(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle);
#endif

int LocalReduce_Count_COMPLEX(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle);

#ifdef HAVE_CCTK_COMPLEX8
int LocalReduce_Count_COMPLEX8(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle);
#endif

#ifdef HAVE_CCTK_COMPLEX16
int LocalReduce_Count_COMPLEX16(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle);
#endif

#ifdef HAVE_CCTK_COMPLEX32
int LocalReduce_Count_COMPLEX32(int i, int weight_on, const void * const weight, CCTK_INT * input_array_offsets, int * indices, int max_iter, int * actual_indices, CCTK_INT * input_array_strides, CCTK_INT * input_array_min_subscripts,const CCTK_INT * input_array_dims, int num_points, int * actual_iters_per_dim, int * iters_per_dim,    int N_dims, const void *const input_arrays[], const CCTK_INT output_number_type_codes[], void * const output_numbers[], int param_table_handle);
#endif

#ifdef __cplusplus
}
#endif

#endif
