/*@@
   @file    CarpetIOBasic_Parameters.c
   @author  Automatically generated by CreateParameterBindings.pl
   @desc
            Creates/extends parameters for this thorn
   @enddesc
 @@*/


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "cctk_Config.h"
#include "cctk_Constants.h"
#include "ParameterBindings.h"
#include "CParameterStructNames.h"
#include "ParameterCRestrictedIO.h"

/* structure containing all private parameters of thorn CarpetIOBasic */
struct
{
  CCTK_REAL outInfo_dt;
  CCTK_REAL real_max;
  CCTK_REAL real_min;
  const char * outInfo_criterion;
  const char * outInfo_reductions;
  const char * outInfo_vars;
  CCTK_INT int_width;
  CCTK_INT iter_width;
  CCTK_INT outHeader_every;
  CCTK_INT outInfo_every;
  CCTK_INT real_prec;
  CCTK_INT real_prec_sci;
  CCTK_INT real_width;
  CCTK_INT time_prec;
  CCTK_INT time_width;
} PRIVATE_CARPETIOBASIC_STRUCT;


/* structure containing all restricted parameters of thorn CarpetIOBasic */
struct
{
  int dummy_parameter;
} RESTRICTED_IOBASIC_STRUCT;


int CCTKi_BindingsCreateCarpetIOBasicParameters(void);
int CCTKi_BindingsCreateCarpetIOBasicParameters(void)
{
  CCTKi_ParameterCreate("int_width",
                        "CarpetIOBasic",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Field width for integer values",
                        "9",
                        &(PRIVATE_CARPETIOBASIC_STRUCT.int_width),
                        0,
                        NULL,
                        1,
                        "1:*", "");

  CCTKi_ParameterCreate("iter_width",
                        "CarpetIOBasic",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Field width for the current iteration",
                        "9",
                        &(PRIVATE_CARPETIOBASIC_STRUCT.iter_width),
                        0,
                        NULL,
                        1,
                        "1:*", "");

  CCTKi_ParameterCreate("outHeader_every",
                        "CarpetIOBasic",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "How often to print the header",
                        "20",
                        &(PRIVATE_CARPETIOBASIC_STRUCT.outHeader_every),
                        0,
                        NULL,
                        2,
                        "1:*", "Output every so many time steps" ,
                        "-1", "No header output");

  CCTKi_ParameterCreate("outInfo_criterion",
                        "CarpetIOBasic",
                        "KEYWORD",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Criterion to select scalar output intervals, overrides out_every",
                        "iteration",
                        &(PRIVATE_CARPETIOBASIC_STRUCT.outInfo_criterion),
                        0,
                        NULL,
                        3,
                        "never", "Never output",
                        "iteration", "Output every so many iterations",
                        "time", "Output every that much coordinate time");

  CCTKi_ParameterCreate("outInfo_dt",
                        "CarpetIOBasic",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "How often to do scalar output, overrides out_dt",
                        "-2",
                        &(PRIVATE_CARPETIOBASIC_STRUCT.outInfo_dt),
                        0,
                        NULL,
                        4,
                        "(0:*", "In intervals of that much coordinate time",
                        "0", "As often as possible",
                        "-1", "Disable output",
                        "-2", "Default to IO::out_dt");

  CCTKi_ParameterCreate("outInfo_every",
                        "CarpetIOBasic",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "How often to do scalar output, overrides IO::out_every",
                        "-2",
                        &(PRIVATE_CARPETIOBASIC_STRUCT.outInfo_every),
                        0,
                        NULL,
                        3,
                        "1:*", "Output every so many time steps",
                        "-1:0", "No output",
                        "-2", "Default to IO::out_every");

  CCTKi_ParameterCreate("outInfo_reductions",
                        "CarpetIOBasic",
                        "STRING",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "List of reductions to output in scalar form",
                        "minimum maximum",
                        &(PRIVATE_CARPETIOBASIC_STRUCT.outInfo_reductions),
                        0,
                        NULL,
                        1,
                        "", "A regex which matches everything");

  CCTKi_ParameterCreate("outInfo_vars",
                        "CarpetIOBasic",
                        "STRING",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Variables to output in scalar form",
                        "",
                        &(PRIVATE_CARPETIOBASIC_STRUCT.outInfo_vars),
                        0,
                        NULL,
                        1,
                        "", "A regex which matches everything");

  CCTKi_ParameterCreate("real_max",
                        "CarpetIOBasic",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Upper bound for numbers that are displayed in fixed notation",
                        "1.0e+3",
                        &(PRIVATE_CARPETIOBASIC_STRUCT.real_max),
                        0,
                        NULL,
                        1,
                        "(0.0:*", "");

  CCTKi_ParameterCreate("real_min",
                        "CarpetIOBasic",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Lower bound for numbers that are displayed in fixed notation",
                        "1.0e-8",
                        &(PRIVATE_CARPETIOBASIC_STRUCT.real_min),
                        0,
                        NULL,
                        1,
                        "(0.0:*", "");

  CCTKi_ParameterCreate("real_prec",
                        "CarpetIOBasic",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Precision for real values",
                        "7",
                        &(PRIVATE_CARPETIOBASIC_STRUCT.real_prec),
                        0,
                        NULL,
                        1,
                        "0:*", "");

  CCTKi_ParameterCreate("real_prec_sci",
                        "CarpetIOBasic",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Precision for real values in scientific notation",
                        "6",
                        &(PRIVATE_CARPETIOBASIC_STRUCT.real_prec_sci),
                        0,
                        NULL,
                        1,
                        "0:*", "");

  CCTKi_ParameterCreate("real_width",
                        "CarpetIOBasic",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Field width for real values",
                        "12",
                        &(PRIVATE_CARPETIOBASIC_STRUCT.real_width),
                        0,
                        NULL,
                        1,
                        "1:*", "");

  CCTKi_ParameterCreate("time_prec",
                        "CarpetIOBasic",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Precision for the simulation time",
                        "3",
                        &(PRIVATE_CARPETIOBASIC_STRUCT.time_prec),
                        0,
                        NULL,
                        1,
                        "0:*", "");

  CCTKi_ParameterCreate("time_width",
                        "CarpetIOBasic",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Field width for the simulation time",
                        "9",
                        &(PRIVATE_CARPETIOBASIC_STRUCT.time_width),
                        0,
                        NULL,
                        1,
                        "1:*", "");

  return 0;
}

int CCTKi_BindingsCarpetIOBasicParameterExtensions(void);
int CCTKi_BindingsCarpetIOBasicParameterExtensions(void)
{
  return 0;
}
