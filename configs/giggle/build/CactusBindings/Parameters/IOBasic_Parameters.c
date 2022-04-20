/*@@
   @file    IOBasic_Parameters.c
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

/* structure containing all private parameters of thorn IOBasic */
struct
{
  CCTK_REAL outInfo_dt;
  CCTK_REAL outScalar_dt;
  const char * outInfo_criterion;
  const char * outInfo_reductions;
  const char * outInfo_vars;
  const char * outScalar_criterion;
  const char * outScalar_reductions;
  const char * outScalar_style;
  const char * outScalar_vars;
  const char * out_dir;
  const char * out_format;
  CCTK_INT outInfo_every;
  CCTK_INT outScalar_every;
} PRIVATE_IOBASIC_STRUCT;


/* structure containing all restricted parameters of thorn IOBasic */
struct
{
  int dummy_parameter;
} RESTRICTED_IOBASIC_STRUCT;


int CCTKi_BindingsCreateIOBasicParameters(void);
int CCTKi_BindingsCreateIOBasicParameters(void)
{
  CCTKi_ParameterCreate("outInfo_criterion",
                        "IOBasic",
                        "KEYWORD",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Criterion to select Info output intervals",
                        "iteration",
                        &(PRIVATE_IOBASIC_STRUCT.outInfo_criterion),
                        0,
                        NULL,
                        3,
                        "never", "Never output",
                        "iteration", "Output every so many iterations",
                        "time", "Output every that much coordinate time");

  CCTKi_ParameterCreate("outInfo_dt",
                        "IOBasic",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "How often to do Info output",
                        "-2",
                        &(PRIVATE_IOBASIC_STRUCT.outInfo_dt),
                        0,
                        NULL,
                        4,
                        "(0:*", "In intervals of that much coordinate time",
                        "0", "As often as possible",
                        "-1", "Disable output",
                        "-2", "Default to IO::out_dt");

  CCTKi_ParameterCreate("outInfo_every",
                        "IOBasic",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "How often to do Info output",
                        "-1",
                        &(PRIVATE_IOBASIC_STRUCT.outInfo_every),
                        0,
                        NULL,
                        3,
                        "1:*", "Every so many iterations",
                        "0:", "Disable Info output",
                        "-1:", "Default to IO::out_every");

  CCTKi_ParameterCreate("outInfo_reductions",
                        "IOBasic",
                        "STRING",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "List of reductions to output as Info to screen",
                        "minimum maximum",
                        &(PRIVATE_IOBASIC_STRUCT.outInfo_reductions),
                        0,
                        NULL,
                        1,
                        ".+", "Space-separated list of reduction operators");

  CCTKi_ParameterCreate("outInfo_vars",
                        "IOBasic",
                        "STRING",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Variables to output as Info to screen",
                        "",
                        &(PRIVATE_IOBASIC_STRUCT.outInfo_vars),
                        0,
                        NULL,
                        2,
                        ".+", "Space-separated list of fully qualified variable/group names",
                        "^$", "An empty string to output nothing");

  CCTKi_ParameterCreate("outScalar_criterion",
                        "IOBasic",
                        "KEYWORD",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Criterion to select Scalar output intervals",
                        "iteration",
                        &(PRIVATE_IOBASIC_STRUCT.outScalar_criterion),
                        0,
                        NULL,
                        3,
                        "never", "Never output",
                        "iteration", "Output every so many iterations",
                        "time", "Output every that much coordinate time");

  CCTKi_ParameterCreate("outScalar_dt",
                        "IOBasic",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "How often to do Scalar output",
                        "-2",
                        &(PRIVATE_IOBASIC_STRUCT.outScalar_dt),
                        0,
                        NULL,
                        4,
                        "(0:*", "In intervals of that much coordinate time",
                        "0", "As often as possible",
                        "-1", "Disable output",
                        "-2", "Default to IO::out_dt");

  CCTKi_ParameterCreate("outScalar_every",
                        "IOBasic",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "How often to do Scalar output",
                        "-1",
                        &(PRIVATE_IOBASIC_STRUCT.outScalar_every),
                        0,
                        NULL,
                        3,
                        "1:*", "Every so many iterations",
                        "0:", "Disable Scalar output",
                        "-1:", "Default to IO::out_every");

  CCTKi_ParameterCreate("outScalar_reductions",
                        "IOBasic",
                        "STRING",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "List of reductions to output into files",
                        "minimum maximum norm1 norm2",
                        &(PRIVATE_IOBASIC_STRUCT.outScalar_reductions),
                        0,
                        NULL,
                        1,
                        ".+", "Space-separated list of reduction operators");

  CCTKi_ParameterCreate("outScalar_style",
                        "IOBasic",
                        "KEYWORD",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Which style for Scalar output",
                        "xgraph",
                        &(PRIVATE_IOBASIC_STRUCT.outScalar_style),
                        0,
                        NULL,
                        2,
                        "gnuplot", "1D output readable by gnuplot",
                        "xgraph", "1D output readable by xgraph");

  CCTKi_ParameterCreate("outScalar_vars",
                        "IOBasic",
                        "STRING",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Variables to output into files",
                        "",
                        &(PRIVATE_IOBASIC_STRUCT.outScalar_vars),
                        0,
                        NULL,
                        2,
                        ".+", "Space-separated list of fully qualified variable/group names",
                        "^$", "An empty string to output nothing");

  CCTKi_ParameterCreate("out_dir",
                        "IOBasic",
                        "STRING",
                        "PRIVATE",
                        CCTK_STEERABLE_RECOVER,
                        "Output directory for IOBasic's scalar files, overrides IO::out_dir",
                        "",
                        &(PRIVATE_IOBASIC_STRUCT.out_dir),
                        0,
                        NULL,
                        2,
                        ".+", "A valid directory name",
                        "^$", "An empty string to choose the default from IO::out_dir");

  CCTKi_ParameterCreate("out_format",
                        "IOBasic",
                        "STRING",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Which format for Scalar floating-point number output",
                        ".13f",
                        &(PRIVATE_IOBASIC_STRUCT.out_format),
                        0,
                        NULL,
                        1,
                        "^(\\.[1]?[0-9])?[EGefg]$", "output with given precision in exponential / floating point notation");

  return 0;
}

int CCTKi_BindingsIOBasicParameterExtensions(void);
int CCTKi_BindingsIOBasicParameterExtensions(void)
{
  return 0;
}
