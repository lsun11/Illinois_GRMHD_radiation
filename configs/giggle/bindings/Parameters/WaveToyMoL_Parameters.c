/*@@
   @file    WaveToyMoL_Parameters.c
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
#include "ParameterCRestrictedMETHODOFLINES.h"

/* structure containing all private parameters of thorn WaveToyMoL */
struct
{
  const char * bound;
  CCTK_INT estimate_error;
  CCTK_INT order;
} PRIVATE_WAVETOYMOL_STRUCT;


/* structure containing all restricted parameters of thorn WaveToyMoL */
struct
{
  CCTK_INT WaveToyMoL_MaxNumEvolvedVars;
} RESTRICTED_WAVETOYMOL_STRUCT;


int CCTKi_BindingsCreateWaveToyMoLParameters(void);
int CCTKi_BindingsCreateWaveToyMoLParameters(void)
{
  CCTKi_ParameterCreate("WaveToyMoL_MaxNumEvolvedVars",
                        "WaveToyMoL",
                        "INT",
                        "RESTRICTED",
                        CCTK_STEERABLE_NEVER,
                        "The maximum number of evolved variables used by WaveMoL",
                        "2",
                        &(RESTRICTED_WAVETOYMOL_STRUCT.WaveToyMoL_MaxNumEvolvedVars),
                        0,
                        NULL,
                        1,
                        "2:2", "");

  CCTKi_ParameterCreate("bound",
                        "WaveToyMoL",
                        "STRING",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Type of boundary condition to use",
                        "None",
                        &(PRIVATE_WAVETOYMOL_STRUCT.bound),
                        0,
                        NULL,
                        1,
                        ".*", "must be a registered boundary condition");

  CCTKi_ParameterCreate("estimate_error",
                        "WaveToyMoL",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Estimate the truncation error",
                        "no",
                        &(PRIVATE_WAVETOYMOL_STRUCT.estimate_error),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("order",
                        "WaveToyMoL",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Finite differencing order",
                        "2",
                        &(PRIVATE_WAVETOYMOL_STRUCT.order),
                        0,
                        NULL,
                        2,
                        "2", "second order",
                        "4", "fourth order");

  return 0;
}

int CCTKi_BindingsWaveToyMoLParameterExtensions(void);
int CCTKi_BindingsWaveToyMoLParameterExtensions(void)
{
  CCTKi_ParameterAccumulatorBase("WaveToyMoL",
                          "WaveToyMoL_MaxNumEvolvedVars",
                          "MethodofLines",
                          "MoL_Num_Evolved_Vars");

  return 0;
}

