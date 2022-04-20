/*@@
   @file    scalarwaveMoL_Parameters.c
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
#include "ParameterCRestrictedGRID.h"
#include "ParameterCRestrictedMETHODOFLINES.h"

/* structure containing all private parameters of thorn scalarwaveMoL */
struct
{
  CCTK_REAL amplitude;
  CCTK_REAL width;
  const char * bound;
  CCTK_INT enable_moving_grid;
  CCTK_INT scalarwave_Symmetry;
} PRIVATE_SCALARWAVEMOL_STRUCT;


/* structure containing all restricted parameters of thorn scalarwaveMoL */
struct
{
  CCTK_INT WaveMoL_MaxNumConstrainedVars;
  CCTK_INT WaveMoL_MaxNumEvolvedVars;
} RESTRICTED_SCALARWAVEMOL_STRUCT;


int CCTKi_BindingsCreatescalarwaveMoLParameters(void);
int CCTKi_BindingsCreatescalarwaveMoLParameters(void)
{
  CCTKi_ParameterCreate("WaveMoL_MaxNumConstrainedVars",
                        "scalarwaveMoL",
                        "INT",
                        "RESTRICTED",
                        CCTK_STEERABLE_NEVER,
                        "The maximum number of constrained variables used by WaveMoL",
                        "1",
                        &(RESTRICTED_SCALARWAVEMOL_STRUCT.WaveMoL_MaxNumConstrainedVars),
                        0,
                        NULL,
                        1,
                        "1:1", "The Analytic-Numerical gridfunction");

  CCTKi_ParameterCreate("WaveMoL_MaxNumEvolvedVars",
                        "scalarwaveMoL",
                        "INT",
                        "RESTRICTED",
                        CCTK_STEERABLE_NEVER,
                        "The maximum number of evolved variables used by WaveMoL",
                        "5",
                        &(RESTRICTED_SCALARWAVEMOL_STRUCT.WaveMoL_MaxNumEvolvedVars),
                        0,
                        NULL,
                        1,
                        "5:5", "Just 5: phi and the four derivatives");

  CCTKi_ParameterCreate("amplitude",
                        "scalarwaveMoL",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "The amplitude of the waves",
                        "1.0",
                        &(PRIVATE_SCALARWAVEMOL_STRUCT.amplitude),
                        0,
                        NULL,
                        1,
                        "*:*", "No restriction");

  CCTKi_ParameterCreate("bound",
                        "scalarwaveMoL",
                        "KEYWORD",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Type of boundary condition to use",
                        "radiation",
                        &(PRIVATE_SCALARWAVEMOL_STRUCT.bound),
                        0,
                        NULL,
                        3,
                        "none", "No boundary condition",
                        "flat", "Flat boundary condition",
                        "radiation", "Radiation boundary condition");

  CCTKi_ParameterCreate("enable_moving_grid",
                        "scalarwaveMoL",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "yes or no",
                        "0",
                        &(PRIVATE_SCALARWAVEMOL_STRUCT.enable_moving_grid),
                        0,
                        NULL,
                        1,
                        "0:1", "1=yes, 0 =no");

  CCTKi_ParameterCreate("scalarwave_Symmetry",
                        "scalarwaveMoL",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Symmetry choice",
                        "4",
                        &(PRIVATE_SCALARWAVEMOL_STRUCT.scalarwave_Symmetry),
                        0,
                        NULL,
                        1,
                        "0:4", "From 0 to 4");

  CCTKi_ParameterCreate("width",
                        "scalarwaveMoL",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "The width of the wave",
                        "1.0",
                        &(PRIVATE_SCALARWAVEMOL_STRUCT.width),
                        0,
                        NULL,
                        1,
                        "0:*", "Positive");

  return 0;
}

int CCTKi_BindingsscalarwaveMoLParameterExtensions(void);
int CCTKi_BindingsscalarwaveMoLParameterExtensions(void)
{
  CCTKi_ParameterAccumulatorBase("scalarwaveMoL",
                          "WaveMoL_MaxNumConstrainedVars",
                          "MethodofLines",
                          "MoL_Num_Constrained_Vars");

  CCTKi_ParameterAccumulatorBase("scalarwaveMoL",
                          "WaveMoL_MaxNumEvolvedVars",
                          "MethodofLines",
                          "MoL_Num_Evolved_Vars");

  return 0;
}
