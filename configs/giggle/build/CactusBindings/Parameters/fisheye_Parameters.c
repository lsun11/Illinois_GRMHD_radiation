/*@@
   @file    fisheye_Parameters.c
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
#include "ParameterCGlobal.h"
#include "ParameterCRestrictedBSSN.h"
#include "ParameterCRestrictedGRID.h"
#include "ParameterCRestrictedIO.h"

/* structure containing all private parameters of thorn fisheye */
struct
{
  int dummy_parameter;
} PRIVATE_FISHEYE_STRUCT;


/* structure containing all restricted parameters of thorn fisheye */
struct
{
  int dummy_parameter;
} RESTRICTED_FISHEYE_STRUCT;


int CCTKi_BindingsCreatefisheyeParameters(void);
int CCTKi_BindingsCreatefisheyeParameters(void)
{
  CCTKi_ParameterCreate("fisheye_enable",
                        "fisheye",
                        "INT",
                        "GLOBAL",
                        CCTK_STEERABLE_NEVER,
                        "Enable Fisheye",
                        "0",
                        &(GLOBAL_PARAMETER_STRUCT.fisheye_enable),
                        0,
                        NULL,
                        1,
                        "0:1", "Either 0 (no) or 1 (yes)");

  return 0;
}

int CCTKi_BindingsfisheyeParameterExtensions(void);
int CCTKi_BindingsfisheyeParameterExtensions(void)
{
  return 0;
}
