/*@@
   @file    excision_Parameters.c
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

/* structure containing all private parameters of thorn excision */
struct
{
  int dummy_parameter;
} PRIVATE_EXCISION_STRUCT;


/* structure containing all restricted parameters of thorn excision */
struct
{
  CCTK_REAL C_ko;
} RESTRICTED_EXCISION_STRUCT;


int CCTKi_BindingsCreateexcisionParameters(void);
int CCTKi_BindingsCreateexcisionParameters(void)
{
  CCTKi_ParameterCreate("excision_enable",
                        "excision",
                        "INT",
                        "GLOBAL",
                        CCTK_STEERABLE_ALWAYS,
                        "disable (0) or enable (1) excision",
                        "0",
                        &(GLOBAL_PARAMETER_STRUCT.excision_enable),
                        0,
                        NULL,
                        1,
                        "*:*", "Any Integer");

  CCTKi_ParameterCreate("excision_radius",
                        "excision",
                        "REAL",
                        "GLOBAL",
                        CCTK_STEERABLE_ALWAYS,
                        "Excision radius.  Set to some fraction of AH radius",
                        "0.1",
                        &(GLOBAL_PARAMETER_STRUCT.excision_radius),
                        0,
                        NULL,
                        1,
                        "0:*", "Any Positive Real");

  CCTKi_ParameterCreate("C_ko",
                        "excision",
                        "REAL",
                        "RESTRICTED",
                        CCTK_STEERABLE_NEVER,
                        "Kreiss-Oliger parameter",
                        "0.2",
                        &(RESTRICTED_EXCISION_STRUCT.C_ko),
                        0,
                        NULL,
                        1,
                        "*:*", "Anything");

  return 0;
}

int CCTKi_BindingsexcisionParameterExtensions(void);
int CCTKi_BindingsexcisionParameterExtensions(void)
{
  return 0;
}
