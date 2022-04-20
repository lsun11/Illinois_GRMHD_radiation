/*@@
   @file    SymBase_Parameters.c
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

/* structure containing all private parameters of thorn SymBase */
struct
{
  CCTK_INT verbose;
} PRIVATE_SYMBASE_STRUCT;


/* structure containing all restricted parameters of thorn SymBase */
struct
{
  int dummy_parameter;
} RESTRICTED_SYMBASE_STRUCT;


int CCTKi_BindingsCreateSymBaseParameters(void);
int CCTKi_BindingsCreateSymBaseParameters(void)
{
  CCTKi_ParameterCreate("verbose",
                        "SymBase",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Output symmetry boundary face descriptions after registration",
                        "yes",
                        &(PRIVATE_SYMBASE_STRUCT.verbose),
                        0,
                        NULL,
                        0);

  return 0;
}

int CCTKi_BindingsSymBaseParameterExtensions(void);
int CCTKi_BindingsSymBaseParameterExtensions(void)
{
  return 0;
}
