/*@@
   @file    CarpetInterp_Parameters.c
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
#include "ParameterCRestrictedCACTUS.h"

/* structure containing all private parameters of thorn CarpetInterp */
struct
{
  CCTK_REAL ipoison;
  CCTK_REAL poison;
  CCTK_INT barriers;
  CCTK_INT check_tree_search;
  CCTK_INT tree_search;
} PRIVATE_CARPETINTERP_STRUCT;


/* structure containing all restricted parameters of thorn CarpetInterp */
struct
{
  int dummy_parameter;
} RESTRICTED_INTERP_STRUCT;


int CCTKi_BindingsCreateCarpetInterpParameters(void);
int CCTKi_BindingsCreateCarpetInterpParameters(void)
{
  CCTKi_ParameterCreate("barriers",
                        "CarpetInterp",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Insert barriers at strategic places for debugging purposes (slows down execution)",
                        "no",
                        &(PRIVATE_CARPETINTERP_STRUCT.barriers),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("check_tree_search",
                        "CarpetInterp",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Cross-check the result of the tree search",
                        "yes",
                        &(PRIVATE_CARPETINTERP_STRUCT.check_tree_search),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("ipoison",
                        "CarpetInterp",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Integer poison value",
                        "-420042",
                        &(PRIVATE_CARPETINTERP_STRUCT.ipoison),
                        0,
                        NULL,
                        1,
                        "*:*", "");

  CCTKi_ParameterCreate("poison",
                        "CarpetInterp",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Poison value",
                        "-4.20042e+30",
                        &(PRIVATE_CARPETINTERP_STRUCT.poison),
                        0,
                        NULL,
                        1,
                        "*:*", "");

  CCTKi_ParameterCreate("tree_search",
                        "CarpetInterp",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Use a tree search to find the source processor",
                        "no",
                        &(PRIVATE_CARPETINTERP_STRUCT.tree_search),
                        0,
                        NULL,
                        0);

  return 0;
}

int CCTKi_BindingsCarpetInterpParameterExtensions(void);
int CCTKi_BindingsCarpetInterpParameterExtensions(void)
{
  return 0;
}

