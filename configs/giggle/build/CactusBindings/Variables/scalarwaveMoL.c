/*@@
   @file    scalarwaveMoL.c
   @author  Automatically generated by GridFuncStuff.pl
   @desc
            Creates the CCTK variables for thorn scalarwaveMoL
   @enddesc
 @@*/


#define THORN_IS_scalarwaveMoL 1

#include <stddef.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameter.h"
#include "cctk_WarnLevel.h"
#include "cctki_Groups.h"
#include "cctki_FortranWrappers.h"

int CactusBindingsVariables_scalarwaveMoL_Initialise(void);
static int CCTKi_BindingsFortranWrapperscalarwaveMoL(void *_GH, void *fpointer);

static int CCTKi_BindingsFortranWrapperscalarwaveMoL(void *_GH, void *fpointer)
{
  cGH *GH = _GH;
  const int _cctk_zero = 0;
#ifndef CCTK_DEBUG
  CCTK_COMPLEX _cctk_dummy_var[4];
#endif
  void (*function)(SCALARWAVEMOL_C2F_PROTO);
  DECLARE_SCALARWAVEMOL_C2F
  INITIALISE_SCALARWAVEMOL_C2F
  (void) (_cctk_zero + 0);
#ifndef CCTK_DEBUG
  (void) (_cctk_dummy_var + 0);
#endif

  function = (void (*) (SCALARWAVEMOL_C2F_PROTO)) fpointer;
  function (PASS_SCALARWAVEMOL_C2F (GH));

  return (0);
}

int CactusBindingsVariables_scalarwaveMoL_Initialise(void)
{
  const char * warn_mixeddim_gfs = "";
  int warn_mixeddim = 0;
  const CCTK_INT *allow_mixeddim_gfs;


  allow_mixeddim_gfs = CCTK_ParameterGet ("allow_mixeddim_gfs", "Cactus", 0);

  if (CCTKi_CreateGroup ("SCALARMOLEVOLVE", "SCALARWAVEMOL", "SCALARWAVEMOL",
                         "GF", "REAL", "PROTECTED",
                         3, 3,
                         "NONE", "DEFAULT",
                         "", "",
                         "tensortypealias=\"Scalar\"",
                         NULL,
                         2,
                         "phi",
                         "phidot") == 1)
  {
    warn_mixeddim_gfs = "scalarMoLevolve";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("SCALARMOLRHS", "SCALARWAVEMOL", "SCALARWAVEMOL",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "",
                         NULL,
                         2,
                         "phi_rhs",
                         "phidot_rhs") == 1)
  {
    warn_mixeddim_gfs = "scalarMoLrhs";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("SCALARMOLSTAGGER", "SCALARWAVEMOL", "SCALARWAVEMOL",
                         "GF", "REAL", "PROTECTED",
                         3, 3,
                         "NONE", "DEFAULT",
                         "", "",
                         "tensortypealias=\"Scalar\" Prolongation=\"STAGGER110\"",
                         NULL,
                         2,
                         "phi_stagger",
                         "phidot_stagger") == 1)
  {
    warn_mixeddim_gfs = "scalarMoLstagger";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("SCALARMOLSTAGGERRHS", "SCALARWAVEMOL", "SCALARWAVEMOL",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "",
                         NULL,
                         2,
                         "phi_stagger_rhs",
                         "phidot_stagger_rhs") == 1)
  {
    warn_mixeddim_gfs = "scalarMoLstaggerrhs";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("SCALARMOLANALYTIC", "SCALARWAVEMOL", "SCALARWAVEMOL",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "",
                         NULL,
                         3,
                         "phi_analytic",
                         "phi_analytic_minus_numeric",
                         "phi_stagger_analytic_minus_numeric") == 1)
  {
    warn_mixeddim_gfs = "scalarMoLanalytic";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("SCAWAVE_HAVE_GLOBAL_BOUNDARY", "SCALARWAVEMOL", "SCALARWAVEMOL",
                         "ARRAY", "INT", "PROTECTED",
                         2, 1,
                         "NONE", "CONSTANT",
                         "16384,3", "",
                         "",
                         NULL,
                         2,
                         "scawave_have_bdry_min",
                         "scawave_have_bdry_max") == 1)
  {
    warn_mixeddim_gfs = "scawave_have_global_boundary";
    warn_mixeddim = 2;
  }

  if (*warn_mixeddim_gfs)
  {
    if (allow_mixeddim_gfs && *allow_mixeddim_gfs)
    {
      CCTK_VWarn (2, __LINE__, __FILE__, "Cactus",
                  "CCTKi_CreateGroup: Working dimension already set, "
                  "creating GF group '%s' with different dimension %d",
                  warn_mixeddim_gfs, warn_mixeddim);
    }
    else
    {
      CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                  "CCTKi_CreateGroup: Working dimension already set, "
                  "cannot create GF group '%s' with dimension %d",
                  warn_mixeddim_gfs, warn_mixeddim);
    }
 }

  CCTKi_RegisterFortranWrapper("scalarwaveMoL", CCTKi_BindingsFortranWrapperscalarwaveMoL);

  return 0;
}
