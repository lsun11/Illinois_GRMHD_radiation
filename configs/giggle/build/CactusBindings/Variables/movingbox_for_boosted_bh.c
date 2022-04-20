/*@@
   @file    movingbox_for_boosted_bh.c
   @author  Automatically generated by GridFuncStuff.pl
   @desc
            Creates the CCTK variables for thorn movingbox_for_boosted_bh
   @enddesc
 @@*/


#define THORN_IS_movingbox_for_boosted_bh 1

#include <stddef.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameter.h"
#include "cctk_WarnLevel.h"
#include "cctki_Groups.h"
#include "cctki_FortranWrappers.h"

int CactusBindingsVariables_movingbox_for_boosted_bh_Initialise(void);
static int CCTKi_BindingsFortranWrappermovingbox_for_boosted_bh(void *_GH, void *fpointer);

static int CCTKi_BindingsFortranWrappermovingbox_for_boosted_bh(void *_GH, void *fpointer)
{
  cGH *GH = _GH;
  const int _cctk_zero = 0;
#ifndef CCTK_DEBUG
  CCTK_COMPLEX _cctk_dummy_var[4];
#endif
  void (*function)(MOVINGBOX_FOR_BOOSTED_BH_C2F_PROTO);
  DECLARE_MOVINGBOX_FOR_BOOSTED_BH_C2F
  INITIALISE_MOVINGBOX_FOR_BOOSTED_BH_C2F
  (void) (_cctk_zero + 0);
#ifndef CCTK_DEBUG
  (void) (_cctk_dummy_var + 0);
#endif

  function = (void (*) (MOVINGBOX_FOR_BOOSTED_BH_C2F_PROTO)) fpointer;
  function (PASS_MOVINGBOX_FOR_BOOSTED_BH_C2F (GH));

  return (0);
}

int CactusBindingsVariables_movingbox_for_boosted_bh_Initialise(void)
{
  const char * warn_mixeddim_gfs = "";
  int warn_mixeddim = 0;
  const CCTK_INT *allow_mixeddim_gfs;


  allow_mixeddim_gfs = CCTK_ParameterGet ("allow_mixeddim_gfs", "Cactus", 0);

  if (CCTKi_CreateGroup ("TMP_GFS", "MOVINGBOX_FOR_BOOSTED_BH", "MOVINGBOX_FOR_BOOSTED_BH",
                         "GF", "REAL", "PRIVATE",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "",
                         NULL,
                         3,
                         "temp0",
                         "tempx",
                         "tempy") == 1)
  {
    warn_mixeddim_gfs = "tmp_gfs";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("CENTER_ORIGIN", "MOVINGBOX_FOR_BOOSTED_BH", "MOVINGBOX_FOR_BOOSTED_BH",
                         "SCALAR", "REAL", "PRIVATE",
                         0, 1,
                         "NONE", "CONSTANT",
                         "", "",
                         "",
                         NULL,
                         4,
                         "xc_prev",
                         "yc_prev",
                         "zc_prev",
                         "max_phi") == 1)
  {
    warn_mixeddim_gfs = "center_origin";
    warn_mixeddim = 0;
  }
  if (CCTKi_CreateGroup ("MAX_PHI_FLAG", "MOVINGBOX_FOR_BOOSTED_BH", "MOVINGBOX_FOR_BOOSTED_BH",
                         "SCALAR", "INT", "PRIVATE",
                         0, 1,
                         "NONE", "CONSTANT",
                         "", "",
                         "",
                         NULL,
                         1,
                         "find_max_phi") == 1)
  {
    warn_mixeddim_gfs = "max_phi_flag";
    warn_mixeddim = 0;
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

  CCTKi_RegisterFortranWrapper("movingbox_for_boosted_bh", CCTKi_BindingsFortranWrappermovingbox_for_boosted_bh);

  return 0;
}
