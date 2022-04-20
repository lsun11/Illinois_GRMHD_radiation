/*@@
   @file    ADMMacros.c
   @author  Automatically generated by GridFuncStuff.pl
   @desc
            Creates the CCTK variables for thorn ADMMacros
   @enddesc
 @@*/


#define THORN_IS_ADMMacros 1

#include <stddef.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameter.h"
#include "cctk_WarnLevel.h"
#include "cctki_Groups.h"
#include "cctki_FortranWrappers.h"

int CactusBindingsVariables_ADMMacros_Initialise(void);
static int CCTKi_BindingsFortranWrapperADMMacros(void *_GH, void *fpointer);

static int CCTKi_BindingsFortranWrapperADMMacros(void *_GH, void *fpointer)
{
  cGH *GH = _GH;
  const int _cctk_zero = 0;
#ifndef CCTK_DEBUG
  CCTK_COMPLEX _cctk_dummy_var[4];
#endif
  void (*function)(ADMMACROS_C2F_PROTO);
  DECLARE_ADMMACROS_C2F
  INITIALISE_ADMMACROS_C2F
  (void) (_cctk_zero + 0);
#ifndef CCTK_DEBUG
  (void) (_cctk_dummy_var + 0);
#endif

  function = (void (*) (ADMMACROS_C2F_PROTO)) fpointer;
  function (PASS_ADMMACROS_C2F (GH));

  return (0);
}

int CactusBindingsVariables_ADMMacros_Initialise(void)
{
  const char * warn_mixeddim_gfs = "";
  int warn_mixeddim = 0;
  const CCTK_INT *allow_mixeddim_gfs;


  allow_mixeddim_gfs = CCTK_ParameterGet ("allow_mixeddim_gfs", "Cactus", 0);

  if (CCTKi_CreateGroup ("LOCAL_SPATIAL_ORDER", "ADMMACROS", "ADMMACROS",
                         "SCALAR", "INT", "PUBLIC",
                         0, 1,
                         "NONE", "CONSTANT",
                         "", "",
                         "",
                         NULL,
                         1,
                         "local_spatial_order") == 1)
  {
    warn_mixeddim_gfs = "local_spatial_order";
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

  CCTKi_RegisterFortranWrapper("ADMMacros", CCTKi_BindingsFortranWrapperADMMacros);

  return 0;
}
