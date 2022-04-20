/*@@
   @file    NaNChecker.c
   @author  Automatically generated by GridFuncStuff.pl
   @desc
            Creates the CCTK variables for thorn NaNChecker
   @enddesc
 @@*/


#define THORN_IS_NaNChecker 1

#include <stddef.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameter.h"
#include "cctk_WarnLevel.h"
#include "cctki_Groups.h"
#include "cctki_FortranWrappers.h"

int CactusBindingsVariables_NaNChecker_Initialise(void);
static int CCTKi_BindingsFortranWrapperNaNChecker(void *_GH, void *fpointer);

static int CCTKi_BindingsFortranWrapperNaNChecker(void *_GH, void *fpointer)
{
  cGH *GH = _GH;
  const int _cctk_zero = 0;
#ifndef CCTK_DEBUG
  CCTK_COMPLEX _cctk_dummy_var[4];
#endif
  void (*function)(NANCHECKER_C2F_PROTO);
  DECLARE_NANCHECKER_C2F
  INITIALISE_NANCHECKER_C2F
  (void) (_cctk_zero + 0);
#ifndef CCTK_DEBUG
  (void) (_cctk_dummy_var + 0);
#endif

  function = (void (*) (NANCHECKER_C2F_PROTO)) fpointer;
  function (PASS_NANCHECKER_C2F (GH));

  return (0);
}

int CactusBindingsVariables_NaNChecker_Initialise(void)
{
  const char * warn_mixeddim_gfs = "";
  int warn_mixeddim = 0;
  const CCTK_INT *allow_mixeddim_gfs;


  allow_mixeddim_gfs = CCTK_ParameterGet ("allow_mixeddim_gfs", "Cactus", 0);

  if (CCTKi_CreateGroup ("NANMASK", "NANCHECKER", "NANCHECKER",
                         "GF", "INT", "PRIVATE",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Prolongation=\"None\" checkpoint=\"no\"",
                         NULL,
                         1,
                         "NaNmask") == 1)
  {
    warn_mixeddim_gfs = "NaNmask";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("NANSFOUND", "NANCHECKER", "NANCHECKER",
                         "SCALAR", "INT", "PRIVATE",
                         0, 1,
                         "NONE", "CONSTANT",
                         "", "",
                         "checkpoint=\"no\"",
                         NULL,
                         1,
                         "NaNsFound") == 1)
  {
    warn_mixeddim_gfs = "NaNsFound";
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

  CCTKi_RegisterFortranWrapper("NaNChecker", CCTKi_BindingsFortranWrapperNaNChecker);

  return 0;
}