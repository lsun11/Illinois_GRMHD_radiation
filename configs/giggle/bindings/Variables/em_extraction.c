/*@@
   @file    em_extraction.c
   @author  Automatically generated by GridFuncStuff.pl
   @desc
            Creates the CCTK variables for thorn em_extraction
   @enddesc
 @@*/


#define THORN_IS_em_extraction 1

#include <stddef.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameter.h"
#include "cctk_WarnLevel.h"
#include "cctki_Groups.h"
#include "cctki_FortranWrappers.h"

int CactusBindingsVariables_em_extraction_Initialise(void);
static int CCTKi_BindingsFortranWrapperem_extraction(void *_GH, void *fpointer);

static int CCTKi_BindingsFortranWrapperem_extraction(void *_GH, void *fpointer)
{
  cGH *GH = _GH;
  const int _cctk_zero = 0;
#ifndef CCTK_DEBUG
  CCTK_COMPLEX _cctk_dummy_var[4];
#endif
  void (*function)(EM_EXTRACTION_C2F_PROTO);
  DECLARE_EM_EXTRACTION_C2F
  INITIALISE_EM_EXTRACTION_C2F
  (void) (_cctk_zero + 0);
#ifndef CCTK_DEBUG
  (void) (_cctk_dummy_var + 0);
#endif

  function = (void (*) (EM_EXTRACTION_C2F_PROTO)) fpointer;
  function (PASS_EM_EXTRACTION_C2F (GH));

  return (0);
}

int CactusBindingsVariables_em_extraction_Initialise(void)
{
  const char * warn_mixeddim_gfs = "";
  int warn_mixeddim = 0;
  const CCTK_INT *allow_mixeddim_gfs;


  allow_mixeddim_gfs = CCTK_ParameterGet ("allow_mixeddim_gfs", "Cactus", 0);

  if (CCTKi_CreateGroup ("NPPHI2_EM", "EM_EXTRACTION", "EM_EXTRACTION",
                         "GF", "REAL", "PROTECTED",
                         3, 3,
                         "NONE", "DEFAULT",
                         "", "",
                         "",
                         NULL,
                         10,
                         "NPphi0re",
                         "NPphi0im",
                         "NPphi2re",
                         "NPphi2im",
                         "SPx",
                         "SPy",
                         "SPz",
                         "T_0x",
                         "T_0y",
                         "T_0z") == 1)
  {
    warn_mixeddim_gfs = "NPphi2_em";
    warn_mixeddim = 3;
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

  CCTKi_RegisterFortranWrapper("em_extraction", CCTKi_BindingsFortranWrapperem_extraction);

  return 0;
}
