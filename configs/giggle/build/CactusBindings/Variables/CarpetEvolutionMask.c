/*@@
   @file    CarpetEvolutionMask.c
   @author  Automatically generated by GridFuncStuff.pl
   @desc
            Creates the CCTK variables for thorn CarpetEvolutionMask
   @enddesc
 @@*/


#define THORN_IS_CarpetEvolutionMask 1

#include <stddef.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameter.h"
#include "cctk_WarnLevel.h"
#include "cctki_Groups.h"
#include "cctki_FortranWrappers.h"

int CactusBindingsVariables_CarpetEvolutionMask_Initialise(void);
static int CCTKi_BindingsFortranWrapperCarpetEvolutionMask(void *_GH, void *fpointer);

static int CCTKi_BindingsFortranWrapperCarpetEvolutionMask(void *_GH, void *fpointer)
{
  cGH *GH = _GH;
  const int _cctk_zero = 0;
#ifndef CCTK_DEBUG
  CCTK_COMPLEX _cctk_dummy_var[4];
#endif
  void (*function)(CARPETEVOLUTIONMASK_C2F_PROTO);
  DECLARE_CARPETEVOLUTIONMASK_C2F
  INITIALISE_CARPETEVOLUTIONMASK_C2F
  (void) (_cctk_zero + 0);
#ifndef CCTK_DEBUG
  (void) (_cctk_dummy_var + 0);
#endif

  function = (void (*) (CARPETEVOLUTIONMASK_C2F_PROTO)) fpointer;
  function (PASS_CARPETEVOLUTIONMASK_C2F (GH));

  return (0);
}

int CactusBindingsVariables_CarpetEvolutionMask_Initialise(void)
{
  const char * warn_mixeddim_gfs = "";
  int warn_mixeddim = 0;
  const CCTK_INT *allow_mixeddim_gfs;


  allow_mixeddim_gfs = CCTK_ParameterGet ("allow_mixeddim_gfs", "Cactus", 0);

  if (CCTKi_CreateGroup ("EVOLUTION_MASK", "CARPETEVOLUTIONMASK", "CARPETEVOLUTIONMASK",
                         "GF", "REAL", "PUBLIC",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "prolongation=\"none\" InterpNumTimelevels=1",
                         NULL,
                         1,
                         "evolution_mask") == 1)
  {
    warn_mixeddim_gfs = "evolution_mask";
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

  CCTKi_RegisterFortranWrapper("CarpetEvolutionMask", CCTKi_BindingsFortranWrapperCarpetEvolutionMask);

  return 0;
}
