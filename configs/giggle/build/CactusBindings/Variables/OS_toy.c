/*@@
   @file    OS_toy.c
   @author  Automatically generated by GridFuncStuff.pl
   @desc
            Creates the CCTK variables for thorn OS_toy
   @enddesc
 @@*/


#define THORN_IS_OS_toy 1

#include <stddef.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameter.h"
#include "cctk_WarnLevel.h"
#include "cctki_Groups.h"
#include "cctki_FortranWrappers.h"

int CactusBindingsVariables_OS_toy_Initialise(void);
static int CCTKi_BindingsFortranWrapperOS_toy(void *_GH, void *fpointer);

static int CCTKi_BindingsFortranWrapperOS_toy(void *_GH, void *fpointer)
{
  cGH *GH = _GH;
  const int _cctk_zero = 0;
#ifndef CCTK_DEBUG
  CCTK_COMPLEX _cctk_dummy_var[4];
#endif
  void (*function)(OS_TOY_C2F_PROTO);
  DECLARE_OS_TOY_C2F
  INITIALISE_OS_TOY_C2F
  (void) (_cctk_zero + 0);
#ifndef CCTK_DEBUG
  (void) (_cctk_dummy_var + 0);
#endif

  function = (void (*) (OS_TOY_C2F_PROTO)) fpointer;
  function (PASS_OS_TOY_C2F (GH));

  return (0);
}

int CactusBindingsVariables_OS_toy_Initialise(void)
{
  const char * warn_mixeddim_gfs = "";
  int warn_mixeddim = 0;
  const CCTK_INT *allow_mixeddim_gfs;


  allow_mixeddim_gfs = CCTK_ParameterGet ("allow_mixeddim_gfs", "Cactus", 0);

  if (CCTKi_CreateGroup ("MHD_OS_VOLINT", "OS_TOY", "OS_TOY",
                         "SCALAR", "REAL", "PRIVATE",
                         0, 1,
                         "NONE", "CONSTANT",
                         "", "",
                         "",
                         NULL,
                         1,
                         "OS_restmass_VolInt") == 1)
  {
    warn_mixeddim_gfs = "mhd_OS_VolInt";
    warn_mixeddim = 0;
  }
  if (CCTKi_CreateGroup ("MHD_OS_PRIVATE", "OS_TOY", "OS_TOY",
                         "GF", "REAL", "PRIVATE",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "prolongation=\"none\" InterpNumTimelevels=1",
                         NULL,
                         1,
                         "OS_rest_mass_integrand") == 1)
  {
    warn_mixeddim_gfs = "mhd_OS_private";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("OS_CENTER_DIAGNOSTICS", "OS_TOY", "OS_TOY",
                         "SCALAR", "REAL", "PRIVATE",
                         0, 1,
                         "NONE", "CONSTANT",
                         "", "",
                         "",
                         NULL,
                         1,
                         "tau_center_OS") == 1)
  {
    warn_mixeddim_gfs = "OS_center_diagnostics";
    warn_mixeddim = 0;
  }
  if (CCTKi_CreateGroup ("ANALYTIC", "OS_TOY", "OS_TOY",
                         "GF", "REAL", "PRIVATE",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "",
                         NULL,
                         1,
                         "rho_b_analytic") == 1)
  {
    warn_mixeddim_gfs = "analytic";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("PARTICLE_TRACER_COORD", "OS_TOY", "OS_TOY",
                         "ARRAY", "REAL", "PRIVATE",
                         2, 1,
                         "NONE", "CONSTANT",
                         "NARR,4", "",
                         "",
                         NULL,
                         5,
                         "coord",
                         "slope",
                         "coordt",
                         "slopet",
                         "slopem") == 1)
  {
    warn_mixeddim_gfs = "particle_tracer_coord";
    warn_mixeddim = 2;
  }
  if (CCTKi_CreateGroup ("V_PREVIOUS", "OS_TOY", "OS_TOY",
                         "GF", "REAL", "PRIVATE",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "InterpNumTimelevels=1 prolongation=\"none\"",
                         NULL,
                         4,
                         "vx_p",
                         "vy_p",
                         "vz_p",
                         "u0_p") == 1)
  {
    warn_mixeddim_gfs = "v_previous";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("PCLE_STUFF", "OS_TOY", "OS_TOY",
                         "ARRAY", "REAL", "PRIVATE",
                         1, 1,
                         "NONE", "CONSTANT",
                         "NARR", "",
                         "",
                         NULL,
                         22,
                         "areal_radius",
                         "gxx_pcle",
                         "gxz_pcle",
                         "gzz_pcle",
                         "phi_pcle",
                         "rho_b_pcle",
                         "P_pcle",
                         "u0_pcle",
                         "E_rad_pcle",
                         "Ec_pcle",
                         "eta",
                         "Q",
                         "F_radx_pcle",
                         "Fc_pcle",
                         "z_anal",
                         "rho_b_anal",
                         "Ec_anal",
                         "Fc_anal",
                         "E_anal",
                         "F_anal",
                         "E_rad_gradient",
                         "F_coeff") == 1)
  {
    warn_mixeddim_gfs = "pcle_stuff";
    warn_mixeddim = 1;
  }
  if (CCTKi_CreateGroup ("MORE_PCLE_STUFF", "OS_TOY", "OS_TOY",
                         "ARRAY", "REAL", "PRIVATE",
                         2, 1,
                         "NONE", "CONSTANT",
                         "NARR,3", "",
                         "",
                         NULL,
                         1,
                         "pos") == 1)
  {
    warn_mixeddim_gfs = "more_pcle_stuff";
    warn_mixeddim = 2;
  }
  if (CCTKi_CreateGroup ("TRACER", "OS_TOY", "OS_TOY",
                         "ARRAY", "REAL", "PRIVATE",
                         1, 1,
                         "NONE", "CONSTANT",
                         "NARR", "",
                         "",
                         NULL,
                         8,
                         "u0_p_int",
                         "vx_p_int",
                         "vy_p_int",
                         "vz_p_int",
                         "u0_int",
                         "vx_int",
                         "vy_int",
                         "vz_int") == 1)
  {
    warn_mixeddim_gfs = "Tracer";
    warn_mixeddim = 1;
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

  CCTKi_RegisterFortranWrapper("OS_toy", CCTKi_BindingsFortranWrapperOS_toy);

  return 0;
}
