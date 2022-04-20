/*@@
   @file    diagnostics_mhd.c
   @author  Automatically generated by GridFuncStuff.pl
   @desc
            Creates the CCTK variables for thorn diagnostics_mhd
   @enddesc
 @@*/


#define THORN_IS_diagnostics_mhd 1

#include <stddef.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameter.h"
#include "cctk_WarnLevel.h"
#include "cctki_Groups.h"
#include "cctki_FortranWrappers.h"

int CactusBindingsVariables_diagnostics_mhd_Initialise(void);
static int CCTKi_BindingsFortranWrapperdiagnostics_mhd(void *_GH, void *fpointer);

static int CCTKi_BindingsFortranWrapperdiagnostics_mhd(void *_GH, void *fpointer)
{
  cGH *GH = _GH;
  const int _cctk_zero = 0;
#ifndef CCTK_DEBUG
  CCTK_COMPLEX _cctk_dummy_var[4];
#endif
  void (*function)(DIAGNOSTICS_MHD_C2F_PROTO);
  DECLARE_DIAGNOSTICS_MHD_C2F
  INITIALISE_DIAGNOSTICS_MHD_C2F
  (void) (_cctk_zero + 0);
#ifndef CCTK_DEBUG
  (void) (_cctk_dummy_var + 0);
#endif

  function = (void (*) (DIAGNOSTICS_MHD_C2F_PROTO)) fpointer;
  function (PASS_DIAGNOSTICS_MHD_C2F (GH));

  return (0);
}

int CactusBindingsVariables_diagnostics_mhd_Initialise(void)
{
  const char * warn_mixeddim_gfs = "";
  int warn_mixeddim = 0;
  const CCTK_INT *allow_mixeddim_gfs;


  allow_mixeddim_gfs = CCTK_ParameterGet ("allow_mixeddim_gfs", "Cactus", 0);

  if (CCTKi_CreateGroup ("VOLINTEGRALS_MHD", "DIAGNOSTICS_MHD", "DIAGNOSTICS_MHD",
                         "SCALAR", "REAL", "PROTECTED",
                         0, 1,
                         "NONE", "CONSTANT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         68,
                         "T_VolInt",
                         "M0_VolInt",
                         "M0_AH_VolInt",
                         "M0_escape30M",
                         "M0_escape50M",
                         "M0_escape70M",
                         "M0_escape100M",
                         "Minternal_VolInt",
                         "Minternal_cold_VolInt",
                         "em_energy_VolInt",
                         "em_energy2_VolInt",
                         "em_energy_outsideBH_VolInt",
                         "b_phi_VolInt",
                         "CoMx_VolInt",
                         "CoMy_VolInt",
                         "CoMz_VolInt",
                         "CoM_VolInt_denominator",
                         "monopole_VolInt",
                         "monopole_outsideBH_VolInt",
                         "brem_qei_VolInt",
                         "brem_qee_VolInt",
                         "synch_VolInt",
                         "M0_horiz_VolInt",
                         "M0_r1_VolInt",
                         "M0_r2_VolInt",
                         "M0_r3_VolInt",
                         "fluid_energy_horiz_VolInt",
                         "fluid_energy_r1_VolInt",
                         "fluid_energy_r2_VolInt",
                         "fluid_energy_r3_VolInt",
                         "em_energy_between_VolInt",
                         "fluid_J_horiz_VolInt",
                         "fluid_J_r1_VolInt",
                         "fluid_J_r2_VolInt",
                         "fluid_J_r3_VolInt",
                         "minternal_horiz_VolInt",
                         "minternal_r1_VolInt",
                         "minternal_r2_VolInt",
                         "minternal_r3_VolInt",
                         "minternal_cold_horiz_VolInt",
                         "minternal_cold_r1_VolInt",
                         "minternal_cold_r2_VolInt",
                         "minternal_cold_r3_VolInt",
                         "em_J_between_VolInt",
                         "half_b2_u0_VolInt",
                         "Tem0_0_VolInt",
                         "half_b2_u0_outsideBH_VolInt",
                         "Tem0_0_outsideBH_VolInt",
                         "Tfluid0_0_VolInt",
                         "Tfluid0_0_outsideBH_VolInt",
                         "em_energy_outsideradius1_VolInt",
                         "em_energy_outsideradius2_VolInt",
                         "density_modes_r0",
                         "density_modes_r1",
                         "density_modes_i1",
                         "density_modes_r2",
                         "density_modes_i2",
                         "density_modes_r3",
                         "density_modes_i3",
                         "density_modes_r4",
                         "density_modes_i4",
                         "density_modes_r5",
                         "density_modes_i5",
                         "density_modes_r6",
                         "density_modes_i6",
                         "rad_energy_VolInt",
                         "rad_energy_nue_VolInt",
                         "rad_energy_nux_VolInt") == 1)
  {
    warn_mixeddim_gfs = "volIntegrals_mhd";
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

  CCTKi_RegisterFortranWrapper("diagnostics_mhd", CCTKi_BindingsFortranWrapperdiagnostics_mhd);

  return 0;
}
