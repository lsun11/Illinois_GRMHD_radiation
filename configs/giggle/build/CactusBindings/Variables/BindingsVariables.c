/*@@
   @file    BindingsVariables.c
   @author  Automatically generated by GridFuncStuff.pl
   @desc
            Calls the variable binding routines for all thorns
   @enddesc
 @@*/


#include "cctk_ActiveThorns.h"

int CactusBindingsVariables_ADMMacros_Initialise(void);
int CactusBindingsVariables_AEILocalInterp_Initialise(void);
int CactusBindingsVariables_Boundary_Initialise(void);
int CactusBindingsVariables_Cactus_Initialise(void);
int CactusBindingsVariables_Carpet_Initialise(void);
int CactusBindingsVariables_CarpetEvolutionMask_Initialise(void);
int CactusBindingsVariables_CarpetIOASCII_Initialise(void);
int CactusBindingsVariables_CarpetIOBasic_Initialise(void);
int CactusBindingsVariables_CarpetIOHDF5_Initialise(void);
int CactusBindingsVariables_CarpetIntegrateTest_Initialise(void);
int CactusBindingsVariables_CarpetInterp_Initialise(void);
int CactusBindingsVariables_CarpetLib_Initialise(void);
int CactusBindingsVariables_CarpetReduce_Initialise(void);
int CactusBindingsVariables_CarpetRegrid_Initialise(void);
int CactusBindingsVariables_CarpetRegrid2_Initialise(void);
int CactusBindingsVariables_CarpetSlab_Initialise(void);
int CactusBindingsVariables_CarpetTracker_Initialise(void);
int CactusBindingsVariables_CartGrid3D_Initialise(void);
int CactusBindingsVariables_Cartoon2D_Initialise(void);
int CactusBindingsVariables_CoordBase_Initialise(void);
int CactusBindingsVariables_Dissipation_Initialise(void);
int CactusBindingsVariables_GSL_Initialise(void);
int CactusBindingsVariables_HDF5_Initialise(void);
int CactusBindingsVariables_IDScalarWaveMoL_Initialise(void);
int CactusBindingsVariables_IDWaveMoL_Initialise(void);
int CactusBindingsVariables_IOASCII_Initialise(void);
int CactusBindingsVariables_IOBasic_Initialise(void);
int CactusBindingsVariables_IOUtil_Initialise(void);
int CactusBindingsVariables_InitBase_Initialise(void);
int CactusBindingsVariables_LocalInterp_Initialise(void);
int CactusBindingsVariables_LocalReduce_Initialise(void);
int CactusBindingsVariables_LoopControl_Initialise(void);
int CactusBindingsVariables_MoL_Initialise(void);
int CactusBindingsVariables_NaNChecker_Initialise(void);
int CactusBindingsVariables_OS_collapse_Initialise(void);
int CactusBindingsVariables_OS_toy_Initialise(void);
int CactusBindingsVariables_SpaceMask_Initialise(void);
int CactusBindingsVariables_SphericalSurface_Initialise(void);
int CactusBindingsVariables_StaticConformal_Initialise(void);
int CactusBindingsVariables_SymBase_Initialise(void);
int CactusBindingsVariables_Time_Initialise(void);
int CactusBindingsVariables_TwoPunctures_Initialise(void);
int CactusBindingsVariables_TwoPuncturesAEI_Initialise(void);
int CactusBindingsVariables_WaveMoL_Initialise(void);
int CactusBindingsVariables_WaveToyMoL_Initialise(void);
int CactusBindingsVariables_ahfinderdirect_Initialise(void);
int CactusBindingsVariables_ang_freq_ofBfield_Initialise(void);
int CactusBindingsVariables_bbh_cookpfeiffer_rot_Initialise(void);
int CactusBindingsVariables_bbhlorene_Initialise(void);
int CactusBindingsVariables_bhns_Initialise(void);
int CactusBindingsVariables_bssn_Initialise(void);
int CactusBindingsVariables_cylindrical2d_Initialise(void);
int CactusBindingsVariables_diagnostics_mhd_Initialise(void);
int CactusBindingsVariables_diagnostics_vacuum_Initialise(void);
int CactusBindingsVariables_disk_powerlaw_Initialise(void);
int CactusBindingsVariables_em_extraction_Initialise(void);
int CactusBindingsVariables_excision_Initialise(void);
int CactusBindingsVariables_fisheye_Initialise(void);
int CactusBindingsVariables_gw_extraction_Initialise(void);
int CactusBindingsVariables_jobsuicide_Initialise(void);
int CactusBindingsVariables_lapse_Initialise(void);
int CactusBindingsVariables_linearized_wave_Initialise(void);
int CactusBindingsVariables_mag_bondi_Initialise(void);
int CactusBindingsVariables_magnetar_Initialise(void);
int CactusBindingsVariables_mhd_evolve_Initialise(void);
int CactusBindingsVariables_movingbox_Initialise(void);
int CactusBindingsVariables_movingbox_for_boosted_bh_Initialise(void);
int CactusBindingsVariables_nonlinear_alfven_wave_Initialise(void);
int CactusBindingsVariables_scalarwaveMoL_Initialise(void);
int CactusBindingsVariables_shift_Initialise(void);
int CactusBindingsVariables_shock_mhd_Initialise(void);
int CactusBindingsVariables_zlib_Initialise(void);

int CCTKi_BindingsVariablesInitialise(void);

int CCTKi_BindingsVariablesInitialise(void)
{
  if (CCTK_IsThornActive("ADMMacros"))
  {
    CactusBindingsVariables_ADMMacros_Initialise();
  }
  if (CCTK_IsThornActive("AEILocalInterp"))
  {
    CactusBindingsVariables_AEILocalInterp_Initialise();
  }
  if (CCTK_IsThornActive("Boundary"))
  {
    CactusBindingsVariables_Boundary_Initialise();
  }
  if (CCTK_IsThornActive("Cactus"))
  {
    CactusBindingsVariables_Cactus_Initialise();
  }
  if (CCTK_IsThornActive("Carpet"))
  {
    CactusBindingsVariables_Carpet_Initialise();
  }
  if (CCTK_IsThornActive("CarpetEvolutionMask"))
  {
    CactusBindingsVariables_CarpetEvolutionMask_Initialise();
  }
  if (CCTK_IsThornActive("CarpetIOASCII"))
  {
    CactusBindingsVariables_CarpetIOASCII_Initialise();
  }
  if (CCTK_IsThornActive("CarpetIOBasic"))
  {
    CactusBindingsVariables_CarpetIOBasic_Initialise();
  }
  if (CCTK_IsThornActive("CarpetIOHDF5"))
  {
    CactusBindingsVariables_CarpetIOHDF5_Initialise();
  }
  if (CCTK_IsThornActive("CarpetIntegrateTest"))
  {
    CactusBindingsVariables_CarpetIntegrateTest_Initialise();
  }
  if (CCTK_IsThornActive("CarpetInterp"))
  {
    CactusBindingsVariables_CarpetInterp_Initialise();
  }
  if (CCTK_IsThornActive("CarpetLib"))
  {
    CactusBindingsVariables_CarpetLib_Initialise();
  }
  if (CCTK_IsThornActive("CarpetReduce"))
  {
    CactusBindingsVariables_CarpetReduce_Initialise();
  }
  if (CCTK_IsThornActive("CarpetRegrid"))
  {
    CactusBindingsVariables_CarpetRegrid_Initialise();
  }
  if (CCTK_IsThornActive("CarpetRegrid2"))
  {
    CactusBindingsVariables_CarpetRegrid2_Initialise();
  }
  if (CCTK_IsThornActive("CarpetSlab"))
  {
    CactusBindingsVariables_CarpetSlab_Initialise();
  }
  if (CCTK_IsThornActive("CarpetTracker"))
  {
    CactusBindingsVariables_CarpetTracker_Initialise();
  }
  if (CCTK_IsThornActive("CartGrid3D"))
  {
    CactusBindingsVariables_CartGrid3D_Initialise();
  }
  if (CCTK_IsThornActive("Cartoon2D"))
  {
    CactusBindingsVariables_Cartoon2D_Initialise();
  }
  if (CCTK_IsThornActive("CoordBase"))
  {
    CactusBindingsVariables_CoordBase_Initialise();
  }
  if (CCTK_IsThornActive("Dissipation"))
  {
    CactusBindingsVariables_Dissipation_Initialise();
  }
  if (CCTK_IsThornActive("GSL"))
  {
    CactusBindingsVariables_GSL_Initialise();
  }
  if (CCTK_IsThornActive("HDF5"))
  {
    CactusBindingsVariables_HDF5_Initialise();
  }
  if (CCTK_IsThornActive("IDScalarWaveMoL"))
  {
    CactusBindingsVariables_IDScalarWaveMoL_Initialise();
  }
  if (CCTK_IsThornActive("IDWaveMoL"))
  {
    CactusBindingsVariables_IDWaveMoL_Initialise();
  }
  if (CCTK_IsThornActive("IOASCII"))
  {
    CactusBindingsVariables_IOASCII_Initialise();
  }
  if (CCTK_IsThornActive("IOBasic"))
  {
    CactusBindingsVariables_IOBasic_Initialise();
  }
  if (CCTK_IsThornActive("IOUtil"))
  {
    CactusBindingsVariables_IOUtil_Initialise();
  }
  if (CCTK_IsThornActive("InitBase"))
  {
    CactusBindingsVariables_InitBase_Initialise();
  }
  if (CCTK_IsThornActive("LocalInterp"))
  {
    CactusBindingsVariables_LocalInterp_Initialise();
  }
  if (CCTK_IsThornActive("LocalReduce"))
  {
    CactusBindingsVariables_LocalReduce_Initialise();
  }
  if (CCTK_IsThornActive("LoopControl"))
  {
    CactusBindingsVariables_LoopControl_Initialise();
  }
  if (CCTK_IsThornActive("MoL"))
  {
    CactusBindingsVariables_MoL_Initialise();
  }
  if (CCTK_IsThornActive("NaNChecker"))
  {
    CactusBindingsVariables_NaNChecker_Initialise();
  }
  if (CCTK_IsThornActive("OS_collapse"))
  {
    CactusBindingsVariables_OS_collapse_Initialise();
  }
  if (CCTK_IsThornActive("OS_toy"))
  {
    CactusBindingsVariables_OS_toy_Initialise();
  }
  if (CCTK_IsThornActive("SpaceMask"))
  {
    CactusBindingsVariables_SpaceMask_Initialise();
  }
  if (CCTK_IsThornActive("SphericalSurface"))
  {
    CactusBindingsVariables_SphericalSurface_Initialise();
  }
  if (CCTK_IsThornActive("StaticConformal"))
  {
    CactusBindingsVariables_StaticConformal_Initialise();
  }
  if (CCTK_IsThornActive("SymBase"))
  {
    CactusBindingsVariables_SymBase_Initialise();
  }
  if (CCTK_IsThornActive("Time"))
  {
    CactusBindingsVariables_Time_Initialise();
  }
  if (CCTK_IsThornActive("TwoPunctures"))
  {
    CactusBindingsVariables_TwoPunctures_Initialise();
  }
  if (CCTK_IsThornActive("TwoPuncturesAEI"))
  {
    CactusBindingsVariables_TwoPuncturesAEI_Initialise();
  }
  if (CCTK_IsThornActive("WaveMoL"))
  {
    CactusBindingsVariables_WaveMoL_Initialise();
  }
  if (CCTK_IsThornActive("WaveToyMoL"))
  {
    CactusBindingsVariables_WaveToyMoL_Initialise();
  }
  if (CCTK_IsThornActive("ahfinderdirect"))
  {
    CactusBindingsVariables_ahfinderdirect_Initialise();
  }
  if (CCTK_IsThornActive("ang_freq_ofBfield"))
  {
    CactusBindingsVariables_ang_freq_ofBfield_Initialise();
  }
  if (CCTK_IsThornActive("bbh_cookpfeiffer_rot"))
  {
    CactusBindingsVariables_bbh_cookpfeiffer_rot_Initialise();
  }
  if (CCTK_IsThornActive("bbhlorene"))
  {
    CactusBindingsVariables_bbhlorene_Initialise();
  }
  if (CCTK_IsThornActive("bhns"))
  {
    CactusBindingsVariables_bhns_Initialise();
  }
  if (CCTK_IsThornActive("bssn"))
  {
    CactusBindingsVariables_bssn_Initialise();
  }
  if (CCTK_IsThornActive("cylindrical2d"))
  {
    CactusBindingsVariables_cylindrical2d_Initialise();
  }
  if (CCTK_IsThornActive("diagnostics_mhd"))
  {
    CactusBindingsVariables_diagnostics_mhd_Initialise();
  }
  if (CCTK_IsThornActive("diagnostics_vacuum"))
  {
    CactusBindingsVariables_diagnostics_vacuum_Initialise();
  }
  if (CCTK_IsThornActive("disk_powerlaw"))
  {
    CactusBindingsVariables_disk_powerlaw_Initialise();
  }
  if (CCTK_IsThornActive("em_extraction"))
  {
    CactusBindingsVariables_em_extraction_Initialise();
  }
  if (CCTK_IsThornActive("excision"))
  {
    CactusBindingsVariables_excision_Initialise();
  }
  if (CCTK_IsThornActive("fisheye"))
  {
    CactusBindingsVariables_fisheye_Initialise();
  }
  if (CCTK_IsThornActive("gw_extraction"))
  {
    CactusBindingsVariables_gw_extraction_Initialise();
  }
  if (CCTK_IsThornActive("jobsuicide"))
  {
    CactusBindingsVariables_jobsuicide_Initialise();
  }
  if (CCTK_IsThornActive("lapse"))
  {
    CactusBindingsVariables_lapse_Initialise();
  }
  if (CCTK_IsThornActive("linearized_wave"))
  {
    CactusBindingsVariables_linearized_wave_Initialise();
  }
  if (CCTK_IsThornActive("mag_bondi"))
  {
    CactusBindingsVariables_mag_bondi_Initialise();
  }
  if (CCTK_IsThornActive("magnetar"))
  {
    CactusBindingsVariables_magnetar_Initialise();
  }
  if (CCTK_IsThornActive("mhd_evolve"))
  {
    CactusBindingsVariables_mhd_evolve_Initialise();
  }
  if (CCTK_IsThornActive("movingbox"))
  {
    CactusBindingsVariables_movingbox_Initialise();
  }
  if (CCTK_IsThornActive("movingbox_for_boosted_bh"))
  {
    CactusBindingsVariables_movingbox_for_boosted_bh_Initialise();
  }
  if (CCTK_IsThornActive("nonlinear_alfven_wave"))
  {
    CactusBindingsVariables_nonlinear_alfven_wave_Initialise();
  }
  if (CCTK_IsThornActive("scalarwaveMoL"))
  {
    CactusBindingsVariables_scalarwaveMoL_Initialise();
  }
  if (CCTK_IsThornActive("shift"))
  {
    CactusBindingsVariables_shift_Initialise();
  }
  if (CCTK_IsThornActive("shock_mhd"))
  {
    CactusBindingsVariables_shock_mhd_Initialise();
  }
  if (CCTK_IsThornActive("zlib"))
  {
    CactusBindingsVariables_zlib_Initialise();
  }
  return 0;
}