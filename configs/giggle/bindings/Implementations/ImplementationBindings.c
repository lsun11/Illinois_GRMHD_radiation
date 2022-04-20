void CCTKi_BindingsThorn_ADMMacros(void);
void CCTKi_BindingsThorn_AEILocalInterp(void);
void CCTKi_BindingsThorn_Boundary(void);
void CCTKi_BindingsThorn_Cactus(void);
void CCTKi_BindingsThorn_Carpet(void);
void CCTKi_BindingsThorn_CarpetEvolutionMask(void);
void CCTKi_BindingsThorn_CarpetIOASCII(void);
void CCTKi_BindingsThorn_CarpetIOBasic(void);
void CCTKi_BindingsThorn_CarpetIOHDF5(void);
void CCTKi_BindingsThorn_CarpetIntegrateTest(void);
void CCTKi_BindingsThorn_CarpetInterp(void);
void CCTKi_BindingsThorn_CarpetLib(void);
void CCTKi_BindingsThorn_CarpetReduce(void);
void CCTKi_BindingsThorn_CarpetRegrid(void);
void CCTKi_BindingsThorn_CarpetRegrid2(void);
void CCTKi_BindingsThorn_CarpetSlab(void);
void CCTKi_BindingsThorn_CarpetTracker(void);
void CCTKi_BindingsThorn_CartGrid3D(void);
void CCTKi_BindingsThorn_Cartoon2D(void);
void CCTKi_BindingsThorn_CoordBase(void);
void CCTKi_BindingsThorn_Dissipation(void);
void CCTKi_BindingsThorn_GSL(void);
void CCTKi_BindingsThorn_HDF5(void);
void CCTKi_BindingsThorn_IDScalarWaveMoL(void);
void CCTKi_BindingsThorn_IDWaveMoL(void);
void CCTKi_BindingsThorn_IOASCII(void);
void CCTKi_BindingsThorn_IOBasic(void);
void CCTKi_BindingsThorn_IOUtil(void);
void CCTKi_BindingsThorn_InitBase(void);
void CCTKi_BindingsThorn_LocalInterp(void);
void CCTKi_BindingsThorn_LocalReduce(void);
void CCTKi_BindingsThorn_LoopControl(void);
void CCTKi_BindingsThorn_MoL(void);
void CCTKi_BindingsThorn_NaNChecker(void);
void CCTKi_BindingsThorn_OS_collapse(void);
void CCTKi_BindingsThorn_OS_toy(void);
void CCTKi_BindingsThorn_SpaceMask(void);
void CCTKi_BindingsThorn_SphericalSurface(void);
void CCTKi_BindingsThorn_StaticConformal(void);
void CCTKi_BindingsThorn_SymBase(void);
void CCTKi_BindingsThorn_Time(void);
void CCTKi_BindingsThorn_TwoPunctures(void);
void CCTKi_BindingsThorn_TwoPuncturesAEI(void);
void CCTKi_BindingsThorn_WaveMoL(void);
void CCTKi_BindingsThorn_WaveToyMoL(void);
void CCTKi_BindingsThorn_ahfinderdirect(void);
void CCTKi_BindingsThorn_ang_freq_ofBfield(void);
void CCTKi_BindingsThorn_bbh_cookpfeiffer_rot(void);
void CCTKi_BindingsThorn_bbhlorene(void);
void CCTKi_BindingsThorn_bhns(void);
void CCTKi_BindingsThorn_bssn(void);
void CCTKi_BindingsThorn_cylindrical2d(void);
void CCTKi_BindingsThorn_diagnostics_mhd(void);
void CCTKi_BindingsThorn_diagnostics_vacuum(void);
void CCTKi_BindingsThorn_disk_powerlaw(void);
void CCTKi_BindingsThorn_em_extraction(void);
void CCTKi_BindingsThorn_excision(void);
void CCTKi_BindingsThorn_fisheye(void);
void CCTKi_BindingsThorn_gw_extraction(void);
void CCTKi_BindingsThorn_jobsuicide(void);
void CCTKi_BindingsThorn_lapse(void);
void CCTKi_BindingsThorn_linearized_wave(void);
void CCTKi_BindingsThorn_mag_bondi(void);
void CCTKi_BindingsThorn_magnetar(void);
void CCTKi_BindingsThorn_mhd_evolve(void);
void CCTKi_BindingsThorn_movingbox(void);
void CCTKi_BindingsThorn_movingbox_for_boosted_bh(void);
void CCTKi_BindingsThorn_nonlinear_alfven_wave(void);
void CCTKi_BindingsThorn_scalarwaveMoL(void);
void CCTKi_BindingsThorn_shift(void);
void CCTKi_BindingsThorn_shock_mhd(void);
void CCTKi_BindingsThorn_zlib(void);

int CCTKi_BindingsImplementationsInitialise(void);
int CCTKi_BindingsImplementationsInitialise(void)
{
  CCTKi_BindingsThorn_ADMMacros();
  CCTKi_BindingsThorn_AEILocalInterp();
  CCTKi_BindingsThorn_Boundary();
  CCTKi_BindingsThorn_Cactus();
  CCTKi_BindingsThorn_Carpet();
  CCTKi_BindingsThorn_CarpetEvolutionMask();
  CCTKi_BindingsThorn_CarpetIOASCII();
  CCTKi_BindingsThorn_CarpetIOBasic();
  CCTKi_BindingsThorn_CarpetIOHDF5();
  CCTKi_BindingsThorn_CarpetIntegrateTest();
  CCTKi_BindingsThorn_CarpetInterp();
  CCTKi_BindingsThorn_CarpetLib();
  CCTKi_BindingsThorn_CarpetReduce();
  CCTKi_BindingsThorn_CarpetRegrid();
  CCTKi_BindingsThorn_CarpetRegrid2();
  CCTKi_BindingsThorn_CarpetSlab();
  CCTKi_BindingsThorn_CarpetTracker();
  CCTKi_BindingsThorn_CartGrid3D();
  CCTKi_BindingsThorn_Cartoon2D();
  CCTKi_BindingsThorn_CoordBase();
  CCTKi_BindingsThorn_Dissipation();
  CCTKi_BindingsThorn_GSL();
  CCTKi_BindingsThorn_HDF5();
  CCTKi_BindingsThorn_IDScalarWaveMoL();
  CCTKi_BindingsThorn_IDWaveMoL();
  CCTKi_BindingsThorn_IOASCII();
  CCTKi_BindingsThorn_IOBasic();
  CCTKi_BindingsThorn_IOUtil();
  CCTKi_BindingsThorn_InitBase();
  CCTKi_BindingsThorn_LocalInterp();
  CCTKi_BindingsThorn_LocalReduce();
  CCTKi_BindingsThorn_LoopControl();
  CCTKi_BindingsThorn_MoL();
  CCTKi_BindingsThorn_NaNChecker();
  CCTKi_BindingsThorn_OS_collapse();
  CCTKi_BindingsThorn_OS_toy();
  CCTKi_BindingsThorn_SpaceMask();
  CCTKi_BindingsThorn_SphericalSurface();
  CCTKi_BindingsThorn_StaticConformal();
  CCTKi_BindingsThorn_SymBase();
  CCTKi_BindingsThorn_Time();
  CCTKi_BindingsThorn_TwoPunctures();
  CCTKi_BindingsThorn_TwoPuncturesAEI();
  CCTKi_BindingsThorn_WaveMoL();
  CCTKi_BindingsThorn_WaveToyMoL();
  CCTKi_BindingsThorn_ahfinderdirect();
  CCTKi_BindingsThorn_ang_freq_ofBfield();
  CCTKi_BindingsThorn_bbh_cookpfeiffer_rot();
  CCTKi_BindingsThorn_bbhlorene();
  CCTKi_BindingsThorn_bhns();
  CCTKi_BindingsThorn_bssn();
  CCTKi_BindingsThorn_cylindrical2d();
  CCTKi_BindingsThorn_diagnostics_mhd();
  CCTKi_BindingsThorn_diagnostics_vacuum();
  CCTKi_BindingsThorn_disk_powerlaw();
  CCTKi_BindingsThorn_em_extraction();
  CCTKi_BindingsThorn_excision();
  CCTKi_BindingsThorn_fisheye();
  CCTKi_BindingsThorn_gw_extraction();
  CCTKi_BindingsThorn_jobsuicide();
  CCTKi_BindingsThorn_lapse();
  CCTKi_BindingsThorn_linearized_wave();
  CCTKi_BindingsThorn_mag_bondi();
  CCTKi_BindingsThorn_magnetar();
  CCTKi_BindingsThorn_mhd_evolve();
  CCTKi_BindingsThorn_movingbox();
  CCTKi_BindingsThorn_movingbox_for_boosted_bh();
  CCTKi_BindingsThorn_nonlinear_alfven_wave();
  CCTKi_BindingsThorn_scalarwaveMoL();
  CCTKi_BindingsThorn_shift();
  CCTKi_BindingsThorn_shock_mhd();
  CCTKi_BindingsThorn_zlib();
  return 0;
}
