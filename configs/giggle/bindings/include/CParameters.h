#include "CParameterStructNames.h"

#ifdef THORN_IS_ADMMacros
#include "ADMMacros_CParameters.h"
#endif

#ifdef THORN_IS_AEILocalInterp
#include "AEILocalInterp_CParameters.h"
#endif

#ifdef THORN_IS_Boundary
#include "Boundary_CParameters.h"
#endif

#ifdef THORN_IS_Cactus
#include "Cactus_CParameters.h"
#endif

#ifdef THORN_IS_Carpet
#include "Carpet_CParameters.h"
#endif

#ifdef THORN_IS_CarpetEvolutionMask
#include "CarpetEvolutionMask_CParameters.h"
#endif

#ifdef THORN_IS_CarpetIOASCII
#include "CarpetIOASCII_CParameters.h"
#endif

#ifdef THORN_IS_CarpetIOBasic
#include "CarpetIOBasic_CParameters.h"
#endif

#ifdef THORN_IS_CarpetIOHDF5
#include "CarpetIOHDF5_CParameters.h"
#endif

#ifdef THORN_IS_CarpetIntegrateTest
#include "CarpetIntegrateTest_CParameters.h"
#endif

#ifdef THORN_IS_CarpetInterp
#include "CarpetInterp_CParameters.h"
#endif

#ifdef THORN_IS_CarpetLib
#include "CarpetLib_CParameters.h"
#endif

#ifdef THORN_IS_CarpetReduce
#include "CarpetReduce_CParameters.h"
#endif

#ifdef THORN_IS_CarpetRegrid
#include "CarpetRegrid_CParameters.h"
#endif

#ifdef THORN_IS_CarpetRegrid2
#include "CarpetRegrid2_CParameters.h"
#endif

#ifdef THORN_IS_CarpetSlab
#include "CarpetSlab_CParameters.h"
#endif

#ifdef THORN_IS_CarpetTracker
#include "CarpetTracker_CParameters.h"
#endif

#ifdef THORN_IS_CartGrid3D
#include "CartGrid3D_CParameters.h"
#endif

#ifdef THORN_IS_Cartoon2D
#include "Cartoon2D_CParameters.h"
#endif

#ifdef THORN_IS_CoordBase
#include "CoordBase_CParameters.h"
#endif

#ifdef THORN_IS_Dissipation
#include "Dissipation_CParameters.h"
#endif

#ifdef THORN_IS_GSL
#include "GSL_CParameters.h"
#endif

#ifdef THORN_IS_HDF5
#include "HDF5_CParameters.h"
#endif

#ifdef THORN_IS_IDScalarWaveMoL
#include "IDScalarWaveMoL_CParameters.h"
#endif

#ifdef THORN_IS_IDWaveMoL
#include "IDWaveMoL_CParameters.h"
#endif

#ifdef THORN_IS_IOASCII
#include "IOASCII_CParameters.h"
#endif

#ifdef THORN_IS_IOBasic
#include "IOBasic_CParameters.h"
#endif

#ifdef THORN_IS_IOUtil
#include "IOUtil_CParameters.h"
#endif

#ifdef THORN_IS_InitBase
#include "InitBase_CParameters.h"
#endif

#ifdef THORN_IS_LocalInterp
#include "LocalInterp_CParameters.h"
#endif

#ifdef THORN_IS_LocalReduce
#include "LocalReduce_CParameters.h"
#endif

#ifdef THORN_IS_LoopControl
#include "LoopControl_CParameters.h"
#endif

#ifdef THORN_IS_MoL
#include "MoL_CParameters.h"
#endif

#ifdef THORN_IS_NaNChecker
#include "NaNChecker_CParameters.h"
#endif

#ifdef THORN_IS_OS_collapse
#include "OS_collapse_CParameters.h"
#endif

#ifdef THORN_IS_OS_toy
#include "OS_toy_CParameters.h"
#endif

#ifdef THORN_IS_SpaceMask
#include "SpaceMask_CParameters.h"
#endif

#ifdef THORN_IS_SphericalSurface
#include "SphericalSurface_CParameters.h"
#endif

#ifdef THORN_IS_StaticConformal
#include "StaticConformal_CParameters.h"
#endif

#ifdef THORN_IS_SymBase
#include "SymBase_CParameters.h"
#endif

#ifdef THORN_IS_Time
#include "Time_CParameters.h"
#endif

#ifdef THORN_IS_TwoPunctures
#include "TwoPunctures_CParameters.h"
#endif

#ifdef THORN_IS_TwoPuncturesAEI
#include "TwoPuncturesAEI_CParameters.h"
#endif

#ifdef THORN_IS_WaveMoL
#include "WaveMoL_CParameters.h"
#endif

#ifdef THORN_IS_WaveToyMoL
#include "WaveToyMoL_CParameters.h"
#endif

#ifdef THORN_IS_ahfinderdirect
#include "ahfinderdirect_CParameters.h"
#endif

#ifdef THORN_IS_ang_freq_ofBfield
#include "ang_freq_ofBfield_CParameters.h"
#endif

#ifdef THORN_IS_bbh_cookpfeiffer_rot
#include "bbh_cookpfeiffer_rot_CParameters.h"
#endif

#ifdef THORN_IS_bbhlorene
#include "bbhlorene_CParameters.h"
#endif

#ifdef THORN_IS_bhns
#include "bhns_CParameters.h"
#endif

#ifdef THORN_IS_bssn
#include "bssn_CParameters.h"
#endif

#ifdef THORN_IS_cylindrical2d
#include "cylindrical2d_CParameters.h"
#endif

#ifdef THORN_IS_diagnostics_mhd
#include "diagnostics_mhd_CParameters.h"
#endif

#ifdef THORN_IS_diagnostics_vacuum
#include "diagnostics_vacuum_CParameters.h"
#endif

#ifdef THORN_IS_disk_powerlaw
#include "disk_powerlaw_CParameters.h"
#endif

#ifdef THORN_IS_em_extraction
#include "em_extraction_CParameters.h"
#endif

#ifdef THORN_IS_excision
#include "excision_CParameters.h"
#endif

#ifdef THORN_IS_fisheye
#include "fisheye_CParameters.h"
#endif

#ifdef THORN_IS_gw_extraction
#include "gw_extraction_CParameters.h"
#endif

#ifdef THORN_IS_jobsuicide
#include "jobsuicide_CParameters.h"
#endif

#ifdef THORN_IS_lapse
#include "lapse_CParameters.h"
#endif

#ifdef THORN_IS_linearized_wave
#include "linearized_wave_CParameters.h"
#endif

#ifdef THORN_IS_mag_bondi
#include "mag_bondi_CParameters.h"
#endif

#ifdef THORN_IS_magnetar
#include "magnetar_CParameters.h"
#endif

#ifdef THORN_IS_mhd_evolve
#include "mhd_evolve_CParameters.h"
#endif

#ifdef THORN_IS_movingbox
#include "movingbox_CParameters.h"
#endif

#ifdef THORN_IS_movingbox_for_boosted_bh
#include "movingbox_for_boosted_bh_CParameters.h"
#endif

#ifdef THORN_IS_nonlinear_alfven_wave
#include "nonlinear_alfven_wave_CParameters.h"
#endif

#ifdef THORN_IS_scalarwaveMoL
#include "scalarwaveMoL_CParameters.h"
#endif

#ifdef THORN_IS_shift
#include "shift_CParameters.h"
#endif

#ifdef THORN_IS_shock_mhd
#include "shock_mhd_CParameters.h"
#endif

#ifdef THORN_IS_zlib
#include "zlib_CParameters.h"
#endif


