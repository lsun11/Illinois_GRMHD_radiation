/* Defines for thorn this file is part of */

#ifdef THORN_IS_ADMMacros
#define CCTK_THORN ADMMacros
#define CCTK_THORNSTRING "ADMMacros"
#define CCTK_ARRANGEMENT CactusEinstein
#define CCTK_ARRANGEMENTSTRING "CactusEinstein"
#endif

#ifdef THORN_IS_AEILocalInterp
#define CCTK_THORN AEILocalInterp
#define CCTK_THORNSTRING "AEILocalInterp"
#define CCTK_ARRANGEMENT CactusBase
#define CCTK_ARRANGEMENTSTRING "CactusBase"
#endif

#ifdef THORN_IS_Boundary
#define CCTK_THORN Boundary
#define CCTK_THORNSTRING "Boundary"
#define CCTK_ARRANGEMENT CactusBase
#define CCTK_ARRANGEMENTSTRING "CactusBase"
#endif

#ifdef THORN_IS_Cactus
#define CCTK_THORN src
#define CCTK_THORNSTRING "src"
#define CCTK_ARRANGEMENT illinois_grmhd_radtest_0903
#define CCTK_ARRANGEMENTSTRING "illinois_grmhd_radtest_0903"
#endif

#ifdef THORN_IS_Carpet
#define CCTK_THORN Carpet
#define CCTK_THORNSTRING "Carpet"
#define CCTK_ARRANGEMENT Carpet
#define CCTK_ARRANGEMENTSTRING "Carpet"
#endif

#ifdef THORN_IS_CarpetEvolutionMask
#define CCTK_THORN CarpetEvolutionMask
#define CCTK_THORNSTRING "CarpetEvolutionMask"
#define CCTK_ARRANGEMENT Carpet
#define CCTK_ARRANGEMENTSTRING "Carpet"
#endif

#ifdef THORN_IS_CarpetIOASCII
#define CCTK_THORN CarpetIOASCII
#define CCTK_THORNSTRING "CarpetIOASCII"
#define CCTK_ARRANGEMENT Carpet
#define CCTK_ARRANGEMENTSTRING "Carpet"
#endif

#ifdef THORN_IS_CarpetIOBasic
#define CCTK_THORN CarpetIOBasic
#define CCTK_THORNSTRING "CarpetIOBasic"
#define CCTK_ARRANGEMENT Carpet
#define CCTK_ARRANGEMENTSTRING "Carpet"
#endif

#ifdef THORN_IS_CarpetIOHDF5
#define CCTK_THORN CarpetIOHDF5
#define CCTK_THORNSTRING "CarpetIOHDF5"
#define CCTK_ARRANGEMENT Carpet
#define CCTK_ARRANGEMENTSTRING "Carpet"
#endif

#ifdef THORN_IS_CarpetIntegrateTest
#define CCTK_THORN CarpetIntegrateTest
#define CCTK_THORNSTRING "CarpetIntegrateTest"
#define CCTK_ARRANGEMENT CarpetExtra
#define CCTK_ARRANGEMENTSTRING "CarpetExtra"
#endif

#ifdef THORN_IS_CarpetInterp
#define CCTK_THORN CarpetInterp
#define CCTK_THORNSTRING "CarpetInterp"
#define CCTK_ARRANGEMENT Carpet
#define CCTK_ARRANGEMENTSTRING "Carpet"
#endif

#ifdef THORN_IS_CarpetLib
#define CCTK_THORN CarpetLib
#define CCTK_THORNSTRING "CarpetLib"
#define CCTK_ARRANGEMENT Carpet
#define CCTK_ARRANGEMENTSTRING "Carpet"
#endif

#ifdef THORN_IS_CarpetReduce
#define CCTK_THORN CarpetReduce
#define CCTK_THORNSTRING "CarpetReduce"
#define CCTK_ARRANGEMENT Carpet
#define CCTK_ARRANGEMENTSTRING "Carpet"
#endif

#ifdef THORN_IS_CarpetRegrid
#define CCTK_THORN CarpetRegrid
#define CCTK_THORNSTRING "CarpetRegrid"
#define CCTK_ARRANGEMENT Carpet
#define CCTK_ARRANGEMENTSTRING "Carpet"
#endif

#ifdef THORN_IS_CarpetRegrid2
#define CCTK_THORN CarpetRegrid2
#define CCTK_THORNSTRING "CarpetRegrid2"
#define CCTK_ARRANGEMENT Carpet
#define CCTK_ARRANGEMENTSTRING "Carpet"
#endif

#ifdef THORN_IS_CarpetSlab
#define CCTK_THORN CarpetSlab
#define CCTK_THORNSTRING "CarpetSlab"
#define CCTK_ARRANGEMENT Carpet
#define CCTK_ARRANGEMENTSTRING "Carpet"
#endif

#ifdef THORN_IS_CarpetTracker
#define CCTK_THORN CarpetTracker
#define CCTK_THORNSTRING "CarpetTracker"
#define CCTK_ARRANGEMENT Carpet
#define CCTK_ARRANGEMENTSTRING "Carpet"
#endif

#ifdef THORN_IS_CartGrid3D
#define CCTK_THORN CartGrid3D
#define CCTK_THORNSTRING "CartGrid3D"
#define CCTK_ARRANGEMENT CactusBase
#define CCTK_ARRANGEMENTSTRING "CactusBase"
#endif

#ifdef THORN_IS_Cartoon2D
#define CCTK_THORN Cartoon2D
#define CCTK_THORNSTRING "Cartoon2D"
#define CCTK_ARRANGEMENT BetaThorns
#define CCTK_ARRANGEMENTSTRING "BetaThorns"
#endif

#ifdef THORN_IS_CoordBase
#define CCTK_THORN CoordBase
#define CCTK_THORNSTRING "CoordBase"
#define CCTK_ARRANGEMENT CactusBase
#define CCTK_ARRANGEMENTSTRING "CactusBase"
#endif

#ifdef THORN_IS_Dissipation
#define CCTK_THORN Dissipation
#define CCTK_THORNSTRING "Dissipation"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_GSL
#define CCTK_THORN GSL
#define CCTK_THORNSTRING "GSL"
#define CCTK_ARRANGEMENT CactusExternal
#define CCTK_ARRANGEMENTSTRING "CactusExternal"
#endif

#ifdef THORN_IS_HDF5
#define CCTK_THORN HDF5
#define CCTK_THORNSTRING "HDF5"
#define CCTK_ARRANGEMENT CactusBase
#define CCTK_ARRANGEMENTSTRING "CactusBase"
#endif

#ifdef THORN_IS_IDScalarWaveMoL
#define CCTK_THORN IDScalarWaveMoL
#define CCTK_THORNSTRING "IDScalarWaveMoL"
#define CCTK_ARRANGEMENT CarpetExtra
#define CCTK_ARRANGEMENTSTRING "CarpetExtra"
#endif

#ifdef THORN_IS_IDWaveMoL
#define CCTK_THORN IDWaveMoL
#define CCTK_THORNSTRING "IDWaveMoL"
#define CCTK_ARRANGEMENT CactusExamples
#define CCTK_ARRANGEMENTSTRING "CactusExamples"
#endif

#ifdef THORN_IS_IOASCII
#define CCTK_THORN IOASCII
#define CCTK_THORNSTRING "IOASCII"
#define CCTK_ARRANGEMENT CactusBase
#define CCTK_ARRANGEMENTSTRING "CactusBase"
#endif

#ifdef THORN_IS_IOBasic
#define CCTK_THORN IOBasic
#define CCTK_THORNSTRING "IOBasic"
#define CCTK_ARRANGEMENT CactusBase
#define CCTK_ARRANGEMENTSTRING "CactusBase"
#endif

#ifdef THORN_IS_IOUtil
#define CCTK_THORN IOUtil
#define CCTK_THORNSTRING "IOUtil"
#define CCTK_ARRANGEMENT CactusBase
#define CCTK_ARRANGEMENTSTRING "CactusBase"
#endif

#ifdef THORN_IS_InitBase
#define CCTK_THORN InitBase
#define CCTK_THORNSTRING "InitBase"
#define CCTK_ARRANGEMENT CactusBase
#define CCTK_ARRANGEMENTSTRING "CactusBase"
#endif

#ifdef THORN_IS_LocalInterp
#define CCTK_THORN LocalInterp
#define CCTK_THORNSTRING "LocalInterp"
#define CCTK_ARRANGEMENT CactusBase
#define CCTK_ARRANGEMENTSTRING "CactusBase"
#endif

#ifdef THORN_IS_LocalReduce
#define CCTK_THORN LocalReduce
#define CCTK_THORNSTRING "LocalReduce"
#define CCTK_ARRANGEMENT CactusBase
#define CCTK_ARRANGEMENTSTRING "CactusBase"
#endif

#ifdef THORN_IS_LoopControl
#define CCTK_THORN LoopControl
#define CCTK_THORNSTRING "LoopControl"
#define CCTK_ARRANGEMENT Carpet
#define CCTK_ARRANGEMENTSTRING "Carpet"
#endif

#ifdef THORN_IS_MoL
#define CCTK_THORN MoL
#define CCTK_THORNSTRING "MoL"
#define CCTK_ARRANGEMENT CactusBase
#define CCTK_ARRANGEMENTSTRING "CactusBase"
#endif

#ifdef THORN_IS_NaNChecker
#define CCTK_THORN NaNChecker
#define CCTK_THORNSTRING "NaNChecker"
#define CCTK_ARRANGEMENT CactusUtils
#define CCTK_ARRANGEMENTSTRING "CactusUtils"
#endif

#ifdef THORN_IS_OS_collapse
#define CCTK_THORN OS_collapse
#define CCTK_THORNSTRING "OS_collapse"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_OS_toy
#define CCTK_THORN OS_toy
#define CCTK_THORNSTRING "OS_toy"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_SpaceMask
#define CCTK_THORN SpaceMask
#define CCTK_THORNSTRING "SpaceMask"
#define CCTK_ARRANGEMENT CactusEinstein
#define CCTK_ARRANGEMENTSTRING "CactusEinstein"
#endif

#ifdef THORN_IS_SphericalSurface
#define CCTK_THORN SphericalSurface
#define CCTK_THORNSTRING "SphericalSurface"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_StaticConformal
#define CCTK_THORN StaticConformal
#define CCTK_THORNSTRING "StaticConformal"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_SymBase
#define CCTK_THORN SymBase
#define CCTK_THORNSTRING "SymBase"
#define CCTK_ARRANGEMENT CactusBase
#define CCTK_ARRANGEMENTSTRING "CactusBase"
#endif

#ifdef THORN_IS_Time
#define CCTK_THORN Time
#define CCTK_THORNSTRING "Time"
#define CCTK_ARRANGEMENT CactusBase
#define CCTK_ARRANGEMENTSTRING "CactusBase"
#endif

#ifdef THORN_IS_TwoPunctures
#define CCTK_THORN TwoPunctures
#define CCTK_THORNSTRING "TwoPunctures"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_TwoPuncturesAEI
#define CCTK_THORN TwoPuncturesAEI
#define CCTK_THORNSTRING "TwoPuncturesAEI"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_WaveMoL
#define CCTK_THORN WaveMoL
#define CCTK_THORNSTRING "WaveMoL"
#define CCTK_ARRANGEMENT CactusExamples
#define CCTK_ARRANGEMENTSTRING "CactusExamples"
#endif

#ifdef THORN_IS_WaveToyMoL
#define CCTK_THORN WaveToyMoL
#define CCTK_THORNSTRING "WaveToyMoL"
#define CCTK_ARRANGEMENT CarpetExtra
#define CCTK_ARRANGEMENTSTRING "CarpetExtra"
#endif

#ifdef THORN_IS_ahfinderdirect
#define CCTK_THORN ahfinderdirect
#define CCTK_THORNSTRING "ahfinderdirect"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_ang_freq_ofBfield
#define CCTK_THORN ang_freq_ofBfield
#define CCTK_THORNSTRING "ang_freq_ofBfield"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_bbh_cookpfeiffer_rot
#define CCTK_THORN bbh_cookpfeiffer_rot
#define CCTK_THORNSTRING "bbh_cookpfeiffer_rot"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_bbhlorene
#define CCTK_THORN bbhlorene
#define CCTK_THORNSTRING "bbhlorene"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_bhns
#define CCTK_THORN bhns
#define CCTK_THORNSTRING "bhns"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_bssn
#define CCTK_THORN bssn
#define CCTK_THORNSTRING "bssn"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_cylindrical2d
#define CCTK_THORN cylindrical2d
#define CCTK_THORNSTRING "cylindrical2d"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_diagnostics_mhd
#define CCTK_THORN diagnostics_mhd
#define CCTK_THORNSTRING "diagnostics_mhd"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_diagnostics_vacuum
#define CCTK_THORN diagnostics_vacuum
#define CCTK_THORNSTRING "diagnostics_vacuum"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_disk_powerlaw
#define CCTK_THORN disk_powerlaw
#define CCTK_THORNSTRING "disk_powerlaw"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_em_extraction
#define CCTK_THORN em_extraction
#define CCTK_THORNSTRING "em_extraction"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_excision
#define CCTK_THORN excision
#define CCTK_THORNSTRING "excision"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_fisheye
#define CCTK_THORN fisheye
#define CCTK_THORNSTRING "fisheye"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_gw_extraction
#define CCTK_THORN gw_extraction
#define CCTK_THORNSTRING "gw_extraction"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_jobsuicide
#define CCTK_THORN jobsuicide
#define CCTK_THORNSTRING "jobsuicide"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_lapse
#define CCTK_THORN lapse
#define CCTK_THORNSTRING "lapse"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_linearized_wave
#define CCTK_THORN linearized_wave
#define CCTK_THORNSTRING "linearized_wave"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_mag_bondi
#define CCTK_THORN mag_bondi
#define CCTK_THORNSTRING "mag_bondi"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_magnetar
#define CCTK_THORN magnetar
#define CCTK_THORNSTRING "magnetar"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_mhd_evolve
#define CCTK_THORN mhd_evolve
#define CCTK_THORNSTRING "mhd_evolve"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_movingbox
#define CCTK_THORN movingbox
#define CCTK_THORNSTRING "movingbox"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_movingbox_for_boosted_bh
#define CCTK_THORN movingbox_for_boosted_bh
#define CCTK_THORNSTRING "movingbox_for_boosted_bh"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_nonlinear_alfven_wave
#define CCTK_THORN nonlinear_alfven_wave
#define CCTK_THORNSTRING "nonlinear_alfven_wave"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_scalarwaveMoL
#define CCTK_THORN scalarwaveMoL
#define CCTK_THORNSTRING "scalarwaveMoL"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_shift
#define CCTK_THORN shift
#define CCTK_THORNSTRING "shift"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_shock_mhd
#define CCTK_THORN shock_mhd
#define CCTK_THORNSTRING "shock_mhd"
#define CCTK_ARRANGEMENT ABE
#define CCTK_ARRANGEMENTSTRING "ABE"
#endif

#ifdef THORN_IS_zlib
#define CCTK_THORN zlib
#define CCTK_THORNSTRING "zlib"
#define CCTK_ARRANGEMENT CactusBase
#define CCTK_ARRANGEMENTSTRING "CactusBase"
#endif

