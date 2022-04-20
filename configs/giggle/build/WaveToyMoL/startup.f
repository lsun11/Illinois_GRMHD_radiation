      subroutine WaveToyMoL_Startup
      implicit none
      integer      CCTK_Equals, CCTK_MyProc, CCTK_nProcs, CCTK_IsThornAc
     &tive
      external     CCTK_Equals, CCTK_MyProc, CCTK_nProcs, CCTK_IsThornAc
     &tive
      integer*8 CCTK_PointerTo, CCTK_NullPointer
      external     CCTK_PointerTo, CCTK_NullPointer
      external Boundary_SelectGroupForBC
      INTEGER*4 Boundary_SelectGroupForBC
      external MoLRegisterEvolvedGroup
      INTEGER*4 MoLRegisterEvolvedGroup
      
      REAL*8  excision_radius
      REAL*8  run_time
      INTEGER*4 Symmetry
      INTEGER*4 bssn_enable
      INTEGER*4 cowling_enable
      INTEGER*4 excision_enable
      INTEGER*4 fisheye_enable
      INTEGER*4 iter_count
      INTEGER*4 number_of_mol_ministeps
      INTEGER*4 rot_metric
      INTEGER*4 trA_detg_enforce
      COMMON /cctk_params_global/excision_radius,run_time,Symmetry,bssn_
     &enable,cowling_enable,excision_enable,fisheye_enable,iter_count,nu
     &mber_of_mol_ministeps,rot_metric,trA_detg_enforce
      INTEGER*4 WaveToyMoL_MaxNumEvolvedVars
      COMMON /WaveToyMoLrest/WaveToyMoL_MaxNumEvolvedVars
      integer*8  bound
      INTEGER*4 estimate_error
      INTEGER*4 order
      COMMON /WaveToyMoLpriv/bound,estimate_error,order
      INTEGER*4 CCTKH0
      INTEGER*4 CCTKH1
      INTEGER*4 CCTKH2
      INTEGER*4 CCTKH3
      INTEGER*4 CCTKH4
      INTEGER*4 CCTKH5
      INTEGER*4 CCTKH6
      INTEGER*4 CCTKH7
      INTEGER*4 CCTKH8
      INTEGER*4 CCTKH9
      INTEGER*4 CCTKH10
      INTEGER*4 CCTKH11
      INTEGER*4 MoL_Num_Evolved_Vars
      INTEGER*4 CCTKH12
      INTEGER*4 CCTKH13
      COMMON /METHODOFLINESrest/CCTKH0,CCTKH1,CCTKH2,CCTKH3,CCTKH4,CCTKH
     &5,CCTKH6,CCTKH7,CCTKH8,CCTKH9,CCTKH10,CCTKH11,MoL_Num_Evolved_Vars
     &,CCTKH12,CCTKH13
      
      integer ierr
      call CCTK_RegisterBanner (ierr, "WaveToyMoL: Evolving a scalar fie
     &ld with MoL")
      if (ierr.ne.0) call CCTK_Warn(0,11,"startup.F77","WaveToyMoL", "in
     &ternal error")
      end
