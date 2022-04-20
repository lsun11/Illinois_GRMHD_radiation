      subroutine WaveToyMoL_EstimateError (cctk_dim,cctk_gsh,cctk_lsh,cc
     &tk_lbnd,cctk_ubnd,cctk_lssh,cctk_from,cctk_to,cctk_bbox,cctk_delta
     &_time,cctk_time,cctk_delta_space,cctk_origin_space,cctk_levfac,cct
     &k_levoff,cctk_levoffdenom,cctk_timefac,cctk_convlevel,cctk_convfac
     &,cctk_nghostzones,cctk_iteration,cctkGH, X0scalarevolvedot,X0scala
     &revolveerrorestimate,X1scalarevolvedot,X1scalarevolveerrorestimate
     &,X2scalarevolvedot,X2scalarevolveerrorestimate,error,phidot,psidot
     &, X0coordinates,X0scalarevolve,X1coordinates,X1scalarevolve,X2coor
     &dinates,X2scalarevolve,coarse_dx,coarse_dy,coarse_dz,phi,phi_p,phi
     &_p_p,psi,psi_p,psi_p_p,r,x,y,z)
      implicit none
      INTEGER cctk_dim
      INTEGER cctk_gsh(cctk_dim),cctk_lsh(cctk_dim)
      INTEGER cctk_lbnd(cctk_dim),cctk_ubnd(cctk_dim)
      INTEGER cctk_lssh(3*cctk_dim)
      INTEGER cctk_from(cctk_dim),cctk_to(cctk_dim)
      INTEGER cctk_bbox(2*cctk_dim)
      REAL*8 cctk_delta_time, cctk_time
      REAL*8 cctk_delta_space(cctk_dim)
      REAL*8 cctk_origin_space(cctk_dim)
      INTEGER cctk_levfac(cctk_dim)
      INTEGER cctk_levoff(cctk_dim)
      INTEGER cctk_levoffdenom(cctk_dim)
      INTEGER cctk_timefac
      INTEGER cctk_convlevel
      INTEGER cctk_convfac
      INTEGER cctk_nghostzones(cctk_dim)
      INTEGER cctk_iteration
      integer*8 cctkGH
      INTEGER X0scalarevolvedot
      INTEGER X0scalarevolveerrorestimate
      INTEGER X1scalarevolvedot
      INTEGER X1scalarevolveerrorestimate
      INTEGER X2scalarevolvedot
      INTEGER X2scalarevolveerrorestimate
      REAL*8 error(X0scalarevolveerrorestimate,X1scalarevolveerrorestima
     &te,X2scalarevolveerrorestimate)
      REAL*8 phidot(X0scalarevolvedot,X1scalarevolvedot,X2scalarevolvedo
     &t)
      REAL*8 psidot(X0scalarevolvedot,X1scalarevolvedot,X2scalarevolvedo
     &t)
      INTEGER X0coordinates
      INTEGER X0scalarevolve
      INTEGER X1coordinates
      INTEGER X1scalarevolve
      INTEGER X2coordinates
      INTEGER X2scalarevolve
      REAL*8 coarse_dx
      REAL*8 coarse_dy
      REAL*8 coarse_dz
      REAL*8 phi(X0scalarevolve,X1scalarevolve,X2scalarevolve)
      REAL*8 phi_p(X0scalarevolve,X1scalarevolve,X2scalarevolve)
      REAL*8 phi_p_p(X0scalarevolve,X1scalarevolve,X2scalarevolve)
      REAL*8 psi(X0scalarevolve,X1scalarevolve,X2scalarevolve)
      REAL*8 psi_p(X0scalarevolve,X1scalarevolve,X2scalarevolve)
      REAL*8 psi_p_p(X0scalarevolve,X1scalarevolve,X2scalarevolve)
      REAL*8 r(X0coordinates,X1coordinates,X2coordinates)
      REAL*8 x(X0coordinates,X1coordinates,X2coordinates)
      REAL*8 y(X0coordinates,X1coordinates,X2coordinates)
      REAL*8 z(X0coordinates,X1coordinates,X2coordinates)
      
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
      
      INTEGER*4 bndwidth
      parameter (bndwidth = 1)
      integer imin(3), imax(3)
      REAL*8 dx(3)
      integer i, j, k
      integer d
      do d=1,3
         imin(d) = 1+cctk_nghostzones(d)
         imax(d) = cctk_lsh(d)-cctk_nghostzones(d)
         if (cctk_bbox(2*d-1).ne.0) imin(d) = 1+bndwidth
         if (cctk_bbox(2*d  ).ne.0) imax(d) = cctk_lsh(d)-bndwidth
      end do
      do d=1,3
         dx(d) = (cctk_delta_space(d)/cctk_levfac(d))
      end do
      do k=imin(3),imax(3)
         do j=imin(2),imax(2)
            do i=imin(1),imax(1)
               error(i,j,k) = max (
     $              abs(phi(i+1,j,k) - phi(i-1,j,k)),
     $              abs(phi(i,j+1,k) - phi(i,j-1,k)),
     $              abs(phi(i,j,k+1) - phi(i,j,k-1)))
            end do
         end do
      end do
      end
      subroutine WaveToyMoL_EstimateErrorBoundaries (cctk_dim,cctk_gsh,c
     &ctk_lsh,cctk_lbnd,cctk_ubnd,cctk_lssh,cctk_from,cctk_to,cctk_bbox,
     &cctk_delta_time,cctk_time,cctk_delta_space,cctk_origin_space,cctk_
     &levfac,cctk_levoff,cctk_levoffdenom,cctk_timefac,cctk_convlevel,cc
     &tk_convfac,cctk_nghostzones,cctk_iteration,cctkGH, X0scalarevolved
     &ot,X0scalarevolveerrorestimate,X1scalarevolvedot,X1scalarevolveerr
     &orestimate,X2scalarevolvedot,X2scalarevolveerrorestimate,error,phi
     &dot,psidot, X0coordinates,X0scalarevolve,X1coordinates,X1scalarevo
     &lve,X2coordinates,X2scalarevolve,coarse_dx,coarse_dy,coarse_dz,phi
     &,phi_p,phi_p_p,psi,psi_p,psi_p_p,r,x,y,z)
      implicit none
      INTEGER cctk_dim
      INTEGER cctk_gsh(cctk_dim),cctk_lsh(cctk_dim)
      INTEGER cctk_lbnd(cctk_dim),cctk_ubnd(cctk_dim)
      INTEGER cctk_lssh(3*cctk_dim)
      INTEGER cctk_from(cctk_dim),cctk_to(cctk_dim)
      INTEGER cctk_bbox(2*cctk_dim)
      REAL*8 cctk_delta_time, cctk_time
      REAL*8 cctk_delta_space(cctk_dim)
      REAL*8 cctk_origin_space(cctk_dim)
      INTEGER cctk_levfac(cctk_dim)
      INTEGER cctk_levoff(cctk_dim)
      INTEGER cctk_levoffdenom(cctk_dim)
      INTEGER cctk_timefac
      INTEGER cctk_convlevel
      INTEGER cctk_convfac
      INTEGER cctk_nghostzones(cctk_dim)
      INTEGER cctk_iteration
      integer*8 cctkGH
      INTEGER X0scalarevolvedot
      INTEGER X0scalarevolveerrorestimate
      INTEGER X1scalarevolvedot
      INTEGER X1scalarevolveerrorestimate
      INTEGER X2scalarevolvedot
      INTEGER X2scalarevolveerrorestimate
      REAL*8 error(X0scalarevolveerrorestimate,X1scalarevolveerrorestima
     &te,X2scalarevolveerrorestimate)
      REAL*8 phidot(X0scalarevolvedot,X1scalarevolvedot,X2scalarevolvedo
     &t)
      REAL*8 psidot(X0scalarevolvedot,X1scalarevolvedot,X2scalarevolvedo
     &t)
      INTEGER X0coordinates
      INTEGER X0scalarevolve
      INTEGER X1coordinates
      INTEGER X1scalarevolve
      INTEGER X2coordinates
      INTEGER X2scalarevolve
      REAL*8 coarse_dx
      REAL*8 coarse_dy
      REAL*8 coarse_dz
      REAL*8 phi(X0scalarevolve,X1scalarevolve,X2scalarevolve)
      REAL*8 phi_p(X0scalarevolve,X1scalarevolve,X2scalarevolve)
      REAL*8 phi_p_p(X0scalarevolve,X1scalarevolve,X2scalarevolve)
      REAL*8 psi(X0scalarevolve,X1scalarevolve,X2scalarevolve)
      REAL*8 psi_p(X0scalarevolve,X1scalarevolve,X2scalarevolve)
      REAL*8 psi_p_p(X0scalarevolve,X1scalarevolve,X2scalarevolve)
      REAL*8 r(X0coordinates,X1coordinates,X2coordinates)
      REAL*8 x(X0coordinates,X1coordinates,X2coordinates)
      REAL*8 y(X0coordinates,X1coordinates,X2coordinates)
      REAL*8 z(X0coordinates,X1coordinates,X2coordinates)
      
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
      
      INTEGER*4 bndwidth
      parameter (bndwidth = 1)
      integer ierr
      ierr = Boundary_SelectGroupForBC(cctkGH, 16383, bndwidth, -1,
     $     "WaveToyMoL::scalarevolveerrorestimate", "scalar");
      if (ierr.ne.0) call CCTK_Warn(0,48,"estimateerror.F77","WaveToyMoL
     &", "internal error")
      end
