#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine stationary_puncture_initial_data_part2(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,global_ext
  real*8 			           :: dT,dX,dY,dZ

  real*8                                   :: phib,lapm1b,rb,PI,actual_excis_radius,excis_radius_fish
  real*8                                   :: shiftxb,shiftyb,shiftzb,shiftrb
  real*8                                   :: gxxb,gxyb,gxzb,gyyb,gyzb,gzzb,Gammaxb,Gammayb,Gammazb
  real*8                                   :: Axxb,Axyb,Axzb,Ayyb,Ayzb,Azzb
  real*8                                   :: gupxxb,gupxyb,gupxzb,gupyyb,gupyzb,gupzzb
  real*8, dimension(1,3)                   :: pointcoords
  ! integer                                  :: N_theta_sm,N_phi_sm,ntot_sm,numtheta_sm,numphi_sm
  !  real*8                                   :: dphi_sm,dcostheta_sm

  integer :: index,dummy,foundit_flag
  integer :: handle,ierr
  real*8  :: ONE,ZERO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(ONE = 1.D0, ZERO = 0.D0)

  INTEGER :: i,j,k !,iter,iterloop,ITERMAX=150,lo_or_hi
  !  REAL*8 :: xlo,xhi,r_iso,outputdiff,outputdiffold,mass,r_s,betar,dbdr,drdbr
  !  REAL*8 :: EPS=1.D-14

  ext = cctk_lsh

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  PI = 3.14159265358979323846D0

  !Convert to fisheye coordinates, if desired
  if(fisheye_enable.eq.1 .and. asymptotic_bh.ne.2) then
     call trans_phys_fish_tensor_flat(cctkGH,cctk_lsh,Symmetry, &
          x, y, z,  &
          PhysicalRadius,RadiusDerivative, &
          gxx,gxy,gxz,gyy,gyz,gzz)
     call trans_phys_fish_tensor_inv(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)
     call trans_phys_fish_tensor(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative, &
          Axx,Axy,Axz,Ayy,Ayz,Azz)
     call trans_phys_fish_phi(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative,phi)
     call trans_phys_fish_gamt_flat(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative,RadiusDerivative2, &
          Gammax,Gammay,Gammaz)
  endif

  Psi = exp(phi)

  shiftxt = 0.D0
  shiftyt = 0.D0
  shiftzt = 0.D0
  lapset = 0.D0

  ! Set matter variables to zero!
  rho = ZERO
  S = ZERO
  Sx = ZERO
  Sy = ZERO
  Sz = ZERO
  Sxx = ZERO
  Sxy = ZERO
  Sxz = ZERO
  Syy = ZERO
  Syz = ZERO
  Szz = ZERO

  ! call find_apparent_horizon
  !  adm::AH_Manager.Find_Horizon(level,0.0,0.0);


  !======================================
  ! Set everything else to Zero! 
  !======================================
  Gammaxxx = ZERO
  Gammaxxy = ZERO
  Gammaxxz = ZERO
  Gammaxyy = ZERO
  Gammaxyz = ZERO
  Gammaxzz = ZERO
  Gammayxx = ZERO
  Gammayxy = ZERO
  Gammayxz = ZERO
  Gammayyy = ZERO
  Gammayyz = ZERO
  Gammayzz = ZERO
  Gammazxx = ZERO
  Gammazxy = ZERO
  Gammazxz = ZERO
  Gammazyy = ZERO
  Gammazyz = ZERO
  Gammazzz = ZERO

!  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call fill_lapse_symmetry_gz(ext,X,Y,Z,Symmetry,lapm1)
!  call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')
  call fill_shift_symmetry_gz(ext,X,Y,Z,Symmetry,shiftx,shifty,shiftz,shiftxt,shiftyt,shiftzt)
!  call CartSymGN(dummy,cctkGH,'shift::shift_vars')

  if(Symmetry==AXISYM) then
     call BndCartoon2DVN(dummy, cctkGH, 0, 'bssn::phi')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'bssn::trK')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'lapse::lapm1')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'lapse::lapset')
     call CCTK_VarIndex(index,'shift::shiftx')
     call BndCartoon2DVI(dummy, cctkGH, 1, index) 
     call CCTK_VarIndex(index,'BSSN::Gammax')
     call BndCartoon2DVI(dummy, cctkGH, 1, index) 
     call CCTK_VarIndex(index,'shift::shiftxt')
     call BndCartoon2DVI(dummy, cctkGH, 1, index) 
     call CCTK_VarIndex(index,'BSSN::Axx')
     call BndCartoon2DVI(dummy, cctkGH, 2, index) 
     call CCTK_VarIndex(index,'BSSN::gxx')
     call BndCartoon2DVI(dummy, cctkGH, 2, index) 
     write(*,*) "inside stationary_puncture_initial_data.F90:1",phi(2,2,2)
     write(*,*) "inside stationary_puncture_initial_data.F90:2",trK(2,2,2)
     write(*,*) "inside stationary_puncture_initial_data.F90:3",lapm1(2,2,2)
     write(*,*) "inside stationary_puncture_initial_data.F90:4",lapset(2,2,2)
     write(*,*) "inside stationary_puncture_initial_data.F90:5",shiftx(2,2,2)
     write(*,*) "inside stationary_puncture_initial_data.F90:6",Gammax(2,2,2)
     write(*,*) "inside stationary_puncture_initial_data.F90:7",shiftxt(2,2,2)
     write(*,*) "inside stationary_puncture_initial_data.F90:8",Axx(2,2,2),gxx(2,2,2)
  end if

  !Following lines are (empirically speaking) ESSENTIAL for excision runs:
  call Derivs(ext,X,Y,Z,dX,dY,dZ,phi,phix,phiy,phiz,Symmetry) 
  call Derivs(ext,X,Y,Z,dX,dY,dZ,lapm1,lapsex,lapsey,lapsez,Symmetry) 

  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,phi,phix,phiy,phiz)
  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,lapm1,lapsex,lapsey,lapsez)

  if(Symmetry==AXISYM) then 
     call CCTK_VarIndex(index,'lapse::lapsex')
     call BndCartoon2DVI(dummy, cctkGH, 1, index) 
     call CCTK_VarIndex(index,'bssn::phix') 
     call BndCartoon2DVI(dummy, cctkGH, 1, index) 
  end if

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call fill_lapsederivs_symmetry_gz(ext,X,Y,Z,Symmetry,lapsex,lapsey,lapsez)
  !  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars') 

end subroutine stationary_puncture_initial_data_part2
