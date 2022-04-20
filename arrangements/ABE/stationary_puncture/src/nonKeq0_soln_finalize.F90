#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine stationary_puncture_id_nonKeq0_soln_finalize(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext
  real*8 			           :: dX,dY,dZ
  real*8 			           :: n1,n2,n3,mf
  logical :: fish_to_phys
  integer :: index,dummy
  real*8  :: ONE,ZERO
  parameter(ONE = 1.D0, ZERO = 0.D0)
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  INTEGER :: i,j,k

  ext = cctk_lsh
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  phi = log(psi)
  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call fill_lapse_symmetry_gz(ext,X,Y,Z,Symmetry,lapm1)
  call fill_shift_symmetry_gz(ext,X,Y,Z,Symmetry,shiftx,shifty,shiftz,shiftxt,shiftyt,shiftzt)
!  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
!  call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')
!  call CartSymGN(dummy,cctkGH,'shift::shift_vars')

  gxx = ONE
  gxy = ZERO
  gxz = ZERO
  gyy = ONE
  gyz = ZERO
  gzz = ONE

  gupxx = ONE
  gupxy = ZERO
  gupxz = ZERO
  gupyy = ONE
  gupyz = ZERO
  gupzz = ONE


  !Convert to fisheye coordinates, if desired
  if(fisheye_enable.eq.1) then
     call trans_phys_fish_tensor_flat(cctkGH,cctk_lsh,Symmetry, &
          x, y, z,  &
          PhysicalRadius,RadiusDerivative, &
          gxx,gxy,gxz,gyy,gyz,gzz)
     call trans_phys_fish_tensor_inv(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)
     call trans_phys_fish_phi(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative,phi)
!     fish_to_phys = .false.
!     call statpunc_trans_fish_phys_vector_fullgridfunction(ext, X, Y, Z, r, PhysicalRadius, RadiusDerivative, &
!          shiftx,shifty,shiftz,fish_to_phys)

     Sz = sqrt(shiftx*shiftx + shifty*shifty + shiftz*shiftz)
     shiftx = Sz * x/r / RadiusDerivative
     shifty = Sz * y/r / RadiusDerivative
     shiftz = Sz * z/r / RadiusDerivative
     
     psi = exp(phi)
  end if

  !======================================
  ! Note:  kset_c uses *physical* metric as input, so convert gij to physical metric
  !======================================
  Sx = psi**(4.D0)
  gxx = gxx*Sx
  gxy = gxy*Sx
  gxz = gxz*Sx
  gyy = gyy*Sx
  gyz = gyz*Sx
  gzz = gzz*Sx

  if(Symmetry==AXISYM) then
     call BndCartoon2DVN(dummy, cctkGH, 0, 'bssn::phi')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'bssn::trK')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'lapse::lapm1')
     call CCTK_VarIndex(index,'shift::shiftx')
     call BndCartoon2DVI(dummy, cctkGH, 1, index) 
     call CCTK_VarIndex(index,'BSSN::gxx')
     call BndCartoon2DVI(dummy, cctkGH, 2, index) 
  end if

  !We use all the Gammaijk's as temporary storage here:
  call kset_c_v2(ext,X,Y,Z, &
       Axx,Axy,Axz,Ayy,Ayz,Azz,trK, &
       gxx,gxy,gxz,gyy,gyz,gzz, &
       shiftx,shifty,shiftz,lapm1, &
       Symmetry, &
       Gammaxxx, Gammaxxy, Gammaxxz, Gammaxyy, Gammaxyz,  &
       Gammayxx, Gammayxy, Gammayxz, Gammayyy, & 
       gxxx,gxxy,gxxz,gxyx,gxyy,gxyz, &
       gxzx,gxzy,gxzz,gyyx,gyyy,gyyz, &
       gyzx,gyzy,gyzz,gzzx,gzzy,gzzz, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
       Gammazzz)
  
!  call statpunc_kset_c(ext,X,Y,Z, &
!       Axx,Axy,Axz,Ayy,Ayz,Azz,trK, &
!       gxx,gxy,gxz,gyy,gyz,gzz, &
!       shiftx,shifty,shiftz,lapm1, &
!       Symmetry)

  !The next function is correct to the order of finite differencing on the grid, but 
  !  does _not_ fill in the boundary ghostzones!  I.e., it only fills the grid interior.
  call statpunc_kset_c_arborder_nogz(cctkGH,cctk_lsh,cctk_nghostzones,Symmetry, &
       dX,dY,dZ, &
       Axx,Axy,Axz,Ayy,Ayz,Azz,trK, &
       gxx,gxy,gxz,gyy,gyz,gzz, &
       shiftx,shifty,shiftz,lapm1)

  !==============================
  ! Convert BACK to tilde metric
  ! use Sx for temporary array allocation...
  !==============================
  Sx = psi**(-4.D0)
  gxx = gxx*Sx
  gxy = gxy*Sx
  gxz = gxz*Sx
  gyy = gyy*Sx
  gyz = gyz*Sx
  gzz = gzz*Sx

  Sx = psi**(-4.D0)
  Axx = Axx*Sx
  Axy = Axy*Sx
  Axz = Axz*Sx
  Ayy = Ayy*Sx
  Ayz = Ayz*Sx
  Azz = Azz*Sx

  Sx =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz &
       - gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / Sx
  gupxy = - ( gxy * gzz - gyz * gxz ) / Sx
  gupxz =   ( gxy * gyz - gyy * gxz ) / Sx
  gupyy =   ( gxx * gzz - gxz * gxz ) / Sx
  gupyz = - ( gxx * gyz - gxy * gxz ) / Sx
  gupzz =   ( gxx * gyy - gxy * gxy ) / Sx  

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
!  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
!  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')
  if(Symmetry .eq. AXISYM) then 
     call CCTK_VarIndex(index,'BSSN::Axx')
     call BndCartoon2DVI(dummy, cctkGH, 2, index)
  end if

  !====================================================
  ! Get Gamma^i, using gxx_p as temporary storage
  !====================================================
  call setgamma_v2(ext,cctk_nghostzones,dX,dY,dZ, &
       phi,gxx,gxy,gxz,gyy,gyz,gzz, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
       Gammax,Gammay,Gammaz,gxx_p, &
       Symmetry)
  gxx_p = gxx

! Old version: did not set up Gamma^i's using arbitrary order in the grid interior.
!  call setgamma(ext,X,Y,Z, &
!       phi,gxx,gxy,gxz,gyy,gyz,gzz, &
!       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
!       Gammax,Gammay,Gammaz, &
!       Symmetry)

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
!  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
!  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')
  if(Symmetry .eq. AXISYM) then 
     call CCTK_VarIndex(index,'BSSN::Gammax')
     call BndCartoon2DVI(dummy, cctkGH, 1, index)
  end if
!  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')

  n1 = 0
  n2 = 0
  n3 = 0
  mf = 1
  !use gxx_p as temporary storage here:
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammax,gxx_p,n1,n2,n3,mf,Symmetry)
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammay,gxx_p,n1,n2,n3,mf,Symmetry)
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammaz,gxx_p,n1,n2,n3,mf,Symmetry)
  gxx_p=gxx

  !FIXME: Check this:
  !Next line NEEDED!
!  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')

end subroutine stationary_puncture_id_nonKeq0_soln_finalize
