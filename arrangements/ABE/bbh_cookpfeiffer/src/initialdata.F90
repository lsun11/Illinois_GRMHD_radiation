#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bbhcp_setup_initialdata(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,global_ext
  real*8 			           :: dT,dX,dY,dZ
  real*8 			           :: n1,n2,n3,mf

  real*8                                   :: PI,dmass,detmin,detmax
  real*8                                   :: px1,py1,pz1,px2,py2,pz2
  real*8                                   :: m1,m2,y1,y2
  real*8, dimension(1,3)                   :: pointcoords

  integer :: index,dummy,foundit_flag
  integer :: handle,ierr
  real*8  :: ONE,ZERO,F1o3,F1o12,TWO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(ONE = 1.D0, ZERO = 0.D0, F1o3 = 1.D0/3.D0, F1o12 = 1.D0/12.D0, TWO = 2.D0)

  INTEGER :: i,j,k

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  PI = 3.14159265358979323846D0

  ext = cctk_lsh

!!  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')
  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')

  psi = exp(phi)

!! call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_AH')
  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_AH')

  ! Okay, now phi and psi are correct everywhere, even ghostzones!

  if(use_cookpfeiffer_lapseshift==0) then
     lapm1 = psi**(-2) - 1.D0
     shiftx = 0.D0
     shifty = 0.D0
     shiftz = 0.D0
  end if



  ! Next, set the conformal metric and its inverse
  gxx = 1.D0
  gxy = 0.D0
  gxz = 0.D0
  gyy = 1.D0
  gyz = 0.D0
  gzz = 1.D0

  gupxx = 1.D0
  gupxy = 0.D0
  gupxz = 0.D0
  gupyy = 1.D0
  gupyz = 0.D0
  gupzz = 1.D0

  trK = 0.D0

  ! Now compute Kzz from trK = 0 and all the other Kij's:
  ! This line is valid since gupij = deltaij
  Kzz = - (Kxx + Kyy)

! Old code snippet:
!  trK =  gupxx * Kxx + gupyy * Kyy + gupzz * Kzz +  &
!       TWO * ( gupxy * Kxy + gupxz * Kxz + gupyz * Kyz )

  !Compute Atilde_ij's
  Axx = Kxx / (Psi*Psi*Psi*Psi) - (1.D0/3.D0) * gxx * trK
  Axy = Kxy / (Psi*Psi*Psi*Psi) - (1.D0/3.D0) * gxy * trK
  Axz = Kxz / (Psi*Psi*Psi*Psi) - (1.D0/3.D0) * gxz * trK
  Ayy = Kyy / (Psi*Psi*Psi*Psi) - (1.D0/3.D0) * gyy * trK
  Ayz = Kyz / (Psi*Psi*Psi*Psi) - (1.D0/3.D0) * gyz * trK
  Azz = Kzz / (Psi*Psi*Psi*Psi) - (1.D0/3.D0) * gzz * trK

  ! Set excision_zone_gf to avoid valgrind memory errors
  excision_zone_gf = 0

  shiftxt = 0.D0
  shiftyt = 0.D0
  shiftzt = 0.D0
  lapset = 0.D0

  Gammax = 0.D0
  Gammay = 0.D0
  Gammaz = 0.D0

  !Look, Gamma^i was just set to zero,
  ! so next lines not strictly needed for brandt-brugmann ID, 
  ! but doesn't hurt in case this code is ever recycled!
  n1 = 0
  n2 = 0
  n3 = 0
  mf = 1
  !use gxx_p as temporary storage here:
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammax,gxx_p,n1,n2,n3,mf,Symmetry)
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammay,gxx_p,n1,n2,n3,mf,Symmetry)
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammaz,gxx_p,n1,n2,n3,mf,Symmetry)
  gxx_p=gxx

  !Next line NEEDED!
!!  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')


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

  if(fisheye_enable .eq. 1) then
     call trans_phys_fish_tensor(cctkGH,cctk_lsh,Symmetry, &
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
  end if
  
  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')
  call CartSymGN(dummy,cctkGH,'shift::shift_vars')
  
  if(Symmetry==AXISYM) then
     write(*,*) "I PITY DA FOO TRIES TO USE AXISYMMETRY FOR ORBITING BLACK HOLES!"
     stop
  end if
  
  !Following lines are (empirically speaking) ESSENTIAL for excision runs:
  call Derivs(ext,X,Y,Z,dX,dY,dZ,phi,phix,phiy,phiz,Symmetry) 
  call Derivs(ext,X,Y,Z,dX,dY,dZ,lapm1,lapsex,lapsey,lapsez,Symmetry) 

  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,phi,phix,phiy,phiz)
  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,lapm1,lapsex,lapsey,lapsez)
  
  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
!!  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars') 
  
end subroutine bbhcp_setup_initialdata
