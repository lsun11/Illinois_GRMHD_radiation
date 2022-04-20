#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine kerr_initial_data_part2(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,global_ext
  real*8 			           :: dT,dX,dY,dZ
  real*8, parameter  :: ZERO = 0.D0
  integer :: n1,n2,n3,mf, dummy
  real*8 :: r0,rh_out,rh_in


  ext = cctk_lsh

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ! Compute Aij and convert to tilde metric
  if (abs(vx_boost) .gt. 0.d0) then 
! *** TEST ***
write(*,*) 'id part 2, gxx = ',gxx(46,30,4),dx
! ************
     ! Storage on input: 
     ! g0tt,x -> gxxy, g0tx,x -> gxyy, g0ty,x -> gxzy g0tz,x -> gyyy.
     ! g0_{mu nu} is the unboosted metric
     ! gij stores the boosted physical 3-metric
     ! After the call, gij -> tilde{g}ij
     call compute_Aij_boost(cctkGH, dx,  dy,  dz, X,Y,Z, & 
	     cctk_nghostzones, ext, lapm1, phi, &
             shiftx,shifty,shiftz,gxx, gxy,  & 
             gxz, gyy, gyz, gzz, gupxx, gupxy, gupxz, & 
             gupyy, gupyz, gupzz, gxxy, gxyy,gxzy,gyyy, &
             gxxx, gxyx, gxzx, gyyx, gyzx, gzzx, &
             !!Gammaxxx, Gammaxxy, Gammaxxz, Gammaxyy, Gammaxyz, Gammaxzz, &
             !!Gammayxx, Gammayxy, Gammayxz, Gammayyy, Gammayyz, Gammayzz, &
             !!Gammazxx, Gammazxy, Gammazxz, Gammazyy, Gammazyz, Gammazzz, & 
             vx_boost, coordinate_type, sam, Axx, Axy, Axz, Ayy, Ayz, Azz, trK)
  end if

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')

  ! Use rho as temporary storage here:
  call setgamma_v2(ext, cctk_nghostzones, dX, dY, dZ, &
     phi, gxx, gxy, gxz, gyy, gyz, gzz, &
     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Gammax, Gammay, Gammaz, rho, &
     Symmetry)

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  
  n1 = 0
  n2 = 0
  n3 = 0
  mf = 1
  !use rho as temporary storage here:
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammax,rho,n1,n2,n3,mf,Symmetry)   
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammay,rho,n1,n2,n3,mf,Symmetry)   
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammaz,rho,n1,n2,n3,mf,Symmetry)

  Psi = exp(phi)

  ! Setup puncture initial lapse and shift
  lapm1 = Psi**(-2) - 1.d0
  shiftx = 0.d0
  shifty = 0.d0
  shiftz = 0.d0

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

  !Following lines are (empirically speaking) ESSENTIAL for excision runs:
  call Derivs(ext,X,Y,Z,dX,dY,dZ,phi,phix,phiy,phiz,Symmetry)
  call Derivs(ext,X,Y,Z,dX,dY,dZ,lapm1,lapsex,lapsey,lapsez,Symmetry)

  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,phi,phix,phiy,phiz)
  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,lapm1,lapsex,lapsey,lapsez)

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')

 end subroutine kerr_initial_data_part2
