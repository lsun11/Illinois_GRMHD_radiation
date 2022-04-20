#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine kerr_initial_data_part1(CCTK_ARGUMENTS)
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

  ! Set up the unboosted metric
  if (coordinate_type==1) then 
     call ks_initial_metric_plus_junk(ext, X, Y, Z, PhysicalRadius,       &
                    gxx, gxy, gxz, gyy, gyz, gzz, trK,             &
                    Axx, Axy, Axz, Ayy, Ayz, Azz, lapm1, shiftx,shifty,shiftz, &
                    gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
                    phi, sam, r_junk, lambda, vx_boost)
  else if (coordinate_type==2) then 
     call ks_initial_metric_bl_junk(ext, X, Y, Z, PhysicalRadius,       &
                    gxx, gxy, gxz, gyy, gyz, gzz, trK,             &
                    Axx, Axy, Axz, Ayy, Ayz, Azz, lapm1, shiftx,shifty,shiftz, &
                    gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
                    phi, sam, r_junk, lambda, vx_boost)
  else
     call kerr_initial_metric_puncture(ext, X, Y, Z, PhysicalRadius,       &
                    gxx, gxy, gxz, gyy, gyz, gzz, trK,             &
                    Axx, Axy, Axz, Ayy, Ayz, Azz, lapm1, shiftx,shifty,shiftz, &
                    gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, phi, sam, vx_boost)
  end if

  ! Compute the x-derivate of the unboosted spacetime metric and boost the metric
  if (abs(vx_boost) .gt. 0.d0) then
     ! temporarily store gtt -> rho, gtx -> Sx, gty -> Sy, gtz -> Sz,
     ! gtt,x -> gxxy, gtx,x -> gxyy, gty,x -> gxzy gtz,x -> gyyy.
     ! Note that gij stores the boosted physical 3-metric after the call
     call metric_derivatives_and_boost(cctkGH, dx,  dy,  dz, &
             cctk_nghostzones, ext, phi, &
             lapm1,shiftx,shifty,shiftz,rho,Sx,Sy,Sz,gxx, gxy,  &
             gxz, gyy, gyz, gzz, gxxy, gxyy,gxzy,gyyy, &
             gxxx, gxyx, gxzx, gyyx, gyzx, gzzx, vx_boost)

     call CartSymGN(dummy,cctkGH,'BSSN::BSSN_aux_restrict2')

  end if

  ! Use rho as temporary storage here:
  !!call setgamma_v2(ext, cctk_nghostzones, dX, dY, dZ, &
  !!   phi, gxx, gxy, gxz, gyy, gyz, gzz, &
  !!   gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Gammax, Gammay, Gammaz, rho, &
  !!   Symmetry)

  !!call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  !!
  !!n1 = 0
  !!n2 = 0
  !!n3 = 0
  !!mf = 1
  !!!use rho as temporary storage here:
  !!call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammax,rho,n1,n2,n3,mf,Symmetry)
  !!call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammay,rho,n1,n2,n3,mf,Symmetry)
  !!call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammaz,rho,n1,n2,n3,mf,Symmetry)

  !!Psi = exp(phi)

  !!! Note: We temporarily store the initial value of trK to lapset_rhs
  !!!lapset_rhs = trK

  !!! Setup puncture initial lapse and shift
  !!lapm1 = Psi**(-2) - 1.d0
  !!shiftx = 0.d0
  !!shifty = 0.d0
  !!shiftz = 0.d0
  !!!!! *** TEST ***
  !!!!rh_out = 1.d0 + sqrt(1.d0 - sam**2)
  !!!!rh_in  = 1.d0 - sqrt(1.d0 - sam**2)
  !!!!! BL radius inside which junk initial data will be filled
  !!!!if (r_junk .lt. 0.d0) then
  !!!!   r0 = 0.5d0*(rh_in + rh_out)
  !!!!   if (abs(sam) .gt. 0.96d0) r0 = 0.6d0
  !!!!else
  !!!!   r0 = r_junk
  !!!!end if

  !!!!shiftx = -X/PhysicalRadius * (PhysicalRadius-0.7d0)* 20.d0* PhysicalRadius**2/(1.d0+PhysicalRadius)**5
  !!!!shifty = -Y/PhysicalRadius * (PhysicalRadius-0.7d0)* 20.d0* PhysicalRadius**2/(1.d0+PhysicalRadius)**5
  !!!!shiftz = -Z/PhysicalRadius * (PhysicalRadius-0.7d0)* 20.d0* PhysicalRadius**2/(1.d0+PhysicalRadius)**5
  !!!!! ************

  !!shiftxt = 0.D0
  !!shiftyt = 0.D0
  !!shiftzt = 0.D0
  !!lapset = 0.D0



  !!! Set matter variables to zero!
  !!rho = ZERO
  !!S = ZERO
  !!Sx = ZERO
  !!Sy = ZERO
  !!Sz = ZERO
  !!Sxx = ZERO
  !!Sxy = ZERO
  !!Sxz = ZERO
  !!Syy = ZERO
  !!Syz = ZERO
  !!Szz = ZERO

  !!!Following lines are (empirically speaking) ESSENTIAL for excision runs:
  !!call Derivs(ext,X,Y,Z,dX,dY,dZ,phi,phix,phiy,phiz,Symmetry)
  !!call Derivs(ext,X,Y,Z,dX,dY,dZ,lapm1,lapsex,lapsey,lapsez,Symmetry)

  !!call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,phi,phix,phiy,phiz)
  !!call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,lapm1,lapsex,lapsey,lapsez)

  !!call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')

end subroutine kerr_initial_data_part1
