#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine kerr_initial_data(CCTK_ARGUMENTS)
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
  if (vx_boost==0.d0 .and. x0_bh==0.d0) then 
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
     else if (coordinate_type==3) then 
        call kerr_initial_metric_puncture(ext, X, Y, Z, PhysicalRadius,       &
                       gxx, gxy, gxz, gyy, gyz, gzz, trK,             &
                       Axx, Axy, Axz, Ayy, Ayz, Azz, lapm1, shiftx,shifty,shiftz, &
                       gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, phi, sam, vx_boost)
     else
        call kerr_initial_metric_quasiisotropic(ext, X, Y, Z, PhysicalRadius,       &
                       gxx, gxy, gxz, gyy, gyz, gzz, trK,             &
                       Axx, Axy, Axz, Ayy, Ayz, Azz, lapm1, shiftx,shifty,shiftz, &
                       gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, phi, sam, vx_boost)
     end if
  else 
     if (coordinate_type .lt. 3) then 
        write(*,*) 'Boost or off center BH for coordinate_type <3 are not supported!'
	stop
     end if

     if (coordinate_type==3) then 
        call compute_boost_metric_Aij(1.d0,sam,x0_bh,vx_boost,ext,X,Y,Z, &
                        lapm1,shiftx,shifty,shiftz,phi,gxx,gxy,gxz, &
                        gyy,gyz,gzz,gupxx,gupxy,gupxz,gupyy, &
                        gupyz,gupzz,Axx,Axy,Axz,Ayy,Ayz,Azz,trK)
     else
        call compute_boost_metric_Aij_quasiisotropic(1.d0,sam,x0_bh,vx_boost,ext,X,Y,Z, &
                        lapm1,shiftx,shifty,shiftz,phi,gxx,gxy,gxz, &
                        gyy,gyz,gzz,gupxx,gupxy,gupxz,gupyy, &
                        gupyz,gupzz,Axx,Axy,Axz,Ayy,Ayz,Azz,trK)
     end if
  end if

  ! Use rho as temporary storage here:
  if (coordinate_type==3) then 
     call setgamma_kerr_puncture(1.d0,sam,x0_bh,vx_boost,ext,X,Y,Z,Gammax,Gammay,Gammaz)
  else if (coordinate_type==4) then 
     call setgamma_kerr_quasiisotropic(1.d0,sam,x0_bh,vx_boost,ext,X,Y,Z,Gammax,Gammay,Gammaz)
  else
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
  end if

  Psi = exp(phi)

  if (gauge_choice==1) then
     ! Setup puncture initial lapse and shift
     lapm1 = Psi**(-2) - 1.d0
     shiftx = 0.0d0
     shifty = 0.0d0
     shiftz = 0.0d0
  end if

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

end subroutine kerr_initial_data
