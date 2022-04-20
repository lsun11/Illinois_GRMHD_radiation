!-----------------------------------------------------
! Driver routine for computing Ricci, and constraints
!-----------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine driver_spatial_derivs(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8                  :: dT, dX,dY,dZ,detmin,detmax
  integer, dimension(3)   :: ext,fake_ext
  integer                 :: dummy,i,j,k,index
  integer, parameter      :: AXISYM = 4
  !
  dT = CCTK_DELTA_TIME
  dX=CCTK_DELTA_SPACE(1)
  dY=CCTK_DELTA_SPACE(2)
  dZ=CCTK_DELTA_SPACE(3)

  ext = cctk_lsh

  gxxx=0.d0
  gxxy=0.d0
  gxxz=0.d0
  gxyx=0.d0
  gxyy=0.d0
  gxyz=0.d0
  gxzx=0.d0
  gxzy=0.d0
  gxzz=0.d0
  gyyx=0.d0
  gyyy=0.d0
  gyyz=0.d0
  gyzx=0.d0
  gyzy=0.d0
  gyzz=0.d0
  gzzx=0.d0
  gzzy=0.d0
  gzzz=0.d0
  lapsex=0.d0
  lapsey=0.d0
  lapsez=0.d0
  shiftxx=0.d0
  shiftxy=0.d0
  shiftxz=0.d0
  shiftyx=0.d0
  shiftyy=0.d0
  shiftyz=0.d0
  shiftzx=0.d0
  shiftzy=0.d0
  shiftzz=0.d0
 call spatial_metric_derivs(cctkGH,  dT,  dx,  dy,  dz, &
       cctk_nghostzones, cctk_lsh, &
       gxx, gxy, gxz, gyy, gyz, gzz, &
       lapm1,shiftx,shifty,shiftz, &
       lapsex,lapsey,lapsez, &
       shiftxx,shiftxy,shiftxz, &
       shiftyx,shiftyy,shiftyz, &
       shiftzx,shiftzy,shiftzz, &
       gxxx, gxxy, gxxz, &
       gxyx, gxyy, gxyz, &
       gxzx, gxzy, gxzz, &
       gyyx, gyyy, gyyz, &
       gyzx, gyzy, gyzz, &
       gzzx, gzzy, gzzz)
  call spatial_metric_derivs_4(cctkGH,  dT,  dx,  dy,  dz, &
       cctk_nghostzones, cctk_lsh, &
       gxx, gxy, gxz, gyy, gyz, gzz, &
       lapm1,shiftx,shifty,shiftz, &
       lapsex,lapsey,lapsez, &
       shiftxx,shiftxy,shiftxz, &
       shiftyx,shiftyy,shiftyz, &
       shiftzx,shiftzy,shiftzz, &
       gxxx, gxxy, gxxz, &
       gxyx, gxyy, gxyz, &
       gxzx, gxzy, gxzz, &
       gyyx, gyyy, gyyz, &
       gyzx, gyzy, gyzz, &
       gzzx, gzzy, gzzz)
  call spatial_metric_derivs_2(cctkGH,  dT,  dx,  dy,  dz, &
       cctk_nghostzones, cctk_lsh, &
       lapm1,gxx, gxy, gxz, gyy, gyz, gzz, &
       lapsex,lapsey,lapsez, &
       shiftx,shifty,shiftz, &
       shiftxx,shiftxy,shiftxz, &
       shiftyx,shiftyy,shiftyz, &
       shiftzx,shiftzy,shiftzz, &
       gxxx, gxxy, gxxz, &
       gxyx, gxyy, gxyz, &
       gxzx, gxzy, gxzz, &
       gyyx, gyyy, gyyz, &
       gyzx, gyzy, gyzz, &
       gzzx, gzzy, gzzz)

  !Following is needed for accurate ADM mass calculation, among other things.
  call Derivs(ext,X,Y,Z,dX,dY,dZ,phi,phix,phiy,phiz,Symmetry)
  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,phi,phix,phiy,phiz)


end subroutine driver_spatial_derivs
