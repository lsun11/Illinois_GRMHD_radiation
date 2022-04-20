!-------------------------------------------------------
! Driver for hyperbolic shift BC's (Spatial_Gauge == 1)
!-------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine hyperbolic_shift_bc_driver(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  interface
     subroutine hyper_shift_bc(ext,fake_ext,dT,dx,dy,dz,X,Y,Z,r_fish, &
          PhysicalRadius, RadiusDerivative,shiftx,shifty,shiftz, &
          shiftxb,shiftyb,shiftzb,lapm1,phi,Symmetry, &
          Xcenter,Ycenter,Zcenter, &
          gxx,gxy,gxz,gyy,gyz,gzz, &
          shift_bc_type,have_global_bdry_min,have_global_bdry_max,use_trans_fish_phys_new,fisheye_enable)
       integer, dimension(3)                   :: ext,fake_ext,have_global_bdry_min,have_global_bdry_max
       real*8                                  :: dT,dx,dy,dz,Xcenter,Ycenter,Zcenter
       real*8, dimension(ext(1),ext(2),ext(3)) :: r_fish,PhysicalRadius, RadiusDerivative
       real*8, dimension(ext(1),ext(2),ext(3)) :: X,Y,Z
       real*8, dimension(ext(1),ext(2),ext(3)) :: shiftx,shifty,shiftz
       real*8, dimension(ext(1),ext(2),ext(3)) :: shiftxb,shiftyb,shiftzb
       real*8, dimension(ext(1),ext(2),ext(3)) :: gxx,gxy,gxz,gyy,gyz,gzz
       real*8, dimension(ext(1),ext(2),ext(3)) :: lapm1,phi
       integer                                 :: Symmetry,shift_bc_type,use_trans_fish_phys_new,fisheye_enable
     end subroutine hyper_shift_bc
  end interface
  integer                                   :: n1,n2,n3,m,i,dummy,vindex
  integer, dimension(3)                     :: ext,fake_ext,thisproc_have_global_bdry_min,thisproc_have_global_bdry_max
  real*8                                    :: dT,dX,dY,dZ,psi6c,phic,levelnumber
  real*8                                    :: HALF,ld_eps,ld_c
  integer                    :: PI_SYMM, AXISYM
  parameter(PI_SYMM = 3, AXISYM = 4)
  parameter ( HALF = 0.5D0 )
  !
  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)
  dT = CCTK_DELTA_TIME

  levelnumber = cctk_levfac(1)
  levelnumber = log(levelnumber)/log(2.D0)+1.D0

  if(CCTK_ITERATION.gt.0 .and. levelnumber.eq.1.D0) then
     thisproc_have_global_bdry_min = have_global_bdry_min(int(levelnumber),:)
     thisproc_have_global_bdry_max = have_global_bdry_max(int(levelnumber),:)

     do i=1,cctk_nghostzones(2)
        fake_ext = cctk_lsh - cctk_nghostzones + i
        if(Symmetry==AXISYM) fake_ext(1) = fake_ext(1) + 1
        call hyper_shift_bc(ext,fake_ext,dT,dx,dy,dz,X,Y,Z,r,PhysicalRadius, RadiusDerivative, &
             shiftx,shifty,shiftz, &
             shiftx_p,shifty_p,shiftz_p, &
             lapm1_p,phi_p,Symmetry, &
             Xcenter,Ycenter,Zcenter, &
             gxx_p,gxy_p,gxz_p,gyy_p,gyz_p,gzz_p, &
             shift_bc_type,thisproc_have_global_bdry_min,thisproc_have_global_bdry_max,use_trans_fish_phys_new,fisheye_enable)
     end do
  end if

  !  call interp_center(cctkGH,cctk_nghostzones,phi_p,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,phic)
  !FIXME: Interpolation in Carpet can only be performed in a GLOBAL context in the schedule.ccl file.
  !  To fix this, you'll need to call the interpolation routine from a GLOBAL context (this function is called
  !  in a LOCAL context only), then compute "psi6c" and "call subtract_radial2" will need to be done in a 
  !  LOCAL context.
!  call CCTK_VarIndex (vindex, "bssn::phi")
!  call interp_center_carp(cctkGH,vindex,phic)
!  psi6c = exp(6.D0*phic)
!  call subtract_radial2(ext,X,Y,Z,PhysicalRadius,shiftx,shifty,shiftz,psi6c,hyper_psi6_init,hyper_Rscale);

  if(Symmetry == AXISYM) then
     !first rotate, then compute values at the boundary. 
     !(needed since BndCartoon2D can only(?) be called inside a (CCTK_ARGUMENTS) function)
     call CCTK_VarIndex(vindex,'shift::shiftx')
     call BndCartoon2DVI(dummy, cctkGH, 1, vindex)
  end if

end subroutine hyperbolic_shift_bc_driver
