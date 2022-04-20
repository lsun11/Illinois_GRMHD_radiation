!----------------------------------------------------------------------------------
! Boundary condition driver for 2nd order Gamma-driving shift (Spatial_Gauge == 5)
!----------------------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine hbpuncture_shift_bc_driver(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  interface
     subroutine hbpuncture_shift_bc(ext,fake_ext,dT,dx,dy,dz, &
          X,Y,Z,r_fish,PhysicalRadius,RadiusDerivative, &
          shiftx,shifty,shiftz, &
          shiftxb,shiftyb,shiftzb,shiftxt,shiftyt,shiftzt, &
          shiftxtb,shiftytb,shiftztb,lapm1,phi,Symmetry, &
          Xcenter,Ycenter,Zcenter, &
          gxx,gxy,gxz,gyy,gyz,gzz, &
          punc_shift_bc_radial_falloff_power,shift_bc_type, &
          have_bdry_min,have_bdry_max,use_trans_fish_phys_new,fisheye_enable)
       integer, dimension(3)                   :: ext,fake_ext,have_bdry_min,have_bdry_max
       real*8                                  :: dT,dx,dy,dz,Xcenter,Ycenter,Zcenter
       real*8, dimension(ext(1),ext(2),ext(3)) :: X,Y,Z,r_fish
       real*8, dimension(ext(1),ext(2),ext(3)) :: PhysicalRadius,RadiusDerivative
       real*8, dimension(ext(1),ext(2),ext(3)) :: shiftx,shifty,shiftz
       real*8, dimension(ext(1),ext(2),ext(3)) :: shiftxb,shiftyb,shiftzb
       real*8, dimension(ext(1),ext(2),ext(3)) :: shiftxt,shiftyt,shiftzt
       real*8, dimension(ext(1),ext(2),ext(3)) :: shiftxtb,shiftytb,shiftztb
       real*8, dimension(ext(1),ext(2),ext(3)) :: gxx,gxy,gxz,gyy,gyz,gzz
       real*8, dimension(ext(1),ext(2),ext(3)) :: lapm1,phi
       integer                                 :: Symmetry,punc_shift_bc_radial_falloff_power,fisheye_enable
       integer                                 :: shift_bc_type,use_trans_fish_phys_new
     end subroutine hbpuncture_shift_bc
  end interface
  integer                                   :: n1,n2,n3,m,i,dummy,index
  integer, dimension(3)                     :: ext,fake_ext,thisproc_have_global_bdry_min,thisproc_have_global_bdry_max
  real*8                                    :: dT,dX,dY,dZ,levelnumber
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
        call hbpuncture_shift_bc(ext,fake_ext,dT,dx,dy,dz, &
             X,Y,Z,r,PhysicalRadius,RadiusDerivative, &
             shiftx,shifty,shiftz, &
             shiftx_p,shifty_p,shiftz_p, &
             shiftxt,shiftyt,shiftzt, &
             shiftxt_p,shiftyt_p,shiftzt_p, &
             lapm1_p,phi_p,Symmetry, &
             Xcenter,Ycenter,Zcenter, &
             gxx_p,gxy_p,gxz_p,gyy_p,gyz_p,gzz_p, &
             punc_shift_bc_radial_falloff_power,shift_bc_type, &
             thisproc_have_global_bdry_min,thisproc_have_global_bdry_max,use_trans_fish_phys_new,fisheye_enable)
     end do

     if(Symmetry == AXISYM) then
        !first rotate, then compute values at the boundary. 
        !(needed since BndCartoon2D can only(?) be called inside a (CCTK_ARGUMENTS) function)
        call CCTK_VarIndex(index,'shift::shiftx')
        call BndCartoon2DVI(dummy, cctkGH, 1, index)
     endif
  end if
end subroutine hbpuncture_shift_bc_driver
