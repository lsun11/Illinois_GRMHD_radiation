!--------------------------------------
! Update boundary driver for opl lapse
!--------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine opl_lapse_bc_driver(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS   

  interface
     subroutine opl_lapse_bc(ext,fake_ext, dT,dx,dy,dz, x, y, z,PhysicalRadius, RadiusDerivative, &
          lapse_new,lapse, phi, &
          Xcenter,Ycenter,Zcenter,Symmetry, &
          gxx, gxy, gxz, gyy, gyz, gzz, &
          have_global_bdry_min,have_global_bdry_max)
       integer, dimension(3)                   :: ext,fake_ext,have_global_bdry_min,have_global_bdry_max
       real*8                                  :: dT,dx,dy,dz,Xcenter,Ycenter,Zcenter
       real*8, dimension(ext(1),ext(2),ext(3)) :: X,Y,Z
       real*8, dimension(ext(1),ext(2),ext(3)) :: PhysicalRadius, RadiusDerivative
       real*8, dimension(ext(1),ext(2),ext(3)) :: lapse_new,lapse, phi
       real*8, dimension(ext(1),ext(2),ext(3)) :: gxx, gxy, gxz, gyy, gyz, gzz
       integer                                 :: Symmetry
     end subroutine opl_lapse_bc
  end interface
  integer               :: dummy,i
  integer, dimension(3) :: ext,fake_ext,thisproc_have_global_bdry_min,thisproc_have_global_bdry_max
  real*8                :: dT,dX,dY,dZ,levelnumber
  integer               :: PI_SYMM, AXISYM
  parameter(PI_SYMM = 3, AXISYM = 4)
  !
  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)
  dT = CCTK_DELTA_TIME

  levelnumber = cctk_levfac(1)
  levelnumber = log(levelnumber)/log(2.D0)+1.D0

  if(CCTK_ITERATION.gt.0 .and. levelnumber==1.D0) then
     thisproc_have_global_bdry_min = have_global_bdry_min(int(levelnumber),:)
     thisproc_have_global_bdry_max = have_global_bdry_max(int(levelnumber),:)

     do i=1,cctk_nghostzones(2)
        fake_ext = cctk_lsh - cctk_nghostzones + i
        if(Symmetry==AXISYM) fake_ext(1) = fake_ext(1) + 1
        call opl_lapse_bc(ext,fake_ext,dT,dx,dy,dz,X,Y,Z,PhysicalRadius, RadiusDerivative, &
             lapm1, lapm1_p,phi_p, &
             Xcenter,Ycenter,Zcenter,Symmetry, &
             gxx_p, gxy_p, gxz_p, gyy_p, gyz_p, gzz_p, &
             thisproc_have_global_bdry_min,thisproc_have_global_bdry_max)
     end do

     if(Symmetry == AXISYM) then
        !first rotate, then compute values at the boundary. 
        !(needed since BndCartoon2D can only(?) be called inside a (CCTK_ARGUMENTS) function)
        call BndCartoon2DVN(dummy, cctkGH, 0, 'lapse::lapm1')
     end if

  end if

  ! TEST: Set a floor on the lapse.
  if(opl_lapse_floor.ne.0.D0) then
     where(lapm1.lt.opl_lapse_floor-1.D0)
        lapm1 = opl_lapse_floor-1.D0
     end where
  end if

end subroutine opl_lapse_bc_driver
