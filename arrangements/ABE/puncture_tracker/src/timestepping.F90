#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!-------------------------------
! Puncture tracker timestepping
!-------------------------------
subroutine puncture_tracker_timestepping(CCTK_ARGUMENTS)

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8                 :: dX,dY,dZ
  integer, dimension(3)  :: ext
  integer                :: vindex,ierr
  real*8, dimension(1)   :: punc1_position_x_rhs,punc1_position_y_rhs,punc1_position_z_rhs
  real*8, dimension(1)   :: punc2_position_x_rhs,punc2_position_y_rhs,punc2_position_z_rhs
  real*8, dimension(1,3) :: pointcoords

  dX=CCTK_DELTA_SPACE(1)
  dY=CCTK_DELTA_SPACE(2)
  dZ=CCTK_DELTA_SPACE(3)
  ext = cctk_lsh

  pointcoords(1,1)=xbh1
  pointcoords(1,2)=ybh1
  pointcoords(1,3)=zbh1

  puncture_position_rhs = 0.D0

  write(*,*) "hi! puncture tracker dxyz",dx,dy,dz
  write(*,*) "hi! puncture tracker",xbh1,ybh1,zbh1

  if(Symmetry.ne.4) then
     call CCTK_VarIndex(vindex,"shift::shiftx")
     call interp_driver_carp(cctkGH,1,pointcoords,vindex,punc1_position_x_rhs)
     call CCTK_VarIndex(vindex,"shift::shifty")
     call interp_driver_carp(cctkGH,1,pointcoords,vindex,punc1_position_y_rhs)
!     call interpolate_pointset(cctkGH,cctk_nghostzones,1,2,ext, &
!          global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d, &
!          Xglobal,Yglobal,Zglobal,pointcoords,shiftx,punc1_position_x_rhs)
     puncture_position_rhs(1,2,1) = -punc1_position_x_rhs(1)
!     call interpolate_pointset(cctkGH,cctk_nghostzones,1,2,ext, &
!          global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d, &
!          Xglobal,Yglobal,Zglobal,pointcoords,shifty,punc1_position_y_rhs)
     puncture_position_rhs(1,2,2) = -punc1_position_y_rhs(1)

     pointcoords(1,1)=xbh2
     pointcoords(1,2)=ybh2
     pointcoords(1,3)=zbh2
     call CCTK_VarIndex(vindex,"shift::shiftx")
     call interp_driver_carp(cctkGH,1,pointcoords,vindex,punc2_position_x_rhs)
     call CCTK_VarIndex(vindex,"shift::shifty")
     call interp_driver_carp(cctkGH,1,pointcoords,vindex,punc2_position_y_rhs)
!     call interpolate_pointset(cctkGH,cctk_nghostzones,1,2,ext, &
!          global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d, &
!          Xglobal,Yglobal,Zglobal,pointcoords,shiftx,punc2_position_x_rhs)
     puncture_position_rhs(2,2,1) = -punc2_position_x_rhs(1)
!     call interpolate_pointset(cctkGH,cctk_nghostzones,1,2,ext, &
!          global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d, &
!          Xglobal,Yglobal,Zglobal,pointcoords,shifty,punc2_position_y_rhs)
     puncture_position_rhs(2,2,2) = -punc2_position_y_rhs(1)

  end if
  if (Symmetry.ne.1) then
     write(*,*) "Inside puncture_tracker_timestepping: WARNING: SYMMETRY != 1 NOT SUPPORTED!"
     call CCTK_VarIndex(vindex,"shift::shiftz")
     call interp_driver_carp(cctkGH,1,pointcoords,vindex,punc1_position_z_rhs)
!     call interpolate_pointset(cctkGH,cctk_nghostzones,1,2,ext, &
!          global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d, &
!          Xglobal,Yglobal,Zglobal,pointcoords,shiftz,punc1_position_z_rhs)
     puncture_position_rhs(1,2,3) = -punc1_position_z_rhs(1)
  end if

end subroutine puncture_tracker_timestepping
