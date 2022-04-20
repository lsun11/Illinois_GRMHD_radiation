!-----------------------------------------------------------------------------
!
!$Id: predict_mhd.F90  $
!
!-----------------------------------------------------------------------------
!
! Predict & Correct RHS routines for mhd variables
!
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
!
! Predictor
!
!-----------------------------------------------------------------------------
subroutine interp_center(cctkGH,cctk_nghostzones,funct,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,finalvalue)
  
  implicit none
  
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_POINTER                             :: cctkGH   
  integer,dimension(3)                     :: ext,cctk_nghostzones,global_ext
  integer                                  :: numpoints,interpolation_order
  real*8, dimension(1,3)                   :: pointcoords
  real*8, dimension(1)                     :: outputdata
  real*8,dimension(ext(1))                 :: Xlocal1d
  real*8,dimension(ext(2))                 :: Ylocal1d
  real*8,dimension(ext(3))                 :: Zlocal1d
  real*8,dimension(global_ext(1))          :: Xglobal
  real*8,dimension(global_ext(2))          :: Yglobal
  real*8,dimension(global_ext(3))          :: Zglobal
  real*8                                   :: dX,dY,dZ,finalvalue
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in):: funct

  numpoints = 1
  interpolation_order = 2
  pointcoords = 0.D0

  call interpolate_pointset(cctkGH,cctk_nghostzones,numpoints,interpolation_order,ext,global_ext, &
       dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal, &
       pointcoords,funct,outputdata) 
  
  finalvalue = outputdata(1)

end subroutine interp_center
