!---------------------------------------------------------------
! Timestepping driver for hyperbolic shift (Spatial_Gauge == 1)
!---------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine hyperbolic_shift_timestepping(CCTK_ARGUMENTS)
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  interface
     subroutine hyperbolic_shift_rhs(cctkGH,cctk_nghostzones,ext, &
          hyper_b1,hyper_b2,hyper_b3, &
          RadiusDerivative,lapm1, &
          shiftx_rhs,shifty_rhs,shiftz_rhs, &
          shiftxt,shiftyt,shiftzt, &
          shiftxt_rhs,shiftyt_rhs,shiftzt_rhs, &
          Gammax_rhs,Gammay_rhs,Gammaz_rhs, &
          Gammax,Gammay,Gammaz, &
          Gammax_drive,Gammay_drive,Gammaz_drive)
       implicit none
       CCTK_POINTER                            :: cctkGH
       integer, dimension(3)                   :: cctk_nghostzones,ext
       real*8, dimension(ext(1),ext(2),ext(3)) :: RadiusDerivative,lapm1
       real*8                                  :: hyper_b1,hyper_b2,hyper_b3
       real*8, dimension(ext(1),ext(2),ext(3)) :: shiftx_rhs,shifty_rhs,shiftz_rhs
       real*8, dimension(ext(1),ext(2),ext(3)) :: shiftxt,shiftyt,shiftzt
       real*8, dimension(ext(1),ext(2),ext(3)) :: shiftxt_rhs,shiftyt_rhs,shiftzt_rhs
       real*8, dimension(ext(1),ext(2),ext(3)) :: Gammax_rhs,Gammay_rhs,Gammaz_rhs
       real*8, dimension(ext(1),ext(2),ext(3)) :: Gammax,Gammay,Gammaz
       real*8, dimension(ext(1),ext(2),ext(3)) :: Gammax_drive,Gammay_drive,Gammaz_drive
     end subroutine hyperbolic_shift_rhs
  end interface
  integer                                   :: n1,n2,n3,m,i,dummy,index
  integer, dimension(3)                     :: ext,fake_ext
  real*8                                    :: dT,dX,dY,dZ,psi6c,phic
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

  if(excision_enable==1) then
     call vector_excision_bc(ext,X,Y,Z, &
          shiftxt,shiftyt,shiftzt, &
          Symmetry,excision_zone_gf)
  end if

  call hyperbolic_shift_rhs(cctkGH,cctk_nghostzones,cctk_lsh, &
       hyper_b1,hyper_b2,hyper_b3, &
       RadiusDerivative,lapm1, &
       shiftx_rhs,shifty_rhs,shiftz_rhs, &
       shiftxt,shiftyt,shiftzt, &
       shiftxt_rhs,shiftyt_rhs,shiftzt_rhs, &
       Gammax_rhs,Gammay_rhs,Gammaz_rhs, &
       Gammax,Gammay,Gammaz, &
       Gammax_drive,Gammay_drive,Gammaz_drive);

end subroutine hyperbolic_shift_timestepping
