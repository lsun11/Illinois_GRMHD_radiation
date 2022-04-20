!-----------------------------------------------------
! Driver routine for computing Ricci, and constraints
!-----------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine driver_h17(CCTK_ARGUMENTS)

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

  write(*,*) "inside driver_h17"


  write(*,*) "X(10,10,10): ",X(10,10,10)
  write(*,*) "Y(10,10,10): ",Y(10,10,10)
  write(*,*) "Z(10,10,10): ",Z(10,10,10)
  write(*,*) "gupxx(10,10,10): ",gupxx(10,10,10)
  write(*,*) "lapm1(10,10,10): ",lapm1(10,10,10)
  write(*,*) "phi(10,10,10): ",phi(10,10,10)
  call compute_h17_rhs(cctkGH,cctk_lsh, &
       			X,Y,Z, &
  			gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,lapm1,phi,h17_rhs)

end subroutine driver_h17
