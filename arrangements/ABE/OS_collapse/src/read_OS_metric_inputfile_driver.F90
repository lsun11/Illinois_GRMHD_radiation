
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine read_OS_metric_inputfile_driver(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8 :: dX,dY,dZ
  integer :: nx,ny,nz
  
  nx=cctk_lsh(1)
  ny=cctk_lsh(2)
  nz=cctk_lsh(3)

  !initialize
  psi=0.d0
  lapm1=0.d0
  write(*,*) "hello"
  write(*,*) "cctk_lsh: ",cctk_lsh
  write(*,*) "X(10,10,10):",X(10,10,10)
  write(*,*) "Y(10,10,10):",Y(10,10,10)
  write(*,*) "Z(10,10,10):",Z(10,10,10)
  write(*,*) "psi(10,10,10):",psi(10,10,10)
  write(*,*) "lapm1(10,10,10):",lapm1(10,10,10)
  coord=10
  call read_OS_metric_inputfile(cctkGH,psi,lapm1,X,Y,Z,nx,ny,nz,coord,narr,R_edge)
  
  if (set_lapse_one .eq. 1) then
     lapm1=0.d0
  endif

  write(*,*) "psi(10,10,10):",psi(10,10,10)
  write(*,*) "lapm1(10,10,10):",lapm1(10,10,10)
end subroutine read_OS_metric_inputfile_driver
