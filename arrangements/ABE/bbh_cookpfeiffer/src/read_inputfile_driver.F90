
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bbhcp_read_inputfile_driver(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8 :: dX,dY,dZ

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  call read_inputfile_bbhcp(cctkGH,genID_cmdline_output_enable,fisheye_enable, &
       phi,kxx,kxy,kxz,kyy,kyz,shiftx,shifty,shiftz,lapm1,X(1,1,1),Y(1,1,1),Z(1,1,1),dX,dY,dZ,cctk_lsh)
  if(genID_cmdline_output_enable==1) then
     phi = 0.D0
     kxx = 0.D0
     kxy = 0.D0
     kxz = 0.D0
     kyy = 0.D0
     kyz = 0.D0
     lapm1 = 0.D0
     shiftx = 0.D0
     shifty = 0.D0
     shiftz = 0.D0
  end if


  write(*,*) "aaMIN:",X(1,1,1),Y(1,1,1),Z(1,1,1)
  write(*,*) "aaMAX:",X(cctk_lsh(1),1,1),Y(1,cctk_lsh(2),1),Z(1,1,cctk_lsh(3))

end subroutine bbhcp_read_inputfile_driver
