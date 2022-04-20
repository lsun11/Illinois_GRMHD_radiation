
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bbhlorene_read_inputfile_driver(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8 :: dX,dY,dZ

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  call read_inputfile_bbhlorene(cctkGH,genID_cmdline_output_enable,fisheye_enable, &
       psi,X(1,1,1),Y(1,1,1),Z(1,1,1),dX,dY,dZ,cctk_lsh, &
       each_black_hole_mass,each_black_hole_momentum,total_binary_separation,each_black_hole_spin)

  if(genID_cmdline_output_enable==1) psi = 1.D0

  write(*,*) "aaMIN:",X(1,1,1),Y(1,1,1),Z(1,1,1)
  write(*,*) "aaMAX:",X(cctk_lsh(1),1,1),Y(1,cctk_lsh(2),1),Z(1,1,cctk_lsh(3))

end subroutine bbhlorene_read_inputfile_driver
