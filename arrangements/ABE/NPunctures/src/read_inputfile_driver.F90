
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine NPunctures_read_inputfile_driver(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8 :: dX,dY,dZ

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  call read_inputfile_NPunctures(cctkGH,genID_cmdline_output_enable,fisheye_enable, &
       psi,X(1,1,1),Y(1,1,1),Z(1,1,1), &
       dx,dy,dz,cctk_lsh, &
       half_binary_separation,  x_offset, npunctures_numbhs, &
       bhmass, &
       bhxpos, bhypos, &
       bh_px, bh_py, bh_pz, &
       bh_sx, bh_sy, bh_sz, &
       BigM)

  if(genID_cmdline_output_enable==1) then 
     psi = 1.D0
  else 
     write(*,*) 'ADM mass of the binary = ',BigM
  end if

end subroutine NPunctures_read_inputfile_driver
