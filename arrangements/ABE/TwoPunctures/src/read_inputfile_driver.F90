
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine TwoPunctures_read_inputfile_driver(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8 :: dX,dY,dZ

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  call read_inputfile_TwoPunctures(cctkGH,genID_cmdline_output_enable,fisheye_enable, &
       psi,X(1,1,1),Y(1,1,1),Z(1,1,1),dX,dY,dZ,cctk_lsh, &
       bh_mass_plus,bh_mass_minus,bh_px_plus,bh_px_minus,bh_py_plus,bh_py_minus, & 
       bh_spin_plus,bh_spin_minus, half_binary_separation, x_offset, BigM)

  if(genID_cmdline_output_enable==1) then 
    psi = 1.D0
  else 
    write(*,*) 'ADM mass of the binary = ',BigM
  end if

end subroutine TwoPunctures_read_inputfile_driver
