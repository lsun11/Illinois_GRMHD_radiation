!-----------------------------------------------------------
! Setup alternate atmosphere (typically used for disk runs)
!-----------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine disk_Setup_alt_Atmosphere(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: i,j,k

  ! testing, initialize
  rho_b_atm_gf = 0.d0
  pfloor_gf = 0.d0

  do  k=1,cctk_lsh(3)
     do  j=1,cctk_lsh(2)
        do  i=1,cctk_lsh(1)
           rho_b_atm_gf(i,j,k) = rho_atm*PhysicalRadius(i,j,k)**rho_atm_index
           pfloor_gf(i,j,k) = p_atm*PhysicalRadius(i,j,k)**p_atm_index
        end do
     end do
  end do

end subroutine disk_Setup_alt_Atmosphere
