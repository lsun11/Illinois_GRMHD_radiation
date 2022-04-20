#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!----------------------------------------------------------------------------------
! Update center of a moving box.  Useful for moving NS grid box in BHNS simulation
!----------------------------------------------------------------------------------
subroutine lw_movingbox_update_center(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  if(MOD(cctk_iteration,out_every)==0) then
     if(Symmetry==EQUATORIAL .or. Symmetry==NO_SYMM) then
        !Set position in CarpetRegrid2
        position_x(1) = 0.5d0*cos(0.5d0*CCTK_TIME)
        position_y(1) = 0.4d0*sin(0.5d0*CCTK_TIME)
        position_z(1) = 0.D0

        write(*,*) "UPDATING CENTER OF NS MOVING BOX TO:",position_x(1),position_y(1),position_z(1)
        write(*,*) "RADIUS OF MOVING BOX1:",radius(1,1)
        write(*,*) "MOVING BOX1: ACTIVE?",active(1)
     else
        write(*,*) "ERROR, UNSUPPORTED: Symmetry=",Symmetry
     end if
  end if
end subroutine lw_movingbox_update_center
