#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!----------------------------------------------------------------------------------
! Update center of a moving box.  Useful for moving NS grid box in BHNS simulation
!----------------------------------------------------------------------------------
subroutine movingbox_update_centers(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: i,j,closest_box_to_bh
  real*8,dimension(3) :: old_position_x,old_position_y,old_position_z
  real*8  :: tempstoragex,tempstoragey
  real*8  :: bh_distance1,bh_distance2
  real*8  :: old_distance1,old_distance2
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  if(MOD(cctk_iteration,out_every)==0) then
     if(Symmetry==EQUATORIAL) then

        old_position_x = position_x
        old_position_y = position_y
        old_position_z = position_z

        !Set position in CarpetRegrid2
        position_x(2) = Box1X_VolInt/Box1denom_VolInt
        position_y(2) = Box1Y_VolInt/Box1denom_VolInt
        position_z(2) = 0.D0

        if(sqrt((position_x(2)-position_x(1))**2 + (position_y(2)-position_y(1))**2 + (position_z(2)-position_z(1))**2) .gt. 3.0*Finest_Refinement_Box_Radius) then
           do i=3,10
              position_x(i) = position_x(2)
              position_y(i) = position_y(2)
              position_z(i) = position_z(2)
           end do
        else
           !Set position in CarpetRegrid2
           position_x(3) = position_x(2)+2.0*Finest_Refinement_Box_Radius
           position_y(3) = position_y(2)
           position_z(3) = 0.D0

           !Set position in CarpetRegrid2
           position_x(4) = position_x(2)+2.0*Finest_Refinement_Box_Radius
           position_y(4) = position_y(2)+2.0*Finest_Refinement_Box_Radius
           position_z(4) = 0.D0

           !Set position in CarpetRegrid2
           position_x(5) = position_x(2)
           position_y(5) = position_y(2)+2.0*Finest_Refinement_Box_Radius
           position_z(5) = 0.D0

           Box2X_VolInt = position_x(3)
           Box2Y_VolInt = position_y(3)
           Box3X_VolInt = position_x(4)
           Box3Y_VolInt = position_y(4)
           Box4X_VolInt = position_x(5)
           Box4Y_VolInt = position_y(5)
           Box2denom_VolInt = 1.D0
           Box3denom_VolInt = 1.D0
           Box4denom_VolInt = 1.D0

           !Set position in CarpetRegrid2
           position_x(6) = position_x(2)-2.0*Finest_Refinement_Box_Radius
           position_y(6) = position_y(2)+2.0*Finest_Refinement_Box_Radius
           position_z(6) = 0.D0

           !Set position in CarpetRegrid2
           position_x(7) = position_x(2)-2.0*Finest_Refinement_Box_Radius
           position_y(7) = position_y(2)
           position_z(7) = 0.D0

           !Set position in CarpetRegrid2
           position_x(8) = position_x(2)-2.0*Finest_Refinement_Box_Radius
           position_y(8) = position_y(2)-2.0*Finest_Refinement_Box_Radius
           position_z(8) = 0.D0

           !Set position in CarpetRegrid2
           position_x(9) = position_x(2)
           position_y(9) = position_y(2)-2.0*Finest_Refinement_Box_Radius
           position_z(9) = 0.D0

           !Set position in CarpetRegrid2
           position_x(10) = position_x(2)+2.0*Finest_Refinement_Box_Radius
           position_y(10) = position_y(2)-2.0*Finest_Refinement_Box_Radius
           position_z(10) = 0.D0

        end if
        write(*,*) "UPDATING CENTER OF NS MOVING BOX TO:",position_x(2),position_y(2),position_z(2)
        write(*,*) "UPDATING CENTER OF 2nd MOVING BOX TO:",position_x(3),position_y(3),position_z(3)
        write(*,*) "UPDATING CENTER OF 3rd MOVING BOX TO:",position_x(4),position_y(4),position_z(4)
        write(*,*) "UPDATING CENTER OF 4th MOVING BOX TO:",position_x(5),position_y(5),position_z(5)
        write(*,*) "RADIUS OF MOVING BOX1:",radius(1,1)
        write(*,*) "RADIUS OF MOVING BOX2:",radius(2,1)
        write(*,*) "MOVING BOX1: ACTIVE (BH)?",active(1)
        write(*,*) "MOVING BOX2: ACTIVE (NS)?",active(2)
        write(*,*) "MOVING BOX3: ACTIVE?",active(3)
        write(*,*) "MOVING BOX4: ACTIVE?",active(4)
        write(*,*) "MOVING BOX5: ACTIVE?",active(5)
     else
        write(*,*) "ERROR, UNSUPPORTED: Symmetry=",Symmetry
     end if
  end if
end subroutine movingbox_update_centers
