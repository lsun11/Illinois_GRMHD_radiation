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

        !Set position in CarpetRegrid2
        position_x(3) = Box2X_VolInt/Box2denom_VolInt
        position_y(3) = Box2Y_VolInt/Box2denom_VolInt
        position_z(3) = 0.D0

        !Set position in CarpetRegrid2
        position_x(4) = Box3X_VolInt/Box3denom_VolInt
        position_y(4) = Box3Y_VolInt/Box3denom_VolInt
        position_z(4) = 0.D0

        !Set position in CarpetRegrid2
        position_x(5) = Box4X_VolInt/Box4denom_VolInt
        position_y(5) = Box4Y_VolInt/Box4denom_VolInt
        position_z(5) = 0.D0

        do i=2,5
           do j=2,5
              if(i.lt.j) then
                 old_distance1 = sqrt((position_x(i)-old_position_x(i))**2 + (position_x(i)-old_position_x(i))**2 + (position_x(i)-old_position_z(i))**2)
                 old_distance2 = sqrt((position_x(j)-old_position_x(i))**2 + (position_x(j)-old_position_x(i))**2 + (position_x(j)-old_position_z(i))**2)
                 if(old_distance1.gt.old_distance2) then
                    tempstoragex = position_x(i)
                    tempstoragey = position_y(i)

                    position_x(i) = position_x(j)
                    position_y(i) = position_y(j)

                    position_x(j) = tempstoragex
                    position_y(j) = tempstoragey
                 end if
              end if
           end do
        end do

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
