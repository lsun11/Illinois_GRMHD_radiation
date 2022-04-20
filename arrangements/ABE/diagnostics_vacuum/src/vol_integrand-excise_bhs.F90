#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine excise_bhs_VolInt(CCTK_ARGUMENTS)
  implicit none

  !~~~~~> Input parameters:
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !~~~~~> Other parameters:
  integer :: i,j,k
  real*8 :: bh1_distance,bh2_distance
  real*8 :: rbh1,rbh2

  rbh1 = BH_Vol_Excise_Radius
  if(num_BHs.gt.1) then
     rbh2 = BH_Vol_Excise_Radius
  else
     rbh2 = -1
  end if

  !~~~~~> Set up integration

  ! IGNORE THE FOLLOWING COMMENT.  WE DIDN'T USE THE PRIVATE() BLOCK IN THE OMP DIRECTIVE WHEN THAT COMMENT WAS WRITTEN, AND THAT'S WHAT CAUSED THE PROBLEM.
  !! Disable OpenMP for this loop, due to an OpenMP bug with the PGI compiler on Ranger.
  !!    Basically, if you uncomment the OpenMP calls around this loop, 
  !!    you'll get the wrong values for HC, MC, etc violations.
!$omp parallel do private(bh1_distance,bh2_distance)
  do k = 1,cctk_lsh(3)
     do j = 1,cctk_lsh(2)
        do i = 1,cctk_lsh(1)
           bh1_distance=sqrt((x(i,j,k)-bh_posn_x(1))**2+(y(i,j,k)-bh_posn_y(1))**2+(z(i,j,k)-bh_posn_z(1))**2)
           bh2_distance=sqrt((x(i,j,k)-bh_posn_x(2))**2+(y(i,j,k)-bh_posn_y(2))**2+(z(i,j,k)-bh_posn_z(2))**2)
           !Here we account for a common horizon:
           if((bh_posn_x(2)==0 .and. bh_posn_y(2)==0 .and. bh_posn_z(2)==0) .and. &
              (bh_posn_x(3).ne.0 .or. bh_posn_y(3).ne.0 .or. bh_posn_z(3).ne.0) ) then
              bh2_distance=sqrt((x(i,j,k)-bh_posn_x(3))**2+(y(i,j,k)-bh_posn_y(3))**2+(z(i,j,k)-bh_posn_z(3))**2)
           end if
           ! if we're inside either BH, set the volume integrand to zero
           if(bh1_distance.lt.rbh1 .or. bh2_distance.lt.rbh2) then
              VolIntegrand(i,j,k) = 0.D0
              VolIntegrand2(i,j,k) = 0.D0
              VolIntegrand3(i,j,k) = 0.D0
              VolIntegrand4(i,j,k) = 0.D0
              !              write(*,*) "HI1.",i,j,k
              !              write(*,*) "HI2.",bh1_distance,bh2_distance
              !           else
           end if
        end do      !~~> i-loop
     end do       !~~> j-loop
  end do        !~~> k-loop
!$omp end parallel do

end subroutine excise_bhs_VolInt
