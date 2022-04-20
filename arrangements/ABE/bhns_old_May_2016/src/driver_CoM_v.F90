#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine find_CoM_velocity(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS


  real*8                                   :: x_NS_CoM_coord,y_NS_CoM_coord
  real*8                                   :: x_NS_CoM_coord2,y_NS_CoM_coord2
  real*8,dimension(1,3)                    :: points_CM
  real*8,dimension(1,3)                    :: points_CM2
  integer                                  :: vindex

!  bhnsinsertBNOW=0

  if(CCTK_ITERATION .eq. ITERATION_TO_INSERT_MAGNETIC_FIELDS .or. (ITERATION_TO_INSERT_MAGNETIC_FIELDS.lt.0 .and. CCTK_ITERATION .eq. 0) ) then

!     bhnsinsertBNOW=1
    
     if (CCTK_ITERATION > 0) then
        if((use_new_bhns_initial_data).ne.3) then
           x_NS_CoM_coord = CoMx_VolInt/CoM_VolInt_denominator !initial_ns_coord_x                                                     
           y_NS_CoM_coord = CoMy_VolInt/CoM_VolInt_denominator !initial_ns_coord_y
        else 
           !Set position in CarpetRegrid2
           x_NS_CoM_coord = Box1X_VolInt1/Box1denom_VolInt1
           y_NS_CoM_coord = Box1Y_VolInt1/Box1denom_VolInt1
           
           !Set position in CarpetRegrid2
           x_NS_CoM_coord2 = Box1X_VolInt2/Box1denom_VolInt2
           y_NS_CoM_coord2 = Box1Y_VolInt2/Box1denom_VolInt2
        end if
     else 
        if((use_new_bhns_initial_data).ne.3) then
           x_NS_CoM_coord = initial_ns_coord_x
           y_NS_CoM_coord = initial_ns_coord_y
        else
           x_NS_CoM_coord  =  initial_ns_coord_x
           y_NS_CoM_coord  =  initial_ns_coord_y
           x_NS_CoM_coord2 =  initial_ns2_coord_x
           y_NS_CoM_coord2 =  initial_ns2_coord_y
        end if
     end if
        
!    Compute the  velocity of the center of mass                         
                                                                                                            
     points_CM(1,1) = x_NS_CoM_coord
     points_CM(1,2) = y_NS_CoM_coord
     points_CM(1,3) = 0.0d0

     call CCTK_VarIndex(vindex,"mhd_evolve::vx")
     call interp_driver_carp(cctkGH,1,points_CM,vindex,bhns_vx_CM)
     call CCTK_VarIndex(vindex,"mhd_evolve::vy")
     call interp_driver_carp(cctkGH,1,points_CM,vindex,bhns_vy_CM)
     call CCTK_VarIndex(vindex,"mhd_evolve::vz")
     call interp_driver_carp(cctkGH,1,points_CM,vindex,bhns_vz_CM)
     
     print *, "reset the velocity of the atmosphere (vx,vy,vz)=",bhns_vx_CM, bhns_vy_CM, bhns_vz_CM
     print *, " (x_NS,y_NS)=",x_NS_CoM_coord,y_NS_CoM_coord

     if (use_new_bhns_initial_data.eq.3 ) then
        points_CM2(1,1) = x_NS_CoM_coord2
        points_CM2(1,2) = y_NS_CoM_coord2
        points_CM2(1,3) = 0.0d0
        
        call CCTK_VarIndex(vindex,"mhd_evolve::vx")
        call interp_driver_carp(cctkGH,1,points_CM2,vindex,bhns_vx_CM2)
        call CCTK_VarIndex(vindex,"mhd_evolve::vy")
        call interp_driver_carp(cctkGH,1,points_CM2,vindex,bhns_vy_CM2)
        call CCTK_VarIndex(vindex,"mhd_evolve::vz")
        call interp_driver_carp(cctkGH,1,points_CM2,vindex,bhns_vz_CM2)
        
        print *, "reset the velocity of the atmosphere (vx2,vy2,vz2)=",bhns_vx_CM2, bhns_vy_CM2, bhns_vz_CM2
        print *, " (x_NS2,y_NS2)=",x_NS_CoM_coord2,y_NS_CoM_coord2
     end if

  end if
end subroutine find_CoM_velocity
