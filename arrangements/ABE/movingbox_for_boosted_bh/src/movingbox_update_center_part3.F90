#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!----------------------------------------------------------------------------------
! Update center of a moving box.  
!----------------------------------------------------------------------------------
subroutine movingbox_update_center_part3(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
 
  real*8 :: num_max, xc,yc
  integer :: handle,vindex,ierr
  
!
  if (mod(cctk_iteration,out_every)==0 .and. cctk_iteration .gt. 0) then
     if (find_max_phi==1) then 
        
        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_VarIndex(vindex,"movingbox_for_boosted_bh::temp0")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, num_max, 1, vindex)

        call CCTK_VarIndex(vindex,"movingbox_for_boosted_bh::tempx")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, xc, 1, vindex)

        call CCTK_VarIndex(vindex,"movingbox_for_boosted_bh::tempy")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, yc, 1, vindex)

        if (num_max == 0.d0) then
	   write(*,*) 'Warning: maximum phi not found! Moving box center is not updated.'
	   xc = position_x(1)
	   yc = position_y(1)
        else
           xc = xc/num_max
           yc = yc/num_max
        end if

        if (abs(xc-position_x(1)) .gt. max_dev .or. abs(yc-position_y(1)) .gt. max_dev) then 
           write(*,*) 'Warning: the maximum phi is located at (x,y,z) = ',xc,yc,0.d0
	   write(*,*) 'This location deviates a lot from the previous moving box center: ',position_x(1),position_y(1),0.d0
	   write(*,*) 'Moving box center is not updated.'
	   xc = position_x(1)
           yc = position_y(1)
	end if

        position_x(1) = xc
        position_y(1) = yc
        position_z(1) = 0.d0
        write(*,*) 'Horizon not found, set the center of moving box to (',xc,yc,0,')'
        initial_guess__coord_ellipsoid__x_center(2) = xc
        initial_guess__coord_ellipsoid__y_center(2) = yc
        initial_guess__coord_ellipsoid__z_center(2) = 0.d0 
        !!initial_guess__coord_ellipsoid__x_center(2) = -1.d0
        !!initial_guess__coord_ellipsoid__y_center(2) = 0.d0
        !!initial_guess__coord_ellipsoid__z_center(2) = 0.d0
	origin_x(2) = xc
        origin_y(2) = yc
        origin_z(2) = 0.d0
        !!initial_guess__coord_ellipsoid__x_radius(2) = (sf_max_x(1)-sf_min_x(1))*0.5d0
        !!initial_guess__coord_ellipsoid__y_radius(2) = (sf_max_y(1)-sf_min_y(1))*0.5d0
        !!initial_guess__coord_ellipsoid__z_radius(2) = sf_max_z(1)
        !!initial_guess__coord_ellipsoid__x_radius(2) = 0.24705d0
        !!initial_guess__coord_ellipsoid__y_radius(2) = 0.28527d0
        !!initial_guess__coord_ellipsoid__z_radius(2) = 0.28527d0
        !!ah_origin_x = xc
	!!ah_origin_y = yc
        !!ah_origin_z = 0.d0
        !!ah_centroid_x = xc
        !!ah_centroid_y = yc
        !!ah_centroid_z = 0.d0
        !!ah_centroid_x_p = xc
        !!ah_centroid_y_p = yc
        !!ah_centroid_z_p = 0.d0
        !!ah_radius = 0.28527d0
        !!ah_found_flag = 1
	!!ah_centroid_valid = 1
        !!sf_mean_radius = 0.28527d0
        !!sf_centroid_x = xc
        !!sf_centroid_y = yc
        !!sf_centroid_z = 0.d0
        !!sf_origin_x = xc
        !!sf_origin_y = yc
        !!sf_origin_z = 0.d0
        !!sf_radius = 0.28527d0
        !!sf_min_x(1) = xc - 0.28527d0
        !!sf_min_y(1) = yc - 0.28527d0
        !!sf_min_z(1) = - 0.28527d0
        !!sf_max_x(1) = xc + 0.28527d0
        !!sf_max_y(1) = yc + 0.28527d0
        !!sf_max_z(1) = 0.28527d0
     else
	position_x(1) = sf_origin_x(1)
        position_y(1) = sf_origin_y(1)
        position_z(1) = 0.d0
        xc_prev = sf_origin_x(1)
        yc_prev = sf_origin_y(1)
        zc_prev = sf_origin_z(1)
        write(*,*) 'New movingbox center: ',position_x(1),position_y(1),position_z(1)
     end if

  end if

end subroutine movingbox_update_center_part3
