#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!----------------------------------------------------------------------------------
! Update center of a moving box.  
!----------------------------------------------------------------------------------
subroutine movingbox_update_centers_part3(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
 
  real*8 :: xc,yc
  integer :: handle,vindex,ierr
  
!
  if (mod(cctk_iteration,out_every)==0 .and. cctk_iteration .gt. 0) then
     if (find_max_phi1==1) then 
        
        call CCTK_ReductionHandle(handle,"maximum")
        call CCTK_VarIndex(vindex,"movingbox_for_bbh::tempx1")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, xc, 1, vindex)

        call CCTK_VarIndex(vindex,"movingbox_for_bbh::tempy1")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, yc, 1, vindex)

        if ( (xc-position_x(1))**2 + (yc-position_y(1))**2 .gt. max_dev1**2) then 
	   write(*,*) 'Warning: maximum phi for BH #1 not found! Moving box center #1 is not updated.'
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
     else
	position_x(1) = sf_origin_x(1)
        position_y(1) = sf_origin_y(1)
        position_z(1) = 0.d0
        xc1_prev = sf_origin_x(1)
        yc1_prev = sf_origin_y(1)
        zc1_prev = sf_origin_z(1)
        write(*,*) 'New movingbox center #1: ',position_x(1),position_y(1),position_z(1)
     end if

     if (find_max_phi2==1) then 
        
        call CCTK_ReductionHandle(handle,"maximum")
        call CCTK_VarIndex(vindex,"movingbox_for_bbh::tempx2")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, xc, 1, vindex)

        call CCTK_VarIndex(vindex,"movingbox_for_bbh::tempy2")
        call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, yc, 1, vindex)

        if ( (xc-position_x(2))**2 + (yc-position_y(2))**2 .gt. max_dev2**2) then 
	   write(*,*) 'Warning: maximum phi for BH #2 not found! Moving box center #2 is not updated.'
	   xc = position_x(2)
           yc = position_y(2)
	end if

        position_x(2) = xc
        position_y(2) = yc
        position_z(2) = 0.d0
        write(*,*) 'Horizon not found, set the center of moving box to (',xc,yc,0,')'
        initial_guess__coord_ellipsoid__x_center(3) = xc
        initial_guess__coord_ellipsoid__y_center(3) = yc
        initial_guess__coord_ellipsoid__z_center(3) = 0.d0 
        !!initial_guess__coord_ellipsoid__x_center(3) = -1.d0
        !!initial_guess__coord_ellipsoid__y_center(3) = 0.d0
        !!initial_guess__coord_ellipsoid__z_center(3) = 0.d0
	origin_x(3) = xc
        origin_y(3) = yc
        origin_z(3) = 0.d0
        !!initial_guess__coord_ellipsoid__x_radius(3) = (sf_max_x(1)-sf_min_x(1))*0.5d0
        !!initial_guess__coord_ellipsoid__y_radius(3) = (sf_max_y(1)-sf_min_y(1))*0.5d0
        !!initial_guess__coord_ellipsoid__z_radius(3) = sf_max_z(1)
        !!initial_guess__coord_ellipsoid__x_radius(3) = 0.24705d0
        !!initial_guess__coord_ellipsoid__y_radius(3) = 0.28527d0
        !!initial_guess__coord_ellipsoid__z_radius(3) = 0.28527d0
     else
	position_x(2) = sf_origin_x(2)
        position_y(2) = sf_origin_y(2)
        position_z(2) = 0.d0
        xc2_prev = sf_origin_x(2)
        yc2_prev = sf_origin_y(2)
        zc2_prev = sf_origin_z(2)
        write(*,*) 'New movingbox center #2: ',position_x(2),position_y(2),position_z(2)
     end if


  end if

end subroutine movingbox_update_centers_part3
