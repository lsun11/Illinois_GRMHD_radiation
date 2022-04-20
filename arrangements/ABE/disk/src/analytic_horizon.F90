#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine analytic_horizon_data(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext
  real*8 			           :: dT,dX,dY,dZ
  real*8                                   :: dphi, dcostheta, m0
  real*8                                   :: w_eh, sym_factor
  real*8                                   :: costheta, sintheta, ph
  real*8                                   :: rhh, rhh_grid, dd
  integer                                  :: n, i,j, k, imin,imax,jmin,jmax,kmin,kmax
  real*8                                   :: xmin,xmax,ymin,ymax,zmin,zmax
  real*8                                   :: grid_rmin, grid_rmax, temp
  real*8                                   :: x_l, y_l, z_l, r_grid, r_phys
  real*8                                   :: PI
  integer :: handle,dummy,proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax
  integer :: proc_kmax,glob_imax,glob_jmax,glob_kmax, index
  integer :: ierr,ONE,ZERO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  real*8  :: find_hdist_ana, find_hdist_ana_bl
  external find_hdist_ana, find_hdist_ana_bl
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  CCTK_REAL reduction_value

  !ZACH SAYS: THERE IS A BUG HERE:
  !PI = acos(1.d0)
  !Fixed line of code:
  PI = acos(-1.d0)
  dphi = 2.d0 * PI / N_phi
  dcostheta = 2.d0 / N_theta
  m0 = 1.d0
  
  ! determine w for the event horizon
  w_eh = m0 + sqrt(m0*m0 - sam*sam)
  
  if (Symmetry==OCTANT) then
     sym_factor = 8.d0
  else if (Symmetry==EQUATORIAL) then
     sym_factor = 2.d0
  else if (Symmetry==NO_SYMM) then
     sym_factor = 1.d0
  else if (Symmetry==PI_SYMM) then
     sym_factor = 4.d0
  else if (Symmetry==AXISYM) then
     sym_factor = 1.d0
     dcostheta = 1.0/N_theta
  end if

  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"grid::X")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xmax,1,index)
  call CCTK_VarIndex(index,"grid::Y")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,ymax,1,index)
  call CCTK_VarIndex(index,"grid::Z")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,zmax,1,index)
  grid_rmin = 0.0
  grid_rmax = sqrt(xmax*xmax + ymax*ymax + zmax*zmax)

  ! fill list of unit normals

  n = 1
  do i = 1, N_theta
    costheta = 1.d0 - (i - 0.5d0)*dcostheta
    sintheta = sqrt(1.d0 - costheta*costheta)
    do j = 1, N_phi
       if(N_phi==1) then
          ph = 0.d0
       else
          ph = (j + 0.5d0)*dphi
       end if

       nx_d(n) = sintheta*cos(ph)
       ny_d(n) = sintheta*sin(ph)
       nz_d(n) = costheta

       if (radial_coordinate_type == 0) then
          rhh = sintheta*sintheta/(w_eh*w_eh + sam*sam)
          rhh = rhh + costheta*costheta/(w_eh*w_eh)
          rhh = 1.0/sqrt(rhh)
       else
          rhh = w_eh
       end if

       ! radially transform the horizon radius:
       call Convert_Physical_Radius(ierr,rhh,grid_rmin,grid_rmax,fisheye_enable,rhh_grid)

!       if (CCTK_MyProc(CCTKGH) == 0) write(*,*) "i, rhh_grid = ", i ," " , rhh_grid 

       ! compute point location
       xn_d(n) = rhh_grid * nx_d(n)
       yn_d(n) = rhh_grid * ny_d(n)
       zn_d(n) = rhh_grid * nz_d(n)

       ! compute point normal:  n = grad(r)
       dd = 0.01*rhh_grid
       if (radial_coordinate_type == 0) then
          nx_d(n) = (find_hdist_ana(xn_d(n)+dd,yn_d(n),zn_d(n),sam,grid_rmin,grid_rmax,fisheye_enable)- & 
               find_hdist_ana(xn_d(n)-dd,yn_d(n),zn_d(n),sam,grid_rmin,grid_rmax,fisheye_enable)) / (2.0*dd)
          ny_d(n) = (find_hdist_ana(xn_d(n),yn_d(n)+dd,zn_d(n),sam,grid_rmin,grid_rmax,fisheye_enable)- &
               find_hdist_ana(xn_d(n),yn_d(n)-dd,zn_d(n),sam,grid_rmin,grid_rmax,fisheye_enable)) / (2.0*dd) 
          nz_d(n) = (find_hdist_ana(xn_d(n),yn_d(n),zn_d(n)+dd,sam,grid_rmin,grid_rmax,fisheye_enable)- &
               find_hdist_ana(xn_d(n),yn_d(n),zn_d(n)-dd,sam,grid_rmin,grid_rmax,fisheye_enable)) / (2.0*dd)
       else
          nx_d(n) = (find_hdist_ana_bl(xn_d(n)+dd,yn_d(n),zn_d(n),sam,grid_rmin,grid_rmax,fisheye_enable)- &
               find_hdist_ana_bl(xn_d(n)-dd,yn_d(n),zn_d(n),sam,grid_rmin,grid_rmax,fisheye_enable)) / (2.0*dd)
          ny_d(n) = (find_hdist_ana_bl(xn_d(n),yn_d(n)+dd,zn_d(n),sam,grid_rmin,grid_rmax,fisheye_enable)- &
               find_hdist_ana_bl(xn_d(n),yn_d(n)-dd,zn_d(n),sam,grid_rmin,grid_rmax,fisheye_enable)) / (2.0*dd) 
          nz_d(n) = (find_hdist_ana_bl(xn_d(n),yn_d(n),zn_d(n)+dd,sam,grid_rmin,grid_rmax,fisheye_enable)- &
               find_hdist_ana_bl(xn_d(n),yn_d(n),zn_d(n)-dd,sam,grid_rmin,grid_rmax,fisheye_enable)) / (2.0*dd)
       end if
       
       nn_d(n) = (nx_d(n)*cos(ph) + ny_d(n)*sin(ph))*sintheta + nz_d(n)*costheta;
       ! normalize with the "Jacobian."
       nx_d(n) = nx_d(n)/nn_d(n);  ny_d(n) = ny_d(n)/nn_d(n);  nz_d(n) = nz_d(n)/nn_d(n);
       
       n = n + 1
    end do
 end do

 ! Finally, you need to calculate the event horizon mask.
 imin = lbound(X,1)
 imax = ubound(X,1)
 jmin = lbound(X,2)
 jmax = ubound(X,2)
 kmin = lbound(X,3)
 kmax = ubound(X,3)

 mskf = 0.d0

 do i = imin+1, imax
    do j = jmin, jmax
       do k = kmin+1, kmax
          x_l = X(i,j,k)
          y_l = Y(i,j,k)
          z_l = Z(i,j,k)
          r_grid = sqrt(x_l**2 + y_l**2 + z_l**2)
          costheta = z_l/r_grid
          sintheta = sqrt(1.d0-costheta**2)

          if (radial_coordinate_type == 0) then
             rhh = sintheta*sintheta/(w_eh*w_eh + sam*sam)
             rhh = rhh + costheta*costheta/(w_eh*w_eh)
             rhh = 1.0/sqrt(rhh)
          else
             rhh = w_eh
          end if
          
          r_phys = PhysicalRadius(i,j,k)
          if (r_phys .lt. rhh) then 
             mskf(i,j,k) = 0.d0
          else
             mskf(i,j,k) = 1.d0
          end if
       end do
    end do
 end do

! Convert some other radii:
 if (excision_radius .gt. 0.0) then
    call Convert_Physical_Radius(ierr,0.8*excision_radius,grid_rmin,grid_rmax,fisheye_enable,temp)
    excision_radius  = temp/0.8
 end if
 call Convert_Physical_Radius(ierr,r_out_flux1,grid_rmin,grid_rmax,fisheye_enable,temp)
 r_out_flux1 = temp
 call Convert_Physical_Radius(ierr,r_out_flux2,grid_rmin,grid_rmax,fisheye_enable,temp)
 r_out_flux2 = temp
 call Convert_Physical_Radius(ierr,r_out_flux3,grid_rmin,grid_rmax,fisheye_enable,temp)
 r_out_flux3 = temp
end subroutine analytic_horizon_data

FUNCTION find_hdist_ana(xx,yy,zz,sam,grid_rmin,grid_rmax,fisheye_enable)
  implicit none
  DECLARE_CCTK_FUNCTIONS

  real*8 :: find_hdist_ana, xx, yy, zz, sam
  real*8 :: grid_rmin, grid_rmax
  real*8 :: m0, w_eh, rad, r_eh, r_eh_grid
  real*8 :: sin2theta, cos2theta
  integer :: fisheye_enable, ierr

  w_eh = m0 + sqrt(m0*m0 - sam*sam)
  rad = sqrt(xx*xx+yy*yy+zz*zz)
  cos2theta = zz*zz/(rad*rad)
  sin2theta = 1.0-cos2theta 
  r_eh = sin2theta/(w_eh*w_eh + sam*sam)
  r_eh = r_eh + cos2theta/(w_eh*w_eh)
  r_eh = 1.0/sqrt(r_eh)

  call Convert_Physical_Radius(ierr,r_eh,grid_rmin,grid_rmax,fisheye_enable,r_eh_grid)

  find_hdist_ana = rad - r_eh_grid

end FUNCTION find_hdist_ana

FUNCTION find_hdist_ana_bl(xx,yy,zz,sam,grid_rmin,grid_rmax,fisheye_enable)
  implicit none
  DECLARE_CCTK_FUNCTIONS

  real*8 :: find_hdist_ana_bl, xx, yy, zz, sam
  real*8 :: grid_rmin,grid_rmax
  real*8 :: m0, w_eh, rad, r_eh, r_eh_grid
  real*8 :: sin2theta, cos2theta
  integer :: fisheye_enable,ierr

  w_eh = m0 + sqrt(m0*m0 - sam*sam)
  rad = sqrt(xx*xx+yy*yy+zz*zz)
  cos2theta = zz*zz/(rad*rad)
  sin2theta = 1.0-cos2theta 
  r_eh = w_eh
  call Convert_Physical_Radius(ierr,r_eh,grid_rmin,grid_rmax,fisheye_enable,r_eh_grid)

  find_hdist_ana_bl = rad - r_eh_grid
end FUNCTION find_hdist_ana_bl

subroutine some_radius_conversions(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext
  real*8 			           :: dT,dX,dY,dZ
  real*8                                   :: dphi, dcostheta, m0
  real*8                                   :: w_eh, sym_factor
  real*8                                   :: costheta, sintheta, ph
  real*8                                   :: rhh, rhh_grid, dd
  integer                                  :: n, i,j, k, imin,imax,jmin,jmax,kmin,kmax
  real*8                                   :: xmin,xmax,ymin,ymax,zmin,zmax
  real*8                                   :: grid_rmin, grid_rmax, temp
  real*8                                   :: x_l, y_l, z_l, r_grid, r_phys
  real*8                                   :: PI
  integer :: handle,dummy,proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax
  integer :: proc_kmax,glob_imax,glob_jmax,glob_kmax, index
  integer :: ierr,ONE,ZERO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  real*8  :: find_hdist_ana, find_hdist_ana_bl
  external find_hdist_ana, find_hdist_ana_bl
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  CCTK_REAL reduction_value

  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"grid::X")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xmax,1,index)
  call CCTK_VarIndex(index,"grid::Y")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,ymax,1,index)
  call CCTK_VarIndex(index,"grid::Z")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,zmax,1,index)
  grid_rmin = 0.0
  grid_rmax = sqrt(xmax*xmax + ymax*ymax + zmax*zmax)

  if (excision_radius .gt. 0.0) then
     call Convert_Physical_Radius(ierr,0.8*excision_radius,grid_rmin,grid_rmax,fisheye_enable,temp)
     excision_radius  = temp/0.8
  end if
  call Convert_Physical_Radius(ierr,r_out_flux1,grid_rmin,grid_rmax,fisheye_enable,temp)
  r_out_flux1 = temp
  call Convert_Physical_Radius(ierr,r_out_flux2,grid_rmin,grid_rmax,fisheye_enable,temp)
  r_out_flux2 = temp
  call Convert_Physical_Radius(ierr,r_out_flux3,grid_rmin,grid_rmax,fisheye_enable,temp)
  r_out_flux3 = temp
  
end subroutine some_radius_conversions
