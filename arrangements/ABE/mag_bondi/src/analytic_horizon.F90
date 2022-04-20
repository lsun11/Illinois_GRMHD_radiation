#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine mag_bondi_analytic_horizon_data(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8 			           :: dT,dX,dY,dZ
  real*8                                   :: dphi, dcostheta, m0
  real*8                                   :: w_eh, sym_factor
  real*8                                   :: costheta, sintheta, ph
  real*8                                   :: rhh, rhh_grid, dd
  integer                                  :: n, i,j, k, imin,imax,jmin,jmax,kmin,kmax
  real*8                                   :: xmin,ymin,ymax_bondi,zmin
  real*8                                   :: grid_rmin, grid_rmax, temp
  real*8                                   :: x_l, y_l, z_l, r_grid, r_phys
  integer :: handle,dummy,proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax
  integer :: proc_kmax,glob_imax,glob_jmax,glob_kmax, index
  integer :: ierr,ONE,ZERO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  integer, parameter :: AXISYM_FULL = 5
  real*8  :: hdist_ana_bl
  real*8, parameter :: PI = 3.14159265358979323846d0
  external hdist_ana_bl
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  CCTK_REAL reduction_value

  dphi = 2.d0 * PI / N_phi
  dcostheta = 2.d0 / N_theta
  m0 = 1.d0

  ! determine w for the event horizon
  w_eh = m0 + sqrt(m0*m0 - sam*sam)
  
  if (Symmetry==OCTANT) then
     sym_factor = 8.d0
     dcostheta = 1.d0/N_theta
     dphi = 0.5d0*PI / N_phi
  else if (Symmetry==EQUATORIAL) then
     sym_factor = 2.d0
      dcostheta = 1.d0/N_theta
  else if (Symmetry==NO_SYMM) then
     sym_factor = 1.d0
  else if (Symmetry==PI_SYMM) then
     sym_factor = 4.d0
  else if (Symmetry==AXISYM .or. Symmetry==AXISYM_FULL) then
     sym_factor = 1.d0
     if (N_phi .ne. 1) then 
	write(*,*) 'N_phi = ',N_phi 
	write(*,*) 'N_phi must be equal to 1 in mag_bh_accretion_analytic_horizon_data when '
	write(*,*) 'evolved in axisymmetry.'
	stop
     end if
     if (Symmetry==AXISYM) dcostheta = 1.d0/N_theta
  end if

  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"grid::X")
  ! Note: xmax has already been computed
  !!call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xmax,1,index)
  call CCTK_VarIndex(index,"grid::Y")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,ymax_bondi,1,index)
  call CCTK_VarIndex(index,"grid::Z")
  ! Note: zmax has already been computed
  !!call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,zmax,1,index)
  grid_rmin = 0.0
  grid_rmax = sqrt(xmax_bondi*xmax_bondi + ymax_bondi*ymax_bondi + zmax_bondi*zmax_bondi)

  ! fill list of unit normals

  n = 1
  do i = 1, N_theta
    costheta = 1.d0 - (i - 0.5d0)*dcostheta
    sintheta = sqrt(1.d0 - costheta*costheta)
    do j = 1, N_phi
       if(N_phi==1) then
          ph = 0.d0
       else
          ph = (j - 0.5d0)*dphi
       end if

       nx_d(n) = sintheta*cos(ph)
       ny_d(n) = sintheta*sin(ph)
       nz_d(n) = costheta

       rhh = w_eh

       ! radially transform the horizon radius:
       if (puncture_id==1) then 
	  rhh_grid = 0.5d0
       else
          rhh_grid = rhh - r0
       end if

       !if (CCTK_MyProc(CCTKGH) == 0) write(*,*) "i, rhh_grid = ", i ," " , rhh_grid 

       ! compute point location
       xn_d(n) = rhh_grid * nx_d(n)
       yn_d(n) = rhh_grid * ny_d(n)
       zn_d(n) = rhh_grid * nz_d(n)

       ! compute point normal:  n = grad(r)
       !!dd = 0.01d0*rhh_grid
       !!nx_d(n) = (hdist_ana_bl(xn_d(n)+dd,yn_d(n),zn_d(n),m0,sam,r0)- &
       !!     hdist_ana_bl(xn_d(n)-dd,yn_d(n),zn_d(n),m0,sam,r0)) / (2.0*dd)
       !!ny_d(n) = (hdist_ana_bl(xn_d(n),yn_d(n)+dd,zn_d(n),m0,sam,r0)- &
       !!     hdist_ana_bl(xn_d(n),yn_d(n)-dd,zn_d(n),m0,sam,r0)) / (2.0*dd) 
       !!nz_d(n) = (hdist_ana_bl(xn_d(n),yn_d(n),zn_d(n)+dd,m0,sam,r0)- &
       !!     hdist_ana_bl(xn_d(n),yn_d(n),zn_d(n)-dd,m0,sam,r0)) / (2.0*dd)
       !!
       !!nn_d(n) = (nx_d(n)*cos(ph) + ny_d(n)*sin(ph))*sintheta + nz_d(n)*costheta;
       ! normalize with the "Jacobian."
       !!nx_d(n) = nx_d(n)/nn_d(n); ny_d(n) = ny_d(n)/nn_d(n);  nz_d(n) = nz_d(n)/nn_d(n);

       n = n + 1
    end do
 end do

! Convert some other radii:
 if (puncture_id==1) then 
    r_out_flux1 = 0.5d0*(r_out_flux1-1.d0 + sqrt(r_out_flux1*(r_out_flux1-2.d0)) )
    r_out_flux2 = 0.5d0*(r_out_flux2-1.d0 + sqrt(r_out_flux2*(r_out_flux2-2.d0)) )
    r_out_flux3 = 0.5d0*(r_out_flux3-1.d0 + sqrt(r_out_flux3*(r_out_flux3-2.d0)) )
 else
    r_out_flux1 = r_out_flux1 - r0
    r_out_flux2 = r_out_flux2 - r0
    r_out_flux3 = r_out_flux3 - r0
 end if

 !!if (excision_radius .gt. 0.0) then
 !!   call Convert_Physical_Radius(ierr,0.8*excision_radius,grid_rmin,grid_rmax,fisheye_enable,temp)
 !!   excision_radius  = temp/0.8
 !!end if
 !!call Convert_Physical_Radius(ierr,r_out_flux1,grid_rmin,grid_rmax,fisheye_enable,temp)
 !!r_out_flux1 = temp
 !!call Convert_Physical_Radius(ierr,r_out_flux2,grid_rmin,grid_rmax,fisheye_enable,temp)
 !!r_out_flux2 = temp
 !!call Convert_Physical_Radius(ierr,r_out_flux3,grid_rmin,grid_rmax,fisheye_enable,temp)
 !!r_out_flux3 = temp

end subroutine mag_bondi_analytic_horizon_data

subroutine mag_bondi_some_radius_conversions(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext
  real*8                                   :: dT,dX,dY,dZ
  real*8                                   :: dphi, dcostheta, m0
  real*8                                   :: w_eh, sym_factor
  real*8                                   :: costheta, sintheta, ph
  real*8                                   :: rhh, rhh_grid, dd
  integer                                  :: n, i,j, k, imin,imax,jmin,jmax,kmin,kmax
  real*8                                   :: xmin,ymin,ymax,zmin
  real*8                                   :: grid_rmin, grid_rmax, temp
  real*8                                   :: x_l, y_l, z_l, r_grid, r_phys
  real*8                                   :: PI
  integer :: handle,dummy,proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax
  integer :: proc_kmax,glob_imax,glob_jmax,glob_kmax, index
  integer :: ierr,ONE,ZERO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  integer, parameter :: AXISYM_FULL = 5
  real*8  :: find_hdist_ana, find_hdist_ana_bl
  external find_hdist_ana, find_hdist_ana_bl
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  CCTK_REAL reduction_value

  !!call CCTK_ReductionHandle(handle,"maximum")
  !!call CCTK_VarIndex(index,"grid::X")
  !!call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xmax,1,index)
  !!call CCTK_VarIndex(index,"grid::Y")
  !!call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,ymax,1,index)
  !!call CCTK_VarIndex(index,"grid::Z")
  !!call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,zmax,1,index)
  !!grid_rmin = 0.0
  !!grid_rmax = sqrt(xmax*xmax + ymax*ymax + zmax*zmax)

  !!if (excision_radius .gt. 0.0) then
  !!   call Convert_Physical_Radius(ierr,0.8*excision_radius,grid_rmin,grid_rmax,fisheye_enable,temp)
  !!   excision_radius  = temp/0.8
  !!end if
  !!call Convert_Physical_Radius(ierr,r_out_flux1,grid_rmin,grid_rmax,fisheye_enable,temp)
  !!r_out_flux1 = temp
  !!call Convert_Physical_Radius(ierr,r_out_flux2,grid_rmin,grid_rmax,fisheye_enable,temp)
  !!r_out_flux2 = temp
  !!call Convert_Physical_Radius(ierr,r_out_flux3,grid_rmin,grid_rmax,fisheye_enable,temp)
  !!r_out_flux3 = temp

  if (excision_radius .gt. 0.0) excision_radius = (0.8d0*excision_radius -r0)/0.8d0
  r_out_flux1 = r_out_flux1 - r0
  r_out_flux2 = r_out_flux2 - r0
  r_out_flux3 = r_out_flux3 - r0

end subroutine mag_bondi_some_radius_conversions

FUNCTION hdist_ana_bl(xx,yy,zz,m0,sam,r0)
  implicit none
  DECLARE_CCTK_FUNCTIONS

  real*8 :: hdist_ana_bl, xx, yy, zz, sam
  real*8 :: m0, rad, r_eh_grid,r0

  r_eh_grid = m0 + sqrt(m0*m0 - sam*sam) - r0
  rad = sqrt(xx*xx+yy*yy+zz*zz)

  hdist_ana_bl = rad - r_eh_grid
end FUNCTION hdist_ana_bl
