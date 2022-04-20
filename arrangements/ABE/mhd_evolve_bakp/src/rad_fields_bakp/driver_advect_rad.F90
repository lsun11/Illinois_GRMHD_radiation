!-------------------------------------------------------
!    :: Driver routine for MHD timestepping, v2.0 ::
! (i.e., computing RHS's of all conservative variables)
!-------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine driver_advect_rad(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  real*8                :: dX,dY,dZ,b2bt
  integer               :: index,ierr,handle,dummy
  CCTK_REAL             :: reduction_value
  integer               :: AXISYM,i,j,k
  parameter(AXISYM = 4)

  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

 write(*,*) "Start driver_advect_rad!!!!!!"

  do k=1, ext(3)
     do j=1, ext(2)
        do i=1, ext(1)

if (i==27.and.j==24.and.k==16) then
   write(*,*) "1.In driver_advect_rad, BEFORE rad_source_cpp, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==25.and.j==14.and.k==19) then
   write(*,*) "2.In driver_advect_rad, BEFORE rad_source_cpp, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_advect_rad, BEFORE rad_source_cpp, tau(i,j,k) is ", tau(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_advect_rad, BEFORE rad_source_cpp, tau_rhs(i,j,k) is ", tau_rhs(i,j,k)
end if


if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_advect_rad, BEFORE rad_source_cpp, tau_rad(i,j,k) is ", tau_rad(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_advect_rad, BEFORE rad_source_cpp, tau_rad_rhs(i,j,k) is ", tau_rad_rhs(i,j,k)
end if

	 end do 
    end do
  end do




!	write(*,*) "Start advect_rad_tau_ct_cpp!!!!!!"
        call advect_rad_tau_ct_cpp(m, cctkGH,ext,X,Y,Z, tau_rad_rhs, &
	     tau_rad_flux, Symmetry, dX, dY, dZ)


  do k=1, ext(3)
     do j=1, ext(2)
        do i=1, ext(1)

if (i==27.and.j==24.and.k==16) then
   write(*,*) "1.In driver_advect_rad, AFTER advect_rad_tau_ct_cpp, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==25.and.j==14.and.k==19) then
   write(*,*) "2.In driver_advect_rad, AFTER advect_rad_tau_ct_cpp, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_advect_rad, AFTER advect_rad_tau_ct_cpp, tau(i,j,k) is ", tau(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_advect_rad, AFTER advect_rad_tau_ct_cpp, tau_rhs(i,j,k) is ", tau_rhs(i,j,k)
end if


if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_advect_rad, AFTER advect_rad_tau_ct_cpp, tau_rad(i,j,k) is ", tau_rad(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_advect_rad, AFTER advect_rad_tau_ct_cpp, tau_rad_rhs(i,j,k) is ", tau_rad_rhs(i,j,k)
end if


    end do
    end do
end do






!        write(*,*) "Start advect_Srad_cpp!!!!!!"
	call advect_Srad_cpp(m, cctkGH, cctk_lsh, cctk_nghostzones, Symmetry, &
             S_rad_x_rhs, S_rad_y_rhs, S_rad_z_rhs, &
	     S_radx_flux, S_rady_flux, S_radz_flux, &
  	     dX,dY,dZ)


!	write(*,*) "Start rad_source_cpp!!!!!!"
        call rad_source_cpp(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry, &
       	     enable_HARM_energyvariable, &
       	     dX, dY,dZ, Z, &
       	     S_rad_x_rhs, S_rad_y_rhs, S_rad_z_rhs,  &
	     tau_rad_rhs, rho_star, &
       	     P, h, u0,  &
       	     vx, vy, vz,  &
       	     sbt, sbx, sby, sbz, &
       	     lapm1, shiftx, shifty, shiftz,  &
       	     phi,  &
       	     gxx, gxy, gxz, gyy, gyz, gzz,  &
       	     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,  &
       	     lapm1_f, shiftx_f, shifty_f, shiftz_f,  &
       	     phi_f,  &
       	     gxx_f, gxy_f, gxz_f, gyy_f, gyz_f, gzz_f, &
       	     temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8, &
       	     E_rad,F_rad0,F_radx,F_rady,F_radz,P_rad, &
       	     rad_closure_scheme,rad_opacity_abs,rad_opacity_sct, rad_const, &
	     T_fluid)

	     

     if(excision_enable == 1) then
        call remove_interior(ext,X,Y,Z,tau_rad_rhs,excision_zone_gf,Symmetry)
        call remove_interior(ext,X,Y,Z,S_rad_x_rhs,excision_zone_gf,Symmetry)
        call remove_interior(ext,X,Y,Z,S_rad_y_rhs,excision_zone_gf,Symmetry)
        call remove_interior(ext,X,Y,Z,S_rad_z_rhs,excision_zone_gf,Symmetry)      
     end if


  do k=1, ext(3)
     do j=1, ext(2)
        do i=1, ext(1)

if (i==27.and.j==24.and.k==16) then
   write(*,*) "1.In driver_advect_rad, AFTER rad_source_cpp, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==25.and.j==14.and.k==19) then
   write(*,*) "2.In driver_advect_rad, AFTER rad_source_cpp, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_advect_rad, AFTER rad_source_cpp, tau(i,j,k) is ", tau(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_advect_rad, AFTER rad_source_cpp, tau_rhs(i,j,k) is ", tau_rhs(i,j,k)
end if


if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_advect_rad, AFTER rad_source_cpp, tau_rad(i,j,k) is ", tau_rad(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_advect_rad, AFTER rad_source_cpp, tau_rad_rhs(i,j,k) is ", tau_rad_rhs(i,j,k)
end if


	end do
    end do
end do

! write(*,*) "End driver_advect_rad!!!!!!"


call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_conservatives')
call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_primitives')


end subroutine driver_advect_rad
