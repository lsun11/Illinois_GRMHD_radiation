!--------------------------------------------------
! Compute B^i fluxes from induction equation: v1.5
!--------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine radfields_fluxes(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !OTHER PARAMETERS: 
  real*8                :: dX,dY,dZ
  integer               :: pow_axi,i,j,k
  integer, dimension(3) :: ext

  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)


!!---------------------construction---------------




  do k=1, ext(3)
     do j=1, ext(2)
        do i=1, ext(1)

if (i==27.and.j==24.and.k==16) then
   write(*,*) "1.In driver_radfields_fluxes, BEFORE cmaxcim, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==25.and.j==14.and.k==19) then
   write(*,*) "2.In driver_radfields_fluxes, BEFORE cmaxcim, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_radfields_fluxes, BEFORE cmaxcim, tau(i,j,k) is ", tau(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_radfields_fluxes, BEFORE cmaxcim, tau_rhs(i,j,k) is ", tau_rhs(i,j,k)
end if


if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_radfields_fluxes, BEFORE cmaxcim, tau_rad(i,j,k) is ", tau_rad(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_radfields_fluxes, BEFORE cmaxcim, tau_rad_rhs(i,j,k) is ", tau_rad_rhs(i,j,k)
end if


         end do
    end do
  end do



  call compute_cmax_cmin_hybrid_rad_cpp(cctkGH,ext,v02r,v02l,cmax,cmin, &
          rho_br,rho_bl, &
          Pr,Pl, vxr,vxl,vyr,vyl,vzr,vzl, &
          Bxr,Bxl,Byr,Byl,Bzr,Bzl,lapm1_f, &
          shiftx_f,shifty_f,shiftz_f, phi_f, &
          gxx_f,gxy_f,gxz_f,gyy_f,gyz_f,gzz_f, &
          gupxx_f,gupyy_f,gupzz_f,m, &
          neos,ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, &
          temp1,temp2,temp3,temp4,temp5)


  do k=1, ext(3)
     do j=1, ext(2)
        do i=1, ext(1)

if (i==27.and.j==24.and.k==16) then
   write(*,*) "1.In driver_radfields_fluxes, BEFORE flux_rad_cpp, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==25.and.j==14.and.k==19) then
   write(*,*) "2.In driver_radfields_fluxes, BEFORE flux_rad_cpp, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_radfields_fluxes, BEFORE flux_rad_cpp, tau(i,j,k) is ", tau(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_radfields_fluxes, BEFORE flux_rad_cpp, tau_rhs(i,j,k) is ", tau_rhs(i,j,k)
end if


if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_radfields_fluxes, BEFORE flux_rad_cpp, tau_rad(i,j,k) is ", tau_rad(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_radfields_fluxes, BEFORE flux_rad_cpp, tau_rad_rhs(i,j,k) is ", tau_rad_rhs(i,j,k)
end if


         end do
    end do
  end do






!!$ This routine computes the flux for tau_rad, S_radx, S_rady, S_radz


call flux_rad_cpp (m, cctkGH, ext, X, &
     		   tau_rad_flux, & 
     		   S_radx_flux, S_rady_flux, S_radz_flux, &
     		   E_radr, E_radl, &
     		   F_rad0r, F_rad0l, &
     		   F_radxr,F_radxl, &
     		   F_radyr, F_radyl, &
     		   F_radzr, F_radzl, &
     		   P_radr,P_radl, &
     		   vxr,vxl,vyr,vyl,vzr,vzl, &
     		   gxx_f, gxy_f, gxz_f, gyy_f, gyz_f, gzz_f, &
     		   cmax,cmin, &
     		   shiftx_f, shifty_f, shiftz_f, &
     		   lapm1_f, phi_f, &
     		   pow_axi,Symmetry)






  do k=1, ext(3)
     do j=1, ext(2)
        do i=1, ext(1)

if (i==27.and.j==24.and.k==16) then
   write(*,*) "1.In driver_radfields_fluxes, AFTER flux_rad_cpp, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==25.and.j==14.and.k==19) then
   write(*,*) "2.In driver_radfields_fluxes, AFTER flux_rad_cpp, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_radfields_fluxes, AFTER flux_rad_cpp, tau(i,j,k) is ", tau(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_radfields_fluxes, AFTER flux_rad_cpp, tau_rhs(i,j,k) is ", tau_rhs(i,j,k)
end if


if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_radfields_fluxes, AFTER flux_rad_cpp, tau_rad(i,j,k) is ", tau_rad(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_radfields_fluxes, AFTER flux_rad_cpp, tau_rad_rhs(i,j,k) is ", tau_rad_rhs(i,j,k)
end if


         end do
    end do
  end do





if(excision_enable == 1) then
	call scalar_excision_bc(ext,X,Y,Z,tau_rad_flux,Symmetry,excision_zone_gf);
	call scalar_excision_bc(ext,X,Y,Z,S_radx_flux,Symmetry,excision_zone_gf);
	call scalar_excision_bc(ext,X,Y,Z,S_rady_flux,Symmetry,excision_zone_gf);
	call scalar_excision_bc(ext,X,Y,Z,S_radz_flux,Symmetry,excision_zone_gf);
end if 



end subroutine radfields_fluxes
