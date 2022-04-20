!-------------------------------------------------------
!    :: Driver routine for MHD timestepping,
! (i.e., computing RHS's of all conservative variables)
!-------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine mhd_check_tau(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  real*8                :: dX,dY,dZ
  integer               :: index,ierr,handle,dummy
  CCTK_REAL             :: reduction_value
  integer               :: AXISYM,i,j,k
  parameter(AXISYM = 4)


ext = cctk_lsh

  do k=1, ext(3)
     do j=1, ext(2)
        do i=1, ext(1)

if (i==27.and.j==24.and.k==16) then
   write(*,*) "1.In Driver check tau, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==25.and.j==14.and.k==19) then
   write(*,*) "2.In Driver check tau, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_check_tau, tau(i,j,k) is ", tau(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_check_tau, tau_rhs(i,j,k) is ", tau_rhs(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_check_tau, tau_rad(i,j,k) is ", tau_rad(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
    write(*,*) " In driver_check_tau, tau_rad_rhs(i,j,k) is ", tau_rad_rhs(i,j,k)
end if

	end do
    end do
  end do




if(1==0) then
do k=1,3
write(*,*) "Inside <<<mhd_check_tau>>>, k, tau(0,0,k) =", k, tau(0,0,k)
write(*,*) "Inside <<<mhd_check_tau>>>, st_x(0,0,k) =", st_x(0,0,k)
write(*,*) "Inside <<<mhd_check_tau>>>, mhd_st_x(0,0,k) =", mhd_st_x(0,0,k) 
write(*,*) "Inside <<<mhd_check_tau>>>, rho_star(0,0,k) =", rho_star(0,0,k) 
write(*,*) "Inside <<<mhd_check_tau>>>, phi(0,0,k) =", phi(0,0,k)
write(*,*) "Inside <<<mhd_check_tau>>>, gxx(0,0,k) =", gxx(0,0,k)
write(*,*) "Inside <<<mhd_check_tau>>>, lapm1(0,0,k) =", lapm1(0,0,k)
write(*,*) "Inside <<<mhd_check_tau>>>, shiftx(0,0,k) =", shiftx(0,0,k)
end do
end if

end  subroutine mhd_check_tau