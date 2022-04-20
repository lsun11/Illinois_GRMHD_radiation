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
          
           if (abs(F_radx(i,j,k)) .gt. 1.0e2) then
              write(*,*) " In driver_check_tau, F_radx too large!!!!"
              write(*,*) "i,j,k X,Y,Z = ", i,j,k, X(i,j,k), Y(i,j,k), Z(i,j,k)
              write(*,*) "E_rad(i,j,k), F_radx(i,j,k), F_rady(i,j,k), F_radz(i,j,k), F_rad0(i,j,k)=",E_rad(i,j,k), F_radx(i,j,k), F_rady(i,j,k), F_radz(i,j,k), F_rad0(i,j,k)
              write(*,*) "tau_rad(i,j,k), S_rad_x(i,j,k), S_rad_y(i,j,k), S_rad_z(i,j,k)=",tau_rad(i,j,k), S_rad_x(i,j,k), S_rad_y(i,j,k), S_rad_z(i,j,k)
              write(*,*) "rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k), u0(i,j,k)=", rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k), u0(i,j,k)
              write(*,*) "Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)=", Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)
           end if

	end do
    end do
  end do

end  subroutine mhd_check_tau



subroutine mhd_check_tau2(CCTK_ARGUMENTS)

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

           if (abs(F_radx(i,j,k)) .gt. 1.0e2) then
              write(*,*) " In driver_check_tau2, F_radx too large!!!!"
              write(*,*) "i,j,k X,Y,Z = ", i,j,k, X(i,j,k), Y(i,j,k), Z(i,j,k)
              write(*,*) "E_rad(i,j,k), F_radx(i,j,k), F_rady(i,j,k), F_radz(i,j,k), F_rad0(i,j,k)=",E_rad(i,j,k), F_radx(i,j,k), F_rady(i,j,k), F_radz(i,j,k), F_rad0(i,j,k)
              write(*,*) "tau_rad(i,j,k), S_rad_x(i,j,k), S_rad_y(i,j,k), S_rad_z(i,j,k)=",tau_rad(i,j,k), S_rad_x(i,j,k), S_rad_y(i,j,k), S_rad_z(i,j,k)
              write(*,*) "rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k), u0(i,j,k)=", rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k), u0(i,j,k)
              write(*,*) "Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)=", Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)
           end if

        end do
    end do
  end do

end  subroutine mhd_check_tau2





subroutine mhd_check_tau3(CCTK_ARGUMENTS)

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

           if (abs(F_radx(i,j,k)) .gt. 1.0e2) then
              write(*,*) " In driver_check_tau3, F_radx too large!!!!"
              write(*,*) "i,j,k X,Y,Z = ", i,j,k, X(i,j,k), Y(i,j,k), Z(i,j,k)
              write(*,*) "E_rad(i,j,k), F_radx(i,j,k), F_rady(i,j,k), F_radz(i,j,k), F_rad0(i,j,k)=",E_rad(i,j,k), F_radx(i,j,k), F_rady(i,j,k), F_radz(i,j,k), F_rad0(i,j,k)
              write(*,*) "tau_rad(i,j,k), S_rad_x(i,j,k), S_rad_y(i,j,k), S_rad_z(i,j,k)=",tau_rad(i,j,k), S_rad_x(i,j,k), S_rad_y(i,j,k), S_rad_z(i,j,k)
              write(*,*) "rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k), u0(i,j,k)=", rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k), u0(i,j,k)
              write(*,*) "Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)=", Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)
           end if

        end do
    end do
  end do

end  subroutine mhd_check_tau3



subroutine mhd_check_tau4(CCTK_ARGUMENTS)

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

           if (abs(F_radx(i,j,k)) .gt. 1.0e2) then
              write(*,*) " In driver_check_tau4, F_radx too large!!!!"
              write(*,*) "i,j,k X,Y,Z = ", i,j,k, X(i,j,k), Y(i,j,k), Z(i,j,k)
              write(*,*) "E_rad(i,j,k), F_radx(i,j,k), F_rady(i,j,k), F_radz(i,j,k), F_rad0(i,j,k)=",E_rad(i,j,k), F_radx(i,j,k), F_rady(i,j,k), F_radz(i,j,k), F_rad0(i,j,k)
              write(*,*) "tau_rad(i,j,k), S_rad_x(i,j,k), S_rad_y(i,j,k), S_rad_z(i,j,k)=",tau_rad(i,j,k), S_rad_x(i,j,k), S_rad_y(i,j,k), S_rad_z(i,j,k)
              write(*,*) "rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k), u0(i,j,k)=", rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k), u0(i,j,k)
              write(*,*) "Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)=", Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)
           end if

        end do
    end do
  end do
end  subroutine mhd_check_tau4

subroutine mhd_check_tau5(CCTK_ARGUMENTS)

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

          if (abs(F_radx(i,j,k)) .gt. 1.0e2) then
              write(*,*) " In driver_check_tau5, F_radx too large!!!!"
              write(*,*) "i,j,k X,Y,Z = ", i,j,k, X(i,j,k), Y(i,j,k), Z(i,j,k)
              write(*,*) "E_rad(i,j,k), F_radx(i,j,k), F_rady(i,j,k), F_radz(i,j,k), F_rad0(i,j,k)=",E_rad(i,j,k), F_radx(i,j,k), F_rady(i,j,k), F_radz(i,j,k), F_rad0(i,j,k)
              write(*,*) "tau_rad(i,j,k), S_rad_x(i,j,k), S_rad_y(i,j,k), S_rad_z(i,j,k)=",tau_rad(i,j,k), S_rad_x(i,j,k), S_rad_y(i,j,k), S_rad_z(i,j,k)
              write(*,*) "rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k), u0(i,j,k)=", rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k), u0(i,j,k)
              write(*,*) "Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)=", Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)
           end if

        end do
    end do
  end do
end  subroutine mhd_check_tau5
