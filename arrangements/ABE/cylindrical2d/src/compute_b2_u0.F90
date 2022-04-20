!------------------------------------------------------------------------
! Compute b^2 (store in Pr) and u0 for data output 
!------------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine cylindrical2d_compute_b2_u0(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  integer :: i,j,k
  real*8 :: s4pi,alpha,fac,sb0_l,sbx_l,sby_l,sbz_l,u0_l,ux_l,uy_l,uz_l
  real*8 :: Bx_l,By_l,Bz_l
  real*8, parameter :: pi = 3.1415926535897932d0
!
  if (mod(cctk_iteration,out_every)==0) then
     s4pi = 2.d0*sqrt(pi)
     do k=1,cctk_lsh(3)
	do j=1,cctk_lsh(2)
	   do i=1,cctk_lsh(1)
	      u0_l = 1.d0/sqrt(1.d0-vx(i,j,k)**2 - vy(i,j,k)**2 - vz(i,j,k)**2)
	      ux_l = u0_l*vx(i,j,k)
              uy_l = u0_l*vy(i,j,k)
              uz_l = u0_l*vz(i,j,k)
	      Bx_l = Bx(i,j,k)/s4pi 
	      By_l = By(i,j,k)/s4pi
	      Bz_l = Bz(i,j,k)/s4pi
	      sb0_l = ux_l*Bx_l + uy_l*By_l + uz_l*Bz_l
	      sbx_l = (Bx_l + ux_l*sb0_l)/u0_l
              sby_l = (By_l + uy_l*sb0_l)/u0_l
              sbz_l = (Bz_l + uz_l*sb0_l)/u0_l
	      u0(i,j,k) = u0_l
	      Pr(i,j,k) = sbx_l**2 + sby_l**2 + sbz_l**2 - sb0_l**2
	   end do
	end do
     end do
  end if
end subroutine cylindrical2d_compute_b2_u0

