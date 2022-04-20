
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bbh_cookpfeiffer_rot_read_inputfile_driver(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8 :: dX,dY,dZ
  integer :: i,j,k
  integer :: gridnumber

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  gridnumber=1
  call read_inputfile_bbh_cookpfeiffer_rot(genID_cmdline_output_enable,gridnumber, &
       Nlograd_rot_interp,Ntheta_rot_interp,Nphi_rot_interp, &
       radmin_rot_interp_1,radmax_rot_interp_1,xbh1_initial, &
       K_rr_rot1,K_rth_rot1,K_rp_rot1,K_thth_rot1,K_thp_rot1,K_pp_rot1, &
       shiftr_rot1,shiftth_rot1,shiftp_rot1,phi_rot1,lapm1_rot1)
  gridnumber=2
  call read_inputfile_bbh_cookpfeiffer_rot(genID_cmdline_output_enable,gridnumber, &
       Nlograd_rot_interp,Ntheta_rot_interp,Nphi_rot_interp, &
       radmin_rot_interp_2,radmax_rot_interp_2,xbh2_initial, &
       K_rr_rot2,K_rth_rot2,K_rp_rot2,K_thth_rot2,K_thp_rot2,K_pp_rot2, &
       shiftr_rot2,shiftth_rot2,shiftp_rot2,phi_rot2,lapm1_rot2)
  gridnumber=3
  call read_inputfile_bbh_cookpfeiffer_rot(genID_cmdline_output_enable,gridnumber, &
       Nlograd_rot_interp,Ntheta_rot_interp,Nphi_rot_interp, &
       radmin_rot_interp_3,radmax_rot_interp_3,0.d0, &
       K_rr_rot3,K_rth_rot3,K_rp_rot3,K_thth_rot3,K_thp_rot3,K_pp_rot3, &
       shiftr_rot3,shiftth_rot3,shiftp_rot3,phi_rot3,lapm1_rot3)
  
!  do k=1,20
!     do j=1,20
!        do i=1,20
!           write(*,*) lograd_arr1(i+20*((j-1)+20*(k-1))),theta_arr1(i+20*((j-1)+20*(k-1))),phi_arr1(i+20*((j-1)+20*(k-1)))
!        end do
!     end do
!  end do
  
  if(genID_cmdline_output_enable==1) then
     phi = 0.D0
     kxx = 0.D0
     kxy = 0.D0
     kxz = 0.D0
     kyy = 0.D0
     kyz = 0.D0
     lapm1 = 0.D0
     shiftx = 0.D0
     shifty = 0.D0
     shiftz = 0.D0
  end if
  
end subroutine bbh_cookpfeiffer_rot_read_inputfile_driver
