
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bhns_initialdata_read_binfile_driver(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  real*8 		:: dX,dY,dZ
  integer               :: i,ierr
  write(*,*) "If Valgrind gives an error right before this line, ignore the error."

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ext = cctk_lsh

  !Set initial BH position to zero.
  bh_posn_x = 0.D0
  bh_posn_y = 0.D0
  bh_posn_z = 0.D0

  write(*,*) "n_ijk:",cctk_lsh
  write(*,*) "ZZZ:",Z(1,1,:)
!  write(*,*) "X:",X(:,1,1)

  if(use_new_bhns_initial_data==0) then
     call read_inputfile_bhns(cctkGH,genID_cmdline_output_enable,fisheye_enable,reset_shift_lapse, &
          X(1,1,1),Y(1,1,1),Z(1,1,1), &
          dX,dY,dZ,cctk_lsh, &
          lapm1,shiftx,shifty,shiftz, &
          phi,Axx,Axy,Axz,Ayy,Ayz, &
          rho_b,vx,vy,vz, &
          gamma_th,bh_posn_x(1),BigM)
  else
     call read_inputfile_bhns_bhspin(cctkGH,genID_cmdline_output_enable,fisheye_enable, &
          X(1,1,1),Y(1,1,1),Z(1,1,1), &
          dX,dY,dZ,cctk_lsh, &
          lapm1,shiftx,shifty,shiftz, &
          phi,Axx,Axy,Axz,Ayy,Ayz,Azz, &
          rho_b,vx,vy,vz, &
          gamma_th,bh_posn_x(1),BigM,initial_ns_coord_x,initial_ns_coord_y)
  end if

  if(genID_cmdline_output_enable.ne.1) then
     lapm1 = lapm1 - 1.D0
  else
     lapm1=0.D0

     shiftx=0.D0
     shifty=0.D0
     shiftz=0.D0

     phi = 1.D0

     gxx = 1.D0
     gxy = 0.D0
     gxz = 0.D0
     gyy = 1.D0
     gyz = 0.D0
     gzz = 1.D0

     Axx = 0.D0
     Axy = 0.D0
     Axz = 0.D0
     Ayy = 0.D0
     Ayz = 0.D0
     Azz = 0.D0

     rho_b = 1.D-7
     vx = 0.D0
     vy = 0.D0
     vz = 0.D0

     gamma_th = 1.3D0
  end if



  write(6,*)'gamma_th,BigM,bh_posn_x(1) = ',gamma_th,BigM,bh_posn_x(1)

  ! K_poly=1 because Keisuke's datafiles are in polytropic units already!!!
  k_poly=1.0d0
  n_poly=1.0d0/(gamma_th-1.0d0)
  write(6,*)'BH irreducible mass scaling factor:',BigM
  write(6,*)'BH position x=',bh_posn_x(1),' xbh/mirr=',bh_posn_x(1)/BigM

  neos=1
  rho_tab(1)=1.0
  P_tab(1)=K_poly
  eps_tab(1)=K_poly/(gamma_th-1.0d0)
  do i=1,2
     k_tab(i)=K_poly
     gamma_tab(i)=gamma_th
  enddo

end subroutine bhns_initialdata_read_binfile_driver
