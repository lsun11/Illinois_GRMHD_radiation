
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine wdns_initialdata_read_binfile_driver(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  real*8 		:: dX,dY,dZ
  real*8                :: gamma1, gamma2, gamma3
  real*8                :: kappa1, kappa2, kappa3
  real*8                :: rhoo1, rhoo2
  real*8                :: an1, an2, an3, c2
  integer               :: i,ierr

  write(*,*) "If Valgrind gives an error right before this line, ignore the error."

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ext = cctk_lsh


  write(*,*) "n_ijk:",cctk_lsh
  write(*,*) "ZZZ:",Z(1,1,:)
!  write(*,*) "X:",X(:,1,1)

      if (lorene_init.eq.0) then
         call read_inputfile_wdns(cctkGH,genID_cmdline_output_enable,fisheye_enable,reset_shift_lapse, &
         X(1,1,1),Y(1,1,1),Z(1,1,1), &
         dX,dY,dZ,cctk_lsh, &
         lapm1,shiftx,shifty,shiftz, &
         phi,Axx,Axy,Axz,Ayy,Ayz, &
         rho_b,vx,vy,vz, &
         gamma1, gamma2, gamma3, kappa1, kappa2, kappa3, rhoo1, rhoo2)
      else if (lorene_init.eq.1) then
         call read_inputfile_wdns_lorene(cctkGH,genID_cmdline_output_enable,fisheye_enable,reset_shift_lapse, &
         X(1,1,1),Y(1,1,1),Z(1,1,1), &
         dX,dY,dZ,cctk_lsh, &
         lapm1,shiftx,shifty,shiftz, &
         phi,Axx,Axy,Axz,Ayy,Ayz, &
         rho_b,vx,vy,vz, &
         gamma1, gamma2, gamma3, kappa1, kappa2, kappa3, rhoo1, rhoo2)
      else	
         call read_inputfile_wdns_lorene2(cctkGH,genID_cmdline_output_enable,fisheye_enable,reset_shift_lapse, &
         X(1,1,1),Y(1,1,1),Z(1,1,1), &
         dX,dY,dZ,cctk_lsh, &
         lapm1,shiftx,shifty,shiftz, &
         phi,Axx,Axy,Axz,Ayy,Ayz, &
         rho_b,vx,vy,vz, &
         gamma1, gamma2, gamma3, kappa1, kappa2, kappa3, rhoo1, rhoo2)
      end if

!!$     write(*,*) "an1,an2,an3 =",gamma1,gamma2,gamma3
!!$     write(*,*) "Gam1,Gam2,Gam3 = ",1.+1./gamma1,1.+1./(gamma2),1.+1./(gamma3)
!!$     write(*,*) "kappa1,kappa2,kappa3 =",kappa1,kappa2,kappa3
!!$     write(*,*) "rhoo1,rhoo2",rhoo1,rhoo2
!!$     stop


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


  end if



!  write(6,*)'gamma1,BigM,bh_posn_x(1) = ',gamma_th,BigM,bh_posn_x(1)


  neos=2
  rho_tab(1)= rhoo1
  rho_tab(2)= rhoo2
  an1 = gamma1
  an2 = gamma2
  an3 = gamma3
  gamma_tab(1) = 1.d0 + 1.d0/an1
  gamma_tab(2) = 1.d0 + 1.d0/an2
  gamma_tab(3) = 1.d0 + 1.d0/an3
  k_tab(1) = kappa1
  k_tab(2) = kappa2
  k_tab(3) = kappa3

  do i=1,2
     P_tab(i)=k_tab(i)*(rho_tab(i)**gamma_tab(i))
  enddo

    eps_tab(1)=an1*P_tab(1)/rhoo1

    c2 = (an1-an2)*kappa1*rhoo1**(1.d0/an1)

    eps_tab(2)=an2*P_tab(2)/rhoo2 + c2


end subroutine wdns_initialdata_read_binfile_driver
