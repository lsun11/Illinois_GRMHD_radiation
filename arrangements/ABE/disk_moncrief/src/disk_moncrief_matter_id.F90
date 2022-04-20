#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine disk_moncrief_matter_id(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  real*8 		:: dX,dY,dZ
  integer               :: i,j,k
  character             :: varname*30
  real*8                :: RoM,costheta,cos2theta,sin2theta,varpi
  real*8                :: Delta,Delta_in,Sigma,Sigma_in,A,A_in
  real*8                :: zeta,zeta_in,xi,xi_in,log_enthalpy,enthalpy,epsilon,rho0,gphiphi,gphit,ut,vphi
!  real*8                :: BigMass
  REAL*8, PARAMETER :: PI_D=3.141592653589793238462643383279502884197
  real*8                     :: HALF, ONE, ZERO, TWO, FOUR
  parameter(HALF = 0.5D0, ONE = 1.D0, ZERO = 0.D0, TWO = 2.D0, FOUR = 4.D0)


  write(*,*) "*****************************"
  write(*,*) "inside disk_moncrief_matter_id"
  write(*,*) "*****************************"

  Delta = 0.1d0

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ext = cctk_lsh

  !commented out temporarily, but something is needed
  ! if (genID_cmdline_output_enable .eq. 1) BigMass=1.d0
  !temporary hack
!  BigMass=1.0
  write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,*) "!!!!!!!!!!!WARNING WARNING !!!!!!!!!!!!!!"
  write(*,*) "!!!!!!!!!get the binary mass right!!!!!!!"
  write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"


  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
           RoM = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)/BigMass
           costheta = (Z(i,j,k)/BigMass)/RoM
           cos2theta = costheta*costheta
           sin2theta = ONE - cos2theta

           ! this is not the physical varpi, but it is only used
           ! for calculating v^y.
           varpi = sqrt(x(i,j,k)**2 + y(i,j,k)**2)/BigMass
            
           if(varpi .gt. 0.9d0*RoM_in) then
              !    if(RoM .gt. 0.9d0*RoM_in) then
              Delta    = RoM**2    - TWO*RoM    + sam_disk*sam_disk
              Delta_in = RoM_in**2 - TWO*RoM_in + sam_disk*sam_disk

              Sigma    = RoM**2    + sam_disk*sam_disk*cos2theta
              Sigma_in = RoM_in**2

              A    = (RoM*RoM + sam_disk*sam_disk)**2 - Delta*sam_disk*sam_disk*sin2theta
              A_in = RoM_in*(RoM_in**3 + RoM_in*sam_disk*sam_disk + TWO*sam_disk*sam_disk)

              zeta = ONE + sqrt(ONE + FOUR*ell**2*Sigma**2*Delta/(A**2*sin2theta))
              zeta = zeta * A / (Sigma * Delta)
              zeta_in = ONE + sqrt(ONE + FOUR*ell**2*Sigma_in**2*Delta_in/A_in**2)
              zeta_in = zeta_in * A_in / (Sigma_in * Delta_in)

              xi = -HALF * sqrt(ONE + FOUR*ell**2*Sigma**2*Delta/(A**2*sin2theta))
              xi = xi - TWO*sam_disk*RoM*ell/A
              xi_in = -HALF * sqrt(ONE + FOUR*ell**2*Sigma_in**2*Delta_in/(A_in**2))
              xi_in = xi_in - TWO*sam_disk*RoM_in*ell/A_in

              log_enthalpy = HALF*log(zeta) + xi - HALF*log(zeta_in) - xi_in 

              if(log_enthalpy.gt.ZERO) then
                 enthalpy = exp(log_enthalpy)
                 epsilon = (enthalpy - ONE)/Gamma_th
                 rho0 = ( (Gamma_th-ONE)*epsilon/(K_poly) )**(ONE/(Gamma_th-ONE))
              else
                 enthalpy = ONE
                 epsilon = ZERO
                 rho0 = ZERO
              end if

              rho_b(i,j,k) = rho0

              !may want to change these
              gphiphi = sin2theta*(Sigma + sam_disk*sam_disk*(1+TWO*RoM/Sigma)*sin2theta)
              gphit = -TWO*sam_disk*RoM*sin2theta/Sigma

              ut = gphiphi+sqrt(gphiphi**2 + FOUR*ell**2*Delta*sin2theta)
              ut = sqrt(ut/(TWO*Delta*sin2theta))

              vphi =  ONE/ut * (ell - gphit*ut**2)/(gphiphi*ut)
              vx(i,j,k) = -Y(i,j,k)*vphi
              vy(i,j,k) = X(i,j,k)*vphi
              vz(i,j,k) = ZERO
              u0(i,j,k) = ut
           else
              rho_b(i,j,k) = 0.d0
              vx(i,j,k) = ZERO
              vy(i,j,k) = ZERO
              vz(i,j,k) = ZERO
             ! vx(i,j,k) = -shiftx(i,j,k)
             ! vy(i,j,k) = -shifty(i,j,k)
             ! vz(i,j,k) = -shiftz(i,j,k)
              
              u0(i,j,k) = 1.d0/(lapm1(i,j,k)+1.d0)
                   
           endif
           if ((i.eq.42).and.(j.eq.42).and.(k.eq.5).and.(1.eq.0)) then
              write(*,*) "X: ",X(i,j,k)
              write(*,*) "Y: ",Y(i,j,k)
              write(*,*) "Z: ",Z(i,j,k)
              write(*,*) "vx: ",vx(i,j,k)
              write(*,*) "vy: ",vy(i,j,k)
              write(*,*) "vz: ",vz(i,j,k)
         
              write(*,*) "shiftx: ",shiftx(i,j,k)
              write(*,*) "shifty: ",shifty(i,j,k)
              write(*,*) "shiftz: ",shiftz(i,j,k)
              
              write(*,*) "gxx: ",gxx(i,j,k)
              write(*,*) "gxy: ",gxy(i,j,k)
              write(*,*) "gxz: ",gxz(i,j,k)
              write(*,*) "gyy: ",gyy(i,j,k)
              write(*,*) "gyz: ",gyz(i,j,k)
              write(*,*) "gzz: ",gzz(i,j,k)
           endif
        end do
     end do
  end do

  neos=1
  rho_tab(1)=1.d0
  P_tab(1)=K_poly
  eps_tab(1)=P_tab(1)/rho_tab(1)/(gamma_th-1.0d0)
  do i=1,2
     k_tab(i)=K_poly
     gamma_tab(i)=gamma_th
  enddo


end subroutine disk_moncrief_matter_id

