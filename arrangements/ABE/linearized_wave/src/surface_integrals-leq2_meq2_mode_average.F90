
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!-------------------------------------------------------------------
!
! Calculate l=2,m=2 decomposition of Psi_4:
!   Psi_4^{l=2,m=2}(t,r) = Integral dOmega Psi_4(t,r,theta,phi) Y^{l=2,m=2}(theta,phi)
!
!-------------------------------------------------------------------

subroutine leq2_meq2_mode_average_lw(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer, dimension(3)                    :: ext,global_ext
  real*8                                   :: dX,dY,dZ
  integer                                  :: N_theta,N_phi,sym_factor,ntot
  real*8                                   :: dcostheta,dphi,hplus_anal,hcross_anal,E_GW_anal,psi4r,psi4i
  real*8                                   :: phiangle,costheta,sintheta,PI,Y22nophi,sintwophi,costwophi
  integer                                  :: i,j,n
  integer                                  :: interp_order
  integer                                  :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ext = cctk_lsh
  global_ext = cctk_gsh

  PI = 3.14159265358979323846D0

  interp_order = 1

  N_theta = 80
  N_phi = 80
  dphi=2.D0 * PI / N_phi
  dcostheta = 2.D0 / N_theta

  if (Symmetry==OCTANT) then 
     N_theta = N_theta/2
     N_phi = N_phi/4
     sym_factor=8
  else if (Symmetry==EQUATORIAL) then
     N_theta = N_theta/2
     sym_factor=2
  else if (Symmetry==NO_SYMM) then 
     sym_factor=1
  else if (Symmetry==PI_SYMM) then 
     N_theta = N_theta/2
     N_phi = N_phi/2
     sym_factor=4
  else if (Symmetry==AXISYM) then
     N_theta = N_theta/2
     N_phi = 1
     dphi=2.0*PI
     sym_factor=2
  end if
   
  ntot = N_theta*N_phi

  Psi4resumlw = 0.D0
  Psi4imsumlw = 0.D0

  n = 1
  do i=1,N_theta
     costheta = 1.D0 - (i - 0.5D0)*dcostheta
     sintheta = sqrt(1.D0 - costheta*costheta)
     do j=1,N_phi
        phiangle = (j - 0.5D0)*dphi
        if(N_phi==1) phiangle = 0.D0

        sintwophi = sin(2.D0*phiangle)
        costwophi = cos(2.D0*phiangle)

        Y22nophi = 0.25D0*sqrt(15.D0/(2.D0*PI))*sintheta*sintheta

        call gw_anal(CCTK_TIME+time_shift,radius_GW_phys,acos(costheta),phiangle, &
             hplus_anal,hcross_anal,E_GW_anal,psi4r,psi4i, & 
             amplitude,width,mode)

        Psi4resumlw = Psi4resumlw + Y22nophi*(psi4r*costwophi + psi4i*sintwophi)
        Psi4imsumlw = Psi4imsumlw + Y22nophi*(-psi4r*sintwophi + psi4i*costwophi)

!        write(*,*) "zy",radius_GW_phys*sintheta*cos(phiangle),radius_GW_phys*sintheta*sin(phiangle),radius_GW_phys*costheta

! Check that integral of Y22 Y22star equals 1.0:
!        Psi4resumlw = Psi4resumlw + Y22nophi*Y22nophi
!        Psi4imsumlw = Psi4imsumlw + Y22nophi*Y22nophi
        
        n = n + 1
     end do
  end do

  Psi4resumlw = Psi4resumlw * sym_factor * dcostheta * dphi
  Psi4imsumlw = Psi4imsumlw * sym_factor * dcostheta * dphi

end subroutine leq2_meq2_mode_average_lw
