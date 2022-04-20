#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!---------------------------------------------------------------------
! This routine computes the areal radius of a spherical surface of 
!   coordinate radius radius_GW and the averaged 1-lapse
!---------------------------------------------------------------------
subroutine areal_radius_averaged_1minus_lapse(cctkGH,Symmetry,N_theta,N_phi,ntot, &
     dphi,dcostheta,radius_GW,phi,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,lapm1, & 
     r_areal,avg_1_lapse)
  implicit none

  CCTK_POINTER :: cctkGH  
  real*8                                   :: radius_GW
  integer                                  :: N_theta,N_phi,sym_factor,ntot,n
  real*8                                   :: dcostheta,dphi

  real*8                                   :: r_areal,avg_1_lapse

  real*8, dimension(ntot)                  :: phi,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,lapm1
  real*8                                   :: phiangle,costheta,sintheta
  real*8				   :: area,x,y,z,dS
  real*8, parameter 			   :: PI = 3.14159265358979323846D0
  integer                                  :: i,j
  integer                                  :: Symmetry
  integer                                  :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  if (Symmetry==OCTANT) then 
     sym_factor=8
  else if (Symmetry==EQUATORIAL) then
     sym_factor=2
  else if (Symmetry==NO_SYMM) then 
     sym_factor=1
  else if (Symmetry==PI_SYMM) then 
     sym_factor=4
  else if (Symmetry==AXISYM) then
     sym_factor=2
  end if

  area = 0.d0
  avg_1_lapse = 0.d0

  n = 1
  do i=1,N_theta
     costheta = 1.D0 - (i - 0.5D0)*dcostheta
     sintheta = sqrt(1.D0 - costheta*costheta)

     do j=1,N_phi
        phiangle = (j - 0.5D0)*dphi
        if(N_phi==1) phiangle = 0.D0

	!-------------------------------------------------------------
	! The surface area of a sphere of coordinate radius r is given 
 	! by 
	!    S = Integrate[ psi^6 r^2 sqrt(gamma^{ij} x_i x_j)/r dcostheta dphi ] 
        !      = Integrate[ psi^4 r sqrt(tilde{gamma}^{ij} x_i x_j) dcostheta dphi ]
        !--------------------------------------------------------------
	x = radius_GW*sintheta*cos(phiangle)
	y = radius_GW*sintheta*sin(phiangle)
	z = radius_GW*costheta
	dS = exp(4.d0*phi(n))*radius_GW*sqrt( gupxx(n)*x*x + &
                 2.d0*gupxy(n)*x*y + 2.d0*gupxz(n)*x*z + & 
                 gupyy(n)*y*y + 2.d0*gupyz(n)*y*z + gupzz(n)*z*z )
	area = area + dS
	avg_1_lapse = avg_1_lapse - lapm1(n)*dS
	
        n = n + 1 
     end do
  end do

  avg_1_lapse = avg_1_lapse / area
  area = area * sym_factor * dcostheta * dphi
  r_areal = sqrt(area/(4.d0*PI))
  
end subroutine areal_radius_averaged_1minus_lapse
