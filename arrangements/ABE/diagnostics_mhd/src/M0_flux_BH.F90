#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine M0_flux_BH(cctkGH,F_M0,N_theta,N_phi,Symmetry,found_horizon)
  implicit none
!  DECLARE_CCTK_ARGUMENTS
!  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_POINTER :: cctkGH
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  real*8 :: F_M0_ijk,F_M0
  real*8 :: vx_i, vy_i, vz_i
  real*8 :: rhos_i
  real*8 :: rr_i, dS, PI, dX, dY, dZ, f1ospi,xh,yh,zh,dRdmu,dRdphi
  real*8 :: cosphi,sinphi,dmudx,dmudy,dmudz,dphidx,dphidy,dphidz,sym_factor
  real*8 :: sintheta,costheta,nn,xn,yn,phiangle,dphi,dcostheta
  integer :: i,j, interp_order, horizon_number, foundflag, n,ntot
  integer :: ind0,ind1,ind2,indm1,indm2,vindex
  integer :: N_theta,N_phi, Symmetry, found_horizon
  real*8, dimension(N_theta*N_phi,3)        :: pointcoords
  real*8, dimension(N_theta*N_phi)          :: nx_d,ny_d,nz_d
  real*8, dimension(N_theta*N_phi)          :: ah_radii, x_ah1,y_ah1,z_ah1
  real*8, dimension(N_theta*N_phi)          :: vxint, vyint, vzint
  real*8, dimension(N_theta*N_phi)          :: rhosint
  CCTK_POINTER,dimension(4)                 :: output_array_pointers

  if (Symmetry .ne. EQUATORIAL .and. Symmetry .ne. NO_SYMM) then 
     write(*,*) 'Symmetry not supported in M0_flux_BH'
     stop
  end if 

  ntot = N_theta*N_phi

  PI = 3.14159265358979323846D0
  dphi = 2.d0 * PI / N_phi
  if (Symmetry==EQUATORIAL) then
     dcostheta = 1.d0 / N_theta
     sym_factor = 2.d0
  else 
     dcostheta = 2.d0 / N_theta
     sym_factor = 1.d0
  end if

! Get the origin of BH
  horizon_number = 1
  foundflag = HorizonLocalCoordinateOrigin(horizon_number,xh,yh,zh)

  n = 1
  do i=1,N_theta
     costheta = 1.d0 - (i - 0.5d0)*dcostheta
     sintheta = sqrt(1.d0 - costheta*costheta)
     do j=1,N_phi
        phiangle = (j - 0.5d0)*dphi
        x_ah1(n) = xh + sintheta*cos(phiangle)
        y_ah1(n) = yh + sintheta*sin(phiangle)
        z_ah1(n) = zh + costheta
        n = n + 1
     end do
  end do

  ! Find horizon radii
  foundflag = HorizonRadiusInDirection(horizon_number,ntot,x_ah1,y_ah1,z_ah1,ah_radii)

  ! Horizon not found, set the flux to 0
  ! Note: If horizon is not found, ahfinderdirect will set ah_radii to -1
  if (ah_radii(1) .lt. 0.d0) then
     F_M0 = 0.d0
     found_horizon = 0
     return
  end if
  found_horizon = 1

  ! Now set the points on the horizon surface 
  n = 1
  do i=1,N_theta
     costheta = 1.d0 - (i - 0.5d0)*dcostheta
     sintheta = sqrt(1.d0 - costheta*costheta)
     do j=1,N_phi
        phiangle = (j - 0.5d0)*dphi
        pointcoords(n,1) = xh + ah_radii(n)*sintheta*cos(phiangle)
        pointcoords(n,2) = yh + ah_radii(n)*sintheta*sin(phiangle)
        pointcoords(n,3) = zh + ah_radii(n)*costheta
        n = n + 1
     end do
  end do

  ! Let f = sqrt[(x-xh)^2+(y-yh)^2+(z-zh)^2] - R(mu,phi), where mu=cos(theta). 
  ! Now compute partial f / partial X^i (X^i = x^i-xh^i) and normalize 
  ! by the "Jacobian", storing them to ni_d. 
  ! 
  n=1
  do i=1,N_theta
     costheta = 1.d0 - (i - 0.5d0)*dcostheta
     sintheta = sqrt(1.d0 - costheta*costheta)
     do j=1,N_phi
        phiangle = (j - 0.5d0)*dphi
        cosphi = cos(phiangle)
        sinphi = sin(phiangle)

        ! dR/dmu 
        if (i==1) then 
   	   ind0 = j
	   ind1 = j + N_phi
	   ind2 = j + 2*N_phi
	   dRdmu = (1.5d0*ah_radii(ind0) - 2.d0*ah_radii(ind1) + 0.5d0*ah_radii(ind2)) / dcostheta
        else if (i==N_theta) then 
           if (Symmetry==EQUATORIAL) then
	      ind1 = j + (i-1)*N_phi  ! Equatorial sym --> ah_radii(ind1) = ah_radii(ind0)
	      indm1 = j + (i-2)*N_phi
	      dRdmu = 0.5d0*(ah_radii(indm1)-ah_radii(ind1))/dcostheta
           else ! This is for NO_SYMM
              ind0 = j + (i-1)*N_phi
              indm1 = j + (i-2)*N_phi
              indm2 = j + (i-3)*N_phi
              dRdmu = (-1.5d0*ah_radii(ind0) + 2.d0*ah_radii(indm1) - 0.5d0*ah_radii(indm2)) / dcostheta
           end if
        else
	   ind1 = j + i*N_phi
	   indm1 = j + (i-2)*N_phi
	   dRdmu = 0.5d0*(ah_radii(indm1)-ah_radii(ind1))/dcostheta
        end if

        ! dR/dphi
        if (j==1) then 
	   ind1 = n+1
 	   indm1 = i*N_phi 
	else if (j==N_phi) then 
           ind1 = 1 + (i-1)*N_phi
	   indm1 = n-1
        else 
           ind1 = n+1
	   indm1 = n-1
	end if
        dRdphi = 0.5d0*(ah_radii(ind1)-ah_radii(indm1))/dphi
        
        dmudx = -costheta*sintheta*cosphi/ah_radii(n)
	dmudy = -costheta*sintheta*sinphi/ah_radii(n)
	dmudz = sintheta*sintheta/ah_radii(n)
        dphidx = -sinphi/(ah_radii(n)*sintheta)
  	dphidy = cosphi/(ah_radii(n)*sintheta)

	nx_d(n) = sintheta*cosphi - dRdmu*dmudx - dRdphi*dphidx
        ny_d(n) = sintheta*sinphi - dRdmu*dmudy - dRdphi*dphidy
	nz_d(n) = costheta - dRdmu*dmudz
	nn = (nx_d(n)*cosphi + ny_d(n)*sinphi)*sintheta + nz_d(n)*costheta
	! normalize by the "Jacobian."
	nx_d(n) = nx_d(n)/nn
        ny_d(n) = ny_d(n)/nn
	nz_d(n) = nz_d(n)/nn

        n = n + 1
     end do
  end do

! Interpolate the grid functions onto the AH
  call CCTK_VarIndex(vindex,"mhd_evolve::vx")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,vxint)
  call CCTK_VarIndex(vindex,"mhd_evolve::vy")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,vyint)
  call CCTK_VarIndex(vindex,"mhd_evolve::vz")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,vzint)
  call CCTK_VarIndex(vindex,"mhd_evolve::rho_star")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,rhosint)

  !
  ! integrate on the horizon
  !

  F_M0 = 0.d0

  do i=1,ntot 
     ! read in data
     vx_i  = vxint(i)
     vy_i  = vyint(i)
     vz_i  = vzint(i)
     rhos_i = rhosint(i)

     F_M0_ijk = (nx_d(i)*vx_i + ny_d(i)*vy_i + nz_d(i)*vz_i)*rhos_i

     dS = dphi * dcostheta * sym_factor * ah_radii(i)**2

     F_M0  = F_M0 + F_M0_ijk*dS

  end do
  
end subroutine M0_flux_BH
