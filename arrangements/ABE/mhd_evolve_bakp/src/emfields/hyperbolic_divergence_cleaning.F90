!---------------------------------------------------------------------
!       :: Driver routine for hyperbolic divergence cleaning ::
!  See e.g., M. Anderson et al. Class.Quant.Grav. 23 (2006) 6503-6524 
!---------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "GenericFD.h"

subroutine hyperbolic_divergence_cleaning(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  real*8                :: dX,dY,dZ,c_h,c_p,Blagrangemultiplierx,Blagrangemultipliery,Blagrangemultiplierz,Bxx,Byy,Bzz,KO_Strength,r_from_BH
  integer               :: index,ierr,handle,dummy
  CCTK_REAL             :: reduction_value
  integer               :: AXISYM,i,j,k
  parameter(AXISYM = 4)
  ! 1st of 2 needed #includes for GenericFD.h:
#include "../../../GenFD_decl_varF90.h"

  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  !     Initialise finite differencing variables
  ! 2nd of 2 needed #includes for GenericFD.h:
#include "../../../GenFD_set_varF90.h"

  !!  write(*,*) "INSIDE hyperbolic_divergence_cleaning.  NUM_BHS = ",num_BHs,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1)

  ! Can't OpenMP this loop, since c_h and c_p are local variables!  TODO: rewrite this routine in C++.
  do k=cctk_nghostzones(3),ext(3)-cctk_nghostzones(3)
     do j=cctk_nghostzones(2),ext(2)-cctk_nghostzones(2)
        do i=cctk_nghostzones(1),ext(1)-cctk_nghostzones(1)

           ! default = 1.0, higher and you'll need to lower the Courant factor
           c_h = c_h_default
           ! between 1 and 12, the higher, the better it is at fixing divergence around larger shocks
           c_p = c_p_default

           if(hyperbolic_divergence_cleaning_centered_differencing==1) then
              Blagrangemultiplierx = D1_c2(Blagrangemultiplier,i,j,k)
              Blagrangemultipliery = D2_c2(Blagrangemultiplier,i,j,k)
              Blagrangemultiplierz = D3_c2(Blagrangemultiplier,i,j,k)
              Bxtilde_or_Ax_rhs(i,j,k) = Bxtilde_or_Ax_rhs(i,j,k) - &
                   exp(4.D0*phi(i,j,k))*( gupxx(i,j,k)*Blagrangemultiplierx + gupxy(i,j,k)*Blagrangemultipliery + gupxz(i,j,k)*Blagrangemultiplierz)
              Bytilde_or_Ay_rhs(i,j,k) = Bytilde_or_Ay_rhs(i,j,k) - &
                   exp(4.D0*phi(i,j,k))*( gupxy(i,j,k)*Blagrangemultiplierx + gupyy(i,j,k)*Blagrangemultipliery + gupyz(i,j,k)*Blagrangemultiplierz)
              Bztilde_or_Az_rhs(i,j,k) = Bztilde_or_Az_rhs(i,j,k) - &
                   exp(4.D0*phi(i,j,k))*( gupxz(i,j,k)*Blagrangemultiplierx + gupyz(i,j,k)*Blagrangemultipliery + gupzz(i,j,k)*Blagrangemultiplierz)

              Bxx = D1_c2(Bxtilde,i,j,k)
              Byy = D2_c2(Bytilde,i,j,k)
              Bzz = D3_c2(Bztilde,i,j,k)
              Blagrangemultiplier_rhs(i,j,k) = -c_h**2/c_p**2*Blagrangemultiplier(i,j,k) - c_h**2*( &
                   Bxx+Byy+Bzz)
           else
              Blagrangemultiplierx = 0.25D0* ( &
                   (Blagrangemultiplier(i+1,j,k) - Blagrangemultiplier(i,j,k))/dx + &
                   (Blagrangemultiplier(i+1,j+1,k) - Blagrangemultiplier(i,j+1,k))/dx + &
                   (Blagrangemultiplier(i+1,j,k+1) - Blagrangemultiplier(i,j,k+1))/dx + &
                   (Blagrangemultiplier(i+1,j+1,k+1) - Blagrangemultiplier(i,j+1,k+1))/dx)
              Blagrangemultipliery = 0.25D0* ( &
                   (Blagrangemultiplier(i,  j+1,k) - Blagrangemultiplier(i,    j,k))/dy + &
                   (Blagrangemultiplier(i+1,j+1,k) - Blagrangemultiplier(i+1,  j,k))/dy + &
                   (Blagrangemultiplier(i,  j+1,k+1) - Blagrangemultiplier(i,  j,k+1))/dy + &
                   (Blagrangemultiplier(i+1,j+1,k+1) - Blagrangemultiplier(i+1,j,k+1))/dy)
              Blagrangemultiplierz = 0.25D0* ( &
                   (Blagrangemultiplier(i,j,    k+1) - Blagrangemultiplier(i,j,    k))/dz + &
                   (Blagrangemultiplier(i,j+1,  k+1) - Blagrangemultiplier(i,j+1,  k))/dz + &
                   (Blagrangemultiplier(i+1,j,  k+1) - Blagrangemultiplier(i+1,j,  k))/dz + &
                   (Blagrangemultiplier(i+1,j+1,k+1) - Blagrangemultiplier(i+1,j+1,k))/dz)

              Bxtilde_or_Ax_rhs(i,j,k) = Bxtilde_or_Ax_rhs(i,j,k) - &
                   exp(4.D0*phi(i,j,k))*( gupxx(i,j,k)*Blagrangemultiplierx + gupxy(i,j,k)*Blagrangemultipliery + gupxz(i,j,k)*Blagrangemultiplierz)
              Bytilde_or_Ay_rhs(i,j,k) = Bytilde_or_Ay_rhs(i,j,k) - &
                   exp(4.D0*phi(i,j,k))*( gupxy(i,j,k)*Blagrangemultiplierx + gupyy(i,j,k)*Blagrangemultipliery + gupyz(i,j,k)*Blagrangemultiplierz)
              Bztilde_or_Az_rhs(i,j,k) = Bztilde_or_Az_rhs(i,j,k) - &
                   exp(4.D0*phi(i,j,k))*( gupxz(i,j,k)*Blagrangemultiplierx + gupyz(i,j,k)*Blagrangemultipliery + gupzz(i,j,k)*Blagrangemultiplierz)

              if(num_BHs.ge.1) then
                 r_from_BH = sqrt((x(i,j,k)-bh_posn_x(1))**2 + (y(i,j,k)-bh_posn_y(1))**2 + (z(i,j,k)-bh_posn_z(1))**2)
                 if(r_from_BH .lt. min_BH_radius) then
!                    c_h = exp(-4.D0/(1.D-10 + r_from_BH))*c_h_default
                    c_h = exp(-4.D0/(1.D-10 + r_from_BH))*c_h_default
                    c_p = c_h*c_p_default/c_h_default
                 end if

                 if(num_BHs.eq.2) then 
                    r_from_BH = sqrt((x(i,j,k)-bh_posn_x(2))**2 + (y(i,j,k)-bh_posn_y(2))**2 + (z(i,j,k)-bh_posn_z(2))**2)
                    if(r_from_BH .lt. min_BH_radius) then
!                       c_h = exp(-4.D0/(1.D-10 + r_from_BH))*c_h_default
                       c_h = exp(-4.D0/(1.D-10 + r_from_BH))*c_h_default
                       c_p = c_h*c_p_default/c_h_default
                    end if
                 end if

                 if(num_BHs.gt.2) then
                    write(*,*) "SORRY, Hyperbolic Divergence Cleaning DOES NOT CURRENTLY SUPPORT num_BHs.gt.2"
                    stop
                 end if
              end if

              ! Gotta be careful!  If c_p is too close to zero or equal to zero, c_h**2/c_p**2 will blow up!
              if(c_p.lt.1.D-20) then
                 c_h = 1.D-20
                 c_p = c_h*c_p_default/c_h_default
              end if

              !FAILED TRIES:
              !                 c_h = exp(-8.D0*phi(i,j,k))
              !                 c_p = exp(-8.D0*phi(i,j,k))
              !                 c_h = exp(-4.D0*phi(i,j,k))
              !                 c_p = exp(-4.D0*phi(i,j,k))
              !                 c_h = exp(-10.D0*phi(i,j,k))
              !                 c_p = exp(-10.D0*phi(i,j,k))
              !FAILED TRIES: GAUSSIAN-LIKE
              !                    c_h = -1.D0*exp(-(r(i,j,k)*2.7D0)**2)+1.D0
              !                    c_h = -1.D0*exp(-(r(i,j,k)*4.0D0)**2)+1.D0
              !                    c_p = 11.D0*exp(-(r(i,j,k)*4.D0)**2)+1.D0
              !                    c_h = (10**(r(i,j,k)**5)-1.D0)/2.89485D0
              !                    c_p = (10**(r(i,j,k)**5)-1.D0)/2.89485D0

              Blagrangemultiplier_rhs(i,j,k) = (-c_h**2/c_p**2*Blagrangemultiplier(i,j,k) - c_h**2*( &
                   0.5*(0.5*(1.0/dx*( &
                   (Bxtilde(i,j,k) - Bxtilde(i-1,j,k)) + &
                   (Bxtilde(i,j-1,k) - Bxtilde(i-1,j-1,k)) + &
                   (Bxtilde(i,j,k-1) - Bxtilde(i-1,j,k-1)) + &
                   (Bxtilde(i,j-1,k-1) - Bxtilde(i-1,j-1,k-1)) &
                   ))) + &
                   0.5*(0.5*(1.0/dy*( &
                   (Bytilde(i,j,k) - Bytilde(i,j-1,k)) +  &
                   (Bytilde(i-1,j,k) - Bytilde(i-1,j-1,k)) +  &
                   (Bytilde(i,j,k-1) - Bytilde(i,j-1,k-1)) +  &
                   (Bytilde(i-1,j,k-1) - Bytilde(i-1,j-1,k-1))  &
                   ))) +  &
                   0.5*(0.5*(1.0/dz*(  &
                   (Bztilde(i,j,k) - Bztilde(i,j,k-1)) +  &
                   (Bztilde(i,j-1,k) - Bztilde(i,j-1,k-1)) +  &
                   (Bztilde(i-1,j,k) - Bztilde(i-1,j,k-1)) +  &
                   (Bztilde(i-1,j-1,k) - Bztilde(i-1,j-1,k-1))  &
                   )))))

           end if
        end do
     end do
  end do

  ! Empirically, we find that the following K-O strength seems to work best:
!  KO_Strength = 0.9D0
!  KO_Strength = 1.8D0
  KO_Strength = 0.9D0
  !$omp parallel do
  do k = cctk_nghostzones(3), ext(3)-cctk_nghostzones(3)
     do j = cctk_nghostzones(2), ext(2)-cctk_nghostzones(2)
        do i = cctk_nghostzones(1), ext(1)-cctk_nghostzones(1)

           Blagrangemultiplier_rhs(i,j,k) = Blagrangemultiplier_rhs(i,j,k) - KO_Strength / 16.D0 &
                * (+ (Blagrangemultiplier(i-2,j,k) - 4.D0*Blagrangemultiplier(i-1,j,k) + 6.D0*Blagrangemultiplier(i,j,k) - 4.D0*Blagrangemultiplier(i+1,j,k) + Blagrangemultiplier(i+2,j,k)) / dx &
                + (Blagrangemultiplier(i,j-2,k) - 4.D0*Blagrangemultiplier(i,j-1,k) + 6.D0*Blagrangemultiplier(i,j,k) - 4.D0*Blagrangemultiplier(i,j+1,k) + Blagrangemultiplier(i,j+2,k)) / dy &
                + (Blagrangemultiplier(i,j,k-2) - 4.D0*Blagrangemultiplier(i,j,k-1) + 6.D0*Blagrangemultiplier(i,j,k) - 4.D0*Blagrangemultiplier(i,j,k+1) + Blagrangemultiplier(i,j,k+2)) / dz)

           if(num_BHs.ge.1) then
              if((sqrt((x(i,j,k)-bh_posn_x(1))**2 + (y(i,j,k)-bh_posn_y(1))**2 + (z(i,j,k)-bh_posn_z(1))**2) .lt. min_BH_radius).or. &
                   (num_BHs.eq.2 .and. sqrt((x(i,j,k)-bh_posn_x(2))**2 + (y(i,j,k)-bh_posn_y(2))**2 + (z(i,j,k)-bh_posn_z(2))**2) .lt. min_BH_radius)) then
                 Bxtilde_or_Ax_rhs(i,j,k) = Bxtilde_or_Ax_rhs(i,j,k) - 1.D0 / 16.D0 &
                      * (+ (Bxtilde(i-2,j,k) - 4.D0*Bxtilde(i-1,j,k) + 6.D0*Bxtilde(i,j,k) - 4.D0*Bxtilde(i+1,j,k) + Bxtilde(i+2,j,k)) / dx &
                      + (Bxtilde(i,j-2,k) - 4.D0*Bxtilde(i,j-1,k) + 6.D0*Bxtilde(i,j,k) - 4.D0*Bxtilde(i,j+1,k) + Bxtilde(i,j+2,k)) / dy &
                      + (Bxtilde(i,j,k-2) - 4.D0*Bxtilde(i,j,k-1) + 6.D0*Bxtilde(i,j,k) - 4.D0*Bxtilde(i,j,k+1) + Bxtilde(i,j,k+2)) / dz)
                 Bytilde_or_Ay_rhs(i,j,k) = Bytilde_or_Ay_rhs(i,j,k) - 1.D0 / 16.D0 &
                      * (+ (Bytilde(i-2,j,k) - 4.D0*Bytilde(i-1,j,k) + 6.D0*Bytilde(i,j,k) - 4.D0*Bytilde(i+1,j,k) + Bytilde(i+2,j,k)) / dx &
                      + (Bytilde(i,j-2,k) - 4.D0*Bytilde(i,j-1,k) + 6.D0*Bytilde(i,j,k) - 4.D0*Bytilde(i,j+1,k) + Bytilde(i,j+2,k)) / dy &
                      + (Bytilde(i,j,k-2) - 4.D0*Bytilde(i,j,k-1) + 6.D0*Bytilde(i,j,k) - 4.D0*Bytilde(i,j,k+1) + Bytilde(i,j,k+2)) / dz)
                 Bztilde_or_Az_rhs(i,j,k) = Bztilde_or_Az_rhs(i,j,k) - 1.D0 / 16.D0 &
                      * (+ (Bztilde(i-2,j,k) - 4.D0*Bztilde(i-1,j,k) + 6.D0*Bztilde(i,j,k) - 4.D0*Bztilde(i+1,j,k) + Bztilde(i+2,j,k)) / dx &
                      + (Bztilde(i,j-2,k) - 4.D0*Bztilde(i,j-1,k) + 6.D0*Bztilde(i,j,k) - 4.D0*Bztilde(i,j+1,k) + Bztilde(i,j+2,k)) / dy &
                      + (Bztilde(i,j,k-2) - 4.D0*Bztilde(i,j,k-1) + 6.D0*Bztilde(i,j,k) - 4.D0*Bztilde(i,j,k+1) + Bztilde(i,j,k+2)) / dz)
              end if
              if(num_BHs.gt.2) then
                 write(*,*) "SORRY, Hyperbolic Divergence Cleaning DOES NOT CURRENTLY SUPPORT num_BHs.gt.2"
                 stop
              end if
           end if
        end do
     end do
  end do
  !$omp end parallel do

  if(1==0) then
     Blagrangemultiplier_rhs(1,:,:) = 0.D0
     Blagrangemultiplier_rhs(:,1,:) = 0.D0
     Blagrangemultiplier_rhs(:,:,1) = 0.D0
     Blagrangemultiplier_rhs(ext(1),:,:) = 0.D0
     Blagrangemultiplier_rhs(:,ext(2),:) = 0.D0
     Blagrangemultiplier_rhs(:,:,ext(3)) = 0.D0
  end if

end subroutine hyperbolic_divergence_cleaning
