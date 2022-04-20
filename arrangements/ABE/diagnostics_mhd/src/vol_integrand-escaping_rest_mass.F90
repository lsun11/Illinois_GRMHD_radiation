!----------------------------
! Set up rest mass integrand
!----------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine escaping_rest_mass_30M_integrand(CCTK_ARGUMENTS)
  implicit none

  !~~~~~> Input parameters:
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !~~~~~> Other variables:
  real*8                             :: dX, dY, dZ, u_0, vrad, rad
  integer                            :: i,j,k
  integer                            :: AXISYM
  parameter(AXISYM = 4)

  WhichIntegral = 119

!  write(*,*) "Start vol_integrand-escaping_rest_mass_30M.F90"
  
  if(MOD(cctk_iteration,Compute_VolIntegrands_Every)==0) then
     !-----------------------------------------------------------------------------
     ! Compute integrand
     !-----------------------------------------------------------------------------
     !Reset VolIntegrand quickly, using OpenMP
     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              VolIntegrand(i,j,k) = 0.D0
           end do
        end do
     end do
     !$omp end parallel do
     
     if (M_ADM.eq.0.0d0) then
        write(*,*) "YOU DIDN'T SET THE PARAMETER diagnostics_mhd::M_ADM IN YOUR PAR FILE AND THE ESCAPING MASS DIAGNOSTIC WILL BE GARBAGE"
     end if

     !$omp parallel do private (u_0,vrad,r)
     do k = cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
        do j = cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
           do i = cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
              ! Have to be in asymptotically flat regime for the diagnostic to be reliable
              ! So, here I set a radius cutoff at 30M. Note that M_ADM must be set in the par file
              rad=sqrt(X(i,j,k)*X(i,j,k)+Y(i,j,k)*Y(i,j,k)+Z(i,j,k)*Z(i,j,k))
              if (rad.gt.(radius_esc1*M_ADM)) then
                 ! compute U=-1-u_0, where u_0 = u^0 * (-alpha^2 + gamma_ij beta^i (beta^j + v^j))
                 u_0 = -1.D0 - u0(i,j,k)*(-(lapm1(i,j,k)+1.D0)**2 + &
                      exp(4.D0*phi(i,j,k)) * ( &
                      gxx(i,j,k)*shiftx(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) + &
                      2.D0*gxy(i,j,k)*shiftx(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) + &
                      2.D0*gxz(i,j,k)*shiftx(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) + &
                      gyy(i,j,k)*shifty(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) + &
                      2.D0*gyz(i,j,k)*shifty(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) + &
                      gzz(i,j,k)*shiftz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k))))
                 ! compute v^r = \partial_i r v^i = Sum_i x^i v^i/r
                 vrad = (X(i,j,k)*vx(i,j,k) + Y(i,j,k)*vy(i,j,k) + Z(i,j,k)*vz(i,j,k))/rad

                 
                 ! -1 - u_0 > 0 and v^r > 0 --> unbound, otherwise bound.                 
                 if ((u_0 .gt. 0.0d0).and.(vrad > 0.0d0)) then                    
                    VolIntegrand(i,j,k) = rho_star(i,j,k) 
                 end if
              end if
           end do ! i-loop
        end do ! j-loop
     end do ! k-loop
     !$omp end parallel do
     if (Symmetry == AXISYM) then
        VolIntegrand = X*VolIntegrand
     end if
  end if

!   write(*,*) "End vol_integrand-escaping_rest_mass_30M.F90"
end subroutine escaping_rest_mass_30M_integrand


subroutine escaping_rest_mass_50M_integrand(CCTK_ARGUMENTS)
  implicit none

  !~~~~~> Input parameters:
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !~~~~~> Other variables:
  real*8                             :: dX, dY, dZ, u_0, vrad, rad
  integer                            :: i,j,k
  integer                            :: AXISYM
  parameter(AXISYM = 4)

  WhichIntegral = 120

!  write(*,*) "Start vol_integrand-escaping_rest_mass_50M.F90"
  
  if(MOD(cctk_iteration,Compute_VolIntegrands_Every)==0) then
     !-----------------------------------------------------------------------------
     ! Compute integrand
     !-----------------------------------------------------------------------------
     !Reset VolIntegrand quickly, using OpenMP
     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              VolIntegrand(i,j,k) = 0.D0
           end do
        end do
     end do
     !$omp end parallel do
     
     if (M_ADM.eq.0.0d0) then
        write(*,*) "YOU DIDN'T SET THE PARAMETER diagnostics_mhd::M_ADM IN YOUR PAR FILE AND THE ESCAPING MASS DIAGNOSTIC WILL BE GARBAGE"
     end if

     !$omp parallel do private (u_0,vrad,r)
     do k = cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
        do j = cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
           do i = cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
              ! Have to be in asymptotically flat regime for the diagnostic to be reliable
              ! So, here I set a radius cutoff at 50M. Note that M_ADM must be set in the par file
              rad=sqrt(X(i,j,k)*X(i,j,k)+Y(i,j,k)*Y(i,j,k)+Z(i,j,k)*Z(i,j,k))
              if (rad.gt.(radius_esc2*M_ADM)) then
                 ! compute U=-1-u_0, where u_0 = u^0 * (-alpha^2 + gamma_ij beta^i (beta^j + v^j))
                 u_0 = -1.D0 - u0(i,j,k)*(-(lapm1(i,j,k)+1.D0)**2 + &
                      exp(4.D0*phi(i,j,k)) * ( &
                      gxx(i,j,k)*shiftx(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) + &
                      2.D0*gxy(i,j,k)*shiftx(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) + &
                      2.D0*gxz(i,j,k)*shiftx(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) + &
                      gyy(i,j,k)*shifty(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) + &
                      2.D0*gyz(i,j,k)*shifty(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) + &
                      gzz(i,j,k)*shiftz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k))))
                 ! compute v^r = \partial_i r v^i = Sum_i x^i v^i/r
                 vrad = (X(i,j,k)*vx(i,j,k) + Y(i,j,k)*vy(i,j,k) + Z(i,j,k)*vz(i,j,k))/rad

                 
                 ! -1 - u_0 > 0 and v^r > 0 --> unbound, otherwise bound.                 
                 if ((u_0 .gt. 0.0d0).and.(vrad > 0.0d0)) then                    
                    VolIntegrand(i,j,k) = rho_star(i,j,k) 
                 end if
              end if
           end do ! i-loop
        end do ! j-loop
     end do ! k-loop
     !$omp end parallel do
     if (Symmetry == AXISYM) then
        VolIntegrand = X*VolIntegrand
     end if
  end if

end subroutine escaping_rest_mass_50M_integrand

subroutine escaping_rest_mass_70M_integrand(CCTK_ARGUMENTS)
  implicit none

  !~~~~~> Input parameters:
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !~~~~~> Other variables:
  real*8                             :: dX, dY, dZ, u_0, vrad, rad
  integer                            :: i,j,k
  integer                            :: AXISYM
  parameter(AXISYM = 4)

  WhichIntegral = 121

  if(MOD(cctk_iteration,Compute_VolIntegrands_Every)==0) then
     !-----------------------------------------------------------------------------
     ! Compute integrand
     !-----------------------------------------------------------------------------
     !Reset VolIntegrand quickly, using OpenMP
     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              VolIntegrand(i,j,k) = 0.D0
           end do
        end do
     end do
     !$omp end parallel do
     
     if (M_ADM.eq.0.0d0) then
        write(*,*) "YOU DIDN'T SET THE PARAMETER diagnostics_mhd::M_ADM IN YOUR PAR FILE AND THE ESCAPING MASS DIAGNOSTIC WILL BE GARBAGE"
     end if

     !$omp parallel do private (u_0,vrad,r)
     do k = cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
        do j = cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
           do i = cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
              ! Have to be in asymptotically flat regime for the diagnostic to be reliable
              ! So, here I set a radius cutoff at 70M. Note that M_ADM must be set in the par file
              rad=sqrt(X(i,j,k)*X(i,j,k)+Y(i,j,k)*Y(i,j,k)+Z(i,j,k)*Z(i,j,k))
              if (rad.gt.(radius_esc3*M_ADM)) then
                 ! compute U=-1-u_0, where u_0 = u^0 * (-alpha^2 + gamma_ij beta^i (beta^j + v^j))
                 u_0 = -1.D0 - u0(i,j,k)*(-(lapm1(i,j,k)+1.D0)**2 + &
                      exp(4.D0*phi(i,j,k)) * ( &
                      gxx(i,j,k)*shiftx(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) + &
                      2.D0*gxy(i,j,k)*shiftx(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) + &
                      2.D0*gxz(i,j,k)*shiftx(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) + &
                      gyy(i,j,k)*shifty(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) + &
                      2.D0*gyz(i,j,k)*shifty(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) + &
                      gzz(i,j,k)*shiftz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k))))
                 ! compute v^r = \partial_i r v^i = Sum_i x^i v^i/r
                 vrad = (X(i,j,k)*vx(i,j,k) + Y(i,j,k)*vy(i,j,k) + Z(i,j,k)*vz(i,j,k))/rad

                 
                 ! -1 - u_0 > 0 and v^r > 0 --> unbound, otherwise bound.                 
                 if ((u_0 .gt. 0.0d0).and.(vrad > 0.0d0)) then                    
                    VolIntegrand(i,j,k) = rho_star(i,j,k) 
                 end if
              end if
           end do ! i-loop
        end do ! j-loop
     end do ! k-loop
     !$omp end parallel do
     if (Symmetry == AXISYM) then
        VolIntegrand = X*VolIntegrand
     end if
  end if

end subroutine escaping_rest_mass_70M_integrand

subroutine escaping_rest_mass_100M_integrand(CCTK_ARGUMENTS)
  implicit none

  !~~~~~> Input parameters:
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !~~~~~> Other variables:
  real*8                             :: dX, dY, dZ, u_0, vrad, rad
  integer                            :: i,j,k
  integer                            :: AXISYM
  parameter(AXISYM = 4)

  WhichIntegral = 122

  if(MOD(cctk_iteration,Compute_VolIntegrands_Every)==0) then
     !-----------------------------------------------------------------------------
     ! Compute integrand
     !-----------------------------------------------------------------------------
     !Reset VolIntegrand quickly, using OpenMP
     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              VolIntegrand(i,j,k) = 0.D0
           end do
        end do
     end do
     !$omp end parallel do
     
     if (M_ADM.eq.0.0d0) then
        write(*,*) "YOU DIDN'T SET THE PARAMETER diagnostics_mhd::M_ADM IN YOUR PAR FILE AND THE ESCAPING MASS DIAGNOSTIC WILL BE GARBAGE"
     end if

     !$omp parallel do private (u_0,vrad,r)
     do k = cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
        do j = cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
           do i = cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
              ! Have to be in asymptotically flat regime for the diagnostic to be reliable
              ! So, here I set a radius cutoff at 100M. Note that M_ADM must be set in the par file
              rad=sqrt(X(i,j,k)*X(i,j,k)+Y(i,j,k)*Y(i,j,k)+Z(i,j,k)*Z(i,j,k))
              if (rad.gt.(radius_esc4*M_ADM)) then
                 ! compute U=-1-u_0, where u_0 = u^0 * (-alpha^2 + gamma_ij beta^i (beta^j + v^j))
                 u_0 = -1.D0 - u0(i,j,k)*(-(lapm1(i,j,k)+1.D0)**2 + &
                      exp(4.D0*phi(i,j,k)) * ( &
                      gxx(i,j,k)*shiftx(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) + &
                      2.D0*gxy(i,j,k)*shiftx(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) + &
                      2.D0*gxz(i,j,k)*shiftx(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) + &
                      gyy(i,j,k)*shifty(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) + &
                      2.D0*gyz(i,j,k)*shifty(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) + &
                      gzz(i,j,k)*shiftz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k))))
                 ! compute v^r = \partial_i r v^i = Sum_i x^i v^i/r
                 vrad = (X(i,j,k)*vx(i,j,k) + Y(i,j,k)*vy(i,j,k) + Z(i,j,k)*vz(i,j,k))/rad

                 
                 ! U=-1 - u_0 > 0 and v^r > 0 --> unbound, otherwise bound.                 
                 if ((u_0 .gt. 0.0d0).and.(vrad > 0.0d0)) then                    
                    VolIntegrand(i,j,k) = rho_star(i,j,k) 
                 end if
              end if
           end do ! i-loop
        end do ! j-loop
     end do ! k-loop
     !$omp end parallel do
     if (Symmetry == AXISYM) then
        VolIntegrand = X*VolIntegrand
     end if
  end if

end subroutine escaping_rest_mass_100M_integrand
