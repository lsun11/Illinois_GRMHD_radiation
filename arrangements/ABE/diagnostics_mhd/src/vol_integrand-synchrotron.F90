!----------------------------
! Set up rest mass integrand
!----------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine synchrotron_integrand(CCTK_ARGUMENTS)
  implicit none

  !~~~~~> Input parameters:
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  INTERFACE
     SUBROUTINE synch_temp_dep(theta,t_dep)
       IMPLICIT NONE
       real*8 :: t_dep,theta
     END SUBROUTINE synch_temp_dep
  END INTERFACE
  !~~~~~> Other variables:
  real*8                             :: dX, dY, dZ,rr1,rr2,rr3
  integer                            :: i,j,k
  integer                            :: AXISYM
  integer                            :: found_horizon1,found_horizon2,found_horizon3,good_point
  real*8                             :: rad_est_1,rad_est_2,rad_est_3
  real*8                             :: q_synch,n_tilde,n_inf,theta,mb_o_me,GM_o_c2
  real*8                             :: lapse,sqrtgamma
  real*8                             :: xh1,xh2,xh3,yh1,yh2,yh3,zh1,zh2,zh3
  real*8                             :: x1,x2,x3,y1,y2,y3,z1,z2,z3  
  real*8                             :: Const,beta
  real*8                             :: e,c,me_c2
    
  parameter(AXISYM = 4)

  call estimate_ah_radius(1,Symmetry,found_horizon1,xh1,yh1,zh1,rad_est_1)
  
  write(*,*) "inside vol_integrand_synchrotron"
  if (rot_metric .ne. 1) then
     if (num_BHs .gt. 1) then
        call estimate_ah_radius(2,Symmetry,found_horizon2,xh2,yh2,zh2,rad_est_2) 
     else 
        !this is specifically for a single schwarz BH at origin.  This should be made more general.
        found_horizon2=0
        xh1=0.d0
        yh1=0.d0
        zh1=0.d0
        rad_est_1=2.d0
     endif
     if (num_BHs .gt. 2) then
        call estimate_ah_radius(3,Symmetry,found_horizon3,xh3,yh3,zh3,rad_est_3) 
     else
        found_horizon3=0
     endif
  else	
     xh1 = xbh1_initial*cos(binary_orb_freq*cctk_time)
     yh1 = xbh1_initial*sin(binary_orb_freq*cctk_time)
     zh1 = 0.d0
     xh2 = xbh2_initial*cos(binary_orb_freq*cctk_time)
     yh2 = xbh2_initial*sin(binary_orb_freq*cctk_time)
     zh2 = 0.d0
     xh3 = 0.d0
     yh3 = 0.d0
     zh3 = 0.d0
     rad_est_1 = rah1_initial
     rad_est_2 = rah2_initial
     rad_est_3 = 0.d0
     found_horizon1 = 1
     found_horizon2 = 1
     found_horizon3 = 0
  endif
  WhichIntegral = 113
  
  mb_o_me = 1.83868d3
  n_inf = 1.d0
  q_synch = 0.d0
  GM_o_c2 = 1.4766d5
  
  Const=2.151889
  beta=10.d0
  me_c2 = 8.187111168d-7   !rest mass energy of electron in ergs
  e = 4.802d-10 !electron charge in esu
  c = 3.0d10 !speed of light in cm/sec

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
     
     !!  !$omp parallel do
     do k = cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
        do j = cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
           do i = cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
              if (found_horizon1 .eq. 1) then
                 x1 = X(i,j,k)-xh1
                 y1 = Y(i,j,k)-yh1
                 z1 = Z(i,j,k)-zh1
              else 
                 x1=0.d0
                 y1=0.d0
                 z1=0.d0
              endif
              
              if (found_horizon2 .eq. 1) then
                 x2 = X(i,j,k)-xh2
                 y2 = Y(i,j,k)-yh2 
                 z2 = Z(i,j,k)-zh2
              else
                 x2=0.d0
                 y2=0.d0
                 z2=0.d0   
              endif
              
              if (found_horizon3 .eq. 1) then 
                 x3 = X(i,j,k)-xh3
                 y3 = Y(i,j,k)-yh3
                 z3 = Z(i,j,k)-zh3
              else
                 x3=0.d0
                 y3=0.d0
                 z3=0.d0     
              endif
              
              rr1=sqrt(x1*x1+y1*y1+z1*z1)
              rr2=sqrt(x2*x2+y2*y2+z2*z2)
              rr3=sqrt(x3*x3+y3*y3+z3*z3)
              
              if (((rr1 .gt. inner_lum_rad_ratio*rad_est_1) &
                   .or. (found_horizon1 .eq. 0)) .and. &
                   ((rr2 .gt. inner_lum_rad_ratio*rad_est_2) &
                   .or. (found_horizon2 .eq. 0)) .and. &
                   ((rr3 .gt. inner_lum_rad_ratio*rad_est_3) &
                   .or. (found_horizon3 .eq. 0)) .and. &
                   (sqrt(X(i,j,k)*X(i,j,k)&
                   +Y(i,j,k)*Y(i,j,k)&
                   +Z(i,j,k)*Z(i,j,k)) .lt. lum_outer_rad)) then
                 lapse = lapm1(i,j,k)+1.d0
                 sqrtgamma = psi(i,j,k)**6
                 !     replace this with correct VolIntegrand
                 theta=0.5d0*P(i,j,k)/rho_b(i,j,k)*mb_o_me
                 n_tilde=rho_b(i,j,k)
                 !     note that I assume rho_b=1 at infinity so rho_b = rho_tilde
                 call synch_temp_dependence(theta,q_synch)
		 !write(*,*) "theta,q_synch: ",theta,q_synch
                 q_synch = q_synch * (4.0505*Const*108.d0*(n_tilde*n_inf)**2*e**4*c)/(sqrt(3.d0)*beta*me_c2)
                 if (theta .lt. 1.d0) q_synch=0.d0
                 !     VolIntegrand(i,j,k) = lapse*sqrtgamma*(q_ei+q_ee) 
                 VolIntegrand(i,j,k) = sqrtgamma*q_synch*GM_o_c2**3
	      endif
           end do               ! i-loop
        end do                  ! j-loop
     end do                    ! k-loop
     !     !  !$omp end parallel do
     if (Symmetry == AXISYM) then
        VolIntegrand = X*VolIntegrand
     end if
  end if
end subroutine synchrotron_integrand
