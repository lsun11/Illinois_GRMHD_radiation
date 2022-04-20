#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
! Particle tracer for BHNS thorn.
! This routine has been optimized, and should not slow down
! your run with 1000 particles.
!
! Read the code & comments carefully to figure out the timestepping.
!-----------------------------------------------------------------------------
  
subroutine bhns_diagnostics_tracer_compute_u_0(CCTK_ARGUMENTS)
  
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  ! L variables used to compute B_poloidal & B_toroidal diagnostics
  real*8  :: horizdirn_x,horizdirn_y,horizdirn_z,max_horiz_radius
  integer :: i,j,k
  real*8  :: phiL,psi4L,psi6L,gxxL,gxyL,gxzL,gyyL,gyzL,gzzL,Bx_L,By_L,Bz_L,x_minus_xBH,y_minus_yBH,partial_phi_magnitude,e_hat_phix,e_hat_phiy,B_x_L,B_y_L,B_toroidal,B2,B_poloidal,Omega,smallb2L,lambda_MRI,dX,rho_starL

  if(cctk_iteration.eq.ITERATION_TO_BOOST_MAGNETIC_FIELDS) then
     Ax = Ax*mythbusters_boost_factor
     Ay = Ay*mythbusters_boost_factor
     Az = Az*mythbusters_boost_factor
     psi6phi = psi6phi*mythbusters_boost_factor

     Ax_p = Ax_p*mythbusters_boost_factor
     Ay_p = Ay_p*mythbusters_boost_factor
     Az_p = Az_p*mythbusters_boost_factor
     psi6phi_p = psi6phi_p*mythbusters_boost_factor

     Ax_p_p = Ax_p_p*mythbusters_boost_factor
     Ay_p_p = Ay_p_p*mythbusters_boost_factor
     Az_p_p = Az_p_p*mythbusters_boost_factor
     psi6phi_p_p = psi6phi_p_p*mythbusters_boost_factor
  end if

  if(MOD(cctk_iteration,out_every)==0 .and. em_evolve_enable==1 .and. cctk_iteration.ge.ITERATION_TO_output_MRI_wavelength) then

     call bhns_compute_b2_cpp(cctkGH,cctk_lsh, phi, lapm1, &
          shiftx,shifty,shiftz,vx,vy,vz,Bx,By,Bz, & 
          gxx, gxy, gxz, gyy, gyz, gzz, temp2)

     !     Get horizon radius in direction where it is likely to be MAXIMIZED:
     horizdirn_x = 0.D0
     horizdirn_y = 100000.D0
     horizdirn_z = 0.D0
     call get_ah_radius_in_dirn(cctkGH,horizdirn_x,horizdirn_y,horizdirn_z,max_horiz_radius)

     dX = CCTK_DELTA_SPACE(1)
     write(*,*) "dX,get_ah_radius_in_dirn = ",dX,max_horiz_radius

     !Here's how we define the toroidal & poloidal components of the magnetic field.
     !
     !  B_{toroidal} \approx |B^{\phi}| = |(-y B^x + x B^y) / (x^2 + y^2)|
     !
     !  B_{poloidal} \approx sqrt(B^2 - (B^{\phi})^2) = sqrt(\gamma_{ij} B^i B^j - \gamma_{\phi\phi} B^{\phi} B^{\phi}),
     !
     ! where 
     !
     !  \gamma_{\phi\phi} = y^2 \gamma_{xx} - 2xy \gamma_{xy} + x^2 \gamma_{yy},
     ! 
     ! and
     ! 
     ! \gamma_{ij} = e^{4 \phi} \tilde \gamma_{ij}, ( \tilde \gamma_{ij} is the evolved variable gij ).

     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              phiL = phi(i,j,k)
              psi4L = exp(4*phiL)
              psi6L = sqrt(psi4L)*psi4L
              gxxL = gxx(i,j,k)
              gxyL = gxy(i,j,k)
              gxzL = gxz(i,j,k)
              gyyL = gyy(i,j,k)
              gyzL = gyz(i,j,k)
              gzzL = gzz(i,j,k)

              Bx_L = Bx(i,j,k)
              By_L = By(i,j,k)
              Bz_L = Bz(i,j,k)

              x_minus_xBH = (X(i,j,k)-bh_posn_x(1))
              y_minus_yBH = (Y(i,j,k)-bh_posn_y(1))

              ! partial_phi_magnitude = | \partial_{\phi} | = length of phi basis vector
              partial_phi_magnitude = sqrt( psi4L * ( gxxL*(-y_minus_yBH)**2 + gyyL*(x_minus_xBH)**2 + 2.D0*gxyL*(x_minus_xBH)*(-y_minus_yBH) ) )

              ! e_hat_phix = e_{\hat \phi}^x = -(y-y_0)/| \partial_{\phi} | = -(y-y_0)/partial_phi_magnitude
              e_hat_phix = -y_minus_yBH / partial_phi_magnitude

              ! e_hat_phiy = e_{\hat \phi}^y =  (x-x_0)/| \partial_{\phi} | =  (x-x_0)/partial_phi_magnitude
              e_hat_phiy  = x_minus_xBH / partial_phi_magnitude

              ! B_x_L = \gamma_{xj} B^j
              B_x_L  = psi4L*( gxxL*Bx_L + gxyL*By_L + gxzL*Bz_L )
              
              ! B_y_L = \gamma_{yj} B^j
              B_y_L  = psi4L*( gxyL*Bx_L + gyyL*By_L + gyzL*Bz_L )

              ! B_toroidal = B^T = | B_x e_{\hat \phi}^x + B_y e_{\hat \phi}^y |
              B_toroidal = abs( B_x_L*e_hat_phix + B_y_L*e_hat_phiy )

              ! B2 = B^2 = \gamma_{ij} B^i B^j , as defined above.
              B2 = psi4L * (gxxL*Bx_L*Bx_L + 2.D0*gxyL*Bx_L*By_L + 2.D0*gxzL*Bx_L*Bz_L + gyyL*By_L*By_L + 2.D0*gyzL*By_L*Bz_L + gzzL*Bz_L*Bz_L)

              ! B_{poloidal} = sqrt(B^2 - B^{T}^2)
              B_poloidal = sqrt( B2 - B_toroidal**2 )

              if((B_toroidal)**2 .gt. B2+1D-20) then
                 write(*,*) "ERROR: Found B_toroidal > B^2 !",B_toroidal**2,B2
!                 stop
              end if
              
              if(rho_star(i,j,k).lt.0.D0) rho_star(i,j,k)=psi6L*rho_b_atm

              rho_starL = rho_star(i,j,k)

              temp10(i,j,k) = B_toroidal
              temp11(i,j,k) = B_poloidal

!              if(cctk_iteration.ge.ITERATION_TO_output_MRI_wavelength .and. MOD(cctk_iteration,out_every)==0) then
              ! Calculate MRI wavelength: McKinney, Tchekhovskoy, & Blandford, 2012, Eq 30
              ! \lambda_{MRI} \approx 2 pi |v_{\theta,A}| / |\Omega(r,\theta)|, where
              !      |v_{\theta,A}| = sqrt(b_{\theta} b^{\theta} / (b^2 + rho_0 h) ), b_{\theta} \approx B_poloidal
              Omega = ((X(i,j,k)-bh_posn_x(1))*vy(i,j,k)-(Y(i,j,k)-bh_posn_y(1))*vx(i,j,k)) / ((X(i,j,k)-bh_posn_x(1))**2+(Y(i,j,k)-bh_posn_y(1))**2) ! = \Omega(r,\theta)
              smallb2L = temp2(i,j,k)
              lambda_MRI = 2.D0*3.14159265358979323846D0*sqrt(B_poloidal*B_poloidal / (smallb2L + rho_b(i,j,k)*h(i,j,k)) ) / (abs(Omega)) ! \approx \lambda_{MRI}
              if(Omega==0.D0 .or. (smallb2L+ rho_b(i,j,k)*h(i,j,k))==0) lambda_MRI=0.D0

              P_thermal(i,j,k) = smallb2L
              temp2(i,j,k) = smallb2L / (2.D0*P(i,j,k)) * rho_starL
              temp3(i,j,k) = B_toroidal*rho_starL
              temp4(i,j,k) = B_poloidal*rho_starL
              temp5(i,j,k) = dX
              temp6(i,j,k) = P(i,j,k)/rho_b(i,j,k)/rho_b(i,j,k)
              temp7(i,j,k) = abs(Omega)
              temp8(i,j,k) = lambda_MRI
              temp9(i,j,k) = rho_starL
              if(rho_b(i,j,k).lt.1.D3*rho_b_atm) then
                 ! temp2 = b^2/(2*P)*rho_star, for rho>1e3 rho_atm
                 ! temp3 = B_toroidal*rho_star, for rho>1e3 rho_atm
                 ! temp4 = B_poloidal*rho_star, for rho>1e3 rho_atm
                 ! temp5 = dX, for rho>1e3 rho_atm
                 ! temp6 = K, for rho>1e3 rho_atm
                 ! temp7 = Omega, for rho>1e3 rho_atm
                 ! temp8 = lambda_MRI, for rho>1e3 rho_atm
                 ! temp9 = rho_star, for rho>1e3 rho_atm
                 ! temp10 = B_toroidal
                 ! temp11 = B_poloidal
                 ! MONOPOLE = rho_star
                 ! P_thermal = b^2

!                 P_thermal(i,j,k)=0.D0
                 temp2(i,j,k)=0.D0
                 temp3(i,j,k)=0.D0
                 temp4(i,j,k)=0.D0
!                 temp5(i,j,k)=0.D0
!                 temp6(i,j,k)=0.D0
!                 temp7(i,j,k)=0.D0
!                 temp8(i,j,k)=0.D0
                 temp9(i,j,k)=0.D0
                 temp10(i,j,k)=0.D0
                 temp11(i,j,k)=0.D0
              end if

              MONOPOLE(i,j,k) = rho_starL

 !             end if

              ! If inside the horizon or at very low density, turn off diagnostic
              !if(rho_b(i,j,k) .lt. 1.D0*rho_b_atm .or. sqrt( (X(i,j,k)-bh_posn_x(1))**2 + (Y(i,j,k)-bh_posn_y(1))**2 + (Z(i,j,k)-bh_posn_z(1))**2 ) .lt. max_horiz_radius) then
              !if(rho_b(i,j,k) .lt. 1.D3*rho_b_atm .or. sqrt( (X(i,j,k)-bh_posn_x(1))**2 + (Y(i,j,k)-bh_posn_y(1))**2 + (Z(i,j,k)-bh_posn_z(1))**2 ) .lt. max_horiz_radius) then
              if(sqrt( (X(i,j,k)-bh_posn_x(1))**2 + (Y(i,j,k)-bh_posn_y(1))**2 + (Z(i,j,k)-bh_posn_z(1))**2 ) .lt. max_horiz_radius) then
                 temp2(i,j,k) = 0.D0
                 temp3(i,j,k) = 0.D0
                 temp4(i,j,k) = 0.D0
!                 temp5(i,j,k) = 0.D0
                 temp6(i,j,k) = 0.D0
                 temp7(i,j,k) = 0.D0
                 temp8(i,j,k) = 0.D0
                 temp9(i,j,k) = 0.D0
                 temp10(i,j,k)= 0.D0
                 temp11(i,j,k)= 0.D0
                 MONOPOLE(i,j,k) = 0.D0
                 P_thermal(i,j,k) = 0.D0
              end if

              if(temp2(i,j,k).lt.0.D0) write(*,*) "BAD temp2!",x(i,j,k),y(i,j,k),z(i,j,k),temp2(i,j,k)
              if(temp3(i,j,k).lt.0.D0) write(*,*) "BAD temp3!",x(i,j,k),y(i,j,k),z(i,j,k),temp3(i,j,k)
              if(temp4(i,j,k).lt.0.D0) write(*,*) "BAD temp4!",x(i,j,k),y(i,j,k),z(i,j,k),temp4(i,j,k)
              if(temp5(i,j,k).lt.0.D0) write(*,*) "BAD temp5!",x(i,j,k),y(i,j,k),z(i,j,k),temp5(i,j,k)
              if(temp6(i,j,k).lt.0.D0) write(*,*) "BAD temp6!",x(i,j,k),y(i,j,k),z(i,j,k),temp6(i,j,k)
              if(temp7(i,j,k).lt.0.D0) write(*,*) "BAD temp7!",x(i,j,k),y(i,j,k),z(i,j,k),temp7(i,j,k)
              if(temp8(i,j,k).lt.0.D0) write(*,*) "BAD temp8!",x(i,j,k),y(i,j,k),z(i,j,k),temp8(i,j,k)
              if(temp9(i,j,k).lt.0.D0) write(*,*) "BAD temp9!",x(i,j,k),y(i,j,k),z(i,j,k),temp9(i,j,k)
              if(temp10(i,j,k).lt.0.D0) write(*,*) "BAD temp10!",x(i,j,k),y(i,j,k),z(i,j,k),temp10(i,j,k)
              if(temp11(i,j,k).lt.0.D0) write(*,*) "BAD temp11!",x(i,j,k),y(i,j,k),z(i,j,k),temp11(i,j,k)

              if(smallb2L .lt. 0.D0 .or. P(i,j,k) .lt.0.D0 .or. rho_star(i,j,k).lt.0.D0) then
                 write(*,*) "ERROR: Found b^2<0 or P<0 or rho_star<0 !",smallb2L,P(i,j,k),rho_star(i,j,k)
                 stop
              end if
              

           end do
        end do
     end do

  end if

  if(cctk_iteration.ge.bhns_particle_tracer_start) then
     
     if(MOD(cctk_iteration,particle_tracer_substep_every)==0 .or. cctk_iteration.eq.bhns_particle_tracer_start) then
        
        ! u_0 = u^0 * (-alpha^2 + gamma_ij beta^i (beta^j + v^j))
        temp1 = -1.D0 - u0*(-(lapm1+1.D0)**2 + &
             exp(4.D0*phi) * ( &
             gxx*shiftx*(shiftx+vx) + &
             2.D0*gxy*shiftx*(shifty+vy) + &
             2.D0*gxz*shiftx*(shiftz+vz) + &
             gyy*shifty*(shifty+vy) + &
             2.D0*gyz*shifty*(shiftz+vz) + &
             gzz*shiftz*(shiftz+vz)))

     end if
  end if

end subroutine bhns_diagnostics_tracer_compute_u_0
