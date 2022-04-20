!-------------------------------------------------------------------------
! Compute the analytic solution for the nonlinear Alfven wave
!-------------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

! Given time t, compute the analytic solution of the nonlinear Alfven wave 
! and store them to several gfs.
!
subroutine nonlinear_alfven_wave_analytic_time(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS 
  real*8 :: yx,yy,yz,zx,zy,zz, s4pi
  real*8 :: ut_l,ut_r,u1_r,u2_r,u3_r,B1_r,B2_r,B3_r,B1_l,B2_l,B3_l
  real*8 :: lambda,gamma,wf_ut_l,wf_u1_l,wf_u2_l,wf_u3_l
  real*8 :: wf_bt_l,wf_b1_l,wf_b2_l,wf_b3_l
  real*8 :: wf_bt_r,wf_b1_r,wf_b2_r,wf_b3_r
  real*8 :: wf_ut_r,wf_u1_r,wf_u2_r,wf_u3_r
  real*8 :: byc,bzc,chi_k,byz_rhs,theta_l, a11_k,a22_k,a12_k,xs
  real*8 :: ut_i,u1_i,u2_i,u3_i,ux_i,uy_i,uz_i
  real*8 :: B1_i,B2_i,B3_i,Bx_i,By_i,Bz_i
  integer :: i,j,k
  real*8, parameter :: pi = 3.1415926535897932d0
! 
  if (mod(cctk_iteration,out_every)==0) then 
     s4pi = 2.d0*sqrt(pi)

     if (abs(nx) .lt. 1.d-300) then
        ! In this case, the relationship between the coordinates in the wave frame
        ! and the lab frame is \hat{xs}=\hat{n}, \hat{ys}=-\hat{x},
        ! \hat{zs} = \hat{xs} cross \hat{ys}
        yx = -1.d0
        yy = 0.d0
        yz = 0.d0
        zx = 0.d0
        zy = -nz
        zz = ny
     else if (abs(ny) .lt. 1.d-300) then
        ! In this case, the relationship between the coordinates in the shock frame     ! and the lab frame is \hat{xs}=\hat{n}, \hat{ys}=\hat{y},
        ! \hat{zs} = \hat{xs} cross \hat{ys}
        yx = 0.d0
        yy = 1.d0
        yz = 0.d0
        zx = -nz
        zy = 0.d0
        zz = nx
     else 
        ! In this case, the relationship between the coordinates in the shock frame     ! and the lab frame is \hat{xs}=\hat{n}, \hat{zs}=\hat{z},
        ! \hat{ys} = \hat{z} cross \hat{n}
        yx = -ny
        yy = nx
        yz = 0.d0
        zx = 1.d0
        zy = 0.d0
        zz = 0.d0
     end if

     B1_l = B1os4pi_l*s4pi
     B2_l = B2os4pi_l*s4pi
     B3_l = 0.d0

     ! Compute right state and wave speed
     call nonlinear_alfven_wave_compute_right_state(u1_l,u2_l,u3_l,rho_b0,P0, &
                   B1os4pi_l,B2os4pi_l,rotations,npolywave,  &
                   ut_l,ut_r,u1_r,u2_r,u3_r,B1_r,B2_r,B3_r, &
                   lambda,gamma,wf_ut_l,wf_u1_l,wf_u2_l,wf_u3_l, &
                   wf_bt_l,wf_b1_l,wf_b2_l,wf_b3_l,  &
                   wf_ut_r,wf_u1_r,wf_u2_r,wf_u3_r, &
                   wf_bt_r,wf_b1_r,wf_b2_r,wf_b3_r,  &
                   byc,bzc,chi_k,byz_rhs,theta_l, a11_k,a22_k,a12_k)

     ! Now compute the analytic solution at time t 
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              xs = X(i,j,k)*nx + Y(i,j,k)*ny + Z(i,j,k)*nz - lambda*cctk_time
              call nonlinear_alfven_wave_analytic(xs,ut_l,u1_l,u2_l,u3_l,&
                     B1_l,B2_l,B3_l, ut_r,u1_r,u2_r,u3_r,B1_r,B2_r,B3_r, &
                     ut_i,u1_i,u2_i,u3_i, B1_i,B2_i,B3_i,wtrans,rotations,theta_l,chi_k, &
                     byz_rhs,a11_k,a12_k,a22_k,byc,bzc, wf_b1_l,wf_b2_l,wf_b3_l, &
                     wf_u1_l,wf_u2_l,wf_u3_l,lambda,gamma)
              call convert_v123_to_vxyz(u1_i,u2_i,u3_i,ux_i,uy_i,uz_i,nx,ny,nz)
              call convert_v123_to_vxyz(B1_i,B2_i,B3_i,Bx_i,By_i,Bz_i,nx,ny,nz)
              ! Store analytic u^i in vxr,vyr and vzr; analytic B^i in Bxr,Byr and Bzr
              vxr(i,j,k) = ux_i
              vyr(i,j,k) = uy_i
              vzr(i,j,k) = uz_i
              Bxr(i,j,k) = Bx_i
              Byr(i,j,k) = By_i
              Bzr(i,j,k) = Bz_i
              ! Store numerical u^i in vxl, vyl and vzl
              ut_i = 1.d0/sqrt(1.d0-vx(i,j,k)**2-vy(i,j,k)**2-vz(i,j,k)**2)
              vxl(i,j,k)=vx(i,j,k)*ut_i 
              vyl(i,j,k)=vy(i,j,k)*ut_i
              vzl(i,j,k)=vz(i,j,k)*ut_i
           end do
        end do
     end do
   end if
end subroutine nonlinear_alfven_wave_analytic_time

!
! Given the left side of the wave u^i_l=(u1_l,u2_l,u3_l), rho_b, P,
! B^i_l = (B1_l,B2_l,0), the wave's width wtrans and "rotations"; compute the 
! right side of the wave u^i_r=(u1_r,u2_r,u3_r), B^i_r = (B1_r,B2_r,B3_r) 
! at the lab frame, as well as quantities at the wave frame (frame at which the 
! wave is stationary), wave speed lambda and a bunch of parameters necessary to 
! compute everything anywhere.
! 
subroutine nonlinear_alfven_wave_compute_right_state(u1_l,u2_l,u3_l,rho_b0,P0, & 
		B1os4pi_l,B2os4pi_l,rotations,npolywave,  & 
                ut_l,ut_r,u1_r,u2_r,u3_r,B1_r,B2_r,B3_r, & 
                lambda,gamma,wf_ut_l,wf_u1_l,wf_u2_l,wf_u3_l, & 
                wf_bt_l,wf_b1_l,wf_b2_l,wf_b3_l,  & 
                wf_ut_r,wf_u1_r,wf_u2_r,wf_u3_r, &
		wf_bt_r,wf_b1_r,wf_b2_r,wf_b3_r,  &
                byc,bzc,chi_k,byz_rhs,theta_l, a11_k,a22_k,a12_k)
  implicit none
  ! Input vars:
  real*8 :: u1_l,u2_l,u3_l,B1os4pi_l,B2os4pi_l,rotations,rho_b0,P0,npolywave
  ! Output vars:
  real*8 :: ut_l,ut_r,u1_r,u2_r,u3_r,B1_r,B2_r,B3_r,byz_rhs,theta_l
  real*8 :: lambda,gamma,wf_ut_l,wf_u1_l,wf_u2_l,wf_u3_l
  real*8 :: wf_bt_l,wf_b1_l,wf_b2_l,wf_b3_l, byc,bzc,chi_k
  real*8 :: wf_ut_r,wf_u1_r,wf_u2_r,wf_u3_r,wf_bt_r,wf_b1_r,wf_b2_r,wf_b3_r
  ! Other vars:
  real*8 :: B1_l,B2_l,B3_l,sbt_l,sb1_l,sb2_l,sb3_l,sbs_l
  real*8 :: sbt_r,sb1_r,sb2_r,sb3_r,sbs_r
  real*8 :: h_l,sq_eps,s4pi,ay_k,az_k,c_k,d_k
  real*8 :: a11_k,a22_k,a33_k,a12_k,a13_k,a23_k,capD_k,coD
  real*8 :: thet,byz
  real*8, parameter :: pi = 3.1415926535897932d0
!
  s4pi = sqrt(4.d0*pi)
  ut_l = sqrt(1.d0 + u1_l**2 + u2_l**2 + u3_l**2)
  B1_l = B1os4pi_l * s4pi
  B2_l = B2os4pi_l * s4pi
  B3_l = 0.d0

  ! Compute b^mu
  sbt_l = u1_l*B1os4pi_l + u2_l*B2os4pi_l
  sb1_l = (B1os4pi_l + sbt_l*u1_l)/ut_l
  sb2_l = (B2os4pi_l + sbt_l*u2_l)/ut_l
  sb3_l = sbt_l*u3_l/ut_l
  sbs_l = sb1_l**2 + sb2_l**2 + sb3_l**2 - sbt_l**2

  h_l = 1.d0 + ( npolywave + 1.d0 ) * P0 / rho_b0

  ! Compute wave speed lambda
  sq_eps = sqrt(rho_b0 * h_l + sbs_l)
  if (u1_l .gt. 0.d0) then
     lambda = (sb1_l/sq_eps - u1_l)/(sbt_l/sq_eps - ut_l)
  else
     lambda = (sb1_l/sq_eps + u1_l)/(sbt_l/sq_eps + ut_l)
  end if

  ! boost to wavefront frame
  gamma = 1.d0/sqrt(1.d0-lambda**2)
  wf_ut_l = gamma*(ut_l - lambda*u1_l)
  wf_u1_l = gamma*(u1_l - lambda*ut_l)
  wf_u2_l = u2_l
  wf_u3_l = u3_l
  wf_bt_l = gamma*(sbt_l - lambda*sb1_l)
  wf_b1_l = gamma*(sb1_l - lambda*sbt_l)
  wf_b2_l = sb2_l
  wf_b3_l = sb3_l

  ! compute Komissarov variables
  chi_k = wf_u1_l / wf_b1_l
  ay_k = (wf_u2_l - chi_k*wf_b2_l)/(wf_ut_l - chi_k*wf_bt_l)
  az_k = (wf_u3_l - chi_k*wf_b3_l)/(wf_ut_l - chi_k*wf_bt_l)
  c_k = chi_k*sbs_l/(wf_ut_l - chi_k*wf_bt_l)
  d_k = sbs_l - wf_b1_l*wf_b1_l
  a11_k = 1.d0 - ay_k*ay_k
  a22_k = 1.d0 - az_k*az_k
  a33_k = -(c_k*c_k + d_k)
  a12_k = -ay_k*az_k
  a13_k = -c_k*ay_k
  a23_k = -c_k*az_k
  capD_k = a11_k*a22_k - a12_k*a12_k
  coD = c_k/capD_k
  ! byz_rhs = d + c^2/D
  byz_rhs = d_k + c_k * coD
  byc = coD * ay_k
  bzc = coD * az_k
  theta_l = atan2(wf_b3_l-bzc,wf_b2_l-byc)

  ! Compute the right side of the wave
  thet = theta_l + rotations * 2.d0*pi
  byz = sqrt( byz_rhs/(a11_k*cos(thet)*cos(thet) &
                         + a22_k*sin(thet)*sin(thet) &
                         + 2.d0*a12_k*cos(thet)*sin(thet)) )
  wf_b1_r = wf_b1_l
  wf_b2_r = byc + byz*cos(thet)
  wf_b3_r = bzc + byz*sin(thet)
  ! [u] = \chi [b]
  wf_u1_r = wf_u1_l
  wf_u2_r = wf_u2_l + chi_k*(wf_b2_r - wf_b2_l)
  wf_u3_r = wf_u3_l + chi_k*(wf_b3_r - wf_b3_l)
  ! u.u = -1
  wf_ut_r = sqrt( 1.d0 + wf_u1_r**2 + wf_u2_r**2 + wf_u3_r**2 )
  ! b.u = 0
  wf_bt_r = (wf_b1_r*wf_u1_r + wf_b2_r*wf_u2_r + wf_b3_r*wf_u3_r)/wf_ut_r
  ! boost back to lab frame
  ut_r = gamma*(wf_ut_r + lambda*wf_u1_r)
  u1_r = gamma*(wf_u1_r + lambda*wf_ut_r)
  u2_r = wf_u2_r
  u3_r = wf_u3_r
  sbt_r = gamma*(wf_bt_r + lambda*wf_b1_r)
  sb1_r = gamma*(wf_b1_r + lambda*wf_bt_r)
  sb2_r = wf_b2_r
  sb3_r = wf_b3_r
  sbs_r = sb1_r*sb1_r + sb2_r*sb2_r + sb3_r*sb3_r - sbt_r*sbt_r
  B1_r = (ut_r*sb1_r - sbt_r*u1_r)*s4pi
  B2_r = (ut_r*sb2_r - sbt_r*u2_r)*s4pi
  B3_r = (ut_r*sb3_r - sbt_r*u3_r)*s4pi
end subroutine nonlinear_alfven_wave_compute_right_state

!
! Given position x1, compute u^mu=(u0,u1,u2,u3) and B^i=(B1,B2,B3).
!
subroutine nonlinear_alfven_wave_analytic(x1,ut_l,u1_l,u2_l,u3_l,& 
          B1_l,B2_l,B3_l, ut_r,u1_r,u2_r,u3_r,B1_r,B2_r,B3_r, &
          u0,u1,u2,u3, B1,B2,B3,wtrans,rotations,theta_l,chi_k, & 
          byz_rhs,a11_k,a12_k,a22_k,byc,bzc, wf_b1_l,wf_b2_l,wf_b3_l, & 
          wf_u1_l,wf_u2_l,wf_u3_l,lambda,gamma)
  implicit none
  ! Input vars:
  real*8 :: x1,ut_l,u1_l,u2_l,u3_l,B1_l,B2_l,B3_l
  real*8 :: ut_r,u1_r,u2_r,u3_r,B1_r,B2_r,B3_r,lambda,gamma
  real*8 :: wtrans,rotations,theta_l,chi_k,byz_rhs,a11_k,a12_k,a22_k
  real*8 :: byc,bzc, wf_b1_l,wf_b2_l,wf_b3_l, wf_u1_l,wf_u2_l,wf_u3_l
  ! Output vars:
  real*8 :: u0,u1,u2,u3,B1,B2,B3
  ! Other vars:
  real*8 :: thet,byz,wf_b1,wf_b2,wf_b3,wf_u1,wf_u2,wf_u3,wf_ut,wf_bt
  real*8 :: sb0,sb1,sb2,sb3,s4pi
  real*8, parameter :: pi = 3.1415926535897932d0
!
  s4pi = sqrt(4.d0*pi)
  if (x1 .le. -0.5d0*wtrans) then 
     ! Left side of the wave
     u0 = ut_l
     u1 = u1_l
     u2 = u2_l
     u3 = u3_l    
     B1 = B1_l
     B2 = B2_l
     B3 = B3_l
  else if (x1 .ge. 0.5d0*wtrans) then 
     ! Right side of the wave
     u0 = ut_r
     u1 = u1_r
     u2 = u2_r
     u3 = u3_r
     B1 = B1_r
     B2 = B2_r
     B3 = B3_r
  else 
     ! "In" side of the wave
     thet = theta_l + rotations * 2.d0*pi * &
                  sin( (x1 + wtrans*0.5d0)/wtrans*pi*0.5d0 )**2
     byz = sqrt( byz_rhs/(a11_k*cos(thet)*cos(thet) &
                            + a22_k*sin(thet)*sin(thet) &
                            + 2.d0*a12_k*cos(thet)*sin(thet)) )
     wf_b1 = wf_b1_l
     wf_b2 = byc + byz*cos(thet)
     wf_b3 = bzc + byz*sin(thet)
     ! [u] = \chi [b]
     wf_u1 = wf_u1_l
     wf_u2 = wf_u2_l + chi_k*(wf_b2 - wf_b2_l)
     wf_u3 = wf_u3_l + chi_k*(wf_b3 - wf_b3_l)
     ! u.u = -1
     wf_ut = sqrt( 1.d0 + wf_u1**2 + wf_u2**2 + wf_u3**2 )
     ! b.u = 0
     wf_bt = (wf_b1*wf_u1 + wf_b2*wf_u2 + wf_b3*wf_u3)/wf_ut
     ! boost back to lab frame
     u0 = gamma*(wf_ut + lambda*wf_u1)
     u1 = gamma*(wf_u1 + lambda*wf_ut)
     u2 = wf_u2
     u3 = wf_u3
     sb0 = gamma*(wf_bt + lambda*wf_b1)
     sb1 = gamma*(wf_b1 + lambda*wf_bt)
     sb2 = wf_b2
     sb3 = wf_b3
     B1 = (u0*sb1 - sb0*u1)*s4pi
     B2 = (u0*sb2 - sb0*u2)*s4pi
     B3 = (u0*sb3 - sb0*u3)*s4pi
  end if
end subroutine nonlinear_alfven_wave_analytic 
