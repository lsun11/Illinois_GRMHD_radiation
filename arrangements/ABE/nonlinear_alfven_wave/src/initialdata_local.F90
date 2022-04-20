!-------------------------------------------------------------------
! Here we set up all required variables for shock tests, including:
!  emfields, BSSN variables, primitives, etc.
!-------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine nonlinear_alfven_wave_initialdata_local(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  real*8                :: dT,dX,dY,dZ,P_max,rho_max,tau_max,rhos_max,rho_fail_max_step,M_fail_step
  real*8                :: xmin,ymin,zmin,xmax,ymax,zmax
  integer               :: proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax
  integer               :: ierr,index,handle,dummy,glob_imax,glob_jmax,glob_kmax
  CCTK_REAL             :: reduction_value

  integer               :: i,j,k
  real*8                :: norm,h_l,ut_l,hm1_l,tau_l,tau_r
  real*8                :: ux_l,ux_r,uy_l,uy_r,uz_l,uz_r
  real*8                :: ut_r,u1_r,u2_r,u3_r
  real*8                :: xs,ys,zs,yx,yy,yz,zx,zy,zz
  real*8                :: Bx_l,Bx_r,By_l,By_r,Bz_l,Bz_r
  real*8                :: B1_l,B2_l,B3_l,B1_r,B2_r,B3_r
  real*8                :: gamma,lambda
  real*8                :: wf_ut_l,wf_u1_l,wf_u2_l,wf_u3_l
  real*8                :: wf_bt_l,wf_b1_l,wf_b2_l,wf_b3_l
  real*8                :: wf_ut_r,wf_u1_r,wf_u2_r,wf_u3_r
  real*8                :: wf_bt_r,wf_b1_r,wf_b2_r,wf_b3_r
  real*8                :: wf_ut,wf_u1,wf_u2,wf_u3
  real*8                :: wf_bt,wf_b1,wf_b2,wf_b3
  real*8 		:: ut_i,u1_i,u2_i,u3_i,sbt_i,sb1_i,sb2_i,sb3_i,sbs_i
  real*8                :: B1_i,B2_i,B3_i,ux_i,uy_i,uz_i
  real*8                :: sbt_l,sbx_l,sby_l,sbz_l,sbt_r,sbx_r,sby_r,sbz_r
  real*8                :: sbx_i,sby_i,sbz_i, Bx_i,By_i,Bz_i
  real*8                :: byc,bzc,bpy,bpz,byz,theta_l,thet
  real*8 		:: A1,A2,A3
  real*8                :: chi_k, ay_k, az_k, c_k, d_k, coD,byz_rhs
  real*8                :: a11_k,a22_k,a33_k,a12_k,a13_k,a23_k,capD_k
  real*8, parameter :: pi = 3.1415926535897932d0, f1o6 = 0.1666666666666666666666667d0
  real*8 :: s4pi,sb1_l,sb1_r,sb2_l,sb2_r,sb3_r,sb3_l,sbs_l,sbs_r
  real*8 :: dx1, B3_im1,B3_imh, int_B3dx_i
  integer, parameter :: num = 1000
  real*8, dimension(num) :: int_B3dx

  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  integer :: ONE,ZERO
  parameter(ONE = 1.D0, ZERO = 0.D0)
!

  if(primitives_solver.eq.0 .or. enable_shocktest_primitive_mode.ne.1) then
     write(*,*) "ERROR.  YOU MUST SET primitives_solver=1 or 2 AND enable_shocktest_primitive_mode=1 FOR THIS TEST"
     stop
  end if
 
  ! sqrt(4 pi) 
  s4pi = 2.d0*sqrt(pi)

  ext = cctk_lsh

  if(cctk_lsh(1).le.2*cctk_nghostzones(1)) then
     write(*,*) "ERROR: MUST SET NX>NGHOSTZONESX*2!"
     stop
  else if(cctk_lsh(2).le.2*cctk_nghostzones(2)) then
     write(*,*) "ERROR: MUST SET NY>NGHOSTZONESY*2!"
     stop
  else if(cctk_lsh(3).le.2*cctk_nghostzones(3)) then
     write(*,*) "ERROR: MUST SET NZ>NGHOSTZONESZ*2!"
     stop
  end if

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ! just to make sure, normalize n:
  !
  norm = sqrt(nx*nx + ny*ny + nz*nz)
  nx = nx/norm
  ny = ny/norm
  nz = nz/norm

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
     ! In this case, the relationship between the coordinates in the shock frame
     ! and the lab frame is \hat{xs}=\hat{n}, \hat{ys}=\hat{y},
     ! \hat{zs} = \hat{xs} cross \hat{ys}
     yx = 0.d0
     yy = 1.d0
     yz = 0.d0
     zx = -nz
     zy = 0.d0
     zz = nx
  else if (abs(nz) .lt. 1.d-300) then 
     ! In this case, the relationship between the coordinates in the shock frame
     ! and the lab frame is \hat{xs}=\hat{n}, \hat{zs}=\hat{z},
     ! \hat{ys} = \hat{z} cross \hat{n}
     yx = -ny
     yy = nx
     yz = 0.d0
     zx = 1.d0
     zy = 0.d0
     zz = 0.d0
  else
     write(*,*) 'At least one of the nx, ny or ny must be set to 0'
     stop
  end if

  B1_l = B1os4pi_l * s4pi
  B2_l = B2os4pi_l * s4pi
  B3_l = 0.d0

  ! Compute b^mu
  sbt_l = u1_l*B1os4pi_l + u2_l*B2os4pi_l 
  sb1_l = (B1os4pi_l + sbt_l*u1_l)/ut_l
  sb2_l = (B2os4pi_l + sbt_l*u2_l)/ut_l
  sb3_l = sbt_l*u3_l/ut_l
  sbs_l = sb1_l**2 + sb2_l**2 + sb3_l**2 - sbt_l**2

  hm1_l = ( npolywave + 1.d0 ) * P0 / rho_b0
  h_l = 1.d0 + hm1_l

  tau_l = rho_b0*ut_l*(h_l*ut_l - 1.d0) - P0 + sbs_l*(ut_l*ut_l-0.5d0) - sbt_l**2

  ! Compute the right side of the wave, wave speed, and a bunch of 
  ! auxiliary parameters used to compute the wave anywhere
  !
  call nonlinear_alfven_wave_compute_right_state(u1_l,u2_l,u3_l,rho_b0,P0, &
                B1os4pi_l,B2os4pi_l,rotations,npolywave,  &
                ut_l,ut_r,u1_r,u2_r,u3_r,B1_r,B2_r,B3_r, &
                lambda,gamma,wf_ut_l,wf_u1_l,wf_u2_l,wf_u3_l, &
                wf_bt_l,wf_b1_l,wf_b2_l,wf_b3_l,  &
                wf_ut_r,wf_u1_r,wf_u2_r,wf_u3_r, &
                wf_bt_r,wf_b1_r,wf_b2_r,wf_b3_r,  &
                byc,bzc,chi_k,byz_rhs,theta_l, a11_k,a22_k,a12_k)

  sbt_r = gamma*(wf_bt_r + lambda*wf_b1_r)
  sb1_r = gamma*(wf_b1_r + lambda*wf_bt_r)
  sb2_r = wf_b2_r
  sb3_r = wf_b3_r
  sbs_r = sb1_r*sb1_r + sb2_r*sb2_r + sb3_r*sb3_r - sbt_r*sbt_r
  ! Note that h_r = h_l since P and rhob are constant everywhere
  tau_r = rho_b0*ut_r*(h_l*ut_r - 1.d0) - P0 + sbs_r*(ut_r*ut_r-0.5d0) - sbt_r**2

  call convert_v123_to_vxyz(u1_l,u2_l,u3_l,ux_l,uy_l,uz_l,nx,ny,nz)
  call convert_v123_to_vxyz(u1_r,u2_r,u3_r,ux_r,uy_r,uz_r,nx,ny,nz)
  call convert_v123_to_vxyz(sb1_l,sb2_l,sb3_l,sbx_l,sby_l,sbz_l,nx,ny,nz)
  call convert_v123_to_vxyz(sb1_r,sb2_r,sb3_r,sbx_r,sby_r,sbz_r,nx,ny,nz)
  call convert_v123_to_vxyz(B1_l,B2_l,B3_l,Bx_l,By_l,Bz_l,nx,ny,nz)
  call convert_v123_to_vxyz(B1_r,B2_r,B3_r,Bx_r,By_r,Bz_r,nx,ny,nz)

  ! Compute Integrate[ B3(x'), {x',-infinity,x1}] for -wtrans/2 < x1 < wtrans/2
  dx1 = wtrans/dble(num-1)
  if ( min(dX,dY,dZ)/dx1 .lt. 1.d0) then 
     write(*,*) 'You need to increase num to ensure accuracy!'
     stop
  end if
  int_B3dx(1) = 0.d0
  ! Integration by Simpson's rule
  B3_im1 = 0.d0
  do i=2,num
     xs = (i-1)*dx1 - 0.5d0*wtrans
     call nonlinear_alfven_wave_analytic(xs-0.5d0*dx1,ut_l,u1_l,u2_l,u3_l,&
          B1_l,B2_l,B3_l, ut_r,u1_r,u2_r,u3_r,B1_r,B2_r,B3_r, &
          ut_i,u1_i,u2_i,u3_i, B1_i,B2_i,B3_imh,wtrans,rotations,theta_l,chi_k, &
          byz_rhs,a11_k,a12_k,a22_k,byc,bzc, wf_b1_l,wf_b2_l,wf_b3_l, &
          wf_u1_l,wf_u2_l,wf_u3_l,lambda,gamma)
     call nonlinear_alfven_wave_analytic(xs,ut_l,u1_l,u2_l,u3_l,&
          B1_l,B2_l,B3_l, ut_r,u1_r,u2_r,u3_r,B1_r,B2_r,B3_r, &
          ut_i,u1_i,u2_i,u3_i, B1_i,B2_i,B3_i,wtrans,rotations,theta_l,chi_k, &
          byz_rhs,a11_k,a12_k,a22_k,byc,bzc, wf_b1_l,wf_b2_l,wf_b3_l, &
          wf_u1_l,wf_u2_l,wf_u3_l,lambda,gamma)
     int_B3dx(i) = int_B3dx(i-1) + dx1*(B3_im1 + 4.d0*B3_imh + B3_i)*f1o6
     B3_im1 = B3_i
  end do

  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
           P(i,j,k) = P0
 	   rho_b(i,j,k) = rho_b0
	   xs = X(i,j,k)*nx + Y(i,j,k)*ny + Z(i,j,k)*nz
           if (xs .le. -0.5d0*wtrans) then
              !
              ! "left" side of partition
              !
              rho_star(i,j,k) = rho_b0 * ut_l
              tau(i,j,k)   = tau_l
              h(i,j,k) = h_l
              w(i,j,k)    = rho_star(i,j,k) * ut_l
              st_x(i,j,k)  = rho_star(i,j,k) *h_l*ux_l
              st_y(i,j,k)  = rho_star(i,j,k) *h_l*uy_l
              st_z(i,j,k)  = rho_star(i,j,k) *h_l*uz_l
	      mhd_st_x(i,j,k) = st_x(i,j,k) + ut_l*sbs_l*ux_l - sbt_l*sbx_l
              mhd_st_y(i,j,k) = st_y(i,j,k) + ut_l*sbs_l*uy_l - sbt_l*sby_l
              mhd_st_z(i,j,k) = st_z(i,j,k) + ut_l*sbs_l*uz_l - sbt_l*sbz_l
              vx(i,j,k)   = ux_l/ut_l
              vy(i,j,k)   = uy_l/ut_l
              vz(i,j,k)   = uz_l/ut_l
	      u0(i,j,k)   = ut_l
   	      Bx(i,j,k)   = Bx_l
              By(i,j,k)   = By_l
              Bz(i,j,k)   = Bz_l
	      Bxtilde(i,j,k) = Bx_l
              Bytilde(i,j,k) = By_l
              Bztilde(i,j,k) = Bz_l
	      sbt(i,j,k) = sbt_l
	      sbx(i,j,k) = sbx_l
	      sby(i,j,k) = sby_l
	      sbz(i,j,k) = sbz_l
           else if (xs .ge. 0.5d0*wtrans) then
              !
              ! "right" side of partition
              !
              rho_star(i,j,k) = rho_b0 * ut_r
              tau(i,j,k)   = tau_r
              h(i,j,k) = h_l
              w(i,j,k)     = rho_star(i,j,k) *ut_r
              st_x(i,j,k)  = rho_star(i,j,k) *h_l*ux_r
              st_y(i,j,k)  = rho_star(i,j,k) *h_l*uy_r
              st_z(i,j,k)  = rho_star(i,j,k) *h_l*uz_r
              mhd_st_x(i,j,k) = st_x(i,j,k) + ut_r*sbs_r*ux_r - sbt_r*sbx_r
              mhd_st_y(i,j,k) = st_y(i,j,k) + ut_r*sbs_r*uy_r - sbt_r*sby_r
              mhd_st_z(i,j,k) = st_z(i,j,k) + ut_r*sbs_r*uz_r - sbt_r*sbz_r
              vx(i,j,k)   = ux_r/ut_r
              vy(i,j,k)   = uy_r/ut_r
              vz(i,j,k)   = uz_r/ut_r
	      u0(i,j,k)   = ut_r
              Bx(i,j,k)   = Bx_r
              By(i,j,k)   = By_r
              Bz(i,j,k)   = Bz_r
              Bxtilde(i,j,k) = Bx_r
              Bytilde(i,j,k) = By_r
              Bztilde(i,j,k) = Bz_r
              sbt(i,j,k) = sbt_r
              sbx(i,j,k) = sbx_r
              sby(i,j,k) = sby_r
              sbz(i,j,k) = sbz_r
  	   else 
              !
              ! "in" side of partition
              !
	      thet = theta_l + rotations * 2.d0*pi * &
                  sin( (xs + wtrans*0.5d0)/wtrans*pi*0.5d0 )**2
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
              ut_i = gamma*(wf_ut + lambda*wf_u1)
              u1_i = gamma*(wf_u1 + lambda*wf_ut)
              u2_i = wf_u2
              u3_i = wf_u3
              sbt_i = gamma*(wf_bt + lambda*wf_b1)
              sb1_i = gamma*(wf_b1 + lambda*wf_bt)
              sb2_i = wf_b2
              sb3_i = wf_b3
              sbs_i = sb1_i*sb1_i + sb2_i*sb2_i + sb3_i*sb3_i - sbt_i*sbt_i
              B1_i = (ut_i*sb1_i - sbt_i*u1_i)*s4pi
              B2_i = (ut_i*sb2_i - sbt_i*u2_i)*s4pi
              B3_i = (ut_i*sb3_i - sbt_i*u3_i)*s4pi

	      call convert_v123_to_vxyz(u1_i,u2_i,u3_i,ux_i,uy_i,uz_i,nx,ny,nz)
	      call convert_v123_to_vxyz(sb1_i,sb2_i,sb3_i,sbx_i,sby_i,sbz_i,nx,ny,nz)
	      call convert_v123_to_vxyz(B1_i,B2_i,B3_i,Bx_i,By_i,Bz_i,nx,ny,nz)

              rho_star(i,j,k) = rho_b0 * ut_i
              tau(i,j,k) = rho_b0*ut_i*(h_l*ut_i - 1.d0) - P0 + sbs_i*(ut_i*ut_i-0.5d0) - sbt_i**2
              h(i,j,k) = h_l
              w(i,j,k)     = rho_star(i,j,k) *ut_i
              st_x(i,j,k)  = rho_star(i,j,k) *h_l*ux_i
              st_y(i,j,k)  = rho_star(i,j,k) *h_l*uy_i
              st_z(i,j,k)  = rho_star(i,j,k) *h_l*uz_i
              mhd_st_x(i,j,k) = st_x(i,j,k) + ut_i*sbs_i*ux_i - sbt_i*sbx_i
              mhd_st_y(i,j,k) = st_y(i,j,k) + ut_i*sbs_i*uy_i - sbt_i*sby_i
              mhd_st_z(i,j,k) = st_z(i,j,k) + ut_i*sbs_i*uz_i - sbt_i*sbz_i
              vx(i,j,k)   = ux_i/ut_i
              vy(i,j,k)   = uy_i/ut_i
              vz(i,j,k)   = uz_i/ut_i
              u0(i,j,k)   = ut_i
              Bx(i,j,k)   = Bx_i
              By(i,j,k)   = By_i
              Bz(i,j,k)   = Bz_i
              Bxtilde(i,j,k) = Bx_i
              Bytilde(i,j,k) = By_i
              Bztilde(i,j,k) = Bz_i
              sbt(i,j,k) = sbt_i
              sbx(i,j,k) = sbx_i
              sby(i,j,k) = sby_i
              sbz(i,j,k) = sbz_i
           end if

           ! Setup vector potential [use the expression of problem 5.33 in
	   ! D.J. Griffth's Introduction to Electrodynamics]
           ! Note that Ai are staggered
           A3 = 0.d0
	   xs = X(i,j,k)*nx + (Y(i,j,k)+0.5d0*dY)*ny + (Z(i,j,k)+0.5d0*dZ)*nz
           ys = X(i,j,k)*yx + (Y(i,j,k)+0.5d0*dY)*yy + (Z(i,j,k)+0.5d0*dZ)*yz
           zs = X(i,j,k)*zx + (Y(i,j,k)+0.5d0*dY)*zy + (Z(i,j,k)+0.5d0*dZ)*zz
     	   call nonlinear_alfven_wave_analytic(xs,ut_l,u1_l,u2_l,u3_l,&
                  B1_l,B2_l,B3_l, ut_r,u1_r,u2_r,u3_r,B1_r,B2_r,B3_r, &
                  ut_i,u1_i,u2_i,u3_i, B1_i,B2_i,B3_i,wtrans,rotations, & 
                  theta_l,chi_k, byz_rhs,a11_k,a12_k,a22_k,byc,bzc,  & 
                  wf_b1_l,wf_b2_l,wf_b3_l, &
                  wf_u1_l,wf_u2_l,wf_u3_l,lambda,gamma)
           call compute_int_B3dx(xs,wtrans,B3_r,int_B3dx_i,int_B3dx,dx1,num)
	   A1 = B2_i*zs
	   A2 = int_B3dx_i - B1_i*zs
	   Ax(i,j,k) = nx*A1 + yx*A2 + zx*A3

           xs = (X(i,j,k)+0.5d0*dX)*nx + Y(i,j,k)*ny + (Z(i,j,k)+0.5d0*dZ)*nz
           ys = (X(i,j,k)+0.5d0*dX)*yx + Y(i,j,k)*yy + (Z(i,j,k)+0.5d0*dZ)*yz
  	   zs = (X(i,j,k)+0.5d0*dX)*zx + Y(i,j,k)*zy + (Z(i,j,k)+0.5d0*dZ)*zz
           call nonlinear_alfven_wave_analytic(xs,ut_l,u1_l,u2_l,u3_l,&
                  B1_l,B2_l,B3_l, ut_r,u1_r,u2_r,u3_r,B1_r,B2_r,B3_r, &
                  ut_i,u1_i,u2_i,u3_i, B1_i,B2_i,B3_i,wtrans,rotations, &
                  theta_l,chi_k, byz_rhs,a11_k,a12_k,a22_k,byc,bzc,  &
                  wf_b1_l,wf_b2_l,wf_b3_l, &
                  wf_u1_l,wf_u2_l,wf_u3_l,lambda,gamma)
           call compute_int_B3dx(xs,wtrans,B3_r,int_B3dx_i,int_B3dx,dx1,num)
           A1 = B2_i*zs
           A2 = int_B3dx_i - B1_i*zs
	   Ay(i,j,k) = ny*A1 + yy*A2 + zy*A3

           xs = (X(i,j,k)+0.5d0*dX)*nx + (Y(i,j,k)+0.5d0*dY)*ny + Z(i,j,k)*nz
           ys = (X(i,j,k)+0.5d0*dX)*yx + (Y(i,j,k)+0.5d0*dY)*yy + Z(i,j,k)*yz
           zs = (X(i,j,k)+0.5d0*dX)*zx + (Y(i,j,k)+0.5d0*dY)*zy + Z(i,j,k)*zz
           call nonlinear_alfven_wave_analytic(xs,ut_l,u1_l,u2_l,u3_l,&
                  B1_l,B2_l,B3_l, ut_r,u1_r,u2_r,u3_r,B1_r,B2_r,B3_r, &
                  ut_i,u1_i,u2_i,u3_i, B1_i,B2_i,B3_i,wtrans,rotations, &
                  theta_l,chi_k, byz_rhs,a11_k,a12_k,a22_k,byc,bzc,  &
                  wf_b1_l,wf_b2_l,wf_b3_l, &
                  wf_u1_l,wf_u2_l,wf_u3_l,lambda,gamma)
           call compute_int_B3dx(xs,wtrans,B3_r,int_B3dx_i,int_B3dx,dx1,num)
           A1 = B2_i*zs
           A2 = int_B3dx_i - B1_i*zs
	   Az(i,j,k) = nz*A1 + yz*A3 + zz*A3
        end do
     end do
  end do

  neos = 1
  write(*,*) "Assuming polytropic EOS..."
  write(*,*) "Polytropic index = ",npolywave
  write(*,*) "Polytropic constant K = ",K_poly
  rho_tab(1) = 0.03d0
  P_tab(1) = K_poly * rho_tab(1)**(1.0+1.0/npolywave)
  eps_tab(1) = npolywave*P_tab(1)/rho_tab(1)
  gamma_tab(1) = 1.d0 + 1.d0/npolywave
  k_tab(1) = K_poly
  gamma_tab(2) = 1.d0 + 1.d0/npolywave 
  k_tab(2) = K_poly

  gamma_th = 1.d0 + 1.d0/npolywave 

  !======================================
  ! Set metric data to Minkowski
  !======================================
  phi = ZERO
  phix = ZERO
  phiy = ZERO
  phiz = ZERO
  psi = exp(phi)
  lapm1 = ZERO
  gxx = ONE
  gyy = ONE
  gzz = ONE
  gxy = ZERO
  gxz = ZERO
  gyz = ZERO
  gupxx = ONE
  gupyy = ONE
  gupzz = ONE
  gupxy = ZERO
  gupxz = ZERO
  gupyz = ZERO
  Gammax = ZERO
  Gammay = ZERO
  Gammaz = ZERO 
  !======================================
  ! Set K_ij to Zero!
  !======================================
  trK = ZERO
  Axx = ZERO
  Axy = ZERO
  Axz = ZERO
  Ayy = ZERO
  Ayz = ZERO
  Azz = ZERO
  !======================================
  ! Set everything else to Zero!
  !======================================
  shiftx = ZERO
  shifty = ZERO
  shiftz = ZERO

  !======================================
  ! Set initial EM fields to zero:
  !======================================

  Ex = By*vz - Bz*vy
  Ey = Bz*vx - Bx*vz
  Ez = Bx*vy - By*vx

  !FIXME: read note:
  ! NOTE: The following lines are not strictly correct, but this should NOT matter, since Cowling is turned on!
  Sxx = ZERO
  Sxy = ZERO
  Sxz = ZERO
  Syy = ZERO
  Syz = ZERO
  Szz = ZERO

  Sx = ZERO
  Sy = ZERO
  Sz = ZERO

  rho = ZERO

end subroutine nonlinear_alfven_wave_initialdata_local

subroutine convert_v123_to_vxyz(v1,v2,v3,vx,vy,vz,nx,ny,nz)
  implicit none 
  real*8 :: v1,v2,v3,vx,vy,vz,nx,ny,nz

  if (abs(nx) .lt. 1.d-300) then 
     ! In this case, the relationship between the coordinates in the shock frame 
     ! and the lab frame is \hat{xs}=\hat{n}, \hat{ys}=-\hat{x}, 
     ! \hat{zs} = \hat{xs} cross \hat{ys}
     vx = -v2
     vy = ny*v1 - nz*v3
     vz = nz*v1 + ny*v3
  else if (abs(ny) .lt. 1.d-300) then 
     ! In this case, the relationship between the coordinates in the shock frame
     ! and the lab frame is \hat{xs}=\hat{n}, \hat{ys}=\hat{y},
     ! \hat{zs} = \hat{xs} cross \hat{ys}
     vx = nx*v1 - nz*v3 
     vy = v2
     vz = nx*v3 + nz*v1 
  else 
     ! In this case, the relationship between the coordinates in the shock frame
     ! and the lab frame is \hat{xs}=\hat{n}, \hat{zs}=\hat{z},
     ! \hat{ys} = \hat{z} cross \hat{n}
     vx = nx*v1 - ny*v2
     vy = nx*v2 + ny*v1
     vz = v3
  end if

end subroutine convert_v123_to_vxyz

subroutine compute_int_B3dx(x,wtrans,B3_r,int_B3dx_i,int_B3dx,dx,num)
  implicit none
  integer :: num
  real*8, dimension(num) :: int_B3dx
  real*8 :: x,int_B3dx_i,wtrans,dx,tmp,B3_r
  integer, parameter :: m = 6
  real*8, dimension(m) :: xa,ya
  integer :: ilo,i
! 
  if (x .le. -0.5d0*wtrans) then 
     int_B3dx_i = 0.d0
  else if (x .ge. 0.5d0*wtrans) then
     int_B3dx_i = int_B3dx(num) + B3_r * (x-0.5d0*wtrans)
  else
     ilo = int( (x+0.5d0*wtrans)/dx ) + 1 - m/2
     if (ilo .lt. 0) ilo = 0
     if (ilo .gt. num-m) ilo = num-m
     do i=1,m
	xa(i) = (ilo+i-1)*dx - 0.5d0*wtrans
	ya(i) = int_B3dx(ilo+i)
     end do
     call polint_nonlinear_alfven_wave(xa,ya,m,x,int_B3dx_i,tmp)
  end if
end subroutine compute_int_B3dx

!-----------------------------------------------------------------------------
!
!  Numerical Recipe interpolation subroutine
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE polint_nonlinear_alfven_wave(xa,ya,n,x,y,dy)
      IMPLICIT NONE
      integer                    :: n,NMAX,i,m,ns
      real*8                     :: dy,x,y
      real*8, dimension(n)       :: xa,ya
      PARAMETER (NMAX=10)
      real*8, dimension(NMAX)    :: c,d
      real*8                     :: den,dif,dift,ho,hp,w
!
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) then
            write(*,*) 'failure in polint'
            stop
          end if
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END SUBROUTINE polint_nonlinear_alfven_wave
