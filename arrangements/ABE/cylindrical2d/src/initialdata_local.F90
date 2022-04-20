!-------------------------------------------------------------------
! Here we set up all required variables for shock tests, including:
!  emfields, BSSN variables, primitives, etc.
!-------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine cylindrical2d_initialdata_local(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8 :: s4pi, Bx_l, u0l,ux_l,uy_l,uz_l,enth,vx_l,vy_l,rho_l,P_l
  real*8 :: sb0l,sbxl,sbyl,sbzl,sbs,dX,dY,dZ
  real*8 :: r_in,r_out,lambda,fact, pomega, decfact
  real*8, parameter :: pi = 3.1415926535897932d0
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM, i,j,k
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
!

  !!if(primitives_solver.eq.0 .or. enable_shocktest_primitive_mode.ne.1) then
  !!   write(*,*) "ERROR.  YOU MUST SET primitives_solver=1 or 2 AND enable_shocktest_primitive_mode=1 FOR THIS TEST"
  !!   stop
  !!end if

  if (disk_outer_radius .lt. disk_inner_radius) then 
     write(*,*) 'ERROR. disk_outer_radius must not be smaller than disk_inner_radius'
     stop
  end if

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)
 
  s4pi = 2.d0*sqrt(pi)
  Bx_l = Bxos4pi * s4pi
  r_in = disk_inner_radius
  r_out = disk_outer_radius
  lambda = 0.5d0/(r_out-r_in)
  if (r_out .gt. r_in) then 
     fact = 1.d0/(1.d0 - exp(-lambda*(r_out-r_in)) )
  else 
     fact = 1.d0
  end if

  if (abs(rot_omega*r_out) .ge. 1.d0) then 
     write(*,*) 'ERROR. The disk is rotating faster than the speed of light!'
     write(*,*) 'Try lower the angular velocity of rotation or decrease the disk radius.' 
     stop
  end if

  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)

	   pomega = sqrt(X(i,j,k)**2 + Y(i,j,k)**2)
           if (pomega .le. r_in) then 
	      rho_l = rho_in
	      P_l = P_in
              vx_l = -rot_omega * Y(i,j,k)
              vy_l =  rot_omega * X(i,j,k)
           else if (pomega .ge. r_out) then 
              rho_l = rho_out
              P_l = P_out
              vx_l = 0.d0
              vy_l = 0.d0
	   else
	      decfact = fact*(1.d0-exp(-lambda*(pomega-r_in)))
	      rho_l = rho_in + (rho_out-rho_in)*decfact
	      P_l = P_in + (P_out-P_in)*decfact
	      vx_l = -rot_omega * Y(i,j,k) * (1.d0-decfact)
	      vy_l =  rot_omega * X(i,j,k) * (1.d0-decfact)
           end if
           rho_b(i,j,k) = rho_l
	   P(i,j,k) = P_l
           vx(i,j,k) = vx_l
	   vy(i,j,k) = vy_l
	   vz(i,j,k) = 0.d0

	   u0l = 1.d0/sqrt(1.d0 - vx_l**2 - vy_l**2) 
	   ux_l = u0l*vx_l
	   uy_l = u0l*vy_l
	   uz_l = 0.d0
	   enth = 1.d0 + (npolyshock+1.d0)*P_l/rho_l
	   sb0l = ux_l*Bxos4pi
	   sbxl = (Bxos4pi + sb0l*ux_l)/u0l
	   sbyl = sb0l*vy_l
	   sbzl = 0.d0
	   sbs = sbxl**2 + sbyl**2 + sbzl**2 - sb0l**2

	   u0(i,j,k) = u0l
	   h(i,j,k) = enth
 	   rho_star(i,j,k) = rho_l * u0l
	   tau(i,j,k) = rho_l*u0l*(enth*u0l-1.d0)-P_l + sbs*(u0l*u0l-0.5d0) -sb0l**2
	   st_x(i,j,k)  = rho_star(i,j,k) * enth * ux_l
           st_y(i,j,k)  = rho_star(i,j,k) * enth * uy_l
           st_z(i,j,k)  = rho_star(i,j,k) * enth * uz_l
	   mhd_st_x(i,j,k) = st_x(i,j,k) + u0l*sbs*ux_l - sb0l*sbxl
           mhd_st_y(i,j,k) = st_y(i,j,k) + u0l*sbs*uy_l - sb0l*sbyl
           mhd_st_z(i,j,k) = st_z(i,j,k) + u0l*sbs*uz_l - sb0l*sbzl
	   Bx(i,j,k) = Bx_l
	   By(i,j,k) = 0.d0
	   Bz(i,j,k) = 0.d0
	   sbt(i,j,k) = sb0l
	   sbx(i,j,k) = sbxl
           sby(i,j,k) = sbyl
           sbz(i,j,k) = sbzl
           Bxtilde(i,j,k) = Bx_l
           Bytilde(i,j,k) = 0.d0
           Bztilde(i,j,k) = 0.d0

           ! Setup vector potential [use the expression of problem 5.33 in
	   ! D.J. Griffth's Introduction to Electrodynamics]
           ! Note that Ai are staggered
	   Ax(i,j,k) = 0.d0
	   Ay(i,j,k) = 0.d0
	   Az(i,j,k) = Bx_l*(Y(i,j,k) + 0.5d0*dY)
        end do
     end do
  end do

  neos = 1
  write(*,*) "Assuming polytropic EOS..."
  write(*,*) "Polytropic index = ",npolyshock
  K_poly = P_out/rho_out**(1.0+1.0/npolyshock)
  write(*,*) "Polytropic constant K = ",K_poly
  rho_tab(1) = rho_out*1.078172662d0
  P_tab(1) = K_poly*rho_tab(1)**(1.0+1.0/npolyshock)
  eps_tab(1) = npolyshock*P_tab(1)/rho_tab(1)
  gamma_tab(1) = 1.d0 + 1.d0/npolyshock
  k_tab(1) = K_poly
  gamma_tab(2) = 1.d0 + 1.d0/npolyshock
  k_tab(2) = K_poly

  gamma_th = 1.d0 + 1.d0/npolyshock

  !======================================
  ! Set metric data to Minkowski
  !======================================
  phi = 0.d0
  phix = 0.d0
  phiy = 0.d0
  phiz = 0.d0
  psi = exp(phi)
  lapm1 = 0.d0
  gxx = 1.d0
  gyy = 1.d0
  gzz = 1.d0
  gxy = 0.d0
  gxz = 0.d0
  gyz = 0.d0
  gupxx = 1.d0
  gupyy = 1.d0
  gupzz = 1.d0
  gupxy = 0.d0
  gupxz = 0.d0
  gupyz = 0.d0
  Gammax = 0.d0
  Gammay = 0.d0
  Gammaz = 0.d0 
  !======================================
  ! Set K_ij to Zero!
  !======================================
  trK = 0.d0
  Axx = 0.d0
  Axy = 0.d0
  Axz = 0.d0
  Ayy = 0.d0
  Ayz = 0.d0
  Azz = 0.d0
  !======================================
  ! Set everything else to Zero!
  !======================================
  shiftx = 0.d0
  shifty = 0.d0
  shiftz = 0.d0

  !======================================
  ! Set initial EM fields to zero:
  !======================================

  Ex = By*vz - Bz*vy
  Ey = Bz*vx - Bx*vz
  Ez = Bx*vy - By*vx

  !FIXME: read note:
  ! NOTE: The following lines are not strictly correct, but this should NOT matter, since Cowling is turned on!
  Sxx = 0.d0
  Sxy = 0.d0
  Sxz = 0.d0
  Syy = 0.d0
  Syz = 0.d0
  Szz = 0.d0

  Sx = 0.d0
  Sy = 0.d0
  Sz = 0.d0

  rho = 0.d0

end subroutine cylindrical2d_initialdata_local
