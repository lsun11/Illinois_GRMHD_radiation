!-------------------------------------------------------------------
! Here we set up all required variables for shock tests, including:
!  emfields, BSSN variables, primitives, radiations, etc.
!-------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine mhd_shock_initialdata_local(CCTK_ARGUMENTS)
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
  real*8                :: norm, w_l,w_r,h_l,h_r,ut_l,ut_r,hm1_l,hm1_r,tau_l,tau_r
  real*8		:: tau_rad_l, tau_rad_r
  real*8                :: ux_l,ux_r,uy_l,uy_r,uz_l,uz_r, F_rad0l, P_radl
  real*8                :: xs,ys,yx,yy,yz,zx,zy,zz
  real*8                :: sbx_l,sbx_r,sby_l,sby_r,sbz_l,sbz_r
  real*8                :: Bx_l,Bx_r,By_l,By_r,Bz_l,Bz_r
  real*8                :: B1_l,B2_l,B3_l,B1_r,B2_r,B3_r
  real*8 		:: A1,A2,A3
  real*8, parameter :: pi = 3.1415926535897932d0
  real*8 :: s4pi,sb0_l,sb0_r,sb1_l,sb1_r,sb2_l,sb2_r,sb3_r,sb3_l,sbs_l,sbs_r
  real*8		:: w_boost, P_save, P_savep1, rho_save, rho_savep1, vx_save, vx_savep1, u0_save, u0_savep1, X_p
  real*8                :: zeta, chi_shock, Fksq, C_A, C_B, P_rad0xl, P_rad0yl, P_rad0zl, P_rad00l

  real*8, dimension(1)  :: xsave, tausave, rhosave,Psave,uxsave,Eradsave,Fradsave
  real*8, dimension(1) :: xread,rho0read,Pread,uxread,Eread,Fxread,tauread

  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  integer :: ONE,ZERO
  parameter(ONE = 1.D0, ZERO = 0.D0)

  

  if(primitives_solver.eq.0 .or. enable_shocktest_primitive_mode.ne.1) then
     write(*,*) "ERROR.  YOU MUST SET primitives_solver=1 or 2 AND enable_shocktest_primitive_mode=1 FOR SHOCK TESTS"
     stop
  end if
 
  ! sqrt(4 pi) 
  s4pi = 2.d0*sqrt(pi)

  ext = cctk_lsh

  write(*,*) ext(1), ext(2), ext(3)



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

  !!write(*,*) "Okay... Looks like you set up the grid correctly EXTENTS:",cctk_lsh
  !!write(*,*) "Okay... Looks like you set up the grid correctly GZS:",cctk_nghostzones

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ! u^0 on the left (w_l) and right (w_r) side of the shock
  w_l = sqrt(1.0 + u1_l**2 + u2_l**2 + u3_l**2)
  w_r = sqrt(1.0 + u1_r**2 + u2_r**2 + u3_r**2)

  B1_l = B1os4pi_l * s4pi
  B2_l = B2os4pi_l * s4pi
  B3_l = B3os4pi_l * s4pi
  B1_r = B1os4pi_r * s4pi
  B2_r = B2os4pi_r * s4pi
  B3_r = B3os4pi_r * s4pi

  ! Compute b^mu
  sb0_l = u1_l*B1os4pi_l + u2_l*B2os4pi_l + u3_l*B3os4pi_l
  sb0_r = u1_r*B1os4pi_r + u2_r*B2os4pi_r + u3_r*B3os4pi_r
  sb1_l = (B1os4pi_l + sb0_l*u1_l)/w_l
  sb1_r = (B1os4pi_r + sb0_r*u1_r)/w_r
  sb2_l = (B2os4pi_l + sb0_l*u2_l)/w_l
  sb2_r = (B2os4pi_r + sb0_r*u2_r)/w_r
  sb3_l = (B3os4pi_l + sb0_l*u3_l)/w_l
  sb3_r = (B3os4pi_r + sb0_r*u3_r)/w_r
  sbs_l = sb1_l**2 + sb2_l**2 + sb3_l**2 - sb0_l**2
  sbs_r = sb1_r**2 + sb2_r**2 + sb3_r**2 - sb0_r**2

  hm1_l = ( npolyshock + 1.d0 ) * P_l / rho_l
  hm1_r = ( npolyshock + 1.d0 ) * P_r / rho_r  
  h_l = 1.d0 + hm1_l
  h_r = 1.d0 + hm1_r

  tau_l = rho_l*w_l*(h_l*w_l - 1.d0) - P_l + sbs_l*(w_l*w_l-0.5d0) - sb0_l**2
  tau_r = rho_r*w_r*(h_r*w_r - 1.d0) - P_r + sbs_r*(w_r*w_r-0.5d0) - sb0_r**2

  ! just to make sure, normalize n:
  !
  norm = sqrt(nx*nx + ny*ny + nz*nz)
  nx = nx/norm
  ny = ny/norm
  nz = nz/norm

  if (nx==0.d0) then 
     ! In this case, the relationship between the coordinates in the shock frame 
     ! and the lab frame is \hat{xs}=\hat{n}, \hat{ys}=-\hat{x}, 
     ! \hat{zs} = \hat{xs} cross \hat{ys}
     yx = -1.d0
     yy = 0.d0
     yz = 0.d0
     zx = 0.d0
     zy = -nz
     zz = ny
     ux_l = -u2_l
     ux_r = -u2_r
     uy_l = ny*u1_l - nz*u3_l
     uy_r = ny*u1_r - nz*u3_r
     uz_l = nz*u1_l + ny*u3_l
     uz_r = nz*u1_r + ny*u3_r
     sbx_l = -sb2_l
     sbx_r = -sb2_r
     sby_l = ny*sb1_l - nz*sb3_l
     sby_r = ny*sb1_r - nz*sb3_r
     sbz_l = nz*sb1_l + ny*sb3_l
     sbz_r = nz*sb1_r + ny*sb3_r
     Bx_l = -B2os4pi_l*s4pi
     Bx_r = -B2os4pi_r*s4pi
     By_l = (ny*B1os4pi_l - nz*B3os4pi_l)*s4pi
     By_r = (ny*B1os4pi_r - nz*B3os4pi_r)*s4pi
     Bz_l = (nz*B1os4pi_l + ny*B3os4pi_l)*s4pi
     Bz_r = (nz*B1os4pi_r + ny*B3os4pi_r)*s4pi
  else if (ny==0.d0) then 
     ! In this case, the relationship between the coordinates in the shock frame
     ! and the lab frame is \hat{xs}=\hat{n}, \hat{ys}=\hat{y},
     ! \hat{zs} = \hat{xs} cross \hat{ys}
     yx = 0.d0
     yy = 1.d0
     yz = 0.d0
     zx = -nz
     zy = 0.d0
     zz = nx
     ux_l = nx*u1_l - nz*u3_l 
     ux_r = nx*u1_r - nz*u3_r
     uy_l = u2_l 
     uy_r = u2_r
     uz_l = nx*u3_l + nz*u1_l 
     uz_r = nx*u3_r + nz*u1_r
     sbx_l = nx*sb1_l - nz*sb3_l
     sbx_r = nx*sb1_r - nz*sb3_r
     sby_l = sb2_l
     sby_r = sb2_r
     sbz_l = nx*sb3_l + nz*sb1_l
     sbz_r = nx*sb3_r + nz*sb1_r
     Bx_l = (nx*B1os4pi_l - nz*B3os4pi_l)*s4pi
     Bx_r = (nx*B1os4pi_r - nz*B3os4pi_r)*s4pi
     By_l = B2os4pi_l*s4pi
     By_r = B2os4pi_r*s4pi
     Bz_l = (nx*B3os4pi_l + nz*B1os4pi_l)*s4pi
     Bz_r = (nx*B3os4pi_r + nz*B1os4pi_r)*s4pi
  else if (nz==0.d0) then 
     ! In this case, the relationship between the coordinates in the shock frame
     ! and the lab frame is \hat{xs}=\hat{n}, \hat{zs}=\hat{z},
     ! \hat{ys} = \hat{z} cross \hat{n}
     yx = -ny
     yy = nx
     yz = 0.d0
     zx = 1.d0
     zy = 0.d0
     zz = 0.d0
     ux_l = nx*u1_l - ny*u2_l 
     ux_r = nx*u1_r - ny*u2_r
     uy_l = nx*u2_l + ny*u1_l
     uy_r = nx*u2_r + ny*u1_r
     uz_l = u3_l
     uz_r = u3_r
     sbx_l = nx*sb1_l - ny*sb2_l
     sbx_r = nx*sb1_r - ny*sb2_r
     sby_l = nx*sb2_l + ny*sb1_l
     sby_r = nx*sb2_r + ny*sb1_r
     sbz_l = sb3_l
     sbz_r = sb3_r
     Bx_l = (nx*B1os4pi_l - ny*B2os4pi_l)*s4pi
     Bx_r = (nx*B1os4pi_r - ny*B2os4pi_r)*s4pi
     By_l = (nx*B2os4pi_l + ny*B1os4pi_l)*s4pi
     By_r = (nx*B2os4pi_r + ny*B1os4pi_r)*s4pi
     Bz_l = B3os4pi_l*s4pi
     Bz_r = B3os4pi_r*s4pi
  else 
     write(*,*) 'At least one of the nx, ny or ny must be set to 0'
     stop
  end if
  

if (read_radiation_shock_data .eq. 1) then
 write(*,*) 'DIRECTLY READ SHOCK DATA FROM INPUT FILE!!'
 write(*,*) 'rad_const =', rad_const
! here we call radiationshock_readfiles to read in rho0, P, E_rad, F_rad, vx, and tau.
! so at this point the 6 GFs should already be filled.
! we just need to fill in others.


        do k=1,ext(3)
          do j=1,ext(2)
            do i=1,ext(1)


	        if (rho_b(i,j,k) == 0.d0) then
		    write(*,*) " At the begining of initialdata_local.F90, zero rho is found!!!!", i,j,k
		end if


	        ! here we need to distinguish wether the shock is boosted!!!!

		w(i,j,k) = 1.d0/sqrt(1.d0 - vx(i,j,k)**2)
		if (boost_shock .eq. 1) then
		! u0(i,j,k) = w(i,j,k)/sqrt(1.d0 - u_boost**2)*(1.d0 - u_boost*vx(i,j,k))
		! The following expression is correct because vx(i,j,k) is already boosted in readinputfile_read.C                                                                
		u0(i,j,k) = 1.d0/sqrt(1.d0 - vx(i,j,k)**2)
 		! F_radx(i,j,k) here is already boosted in readinputfile_read.C so only a factor of -u_boost is needed.
		!F_rad0(i,j,k) = u_boost* F_radx(i,j,k)
		F_rad0(i,j,k) = F_radx(i,j,k)*vx(i,j,k)
                else		
                u0(i,j,k) = w(i,j,k)
                ! F_rad0 = - F_radx*u_x/u_0 = - F_radx*v_x = F_radx*vx (u_0 = - u0 / v_x = -vx) 
		F_rad0(i,j,k) = F_radx(i,j,k)*vx(i,j,k)
		end if 


if (i==4.and.j==4.and.k==4) then
    write(*,*) " Inside shock_mhd initial_data_local  u0(i,j,k) is ", u0(i,j,k)
    write(*,*) " w(i,j,k), u_boost, vx(i,j,k) are",  w(i,j,k), u_boost, vx(i,j,k)
    write(*,*) "E_rad and F_radx are", E_rad(i,j,k), F_radx(i,j,k)
end if


		!vx(i,j,k) = ux(i,j,k)/w(i,j,k)
                vy(i,j,k) = uy_l/u0(i,j,k)
                vz(i,j,k) = uz_l/u0(i,j,k)

                F_rady(i,j,k) = ZERO
                F_radz(i,j,k) = ZERO
                P_radl = E_rad(i,j,k)/3.0

                rho_star(i,j,k) = rho_b(i,j,k) * u0(i,j,k)

		! Lunan : add this if statement to aviod dividing by zero.

		if (rho_b(i,j,k) == 0.d0 .and. P(i,j,k) == 0.d0) then
		h(i,j,k) = 1.d0 + ( npolyshock + 1.d0 )
		else
                h(i,j,k) = 1.d0 + ( npolyshock + 1.d0 ) * P(i,j,k) / rho_b(i,j,k)
		end if



		tau(i,j,k) = rho_star(i,j,k)*(h(i,j,k)*u0(i,j,k)-1.d0)-P(i,j,k)+sbs_l*(u0(i,j,k)*u0(i,j,k)-0.5d0)-sb0_l**2


                st_x(i,j,k)  = rho_star(i,j,k)*h(i,j,k)*vx(i,j,k)/sqrt(1.d0-vx(i,j,k)**2)
                !st_x(i,j,k)  = rho_star(i,j,k)*h(i,j,k)*vx(i,j,k)/(1.d0-vx(i,j,k)**2)
                st_y(i,j,k)  = rho_star(i,j,k)*h(i,j,k)*uy_l
                st_z(i,j,k)  = rho_star(i,j,k)*h(i,j,k)*uz_l
                
                mhd_st_x(i,j,k) = st_x(i,j,k) + u0(i,j,k)*sbs_l*vx(i,j,k)/sqrt(1-vx(i,j,k)**2) - sb0_l*sbx_l
                mhd_st_y(i,j,k) = st_y(i,j,k) + u0(i,j,k)*sbs_l*uy_l - sb0_l*sby_l
                mhd_st_z(i,j,k) = st_z(i,j,k) + u0(i,j,k)*sbs_l*uz_l - sb0_l*sbz_l
                

                Bx(i,j,k)   = Bx_l
                By(i,j,k)   = By_l
                Bz(i,j,k)   = Bz_l
                Bxtilde(i,j,k) = Bx_l
                Bytilde(i,j,k) = By_l
                Bztilde(i,j,k) = Bz_l
                sbt(i,j,k) = sb0_l
                sbx(i,j,k) = sbx_l
                sby(i,j,k) = sby_l
                sbz(i,j,k) = sbz_l

		!T_fluid(i,j,k) = (E_rad(i,j,k)/rad_const)**(0.25)
		T_fluid(i,j,k) = P(i,j,k)/rho(i,j,k)
                if (rad_evolve_enable .eq. 1) then
                   if(rad_closure_scheme .eq. 0) then
                      tau_rad(i,j,k) = (E_rad(i,j,k) + P_radl)*u0(i,j,k)*u0(i,j,k) + 2.0*u0(i,j,k)*F_rad0(i,j,k)-P_radl
                      ! In flat space, F_rad_alpha = g_{alpha beta} F_rad^alpha
                      ! Notice that: F_rad_0 = -F_rad0, F_rad_i = F_radi, u_0 = -u0, v_i = -vi, u_i = ui
                      S_rad_x(i,j,k) = (E_rad(i,j,k) + P_radl)*u0(i,j,k)*u0(i,j,k)*vx(i,j,k) + F_rad0(i,j,k)*u0(i,j,k)*vx(i,j,k) + F_radx(i,j,k)*u0(i,j,k)
                      S_rad_y(i,j,k) = (E_rad(i,j,k) + P_radl)*u0(i,j,k)*uy_l + F_rad0(i,j,k)*uy_l + F_rady(i,j,k)*u0(i,j,k)
                      S_rad_z(i,j,k) = (E_rad(i,j,k) + P_radl)*u0(i,j,k)*uz_l + F_rad0(i,j,k)*uz_l + F_radz(i,j,k)*u0(i,j,k)
                   else
                      zeta = sqrt((-F_rad0(i,j,k)**2 + F_radx(i,j,k)**2 + F_rady(i,j,k)**2 + F_radz(i,j,k)**2)/E_rad(i,j,k)**2)
                      chi_shock  = 1.0/3.0 + zeta**2*(6.0-2.0*zeta+6.0*zeta**2)/15.0
                      chi_rad(i,j,k)=chi_shock
                      Fksq = F_radx(i,j,k)**2 + F_rady(i,j,k)**2 + F_radz(i,j,k)**2
                      C_A = (3.0*chi_shock -1.0)/2.0
                      C_B = 1.5*(1.0-chi_shock)
                      P_radxx(i,j,k) = F_radx(i,j,k)**2*E_rad(i,j,k)/Fksq*C_A + E_rad(i,j,k)*(1.0+(u0(i,j,k)*vx(i,j,k))**2)/3.0*C_B
                      P_radyy(i,j,k) = F_rady(i,j,k)**2*E_rad(i,j,k)/Fksq*C_A + E_rad(i,j,k)*(1.0+uy_l**2)/3.0*C_B
                      P_radzz(i,j,k) = F_radz(i,j,k)**2*E_rad(i,j,k)/Fksq*C_A + E_rad(i,j,k)*(1.0+uz_l**2)/3.0*C_B
                      P_radxy(i,j,k) = F_radx(i,j,k)*F_rady(i,j,k)/Fksq*C_A + E_rad(i,j,k)* u0(i,j,k)*vx(i,j,k)*uy_l/3.0*C_B
                      P_radxz(i,j,k) = F_radx(i,j,k)*F_radz(i,j,k)/Fksq*C_A + E_rad(i,j,k)* u0(i,j,k)*vz(i,j,k)*uz_l/3.0*C_B
                      P_radyz(i,j,k) = F_rady(i,j,k)*F_radz(i,j,k)/Fksq*C_A + E_rad(i,j,k)* uy_l*uz_l/3.0*C_B
                    ! In flat space, P_rad^i_j = P_rad^{ij} 
                      P_rad0xl = (P_radxx(i,j,k)*u0(i,j,k)*vx(i,j,k) + P_radxy(i,j,k) * uy_l + P_radxz(i,j,k) * uz_l)/u0(i,j,k)
                      P_rad0yl = (P_radxy(i,j,k)*u0(i,j,k)*vx(i,j,k) + P_radyy(i,j,k) * uy_l + P_radyz(i,j,k) * uz_l)/u0(i,j,k)
                      P_rad0yl = (P_radxz(i,j,k)*u0(i,j,k)*vx(i,j,k) + P_radyz(i,j,k) * uy_l + P_radzz(i,j,k) * uz_l)/u0(i,j,k)
                      P_rad00l = (P_rad0xl*u0(i,j,k)*vx(i,j,k) + P_rad0yl * uy_l + P_rad0zl * uz_l)/u0(i,j,k)

                      tau_rad(i,j,k) = E_rad(i,j,k)*u0(i,j,k)*u0(i,j,k) + 2.0*u0(i,j,k)*F_rad0(i,j,k) + P_rad00l
                      S_rad_x(i,j,k) = E_rad(i,j,k)*u0(i,j,k)*u0(i,j,k)*vx(i,j,k) + F_rad0(i,j,k)*u0(i,j,k)*vx(i,j,k) &
                                       + F_radx(i,j,k)*u0(i,j,k) + P_rad0xl
                      S_rad_y(i,j,k) = E_rad(i,j,k)*u0(i,j,k)*uy_l + F_rad0(i,j,k)*uy_l + F_rady(i,j,k)*u0(i,j,k) + P_rad0yl
                      S_rad_z(i,j,k) = E_rad(i,j,k)*u0(i,j,k)*uz_l + F_rad0(i,j,k)*uz_l + F_radz(i,j,k)*u0(i,j,k) + P_rad0zl
                   end if
                else 
		   tau_rad(i,j,k) = 0.d0   
		   S_rad_x(i,j,k) = 0.d0
		   S_rad_y(i,j,k) = 0.d0
		   S_rad_z(i,j,k) = 0.d0
		end if




           ! D.J. Griffth's Introduction to Electrodynamics]
           ! Note that Ai are staggered
                A1 = 0.d0
                xs = X(i,j,k)*nx + (Y(i,j,k)+0.5d0*dY)*ny + (Z(i,j,k)+0.5d0*dZ)*nz
                ys = X(i,j,k)*yx + (Y(i,j,k)+0.5d0*dY)*yy + (Z(i,j,k)+0.5d0*dZ)*yz
                if (xs .le. 0.d0) then
                   A2 = B3_l*xs
	                 A3 = B1_l*ys - B2_l*xs
                else
                   A2 = B3_r*xs
	                 A3 = B1_r*ys - B2_r*xs
                end if
           	Ax(i,j,k) = nx*A1 + yx*A2 + zx*A3
           	xs = (X(i,j,k)+0.5d0*dX)*nx + Y(i,j,k)*ny + (Z(i,j,k)+0.5d0*dZ)*nz
           	ys = (X(i,j,k)+0.5d0*dX)*yx + Y(i,j,k)*yy + (Z(i,j,k)+0.5d0*dZ)*yz
           	if (xs .le. 0.d0) then
              	   A2 = B3_l*xs
              	   A3 = B1_l*ys - B2_l*xs
           	else
		   A2 = B3_r*xs
              	   A3 = B1_r*ys - B2_r*xs
           	end if
           	   Ay(i,j,k) = ny*A1 + yy*A2 + zy*A3
           	   xs = (X(i,j,k)+0.5d0*dX)*nx + (Y(i,j,k)+0.5d0*dY)*ny + Z(i,j,k)*nz
          	   ys = (X(i,j,k)+0.5d0*dX)*yx + (Y(i,j,k)+0.5d0*dY)*yy + Z(i,j,k)*yz
		   if (xs .le. 0.d0) then
                   A2 = B3_l*xs
	                 A3 = B1_l*ys - B2_l*xs
                else
                   A2 = B3_r*xs
	                 A3 = B1_r*ys - B2_r*xs
		end if
                Az(i,j,k) = nz*A1 + yz*A3 + zz*A3

	end do
     end do
  end do


else

 write(*,*) 'CALCULATE RAD SHOCK DATA USING PARFILE!'
	 
 if (rad_evolve_enable .eq. 1) then
       write(*,*) 'integrate_flag is', integrate_flag
	        
  call  initialdata_radiation (ext, x, npolyshock, &
	rad_opacity_abs, rad_opacity_sct, &
	rho_l, P_l, u1_l, Erad_l, Fradx_l, Fradx_r, & 
	cctkGH, cctk_lsh, integrate_flag, xsave, tausave, &
	rhosave, Psave, uxsave, Eradsave, Fradsave)

  	do k=1,ext(3)
     	  do j=1,ext(2)
            do i=1,ext(1)

              xs = X(i,j,k)*nx + Y(i,j,k)*ny + Z(i,j,k)*nz
              if (xs .le. 0.D0) then
             !
             ! "left" side of partition

	   	w(i,j,k) = sqrt(1.d0 + uxsave(i)**2)
		rho_b(i,j,k) = rhosave(i)
	   	P(i,j,k) = Psave(i)
	   	vx(i,j,k) = uxsave(i)/w(i,j,k)
	   	vy(i,j,k) = uy_l/w(i,j,k)
	   	vz(i,j,k) = uz_l/w(i,j,k)

	   	E_rad(i,j,k) = Eradsave(i)
	   	F_radx(i,j,k) = Fradsave(i)
		
	        F_rad0(i,j,k) = ZERO
                F_rady(i,j,k) = ZERO
                F_radz(i,j,k) = ZERO
                P_radl = E_rad(i,j,k)/3.0

	   	tau (i,j,k) = tausave(i)
	   	rho_star(i,j,k) = rho_b(i,j,k) * w(i,j,k)
	   	h(i,j,k) = 1.d0 + ( npolyshock + 1.d0 ) * P(i,j,k) / rho_b(i,j,k)

	   	st_x(i,j,k)  = rho_star(i,j,k)*h(i,j,k)*uxsave(i)
	   	st_y(i,j,k)  = rho_star(i,j,k)*h(i,j,k)*uy_l
          	st_z(i,j,k)  = rho_star(i,j,k)*h(i,j,k)*uz_l
           	mhd_st_x(i,j,k) = st_x(i,j,k) + w(i,j,k)*sbs_l*uxsave(i) - sb0_l*sbx_l
           	mhd_st_y(i,j,k) = st_y(i,j,k) + w(i,j,k)*sbs_l*uy_l - sb0_l*sby_l
           	mhd_st_z(i,j,k) = st_z(i,j,k) + w(i,j,k)*sbs_l*uz_l - sb0_l*sbz_l
	   	u0(i,j,k) = w (i,j,k)
	   
		Bx(i,j,k)   = Bx_l
        	By(i,j,k)   = By_l
           	Bz(i,j,k)   = Bz_l
           	Bxtilde(i,j,k) = Bx_l
           	Bytilde(i,j,k) = By_l
           	Bztilde(i,j,k) = Bz_l
           	sbt(i,j,k) = sb0_l
           	sbx(i,j,k) = sbx_l
           	sby(i,j,k) = sby_l
           	sbz(i,j,k) = sbz_l

		!T_fluid(i,j,k) = (E_rad(i,j,k)/rad_const)**(0.25)
                T_fluid(i,j,k) = P(i,j,k)/rho(i,j,k)
                tau_rad(i,j,k) = w(i,j,k)*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*u0(i,j,k) + 2.0*u0(i,j,k)*F_rad0(i,j,k)-P_radl)
                S_rad_x(i,j,k) = w(i,j,k)*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*uxsave(i) + F_rad0(i,j,k)*uxsave(i) + F_radx(i,j,k)*u0(i,j,k))
                S_rad_y(i,j,k) = w(i,j,k)*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*uy_l + F_rad0(i,j,k)*uy_l + F_rady(i,j,k)*u0(i,j,k))
                S_rad_z(i,j,k) = w(i,j,k)*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*uz_l + F_rad0(i,j,k)*uz_l + F_radz(i,j,k)*u0(i,j,k))

	   else
              !
              ! "right" side of partition
              !
              P(i,j,k)    = P_r
              rho_b(i,j,k) = rho_r
              rho_star(i,j,k) = rho_r * w_r
              tau(i,j,k)   = tau_r
              h(i,j,k) = h_r
              w(i,j,k)    = rho_star(i,j,k) *w_r
              st_x(i,j,k)  = rho_star(i,j,k) *h_r*ux_r
              st_y(i,j,k)  = rho_star(i,j,k) *h_r*uy_r
              st_z(i,j,k)  = rho_star(i,j,k) *h_r*uz_r
              mhd_st_x(i,j,k) = st_x(i,j,k) + w_r*sbs_r*ux_r - sb0_r*sbx_r
              mhd_st_y(i,j,k) = st_y(i,j,k) + w_r*sbs_r*uy_r - sb0_r*sby_r
              mhd_st_z(i,j,k) = st_z(i,j,k) + w_r*sbs_r*uz_r - sb0_r*sbz_r
              vx(i,j,k)   = ux_r/w_r
              vy(i,j,k)   = uy_r/w_r
              vz(i,j,k)   = uz_r/w_r
              u0(i,j,k)   = w_r
	      Bx(i,j,k)   = Bx_r
              By(i,j,k)   = By_r
              Bz(i,j,k)   = Bz_r
              Bxtilde(i,j,k) = Bx_r
              Bytilde(i,j,k) = By_r
              Bztilde(i,j,k) = Bz_r
              sbt(i,j,k) = sb0_r
              sbx(i,j,k) = sbx_r
              sby(i,j,k) = sby_r
              sbz(i,j,k) = sbz_r
              E_rad(i,j,k) = Erad_r
              F_rad0(i,j,k) = ZERO
              F_radx(i,j,k) = Fradx_r
              F_rady(i,j,k) = ZERO
              F_radz(i,j,k) = ZERO
              P_radl = Erad_r/3.0
              !T_fluid(i,j,k) = P(i,j,k)*M_B/rho_star(i,j,k)
              T_fluid(i,j,k) = (E_rad(i,j,k)/rad_const)**(0.25)
              tau_rad(i,j,k) = w_r*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*u0(i,j,k) + 2.0*u0(i,j,k)*F_rad0(i,j,k)-P_radl)
              S_rad_x(i,j,k) = w_r*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*ux_r + F_rad0(i,j,k)*ux_r + F_radx(i,j,k)*u0(i,j,k))
              S_rad_y(i,j,k) = w_r*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*uy_r + F_rad0(i,j,k)*uy_r + F_rady(i,j,k)*u0(i,j,k))
              S_rad_z(i,j,k) = w_r*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*uz_r + F_rad0(i,j,k)*uz_r + F_radz(i,j,k)*u0(i,j,k))
	   end if 

           ! D.J. Griffth's Introduction to Electrodynamics]
           ! Note that Ai are staggered
           A1 = 0.d0
           xs = X(i,j,k)*nx + (Y(i,j,k)+0.5d0*dY)*ny + (Z(i,j,k)+0.5d0*dZ)*nz
           ys = X(i,j,k)*yx + (Y(i,j,k)+0.5d0*dY)*yy + (Z(i,j,k)+0.5d0*dZ)*yz
           if (xs .le. 0.d0) then
              A2 = B3_l*xs
              A3 = B1_l*ys - B2_l*xs
           else
              A2 = B3_r*xs
              A3 = B1_r*ys - B2_r*xs
           end if
           Ax(i,j,k) = nx*A1 + yx*A2 + zx*A3
           xs = (X(i,j,k)+0.5d0*dX)*nx + Y(i,j,k)*ny + (Z(i,j,k)+0.5d0*dZ)*nz
           ys = (X(i,j,k)+0.5d0*dX)*yx + Y(i,j,k)*yy + (Z(i,j,k)+0.5d0*dZ)*yz
           if (xs .le. 0.d0) then
              A2 = B3_l*xs
              A3 = B1_l*ys - B2_l*xs
           else
              A2 = B3_r*xs
              A3 = B1_r*ys - B2_r*xs
           end if
           Ay(i,j,k) = ny*A1 + yy*A2 + zy*A3
           xs = (X(i,j,k)+0.5d0*dX)*nx + (Y(i,j,k)+0.5d0*dY)*ny + Z(i,j,k)*nz
           ys = (X(i,j,k)+0.5d0*dX)*yx + (Y(i,j,k)+0.5d0*dY)*yy + Z(i,j,k)*yz
           if (xs .le. 0.d0) then
              A2 = B3_l*xs
              A3 = B1_l*ys - B2_l*xs
           else
              A2 = B3_r*xs
              A3 = B1_r*ys - B2_r*xs
           end if
           Az(i,j,k) = nz*A1 + yz*A3 + zz*A3
   	end do
     end do
  end do


  else ! no radiation 

  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
	      xs = X(i,j,k)*nx + Y(i,j,k)*ny + Z(i,j,k)*nz
              if (xs .le. 0.D0) then
             !
             ! "left" side of partition

	     ! If radiation is involved, we calcuate the 1D profile (in x-direction) using initialdata_radiation and fill
	     !the 3D arrays with this 1D array.
	   
              P(i,j,k)    = P_l
              rho_b(i,j,k) = rho_l
              rho_star(i,j,k) = rho_l * w_l
              tau(i,j,k)   = tau_l
              h(i,j,k) = h_l
              w(i,j,k)    = rho_star(i,j,k) *w_l
              st_x(i,j,k)  = rho_star(i,j,k) *h_l*ux_l
              st_y(i,j,k)  = rho_star(i,j,k) *h_l*uy_l
              st_z(i,j,k)  = rho_star(i,j,k) *h_l*uz_l
	      mhd_st_x(i,j,k) = st_x(i,j,k) + w_l*sbs_l*ux_l - sb0_l*sbx_l
              mhd_st_y(i,j,k) = st_y(i,j,k) + w_l*sbs_l*uy_l - sb0_l*sby_l
              mhd_st_z(i,j,k) = st_z(i,j,k) + w_l*sbs_l*uz_l - sb0_l*sbz_l
              vx(i,j,k)   = ux_l/w_l
              vy(i,j,k)   = uy_l/w_l
              vz(i,j,k)   = uz_l/w_l
	      u0(i,j,k)   = w_l
   	      Bx(i,j,k)   = Bx_l
              By(i,j,k)   = By_l
              Bz(i,j,k)   = Bz_l
	      Bxtilde(i,j,k) = Bx_l
              Bytilde(i,j,k) = By_l
              Bztilde(i,j,k) = Bz_l
	      sbt(i,j,k) = sb0_l
	      sbx(i,j,k) = sbx_l
	      sby(i,j,k) = sby_l
	      sbz(i,j,k) = sbz_l
	      E_rad(i,j,k) = Erad_l
	      F_rad0(i,j,k) = ZERO
	      F_radx(i,j,k) = Fradx_l
	      F_rady(i,j,k) = ZERO
	      F_radz(i,j,k) = ZERO
	      P_radl = Erad_l/3.0
	      
	      !T_fluid(i,j,k) = (E_rad(i,j,k)/rad_const)**(0.25)
	      T_fluid(i,j,k) = P(i,j,k)/rho(i,j,k)
              tau_rad(i,j,k) = w_l*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*u0(i,j,k) + 2.0*u0(i,j,k)*F_rad0(i,j,k)-P_radl)
              S_rad_x(i,j,k) = w_l*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*ux_l + F_rad0(i,j,k)*ux_l + F_radx(i,j,k)*u0(i,j,k))
              S_rad_y(i,j,k) = w_l*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*uy_l + F_rad0(i,j,k)*uy_l + F_rady(i,j,k)*u0(i,j,k))
              S_rad_z(i,j,k) = w_l*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*uz_l + F_rad0(i,j,k)*uz_l + F_radz(i,j,k)*u0(i,j,k))
           else
              !
              ! "right" side of partition
              !
              P(i,j,k)    = P_r
              rho_b(i,j,k) = rho_r
              rho_star(i,j,k) = rho_r * w_r
              tau(i,j,k)   = tau_r
              h(i,j,k) = h_r
              w(i,j,k)    = rho_star(i,j,k) *w_r
              st_x(i,j,k)  = rho_star(i,j,k) *h_r*ux_r
              st_y(i,j,k)  = rho_star(i,j,k) *h_r*uy_r
              st_z(i,j,k)  = rho_star(i,j,k) *h_r*uz_r
              mhd_st_x(i,j,k) = st_x(i,j,k) + w_r*sbs_r*ux_r - sb0_r*sbx_r
              mhd_st_y(i,j,k) = st_y(i,j,k) + w_r*sbs_r*uy_r - sb0_r*sby_r
              mhd_st_z(i,j,k) = st_z(i,j,k) + w_r*sbs_r*uz_r - sb0_r*sbz_r
              vx(i,j,k)   = ux_r/w_r
              vy(i,j,k)   = uy_r/w_r
              vz(i,j,k)   = uz_r/w_r
	      u0(i,j,k)   = w_r
              Bx(i,j,k)   = Bx_r
              By(i,j,k)   = By_r
              Bz(i,j,k)   = Bz_r
              Bxtilde(i,j,k) = Bx_r
              Bytilde(i,j,k) = By_r
              Bztilde(i,j,k) = Bz_r
              sbt(i,j,k) = sb0_r
              sbx(i,j,k) = sbx_r
              sby(i,j,k) = sby_r
              sbz(i,j,k) = sbz_r
	      E_rad(i,j,k) = Erad_r	   
              F_rad0(i,j,k) = ZERO 
	      F_radx(i,j,k) = Fradx_r
  	      F_rady(i,j,k) = ZERO
	      F_radz(i,j,k) = ZERO
	      P_radl = Erad_r/3.0
	      
	      !T_fluid(i,j,k) = (E_rad(i,j,k)/rad_const)**(0.25)
	      T_fluid(i,j,k) = P(i,j,k)/rho(i,j,k)
              tau_rad(i,j,k) = w_r*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*u0(i,j,k) + 2.0*u0(i,j,k)*F_rad0(i,j,k)-P_radl)
              S_rad_x(i,j,k) = w_r*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*ux_r + F_rad0(i,j,k)*ux_r + F_radx(i,j,k)*u0(i,j,k))
              S_rad_y(i,j,k) = w_r*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*uy_r + F_rad0(i,j,k)*uy_r + F_rady(i,j,k)*u0(i,j,k))
              S_rad_z(i,j,k) = w_r*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*uz_r + F_rad0(i,j,k)*uz_r + F_radz(i,j,k)*u0(i,j,k))          
	end if

           ! Setup vector potential [use the expression of problem 5.33 in
	   ! D.J. Griffth's Introduction to Electrodynamics]
           ! Note that Ai are staggered
           A1 = 0.d0
	   xs = X(i,j,k)*nx + (Y(i,j,k)+0.5d0*dY)*ny + (Z(i,j,k)+0.5d0*dZ)*nz
           ys = X(i,j,k)*yx + (Y(i,j,k)+0.5d0*dY)*yy + (Z(i,j,k)+0.5d0*dZ)*yz
           if (xs .le. 0.d0) then
	      A2 = B3_l*xs
	      A3 = B1_l*ys - B2_l*xs
	   else
              A2 = B3_r*xs
              A3 = B1_r*ys - B2_r*xs
	   end if
	   Ax(i,j,k) = nx*A1 + yx*A2 + zx*A3
           xs = (X(i,j,k)+0.5d0*dX)*nx + Y(i,j,k)*ny + (Z(i,j,k)+0.5d0*dZ)*nz
           ys = (X(i,j,k)+0.5d0*dX)*yx + Y(i,j,k)*yy + (Z(i,j,k)+0.5d0*dZ)*yz
           if (xs .le. 0.d0) then
              A2 = B3_l*xs
              A3 = B1_l*ys - B2_l*xs
           else
              A2 = B3_r*xs
              A3 = B1_r*ys - B2_r*xs
           end if
	   Ay(i,j,k) = ny*A1 + yy*A2 + zy*A3
           xs = (X(i,j,k)+0.5d0*dX)*nx + (Y(i,j,k)+0.5d0*dY)*ny + Z(i,j,k)*nz
           ys = (X(i,j,k)+0.5d0*dX)*yx + (Y(i,j,k)+0.5d0*dY)*yy + Z(i,j,k)*yz
           if (xs .le. 0.d0) then
              A2 = B3_l*xs
              A3 = B1_l*ys - B2_l*xs
           else
              A2 = B3_r*xs
              A3 = B1_r*ys - B2_r*xs
           end if
	   Az(i,j,k) = nz*A1 + yz*A3 + zz*A3

           end do
     	end do
      end do
  end if
end if

  !tau_p = tau
  !h_p = h

  neos = 1
  write(*,*) "Assuming polytropic EOS..."
  write(*,*) "Polytropic index = ",npolyshock
  write(*,*) "Polytropic constant K = ",K_poly
  rho_tab(1) = 0.03d0
  P_tab(1) = K_poly * rho_tab(1)**(1.0+1.0/npolyshock)
  eps_tab(1) = npolyshock*P_tab(1)/rho_tab(1)
  gamma_tab(1) = 1.d0 + 1.d0/npolyshock
  k_tab(1) = K_poly
  gamma_tab(2) = 1.d0 + 1.d0/npolyshock 
  k_tab(2) = K_poly

  gamma_th = 1.d0 + 1.d0/npolyshock 


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
  !!Ex = ZERO
  !!Ey = ZERO
  !!Ez = ZERO

  !FIXME: read note:
  ! NOTE: The following lines are not strictly correct, but this should NOT matter, since Cowling is turned on!
  Sxx = ZERO
  Sxy = ZERO
  Sxz = ZERO
  Syy = ZERO
  Syz = ZERO
  Szz = ZERO

  write(*,*) "checkpoint5"
  
  Sx = ZERO
  Sy = ZERO
  Sz = ZERO

  rho = ZERO

  write(*,*) "LOCAL INITIALDATA FINISHED!",dx,mhd_st_x(4,4,4), rho_star(4,4,4), h(4,4,4), vx(4,4,4), npolyshock, P(4,4,4),rho_b(4,4,4)
  write(*,*) "LOCAL INITIALDATA FINISHED! y=",y(1,:,1)







end subroutine mhd_shock_initialdata_local

