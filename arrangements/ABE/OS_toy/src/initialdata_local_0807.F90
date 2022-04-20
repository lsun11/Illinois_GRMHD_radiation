!-------------------------------------------------------------------
! Here we set up all required variables for OS tests, including:
!  emfields, BSSN variables, primitives, etc.
!-------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine mhd_OS_initialdata_local(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  real*8                :: dT,dX,dY,dZ,P_max,rho_max,tau_max,rhos_max,rho_fail_max_step,M_fail_step,rhobatm
  real*8                :: xmin,ymin,zmin,xmax,ymax,zmax
  integer               :: proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax
  integer               :: ierr,index,handle,dummy,glob_imax,glob_jmax,glob_kmax
  CCTK_REAL             :: reduction_value

  integer               :: i,j,k
  real*8                :: xs,ys, xspsi2
  real*8, parameter :: pi = 3.1415926535897932d0
  real*8  	        :: s4pi,sb0_l,sb0_r,sb1_l,sb1_r,sb2_l,sb2_r,sb3_r,sb3_l,sbs_l,sbs_r
  real*8		:: rho_OS, r_iso, tempOS1, tempOS2, fac, fac2, psi2, psi4, psim4,F_rad0l, P_radl, u_x, u_y, u_z
  real*8                :: zeta, chi_shock, Fksq, C_A, C_B, P_rad0xl, P_rad0yl, P_rad0zl, P_rad00l, P_rad0_xl, P_rad0_yl, P_rad0_zl




  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  integer :: ONE,ZERO
  parameter(ONE = 1.D0, ZERO = 0.D0)

  write(*,*) "Start OS_toy initialdata_local!!!!!"

  if(primitives_solver.eq.0 .or. enable_OS_collapse.ne.1) then
     write(*,*) "ERROR.  YOU MUST SET primitives_solver=1 or 2 AND enable_OS_collapse=1 FOR OS_toy TESTS"
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

  write(*,*) "Okay... Looks like you set up the grid correctly EXTENTS:",cctk_lsh
  write(*,*) "Okay... Looks like you set up the grid correctly GZS:",cctk_nghostzones

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  !initialize all this stuff                                                                                             
  P=ZERO
  rho_star=ZERO
  rhoYe=ZERO
  h=ZERO
  w=ZERO
  tau=ZERO
  st_x=ZERO
  st_y=ZERO
  st_z=ZERO
  tau_rad=ZERO
  S_rad_x=ZERO
  S_rad_y=ZERO
  S_rad_z=ZERO

  sbt=ZERO
  sbx=ZERO
  sby=ZERO
  sbz=ZERO
  Bx=ZERO
  By=ZERO
  Bz=ZERO
  Ex=ZERO
  Ey=ZERO
  Ez=ZERO
  mhd_st_x=ZERO
  mhd_st_y=ZERO
  mhd_st_z=ZERO

  ! Set excision_zone_gf to avoid valgrind memory errors                                                                 
  excision_zone_gf = 0


  write(*,*) "Start generating OS ID!!!! M=",M_OS
  write(*,*) "Start generating OS ID!!!! R=",R_OS
  write(*,*) "Start generating OS ID!!!! x,y,z=",xc_OS, yc_OS, zc_OS
  write(*,*) "Start generating OS ID!!!! P_over_rho=",P_over_rho
  write(*,*) "Start generating OS ID!!!! E_over_rho=",E_over_rho



  rho_OS = 3.0*M_OS/(4.0*pi*R_OS*R_OS*R_OS)

  r_iso = R_OS * (1.0 - M_OS/R_OS + sqrt(1.0 - 2.0*M_OS/R_OS))/2.0

  write(*,*) "Stellar density and isotropic radius are", rho_OS, r_iso


  !!!!First, set the BSSN variables!!!!


  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
           xs = sqrt((X(i,j,k)-xc_OS)**2 + (Y(i,j,k)-yc_OS)**2 + (Z(i,j,k)-zc_OS)**2)
	   if ( xs .gt. r_iso) then
	      psi(i,j,k) = 1.0 + M_OS/(2.0 *xs) 
	   else
	      tempOS1 = (1.0 + sqrt (1.0 - 2.0 * M_OS /R_OS)) * r_iso * R_OS**2
	      tempOS2 = 2.0 * r_iso**3 + M_OS * xs**2 	
	      psi(i,j,k) = sqrt (tempOS1/tempOS2)
	   end if
	   
	   phi(i,j,k) = log(psi(i,j,k))
	end do   
     end do
  end do
  !======================================
  ! Set other metric data to Minkowski
  !======================================
  phix = ZERO
  phiy = ZERO
  phiz = ZERO
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
  ! Set derivatives to zero                                                                                                                                            
  !======================================                                
  gxxx = 0.D0
  gxxy = 0.D0
  gxxz = 0.D0
  gxyx = 0.D0
  gxyy = 0.D0
  gxyz = 0.D0
  gxzx = 0.D0
  gxzy = 0.D0
  gxzz = 0.D0
  gyyx = 0.D0
  gyyy = 0.D0
  gyyz = 0.D0
  gyzx = 0.D0
  gyzy = 0.D0
  gyzz = 0.D0
  gzzx = 0.D0
  gzzy = 0.D0
  gzzz = 0.D0

  Gammaxxx = 0.D0
  Gammaxxy = 0.D0
  Gammaxxz = 0.D0
  Gammaxyy = 0.D0
  Gammaxyz = 0.D0
  Gammaxzz = 0.D0
  Gammayxx = 0.D0
  Gammayxy = 0.D0
  Gammayxz = 0.D0
  Gammayyy = 0.D0
  Gammayyz = 0.D0
  Gammayzz = 0.D0
  Gammazxx = 0.D0
  Gammazxy = 0.D0
  Gammazxz = 0.D0
  Gammazyy = 0.D0
  Gammazyz = 0.D0
  Gammazzz = 0.D0
  !======================================   
  ! Set BSSN matter sources to zero!                                                                                                                                                          
  !======================================   
  rho = ZERO
  S = ZERO
  Sx = ZERO
  Sy = ZERO
  Sz = ZERO
  Sxx = ZERO
  Sxy = ZERO
  Sxz = ZERO
  Syy = ZERO
  Syz = ZERO
  Szz = ZERO
  !======================================
  ! Set everything else to Zero!
  !======================================
  shiftx = ZERO
  shifty = ZERO
  shiftz = ZERO
  shiftxt = 0.D0
  shiftyt = 0.D0
  shiftzt = 0.D0
  lapset = 0.D0


  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
           Kxx(i,j,k) = 0.d0
           Kxy(i,j,k) = 0.d0
           Kxz(i,j,k) = 0.d0
           Kyy(i,j,k) = 0.d0
           Kyz(i,j,k) = 0.d0
           Kzz(i,j,k) = 0.d0
        end do
     end do
  end do

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')
  call CartSymGN(dummy,cctkGH,'shift::shift_vars')




  !!!!Second, set rho_OS, P, and E_rad (only non-zero matter varibales)!!!!
  !!! Here rhobatm = rho_fact * rho_OS, since rho_OS = rho_max at t = 0
  rhobatm = rho_fact * rho_OS

  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
           psi2=psi(i,j,k)*psi(i,j,k) 
           psi4=psi2*psi2
           xs = sqrt((X(i,j,k)-xc_OS)**2 + (Y(i,j,k)-yc_OS)**2 + (Z(i,j,k)-zc_OS)**2)          
           xspsi2 = xs*psi2

!!! Set other to zero (or one)
           vx(i,j,k)   = 0.d0
           vy(i,j,k)   = 0.d0
           vz(i,j,k)   = 0.d0
           u0(i,j,k)   = 1.d0
           Bx(i,j,k)   = 0.d0
           By(i,j,k)   = 0.d0
           Bz(i,j,k)   = 0.d0
           Bxtilde(i,j,k)   = 0.d0
           Bytilde(i,j,k)   = 0.d0
           Bztilde(i,j,k)   = 0.d0
           sbt(i,j,k) = 0.d0
           sbx(i,j,k) = 0.d0
           sby(i,j,k) = 0.d0
           sbz(i,j,k) = 0.d0
              

           if ( xs .lt. r_iso) then
              rho_b(i,j,k) = rho_OS
              P(i,j,k) = P_over_rho * rho_b(i,j,k)**(gamma_th)
              
              if (rad_evolve_enable .eq. 1) then               
                 E_rad(i,j,k) = P(i,j,k)/P_over_rho*E_over_rho
              else
                 E_rad(i,j,k) = 0.d0
              end if
              
           else
              !rho_b(i,j,k) = rho_OS/(exp((xs-r_iso)/rounding)+1.0)+rhobatm*(1.0-1.0/(exp((xs-r_iso)/rounding)+1.0))
              rho_b(i,j,k) = rhobatm
              P(i,j,k) = 0.d0              
              if (rad_closure_scheme .eq. 0) then
                 !   P(i,j,k) = 0.d0
                 E_rad(i,j,k) = 0.d0
              else
                 !   P(i,j,k) = P_over_rho * rho_b(i,j,k)**(gamma_th)
                 E_rad(i,j,k) = E_over_rho*rho_b(i,j,k)**(gamma_th)
              end if
           end if
           
           
           Y_e(i,j,k)=0.5d0
           !              P(i,j,k) = P_over_rho * rho_b(i,j,k)**(gamma_th)
           rho_star(i,j,k) = psi(i,j,k)**6*rho_b(i,j,k)*u0(i,j,k)
           w(i,j,k)    = rho_star(i,j,k) * u0(i,j,k)
           rhoYe(i,j,k) = rho_star(i,j,k)*Y_e(i,j,k)
           
           !             New setup             
           !h(i,j,k) = 1.0 + gamma_th/(gamma_th-1.0)*P(i,j,k)/rho_b(i,j,k) 
           h(i,j,k) = 1.0 + P_over_rho * rho_b(i,j,k)**(gamma_th-1.0)
           tau(i,j,k)   = w(i,j,k)*h(i,j,k)-psi(i,j,k)**6*P(i,j,k)-rho_star(i,j,k)
           
           st_x(i,j,k)  = 0.d0
           st_y(i,j,k)  = 0.d0
           st_z(i,j,k)  = 0.d0
           mhd_st_x(i,j,k) = 0.d0
           mhd_st_y(i,j,k) = 0.d0
           mhd_st_z(i,j,k) = 0.d0
           
           
           Ax(i,j,k) = 0.d0
           Ay(i,j,k) = 0.d0
           Az(i,j,k) = 0.d0
           
           T_fluid(i,j,k) = P(i,j,k)/rho_b(i,j,k)
           
           u_x = u0(i,j,k)*psi(i,j,k)**4*(gxx(i,j,k)*(shiftx(i,j,k)+vx(i,j,k))+gxy(i,j,k)*(shifty(i,j,k)+vx(i,j,k))+gxz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)))
           u_y = u0(i,j,k)*psi(i,j,k)**4*(gxy(i,j,k)*(shiftx(i,j,k)+vx(i,j,k))+gyy(i,j,k)*(shifty(i,j,k)+vx(i,j,k))+gyz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)))
           u_z = u0(i,j,k)*psi(i,j,k)**4*(gxz(i,j,k)*(shiftx(i,j,k)+vx(i,j,k))+gyz(i,j,k)*(shifty(i,j,k)+vx(i,j,k))+gzz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)))
           
           if (rad_evolve_enable .eq. 1) then
             P_radl = E_rad(i,j,k)/3.d0
             F_rad0(i,j,k) = 0.d0
             F_radx(i,j,k) = 0.d0
             F_rady(i,j,k) = 0.d0
             F_radz(i,j,k) = 0.d0
             psim4 = 1.0/psi4 
             
             F_rad_scalar(i,j,k) =  sqrt(psi4*(gxx(i,j,k)*F_radx(i,j,k)**2 + gyy(i,j,k)*F_rady(i,j,k)**2 + gzz(i,j,k)*F_radz(i,j,k)**2 + &
                  2.0*(gxy(i,j,k)*F_radx(i,j,k)*F_rady(i,j,k) + gxz(i,j,k)*F_radx(i,j,k)*F_radz(i,j,k) + gyz(i,j,k)*F_rady(i,j,k)*F_radz(i,j,k) ) ) )
             
             
             if(rad_closure_scheme .eq. 0) then
                tau_rad(i,j,k) = psi(i,j,k)**6*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*u0(i,j,k) + 2.0*u0(i,j,k)*F_rad0(i,j,k)-P_radl)
                S_rad_x(i,j,k) = psi(i,j,k)**6*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*u0(i,j,k)*vx(i,j,k) + F_rad0(i,j,k)*u0(i,j,k)*vx(i,j,k) + F_radx(i,j,k)*u0(i,j,k))
                S_rad_y(i,j,k) = psi(i,j,k)**6*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*u0(i,j,k)*vy(i,j,k) + F_rad0(i,j,k)*u0(i,j,k)*vy(i,j,k) + F_rady(i,j,k)*u0(i,j,k))
                S_rad_z(i,j,k) = psi(i,j,k)**6*((E_rad(i,j,k) + P_radl)*u0(i,j,k)*u0(i,j,k)*vz(i,j,k) + F_rad0(i,j,k)*u0(i,j,k)*vz(i,j,k) + F_radz(i,j,k)*u0(i,j,k))
             else
                ! At t=0, we set metric to Minkowsiki. So F^0=-F_0, F^i=F_i 
                F_radx(i,j,k) = 0.0e-18
                F_rady(i,j,k) = 0.0e-18
                F_radz(i,j,k) = 0.0e-18
                F_rad0(i,j,k) = F_radx(i,j,k)*u_x + F_rady(i,j,k)*u_y + F_radz(i,j,k)*u_z 

               if (E_rad(i,j,k).le. 0.d0) then
                   zeta=0.d0
                else
                   zeta = sqrt((-F_rad0(i,j,k)**2 + F_radx(i,j,k)**2 + F_rady(i,j,k)**2 + F_radz(i,j,k)**2)/E_rad(i,j,k)**2)
                end if
                chi_shock  = 1.0/3.0 + zeta**2*(6.0-2.0*zeta+6.0*zeta**2)/15.0
                !chi_shock = 1.0/3.0
                chi_rad(i,j,k) = chi_shock
                Fksq = F_radx(i,j,k)**2 + F_rady(i,j,k)**2 + F_radz(i,j,k)**2
                C_A = (3.0*chi_shock -1.0)/2.0
                C_B = 1.5*(1.0-chi_shock)
                if (F_rad_scalar(i,j,k) .eq. 0.0) then
                   P_radxx(i,j,k) = E_rad(i,j,k)*(psim4*gupxx(i,j,k)+(u0(i,j,k)*vx(i,j,k))**2)/3.0*C_B
                   P_radyy(i,j,k) = E_rad(i,j,k)*(psim4*gupyy(i,j,k)+(u0(i,j,k)*vy(i,j,k))**2)/3.0*C_B
                   P_radzz(i,j,k) = E_rad(i,j,k)*(psim4*gupzz(i,j,k)+(u0(i,j,k)*vz(i,j,k))**2)/3.0*C_B
                   P_radxy(i,j,k) = E_rad(i,j,k)*(psim4*gupxy(i,j,k)+u0(i,j,k)**2*vx(i,j,k)*vy(i,j,k))/3.0*C_B
                   P_radxz(i,j,k) = E_rad(i,j,k)*(psim4*gupxz(i,j,k)+u0(i,j,k)**2*vx(i,j,k)*vz(i,j,k))/3.0*C_B
                   P_radyz(i,j,k) = E_rad(i,j,k)*(psim4*gupyz(i,j,k)+u0(i,j,k)**2*vy(i,j,k)*vz(i,j,k))/3.0*C_B
                else
                   P_radxx(i,j,k) = E_rad(i,j,k)*(psim4*gupxx(i,j,k)+(u0(i,j,k)*vx(i,j,k))**2)*( C_A*psi4*gxx(i,j,k)*F_radx(i,j,k)**2/F_rad_scalar(i,j,k)**2 + C_B/3.0 )
                   P_radyy(i,j,k) = E_rad(i,j,k)*(psim4*gupyy(i,j,k)+(u0(i,j,k)*vy(i,j,k))**2)*( C_A*psi4*gyy(i,j,k)*F_rady(i,j,k)**2/F_rad_scalar(i,j,k)**2 + C_B/3.0 )
                   P_radzz(i,j,k) = E_rad(i,j,k)*(psim4*gupzz(i,j,k)+(u0(i,j,k)*vz(i,j,k))**2)*( C_A*psi4*gzz(i,j,k)*F_radz(i,j,k)**2/F_rad_scalar(i,j,k)**2 + C_B/3.0 )
                   P_radxy(i,j,k) = E_rad(i,j,k)*(psim4*gupxy(i,j,k)+(u0(i,j,k)**2*vx(i,j,k)*vy(i,j,k)))*( C_A*psi4*gxy(i,j,k)*F_radx(i,j,k)*F_rady(i,j,k)/F_rad_scalar(i,j,k)**2 + C_B/3.0 )
                   P_radxz(i,j,k) = E_rad(i,j,k)*(psim4*gupxz(i,j,k)+(u0(i,j,k)**2*vx(i,j,k)*vz(i,j,k)))*( C_A*psi4*gxz(i,j,k)*F_radx(i,j,k)*F_radz(i,j,k)/F_rad_scalar(i,j,k)**2 + C_B/3.0 )
                   P_radyz(i,j,k) = E_rad(i,j,k)*(psim4*gupyz(i,j,k)+(u0(i,j,k)**2*vy(i,j,k)*vz(i,j,k)))*( C_A*psi4*gyz(i,j,k)*F_rady(i,j,k)*F_radz(i,j,k)/F_rad_scalar(i,j,k)**2 + C_B/3.0 )
                end if
                
                P_rad0xl = P_radxx(i,j,k)*vx(i,j,k) + P_radxy(i,j,k) * vy(i,j,k) + P_radxz(i,j,k) * vz(i,j,k)
                P_rad0yl = P_radxy(i,j,k)*vx(i,j,k) + P_radyy(i,j,k) * vy(i,j,k) + P_radyz(i,j,k) * vz(i,j,k)
                P_rad0yl = P_radxz(i,j,k)*vx(i,j,k) + P_radyz(i,j,k) * vy(i,j,k) + P_radzz(i,j,k) * vz(i,j,k)
                P_rad00l = P_rad0xl*vx(i,j,k) + P_rad0yl*vy(i,j,k)  + P_rad0zl*vz(i,j,k)
                
                ! shift_i = 0.0
                P_rad0_xl = P_rad00l*0.0 + psi4 * (P_rad0xl*gxx(i,j,k) + P_rad0yl*gxy(i,j,k) + P_rad0zl*gxz(i,j,k))
                P_rad0_yl = P_rad00l*0.0 + psi4 * (P_rad0xl*gxy(i,j,k) + P_rad0yl*gyy(i,j,k) + P_rad0zl*gyz(i,j,k))
                P_rad0_zl = P_rad00l*0.0 + psi4 * (P_rad0xl*gxz(i,j,k) + P_rad0yl*gyz(i,j,k) + P_rad0zl*gzz(i,j,k))
                
                tau_rad(i,j,k) = psi(i,j,k)**6*(E_rad(i,j,k)*u0(i,j,k)*u0(i,j,k) + 2.0*u0(i,j,k)*F_rad0(i,j,k) + P_rad00l)
                S_rad_x(i,j,k) = psi(i,j,k)**6*(E_rad(i,j,k)*u0(i,j,k)*u0(i,j,k)*vx(i,j,k) + F_rad0(i,j,k)*u0(i,j,k)*vx(i,j,k) &
                     + F_radx(i,j,k)*u0(i,j,k) + P_rad0_xl)
                S_rad_y(i,j,k) = psi(i,j,k)**6*(E_rad(i,j,k)*u0(i,j,k)*u0(i,j,k)*vy(i,j,k) + F_rad0(i,j,k)*u0(i,j,k)*vy(i,j,k) &
                     + F_rady(i,j,k)*u0(i,j,k) + P_rad0_yl)
                S_rad_z(i,j,k) = psi(i,j,k)**6*(E_rad(i,j,k)*u0(i,j,k)*u0(i,j,k)*vz(i,j,k) + F_rad0(i,j,k)*u0(i,j,k)*vz(i,j,k) &
                     + F_radz(i,j,k)*u0(i,j,k) + P_rad0_zl)
                
                if (isnan(tau_rad(i,j,k))) then
                   write(*,*) "In OS_initialdata_local, tau_rad is nan!!!! i,j,k are", i,j,k
                   write(*,*) "P_rad00l, P_radxx(i,j,k), zeta, E_rad(i,j,k) = ",P_rad00l, P_radxx(i,j,k), zeta, E_rad(i,j,k)
                end if
             end if
          else
             tau_rad(i,j,k) = 0.d0
             S_rad_x(i,j,k) = 0.d0
             S_rad_y(i,j,k) = 0.d0
             S_rad_z(i,j,k) = 0.d0	  
          end if
          
          rho (i,j,k) = h(i,j,k)*w(i,j,k)*psi(i,j,k)**(-6.0)-P(i,j,k) + E_rad(i,j,k)
          
          Sxx(i,j,k) = psi(i,j,k)**4 * (P(i,j,k) + P_radl)
          Syy(i,j,k) = psi(i,j,k)**4 * (P(i,j,k) + P_radl)
          Szz(i,j,k) = psi(i,j,k)**4 * (P(i,j,k) + P_radl)
          
          
          if(isnan(tau_rad(i,j,k))) write(*,*) "Inside OS_toy ID. tau_rad is NAN, i,j,k=",i,j,k
          if(isnan(S_rad_x(i,j,k))) write(*,*) "Inside OS_toy ID. S_rad_x is NAN, i,j,k=",i,j,k
          if(isnan(P_radxx(i,j,k))) write(*,*) "Inside OS_toy ID. P_radxx is NAN, i,j,k=",i,j,k
          if(isnan(Sxx(i,j,k))) write(*,*) "Inside OS_toy ID. Sxx(i,j,k) is NAN, i,j,k=",i,j,k
          
       end do
    end do
 end do
 
!!!Since Everywhere is in LTE. We can compute rad_const anywhere
!!! To make sure rad_const is constant everywhere, rad_const = E/T^4 = E_over_rho * rho**gamma / (P_over_rho* rho**(gamma-1))**4
!!! = E_over_rho/P_over_rho**4  * rho**(4-3*gamma) = const. --> gamma = 4/3

!  rad_const = E_over_rho/(P_over_rho)**4*M_OS/(4.d0*pi*R_OS**3/3.d0)
  rad_const = E_over_rho/P_over_rho**4 * (rho_OS)**(4.0-3.0*gamma_th)
!  rad_const = E_over_rho/P_over_rho**4 * (rho_OS)**(5.0-4.0*gamma_th)
  write(*,*) "rad_const is", rad_const
!  write(*,*) "Check: E_over_rho/P_over_rho * rho_b**(4.0-3*gamma_th) =", E_over_rho/P_over_rho**4 * (M_OS/(4.d0*pi*R_OS**3/3.d0) )**(4.0-3.0*gamma_th)

  neos = 1
  write(*,*) "Assuming polytropic EOS..."
  write(*,*) "Polytropic constant K = ",P_over_rho
  rho_tab(1) = rho_OS*100.0

!  P_tab(1) = P_over_rho*rho_tab(1)**(gamma_th)
!  eps_tab(1) = 1.0/(gamma_th - 1.0) * P_tab(1)/rho_tab(1)
!  gamma_tab(1) = gamma_th
!  k_tab(1) = P_over_rho
!  gamma_tab(2) = gamma_th
!  k_tab(2) = P_over_rho

  P_tab(1) = P_over_rho*rho_tab(1)
  eps_tab(1) = 0.0
  gamma_tab(1) = gamma_th
  k_tab(1) = P_over_rho
  gamma_tab(2) = gamma_th
  k_tab(2) = P_over_rho

  write(*,*) "rho_tab(1), P_tab(1), gamma_tab(1), k_tab(1) = ", rho_tab(1), P_tab(1), gamma_tab(1), k_tab(1)

  !gamma_th = 0.d0 




  !======================================
  ! Set initial EM fields to zero:
  !======================================

  Ex = By*vz - Bz*vy
  Ey = Bz*vx - Bx*vz
  Ez = Bx*vy - By*vx
  !!Ex = ZERO
  !!Ey = ZERO
  !!Ez = ZERO


  ! With P = 0, Sx=Sy=Sz = 0, Sxy=Syz=Sxz=0, however, Sxx, Syy, Szz and rho are not zero!!!!
  
!  Sxx = ZERO
  Sxy = ZERO
  Sxz = ZERO
!  Syy = ZERO
  Syz = ZERO
!  Szz = ZERO

  Sx = ZERO
  Sy = ZERO
  Sz = ZERO

!    !particle tracer stuff                                                                                                 
  do index=1,narr
     coord(index,1)=0.d0  
     coord(index,2)=index*(particle_rad_cut*r_iso/narr)
     coord(index,3)=0.d0                                                                                 
     coord(index,4)=0.d0                                                 
  end do    




  write(*,*) "LOCAL INITIALDATA FINISHED!",dx,mhd_st_x(4,4,4)
  write(*,*) "LOCAL INITIALDATA FINISHED! y=",y(1,:,1)


end subroutine mhd_OS_initialdata_local
