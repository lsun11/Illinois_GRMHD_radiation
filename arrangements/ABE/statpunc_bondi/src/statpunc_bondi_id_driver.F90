#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine statpunc_bondi_initialdata_driver(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  INTERFACE
     SUBROUTINE mb_get_rho0(rho0_derivs_s,Rs_s,a2_s,rho0_s,M,Rs,Mdot,K_poly,Gamma,G29_rhs,rho0,a2)
       USE nrtype
       IMPLICIT NONE
       real(DP), INTENT(IN) :: Rs_s,a2_s,rho0_s,M,Rs,Mdot,K_poly,Gamma,G29_rhs,rho0_derivs_s
       real(DP), INTENT(OUT)             :: rho0,a2
     END SUBROUTINE mb_get_rho0
  END INTERFACE
  INTERFACE
     SUBROUTINE sonic_statpunc(M,a2_inf,gamma_th,a2_s,Rs_s,u2_s)
       IMPLICIT NONE
       real*8 :: M,a2_inf,gamma_th,a2_s,Rs_s,u2_s,a2_s_prev
     END SUBROUTINE sonic_statpunc
  END INTERFACE
  INTERFACE
     FUNCTION r_star_func_statpunc(x,M,Mdot,rho0_star,a2_inf,a2_star,gamma_th)
       IMPLICIT NONE
       real*8 :: x,M,Mdot,rho0_star,a2_inf,a2_star,gamma_th,r_star_func_statpunc
     END FUNCTION r_star_func_statpunc
  END INTERFACE

  integer, dimension(3) :: ext
  real*8 		:: dX,dY,dZ
  real*8                :: X_stat_frame
  integer               :: i,j,k
  character             :: varname*30
  real*8 :: u2_s,a2_s,a2_s_prev,G29_rhs,Rs,Ri,rho0,rho0_s,a2,Rs_a,Ri_s,Rs_s,a2_a,rho0_a,u2_a,rho0_derivi_a,rho0_derivs_s,Delta,us,nx,ny,nz,psi2,psi4,alph,ui,u0L,drhosurf,rr,gamv,sol
  real*8                :: K_poly_outer,G29_rhs_outer,a2_star_outer,a2_star_inner
  real*8                :: rho0_star,P_star,R_star,u2_star
  real*8                :: f_L,f_R,f_mid
  real*8                :: mb_o_me
  real*8                :: bh_mass_rescale
  REAL*8, PARAMETER :: PI_D=3.141592653589793238462643383279502884197

  bh_mass_rescale = bh_mass * mass_rescale_factor
  mb_o_me = 1.837152755d3 !ratio of baryon mass to electron mass

  !Delta = 0.1d0
  Delta = 0.1
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ext = cctk_lsh

  !set boost velocity
  v_boost = 0.d0
  !if (use_adm .eq. 1) then
  !   v_boost = bh_px_plus/bh_mass_rescale
  !endif

  !add wind
  v_boost=v_boost+wind

  ! write(*,*) "v_boost: ",v_boost

  !Set initial BH position to zero.
  bh_posn_x = xh0
  bh_posn_y = yh0
  bh_posn_z = zh0

  !Note that a capital Rs denotes schwarzchild radius, and Ri denotes isotropic radius
  !a2_inf and rho_inf are given in .par file

  !use eq G.28 to get K_poly
  K_poly_outer = a2_inf*rho0_inf**(1.d0-gamma_th_outer)/(gamma_th_outer - gamma_th_outer/(gamma_th_outer-1.d0)*a2_inf)

  !set this just in case
  n_poly = 1.d0/(gamma_th-1.d0)

  !get info about sonic point
  call sonic_statpunc(bh_mass_rescale,a2_inf,gamma_th_outer,a2_s,Rs_s,u2_s) !check this

  !now we can set density at sonic radius
  rho0_s=(1.d0/K_poly_outer*a2_s/(gamma_th_outer - gamma_th_outer/(gamma_th_outer-1.d0)*a2_s))**(1.d0/(gamma_th_outer-1.d0))

  !find derivative of rho0 with respect to Schwarzschild radius, evaluated at sonic radius
  rho0_derivs_s = -2.d0 *rho0_s/Rs_s 

  !now we can calculate Mdot
  Mdot=4.d0*PI_D*rho0_s*sqrt(u2_s)*Rs_s*Rs_s

  !this conversion between coords concerns me.  How can we do better?
  !compute schwarszhild adjustment radius
  Rs_a=Ri_a*(1.d0+bh_mass_rescale/(2.d0*Ri_a))**2

  !set constant on rhs of G.29
  G29_rhs_outer = (1.d0+a2_inf/(gamma_th_outer-1.d0-a2_inf))**2

  a2_star_outer = gamma_th_outer*(4.d0/3.d0)/mb_o_me/(1.d0+gamma_th_outer/(gamma_th_outer-1.d0)*(4.d0/3.d0)/mb_o_me)

  if (const_gamma .eq. 1) then
     if ((gamma_th_outer-gamma_th)/gamma_th .gt. 1.d-4) then
        write(*,*) "you set const_gamma=1, so make sure you actually set gamma_th=gamma_th_outer"
        stop
     endif
     G29_rhs = G29_rhs_outer
     K_poly = K_poly_outer
  else
     rho0_star = (1.d0/K_poly_outer*a2_star_outer/(gamma_th_outer - gamma_th_outer/(gamma_th_outer-1.d0)*a2_star_outer))**(1.d0/(gamma_th_outer-1.d0))

     P_star = K_poly_outer*rho0_star**gamma_th_outer

     R_star = 300.d0
     dx=200.d0

     f_L=r_star_func_statpunc(R_star,bh_mass_rescale,Mdot,rho0_star,a2_inf,a2_star_outer,gamma_th_outer)
     f_R=r_star_func_statpunc(R_star+dx,bh_mass_rescale,Mdot,rho0_star,a2_inf,a2_star_outer,gamma_th_outer)

     do 
        dx=dx/2.d0
        f_mid = r_star_func_statpunc(R_star+dx,bh_mass_rescale,Mdot,rho0_star,a2_inf,a2_star_outer,gamma_th_outer)
        if (f_mid*f_L .gt. 0.d0) R_star=R_star+dx
        if (dx .lt. 1.e-16) exit
     end do
     write(*,*) "R_star: ",R_star
     u2_star =  (Mdot/(4.d0*PI_D*rho0_star*R_star*R_star))**2
     write(*,*) "u2_star: ",u2_star

     K_poly = P_star/rho0_star**gamma_th
     a2_star_inner = gamma_th*K_poly*rho0_star**(gamma_th-1.d0)&
          /(1.d0 + gamma_th/(gamma_th-1.d0)*K_poly*rho0_star**(gamma_th-1.d0))
     write(*,*) "a2_star_inner: ",a2_star_inner
     
     G29_rhs = (1.d0-2.d0*bh_mass_rescale/R_star+u2_star)*(1.d0+a2_star_inner/(gamma_th-1.d0-a2_star_inner))**2


  endif
  

  !find sound speed at adjustment radius
  call mb_get_rho0(rho0_derivs_s,Rs_s,a2_s,rho0_s,bh_mass_rescale,Rs_a,Mdot,K_poly,gamma_th,G29_rhs,rho0_a,a2_a)
write(*,*) "rho0_a: ",rho0_a
  a2_a = gamma_th*K_poly*rho0_a**(gamma_th-1.d0)/(1.d0+gamma_th/(gamma_th-1.d0)*K_poly*rho0_a**(gamma_th-1.d0))
write(*,*) "a2_a: ",a2_a
  !find u^2 at adjustment radius
  u2_a=(Mdot/(4.d0*PI_D*rho0_a*Rs_a*Rs_a))**2

  !find derivative of rho0 with respect to isotropic radius, evaluated at adjustment radius
  rho0_derivi_a=-(1.d0+bh_mass_rescale/(2.d0*Ri_a))*(1.d0-bh_mass_rescale/(2.d0*Ri_a))*(rho0_a/Rs_a)*(2.d0*u2_a-bh_mass_rescale/Rs_a)/(u2_a-(1.d0-2.d0*bh_mass_rescale/Rs_a+u2_a)*a2_a)
  
  i=0
  j=0
  k=0
  do k=1,ext(3)
     !write(*,*) "i: ",i,"j: ",j,"k: ",k,"X: ",X(i,j,k),"Y: ",Y(i,j,k),"Z: ",Z(i,j,k) 
     do j=1,ext(2)
        do i=1,ext(1)
           !decide whether or not to try and correct for lortentz contraction
           if (lorentz_corr .eq. 1) then
              x_stat_frame = 1.d0/sqrt(1.d0-v_boost**2) * (X(i,j,k)-xh0)
           else
              x_stat_frame = X(i,j,k)
           endif

           !set psi2,psi4    
           psi2 = exp(2.d0*phi(i,j,k))
           psi4 = psi2*psi2
           !set lapse
           alph = lapm1(i,j,k)+1.0

           
           !dist from BH 1
           Ri=sqrt(x_stat_frame*x_stat_frame+Y(i,j,k)*Y(i,j,k)+Z(i,j,k)*Z(i,j,k))
           Rs=Ri*(1.d0+bh_mass_rescale/(2.d0*Ri))**2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           if (Ri .gt. Ri_a) then
           !   write(*,*) "test0"
           !   write(*,*) "rho0_derivs_s: ",rho0_derivs_s
           !   write(*,*) "Rs_s: ",Rs_s
           !   write(*,*) "a2_s: ",a2_s
           !   write(*,*) "rho0_s: ",rho0_s
           !   write(*,*) "bh_mass_rescale: ",bh_mass_rescale
           !   write(*,*) "Rs: ",Rs
           !   write(*,*) "Mdot: ",Mdot
           !   write(*,*) "K_poly: ",K_poly
           !   write(*,*) "gamma_th: ",gamma_th
           !   write(*,*) "G29_rhs: ",G29_rhs
              call mb_get_rho0(rho0_derivs_s,Rs_s,a2_s,rho0_s,bh_mass_rescale,Rs,Mdot,K_poly,gamma_th,G29_rhs,rho0,a2)
              
              !radial component of 4-vel in schwarzschild coords
              us=Mdot/(4*PI_D*rho0*Rs*Rs)
              !radial component of 4-vel in isotropic coords
              ui=us/((1.d0-bh_mass_rescale/(2.d0*Ri))*(1.d0+bh_mass_rescale/(2.d0*Ri)))
              !set u0L
              u0L = sqrt(1.D0+ui**2*psi4)/alph
              gamv=1.0/sqrt(1.0-v_boost**2)
              sol=alph/psi2
           endif

           !if ((Ri .le. Ri_a) .and. (Ri .ge. bh_mass_rescale/2.d0)) then
           if (Ri .le. Ri_a) then
           !rho0=rho0_a+rho0_derivi_a*(Ri_a-Ri-(Ri_a*Ri_a-Ri*Ri)/Ri_a)
              rho0=rho0_a+rho0_derivi_a*(Ri*Ri-Ri_a*Ri_a)/(2.d0*Ri_a)
              a2=gamma_th*K_poly*rho0**(gamma_th-1.d0)/&
                   (1.d0+gamma_th/(gamma_th-1.d0)*K_poly*rho0**(gamma_th-1.d0))
              !radial component of 4-vel in schwarzschild coords
              us=Mdot/(4.d0*PI_D*Rs_a*Rs_a*rho0_a)*Ri/Ri_a
              !radial component of 4-vel in isotropic coords
              ui=us/((1.d0-bh_mass_rescale/(2.d0*Ri_a))*(1.d0+bh_mass_rescale/(2.d0*Ri_a)))
              !set u0 
              u0L = sqrt(1.D0+ui**2*psi4)/alph  
              gamv=1.0/sqrt(1.0-v_boost**2)
              sol=alph/psi2
           endif
           !I don't think we need this stuff
           if (1 .eq. 0) then
              if (Ri .lt. bh_mass_rescale/2.d0) then
                 rho0=(rho0_a+rho0_derivi_a*Ri_a/4.d0)*0.5d0*(1.d0+Delta-(1.d0-Delta)*cos(2.d0*PI_D*Ri/Ri_a))
                 a2=gamma_th*K_poly*rho0**(gamma_th-1.d0)/(1.d0+gamma_th/(gamma_th-1.d0)*K_poly*rho0**(gamma_th-1.d0))
                 !radial component of 4-vel in schwarzschild coords
                 us=Mdot/(4.d0*PI_D*rs_a*rs_a*rho0_a)*Ri/Ri_a
                 !radial component of 4-vel in isotropic coords
                 ui=us/((1.d0-bh_mass_rescale/(2.d0*Ri_a))*(1.d0+bh_mass_rescale/(2.d0*Ri_a)))
                 !set psi2,psi4
                 psi2 = exp(2.d0*phi(i,j,k))
                 psi4 = psi2*psi2
                 !set lapse
                 alph = lapm1(i,j,k)+1.0
                 !set u0
                 u0L = sqrt(1.D0+ui**2*psi4)/alph
                 gamv=1.0/sqrt(1.0-v_boost**2)
                 sol=alph/psi2
              endif
           endif

           !components of radial unit vector
           if (Ri .gt. 0) then
              nx=x_stat_frame/Ri
              ny=Y(i,j,k)/Ri
              nz=Z(i,j,k)/Ri
           else
              nx=0
              ny=0
              nz=0
           endif
           !josh's method
           if (sol_method .eq. 1) then
              !components of 3-vel in isotropic coords
              vx(i,j,k)=(sol*v_boost-nx*ui/u0L)/(1.d0-v_boost*ui*nx/u0L/sol)
              vy(i,j,k)=-ny*ui/u0L/gamv/(1.d0-v_boost*ui*nx/u0L/sol)
              vz(i,j,k)=-nz*ui/u0L/gamv/(1.d0-v_boost*ui*nx/u0L/sol)
              u0(i,j,k)=(alph**2-psi4*(vx(i,j,k)*vx(i,j,k)+vy(i,j,k)*vy(i,j,k)+vz(i,j,k)*vz(i,j,k)))**-0.5
           endif
           !try setting speed of light to 1
           if (sol_method .eq. 2) then
              sol=1.d0
              !components of 3-vel in isotropic coords
              vx(i,j,k) = (v_boost - nx*ui/u0L)/(1.d0 - v_boost*ui*nx/u0L/sol/sol)
              vy(i,j,k)=-ny*ui/u0L/gamv/(1.d0-v_boost*ui*nx/u0L/sol/sol)   
              vz(i,j,k)=-nz*ui/u0L/gamv/(1.d0-v_boost*ui*nx/u0L/sol/sol)
              u0(i,j,k)=(alph**2-psi4*(vx(i,j,k)*vx(i,j,k)+vy(i,j,k)*vy(i,j,k)+vz(i,j,k)*vz(i,j,k)))**-0.5
           endif

           !try something else
           if (sol_method .eq. 3) then
              !components of 3-vel in isotropic coords
              vx(i,j,k) = (v_boost - nx*ui/u0L)/(1.d0 - v_boost*ui*nx/u0L/sol/sol)
              vy(i,j,k)=-ny*ui/u0L/gamv/(1.d0-v_boost*ui*nx/u0L/sol/sol)
              vz(i,j,k)=-nz*ui/u0L/gamv/(1.d0-v_boost*ui*nx/u0L/sol/sol) 
              u0(i,j,k)=(alph**2-psi4*(vx(i,j,k)*vx(i,j,k)+vy(i,j,k)*vy(i,j,k)+vz(i,j,k)*vz(i,j,k)))**-0.5 
           endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           
           !use id from closest BH
           rho_b(i,j,k) = rho0

           !is this redundant?
           u0(i,j,k)=(alph**2-psi4*(vx(i,j,k)*vx(i,j,k)+vy(i,j,k)*vy(i,j,k)+vz(i,j,k)*vz(i,j,k)))**(-0.5)
          
        end do
     end do
  end do

  write(*,*) "outside loop"

  neos=1
  rho_tab(1)=1.0
  P_tab(1)=K_poly
  eps_tab(1)=P_tab(1)/rho_tab(1)/(gamma_th-1.0d0)
  !!write(*,*) "P_tab: ",P_tab
  !!write(*,*) "eps_tab: ",eps_tab

  do i=1,2
     k_tab(i)=K_poly
     gamma_tab(i)=gamma_th
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  if(CCTK_MyProc(CCTKGH)==0) then
     open(UNIT=15,FILE="run_info.txt",STATUS="REPLACE")
     write(15,*) "************************************************"
     write(15,*) "hydro info**************************************"
     write(15,*) "************************************************"
     write(15,*) "isotropic sonic radius: ",Ri_s
     write(15,*) "Mdot: ",Mdot
     write(15,*) "Gamma: ",gamma_th
!     write(15,*) "BH momentum: ",bh_px_plus
     write(15,*) "a_inf: ",sqrt((gamma_th-1.d0)*(sqrt(G29_rhs)-1)/sqrt(G29_rhs))
     write(15,*) "a_s: ",sqrt(a2_s)
     write(15,*) "************************************************"
     write(15,*) "grid info***************************************"
     write(15,*) "************************************************"
     write(15,*) "outer boundary: "
     write(15,*) "max resolution: "
     write(15,*) "min resolution: "
     write(15,*) "# refinement levels: "
     write(15,*) "************************************************"
     write(15,*) "shift info**************************************"
     write(15,*) "************************************************"
     write(15,*) "shift spatial_gauge: ",spatial_gauge
     write(15,*) "hbpunc_advect_enable: "
     write(15,*) "eta: "
     write(15,*) "************************************************"
     write(15,*) "lapse info**************************************"
     write(15,*) "************************************************"
     write(15,*) "slicing_type: "
     write(15,*) "opl_lapse_floor: "
     write(15,*) "opl_advect_enable: ",opl_advect_enable
     write(15,*) "*************************************************"
     close(15)
  endif
end subroutine statpunc_bondi_initialdata_driver

FUNCTION r_star_func_statpunc(x,M,Mdot,rho0_star,a2_inf,a2_star,gamma_th)
  IMPLICIT NONE
  real*8 :: x,M,Mdot,rho0_star,a2_inf,a2_star,gamma_th,r_star_func_statpunc
  real*8 :: pi
  pi=acos(-1.d0)

  r_star_func_statpunc = (1.d0-2.d0*M/x+(Mdot/(4.d0*pi*rho0_star*x**2))**2)&
       *(1.d0+a2_star/(gamma_th-1.d0-a2_star))**2&
       -(1.d0+a2_inf/(gamma_th-1.d0-a2_inf))**2
END FUNCTION r_star_func_statpunc

SUBROUTINE sonic_statpunc(M,a2_inf,gamma_th,a2_s,Rs_s,u2_s)
  IMPLICIT NONE
  real*8 :: M,a2_inf,gamma_th,a2_s,Rs_s,u2_s,a2_s_prev
  a2_s = a2_inf
  a2_s_prev=a2_s

  if (abs(gamma_th - 5.d0/3.d0) .gt. 1.d-4) then
     do
        a2_s = (2.d0*(gamma_th-1.d0)*a2_inf+&
             (1.d0-6.d0*(gamma_th-1.d0))*a2_s*a2_s+&
             3.d0*a2_s*a2_s*a2_s-a2_inf*a2_inf)/(2.d0*(gamma_th-1.d0)-&
             3.d0*(gamma_th-1.d0)**2)
        if (abs((a2_s-a2_s_prev)/a2_s_prev) .lt. 1.d-12) exit
        a2_s_prev = a2_s
     end do
  else
     do 
        a2_s = sqrt(a2_s*a2_s*a2_s+4.d0/9.d0*a2_inf-1.d0/3.d0*a2_inf*a2_inf)
        if (abs((a2_s-a2_s_prev)/a2_s_prev) .lt. 1.d-14) exit
        a2_s_prev = a2_s
     end do
  endif

  !inward 4-vel squared at sonic radius
  u2_s = a2_s/(1.d0+3.d0*a2_s)

  !sonic radius in schwarzchild coords
  Rs_s = (1.d0+3.d0*a2_s)/(2.d0*a2_s)*M 

END SUBROUTINE sonic_statpunc
