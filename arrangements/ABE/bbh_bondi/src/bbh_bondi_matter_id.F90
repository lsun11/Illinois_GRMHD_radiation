#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bbh_bondi_matter_id(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  INTERFACE
     FUNCTION get_rho0(rho0_derivs_s,Rs_s,a2_s,rho0_s,M,Rs,Mdot,K_poly,Gamma,G29_rhs)
       USE nrtype
       IMPLICIT NONE
       real(DP), INTENT(IN) :: Rs_s,a2_s,rho0_s,M,Rs,Mdot,K_poly,Gamma,G29_rhs,rho0_derivs_s
       real(DP)             :: get_rho0
     END FUNCTION get_rho0
  END INTERFACE
  INTERFACE
     SUBROUTINE bbh_sonic(M,a2_inf,gamma_th,a2_s,Rs_s,u2_s)
       IMPLICIT NONE
       real*8 :: M,a2_inf,gamma_th,a2_s,Rs_s,u2_s,a2_s_prev
     END SUBROUTINE bbh_sonic
  END INTERFACE
  INTERFACE
     FUNCTION bbh_r_star_func(x,M,Mdot,rho0_star,a2_inf,a2_star,gamma_th)
       IMPLICIT NONE
       real*8 :: x,M,Mdot,rho0_star,a2_inf,a2_star,gamma_th,bbh_r_star_func
     END FUNCTION bbh_r_star_func
  END INTERFACE

  integer, dimension(3) :: ext
  real*8 		:: dX,dY,dZ
  integer               :: i,j,k
  character             :: varname*30
  real*8 :: a2_s_prev
  real*8 :: Mdot,Rs_s
  real*8 :: xbh1,xbh2,ybh1,ybh2,zbh1,zbh2,u2_s,a2_s,G29_rhs,Rs_1,Ri_1,Rs_2,Ri_2,rho0_1,rho0_2,rho0_s,a2,Rs_a,a2_a,rho0_a,u2_a,rho0_derivi_a,rho0_derivs_s,Delta,us,nx,ny,nz,psi2,psi4,alph,ui,ui_1,ui_2,u0_1,u0_2,u0L,vx_1,vx_2,vy_1,vy_2,vz_1,vz_2,ux_1,uy_1,uz_1,ux_2,uy_2,uz_2,uxL,uyL,uzL,drhosurf,rr,gamv,sol
  real*8                :: K_poly_outer,G29_rhs_outer,a2_star_outer,a2_star_inner
  real*8                :: rho0_star,P_star,R_star,u2_star
  real*8                :: f_L,f_R,f_mid
  real*8                :: mb_o_me
  real*8                :: R_star_L,R_star_R,x1,x2
  real*8                :: BigMass
  integer               :: num,found,foundflag
  REAL*8, PARAMETER :: PI_D=3.141592653589793238462643383279502884197
  
  mb_o_me = 1.837152755d3 !ratio of baryon mass to electron mass           

  write(*,*) "*****************************"
  write(*,*) "inside bbh_bondi_matter_id"
  write(*,*) "*****************************"

  Delta = 0.1d0

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ext = cctk_lsh

  !commented out temporarily, but something is needed
  ! if (genID_cmdline_output_enable .eq. 1) BigMass=1.d0
  !temporary hack
  BigMass=1.0
  write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,*) "!!!!!!!!!!!WARNING WARNING !!!!!!!!!!!!!!"
  write(*,*) "!!!!!!!!!get the binary mass right!!!!!!!"
  write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  
!    I DON'T THINK THIS IS NEEDED FOR BBH
  !actually, I think it is needed
!!$  !Set initial BH position to zero.
  bh_posn_x = 0.D0
  bh_posn_y = 0.D0
  bh_posn_z = 0.D0


  foundflag = HorizonLocalCoordinateOrigin(1,xbh1,ybh1,zbh1)        
  foundflag = HorizonLocalCoordinateOrigin(2,xbh2,ybh2,zbh2)
  write(*,*) "xbh1: ",xbh1
  write(*,*) "xbh2: ",xbh2
  
   !Note that a capital Rs denotes schwarzchild radius, and Ri denotes isotropic radius

  !use eq G.28 to get K_poly
  K_poly_outer = a2_inf*rho0_inf**(1.d0-gamma_th_outer)/(gamma_th_outer - gamma_th_outer/(gamma_th_outer-1.d0)*a2_inf)

  !set this just in case
  n_poly = 1.d0/(gamma_th-1.d0)


  if (bbh_id_method .eq. 1) then
     write(*,*) "this method isn't currently supported"
     stop


!!$     !get info about sonic point
!!$     call bbh_sonic(bh_mass_plus,a2_inf,gamma_th_outer,a2_s,Rs_s,u2_s) !check this
!!$     
!!$     !this didn't work with checkpointing so I will do it by hand
!!$     !lum_outer_rad=4.d0*Rs_s
!!$
!!$     !now we can set density at sonic radius
!!$     rho0_s=(1.d0/K_poly_outer*a2_s/(gamma_th_outer - gamma_th_outer/(gamma_th_outer-1.d0)*a2_s))**(1.d0/(gamma_th_outer-1.d0))
!!$
!!$     !find derivative of rho0 with respect to Schwarzschild radius, evaluated at sonic radius
!!$     rho0_derivs_s = -2.d0 *rho0_s/Rs_s 
!!$
!!$     !now we can calculate Mdot
!!$     Mdot=4.d0*PI_D*rho0_s*sqrt(u2_s)*Rs_s*Rs_s
!!$
!!$     write(*,*) "a2_inf: ",a2_inf
!!$     write(*,*) "a2_s: ",a2_s
!!$     write(*,*) "u2_s: ",u2_s
!!$     write(*,*) "Rs_s: ",Rs_s
!!$     write(*,*) "Mdot: ",Mdot
!!$     write(*,*) "rho0_s: ",rho0_s
!!$
!!$
!!$     !this conversion between coords concerns me.  How can we do better?
!!$     !compute schwarszhild adjustment radius
!!$     Rs_a=Ri_a*(1.d0+bh_mass_plus/(2.d0*Ri_a))**2
!!$
!!$     !set constant on rhs of G.29
!!$     G29_rhs_outer = (1.d0+a2_inf/(gamma_th_outer-1.d0-a2_inf))**2
!!$
!!$     a2_star_outer = gamma_th_outer*(4.d0/3.d0)/mb_o_me/(1.d0+gamma_th_outer/(gamma_th_outer-1.d0)*(4.d0/3.d0)/mb_o_me)
!!$
!!$     if (const_gamma .eq. 1) then
!!$        if ((gamma_th_outer-gamma_th)/gamma_th .gt. 1.d-4) then
!!$           write(*,*) "you set const_gamma=1, so make sure you actually set gamma_th=gamma_th_outer"
!!$           stop
!!$        endif
!!$        G29_rhs = G29_rhs_outer
!!$        K_poly = K_poly_outer
!!$     else
!!$        
!!$        rho0_star = (1.d0/K_poly_outer*a2_star_outer/(gamma_th_outer - gamma_th_outer/(gamma_th_outer-1.d0)*a2_star_outer))**(1.d0/(gamma_th_outer-1.d0))
!!$        
!!$        P_star = K_poly_outer*rho0_star**gamma_th_outer
!!$        
!!$        R_star_L = 300.d0
!!$        R_star_R = 500.d0
!!$        
!!$        num=100
!!$        found=0
!!$        do while (found .eq. 0)
!!$           num=num*10
!!$           dx=(R_star_R-R_star_L)/num
!!$           
!!$           do i=1,num
!!$              x1=R_star_L+i*dx
!!$              x2=R_star_L+(i+1)*dx
!!$              if (bbh_r_star_func(x1,bh_mass_plus,Mdot,rho0_star,a2_inf,a2_star_outer,gamma_th_outer)*&
!!$                   bbh_r_star_func(x2,bh_mass_plus,Mdot,rho0_star,a2_inf,a2_star_outer,gamma_th_outer) .lt. 0) then
!!$                 found=1
!!$                 exit
!!$              endif
!!$           end do
!!$        end do
!!$        
!!$        R_star=x1
!!$        dx=x2-x1
!!$        
!!$        f_L=bbh_r_star_func(R_star,bh_mass_plus,Mdot,rho0_star,a2_inf,a2_star_outer,gamma_th_outer)
!!$        f_R=bbh_r_star_func(R_star+dx,bh_mass_plus,Mdot,rho0_star,a2_inf,a2_star_outer,gamma_th_outer)
!!$        
!!$        do 
!!$           dx=dx/2.d0
!!$           f_mid = bbh_r_star_func(R_star+dx,bh_mass_plus,Mdot,rho0_star,a2_inf,a2_star_outer,gamma_th_outer)
!!$           if (f_mid*f_L .gt. 0.d0) R_star=R_star+dx
!!$           if (dx .lt. 1.e-16) exit
!!$        end do
!!$        
!!$        u2_star =  (Mdot/(4.d0*PI_D*rho0_star*R_star*R_star))**2
!!$        
!!$        K_poly = P_star/rho0_star**gamma_th
!!$        a2_star_inner = gamma_th*K_poly*rho0_star**(gamma_th-1.d0)&
!!$             /(1.d0 + gamma_th/(gamma_th-1.d0)*K_poly*rho0_star**(gamma_th-1.d0))
!!$        
!!$        G29_rhs = (1.d0-2.d0*bh_mass_plus/R_star+u2_star)*(1.d0+a2_star_inner/(gamma_th-1.d0-a2_star_inner))**2
!!$        
!!$     endif
!!$
!!$     
!!$     
!!$     !find sound speed at adjustment radius
!!$     rho0_a=get_rho0(rho0_derivs_s,Rs_s,a2_s,rho0_s,bh_mass_plus,Rs_a,Mdot,K_poly,gamma_th,G29_rhs)
!!$
!!$     a2_a = gamma_th*K_poly*rho0_a**(gamma_th-1.d0)/(1.d0+gamma_th/(gamma_th-1.d0)*K_poly*rho0_a**(gamma_th-1.d0))
!!$
!!$     !find u^2 at adjustment radius
!!$     u2_a=(Mdot/(4.d0*PI_D*rho0_a*Rs_a*Rs_a))**2
!!$
!!$     !find derivative of rho0 with respect to isotropic radius, evaluated at adjustment radius
!!$     rho0_derivi_a=-(1.d0+bh_mass_plus/(2.d0*Ri_a))*(1.d0-bh_mass_plus/(2.d0*Ri_a))*(rho0_a/Rs_a)*(2.d0*u2_a-bh_mass_plus/Rs_a)/(u2_a-(1.d0-2.d0*bh_mass_plus/Rs_a+u2_a)*a2_a)
!!$
!!$
!!$     do k=1,ext(3)
!!$        do j=1,ext(2)
!!$           do i=1,ext(1)
!!$
!!$              !set psi2,psi4    
!!$              psi2 = exp(2.d0*phi(i,j,k))
!!$              psi4 = psi2*psi2
!!$              !set lapse
!!$              alph = lapm1(i,j,k)+1.0
!!$
!!$
!!$              !dist from BH 1
!!$              Ri_1=sqrt((X(i,j,k)-xbh1)*(X(i,j,k)-xbh1)+Y(i,j,k)*Y(i,j,k)+Z(i,j,k)*Z(i,j,k))
!!$              Rs_1=Ri_1*(1.d0+bh_mass_plus/(2.d0*Ri_1))**2
!!$              !dist from BH 2
!!$              Ri_2=sqrt((X(i,j,k)-xbh2)*(X(i,j,k)-xbh2)+Y(i,j,k)*Y(i,j,k)+Z(i,j,k)*Z(i,j,k))
!!$              Rs_2=Ri_2*(1.d0+bh_mass_plus/(2.d0*Ri_2))**2
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$              if (Ri_1 .gt. Ri_a) then
!!$                 rho0_1=get_rho0(rho0_derivs_s,Rs_s,a2_s,rho0_s,bh_mass_plus,Rs_1,Mdot,K_poly,gamma_th,G29_rhs)
!!$                 !radial component of 4-vel in schwarzschild coords
!!$                 us=Mdot/(4*PI_D*rho0_1*Rs_1*Rs_1)
!!$                 !radial component of 4-vel in isotropic coords
!!$                 ui_1=us/((1.d0-bh_mass_plus/(2.d0*Ri_1))*(1.d0+bh_mass_plus/(2.d0*Ri_1)))
!!$              endif
!!$
!!$              if (Ri_1 .le. Ri_a) then
!!$                 !inside adjustment radius, set rho0 to be constant value
!!$                 ! rho0=rho0_a
!!$                 rho0_1=rho0_a + rho0_derivi_a*(Ri_1*Ri_1-Ri_a*Ri_a)/(2.d0*Ri_a)
!!$                 a2=gamma_th*K_poly*rho0_1**(gamma_th-1.d0)/&
!!$                      (1.d0+gamma_th/(gamma_th-1.d0)*K_poly*rho0_1**(gamma_th-1.d0))
!!$                 !radial component of 4-vel in schwarzschild coords
!!$                 us=Mdot/(4.d0*PI_D*Rs_a*Rs_a*rho0_a)*Ri_1/Ri_a
!!$                 !radial component of 4-vel in isotropic coords
!!$                 ui_1=us/((1.d0-bh_mass_plus/(2.d0*Ri_a))*(1.d0+bh_mass_plus/(2.d0*Ri_a)))
!!$              endif
!!$
!!$	      !set u0,etc 
!!$              u0_1 = sqrt(1.D0+ui_1**2*psi4)/alph
!!$
!!$              gamv=1.0/sqrt(1.0-vxboost**2)
!!$              sol=alph/psi2
!!$
!!$              !components of radial unit vector
!!$              if (Ri_1 .gt. 0) then
!!$                 nx=(X(i,j,k)-xbh1)/Ri_1
!!$                 ny=Y(i,j,k)/Ri_1
!!$                 nz=Z(i,j,k)/Ri_1
!!$              else
!!$                 nx=0
!!$                 ny=0
!!$                 nz=0
!!$              endif
!!$
!!$	      !components of 4-vel in isotropic coords
!!$	      ux_1=-nx*ui_1
!!$	      uy_1=-ny*ui_1
!!$	      uz_1=-nz*ui_1
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$              if (Ri_2 .gt. Ri_a) then
!!$                 rho0_2=get_rho0(rho0_derivs_s,Rs_s,a2_s,rho0_s,bh_mass_plus,Rs_2,Mdot,K_poly,gamma_th,G29_rhs)
!!$                 !radial component of 4-vel in schwarzschild coords
!!$                 us=Mdot/(4*PI_D*rho0_2*Rs_2*Rs_2)
!!$                 !radial component of 4-vel in isotropic coords
!!$                 ui_2=us/((1.d0-bh_mass_plus/(2.d0*Ri_2))*(1.d0+bh_mass_plus/(2.d0*Ri_2)))
!!$              endif
!!$
!!$              if (Ri_2 .le. Ri_a) then
!!$                 !inside adjustment radius, set rho0 to be constant value
!!$                 ! rho0=rho0_a
!!$                 rho0_2=rho0_a + rho0_derivi_a*(Ri_2*Ri_2-Ri_a*Ri_a)/(2.d0*Ri_a)
!!$                 a2=gamma_th*K_poly*rho0_2**(gamma_th-1.d0)/&
!!$                      (1.d0+gamma_th/(gamma_th-1.d0)*K_poly*rho0_2**(gamma_th-1.d0))
!!$                 !radial component of 4-vel in schwarzschild coords
!!$                 us=Mdot/(4.d0*PI_D*Rs_a*Rs_a*rho0_a)*Ri_2/Ri_a
!!$                 !radial component of 4-vel in isotropic coords
!!$                 ui_2=us/((1.d0-bh_mass_plus/(2.d0*Ri_a))*(1.d0+bh_mass_plus/(2.d0*Ri_a)))
!!$              endif
!!$
!!$              !components of radial unit vector                                                                                                                                                                     
!!$              if (Ri_2 .gt. 0) then         
!!$                 nx=(X(i,j,k)-xbh2)/Ri_2
!!$                 ny=Y(i,j,k)/Ri_2
!!$                 nz=Z(i,j,k)/Ri_2
!!$              else
!!$                 nx=0
!!$                 ny=0
!!$                 nz=0
!!$              endif
!!$
!!$              !components of 4-vel in isotropic coords
!!$	      ux_2=-nx*ui_2
!!$	      uy_2=-ny*ui_2
!!$	      uz_2=-nz*ui_2
!!$
!!$
!!$              if (Ri_1 .lt. Ri_2) then
!!$	      	 rho_b(i,j,k) = rho0_1
!!$                 uxL=ux_1
!!$                 uyL=uy_1
!!$                 uzL=uz_1
!!$              else
!!$	         rho_b(i,j,k) = rho0_2
!!$                 uxL=ux_2
!!$                 uyL=uy_2
!!$                 uzL=uz_2
!!$              endif
!!$
!!$	      u0L = sqrt(1.D0+(uxL*uxL+uyL*uyL+uzL*uzL)**2*psi4)/alph
!!$              gamv=1.0/sqrt(1.0-vxboost**2)
!!$
!!$
!!$
!!$              uxL=uxL*gamv+u0L*gamv*vxboost
!!$
!!$              u0L = sqrt(1.D0+(uxL*uxL+uyL*uyL+uzL*uzL)*psi4)/alph
!!$
!!$              vx(i,j,k)=uxL/u0L
!!$              vy(i,j,k)=uyL/u0L
!!$              vz(i,j,k)=uzL/u0L
!!$
!!$
!!$              !is this redundant?
!!$              u0(i,j,k)=(alph**2-psi4*(vx(i,j,k)*vx(i,j,k)+vy(i,j,k)*vy(i,j,k)+vz(i,j,k)*vz(i,j,k)))**(-0.5)
!!$
!!$           end do
!!$        end do
!!$     end do
!!$
!!$     neos=1
!!$     rho_tab(1)=1.0
!!$     P_tab(1)=K_poly
!!$     eps_tab(1)=P_tab(1)/rho_tab(1)/(gamma_th-1.0d0)
!!$     write(*,*) "P_tab: ",P_tab
!!$     write(*,*) "eps_tab: ",eps_tab
!!$
!!$     do i=1,2
!!$        k_tab(i)=K_poly
!!$        gamma_tab(i)=gamma_th
!!$     enddo
!!$

  else if (bbh_id_method .eq. 3) then
     !get info about sonic point
     call bbh_sonic(BigMass,a2_inf,gamma_th_outer,a2_s,Rs_s,u2_s) !check this
     
     write(*,*) "BigMass: ",BigMass
     write(*,*) "a2_inf: ",a2_inf
     write(*,*) "gamma_th_outer: ",gamma_th_outer
     write(*,*) "a2_s: ",a2_s
     write(*,*) "Rs_s: ",Rs_s
     write(*,*) "u2_s: ",u2_s
     !this didn't work with checkpointing so I will do it by hand
     !lum_outer_rad=4.d0*Rs_s

     !now we can set density at sonic radius
     rho0_s=(1.d0/K_poly_outer*a2_s/(gamma_th_outer - gamma_th_outer/(gamma_th_outer-1.d0)*a2_s))**(1.d0/(gamma_th_outer-1.d0))

     !find derivative of rho0 with respect to Schwarzschild radius, evaluated at sonic radius
     rho0_derivs_s = -2.d0 *rho0_s/Rs_s 

     !now we can calculate Mdot
     Mdot=4.d0*PI_D*rho0_s*sqrt(u2_s)*Rs_s*Rs_s

     write(*,*) "a2_inf: ",a2_inf
     write(*,*) "a2_s: ",a2_s
     write(*,*) "u2_s: ",u2_s
     write(*,*) "Rs_s: ",Rs_s
     write(*,*) "Mdot: ",Mdot
     write(*,*) "BigMass: ",BigMass
     write(*,*) "rho0_s: ",rho0_s
     write(*,*) "Ri_a: ",Ri_a
     !compute schwarszhild adjustment radius
     Rs_a=Ri_a*(1.d0+BigMass/(2.d0*Ri_a))**2

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

        R_star_L = 300.d0
        R_star_R = 500.d0
        
        num=100
        found=0
        do while (found .eq. 0)
           num=num*10
           dx=(R_star_R-R_star_L)/num
           do i=1,num
              x1=R_star_L+i*dx
              x2=R_star_L+(i+1)*dx
              if (bbh_r_star_func(x1,BigMass,Mdot,rho0_star,a2_inf,a2_star_outer,gamma_th_outer)*&
                   bbh_r_star_func(x2,BigMass,Mdot,rho0_star,a2_inf,a2_star_outer,gamma_th_outer) .lt. 0) then
                 found=1
                 exit
              endif
            end do
         end do
         R_star=x1
         dx=x2-x1

        
        f_L=bbh_r_star_func(R_star,BigMass,Mdot,rho0_star,a2_inf,a2_star_outer,gamma_th_outer)
        f_R=bbh_r_star_func(R_star+dx,BigMass,Mdot,rho0_star,a2_inf,a2_star_outer,gamma_th_outer)
        do 
           dx=dx/2.d0
           f_mid = bbh_r_star_func(R_star+dx,BigMass,Mdot,rho0_star,a2_inf,a2_star_outer,gamma_th_outer)
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
        G29_rhs = (1.d0-2.d0*BigMass/R_star+u2_star)*(1.d0+a2_star_inner/(gamma_th-1.d0-a2_star_inner))**2

     endif


     write(*,*) "test1"
     !find sound speed at adjustment radius
     rho0_a=get_rho0(rho0_derivs_s,Rs_s,a2_s,rho0_s,BigMass,Rs_a,Mdot,K_poly,gamma_th,G29_rhs)
     write(*,*) "test2"
     a2_a = gamma_th*K_poly*rho0_a**(gamma_th-1.d0)/(1.d0+gamma_th/(gamma_th-1.d0)*K_poly*rho0_a**(gamma_th-1.d0))

     !find u^2 at adjustment radius
     u2_a=(Mdot/(4.d0*PI_D*rho0_a*Rs_a*Rs_a))**2

     !find derivative of rho0 with respect to isotropic radius, evaluated at adjustment radius
     rho0_derivi_a=-(1.d0+BigMass/(2.d0*Ri_a))*(1.d0-BigMass/(2.d0*Ri_a))*(rho0_a/Rs_a)*(2.d0*u2_a-BigMass/Rs_a)/(u2_a-(1.d0-2.d0*BigMass/Rs_a+u2_a)*a2_a)


     do k=1,ext(3)
        do j=1,ext(2)
           do i=1,ext(1)
              !set psi2,psi4                                                                                                                                                                                        
              psi2 = exp(2.d0*phi(i,j,k))
              psi4 = psi2*psi2
              !set lapse
              alph = lapm1(i,j,k)+1.0


              !dist from BH 1
              Ri_1=sqrt(X(i,j,k)*X(i,j,k)+Y(i,j,k)*Y(i,j,k)+Z(i,j,k)*Z(i,j,k))
              Rs_1=Ri_1*(1.d0+BigMass/(2.d0*Ri_1))**2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              if (Ri_1 .gt. Ri_a) then
                 !write(*,*) "rho0_derivs_s: ",rho0_derivs_s
                 !write(*,*) "Rs_s: ",Rs_s
                 !write(*,*) "a2_s: ",a2_s
                 !write(*,*) "rho0_s: ",rho0_s
                 !write(*,*) "BigMass: ",BigMass
                 !write(*,*) "Rs_1: ",Rs_1
                 !write(*,*) "Mdot: ",Mdot
                 !write(*,*) "K_poly: ",K_poly
                 !write(*,*) "gamma_th: ",gamma_th
                 !write(*,*) "G29_rhs: ",G29_rhs
                 rho0_1=get_rho0(rho0_derivs_s,Rs_s,a2_s,rho0_s,BigMass,Rs_1,Mdot,K_poly,gamma_th,G29_rhs)
                 !write(*,*) "rho0_1: ",rho0_1
                 !stop
                 !radial component of 4-vel in schwarzschild coords
                 us=Mdot/(4*PI_D*rho0_1*Rs_1*Rs_1)
                 !radial component of 4-vel in isotropic coords
                 ui_1=us/((1.d0-BigMass/(2.d0*Ri_1))*(1.d0+BigMass/(2.d0*Ri_1)))
                 !set u0 
                 u0_1 = sqrt(1.D0+ui_1**2*psi4)/alph
              endif

              if (Ri_1 .le. Ri_a) then
                 !inside adjustment radius, set rho0 to be constant value
                 rho0_1=rho0_a + rho0_derivi_a*(Ri_1*Ri_1-Ri_a*Ri_a)/(2.d0*Ri_a)
                 !rho0_1=rho0_a + rho0_derivi_a*(Ri_1-Ri_a)
                 a2=gamma_th*K_poly*rho0_1**(gamma_th-1.d0)/&
                      (1.d0+gamma_th/(gamma_th-1.d0)*K_poly*rho0_1**(gamma_th-1.d0))
                 !radial component of 4-vel in schwarzschild coords
                 us=Mdot/(4.d0*PI_D*Rs_a*Rs_a*rho0_a)*Ri_1/Ri_a
                 !radial component of 4-vel in isotropic coords
                 ui_1=us/((1.d0-BigMass/(2.d0*Ri_a))*(1.d0+BigMass/(2.d0*Ri_a)))
                 !set u0 
                 u0_1 = sqrt(1.D0+ui_1**2*psi4)/alph  
              endif

	      gamv=1.0/sqrt(1.0-vxboost**2)
              sol=alph/psi2

              !components of radial unit vector
              if (Ri_1 .gt. 0) then
                 nx=X(i,j,k)/Ri_1
                 ny=Y(i,j,k)/Ri_1
                 nz=Z(i,j,k)/Ri_1
              else
                 nx=0
                 ny=0
                 nz=0
              endif
              !components of 3-vel in isotropic coords
	      vx(i,j,k)=(sol*vxboost-nx*ui_1/u0_1)/(1.d0-vxboost*ui_1*nx/u0_1/sol)
              vy(i,j,k)=-ny*ui_1/u0_1/gamv/(1.d0-vxboost*ui_1*nx/u0_1/sol)
              vz(i,j,k)=-nz*ui_1/u0_1/gamv/(1.d0-vxboost*ui_1*nx/u0_1/sol)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              rho_b(i,j,k) = rho0_1
              !is this redundant?
              u0(i,j,k)=(alph**2-psi4*(vx(i,j,k)*vx(i,j,k)+vy(i,j,k)*vy(i,j,k)+vz(i,j,k)*vz(i,j,k)))**(-0.5)

              !this is wrong, but we are just initializing.  Brian only added this for debugging purposes
              !!u_0_check(i,j,k)=1.d0
           end do
        end do
     end do
     write(*,*) "test3"
     neos=1
     rho_tab(1)=1.0
     P_tab(1)=K_poly
     eps_tab(1)=P_tab(1)/rho_tab(1)/(gamma_th-1.0d0)
     write(*,*) "K_poly: ",K_poly
     do i=1,2
        k_tab(i)=K_poly
        gamma_tab(i)=gamma_th
     enddo

  else 
     write(*,*) "choose an appropriate initial data method"
     stop
  endif

end subroutine bbh_bondi_matter_id

FUNCTION bbh_r_star_func(x,M,Mdot,rho0_star,a2_inf,a2_star,gamma_th)
  IMPLICIT NONE
  real*8 :: x,M,Mdot,rho0_star,a2_inf,a2_star,gamma_th,bbh_r_star_func
  real*8 :: pi
  pi=acos(-1.d0)

  bbh_r_star_func = (1.d0-2.d0*M/x+(Mdot/(4.d0*pi*rho0_star*x**2))**2)&
       *(1.d0+a2_star/(gamma_th-1.d0-a2_star))**2&
       -(1.d0+a2_inf/(gamma_th-1.d0-a2_inf))**2
END FUNCTION bbh_r_star_func

SUBROUTINE bbh_sonic(M,a2_inf,gamma_th,a2_s,Rs_s,u2_s)
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
        if (abs((a2_s-a2_s_prev)/a2_s_prev) .lt. 1.d-14) exit
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

END SUBROUTINE bbh_sonic
