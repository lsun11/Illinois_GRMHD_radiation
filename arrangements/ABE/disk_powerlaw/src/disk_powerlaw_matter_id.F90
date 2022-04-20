#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine disk_powerlaw_matter_id(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  real*8 		:: dX,dY,dZ
  integer               :: i,j,k
  character             :: varname*30
  real*8                :: RoM,costheta,cos2theta,sin2theta,varpi2
  real*8                :: Delta,Sigma,A
  real*8		:: accuracy
  real*8                :: Xks,Yks,Zks
  real*8                :: zeta,zeta_in,xi,xi_in,log_enthalpy,enthalpy,epsilon,rho0,g_phiphi,g_tphi,vphi
!  real*8                :: BigMass
  REAL*8, PARAMETER :: PI_D=3.141592653589793238462643383279502884197
  real*8                     :: HALF, ONE, ZERO, TWO, FOUR
  real*8      :: lambda,lambda_inner,eta,ell,alpha_exp,f_inner,f,u_t,ut_bl,uphi_bl,ur_bl,uphi_ks,ut_ks,vphi_ks,k_const
  real*8      :: u_xL,u_yL,Psi4,lapse,ut,gupijuiuj
  real*8      :: uxL,uyL,uzL
  real*8      :: AA,BB,CC
  real*8      :: shift_x,shift_y,shift_z,beta2
  real*8      :: g_tt,g_tr,g_tth,g_tph,g_rr,g_rth,g_rph,g_thth,g_thph,g_phph
  real*8      :: gup_tt,gup_tph,gup_phph
  real*8      :: gup_tt_inner,gup_tph_inner,gup_phph_inner
  real*8      :: Sigma_inner,Delta_inner,A_inner,g_phph_inner,g_tt_inner,g_tph_inner,u_t_inner
  real*8      :: shiftph
  parameter(HALF = 0.5D0, ONE = 1.D0, ZERO = 0.D0, TWO = 2.D0, FOUR = 4.D0)

  
  write(*,*) "*****************************"
  write(*,*) "inside disk_powerlaw_matter_id"
  write(*,*) "*****************************"

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ext = cctk_lsh

  !commented out temporarily, but something is needed
  ! if (genID_cmdline_output_enable .eq. 1) BigMass=1.d0
  !temporary hack
!  BigMass=1.0
  write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,*) "!!!!!!!!!!!WARNING WARNING !!!!!!!!!!!!!!"
  write(*,*) "!!!!!!!!!get the binary mass right!!!!!!!"
  write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"


!  RoM_inner=15.d0
  Sigma_inner=RoM_inner*RoM_inner
  Delta_inner=RoM_inner*RoM_inner-2.d0*RoM_inner+sam_disk*sam_disk
  A_inner = (RoM_inner*RoM_inner + sam_disk*sam_disk)**2 - sam_disk*sam_disk*Delta_inner
  g_phph_inner = A_inner/Sigma_inner
  g_tt_inner = -(ONE-TWO*RoM_inner/Sigma_inner)
  g_tph_inner = -TWO*sam_disk*RoM_inner/Sigma_inner
  lambda_inner = sqrt(-g_phph_inner/g_tt_inner)

  alpha_exp = q/(q-TWO)
  eta = ell_inner/(lambda_inner**(TWO-q))
  k_const = eta**(-TWO/(q-TWO))
  
!  write(*,*) "eta: ",eta
!  gup_tt_inner = -(RoM_inner*RoM_inner + sam_disk*sam_disk + 2*RoM_inner*sam_disk*sam_disk/Sigma_inner)/Delta_inner
!  gup_tph_inner = -TWO*sam_disk*RoM_inner/Sigma_inner/Delta_inner
!  gup_phph_inner = (1-TWO*RoM_inner/Sigma_inner)/Delta_inner  

 ! u_t_inner = -1.d0/sqrt(-gup_tt_inner+2.d0*ell_inner*gup_tph_inner-ell_inner*ell_inner*gup_phph_inner)
  u_t_inner = -sqrt(Delta_inner/(g_tt_inner*ell_inner*ell_inner+TWO*g_tph_inner*ell_inner+g_phph_inner)) 
  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
           Xks = X(i,j,k)
           Yks = Y(i,j,k)
           Zks = Z(i,j,k)
           
           RoM = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)/BigMass
           
           costheta = (Z(i,j,k)/BigMass)/RoM
           cos2theta = costheta*costheta
           sin2theta = ONE - cos2theta
           
           Sigma = RoM*RoM + sam_disk*sam_disk*cos2theta
           Delta = RoM*RoM - 2.d0*RoM + sam_disk*sam_disk 
           A = (RoM*RoM + sam_disk*sam_disk)**2 - sam_disk*sam_disk*Delta*sin2theta
           
           !not needed
           !Xbl = Xks*cos(-sam_disk/Delta*RoM)-Yks*sin(-sam_disk/Delta*RoM)
           !Ybl = Xks*sin(-sam_disk/Delta*RoM)+Yks*cos(-sam_disk/Delta*RoM)
           !Zbl = Xks

           !"distance to the rotation axis" squared
           varpi2 = Delta*sin2theta
           
           if(sqrt(varpi2) .gt. 0.9d0*RoM_inner) then
              !compute metric in standard boyer-lindquist coordinates
              g_tt = -(ONE-TWO*RoM/Sigma)
              g_tr = 0.d0
              g_tth = 0.d0
              g_tph = -TWO*sam_disk*RoM*sin2theta/Sigma
              g_rr = Sigma/Delta
              g_rth = 0.d0
              g_rph = 0.d0
              g_thth = Sigma
              g_thph = 0.d0
              g_phph = A*sin2theta/Sigma
              
              gup_tt = -(RoM*RoM + sam_disk*sam_disk + 2*RoM*sam_disk*sam_disk/Sigma*sin2theta)/Delta
              !gup_tr = 0.d0
              !gup_tth = 0.d0
              gup_tph = -TWO*sam_disk*RoM/Sigma/Delta
              !gup_rr = Delta/Sigma
              !gup_rth = 0.d0
              !gup_rph = 0.d0
              !gup_thth = ONE/Sigma
              !gup_thph = 0.do
              gup_phph = (1-TWO*RoM/Sigma)/Delta/sin2theta
              
              lambda = sqrt(-g_phph/g_tt)
              
              ell = eta*lambda**(TWO-q)
!     Hassan : I add the following if to compute lambda accuratly for kerr black hole . current lambda in the code is for schwartzchild black hole. March 21
              if (abs(sam_disk)>1.d-12 .and. (1.eq.1)) then
                 accuracy = 1.d0
                 !write(*,*) "lambda" , lambda , "l" , ell , "accuracy" , accuracy
                 Do while (abs(accuracy) .gt. 1.0E-14)
                    accuracy = (lambda - sqrt(ell*(gup_tt - ell*gup_tph) / (gup_tph - ell*gup_phph)))/lambda
                    lambda = sqrt(ell*(gup_tt - ell*gup_tph) / (gup_tph - ell*gup_phph))
                    ell = eta*lambda**(TWO-q)
                    !write(*,*) "lambda" , lambda , "l" , ell , "accuracy" , accuracy 
                 end do
!                 write (*,*) "yeap this point is gone now I take care of the next point"
              end if

              !Hassan : end of change March 21.
                            

         
              u_t = -sqrt(varpi2/(g_tt*ell*ell+TWO*g_tph*ell+g_phph))
              ut_bl  = gup_tt*u_t - gup_tph * ell * u_t
             
              !f_inner = abs(1-k_const*ell_inner**(alpha_exp+ONE))**(ONE/(alpha_exp+ONE))
              !f       = abs(1-k_const*ell**(alpha_exp+ONE))**(ONE/(alpha_exp+ONE))
              
	      f_inner = abs(1-eta*eta*lambda_inner**(TWO-TWO*q))**(ONE/(alpha_exp+ONE))
	      f = abs(1-eta*eta*lambda**(TWO-TWO*q))**(ONE/(alpha_exp+ONE))


              enthalpy = u_t_inner*f_inner/u_t/f
              !if ((i.eq. 45) .and. (j.eq.41) .and. (k.eq.4)) then
              !   write(*,*) "X: ",X(i,j,k)
              !   write(*,*) "Y: ",Y(i,j,k) 
              !   write(*,*) "Z: ",Z(i,j,k) 
              !   write(*,*) "lambda: ",lambda
              !   write(*,*) "enthalpy: ",enthalpy
              !   write(*,*) "f: ",f
              !   write(*,*) "g_phph: ",g_phph
              !   write(*,*) "g_tt: ",g_tt
              !   write(*,*) "u_t: ",u_t
              !   stop
              !endif
              if(enthalpy.gt.ONE) then
                 epsilon = (enthalpy - ONE)/Gamma_th
                 rho0 = ( (Gamma_th-ONE)*epsilon/(K_poly) )**(ONE/(Gamma_th-ONE))
              else
                 enthalpy = ONE
                 epsilon = ZERO
                 rho0 = ZERO
              end if
              rho_b(i,j,k) = rho0
              uphi_bl =  (gup_tph-ell*gup_phph)*u_t
              ur_bl = 0.d0
              uphi_ks = sam_disk/Delta*ur_bl + uphi_bl
              ut_ks = TWO*RoM/Delta*ur_bl + ut_bl
              vphi_ks = uphi_ks / ut_ks
              
              uxL = -Yks/BigMass*uphi_ks
              uyL = Xks/BigMass*uphi_ks
              uzL = 0.d0
!!$              
!!$              vx(i,j,k) = uxL/ut_ks
!!$              vy(i,j,k) = uyL/ut_ks
!!$              vz(i,j,k) = uzL/ut_ks

             ! !recall u_z=0
             ! u_xL = -Yks/(Xks*Xks+Yks*Yks)*(-ell*u_t)
             ! u_yL =  Xks/(Xks*Xks+Yks*Yks)*(-ell*u_t)
              
              Psi4 = exp(4.d0*phi(i,j,k))
              lapse = 1.d0+lapm1(i,j,k)

              shift_x = Psi4*(gxx(i,j,k)*shiftx(i,j,k) + & 
                   gxy(i,j,k)*shifty(i,j,k) + gxz(i,j,k)*shiftz(i,j,k))
              shift_y = Psi4*(gxy(i,j,k)*shiftx(i,j,k) + & 
                   gyy(i,j,k)*shifty(i,j,k) + gyz(i,j,k)*shiftz(i,j,k)) 
              shift_z = Psi4*(gxz(i,j,k)*shiftx(i,j,k) + & 
                   gyz(i,j,k)*shifty(i,j,k) + gzz(i,j,k)*shiftz(i,j,k)) 
              
              beta2 = shiftx(i,j,k)*shift_x + shifty(i,j,k)*shift_y + shiftz(i,j,k)*shift_z

             ! gijuiuj = Psi4*(gxx(i,j,k)*uxL*uxL + 2.d0*gxy(i,j,k)*uxL*uyL + gyy(i,j,k)*uyL*uyL)

!!$              ut = sqrt(-1.d0 / (-lapse*lapse+beta2+&
!!$                   2.d0*(shift_x*vx(i,j,k)+shift_y*vy(i,j,k)+shift_z*vz(i,j,k)) + &
!!$                   Psi4*(gxx(i,j,k)*vx(i,j,k)*vx(i,j,k) + &
!!$                   gyy(i,j,k)*vy(i,j,k)*vy(i,j,k) + &
!!$                   gzz(i,j,k)*vz(i,j,k)*vz(i,j,k) + &
!!$                   2.d0*gxy(i,j,k)*vx(i,j,k)*vy(i,j,k) + &
!!$                   2.d0*gxz(i,j,k)*vx(i,j,k)*vz(i,j,k) + &
!!$                   2.d0*gyz(i,j,k)*vy(i,j,k)*vz(i,j,k))))
!!$              

              AA = -lapse*lapse+beta2
              BB = 2.d0 * (shift_x*uxL + &
                   shift_y*uyL + &
                   shift_z*uzL)
              CC =  1.d0+ Psi4*(gxx(i,j,k)*uxL*uxL + &
                   gyy(i,j,k)*uyL*uyL + &
                   gzz(i,j,k)*uzL*uzL + &
                   2.d0*gxy(i,j,k)*uxL*uyL + &
                   2.d0*gxz(i,j,k)*uxL*uzL + &
                   2.d0*gyz(i,j,k)*uyL*uzL)
              
              ut = (-BB - sqrt(BB*BB-4.d0*AA*CC))/(2.d0*AA)
              
              vx(i,j,k) = uxL/ut
              vy(i,j,k) = uyL/ut
              vz(i,j,k) = uzL/ut

             ! gupijuiuj = (gupxx(i,j,k)*u_xL*u_xL + 2.d0*gupxy(i,j,k)*u_xL*u_yL + gupyy(i,j,k)*u_yL*u_yL)/Psi4
            
             ! ut = 1.d0/lapse * sqrt(1.d0 + gupijuiuj)
              
        !      vx(i,j,k) = 1.d0/ut * (gupxx(i,j,k)*u_xL + gupxy(i,j,k)*u_yL)/Psi4 - shiftx(i,j,k)
        !      vy(i,j,k) = 1.d0/ut * (gupxy(i,j,k)*u_xL + gupyy(i,j,k)*u_yL)/Psi4 - shifty(i,j,k)
        !      vz(i,j,k) = 1.d0/ut * (gupxz(i,j,k)*u_xL + gupyz(i,j,k)*u_yL)/Psi4 - shiftz(i,j,k)

              u0(i,j,k) = ut
              
              !vx(i,j,k) = -Yks*vphi_ks
              !vy(i,j,k) = Xks*vphi_ks
              !vz(i,j,k) = ZERO
              !u0(i,j,k) = ut_ks
           else
              rho_b(i,j,k) = 1.d-20
              !vx(i,j,k) = ZERO
              !vy(i,j,k) = ZERO
              !vz(i,j,k) = ZERO
               vx(i,j,k) = -shiftx(i,j,k)
               vy(i,j,k) = -shifty(i,j,k)
               vz(i,j,k) = -shiftz(i,j,k)
              
              u0(i,j,k) = 1.d0/(lapm1(i,j,k)+1.d0)
           endif
        end do
     end do
  end do
 
  neos=1
  rho_tab(1)=1.d0
  P_tab(1)=K_poly
  eps_tab(1)=P_tab(1)/rho_tab(1)/(gamma_th-1.0d0)
  do i=1,2
     k_tab(i)=K_poly
     gamma_tab(i)=gamma_th
  enddo


end subroutine disk_powerlaw_matter_id

