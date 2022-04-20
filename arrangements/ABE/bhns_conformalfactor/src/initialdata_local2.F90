#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bhns_initialdata_local2(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext

  real*8 			           :: dT,dX,dY,dZ
  real*8 			           :: n1,n2,n3,mf

  integer :: vindex,handle,ierr
  real*8  :: rho_max,tau_max,valuetosetvz
  real*8  :: ONE,ZERO,W_L,U2, gxxL, gxyL, gxzL, gyyL, gyzL, gzzL,psiL,psiL4
  real*8  :: cn2, cn3, an1, an2, an3
  real*8  :: gupxxL, gupxyL, gupxzL, gupyyL, gupyzL, gupzzL, u0L,lapseL
  real*8  :: P_deplete, RADEQUAT, delta_bar, Omega_Frame
  real*8  :: rr, re2, ww, phia, pai=3.141592653589793d0
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(ONE = 1.D0, ZERO = 0.D0)

  real*8 :: Omega_r

  INTEGER :: i,j,k, ia

  character :: filename*30,c2*2,c1
  CCTK_REAL :: reduction_value

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ext = cctk_lsh



write(*,*) 'In bhns local2, rho_b_atm = ',rho_b_atm,dX

if ((use_new_bhns_initial_data/=5)) then


   ! the v's are really the u_i's when use_new_bhns_initial_data<2
   if ((use_new_bhns_initial_data==1).or. &
         (use_new_bhns_initial_data==3).or. & 
         (use_new_bhns_initial_data==4)) then
      u0=sqrt(1.0+psi**(-4)*(vx*vx+vy*vy+vz*vz))/(lapm1+1.0d0)
      !
   else if(use_new_bhns_initial_data==2) then
      ! The v's are the U's determined through u^a = W(n^a + U^a) when use_new_bhns_initial_data==2
      ! The Lorentz factor as measured by a normal observer is W = 1/Sqrt(1-g_ij U^i U^j)
      ! Knowing W we can find u^0 = W/lapse
      do k=1,ext(3)
         do j=1,ext(2)
            do i=1,ext(1)
               psiL=psi(i,j,k)
               psiL4=psiL*psiL*psiL*psiL
               ! Full 3-metric
               gxxL=psiL4*gxx(i,j,k)
               gxyL=psiL4*gxy(i,j,k)
               gxzL=psiL4*gxz(i,j,k)
               gyyL=psiL4*gyy(i,j,k)
               gyzL=psiL4*gyz(i,j,k)
               gzzL=psiL4*gzz(i,j,k)
               if(rho_b(i,j,k) .le. 0.d0) then
                  vx(i,j,k) = 0.d0
                  vy(i,j,k) = 0.d0
                  vz(i,j,k) = 0.d0
               end if
               U2 = gxxL*vx(i,j,k)*vx(i,j,k) + gyyL*vy(i,j,k)*vy(i,j,k) + gzzL*vz(i,j,k)*vz(i,j,k) + &
                     2.d0*(gxyL*vx(i,j,k)*vy(i,j,k) + gxzL*vx(i,j,k)*vz(i,j,k) + gyzL*vy(i,j,k)*vz(i,j,k) )
               if (abs(U2-1.d0).le.1.0d-3) then
                  write(*,*) "Stopping: can't have near luminal initial data: check your data"
                  stop
               end if
               W_L = 1.d0/sqrt(1.d0-U2)
               ! Compute u^0
               u0(i,j,k) =W_L/(lapm1(i,j,k)+1.0d0)
               ! And now redefine U^i to correspond to u_i = W * U_i
               vx(i,j,k) = W_L*( gxxL*vx(i,j,k)+gxyL*vy(i,j,k)+gxzL*vz(i,j,k) )
               vy(i,j,k) = W_L*( gxyL*vx(i,j,k)+gyyL*vy(i,j,k)+gyzL*vz(i,j,k) )
               vz(i,j,k) = W_L*( gxzL*vx(i,j,k)+gyzL*vy(i,j,k)+gzzL*vz(i,j,k) )
            enddo
         enddo
      enddo
   end if
   
        
  rho_b= max(rho_b,rho_b_atm)
  
  if (nperturb.ne.0) then
     write(6,*) "****************************** ADDING A PERTURBATION **************************************"
     write(6,*) "nperturb=", nperturb
     write(6,*) "ampl_perturb=", ampl_perturb
     write(6,*) "radi_perturb=", radi_perturb
     
     do k=1,ext(3)
        do j=1,ext(2)
           do i=1,ext(1)
              rr = x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k)
              ww = x(i,j,k)*x(i,j,k) - y(i,j,k)*y(i,j,k)
              re2= radi_perturb*radi_perturb
              if (nperturb==4 .or. nperturb==3)  rho_b(i,j,k) = rho_b(i,j,k)*(1.0d0 + ampl_perturb*ww/re2)
              
              if (rr .ne. 0) then
                 phia = atan2(y(i,j,k), x(i,j,k))
                 if (y(i,j,k)<0)    phia = 2.0d0*pai + phia
                 if (nperturb==1)  rho_b(i,j,k) = rho_b(i,j,k)*(1.0d0-ampl_perturb + ampl_perturb*cos(2.0d0*phia))
                 if (nperturb==2)  rho_b(i,j,k) = rho_b(i,j,k)*(1.0d0-ampl_perturb + ampl_perturb*cos(Omega_value*phia))
              end if
           end do
        end do
     end do
  end if


  !compute u0 before rho_star
  if(use_new_bhns_initial_data==6) then
      ! Remember now v^i are actually u_i  and we want v^i = u^i/u^0       

     print *, "Compute u0 for Idata =6"
     do k=1,ext(3)
        do j=1,ext(2)
           do i=1,ext(1)
              psiL  =   psi(i,j,k)
              psiL4 =   psiL*psiL*psiL*psiL
              lapseL=   lapm1(i,j,k) + 1.0d0
              
              ! Full inverse metric
              gupxxL = gupxx(i,j,k)/psiL4
              gupxyL = gupxy(i,j,k)/psiL4
              gupxzL = gupxz(i,j,k)/psiL4
              gupyyL = gupyy(i,j,k)/psiL4
              gupyzL = gupyz(i,j,k)/psiL4
              gupzzL = gupzz(i,j,k)/psiL4
              lapseL = lapm1(i,j,k)+1.d0
              
              u0(i,j,k) = sqrt(1.d0 + gupxxL*vx(i,j,k)*vx(i,j,k) + gupyyL*vy(i,j,k)*vy(i,j,k) + &
                    gupzzL*vz(i,j,k)*vz(i,j,k) + 2.d0*gupxyL*vx(i,j,k)*vy(i,j,k) +              &
                    2.d0*gupxzL*vx(i,j,k)*vz(i,j,k) + 2.d0*gupyzL*vz(i,j,k)*vy(i,j,k))/lapseL
           end do
        end do
     end do
  end if
 

  rho_star=rho_b*(lapm1+1.0d0)*u0*psi**6

  ! compute  pressure
  if(piecewise.eq.0) then  
     
     P=K_poly*rho_b**(gamma_th)
     h=1.0+gamma_th/(gamma_th-1.0)*P/rho_b
     
  else     
     
     print *, "^^^^^^^^^^^^^^^^^^^^^^ INSIDE initialdata_local2.F90 ^^^^^^^^^^^^^^^^^^^^^^^^ "
     write(6,'(2a18)') "Gamma", "Kappa"
     do ia=1, 11
        write(6,'(1p,2e18.10)') gamma_tab(ia), k_tab(ia)
     enddo
     write(6,'(3a18)') "rho0", "P", "epsilon"
     do ia=1, 10
        write(6,'(1p,3e18.10)') rho_tab(ia), P_tab(ia), eps_tab(ia)
     enddo
     
     
     where((rho_b - rho_tab(1))/rho_tab(1).le.1.0d-15)     ! atmosphere
        P    = k_tab(1)*rho_b**gamma_tab(1)
        
        h    = 1.0d0 + P/rho_b*gamma_tab(1)/(gamma_tab(1) - 1.0d0)
     elsewhere((rho_b - rho_tab(neos))/rho_tab(neos).ge.1.0d-15)   ! inner core
        P    = k_tab(neos+1)*rho_b**gamma_tab(neos+1)
        
        h    = 1.0d0+ eps_tab(neos) + P_tab(neos)/rho_tab(neos) &
              + (P/rho_b - P_tab(neos)/rho_tab(neos) )*gamma_tab(neos+1)/(gamma_tab(neos+1) - 1.0d0)
     endwhere
     do ia=1,neos-1
!!       where(((rho_b - rho_tab(ia))*(rho_b - rho_tab(ia+1))).le.0.0d0) 
        where ( ( (rho_b - rho_tab(ia))/rho_tab(ia).gt. 1.d-15).and.(  (rho_b - rho_tab(ia+1))/rho_tab(ia+1).le.1.d-15 ) )
           
           P    = k_tab(ia+1)*rho_b**gamma_tab(ia+1)
           
           h    = 1.0d0 + eps_tab(ia) + P_tab(ia)/rho_tab(ia) &
                 + (P/rho_b - P_tab(ia)/rho_tab(ia) )*gamma_tab(ia+1)/(gamma_tab(ia+1) - 1.0d0)
           
        endwhere
     end do
  end if
  
  w = rho_star*(lapm1 + 1.0d0)*u0
  tau = w*h-psi**(6)*P-rho_star
  

  st_x = rho_star*h*vx
  st_y = rho_star*h*vy
  st_z = rho_star*h*vz
  mhd_st_x=st_x
  mhd_st_y=st_y
  mhd_st_z=st_z
  
  ! and really the v's: v^i = u^i/u^0
  if ((use_new_bhns_initial_data==1).or. &
        (use_new_bhns_initial_data==3).or. & 
        (use_new_bhns_initial_data==4)) then
     vx=psi**(-4)*vx/u0-shiftx
     vy=psi**(-4)*vy/u0-shifty
     vz=psi**(-4)*vz/u0-shiftz
  else if ((use_new_bhns_initial_data==2).or.(use_new_bhns_initial_data==6)) then
     ! Remember now v^i are actually u_i
     ! and we want v^i = u^i/u^0
     do k=1,ext(3)
        do j=1,ext(2)
           do i=1,ext(1)
              psiL=psi(i,j,k)
              psiL4=psiL*psiL*psiL*psiL
              ! Full inverse metric
              gupxxL=gupxx(i,j,k)/psiL4
              gupxyL=gupxy(i,j,k)/psiL4
              gupxzL=gupxz(i,j,k)/psiL4
              gupyyL=gupyy(i,j,k)/psiL4
              gupyzL=gupyz(i,j,k)/psiL4
              gupzzL=gupzz(i,j,k)/psiL4
              u0L = u0(i,j,k)
              ! And now set v^i = g^{ij}u_i/u^0 - beta^i
              vx(i,j,k) = ( gupxxL*vx(i,j,k)+gupxyL*vy(i,j,k)+gupxzL*vz(i,j,k) )/u0L - shiftx(i,j,k)
              vy(i,j,k) = ( gupxyL*vx(i,j,k)+gupyyL*vy(i,j,k)+gupyzL*vz(i,j,k) )/u0L - shifty(i,j,k)
              vz(i,j,k) = ( gupxzL*vx(i,j,k)+gupyzL*vy(i,j,k)+gupzzL*vz(i,j,k) )/u0L - shiftz(i,j,k)
           enddo
        enddo
     enddo
  end if

  if(fisheye_enable .eq. 1) then
     call trans_phys_fish_tensor_flat(cctkGH,cctk_lsh,Symmetry, &
          x, y, z,  &
          PhysicalRadius,RadiusDerivative, &
          gxx,gxy,gxz,gyy,gyz,gzz)
     call trans_phys_fish_tensor_inv(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)
     call trans_phys_fish_tensor(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative, &
          Axx,Axy,Axz,Ayy,Ayz,Azz)
     call trans_phys_fish_phi(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative,phi)
     call trans_phys_fish_gamt_flat(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative,RadiusDerivative2, &
          Gammax,Gammay,Gammaz)
     call trans_phys_fish_matter(ext,X,Y,Z, &
          PhysicalRadius,RadiusDerivative, &
          rho_star,tau,st_x,st_y,st_z,w,vx,vy,vz,Symmetry)
  end if

  !still need to calculate basic matter variables

  rho = h*w*exp(-6.0*phi)-P
  Sx = st_x*exp(-6.0*phi)
  Sy = st_y*exp(-6.0*phi)
  Sz = st_z*exp(-6.0*phi)
  Sxx = st_x*st_x/w/h*exp(-6.0*phi) + P*gxx*exp(4.0*phi)
  Sxy = st_x*st_y/w/h*exp(-6.0*phi) + P*gxy*exp(4.0*phi)
  Sxz = st_x*st_z/w/h*exp(-6.0*phi) + P*gxz*exp(4.0*phi)
  Syy = st_y*st_y/w/h*exp(-6.0*phi) + P*gyy*exp(4.0*phi)
  Syz = st_y*st_z/w/h*exp(-6.0*phi) + P*gyz*exp(4.0*phi)
  Szz = st_z*st_z/w/h*exp(-6.0*phi) + P*gzz*exp(4.0*phi)
  S = exp(-4.0*phi)*(gupxx*Sxx+gupyy*Syy+gupzz*Szz+ &
       2.0*(gupxy*Sxy+gupxz*Sxz+gupyz*Syz))

  mhd_st_x=st_x
  mhd_st_y=st_y
  mhd_st_z=st_z

end if

! *********************
! ****   Cook ID   ****
! *********************
if (use_new_bhns_initial_data==5) then


   neos=1
   rho_tab(1)=1.0
   P_tab(1)=K_poly
   eps_tab(1)=K_poly/(gamma_th-1.0d0)
   do i=1,2
      k_tab(i)=K_poly
      gamma_tab(i)=gamma_th
   enddo
   
   n_poly = 1.0d0/(gamma_th-1.d0)
   delta_bar = 0.d0
   RADEQUAT = 1.d0


   ! Notice that the code has already read the rns ID   
   if (nperturb.ne.0) then
      write(6,*) "****************************** ADDING A PERTURBATION **************************************"
      write(6,*) "nperturb=", nperturb
      write(6,*) "ampl_perturb=", ampl_perturb
      write(6,*) "radi_perturb=", radi_perturb
      
      do k=1,ext(3)
         do j=1,ext(2)
            do i=1,ext(1)
               rr = x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k)
               ww = x(i,j,k)*x(i,j,k) - y(i,j,k)*y(i,j,k)
               re2= radi_perturb*radi_perturb
               if (nperturb==4 .or. nperturb==3)  rho_b(i,j,k) = rho_b(i,j,k)*(1.0d0 + ampl_perturb*ww/re2)
               
               if (rr .ne. 0) then
                  phia = atan2(y(i,j,k), x(i,j,k))
                  if (y(i,j,k)<0)    phia = 2.0d0*pai + phia
                  
                  !            if ( dabs(x(i,j,k))<=2.0 .and. dabs(y(i,j,k))<=2.0 .and. dabs(z(i,j,k))<=1.0   )  write(6,'(a6,1p,10e20.10)')  "x,y,z=", x(i,j,k), y(i,j,k), z(i,j,k), phia, ww, re2
                  
                  if (nperturb==1)  rho_b(i,j,k) = rho_b(i,j,k)*(1.0d0-ampl_perturb + ampl_perturb*cos(2.0d0*phia))
                  if (nperturb==2)  rho_b(i,j,k) = rho_b(i,j,k)*(1.0d0-ampl_perturb + ampl_perturb*cos(Omega_value*phia))
               end if
            end do
         end do
      end do
      P=K_poly*rho_b**(gamma_th)
   end if

   call compute_magnetar_hybrid_bhns(ext,P_deplete_bhns, RADEQUAT, delta_bar, &
        X,Y,Z, &
        neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th, &
        lapm1,phi,shiftx,shifty, &
        gxx,gxy,gxz,gyy,gyz,gzz, &
        rho,S,Sx,Sy,Sz, &
        Sxx,Sxy,Sxz,Syy,Syz,Szz, &
        rho_star,tau,st_x,st_y,st_z, &
        P,w,vx,vy,vz,rho_b, &
        u0,h,rho_b_atm,Omega_Frame,PhysicalRadius,eps_flag,K_poly,n_poly,RESET_RHO_B_ATM,rhob_fac2)
 


   if(reset_shift_lapse==1) then
       shiftx=0.
      shifty=0.
      shiftz=0.
   endif

   mhd_st_x=st_x
   mhd_st_y=st_y
   mhd_st_z=st_z

  if(fisheye_enable .eq. 1) then
     call trans_phys_fish_tensor_flat(cctkGH,cctk_lsh,Symmetry, &
          x, y, z,  &
          PhysicalRadius,RadiusDerivative, &
          gxx,gxy,gxz,gyy,gyz,gzz)
     call trans_phys_fish_tensor_inv(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)
     call trans_phys_fish_tensor(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative, &
          Axx,Axy,Axz,Ayy,Ayz,Azz)
     call trans_phys_fish_phi(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative,phi)
     call trans_phys_fish_gamt_flat(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative,RadiusDerivative2, &
          Gammax,Gammay,Gammaz)
     call trans_phys_fish_matter(ext,X,Y,Z, &
          PhysicalRadius,RadiusDerivative, &
          rho_star,tau,st_x,st_y,st_z,w,vx,vy,vz,Symmetry)
  end if

end if

if(1==0) then
   if ((use_new_bhns_initial_data.ne.5).and.((nperturb==5 .or. nperturb==3))) then
      write(6,*) "****************************** ADDING A PERTURBATION IN VELOCITY ***********************"
      write(6,*) "nperturb=", nperturb
      write(6,*) "lambda_perturb=", lambda_perturb
      write(6,*) "a2oa1_perturb=", a2oa1_perturb
      
      do k=1,ext(3)
         do j=1,ext(2)
            do i=1,ext(1)
               vx(i,j,k) = vx(i,j,k) + y(i,j,k)*lambda_perturb/a2oa1_perturb
               vy(i,j,k) = vy(i,j,k) - x(i,j,k)*lambda_perturb*a2oa1_perturb
            end do
         end do
      end do
   end if
end if

  if(Symmetry==AXISYM) then
     write(*,*) "I PITY DA FOO TRIES TO USE AXISYMMETRY FOR ORBITING BLACK HOLES!"
     stop
  end if

  !Following lines are (empirically speaking) ESSENTIAL for excision runs:
  call Derivs(ext,X,Y,Z,dX,dY,dZ,phi,phix,phiy,phiz,Symmetry) 
  call Derivs(ext,X,Y,Z,dX,dY,dZ,lapm1,lapsex,lapsey,lapsez,Symmetry) 
  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,phi,phix,phiy,phiz)
  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,lapm1,lapsex,lapsey,lapsez)

  !call CartSymGN(ierr,cctkGH,'BSSN::BSSN_vars')
  !  call CCTK_SyncGroup(ierr,cctkGH,'BSSN::BSSN_vars') 

  ! Initialize t_last, M0dot_last and int_M0dot: variables necessary to compute the
  !  time integrated M0 flux across BH
  t_last = 0.d0
  M0dot_last = 0.d0
  int_M0dot = 0.d0

  ! Initialize the horizon position and radius
  xh_last = xh0
  yh_last = yh0
  zh_last = zh0
  ah_radii_last = r0_ah

  if(subtract_off_Omega_r==1) then
     vx = vx + Omega_value*y
     vy = vy - Omega_value*x
  end if

end subroutine bhns_initialdata_local2






