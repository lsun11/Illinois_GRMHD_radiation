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
  real*8  :: betaxL,betayL,betazL
  real*8  :: P_deplete, RADEQUAT, delta_bar, Omega_Frame, P_radl
  real*8  :: P_rad0xl, P_rad0yl, P_rad0zl, P_rad00l  
  real*8  :: F_rad_xl, F_rad_yl, F_rad_zl, F_rad_0l, temp_rad, temp_rad1
  real*8  :: beta2, udotbeta, g_00l, u_0l
  real*8  :: rr, re2, ww, phia, pai=3.141592653589793d0
  real*8  :: u_x, u_y, u_z, uxl, uyl, uzl, v_xl, v_yl, v_zl
  real*8  :: shift_x, shift_y, shift_z
  real*8  :: Fasq, zeta_temp, zeta_cut, zeta, chil, eps_l
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(ONE = 1.D0, ZERO = 0.D0)





  real*8 :: Omega_r


  real*8 :: xNS1, yNS1, xNS2, yNS2




  INTEGER :: i,j,k, ia

  character :: filename*30,c2*2,c1
  CCTK_REAL :: reduction_value

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ext = cctk_lsh

! Calculate CoM for initialize temperature, add a small value of the denominator.



 if(CCTK_ITERATION==0) then
      
        ! compute the position of the NSs (iter = 0)
        xNS1 = initial_ns_coord_x
        yNS1 = initial_ns_coord_y

        if((use_new_bhns_initial_data.eq.3).or.(use_new_bhns_initial_data.eq.4)) then
        xNS2 = initial_ns2_coord_x
        yNS2 = initial_ns2_coord_y          
        end if

 else
        if (Box1denom_VolInt == 0.0) then
   	Box1denom_VolInt = 1.0e-10
	end if
	if (Box1denom_VolInt1 == 0.0) then
   	Box1denom_VolInt2 = 1.0e-10
	end if
	if (Box1denom_VolInt1 == 0.0) then
   	Box1denom_VolInt2 = 1.0e-10
	end if

        ! compute the position of the NSs (iter > 0)
     if((use_new_bhns_initial_data.eq.3).or.(use_new_bhns_initial_data.eq.4)) then
        xNS1 = Box1X_VolInt1/Box1denom_VolInt1
        yNS1 = Box1Y_VolInt1/Box1denom_VolInt1

        xNS2 = Box1X_VolInt2/Box1denom_VolInt2
        yNS2 = Box1Y_VolInt2/Box1denom_VolInt2
     else
        ! compute the position of the NS
	xNS1 = Box1X_VolInt/Box1denom_VolInt
	yNS1 = Box1Y_VolInt/Box1denom_VolInt
     end if
 
 end if



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
  if((use_new_bhns_initial_data==6).or. & 
         (use_new_bhns_initial_data==8)) then
      ! Remember now v^i are actually u_i  and we want v^i = u^i/u^0       

     print *, "Compute u0 for Idata =6 or Idata =8"
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
 
  if(ergo_star.eq.0) then
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
        where ( ( (rho_b - rho_tab(ia))/rho_tab(ia) .gt. 1.d-15).and.(  (rho_b - rho_tab(ia+1))/rho_tab(ia+1).le.1.d-15 ) )
           
           P    = k_tab(ia+1)*rho_b**gamma_tab(ia+1)
          
           h    = 1.0d0 + eps_tab(ia) + P_tab(ia)/rho_tab(ia) &
	        + (P/rho_b - P_tab(ia)/rho_tab(ia) )*gamma_tab(ia+1)/(gamma_tab(ia+1) - 1.0d0)
           
        endwhere
     end do
  
  else

     print *, "^^^^^^^^^^^^^^^^^^^^^^ INSIDE initialdata_local2.F90- Ergo Star ^^^^^^^^^^^^^^^^^^^^^^^^ "
     print *, "neos = ", neos 
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
      P = (ergo_sigma*((1.0d0+eps_tab(neos)+P_tab(neos)/rho_tab(neos))/rho_tab(neos)**ergo_sigma)*rho_b**(ergo_sigma+1) + P_tab(neos) - ergo_sigma*(rho_tab(neos)*(eps_tab(neos)+1)))/(ergo_sigma+1)
      h = ((1.0d0+eps_tab(neos)+P_tab(neos)/rho_tab(neos))/rho_tab(neos)**ergo_sigma)*rho_b**ergo_sigma
     endwhere

     do ia=1,neos-1

	 where ( ( (rho_b - rho_tab(ia))/rho_tab(ia).gt. 1.d-15).and.(  (rho_b - rho_tab(ia+1))/rho_tab(ia+1).le.1.d-15 ) )
           P    = k_tab(ia+1)*rho_b**gamma_tab(ia+1)
           h    = 1.0d0 + eps_tab(ia) + P_tab(ia)/rho_tab(ia) &
                 + (P/rho_b - P_tab(ia)/rho_tab(ia) )*gamma_tab(ia+1)/(gamma_tab(ia+1) - 1.0d0)
        endwhere

     end do     
  end if
  
 
  if (compute_microphysics .eq. 1) then 
       P_cld = P
       eps_cld = h - 1.d0 - P/rho_b

!!!! Experiment: Add a little addiional eps (and pressure) as the thermal part to make the temperature finite
	     do k=1,ext(3)
		 do j=1,ext(2)
           	     do i=1,ext(1)
		     eps_thermal(i,j,k) = (eps_thermal_bhns -1.0) * eps_cld(i,j,k) !!! This should be zero since eps_thermal_bhns = 1.0
		     eps_tot(i,j,k) = eps_thermal(i,j,k) + eps_cld(i,j,k)

            	     call compute_P_T_microphys(P(i,j,k), T_fluid(i,j,k), P_cld(i,j,k), eps_tot(i,j,k), eps_cld(i,j,k), rho_b(i,j,k))
	    	    end do
	         end do
	    end do	    

!	   write(*,*) "Inside initialdata_local2.F90, P", P
!	   write(*,*) "Inside initialdata_local2.F90, T_fluid", T_fluid
!	   write(*,*) "Inside initialdata_local2.F90, P_cld", P_cld
!	   write(*,*) "Inside initialdata_local2.F90, eps_tot", eps_tot
!	   write(*,*) "Inside initialdata_local2.F90, eps_cld", eps_cld	   

!	   P = P + (gamma_th-1.d0) * rho_b * eps_thermal
	   h = 1.d0 + P/rho_b + eps_tot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end if

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
  else if ((use_new_bhns_initial_data==2).or.(use_new_bhns_initial_data==6).or. & 
        (use_new_bhns_initial_data==8)) then
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

   
   ! single Polytrope
   if(piecewise.eq.0) then  
   print *, "assuming NOT piecewise EOS for Cook ID"   
   
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
      
      print *, "Now call:compute_magnetar_hybrid_bhns"
      
      call compute_magnetar_hybrid_bhns(ext,P_deplete_bhns, RADEQUAT, delta_bar, &
            X,Y,Z, &
            neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th, &
            lapm1,phi,shiftx,shifty, &
            gxx,gxy,gxz,gyy,gyz,gzz, &
            rho,S,Sx,Sy,Sz, &
            Sxx,Sxy,Sxz,Syy,Syz,Szz, &
            rho_star,tau,st_x,st_y,st_z, &
            P,w,vx,vy,vz,rho_b, &
            u0,h,rho_b_atm, &
            Omega_Frame,PhysicalRadius,eps_flag,K_poly,n_poly,RESET_RHO_B_ATM,rhob_fac2,enable_OS_collapse)
      
   else			
      ! piecewise polytrope
      print *, "assuming piecewise EOS for Cook ID"
      gamma_th= 1.66d0

      call compute_magnetar_hybrid_bhnsII(ext,P_deplete_bhns,X,Y,Z, &
            neos,ergo_star, ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th, &
            lapm1,phi,shiftx,shifty, &
            gxx,gxy,gxz,gyy,gyz,gzz, &
            rho,S,Sx,Sy,Sz, &
            Sxx,Sxy,Sxz,Syy,Syz,Szz, &
            rho_star,tau,st_x,st_y,st_z, &
            P,w,vx,vy,vz,&
	    P_cld, eps_cld, eps_thermal,&
	    initial_ns_coord_x,initial_ns_coord_y,&
	    rho_b,u0,h,rho_b_atm,Omega_Frame,PhysicalRadius,eps_flag,K_poly,n_poly,RESET_RHO_B_ATM,rhob_fac2, eps_thermal_bhns)


         do k=1,ext(3)
            do j=1,ext(2)
               do i=1,ext(1)

               psiL=psi(i,j,k)
               psiL4=psiL*psiL*psiL*psiL

               lapseL=   lapm1(i,j,k) + 1.0d0
               betaxL  =  shiftx(i,j,k)
               betayL  =  shifty(i,j,k)
               betazL  =  shiftz(i,j,k)

               ! Full 3-metric
               gxxL=psiL4*gxx(i,j,k)
               gxyL=psiL4*gxy(i,j,k)
               gxzL=psiL4*gxz(i,j,k)
               gyyL=psiL4*gyy(i,j,k)
               gyzL=psiL4*gyz(i,j,k)
               gzzL=psiL4*gzz(i,j,k)

               ! g00
               temp_g00(i,j,k) = psiL4*(gxxL*betaxL*betaxL  +         &
                                 gyyL*betayL*betayL  +         &
                                 gzzL*betazL*betazL  +         &
                          2.0d0*(gxyL*betaxL*betayL  +         &
                                 gxzL*betaxL*betazL  +         &
                                 gyzL*betayL*betazL)) - (lapseL*lapseL)

               end do
            end do
         end do


      end if



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


     call CCTK_ReductionHandle(handle,"maximum") 
     call CCTK_VarIndex(vindex,"mhd_evolve::rho_b")	
     
     if (handle .gt. 0) then
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,reduction_value,1,vindex)
     if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH).eq.0) then
     print *,"INSIDE LOCAL2: Maximum value of rho_b is ",reduction_value    
     end if 
     else 
     call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
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


if (compute_microphysics .eq. 1) then

  do k=1,ext(3)
    do j=1,ext(2)
      do i=1,ext(1)

         Y_e(i,j,k) = 0.0d0
         rhoYe(i,j,k) = Y_e(i,j,k) * rho_star(i,j,k)
         
	 eps_l = eps_thermal(i,j,k) + eps_cld(i,j,k)
!	 call compute_P_T_microphys(P(i,j,k), T_fluid(i,j,k), P_cld(i,j,k), eps_l, eps_cld(i,j,k), rho_b(i,j,k))

      end do
   end do
end do

end if


   if (rad_evolve_enable .eq. 1  .and. iteration_to_insert_rad .eq. 0) then

write(*,*) "Inside BHNS initialdata_local2.F90, start imposed radiation!!!!"

  do k=1,ext(3)
    do j=1,ext(2)
      do i=1,ext(1)     
           
           psiL=psi(i,j,k)
           psiL4=psiL*psiL*psiL*psiL

           shift_x = psiL4 *(gxx(i,j,k)*shiftx(i,j,k) + gxy(i,j,k)*shifty(i,j,k) + gxz(i,j,k)*shiftz(i,j,k))
           shift_y = psiL4 *(gxy(i,j,k)*shiftx(i,j,k) + gyy(i,j,k)*shifty(i,j,k) + gyz(i,j,k)*shiftz(i,j,k))
           shift_z = psiL4 *(gxz(i,j,k)*shiftx(i,j,k) + gyz(i,j,k)*shifty(i,j,k) + gzz(i,j,k)*shiftz(i,j,k))
             
           v_xl = psiL4 * (gxx(i,j,k)*vx(i,j,k) + gxy(i,j,k)*vy(i,j,k) + gxz(i,j,k)*vz(i,j,k))
           v_yl = psiL4 * (gxy(i,j,k)*vx(i,j,k) + gyy(i,j,k)*vy(i,j,k) + gyz(i,j,k)*vz(i,j,k))
           v_zl = psiL4 * (gxz(i,j,k)*vx(i,j,k) + gyz(i,j,k)*vy(i,j,k) + gzz(i,j,k)*vz(i,j,k))

           beta2 = shiftx(i,j,k)*shift_x + shifty(i,j,k)*shift_y + shiftz(i,j,k)*shift_z
           udotbeta = u0(i,j,k)*(vx(i,j,k)*shift_x + vy(i,j,k)*shift_y + vz(i,j,k)*shift_z)
           g_00l =beta2-(lapm1(i,j,k)+1.0)*(lapm1(i,j,k)+1.0)
           u_0l = g_00l*u0(i,j,k) + udotbeta
             
           if(u_0l .eq. 0.0 .or. isnan(u_0l)) then
              write(*,*) "Inside BHNS initial2, u_0l is zero or nan!!!!, u_0l, g_00l, u0(i,j,k), udotbeta are", u_0l, g_00l, u0(i,j,k), udotbeta
           end if

           uxl = u0(i,j,k)*vx(i,j,k)
           uyl = u0(i,j,k)*vy(i,j,k)
           uzl = u0(i,j,k)*vz(i,j,k)


           E_rad(i,j,k) = rho_b(i,j,k) * Erad_over_rho
           F_radx(i,j,k) = 0.0
           F_rady(i,j,k) = 0.0
           F_radz(i,j,k) = 0.0
             
           if (E_rad(i,j,k) .lt. Erad_atm_cut) then
              E_rad(i,j,k) = Erad_atm_cut
           end if

           u_x = u0(i,j,k)*psiL4*(gxx(i,j,k)*(vx(i,j,k) + shiftx(i,j,k)) + gxy(i,j,k)*(vy(i,j,k) + shifty(i,j,k)) + gxz(i,j,k)*(vz(i,j,k) + shiftz(i,j,k)))
           u_y = u0(i,j,k)*psiL4*(gxy(i,j,k)*(vx(i,j,k) + shiftx(i,j,k)) + gyy(i,j,k)*(vy(i,j,k) + shifty(i,j,k)) + gyz(i,j,k)*(vz(i,j,k) + shiftz(i,j,k)))
           u_z = u0(i,j,k)*psiL4*(gxz(i,j,k)*(vx(i,j,k) + shiftx(i,j,k)) + gyz(i,j,k)*(vy(i,j,k) + shifty(i,j,k)) + gzz(i,j,k)*(vz(i,j,k) + shiftz(i,j,k)))

           F_rad0(i,j,k) = - (F_radx(i,j,k)*u_x + F_rady(i,j,k)*u_y + F_radz(i,j,k)*u_z)/u_0l

           if (isnan(F_rad0(i,j,k))) then
              write(*,*) "Inside BHNS initial2, F_rad0(i,j,k) is nan!!!!, F_radx(i,j,k), u_x, F_rady(i,j,k), u_y, F_radz(i,j,k), u_z, u_0l are",  F_radx(i,j,k), u_x, F_rady(i,j,k), u_y, F_radz(i,j,k), u_z, u_0l
           end if


           F_rad_xl = psiL4 * (gxx(i,j,k) * F_radx(i,j,k) + gxy(i,j,k) * F_rady(i,j,k) + gxz(i,j,k) * F_radz(i,j,k)) + shift_x* F_rad0(i,j,k)
           F_rad_yl = psiL4 * (gxy(i,j,k) * F_radx(i,j,k) + gyy(i,j,k) * F_rady(i,j,k) + gyz(i,j,k) * F_radz(i,j,k)) + shift_y* F_rad0(i,j,k)
           F_rad_zl = psiL4 * (gxz(i,j,k) * F_radx(i,j,k) + gyz(i,j,k) * F_rady(i,j,k) + gzz(i,j,k) * F_radz(i,j,k)) + shift_z* F_rad0(i,j,k)
           F_rad_0l = - (F_rad_xl*uxl + F_rad_yl*uyl + F_rad_zl*uzl)/u0(i,j,k)
                       
           if (rad_closure_scheme.eq.0) then
              P_radl = E_rad(i,j,k)/3.d0
              
              temp_rad = (lapm1(i,j,k)+1.0)*u0l
              temp_rad1 = temp_rad*temp_rad*(E_rad(i,j,k)+P_radl) - P_radl + 2.d0*(lapm1(i,j,k)+1.0)*u0(i,j,k)*F_rad0(i,j,k)
              
              tau_rad(i,j,k) = (lapm1(i,j,k)+1.0)*(lapm1(i,j,k)+1.0)*psiL**6*(E_rad(i,j,k)*u0(i,j,k)*u0(i,j,k)+2.0*F_rad0(i,j,k)*u0(i,j,k)+P_radl*u0(i,j,k)*u0(i,j,k))-psiL**6*P_radl
              S_rad_x(i,j,k) = (lapm1(i,j,k)+1.0)*psiL**6*((E_rad(i,j,k)+P_radl)*u0(i,j,k)*u_x + F_rad0(i,j,k)*u_x + F_rad_xl * u0(i,j,k))
              S_rad_y(i,j,k) = (lapm1(i,j,k)+1.0)*psiL**6*((E_rad(i,j,k)+P_radl)*u0(i,j,k)*u_y + F_rad0(i,j,k)*u_y + F_rad_yl * u0(i,j,k))
              S_rad_z(i,j,k) = (lapm1(i,j,k)+1.0)*psiL**6*((E_rad(i,j,k)+P_radl)*u0(i,j,k)*u_z + F_rad0(i,j,k)*u_z + F_rad_zl * u0(i,j,k))
              
              rho(i,j,k) = rho(i,j,k) + temp_rad1
              Sx(i,j,k) = Sx(i,j,k) + temp_rad*( ( (E_rad(i,j,k) + P_radl) * u0(i,j,k) + F_rad0(i,j,k)) * (shift_x + v_xl) + F_rad_xl)
              Sy(i,j,k) = Sy(i,j,k) + temp_rad*( ( (E_rad(i,j,k) + P_radl) * u0(i,j,k) + F_rad0(i,j,k)) * (shift_y + v_yl) + F_rad_yl)
              Sz(i,j,k) = Sz(i,j,k) + temp_rad*( ( (E_rad(i,j,k) + P_radl) * u0(i,j,k) + F_rad0(i,j,k)) * (shift_z + v_zl) + F_rad_zl)
              Sxx(i,j,k) = Sxx(i,j,k) + (E_rad(i,j,k)+P_radl)*(u0(i,j,k)*(shift_x + v_xl))**2 + 2.0*F_rad_xl*u0(i,j,k)*(shift_x + v_xl) + psiL4 * P_radl * gxx(i,j,k)
              Syy(i,j,k) = Syy(i,j,k) + (E_rad(i,j,k)+P_radl)*(u0(i,j,k)*(shift_y + v_yl))**2 + 2.0*F_rad_yl*u0(i,j,k)*(shift_y + v_yl) + psiL4 * P_radl * gyy(i,j,k)
              Szz(i,j,k) = Szz(i,j,k) + (E_rad(i,j,k)+P_radl)*(u0(i,j,k)*(shift_z + v_zl))**2 + 2.0*F_rad_zl*u0(i,j,k)*(shift_z + v_zl) + psiL4 * P_radl * gzz(i,j,k)
              Sxy(i,j,k) = Sxy(i,j,k) + (E_rad(i,j,k)+P_radl)*u0(i,j,k)**2*(shift_x + v_xl)*(shift_y + v_yl) + u0(i,j,k)*(F_rad_xl*(shift_y + v_yl)+F_rad_yl*(shift_x + v_xl)) + psiL4 * P_radl * gxy(i,j,k)
              Sxz(i,j,k) = Sxz(i,j,k) + (E_rad(i,j,k)+P_radl)*u0(i,j,k)**2*(shift_x + v_xl)*(shift_z + v_zl) + u0(i,j,k)*(F_rad_xl*(shift_z + v_zl)+F_rad_zl*(shift_x + v_xl)) + psiL4 * P_radl * gxz(i,j,k)
              Syz(i,j,k) = Syz(i,j,k) + (E_rad(i,j,k)+P_radl)*u0(i,j,k)**2*(shift_y + v_yl)*(shift_z + v_zl) + u0(i,j,k)*(F_rad_yl*(shift_z + v_zl)+F_rad_zl*(shift_y + v_yl)) + psiL4 * P_radl * gyz(i,j,k)
              
           else

              Fasq = F_rad_0l*F_rad0(i,j,k) + F_rad_xl*F_radx(i,j,k) +  F_rad_yl*F_rady(i,j,k) +  F_rad_zl*F_radz(i,j,k)

              zeta_temp = sqrt(abs(F_rad_0l*F_rad0(i,j,k) + F_rad_xl*F_radx(i,j,k) +  F_rad_yl*F_rady(i,j,k) +  F_rad_zl*F_radz(i,j,k))/E_rad(i,j,k)**2)
              
!              zeta_cut = Erad_atm_cut*1.5
              zeta_cut = 1.0e-40
              if (E_rad(i,j,k).le.zeta_cut) then
                 zeta = 1.0
              else
                 zeta = zeta_temp
              end if
              
              if (zeta .gt. 1.0) then
                 zeta = 1.0;
              end if
              chil = 1/3.0 + zeta**2*(6.0-2.0*zeta+6.0*zeta**2)/15.0
              
              if (chil .gt. 1.0) then
                 chil = 1.0;
              end if

              zeta_rad(i,j,k) = zeta
              chi_rad(i,j,k) = chil
              
              if (E_rad(i,j,k) .lt. Erad_atm_cut) then
                 P_radxx(i,j,k) = 0.0
                 P_radyy(i,j,k) = 0.0
                 P_radzz(i,j,k) = 0.0
                 P_radxy(i,j,k) = 0.0
                 P_radxz(i,j,k) = 0.0
                 P_radyz(i,j,k) = 0.0
              else
                 if (Fasq .le. 0) then
                    P_radxx(i,j,k) = E_rad(i,j,k)*(gupxx(i,j,k)/psiL4 - shiftx(i,j,k)*shiftx(i,j,k)/(lapm1(i,j,k)+1.0)**2 + uxl**2)/2.0*(1.0-chil)
                    P_radyy(i,j,k) = E_rad(i,j,k)*(gupyy(i,j,k)/psiL4 - shifty(i,j,k)*shifty(i,j,k)/(lapm1(i,j,k)+1.0)**2 + uyl**2)/2.0*(1.0-chil)
                    P_radzz(i,j,k) = E_rad(i,j,k)*(gupzz(i,j,k)/psiL4 - shiftz(i,j,k)*shiftz(i,j,k)/(lapm1(i,j,k)+1.0)**2 + uzl**2)/2.0*(1.0-chil)
                    P_radxy(i,j,k) = E_rad(i,j,k)*(gupxy(i,j,k)/psiL4 - shiftx(i,j,k)*shifty(i,j,k)/(lapm1(i,j,k)+1.0)**2 + uxl*uyl)/2.0*(1.0-chil)
                    P_radxz(i,j,k) = E_rad(i,j,k)*(gupxz(i,j,k)/psiL4 - shiftx(i,j,k)*shiftz(i,j,k)/(lapm1(i,j,k)+1.0)**2 + uxl*uzl)/2.0*(1.0-chil)
                    P_radyz(i,j,k) = E_rad(i,j,k)*(gupyz(i,j,k)/psiL4 - shifty(i,j,k)*shiftz(i,j,k)/(lapm1(i,j,k)+1.0)**2 + uyl*uzl)/2.0*(1.0-chil)
                 else
                    P_radxx(i,j,k) = E_rad(i,j,k)*((F_radx(i,j,k)**2/Fasq)*(3.0*chil -1.0)/2.0 + (gupxx(i,j,k)/psiL4 - shiftx(i,j,k)*shiftx(i,j,k)/(lapm1(i,j,k)+1.0)**2 + uxl**2)/2.0*(1.0-chil))
                    P_radyy(i,j,k) = E_rad(i,j,k)*((F_rady(i,j,k)**2/Fasq)*(3.0*chil -1.0)/2.0 + (gupyy(i,j,k)/psiL4 - shifty(i,j,k)*shifty(i,j,k)/(lapm1(i,j,k)+1.0)**2 + uyl**2)/2.0*(1.0-chil))
                    P_radzz(i,j,k) = E_rad(i,j,k)*((F_radz(i,j,k)**2/Fasq)*(3.0*chil -1.0)/2.0 + (gupzz(i,j,k)/psiL4 - shiftz(i,j,k)*shiftz(i,j,k)/(lapm1(i,j,k)+1.0)**2 + uzl**2)/2.0*(1.0-chil))
                    P_radxy(i,j,k) = E_rad(i,j,k)*((F_radx(i,j,k)*F_rady(i,j,k)/Fasq)*(3.0*chil -1.0)/2.0 + (gupxy(i,j,k)/psiL4 - shiftx(i,j,k)*shifty(i,j,k)/(lapm1(i,j,k)+1.0)**2 + uxl*uyl)/2.0*(1.0-chil))
                    P_radxz(i,j,k) = E_rad(i,j,k)*((F_radx(i,j,k)*F_radz(i,j,k)/Fasq)*(3.0*chil -1.0)/2.0 + (gupxz(i,j,k)/psiL4 - shiftx(i,j,k)*shiftz(i,j,k)/(lapm1(i,j,k)+1.0)**2 + uxl*uzl)/2.0*(1.0-chil))
                    P_radyz(i,j,k) = E_rad(i,j,k)*((F_rady(i,j,k)*F_radz(i,j,k)/Fasq)*(3.0*chil -1.0)/2.0 + (gupyz(i,j,k)/psiL4 - shifty(i,j,k)*shiftz(i,j,k)/(lapm1(i,j,k)+1.0)**2 ++ uyl*uzl)/2.0*(1.0-chil))
                 end if
              end if
              
              P_rad0xl = - (P_radxx(i,j,k) * u_x + P_radxy(i,j,k) * u_y + P_radxz(i,j,k) * u_z)/u_0l
              P_rad0yl = - (P_radxy(i,j,k) * u_x + P_radyy(i,j,k) * u_y + P_radyz(i,j,k) * u_z)/u_0l
              P_rad0zl = - (P_radxz(i,j,k) * u_x + P_radyz(i,j,k) * u_y + P_radzz(i,j,k) * u_z)/u_0l
              P_rad00l = - (P_rad0xl * u_x + P_rad0yl * u_y + P_rad0zl * u_z)/u_0l
              
              
              tau_rad(i,j,k) = (lapm1(i,j,k)+1.0)*(lapm1(i,j,k)+1.0)*psiL**6*(E_rad(i,j,k)*u0(i,j,k)*u0(i,j,k)+2.0*F_rad0(i,j,k)*u0(i,j,k)+P_rad00l)
              S_rad_x(i,j,k) = (lapm1(i,j,k)+1.0)*psiL**6*(E_rad(i,j,k)*u0(i,j,k)*u_x + F_rad0(i,j,k)*u_x + F_rad_xl * u0(i,j,k) + P_rad00l*shift_x + psiL4*(P_rad0xl*gxx(i,j,k) + P_rad0yl*gxy(i,j,k) + P_rad0zl*gxz(i,j,k)))
              S_rad_y(i,j,k) = (lapm1(i,j,k)+1.0)*psiL**6*(E_rad(i,j,k)*u0(i,j,k)*u_y + F_rad0(i,j,k)*u_y + F_rad_yl * u0(i,j,k) + P_rad00l*shift_y + psiL4*(P_rad0xl*gxy(i,j,k) + P_rad0yl*gyy(i,j,k) + P_rad0zl*gyz(i,j,k)))
              S_rad_z(i,j,k) = (lapm1(i,j,k)+1.0)*psiL**6*(E_rad(i,j,k)*u0(i,j,k)*u_z + F_rad0(i,j,k)*u_z + F_rad_zl * u0(i,j,k) + P_rad00l*shift_z + psiL4*(P_rad0xl*gxz(i,j,k) + P_rad0yl*gyz(i,j,k) + P_rad0zl*gzz(i,j,k)))
              
              
              rho(i,j,k) = rho(i,j,k) + (lapm1(i,j,k)+1.0)**2.0*(E_rad(i,j,k)*u0(i,j,k)**2 + 2.0 * F_rad0(i,j,k) * u0(i,j,k) + P_rad00l)
              
              Sx(i,j,k) = Sx(i,j,k) + (lapm1(i,j,k)+1.0)*(u0(i,j,k)*E_rad(i,j,k)*u_x + F_rad0(i,j,k)*u_x + u0(i,j,k)*F_rad_xl +&
                   P_rad00l*shift_x + psiL4*(P_rad0xl*gxx(i,j,k) + P_rad0yl*gxy(i,j,k) + P_rad0zl*gxz(i,j,k)))
              
              Sy(i,j,k) = Sy(i,j,k) + (lapm1(i,j,k)+1.0)*(u0(i,j,k)*E_rad(i,j,k)*u_y + F_rad0(i,j,k)*u_y + u0(i,j,k)*F_rad_yl +&
                   P_rad00l*shift_y + psiL4*(P_rad0xl*gxy(i,j,k) + P_rad0yl*gyy(i,j,k) + P_rad0zl*gyz(i,j,k)))
              
              Sz(i,j,k) = Sz(i,j,k) + (lapm1(i,j,k)+1.0)*(u0(i,j,k)*E_rad(i,j,k)*u_z + F_rad0(i,j,k)*u_z + u0(i,j,k)*F_rad_zl +&
                   P_rad00l*shift_z + psiL4*(P_rad0xl*gxz(i,j,k) + P_rad0yl*gyz(i,j,k) + P_rad0zl*gzz(i,j,k)))
              
              Sxx(i,j,k) = Sxx(i,j,k) + E_rad(i,j,k)*u_x*u_x + 2.0*F_rad_xl*u_x +&
                   shift_x**2.0*P_rad00l + shift_x*2.0*(gxx(i,j,k)*P_rad0xl+gxy(i,j,k)*P_rad0yl+gxz(i,j,k)*P_rad0zl) +&
                   psiL4**2.0*( gxx(i,j,k)**2.0*P_radxx(i,j,k) + gxy(i,j,k)**2.0*P_radyy(i,j,k) + gxz(i,j,k)**2.0*P_radzz(i,j,k) +&
                   2.0*(gxx(i,j,k)*gxy(i,j,k)*P_radxy(i,j,k)+gxx(i,j,k)*gxz(i,j,k)*P_radxz(i,j,k)+gxy(i,j,k)*gxz(i,j,k)*P_radyz(i,j,k)) )
              
              Syy(i,j,k) = Syy(i,j,k) + E_rad(i,j,k)*u_y*u_y + 2.0*F_rad_yl*u_y +&
                   shift_y**2.0*P_rad00l + shift_y*2.0*(gxy(i,j,k)*P_rad0xl+gyy(i,j,k)*P_rad0yl+gyz(i,j,k)*P_rad0zl) +&
                   psiL4**2.0*( gxy(i,j,k)**2.0*P_radxx(i,j,k) + gyy(i,j,k)**2.0*P_radyy(i,j,k) + gyz(i,j,k)**2.0*P_radzz(i,j,k) +&
                   2.0*(gxy(i,j,k)*gyy(i,j,k)*P_radxy(i,j,k)+gxy(i,j,k)*gyz(i,j,k)*P_radxz(i,j,k)+gyy(i,j,k)*gyz(i,j,k)*P_radyz(i,j,k)) )
              
              Szz(i,j,k) = Szz(i,j,k) + E_rad(i,j,k)*u_z*u_z + 2.0*F_rad_zl*u_z +&
                   shift_z**2.0*P_rad00l + shift_z*2.0*(gxz(i,j,k)*P_rad0xl+gyz(i,j,k)*P_rad0yl+gzz(i,j,k)*P_rad0zl) +&
                   psiL4**2.0*( gxz(i,j,k)**2.0*P_radxx(i,j,k) + gyz(i,j,k)**2.0*P_radyy(i,j,k) + gzz(i,j,k)**2.0*P_radzz(i,j,k) +&
                   2.0*(gxz(i,j,k)*gyz(i,j,k)*P_radxy(i,j,k)+gxz(i,j,k)*gzz(i,j,k)*P_radxz(i,j,k)+gyz(i,j,k)*gzz(i,j,k)*P_radyz(i,j,k)) )
              
              Sxy(i,j,k) = Sxy(i,j,k) + E_rad(i,j,k)*u_x*u_y + F_rad_xl*u_y + F_rad_yl*u_x +&
                   shift_x*shift_y*P_rad00l +&
                   psiL4*(shift_x*(gxy(i,j,k)*P_rad0xl + gyy(i,j,k)*P_rad0yl + gyz(i,j,k)*P_rad0zl) + (shift_y*(gxx(i,j,k)*P_rad0xl + gxy(i,j,k)*P_rad0yl + gxz(i,j,k)*P_rad0zl)) )+&
                   psiL4**2.0*(gxx(i,j,k)*gxy(i,j,k)*P_radxx(i,j,k) + gxy(i,j,k)*gyy(i,j,k)*P_radyy(i,j,k) + gxz(i,j,k)*gyz(i,j,k)*P_radzz(i,j,k) +&
                   (gxx(i,j,k)*gyy(i,j,k) + gxy(i,j,k)*gxy(i,j,k))*P_radxy(i,j,k) + (gxx(i,j,k)*gyz(i,j,k) + gxz(i,j,k)*gxy(i,j,k))*P_radxz(i,j,k) + (gxy(i,j,k)*gyz(i,j,k) + gxz(i,j,k)*gyy(i,j,k))*P_radyz(i,j,k))
              
              Sxz(i,j,k) = Sxz(i,j,k) + E_rad(i,j,k)*u_x*u_z + F_rad_xl*u_z + F_rad_zl*u_x +&
                   shift_x*shift_z*P_rad00l +&
                   psiL4*(shift_x*(gxz(i,j,k)*P_rad0xl + gyz(i,j,k)*P_rad0yl + gzz(i,j,k)*P_rad0zl) + (shift_z*(gxx(i,j,k)*P_rad0xl + gxy(i,j,k)*P_rad0yl + gxz(i,j,k)*P_rad0zl)) )+&
                   psiL4**2.0*(gxx(i,j,k)*gxz(i,j,k)*P_radxx(i,j,k) + gxy(i,j,k)*gyz(i,j,k)*P_radyy(i,j,k) + gxz(i,j,k)*gzz(i,j,k)*P_radzz(i,j,k) +&
                   (gxx(i,j,k)*gyz(i,j,k) + gxy(i,j,k)*gxz(i,j,k))*P_radxy(i,j,k) + (gxx(i,j,k)*gzz(i,j,k) + gxz(i,j,k)*gxz(i,j,k))*P_radxz(i,j,k) + (gxy(i,j,k)*gzz(i,j,k) + gxz(i,j,k)*gyz(i,j,k))*P_radyz(i,j,k))

              Syz(i,j,k) = Syz(i,j,k) + E_rad(i,j,k)*u_y*u_z + F_rad_yl*u_z + F_rad_zl*u_y +&
                   shift_y*shift_z*P_rad00l +&
                   psiL4*(shift_y*(gxz(i,j,k)*P_rad0xl + gyz(i,j,k)*P_rad0yl + gzz(i,j,k)*P_rad0zl) + (shift_z*(gxy(i,j,k)*P_rad0xl + gyy(i,j,k)*P_rad0yl + gyz(i,j,k)*P_rad0zl)) )+&
                   psiL4**2.0*(gxy(i,j,k)*gxz(i,j,k)*P_radxx(i,j,k) + gyy(i,j,k)*gyz(i,j,k)*P_radyy(i,j,k) + gyz(i,j,k)*gzz(i,j,k)*P_radzz(i,j,k) +&
                   (gxy(i,j,k)*gyz(i,j,k) + gxz(i,j,k)*gyy(i,j,k))*P_radxy(i,j,k) + (gxy(i,j,k)*gzz(i,j,k) + gyz(i,j,k)*gxz(i,j,k))*P_radxz(i,j,k) + (gyy(i,j,k)*gzz(i,j,k) + gyz(i,j,k)*gyz(i,j,k))*P_radyz(i,j,k))
              

              if(isnan(tau_rad(i,j,k)) .and. isnan(S_rad_x(i,j,k))) then
                 write(*,*) "after BHNS write(*,*) initialdata_local2.F90, tau_rad or S_rad_x is nan", tau_rad(i,j,k), S_rad_x(i,j,k)
                 write(*,*) "E_rad, F_radx, P_radxx, zeta_rad, chi_rad, u_0l", E_rad(i,j,k), F_radx(i,j,k), P_radxx(i,j,k), zeta_rad(i,j,k), chi_rad(i,j,k), u_0l
                 write(*,*) "lapm1, uxl, uyl, uzl, Fasq", lapm1(i,j,k), uxl,uyl,uzl, Fasq
              end if
            
           end if
    
           
        end do
     end do
  end do

end if


end subroutine bhns_initialdata_local2






