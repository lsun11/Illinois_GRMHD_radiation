!--------------------------------------------------------------           
! Driver for primitives and boundary conditions routines, v2.0
!--------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!!$#include "GenericFD.h"

!
subroutine mhd_primitives_and_bcs_v2(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext,fake_ext,thisproc_have_global_bdry_min,thisproc_have_global_bdry_max
  real*8                :: dT,dX,dY,dZ,b2bt,inside_horizon_boundary_psi6, rho_cent
  real*8                :: rho_fail_max_step,M_fail_step,randomnumber
  real*8                :: xcenter,ycenter,zcenter,horizdirn_x,horizdirn_y,horizdirn_z,distfrombhcenter,horiz_radius
  real*8		:: xNS1,yNS1,xNS2,yNS2, rad
  integer               :: proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax
  integer               :: glob_imax,glob_jmax,glob_kmax
  integer               :: i,j,k,vindex,NFont,NFont_L
  integer               :: handle,dummy,ierr
  integer               :: repairs_needed,repairs_rad_needed,ignore_ghostzones
  CCTK_REAL             :: reduction_value
  integer               :: AXISYM, levelnumber
  parameter(AXISYM = 4)

! Storage variables:
  real*8                             :: Mx,My,Mz,tu,Psi6,AH_radius_minimum
  real*8                             :: M_PI

  M_PI = acos(-1.D0)


  print *, "^^^^^^^^^^^^^^^^^^^^^^ Start driver_primitives_and_bcs_v2 ^^^^^^^^^^^^^^^^^^^^^^^^"
  ext = cctk_lsh

  !!levelnumber = cctk_levfac(1)
  !!levelnumber = log(levelnumber)/log(2.D0)+1.D0
  levelnumber = int( log(real(cctk_levfac(1)))/log(2.0) + 1.e-5 ) + 1

  glob_imax = ext(1)
  glob_jmax = ext(2)
  glob_kmax = ext(3)
  proc_imin = -100
  proc_jmin = -100
  proc_kmin = -100
  proc_imax = -1
  proc_jmax = -1
  proc_kmax = -1


! calcualte CoM for tempearture computation.

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



  if(cctk_iteration.gt.0) then
     dT = CCTK_DELTA_TIME
     dX = CCTK_DELTA_SPACE(1)
     dY = CCTK_DELTA_SPACE(2)
     dZ = CCTK_DELTA_SPACE(3)

     thisproc_have_global_bdry_min = have_global_bdry_min(levelnumber,:)
     thisproc_have_global_bdry_max = have_global_bdry_max(levelnumber,:)

     if(thisproc_have_global_bdry_min(1)==1) proc_imin = 0
     if(thisproc_have_global_bdry_min(2)==1) proc_jmin = 0
     if(thisproc_have_global_bdry_min(3)==1) proc_kmin = 0
     if(thisproc_have_global_bdry_max(1)==1) proc_imax = glob_imax
     if(thisproc_have_global_bdry_max(2)==1) proc_jmax = glob_jmax
     if(thisproc_have_global_bdry_max(3)==1) proc_kmax = glob_kmax

     !     if(levelnumber==1 .or. primitives_solver.ne.1) then
     ! Should check if there are NaNs first...

     if(em_evolve_enable==1) then
        if(levelnumber==1) then
           ! First apply BC's to Bitilde's
           ! This is a quick and dirty fix: set thisproc_have_global_bdry_min/max so that
           !    ALL boundaries are updated, EXCEPT in symmetry ghostzones!
           !    We must apply this fix, since have_global_bdry_min/max does not save all the
           !    necessary information when being saved to checkpoint.
           !    In other words, have_global_bdry_min/max is BROKEN.
           call global_edges(dX,dY,dZ, &
                X(1,1,1),Y(1,1,1),Z(1,1,1), &
                X(1,1,1),Y(1,1,1),Z(1,1,1), &
                X(cctk_lsh(1),1,1),Y(1,cctk_lsh(2),1),Z(1,1,cctk_lsh(3)), &
                X(cctk_lsh(1),1,1),Y(1,cctk_lsh(2),1),Z(1,1,cctk_lsh(3)), &
                thisproc_have_global_bdry_min,thisproc_have_global_bdry_max,Symmetry)
           do i=1,cctk_nghostzones(1)
              fake_ext = cctk_lsh - cctk_nghostzones + i
              if(Symmetry==AXISYM) fake_ext(1) = fake_ext(1) + 1

	      if (constrained_transport_scheme .ne. 3) then 
                 call emfields_bc_newv2(ext,fake_ext,X,Y,Z, &
                      Bxtilde_p,Bytilde_p,Bztilde_p,Bxtilde,Bytilde,Bztilde,Symmetry, &
                      thisproc_have_global_bdry_min,thisproc_have_global_bdry_max,EM_BC)
	      end if
           end do
 
           if(hyperbolic_divergence_cleaning_enable==1) then
              if(1==1) then
                 do i=1,cctk_nghostzones(1)
                    fake_ext = cctk_lsh - cctk_nghostzones + i
                    if(Symmetry==AXISYM) fake_ext(1) = fake_ext(1) + 1
                    call update_boundary_Blagrangemultiplier(ext,fake_ext,X,Y,Z, &
                         PhysicalRadius, RadiusDerivative, &
                         dT,dX,dY,dZ, lapm1_p, &
                         phi_p,gxx_p,gxy_p,gxz_p,gyy_p,gyz_p,gzz_p, &
                         Blagrangemultiplier_p,Blagrangemultiplier, &
                         Symmetry, &
                         thisproc_have_global_bdry_min,thisproc_have_global_bdry_max)
                 end do
              else
                 Blagrangemultiplier(1:3,:,:) = 0.D0
                 Blagrangemultiplier(:,1:3,:) = 0.D0
                 Blagrangemultiplier(:,:,1:3) = 0.D0
                 Blagrangemultiplier(ext(1)-2:ext(1),:,:) = 0.D0
                 Blagrangemultiplier(:,ext(2)-2:ext(2),:) = 0.D0
                 Blagrangemultiplier(:,:,ext(3)-2:ext(3)) = 0.D0
              end if
           end if
        end if

        if (enable_trace_field_line==1) then
           do i=1,cctk_nghostzones(1)
              fake_ext = cctk_lsh - cctk_nghostzones + i
              call update_boundary_field_line_tracer(ext,fake_ext, Symmetry, &
                   mhd_psi_line, mhd_u_psi, mhd_chi_line, mhd_u_chi, &
                   mhd_psi_line_p, mhd_u_psi_p, mhd_chi_line_p, mhd_u_chi_p, &
                   thisproc_have_global_bdry_min,thisproc_have_global_bdry_max,EM_BC)
           end do
        end if

        call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativex')
        call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativey')
        call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativez')
        call CartSymGN(dummy,cctkGH,'mhd_evolve::em_Blagrangemultiplier')
        call CartSymGN(dummy,cctkGH,'mhd_evolve::field_line_variables')


   	!----------------------------------------
	! Impose frozen BC on the "conservative" variables 
	! before primitive recovery
        if (Matter_BC==2) then 
           if(levelnumber==1) then
              call global_edges(dX,dY,dZ, &
                   X(1,1,1),Y(1,1,1),Z(1,1,1), &
                   X(1,1,1),Y(1,1,1),Z(1,1,1), &
                   X(cctk_lsh(1),1,1),Y(1,cctk_lsh(2),1),Z(1,1,cctk_lsh(3)), &
                   X(cctk_lsh(1),1,1),Y(1,cctk_lsh(2),1),Z(1,1,cctk_lsh(3)), &
                   thisproc_have_global_bdry_min,thisproc_have_global_bdry_max,Symmetry)

              do i=1,cctk_nghostzones(1)
                 fake_ext = cctk_lsh - cctk_nghostzones + i
                 if(Symmetry==AXISYM) fake_ext(1) = fake_ext(1) + 1

                 if (constrained_transport_scheme==3) then 	
                    call apply_matter_frozen_bc(ext,fake_ext, &
                         rho_star,tau,mhd_st_x,mhd_st_y,mhd_st_z, &
                         Ax,Ay,Az, &
                         rho_star_p,tau_p,mhd_st_x_p,mhd_st_y_p,mhd_st_z_p, &
                         Ax_p,Ay_p,Az_p,thisproc_have_global_bdry_min, & 
                         thisproc_have_global_bdry_max)
	   	 else
                    call apply_matter_frozen_bc(ext,fake_ext, &
                         rho_star,tau,mhd_st_x,mhd_st_y,mhd_st_z, &
                         Bxtilde,Bytilde,Bztilde, &
                         rho_star_p,tau_p,mhd_st_x_p,mhd_st_y_p,mhd_st_z_p, &
                         Bxtilde_p,Bytilde_p,Bztilde_p, & 
                         thisproc_have_global_bdry_min, & 
                         thisproc_have_global_bdry_max)
	         end if
              end do
           end if
        end if
	!----------------------------------------

        !---------------------------------------
        !Next we recompute B^i from \tilde{B^i} or compute B^i from Ai
        !(This is the easiest primitives recovery!)
        ! Note that when constrained_transport_scheme==3, B has 
        ! already been computed.
	if (constrained_transport_scheme .ne. 3) then
	   ! Compute B^i from tilde{B}^i
           !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
                    !Note that the primitives solver requires Bx,By,Bz as input.
                    Bx(i,j,k) = Bxtilde(i,j,k)
                    By(i,j,k) = Bytilde(i,j,k)
                    Bz(i,j,k) = Bztilde(i,j,k)
                 end do
              end do
           end do
           !$omp end parallel do
           b2bt = -1.D0
           call convert_b(ext,Bx,By,Bz,phi,b2bt) 
	end if
        !---------------------------------------
     end if

     !-----------------------------------------------------------------------------------------
     ! PRIMITIVES SOLVER STEP:
     call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativex')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativey')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativez')
     !     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_conservatives')


  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
           if (compute_microphysics ==1 .and. T_fluid(i,j,k) .gt. 1.0e-30) then
              write(*,*) " In driver_primitives_and_bcs_v2.F90 (1), T_fluid too large!!!!"
              write(*,*) "i,j,k X,Y,Z, T_fluid(i,j,k) = ", i,j,k, X(i,j,k), Y(i,j,k), Z(i,j,k), T_fluid(i,j,k)
              write(*,*) "rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k)=",rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k)
              write(*,*) "eps_tot(i,j,k), eps_cld(i,j,k), eps_thermal(i,j,k), P_thermal(i,j,k)=",eps_tot(i,j,k), eps_cld(i,j,k), eps_thermal(i,j,k), P_thermal(i,j,k)
              write(*,*) "tau(i,j,k), mhd_st_x(i,j,k), mhd_st_y(i,j,k), mhd_st_z(i,j,k), u0(i,j,k)=", tau(i,j,k), mhd_st_x(i,j,k), mhd_st_y(i,j,k), mhd_st_z(i,j,k), u0(i,j,k)
              write(*,*) "Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)=", Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)
!              call EXIT(STATUS)                                                                                                                                                                                          
           end if
        end do
     end do
  end do

     if (enable_primitives_disk == 0) then
        if(primitives_solver==11) then
           if(rho_star(1,1,1) == rho_star(1,1,1) ) then
           else
              write(*,*) "BAD"
              stop
           end if
           !Here we use the HARM primitives solver, with a new prescription that minimizes changes to
           !  conservatives without applying the "Font" fix.
           !We intend to make the below function truly generic, but currently it only supports Gamma-law EOS.
           !We hope the below function will eventually be used in place of other primitives solver options,
           !  since it can be easily extended to use your desired technique for solving primitives.

  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
           if (isnan(tau(i,j,k)) .or. isnan(mhd_st_x(i,j,k))) then
              write(*,*) "inside driver_primitives_and_bcs_v2, BEFORE generic, tau and mhd_st_x=", tau(i,j,k), mhd_st_x(i,j,k)
              write(*,*) "i,j,k =", i,j,k
              write(*,*) "tau_rhs, mhd_st_x_rhs, mhd_st_y_rhs, mhd_st_z_rhs", tau_rhs(i,j,k), mhd_st_x_rhs(i,j,k), mhd_st_y_rhs(i,j,k), mhd_st_z_rhs(i,j,k)
              write(*,*) "rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k), u0(i,j,k)=", rho_b(i,j,k), P(i,j,k), vx(i,j,k), u0(i,j,k)
              write(*,*) "lapm1(i,j,k), shiftx(i,j,k), E_rad(i,j,k), F_radx(i,j,k), tau_rad(i,j,k), h(i,j,k)", lapm1(i,j,k), shiftx(i,j,k), E_rad(i,j,k),F_radx(i,j,k)\
              , tau_rad(i,j,k), h(i,j,k)
              write(*,*) "eps_cld(i,j,k), eps_thermal(i,j,k), eps_tot(i,j,k), P_cld(i,j,k), P_thermal(i,j,k), T_fluid(i,j,k)=", eps_cld(i,j,k), eps_thermal(i,j,k), eps_tot(i,j,k), P_cld(i,j,k), P_thermal(i,j,k),T_fluid(i,j,k)
           end if
        end do
     end do
  end do




!	  print *, "^^^^^^^^^^^^^^^^^^^^^^ Inside driver_v2, start primsol_generic + metric source ^^^^^^^^^^^^^^^^^^^^^^^^"        
  	  call primitives_generic(cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z, &
                phi,gxx,gxy,gxz,gyy,gyz,gzz, &
                gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                lapm1,shiftx,shifty,shiftz, &
                Bx,By,Bz, &
                eps_tot, eps_thermal, eps_cld, P_cld,&
                ka_gf, ks_gf, emission_gf, chi_rad, chi_rad_nue, Y_e, optd, eta_nue,&
                ka_gf_nue, ks_gf_nue, emission_gf_nue, ka_gf_nux,ks_gf_nux, emission_gf_nux,&
                mhd_st_x,mhd_st_y,mhd_st_z,tau,rho_star, rhoYe, &
                vx,vy,vz,P,rho_b,h,u0, &
                rho_b_atm, tau_atm, rho_b_max,&
		xNS1, yNS1, xNS2, yNS2,  M_B, rad_T_fac, rad_T_cutoff, rad_T_pow, rad_T_floor, T_fluid,&
                neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
                temp1,temp2,temp3,temp4,temp5, &
                primitives_debug,Psi6threshold,rad_const,horizon_enforce_rho_profile,enable_OS_collapse, &
		rad_evolve_enable, compute_microphysics, &
                microphysics_scheme, T_fluid_cgs_atm, rad_fix)

           !The goal of the below function is to just update the metric source terms.  We should
           !  really get rid of the old h,w,st_i,Ei gridfunctions...
           !Currently this only works for Gamma-law EOS's but could be extended for arbitary EOSs
           !  with some work.
           call metric_source_terms_and_misc_vars(cctkGH,cctk_lsh, &
                rho,Sx,Sy,Sz, &
                Sxx,Sxy,Sxz,Syy,Syz,Szz, &
		tau_rad, S_rad_x,S_rad_y, S_rad_z,&
                rhoYe, Y_e,&
                E_rad, F_radx, F_rady, F_radz, F_rad0, F_rad_scalar,&
                h,w,st_x,st_y,st_z, &
                Ex,Ey,Ez, &
                sbt,sbx,sby,sbz, &
                lapm1,shiftx,shifty,shiftz, &
                phi,gxx,gxy,gxz,gyy,gyz,gzz, &
                gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                P,rho_b,u0,vx,vy,vz, &
                Bx,By,Bz, &
                rho_star,mhd_st_x,mhd_st_y,mhd_st_z,&
                neos,ergo_star, ergo_sigma,&
                rho_tab, P_tab,eps_tab,k_tab,gamma_tab,&
                rho_b_atm, enable_OS_collapse, rad_evolve_enable, rad_closure_scheme, Psi6threshold, Erad_atm_cut)
          
!	print *, "^^^^^^^^^^^^^^^^^^^^^^ Inside driver_v2, END primsol_generic + metric source ^^^^^^^^^^^^^^^^^^^^^^^^"
           ! After getting ka_gf/ks_gf, compute optical depth.
!           call compute_optical_depth(ext, X,Y,Z, dX,dY,dZ, optd, optd_x, optd_y, optd_z, ka_gf, ks_gf)



  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
           if (isnan(tau(i,j,k)) .or. isnan(mhd_st_x(i,j,k))) then
              write(*,*) "inside driver_primitives_and_bcs_v2, BEFORE rad, tau and mhd_st_x=", tau(i,j,k), mhd_st_x(i,j,k)
              write(*,*) "i,j,k =", i,j,k
              write(*,*) "tau_rhs, mhd_st_x_rhs, mhd_st_y_rhs, mhd_st_z_rhs", tau_rhs(i,j,k), mhd_st_x_rhs(i,j,k), mhd_st_y_rhs(i,j,k), mhd_st_z_rhs(i,j,k)
              write(*,*) "rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k), u0(i,j,k)=", rho_b(i,j,k), P(i,j,k), vx(i,j,k), u0(i,j,k)
              write(*,*) "lapm1(i,j,k), shiftx(i,j,k), E_rad(i,j,k), F_radx(i,j,k), tau_rad(i,j,k), h(i,j,k)", lapm1(i,j,k), shiftx(i,j,k), E_rad(i,j,k),F_radx(i,j,k)\
              , tau_rad(i,j,k), h(i,j,k)
              write(*,*) "eps_cld(i,j,k), eps_thermal(i,j,k), eps_tot(i,j,k), P_cld(i,j,k), P_thermal(i,j,k), T_fluid(i,j,k)=", eps_cld(i,j,k), eps_thermal(i,j,k), eps_tot(i,j,k), P_cld(i,j,k), P_thermal(i,j,k),T_fluid(i,j,k)
           end if
        end do
     end do
  end do



if(rad_evolve_enable==1) then
   do k=1, ext(3)
      do j=1, ext(2)
         do i=1, ext(1)
            if ( tau_rad(i,j,k) .lt. 0.d0) then
               tau_rad(i,j,k) = Erad_atm_cut*psi(i,j,k)**6*(lapm1(i,j,k)+1.0)*u0(i,j,k)*u0(i,j,k)
               S_rad_x(i,j,k) = 0.d0
               S_rad_y(i,j,k) = 0.d0
               S_rad_z(i,j,k) = 0.d0
            end if
         end do
      end do
   end do
end if


  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
           if (compute_microphysics ==1 .and. T_fluid(i,j,k) .gt. 1.0e-30) then
              write(*,*) " In driver_primitives_and_bcs_v2.F90 (2), T_fluid too large!!!!"
              write(*,*) "i,j,k X,Y,Z, T_fluid(i,j,k) = ", i,j,k, X(i,j,k), Y(i,j,k), Z(i,j,k), T_fluid(i,j,k)
              write(*,*) "rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k)=",rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k)
              write(*,*) "eps_tot(i,j,k), eps_cld(i,j,k), eps_thermal(i,j,k), P_thermal(i,j,k)=",eps_tot(i,j,k), eps_cld(i,j,k), eps_thermal(i,j,k), P_thermal(i,j,k)
              write(*,*) "tau(i,j,k), mhd_st_x(i,j,k), mhd_st_y(i,j,k), mhd_st_z(i,j,k), u0(i,j,k)=", tau(i,j,k), mhd_st_x(i,j,k), mhd_st_y(i,j,k), mhd_st_z(i,j,k), u0(i,j,k)
              write(*,*) "Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)=", Bx(i,j,k), By(i,j,k), Bz(i,j,k), smallb2(i,j,k)
!              call EXIT(STATUS)     
           end if
        end do
     end do
  end do



        ignore_ghostzones = 0
   if(rad_evolve_enable==1 .and. rad_closure_scheme == 1) then
           call  primitive_vars_rad_cpp(ext,cctk_nghostzones, X, Y, Z, &
                rho_star, &
                tau_rad,S_rad_x,S_rad_y,S_rad_z, &
                Sx, Sy, Sz, rho, &
                Sxx, Sxy, Sxz, Syy, Syz, Szz,&
                E_rad,F_radx,F_rady,F_radz, F_rad0, F_rad_scalar,&
                P_radxx, P_radyy, P_radzz, P_radxy, P_radxz, P_radyz,&
                phi, lapm1, shiftx, shifty, shiftz, gxx, gxy, gxz, gyy, gyz, gzz, &
                gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
                vx, vy, vz, u0, chi_rad, zeta_rad,&
                cctkGH, temp6, ignore_ghostzones, repairs_rad_needed, Psi6threshold, Erad_atm_cut)
     
! Here we only use this function to repair radiation!!!
        if (repairs_rad_needed==1) then
           write(*,*) "Repair Electron Anti!!!!!!!!!!!!!!"
           call repair_failures_mhd_hybrid(ext,Z,gamma_th, &
                temp10, temp6, rho_b, P, &
                vx, vy, vz, u0, w, h, rho_star, tau, &
                st_x, st_y, st_z, mhd_st_x, mhd_st_y, mhd_st_z, &
                rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
                E_rad,F_radx,F_rady,F_radz,F_rad0, F_rad_scalar,&
                P_radxx, P_radyy, P_radzz, P_radxy, P_radxz, P_radyz,&
                tau_rad, S_rad_x, S_rad_y, S_rad_z, &
                Erad_atm_cut, rad_closure_scheme, enable_OS_collapse, &
                lapm1, shiftx, shifty, shiftz, phi, &
                gxx, gxy, gxz, gyy, gyz, gzz, &
                gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Bx, By, Bz, &
                sbt, sbx, sby, sbz, rho_b_atm, &
                neos,ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
                proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
                glob_imax,glob_jmax,glob_kmax,Symmetry, primitives_solver)
	end if	
     


!!! Evolve other specices of neutrino if microphysics_scheme == 1

     if(microphysics_scheme==1) then
     	write(*,*) "Inside driver_primitives_and_bcs_v2.F90, microphysics_scheme, rad_evolve_enable, and rad_closure_scheme=", microphysics_scheme, rad_evolve_enable, rad_closure_scheme
        call  primitive_vars_rad_cpp(ext,cctk_nghostzones, X, Y, Z, &
             rho_star, &
             tau_rad_nue,S_rad_x_nue,S_rad_y_nue,S_rad_z_nue, &
             Sx, Sy, Sz, rho, &
             Sxx, Sxy, Sxz, Syy, Syz, Szz,&
             E_rad_nue,F_radx_nue,F_rady_nue,F_radz_nue, F_rad0_nue, F_rad_scalar_nue,&
             P_radxx_nue, P_radyy_nue, P_radzz_nue, P_radxy_nue, P_radxz_nue, P_radyz_nue,&
             phi, lapm1, shiftx, shifty, shiftz, gxx, gxy, gxz, gyy, gyz, gzz, &
             gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
             vx, vy, vz, u0, chi_rad_nue, zeta_rad_nue,&
             cctkGH, temp6, ignore_ghostzones, repairs_rad_needed, Psi6threshold, Erad_atm_cut)
        
! Here we only use this function to repair radiation!!!
        if (repairs_rad_needed==1) then
           write(*,*) "Repair NUE!!!!!!!!!!!!!!"
           call repair_failures_mhd_hybrid(ext,Z,gamma_th, &
                temp10, temp6, rho_b, P, &
                vx, vy, vz, u0, w, h, rho_star, tau, &
                st_x, st_y, st_z, mhd_st_x, mhd_st_y, mhd_st_z, &
                rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
                E_rad_nue,F_radx_nue,F_rady_nue,F_radz_nue,F_rad0_nue, F_rad_scalar_nue,&
                P_radxx_nue, P_radyy_nue, P_radzz_nue, P_radxy_nue, P_radxz_nue, P_radyz_nue,&
                tau_rad_nue, S_rad_x_nue, S_rad_y_nue, S_rad_z_nue, &
                Erad_atm_cut, rad_closure_scheme, enable_OS_collapse, &
                lapm1, shiftx, shifty, shiftz, phi, &
                gxx, gxy, gxz, gyy, gyz, gzz, &
                gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Bx, By, Bz, &
                sbt, sbx, sby, sbz, rho_b_atm, &
                neos,ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
                proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
                glob_imax,glob_jmax,glob_kmax,Symmetry, primitives_solver)
        end if


        call  primitive_vars_rad_cpp(ext,cctk_nghostzones, X, Y, Z, &
             rho_star, &
             tau_rad_nux,S_rad_x_nux,S_rad_y_nux,S_rad_z_nux, &
             Sx, Sy, Sz, rho, &
             Sxx, Sxy, Sxz, Syy, Syz, Szz,&
             E_rad_nux,F_radx_nux,F_rady_nux,F_radz_nux, F_rad0_nux, F_rad_scalar_nux,&
             P_radxx_nux, P_radyy_nux, P_radzz_nux, P_radxy_nux, P_radxz_nux, P_radyz_nux,&
             phi, lapm1, shiftx, shifty, shiftz, gxx, gxy, gxz, gyy, gyz, gzz, &
             gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
             vx, vy, vz, u0, chi_rad_nux, zeta_rad_nux,&
             cctkGH, temp6, ignore_ghostzones, repairs_rad_needed, Psi6threshold, Erad_atm_cut)
        
        ! Here we only use this function to repair radiation!!!   
        if (repairs_rad_needed==1) then
           write(*,*) "Repair NUX!!!!!!!!!!!!!!"
           call repair_failures_mhd_hybrid(ext,Z,gamma_th, &
                temp10, temp6, rho_b, P, &
                vx, vy, vz, u0, w, h, rho_star, tau, &
                st_x, st_y, st_z, mhd_st_x, mhd_st_y, mhd_st_z, &
                rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
                E_rad_nux,F_radx_nux,F_rady_nux,F_radz_nux,F_rad0_nux, F_rad_scalar_nux,&
                P_radxx_nux, P_radyy_nux, P_radzz_nux, P_radxy_nux, P_radxz_nux, P_radyz_nux,&
                tau_rad_nux, S_rad_x_nux, S_rad_y_nux, S_rad_z_nux, &
                Erad_atm_cut, rad_closure_scheme, enable_OS_collapse, &
                lapm1, shiftx, shifty, shiftz, phi, &
                gxx, gxy, gxz, gyy, gyz, gzz, &
                gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Bx, By, Bz, &
                sbt, sbx, sby, sbz, rho_b_atm, &
                neos,ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
                proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
                glob_imax,glob_jmax,glob_kmax,Symmetry, primitives_solver)
        end if

     end if
  end if


  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
           if (isnan(tau(i,j,k)) .or. isnan(mhd_st_x(i,j,k))) then
              write(*,*) "inside driver_primitives_and_bcs_v2, AFTER rad, tau and mhd_st_x=", tau(i,j,k), mhd_st_x(i,j,k)
              write(*,*) "i,j,k =", i,j,k
              write(*,*) "tau_rhs, mhd_st_x_rhs, mhd_st_y_rhs, mhd_st_z_rhs", tau_rhs(i,j,k), mhd_st_x_rhs(i,j,k), mhd_st_y_rhs(i,j,k), mhd_st_z_rhs(i,j,k)
              write(*,*) "rho_b(i,j,k), P(i,j,k), vx(i,j,k), vy(i,j,k), vz(i,j,k), u0(i,j,k)=", rho_b(i,j,k), P(i,j,k), vx(i,j,k), u0(i,j,k)
              write(*,*) "lapm1(i,j,k), shiftx(i,j,k), E_rad(i,j,k), F_radx(i,j,k), tau_rad(i,j,k), h(i,j,k)", lapm1(i,j,k), shiftx(i,j,k), E_rad(i,j,k),F_radx(i,j,k)\
              , tau_rad(i,j,k), h(i,j,k)
              write(*,*) "eps_cld(i,j,k), eps_thermal(i,j,k), eps_tot(i,j,k), P_cld(i,j,k), P_thermal(i,j,k), T_fluid(i,j,k)=", eps_cld(i,j,k), eps_thermal(i,j,k), eps_tot(i,j,k), P_cld(i,j,k), P_thermal(i,j,k),T_fluid(i,j,k)
           end if
        end do
     end do
  end do
  





	else if(primitives_solver==1) then
           !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)

                    !Here, temp2 is supposed to be w_p
                    temp2(i,j,k) = u0(i,j,k)*(lapm1(i,j,k)+1.D0)*rho_star(i,j,k)
                 end do
              end do
           end do
           !$omp end parallel do

           !$omp parallel do
           do k=cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
              do j=cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
                 do i=cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
                    if(temp2(i,j,k).ne.temp2(i,j,k)) then
                       do vindex=1,100
                          write(*,*) "ERROR AT i,j,k:",i,j,k
                          write(*,*) "THERE IS A NAN!  u0,lapm1,rho_star=",u0(i,j,k),lapm1(i,j,k),rho_star(i,j,k)
                       end do
                    end if

                    if(psi(i,j,k).eq.0.D0) then
                       write(*,*) "BAD PSI at:",x(i,j,k),y(i,j,k),z(i,j,k),psi(i,j,k),phi(i,j,k)
                       write(*,*) "BAD PSI gij:",gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
                    end if
                 end do
              end do
           end do
           !$omp end parallel do

           call hydro_primitives(cctkGH,ext,cctk_nghostzones, X, Y, Z, rho_star, tau,&
                mhd_st_x, mhd_st_y, mhd_st_z,&
                u0,vx,vy,vz,&
                w, temp2, rho_b,rho, P, h, &
                Sx, Sy, Sz, &
                Sxx, Sxy, Sxz, Syy, Syz, Szz, &
                phi, lapm1, shiftx,shifty,shiftz, &
                gxx, gxy, gxz, gyy, gyz, gzz, &
                gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,&
                tau_stildefix_enable,tau_atm,enable_shocktest_primitive_mode,&
                rho_b_max, rho_fail_max_step, M_fail_step, rho_b_atm,&
                gamma_th,K_poly,sdots_o_rhot,Symmetry,0)

           ! Following lines are necessary
           !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
                    st_x(i,j,k) = mhd_st_x(i,j,k)
                    st_y(i,j,k) = mhd_st_y(i,j,k)
                    st_z(i,j,k) = mhd_st_z(i,j,k)
                 end do
              end do
           end do
           !$omp end parallel do

        else if(primitives_solver==0) then

           ! Initial guess for h:
           !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
	            h_p(i,j,k) = h(i,j,k)
                 end do
              end do
           end do
           !$omp end parallel do

           !TODO: remove temp1 (unused) parameter.
           call primitive_vars_hybrid2(ext,cctk_nghostzones,X,Y,Z, &
                rho_star,tau,st_x,st_y,st_z, &
                mhd_st_x,mhd_st_y,mhd_st_z,&
		neos, ergo_star, ergo_sigma, &
                rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, &
                w,temp1,rho_b,rho,P,h,Sx,Sy,Sz, &
                Sxx,Sxy,Sxz,Syy,Syz,Szz, &
                phi,lapm1,shiftx,shifty,shiftz, &
                gxx,gxy,gxz,gyy,gyz,gzz, &
                gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                h_p,u0,rho_b_max,rho_b_atm, &
                rho_fail_max_step,M_fail_step,rhos_max, &
                Bx,By,Bz,Ex,Ey,Ez, &
                vx,vy,vz, &
                sbt,sbx,sby,sbz, temp4, &
                proc_imin,proc_jmin,proc_kmin, &
                proc_imax,proc_jmax,proc_kmax, &
                glob_imax,glob_jmax,glob_kmax, &
                Symmetry,pfloor,excision_enable, &
                excision_zone_gf, tau_stildefix_enable,tau_atm,0)
           
                
        else if(primitives_solver==2) then
	   !!print *, "^^^^^^^^^^^^^^^^^^^^^^ Use primitives solver 2!!!!  ^^^^^^^^^^^^^^^^^^^^^^^^"

           !primitive solver (next function call) overwrites h_p.  We want to keep _p values.
           !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
	 	    h_p(i,j,k) = h(i,j,k)
                 end do
              end do
           end do
           !$omp end parallel do

           ignore_ghostzones = 0

           if((MOD(cctk_iteration,64)==0) .and. 1==0) then
              xcenter = bh_posn_x(1)
              ycenter = bh_posn_y(1)
              zcenter = bh_posn_z(1)

              !     Get horizon radius in direction where it is likely to be minimized:
              horizdirn_x = 0.D0
              horizdirn_y = 0.D0
              horizdirn_z = 100000.D0
              call get_ah_radius_in_dirn(cctkGH,horizdirn_x,horizdirn_y,horizdirn_z,horiz_radius)

              do k=1,cctk_lsh(3)
                 do j=1,cctk_lsh(2)
                    do i=1,cctk_lsh(1)
                       distfrombhcenter = sqrt((x(i,j,k)-xcenter)**2 + (y(i,j,k)-ycenter)**2 + (z(i,j,k)-zcenter)**2)
                       if(distfrombhcenter.lt.horiz_radius*0.5D0) then
                          rho_star(i,j,k) = -1.D0
                       end if
                    end do
                 end do
              end do
           end if 

!           if (enable_OS_collapse==1 .and. rad_evolve_enable==1 .and. rad_closure_scheme==0) then
           if (enable_OS_collapse==1 .and. rad_evolve_enable==1.and.rad_closure_scheme==0) then
              write(*,*) "Start to cutoff rad variables outside OS star, at radius grater than", OS_surf_rad
              do k=1,cctk_lsh(3)
                 do j=1,cctk_lsh(2)
                    do i=1,cctk_lsh(1)
                       if(isnan(tau_rad(i,j,k)).or.isnan(S_rad_x(i,j,k))) then
                          write(*,*) "In driver_primitives_and_bcs_v2, before primitive_vars_hybrid2_cpp, tau_rad and S_rad_x=", tau_rad(i,j,k), S_rad_x(i,j,k)
                          write(*,*) "E_rad, F_radx, chi_rad, u0, vx=", E_rad(i,j,k), F_radx(i,j,k), chi_rad(i,j,k), u0(i,j,k), vx(i,j,k)
                          write(*,*) "X,Y,Z=", X(i,j,k), Y(i,j,k), Z(i,j,k)
                       end if
                       rad = sqrt(X(i,j,k)**2+Y(i,j,k)**2+Z(i,j,k)**2)
                       !if (rho_b(i,j,k) .lt. rho_cent*Erad_cut) then                                                                                                            
                       if (rad .gt. OS_surf_rad) then
                          tau_rad(i,j,k) = 0.d0
!                          T_fluid(i,j,k) = 0.d0                                                                                                                                 
                          S_rad_x(i,j,k) = 0.d0
                          S_rad_y(i,j,k) = 0.d0
                          S_rad_z(i,j,k) = 0.d0
                       end if
                       end do
                    end do
                 end do
                 write(*,*) "end of cutoff rad variables outside OS star"
           end if

              do k=1,cctk_lsh(3)
                 do j=1,cctk_lsh(2)
                    do i=1,cctk_lsh(1)
                       if(isnan(tau_rad(i,j,k)).or.isnan(S_rad_x(i,j,k))) then
                          write(*,*) "In driver_primitives_and_bcs_v2, before primitive_vars_hybrid2_cpp, tau_rad and S_rad_x=", tau_rad(i,j,k), S_rad_x(i,j,k)
                          write(*,*) "E_rad, F_radx, P_radxx, chi_rad, u0, vx=", E_rad(i,j,k), F_radx(i,j,k), P_radxx(i,j,k), chi_rad(i,j,k), u0(i,j,k),vx(i,j,k)
                          write(*,*) "X,Y,Z=", X(i,j,k), Y(i,j,k), Z(i,j,k)
                       end if
                    end do
                 end do
              end do

             
           call primitive_vars_hybrid2_cpp(ext,cctk_nghostzones,X,Y,Z, &
                rho_star,tau,st_x,st_y,st_z, &
                mhd_st_x,mhd_st_y,mhd_st_z, rhoYe,&
	        tau_rad,S_rad_x,S_rad_y,S_rad_z,&
		neos, ergo_star, ergo_sigma,&
                rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, &
                w,temp1,rho_b,rho,P,h,Sx,Sy,Sz, &
                Sxx,Sxy,Sxz,Syy,Syz,Szz, &
		E_rad,F_radx,F_rady,F_radz, F_rad0, F_rad_scalar, Y_e, eps_thermal, &
                phi,lapm1,shiftx,shifty,shiftz, &
                gxx,gxy,gxz,gyy,gyz,gzz, &
                gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                h_p,u0,rho_b_max,rho_b_atm, &
                rho_fail_max_step,M_fail_step,rhos_max, &
                Bx,By,Bz,Ex,Ey,Ez, &
                vx,vy,vz, &
                sbt,sbx,sby,sbz,&
		xNS1,yNS2, xNS2, yNS2,  M_B, rad_T_fac,rad_T_cutoff, rad_T_pow, rad_T_floor, T_fluid,&
                proc_imin,proc_jmin,proc_kmin, &
                proc_imax,proc_jmax,proc_kmax, &
                glob_imax,glob_jmax,glob_kmax, &
                Symmetry,pfloor,excision_enable, &
                excision_zone_gf, tau_stildefix_enable,tau_atm,temp4,cctkGH,ignore_ghostzones, &
                enable_shocktest_primitive_mode,repairs_needed,rad_closure_scheme, rad_const, &
                enable_OS_collapse, compute_microphysics, Psi6threshold, Erad_atm_cut)

           !impose flux-less BC for OS_stars (maybe good for other models as well)    
           call CCTK_VarIndex (vindex, "mhd_evolve::rho_b")
           call interp_center_carp(cctkGH,vindex,rho_cent)


           if(rad_evolve_enable==1 .and. rad_closure_scheme == 1) then 
              call   primitive_vars_rad_cpp(ext,cctk_nghostzones, X, Y, Z, &
                     rho_star, &
                     tau_rad,S_rad_x,S_rad_y,S_rad_z, &
                     Sx, Sy, Sz, rho, &
                     Sxx, Sxy, Sxz, Syy, Syz, Szz,& 
                     E_rad,F_radx,F_rady,F_radz, F_rad0, F_rad_scalar,&
                     P_radxx, P_radyy, P_radzz, P_radxy, P_radxz, P_radyz,&
                     phi, lapm1, shiftx, shifty, shiftz, gxx, gxy, gxz, gyy, gyz, gzz, &
                     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
                     vx, vy, vz, u0, chi_rad, zeta_rad,&
                     cctkGH,temp5,ignore_ghostzones, repairs_rad_needed, Psi6threshold, Erad_atm_cut)
           end if 


           if(1==0) then
              inside_horizon_boundary_psi6 = 10.D0
              call hydro_primitives_inside_horizon(cctkGH,ext,cctk_nghostzones, X, Y, Z, rho_star, tau,&
                   mhd_st_x, mhd_st_y, mhd_st_z,&
                   u0,vx,vy,vz,&
                   w, temp2, rho_b,rho, P, h, &
                   Sx, Sy, Sz, &
                   Sxx, Sxy, Sxz, Syy, Syz, Szz, &
                   phi, lapm1, shiftx,shifty,shiftz, &
                   gxx, gxy, gxz, gyy, gyz, gzz, &
                   gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,&
                   tau_stildefix_enable,tau_atm,enable_shocktest_primitive_mode,&
                   rho_b_max, rho_fail_max_step, M_fail_step, rho_b_atm,&
                   gamma_th,K_poly,sdots_o_rhot,Symmetry,0,inside_horizon_boundary_psi6)

              ! Following lines are necessary
              !$omp parallel do
              do k=1,cctk_lsh(3)
                 do j=1,cctk_lsh(2)
                    do i=1,cctk_lsh(1)
                       if(exp(phi(i,j,k)*6.D0).gt.inside_horizon_boundary_psi6) then
                          st_x(i,j,k) = mhd_st_x(i,j,k)
                          st_y(i,j,k) = mhd_st_y(i,j,k)
                          st_z(i,j,k) = mhd_st_z(i,j,k)
                       end if
                    end do
                 end do
              end do
              !$omp end parallel do
           end if

           if (repairs_needed==1 .or. repairs_rad_needed==1) then 
	      call repair_failures_mhd_hybrid(ext,Z,gamma_th, & 
                   temp4, temp5, rho_b, P, &
                   vx, vy, vz, u0, w, h, rho_star, tau, &
                   st_x, st_y, st_z, mhd_st_x, mhd_st_y, mhd_st_z, &
                   rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
                   E_rad,F_radx,F_rady,F_radz,F_rad0, F_rad_scalar,&
                   P_radxx, P_radyy, P_radzz, P_radxy, P_radxz, P_radyz,&
                   tau_rad, S_rad_x, S_rad_y, S_rad_z, &
                   Erad_atm_cut, rad_closure_scheme, enable_OS_collapse, &
                   lapm1, shiftx, shifty, shiftz, phi, &
                   gxx, gxy, gxz, gyy, gyz, gzz, &
                   gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Bx, By, Bz, &
                   sbt, sbx, sby, sbz, rho_b_atm, &
                   neos,ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
                   proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
                   glob_imax,glob_jmax,glob_kmax,Symmetry,primitives_solver)
           end if

           if(excision_enable==1) then
              call hydro_ezbc_hybrid(ext,X,Y,Z,rho_star,tau, &
                   mhd_st_x,mhd_st_y,mhd_st_z,st_x,st_y,st_z, &
                   rho_b,P,h,vx,vy,vz,w,&
                   sbt,sbx,sby,sbz,Bx,By,Bz, &
                   lapm1,shiftx,shifty,shiftz,phi,&
                   gxx,gxy,gxz,gyy,gyz,gzz, &
                   gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,&
                   Symmetry,excision_zone_gf,gamma_th,neos,ergo_star,ergo_sigma,rho_tab,&
                   P_tab,eps_tab,k_tab,gamma_tab)
              if (hyperbolic_divergence_cleaning_enable==1) then 
                 call scalar_excision_bc(ext,X,Y,Z,Blagrangemultiplier,Symmetry,excision_zone_gf)
              end if
              call remove_interior2(ext,X,Y,Z,sbt,excision_zone_gf,Symmetry)
              call remove_interior2(ext,X,Y,Z,sbx,excision_zone_gf,Symmetry)
              call remove_interior2(ext,X,Y,Z,sby,excision_zone_gf,Symmetry)
              call remove_interior2(ext,X,Y,Z,sbz,excision_zone_gf,Symmetry)
           end if

  	   ! Note: No need to apply sym BC on Bitilde since they are 
  	   ! not modified. 
!!!!call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_conservatives')
           !!call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativex')
           !!call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativey')
           !!call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativez')
           call CartSymGN(dummy,cctkGH,'mhd_evolve::em_Blagrangemultiplier')
           call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_primitives')
           call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vs')
           call CartSymGN(dummy,cctkGH,'bssn::BSSN_matter')
           call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_primitives')
           call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_pressure')
           call CartSymGN(dummy,cctkGH,'mhd_evolve::microphys_primitives')
!           call CartSymGN(dummy,cctkGH,'mhd_evolve::temperatures')
           call CartSymGN(dummy,cctkGH,'mhd_evolve::micphys_conservatives') 
           call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_conservatives')
           call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_conservatives')
        end if
     else if(enable_primitives_disk==1 .and. use_harm_primitives==0) then

        call primitive_vars_alt_disk(ext,X,Y,Z, &
             rho_star,tau,st_x,st_y,st_z, &
             mhd_st_x,mhd_st_y,mhd_st_z,neos,ergo_star,ergo_sigma, &
             rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, &
             w,rho_b,rho,P,h,Sx,Sy,Sz, &
             Sxx,Sxy,Sxz,Syy,Syz,Szz, &
             phi,lapm1,shiftx,shifty,shiftz, &
             gxx,gxy,gxz,gyy,gyz,gzz, &
             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
             h_p,u0,rho_b_max,rho_b_atm_gf, &
             rho_fail_max_step,M_fail_step, &
             Bx,By,Bz,Ex,Ey,Ez, &
             vx,vy,vz, &
             sbt,sbx,sby,sbz, &
             proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
             glob_imax,glob_jmax,glob_kmax,Symmetry, &
             Fontfix_tracker_gf,pfloor_gf, NFont_L, &
             enable_HARM_energyvariable, excision_zone_gf, force_font_fix_fail, excision_enable)

        call CCTK_ReductionHandle(handle,"sum")
        call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,NFont_L,NFont,CCTK_VARIABLE_REAL)

     end if

     !-----------------------------------------------------------------------------------------
     ! OUTER BOUNDARY UPDATE STEP:
     !     if(levelnumber==1) then

     if(1==0) then
        call random_number( randomnumber )
        call random_number( randomnumber )
        call random_number( randomnumber )
        write(*,*) "howdy.", randomnumber
        Bx = (1.D0 + 1.D-15*randomnumber)*Bx
        call random_number( randomnumber )
        write(*,*) "howdy.", randomnumber
        By = (1.D0 + 1.D-15*randomnumber)*By
        call random_number( randomnumber )
        write(*,*) "howdy.", randomnumber
        Bz = (1.D0 + 1.D-15*randomnumber)*Bz
        call random_number( randomnumber )
        write(*,*) "howdy.", randomnumber
        vx = (1.D0 + 1.D-15*randomnumber)*vx
        call random_number( randomnumber )
        write(*,*) "howdy.", randomnumber
        vy = (1.D0 + 1.D-15*randomnumber)*vy
        call random_number( randomnumber )
        write(*,*) "howdy.", randomnumber
        vz = (1.D0 + 1.D-15*randomnumber)*vz
        call random_number( randomnumber )
        write(*,*) "howdy.", randomnumber
        rho_b = (1.D0 + 1.D-15*randomnumber)*rho_b
        call random_number( randomnumber )
        write(*,*) "howdy.", randomnumber
        u0    = (1.D0 + 1.D-15*randomnumber)*u0
     end if



     if(levelnumber==1 .and. Matter_BC .ne. 2) then
        !write(*,*) "EDGESMHD:",thisproc_have_global_bdry_min,thisproc_have_global_bdry_max
        !thisproc_have_global_bdry_min = have_global_bdry_min(levelnumber,:)
        !thisproc_have_global_bdry_max = have_global_bdry_max(levelnumber,:)

        ! This is a quick and dirty fix: set thisproc_have_global_bdry_min/max so that
        !    ALL boundaries are updated, EXCEPT in symmetry ghostzones!
        call global_edges(dX,dY,dZ, &
             X(1,1,1),Y(1,1,1),Z(1,1,1), &
             X(1,1,1),Y(1,1,1),Z(1,1,1), &
             X(cctk_lsh(1),1,1),Y(1,cctk_lsh(2),1),Z(1,1,cctk_lsh(3)), &
             X(cctk_lsh(1),1,1),Y(1,cctk_lsh(2),1),Z(1,1,cctk_lsh(3)), &
             thisproc_have_global_bdry_min,thisproc_have_global_bdry_max,Symmetry)


              do k=1,cctk_lsh(3)
                 do j=1,cctk_lsh(2)
                    do i=1,cctk_lsh(1)
                       if(isnan(tau_rad(i,j,k)).or.isnan(S_rad_x(i,j,k))) then
                          write(*,*) "In driver_primitives_and_bcs_v2, before apply_matter_bc_mhd_new, tau_rad and S_rad_x=", tau_rad(i,j,k), S_rad_x(i,j,k)
                          write(*,*) "E_rad, F_radx, P_radxx, chi_rad, u0, vx=", E_rad(i,j,k), F_radx(i,j,k), P_radxx(i,j,k), chi_rad(i,j,k), u0(i,j,k),vx(i,j,k)
                          write(*,*) "X,Y,Z=", X(i,j,k), Y(i,j,k), Z(i,j,k)
                       end if
                    end do
                 end do
              end do

              if(use_disk_bcs==0) then
                 if(Matter_BC.le.6) then
                    do i=1,cctk_nghostzones(1)
                       fake_ext = cctk_lsh - cctk_nghostzones + i
                       if(Symmetry==AXISYM) fake_ext(1) = fake_ext(1) + 1
                       call apply_matter_bc_mhd_new(ext,fake_ext,X,Y,Z, &
                            rho_star,tau,mhd_st_x,mhd_st_y,mhd_st_z, &
                            st_x,st_y,st_z, &
                            tau_rad, S_rad_x, S_rad_y, S_rad_z, rhoYe,&
                            tau_rad_nue, S_rad_x_nue, S_rad_y_nue, S_rad_z_nue, tau_rad_nux, S_rad_x_nux, S_rad_y_nux, S_rad_z_nux, &
                            rho_b,P,vx,vy,vz,Bx,By,Bz, &
                            E_rad, F_radx, F_rady, F_radz, Y_e,&
                            E_rad_nue, F_radx_nue, F_rady_nue, F_radz_nue, E_rad_nux, F_radx_nux, F_rady_nux, F_radz_nux,&
                            F_rad0, P_radxx, P_radyy, P_radzz, P_radxy, P_radxz, P_radyz, & 
                            F_rad0_nue, P_radxx_nue, P_radyy_nue, P_radzz_nue, P_radxy_nue, P_radxz_nue, P_radyz_nue, &
                            F_rad0_nux, P_radxx_nux, P_radyy_nux, P_radzz_nux, P_radxy_nux, P_radxz_nux, P_radyz_nux, &
			    eps_tot, eps_cld, eps_thermal, P_cld, &
                            w,h,u0,Ex,Ey,Ez,sbt,sbx,sby,sbz, T_fluid,&
                            rho,Sx,Sy,Sz, &
                            Sxx,Sxy,Sxz,Syy,Syz,Szz, &
                            phi,lapm1,shiftx,shifty,shiftz, &
                            gxx,gxy,gxz,gyy,gyz,gzz, &
                            gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                            Symmetry,Matter_BC, &
                            neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th,rho_b_atm, rad_const, &
                            thisproc_have_global_bdry_min,thisproc_have_global_bdry_max,compute_microphysics,enable_OS_collapse,rad_evolve_enable,rad_closure_scheme,microphysics_scheme,Erad_atm_cut)
                    end do
	   else
              !FIXME: this BC does not support arbitrary number of ghostzones... yet.
              call matter_bc_mhd_movbh(ext,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1),mhd_bigP, &
                   PhysicalRadius,RadiusDerivative,X,Y,Z, &
                   rho_star,tau,mhd_st_x,mhd_st_y,mhd_st_z, &
		   tau_rad, S_rad_x, S_rad_y,S_rad_z, rhoYe,&
                   st_x,st_y,st_z, &
                   rho_b,P,vx,vy,vz,Bx,By,Bz, &
    		   E_rad, F_radx, F_rady, F_radz, Y_e, &
     		   F_rad0, P_radxx, P_radyy, P_radzz, P_radxy, P_radxz, P_radyz, &
     		   eps_tot, eps_cld, eps_thermal, P_cld, &
                   w,h,u0,Ex,Ey,Ez,sbt,sbx,sby,sbz, T_fluid,&
                   rho,Sx,Sy,Sz, &
                   Sxx,Sxy,Sxz,Syy,Syz,Szz, &
                   phi,lapm1,shiftx,shifty,shiftz, &
                   gxx,gxy,gxz,gyy,gyz,gzz, &
                   gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                   Symmetry,Matter_BC, &
                   neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th,rho_b_atm, rad_const,&
                   nrhovec,dxvec,rhovec,Pvec,vvec, &
                   thisproc_have_global_bdry_min,thisproc_have_global_bdry_max, compute_microphysics,enable_OS_collapse,rad_evolve_enable, rad_closure_scheme) 
           end if
        else if(use_disk_bcs==1) then
           !FIXME: this BC does not support arbitrary number of ghostzones... yet.
           call matter_bc_mhd_disk(ext,Z, &
                rho_star,tau,mhd_st_x,mhd_st_y,mhd_st_z, &
                st_x,st_y,st_z, &
                rho_b,P,vx,vy,vz,Bx,By,Bz, &
                w,h,u0,Ex,Ey,Ez,sbt,sbx,sby,sbz, T_fluid,&
                rho,Sx,Sy,Sz, &
                Sxx,Sxy,Sxz,Syy,Syz,Szz, &
                phi,lapm1,shiftx,shifty,shiftz, &
                gxx,gxy,gxz,gyy,gyz,gzz, &
                gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                Symmetry,Matter_BC, &
                neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th,rho_b_atm_gf, &
                pfloor_gf, enable_HARM_energyvariable, &
                thisproc_have_global_bdry_min,thisproc_have_global_bdry_max, enable_OS_collapse,rad_closure_scheme)
        end if
     end if





     !-----------------------------------------------------------------------------------------

     !====================================================
     ! Get the matter sources on the tilde level!
     !====================================================
!!$     call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_vs')
     !     call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_conservatives')
!!$     call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_primitives')

     ! Following line is needed since we ignore ghostzones inside the primitives solver, so, e.g. Sy(i,j,1) is incorrect:
     call CartSymGN(dummy,cctkGH,'bssn::BSSN_matter')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_primitives')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vs')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_primitives')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_pressure')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::microphys_primitives')
!     call CartSymGN(dummy,cctkGH,'mhd_evolve::temperatures')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_conservatives')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_conservatives')
  end if
!  print *, "^^^^^^^^^^^^^^^^^^^^^^ Done with mhd_primitives_and_bcs_v2 ^^^^^^^^^^^^^^^^^^^^^^^^"

end subroutine mhd_primitives_and_bcs_v2


