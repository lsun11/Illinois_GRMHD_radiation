!--------------------------------------------------------------
! Driver for primitives and boundary conditions routines, v2.0
!--------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
subroutine mhd_primitives_postmolstep(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext,total_nghostzones
  integer, dimension(3) :: thisproc_have_global_bdry_min,thisproc_have_global_bdry_max
  real*8                :: dT,dX,dY,dZ,b2bt
  real*8                :: rho_fail_max_step,M_fail_step,randomnumber
  real*8		:: xNS1,yNS1, xNS2, yNS2
  integer               :: proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax
  integer               :: glob_imax,glob_jmax,glob_kmax
  integer               :: i,j,k,vindex,NFont,NFont_L
  integer               :: handle,dummy,ierr
  integer               :: repairs_needed,repairs_rad_needed,ignore_ghostzones
  CCTK_REAL             :: reduction_value
  integer               :: AXISYM, levelnumber
  real*8  :: alpn1, psiL, psiL4, psiL6, P_radl, u0l, T_fluid_cgs
  real*8  :: P_rad0xl, P_rad0yl, P_rad0zl, P_rad00l
  real*8  :: F_rad_xl, F_rad_yl, F_rad_zl, F_rad_0l, temp_rad, temp_rad1
  real*8  :: beta2, udotbeta, g_00l, u_0l
  real*8  :: u_x, u_y, u_z, uxl, uyl, uzl, v_xl, v_yl, v_zl
  real*8  :: shift_x, shift_y, shift_z
  real*8  :: Fasq, zeta_temp, zeta_cut, zeta, chil, eps_l
  real*8  :: rad_const_cgs, G_cgs, c_cgs, kb_cgs, c2k, k2c
  real*8  :: er, el, au0r1, rho_sl, wl
  real*8  :: B_x, B_y, B_z, B2, sbtl, sb2, sb_x, sb_y, sb_z
  real*8  :: fac, E_x, E_y, E_z, Exl, Eyl, Ezl, temp_a, f1o4pa, f1o4p, f1o8p
  real*8, parameter        :: fac2 = 0.99d0

  parameter(AXISYM = 4)

  f1o4p = 1.d0/(4.d0*acos(-1.d0))
  f1o8p = 1.d0/(8.d0*acos(-1.d0))

  ext = cctk_lsh


  rad_const_cgs = 7.5646e-15
  G_cgs = 6.6743015e-8
  c_cgs = 29979245800.0
  kb_cgs = 1.380649e-16 
  c2k = 1.0e-5 
  k2c = 1.0e5


  print *, "^^^^^^^^^^^^^^^^^^^^^^ Start driver_primitives_postmolstep ^^^^^^^^^^^^^^^^^^^^^^^^"
  print *, " iteration_to_insert_rad = ", iteration_to_insert_rad
!  print *, "primitives_solver = ", primitives_solver
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



  if(levelnumber .gt. 1 .and. cctk_iteration.gt.0) then
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

     !---------------------------------------
     !Recompute B^i from \tilde{B^i} or compute B^i from Ai
     !(This is the easiest primitives recovery!)
     ! Note that when constrained_transport_scheme==3, B has 
     ! already been computed.
     if (em_evolve_enable==1 .and. constrained_transport_scheme .ne. 3) then
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

     !-----------------------------------------------------------------------------------------
     ! PRIMITIVE RECONSTRUCTION STEP:
     call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativex')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativey')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativez')
     !     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_conservatives')

     if (enable_primitives_disk == 0) then

        if(primitives_solver==11) then
           !Here we use the HARM primitives solver, with a new prescription that minimizes changes to
           !  conservatives without applying the "Font" fix.
           !We intend to make the below function truly generic, but currently it only supports Gamma=2.
           !  It can be trivially extended for arbitrary Gamma-law EOS's, but will require work for more
           !  generic EOS's.
           !We hope the below function will eventually be used in place of other primitives solver options,
           !  since it can be easily extended to use your desired technique for solving primitives.
           
	   !write(6,*) "(1)Inside driver_primitives_postmolstep, test gamma_tab =", gamma_tab

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
                rho_b_atm,tau_atm, rho_b_max,&
		xNS1,yNS1, xNS2, yNS2, M_B, rad_T_fac, rad_T_cutoff, rad_T_pow, rad_T_floor, T_fluid,&
                neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
                temp1,temp2,temp3,temp4,temp5, &
                primitives_debug,Psi6threshold,rad_const,horizon_enforce_rho_profile,enable_OS_collapse, &
		rad_evolve_enable, compute_microphysics, &
                microphysics_scheme, T_fluid_cgs_atm, rad_fix)

           !The goal of the below function is to just update the metric source terms.  We should
           !  really get rid of the old h,w,st_i,Ei gridfunctions...
           !Currently this only works for Gamma-law EOS's but could be extended for arbitary EOSs
           !  with some work.
!	   write(*,*) "inside driver_primitives_postmolstep.f90"
           call metric_source_terms_and_misc_vars(cctkGH,cctk_lsh, &
                rho,Sx,Sy,Sz, &
                Sxx,Sxy,Sxz,Syy,Syz,Szz, &
        	tau_rad, S_rad_x,S_rad_y, S_rad_z,&
                rhoYe,Y_e,&
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


  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
           if (T_fluid(i,j,k) .gt. 1.0e-20) then
              write(*,*) " In driver_primitives_postmolsetp, T_fluid too large!!!!"
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
           if(cctk_iteration .gt. iteration_to_insert_rad .and. rad_evolve_enable==1 .and. rad_closure_scheme == 1) then
              call primitive_vars_rad_cpp(ext,cctk_nghostzones, X, Y, Z, &
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

           

           if(microphysics_scheme==1) then
	   write(*,*) "Inside driver_primitives_postmolstep.F90, microphysics_scheme,rad_evolve_enable,rad_closure_scheme=",microphysics_scheme,rad_evolve_enable,rad_closure_scheme
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

           total_nghostzones = number_of_mol_ministeps*cctk_nghostzones
           call hydro_primitives(cctkGH,ext,total_nghostzones, X, Y, Z, rho_star, tau,&
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
                gamma_th,K_poly,sdots_o_rhot,Symmetry,-1)

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

           total_nghostzones = number_of_mol_ministeps*cctk_nghostzones
           !TODO: remove temp1 (unused) parameter.
           call primitive_vars_hybrid2(ext,cctk_nghostzones,X,Y,Z, &
                rho_star,tau,st_x,st_y,st_z, &
                mhd_st_x,mhd_st_y,mhd_st_z,neos, ergo_star, ergo_sigma,&
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
                excision_zone_gf, tau_stildefix_enable,tau_atm,-1)

        else if(primitives_solver==2) then
           print *, "^^^^^^^^^^^^^^^^^^^^^^ primitives_solver==2 ^^^^^^^^^^^^^^^^^^^^^^"
           !primitive solver (next function call) overwrites h_p.  We want to keep _p values.
           !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
		    h_p(i,j,k) = h(i,j,k)      
                    if(isnan(tau_rad(i,j,k)).or.isnan(S_rad_x(i,j,k))) then
                       write(*,*) "In driver_primitives_postmolstep, before primitive_vars_hybrid2_cpp, tau_rad and S_rad_x=", tau_rad(i,j,k), S_rad_x(i,j,k)
                       write(*,*) "X,Y,Z=", X(i,j,k), Y(i,j,k), Z(i,j,k)
                    end if
                 end do
              end do
          end do

          ignore_ghostzones = 0
          total_nghostzones = number_of_mol_ministeps*cctk_nghostzones
          
           !TODO: remove temp1 (unused) parameter.
	   print *, "^^^^^^^^^^^^^^^^^^^^^^ Inside driver_primitive_postmolstep.F90 start primitive_vars_hybrid2_cpp ^^^^^^^^^^^^^^^^^^^^^^^^"
           call primitive_vars_hybrid2_cpp(ext,cctk_nghostzones,X,Y,Z, &
                rho_star,tau,st_x,st_y,st_z, &
		mhd_st_x,mhd_st_y,mhd_st_z, rhoYe,&
                tau_rad,S_rad_x,S_rad_y,S_rad_z,&
                neos, ergo_star, ergo_sigma,&
                rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, &
                w,temp1,rho_b,rho,P,h,Sx,Sy,Sz, &
                Sxx,Sxy,Sxz,Syy,Syz,Szz, &
                E_rad,F_radx,F_rady,F_radz,F_rad0, F_rad_scalar, Y_e, eps_thermal,&
                phi,lapm1,shiftx,shifty,shiftz, &
                gxx,gxy,gxz,gyy,gyz,gzz, &
                gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                h_p,u0,rho_b_max,rho_b_atm, &
                rho_fail_max_step,M_fail_step,rhos_max, &
                Bx,By,Bz,Ex,Ey,Ez, &
                vx,vy,vz, &
                sbt,sbx,sby,sbz, &
         	xNS1,yNS1, xNS2, yNS2, M_B, rad_T_fac, rad_T_cutoff, rad_T_pow, rad_T_floor, T_fluid,&
                proc_imin,proc_jmin,proc_kmin, &
                proc_imax,proc_jmax,proc_kmax, &
                glob_imax,glob_jmax,glob_kmax, &
                Symmetry,pfloor,excision_enable, &
                excision_zone_gf, tau_stildefix_enable,tau_atm,temp4,cctkGH,ignore_ghostzones, &
                enable_shocktest_primitive_mode,repairs_needed, rad_closure_scheme, &
                rad_const,enable_OS_collapse, compute_microphysics, Psi6threshold, Erad_atm_cut)

           if(cctk_iteration .gt. iteration_to_insert_rad .and. rad_evolve_enable==1 .and. rad_closure_scheme == 1) then
              call primitive_vars_rad_cpp(ext,cctk_nghostzones, X, Y, Z, &
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
                   neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
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
           call CartSymGN(dummy,cctkGH,'mhd_evolve::em_Blagrangemultiplier')
           call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_primitives')
           call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vs')
           call CartSymGN(dummy,cctkGH,'bssn::BSSN_matter')
           call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_primitives')
           call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_pressure')
           call CartSymGN(dummy,cctkGH,'mhd_evolve::microphys_primitives')
!           call CartSymGN(dummy,cctkGH,'mhd_evolve::temperatures')
        end if
     else if(enable_primitives_disk==1 .and. use_harm_primitives==0) then

        call primitive_vars_alt_disk(ext,X,Y,Z, &
             rho_star,tau,st_x,st_y,st_z, &
             mhd_st_x,mhd_st_y,mhd_st_z,neos, ergo_star, ergo_sigma,&
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

  end if

write(*,*) "end subroutine mhd_primitives_postmolstep"
end subroutine mhd_primitives_postmolstep
