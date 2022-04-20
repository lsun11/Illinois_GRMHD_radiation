!-------------------------------------------------
! Stuff to run right after initial data is set up
!-------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine mhd_PostInitialData(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  real*8                :: dX,dY,dZ,b2bt, xs, r_bound
  integer               :: dummy,vindex,handle,ierr,i,j,k,cell_centering_enabled

  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  write(*,*) "Start mhd_postid!!!!!!!"
  write(*,*) "check: compute_microphysics, enable_OS_collapse, rad_opacity_abs = ", compute_microphysics, enable_OS_collapse, rad_opacity_abs 

  if(really_want_gamma_th_of_1_point_3==0 .and. abs(gamma_th-1.3).lt.1e-9) then
     write(*,*) "ERROR.  Default gamma_th (=1.3) detected.  You must set gamma_th to the appropriate value in your initial data thorn, or your .par file!"
     write(*,*) ""
     write(*,*) "If you really want to use gamma_th==1.3, please set mhd_evolve::really_want_gamma_th_of_1_point_3=1 in your .par file!"
     stop
  end if

  if(Symmetry.ne.4 .and. Symmetry.ne.0) then
     call cell_centered_or_not(cctk_lsh,X,Y,Z,cell_centering_enabled)
  else
     cell_centering_enabled=1
  end if

  !For emfields, we assume that you've set Bx, By, Bz (the UN-tilded B^i's)
  ! or Ax, Ay, Az (if using constrained transport scheme of Del Zanna)

  write(*,*) "------------------inside mhd_postinitialdata.F90-------------------"

  ! Do CartSymGN's on all the MHD quantities:
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_conservatives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativex')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativey')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativez')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::em_Blagrangemultiplier')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_primitives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vs')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_conservatives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_primitives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_pressure')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::micphys_conservatives') 
  call CartSymGN(dummy,cctkGH,'mhd_evolve::microphys_primitives')
!  call CartSymGN(dummy,cctkGH,'mhd_evolve::temperatures')

  if (constrained_transport_scheme .ne. 3) then 
     !$omp parallel do
     do k=1,ext(3)
	do j=1,ext(2)
	   do i=1,ext(1)
              Bxtilde(i,j,k) = Bx(i,j,k) 
              Bytilde(i,j,k) = By(i,j,k)
              Bztilde(i,j,k) = Bz(i,j,k)
	   end do
        end do
     end do
     !$omp end parallel do

     ! Here, we convert B^i to tilde B^i
     b2bt = 1.D0
     call convert_b(ext,Bxtilde,Bytilde,Bztilde,phi,b2bt)
  end if

  !------------------------------------------------------------------
  ! FILL _p AND _p_p TIMELEVELS:

!  write(*,*) "------------------test 3-------------------"


  !$omp parallel do
  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
           psi6phi_p(i,j,k) = psi6phi(i,j,k)
           Ax_p(i,j,k) = Ax(i,j,k) 
           Ay_p(i,j,k) = Ay(i,j,k)
           Az_p(i,j,k) = Az(i,j,k)
	   psi6phi_p_p(i,j,k) = psi6phi(i,j,k)
           Ax_p_p(i,j,k) = Ax(i,j,k) 
           Ay_p_p(i,j,k) = Ay(i,j,k)
           Az_p_p(i,j,k) = Az(i,j,k)

	   Bxtilde_p(i,j,k) = Bxtilde(i,j,k) 
           Bytilde_p(i,j,k) = Bytilde(i,j,k)
           Bztilde_p(i,j,k) = Bztilde(i,j,k)
           Bxtilde_p_p(i,j,k) = Bxtilde(i,j,k) 
           Bytilde_p_p(i,j,k) = Bytilde(i,j,k)
           Bztilde_p_p(i,j,k) = Bztilde(i,j,k)

           if (isnan(Ax(i,j,k)) .or. isnan(Ay(i,j,k)) .or. isnan(Az(i,j,k))) then
              write (*,*) "Inside mhd_postinitialdata.F90, Ax(i,j,k) is nan!!!!"
           end if

        end do
     end do
  end do
  !$omp end parallel do

  !mhd_evolve thorn assumes that you've set the conservatives & primitives

  !$omp parallel do
  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
           rho_star_p(i,j,k) = rho_star(i,j,k)
           rhoYe_p(i,j,k) = rhoYe(i,j,k)
           tau_p(i,j,k) = tau(i,j,k)
           mhd_st_x_p(i,j,k) = mhd_st_x(i,j,k) 
           mhd_st_y_p(i,j,k) = mhd_st_y(i,j,k)
           mhd_st_z_p(i,j,k) = mhd_st_z(i,j,k)
	   tau_rad_p(i,j,k) = tau_rad(i,j,k)
	   S_rad_x_p(i,j,k) = S_rad_x(i,j,k)
	   S_rad_y_p(i,j,k) = S_rad_y(i,j,k)
	   S_rad_z_p(i,j,k) = S_rad_z(i,j,k)
           tau_rad_nue_p(i,j,k) = tau_rad_nue(i,j,k)
           S_rad_x_nue_p(i,j,k) = S_rad_x_nue(i,j,k)
           S_rad_y_nue_p(i,j,k) = S_rad_y_nue(i,j,k)
           S_rad_z_nue_p(i,j,k) = S_rad_z_nue(i,j,k)
           tau_rad_nux_p(i,j,k) = tau_rad_nux(i,j,k)
           S_rad_x_nux_p(i,j,k) = S_rad_x_nux(i,j,k)
           S_rad_y_nux_p(i,j,k) = S_rad_y_nux(i,j,k)
           S_rad_z_nux_p(i,j,k) = S_rad_z_nux(i,j,k)

           rho_star_p_p(i,j,k) = rho_star(i,j,k)
           rhoYe_p_p(i,j,k) = rhoYe(i,j,k)
           tau_p_p(i,j,k) = tau(i,j,k)
           mhd_st_x_p_p(i,j,k) = mhd_st_x(i,j,k) 
           mhd_st_y_p_p(i,j,k) = mhd_st_y(i,j,k)
           mhd_st_z_p_p(i,j,k) = mhd_st_z(i,j,k)
	   tau_rad_p_p(i,j,k) = tau_rad(i,j,k)
           S_rad_x_p_p(i,j,k) = S_rad_x(i,j,k)
           S_rad_y_p_p(i,j,k) = S_rad_y(i,j,k)
           S_rad_z_p_p(i,j,k) = S_rad_z(i,j,k)
           tau_rad_nue_p_p(i,j,k) = tau_rad_nue(i,j,k)
           S_rad_x_nue_p_p(i,j,k) = S_rad_x_nue(i,j,k)
           S_rad_y_nue_p_p(i,j,k) = S_rad_y_nue(i,j,k)
           S_rad_z_nue_p_p(i,j,k) = S_rad_z_nue(i,j,k)
           tau_rad_nux_p_p(i,j,k) = tau_rad_nux(i,j,k)
           S_rad_x_nux_p_p(i,j,k) = S_rad_x_nux(i,j,k)
           S_rad_y_nux_p_p(i,j,k) = S_rad_y_nux(i,j,k)
           S_rad_z_nux_p_p(i,j,k) = S_rad_z_nux(i,j,k)
        end do
     end do
  end do
  !$omp end parallel do
  !-----------------------------------------
  !Set aux. MHD gridfunctions to zero, so that valgrind is happy.
  !Should probably comment these out if you're debugging the code
  !$omp parallel do
  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
           gupxx_f(i,j,k) = 0.D0
           gupyy_f(i,j,k) = 0.D0
           gupzz_f(i,j,k) = 0.D0 
           cmax(i,j,k) = 0.D0
           cmin(i,j,k) = 0.D0
           cmax_rad(i,j,k) = 0.D0
           cmin_rad(i,j,k) = 0.D0
           cmax_rad_nue(i,j,k) = 0.D0
           cmin_rad_nue(i,j,k) = 0.D0
           cmax_rad_nux(i,j,k) = 0.D0
           cmin_rad_nux(i,j,k) = 0.D0

           lapm1_f(i,j,k) = 0.D0
           shiftx_f(i,j,k) = 0.D0
           shifty_f(i,j,k) = 0.D0
           shiftz_f(i,j,k) = 0.D0 
           gxx_f(i,j,k) = 0.D0
           gxy_f(i,j,k) = 0.D0
           gxz_f(i,j,k) = 0.D0
           gyy_f(i,j,k) = 0.D0
           gyz_f(i,j,k) = 0.D0 
           gzz_f(i,j,k) = 0.D0 
           phi_f(i,j,k) = 1.D0

           dP_m(i,j,k) = 0.D0
           Pr(i,j,k) = 0.D0
           Pl(i,j,k) = 0.D0
           drho_b_m(i,j,k) = 0.D0
           rho_br(i,j,k) = 0.D0
           rho_bl(i,j,k) = 0.D0
           dvx_m(i,j,k) = 0.D0
           dvy_m(i,j,k) = 0.D0
           dvz_m(i,j,k) = 0.D0
           vxr(i,j,k) = 0.D0
           vxl(i,j,k) = 0.D0
           vyr(i,j,k) = 0.D0
           vyl(i,j,k) = 0.D0
           vzr(i,j,k) = 0.D0
           vzl(i,j,k) = 0.D0
           Bxr(i,j,k) = 0.D0
           Bxl(i,j,k) = 0.D0
           Byr(i,j,k) = 0.D0
           Byl(i,j,k) = 0.D0
           Bzr(i,j,k) = 0.D0
           Bzl(i,j,k) = 0.D0
           v02r(i,j,k) = 0.D0
           v02l(i,j,k) = 0.D0
           v02_radr(i,j,k) = 0.D0
           v02_radl(i,j,k) = 0.D0
           v02_rad_nuer(i,j,k) = 0.D0
           v02_rad_nuel(i,j,k) = 0.D0
           v02_rad_nuxr(i,j,k) = 0.D0
           v02_rad_nuxl(i,j,k) = 0.D0
           dBx_m(i,j,k) = 0.D0
           dBy_m(i,j,k) = 0.D0
           dBz_m(i,j,k) = 0.D0
           Y_el(i,j,k) = 0.D0
           Y_er(i,j,k) = 0.D0
	   T_fluidle(i,j,k) = 0.D0
           T_fluidr(i,j,k) = 0.D0

	   E_radr(i,j,k) = 0.D0
	   E_radl(i,j,k) = 0.D0
	   F_radxr(i,j,k) = 0.D0
           F_radxle(i,j,k) = 0.D0
           F_radyr(i,j,k) = 0.D0
           F_radyle(i,j,k) = 0.D0
	   F_radzr(i,j,k) = 0.D0 
           F_radzle(i,j,k) = 0.D0
           FaFar(i,j,k) = 0.D0
           FaFal(i,j,k) = 0.D0
           E_rad_nuer(i,j,k) = 0.D0
           E_rad_nuel(i,j,k) = 0.D0
           F_radx_nuer(i,j,k) = 0.D0
           F_radx_nuele(i,j,k) = 0.D0
           F_rady_nuer(i,j,k) = 0.D0
           F_rady_nuele(i,j,k) = 0.D0
           F_radz_nuer(i,j,k) = 0.D0
           F_radz_nuele(i,j,k) = 0.D0
           FaFar_nue(i,j,k) = 0.D0
           FaFal_nue(i,j,k) = 0.D0
           E_rad_nuxr(i,j,k) = 0.D0
           E_rad_nuxl(i,j,k) = 0.D0
           F_radx_nuxr(i,j,k) = 0.D0
           F_radx_nuxle(i,j,k) = 0.D0
           F_rady_nuxr(i,j,k) = 0.D0
           F_rady_nuxle(i,j,k) = 0.D0
           F_radz_nuxr(i,j,k) = 0.D0
           F_radz_nuxle(i,j,k) = 0.D0
           FaFar_nux(i,j,k) = 0.D0
           FaFal_nux(i,j,k) = 0.D0
!           eps_thermal(i,j,k) = 0.D0
!           eps_cld(i,j,k) = 0.D0

           rho_star_rhs(i,j,k) = 0.D0
           rhoYe_rhs(i,j,k) = 0.D0
           tau_rhs(i,j,k) = 0.D0 
           mhd_st_x_rhs(i,j,k) = 0.D0 
           mhd_st_y_rhs(i,j,k) = 0.D0
           mhd_st_z_rhs(i,j,k) = 0.D0 
	   
	   tau_rad_rhs(i,j,k) = 0.D0
	   S_rad_x_rhs(i,j,k) = 0.D0
	   S_rad_y_rhs(i,j,k) = 0.D0
	   S_rad_z_rhs(i,j,k) = 0.D0
           tau_rad_nue_rhs(i,j,k) = 0.D0
           S_rad_x_nue_rhs(i,j,k) = 0.D0
           S_rad_y_nue_rhs(i,j,k) = 0.D0
           S_rad_z_nue_rhs(i,j,k) = 0.D0
           tau_rad_nux_rhs(i,j,k) = 0.D0
           S_rad_x_nux_rhs(i,j,k) = 0.D0
           S_rad_y_nux_rhs(i,j,k) = 0.D0
           S_rad_z_nux_rhs(i,j,k) = 0.D0

	   tau_rad_flux(i,j,k) = 0.D0
	   S_radx_flux(i,j,k) = 0.D0
	   S_rady_flux(i,j,k) = 0.D0
	   S_radz_flux(i,j,k) = 0.D0
           tau_rad_nue_flux(i,j,k) = 0.D0
           S_radx_nue_flux(i,j,k) = 0.D0
           S_rady_nue_flux(i,j,k) = 0.D0
           S_radz_nue_flux(i,j,k) = 0.D0
           tau_rad_nux_flux(i,j,k) = 0.D0
           S_radx_nux_flux(i,j,k) = 0.D0
           S_rady_nux_flux(i,j,k) = 0.D0
           S_radz_nux_flux(i,j,k) = 0.D0
           drhoYe_m(i,j,k) = 0.D0
	   dT_fluid_m(i,j,k) = 0.D0

           Blagrangemultiplier_rhs(i,j,k) = 0.D0 
           Blagrangemultiplier(i,j,k) = 0.D0
           Blagrangemultiplier_p(i,j,k) = 0.D0
           Blagrangemultiplier_p_p(i,j,k) = 0.D0
           rhob_floor(i,j,k) = 0.D0
           P_floor(i,j,k) = 0.D0

!           T_fluid(i,j,k) = 0.D0
        end do
     end do
  end do
  !$omp end parallel do

  ! The following lines are needed for the mhd_post_regrid that occurs right after timestep zero.
  !    Otherwise, we get valgrind errors when using the primitive_vars_hybrid2() primitives solver.
  !$omp parallel do
  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)  
           h_p(i,j,k) = h(i,j,k)
        end do
     end do
  end do
  !$omp end parallel do


  !Here we compute del dot B all over the grid and store result in Pl:
  if(constrained_transport_scheme==2) then
     !here we use vxr,vxl,vyl,vyr, and vzl as temporary storage, and store the measured monopole violation inside Pl
     call compute_absdivb2(ext,dX,dY,dZ,X,Z,Bx,By,Bz,phi,Pl,Sym_Bz,Symmetry,vxr,vxl,vyl,vyr,vzl)
  elseif (constrained_transport_scheme==3) then 
     ! Note that Bitilde store B^i on staggered grid
     call compute_absdivb_ct3(ext,dX,dY,dZ,phi,Bx_stagger,By_stagger,Bz_stagger, & 
          temp1,temp2,temp3,Pl)
  else
     !here we use vxr,vxl,vyl,vyr, and vzl as temporary storage, and store the measured monopole violation inside Pl
     call compute_absdivb(ext,dX,dY,dZ,X,Z,Bx,By,Bz,phi,Pl,Sym_Bz,Symmetry,vxr,vxl,vyl,vyr,vzl)
  end if
  if(Symmetry==4) then
     !$omp parallel do
     do k=1,ext(3)
        do j=1,ext(2)
           Pl(1,j,k) = Pl(2,j,k)
        end do
     end do
     !$omp end parallel do
     !$omp parallel do
     do j=1,ext(2)
        do i=1,ext(1)
           Pl(i,j,1) = Pl(i,j,2)
        end do
     end do
     !$omp end parallel do
  else if(Symmetry==1) then
     !FIXME: incorrect for arbitrary number of ghostzones!
     !Also incorrect if local processor does not complain z<0!
     !$omp parallel do
     do j=1,ext(2)
        do i=1,ext(1)
           Pl(i,j,1) = Pl(i,j,2)
        end do
     end do
     !$omp end parallel do
  end if


    do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
           if(isnan(tau_rad(i,j,k)) .or. isnan(S_rad_x(i,j,k))) then
              write(*,*) "Inside mhd_postinitialdata.F90, before metric_source_terms_and_misc_vars_initial, tau_rad, S_rad_x, E_rad, F_rad_x=", tau_rad, S_rad_x, E_rad, F_radx
           end if
        end do
     end do
  end do

  if(zero_out_matter_source_terms==0) then
     call metric_source_terms_and_misc_vars_initial(cctkGH,cctk_lsh, &
          rho,Sx,Sy,Sz, &
          Sxx,Sxy,Sxz,Syy,Syz,Szz, &
	  tau_rad, S_rad_x,S_rad_y, S_rad_z,&
          E_rad, F_radx, F_rady, F_radz, F_rad0, F_rad_scalar,&
          h,w,&
	  st_x,st_y,st_z, &
          Ex,Ey,Ez, &
          sbt,sbx,sby,sbz, &
          lapm1,shiftx,shifty,shiftz, &
          phi,gxx,gxy,gxz,gyy,gyz,gzz, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
          P,rho_b,u0,vx,vy,vz, &
          Bx,By,Bz, &
          rho_star,mhd_st_x,mhd_st_y,mhd_st_z, &
	  neos,ergo_star, ergo_sigma,&
          rho_tab, P_tab,eps_tab,k_tab,gamma_tab,&
          rho_b_atm, enable_OS_collapse, rad_evolve_enable, rad_closure_scheme, Psi6threshold, Erad_atm_cut)
  end if



  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
           if(isnan(tau_rad(i,j,k)).or.isnan(S_rad_x(i,j,k))) then
              write(*,*) "In mhd_postinitialdata.F90,  tau_rad and S_rad_x=", tau_rad(i,j,k), S_rad_x(i,j,k)
              write(*,*) "X,Y,Z=", X(i,j,k), Y(i,j,k), Z(i,j,k)
           end if
        end do
     end do
  end do
  

 write(*,*) "------------------end mhd_PostInitialData-------------------"

end subroutine mhd_PostInitialData
