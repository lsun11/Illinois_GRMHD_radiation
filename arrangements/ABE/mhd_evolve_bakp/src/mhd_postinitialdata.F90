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
  real*8                :: dX,dY,dZ,b2bt
  integer               :: dummy,vindex,handle,ierr,i,j,k,cell_centering_enabled

  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  write(*,*) "Start mhd_postid!!!!!!!"

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

!  write(*,*) "------------------test 1-------------------"

  call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_conservatives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::rad_primitives')
  
!  write(*,*) "------------------test 2-------------------"

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

!	   E_rad_p(i,j,k) = E_rad(i,j,k)
!	   F_radx_p(i,j,k) = F_radx(i,j,k)
!	   F_rady_p(i,j,k) = F_rady(i,j,k)
!	   F_radz_p(i,j,k) = F_radz(i,j,k)
           
!           E_rad_p_p(i,j,k) = E_rad(i,j,k)
!           F_radx_p_p(i,j,k) = F_radx(i,j,k)
!           F_rady_p_p(i,j,k) = F_rady(i,j,k)
!           F_radz_p_p(i,j,k) = F_radz(i,j,k)

	   Bxtilde_p(i,j,k) = Bxtilde(i,j,k) 
           Bytilde_p(i,j,k) = Bytilde(i,j,k)
           Bztilde_p(i,j,k) = Bztilde(i,j,k)
           Bxtilde_p_p(i,j,k) = Bxtilde(i,j,k) 
           Bytilde_p_p(i,j,k) = Bytilde(i,j,k)
           Bztilde_p_p(i,j,k) = Bztilde(i,j,k)
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
           tau_p(i,j,k) = tau(i,j,k)
           mhd_st_x_p(i,j,k) = mhd_st_x(i,j,k) 
           mhd_st_y_p(i,j,k) = mhd_st_y(i,j,k)
           mhd_st_z_p(i,j,k) = mhd_st_z(i,j,k)
	   tau_rad_p(i,j,k) = tau_rad(i,j,k)
	   S_rad_x_p(i,j,k) = S_rad_x(i,j,k)
	   S_rad_y_p(i,j,k) = S_rad_y(i,j,k)
	   S_rad_z_p(i,j,k) = S_rad_z(i,j,k)

           rho_star_p_p(i,j,k) = rho_star(i,j,k)
           tau_p_p(i,j,k) = tau(i,j,k)
           mhd_st_x_p_p(i,j,k) = mhd_st_x(i,j,k) 
           mhd_st_y_p_p(i,j,k) = mhd_st_y(i,j,k)
           mhd_st_z_p_p(i,j,k) = mhd_st_z(i,j,k)
	   tau_rad_p_p(i,j,k) = tau_rad(i,j,k)
           S_rad_x_p_p(i,j,k) = S_rad_x(i,j,k)
           S_rad_y_p_p(i,j,k) = S_rad_y(i,j,k)
           S_rad_z_p_p(i,j,k) = S_rad_z(i,j,k)
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
           dBx_m(i,j,k) = 0.D0
           dBy_m(i,j,k) = 0.D0
           dBz_m(i,j,k) = 0.D0

	   E_radr(i,j,k) = 0.D0
	   E_radl(i,j,k) = 0.D0
	   F_rad0r(i,j,k) = 0.D0
	   F_rad0l(i,j,k) = 0.D0
	   F_radxr(i,j,k) = 0.D0
           F_radxl(i,j,k) = 0.D0
           F_radyr(i,j,k) = 0.D0
           F_radyl(i,j,k) = 0.D0
	   F_radzr(i,j,k) = 0.D0 
           F_radzl(i,j,k) = 0.D0
	   P_radr(i,j,k) = 0.D0
           P_radl(i,j,k) = 0.D0

           rho_star_rhs(i,j,k) = 0.D0
           tau_rhs(i,j,k) = 0.D0 
           mhd_st_x_rhs(i,j,k) = 0.D0 
           mhd_st_y_rhs(i,j,k) = 0.D0
           mhd_st_z_rhs(i,j,k) = 0.D0 
	   
	   tau_rad_rhs(i,j,k) = 0.D0
	   S_rad_x_rhs(i,j,k) = 0.D0
	   S_rad_y_rhs(i,j,k) = 0.D0
	   S_rad_z_rhs(i,j,k) = 0.D0


           Blagrangemultiplier_rhs(i,j,k) = 0.D0 
           Blagrangemultiplier(i,j,k) = 0.D0
           Blagrangemultiplier_p(i,j,k) = 0.D0
           Blagrangemultiplier_p_p(i,j,k) = 0.D0
           rhob_floor(i,j,k) = 0.D0
           P_floor(i,j,k) = 0.D0
        end do
     end do
  end do
  !$omp end parallel do


! write(*,*) "------------------test 4-------------------"

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

  !  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_conservatives')
  !  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_primitives')
  !  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_vs')

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


! write(*,*) "------------------test 5-------------------"

  if(zero_out_matter_source_terms==0) then
     call metric_source_terms_and_misc_vars(cctkGH,cctk_lsh, &
          rho,Sx,Sy,Sz, &
          Sxx,Sxy,Sxz,Syy,Syz,Szz, &
	  tau_rad, S_rad_x,S_rad_y, S_rad_z,&
          E_rad, F_rad0, F_radx, F_rady, F_radz, P_rad,&
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
          rho_b_atm)
  end if

 write(*,*) "------------------end mhd_PostInitialData-------------------"

end subroutine mhd_PostInitialData
