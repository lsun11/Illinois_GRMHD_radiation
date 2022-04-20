!------------------------------------------------------------
! Update all auxiliary MHD variables after an AMR regridding
!------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
subroutine mhd_post_regrid(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext,fake_ext,thisproc_have_global_bdry_min,thisproc_have_global_bdry_max
  real*8                :: dX,dY,dZ,levelnumber
  real*8                :: Xglobmin,Yglobmin,Zglobmin,Xglobmax,Yglobmax,Zglobmax
  real*8                :: rho_fail_max_step,M_fail_step,b2bt
  integer               :: proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax
  integer               :: glob_imax,glob_jmax,glob_kmax
  integer               :: i,j,k,vindex,NFont,NFont_L
  integer               :: handle,dummy,ierr,repairs_needed
  CCTK_REAL             :: reduction_value

  integer               :: AXISYM
  parameter(AXISYM = 4)

  ext = cctk_lsh
  levelnumber = cctk_levfac(1)
  levelnumber = log(levelnumber)/log(2.D0)+1.D0
  write(*,*) "At the beginning of driver_post_regrid. levelnumber is ", levelnumber
  write(*,*) "cctk_levfac(1) is ", cctk_levfac(1)
  write(*,*) "cctk_levfac is ", cctk_levfac 

  ! IMPORTANT NOTE:
  ! After an AMR regridding (i.e., AMR grid movement), Carpet sets only the evolved variables (e.g., gij's, Aij's) 
  !   at all points on all grids.
  ! However, auxiliary variables (such as, e.g., matter source terms, Rij, etc.) are set only in regions where the grid
  !   _was_, and not in new regions.  
  !
  ! To properly explain this, it's best to draw a 1D picture of a gridfunction defined along a line parallel to the x-axis.
  ! Let '+' denote the gridpoints where the gridfunction is defined, 'X' denote the gridpoints where the gridfunction is
  !   not defined, and '|' the grid boundary
  ! Before a regrid, here's what our grid looks like:
  ! |++++++++++++++++++++++++++|
  !  <---      x-axis      --->  
  ! After a regrid, the grid has moved a little to the right:
  !      |+++++++++++++++++++++XXXXX|
  !  <---      x-axis      --->
  !  Attempting to access a gridfunction at the 'X' points will result in memory errors, and ultimately the simulation may crash.
  ! This function is designed to fill in all the gridpoints after a regrid, including the 'X' points.

  if(cctk_iteration.ge.0) then
     dX = CCTK_DELTA_SPACE(1)
     dY = CCTK_DELTA_SPACE(2)
     dZ = CCTK_DELTA_SPACE(3)
     
     !write(*,*) "checkpoint (1). levelnumber is ", levelnumber

     ! Should check if there are NaNs first...
     !call fill_bssn_symmetry_gz_bssn_vars(ext,X,Y,Z,Symmetry,phi,chi,trK,Gammax,Gammay,Gammaz,gxx,gxy,gxz,gyy,gyz,gzz,Axx,Axy,Axz,Ayy,Ayz,Azz)
     !     call fill_shift_symmetry_gz(ext,X,Y,Z,Symmetry,shiftx,shifty,shiftz,shiftxt,shiftyt,shiftzt)

     ! FOLLOWING CARTSYM'S ARE NEEDED.  Otherwise, we'll end up with memory errors inside ghostzones (even BSSN quantities)!
     call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')
     call CartSymGN(dummy,cctkGH,'bssn::BSSN_vars')
     call CartSymGN(dummy,cctkGH,'shift::shift_vars')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_conservatives')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativex')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativey') 
     call CartSymGN(dummy,cctkGH,'mhd_evolve::em_conservativez')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::em_Blagrangemultiplier')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_primitives')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vs')


     !write(*,*) "checkpoint (2). levelnumber is ", levelnumber
     ! gupij is not evolved, and so is not defined anywhere that the grid has moved.
     !  Here we recompute gupij from the newly updated gij's:
     call BSSN_compute_gupij(cctkGH,cctk_lsh, &
          gxx,gxy,gxz,gyy,gyz,gzz, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)


     !write(*,*) "checkpoint (3). levelnumber is ", levelnumber
     if (constrained_transport_scheme .ne. 3) then
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

        ! Notice that we set B^i = tildeB^i above.
        b2bt = -1.D0
        call convert_b(ext,Bx,By,Bz,phi,b2bt)
     end if

     !write(*,*) "checkpoint (4). levelnumber is ", levelnumber
     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              !Set face values of field values to zero, so that valgrind is happy.
              !Should probably comment these out if you're debugging the code
              gupxx_f(i,j,k) = 0.D0
              gupyy_f(i,j,k) = 0.D0
              gupzz_f(i,j,k) = 0.D0 
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
           end do
        end do
     end do
     !$omp end parallel do

     !-----------------------------------------------------------------------------------------
     ! PRIMITIVE RECONSTRUCTION STEP:

     !write(*,*) "checkpoint (5). levelnumber is ", levelnumber
     glob_imax = ext(1)
     glob_jmax = ext(2)
     glob_kmax = ext(3)
     proc_imin = -100
     proc_jmin = -100
     proc_kmin = -100
     proc_imax = -1
     proc_jmax = -1
     proc_kmax = -1

     thisproc_have_global_bdry_min = have_global_bdry_min(int(levelnumber),:)
     thisproc_have_global_bdry_max = have_global_bdry_max(int(levelnumber),:)


     !write(*,*) "checkpoint (6). levelnumber is ", levelnumber

     if(thisproc_have_global_bdry_min(1)==1) proc_imin = 0
     if(thisproc_have_global_bdry_min(2)==1) proc_jmin = 0
     if(thisproc_have_global_bdry_min(3)==1) proc_kmin = 0
     if(thisproc_have_global_bdry_max(1)==1) proc_imax = glob_imax
     if(thisproc_have_global_bdry_max(2)==1) proc_jmax = glob_jmax
     if(thisproc_have_global_bdry_max(3)==1) proc_kmax = glob_kmax

     write(*,*) "HELLO in MHD POST_REGRID",levelnumber,thisproc_have_global_bdry_min,thisproc_have_global_bdry_max

     if (enable_primitives_disk == 0) then
        if(primitives_solver==11) then

           !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
                    ! Here we provide approximations to rho_b and P, as first guesses for the primitives solver
                    rho_b(i,j,k) = max(rho_star(i,j,k) * exp(-6.D0*phi(i,j,k)), rho_b_atm)
	            call compute_pcold_epscold(rho_b(i,j,k),P(i,j,k),temp1(i,j,k), & 
                         neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)

                    vx(i,j,k) = -shiftx(i,j,k)
                    vy(i,j,k) = -shifty(i,j,k)
                    vz(i,j,k) = -shiftz(i,j,k)

                    u0(i,j,k) = 1.D0/(lapm1(i,j,k)+1.d0)
                 end do
              end do
           end do
           !$omp end parallel do

           !Here we use the HARM primitives solver, with a new prescription that minimizes changes to
           !  conservatives without applying the "Font" fix.
           !We intend to make the below function truly generic, but currently it only supports Gamma=2.
           !  It can be trivially extended for arbitrary Gamma-law EOS's, but will require work for more
           !  generic EOS's.
           !We hope the below function will eventually be used in place of other primitives solver options,
           !  since it can be easily extended to use your desired technique for solving primitives.
           
	    !write(6,*) "(1)Inside driver_post_regrid, test gamma_tab =", gamma_tab
	    
	   if (levelnumber == max_refinement_level) then
	   write(*,*) "(1)levelnumber is at max refinement level!!!"
	   else if (levelnumber.le.max_refinement_level-1) then
	   call primitives_generic(cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z, &
                phi,gxx,gxy,gxz,gyy,gyz,gzz, &
                gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                lapm1,shiftx,shifty,shiftz, &
                Bx,By,Bz, &
                mhd_st_x,mhd_st_y,mhd_st_z,tau,rho_star, &
                vx,vy,vz,P,rho_b,h,u0, &
                rho_b_atm,tau_atm, rho_b_max,&
                neos,ergo_star, ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
                temp1,temp2,temp3,temp4,temp5, &
                primitives_debug,Psi6threshold,horizon_enforce_rho_profile);
	   write(*,*) "levelnumber is less than max_refinement level, use prim solver" 
	   end if		
           print *,"after primitives_generic ***************"

           !The goal of the below function is to just update the metric source terms.  We should
           !  really get rid of the old h,w,st_i,Ei gridfunctions...
           !Currently this only works for Gamma-law EOS's but could be extended for arbitary EOSs
           !  with some work.
	   if (levelnumber == max_refinement_level) then
           write(*,*) "(2)levelnumber is 7!!!"
           else if (levelnumber.le.max_refinement_level-1) then
           call metric_source_terms_and_misc_vars(cctkGH,cctk_lsh, &
                rho,Sx,Sy,Sz, &
                Sxx,Sxy,Sxz,Syy,Syz,Szz, &
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
		rho_b_atm);
		write(*,*) "call metric_source_terms_and_misc_vars, level number is",levelnumber
	   end if
           write(*,*) "HI tau=",tau(1,1,1),levelnumber
	  
        else if(primitives_solver==1) then
           !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
                    !Note that the correct eqn for w is: 
                    ! w=rho_star*(lapm1+1.0)*u0
                    ! We must fill in w everywhere here since, for primitive_vars_hybrid2(), we need the old value of w
                    !    everywhere on the grid.  This is a problem since w is not prolongated and hence is not set
                    !    everywhere after regridding.  The below line just provides an initial guess for primitive_vars_hybrid2()
                    !Here, temp2 is supposed to be w_p
                    temp2(i,j,k) = rho_star(i,j,k) !*(lapm1(i,j,k)+1.0)*1.D0 <-- Yuk Tung says it's better to just set w_p to rho_star in this case.
                 end do
              end do
           end do
           !$omp end parallel do

           !write(*,*) "PHI=",phi(79,80,3),ext
           !  Zach says: Warning: if ignore_ghostzones=1, we get undefined values post-regridding
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

           !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
                    ! Following lines are necessary
                    st_x(i,j,k) = mhd_st_x(i,j,k)
                    st_y(i,j,k) = mhd_st_y(i,j,k)
                    st_z(i,j,k) = mhd_st_z(i,j,k)
                    !MAYBE FIXME: Put in a warning if B!=0.  However, if B!=0, you shouldn't be using this primitives solver anyway!
                    sbt(i,j,k) = 0.D0
                    sbx(i,j,k) = 0.D0
                    sby(i,j,k) = 0.D0
                    sbz(i,j,k) = 0.D0
                    Bx(i,j,k) = 0.D0
                    By(i,j,k) = 0.D0
                    Bz(i,j,k) = 0.D0
                    Ex(i,j,k) = 0.D0
                    Ey(i,j,k) = 0.D0
                    Ez(i,j,k) = 0.D0
                 end do
              end do
           end do
           !$omp end parallel do

        else if(primitives_solver==2) then
           !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)                  
                    ! Here we provide approximations to rho_b and P, as first guesses for the primitives solver
                    rho_b(i,j,k) = max(rho_star(i,j,k) * exp(-6.D0*phi(i,j,k)), rho_b_atm)
	            call compute_pcold_epscold(rho_b(i,j,k),P(i,j,k),temp1(i,j,k), & 
                         neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
                    ! Finally, an approximation of h_p, so that eps = 1.D0:
                    !   (Recall that eps = h_old[index]-1.0-P_old[index]/rho_b_old)
	            h_p(i,j,k) = 1.d0 + P(i,j,k)/rho_b(i,j,k) + temp1(i,j,k)

                    vx(i,j,k) = -shiftx(i,j,k)
                    vy(i,j,k) = -shifty(i,j,k)
                    vz(i,j,k) = -shiftz(i,j,k)

                    u0(i,j,k) = 1.D0/(lapm1(i,j,k)+1.d0)
                 end do
              end do
           end do
           !$omp end parallel do

           !TODO: remove temp1 (unused) parameter.
           !  Zach says: Warning: if ignore_ghostzones=1, we get undefined values post-regridding
           call primitive_vars_hybrid2_cpp(ext,cctk_nghostzones,X,Y,Z, &
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
                sbt,sbx,sby,sbz, &
                proc_imin,proc_jmin,proc_kmin, &
                proc_imax,proc_jmax,proc_kmax, &
                glob_imax,glob_jmax,glob_kmax, &
                Symmetry,pfloor,excision_enable, &
                excision_zone_gf, tau_stildefix_enable,tau_atm,temp4,cctkGH,0, &
                enable_shocktest_primitive_mode,repairs_needed)

           if (repairs_needed==1) then
              call repair_failures_mhd_hybrid(ext,Z,gamma_th, &
                  temp4, rho_b, P, &
                  vx, vy, vz, u0, w, h, rho_star, tau, &
                  st_x, st_y, st_z, mhd_st_x, mhd_st_y, mhd_st_z, &
                  rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
                  lapm1, shiftx, shifty, shiftz, phi, &
                  gxx, gxy, gxz, gyy, gyz, gzz, &
                  gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Bx, By, Bz, &
                  sbt, sbx, sby, sbz, rho_b_atm, &
                  neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
                  proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
                  glob_imax,glob_jmax,glob_kmax,Symmetry)
           end if

           if(excision_enable==1) then
             call hydro_ezbc_hybrid(ext,X,Y,Z,rho_star,tau, &
                  mhd_st_x,mhd_st_y,mhd_st_z,st_x,st_y,st_z, &
                  rho_b,P,h,vx,vy,vz,w,&
                  sbt,sbx,sby,sbz,Bx,By,Bz, &
                  lapm1,shiftx,shifty,shiftz,phi,&
                  gxx,gxy,gxz,gyy,gyz,gzz, &
                  gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,&
                  Symmetry,excision_zone_gf,gamma_th,neos,ergo_star, ergo_sigma, rho_tab,&
                  P_tab,eps_tab,k_tab,gamma_tab)
             if (hyperbolic_divergence_cleaning_enable==1) then
                call scalar_excision_bc(ext,X,Y,Z,Blagrangemultiplier,Symmetry,excision_zone_gf)
             end if
             call remove_interior2(ext,X,Y,Z,sbt,excision_zone_gf,Symmetry)
             call remove_interior2(ext,X,Y,Z,sbx,excision_zone_gf,Symmetry)
             call remove_interior2(ext,X,Y,Z,sby,excision_zone_gf,Symmetry)
             call remove_interior2(ext,X,Y,Z,sbz,excision_zone_gf,Symmetry)
           end if

        else if(primitives_solver==0) then
           !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
                    !if(i.le.cctk_nghostzones(1) .or. j.le.cctk_nghostzones(2) .or. k.le.cctk_nghostzones(3)) then

                    ! Here we provide approximations to rho_b and P, as first guesses for the primitives solver
                    rho_b(i,j,k) = max(rho_star(i,j,k) * exp(-6.D0*phi(i,j,k)), rho_b_atm)
                    call compute_pcold_epscold(rho_b(i,j,k),P(i,j,k),temp1(i,j,k), &
                                                   neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
                    ! Finally, an approximation of h_p, so that eps = 1.D0:
                    !   (Recall that eps = h_old[index]-1.0-P_old[index]/rho_b_old)
                    h_p(i,j,k) = 1.d0 + P(i,j,k)/rho_b(i,j,k) + temp1(i,j,k)

                    vx(i,j,k) = 0.D0
                    vy(i,j,k) = 0.D0
                    vz(i,j,k) = 0.D0
                 end do
              end do
           end do
            !$omp end parallel do

           !TODO: remove temp1 (unused) parameter.
           !  Zach says: Warning: if ignore_ghostzones=1, we get undefined values post-regridding
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
                excision_zone_gf, tau_stildefix_enable,tau_atm,0)
           
        end if
     end if

     write(*,*) "HELLO postregrid:tau,dX=",tau(1,1,1),dX
  end if
end subroutine mhd_post_regrid
