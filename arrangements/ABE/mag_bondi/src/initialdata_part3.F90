#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine mag_bondi_initialdata_part3(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext
  real*8 			           :: dT,dX,dY,dZ,rhos_max,rho_max
  real*8 			           :: xmin,ymin,zmin,ymax
  real*8                                   :: hc_mask_radius, fac, max_b2oP, maxb2, maxP, tau_max
  real*8                                   :: detmin_l, detmax_l
  real*8                                   :: rho_fail_max_step,M_fail_step
  real*8                                   :: Xglobmin,Yglobmin,Zglobmin
  real*8                                   :: Xglobmax,Yglobmax,Zglobmax
  logical :: fish_to_phys
  integer :: n1, n2, n3, mf
  integer :: handle,dummy,proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax
  integer :: proc_kmax,glob_imax,glob_jmax,glob_kmax
  integer :: index, Nfont, Nfont_l, vindex
  integer :: ignore_ghostzones,repairs_needed
  integer :: ierr,ONE,ZERO, i,j,k
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  integer, parameter :: AXISYM_FULL = 5
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  CCTK_REAL reduction_value

  parameter(ONE = 1.D0, ZERO = 0.D0)
!
  ext = cctk_lsh

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  glob_imax = ext(1)
  glob_jmax = ext(2)
  glob_kmax = ext(3)

  !minimum of our grid:
  call CCTK_ReductionHandle(handle,"minimum")
  call CCTK_VarIndex(vindex,"grid::X")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,Xglobmin,1,vindex)  
  call CCTK_VarIndex(vindex,"grid::Y")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,Yglobmin,1,vindex)
  call CCTK_VarIndex(vindex,"grid::Z")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,Zglobmin,1,vindex)
  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(vindex,"grid::X")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,Xglobmax,1,vindex)
  call CCTK_VarIndex(vindex,"grid::Y")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,Yglobmax,1,vindex)
  call CCTK_VarIndex(vindex,"grid::Z")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,Zglobmax,1,vindex)
  if(X(1,1,1) .eq. Xglobmin) then
     proc_imin = 0
  else
     proc_imin = -100
  end if
  if(Y(1,1,1) .eq. Yglobmin) then
     proc_jmin = 0
  else
     proc_jmin = -100
  end if
  if(Z(1,1,1) .eq. Zglobmin) then
     proc_kmin = 0
  else
     proc_kmin = -100
  end if

  if(X(cctk_lsh(1),1,1) .eq. Xglobmax) then
     proc_imax = glob_imax
  else
     proc_imax = -1
  end if
  if(Y(1,cctk_lsh(2),1) .eq. Yglobmax) then
     proc_jmax = glob_jmax
  else
     proc_jmax = -1
  end if
  if(Z(1,1,cctk_lsh(3)) .eq. Zglobmax) then
     proc_kmax = glob_kmax
  else
     proc_kmax = -1
  end if
  
  pfloor = P_fact
  rho_b_atm_gf = rho_b_atm 
  pfloor_gf = pfloor 

  !!
  !! Setup matter and EM profile   !!
  ! NOTE THAT Xglobmax,Zglobmax is NOT set properly above! (This doesn't matter, since the primitives solver doesn't really need it.)

  call mag_bondi_accreting_matter(ext, X, Y, Z, &
       neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, gamma_th, &
     rho_star, tau, st_x, st_y, st_z, mhd_st_x,mhd_st_y,mhd_st_z, &
     w, h, rho_b, P, u0, vx,vy,vz, npoly, r_sonic, constrained_transport_scheme, & 
     Ax,Ay,Az,Bx_stagger,By_stagger,Bz_stagger, &
     Bx,By,Bz,sbt,sbx,sby,sbz,Ex,Ey,Ez,Bxtilde,Bytilde,Bztilde, &
     Bstrength,r0, puncture_id, excision_enable, excision_radius, xmax_bondi,zmax_bondi)

  write(*,*) 'XMAX? ZMAX?', xmax_bondi,zmax_bondi

  tau_atm = tau_fact

  ! Let's initialize the Fontfix_tracker_gf array
  Fontfix_tracker_gf = 0.d0
  rhos_max = 1.d0
  rho_max = 1.d0

  h_p = h

  if (enable_primitives_disk == 0) then
     if(primitives_solver==1) then
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
             rhos_max, rho_fail_max_step, M_fail_step, rho_b_atm, &
             gamma_th,K_poly,sdots_o_rhot,Symmetry,0)

        ! Following lines are necessary
        !$omp parallel do
        do k=1,cctk_lsh(3)
           do j=1,cctk_lsh(2)
              do i=1,cctk_lsh(1)
                 st_x(i,j,k) = mhd_st_x(i,j,k)
                 st_y(i,j,k) = mhd_st_y(i,j,k)
                 st_z(i,j,k) = mhd_st_z(i,j,k)

                 ! Duh, the B fields are zero when using the quartic primitives solver, so we must
                 !  set sb's to zero:
                 sbt(i,j,k) = 0.D0
                 sbx(i,j,k) = 0.D0
                 sby(i,j,k) = 0.D0
                 sbz(i,j,k) = 0.D0
              end do
           end do
        end do
        !$omp end parallel do
     else if(primitives_solver==0) then

        !For initial data, h_p HAS NOT yet been set:
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
             mhd_st_x,mhd_st_y,mhd_st_z,neos, &
             rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, &
             w,temp1,rho_b,rho,P,h,Sx,Sy,Sz, &
             Sxx,Sxy,Sxz,Syy,Syz,Szz, &
             phi,lapm1,shiftx,shifty,shiftz, &
             gxx,gxy,gxz,gyy,gyz,gzz, &
             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
             h_p,u0,rhos_max,rho_b_atm, &
             rho_fail_max_step,M_fail_step,1.d0, &
             Bx,By,Bz,Ex,Ey,Ez, &
             vx,vy,vz, &
             sbt,sbx,sby,sbz, temp4, &
             proc_imin,proc_jmin,proc_kmin, &
             proc_imax,proc_jmax,proc_kmax, &
             glob_imax,glob_jmax,glob_kmax, &
             Symmetry,pfloor,excision_enable, &
             excision_zone_gf, tau_stildefix_enable,tau_atm,0)

     else if(primitives_solver==2) then

        !For initial data, h_p HAS NOT yet been set:
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
        ignore_ghostzones = 0
        call primitive_vars_hybrid2_cpp(ext,cctk_nghostzones,X,Y,Z, &
             rho_star,tau,st_x,st_y,st_z, &
             mhd_st_x,mhd_st_y,mhd_st_z,neos, &
             rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, &
             w,temp1,rho_b,rho,P,h,Sx,Sy,Sz, &
             Sxx,Sxy,Sxz,Syy,Syz,Szz, &
             phi,lapm1,shiftx,shifty,shiftz, &
             gxx,gxy,gxz,gyy,gyz,gzz, &
             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
             h_p,u0,rho_max,rho_b_atm, &
             rho_fail_max_step,M_fail_step,rhos_max, &
             Bx,By,Bz,Ex,Ey,Ez, &
             vx,vy,vz, &
             sbt,sbx,sby,sbz, &
             proc_imin,proc_jmin,proc_kmin, &
             proc_imax,proc_jmax,proc_kmax, &
             glob_imax,glob_jmax,glob_kmax, &
             Symmetry,pfloor,excision_enable, &
             excision_zone_gf, tau_stildefix_enable,tau_atm,temp1,cctkgh,ignore_ghostzones, &
             enable_shocktest_primitive_mode,repairs_needed)

     end if
  else if(enable_primitives_disk==1 .and. use_harm_primitives==0) then
     call primitive_vars_alt_disk(ext,X,Y,Z, &
          rho_star,tau,st_x,st_y,st_z, &
          mhd_st_x,mhd_st_y,mhd_st_z,neos, &
          rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, &
          w,rho_b,rho,P,h,Sx,Sy,Sz, &
          Sxx,Sxy,Sxz,Syy,Syz,Szz, &
          phi,lapm1,shiftx,shifty,shiftz, &
          gxx,gxy,gxz,gyy,gyz,gzz, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
          h_p,u0,rhos_max,rho_b_atm_gf, &
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
!!$  else if(use_harm_primitives==1) then
!!$     ! call the C function--a wrapper for the Noble et al code.
!!$
!!$     ! At one time, the second line of the function call was:
!!$     ! rho_b, rho_b_p, P, P_p, vx, vy, vz, &
!!$     call primitive_vars_harm(NFont_L, CCTKGH, X, Y, Z, ext, &
!!$          rho_star, tau, mhd_st_x, mhd_st_y, mhd_st_z, &
!!$          rho_b, rho_b, P, P, vx, vy, vz, &
!!$          u0, Bx, By, Bz, lapm1, phi, &
!!$          shiftx, shifty, shiftz, Fontfix_tracker_gf, rho_b_atm_gf, &
!!$          gxx, gxy, gxz, gyy, gyz, gzz, &
!!$          gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
!!$          Symmetry,  excision_radius, &
!!$          glob_imax, glob_jmax, glob_kmax, &
!!$          proc_imax, proc_jmax, proc_kmax, &
!!$          proc_imin, proc_jmin, proc_kmin)
!!$
!!$     call harm_impose_floors(ext,rho_b_atm_gf,pfloor_gf,rho_b,P, &
!!$          vx,vy,vz,Bx,By,Bz, &
!!$          rho_star,tau,mhd_st_x,mhd_st_y,mhd_st_z, &
!!$          lapm1,phi,shiftx,shifty,shiftz, &
!!$          gxx,gxy,gxz,gyy,gyz,gzz, &
!!$          neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th)
!!$
!!$     ! make sure you have the failures array stored in Fontfix_tracker_gf.
!!$     call repair_failures_mhd_disk(ext,X,Z,gamma_th,Fontfix_tracker_gf,rho_b,P, &
!!$          vx,vy,vz,u0,w,h,rho_star,tau, &
!!$          st_x,st_y,st_z,mhd_st_x,mhd_st_y,mhd_st_z, &
!!$          rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz, &
!!$          lapm1,shiftx,shifty,shiftz,phi, &
!!$          gxx,gxy,gxz,gyy,gyz,gzz, &
!!$          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,Bx,By,Bz, &
!!$          sbt,sbx,sby,sbz,rho_b_atm_gf, &
!!$          neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
!!$          proc_imin,proc_jmin,proc_kmin, &
!!$          proc_imax,proc_jmax,proc_kmax, &
!!$          glob_imax,glob_jmax,glob_kmax, &
!!$          Symmetry,pfloor_gf,enable_HARM_energyvariable,excision_radius)
!!$
!!$     call calculate_aux_harm(ext,gamma_th,rho_b,P, &
!!$          vx,vy,vz,u0,w,h,rho_star,tau, &
!!$          st_x,st_y,st_z,mhd_st_x,mhd_st_y,mhd_st_z, &
!!$          rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz, &
!!$          lapm1,shiftx,shifty,shiftz,phi, &
!!$          gxx,gxy,gxz,gyy,gyz,gzz, &
!!$          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,Bx,By,Bz, &
!!$          sbt,sbx,sby,sbz, &
!!$          neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,proc_kmin)
!!$  end if

  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_conservatives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_primitives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vs')

  position_x(1) = 0.d0
  position_y(1) = 0.d0
  position_z(1) = 0.d0

  ! You should probably convert these diagnostics
!!$  adm::AH_Manager.Find_Horizon(level,0.0,PhysTime);
!!$  double J_dyn = 0.0;
!!$//   adm::AH_Manager.dynamical_horizon_J(J_dyn,phi,
!!$// 				      gupxx,gupxy,gupxz,
!!$// 				      gupyy,gupyz,gupzz,
!!$// 				      Axx,Axy,Axz,
!!$// 				      Ayy,Ayz,Azz,trK,
!!$// 				      level,IDENT);
!!$  adm::AH_Manager.dynamical_horizon_J_ana(J_dyn,
!!$					  n_tot_horizon, N_theta, N_phi,
!!$					  xn, yn, zn, nx_g, ny_g, nz_g, nn,
!!$					  phi,
!!$					  gupxx,gupxy,gupxz,
!!$					  gupyy,gupyz,gupzz,
!!$					  Axx,Axy,Axz,
!!$					  Ayy,Ayz,Azz,trK,
!!$					  level,IDENT);
!!$  if (I_Print) cout << "J_dynamical_horizon = " << J_dyn << endl;

end subroutine mag_bondi_initialdata_part3
