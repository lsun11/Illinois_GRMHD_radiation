!-----------------------------------------------------------------------------
! Update outer boundaries
!-----------------------------------------------------------------------------
subroutine apply_matter_bc_mhd_new(ext, fake_ext, X,Y,Z, &
     rho_star, tau, mhd_st_x, mhd_st_y, mhd_st_z, &
     st_x, st_y, st_z, &
     tau_rad, S_rad_x, S_rad_y, S_rad_z, rho_Ye,&
     tau_rad_nue, S_rad_x_nue, S_rad_y_nue, S_rad_z_nue, tau_rad_nux, S_rad_x_nux, S_rad_y_nux, S_rad_z_nux, &
     rho_b, P, vx, vy, vz, Bx, By, Bz, &
     E_rad, F_radx, F_rady, F_radz, Y_e, &
     E_rad_nue, F_radx_nue, F_rady_nue, F_radz_nue, E_rad_nux, F_radx_nux, F_rady_nux, F_radz_nux,&
     F_rad0, P_radxx, P_radyy, P_radzz, P_radxy, P_radxz, P_radyz, &
     F_rad0_nue, P_radxx_nue, P_radyy_nue, P_radzz_nue, P_radxy_nue, P_radxz_nue, P_radyz_nue, &
     F_rad0_nux, P_radxx_nux, P_radyy_nux, P_radzz_nux, P_radxy_nux, P_radxz_nux, P_radyz_nux, &
     eps_tot, eps_cld, eps_thermal, P_cld, &
     w, h, u0, Ex, Ey, Ez, sbt, sbx, sby, sbz, T_fluid,&
     rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
     phi, lapse, shiftx, shifty, shiftz, &
     gxx, gxy, gxz, gyy, gyz, gzz, &
     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
     Symmetry,bc, &
     neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, rho_b_atm, rad_const,  &
     have_bdry_min,have_bdry_max, compute_microphysics, enable_OS_collapse, rad_evolve_enable, rad_closure_scheme, microphysics_scheme, Erad_atm_cut)
  implicit none
  integer, dimension(3)                   :: ext,fake_ext,have_bdry_min,have_bdry_max
  real*8, dimension(ext(1),ext(2),ext(3)) :: X,Y,Z
  real*8, dimension(ext(1),ext(2),ext(3)) :: rho_star, tau, mhd_st_x, mhd_st_y, mhd_st_z
  real*8, dimension(ext(1),ext(2),ext(3)) :: st_x, st_y, st_z
  real*8, dimension(ext(1),ext(2),ext(3)) :: tau_rad, S_rad_x, S_rad_y, S_rad_z, rho_Ye
  real*8, dimension(ext(1),ext(2),ext(3)) :: tau_rad_nue, S_rad_x_nue, S_rad_y_nue, S_rad_z_nue, tau_rad_nux, S_rad_x_nux, S_rad_y_nux, S_rad_z_nux
  real*8, dimension(ext(1),ext(2),ext(3)) :: rho_b, P, vx, vy, vz, Bx, By, Bz
  real*8, dimension(ext(1),ext(2),ext(3)) :: E_rad, F_radx, F_rady, F_radz, Y_e
  real*8, dimension(ext(1),ext(2),ext(3)) :: E_rad_nue, F_radx_nue, F_rady_nue, F_radz_nue, E_rad_nux, F_radx_nux, F_rady_nux, F_radz_nux
  real*8, dimension(ext(1),ext(2),ext(3)) :: F_rad0, P_radxx, P_radyy, P_radzz, P_radxy, P_radxz, P_radyz
  real*8, dimension(ext(1),ext(2),ext(3)) :: F_rad0_nue, P_radxx_nue, P_radyy_nue, P_radzz_nue, P_radxy_nue, P_radxz_nue, P_radyz_nue
  real*8, dimension(ext(1),ext(2),ext(3)) :: F_rad0_nux, P_radxx_nux, P_radyy_nux, P_radzz_nux, P_radxy_nux, P_radxz_nux, P_radyz_nux
  real*8, dimension(ext(1),ext(2),ext(3)) :: eps_tot, eps_cld, eps_thermal, P_cld
  real*8, dimension(ext(1),ext(2),ext(3)) :: w, h, u0, Ex, Ey, Ez, sbt, sbx, sby, sbz, T_fluid
  real*8, dimension(ext(1),ext(2),ext(3)) :: rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz
  real*8, dimension(ext(1),ext(2),ext(3)) :: phi, lapse, shiftx, shifty, shiftz
  real*8, dimension(ext(1),ext(2),ext(3)) :: gxx, gxy, gxz, gyy, gyz, gzz
  real*8, dimension(ext(1),ext(2),ext(3)) :: gupxx, gupxy, gupxz, gupyy, gupyz, gupzz
  integer                                 :: Symmetry
  integer                                 :: bc
  integer                                 :: neos,ergo_star,enable_OS_collapse,rad_evolve_enable,rad_closure_scheme,microphysics_scheme,compute_microphysics
  real*8, dimension(10)                 :: rho_tab,P_tab,eps_tab
  real*8, dimension(11)               :: k_tab, gamma_tab
  real*8                                  :: rho_b_atm, rad_const,  gamma_th,ergo_sigma, Erad_atm_cut
!
! Other variables:
! 
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin, imax, jmax, kmax
  real*8, parameter                  :: ZERO = 0.D0
  real*8, parameter                  :: TWO = 2.D0, THREE = 3.D0
  real*8                             :: nx,ny,nz
  integer                            :: NO_SYMM, EQUATORIAL, OCTANT
  integer                            :: PI_SYMM, AXISYM
  integer                            :: FREEZE, EXTRAP, PERIODIC, OUTF, COPY, QUAD
  integer                            :: PLANAR
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(OUTF = 1, FREEZE = 2, COPY = 3, EXTRAP = 4, QUAD = 5, PLANAR = 6)
! Storage variables:
  integer, parameter                 :: m = 17
  real*8, dimension(m)               :: AUX
  real*8                             :: rho_bl, Pl, vxl, vyl, vzl, Bxl, Byl, Bzl
  real*8                             :: rho_sl, taul, mhd_st_xl, mhd_st_yl, mhd_st_zl
  real*8                             :: st_xl, st_yl, st_zl, wl, hl, u0l, Exl, Eyl, Ezl
  real*8                             :: E_radl,F_radxl,F_radyl,F_radzl,Y_el, T_fluidl
  real*8                             :: E_rad_nuel,F_radx_nuel,F_rady_nuel,F_radz_nuel,E_rad_nuxl,F_radx_nuxl,F_rady_nuxl,F_radz_nuxl
  real*8                             :: F_rad0l, P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl
  real*8                             :: F_rad0_nuel, P_radxx_nuel, P_radyy_nuel, P_radzz_nuel, P_radxy_nuel, P_radxz_nuel, P_radyz_nuel
  real*8                             :: F_rad0_nuxl, P_radxx_nuxl, P_radyy_nuxl, P_radzz_nuxl, P_radxy_nuxl, P_radxz_nuxl, P_radyz_nuxl
  real*8                             :: eps_totl, eps_cldl, eps_thermall, P_cldl
  real*8                             :: tau_radl, S_rad_xl, S_rad_yl, S_rad_zl, rho_Yel
  real*8                             :: tau_rad_nuel, S_rad_x_nuel, S_rad_y_nuel, S_rad_z_nuel, tau_rad_nuxl, S_rad_x_nuxl, S_rad_y_nuxl, S_rad_z_nuxl
  real*8                             :: sbtl, sbxl, sbyl, sbzl, rhol, Sxl, Syl, Szl
  real*8                             :: Sxxl, Sxyl, Sxzl, Syyl, Syzl, Szzl
!
! Input translation
!
  imin = ext(1)-fake_ext(1)+1
  jmin = ext(2)-fake_ext(2)+1
  kmin = ext(3)-fake_ext(3)+1
  imax = fake_ext(1)
  jmax = fake_ext(2)
  kmax = fake_ext(3)
  !first apply boundary conditions on the hydro primitive variables
  if(bc==FREEZE) then
     write(*,*) 'Warning:  Freeze boundary conditions are not supported!'
     stop
  else if(bc==OUTF) then
     if (have_bdry_max(1)==1) then
        rho_b(imax,:,:) = rho_b(imax-1,:,:)
        P(imax,:,:) = P(imax-1,:,:)
        vx(imax,:,:) = vx(imax-1,:,:)
        vy(imax,:,:) = vy(imax-1,:,:)
        vz(imax,:,:) = vz(imax-1,:,:)
        where (vx(imax,:,:) .lt. ZERO)
              vx(imax,:,:) = ZERO
        end where
        if (rad_evolve_enable==1) then
           E_rad(imax,:,:) = E_rad(imax-1,:,:)
           F_radx(imax,:,:) = F_radx(imax-1,:,:)
           F_rady(imax,:,:) = F_rady(imax-1,:,:)
           F_radz(imax,:,:) = F_radz(imax-1,:,:)
           where (F_radx(imax,:,:) .lt. ZERO)
              F_radx(imax,:,:) = ZERO
           end where
           if (microphysics_scheme==1) then
              E_rad_nue(imax,:,:) = E_rad_nue(imax-1,:,:)
              F_radx_nue(imax,:,:) = F_radx_nue(imax-1,:,:)
              F_rady_nue(imax,:,:) = F_rady_nue(imax-1,:,:)
              F_radz_nue(imax,:,:) = F_radz_nue(imax-1,:,:)
              where (F_radx_nue(imax,:,:) .lt. ZERO)
                 F_radx_nue(imax,:,:) = ZERO
              end where
              E_rad_nux(imax,:,:) = E_rad_nux(imax-1,:,:)
              F_radx_nux(imax,:,:) = F_radx_nux(imax-1,:,:)
              F_rady_nux(imax,:,:) = F_rady_nux(imax-1,:,:)
              F_radz_nux(imax,:,:) = F_radz_nux(imax-1,:,:)
              where (F_radx_nux(imax,:,:) .lt. ZERO)
                 F_radx_nux(imax,:,:) = ZERO
              end where
           end if
        end if
        Y_e(imax,:,:) = Y_e(imax-1,:,:)
        T_fluid(imax,:,:) = T_fluid(imax-1,:,:)
     end if
     if (have_bdry_min(1)==1) then
        rho_b(imin,:,:) = rho_b(imin+1,:,:)
        P(imin,:,:) = P(imin+1,:,:)
        vx(imin,:,:) = vx(imin+1,:,:)
        vy(imin,:,:) = vy(imin+1,:,:)
        vz(imin,:,:) = vz(imin+1,:,:)
        where (vx(imin,:,:) .gt. ZERO)
           vx(imin,:,:) = ZERO
        end where
        if (rad_evolve_enable==1) then
           E_rad(imin,:,:) = E_rad(imin+1,:,:)
           F_radx(imin,:,:) = F_radx(imin+1,:,:)
           F_rady(imin,:,:) = F_rady(imin+1,:,:)
           F_radz(imin,:,:) = F_radz(imin+1,:,:)
           where (F_radx(imin,:,:) .gt. ZERO)
              F_radx(imin,:,:) = ZERO
           end where
           if (microphysics_scheme==1) then
              E_rad_nue(imin,:,:) = E_rad_nue(imin+1,:,:)
              F_radx_nue(imin,:,:) = F_radx_nue(imin+1,:,:)
              F_rady_nue(imin,:,:) = F_rady_nue(imin+1,:,:)
              F_radz_nue(imin,:,:) = F_radz_nue(imin+1,:,:)
              where (F_radx_nue(imin,:,:) .gt. ZERO)
                 F_radx_nue(imin,:,:) = ZERO
              end where
              E_rad_nux(imin,:,:) = E_rad_nux(imin+1,:,:)
              F_radx_nux(imin,:,:) = F_radx_nux(imin+1,:,:)
              F_rady_nux(imin,:,:) = F_rady_nux(imin+1,:,:)
              F_radz_nux(imin,:,:) = F_radz_nux(imin+1,:,:)
              where (F_radx_nux(imin,:,:) .gt. ZERO)
                 F_radx_nux(imin,:,:) = ZERO
              end where
           end if
        end if
        Y_e(imin,:,:) = Y_e(imin+1,:,:)
        T_fluid(imin,:,:) = T_fluid(imin+1,:,:)
     end if
     if (have_bdry_max(2)==1) then
        rho_b(:,jmax,:) = rho_b(:,jmax-1,:)
        P(:,jmax,:) = P(:,jmax-1,:)
        vx(:,jmax,:) = vx(:,jmax-1,:)
        vy(:,jmax,:) = vy(:,jmax-1,:)
        vz(:,jmax,:) = vz(:,jmax-1,:)
        where (vy(:,jmax,:) .lt. ZERO)
               vy(:,jmax,:) = ZERO
        end where
        if (rad_evolve_enable==1) then
           E_rad(:,jmax,:) = E_rad(:,jmax-1,:)
           F_radx(:,jmax,:) = F_radx(:,jmax-1,:)
           F_rady(:,jmax,:) = F_rady(:,jmax-1,:)
           F_radz(:,jmax,:) = F_radz(:,jmax-1,:)
           where (F_rady(:,jmax,:) .lt. ZERO)
              F_rady(:,jmax,:) = ZERO
           end where
           if (microphysics_scheme==1) then
              E_rad_nue(:,jmax,:) = E_rad_nue(:,jmax-1,:)
              F_radx_nue(:,jmax,:) = F_radx_nue(:,jmax-1,:)
              F_rady_nue(:,jmax,:) = F_rady_nue(:,jmax-1,:)
              F_radz_nue(:,jmax,:) = F_radz_nue(:,jmax-1,:)
              where (F_rady_nue(:,jmax,:) .lt. ZERO)
                 F_rady_nue(:,jmax,:) = ZERO
              end where
              E_rad_nux(:,jmax,:) = E_rad_nux(:,jmax-1,:)
              F_radx_nux(:,jmax,:) = F_radx_nux(:,jmax-1,:)
              F_rady_nux(:,jmax,:) = F_rady_nux(:,jmax-1,:)
              F_radz_nux(:,jmax,:) = F_radz_nux(:,jmax-1,:)
              where (F_rady_nux(:,jmax,:) .lt. ZERO)
                 F_rady_nux(:,jmax,:) = ZERO
              end where
           end if
        end if
        Y_e(:,jmax,:) = Y_e(:,jmax-1,:)
        T_fluid(:,jmax,:) = T_fluid(:,jmax-1,:)
     end if
     if (have_bdry_min(2)==1) then
        rho_b(:,jmin,:) = rho_b(:,jmin+1,:)
        P(:,jmin,:) = P(:,jmin+1,:)
        vx(:,jmin,:) = vx(:,jmin+1,:)
        vy(:,jmin,:) = vy(:,jmin+1,:)
        vz(:,jmin,:) = vz(:,jmin+1,:)
        where (vy(:,jmin,:) .gt. ZERO)
           vy(:,jmin,:) = ZERO
        end where
        if (rad_evolve_enable==1) then
           E_rad(:,jmin,:) = E_rad(:,jmin+1,:)
           F_radx(:,jmin,:) = F_radx(:,jmin+1,:)
           F_rady(:,jmin,:) = F_rady(:,jmin+1,:)
           F_radz(:,jmin,:) = F_radz(:,jmin+1,:)
           where (F_rady(:,jmin,:) .gt. ZERO)
              F_rady(:,jmin,:) = ZERO
           end where
           if (microphysics_scheme==1) then
              E_rad_nue(:,jmin,:) = E_rad_nue(:,jmin+1,:)
              F_radx_nue(:,jmin,:) = F_radx_nue(:,jmin+1,:)
              F_rady_nue(:,jmin,:) = F_rady_nue(:,jmin+1,:)
              F_radz_nue(:,jmin,:) = F_radz_nue(:,jmin+1,:)
              where (F_rady_nue(:,jmin,:) .gt. ZERO)
                 F_rady_nue(:,jmin,:) = ZERO
              end where
              E_rad_nux(:,jmin,:) = E_rad_nux(:,jmin+1,:)
              F_radx_nux(:,jmin,:) = F_radx_nux(:,jmin+1,:)
              F_rady_nux(:,jmin,:) = F_rady_nux(:,jmin+1,:)
              F_radz_nux(:,jmin,:) = F_radz_nux(:,jmin+1,:)
              where (F_rady_nux(:,jmin,:) .gt. ZERO)
                 F_rady_nux(:,jmin,:) = ZERO
              end where
           end if
        end if
        Y_e(:,jmin,:) = Y_e(:,jmin+1,:)
        T_fluid(:,jmin,:) = T_fluid(:,jmin+1,:)
     end if
     if (have_bdry_max(3)==1) then
        rho_b(:,:,kmax) = rho_b(:,:,kmax-1)
        P(:,:,kmax) = P(:,:,kmax-1)
        vx(:,:,kmax) = vx(:,:,kmax-1)
        vy(:,:,kmax) = vy(:,:,kmax-1)
        vz(:,:,kmax) = vz(:,:,kmax-1)
        where (vz(:,:,kmax) .lt. ZERO)
               vz(:,:,kmax) = ZERO
        end where
        if (rad_evolve_enable==1) then
           E_rad(:,:,kmax) = E_rad(:,:,kmax-1)
           F_radx(:,:,kmax) = F_radx(:,:,kmax-1)
           F_rady(:,:,kmax) = F_rady(:,:,kmax-1)
           F_radz(:,:,kmax) = F_radz(:,:,kmax-1)
           where (F_radz(:,:,kmax) .lt. ZERO)
              F_radz(:,:,kmax) = ZERO
           end where
           if (microphysics_scheme==1) then
              E_rad_nue(:,:,kmax) = E_rad_nue(:,:,kmax-1)
              F_radx_nue(:,:,kmax) = F_radx_nue(:,:,kmax-1)
              F_rady_nue(:,:,kmax) = F_rady_nue(:,:,kmax-1)
              F_radz_nue(:,:,kmax) = F_radz_nue(:,:,kmax-1)
              where (F_radz_nue(:,:,kmax) .lt. ZERO)
                 F_radz_nue(:,:,kmax) = ZERO
              end where
              E_rad_nux(:,:,kmax) = E_rad_nux(:,:,kmax-1)
              F_radx_nux(:,:,kmax) = F_radx_nux(:,:,kmax-1)
              F_rady_nux(:,:,kmax) = F_rady_nux(:,:,kmax-1)
              F_radz_nux(:,:,kmax) = F_radz_nux(:,:,kmax-1)
              where (F_radz_nux(:,:,kmax) .lt. ZERO)
                 F_radz_nux(:,:,kmax) = ZERO
              end where
           end if
        end if
        Y_e(:,:,kmax) = Y_e(:,:,kmax-1)
        T_fluid(:,:,kmax) = T_fluid(:,:,kmax-1)

     end if
     if (have_bdry_min(3)==1) then
        rho_b(:,:,kmin) = rho_b(:,:,kmin+1)
        P(:,:,kmin) = P(:,:,kmin+1)
        vx(:,:,kmin) = vx(:,:,kmin+1)
        vy(:,:,kmin) = vy(:,:,kmin+1)
        vz(:,:,kmin) = vz(:,:,kmin+1)
        where (vz(:,:,kmin) .gt. ZERO)
              vz(:,:,kmin) = ZERO
        end where
        if (rad_evolve_enable==1) then
           E_rad(:,:,kmin) = E_rad(:,:,kmin+1)
           F_radx(:,:,kmin) = F_radx(:,:,kmin+1)
           F_rady(:,:,kmin) = F_rady(:,:,kmin+1)
           F_radz(:,:,kmin) = F_radz(:,:,kmin+1)
           where (F_radz(:,:,kmin) .gt. ZERO)
              F_radz(:,:,kmin) = ZERO
           end where
           if (microphysics_scheme==1) then
              E_rad_nue(:,:,kmin) = E_rad_nue(:,:,kmin+1)
              F_radx_nue(:,:,kmin) = F_radx_nue(:,:,kmin+1)
              F_rady_nue(:,:,kmin) = F_rady_nue(:,:,kmin+1)
              F_radz_nue(:,:,kmin) = F_radz_nue(:,:,kmin+1)
              where (F_radz_nue(:,:,kmin) .gt. ZERO)
                 F_radz_nue(:,:,kmin) = ZERO
              end where
              E_rad_nux(:,:,kmin) = E_rad_nux(:,:,kmin+1)
              F_radx_nux(:,:,kmin) = F_radx_nux(:,:,kmin+1)
              F_rady_nux(:,:,kmin) = F_rady_nux(:,:,kmin+1)
              F_radz_nux(:,:,kmin) = F_radz_nux(:,:,kmin+1)
              where (F_radz_nux(:,:,kmin) .gt. ZERO)
                 F_radz_nux(:,:,kmin) = ZERO
              end where
           end if
        end if
        Y_e(:,:,kmin) = Y_e(:,:,kmin+1)
        T_fluid(:,:,kmin) = T_fluid(:,:,kmin+1)
     end if
  else if(bc==EXTRAP) then
     if (have_bdry_max(1)==1) then
        rho_b(imax,:,:) = TWO*rho_b(imax-1,:,:) &
             - rho_b(imax-2,:,:)
        P(imax,:,:) = TWO*P(imax-1,:,:) &
             - P(imax-2,:,:)
        vx(imax,:,:) = TWO*vx(imax-1,:,:) &
             - vx(imax-2,:,:)
        vy(imax,:,:) = TWO*vy(imax-1,:,:) &
             - vy(imax-2,:,:)
        vz(imax,:,:) = TWO*vz(imax-1,:,:) &
             - vz(imax-2,:,:)
     end if
     if (have_bdry_max(2)==1) then
        rho_b(:,jmax,:) = TWO*rho_b(:,jmax-1,:) &
             - rho_b(:,jmax-2,:)
        P(:,jmax,:) = TWO*P(:,jmax-1,:) &
             - P(:,jmax-2,:)
        vx(:,jmax,:) = TWO*vx(:,jmax-1,:) &
             - vx(:,jmax-2,:)
        vy(:,jmax,:) = TWO*vy(:,jmax-1,:) &
             - vy(:,jmax-2,:)
        vz(:,jmax,:) = TWO*vz(:,jmax-1,:) &
             - vz(:,jmax-2,:)
     end if
     if (have_bdry_max(3)==1) then
        rho_b(:,:,kmax) = TWO*rho_b(:,:,kmax-1) &
             - rho_b(:,:,kmax-2)
        P(:,:,kmax) = TWO*P(:,:,kmax-1) &
             - P(:,:,kmax-2)
        vx(:,:,kmax) = TWO*vx(:,:,kmax-1) &
             - vx(:,:,kmax-2)
        vy(:,:,kmax) = TWO*vy(:,:,kmax-1) &
             - vy(:,:,kmax-2)
        vz(:,:,kmax) = TWO*vz(:,:,kmax-1) &
             - vz(:,:,kmax-2)
     end if
     if (have_bdry_min(1)==1) then
        rho_b(imin,:,:) = TWO*rho_b(imin+1,:,:) &
             - rho_b(imin+2,:,:)
        P(imin,:,:) = TWO*P(imin+1,:,:) &
             - P(imin+2,:,:)
        vx(imin,:,:) = TWO*vx(imin+1,:,:) &
             - vx(imin+2,:,:)
        vy(imin,:,:) = TWO*vy(imin+1,:,:) &
             - vy(imin+2,:,:)
        vz(imin,:,:) = TWO*vz(imin+1,:,:) &
             - vz(imin+2,:,:)
     end if
     if (have_bdry_min(2)==1) then
        rho_b(:,jmin,:) = TWO*rho_b(:,jmin+1,:) &
             - rho_b(:,jmin+2,:)
        P(:,jmin,:) = TWO*P(:,jmin+1,:) &
             - P(:,jmin+2,:)
        vx(:,jmin,:) = TWO*vx(:,jmin+1,:) &
             - vx(:,jmin+2,:)
        vy(:,jmin,:) = TWO*vy(:,jmin+1,:) &
             - vy(:,jmin+2,:)
        vz(:,jmin,:) = TWO*vz(:,jmin+1,:) &
             - vz(:,jmin+2,:)
     end if
     if (have_bdry_min(3)==1) then
        rho_b(:,:,kmin) = TWO*rho_b(:,:,kmin+1) &
             - rho_b(:,:,kmin+2)
        P(:,:,kmin) = TWO*P(:,:,kmin+1) &
             - P(:,:,kmin+2)
        vx(:,:,kmin) = TWO*vx(:,:,kmin+1) &
             - vx(:,:,kmin+2)
        vy(:,:,kmin) = TWO*vy(:,:,kmin+1) &
             - vy(:,:,kmin+2)
        vz(:,:,kmin) = TWO*vz(:,:,kmin+1) &
             - vz(:,:,kmin+2)
     end if

     
     if (have_bdry_max(1)==1) then
        if (rad_evolve_enable==1) then
           E_rad(imax,:,:) = TWO*E_rad(imax-1,:,:) &
                - E_rad(imax-2,:,:)
           F_radx(imax,:,:) = TWO*F_radx(imax-1,:,:) &
                - F_radx(imax-2,:,:)
           F_rady(imax,:,:) = TWO*F_rady(imax-1,:,:) &
                - F_rady(imax-2,:,:)
           F_radz(imax,:,:) = TWO*F_radz(imax-1,:,:) &
                - F_radz(imax-2,:,:)
           if (microphysics_scheme==1) then
              E_rad_nue(imax,:,:) = TWO*E_rad_nue(imax-1,:,:) &
                   - E_rad_nue(imax-2,:,:)
              F_radx_nue(imax,:,:) = TWO*F_radx_nue(imax-1,:,:) &
                   - F_radx_nue(imax-2,:,:)
              F_rady_nue(imax,:,:) = TWO*F_rady_nue(imax-1,:,:) &
                   - F_rady_nue(imax-2,:,:)
              F_radz_nue(imax,:,:) = TWO*F_radz_nue(imax-1,:,:) &
                   - F_radz_nue(imax-2,:,:)

              E_rad_nux(imax,:,:) = TWO*E_rad_nux(imax-1,:,:) &
                   - E_rad_nux(imax-2,:,:)
              F_radx_nux(imax,:,:) = TWO*F_radx_nux(imax-1,:,:) &
                   - F_radx_nux(imax-2,:,:)
              F_rady_nux(imax,:,:) = TWO*F_rady_nux(imax-1,:,:) &
                   - F_rady_nux(imax-2,:,:)
              F_radz_nux(imax,:,:) = TWO*F_radz_nux(imax-1,:,:) &
                   - F_radz_nux(imax-2,:,:)
           end if
        end if
        Y_e(imax,:,:) = TWO*Y_e(imax-1,:,:) &
             - Y_e(imax-2,:,:)
        T_fluid(imax,:,:) = TWO*T_fluid(imax-1,:,:) &
             - T_fluid(imax-2,:,:)
     end if


     if (have_bdry_max(2)==1) then
        if (rad_evolve_enable==1) then
           E_rad(:,jmax,:) = TWO*E_rad(:,jmax-1,:) &
                - E_rad(:,jmax-2,:)
           F_radx(:,jmax,:) = TWO*F_radx(:,jmax-1,:) &
                - F_radx(:,jmax-2,:)
           F_rady(:,jmax,:) = TWO*F_rady(:,jmax-1,:) &
                - F_rady(:,jmax-2,:)
           F_radz(:,jmax,:) = TWO*F_radz(:,jmax-1,:) &
                - F_radz(:,jmax-2,:)

           if (microphysics_scheme==1) then
              E_rad_nue(:,jmax,:) = TWO*E_rad_nue(:,jmax-1,:) &
                   - E_rad_nue(:,jmax-2,:)
              F_radx_nue(:,jmax,:) = TWO*F_radx_nue(:,jmax-1,:) &
                   - F_radx_nue(:,jmax-2,:)
              F_rady_nue(:,jmax,:) = TWO*F_rady_nue(:,jmax-1,:) &
                   - F_rady_nue(:,jmax-2,:)
              F_radz_nue(:,jmax,:) = TWO*F_radz_nue(:,jmax-1,:) &
                   - F_radz_nue(:,jmax-2,:)

              E_rad_nux(:,jmax,:) = TWO*E_rad_nux(:,jmax-1,:) &
                   - E_rad_nux(:,jmax-2,:)
              F_radx_nux(:,jmax,:) = TWO*F_radx_nux(:,jmax-1,:) &
                   - F_radx_nux(:,jmax-2,:)
              F_rady_nux(:,jmax,:) = TWO*F_rady_nux(:,jmax-1,:) &
                   - F_rady_nux(:,jmax-2,:)
              F_radz_nux(:,jmax,:) = TWO*F_radz_nux(:,jmax-1,:) &
                   - F_radz_nux(:,jmax-2,:)
           end if
        end if
        Y_e(:,jmax,:) = TWO*Y_e(:,jmax-1,:) &
             - Y_e(:,jmax-2,:)
        T_fluid(:,jmax,:) = TWO*T_fluid(:,jmax-1,:) &
             - T_fluid(:,jmax-2,:)
     end if


     if (have_bdry_max(3)==1) then
        if (rad_evolve_enable==1) then
           E_rad(:,:,kmax) = TWO*E_rad(:,:,kmax-1) &
                - E_rad(:,:,kmax-2)
           F_radx(:,:,kmax) = TWO*F_radx(:,:,kmax-1) &
                - F_radx(:,:,kmax-2)
           F_rady(:,:,kmax) = TWO*F_rady(:,:,kmax-1) &
                - F_rady(:,:,kmax-2)
           F_radz(:,:,kmax) = TWO*F_radz(:,:,kmax-1) &
                - F_radz(:,:,kmax-2)

           if (microphysics_scheme==1) then
              E_rad_nue(:,:,kmax) = TWO*E_rad_nue(:,:,kmax-1) &
                   - E_rad_nue(:,:,kmax-2)
              F_radx_nue(:,:,kmax) = TWO*F_radx_nue(:,:,kmax-1) &
                   - F_radx_nue(:,:,kmax-2)
              F_rady_nue(:,:,kmax) = TWO*F_rady_nue(:,:,kmax-1) &
                   - F_rady_nue(:,:,kmax-2)
              F_radz_nue(:,:,kmax) = TWO*F_radz_nue(:,:,kmax-1) &
                   - F_radz_nue(:,:,kmax-2)

              E_rad_nux(:,:,kmax) = TWO*E_rad_nux(:,:,kmax-1) &
                   - E_rad_nux(:,:,kmax-2)
              F_radx_nux(:,:,kmax) = TWO*F_radx_nux(:,:,kmax-1) &
                   - F_radx_nux(:,:,kmax-2)
              F_rady_nux(:,:,kmax) = TWO*F_rady_nux(:,:,kmax-1) &
                   - F_rady_nux(:,:,kmax-2)
              F_radz_nux(:,:,kmax) = TWO*F_radz_nux(:,:,kmax-1) &
                   - F_radz_nux(:,:,kmax-2)
           end if
        end if
        Y_e(:,:,kmax) = TWO*Y_e(:,:,kmax-1) &
             - Y_e(:,:,kmax-2)
        T_fluid(:,:,kmax) = TWO*T_fluid(:,:,kmax-1) &
             - T_fluid(:,:,kmax-2)
     end if


     if (have_bdry_min(1)==1) then
        if (rad_evolve_enable==1) then
           E_rad(imin,:,:) = TWO*E_rad(imin+1,:,:) &
                - E_rad(imin+2,:,:)
           F_radx(imin,:,:) = TWO*F_radx(imin+1,:,:) &
                - F_radx(imin+2,:,:)
           F_rady(imin,:,:) = TWO*F_rady(imin+1,:,:) &
                - F_rady(imin+2,:,:)
           F_radz(imin,:,:) = TWO*F_radz(imin+1,:,:) &
                - F_radz(imin+2,:,:)

           if (microphysics_scheme==1) then
              E_rad_nue(imin,:,:) = TWO*E_rad_nue(imin+1,:,:) &
                   - E_rad_nue(imin+2,:,:)
              F_radx_nue(imin,:,:) = TWO*F_radx_nue(imin+1,:,:) &
                   - F_radx_nue(imin+2,:,:)
              F_rady_nue(imin,:,:) = TWO*F_rady_nue(imin+1,:,:) &
                   - F_rady_nue(imin+2,:,:)
              F_radz_nue(imin,:,:) = TWO*F_radz_nue(imin+1,:,:) &
                   - F_radz_nue(imin+2,:,:)
              
              E_rad_nux(imin,:,:) = TWO*E_rad_nux(imin+1,:,:) &
                   - E_rad_nux(imin+2,:,:)
              F_radx_nux(imin,:,:) = TWO*F_radx_nux(imin+1,:,:) &
                   - F_radx_nux(imin+2,:,:)
              F_rady_nux(imin,:,:) = TWO*F_rady_nux(imin+1,:,:) &
                   - F_rady_nux(imin+2,:,:)
              F_radz_nux(imin,:,:) = TWO*F_radz_nux(imin+1,:,:) &
                   - F_radz_nux(imin+2,:,:)
            end if
         end if
        Y_e(imin,:,:) = TWO*Y_e(imin+1,:,:) &
             - Y_e(imin+2,:,:)
        T_fluid(imin,:,:) = TWO*T_fluid(imin+1,:,:) &
             - T_fluid(imin+2,:,:)
     end if

     if (have_bdry_min(2)==1) then
       if (rad_evolve_enable==1) then
          E_rad(:,jmin,:) = TWO*E_rad(:,jmin+1,:) &
               - E_rad(:,jmin+2,:)
          F_radx(:,jmin,:) = TWO*F_radx(:,jmin+1,:) &
               - F_radx(:,jmin+2,:)
          F_rady(:,jmin,:) = TWO*F_rady(:,jmin+1,:) &
               - F_rady(:,jmin+2,:)
          F_radz(:,jmin,:) = TWO*F_radz(:,jmin+1,:) &
               - F_radz(:,jmin+2,:)

          if (microphysics_scheme==1) then
             E_rad_nue(:,jmin,:) = TWO*E_rad_nue(:,jmin+1,:) &
                  - E_rad_nue(:,jmin+2,:)
             F_radx_nue(:,jmin,:) = TWO*F_radx_nue(:,jmin+1,:) &
                  - F_radx_nue(:,jmin+2,:)
             F_rady_nue(:,jmin,:) = TWO*F_rady_nue(:,jmin+1,:) &
                  - F_rady_nue(:,jmin+2,:)
             F_radz_nue(:,jmin,:) = TWO*F_radz_nue(:,jmin+1,:) &
                  - F_radz_nue(:,jmin+2,:)

             E_rad_nux(:,jmin,:) = TWO*E_rad_nux(:,jmin+1,:) &
                  - E_rad_nux(:,jmin+2,:)
             F_radx_nux(:,jmin,:) = TWO*F_radx_nux(:,jmin+1,:) &
                  - F_radx_nux(:,jmin+2,:)
             F_rady_nux(:,jmin,:) = TWO*F_rady_nux(:,jmin+1,:) &
                  - F_rady_nux(:,jmin+2,:)
             F_radz_nux(:,jmin,:) = TWO*F_radz_nux(:,jmin+1,:) &
                  - F_radz_nux(:,jmin+2,:)
          end if
       end if
       Y_e(:,jmin,:) = TWO*Y_e(:,jmin+1,:) &
            - Y_e(:,jmin+2,:)
       T_fluid(:,jmin,:) = TWO*T_fluid(:,jmin+1,:) &
            - T_fluid(:,jmin+2,:)
    end if
        
    if (have_bdry_min(3)==1) then
       if (rad_evolve_enable==1) then    
          E_rad(:,:,kmin) = TWO*E_rad(:,:,kmin+1) &
               - E_rad(:,:,kmin+2)
          F_radx(:,:,kmin) = TWO*F_radx(:,:,kmin+1) &
               - F_radx(:,:,kmin+2)
          F_rady(:,:,kmin) = TWO*F_rady(:,:,kmin+1) &
               - F_rady(:,:,kmin+2)
          F_radz(:,:,kmin) = TWO*F_radz(:,:,kmin+1) &
               - F_radz(:,:,kmin+2)

          if (microphysics_scheme==1) then
             E_rad_nue(:,:,kmin) = TWO*E_rad_nue(:,:,kmin+1) &
                  - E_rad_nue(:,:,kmin+2)
             F_radx_nue(:,:,kmin) = TWO*F_radx_nue(:,:,kmin+1) &
                  - F_radx_nue(:,:,kmin+2)
             F_rady_nue(:,:,kmin) = TWO*F_rady_nue(:,:,kmin+1) &
                  - F_rady_nue(:,:,kmin+2)
             F_radz_nue(:,:,kmin) = TWO*F_radz_nue(:,:,kmin+1) &
                  - F_radz_nue(:,:,kmin+2)

             E_rad_nux(:,:,kmin) = TWO*E_rad_nux(:,:,kmin+1) &
                  - E_rad_nux(:,:,kmin+2)
             F_radx_nux(:,:,kmin) = TWO*F_radx_nux(:,:,kmin+1) &
                  - F_radx_nux(:,:,kmin+2)
             F_rady_nux(:,:,kmin) = TWO*F_rady_nux(:,:,kmin+1) &
                  - F_rady_nux(:,:,kmin+2)
             F_radz_nux(:,:,kmin) = TWO*F_radz_nux(:,:,kmin+1) &
                  - F_radz_nux(:,:,kmin+2) 
          end if
       end if
       Y_e(:,:,kmin) = TWO*Y_e(:,:,kmin+1) &
            - Y_e(:,:,kmin+2)
       T_fluid(:,:,kmin) = TWO*T_fluid(:,:,kmin+1) &
            - T_fluid(:,:,kmin+2)
    end if

     else if(bc==PLANAR) then
        do k=kmin+1,kmax-1
           do j=jmin,jmax
              do i=imin,imax
                 if(i==imax) then
                    if(j==jmax) then
                       rho_b(i,j,k) = rho_b(i-1,j-1,k)
                       P(i,j,k) = P(i-1,j-1,k)
                       vx(i,j,k) = vx(i-1,j-1,k)
                       vy(i,j,k) = vy(i-1,j-1,k)
                       vz(i,j,k) = vz(i-1,j-1,k)
                    else if(j==jmax-1) then
                       rho_b(i,j,k) = rho_b(i-1,j,k)
                       P(i,j,k) = P(i-1,j,k)
                       vx(i,j,k) = vx(i-1,j,k)
                       vy(i,j,k) = vy(i-1,j,k)
                       vz(i,j,k) = vz(i-1,j,k)
                    else
                       rho_b(i,j,k) = rho_b(i-1,j+1,k)
                       P(i,j,k) = P(i-1,j+1,k)
                       vx(i,j,k) = vx(i-1,j+1,k)
                       vy(i,j,k) = vy(i-1,j+1,k)
                       vz(i,j,k) = vz(i-1,j+1,k)
                    end if
                 else if(i==imin) then
                    if(j==jmin) then
                       rho_b(i,j,k) = rho_b(i+1,j+1,k)
                       P(i,j,k) = P(i+1,j+1,k)
                       vx(i,j,k) = vx(i+1,j+1,k)
                       vy(i,j,k) = vy(i+1,j+1,k)
                       vz(i,j,k) = vz(i+1,j+1,k)
                    else if(j==jmin+1) then
                       rho_b(i,j,k) = rho_b(i+1,j,k)
                       P(i,j,k) = P(i+1,j,k)
                       vx(i,j,k) = vx(i+1,j,k)
                       vy(i,j,k) = vy(i+1,j,k)
                       vz(i,j,k) = vz(i+1,j,k)
                    else
                       rho_b(i,j,k) = rho_b(i+1,j-1,k)
                       P(i,j,k) = P(i+1,j-1,k)
                       vx(i,j,k) = vx(i+1,j-1,k)
                       vy(i,j,k) = vy(i+1,j-1,k)
                       vz(i,j,k) = vz(i+1,j-1,k)
                    end if
                 else if(j==jmin) then
                    if(i==imin+1) then
                       rho_b(i,j,k) = rho_b(i,j+1,k)
                       P(i,j,k) = P(i,j+1,k)
                       vx(i,j,k) = vx(i,j+1,k)
                       vy(i,j,k) = vy(i,j+1,k)
                       vz(i,j,k) = vz(i,j+1,k)
                    else
                       rho_b(i,j,k) = rho_b(i-1,j+1,k)
                       P(i,j,k) = P(i-1,j+1,k)
                       vx(i,j,k) = vx(i-1,j+1,k)
                       vy(i,j,k) = vy(i-1,j+1,k)
                       vz(i,j,k) = vz(i-1,j+1,k)
                    end if
                 else if(j==jmax) then
                    if(i==imax-1) then
                       rho_b(i,j,k) = rho_b(i,j-1,k)
                       P(i,j,k) = P(i,j-1,k)
                       vx(i,j,k) = vx(i,j-1,k)
                       vy(i,j,k) = vy(i,j-1,k)
                       vz(i,j,k) = vz(i,j-1,k)
                    else
                       rho_b(i,j,k) = rho_b(i+1,j-1,k)
                       P(i,j,k) = P(i+1,j-1,k)
                       vx(i,j,k) = vx(i+1,j-1,k)
                       vy(i,j,k) = vy(i+1,j-1,k)
                       vz(i,j,k) = vz(i+1,j-1,k)
                    end if
              end if
              if(i==imax) then
                 if(j==jmax) then
                    if (rad_evolve_enable==1) then
                       E_rad(i,j,k) = E_rad(i-1,j-1,k)
                       F_radx(i,j,k) = F_radx(i-1,j-1,k)
                       F_rady(i,j,k) = F_rady(i-1,j-1,k)
                       F_radz(i,j,k) = F_radz(i-1,j-1,k)
                       if (microphysics_scheme==1) then
                          E_rad_nue(i,j,k) = E_rad_nue(i-1,j-1,k)
                          F_radx_nue(i,j,k) = F_radx_nue(i-1,j-1,k)
                          F_rady_nue(i,j,k) = F_rady_nue(i-1,j-1,k)
                          F_radz_nue(i,j,k) = F_radz_nue(i-1,j-1,k)

                          E_rad_nux(i,j,k) = E_rad_nux(i-1,j-1,k)
                          F_radx_nux(i,j,k) = F_radx_nux(i-1,j-1,k)
                          F_rady_nux(i,j,k) = F_rady_nux(i-1,j-1,k)
                          F_radz_nux(i,j,k) = F_radz_nux(i-1,j-1,k)
                       end if
                    end if
                    Y_e(i,j,k) = Y_e(i-1,j-1,k)
                    T_fluid(i,j,k) = T_fluid(i-1,j-1,k)
                 else if(j==jmax-1) then
                    if (rad_evolve_enable==1) then
                       E_rad(i,j,k) = E_rad(i-1,j,k)
                       F_radx(i,j,k) = F_radx(i-1,j,k)
                       F_rady(i,j,k) = F_rady(i-1,j,k)
                       F_radz(i,j,k) = F_radz(i-1,j,k)
                       if (microphysics_scheme==1) then
                          E_rad_nue(i,j,k) = E_rad_nue(i-1,j,k)
                          F_radx_nue(i,j,k) = F_radx_nue(i-1,j,k)
                          F_rady_nue(i,j,k) = F_rady_nue(i-1,j,k)
                          F_radz_nue(i,j,k) = F_radz_nue(i-1,j,k)

                          E_rad_nux(i,j,k) = E_rad_nux(i-1,j,k)
                          F_radx_nux(i,j,k) = F_radx_nux(i-1,j,k)
                          F_rady_nux(i,j,k) = F_rady_nux(i-1,j,k)
                          F_radz_nux(i,j,k) = F_radz_nux(i-1,j,k)
                       end if
                    end if
                    Y_e(i,j,k) = Y_e(i-1,j,k)
                    T_fluid(i,j,k) = T_fluid(i-1,j,k)
                 else
                    if (rad_evolve_enable==1) then
                       E_rad(i,j,k) = E_rad(i-1,j+1,k)
                       F_radx(i,j,k) = F_radx(i-1,j+1,k)
                       F_rady(i,j,k) = F_rady(i-1,j+1,k)
                       F_radz(i,j,k) = F_radz(i-1,j+1,k)
                       if (microphysics_scheme==1) then
                          E_rad_nue(i,j,k) = E_rad_nue(i-1,j+1,k)
                          F_radx_nue(i,j,k) = F_radx_nue(i-1,j+1,k)
                          F_rady_nue(i,j,k) = F_rady_nue(i-1,j+1,k)
                          F_radz_nue(i,j,k) = F_radz_nue(i-1,j+1,k)

                          E_rad_nux(i,j,k) = E_rad_nux(i-1,j+1,k)
                          F_radx_nux(i,j,k) = F_radx_nux(i-1,j+1,k)
                          F_rady_nux(i,j,k) = F_rady_nux(i-1,j+1,k)
                          F_radz_nux(i,j,k) = F_radz_nux(i-1,j+1,k)
                       end if
                    end if
                    Y_e(i,j,k) = Y_e(i-1,j+1,k)
                    T_fluid(i,j,k) = T_fluid(i-1,j+1,k)
                 end if
             else if(i==imin) then
                 if(j==jmin) then
                    if (rad_evolve_enable==1) then
                       E_rad(i,j,k) = E_rad(i+1,j+1,k)
                       F_radx(i,j,k) = F_radx(i+1,j+1,k)
                       F_rady(i,j,k) = F_rady(i+1,j+1,k)
                       F_radz(i,j,k) = F_radz(i+1,j+1,k)
                       if (microphysics_scheme==1) then
                          E_rad_nue(i,j,k) = E_rad_nue(i+1,j+1,k)
                          F_radx_nue(i,j,k) = F_radx_nue(i+1,j+1,k)
                          F_rady_nue(i,j,k) = F_rady_nue(i+1,j+1,k)
                          F_radz_nue(i,j,k) = F_radz_nue(i+1,j+1,k)
                          
                          E_rad_nux(i,j,k) = E_rad_nux(i+1,j+1,k)
                          F_radx_nux(i,j,k) = F_radx_nux(i+1,j+1,k)
                          F_rady_nux(i,j,k) = F_rady_nux(i+1,j+1,k)
                          F_radz_nux(i,j,k) = F_radz_nux(i+1,j+1,k)
                       end if
                    end if
                    Y_e(i,j,k) = Y_e(i+1,j+1,k)
                    T_fluid(i,j,k) = T_fluid(i+1,j+1,k)
                 else if(j==jmin+1) then
                    if (rad_evolve_enable==1) then
                       E_rad(i,j,k) = E_rad(i+1,j,k)
                       F_radx(i,j,k) = F_radx(i+1,j,k)
                       F_rady(i,j,k) = F_rady(i+1,j,k)
                       F_radz(i,j,k) = F_radz(i+1,j,k)
                       if (microphysics_scheme==1) then
                          E_rad_nue(i,j,k) = E_rad_nue(i+1,j,k)
                          F_radx_nue(i,j,k) = F_radx_nue(i+1,j,k)
                          F_rady_nue(i,j,k) = F_rady_nue(i+1,j,k)
                          F_radz_nue(i,j,k) = F_radz_nue(i+1,j,k)

                          E_rad_nux(i,j,k) = E_rad_nux(i+1,j,k)
                          F_radx_nux(i,j,k) = F_radx_nux(i+1,j,k)
                          F_rady_nux(i,j,k) = F_rady_nux(i+1,j,k)
                          F_radz_nux(i,j,k) = F_radz_nux(i+1,j,k)
                       end if
                    end if
                    Y_e(i,j,k) = Y_e(i+1,j,k)
                    T_fluid(i,j,k) = T_fluid(i+1,j,k)
                 else
                    if (rad_evolve_enable==1) then
                       E_rad(i,j,k) = E_rad(i+1,j-1,k)
                       F_radx(i,j,k) = F_radx(i+1,j-1,k)
                       F_rady(i,j,k) = F_rady(i+1,j-1,k)
                       F_radz(i,j,k) = F_radz(i+1,j-1,k)
                       if (microphysics_scheme==1) then
                          E_rad_nue(i,j,k) = E_rad_nue(i+1,j-1,k)
                          F_radx_nue(i,j,k) = F_radx_nue(i+1,j-1,k)
                          F_rady_nue(i,j,k) = F_rady_nue(i+1,j-1,k)
                          F_radz_nue(i,j,k) = F_radz_nue(i+1,j-1,k)
                          
                          E_rad_nux(i,j,k) = E_rad_nux(i+1,j-1,k)
                          F_radx_nux(i,j,k) = F_radx_nux(i+1,j-1,k)
                          F_rady_nux(i,j,k) = F_rady_nux(i+1,j-1,k)
                          F_radz_nux(i,j,k) = F_radz_nux(i+1,j-1,k)
                       end if
                    end if
                    Y_e(i,j,k) = Y_e(i+1,j-1,k)
                    T_fluid(i,j,k) = T_fluid(i+1,j-1,k)
                 end if
              else if(j==jmin) then
                 if(i==imin+1) then
                    if (rad_evolve_enable==1) then
                       E_rad(i,j,k) = E_rad(i,j+1,k)
                       F_radx(i,j,k) = F_radx(i,j+1,k)
                       F_rady(i,j,k) = F_rady(i,j+1,k)
                       F_radz(i,j,k) = F_radz(i,j+1,k)
                       if (microphysics_scheme==1) then
                          E_rad_nue(i,j,k) = E_rad_nue(i,j+1,k)
                          F_radx_nue(i,j,k) = F_radx_nue(i,j+1,k)
                          F_rady_nue(i,j,k) = F_rady_nue(i,j+1,k)
                          F_radz_nue(i,j,k) = F_radz_nue(i,j+1,k)
                          
                          E_rad_nux(i,j,k) = E_rad_nux(i,j+1,k)
                          F_radx_nux(i,j,k) = F_radx_nux(i,j+1,k)
                          F_rady_nux(i,j,k) = F_rady_nux(i,j+1,k)
                          F_radz_nux(i,j,k) = F_radz_nux(i,j+1,k)
                       end if
                    end if
                    Y_e(i,j,k) = Y_e(i,j+1,k)
                    T_fluid(i,j,k) = T_fluid(i,j+1,k)
                 else
                    if (rad_evolve_enable==1) then
                       E_rad(i,j,k) = E_rad(i-1,j+1,k)
                       F_radx(i,j,k) = F_radx(i-1,j+1,k)
                       F_rady(i,j,k) = F_rady(i-1,j+1,k)
                       F_radz(i,j,k) = F_radz(i-1,j+1,k)
                       if (microphysics_scheme==1) then
                          E_rad_nue(i,j,k) = E_rad_nue(i-1,j+1,k)
                          F_radx_nue(i,j,k) = F_radx_nue(i-1,j+1,k)
                          F_rady_nue(i,j,k) = F_rady_nue(i-1,j+1,k)
                          F_radz_nue(i,j,k) = F_radz_nue(i-1,j+1,k)

                          E_rad_nux(i,j,k) = E_rad_nux(i-1,j+1,k)
                          F_radx_nux(i,j,k) = F_radx_nux(i-1,j+1,k)
                          F_rady_nux(i,j,k) = F_rady_nux(i-1,j+1,k)
                          F_radz_nux(i,j,k) = F_radz_nux(i-1,j+1,k) 
                       end if
                    end if
                    Y_e(i,j,k) = Y_e(i-1,j+1,k)
                    T_fluid(i,j,k) = T_fluid(i-1,j+1,k)
                 end if
              else if(j==jmax) then
                 if(i==imax-1) then
                    if (rad_evolve_enable==1) then
                       E_rad(i,j,k) = E_rad(i,j-1,k)
                       F_radx(i,j,k) = F_radx(i,j-1,k)
                       F_rady(i,j,k) = F_rady(i,j-1,k)
                       F_radz(i,j,k) = F_radz(i,j-1,k)
                       if (microphysics_scheme==1) then
                          E_rad_nue(i,j,k) = E_rad_nue(i,j-1,k)
                          F_radx_nue(i,j,k) = F_radx_nue(i,j-1,k)
                          F_rady_nue(i,j,k) = F_rady_nue(i,j-1,k)
                          F_radz_nue(i,j,k) = F_radz_nue(i,j-1,k)

                          E_rad_nux(i,j,k) = E_rad_nux(i,j-1,k)
                          F_radx_nux(i,j,k) = F_radx_nux(i,j-1,k)
                          F_rady_nux(i,j,k) = F_rady_nux(i,j-1,k)
                          F_radz_nux(i,j,k) = F_radz_nux(i,j-1,k)
                       end if
                    end if
                    Y_e(i,j,k) = Y_e(i,j-1,k)
                    T_fluid(i,j,k) = T_fluid(i,j-1,k)
                 else
                    if (rad_evolve_enable==1) then
                       E_rad(i,j,k) = E_rad(i+1,j-1,k)
                       F_radx(i,j,k) = F_radx(i+1,j-1,k)
                       F_rady(i,j,k) = F_rady(i+1,j-1,k)
                       F_radz(i,j,k) = F_radz(i+1,j-1,k)
                       if (microphysics_scheme==1) then
                          E_rad_nue(i,j,k) = E_rad_nue(i+1,j-1,k)
                          F_radx_nue(i,j,k) = F_radx_nue(i+1,j-1,k)
                          F_rady_nue(i,j,k) = F_rady_nue(i+1,j-1,k)
                          F_radz_nue(i,j,k) = F_radz_nue(i+1,j-1,k)

                          E_rad_nux(i,j,k) = E_rad_nux(i+1,j-1,k)
                          F_radx_nux(i,j,k) = F_radx_nux(i+1,j-1,k)
                          F_rady_nux(i,j,k) = F_rady_nux(i+1,j-1,k)
                          F_radz_nux(i,j,k) = F_radz_nux(i+1,j-1,k)
                       end if
                    end if
                    Y_e(i,j,k) = Y_e(i+1,j-1,k)
                    T_fluid(i,j,k) = T_fluid(i+1,j-1,k)
                 end if
              end if
           end do
        end do
     end do
  else if (bc==COPY) then !!!! TODO: Add microphysics_scheme==1 part of COPY and QUAD bc.
     if (have_bdry_max(1)==1) then
        rho_b(imax,:,:) = rho_b(imax-1,:,:)
        P(imax,:,:) = P(imax-1,:,:)
        vx(imax,:,:) = vx(imax-1,:,:)
        vy(imax,:,:) = vy(imax-1,:,:)
        vz(imax,:,:) = vz(imax-1,:,:)
     end if
     if (have_bdry_min(1)==1) then
        rho_b(imin,:,:) = rho_b(imin+1,:,:)
        P(imin,:,:) = P(imin+1,:,:)
        vx(imin,:,:) = vx(imin+1,:,:)
        vy(imin,:,:) = vy(imin+1,:,:)
        vz(imin,:,:) = vz(imin+1,:,:)
     end if
     if (have_bdry_max(2)==1) then
        rho_b(:,jmax,:) = rho_b(:,jmax-1,:)
        P(:,jmax,:) = P(:,jmax-1,:)
        vx(:,jmax,:) = vx(:,jmax-1,:)
        vy(:,jmax,:) = vy(:,jmax-1,:)
        vz(:,jmax,:) = vz(:,jmax-1,:)
     end if
     if (have_bdry_min(2)==1) then
        rho_b(:,jmin,:) = rho_b(:,jmin+1,:)
        P(:,jmin,:) = P(:,jmin+1,:)
        vx(:,jmin,:) = vx(:,jmin+1,:)
        vy(:,jmin,:) = vy(:,jmin+1,:)
        vz(:,jmin,:) = vz(:,jmin+1,:)
     end if
     if (have_bdry_max(3)==1) then
        rho_b(:,:,kmax) = rho_b(:,:,kmax-1)
        P(:,:,kmax) = P(:,:,kmax-1)
        vx(:,:,kmax) = vx(:,:,kmax-1)
        vy(:,:,kmax) = vy(:,:,kmax-1)
        vz(:,:,kmax) = vz(:,:,kmax-1)
     end if
     if (have_bdry_min(3)==1) then
        rho_b(:,:,kmin) = rho_b(:,:,kmin+1)
        P(:,:,kmin) = P(:,:,kmin+1)
        vx(:,:,kmin) = vx(:,:,kmin+1)
        vy(:,:,kmin) = vy(:,:,kmin+1)
        vz(:,:,kmin) = vz(:,:,kmin+1)
     end if
     if (rad_evolve_enable==1) then
        if (have_bdry_max(1)==1) then
           E_rad(imax,:,:) = E_rad(imax-1,:,:)
           F_radx(imax,:,:) = F_radx(imax-1,:,:)
           F_rady(imax,:,:) = F_rady(imax-1,:,:)
           F_radz(imax,:,:) = F_radz(imax-1,:,:)
           Y_e(imax,:,:) = Y_e(imax-1,:,:)
           T_fluid(imax,:,:) = T_fluid(imax-1,:,:)
        end if
        if (have_bdry_min(1)==1) then
           E_rad(imin,:,:) = E_rad(imin+1,:,:)
           F_radx(imin,:,:) = F_radx(imin+1,:,:)
           F_rady(imin,:,:) = F_rady(imin+1,:,:)
           F_radz(imin,:,:) = F_radz(imin+1,:,:)
           Y_e(imin,:,:) = Y_e(imin+1,:,:)
           T_fluid(imin,:,:) = T_fluid(imin+1,:,:)
        end if
        if (have_bdry_max(2)==1) then
           E_rad(:,jmax,:) = E_rad(:,jmax-1,:)
           F_radx(:,jmax,:) = F_radx(:,jmax-1,:)
           F_rady(:,jmax,:) = F_rady(:,jmax-1,:)
           F_radz(:,jmax,:) = F_radz(:,jmax-1,:)
           Y_e(:,jmax,:) = Y_e(:,jmax-1,:)
           T_fluid(:,jmax,:) = T_fluid(:,jmax-1,:)
        end if
        if (have_bdry_min(2)==1) then
           E_rad(:,jmin,:) = E_rad(:,jmin+1,:)
           F_radx(:,jmin,:) = F_radx(:,jmin+1,:)
           F_rady(:,jmin,:) = F_rady(:,jmin+1,:)
           F_radz(:,jmin,:) = F_radz(:,jmin+1,:)
           Y_e(:,jmin,:) = Y_e(:,jmin+1,:)
           T_fluid(:,jmin,:) = T_fluid(:,jmin+1,:)
        end if
        if (have_bdry_max(3)==1) then
           E_rad(:,:,kmax) = E_rad(:,:,kmax-1)
           F_radx(:,:,kmax) = F_radx(:,:,kmax-1)
           F_rady(:,:,kmax) = F_rady(:,:,kmax-1)
           F_radz(:,:,kmax) = F_radz(:,:,kmax-1)
           Y_e(:,:,kmax) = Y_e(:,:,kmax-1)
           T_fluid(:,:,kmax) = T_fluid(:,:,kmax-1)
        end if
        if (have_bdry_min(3)==1) then
           E_rad(:,:,kmin) = E_rad(:,:,kmin+1)
           F_radx(:,:,kmin) = F_radx(:,:,kmin+1)
           F_rady(:,:,kmin) = F_rady(:,:,kmin+1)
           F_radz(:,:,kmin) = F_radz(:,:,kmin+1)
           Y_e(:,:,kmin) = Y_e(:,:,kmin+1)
           T_fluid(:,:,kmin) = T_fluid(:,:,kmin+1)
        end if
     end if
  else if(bc==QUAD) then
     if (have_bdry_max(1)==1) then
        rho_b(imax,:,:) = THREE*rho_b(imax-1,:,:) &
             - THREE*rho_b(imax-2,:,:) + rho_b(imax-3,:,:)
        P(imax,:,:) = THREE*P(imax-1,:,:) &
             - THREE*P(imax-2,:,:) + P(imax-3,:,:)
        vx(imax,:,:) = THREE*vx(imax-1,:,:) &
             - THREE*vx(imax-2,:,:) + vx(imax-3,:,:)
        vy(imax,:,:) = THREE*vy(imax-1,:,:) &
             - THREE*vy(imax-2,:,:) + vy(imax-3,:,:)
        vz(imax,:,:) = THREE*vz(imax-1,:,:) &
             - THREE*vz(imax-2,:,:) + vz(imax-3,:,:)
     end if
     if (have_bdry_min(1)==1) then
        rho_b(imin,:,:) = THREE*rho_b(imin+1,:,:) &
             - THREE*rho_b(imin+2,:,:) + rho_b(imin+3,:,:)
        P(imin,:,:) = THREE*P(imin+1,:,:) &
             - THREE*P(imin+2,:,:) + P(imin+3,:,:)
        vx(imin,:,:) = THREE*vx(imin+1,:,:) &
             - THREE*vx(imin+2,:,:) + vx(imin+3,:,:)
        vy(imin,:,:) = THREE*vy(imin+1,:,:) &
             - THREE*vy(imin+2,:,:) + vy(imin+3,:,:)
        vz(imin,:,:) = THREE*vz(imin+1,:,:) &
             - THREE*vz(imin+2,:,:) + vz(imin+3,:,:)
     end if
     if (have_bdry_max(2)==1) then
        rho_b(:,jmax,:) = THREE*rho_b(:,jmax-1,:) &
             - THREE*rho_b(:,jmax-2,:) + rho_b(:,jmax-3,:)
        P(:,jmax,:) = THREE*P(:,jmax-1,:) &
             - THREE*P(:,jmax-2,:) + P(:,jmax-3,:)
        vx(:,jmax,:) = THREE*vx(:,jmax-1,:) &
             - THREE*vx(:,jmax-2,:) + vx(:,jmax-3,:)
        vy(:,jmax,:) = THREE*vy(:,jmax-1,:) &
             - THREE*vy(:,jmax-2,:) + vy(:,jmax-3,:)
        vz(:,jmax,:) = THREE*vz(:,jmax-1,:) &
             - THREE*vz(:,jmax-2,:) + vz(:,jmax-3,:)
     end if
     if (have_bdry_min(2)==1) then
        rho_b(:,jmin,:) = THREE*rho_b(:,jmin+1,:) &
             - THREE*rho_b(:,jmin+2,:) + rho_b(:,jmin+3,:)
        P(:,jmin,:) = THREE*P(:,jmin+1,:) &
             - THREE*P(:,jmin+2,:) + P(:,jmin+3,:)
        vx(:,jmin,:) = THREE*vx(:,jmin+1,:) &
             - THREE*vx(:,jmin+2,:) + vx(:,jmin+3,:)
        vy(:,jmin,:) = THREE*vy(:,jmin+1,:) &
             - THREE*vy(:,jmin+2,:) + vy(:,jmin+3,:)
        vz(:,jmin,:) = THREE*vz(:,jmin+1,:) &
             - THREE*vz(:,jmin+2,:) + vz(:,jmin+3,:)
     end if
     if (have_bdry_max(3)==1) then
        rho_b(:,:,kmax) = THREE*rho_b(:,:,kmax-1) &
             - THREE*rho_b(:,:,kmax-2) + rho_b(:,:,kmax-3)
        P(:,:,kmax) = THREE*P(:,:,kmax-1) &
             - THREE*P(:,:,kmax-2) + P(:,:,kmax-3)
        vx(:,:,kmax) = THREE*vx(:,:,kmax-1) &
             - THREE*vx(:,:,kmax-2) + vx(:,:,kmax-3)
        vy(:,:,kmax) = THREE*vy(:,:,kmax-1) &
             - THREE*vy(:,:,kmax-2) + vy(:,:,kmax-3)
        vz(:,:,kmax) = THREE*vz(:,:,kmax-1) &
             - THREE*vz(:,:,kmax-2) + vz(:,:,kmax-3)
     end if
     if (have_bdry_min(3)==1) then
        rho_b(:,:,kmin) = THREE*rho_b(:,:,kmin+1) &
             - THREE*rho_b(:,:,kmin+2) + rho_b(:,:,kmin+3)
        P(:,:,kmin) = THREE*P(:,:,kmin+1) &
             - THREE*P(:,:,kmin+2) + P(:,:,kmin+3)
        vx(:,:,kmin) = THREE*vx(:,:,kmin+1) &
             - THREE*vx(:,:,kmin+2) + vx(:,:,kmin+3)
        vy(:,:,kmin) = THREE*vy(:,:,kmin+1) &
             - THREE*vy(:,:,kmin+2) + vy(:,:,kmin+3)
        vz(:,:,kmin) = THREE*vz(:,:,kmin+1) &
             - THREE*vz(:,:,kmin+2) + vz(:,:,kmin+3)
     end if
     if (rad_evolve_enable==1) then
        if (have_bdry_max(1)==1) then
           E_rad(imax,:,:) = THREE*E_rad(imax-1,:,:) &
                - THREE*E_rad(imax-2,:,:) + E_rad(imax-3,:,:)
           F_radx(imax,:,:) = THREE*F_radx(imax-1,:,:) &
                - THREE*F_radx(imax-2,:,:) + F_radx(imax-3,:,:)
           F_rady(imax,:,:) = THREE*F_rady(imax-1,:,:) &
                - THREE*F_rady(imax-2,:,:) + F_rady(imax-3,:,:)
           F_radz(imax,:,:) = THREE*F_radz(imax-1,:,:) &
                - THREE*F_radz(imax-2,:,:) + F_radz(imax-3,:,:)
           Y_e(imax,:,:) = THREE*Y_e(imax-1,:,:) &
                - THREE*Y_e(imax-2,:,:) + Y_e(imax-3,:,:)
           T_fluid(imax,:,:) = THREE*T_fluid(imax-1,:,:) &
                - THREE*T_fluid(imax-2,:,:) + T_fluid(imax-3,:,:)
        end if
        
        if (have_bdry_min(1)==1) then
           E_rad(imin,:,:) = THREE*E_rad(imin+1,:,:) &
                - THREE*E_rad(imin+2,:,:) + E_rad(imin+3,:,:)
           F_radx(imin,:,:) = THREE*F_radx(imin+1,:,:) &
                - THREE*F_radx(imin+2,:,:) + F_radx(imin+3,:,:)
           F_rady(imin,:,:) = THREE*F_rady(imin+1,:,:) &
                - THREE*F_rady(imin+2,:,:) + F_rady(imin+3,:,:)
           F_radz(imin,:,:) = THREE*F_radz(imin+1,:,:) &
                - THREE*F_radz(imin+2,:,:) + F_radz(imin+3,:,:)
           Y_e(imin,:,:) = THREE*Y_e(imin+1,:,:) &
                - THREE*Y_e(imin+2,:,:) + Y_e(imin+3,:,:)
           T_fluid(imin,:,:) = THREE*T_fluid(imin+1,:,:) &
                - THREE*T_fluid(imin+2,:,:) + T_fluid(imin+3,:,:)
        end if
        if (have_bdry_max(2)==1) then
           E_rad(:,jmax,:) = THREE*E_rad(:,jmax-1,:) &
                - THREE*E_rad(:,jmax-2,:) + E_rad(:,jmax-3,:)
           F_radx(:,jmax,:) = THREE*F_radx(:,jmax-1,:) &
                - THREE*F_radx(:,jmax-2,:) + F_radx(:,jmax-3,:)
           F_rady(:,jmax,:) = THREE*F_rady(:,jmax-1,:) &
                - THREE*F_rady(:,jmax-2,:) + F_rady(:,jmax-3,:)
           F_radz(:,jmax,:) = THREE*F_radz(:,jmax-1,:) &
                - THREE*F_radz(:,jmax-2,:) + F_radz(:,jmax-3,:)
           Y_e(:,jmax,:) = THREE*Y_e(:,jmax-1,:) &
                - THREE*Y_e(:,jmax-2,:) + Y_e(:,jmax-3,:)
           T_fluid(:,jmax,:) = THREE*T_fluid(:,jmax-1,:) &
                - THREE*T_fluid(:,jmax-2,:) + T_fluid(:,jmax-3,:)
        end if
        if (have_bdry_min(2)==1) then
           E_rad(:,jmin,:) = THREE*E_rad(:,jmin+1,:) &
                - THREE*E_rad(:,jmin+2,:) + E_rad(:,jmin+3,:)
           F_radx(:,jmin,:) = THREE*F_radx(:,jmin+1,:) &
                - THREE*F_radx(:,jmin+2,:) + F_radx(:,jmin+3,:)
           F_rady(:,jmin,:) = THREE*F_rady(:,jmin+1,:) &
                - THREE*F_rady(:,jmin+2,:) + F_rady(:,jmin+3,:)
           F_radz(:,jmin,:) = THREE*F_radz(:,jmin+1,:) &
                - THREE*F_radz(:,jmin+2,:) + F_radz(:,jmin+3,:)
           Y_e(:,jmin,:) = THREE*Y_e(:,jmin+1,:) &
                - THREE*Y_e(:,jmin+2,:) + Y_e(:,jmin+3,:)
           T_fluid(:,jmin,:) = THREE*T_fluid(:,jmin+1,:) &
                - THREE*T_fluid(:,jmin+2,:) + T_fluid(:,jmin+3,:)
        end if
        if (have_bdry_max(3)==1) then
           E_rad(:,:,kmax) = THREE*E_rad(:,:,kmax-1) &
                - THREE*E_rad(:,:,kmax-2) + E_rad(:,:,kmax-3)
           F_radx(:,:,kmax) = THREE*F_radx(:,:,kmax-1) &
                - THREE*F_radx(:,:,kmax-2) + F_radx(:,:,kmax-3)
           F_rady(:,:,kmax) = THREE*F_rady(:,:,kmax-1) &
                - THREE*F_rady(:,:,kmax-2) + F_rady(:,:,kmax-3)
           F_radz(:,:,kmax) = THREE*F_radz(:,:,kmax-1) &
                - THREE*F_radz(:,:,kmax-2) + F_radz(:,:,kmax-3)
           Y_e(:,:,kmax) = THREE*Y_e(:,:,kmax-1) &
                - THREE*Y_e(:,:,kmax-2) + Y_e(:,:,kmax-3)
           T_fluid(:,:,kmax) = THREE*T_fluid(:,:,kmax-1) &
                - THREE*T_fluid(:,:,kmax-2) + T_fluid(:,:,kmax-3)
        end if
        if (have_bdry_min(3)==1) then
           E_rad(:,:,kmin) = THREE*E_rad(:,:,kmin+1) &
                - THREE*E_rad(:,:,kmin+2) + E_rad(:,:,kmin+3)
           F_radx(:,:,kmin) = THREE*F_radx(:,:,kmin+1) &
                - THREE*F_radx(:,:,kmin+2) + F_radx(:,:,kmin+3)
           F_rady(:,:,kmin) = THREE*F_rady(:,:,kmin+1) &
                - THREE*F_rady(:,:,kmin+2) + F_rady(:,:,kmin+3)
           F_radz(:,:,kmin) = THREE*F_radz(:,:,kmin+1) &
                - THREE*F_radz(:,:,kmin+2) + F_radz(:,:,kmin+3)
           Y_e(:,:,kmin) = THREE*Y_e(:,:,kmin+1) &
                - THREE*Y_e(:,:,kmin+2) + Y_e(:,:,kmin+3)
           T_fluid(:,:,kmin) = THREE*T_fluid(:,:,kmin+1) &
                - THREE*T_fluid(:,:,kmin+2) + T_fluid(:,:,kmin+3)
        end if
     end if
  end if
     !write(*,*) midway bc rho_star = ,rho_b(15,15,5),P(15,15,5),vx(15,15,5),vy(15,15,5),vz(15,15,5),have_bdry_max(3)==1,proc_kmin,Z(1,
     ! now recompute the conserved variables and the auxiliary variables.
  ! i=imax, k=kmax
  !!!! Lunan: before this step, F_rad0 should also be updated at BC? otherwise it would be pass as zero at boundaries?


  if (Symmetry == AXISYM) then
     if (have_bdry_max(1)==1) then
        i=imax
        j=2
        do k=kmin,kmax
           rho_bl = rho_b(i,j,k)
           Pl = P(i,j,k)
           vxl = vx(i,j,k)
           vyl = vy(i,j,k)
           vzl = vz(i,j,k)
           Bxl = Bx(i,j,k)
           Byl = By(i,j,k)
           Bzl = Bz(i,j,k)
           E_radl = E_rad(i,j,k)
           F_radxl = F_radx(i,j,k)
           F_radyl = F_rady(i,j,k)
           F_radzl = F_radz(i,j,k)
           E_rad_nuel = E_rad_nue(i,j,k)
           F_radx_nuel = F_radx_nue(i,j,k)
           F_rady_nuel = F_rady_nue(i,j,k)
           F_radz_nuel = F_radz_nue(i,j,k)
           E_rad_nuxl = E_rad_nux(i,j,k)
           F_radx_nuxl = F_radx_nux(i,j,k)
           F_rady_nuxl = F_rady_nux(i,j,k)
           F_radz_nuxl = F_radz_nux(i,j,k)
           Y_el    = Y_e(i,j,k)
           AUX(1)  = phi(i,j,k)
           AUX(2)  = lapse(i,j,k)
           AUX(3)  = shiftx(i,j,k)
           AUX(4)  = shifty(i,j,k)
           AUX(5)  = shiftz(i,j,k)
           AUX(6)  = gxx(i,j,k)
           AUX(7)  = gxy(i,j,k)
           AUX(8)  = gxz(i,j,k)
           AUX(9)  = gyy(i,j,k)
           AUX(10) = gyz(i,j,k)
           AUX(11) = gzz(i,j,k)
           AUX(12) = gupxx(i,j,k)
           AUX(13) = gupxy(i,j,k)
           AUX(14) = gupxz(i,j,k)
           AUX(15) = gupyy(i,j,k)
           AUX(16) = gupyz(i,j,k)
           AUX(17) = gupzz(i,j,k)
           call recompute_vars(rho_bl,Pl,vxl,vyl,vzl,Bxl,Byl,Bzl,AUX,m,rho_sl,taul, &
                mhd_st_xl,mhd_st_yl,mhd_st_zl,st_xl,st_yl,st_zl,wl,hl,u0l,Exl,Eyl,Ezl,&
                E_radl,F_radxl,F_radyl,F_radzl,Y_el, T_fluidl,tau_radl, S_rad_xl, S_rad_yl, S_rad_zl,rho_Yel,&
                F_rad0l, P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl, &
                E_rad_nuel,F_radx_nuel,F_rady_nuel,F_radz_nuel, tau_rad_nuel, S_rad_x_nuel, S_rad_y_nuel, S_rad_z_nuel, &
                F_rad0_nuel, P_radxx_nuel, P_radyy_nuel, P_radzz_nuel, P_radxy_nuel, P_radxz_nuel, P_radyz_nuel, &
                E_rad_nuxl,F_radx_nuxl,F_rady_nuxl,F_radz_nuxl, tau_rad_nuxl, S_rad_x_nuxl, S_rad_y_nuxl, S_rad_z_nuxl, &
                F_rad0_nuxl, P_radxx_nuxl, P_radyy_nuxl, P_radzz_nuxl, P_radxy_nuxl, P_radxz_nuxl, P_radyz_nuxl, &
		eps_totl, eps_cldl, eps_thermall, P_cldl, &
                sbtl,sbxl,sbyl,sbzl,rhol,Sxl,Syl,Szl,Sxxl,Sxyl,Sxzl,Syyl,Syzl,Szzl, &
                neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, &
		gamma_tab, gamma_th,rho_b_atm, rad_const, &
     		compute_microphysics, enable_OS_collapse,rad_evolve_enable, rad_closure_scheme, microphysics_scheme, Erad_atm_cut)
           ! Assign primitive variables which may have been modified
           vx(i,j,k) = vxl
           vy(i,j,k) = vyl
           vz(i,j,k) = vzl
           E_rad(i,j,k) = E_radl
           E_rad_nue(i,j,k) = E_rad_nuel
           E_rad_nux(i,j,k) = E_rad_nuxl
           ! Assign conserved variables
           rho_star(i,j,k) = rho_sl
           tau(i,j,k) = taul
           mhd_st_x(i,j,k) = mhd_st_xl
           mhd_st_y(i,j,k) = mhd_st_yl
           mhd_st_z(i,j,k) = mhd_st_zl
           tau_rad(i,j,k) = tau_radl
           S_rad_x(i,j,k) = S_rad_xl
           S_rad_y(i,j,k) = S_rad_yl
           S_rad_z(i,j,k) = S_rad_zl
           tau_rad_nue(i,j,k) = tau_rad_nuel
           S_rad_x_nue(i,j,k) = S_rad_x_nuel
           S_rad_y_nue(i,j,k) = S_rad_y_nuel
           S_rad_z_nue(i,j,k) = S_rad_z_nuel
           tau_rad_nux(i,j,k) = tau_rad_nuxl
           S_rad_x_nux(i,j,k) = S_rad_x_nuxl
           S_rad_y_nux(i,j,k) = S_rad_y_nuxl
           S_rad_z_nux(i,j,k) = S_rad_z_nuxl
           rho_Ye(i,j,k) = rho_Yel
           ! Assign other variables
           st_x(i,j,k) = st_xl
           st_y(i,j,k) = st_yl
           st_z(i,j,k) = st_zl
           w(i,j,k) = wl
           h(i,j,k) = hl
           u0(i,j,k) = u0l
           Ex(i,j,k) = Exl
           Ey(i,j,k) = Eyl
           Ez(i,j,k) = Ezl
           sbt(i,j,k) = sbtl
           sbx(i,j,k) = sbxl
           sby(i,j,k) = sbyl
           sbz(i,j,k) = sbzl
           rho(i,j,k) = rhol
           F_rad0(i,j,k) = F_rad0l
           P_radxx(i,j,k) = P_radxxl
           P_radyy(i,j,k) = P_radyyl
           P_radzz(i,j,k) = P_radzzl
           P_radxy(i,j,k) = P_radxyl
           P_radxz(i,j,k) = P_radxzl
           P_radyz(i,j,k) = P_radyzl
           F_rad0_nue(i,j,k) = F_rad0_nuel
           P_radxx_nue(i,j,k) = P_radxx_nuel
           P_radyy_nue(i,j,k) = P_radyy_nuel
           P_radzz_nue(i,j,k) = P_radzz_nuel
           P_radxy_nue(i,j,k) = P_radxy_nuel
           P_radxz_nue(i,j,k) = P_radxz_nuel
           P_radyz_nue(i,j,k) = P_radyz_nuel
           F_rad0_nux(i,j,k) = F_rad0_nuxl
           P_radxx_nux(i,j,k) = P_radxx_nuxl
           P_radyy_nux(i,j,k) = P_radyy_nuxl
           P_radzz_nux(i,j,k) = P_radzz_nuxl
           P_radxy_nux(i,j,k) = P_radxy_nuxl
           P_radxz_nux(i,j,k) = P_radxz_nuxl
           P_radyz_nux(i,j,k) = P_radyz_nuxl
	   eps_tot(i,j,k) = eps_totl
	   eps_cld(i,j,k) = eps_cldl
	   eps_thermal(i,j,k) = eps_thermall
	   P_cld(i,j,k) = P_cldl
           Sx(i,j,k) = Sxl
           Sy(i,j,k) = Syl
           Sz(i,j,k) = Szl
           Sxx(i,j,k) = Sxxl
           Sxy(i,j,k) = Sxyl
           Sxz(i,j,k) = Sxzl
           Syy(i,j,k) = Syyl
           Syz(i,j,k) = Syzl
           Szz(i,j,k) = Szzl
           T_fluid(i,j,k)=T_fluidl
        end do
     end if
     if (have_bdry_max(3)==1) then
        k=kmax
        j=2
        do i=imin,imax
           rho_bl = rho_b(i,j,k)
           Pl = P(i,j,k)
           vxl = vx(i,j,k)
           vyl = vy(i,j,k)
           vzl = vz(i,j,k)
           Bxl = Bx(i,j,k)
           Byl = By(i,j,k)
           Bzl = Bz(i,j,k)
           E_radl = E_rad(i,j,k)
           F_radxl = F_radx(i,j,k)
           F_radyl = F_rady(i,j,k)
           F_radzl = F_radz(i,j,k)
           E_rad_nuxl = E_rad_nux(i,j,k)
           F_radx_nuxl = F_radx_nux(i,j,k)
           F_rady_nuxl = F_rady_nux(i,j,k)
           F_radz_nuxl = F_radz_nux(i,j,k)
           Y_el   = Y_e(i,j,k)
           AUX(1)  = phi(i,j,k)
           AUX(2)  = lapse(i,j,k)
           AUX(3)  = shiftx(i,j,k)
           AUX(4)  = shifty(i,j,k)
           AUX(5)  = shiftz(i,j,k)
           AUX(6)  = gxx(i,j,k)
           AUX(7)  = gxy(i,j,k)
           AUX(8)  = gxz(i,j,k)
           AUX(9)  = gyy(i,j,k)
           AUX(10) = gyz(i,j,k)
           AUX(11) = gzz(i,j,k)
           AUX(12) = gupxx(i,j,k)
           AUX(13) = gupxy(i,j,k)
           AUX(14) = gupxz(i,j,k)
           AUX(15) = gupyy(i,j,k)
           AUX(16) = gupyz(i,j,k)
           AUX(17) = gupzz(i,j,k)
           call recompute_vars(rho_bl,Pl,vxl,vyl,vzl,Bxl,Byl,Bzl,AUX,m,rho_sl,taul, &
                mhd_st_xl,mhd_st_yl,mhd_st_zl,st_xl,st_yl,st_zl,wl,hl,u0l,Exl,Eyl,Ezl,&
                E_radl,F_radxl,F_radyl,F_radzl,Y_el, T_fluidl,tau_radl, S_rad_xl, S_rad_yl, S_rad_zl,rho_Yel,&
                F_rad0l, P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl, &
                E_rad_nuel,F_radx_nuel,F_rady_nuel,F_radz_nuel, tau_rad_nuel, S_rad_x_nuel, S_rad_y_nuel, S_rad_z_nuel, &
                F_rad0_nuel, P_radxx_nuel, P_radyy_nuel, P_radzz_nuel, P_radxy_nuel, P_radxz_nuel, P_radyz_nuel, &
                E_rad_nuxl,F_radx_nuxl,F_rady_nuxl,F_radz_nuxl, tau_rad_nuxl, S_rad_x_nuxl, S_rad_y_nuxl, S_rad_z_nuxl, &
                F_rad0_nuxl, P_radxx_nuxl, P_radyy_nuxl, P_radzz_nuxl, P_radxy_nuxl, P_radxz_nuxl, P_radyz_nuxl, &
		eps_totl, eps_cldl, eps_thermall, P_cldl, &
                sbtl,sbxl,sbyl,sbzl,rhol,Sxl,Syl,Szl,Sxxl,Sxyl,Sxzl,Syyl,Syzl,Szzl, &
                neos, ergo_star, ergo_sigma,rho_tab, P_tab, eps_tab, k_tab, &
     		gamma_tab, gamma_th,rho_b_atm, rad_const, &
     		compute_microphysics, enable_OS_collapse,rad_evolve_enable, rad_closure_scheme, microphysics_scheme, Erad_atm_cut)
           ! Assign primitive variables which may have been modified
           vx(i,j,k) = vxl
           vy(i,j,k) = vyl
           vz(i,j,k) = vzl
           E_rad(i,j,k) = E_radl
           E_rad_nue(i,j,k) = E_rad_nuel
           E_rad_nux(i,j,k) = E_rad_nuxl
           ! Assign conserved variables
           rho_star(i,j,k) = rho_sl
           tau(i,j,k) = taul
           mhd_st_x(i,j,k) = mhd_st_xl
           mhd_st_y(i,j,k) = mhd_st_yl
           mhd_st_z(i,j,k) = mhd_st_zl
           tau_rad(i,j,k) = tau_radl
           S_rad_x(i,j,k) = S_rad_xl
           S_rad_y(i,j,k) = S_rad_yl
           S_rad_z(i,j,k) = S_rad_zl
           tau_rad_nue(i,j,k) = tau_rad_nuel
           S_rad_x_nue(i,j,k) = S_rad_x_nuel
           S_rad_y_nue(i,j,k) = S_rad_y_nuel
           S_rad_z_nue(i,j,k) = S_rad_z_nuel
           tau_rad_nux(i,j,k) = tau_rad_nuxl
           S_rad_x_nux(i,j,k) = S_rad_x_nuxl
           S_rad_y_nux(i,j,k) = S_rad_y_nuxl
           S_rad_z_nux(i,j,k) = S_rad_z_nuxl
           rho_Ye(i,j,k) = rho_Yel
           ! Assign other variables
           st_x(i,j,k) = st_xl
           st_y(i,j,k) = st_yl
           st_z(i,j,k) = st_zl
           w(i,j,k) = wl
           h(i,j,k) = hl
           u0(i,j,k) = u0l
           Ex(i,j,k) = Exl
           Ey(i,j,k) = Eyl
           Ez(i,j,k) = Ezl
           sbt(i,j,k) = sbtl
           sbx(i,j,k) = sbxl
           sby(i,j,k) = sbyl
           sbz(i,j,k) = sbzl
           rho(i,j,k) = rhol
           F_rad0(i,j,k) = F_rad0l
           P_radxx(i,j,k) = P_radxxl
           P_radyy(i,j,k) = P_radyyl
           P_radzz(i,j,k) = P_radzzl
           P_radxy(i,j,k) = P_radxyl
           P_radxz(i,j,k) = P_radxzl
           P_radyz(i,j,k) = P_radyzl
           F_rad0_nue(i,j,k) = F_rad0_nuel
           P_radxx_nue(i,j,k) = P_radxx_nuel
           P_radyy_nue(i,j,k) = P_radyy_nuel
           P_radzz_nue(i,j,k) = P_radzz_nuel
           P_radxy_nue(i,j,k) = P_radxy_nuel
           P_radxz_nue(i,j,k) = P_radxz_nuel
           P_radyz_nue(i,j,k) = P_radyz_nuel
           F_rad0_nux(i,j,k) = F_rad0_nuxl
           P_radxx_nux(i,j,k) = P_radxx_nuxl
           P_radyy_nux(i,j,k) = P_radyy_nuxl
           P_radzz_nux(i,j,k) = P_radzz_nuxl
           P_radxy_nux(i,j,k) = P_radxy_nuxl
           P_radxz_nux(i,j,k) = P_radxz_nuxl
           P_radyz_nux(i,j,k) = P_radyz_nuxl
           eps_tot(i,j,k) = eps_totl
           eps_cld(i,j,k) = eps_cldl
           eps_thermal(i,j,k) = eps_thermall
           P_cld(i,j,k) = P_cldl
           Sx(i,j,k) = Sxl
           Sy(i,j,k) = Syl
           Sz(i,j,k) = Szl
           Sxx(i,j,k) = Sxxl
           Sxy(i,j,k) = Sxyl
           Sxz(i,j,k) = Sxzl
           Syy(i,j,k) = Syyl
           Syz(i,j,k) = Syzl
           Szz(i,j,k) = Szzl
           T_fluid(i,j,k)=T_fluidl
        end do
     end if
     if (have_bdry_min(3)==1) then
        ! k = kmin
        k=kmin
        j=2
        do i=imin,imax
           rho_bl = rho_b(i,j,k)
           Pl = P(i,j,k)
           vxl = vx(i,j,k)
           vyl = vy(i,j,k)
           vzl = vz(i,j,k)
           Bxl = Bx(i,j,k)
           Byl = By(i,j,k)
           Bzl = Bz(i,j,k)
           E_radl = E_rad(i,j,k)
           F_radxl = F_radx(i,j,k)
           F_radyl = F_rady(i,j,k)
           F_radzl = F_radz(i,j,k)
           E_rad_nuel = E_rad_nue(i,j,k)
           F_radx_nuel = F_radx_nue(i,j,k)
           F_rady_nuel = F_rady_nue(i,j,k)
           F_radz_nuel = F_radz_nue(i,j,k)
           E_rad_nuxl = E_rad_nux(i,j,k)
           F_radx_nuxl = F_radx_nux(i,j,k)
           F_rady_nuxl = F_rady_nux(i,j,k)
           F_radz_nuxl = F_radz_nux(i,j,k)
           Y_el   = Y_e(i,j,k)
           AUX(1)  = phi(i,j,k)
           AUX(2)  = lapse(i,j,k)
           AUX(3)  = shiftx(i,j,k)
           AUX(4)  = shifty(i,j,k)
           AUX(5)  = shiftz(i,j,k)
           AUX(6)  = gxx(i,j,k)
           AUX(7)  = gxy(i,j,k)
           AUX(8)  = gxz(i,j,k)
           AUX(9)  = gyy(i,j,k)
           AUX(10) = gyz(i,j,k)
           AUX(11) = gzz(i,j,k)
           AUX(12) = gupxx(i,j,k)
           AUX(13) = gupxy(i,j,k)
           AUX(14) = gupxz(i,j,k)
           AUX(15) = gupyy(i,j,k)
           AUX(16) = gupyz(i,j,k)
           AUX(17) = gupzz(i,j,k)
           call recompute_vars(rho_bl,Pl,vxl,vyl,vzl,Bxl,Byl,Bzl,AUX,m,rho_sl,taul, &
                mhd_st_xl,mhd_st_yl,mhd_st_zl,st_xl,st_yl,st_zl,wl,hl,u0l,Exl,Eyl,Ezl,&
                E_radl,F_radxl,F_radyl,F_radzl,Y_el, T_fluidl,tau_radl, S_rad_xl, S_rad_yl, S_rad_zl,rho_Yel,&
                F_rad0l, P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl, &
                E_rad_nuel,F_radx_nuel,F_rady_nuel,F_radz_nuel, tau_rad_nuel, S_rad_x_nuel, S_rad_y_nuel, S_rad_z_nuel, &
                F_rad0_nuel, P_radxx_nuel, P_radyy_nuel, P_radzz_nuel, P_radxy_nuel, P_radxz_nuel, P_radyz_nuel, &
                E_rad_nuxl,F_radx_nuxl,F_rady_nuxl,F_radz_nuxl, tau_rad_nuxl, S_rad_x_nuxl, S_rad_y_nuxl, S_rad_z_nuxl, &
                F_rad0_nuxl, P_radxx_nuxl, P_radyy_nuxl, P_radzz_nuxl, P_radxy_nuxl, P_radxz_nuxl, P_radyz_nuxl, &
		eps_totl, eps_cldl, eps_thermall, P_cldl, &                
		sbtl,sbxl,sbyl,sbzl,rhol,Sxl,Syl,Szl,Sxxl,Sxyl,Sxzl,Syyl,Syzl,Szzl, &
                neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, &
		gamma_tab, gamma_th,rho_b_atm, rad_const,&
     		compute_microphysics, enable_OS_collapse,rad_evolve_enable, rad_closure_scheme, microphysics_scheme, Erad_atm_cut)
           ! Assign primitive variables which may have been modified
           vx(i,j,k) = vxl
           vy(i,j,k) = vyl
           vz(i,j,k) = vzl
           E_rad(i,j,k) = E_radl
           E_rad_nue(i,j,k) = E_rad_nuel
           E_rad_nux(i,j,k) = E_rad_nuxl
           ! Assign conserved variables
           rho_star(i,j,k) = rho_sl
           tau(i,j,k) = taul
           mhd_st_x(i,j,k) = mhd_st_xl
           mhd_st_y(i,j,k) = mhd_st_yl
           mhd_st_z(i,j,k) = mhd_st_zl
           tau_rad(i,j,k) = tau_radl
           S_rad_x(i,j,k) = S_rad_xl
           S_rad_y(i,j,k) = S_rad_yl
           S_rad_z(i,j,k) = S_rad_zl
           tau_rad_nue(i,j,k) = tau_rad_nuel
           S_rad_x_nue(i,j,k) = S_rad_x_nuel
           S_rad_y_nue(i,j,k) = S_rad_y_nuel
           S_rad_z_nue(i,j,k) = S_rad_z_nuel
           tau_rad_nux(i,j,k) = tau_rad_nuxl
           S_rad_x_nux(i,j,k) = S_rad_x_nuxl
           S_rad_y_nux(i,j,k) = S_rad_y_nuxl
           S_rad_z_nux(i,j,k) = S_rad_z_nuxl
           rho_Ye(i,j,k) = rho_Yel
           ! Assign other variables
           st_x(i,j,k) = st_xl
           st_y(i,j,k) = st_yl
           st_z(i,j,k) = st_zl
           w(i,j,k) = wl
           h(i,j,k) = hl
           u0(i,j,k) = u0l
           Ex(i,j,k) = Exl
           Ey(i,j,k) = Eyl
           Ez(i,j,k) = Ezl
           sbt(i,j,k) = sbtl
           sbx(i,j,k) = sbxl
           sby(i,j,k) = sbyl
           sbz(i,j,k) = sbzl
           rho(i,j,k) = rhol
           F_rad0(i,j,k) = F_rad0l
           P_radxx(i,j,k) = P_radxxl
           P_radyy(i,j,k) = P_radyyl
           P_radzz(i,j,k) = P_radzzl
           P_radxy(i,j,k) = P_radxyl
           P_radxz(i,j,k) = P_radxzl
           P_radyz(i,j,k) = P_radyzl
           F_rad0_nue(i,j,k) = F_rad0_nuel
           P_radxx_nue(i,j,k) = P_radxx_nuel
           P_radyy_nue(i,j,k) = P_radyy_nuel
           P_radzz_nue(i,j,k) = P_radzz_nuel
           P_radxy_nue(i,j,k) = P_radxy_nuel
           P_radxz_nue(i,j,k) = P_radxz_nuel
           P_radyz_nue(i,j,k) = P_radyz_nuel
           F_rad0_nux(i,j,k) = F_rad0_nuxl
           P_radxx_nux(i,j,k) = P_radxx_nuxl
           P_radyy_nux(i,j,k) = P_radyy_nuxl
           P_radzz_nux(i,j,k) = P_radzz_nuxl
           P_radxy_nux(i,j,k) = P_radxy_nuxl
           P_radxz_nux(i,j,k) = P_radxz_nuxl
           P_radyz_nux(i,j,k) = P_radyz_nuxl
           eps_tot(i,j,k) = eps_totl
           eps_cld(i,j,k) = eps_cldl
           eps_thermal(i,j,k) = eps_thermall
           P_cld(i,j,k) = P_cldl
           Sx(i,j,k) = Sxl
           Sy(i,j,k) = Syl
           Sz(i,j,k) = Szl
           Sxx(i,j,k) = Sxxl
           Sxy(i,j,k) = Sxyl
           Sxz(i,j,k) = Sxzl
           Syy(i,j,k) = Syyl
           Syz(i,j,k) = Syzl
           Szz(i,j,k) = Szzl
           T_fluid(i,j,k)=T_fluidl
        end do
     end if
  else
     if (have_bdry_max(1)==1) then
        i=imax
        do k=kmin,kmax
           do j=jmin,jmax
              rho_bl = rho_b(i,j,k)
              Pl = P(i,j,k)
              vxl = vx(i,j,k)
              vyl = vy(i,j,k)
              vzl = vz(i,j,k)
              Bxl = Bx(i,j,k)
              Byl = By(i,j,k)
              Bzl = Bz(i,j,k)
              E_radl = E_rad(i,j,k)
              F_radxl = F_radx(i,j,k)
              F_radyl = F_rady(i,j,k)
              F_radzl = F_radz(i,j,k)
              E_rad_nuel = E_rad_nue(i,j,k)
              F_radx_nuel = F_radx_nue(i,j,k)
              F_rady_nuel = F_rady_nue(i,j,k)
              F_radz_nuel = F_radz_nue(i,j,k)
              E_rad_nuxl = E_rad_nux(i,j,k)
              F_radx_nuxl = F_radx_nux(i,j,k)
              F_rady_nuxl = F_rady_nux(i,j,k)
              F_radz_nuxl = F_radz_nux(i,j,k)
              Y_el   = Y_e(i,j,k)
              AUX(1)  = phi(i,j,k)
              AUX(2)  = lapse(i,j,k)
              AUX(3)  = shiftx(i,j,k)
              AUX(4)  = shifty(i,j,k)
              AUX(5)  = shiftz(i,j,k)
              AUX(6)  = gxx(i,j,k)
              AUX(7)  = gxy(i,j,k)
              AUX(8)  = gxz(i,j,k)
              AUX(9)  = gyy(i,j,k)
              AUX(10) = gyz(i,j,k)
              AUX(11) = gzz(i,j,k)
              AUX(12) = gupxx(i,j,k)
              AUX(13) = gupxy(i,j,k)
              AUX(14) = gupxz(i,j,k)
              AUX(15) = gupyy(i,j,k)
              AUX(16) = gupyz(i,j,k)
              AUX(17) = gupzz(i,j,k)
              call recompute_vars(rho_bl,Pl,vxl,vyl,vzl,Bxl,Byl,Bzl,AUX,m,rho_sl,taul, &
                   mhd_st_xl,mhd_st_yl,mhd_st_zl,st_xl,st_yl,st_zl,wl,hl,u0l,Exl,Eyl,Ezl,&
                   E_radl,F_radxl,F_radyl,F_radzl,Y_el, T_fluidl,tau_radl, S_rad_xl, S_rad_yl, S_rad_zl,rho_Yel,&
                   F_rad0l, P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl, &
                   E_rad_nuel,F_radx_nuel,F_rady_nuel,F_radz_nuel, tau_rad_nuel, S_rad_x_nuel, S_rad_y_nuel, S_rad_z_nuel, &
                   F_rad0_nuel, P_radxx_nuel, P_radyy_nuel, P_radzz_nuel, P_radxy_nuel, P_radxz_nuel, P_radyz_nuel, &
                   E_rad_nuxl,F_radx_nuxl,F_rady_nuxl,F_radz_nuxl, tau_rad_nuxl, S_rad_x_nuxl, S_rad_y_nuxl, S_rad_z_nuxl, &
                   F_rad0_nuxl, P_radxx_nuxl, P_radyy_nuxl, P_radzz_nuxl, P_radxy_nuxl, P_radxz_nuxl, P_radyz_nuxl, &
		   eps_totl, eps_cldl, eps_thermall, P_cldl, &
                   sbtl,sbxl,sbyl,sbzl,rhol,Sxl,Syl,Szl,Sxxl,Sxyl,Sxzl,Syyl,Syzl,Szzl, &
                   neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab,k_tab, & 
		   gamma_tab, gamma_th,rho_b_atm, rad_const,& 
     	    	   compute_microphysics, enable_OS_collapse,rad_evolve_enable, rad_closure_scheme, microphysics_scheme, Erad_atm_cut)
              ! Assign primitive variables which may have been modified
              vx(i,j,k) = vxl
              vy(i,j,k) = vyl
              vz(i,j,k) = vzl
              E_rad(i,j,k) = E_radl
              E_rad_nue(i,j,k) = E_rad_nuel
              E_rad_nux(i,j,k) = E_rad_nuxl
              ! Assign conserved variables
              rho_star(i,j,k) = rho_sl
              tau(i,j,k) = taul
              mhd_st_x(i,j,k) = mhd_st_xl
              mhd_st_y(i,j,k) = mhd_st_yl
              mhd_st_z(i,j,k) = mhd_st_zl
              tau_rad(i,j,k) = tau_radl
              S_rad_x(i,j,k) = S_rad_xl
              S_rad_y(i,j,k) = S_rad_yl
              S_rad_z(i,j,k) = S_rad_zl
              tau_rad_nue(i,j,k) = tau_rad_nuel
              S_rad_x_nue(i,j,k) = S_rad_x_nuel
              S_rad_y_nue(i,j,k) = S_rad_y_nuel
              S_rad_z_nue(i,j,k) = S_rad_z_nuel
              tau_rad_nux(i,j,k) = tau_rad_nuxl
              S_rad_x_nux(i,j,k) = S_rad_x_nuxl
              S_rad_y_nux(i,j,k) = S_rad_y_nuxl
              S_rad_z_nux(i,j,k) = S_rad_z_nuxl
              rho_Ye(i,j,k) = rho_Yel
              ! Assign other variables
              st_x(i,j,k) = st_xl
              st_y(i,j,k) = st_yl
              st_z(i,j,k) = st_zl
              w(i,j,k) = wl
              h(i,j,k) = hl
              u0(i,j,k) = u0l
              Ex(i,j,k) = Exl
              Ey(i,j,k) = Eyl
              Ez(i,j,k) = Ezl
              sbt(i,j,k) = sbtl
              sbx(i,j,k) = sbxl
              sby(i,j,k) = sbyl
              sbz(i,j,k) = sbzl
              rho(i,j,k) = rhol
              F_rad0(i,j,k) = F_rad0l
              P_radxx(i,j,k) = P_radxxl
              P_radyy(i,j,k) = P_radyyl
              P_radzz(i,j,k) = P_radzzl
              P_radxy(i,j,k) = P_radxyl
              P_radxz(i,j,k) = P_radxzl
              P_radyz(i,j,k) = P_radyzl
              F_rad0_nue(i,j,k) = F_rad0_nuel
              P_radxx_nue(i,j,k) = P_radxx_nuel
              P_radyy_nue(i,j,k) = P_radyy_nuel
              P_radzz_nue(i,j,k) = P_radzz_nuel
              P_radxy_nue(i,j,k) = P_radxy_nuel
              P_radxz_nue(i,j,k) = P_radxz_nuel
              P_radyz_nue(i,j,k) = P_radyz_nuel
              F_rad0_nux(i,j,k) = F_rad0_nuxl
              P_radxx_nux(i,j,k) = P_radxx_nuxl
              P_radyy_nux(i,j,k) = P_radyy_nuxl
              P_radzz_nux(i,j,k) = P_radzz_nuxl
              P_radxy_nux(i,j,k) = P_radxy_nuxl
              P_radxz_nux(i,j,k) = P_radxz_nuxl
              P_radyz_nux(i,j,k) = P_radyz_nuxl
              eps_tot(i,j,k) = eps_totl
              eps_cld(i,j,k) = eps_cldl
              eps_thermal(i,j,k) = eps_thermall
              P_cld(i,j,k) = P_cldl
              Sx(i,j,k) = Sxl
              Sy(i,j,k) = Syl
              Sz(i,j,k) = Szl
              Sxx(i,j,k) = Sxxl
              Sxy(i,j,k) = Sxyl
              Sxz(i,j,k) = Sxzl
              Syy(i,j,k) = Syyl
              Syz(i,j,k) = Syzl
              Szz(i,j,k) = Szzl
              T_fluid(i,j,k)=T_fluidl
           end do
        end do
     end if
     if (have_bdry_min(1)==1) then
        i=imin
        do k=kmin,kmax
           do j=jmin,jmax
              rho_bl = rho_b(i,j,k)
              Pl = P(i,j,k)
              vxl = vx(i,j,k)
              vyl = vy(i,j,k)
              vzl = vz(i,j,k)
              Bxl = Bx(i,j,k)
              Byl = By(i,j,k)
              Bzl = Bz(i,j,k)
              E_radl = E_rad(i,j,k)
              F_radxl = F_radx(i,j,k)
              F_radyl = F_rady(i,j,k)
              F_radzl = F_radz(i,j,k)
              E_rad_nuel = E_rad_nue(i,j,k)
              F_radx_nuel = F_radx_nue(i,j,k)
              F_rady_nuel = F_rady_nue(i,j,k)
              F_radz_nuel = F_radz_nue(i,j,k)
              E_rad_nuxl = E_rad_nux(i,j,k)
              F_radx_nuxl = F_radx_nux(i,j,k)
              F_rady_nuxl = F_rady_nux(i,j,k)
              F_radz_nuxl = F_radz_nux(i,j,k)
              Y_el   = Y_e(i,j,k)
              AUX(1)  = phi(i,j,k)
              AUX(2)  = lapse(i,j,k)
              AUX(3)  = shiftx(i,j,k)
              AUX(4)  = shifty(i,j,k)
              AUX(5)  = shiftz(i,j,k)
              AUX(6)  = gxx(i,j,k)
              AUX(7)  = gxy(i,j,k)
              AUX(8)  = gxz(i,j,k)
              AUX(9)  = gyy(i,j,k)
              AUX(10) = gyz(i,j,k)
              AUX(11) = gzz(i,j,k)
              AUX(12) = gupxx(i,j,k)
              AUX(13) = gupxy(i,j,k)
              AUX(14) = gupxz(i,j,k)
              AUX(15) = gupyy(i,j,k)
              AUX(16) = gupyz(i,j,k)
              AUX(17) = gupzz(i,j,k)
              call recompute_vars(rho_bl,Pl,vxl,vyl,vzl,Bxl,Byl,Bzl,AUX,m,rho_sl,taul, &
                   mhd_st_xl,mhd_st_yl,mhd_st_zl,st_xl,st_yl,st_zl,wl,hl,u0l,Exl,Eyl,Ezl,&
                   E_radl,F_radxl,F_radyl,F_radzl,Y_el, T_fluidl,tau_radl, S_rad_xl, S_rad_yl, S_rad_zl,rho_Yel,&
                   F_rad0l, P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl, &
                   E_rad_nuel,F_radx_nuel,F_rady_nuel,F_radz_nuel, tau_rad_nuel, S_rad_x_nuel, S_rad_y_nuel, S_rad_z_nuel, &
                   F_rad0_nuel, P_radxx_nuel, P_radyy_nuel, P_radzz_nuel, P_radxy_nuel, P_radxz_nuel, P_radyz_nuel, &
                   E_rad_nuxl,F_radx_nuxl,F_rady_nuxl,F_radz_nuxl, tau_rad_nuxl, S_rad_x_nuxl, S_rad_y_nuxl, S_rad_z_nuxl, &
                   F_rad0_nuxl, P_radxx_nuxl, P_radyy_nuxl, P_radzz_nuxl, P_radxy_nuxl, P_radxz_nuxl, P_radyz_nuxl, &
		   eps_totl, eps_cldl, eps_thermall, P_cldl, &
                   sbtl,sbxl,sbyl,sbzl,rhol,Sxl,Syl,Szl,Sxxl,Sxyl,Sxzl,Syyl,Syzl,Szzl, &
                   neos, ergo_star, ergo_sigma,rho_tab, P_tab, eps_tab, k_tab, &
     		   gamma_tab, gamma_th,rho_b_atm, rad_const,&
		   compute_microphysics, enable_OS_collapse,rad_evolve_enable, rad_closure_scheme, microphysics_scheme, Erad_atm_cut)
              ! Assign primitive variables which may have been modified
              vx(i,j,k) = vxl
              vy(i,j,k) = vyl
              vz(i,j,k) = vzl
              E_rad(i,j,k) = E_radl
              E_rad_nue(i,j,k) = E_rad_nuel
              E_rad_nux(i,j,k) = E_rad_nuxl
              ! Assign conserved variables
              rho_star(i,j,k) = rho_sl
              tau(i,j,k) = taul
              mhd_st_x(i,j,k) = mhd_st_xl
              mhd_st_y(i,j,k) = mhd_st_yl
              mhd_st_z(i,j,k) = mhd_st_zl
              tau_rad(i,j,k) = tau_radl
              S_rad_x(i,j,k) = S_rad_xl
              S_rad_y(i,j,k) = S_rad_yl
              S_rad_z(i,j,k) = S_rad_zl
              tau_rad_nue(i,j,k) = tau_rad_nuel
              S_rad_x_nue(i,j,k) = S_rad_x_nuel
              S_rad_y_nue(i,j,k) = S_rad_y_nuel
              S_rad_z_nue(i,j,k) = S_rad_z_nuel
              tau_rad_nux(i,j,k) = tau_rad_nuxl
              S_rad_x_nux(i,j,k) = S_rad_x_nuxl
              S_rad_y_nux(i,j,k) = S_rad_y_nuxl
              S_rad_z_nux(i,j,k) = S_rad_z_nuxl
              rho_Ye(i,j,k) = rho_Yel
              ! Assign other variables
              st_x(i,j,k) = st_xl
              st_y(i,j,k) = st_yl
              st_z(i,j,k) = st_zl
              w(i,j,k) = wl
              h(i,j,k) = hl
              u0(i,j,k) = u0l
              Ex(i,j,k) = Exl
              Ey(i,j,k) = Eyl
              Ez(i,j,k) = Ezl
              sbt(i,j,k) = sbtl
              sbx(i,j,k) = sbxl
              sby(i,j,k) = sbyl
              sbz(i,j,k) = sbzl
              rho(i,j,k) = rhol
              F_rad0(i,j,k) = F_rad0l
              P_radxx(i,j,k) = P_radxxl
              P_radyy(i,j,k) = P_radyyl
              P_radzz(i,j,k) = P_radzzl
              P_radxy(i,j,k) = P_radxyl
              P_radxz(i,j,k) = P_radxzl
              P_radyz(i,j,k) = P_radyzl
              F_rad0_nue(i,j,k) = F_rad0_nuel
              P_radxx_nue(i,j,k) = P_radxx_nuel
              P_radyy_nue(i,j,k) = P_radyy_nuel
              P_radzz_nue(i,j,k) = P_radzz_nuel
              P_radxy_nue(i,j,k) = P_radxy_nuel
              P_radxz_nue(i,j,k) = P_radxz_nuel
              P_radyz_nue(i,j,k) = P_radyz_nuel
              F_rad0_nux(i,j,k) = F_rad0_nuxl
              P_radxx_nux(i,j,k) = P_radxx_nuxl
              P_radyy_nux(i,j,k) = P_radyy_nuxl
              P_radzz_nux(i,j,k) = P_radzz_nuxl
              P_radxy_nux(i,j,k) = P_radxy_nuxl
              P_radxz_nux(i,j,k) = P_radxz_nuxl
              P_radyz_nux(i,j,k) = P_radyz_nuxl
              eps_tot(i,j,k) = eps_totl
              eps_cld(i,j,k) = eps_cldl
              eps_thermal(i,j,k) = eps_thermall
              P_cld(i,j,k) = P_cldl
              Sx(i,j,k) = Sxl
              Sy(i,j,k) = Syl
              Sz(i,j,k) = Szl
              Sxx(i,j,k) = Sxxl
              Sxy(i,j,k) = Sxyl
              Sxz(i,j,k) = Sxzl
              Syy(i,j,k) = Syyl
              Syz(i,j,k) = Syzl
              Szz(i,j,k) = Szzl
              T_fluid(i,j,k)=T_fluidl
           end do
        end do
     end if
     if (have_bdry_max(2)==1) then
        j=jmax
        do k=kmin,kmax
           do i=imin,imax
              rho_bl = rho_b(i,j,k)
              Pl = P(i,j,k)
              vxl = vx(i,j,k)
              vyl = vy(i,j,k)
              vzl = vz(i,j,k)
              Bxl = Bx(i,j,k)
              Byl = By(i,j,k)
              Bzl = Bz(i,j,k)
              E_radl = E_rad(i,j,k)
              F_radxl = F_radx(i,j,k)
              F_radyl = F_rady(i,j,k)
              F_radzl = F_radz(i,j,k)
              E_rad_nuel = E_rad_nue(i,j,k)
              F_radx_nuel = F_radx_nue(i,j,k)
              F_rady_nuel = F_rady_nue(i,j,k)
              F_radz_nuel = F_radz_nue(i,j,k)
              E_rad_nuxl = E_rad_nux(i,j,k)
              F_radx_nuxl = F_radx_nux(i,j,k)
              F_rady_nuxl = F_rady_nux(i,j,k)
              F_radz_nuxl = F_radz_nux(i,j,k)
              Y_el   = Y_e(i,j,k)
              AUX(1)  = phi(i,j,k)
              AUX(2)  = lapse(i,j,k)
              AUX(3)  = shiftx(i,j,k)
              AUX(4)  = shifty(i,j,k)
              AUX(5)  = shiftz(i,j,k)
              AUX(6)  = gxx(i,j,k)
              AUX(7)  = gxy(i,j,k)
              AUX(8)  = gxz(i,j,k)
              AUX(9)  = gyy(i,j,k)
              AUX(10) = gyz(i,j,k)
              AUX(11) = gzz(i,j,k)
              AUX(12) = gupxx(i,j,k)
              AUX(13) = gupxy(i,j,k)
              AUX(14) = gupxz(i,j,k)
              AUX(15) = gupyy(i,j,k)
              AUX(16) = gupyz(i,j,k)
              AUX(17) = gupzz(i,j,k)
              call recompute_vars(rho_bl,Pl,vxl,vyl,vzl,Bxl,Byl,Bzl,AUX,m,rho_sl,taul, &
                   mhd_st_xl,mhd_st_yl,mhd_st_zl,st_xl,st_yl,st_zl,wl,hl,u0l,Exl,Eyl,Ezl,&
                   E_radl,F_radxl,F_radyl,F_radzl,Y_el, T_fluidl,tau_radl, S_rad_xl, S_rad_yl, S_rad_zl,rho_Yel,&
                   F_rad0l, P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl, &
                   E_rad_nuel,F_radx_nuel,F_rady_nuel,F_radz_nuel, tau_rad_nuel, S_rad_x_nuel, S_rad_y_nuel, S_rad_z_nuel, &
                   F_rad0_nuel, P_radxx_nuel, P_radyy_nuel, P_radzz_nuel, P_radxy_nuel, P_radxz_nuel, P_radyz_nuel, &
                   E_rad_nuxl,F_radx_nuxl,F_rady_nuxl,F_radz_nuxl, tau_rad_nuxl, S_rad_x_nuxl, S_rad_y_nuxl, S_rad_z_nuxl, &
                   F_rad0_nuxl, P_radxx_nuxl, P_radyy_nuxl, P_radzz_nuxl, P_radxy_nuxl, P_radxz_nuxl, P_radyz_nuxl, &
		   eps_totl, eps_cldl, eps_thermall, P_cldl, &
                   sbtl,sbxl,sbyl,sbzl,rhol,Sxl,Syl,Szl,Sxxl,Sxyl,Sxzl,Syyl,Syzl,Szzl, &
                   neos,ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, &
     		   gamma_tab, gamma_th,rho_b_atm, rad_const,&
     		   compute_microphysics, enable_OS_collapse,rad_evolve_enable, rad_closure_scheme, microphysics_scheme, Erad_atm_cut)
              ! Assign primitive variables which may have been modified
              vx(i,j,k) = vxl
              vy(i,j,k) = vyl
              vz(i,j,k) = vzl
              E_rad(i,j,k) = E_radl
              E_rad_nue(i,j,k) = E_rad_nuel
              E_rad_nux(i,j,k) = E_rad_nuxl
              ! Assign conserved variables
              rho_star(i,j,k) = rho_sl
              tau(i,j,k) = taul
              mhd_st_x(i,j,k) = mhd_st_xl
              mhd_st_y(i,j,k) = mhd_st_yl
              mhd_st_z(i,j,k) = mhd_st_zl
              tau_rad(i,j,k) = tau_radl
              S_rad_x(i,j,k) = S_rad_xl
              S_rad_y(i,j,k) = S_rad_yl
              S_rad_z(i,j,k) = S_rad_zl
              tau_rad_nue(i,j,k) = tau_rad_nuel
              S_rad_x_nue(i,j,k) = S_rad_x_nuel
              S_rad_y_nue(i,j,k) = S_rad_y_nuel
              S_rad_z_nue(i,j,k) = S_rad_z_nuel
              tau_rad_nux(i,j,k) = tau_rad_nuxl
              S_rad_x_nux(i,j,k) = S_rad_x_nuxl
              S_rad_y_nux(i,j,k) = S_rad_y_nuxl
              S_rad_z_nux(i,j,k) = S_rad_z_nuxl
              rho_Ye(i,j,k) = rho_Yel
              ! Assign other variables
              st_x(i,j,k) = st_xl
              st_y(i,j,k) = st_yl
              st_z(i,j,k) = st_zl
              w(i,j,k) = wl
              h(i,j,k) = hl
              u0(i,j,k) = u0l
              Ex(i,j,k) = Exl
              Ey(i,j,k) = Eyl
              Ez(i,j,k) = Ezl
              sbt(i,j,k) = sbtl
              sbx(i,j,k) = sbxl
              sby(i,j,k) = sbyl
              sbz(i,j,k) = sbzl
              rho(i,j,k) = rhol
              F_rad0(i,j,k) = F_rad0l
              P_radxx(i,j,k) = P_radxxl
              P_radyy(i,j,k) = P_radyyl
              P_radzz(i,j,k) = P_radzzl
              P_radxy(i,j,k) = P_radxyl
              P_radxz(i,j,k) = P_radxzl
              P_radyz(i,j,k) = P_radyzl
              F_rad0_nue(i,j,k) = F_rad0_nuel
              P_radxx_nue(i,j,k) = P_radxx_nuel
              P_radyy_nue(i,j,k) = P_radyy_nuel
              P_radzz_nue(i,j,k) = P_radzz_nuel
              P_radxy_nue(i,j,k) = P_radxy_nuel
              P_radxz_nue(i,j,k) = P_radxz_nuel
              P_radyz_nue(i,j,k) = P_radyz_nuel
              F_rad0_nux(i,j,k) = F_rad0_nuxl
              P_radxx_nux(i,j,k) = P_radxx_nuxl
              P_radyy_nux(i,j,k) = P_radyy_nuxl
              P_radzz_nux(i,j,k) = P_radzz_nuxl
              P_radxy_nux(i,j,k) = P_radxy_nuxl
              P_radxz_nux(i,j,k) = P_radxz_nuxl
              P_radyz_nux(i,j,k) = P_radyz_nuxl
              eps_tot(i,j,k) = eps_totl
              eps_cld(i,j,k) = eps_cldl
              eps_thermal(i,j,k) = eps_thermall
               P_cld(i,j,k) = P_cldl
              Sx(i,j,k) = Sxl
              Sy(i,j,k) = Syl
              Sz(i,j,k) = Szl
              Sxx(i,j,k) = Sxxl
              Sxy(i,j,k) = Sxyl
              Sxz(i,j,k) = Sxzl
              Syy(i,j,k) = Syyl
              Syz(i,j,k) = Syzl
              Szz(i,j,k) = Szzl
              T_fluid(i,j,k)=T_fluidl
           end do
        end do
     end if
     if (have_bdry_min(2)==1) then
        j=jmin
        do k=kmin,kmax
           do i=imin,imax
              rho_bl = rho_b(i,j,k)
              Pl = P(i,j,k)
              vxl = vx(i,j,k)
              vyl = vy(i,j,k)
              vzl = vz(i,j,k)
              Bxl = Bx(i,j,k)
              Byl = By(i,j,k)
              Bzl = Bz(i,j,k)
              E_radl = E_rad(i,j,k)
              F_radxl = F_radx(i,j,k)
              F_radyl = F_rady(i,j,k)
              F_radzl = F_radz(i,j,k)
              Y_el   = Y_e(i,j,k)
              AUX(1)  = phi(i,j,k)
              AUX(2)  = lapse(i,j,k)
              AUX(3)  = shiftx(i,j,k)
              AUX(4)  = shifty(i,j,k)
              AUX(5)  = shiftz(i,j,k)
              AUX(6)  = gxx(i,j,k)
              AUX(7)  = gxy(i,j,k)
              AUX(8)  = gxz(i,j,k)
              AUX(9)  = gyy(i,j,k)
              AUX(10) = gyz(i,j,k)
              AUX(11) = gzz(i,j,k)
              AUX(12) = gupxx(i,j,k)
              AUX(13) = gupxy(i,j,k)
              AUX(14) = gupxz(i,j,k)
              AUX(15) = gupyy(i,j,k)
              AUX(16) = gupyz(i,j,k)
              AUX(17) = gupzz(i,j,k)
              call recompute_vars(rho_bl,Pl,vxl,vyl,vzl,Bxl,Byl,Bzl,AUX,m,rho_sl,taul, &
                   mhd_st_xl,mhd_st_yl,mhd_st_zl,st_xl,st_yl,st_zl,wl,hl,u0l,Exl,Eyl,Ezl,&
                   E_radl,F_radxl,F_radyl,F_radzl,Y_el, T_fluidl,tau_radl, S_rad_xl, S_rad_yl, S_rad_zl,rho_Yel,&
                   F_rad0l, P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl, &
                   E_rad_nuel,F_radx_nuel,F_rady_nuel,F_radz_nuel, tau_rad_nuel, S_rad_x_nuel, S_rad_y_nuel, S_rad_z_nuel, &
                   F_rad0_nuel, P_radxx_nuel, P_radyy_nuel, P_radzz_nuel, P_radxy_nuel, P_radxz_nuel, P_radyz_nuel, &
                   E_rad_nuxl,F_radx_nuxl,F_rady_nuxl,F_radz_nuxl, tau_rad_nuxl, S_rad_x_nuxl, S_rad_y_nuxl, S_rad_z_nuxl, &
                   F_rad0_nuxl, P_radxx_nuxl, P_radyy_nuxl, P_radzz_nuxl, P_radxy_nuxl, P_radxz_nuxl, P_radyz_nuxl, &
		   eps_totl, eps_cldl, eps_thermall, P_cldl, &
                   sbtl,sbxl,sbyl,sbzl,rhol,Sxl,Syl,Szl,Sxxl,Sxyl,Sxzl,Syyl,Syzl,Szzl, &
                   neos, ergo_star, ergo_sigma,rho_tab, P_tab, eps_tab, k_tab, &
     		   gamma_tab, gamma_th,rho_b_atm, rad_const,&
     		   compute_microphysics, enable_OS_collapse,rad_evolve_enable, rad_closure_scheme, microphysics_scheme, Erad_atm_cut)
              ! Assign primitive variables which may have been modified
              vx(i,j,k) = vxl
              vy(i,j,k) = vyl
              vz(i,j,k) = vzl
              E_rad(i,j,k) = E_radl
              E_rad_nue(i,j,k) = E_rad_nuel
              E_rad_nux(i,j,k) = E_rad_nuxl
              ! Assign conserved variables
              rho_star(i,j,k) = rho_sl
              tau(i,j,k) = taul
              mhd_st_x(i,j,k) = mhd_st_xl
              mhd_st_y(i,j,k) = mhd_st_yl
              mhd_st_z(i,j,k) = mhd_st_zl
              tau_rad(i,j,k) = tau_radl
              S_rad_x(i,j,k) = S_rad_xl
              S_rad_y(i,j,k) = S_rad_yl
              S_rad_z(i,j,k) = S_rad_zl
              tau_rad_nue(i,j,k) = tau_rad_nuel
              S_rad_x_nue(i,j,k) = S_rad_x_nuel
              S_rad_y_nue(i,j,k) = S_rad_y_nuel
              S_rad_z_nue(i,j,k) = S_rad_z_nuel
              tau_rad_nux(i,j,k) = tau_rad_nuxl
              S_rad_x_nux(i,j,k) = S_rad_x_nuxl
              S_rad_y_nux(i,j,k) = S_rad_y_nuxl
              S_rad_z_nux(i,j,k) = S_rad_z_nuxl
              rho_Ye(i,j,k) = rho_Yel
              ! Assign other variables
              st_x(i,j,k) = st_xl
              st_y(i,j,k) = st_yl
              st_z(i,j,k) = st_zl
              w(i,j,k) = wl
              h(i,j,k) = hl
              u0(i,j,k) = u0l
              Ex(i,j,k) = Exl
              Ey(i,j,k) = Eyl
              Ez(i,j,k) = Ezl
              sbt(i,j,k) = sbtl
              sbx(i,j,k) = sbxl
              sby(i,j,k) = sbyl
              sbz(i,j,k) = sbzl
              rho(i,j,k) = rhol
              F_rad0(i,j,k) = F_rad0l
              P_radxx(i,j,k) = P_radxxl
              P_radyy(i,j,k) = P_radyyl
              P_radzz(i,j,k) = P_radzzl
              P_radxy(i,j,k) = P_radxyl
              P_radxz(i,j,k) = P_radxzl
              P_radyz(i,j,k) = P_radyzl
              F_rad0_nue(i,j,k) = F_rad0_nuel
              P_radxx_nue(i,j,k) = P_radxx_nuel
              P_radyy_nue(i,j,k) = P_radyy_nuel
              P_radzz_nue(i,j,k) = P_radzz_nuel
              P_radxy_nue(i,j,k) = P_radxy_nuel
              P_radxz_nue(i,j,k) = P_radxz_nuel
              P_radyz_nue(i,j,k) = P_radyz_nuel
              F_rad0_nux(i,j,k) = F_rad0_nuxl
              P_radxx_nux(i,j,k) = P_radxx_nuxl
              P_radyy_nux(i,j,k) = P_radyy_nuxl
              P_radzz_nux(i,j,k) = P_radzz_nuxl
              P_radxy_nux(i,j,k) = P_radxy_nuxl
              P_radxz_nux(i,j,k) = P_radxz_nuxl
              P_radyz_nux(i,j,k) = P_radyz_nuxl
              eps_tot(i,j,k) = eps_totl
              eps_cld(i,j,k) = eps_cldl
              eps_thermal(i,j,k) = eps_thermall
              P_cld(i,j,k) = P_cldl
              Sx(i,j,k) = Sxl
              Sy(i,j,k) = Syl
              Sz(i,j,k) = Szl
              Sxx(i,j,k) = Sxxl
              Sxy(i,j,k) = Sxyl
              Sxz(i,j,k) = Sxzl
              Syy(i,j,k) = Syyl
              Syz(i,j,k) = Syzl
              Szz(i,j,k) = Szzl
              T_fluid(i,j,k)=T_fluidl
           end do
        end do
     end if
     if (have_bdry_max(3)==1) then
        k=kmax
        do j=jmin,jmax
           do i=imin,imax
              rho_bl = rho_b(i,j,k)
              Pl = P(i,j,k)
              vxl = vx(i,j,k)
              vyl = vy(i,j,k)
              vzl = vz(i,j,k)
              Bxl = Bx(i,j,k)
              Byl = By(i,j,k)
              Bzl = Bz(i,j,k)
              E_radl = E_rad(i,j,k)
              F_radxl = F_radx(i,j,k)
              F_radyl = F_rady(i,j,k)
              F_radzl = F_radz(i,j,k)
              E_rad_nuel = E_rad_nue(i,j,k)
              F_radx_nuel = F_radx_nue(i,j,k)
              F_rady_nuel = F_rady_nue(i,j,k)
              F_radz_nuel = F_radz_nue(i,j,k)
              E_rad_nuxl = E_rad_nux(i,j,k)
              F_radx_nuxl = F_radx_nux(i,j,k)
              F_rady_nuxl = F_rady_nux(i,j,k)
              F_radz_nuxl = F_radz_nux(i,j,k)
              Y_el   = Y_e(i,j,k)
              AUX(1)  = phi(i,j,k)
              AUX(2)  = lapse(i,j,k)
              AUX(3)  = shiftx(i,j,k)
              AUX(4)  = shifty(i,j,k)
              AUX(5)  = shiftz(i,j,k)
              AUX(6)  = gxx(i,j,k)
              AUX(7)  = gxy(i,j,k)
              AUX(8)  = gxz(i,j,k)
              AUX(9)  = gyy(i,j,k)
              AUX(10) = gyz(i,j,k)
              AUX(11) = gzz(i,j,k)
              AUX(12) = gupxx(i,j,k)
              AUX(13) = gupxy(i,j,k)
              AUX(14) = gupxz(i,j,k)
              AUX(15) = gupyy(i,j,k)
              AUX(16) = gupyz(i,j,k)
              AUX(17) = gupzz(i,j,k)
              !if(i.eq.24 .and. j.eq.24) then
              !   write(*,*) midway2 bc rho_star = ,AUX,rho_bl,Pl,vxl,vyl,vzl,Bxl,Byl,Bzl
              !end if
              call recompute_vars(rho_bl,Pl,vxl,vyl,vzl,Bxl,Byl,Bzl,AUX,m,rho_sl,taul, &
                   mhd_st_xl,mhd_st_yl,mhd_st_zl,st_xl,st_yl,st_zl,wl,hl,u0l,Exl,Eyl,Ezl,&
                   E_radl,F_radxl,F_radyl,F_radzl,Y_el, T_fluidl,tau_radl, S_rad_xl, S_rad_yl, S_rad_zl,rho_Yel,&
                   F_rad0l, P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl, &
                   E_rad_nuel,F_radx_nuel,F_rady_nuel,F_radz_nuel, tau_rad_nuel, S_rad_x_nuel, S_rad_y_nuel, S_rad_z_nuel, &
                   F_rad0_nuel, P_radxx_nuel, P_radyy_nuel, P_radzz_nuel, P_radxy_nuel, P_radxz_nuel, P_radyz_nuel, &
                   E_rad_nuxl,F_radx_nuxl,F_rady_nuxl,F_radz_nuxl, tau_rad_nuxl, S_rad_x_nuxl, S_rad_y_nuxl, S_rad_z_nuxl, &
                   F_rad0_nuxl, P_radxx_nuxl, P_radyy_nuxl, P_radzz_nuxl, P_radxy_nuxl, P_radxz_nuxl, P_radyz_nuxl, &
		   eps_totl, eps_cldl, eps_thermall, P_cldl, &
                   sbtl,sbxl,sbyl,sbzl,rhol,Sxl,Syl,Szl,Sxxl,Sxyl,Sxzl,Syyl,Syzl,Szzl, &
                   neos, ergo_star, ergo_sigma,rho_tab, P_tab, eps_tab, k_tab, &
     		   gamma_tab, gamma_th,rho_b_atm, rad_const,& 
     		   compute_microphysics, enable_OS_collapse,rad_evolve_enable, rad_closure_scheme, microphysics_scheme, Erad_atm_cut)
              ! Assign primitive variables which may have been modified
              vx(i,j,k) = vxl
              vy(i,j,k) = vyl
              vz(i,j,k) = vzl
              E_rad(i,j,k) = E_radl
              E_rad_nue(i,j,k) = E_rad_nuel
              E_rad_nux(i,j,k) = E_rad_nuxl
              ! Assign conserved variables
              rho_star(i,j,k) = rho_sl
              tau(i,j,k) = taul
              mhd_st_x(i,j,k) = mhd_st_xl
              mhd_st_y(i,j,k) = mhd_st_yl
              mhd_st_z(i,j,k) = mhd_st_zl
              tau_rad(i,j,k) = tau_radl
              S_rad_x(i,j,k) = S_rad_xl
              S_rad_y(i,j,k) = S_rad_yl
              S_rad_z(i,j,k) = S_rad_zl
              tau_rad_nue(i,j,k) = tau_rad_nuel
              S_rad_x_nue(i,j,k) = S_rad_x_nuel
              S_rad_y_nue(i,j,k) = S_rad_y_nuel
              S_rad_z_nue(i,j,k) = S_rad_z_nuel
              tau_rad_nux(i,j,k) = tau_rad_nuxl
              S_rad_x_nux(i,j,k) = S_rad_x_nuxl
              S_rad_y_nux(i,j,k) = S_rad_y_nuxl
              S_rad_z_nux(i,j,k) = S_rad_z_nuxl
              rho_Ye(i,j,k) = rho_Yel
              ! Assign other variables
              st_x(i,j,k) = st_xl
              st_y(i,j,k) = st_yl
              st_z(i,j,k) = st_zl
              w(i,j,k) = wl
              h(i,j,k) = hl
              u0(i,j,k) = u0l
              Ex(i,j,k) = Exl
              Ey(i,j,k) = Eyl
              Ez(i,j,k) = Ezl
              sbt(i,j,k) = sbtl
              sbx(i,j,k) = sbxl
              sby(i,j,k) = sbyl
              sbz(i,j,k) = sbzl
              rho(i,j,k) = rhol
              F_rad0(i,j,k) = F_rad0l
              P_radxx(i,j,k) = P_radxxl
              P_radyy(i,j,k) = P_radyyl
              P_radzz(i,j,k) = P_radzzl
              P_radxy(i,j,k) = P_radxyl
              P_radxz(i,j,k) = P_radxzl
              P_radyz(i,j,k) = P_radyzl
              F_rad0_nue(i,j,k) = F_rad0_nuel
              P_radxx_nue(i,j,k) = P_radxx_nuel
              P_radyy_nue(i,j,k) = P_radyy_nuel
              P_radzz_nue(i,j,k) = P_radzz_nuel
              P_radxy_nue(i,j,k) = P_radxy_nuel
              P_radxz_nue(i,j,k) = P_radxz_nuel
              P_radyz_nue(i,j,k) = P_radyz_nuel
              F_rad0_nux(i,j,k) = F_rad0_nuxl
              P_radxx_nux(i,j,k) = P_radxx_nuxl
              P_radyy_nux(i,j,k) = P_radyy_nuxl
              P_radzz_nux(i,j,k) = P_radzz_nuxl
              P_radxy_nux(i,j,k) = P_radxy_nuxl
              P_radxz_nux(i,j,k) = P_radxz_nuxl
              P_radyz_nux(i,j,k) = P_radyz_nuxl
              eps_tot(i,j,k) = eps_totl
              eps_cld(i,j,k) = eps_cldl
              eps_thermal(i,j,k) = eps_thermall
              P_cld(i,j,k) = P_cldl
              Sx(i,j,k) = Sxl
              Sy(i,j,k) = Syl
              Sz(i,j,k) = Szl
              Sxx(i,j,k) = Sxxl
              Sxy(i,j,k) = Sxyl
              Sxz(i,j,k) = Sxzl
              Syy(i,j,k) = Syyl
              Syz(i,j,k) = Syzl
              Szz(i,j,k) = Szzl
              T_fluid(i,j,k)=T_fluidl
           end do
        end do
     end if
     if (have_bdry_min(3)==1) then
        k=kmin
        do j=jmin,jmax
           do i=imin,imax
              rho_bl = rho_b(i,j,k)
              Pl = P(i,j,k)
              vxl = vx(i,j,k)
              vyl = vy(i,j,k)
              vzl = vz(i,j,k)
              Bxl = Bx(i,j,k)
              Byl = By(i,j,k)
              Bzl = Bz(i,j,k)
              E_radl = E_rad(i,j,k)
              F_radxl = F_radx(i,j,k)
              F_radyl = F_rady(i,j,k)
              F_radzl = F_radz(i,j,k)
              E_rad_nuel = E_rad_nue(i,j,k)
              F_radx_nuel = F_radx_nue(i,j,k)
              F_rady_nuel = F_rady_nue(i,j,k)
              F_radz_nuel = F_radz_nue(i,j,k)
              E_rad_nuxl = E_rad_nux(i,j,k)
              F_radx_nuxl = F_radx_nux(i,j,k)
              F_rady_nuxl = F_rady_nux(i,j,k)
              F_radz_nuxl = F_radz_nux(i,j,k)
              Y_el   = Y_e(i,j,k)
              AUX(1)  = phi(i,j,k)
              AUX(2)  = lapse(i,j,k)
              AUX(3)  = shiftx(i,j,k)
              AUX(4)  = shifty(i,j,k)
              AUX(5)  = shiftz(i,j,k)
              AUX(6)  = gxx(i,j,k)
              AUX(7)  = gxy(i,j,k)
              AUX(8)  = gxz(i,j,k)
              AUX(9)  = gyy(i,j,k)
              AUX(10) = gyz(i,j,k)
              AUX(11) = gzz(i,j,k)
              AUX(12) = gupxx(i,j,k)
              AUX(13) = gupxy(i,j,k)
              AUX(14) = gupxz(i,j,k)
              AUX(15) = gupyy(i,j,k)
              AUX(16) = gupyz(i,j,k)
              AUX(17) = gupzz(i,j,k)
              call recompute_vars(rho_bl,Pl,vxl,vyl,vzl,Bxl,Byl,Bzl,AUX,m,rho_sl,taul, &
                   mhd_st_xl,mhd_st_yl,mhd_st_zl,st_xl,st_yl,st_zl,wl,hl,u0l,Exl,Eyl,Ezl,&
                   E_radl,F_radxl,F_radyl,F_radzl,Y_el, T_fluidl,tau_radl,S_rad_xl, S_rad_yl, S_rad_zl,rho_Yel,&
                   F_rad0l, P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl, &
                   E_rad_nuel,F_radx_nuel,F_rady_nuel,F_radz_nuel, tau_rad_nuel, S_rad_x_nuel, S_rad_y_nuel, S_rad_z_nuel, &
                   F_rad0_nuel, P_radxx_nuel, P_radyy_nuel, P_radzz_nuel, P_radxy_nuel, P_radxz_nuel, P_radyz_nuel, &
                   E_rad_nuxl,F_radx_nuxl,F_rady_nuxl,F_radz_nuxl, tau_rad_nuxl, S_rad_x_nuxl, S_rad_y_nuxl, S_rad_z_nuxl, &
                   F_rad0_nuxl, P_radxx_nuxl, P_radyy_nuxl, P_radzz_nuxl, P_radxy_nuxl, P_radxz_nuxl, P_radyz_nuxl, &
		   eps_totl, eps_cldl, eps_thermall, P_cldl, &
                   sbtl,sbxl,sbyl,sbzl,rhol,Sxl,Syl,Szl,Sxxl,Sxyl,Sxzl,Syyl,Syzl,Szzl, &
                   neos, ergo_star, ergo_sigma,rho_tab, P_tab, eps_tab, k_tab, &
     		   gamma_tab, gamma_th,rho_b_atm, rad_const,& 
     		   compute_microphysics, enable_OS_collapse,rad_evolve_enable, rad_closure_scheme, microphysics_scheme, Erad_atm_cut)
              ! Assign primitive variables which may have been modified
              vx(i,j,k) = vxl
              vy(i,j,k) = vyl
              vz(i,j,k) = vzl
              E_rad(i,j,k) = E_radl
              E_rad_nue(i,j,k) = E_rad_nuel
              E_rad_nux(i,j,k) = E_rad_nuxl
              ! Assign conserved variables
              rho_star(i,j,k) = rho_sl
              tau(i,j,k) = taul
              mhd_st_x(i,j,k) = mhd_st_xl
              mhd_st_y(i,j,k) = mhd_st_yl
              mhd_st_z(i,j,k) = mhd_st_zl
              tau_rad(i,j,k) = tau_radl
              S_rad_x(i,j,k) = S_rad_xl
              S_rad_y(i,j,k) = S_rad_yl
              S_rad_z(i,j,k) = S_rad_zl
              tau_rad_nue(i,j,k) = tau_rad_nuel
              S_rad_x_nue(i,j,k) = S_rad_x_nuel
              S_rad_y_nue(i,j,k) = S_rad_y_nuel
              S_rad_z_nue(i,j,k) = S_rad_z_nuel
              tau_rad_nux(i,j,k) = tau_rad_nuxl
              S_rad_x_nux(i,j,k) = S_rad_x_nuxl
              S_rad_y_nux(i,j,k) = S_rad_y_nuxl
              S_rad_z_nux(i,j,k) = S_rad_z_nuxl
              rho_Ye(i,j,k) = rho_Yel
              ! Assign other variables
              st_x(i,j,k) = st_xl
              st_y(i,j,k) = st_yl
              st_z(i,j,k) = st_zl
              w(i,j,k) = wl
              h(i,j,k) = hl
              u0(i,j,k) = u0l
              Ex(i,j,k) = Exl
              Ey(i,j,k) = Eyl
              Ez(i,j,k) = Ezl
              sbt(i,j,k) = sbtl
              sbx(i,j,k) = sbxl
              sby(i,j,k) = sbyl
              sbz(i,j,k) = sbzl
              rho(i,j,k) = rhol
              F_rad0(i,j,k) = F_rad0l
              P_radxx(i,j,k) = P_radxxl
              P_radyy(i,j,k) = P_radyyl
              P_radzz(i,j,k) = P_radzzl
              P_radxy(i,j,k) = P_radxyl
              P_radxz(i,j,k) = P_radxzl
              P_radyz(i,j,k) = P_radyzl
              F_rad0_nue(i,j,k) = F_rad0_nuel
              P_radxx_nue(i,j,k) = P_radxx_nuel
              P_radyy_nue(i,j,k) = P_radyy_nuel
              P_radzz_nue(i,j,k) = P_radzz_nuel
              P_radxy_nue(i,j,k) = P_radxy_nuel
              P_radxz_nue(i,j,k) = P_radxz_nuel
              P_radyz_nue(i,j,k) = P_radyz_nuel
              F_rad0_nux(i,j,k) = F_rad0_nuxl
              P_radxx_nux(i,j,k) = P_radxx_nuxl
              P_radyy_nux(i,j,k) = P_radyy_nuxl
              P_radzz_nux(i,j,k) = P_radzz_nuxl
              P_radxy_nux(i,j,k) = P_radxy_nuxl
              P_radxz_nux(i,j,k) = P_radxz_nuxl
              P_radyz_nux(i,j,k) = P_radyz_nuxl
              eps_tot(i,j,k) = eps_totl
              eps_cld(i,j,k) = eps_cldl
              eps_thermal(i,j,k) = eps_thermall
              P_cld(i,j,k) = P_cldl
              Sx(i,j,k) = Sxl
              Sy(i,j,k) = Syl
              Sz(i,j,k) = Szl
              Sxx(i,j,k) = Sxxl
              Sxy(i,j,k) = Sxyl
              Sxz(i,j,k) = Sxzl
              Syy(i,j,k) = Syyl
              Syz(i,j,k) = Syzl
              Szz(i,j,k) = Szzl
              T_fluid(i,j,k)=T_fluidl
           end do
        end do
     end if
  end if

end subroutine apply_matter_bc_mhd_new





subroutine apply_matter_frozen_bc(ext, fake_ext, &
                rho_star,tau,mhd_st_x,mhd_st_y,mhd_st_z, &
                Bxtilde,Bytilde,Bztilde, &
                rho_star_old,tau_old,mhd_st_x_old,mhd_st_y_old,mhd_st_z_old, &
                Bxtilde_old,Bytilde_old,Bztilde_old, &
                have_bdry_min,have_bdry_max)
  implicit none
  integer, dimension(3)                   :: ext,fake_ext,have_bdry_min,have_bdry_max
  real*8, dimension(ext(1),ext(2),ext(3)) :: rho_star, tau, mhd_st_x, mhd_st_y, mhd_st_z
  real*8, dimension(ext(1),ext(2),ext(3)) :: Bxtilde,Bytilde,Bztilde
  real*8, dimension(ext(1),ext(2),ext(3)) :: rho_star_old, tau_old, mhd_st_x_old, mhd_st_y_old, mhd_st_z_old
  real*8, dimension(ext(1),ext(2),ext(3)) :: Bxtilde_old,Bytilde_old,Bztilde_old
  integer :: imin,jmin,kmin,imax,jmax,kmax, i,j,k
!
  imin = ext(1)-fake_ext(1)+1
  jmin = ext(2)-fake_ext(2)+1
  kmin = ext(3)-fake_ext(3)+1
  imax = fake_ext(1)
  jmax = fake_ext(2)
  kmax = fake_ext(3)
  if (have_bdry_max(1)==1) then
     !$omp parallel do
     do k=1,ext(3)
        do j=1,ext(2)
           rho_star(imax,j,k) = rho_star_old(imax,j,k)
           tau(imax,j,k) = tau_old(imax,j,k)
           mhd_st_x(imax,j,k) = mhd_st_x_old(imax,j,k)
           mhd_st_y(imax,j,k) = mhd_st_y_old(imax,j,k)
           mhd_st_z(imax,j,k) = mhd_st_z_old(imax,j,k)
           Bxtilde(imax,j,k) = Bxtilde_old(imax,j,k)
           Bytilde(imax,j,k) = Bytilde_old(imax,j,k)
           Bztilde(imax,j,k) = Bztilde_old(imax,j,k)
        end do
     end do
     !$omp end parallel do
  end if
  if (have_bdry_min(1)==1) then
     !$omp parallel do
     do k=1,ext(3)
        do j=1,ext(2)
           rho_star(imin,j,k) = rho_star_old(imin,j,k)
           tau(imin,j,k) = tau_old(imin,j,k)
           mhd_st_x(imin,j,k) = mhd_st_x_old(imin,j,k)
           mhd_st_y(imin,j,k) = mhd_st_y_old(imin,j,k)
           mhd_st_z(imin,j,k) = mhd_st_z_old(imin,j,k)
           Bxtilde(imin,j,k) = Bxtilde_old(imin,j,k)
           Bytilde(imin,j,k) = Bytilde_old(imin,j,k)
           Bztilde(imin,j,k) = Bztilde_old(imin,j,k)
        end do
     end do
     !$omp end parallel do
  end if
  if (have_bdry_max(2)==1) then
     !$omp parallel do
     do k=1,ext(3)
        do i=1,ext(1)
           rho_star(i,jmax,k) = rho_star_old(i,jmax,k)
           tau(i,jmax,k) = tau_old(i,jmax,k)
           mhd_st_x(i,jmax,k) = mhd_st_x_old(i,jmax,k)
           mhd_st_y(i,jmax,k) = mhd_st_y_old(i,jmax,k)
           mhd_st_z(i,jmax,k) = mhd_st_z_old(i,jmax,k)
           Bxtilde(i,jmax,k) = Bxtilde_old(i,jmax,k)
           Bytilde(i,jmax,k) = Bytilde_old(i,jmax,k)
           Bztilde(i,jmax,k) = Bztilde_old(i,jmax,k)
        end do
     end do
     !$omp end parallel do
  end if
  if (have_bdry_min(2)==1) then
     !$omp parallel do
     do k=1,ext(3)
        do i=1,ext(1)
           rho_star(i,jmin,k) = rho_star_old(i,jmin,k)
           tau(i,jmin,k) = tau_old(i,jmin,k)
           mhd_st_x(i,jmin,k) = mhd_st_x_old(i,jmin,k)
           mhd_st_y(i,jmin,k) = mhd_st_y_old(i,jmin,k)
           mhd_st_z(i,jmin,k) = mhd_st_z_old(i,jmin,k)
           Bxtilde(i,jmin,k) = Bxtilde_old(i,jmin,k)
           Bytilde(i,jmin,k) = Bytilde_old(i,jmin,k)
           Bztilde(i,jmin,k) = Bztilde_old(i,jmin,k)
        end do
     end do
     !$omp end parallel do
  end if
  if (have_bdry_max(3)==1) then
     !$omp parallel do
     do j=1,ext(2)
        do i=1,ext(1)
           rho_star(i,j,kmax) = rho_star_old(i,j,kmax)
           tau(i,j,kmax) = tau_old(i,j,kmax)
           mhd_st_x(i,j,kmax) = mhd_st_x_old(i,j,kmax)
           mhd_st_y(i,j,kmax) = mhd_st_y_old(i,j,kmax)
           mhd_st_z(i,j,kmax) = mhd_st_z_old(i,j,kmax)
           Bxtilde(i,j,kmax) = Bxtilde_old(i,j,kmax)
           Bytilde(i,j,kmax) = Bytilde_old(i,j,kmax)
           Bztilde(i,j,kmax) = Bztilde_old(i,j,kmax)
        end do
     end do
     !$omp end parallel do
  end if
  if (have_bdry_min(3)==1) then
     !$omp parallel do
     do j=1,ext(2)
        do i=1,ext(1)
           rho_star(i,j,kmin) = rho_star_old(i,j,kmin)
           tau(i,j,kmin) = tau_old(i,j,kmin)
           mhd_st_x(i,j,kmin) = mhd_st_x_old(i,j,kmin)
           mhd_st_y(i,j,kmin) = mhd_st_y_old(i,j,kmin)
           mhd_st_z(i,j,kmin) = mhd_st_z_old(i,j,kmin)
           Bxtilde(i,j,kmin) = Bxtilde_old(i,j,kmin)
           Bytilde(i,j,kmin) = Bytilde_old(i,j,kmin)
           Bztilde(i,j,kmin) = Bztilde_old(i,j,kmin)
        end do
     end do
     !$omp end parallel do
  end if
end subroutine apply_matter_frozen_bc
