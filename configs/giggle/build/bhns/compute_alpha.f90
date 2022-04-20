subroutine bhns_compute_alpha_1d(cctkGH,time,bh_posn_x,bh_posn_y,bh_posn_z,&
     N_varpi,N_phi,N_Z,Z_max,varpi_min,varpi_max,Symmetry,myproc,alpha_rho_cut_off)
  implicit none
!  INTEGER cctk_dim&&                           INTEGER cctk_gsh(cctk_dim),cctk_lsh(cctk_dim)&&                           INTEGER cc
!  REAL*8  excision_radius&&REAL*8  run_time&&INTEGER*4 Symmetry&&INTEGER*4 bssn_enable&&INTEGER*4 cowling_enable&&INTEGER*4 excisio
  external     CCTK_PointerTo
  integer*8 CCTK_PointerTo
  interface
  integer function CCTK_Equals (arg1, arg2)
  implicit none
  integer*8 arg1
  character(*) arg2
  end function CCTK_Equals
  integer function CCTK_MyProc (cctkGH)
  implicit none
  integer*8 cctkGH
  end function CCTK_MyProc
  integer function CCTK_nProcs (cctkGH)
  implicit none
  integer*8 cctkGH
  end function CCTK_nProcs
  integer function CCTK_IsThornActive (name)
  implicit none
  character(*) name
  end function CCTK_IsThornActive
  integer*8 function CCTK_NullPointer ()
  implicit none
  end function CCTK_NullPointer
  end interface
  interface
  INTEGER*4 function HorizonRadiusInDirection (horizon_number, N_points, x, y, z, radius)
  implicit none
  INTEGER*4 horizon_number
  INTEGER*4 N_points
  REAL*8 x(*)
  REAL*8 y(*)
  REAL*8 z(*)
  REAL*8 radius(*)
  end function HorizonRadiusInDirection
  end interface
  
  integer*8 :: cctkGH
  character                                :: filename_1d*50,filename_polar*50,filename_polar_avg*50
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  real*8 :: time, alpha_rho_cut_off
  real*8 :: rr_i, dS, PI, f1ospi
  real*8 :: sintheta,costheta,nn,xn,yn,phiangle,phiangle_BH,dphi,dcostheta
  integer :: i,j,k,interp_order,n,ntot
  integer :: vindex,myproc
  integer :: N_Z,N_phi,N_varpi, Symmetry
  real*8 :: Z_max,varpi_min,varpi_max
  real*8 :: varpi,ZL,dvarpi,dZ,integrand
  real*8 :: bh_posn_x,bh_posn_y,bh_posn_z
  real*8, allocatable, dimension(:,:)        :: pointcoords
  real*8, allocatable, dimension(:)          :: rho_b_int,P_int,h_int,vx_int,vy_int,vz_int,u0_int,lapm1_int
  real*8, allocatable, dimension(:)          :: smallb2_int,sbt_int,sbx_int,sby_int,sbz_int
  real*8, allocatable, dimension(:)          :: phi_int,shiftx_int,shifty_int,shiftz_int
  real*8, allocatable, dimension(:)          :: gxx_int,gxy_int,gxz_int,gyy_int,gyz_int,gzz_int
  real*8, allocatable, dimension(:)          :: gupxx_int,gupxy_int,gupxz_int,gupyy_int,gupyz_int,gupzz_int
!  real*8, allocatable, dimension(:)          :: sigma_xx_int,sigma_xy_int,sigma_xz_int,sigma_yy_int,sigma_yz_int,sigma_zz_int
!  real*8, allocatable, dimension(:)          :: eta_vis_int
  if (Symmetry .ne. EQUATORIAL) then
     write(*,*) 'Symmetry not supported in surface_density_profile'
     stop
  end if
  write(*,*) "inside compute_alpha"
  ntot = N_Z*N_phi*N_varpi
  ! allocate memory  
  allocate(pointcoords(ntot,3))
  allocate(rho_b_int(ntot))
  allocate(P_int(ntot))
  allocate(h_int(ntot))
  allocate(vx_int(ntot))
  allocate(vy_int(ntot))
  allocate(vz_int(ntot))
  allocate(u0_int(ntot))
  allocate(smallb2_int(ntot))
  allocate(sbt_int(ntot))
  allocate(sbx_int(ntot))
  allocate(sby_int(ntot))
  allocate(sbz_int(ntot))
  allocate(lapm1_int(ntot))
  allocate(phi_int(ntot))
  allocate(shiftx_int(ntot))
  allocate(shifty_int(ntot))
  allocate(shiftz_int(ntot))
  allocate(gxx_int(ntot))
  allocate(gxy_int(ntot))
  allocate(gxz_int(ntot))
  allocate(gyy_int(ntot))
  allocate(gyz_int(ntot))
  allocate(gzz_int(ntot))
  allocate(gupxx_int(ntot))
  allocate(gupxy_int(ntot))
  allocate(gupxz_int(ntot))
  allocate(gupyy_int(ntot))
  allocate(gupyz_int(ntot))
  allocate(gupzz_int(ntot))
!  allocate(sigma_xx_int(ntot))
!  allocate(sigma_xy_int(ntot))
!  allocate(sigma_xz_int(ntot))
!  allocate(sigma_yy_int(ntot))
!  allocate(sigma_yz_int(ntot))
!  allocate(sigma_zz_int(ntot))
!  allocate(eta_vis_int(ntot))
  PI = 3.14159265358979323844D0
  dphi = 2.0 * PI / N_phi
  dZ = Z_max / N_Z
  dvarpi = (varpi_max - varpi_min)/N_varpi
  n = 1
  do i=1,N_varpi
     varpi = varpi_min+(i-0.5)*dvarpi
     do j=1,N_phi
        phiangle = (j - 0.5)*dphi
        do k=1,N_Z
           ZL = dZ*(k-0.5)
           pointcoords(n,1) = varpi*cos(phiangle) + bh_posn_x
           pointcoords(n,2) = varpi*sin(phiangle) + bh_posn_y
           pointcoords(n,3) = ZL + bh_posn_z
           n = n + 1
        end do
     end do
  end do
!we compute b^2, for diagnostic purposes
!  call mhd_evolve_compute_b2_cpp(cctkGH,cctk_lsh, phi, lapm1, &
!       shiftx,shifty,shiftz,vx,vy,vz,Bx,By,Bz, &
!       gxx, gxy, gxz, gyy, gyz, gzz, temp10)
  ! Interpolate the grid function
  call CCTK_VarIndex(vindex,"mhd_evolve::rho_b")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,rho_b_int)
  call CCTK_VarIndex(vindex,"mhd_evolve::P")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,P_int)
  call CCTK_VarIndex(vindex,"mhd_evolve::h")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,h_int)
  call CCTK_VarIndex(vindex,"mhd_evolve::temp10")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,smallb2_int)
  call CCTK_VarIndex(vindex,"mhd_evolve::sbt")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,sbt_int)
  call CCTK_VarIndex(vindex,"mhd_evolve::sbx")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,sbx_int)
  call CCTK_VarIndex(vindex,"mhd_evolve::sby")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,sby_int)
  call CCTK_VarIndex(vindex,"mhd_evolve::sbz")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,sbz_int)
  call CCTK_VarIndex(vindex,"mhd_evolve::vx")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,vx_int)
  call CCTK_VarIndex(vindex,"mhd_evolve::vy")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,vy_int)
  call CCTK_VarIndex(vindex,"mhd_evolve::vz")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,vz_int)
  call CCTK_VarIndex(vindex,"mhd_evolve::u0")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,u0_int)
  call CCTK_VarIndex(vindex,"lapse::lapm1")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,lapm1_int)
  call CCTK_VarIndex(vindex,"bssn::phi")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,phi_int)
  call CCTK_VarIndex(vindex,"shift::shiftx")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,shiftx_int)
  call CCTK_VarIndex(vindex,"shift::shifty")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,shifty_int)
  call CCTK_VarIndex(vindex,"shift::shiftz")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,shiftz_int)
  call CCTK_VarIndex(vindex,"bssn::gxx")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gxx_int)
  call CCTK_VarIndex(vindex,"bssn::gxy")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gxy_int)
  call CCTK_VarIndex(vindex,"bssn::gxz")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gxz_int)
  call CCTK_VarIndex(vindex,"bssn::gyy")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gyy_int)
  call CCTK_VarIndex(vindex,"bssn::gyz")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gyz_int)
  call CCTK_VarIndex(vindex,"bssn::gzz")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gzz_int)
  call CCTK_VarIndex(vindex,"bssn::gupxx")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gupxx_int)
  call CCTK_VarIndex(vindex,"bssn::gupxy")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gupxy_int)
  call CCTK_VarIndex(vindex,"bssn::gupxz")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gupxz_int)
  call CCTK_VarIndex(vindex,"bssn::gupyy")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gupyy_int)
  call CCTK_VarIndex(vindex,"bssn::gupyz")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gupyz_int)
  call CCTK_VarIndex(vindex,"bssn::gupzz")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gupzz_int)
!  call CCTK_VarIndex(vindex,mhd_evolve::sigma_xx)
!  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,sigma_xx_int)
!  call CCTK_VarIndex(vindex,mhd_evolve::sigma_xy)
!  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,sigma_xy_int)
!  call CCTK_VarIndex(vindex,mhd_evolve::sigma_xz)
!  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,sigma_xz_int)
!  call CCTK_VarIndex(vindex,mhd_evolve::sigma_yy)
!  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,sigma_yy_int)
!  call CCTK_VarIndex(vindex,mhd_evolve::sigma_yz)
!  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,sigma_yz_int)
!  call CCTK_VarIndex(vindex,mhd_evolve::sigma_zz)
!  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,sigma_zz_int)
!  call CCTK_VarIndex(vindex,mhd_evolve::eta_viscosity)
!  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,eta_vis_int)
  if (myproc.eq.0) then
  call bhns_compute_alpha_IO(cctkGH,time,&
          varpi_min,dvarpi,dZ, &
          N_varpi,N_phi,N_Z,&
          rho_b_int,P_int,h_int,vx_int,vy_int,vz_int,u0_int,&
          sbt_int,sbx_int,sby_int,sbz_int,&
          lapm1_int,shiftx_int,shifty_int,shiftz_int,&
          phi_int,gxx_int,gxy_int,gxz_int,&
          gyy_int,gyz_int, gzz_int,&
          gupxx_int,gupxy_int,gupxz_int,&
          gupyy_int,gupyz_int, gupzz_int,bh_posn_x,bh_posn_y,bh_posn_z,alpha_rho_cut_off)
!          sigma_xx_int,sigma_xy_int,sigma_xz_int,&
!          sigma_yy_int,sigma_yz_int,sigma_zz_int,&
!          eta_vis_int)
  endif
 ! deallocate memory
 deallocate(pointcoords)
 deallocate(rho_b_int,P_int,h_int,vx_int,vy_int,vz_int,u0_int,smallb2_int,sbt_int,sbx_int,sby_int,sbz_int,lapm1_int,&
  phi_int,shiftx_int,shifty_int,shiftz_int,gxx_int,gxy_int,gxz_int,&
  gyy_int,gyz_int,gzz_int,gupxx_int,gupxy_int,gupxz_int,gupyy_int,gupyz_int,gupzz_int)
!  deallocate(sigma_xx_int,sigma_xy_int,sigma_xz_int,sigma_yy_int,&
!             sigma_yz_int,sigma_zz_int,eta_vis_int)
end subroutine bhns_compute_alpha_1d
