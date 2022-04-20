#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bhns_surface_dens_1d(cctkGH,time,bh_posn_x,bh_posn_y,bh_posn_z,&
      N_varpi,N_phi,N_Z,Z_max,varpi_min,varpi_max,Symmetry,myproc)
  implicit none
!  DECLARE_CCTK_ARGUMENTS
!  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_POINTER :: cctkGH
  character                                :: filename_1d*50,filename_polar*50,filename_polar_avg*50
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  real*8 :: time
  real*8 :: Sigma_avgL,SigmaL,torque_densL
  real*8 :: rr_i, dS, PI, f1ospi
  real*8 :: cosphi,sinphi,sym_factor,Lxph,Lyph,dPhi_dphi
  real*8 :: sintheta,costheta,nn,xn,yn,phiangle,phiangle_BH,dphi,dcostheta
  integer :: i,j,k,interp_order,n,ntot
  integer :: ind0,ind1,ind2,indm1,indm2,vindex,myproc
  integer :: N_Z,N_phi,N_varpi, Symmetry
  real*8 :: Z_max,varpi_min,varpi_max
  real*8 :: varpi,ZL,dvarpi,dZ,integrand
  real*8 :: bh_posn_x,bh_posn_y,bh_posn_z
  real*8, allocatable, dimension(:,:)        :: pointcoords
  !real*8, dimension(N_Z*N_phi*N_varpi)          :: nx_d,ny_d,nz_d
  !real*8, dimension(N_Z*N_phi*N_varpi)          :: ah_radii, x_ah1,y_ah1,z_ah1
  !real*8, dimension(N_Z*N_phi*N_varpi)          :: vxint, vyint, vzint
  real*8  :: u_xL,u_yL
  real*8  :: shiftx_intL,shifty_intL,shiftz_intL
  real*8  :: shift_xL,shift_yL,shift_zL
  real*8  :: g4xxx,g4xxy,g4xyx,g4xyy,g4xzx,g4xzy,g4yyx,g4yyy,g4yzx,g4yzy,g4zzx,g4zzy
  real*8  :: g4tt,g4tx,g4ty,g4tz,g4xx,g4xy,g4xz,g4yy,g4yz,g4zz
  real*8  :: dx_gxj_shiftj,dx_gyj_shiftj,dx_gzj_shiftj,dy_gxj_shiftj,dy_gyj_shiftj,dy_gzj_shiftj
  real*8  :: dx_shiftj_gxj,dx_shiftj_gyj,dx_shiftj_gzj,dy_shiftj_gxj,dy_shiftj_gyj,dy_shiftj_gzj
  real*8  :: g4txx,g4tyx,g4tzx,g4txy,g4tyy,g4tzy
  real*8  :: g4ttx,g4tty,g4ttz
  
  real*8, allocatable, dimension(:)          :: rho_b_int,P_int,h_int,vx_int,vy_int,vz_int,u0_int,lapm1_int,lapsex_int
  real*8, allocatable, dimension(:)          :: lapsey_int,phi_int,phix_int,phiy_int,shiftx_int,shifty_int,shiftz_int,shiftxx_int
  real*8, allocatable, dimension(:)          :: shiftyx_int,shiftzx_int,shiftxy_int,shiftyy_int,shiftzy_int,gxx_int,gxy_int,gxz_int
  real*8, allocatable, dimension(:)          :: gyy_int,gyz_int,gzz_int,gupxx_int,gupxy_int,gupxz_int,gupyy_int,gupyz_int
  real*8, allocatable, dimension(:)          :: gupzz_int,gxxx_int,gxyx_int,gxzx_int,gyyx_int,gyzx_int,gzzx_int,gxxy_int
  real*8, allocatable, dimension(:)          :: gxyy_int,gxzy_int,gyyy_int,gyzy_int,gzzy_int
  if (Symmetry .ne. EQUATORIAL) then 
     write(*,*) 'Symmetry not supported in surface_density_profile'
     stop
  end if 
  write(*,*) "inside 1d"
  sym_factor = 2.d0
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
  allocate(lapm1_int(ntot))
  allocate(lapsex_int(ntot))
  allocate(lapsey_int(ntot))
  allocate(phi_int(ntot))
  allocate(phix_int(ntot))
  allocate(phiy_int(ntot))  
  allocate(shiftx_int(ntot))
  allocate(shifty_int(ntot))
  allocate(shiftz_int(ntot))
  allocate(shiftxx_int(ntot))
  allocate(shiftyx_int(ntot))
  allocate(shiftzx_int(ntot))
  allocate(shiftxy_int(ntot))
  allocate(shiftyy_int(ntot))
  allocate(shiftzy_int(ntot))
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
  allocate(gxxx_int(ntot))
  allocate(gxyx_int(ntot))
  allocate(gxzx_int(ntot))
  allocate(gyyx_int(ntot))
  allocate(gyzx_int(ntot))
  allocate(gzzx_int(ntot))
  allocate(gxxy_int(ntot))
  allocate(gxyy_int(ntot))
  allocate(gxzy_int(ntot))
  allocate(gyyy_int(ntot))
  allocate(gyzy_int(ntot))
  allocate(gzzy_int(ntot))
 
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
  ! Interpolate the grid function
  call CCTK_VarIndex(vindex,"mhd_evolve::rho_b")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,rho_b_int)
  call CCTK_VarIndex(vindex,"mhd_evolve::P")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,P_int)
  call CCTK_VarIndex(vindex,"mhd_evolve::h")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,h_int)


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
  call CCTK_VarIndex(vindex,"lapse::lapsex")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,lapsex_int)
  call CCTK_VarIndex(vindex,"lapse::lapsey")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,lapsey_int)

  call CCTK_VarIndex(vindex,"bssn::phi")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,phi_int)
  call CCTK_VarIndex(vindex,"bssn::phix")	
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,phix_int)
  call CCTK_VarIndex(vindex,"bssn::phiy")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,phiy_int)
  
  call CCTK_VarIndex(vindex,"shift::shiftx")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,shiftx_int)
  call CCTK_VarIndex(vindex,"shift::shifty")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,shifty_int)
  call CCTK_VarIndex(vindex,"shift::shiftz")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,shiftz_int)

  call CCTK_VarIndex(vindex,"mhd_evolve::temp1")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,shiftxx_int)
  call CCTK_VarIndex(vindex,"mhd_evolve::temp4")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,shiftyx_int)
  call CCTK_VarIndex(vindex,"mhd_evolve::temp7")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,shiftzx_int)
  call CCTK_VarIndex(vindex,"mhd_evolve::temp2")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,shiftxy_int)
  call CCTK_VarIndex(vindex,"mhd_evolve::temp5")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,shiftyy_int)
  call CCTK_VarIndex(vindex,"mhd_evolve::temp8")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,shiftzy_int)
     
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


  call CCTK_VarIndex(vindex,"bssn::gxxx")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gxxx_int)
  call CCTK_VarIndex(vindex,"bssn::gxxy")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gxxy_int)
  call CCTK_VarIndex(vindex,"bssn::gxyx")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gxyx_int)
  call CCTK_VarIndex(vindex,"bssn::gxyy")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gxyy_int)
  call CCTK_VarIndex(vindex,"bssn::gxzx")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gxzx_int)
  call CCTK_VarIndex(vindex,"bssn::gxzy")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gxzy_int)
  call CCTK_VarIndex(vindex,"bssn::gyyx")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gyyx_int)
  call CCTK_VarIndex(vindex,"bssn::gyyy")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gyyy_int)
  call CCTK_VarIndex(vindex,"bssn::gyzx")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gyzx_int)
  call CCTK_VarIndex(vindex,"bssn::gyzy")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gyzy_int)
  call CCTK_VarIndex(vindex,"bssn::gzzx")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gzzx_int)
  call CCTK_VarIndex(vindex,"bssn::gzzy")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gzzy_int)

  if (myproc.eq.0) then
     call bhns_surf_dens_1d_IO(cctkGH,time,&
          varpi_min,dvarpi,dZ, &
          N_varpi,N_phi,N_Z,&
          rho_b_int,P_int,h_int,vx_int,vy_int,vz_int,u0_int,&
          lapm1_int,shiftx_int,shifty_int,shiftz_int,&
          lapsex_int,shiftxx_int,shiftyx_int,shiftzx_int,&
          lapsey_int,shiftxy_int,shiftyy_int,shiftzy_int,&
          phi_int,gxx_int,gxy_int,gxz_int,&
          gyy_int,gyz_int, gzz_int,&
          gupxx_int,gupxy_int,gupxz_int,&
          gupyy_int,gupyz_int, gupzz_int,&
          phix_int,gxxx_int,gxyx_int,gxzx_int,&
          gyyx_int,gyzx_int, gzzx_int,&
          phiy_int,gxxy_int,gxyy_int,gxzy_int,&
          gyyy_int,gyzy_int, gzzy_int,bh_posn_x,bh_posn_y,bh_posn_z)
  endif

 ! deallocate memory
 deallocate(pointcoords)
 deallocate(rho_b_int,P_int,h_int,vx_int,vy_int,vz_int,u0_int,lapm1_int,lapsex_int,&
  lapsey_int,phi_int,phix_int,phiy_int,shiftx_int,shifty_int,shiftz_int,shiftxx_int,&
  shiftyx_int,shiftzx_int,shiftxy_int,shiftyy_int,shiftzy_int,gxx_int,gxy_int,gxz_int,&
  gyy_int,gyz_int,gzz_int,gupxx_int,gupxy_int,gupxz_int,gupyy_int,gupyz_int,&
  gupzz_int,gxxx_int,gxyx_int,gxzx_int,gyyx_int,gyzx_int,gzzx_int,gxxy_int,&
  gxyy_int,gxzy_int,gyyy_int,gyzy_int,gzzy_int)
end subroutine bhns_surface_dens_1d
