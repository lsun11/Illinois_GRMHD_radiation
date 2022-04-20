#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

! This version of fluxes uses an analytic determination of the event
! horizon radius.  For use with Cowling runs--such as disks around black
! holes a la McKinney and Gammie. 

! The strange thing about this is that the output is placed in
! steerable parameters instead of grid variables.
subroutine mag_bondi_accretion_fluxes_ana(cctkGH, Symmetry, N_theta, N_phi, & 
					     xn_d, yn_d, zn_d, nx_d,ny_d,nz_d, &
                                             F_M0,F_E_fluid,F_E_em, & 
                                             F_J_fluid,F_J_em)
  implicit none
  DECLARE_CCTK_FUNCTIONS
  CCTK_POINTER :: cctkGH

  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  integer :: Symmetry, N_theta, N_phi
  real*8, dimension(N_theta*N_phi) :: xn_d, yn_d, zn_d, nx_d,ny_d,nz_d
  integer, parameter :: AXISYM_FULL = 5
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  real*8 :: dphi, dcostheta, dOmega
  real*8 :: F_M0,F_E_fluid,F_E_em,F_J_fluid,F_J_em
  real*8 :: F_M0_ijk, F_E_fluid_ijk, F_E_em_ijk, F_J_fluid_ijk, F_J_em_ijk
  real*8 :: psi6_i, psi4_i, lap_i, fact, btx_i, bty_i, btz_i,psi2_i
  real*8 :: gxx_i, gxy_i, gxz_i, gyy_i, gyz_i, gzz_i
  real*8 :: bt_i, bx_i, by_i, bz_i, ut_i, vx_i, vy_i, vz_i
  real*8 :: rhos_i, h_i, P_i, bt_x_i, bt_y_i, bt_z_i
  real*8 :: bb, bv, Bbeta, u_t_i, b_t_i, u_x_i, u_y_i, u_z_i, u_phi_i
  real*8 :: b_x_i, b_y_i, b_phi_i, BBx_i,BBy_i,BBz_i
  real*8 :: F_E_emx, F_E_emy, F_E_emz, F_J_fx, F_J_fy, F_J_fz
  real*8 :: F_J_emx, F_J_emy, F_J_emz, rr_i, dS, PI, f1ospi
  real*8 :: F_E_emxt, F_E_emyt, F_E_emzt, v2
  real*8, parameter :: max_gamma = 60.d0
  real*8 :: v2m
  integer :: i, interp_order, n_tot
  integer, dimension(3)                     :: ext,global_ext
  real*8, allocatable, dimension(:)         :: phiint, vxint, vyint, vzint
  real*8, allocatable, dimension(:)         :: gintxx,gintxy,gintxz,gintyy,gintyz,gintzz
  real*8, allocatable, dimension(:)         :: lapseint,shiftintx,shiftinty,shiftintz
  real*8, allocatable, dimension(:)         :: rhosint, hint, Pint
  real*8, allocatable, dimension(:)         :: Bxint, Byint, Bzint


  PI = 3.14159265358979323846D0
  dphi = 2.d0 * PI / N_phi
  if (Symmetry==AXISYM_FULL .or. Symmetry==NO_SYMM) then 
     dcostheta = 2.d0 / N_theta
  else
     dcostheta = 1.d0 / N_theta
  end if
  interp_order = 2

  v2m = 1.d0 - 1.d0/max_gamma**2

  !!if (Symmetry==EQUATORIAL) then 
  !!   sym_factor = 2.d0
  !!else if (Symmetry==OCTANT) then 
  !!   sym_factor = 4.d0
  !!else 
  !!   sym_factor = 1.d0
  !!end if

  if (Symmetry==EQUATORIAL) then
     dOmega = 2.d0*dcostheta*dphi
  else if (Symmetry==OCTANT) then
     dOmega = 4.d0*dcostheta*dphi
  else if (Symmetry==AXISYM) then
     dOmega = 4.d0*PI*dcostheta
  else if (Symmetry==AXISYM_FULL) then
     dOmega = 2.d0*PI*dcostheta
  else
     dOmega = dcostheta*dphi
  end if

  n_tot = N_theta*N_phi

  ! Allocate memory
  allocate(phiint(n_tot))
  allocate(vxint(n_tot))
  allocate(vyint(n_tot))
  allocate(vzint(n_tot))
  allocate(gintxx(n_tot))
  allocate(gintxy(n_tot))
  allocate(gintxz(n_tot))
  allocate(gintyy(n_tot))
  allocate(gintyz(n_tot))
  allocate(gintzz(n_tot))
  allocate(lapseint(n_tot))
  allocate(shiftintx(n_tot))
  allocate(shiftinty(n_tot))
  allocate(shiftintz(n_tot))
  allocate(rhosint(n_tot))
  allocate(hint(n_tot))
  allocate(Pint(n_tot))
  allocate(Bxint(n_tot))
  allocate(Byint(n_tot))
  allocate(Bzint(n_tot))

  ! call interpolation routine.
  call interpolation_for_mag_bondi_diagnostics(cctkGH, xn_d,yn_d,zn_d, &
                        phiint,lapseint,vxint,vyint,vzint, gintxx,gintxy,gintxz, &
                        gintyy,gintyz,gintzz,shiftintx,shiftinty,shiftintz, &
                        rhosint,hint,Pint,Bxint,Byint,Bzint, n_tot)
  !
  ! integrate on the horizon
  !

  F_M0 = 0.d0
  F_E_fluid = 0.d0
  F_E_em = 0.d0
  F_J_fluid = 0.d0
  F_J_em = 0.d0

  F_E_emxt = 0.d0
  F_E_emyt = 0.d0
  F_E_emzt = 0.d0

  f1ospi = 1.d0/sqrt(4.d0*PI)
  do i=1,n_tot 
     psi2_i = exp(2.d0*phiint(i))
     psi4_i = psi2_i*psi2_i
     psi6_i = psi2_i*psi4_i
     lap_i = lapseint(i) + 1.d0
     btx_i = shiftintx(i)
     bty_i = shiftinty(i)
     btz_i = shiftintz(i)
     gxx_i = gintxx(i)*psi4_i
     gxy_i = gintxy(i)*psi4_i
     gxz_i = gintxz(i)*psi4_i
     gyy_i = gintyy(i)*psi4_i
     gyz_i = gintyz(i)*psi4_i
     gzz_i = gintzz(i)*psi4_i
     vx_i  = vxint(i)
     vy_i  = vyint(i)
     vz_i  = vzint(i)
     rhos_i = rhosint(i)
     h_i  = hint(i)
     P_i  = Pint(i)
     BBx_i = Bxint(i)
     BBy_i = Byint(i)
     BBz_i = Bzint(i)

     v2 = (gxx_i*(vx_i+btx_i)**2 + 2.d0*gxy_i*(vx_i+btx_i)*(vy_i+bty_i) + & 
		2.d0*gxz_i*(vx_i+btx_i)*(vz_i+btz_i) + gyy_i*(vy_i+bty_i)**2 + &
		2.d0*gyz_i*(vy_i+bty_i)*(vz_i+btz_i) + gzz_i*(vz_i+btz_i)**2 )/lap_i**2
     ! Restrict v2 to prevent superluminal velocity
     if (v2 .gt. v2m) then 
        fact = sqrt(v2m/v2) 
	vx_i = (vx_i+btx_i)*fact - btx_i
        vy_i = (vy_i+bty_i)*fact - bty_i
        vz_i = (vz_i+btz_i)*fact - btz_i
	v2 = v2m
     end if
     ut_i = 1.d0/sqrt(1.d0-v2)/lap_i

     ! \beta_i
     bt_x_i = gxx_i*btx_i + gxy_i*bty_i + gxz_i*btz_i
     bt_y_i = gxy_i*btx_i + gyy_i*bty_i + gyz_i*btz_i
     bt_z_i = gxz_i*btx_i + gyz_i*bty_i + gzz_i*btz_i
     ! \beta.\beta
     bb = btx_i*bt_x_i + bty_i*bt_y_i + btz_i*bt_z_i
     ! \beta.v
     bv = vx_i*bt_x_i + vy_i*bt_y_i + vz_i*bt_z_i
     ! \beta.b
     Bbeta = bx_i*bt_x_i + by_i*bt_y_i + bz_i*bt_z_i

     ! u_t
     u_t_i = ut_i*(-lap_i*lap_i + bb + bv)

     ! b_t
     b_t_i = bt_i*(-lap_i*lap_i + bb) + Bbeta

     ! u_x, u_y, and u_phi
     u_x_i = ut_i*( gxx_i*(vx_i + btx_i) + gxy_i*(vy_i + bty_i) &
          + gxz_i*(vz_i + btz_i) )
     u_y_i = ut_i*( gxy_i*(vx_i + btx_i) + gyy_i*(vy_i + bty_i) &
          + gyz_i*(vz_i + btz_i) )
     u_z_i = ut_i*( gxz_i*(vx_i + btx_i) + gyz_i*(vy_i + bty_i) &
          + gzz_i*(vz_i + btz_i) )
     u_phi_i = xn_d(i)*u_y_i - yn_d(i)*u_x_i

     ! b^mu
     fact = f1ospi/lap_i
     bt_i = (u_x_i*BBx_i + u_y_i*BBy_i + u_z_i*BBz_i) * fact
     bx_i = BBx_i*fact/ut_i + bt_i*vx_i
     by_i = BBy_i*fact/ut_i + bt_i*vy_i
     bz_i = BBz_i*fact/ut_i + bt_i*vz_i

     ! b_x, b_y, and b_phi 
     b_x_i = bt_x_i*bt_i + gxx_i*bx_i + gxy_i*by_i + gxz_i*bz_i
     b_y_i = bt_y_i*bt_i + gxy_i*bx_i + gyy_i*by_i + gyz_i*bz_i
     b_phi_i = xn_d(i)*b_y_i - yn_d(i)*b_x_i

     ! b.b
     bb = -lap_i*lap_i*bt_i*bt_i &
          + gxx_i*(bx_i+btx_i*bt_i)*(bx_i+btx_i*bt_i) &
          + gyy_i*(by_i+bty_i*bt_i)*(by_i+bty_i*bt_i) &
          + gzz_i*(bz_i+btz_i*bt_i)*(bz_i+btz_i*bt_i) &
          + 2.0*( gxy_i*(bx_i+btx_i*bt_i)*(by_i+bty_i*bt_i) &
          + gxz_i*(bx_i+btx_i*bt_i)*(bz_i+btz_i*bt_i) &
          + gyz_i*(by_i+bty_i*bt_i)*(bz_i+btz_i*bt_i) )

     F_E_emx = bb*ut_i*vx_i*u_t_i - bx_i*b_t_i
     F_E_emy = bb*ut_i*vy_i*u_t_i - by_i*b_t_i
     F_E_emz = bb*ut_i*vz_i*u_t_i - bz_i*b_t_i

     F_J_fx = rhos_i*h_i*vx_i*u_phi_i - lap_i*psi6_i*P_i*yn_d(i)
     F_J_fy = rhos_i*h_i*vy_i*u_phi_i + lap_i*psi6_i*P_i*xn_d(i)
     F_J_fz = rhos_i*h_i*vz_i*u_phi_i

     F_J_emx = (ut_i*vx_i*u_phi_i - 0.5*yn_d(i))*bb - bx_i*b_phi_i
     F_J_emy = (ut_i*vy_i*u_phi_i + 0.5*xn_d(i))*bb - by_i*b_phi_i
     F_J_emz = ut_i*vz_i*u_phi_i*bb - bz_i*b_phi_i

     F_M0_ijk = (nx_d(i)*vx_i + ny_d(i)*vy_i + nz_d(i)*vz_i)*rhos_i
     F_E_fluid_ijk = -F_M0_ijk*h_i*u_t_i
     F_E_em_ijk = -(nx_d(i)*F_E_emx + ny_d(i)*F_E_emy + nz_d(i)*F_E_emz)*lap_i*psi6_i
     F_J_fluid_ijk = (nx_d(i)*F_J_fx + ny_d(i)*F_J_fy + nz_d(i)*F_J_fz)
     F_J_em_ijk = (nx_d(i)*F_J_emx + ny_d(i)*F_J_emy + nz_d(i)*F_J_emz)*lap_i*psi6_i

     rr_i = xn_d(i)*xn_d(i) + yn_d(i)*yn_d(i) + zn_d(i)*zn_d(i)

     !!if (Symmetry==AXISYM) then
     !!   dS = rr_i * 4.d0*PI * dcostheta
     !!elseif (Symmetry==AXISYM_FULL) then 
     !!   dS = rr_i * 2.d0*PI * dcostheta
     !!else
     !!   dS = rr_i * dphi * dcostheta * sym_factor
     !!end if

     dS = rr_i * dOmega

     F_M0      = F_M0      + F_M0_ijk      * dS
     F_E_fluid = F_E_fluid + F_E_fluid_ijk * dS
     F_E_em    = F_E_em    + F_E_em_ijk    * dS
     F_J_fluid = F_J_fluid + F_J_fluid_ijk * dS
     F_J_em    = F_J_em    + F_J_em_ijk    * dS

     F_E_emxt = F_E_emxt -nx_d(i)*F_E_emx*lap_i*psi6_i*dS
     F_E_emyt = F_E_emyt -ny_d(i)*F_E_emy*lap_i*psi6_i*dS
     F_E_emzt = F_E_emzt -nz_d(i)*F_E_emz*lap_i*psi6_i*dS

     !      write(*,100) bb, ut_i, vz_i, u_t_i, bz_i, b_t_i
     ! Uncomment below line to enable velocity output
     !      write(*,100) ut_i, vx_i, vy_i, vz_i

  end do

  deallocate(phiint, vxint, vyint, vzint)
  deallocate(gintxx,gintxy,gintxz,gintyy,gintyz,gintzz)
  deallocate(lapseint,shiftintx,shiftinty,shiftintz)
  deallocate(rhosint, hint, Pint)
  deallocate(Bxint, Byint, Bzint)

end subroutine mag_bondi_accretion_fluxes_ana

subroutine interpolation_for_mag_bondi_diagnostics(cctkGH, xinterp,yinterp,zinterp, & 
			phiint,lapseint,vxint,vyint,vzint, gintxx,gintxy,gintxz, & 
			gintyy,gintyz,gintzz,shiftintx,shiftinty,shiftintz, &
			rhosint,hint,Pint,Bxint,Byint,Bzint, ntot)
  implicit none
  DECLARE_CCTK_FUNCTIONS
  CCTK_POINTER                             :: cctkGH
  integer :: ntot,i
  real*8, dimension(ntot) :: xinterp,yinterp,zinterp
  real*8, dimension(ntot) :: phiint,lapseint,vxint,vyint,vzint
  real*8, dimension(ntot) :: gintxx,gintxy,gintxz,gintyy,gintyz,gintzz
  real*8, dimension(ntot) :: shiftintx,shiftinty,shiftintz,rhosint,hint,Pint
  real*8, dimension(ntot) :: Bxint,Byint,Bzint
  integer, parameter :: N_dims = 3, interpolation_order = 2
  integer, parameter :: N_input_arrays = 20, N_output_arrays = 20
  integer,dimension(N_input_arrays)  :: input_array_type_codes,input_array_varindices
  integer,dimension(N_output_arrays) :: output_array_type_codes
  character(60)                             :: options_string
  CCTK_POINTER,dimension(N_output_arrays)   :: output_array_pointers
  CCTK_POINTER, dimension(3)                :: interp_coords
  integer :: ierr,interp_handle,param_table_handle,coord_system_handle
!
  interp_handle = -1
  call CCTK_InterpHandle (interp_handle, "Lagrange polynomial interpolation")
  !  call CCTK_InterpHandle (interp_handle, "uniform cartesian")
  if (interp_handle .lt. 0) then
     call CCTK_WARN(0,"Cannot get handle for interpolation ! Forgot to activate an implementation providing interpolation operators ??")
  endif

  param_table_handle = -1
  options_string = "order = " // char(ichar('0') + interpolation_order)
  call Util_TableCreateFromString (param_table_handle, options_string)
  if (param_table_handle .lt. 0) then
     call CCTK_WARN(0,"Cannot create parameter table for interpolator")
  endif

  coord_system_handle = -1
  call CCTK_CoordSystemHandle (coord_system_handle, "cart3d")
  if (coord_system_handle .lt. 0) then
     call CCTK_WARN(0,"Cannot get handle for cart3d coordinate system ! Forgot to activate an implementation providing coordinates ??")
  endif

  input_array_type_codes = CCTK_VARIABLE_REAL
  output_array_type_codes = CCTK_VARIABLE_REAL

  ! Specify interpolation input arrays, output arrays:
  call CCTK_VarIndex (input_array_varindices(1), "bssn::phi")
  call CCTK_VarIndex (input_array_varindices(2), "lapse::lapm1")
  call CCTK_VarIndex (input_array_varindices(3), "mhd_evolve::vx")
  call CCTK_VarIndex (input_array_varindices(4), "mhd_evolve::vy")
  call CCTK_VarIndex (input_array_varindices(5), "mhd_evolve::vz")
  call CCTK_VarIndex (input_array_varindices(6), "bssn::gxx")
  call CCTK_VarIndex (input_array_varindices(7), "bssn::gxy")
  call CCTK_VarIndex (input_array_varindices(8), "bssn::gxz")
  call CCTK_VarIndex (input_array_varindices(9), "bssn::gyy")
  call CCTK_VarIndex (input_array_varindices(10),"bssn::gyz")
  call CCTK_VarIndex (input_array_varindices(11),"bssn::gzz")
  call CCTK_VarIndex (input_array_varindices(12),"shift::shiftx")
  call CCTK_VarIndex (input_array_varindices(13),"shift::shifty")
  call CCTK_VarIndex (input_array_varindices(14),"shift::shiftz")
  call CCTK_VarIndex (input_array_varindices(15),"mhd_evolve::rho_star")
  call CCTK_VarIndex (input_array_varindices(16),"mhd_evolve::h")
  call CCTK_VarIndex (input_array_varindices(17),"mhd_evolve::P")
  call CCTK_VarIndex (input_array_varindices(18),"mhd_evolve::Bx")
  call CCTK_VarIndex (input_array_varindices(19),"mhd_evolve::By")
  call CCTK_VarIndex (input_array_varindices(20),"mhd_evolve::Bz")

  output_array_pointers(1) = CCTK_PointerTo(phiint)
  output_array_pointers(2) = CCTK_PointerTo(lapseint)
  output_array_pointers(3) = CCTK_PointerTo(vxint)
  output_array_pointers(4) = CCTK_PointerTo(vyint)
  output_array_pointers(5) = CCTK_PointerTo(vzint)
  output_array_pointers(6) = CCTK_PointerTo(gintxx)
  output_array_pointers(7) = CCTK_PointerTo(gintxy)
  output_array_pointers(8) = CCTK_PointerTo(gintxz)
  output_array_pointers(9) = CCTK_PointerTo(gintyy)
  output_array_pointers(10)= CCTK_PointerTo(gintyz)
  output_array_pointers(11)= CCTK_PointerTo(gintzz)
  output_array_pointers(12)= CCTK_PointerTo(shiftintx)
  output_array_pointers(13)= CCTK_PointerTo(shiftinty)
  output_array_pointers(14)= CCTK_PointerTo(shiftintz)
  output_array_pointers(15)= CCTK_PointerTo(rhosint)
  output_array_pointers(16)= CCTK_PointerTo(hint)
  output_array_pointers(17)= CCTK_PointerTo(Pint)
  output_array_pointers(18)= CCTK_PointerTo(Bxint)
  output_array_pointers(19)= CCTK_PointerTo(Byint)
  output_array_pointers(20)= CCTK_PointerTo(Bzint)

  interp_coords(1) = CCTK_PointerTo(xinterp)
  interp_coords(2) = CCTK_PointerTo(yinterp)
  interp_coords(3) = CCTK_PointerTo(zinterp)

  ! Perform interpolation:
  call CCTK_InterpGridArrays(ierr,cctkGH,N_dims,interp_handle, &
       param_table_handle,coord_system_handle, &
       ntot,input_array_type_codes,interp_coords, &
       N_input_arrays,input_array_varindices, &
       N_output_arrays, output_array_type_codes, output_array_pointers)

end subroutine interpolation_for_mag_bondi_diagnostics
