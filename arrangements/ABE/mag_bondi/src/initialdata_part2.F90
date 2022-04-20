#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine mag_bondi_initialdata_part2(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext
  real*8 			           :: dT,dX,dY,dZ,rhos_max
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
  integer :: ierr,ONE,ZERO, i,j,k
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  integer, parameter :: AXISYM_FULL = 5
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  CCTK_REAL reduction_value

  parameter(ONE = 1.D0, ZERO = 0.D0)
!
  if (abs(Sym_Bz +1.d0) .gt. 1.d-10) then 
     write(*,*) 'Sym_Bz = ',Sym_Bz 
     write(*,*) 'Sym_Bz must be set to -1 to run the magnetized Bondi accretion!'
     stop
  end if

  !!if (puncture_id==1 .and. cowling_enable==1) then 
  !!   write(*,*) 'You should not turn on Cowling when using puncture ID.'
  !!   stop
  !!end if

  !!if (puncture_id==1 .and. self_gravity==1) then 
  !!   write(*,*) 'The parameter self_gravity should be set to 0'
  !!   stop
  !!end if

  if (abs(sam) .gt. 1.d-10) write(*,*) 'Warning: BH spin = ',sam

  ext = cctk_lsh

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  !|~~~~~> Set up *untilded* metric, extrinsic curvature, lapse,shifts
  if (puncture_id==1) then 
     psi = 1.d0+0.5d0/sqrt(X**2+Y**2+Z**2)
     phi = log(psi)
     gxx = 1.d0
     gxy = 0.d0
     gxz = 0.d0
     gyy = 1.d0
     gyz = 0.d0
     gzz = 1.d0
     gupxx = 1.d0
     gupxy = 0.d0
     gupxz = 0.d0
     gupyy = 1.d0
     gupyz = 0.d0
     gupzz = 1.d0
     trK = 0.d0
     Axx = 0.d0
     Axy = 0.d0
     Axz = 0.d0
     Ayy = 0.d0
     Ayz = 0.d0
     Azz = 0.d0
     lapm1 = 1.d0/psi**2 - 1.d0
     shiftx = 0.d0
     shifty = 0.d0
     shiftz = 0.d0
     lapset = 0.d0
     shiftxt = 0.d0
     shiftyt = 0.d0
     shiftzt = 0.d0
     Gammax = 0.d0
     Gammay = 0.d0
     Gammaz = 0.d0
     call Derivs(ext,X,Y,Z,dX,dY,dZ,phi,phix,phiy,phiz,Symmetry)
     call Derivs(ext,X,Y,Z,dX,dY,dZ,lapm1,lapsex,lapsey,lapsez,Symmetry)
     call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
     call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')
  else
     call ks_initial_metric_shift_origin_bl(ext,X,Y,Z, &
          gxx,gxy,gxz,gyy,gyz,gzz,trK,Axx,Axy,Axz,Ayy,Ayz,Azz, &
          lapm1,shiftx,shifty,shiftz,sam,r0,Symmetry)

     call convert(ext,phi,gxx,gxy,gxz,gyy,gyz,gzz, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
          psi,detmin_l,detmax_l)
     psi = psi**(1.d0/12.d0)

    call setgamma_v2(ext, cctk_nghostzones, dX, dY, dZ, &
        phi, gxx, gxy, gxz, gyy, gyz, gzz, &
        gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Gammax, Gammay, Gammaz, rho, &
        Symmetry)

     call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')

     n1 = 0 
     n2 = 0 
     n3 = 0 
     mf = 1    
     ! Use rho as temporary storage here:
     call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammax,rho,n1,n2,n3,mf,Symmetry)
     call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammay,rho,n1,n2,n3,mf,Symmetry)
     call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammaz,rho,n1,n2,n3,mf,Symmetry)

     call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')

     call Derivs(ext,X,Y,Z,dX,dY,dZ,phi,phix,phiy,phiz,Symmetry)
     call Derivs(ext,X,Y,Z,dX,dY,dZ,lapm1,lapsex,lapsey,lapsez,Symmetry)
     call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
     call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')

     !|~~~~~> convert traceless K_{ij} to tilded traceless K_{ij}
     Axx = Axx*exp(-4.D0*phi)
     Axy = Axy*exp(-4.D0*phi)
     Axz = Axz*exp(-4.D0*phi)
     Ayy = Ayy*exp(-4.D0*phi)
     Ayz = Ayz*exp(-4.D0*phi)
     Azz = Azz*exp(-4.D0*phi)
  end if

  !|~~~~~> Set everything else to Zero!
  rho = 0.d0
  S = 0.d0
  Sx = 0.d0
  Sy = 0.d0
  Sz = 0.d0
  Sxx = 0.d0
  Sxy = 0.d0
  Sxz = 0.d0
  Syy = 0.d0
  Syz = 0.d0
  Szz = 0.d0

  Pl = 0.d0

  rho_b_atm = rho_fact
  tau_atm = tau_fact
  rho_b_max = 1.d0

  ! Setup excision_zone_gf 
  if (excision_enable==1) then
     write(*,*) 'excision_radius = ',excision_radius
     call find_excision_zone(ext,x,y,z,excision_zone_gf,excision_radius,Symmetry)
  end if

  position_x(1) = 0.d0
  position_y(1) = 0.d0
  position_z(1) = 0.d0

end subroutine mag_bondi_initialdata_part2
