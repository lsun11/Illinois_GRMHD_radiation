#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine driver_mag_bondi_post_regrid(CCTK_ARGUMENTS)
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
  if (puncture_id==0) then 
     ext = cctk_lsh

     dT = CCTK_DELTA_TIME
     dX = CCTK_DELTA_SPACE(1)
     dY = CCTK_DELTA_SPACE(2)
     dZ = CCTK_DELTA_SPACE(3)

     !|~~~~~> Set up *untilded* metric, extrinsic curvature, lapse,shifts
     call ks_initial_metric_shift_origin_bl(ext,X,Y,Z, &
             gxx,gxy,gxz,gyy,gyz,gzz,trK,Axx,Axy,Axz,Ayy,Ayz,Azz, &
             lapm1,shiftx,shifty,shiftz,sam,r0,Symmetry)

     ! Sync before you call setgamma (which takes derivatives)
     !!call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
     !!call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')
     !!call CartSymGN(dummy,cctkGH,'shift::shift_vars')
     !!call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')

     !!if(Symmetry .eq. AXISYM .or. Symmetry .eq. AXISYM_FULL) then 
     !!   call CCTK_VarIndex(index,'BSSN::Axx')
     !!   call BndCartoon2DVI(dummy, cctkGH, 2, index)
     !!end if

     !|~~~~~> convert to tilded metric and initialize \phi and \Gamma^a
     !Convert to tilde metric and invert tilde metric:
     !  call convert_cpp(cctkGH,cctk_lsh,phi, &
     !       gxx,gxy,gxz,gyy,gyz,gzz, &
     !       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)
     call convert(ext,phi,gxx,gxy,gxz,gyy,gyz,gzz, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
          psi,detmin_l,detmax_l)


     !!call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
     !!call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')

     !Get Gamma^i:
     !  call setgamma_cpp(cctkGH,cctk_lsh, &
     !       dx,dy,dz, &
     !       gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
     !       Gammax, Gammay, Gammaz)
     !!call setgamma(ext,X,Y,Z, &
     !!     phi,gxx,gxy,gxz,gyy,gyz,gzz, &
     !!     gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
     !!     Gammax,Gammay,Gammaz, &
     !!     Symmetry)

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
     !$omp parallel do
     do k=1,ext(3)
        do j=1,ext(2)
           do i=1,ext(1)
              psi(i,j,k) = psi(i,j,k)**(1.d0/12.d0)
              Axx(i,j,k) = Axx(i,j,k)/psi(i,j,k)**4
              Axy(i,j,k) = Axy(i,j,k)/psi(i,j,k)**4
              Axz(i,j,k) = Axz(i,j,k)/psi(i,j,k)**4
              Ayy(i,j,k) = Ayy(i,j,k)/psi(i,j,k)**4
              Ayz(i,j,k) = Ayz(i,j,k)/psi(i,j,k)**4
              Azz(i,j,k) = Azz(i,j,k)/psi(i,j,k)**4
              rho(i,j,k) = 0.d0
              S(i,j,k) = 0.d0
              Sx(i,j,k) = 0.d0
              Sy(i,j,k) = 0.d0
              Sz(i,j,k) = 0.d0
              Sxx(i,j,k) = 0.d0
              Sxy(i,j,k) = 0.d0
              Sxz(i,j,k) = 0.d0
              Syy(i,j,k) = 0.d0
              Syz(i,j,k) = 0.d0
              Szz(i,j,k) = 0.d0
           end do
        end do
     end do
     !$omp end parallel do
  end if

end subroutine driver_mag_bondi_post_regrid
