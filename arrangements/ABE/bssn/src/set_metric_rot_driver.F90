#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine set_metric_rotation(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,global_ext
  real*8 			           :: dT,dX,dY,dZ
  real*8 			           :: n1,n2,n3,mf
  integer :: index,dummy,foundit_flag
  integer :: handle,ierr
  real*8  :: ONE,ZERO,F1o3,F1o12,TWO
  real*8  :: phi_rot_ang
  real*8  :: Psi4L
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  real*8  :: Mass_adm,Mass_1,Mass_2,eta,omega

  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(ONE = 1.D0, ZERO = 0.D0, F1o3 = 1.D0/3.D0, F1o12 = 1.D0/12.D0, TWO = 2.D0)

  INTEGER :: i,j,k
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! NOTE THAT THIS SUBROUTINE IS JUST LIKE THE ONE IN THE BBH_COOKPFEIFFER_NEW THORN.  MAKE SURE THAT ANY CHANGES ARE MADE IN BOTH PLACES.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)


  ext = cctk_lsh
  
!  if (set_omega_test .eq. 1) then
!     phi_rot_ang=omega_test*cctk_time
!  else
     phi_rot_ang=binary_orb_freq*cctk_time
!  endif
  
  call set_metric_rot_lowlevel(cctkGH,phi_rot_ang,xbh1_initial,xbh2_initial, &
       radmin_rot_interp_1,radmax_rot_interp_1,radmin_rot_interp_2, &
       radmax_rot_interp_2,radmin_rot_interp_3,radmax_rot_interp_3,&
       Nlograd_rot_interp,Ntheta_rot_interp,Nphi_rot_interp,&
       K_rr_rot1,K_rth_rot1,K_rp_rot1,K_thth_rot1,K_thp_rot1,K_pp_rot1, &
       shiftr_rot1,shiftth_rot1,shiftp_rot1,phi_rot1,lapm1_rot1, &
       K_rr_rot2,K_rth_rot2,K_rp_rot2,K_thth_rot2,K_thp_rot2,K_pp_rot2, &
       shiftr_rot2,shiftth_rot2,shiftp_rot2,phi_rot2,lapm1_rot2, &
       K_rr_rot3,K_rth_rot3,K_rp_rot3,K_thth_rot3,K_thp_rot3,K_pp_rot3, &
       shiftr_rot3,shiftth_rot3,shiftp_rot3,phi_rot3,lapm1_rot3, &
       kxx,kxy,kxz,kyy,kyz,kzz,shiftx,shifty,shiftz,phi,lapm1, &
       X,Y,Z,X(1,1,1),Y(1,1,1),Z(1,1,1),dX,dY,dZ,cctk_lsh)
  
  !!  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')
  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')


  psi = exp(phi)

!! call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_AH')
  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_AH')

  ! Okay, now phi and psi are correct everywhere, even ghostzones!

!we don't need this option when we are rotating metric
!  if(use_cookpfeiffer_lapseshift==0) then
!     lapm1 = psi**(-2) - 1.D0
!     shiftx = 0.D0
!     shifty = 0.D0
!     shiftz = 0.D0
!  end if

  ! Next, set the conformal metric and its inverse
  !$omp parallel do private (Psi4L)
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
           Psi4L=Psi(i,j,k)*Psi(i,j,k)*Psi(i,j,k)*Psi(i,j,k)
           
           gxx(i,j,k) = 1.D0
           gxy(i,j,k) = 0.D0
           gxz(i,j,k) = 0.D0
           gyy(i,j,k) = 1.D0
           gyz(i,j,k) = 0.D0
           gzz(i,j,k) = 1.D0

           gupxx(i,j,k) = 1.D0
           gupxy(i,j,k) = 0.D0
           gupxz(i,j,k) = 0.D0
           gupyy(i,j,k) = 1.D0
           gupyz(i,j,k) = 0.D0
           gupzz(i,j,k) = 1.D0

           trK(i,j,k) = 0.D0

           ! Now compute Kzz from trK = 0 and all the other Kij's:
           ! This line is valid since gupij = deltaij
           Kzz(i,j,k) = - (Kxx(i,j,k) + Kyy(i,j,k))

           ! Old code snippet:
           !  trK =  gupxx * Kxx + gupyy * Kyy + gupzz * Kzz +  &
           !       TWO * ( gupxy * Kxy + gupxz * Kxz + gupyz * Kyz )

           !Compute Atilde_ij's.  Use fact that Kij is already traceless
           Axx(i,j,k) = Kxx(i,j,k) / Psi4L
           Axy(i,j,k) = Kxy(i,j,k) / Psi4L 
           Axz(i,j,k) = Kxz(i,j,k) / Psi4L 
           Ayy(i,j,k) = Kyy(i,j,k) / Psi4L
           Ayz(i,j,k) = Kyz(i,j,k) / Psi4L 
           Azz(i,j,k) = Kzz(i,j,k) / Psi4L 

           shiftxt(i,j,k) = 0.D0
           shiftyt(i,j,k) = 0.D0
           shiftzt(i,j,k) = 0.D0
           lapset(i,j,k) = 0.D0

           Gammax(i,j,k) = 0.D0
           Gammay(i,j,k) = 0.D0
           Gammaz(i,j,k) = 0.D0

           ! Set matter variables to zero!
           rho(i,j,k) = ZERO
           S(i,j,k) = ZERO
           Sx(i,j,k) = ZERO
           Sy(i,j,k) = ZERO
           Sz(i,j,k) = ZERO
           Sxx(i,j,k) = ZERO
           Sxy(i,j,k) = ZERO
           Sxz(i,j,k) = ZERO
           Syy(i,j,k) = ZERO
           Syz(i,j,k) = ZERO
           Szz(i,j,k) = ZERO

           !======================================
           ! Set everything else to Zero! 
           !======================================
           Gammaxxx(i,j,k) = ZERO
           Gammaxxy(i,j,k) = ZERO
           Gammaxxz(i,j,k) = ZERO
           Gammaxyy(i,j,k) = ZERO
           Gammaxyz(i,j,k) = ZERO
           Gammaxzz(i,j,k) = ZERO
           Gammayxx(i,j,k) = ZERO
           Gammayxy(i,j,k) = ZERO
           Gammayxz(i,j,k) = ZERO
           Gammayyy(i,j,k) = ZERO
           Gammayyz(i,j,k) = ZERO
           Gammayzz(i,j,k) = ZERO
           Gammazxx(i,j,k) = ZERO
           Gammazxy(i,j,k) = ZERO
           Gammazxz(i,j,k) = ZERO
           Gammazyy(i,j,k) = ZERO
           Gammazyz(i,j,k) = ZERO
           Gammazzz(i,j,k) = ZERO
        end do
     end do
  end do
  !$omp end parallel do

  
 
  ! Set excision_zone_gf to avoid valgrind memory errors
  excision_zone_gf = 0

  shiftxt = 0.D0
  shiftyt = 0.D0
  shiftzt = 0.D0
  lapset = 0.D0

  Gammax = 0.D0
  Gammay = 0.D0
  Gammaz = 0.D0

  !Look, Gamma^i was just set to zero,
  ! so next lines not strictly needed for brandt-brugmann ID, 
  ! but doesn't hurt in case this code is ever recycled!
  n1 = 0
  n2 = 0
  n3 = 0
  mf = 1
  !use gxx_p as temporary storage here:
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammax,gxx_p,n1,n2,n3,mf,Symmetry)
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammay,gxx_p,n1,n2,n3,mf,Symmetry)
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammaz,gxx_p,n1,n2,n3,mf,Symmetry)
  gxx_p=gxx

  !Next line NEEDED!
!!  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')



  if(fisheye_enable .eq. 1) then
     call trans_phys_fish_tensor(cctkGH,cctk_lsh,Symmetry, &
          x, y, z,  &
          PhysicalRadius,RadiusDerivative, &
          gxx,gxy,gxz,gyy,gyz,gzz)
     call trans_phys_fish_tensor_inv(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)
     call trans_phys_fish_tensor(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative, &
          Axx,Axy,Axz,Ayy,Ayz,Azz)
     call trans_phys_fish_phi(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative,phi)
     call trans_phys_fish_gamt_flat(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative,RadiusDerivative2, &
          Gammax,Gammay,Gammaz)
  end if
  
  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')
  call CartSymGN(dummy,cctkGH,'shift::shift_vars')
  
  if(Symmetry==AXISYM) then
     write(*,*) "I PITY DA FOO TRIES TO USE AXISYMMETRY FOR ORBITING BLACK HOLES!"
     stop
  end if
  
  !Following lines are (empirically speaking) ESSENTIAL for excision runs:
  call Derivs(ext,X,Y,Z,dX,dY,dZ,phi,phix,phiy,phiz,Symmetry) 
  call Derivs(ext,X,Y,Z,dX,dY,dZ,lapm1,lapsex,lapsey,lapsez,Symmetry) 

  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,phi,phix,phiy,phiz)
  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,lapm1,lapsex,lapsey,lapsez)
  
  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
!!  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars') 


  
end subroutine set_metric_rotation
