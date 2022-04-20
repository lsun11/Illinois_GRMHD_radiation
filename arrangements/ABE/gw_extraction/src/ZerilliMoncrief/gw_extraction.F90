#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!===========================================================
!  Compute the total GW luminosity through a shell
!   by decomposing gij into its L=2,|M|=2 
!   moments, and getting L from them.
!===========================================================
subroutine gw_wave_flux(cctkGH,GW_measurement_radius,Mass, &
     dT,momentsr,momentsi, &
     momentsr_old,momentsi_old,odd_momentsr,odd_momentsi, & 
     odd_momentsr_old,odd_momentsi_old,int_momentsr,int_momentsi, & 
     hplus,hcross,th,ph,Zmin,dZ)
  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_POINTER :: cctkgh
  real*8                                   :: GW_measurement_radius,Mass
  real*8                                   :: dT,Zmin,dZ
  real*8,dimension(nmodes_ZM)                 :: momentsr,momentsi
  real*8,dimension(nmodes_ZM)                 :: momentsr_old,momentsi_old
  real*8,dimension(nmodes_ZM)                 :: odd_momentsr,odd_momentsi
  real*8,dimension(nmodes_ZM)                 :: odd_momentsr_old,odd_momentsi_old
  real*8                                   :: PI,PhysR,dRdr,th,ph
  real*8                                   :: symfactor,phiangle,costheta,sintheta,dphi,dcostheta
  real*8                                   :: r1,r2,rs_1,rs_2,grr_1,grr_2
  real*8                                   :: sint
  real*8,allocatable,dimension(:)          :: phiint,gintxx,gintxy,gintxz,gintyy,gintyz,gintzz
  real*8,allocatable,dimension(:)          :: grr,gtt,gpp,gtp,grt,grp
  real*8,allocatable,dimension(:,:)        :: pointcoords1,pointcoords2
  real*8,dimension(nmodes_ZM)                 :: H2r_1,h1r_1,Gr_1,Kr_1
  real*8,dimension(nmodes_ZM)                 :: H2i_1,h1i_1,Gi_1,Ki_1
  real*8,dimension(nmodes_ZM)                 :: H2r_2,h1r_2,Gr_2,Kr_2
  real*8,dimension(nmodes_ZM)                 :: H2i_2,h1i_2,Gi_2,Ki_2
  real*8,dimension(nmodes_ZM)                 :: k1r,k2r,k1i,k2i,Rr,Ri
  real*8,dimension(nmodes_ZM)                 :: dtmomentsr,dtmomentsi
  real*8,dimension(nmodes_ZM)                 :: int_momentsr,int_momentsi
  real*8,dimension(nmodes_ZM)                 :: Cr_1,Cr_2,Ci_1,Ci_2
  real*8,dimension(nmodes_ZM)                 :: Dr_1,Dr_2,Di_1,Di_2
  real*8                                   :: rs,S,A2,S_1,S_2,dr,Lambda,fac
  real*8                                   :: hplus, hcross
  real*8                                   :: Wlmr,Wlmi,Xlmr,Xlmi
  real*8, parameter			   :: sqrt2 = 1.41421356237309505d0
  integer                                  :: i,j,k,l,m,n
  integer,dimension(nmodes_ZM,2)              :: mode_array,odd_mode_array
  integer                                  :: interpolation_order,vindex
  integer                                  :: N_theta,N_phi,n_tot
  integer                                  :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  r1 = GW_measurement_radius - 0.02D0
  r2 = GW_measurement_radius + 0.02D0

  PI = 3.14159265358979323846D0

  sint = sin(th)

  dphi = 2.D0 * PI / numphi
  dcostheta = 2.D0 / numtheta
  N_theta = numtheta
  N_phi = numphi

  interpolation_order = 2

  if(Symmetry==AXISYM) dphi = 4.d0*PI

  if (Symmetry==EQUATORIAL) then
     symfactor = 2.D0
     N_theta = numtheta/2
  else if (Symmetry==NO_SYMM) then
     symfactor = 1.D0
  else if (Symmetry==PI_SYMM) then
     symfactor = 4.D0
     N_theta = numtheta/2
     N_phi   = numphi  /2    
  else if (Symmetry==AXISYM) then
     symfactor = 1.D0
     N_theta = numtheta/2
     N_phi   = 1.D0
     if (Zmin+dZ .lt. 0.d0) then
        dphi = 2.D0*PI
     else
        dphi = 4.D0*PI
     end if
  else if (Symmetry==OCTANT) then
     symfactor = 8.D0
     N_theta = numtheta/2
     N_phi   = numphi  /4
  end if

  n_tot = N_theta*N_phi

  allocate(pointcoords1(1:n_tot,1:3))
  allocate(pointcoords2(1:n_tot,1:3))

  allocate(phiint(1:n_tot))
  allocate(gintxx(1:n_tot))
  allocate(gintxy(1:n_tot))
  allocate(gintxz(1:n_tot))
  allocate(gintyy(1:n_tot))
  allocate(gintyz(1:n_tot))
  allocate(gintzz(1:n_tot))

  allocate(grr(1:n_tot))
  allocate(gtt(1:n_tot))
  allocate(gpp(1:n_tot))
  allocate(gtp(1:n_tot))
  allocate(grt(1:n_tot))
  allocate(grp(1:n_tot))

  n = 1
  do i=1,N_theta
     costheta = 1.D0 - (i - 0.5D0)*dcostheta
     sintheta = sqrt(1.D0 - costheta*costheta)
     do j=1,N_phi
        phiangle = (j - 0.5D0)*dphi
        if(N_phi==1) phiangle = 0.D0
        pointcoords1(n,1) = r1*sintheta*cos(phiangle)
        pointcoords1(n,2) = r1*sintheta*sin(phiangle)
        pointcoords1(n,3) = r1*costheta

        pointcoords2(n,1) = pointcoords1(n,1)*r2/r1
        pointcoords2(n,2) = pointcoords1(n,2)*r2/r1
        pointcoords2(n,3) = pointcoords1(n,3)*r2/r1
        n = n + 1
     end do
  end do

  !FIXME: This doesn't seem correct in general....
  if (Symmetry==AXISYM) then
     m=2
     if (Zmin+dZ .lt. 0.d0) m=1
     do i=1,nmodes_ZM 
	mode_array(i,1) = m*i
	mode_array(i,2) = 0
	odd_mode_array(i,1) = m*i+1
	odd_mode_array(i,2) = 0
     end do
  else if (Symmetry==EQUATORIAL) then
     i = 1
     l = 2
     do 
        m = l
        do 
           mode_array(i,1) = l
           mode_array(i,2) = m 
	   odd_mode_array(i,1) = l
	   odd_mode_array(i,2) = m-1
	   m = m - 2
	   i = i + 1
           if (m .lt. 0 .or. i .gt. nmodes_ZM) exit
        end do
        if (i .gt. nmodes_ZM) exit
        l = l+1
     end do
  else if (Symmetry==NO_SYMM) then 
     i = 1
     l = 2
     do 
        m = l
        do
           mode_array(i,1) = l
           mode_array(i,2) = m
           odd_mode_array(i,1) = l
           odd_mode_array(i,2) = m
           m = m - 1
           i = i + 1
           if (m .lt. 0 .or. i .gt. nmodes_ZM) exit
        end do
        if (i .gt. nmodes_ZM) exit
        l = l+1
     end do
  else
     write(*,*) 'Symmetry not supported in gw_extraction' 
  end if

  !==================================
  ! First compute variables on r1
  !==================================  
  if (fisheye_enable==1) then
     call CCTK_VarIndex(vindex,"fisheye::PhysicalRadius")
     call interp_driver_carp(cctkGH,n_tot,pointcoords1,vindex,phiint)
     PhysR = phiint(1)
     call CCTK_VarIndex(vindex,"fisheye::RadiusDerivative")
     call interp_driver_carp(cctkGH,n_tot,pointcoords1,vindex,phiint)
     dRdr = phiint(1)
  else
     PhysR = r1
     dRdr = 1.D0
  end if
  call CCTK_VarIndex(vindex,"bssn::phi")
  call interp_driver_carp(cctkGH,n_tot,pointcoords1,vindex,phiint)
  phiint = exp(4.D0*phiint)*(PhysR/r1)**(-4.d0/3.d0) * dRdr**(-2.d0/3.d0)
  call CCTK_VarIndex(vindex,"bssn::gxx")
  call interp_driver_carp(cctkGH,n_tot,pointcoords1,vindex,gintxx)
  call CCTK_VarIndex(vindex,"bssn::gxy")
  call interp_driver_carp(cctkGH,n_tot,pointcoords1,vindex,gintxy)
  call CCTK_VarIndex(vindex,"bssn::gxz")
  call interp_driver_carp(cctkGH,n_tot,pointcoords1,vindex,gintxz)
  call CCTK_VarIndex(vindex,"bssn::gyy")
  call interp_driver_carp(cctkGH,n_tot,pointcoords1,vindex,gintyy)
  call CCTK_VarIndex(vindex,"bssn::gyz")
  call interp_driver_carp(cctkGH,n_tot,pointcoords1,vindex,gintyz)
  call CCTK_VarIndex(vindex,"bssn::gzz")
  call interp_driver_carp(cctkGH,n_tot,pointcoords1,vindex,gintzz)

  ! transform to spherical coordinates
  call cart_to_spher(n_tot,pointcoords1,r1,gintxx,gintxy,gintxz,gintyy,gintyz,gintzz,grr,gtt,gpp,gtp,grt,grp,PhysR,dRdr)
  call Compute_Areal_Radius(n_tot,pointcoords1,phiint,gtt,gpp,rs_1)
  call Compute_grr_Average(n_tot,phiint,grr,grr_1)

  do k=1,nmodes_ZM
     l = mode_array(k,1)
     m = mode_array(k,2)
     ! Calculate H_2, h_1, G, and K from quadrature [Eqs. (46)--(49) in the GW 
     !						    extraction note]
     ! Here H2r is the real part of H_2 and H2i is the imaginary part of H_2, similarly
     ! for the other variables
     call Get_Efuncslm(n_tot,nmodes_ZM,dphi,symfactor,dcostheta, & 
                       H2r_1(k),h1r_1(k),Gr_1(k),Kr_1(k), & 
                       H2i_1(k),h1i_1(k),Gi_1(k),Ki_1(k), & 
                       r1,l,m,rs_1,grr_1,pointcoords1, & 
                       phiint,grr,gtt,gpp,gtp,grt,grp,  & 
                       Mass)
     ! Calculate C and D [Eqs. (37) & (38) of the GW note]
     l = odd_mode_array(k,1)
     m = odd_mode_array(k,2)
     if (m .ge. 0) then
	call Get_Ofuncslm(n_tot,nmodes_ZM,dphi,symfactor,dcostheta, &
                         Cr_1(k),Ci_1(k),Dr_1(k),Di_1(k), &
                         r1,l,m,rs_1,pointcoords1, &
                         phiint,grr,gtt,gpp,gtp,grt,grp)
     else
	Cr_1(k) = 0.d0
	Ci_1(k) = 0.d0
	Dr_1(k) = 0.d0
	Di_1(k) = 0.d0
     end if
  end do

  !==================================
  ! Next compute variables on r2
  !==================================  
  if (fisheye_enable==1) then
     call CCTK_VarIndex(vindex,"fisheye::PhysicalRadius")
     call interp_driver_carp(cctkGH,n_tot,pointcoords2,vindex,phiint)
     PhysR = phiint(1)
     call CCTK_VarIndex(vindex,"fisheye::RadiusDerivative")
     call interp_driver_carp(cctkGH,n_tot,pointcoords2,vindex,phiint)
     dRdr = phiint(1)
  else
     PhysR = r2
     dRdr = 1.d0
  end if
  call CCTK_VarIndex(vindex,"bssn::phi")
  call interp_driver_carp(cctkGH,n_tot,pointcoords2,vindex,phiint)
  phiint = exp(4.D0*phiint)*(PhysR/r2)**(-4.d0/3.d0) * dRdr**(-2.d0/3.d0)
  call CCTK_VarIndex(vindex,"bssn::gxx")
  call interp_driver_carp(cctkGH,n_tot,pointcoords2,vindex,gintxx)
  call CCTK_VarIndex(vindex,"bssn::gxy")
  call interp_driver_carp(cctkGH,n_tot,pointcoords2,vindex,gintxy)
  call CCTK_VarIndex(vindex,"bssn::gxz")
  call interp_driver_carp(cctkGH,n_tot,pointcoords2,vindex,gintxz)
  call CCTK_VarIndex(vindex,"bssn::gyy")
  call interp_driver_carp(cctkGH,n_tot,pointcoords2,vindex,gintyy)
  call CCTK_VarIndex(vindex,"bssn::gyz")
  call interp_driver_carp(cctkGH,n_tot,pointcoords2,vindex,gintyz)
  call CCTK_VarIndex(vindex,"bssn::gzz")
  call interp_driver_carp(cctkGH,n_tot,pointcoords2,vindex,gintzz)

  ! transform to spherical coordinates
  call cart_to_spher(n_tot,pointcoords2,r2,gintxx,gintxy,gintxz,gintyy,gintyz,gintzz,grr,gtt,gpp,gtp,grt,grp,PhysR,dRdr)

  call Compute_Areal_Radius(n_tot,pointcoords2,phiint,gtt,gpp,rs_2)
  call Compute_grr_Average(n_tot,phiint,grr,grr_2)

  do k=1,nmodes_ZM
     l = mode_array(k,1)
     m = mode_array(k,2)
     ! Calculate H2, h1, G, and K from quadrature [Eqs. (46)--(49) in the GW
     !                                              extraction note]
     call Get_Efuncslm(n_tot,nmodes_ZM,dphi,symfactor,dcostheta, &
                       H2r_2(k),h1r_2(k),Gr_2(k),Kr_2(k), &
                       H2i_2(k),h1i_2(k),Gi_2(k),Ki_2(k), &
                       r2,l,m,rs_2,grr_2,pointcoords2, &
                       phiint,grr,gtt,gpp,gtp,grt,grp,  &
                       Mass)
     ! Calculate C and D [Eqs. (37) & (38) of the GW note]
     l = odd_mode_array(k,1)
     m = odd_mode_array(k,2)
     if (m .ge. 0) then
        call Get_Ofuncslm(n_tot,nmodes_ZM,dphi,symfactor,dcostheta, &
                         Cr_2(k),Ci_2(k),Dr_2(k),Di_2(k), &
                         r2,l,m,rs_2,pointcoords2, &
                         phiint,grr,gtt,gpp,gtp,grt,grp)
     else
        Cr_2(k) = 0.d0
        Ci_2(k) = 0.d0
        Dr_2(k) = 0.d0
        Di_2(k) = 0.d0
     end if
  end do

  !=======================================
  ! Finally, compute R_lm and Q_lm [Eqs. (54) && (50) of the GW note]
  !=======================================
  rs = 0.5d0*(rs_1 + rs_2)
  S = (1.d0 - 2.d0*Mass/rs)
  A2 = 1.d0/S
  S_1 = (1.d0 - 2.d0*Mass/rs_1)
  S_2 = (1.d0 - 2.d0*Mass/rs_2)
  dr = rs_2 - rs_1

  hplus = 0.D0
  hcross = 0.D0
  do k=1,nmodes_ZM
     l = mode_array(k,1)
     m = mode_array(k,2)
     !
     ! First calculate k1 and k2 from Eqs. (55) & (56) of the GW note.
     !
     k1r(k) = 0.5D0*(Kr_1(k)+Kr_2(k)) + l*(l+1.D0)*0.5D0*(Gr_1(k)+Gr_2(k)) &
          + 2.D0*rs*S*(Gr_2(k)-Gr_1(k))/dr - S*(h1r_1(k)+h1r_2(k))/rs
     k1i(k) = 0.5D0*(Ki_1(k)+Ki_2(k)) + l*(l+1.D0)*0.5D0*(Gi_1(k)+Gi_2(k)) &
          + 2.D0*rs*S*(Gi_2(k)-Gi_1(k))/dr - S*(h1i_1(k)+h1i_2(k))/rs
     k2r(k) = 0.25d0*(H2r_1(k)+H2r_2(k))/S &
          - (0.5D0/dr)*(rs_2/sqrt(S_2)*(Kr_2(k)+l*(l+1.D0)*Gr_2(k)) &
          - rs_1/sqrt(S_1)*(Kr_1(k)+l*(l+1.D0)*Gr_1(k)))
     k2i(k) = 0.25d0*(H2i_1(k)+H2i_2(k))/S &
          - (0.5D0/dr)*(rs_2/sqrt(S_2)*(Ki_2(k)+l*(l+1.D0)*Gi_2(k)) &
          - rs_1/sqrt(S_1)*(Ki_1(k)+l*(l+1.D0)*Gi_1(k)))
     Lambda = ( (l-1.d0)*(l+2.D0) + 6.D0*Mass/rs ) * l*(l+1.d0)
     ! Now compute R using Eq. (54)
     Rr(k) = (4.D0*S*S*k2r(k)+l*(l+1.D0)*k1r(k)) * rs / Lambda
     Ri(k) = (4.D0*S*S*k2i(k)+l*(l+1.D0)*k1i(k)) * rs / Lambda
     momentsr(k) = Rr(k)
     momentsi(k) = Ri(k)

     ! Now add h_+ and h_x from the even modes [Eqs. (67) & (68) of the GW note]
     ! Note that the factor of 2 for the m/=0 modes has already been taken into account 
     ! when computing the Wlm and Xlm in compute_Wlm_Xlm
     !
     if (m .ge. 0) then
        call compute_Wlm_Xlm(l,m,th,ph,Wlmr,Wlmi,Xlmr,Xlmi)
        hplus = hplus + (Rr(k)*Wlmr - Ri(k)*Wlmi)/rs
        hcross = hcross + (Rr(k)*Xlmr - Ri(k)*Xlmi)/sint/rs
     end if

   ! The following 4 lines are useless if we don't compute E_gw and J_gw.
     !dtmomentsr(i) = (momentsr(i) - momentsr_old(i))/dT
     !dtmomentsi(i) = (momentsi(i) - momentsi_old(i))/dT
     !momentsr_old(i) = momentsr(i)
     !momentsi_old(i) = momentsi(i)

     ! 
     ! Now compute the odd modes
     !
     l = odd_mode_array(k,1)
     m = odd_mode_array(k,2)
     !
     ! Calculate Q from Eq. (50).
     ! Note that here R is really the Q in the GW note. I recycle the R variables 
     ! here
     !
     Rr(k) = (0.5d0*(Cr_1(k)+Cr_2(k)) + rs*rs*(Dr_2(k)-Dr_1(k))/dr)*S/rs
     Ri(k) = (0.5d0*(Ci_1(k)+Ci_2(k)) + rs*rs*(Di_2(k)-Di_1(k))/dr)*S/rs
     odd_momentsr(k) = Rr(k)
     odd_momentsi(k) = Ri(k)
     !
     ! Add h_+ and h_x from the odd modes [Eqs. (65) & (66) of the GW note]
     ! Note that the factor of 2 for the m/=0 modes has already been taken into account
     ! when computing the Wlm and Xlm in compute_Wlm_Xlm
     !
     int_momentsr(k) = int_momentsr(k) + &
                       0.5d0*dT*(odd_momentsr(k) + odd_momentsr_old(k))
     int_momentsi(k) = int_momentsi(k) + &
                       0.5d0*dT*(odd_momentsi(k) + odd_momentsi_old(k))
     odd_momentsr_old(k) = odd_momentsr(k)
     odd_momentsi_old(k) = odd_momentsi(k)
     if (m .ge. 0) then
        call compute_Wlm_Xlm(l,m,th,ph,Wlmr,Wlmi,Xlmr,Xlmi)
        hplus = hplus - (int_momentsr(k)*Xlmr - int_momentsi(k)*Xlmi)/sint/rs
        hcross = hcross + (int_momentsr(k)*Wlmr - int_momentsi(k)*Wlmi)/rs
     end if

  end do

end subroutine gw_wave_flux

subroutine cart_to_spher(n_tot,pointcoords,ri,sxx,sxy,sxz,syy,syz,szz,srr,stt,spp,stp,srt,srp,PhysR,dRdr)
  implicit none

  real*8,dimension(3,3)                    :: L,g
  integer                                  :: n_tot
  real*8,dimension(n_tot,3)                :: pointcoords
  real*8,dimension(n_tot)                  :: sxx,sxy,sxz,syy,syz,szz
  real*8,dimension(n_tot)                  :: srr,stt,spp,stp,srt,srp
  integer                                  :: i,ll,m
  real*8                                   :: ri,cst,snt,csp,snp,PhysR,dRdr,fac

  do i=1,n_tot
     cst = pointcoords(i,3)/ri
     snt = sqrt(1.D0-cst*cst)
     csp = pointcoords(i,1)/(ri*snt)
     snp = pointcoords(i,2)/(ri*snt)
     !===========================================
     ! transformation matrix
     ! Notation:  L(cartesian index,polar index)
     !===========================================
     L(1,1) = pointcoords(i,1)/ri
     L(2,1) = pointcoords(i,2)/ri
     L(3,1) = pointcoords(i,3)/ri
     L(1,2) = ri*csp*cst
     L(2,2) = ri*snp*cst
     L(3,2) = -ri*snt
     L(1,3) = -pointcoords(i,2)
     L(2,3) = pointcoords(i,1)
     L(3,3) = 0.D0

     g(1,1) = sxx(i)
     g(1,2) = sxy(i)
     g(2,1) = sxy(i)
     g(1,3) = sxz(i)
     g(3,1) = sxz(i)
     g(2,2) = syy(i)
     g(2,3) = syz(i)
     g(3,2) = syz(i)
     g(3,3) = szz(i)

     srr(i) = 0.D0
     stt(i) = 0.D0
     spp(i) = 0.D0
     stp(i) = 0.D0
     srt(i) = 0.D0
     srp(i) = 0.D0

     do ll=1,3
        do m=1,3
           srr(i) = srr(i) + L(ll,1)*L(m,1)*g(ll,m)
           srt(i) = srt(i) + L(ll,1)*L(m,2)*g(ll,m)
           srp(i) = srp(i) + L(ll,1)*L(m,3)*g(ll,m)
           stt(i) = stt(i) + L(ll,2)*L(m,2)*g(ll,m)
           stp(i) = stp(i) + L(ll,2)*L(m,3)*g(ll,m)
           spp(i) = spp(i) + L(ll,3)*L(m,3)*g(ll,m)
        end do
     end do
  end do

  srr = srr * (PhysR/ri/dRdr)**(4.d0/3.d0)
  fac = dRdr**(2.d0/3.d0) * (PhysR/ri)**(4.d0/3.d0)
  stt = stt * fac
  spp = spp * fac
  stp = stp * fac
  fac = dRdr**(-1.d0/3.d0) * (PhysR/ri)**(4.d0/3.d0)
  srt = srt * fac
  srp = srp * fac
  
end subroutine cart_to_spher

! Compute the Moncrief functions (for "even" modes)
!
subroutine Get_Efuncslm(n_tot,nmodes_ZM,dphi,symfactor,dcostheta, & 
                         H2r_i,h1r_i,Gr_i,Kr_i,H2i_i,h1i_i,Gi_i,Ki_i, & 
                         ri,l,msign,rs,grr_avg,pointcoords, & 
                         Psi4,srr,stt,spp,stp,srt,srp, &
                         Mass)
  use gw_functions
  implicit none
  integer                                  :: n_tot,nmodes_ZM,l,msign
  real*8                                   :: grr_avg
  real*8                                   :: H2r_i,h1r_i,Gr_i,Kr_i
  real*8                                   :: H2i_i,h1i_i,Gi_i,Ki_i
  real*8                                   :: ri,rs,PI,dphi,symfactor,dcostheta
  real*8                                   :: fac1,fac2,fac1_0,fac2_0,Plm,Plm1,Pl1
  real*8                                   :: Ylmr,Ylmtr,Ylmpr,Wlmr,Xlmr
  real*8                                   :: Ylmi,Ylmti,Ylmpi,Wlmi,Xlmi
  real*8                                   :: cost,sint,cott,p,cosmp,sinmp
  real*8                                   :: Mass,facM, pfac
  real*8,dimension(n_tot,3)                :: pointcoords
  real*8,dimension(n_tot)                  :: Psi4,srr,stt,spp,stp,srt,srp

  integer m,i

  PI = 3.14159265358979323846D0

  facM = sqrt(1.d0-2.d0*Mass/rs)

  H2r_i = 0.D0
  h1r_i = 0.D0
  Gr_i = 0.D0
  Kr_i = 0.D0
  H2i_i = 0.D0
  h1i_i = 0.D0
  Gi_i = 0.D0
  Ki_i = 0.D0
  m = abs(msign)
  !fac2 = sqrt( (2.D0*l+1.D0)*factorial(l-m) / ( 4.D0*PI*factorial(l+m) ) )
  fac1_0 = sqrt( (2.D0*l+1.D0)/(4.d0*PI) )
  fac2 = (2.D0*l+1.D0)/(4.d0*PI)
  do i=l-m+1,l+m
     fac2 = fac2/i
  end do
  fac2 = sqrt(fac2)
  fac1 = fac2 * (l+m)*(l-m+1.d0)
  do i=1,n_tot
     cost = pointcoords(i,3)/ri
     sint = sqrt(1.D0 - cost*cost)
     cott = cost/sint
     ! Note: The function atan2(y,x) returns a real number in (-pi,pi]. 
     !       It already takes care of the quadrant the angle is in.
     p = atan2(pointcoords(i,2),pointcoords(i,1))
     cosmp = cos( msign * p )
     sinmp = sin( msign * p )
     Plm = plgndr(l,m,cost)
     if(m==0) then
        Pl1 = plgndr(l,1,cost)
        Ylmr = fac2*Plm
        Ylmtr = fac1_0 * Pl1
        Ylmpr = 0.D0
        Wlmr = -l*(l+1.D0)*Ylmr - 2.D0*cott*Ylmtr
        Xlmr = 0.D0
        Ylmi = 0.d0
	Ylmti = 0.d0
	Ylmpi = 0.d0
	Wlmi = 0.d0
	Xlmi = 0.d0
     else
        Plm1 = plgndr(l,m-1,cost)
        Ylmr = fac2*Plm     
        Ylmtr = -fac1*Plm1 - m*cott*fac2*Plm  
        if (msign .gt. 0) then 
	   pfac = 1.d0
	else
	   pfac = (-1.d0)**(mod(m,2))
	end if
        Ylmr = Ylmr*pfac      ! Y_{lm} e^{-i m phi}
        Ylmtr = Ylmtr *pfac   ! \partial_{theta} Y_{lm} e^{-i m phi}
        Ylmpr = msign*Ylmr    ! \partial_{phi} Y_{lm} e^{-i m phi} / i
        Wlmr = -l*(l+1.d0)*Ylmr -2.d0*cott*Ylmtr + 2.d0*Ylmr*(m/sint)**2 ! W_{lm} e^{-i m phi}
    	Xlmr = 2.d0*(Ylmtr - cott*Ylmr)*msign ! X_{lm} e^{-i m phi} / i
	Ylmi = Ylmr*sinmp    ! Im(Y_{lm})
	Ylmti = Ylmtr*sinmp  ! Im(\partial_{theta} Y_{lm})
	Ylmpi = Ylmpr*cosmp  ! Im(\partial_{phi} Y_{lm})
	Wlmi = Wlmr * sinmp  ! Im(W_{lm})
	Xlmi = Xlmr*cosmp    ! Im(X_{lm})
	Ylmr = Ylmr*cosmp    ! Re(Y_{lm})
	Ylmtr = Ylmtr*cosmp  ! Re(\partial_{theta} Y_{lm})
	Ylmpr = -Ylmpr*sinmp ! Re(\partial_{phi} Y_{lm})
	Wlmr = Wlmr*cosmp    ! Re(W_{lm})
	Xlmr = -Xlmr*sinmp   ! Re(X_{lm})
     end if
     H2r_i = H2r_i + Ylmr * (Psi4(i) * srr(i) - grr_avg) 
     h1r_i = h1r_i + (1.D0/((l+1.d0)*l)) * Psi4(i) * &
          (Ylmtr * srt(i) + Ylmpr * srp(i)/(sint*sint))
     Gr_i = Gr_i + (0.5D0/((l-1.d0)*l*(l+1)*(l+2)))*(Psi4(i)/(rs*rs)) * &
          (Wlmr * (stt(i) - spp(i)/(sint*sint)) + 2.D0*Xlmr*stp(i)/(sint*sint))
     Kr_i = Kr_i + 0.5D0 * Ylmr * ( Psi4(i)/(rs*rs) * &
                  ( stt(i) + spp(i)/(sint*sint) ) - 2.d0 ) 
     H2i_i = H2i_i - Ylmi * (Psi4(i) * srr(i) - grr_avg)
     h1i_i = h1i_i - (1.D0/((l+1)*l)) * Psi4(i) * &
          (Ylmti * srt(i) + Ylmpi * srp(i)/(sint*sint))
     Gi_i = Gi_i - (0.5D0/((l-1)*l*(l+1)*(l+2)))*Psi4(i)/(rs*rs) * &
          (Wlmi * (stt(i) - spp(i)/(sint*sint)) + 2.D0*Xlmi*stp(i)/(sint*sint))
     Ki_i = Ki_i - 0.5D0 * Ylmi * ( Psi4(i) * ( stt(i) + spp(i)/(sint*sint) ) &
                     - 2.d0 ) / (rs*rs)
  end do
  H2r_i = H2r_i * dphi * dcostheta * symfactor * facM
  h1r_i = h1r_i * dphi * dcostheta * symfactor 
  Gr_i  = Gr_i  * dphi * dcostheta * symfactor 
  Kr_i  = Kr_i  * dphi * dcostheta * symfactor 
  H2i_i = H2i_i * dphi * dcostheta * symfactor * facM
  h1i_i = h1i_i * dphi * dcostheta * symfactor
  Gi_i  = Gi_i  * dphi * dcostheta * symfactor
  Ki_i  = Ki_i  * dphi * dcostheta * symfactor

end subroutine Get_Efuncslm

! Compute the Moncrief functions (for "odd" modes)
!
subroutine Get_Ofuncslm(n_tot,nmodes_ZM,dphi,symfactor,dcostheta, &
                         Cr_i,Ci_i,Dr_i,Di_i, &
                         ri,l,msign,rs,pointcoords, &
                         Psi4,srr,stt,spp,stp,srt,srp)
  use gw_functions
  implicit none
  integer                                  :: n_tot,nmodes_ZM,l,msign
  real*8                                   :: Cr_i,Ci_i,Dr_i,Di_i
  real*8                                   :: ri,rs,PI,dphi,symfactor,dcostheta
  real*8                                   :: fac1,fac2,fac1_0,fac2_0,fac,fad
  real*8                                   :: Plm,Plm1,Pl1, pfac
  real*8                                   :: Ylmr,Ylmtr,Ylmpr,Wlmr,Xlmr
  real*8                                   :: Ylmi,Ylmti,Ylmpi,Wlmi,Xlmi
  real*8                                   :: cost,sint,cott,p,cosmp,sinmp
  real*8,dimension(n_tot,3)                :: pointcoords
  real*8,dimension(n_tot)                  :: Psi4,srr,stt,spp,stp,srt,srp
  integer m,i

  PI = 3.14159265358979323846D0

  Cr_i = 0.d0
  Ci_i = 0.d0
  Dr_i = 0.d0
  Di_i = 0.d0
  m = abs(msign)
  !fac2 = sqrt( (2.D0*l+1.D0)*factorial(l-m) / ( 4.D0*PI*factorial(l+m) ) )
  fac1_0 = sqrt( (2.D0*l+1.D0)/(4.d0*PI) )
  fac2 = (2.D0*l+1.D0)/(4.d0*PI)
  do i=l-m+1,l+m
     fac2 = fac2/i
  end do
  fac2 = sqrt(fac2)
  fac1 = fac2 * (l+m)*(l-m+1.d0)
  do i=1,n_tot
     cost = pointcoords(i,3)/ri
     sint = sqrt(1.D0 - cost*cost)
     cott = cost/sint
     ! Note: The function atan2(y,x) returns a real number in (-pi,pi].
     !       It already takes care of the quadrant the angle is in.
     p = atan2(pointcoords(i,2),pointcoords(i,1))
     cosmp = cos( msign * p )
     sinmp = sin( msign * p )
     Plm = plgndr(l,m,cost)
     if(m==0) then
        Pl1 = plgndr(l,1,cost)
        Ylmr = fac2*Plm
        Ylmtr = fac1_0 * Pl1
        Ylmpr = 0.D0
        Wlmr = -l*(l+1.D0)*Ylmr - 2.D0*cott*Ylmtr
        Xlmr = 0.D0
        Ylmi = 0.d0
        Ylmti = 0.d0
        Ylmpi = 0.d0
        Wlmi = 0.d0
        Xlmi = 0.d0
     else
        Plm1 = plgndr(l,m-1,cost)
        Ylmr = fac2*Plm
        Ylmtr = -fac1*Plm1 - m*cott*fac2*Plm
	if (msign .gt. 0) then 
	   pfac = 1.d0
	else
	   pfac = (-1.d0)**(mod(m,2))
	end if
	Ylmr = pfac*Ylmr
	Ylmtr = pfac*Ylmtr
        Ylmpr = msign*Ylmr
        Wlmr = -l*(l+1.d0)*Ylmr -2.d0*cott*Ylmtr + 2.d0*Ylmr*(m/sint)**2
        Xlmr = 2.d0*(Ylmtr - cott*Ylmr)*msign
        Ylmi = Ylmr*sinmp
        Ylmti = Ylmtr*sinmp
        Ylmpi = Ylmpr*cosmp
        Wlmi = Wlmr * sinmp
        Xlmi = Xlmr*cosmp
        Ylmr = Ylmr*cosmp
        Ylmtr = Ylmtr*cosmp
        Ylmpr = -Ylmpr*sinmp
        Wlmr = Wlmr*cosmp
        Xlmr = -Xlmr*sinmp
     end if

     Cr_i = Cr_i + Psi4(i)/sint*(-Ylmpr*srt(i) + Ylmtr*srp(i))
     Ci_i = Ci_i - Psi4(i)/sint*(-Ylmpi*srt(i) + Ylmti*srp(i))
     Dr_i = Dr_i + Psi4(i) * ( (stt(i) - spp(i)/sint**2)/sint*Xlmr*0.5d0 & 
				- stp(i)*Wlmr/sint )
     Di_i = Di_i - Psi4(i) * ( (stt(i) - spp(i)/sint**2)/sint*Xlmi*0.5d0 & 
                                - stp(i)*Wlmi/sint )
  end do

  fac = dphi*dcostheta*symfactor/(l*(l+1.d0))
  fad = fac/( (l-1.d0)*(l+2.d0)*rs**2 )
  Cr_i = Cr_i*fac
  Ci_i = Ci_i*fac
  Dr_i = Dr_i*fad
  Di_i = Di_i*fad

end subroutine Get_Ofuncslm

subroutine compute_Wlm_Xlm(l,m,th,ph,Wlmr,Wlmi,Xlmr,Xlmi)
  use gw_functions
  implicit none
  integer :: l,m,abm,i
  real*8 :: th,ph
  real*8 :: cost,sint,cott,cosmp,sinmp,Plm,fac1,fac1_0,fac2,Pl1
  real*8 :: Plm1,Ylmr,Ylmi,Ylmtr,Ylmti,Wlmr,Wlmi,Xlmr,Xlmi
  real*8, parameter :: PI = 3.14159265358979323846d0
!
  abm = abs(m)
  fac1_0 = sqrt( (2.D0*l+1.D0)/(4.d0*PI) )
  !fac2 = sqrt( (2.D0*l+1.D0)*factorial(l-abm) / ( 4.D0*PI*factorial(l+abm) ) )
  fac2 = (2.D0*l+1.D0)/(4.D0*PI)
  do i=l-abm+1,l+abm
     fac2 = fac2/i
  end do
  fac2 = sqrt(fac2)
  fac1 = fac2 * (l+m)*(l-m+1.d0)
  cost = cos(th)
  sint = sin(th)
  cott = cost/sint
  cosmp = cos( abm * ph )
  sinmp = sin( abm * ph )
  Plm = plgndr(l,abm,cost)
  if(m==0) then
     Pl1 = plgndr(l,1,cost)
     Ylmr = fac2*Plm
     Ylmtr = fac1_0 * Pl1
     Wlmr = -l*(l+1.D0)*Ylmr - 2.D0*cott*Ylmtr
     Xlmr = 0.D0
     Wlmi = 0.d0
     Xlmi = 0.d0
  else
     Plm1 = plgndr(l,abm-1,cost)
     Ylmr = fac2*Plm
     Ylmtr = -fac1*Plm1 - abm*cott*fac2*Plm
     ! Note the factor of 2 is to account for the +m and -m modes
     Wlmr = 2.d0*( -l*(l+1.d0)*Ylmr -2.d0*cott*Ylmtr + 2.d0*Ylmr*(m/sint)**2 )
     Xlmr = 4.d0*(Ylmtr - cott*Ylmr)*abm
     Wlmi = Wlmr * sinmp
     Xlmi = Xlmr*cosmp
     Wlmr = Wlmr*cosmp
     Xlmr = -Xlmr*sinmp
  end if
end subroutine compute_Wlm_Xlm

subroutine Compute_Areal_Radius(n_tot,pointcoords,Psi4,stt,spp,rtot)
  implicit none

  integer n_tot
  real*8,dimension(n_tot,3) :: pointcoords
  real*8,dimension(n_tot) :: Psi4,stt,spp
  real*8 rtot,sint

  integer i
  
  rtot = 0.D0
  do i=1,n_tot
     sint = sqrt( (pointcoords(i,1)*pointcoords(i,1) + pointcoords(i,2)*pointcoords(i,2)) &
          /(pointcoords(i,1)*pointcoords(i,1) + pointcoords(i,2)*pointcoords(i,2) + pointcoords(i,3)*pointcoords(i,3)) ) 
     rtot = rtot + Psi4(i)*(stt(i) + spp(i)/(sint*sint))/2.D0
  end do
  rtot = sqrt( rtot / n_tot )
end subroutine Compute_Areal_Radius

subroutine Compute_grr_Average(n_tot,Psi4,srr,Gtot)
  implicit none

  integer n_tot
  real*8,dimension(n_tot) :: Psi4,srr
  real*8 Gtot,sint

  integer i  

  Gtot = 0.D0
  do i=1,n_tot
     Gtot = Gtot + Psi4(i)*srr(i)
  end do
  Gtot =  Gtot/dble(n_tot)
end subroutine Compute_grr_Average
