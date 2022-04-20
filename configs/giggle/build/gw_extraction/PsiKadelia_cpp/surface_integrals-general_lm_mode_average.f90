!---------------------------------------------------------------------
! This routine interpolates Psi_4 onto a spherical shell of radius r,
!   then it decomposes Psi_4 into spin weight -2 spherical harmonics
!---------------------------------------------------------------------
subroutine lm_modes(cctkGH,Symmetry,l,m,N_theta,N_phi,ntot, &
     dphi,dcostheta,radius_GW,Psi4reint,Psi4imint,Psi4relm,Psi4imlm)
  use psikad_lm_functions
  implicit none
  integer*8 :: cctkGH
  real*8                                   :: radius_GW
  integer                                  :: N_theta,N_phi,sym_factor,ntot
  real*8                                   :: dcostheta,dphi
  real*8                                   :: Psi4relm,Psi4imlm
  real*8, dimension(ntot)                  :: Psi4reint,Psi4imint
  real*8                                   :: phiangle,costheta,sintheta,PI,Y22nophi,sintwophi,costwophi
  real*8                                   :: Wlmr,Wlmi,Xlmr,Xlmi
  integer                                  :: i,j,n, odd,even,vindex
  integer                                  :: l,m,mmin,mmax
  real*8                                   :: l_normalization,tempXr,tempXi,tempWr,tempWi
  integer                                  :: Symmetry
  integer                                  :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  PI = 3.14159265358979323846D0
  if (Symmetry==OCTANT) then
     sym_factor=8
  else if (Symmetry==EQUATORIAL) then
     sym_factor=2
  else if (Symmetry==NO_SYMM) then
     sym_factor=1
  else if (Symmetry==PI_SYMM) then
     sym_factor=4
  else if (Symmetry==AXISYM) then
     sym_factor=2
  end if
  odd = mod(l+m,2)
  even = 1-odd
  Psi4relm = 0.D0
  Psi4imlm = 0.D0
  n = 1
  do i=1,N_theta
     costheta = 1.D0 - (i - 0.5D0)*dcostheta
     sintheta = sqrt(1.D0 - costheta*costheta)
     do j=1,N_phi
        phiangle = (j - 0.5D0)*dphi
        if(N_phi==1) phiangle = 0.D0
        call psikad_compute_Wlm_Xlm(l,m,costheta,sintheta,phiangle,Wlmr,Wlmi,Xlmr,Xlmi)
        Xlmr = Xlmr/sintheta
        Xlmi = Xlmi/sintheta
        l_normalization = sqrt( 1.d0/((l+2.d0)*(l+1.d0)*l*(l-1.d0)) )
        !==============================================================================
        ! IMPORTANT DERIVATION DONT DELETE:
        ! Notice that weve already divided Xlmr, Xlmi by sintheta!
        !    Integrand = Psi4 * Y^{star}lm(s=-2) 
        !              = Psi4 * lnorm * ((Wlmr + i*Wlmi) - i*(Xlmr + i*Xlmi))^{star}
        !              = Psi4 * lnorm * (Wlmr + i*Wlmi - i*Xlmr + Xlmi)^{star}
        !              = Psi4 * lnorm * (Wlmr - i*Wlmi + i*Xlmr + Xlmi)
        !              = Psi4 * lnorm * (Wlmr + Xlmi + i(Xlmr - Wlmi))
        !              = lnorm * (Psi4r + i*Psi4i) * (Wlmr + Xlmi + i(Xlmr - Wlmi))
        !              = lnorm * (Psi4r * (Wlmr + Xlmi) - Psi4i * (Xlmr - Wlmi) + 
        !                       i(Psi4i * (Wlmr + Xlmi) + Psi4r * (Xlmr - Wlmi)))
        !==============================================================================
        ! Okay, Xlmi and Xlmr are antisymmetric across the equatorial plane, and the Wlms are symmetric
        !  Also, the Psi4re is symmetric and Psi4im is antisymmetric.
        !  We must throw away all equatorially antisymmetric products manually, or well get a wrong answer.
        if (Symmetry .ne. NO_SYMM) then
           Psi4relm = Psi4relm + l_normalization * (Psi4reint(n) * (Wlmr*even + Xlmi*odd) &
                - Psi4imint(n) * (Xlmr*even - Wlmi*odd) )
           Psi4imlm = Psi4imlm + l_normalization * (Psi4imint(n) * (Xlmi*even + Wlmr*odd) &
                + Psi4reint(n) * (Xlmr*odd - Wlmi*even) )
        else
           Psi4relm = Psi4relm + l_normalization * (Psi4reint(n) * (Wlmr + Xlmi) &
                - Psi4imint(n) * (Xlmr - Wlmi) )
           Psi4imlm = Psi4imlm + l_normalization * (Psi4imint(n) * (Xlmi + Wlmr) &
                + Psi4reint(n) * (Xlmr - Wlmi) )
         end if
        !if(isnan(Psi4relm)) then
        !   write(*,*) OUCH. NaN found Re.,Psi4reint(n),Psi4imint(n),radius_GW,costheta,phiangle
        !   write(*,*) coords,radius_GW*sintheta*cos(phiangle),radius_GW*sintheta*sin(phiangle),radius_GW*costheta
        !   stop
        !end if
        !              call psikad_compute_Wlm_Xlm(l,m-2,costheta,sintheta,phiangle,tempWr,tempWi,tempXr,tempXi)
        !              tempXr = tempXr/sintheta
        !              tempXi = tempXi/sintheta
        !==============================================================================
        ! CHECK THAT INTEGRAL[Y^{star}lm(s=-2) Ylm(s=-2) dOmega] == 0 for m!=m:
        !    Integrand = Y^{star}lm(s=-2) Ylm(s=-2)
        !              = lnorm^2 * (tempWr + tempXi + i(tempXr - tempWi)) * (Wlmr + Xlmi - i(Xlmr - Wlmi))
        !              = lnorm^2  *[ (tempWr + tempXi)*(Wlmr + Xlmi) + (tempXr - tempWi)*(Xlmr - Wlmi)
        !                         *i((tempXr - tempWi)*(Wlmr + Xlmi) - (tempXr - tempWi)*(Wlmr + Xlmi))]
        !              Psi4resum = Psi4resum + l_normalization**2 * (tempWr*Wlmr + tempXi*Xlmi + tempXr*Xlmr + tempWi*Wlmi)
        !              Psi4imsum = Psi4imsum + l_normalization**2 * (tempXr*Wlmr - tempWi*Xlmi - tempXr*Wlmr + tempWi*Xlmi)
        !==============================================================================
        ! CHECK THAT INTEGRAL[Y^{star}lm(s=-2) Ylm(s=-2) dOmega] == 1:
        !    Integrand = Y^{star}lm(s=-2) Ylm(s=-2)
        !              = lnorm^2 * (Wlmr + Xlmi + i(Xlmr - Wlmi)) * (Wlmr + Xlmi - i(Xlmr - Wlmi))
        !              = lnorm^2 * [((Wlmr + Xlmi)*(Wlmr + Xlmi) + (Xlmr - Wlmi)*(Xlmr - Wlmi)) 
        !                         +i((Xlmr - Wlmi)*(Wlmr + Xlmi) - (Wlmr + Xlmi)*(Xlmr - Wlmi))]
        !              = lnorm^2 * ((Wlmr + Xlmi)**2 + (Xlmr - Wlmi)**2)
        !==============================================================================
        !              Psi4resum = Psi4resum + l_normalization**2 * ((Wlmr**2 + Xlmi**2) + (Xlmr**2 + Wlmi**2))
        !              Psi4imsum = Psi4imsum + 0.D0
        !              write(*,*) Wlmr**2,Xlmi**2,Xlmr**2,Wlmi**2
        !              Y22nophi = sqrt(5.D0/(64.D0*PI))*(1.D0+costheta)**2 
        !              Psi4resum = Psi4resum + Y22nophi*Y22nophi
        !              Psi4imsum = Psi4imsum + Y22nophi*Y22nophi
        !              write(*,*) HIIIII,i,j,Y22nophi*Y22nophi
        !write(*,*) hi.,Y22nophi*Y22nophi-l_normalization**2 * ((Wlmr + Xlmi)**2 + (Xlmr - Wlmi)**2)
        n = n + 1
     end do
  end do
  Psi4relm = Psi4relm * sym_factor * dcostheta * dphi
  Psi4imlm = Psi4imlm * sym_factor * dcostheta * dphi
end subroutine lm_modes
subroutine psikad_compute_Wlm_Xlm(l,m,cost,sint,ph,Wlmr,Wlmi,Xlmr,Xlmi)
  use psikad_lm_functions
  implicit none
  integer :: l,m
  real*8 :: th,ph
  real*8 :: cost,sint,cott,cosmp,sinmp,Plm,fac1,fac1_0,fac2,Pl1
  real*8 :: Plm1,Ylmr,Ylmi,Ylmtr,Ylmti,Wlmr,Wlmi,Xlmr,Xlmi
  real*8, parameter :: PI = 3.14159265358979323846d0
!
  cott = cost/sint
  fac2 = sqrt( (2.D0*l+1.D0)*factorial(l-abs(m)) / ( 4.D0*PI*factorial(l+abs(m)) ) )
  fac1 = fac2 * (l+abs(m))*(l-abs(m)+1.d0)
  fac1_0 = sqrt( (2.D0*l+1.D0)/(4.d0*PI) )
  cosmp = cos( m * ph )
  sinmp = sin( m * ph )
  Plm = plgndr(l,abs(m),cost)
  if(m==0) then
     Pl1 = plgndr(l,1,cost)
     Ylmr = fac2*Plm
     Ylmtr = fac1_0 * Pl1
     Wlmr = -l*(l+1.D0)*Ylmr - 2.D0*cott*Ylmtr
     Xlmr = 0.D0
     Wlmi = 0.d0
     Xlmi = 0.d0
  else
     Plm1 = plgndr(l,abs(m)-1,cost)
     Ylmr = fac2*Plm
     Ylmtr = -fac1*Plm1 - abs(m)*cott*fac2*Plm
     if(m.lt.0) then
        Ylmr = Ylmr*(-1.D0)**(-m)
        Ylmtr = Ylmtr*(-1.D0)**(-m)
     end if
     Wlmr = ( -l*(l+1.d0)*Ylmr -2.d0*cott*Ylmtr + 2.d0*Ylmr*(m/sint)**2 )
     Xlmr = 2.d0*(Ylmtr - cott*Ylmr)*m
     Wlmi = Wlmr * sinmp
     Xlmi = Xlmr*cosmp
     Wlmr = Wlmr*cosmp
     Xlmr = -Xlmr*sinmp
  end if
end subroutine psikad_compute_Wlm_Xlm
