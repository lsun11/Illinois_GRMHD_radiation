!---------------------------------------------------------------------
! This routine computes the areal radius of a spherical surface of 
!   coordinate radius radius_GW and the averaged values of g^tt, g^tr and g^rr
!---------------------------------------------------------------------
subroutine areal_radius_averaged_gutt_gutr_guprr(cctkGH,Symmetry,N_theta,N_phi,ntot, &
     dphi,dcostheta,radius_GW,phi,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,lapm1, &
     betar,r_areal,avg_gutt,avg_gutr,avg_gurr)
  implicit none
  integer*8 :: cctkGH
  real*8                                   :: radius_GW
  integer                                  :: N_theta,N_phi,sym_factor,ntot,n
  real*8                                   :: dcostheta,dphi
  real*8                                   :: r_areal,avg_gutt,avg_gutr,avg_gurr
  real*8, dimension(ntot)                  :: phi,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,lapm1,betar
  real*8                                   :: phiangle,costheta,sintheta
  real*8                                   :: area,x,y,z,dS,xor,yor,zor,f1olap2,psi4
  real*8, parameter                            :: PI = 3.14159265358979323846D0
  integer                                  :: i,j
  integer                                  :: Symmetry
  integer                                  :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
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
  area = 0.d0
  avg_gutt = 0.d0
  avg_gutr = 0.d0
  avg_gurr = 0.d0
  n = 1
  do i=1,N_theta
     costheta = 1.D0 - (i - 0.5D0)*dcostheta
     sintheta = sqrt(1.D0 - costheta*costheta)
     do j=1,N_phi
        phiangle = (j - 0.5D0)*dphi
        if(N_phi==1) phiangle = 0.D0
        !-------------------------------------------------------------
        ! The surface area of a sphere of coordinate radius r is given 
         ! by 
        !    S = Integrate[ psi^6 r^2 sqrt(gamma^{ij} x_i x_j)/r dcostheta dphi ] 
        !      = Integrate[ psi^4 r sqrt(tilde{gamma}^{ij} x_i x_j) dcostheta dphi ]
        !--------------------------------------------------------------
        xor = sintheta*cos(phiangle)
        yor = sintheta*sin(phiangle)
        zor = costheta
        x = radius_GW*xor
        y = radius_GW*yor
        z = radius_GW*zor
        psi4 = exp(4.d0*phi(n))
        dS = psi4*radius_GW*sqrt( gupxx(n)*x*x + &
                 2.d0*gupxy(n)*x*y + 2.d0*gupxz(n)*x*z + &
                 gupyy(n)*y*y + 2.d0*gupyz(n)*y*z + gupzz(n)*z*z )
        area = area + dS
        f1olap2 = 1.d0/(1.d0+lapm1(n))**2
        avg_gutt = avg_gutt - dS*f1olap2
        avg_gutr = avg_gutr + dS*betar(n)*f1olap2
        avg_gurr = avg_gurr + dS*( (gupxx(n)*xor*xor + 2.d0*gupxy(n)*xor*yor + &
                        2.d0*gupxz(n)*xor*zor + gupyy(n)*yor*yor + &
                        2.d0*gupyz(n)*yor*zor + gupzz(n)*zor*zor)/psi4 &
                        - betar(n)*betar(n)*f1olap2)
        n = n + 1
     end do
  end do
  avg_gutt = avg_gutt / area
  avg_gutr = avg_gutr / area
  avg_gurr = avg_gurr / area
  area = area * sym_factor * dcostheta * dphi
  r_areal = sqrt(area/(4.d0*PI))
end subroutine areal_radius_averaged_gutt_gutr_guprr
