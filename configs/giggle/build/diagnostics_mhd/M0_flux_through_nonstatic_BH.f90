subroutine M0_flux_through_nonstatic_BH(cctkGH,F_M0, horizon_number, &
                       t,t_last,xh_last,yh_last,zh_last, &
                       ah_radii_last,N_theta,N_phi,Symmetry,found_horizon)
  implicit none
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
  INTEGER*4 function HorizonCentroid (horizon_number, centroid_x, centroid_y, centroid_z)
  implicit none
  INTEGER*4 horizon_number
  REAL*8 centroid_x
  REAL*8 centroid_y
  REAL*8 centroid_z
  end function HorizonCentroid
  end interface
  interface
  INTEGER*4 function HorizonLocalCoordinateOrigin (horizon_number, origin_x, origin_y, origin_z)
  implicit none
  INTEGER*4 horizon_number
  REAL*8 origin_x
  REAL*8 origin_y
  REAL*8 origin_z
  end function HorizonLocalCoordinateOrigin
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
  interface
  INTEGER*4 function HorizonWasFound (horizon_number)
  implicit none
  INTEGER*4 horizon_number
  end function HorizonWasFound
  end interface
  
  integer*8 :: cctkGH
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  real*8 :: F_M0_ijk,F_M0,t
  real*8 :: vx_i, vy_i, vz_i
  real*8 :: rhos_i
  real*8 :: rr_i, dS, PI, dX, dY, dZ, f1ospi,xh,yh,zh,dRdmu,dRdphi,dRdt,dt
  real*8 :: cosphi,sinphi,dmudx,dmudy,dmudz,dphidx,dphidy,dphidz,sym_factor
  real*8 :: dmudt,dphidt
  real*8 :: sintheta,costheta,nn,xn,yn,phiangle,dphi,dcostheta
  real*8 :: xh_last,yh_last,zh_last,t_last,xh_dot,yh_dot,zh_dot
  integer :: i,j, interp_order, horizon_number, foundflag, n,ntot
  integer :: ind0,ind1,ind2,indm1,indm2,vindex
  integer :: N_theta,N_phi, Symmetry, found_horizon
  real*8, dimension(N_theta*N_phi)          :: ah_radii_last
  real*8, allocatable, dimension(:,:)       :: pointcoords
  real*8, allocatable, dimension(:)            :: nt_d,nx_d,ny_d,nz_d
  real*8, allocatable, dimension(:)         :: ah_radii, x_ah1,y_ah1,z_ah1
  real*8, allocatable, dimension(:)         :: vxint, vyint, vzint
  real*8, allocatable, dimension(:)         :: rhosint
!
  foundflag = HorizonWasFound(horizon_number)
  ! Horizon not found, set the flux to 0
  if (foundflag .ne. 1) then
     F_M0 = 0.d0
     found_horizon = 0
     return
  end if
  found_horizon = 1
  if (Symmetry .ne. EQUATORIAL .and. Symmetry .ne. NO_SYMM) then
     write(*,*) 'Symmetry not supported in M0_flux_through_nonstatic_BH'
     stop
  end if
  ntot = N_theta*N_phi
  dt = t-t_last
  ! allocate memory
  allocate(pointcoords(ntot,3))
  allocate(nt_d(ntot))
  allocate(nx_d(ntot))
  allocate(ny_d(ntot))
  allocate(nz_d(ntot))
  allocate(ah_radii(ntot))
  allocate(x_ah1(ntot))
  allocate(y_ah1(ntot))
  allocate(z_ah1(ntot))
  allocate(vxint(ntot))
  allocate(vyint(ntot))
  allocate(vzint(ntot))
  allocate(rhosint(ntot))
  PI = 3.14159265358979323846D0
  dphi = 2.d0 * PI / N_phi
  if (Symmetry==EQUATORIAL) then
     dcostheta = 1.d0 / N_theta
     sym_factor = 2.d0
  else
     dcostheta = 2.d0 / N_theta
     sym_factor = 1.d0
  end if
! Get the origin of BH
  foundflag = HorizonLocalCoordinateOrigin(horizon_number,xh,yh,zh)
  xh_dot = (xh - xh_last)/dt
  yh_dot = (yh - yh_last)/dt
  zh_dot = (zh - zh_last)/dt
  n = 1
  do i=1,N_theta
     costheta = 1.d0 - (i - 0.5d0)*dcostheta
     sintheta = sqrt(1.d0 - costheta*costheta)
     do j=1,N_phi
        phiangle = (j - 0.5d0)*dphi
        x_ah1(n) = xh + sintheta*cos(phiangle)
        y_ah1(n) = yh + sintheta*sin(phiangle)
        z_ah1(n) = zh + costheta
        n = n + 1
     end do
  end do
  ! Find horizon radii
  foundflag = HorizonRadiusInDirection(horizon_number,ntot,x_ah1,y_ah1,z_ah1,ah_radii)
  ! Now set the points on the horizon surface 
  n = 1
  do i=1,N_theta
     costheta = 1.d0 - (i - 0.5d0)*dcostheta
     sintheta = sqrt(1.d0 - costheta*costheta)
     do j=1,N_phi
        phiangle = (j - 0.5d0)*dphi
        pointcoords(n,1) = xh + ah_radii(n)*sintheta*cos(phiangle)
        pointcoords(n,2) = yh + ah_radii(n)*sintheta*sin(phiangle)
        pointcoords(n,3) = zh + ah_radii(n)*costheta
        n = n + 1
     end do
  end do
  ! Let f = sqrt[(x-xh)^2+(y-yh)^2+(z-zh)^2] - R(mu,phi), where mu=cos(theta). 
  ! Now compute partial f / partial x^\mu and normalize 
  ! by the Jacobian.
  ! 
  n=1
  do i=1,N_theta
     costheta = 1.d0 - (i - 0.5d0)*dcostheta
     sintheta = sqrt(1.d0 - costheta*costheta)
     do j=1,N_phi
        phiangle = (j - 0.5d0)*dphi
        cosphi = cos(phiangle)
        sinphi = sin(phiangle)
        ! dR/dmu
        if (i==1) then
              ind0 = j
           ind1 = j + N_phi
           ind2 = j + 2*N_phi
           dRdmu = (1.5d0*ah_radii(ind0) - 2.d0*ah_radii(ind1) + 0.5d0*ah_radii(ind2)) / dcostheta
        else if (i==N_theta) then
           if (Symmetry==EQUATORIAL) then
              ind1 = j + (i-1)*N_phi
              indm1 = j + (i-2)*N_phi
              dRdmu = 0.5d0*(ah_radii(indm1)-ah_radii(ind1))/dcostheta
           else
              ind0 = j + (i-1)*N_phi
              indm1 = j + (i-2)*N_phi
              indm2 = j + (i-3)*N_phi
              dRdmu = (-1.5d0*ah_radii(ind0) + 2.d0*ah_radii(indm1) - 0.5d0*ah_radii(indm2)) / dcostheta
           end if
        else
           ind1 = j + i*N_phi
           indm1 = j + (i-2)*N_phi
           dRdmu = 0.5d0*(ah_radii(indm1)-ah_radii(ind1))/dcostheta
        end if
        ! dR/dphi
        if (j==1) then
           ind1 = n+1
            indm1 = i*N_phi
        else if (j==N_phi) then
           ind1 = 1 + (i-1)*N_phi
           indm1 = n-1
        else
           ind1 = n+1
           indm1 = n-1
        end if
        dRdphi = 0.5d0*(ah_radii(ind1)-ah_radii(indm1))/dphi
        ! dR/dt
            dRdt = (ah_radii(n) - ah_radii_last(n))/dt
         dmudt = (costheta*sintheta*cosphi*xh_dot   &
                 + costheta*sintheta*sinphi*yh_dot &
                 - sintheta*sintheta*zh_dot) / ah_radii(n)
        dmudx = -costheta*sintheta*cosphi/ah_radii(n)
        dmudy = -costheta*sintheta*sinphi/ah_radii(n)
        dmudz = sintheta*sintheta/ah_radii(n)
            dphidt = (sinphi*xh_dot - cosphi*yh_dot)/(ah_radii(n)*sintheta)
        dphidx = -sinphi/(ah_radii(n)*sintheta)
          dphidy = cosphi/(ah_radii(n)*sintheta)
         nt_d(n) = -sintheta*cosphi*xh_dot - sintheta*sinphi*yh_dot &
                  - costheta*zh_dot - dRdt - dRdmu*dmudt - dRdphi*dphidt
        nx_d(n) = sintheta*cosphi - dRdmu*dmudx - dRdphi*dphidx
        ny_d(n) = sintheta*sinphi - dRdmu*dmudy - dRdphi*dphidy
        nz_d(n) = costheta - dRdmu*dmudz
        nn = (nx_d(n)*cosphi + ny_d(n)*sinphi)*sintheta + nz_d(n)*costheta
        ! normalize by the Jacobian.
        nt_d(n) = nt_d(n)/nn
        nx_d(n) = nx_d(n)/nn
        ny_d(n) = ny_d(n)/nn
        nz_d(n) = nz_d(n)/nn
        ! update ah_radii_last
        ah_radii_last(n) = ah_radii(n)
        n = n + 1
     end do
  end do
  ! update t_last, xh_last, yh_last and zh_last
  t_last = t
  xh_last = xh
  yh_last = yh
  zh_last = zh
! Interpolate the grid functions onto the AH
  call CCTK_VarIndex(vindex,"mhd_evolve::vx")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,vxint)
  call CCTK_VarIndex(vindex,"mhd_evolve::vy")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,vyint)
  call CCTK_VarIndex(vindex,"mhd_evolve::vz")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,vzint)
  call CCTK_VarIndex(vindex,"mhd_evolve::rho_star")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,rhosint)
  !
  ! integrate on the horizon
  !
  F_M0 = 0.d0
  do i=1,ntot
     ! read in data
     vx_i  = vxint(i)
     vy_i  = vyint(i)
     vz_i  = vzint(i)
     rhos_i = rhosint(i)
     F_M0_ijk = (nt_d(i) + nx_d(i)*vx_i + ny_d(i)*vy_i + nz_d(i)*vz_i)*rhos_i
     dS = dphi * dcostheta * sym_factor * ah_radii(i)**2
     F_M0  = F_M0 + F_M0_ijk*dS
  end do
  ! deallocate memory
  deallocate(pointcoords)
  deallocate(nt_d, nx_d, ny_d, nz_d, ah_radii, x_ah1, y_ah1, z_ah1)
  deallocate(vxint, vyint, vzint, rhosint)
end subroutine M0_flux_through_nonstatic_BH
