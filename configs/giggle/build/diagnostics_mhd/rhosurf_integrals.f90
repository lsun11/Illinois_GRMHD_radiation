subroutine rhosurf_integrals(cctkGH,F_M0,rho_const_surf,surf_area,ave_radius,t,t_last,xh_last,yh_last,zh_last,surf_radii_last,N_rad&
  &_ray,N_theta,N_phi,rhosurf_rmin,rhosurf_rmax,Symmetry,use_origin)
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
  character                                :: filename*50
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  real*8 :: F_M0_ijk,F_M0,rho_const_surf,surf_area,ave_radius,t
  real*8 :: vx_i, vy_i, vz_i
  real*8 :: rhos_i
  real*8 :: rr_i, dsurf_area,dS, PI, dX, dY, dZ, f1ospi,xh,yh,zh,dRdmu,dRdphi,dRdt,dt
  real*8 :: cosphi,sinphi,dmudx,dmudy,dmudz,dmudt,dphidx,dphidy,dphidz,dphidt,sym_factor
  real*8 :: sintheta,costheta,nn,xn,yn,phiangle,dphi,dcostheta
  integer :: i,j, interp_order, horizon_number, foundflag, n,ntot,ii
  integer :: ind0,ind1,ind2,indm1,indm2,vindex
  integer :: N_theta,N_phi,N_rad_ray,Symmetry,use_origin
  real*8 :: rhosurf_rmin,rhosurf_rmax
  real*8 :: xh_last,yh_last,zh_last,t_last,xh_dot,yh_dot,zh_dot
  real*8, allocatable,dimension(:,:)         :: pointcoordsray
  real*8, allocatable,dimension(:)              :: rho_b_interp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real*8, dimension(N_theta*N_phi)          :: surf_radii_last
  real*8, allocatable,dimension(:,:)        :: pointcoords
  real*8, allocatable,dimension(:)          :: nx_d,ny_d,nz_d,nt_d
  real*8, allocatable,dimension(:)          :: surf_radii, phiangle_vec,costheta_vec
  real*8, allocatable,dimension(:)          :: gupxx_int,gupxy_int,gupxz_int
  real*8, allocatable,dimension(:)          :: gupyy_int,gupyz_int,gupzz_int
  real*8, allocatable,dimension(:)          :: vxint,vyint,vzint
  real*8, allocatable,dimension(:)          :: rhosint,rho_b_int
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real*8                                    :: gupxx_i,gupxy_i,gupxz_i
  real*8                                    :: gupyy_i,gupyz_i,gupzz_i
  real*8                                    :: xp1,yp1,zp1
  real*8                                    :: dxray,rhoold,rhonew,frac,rsurf
  if (Symmetry .ne. EQUATORIAL) then
     write(*,*) 'Symmetry not supported in rhoflux_brian'
     stop
  end if
  sym_factor = 2.d0
  ntot = N_theta*N_phi
  dt = t-t_last
  ! allocate memory
  allocate(pointcoordsray(N_rad_ray,3))
  allocate(rho_b_interp(N_rad_ray))
  allocate(pointcoords(ntot,3))
  allocate(nt_d(ntot))
  allocate(nx_d(ntot))
  allocate(ny_d(ntot))
  allocate(nz_d(ntot))
  allocate(phiangle_vec(ntot))
  allocate(costheta_vec(ntot))
  allocate(surf_radii(ntot))
  allocate(vxint(ntot))
  allocate(vyint(ntot))
  allocate(vzint(ntot))
  allocate(gupxx_int(ntot))
  allocate(gupxy_int(ntot))
  allocate(gupxz_int(ntot))
  allocate(gupyy_int(ntot))
  allocate(gupyz_int(ntot))
  allocate(gupzz_int(ntot))
  allocate(rhosint(ntot))
  allocate(rho_b_int(ntot))
  PI = 3.14159265358979323844D0
  dphi = 2.0 * PI / N_phi
  dcostheta = 1.0 / N_theta
  horizon_number = 1
  if (use_origin .eq. 1) then
     xh=0.d0
     yh=0.d0
     zh=0.d0
  else
     ! Get the origin of BH
     foundflag = HorizonLocalCoordinateOrigin(horizon_number,xh,yh,zh)
  endif
  xh_dot = (xh - xh_last)/dt
  yh_dot = (yh - yh_last)/dt
  zh_dot = (zh - zh_last)/dt
  !THIS SHOULD NOT BE HARDCODED
  dxray = 0.1
  n=1
  do i=1,N_theta
     costheta = 1.0 - (i - 0.5)*dcostheta
     sintheta = sqrt(1.0 - costheta*costheta)
     do j=1,N_phi
        phiangle = (j - 0.5)*dphi
        phiangle_vec(n) = phiangle
        costheta_vec(n) = costheta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Find surface radii
        do ii=1,N_rad_ray
           xp1 = ((ii-1)*dxray+rhosurf_rmin)*sintheta*cos(phiangle) + xh
           yp1 = ((ii-1)*dxray+rhosurf_rmin)*sintheta*sin(phiangle) + yh
           zp1 = ((ii-1)*dxray+rhosurf_rmin)*costheta + zh
           if(abs(xp1-xh).lt.rhosurf_rmax.and.abs(yp1-yh).lt.rhosurf_rmax.and.abs(zp1-zh).lt.rhosurf_rmax) then
              pointcoordsray(ii,1) = xp1
              pointcoordsray(ii,2) = yp1
              pointcoordsray(ii,3) = zp1
           else
              if(ii.eq.1)then
                 write(6,*)'impossible!',xp1,yp1,zp1,rhosurf_rmax
                 stop
              endif
              pointcoordsray(ii,1) = pointcoordsray(ii-1,1)
              pointcoordsray(ii,2) = pointcoordsray(ii-1,2)
              pointcoordsray(ii,3) = pointcoordsray(ii-1,3)
           endif
        end do
        call CCTK_VarIndex(vindex,"mhd_evolve::rho_b")
        call interp_driver_carp(cctkGH,N_rad_ray,pointcoordsray,vindex,rho_b_interp)
        foundflag=0
        rsurf=-1.D0
              rhoold=rho_b_interp(1)
        do ii=2,N_rad_ray
           rhonew = rho_b_interp(ii)
           if (rhonew.lt.rho_const_surf .and. rhoold.gt.rho_const_surf) then
              frac=(rho_const_surf-rhoold)/(rhonew-rhoold)
              surf_radii(n)=(ii+frac-2)*dxray+rhosurf_rmin
              foundflag = 1
              EXIT
           else
              rhoold=rhonew
           end if
        end do
        if (foundflag .eq. 0) then
           F_M0 = -10.0
           return
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             n = n + 1
     end do
  end do
  ! Now set the points on the horizon surface 
  n = 1
  do i=1,N_theta
     costheta = 1.0 - (i - 0.5)*dcostheta
     sintheta = sqrt(1.0 - costheta*costheta)
    do j=1,N_phi
        phiangle = (j - 0.5)*dphi
        pointcoords(n,1) = xh + surf_radii(n)*sintheta*cos(phiangle)
        pointcoords(n,2) = yh + surf_radii(n)*sintheta*sin(phiangle)
        pointcoords(n,3) = zh + surf_radii(n)*costheta
        n = n + 1
     end do
  end do
  ! Let f = sqrt[(x-xh)^2+(y-yh)^2+(z-zh)^2] - R(mu,phi), where mu=cos(theta). 
  ! Now compute partial f / partial X^i (X^i = x^i-xh^i) and normalize 
  ! by the Jacobian, storing them to ni_d. 
  ! 
  n=1
  do i=1,N_theta
     costheta = 1.0 - (i - 0.5)*dcostheta
     sintheta = sqrt(1.0 - costheta*costheta)
     do j=1,N_phi
        phiangle = (j - 0.5)*dphi
        cosphi = cos(phiangle)
        sinphi = sin(phiangle)
        ! dR/dmu 
        if (i==1) then
              ind0 = j
           ind1 = j + N_phi
           ind2 = j + 2*N_phi
           dRdmu = (1.5d0*surf_radii(ind0) - 2.d0*surf_radii(ind1) + 0.5d0*surf_radii(ind2)) / dcostheta
        else if (i==N_theta) then
           ind1 = j + (i-1)*N_phi
           indm1 = j + (i-2)*N_phi
           dRdmu = 0.5d0*(surf_radii(indm1)-surf_radii(ind1))/dcostheta
        else
           ind1 = j + i*N_phi
           indm1 = j + (i-2)*N_phi
           dRdmu = 0.5d0*(surf_radii(indm1)-surf_radii(ind1))/dcostheta
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
        dRdphi = 0.5d0*(surf_radii(ind1)-surf_radii(indm1))/dphi
        ! dR/dt
        dRdt = (surf_radii(n) - surf_radii_last(n))/dt
        dmudt = (costheta*sintheta*cosphi*xh_dot   &
             + costheta*sintheta*sinphi*yh_dot &
             - sintheta*sintheta*zh_dot) / surf_radii(n)
        dmudx = -costheta*sintheta*cosphi/surf_radii(n)
        dmudy = -costheta*sintheta*sinphi/surf_radii(n)
        dmudz = sintheta*sintheta/surf_radii(n)
        dphidx = -sinphi/(surf_radii(n)*sintheta)
          dphidy = cosphi/(surf_radii(n)*sintheta)
         nt_d(n) = -sintheta*cosphi*xh_dot - sintheta*sinphi*yh_dot &
             - costheta*zh_dot - dRdt - dRdmu*dmudt - dRdphi*dphidt
        nx_d(n) = sintheta*cosphi - dRdmu*dmudx - dRdphi*dphidx
        ny_d(n) = sintheta*sinphi - dRdmu*dmudy - dRdphi*dphidy
        nz_d(n) = costheta - dRdmu*dmudz
        nn = (nx_d(n)*cosphi + ny_d(n)*sinphi)*sintheta + nz_d(n)*costheta
        !        ! normalize by the Jacobian.
        nx_d(n) = nx_d(n)/nn
        ny_d(n) = ny_d(n)/nn
        nz_d(n) = nz_d(n)/nn
        ! update surf_radii_last
        surf_radii_last(n) = surf_radii(n)
        n = n + 1
     end do
  end do
!!$  ! update t_last, xh_last, yh_last and zh_last
!!$  if (reset .eq. 1) then
!!$     t_last = t
!!$     xh_last = xh
!!$     yh_last = yh
!!$     zh_last = zh
!!$  endif
  ! Interpolate the grid functions onto the Surface
  call CCTK_VarIndex(vindex,"mhd_evolve::vx")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,vxint)
  call CCTK_VarIndex(vindex,"mhd_evolve::vy")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,vyint)
  call CCTK_VarIndex(vindex,"mhd_evolve::vz")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,vzint)
  call CCTK_VarIndex(vindex,"mhd_evolve::rho_star")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,rhosint)
  call CCTK_VarIndex(vindex,"mhd_evolve::rho_b")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,rho_b_int)
!!$  write(filename,31)t
!!$31 FORMAT(check_surf_,E)
!!$  open(UNIT=16,FILE=filename,STATUS=REPLACE)
!!$  do i=1,ntot
!!$     write(16,(5e18.10)) pointcoords(i,1),pointcoords(i,2),pointcoords(i,3),rho_b_int(i),sqrt((pointcoords(i,1)-xh)**2+(pointcoor
!!$  end do
!!$  close(16)
!!$  
  !! Interpolate the grid functions onto the Surface
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
  !  F_M0 = 0.d0
  surf_area = 0.d0
  F_M0 = 0.d0
  ave_radius=0.d0
  do i=1,ntot
     ! read in data
     vx_i  = vxint(i)
     vy_i  = vyint(i)
     vz_i  = vzint(i)
     rhos_i = rhosint(i)
     F_M0_ijk = (nt_d(i) + nx_d(i)*vx_i + ny_d(i)*vy_i + nz_d(i)*vz_i)*rhos_i
     gupxx_i = gupxx_int(i)
     gupxy_i = gupxy_int(i)
     gupxz_i = gupxz_int(i)
     gupyy_i = gupyy_int(i)
     gupyz_i = gupyz_int(i)
     gupzz_i = gupzz_int(i)
     dsurf_area = sqrt(gupxx_i*nx_d(i)*nx_d(i) + &
          gupyy_i*ny_d(i)*ny_d(i) + &
          gupzz_i*nz_d(i)*nz_d(i) + &
          2.d0*gupxy_i*nx_d(i)*ny_d(i) + &
          2.d0*gupxz_i*nx_d(i)*nz_d(i) + &
          2.d0*gupyz_i*ny_d(i)*nz_d(i))* &
          dphi * dcostheta * sym_factor * surf_radii(i)**2
     dS = dphi * dcostheta * sym_factor * surf_radii(i)**2
     F_M0  = F_M0 + F_M0_ijk*dS
     surf_area = surf_area + dsurf_area
     ave_radius = ave_radius + dsurf_area*surf_radii(i)
  end do
  ave_radius = ave_radius / surf_area
  ! deallocate memory
  deallocate(pointcoords)
  deallocate(pointcoordsray,rho_b_interp)
  deallocate(nt_d, nx_d, ny_d, nz_d, surf_radii)
  deallocate(vxint, vyint, vzint, rhosint,rho_b_int)
  deallocate(gupxx_int,gupxy_int,gupxz_int,gupyy_int,gupyz_int,gupzz_int)
  deallocate(phiangle_vec,costheta_vec)
end subroutine rhosurf_integrals
