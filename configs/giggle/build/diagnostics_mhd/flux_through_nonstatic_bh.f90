subroutine flux_through_nonstatic_BH(cctkGH,F_E_fluid,F_E_em,F_M0,F_J_fluid,F_J_em, horizon_number, &
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
  real*8 :: F_M0,f_E_fluid,F_E_em,f_J_fluid,F_J_em
  real*8 :: F_M0_ijk,F_E_fluid_ijk,F_E_em_ijk,F_J_fluid_ijk,F_J_em_ijk
  real*8 :: F_E_emx,F_E_emy,F_E_emz,F_E_emt
  real*8 :: F_E_fluidx,F_E_fluidy,F_E_fluidz,F_E_fluidt
  real*8 :: F_J_emx,F_J_emy,F_J_emz,F_J_emt
  real*8 :: F_J_fluidx,F_J_fluidy,F_J_fluidz,F_J_fluidt
  real*8 :: t
  real*8 :: vx_i, vy_i, vz_i,u0_i,u_0_i,u_x_i,u_y_i,u_phi_i
  real*8 :: lapse_i,shiftx_i,shifty_i,shiftz_i,shift_x_i,shift_y_i,shift_z_i
  real*8 :: shift2_i,shiftv_i
  real*8 :: gxx_i,gxy_i,gxz_i,gyy_i,gyz_i,gzz_i,psi4_i,psi6_i
  real*8 :: rhos_i,h_i,P_i
  real*8 :: bx_i,by_i,bz_i,bt_i,b2,shiftb_i
  real*8 :: b_x_i,b_y_i,b_z_i,b_t_i,b_phi_i
  real*8 :: rr_i, dS, PI, dX, dY, dZ, f1ospi,xh,yh,zh,dRdmu,dRdphi,dRdt,dt
  real*8 :: xh_local,yh_local,zh_local,xh_centroid,yh_centroid,zh_centroid
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
  real*8, allocatable, dimension(:)         :: vxint, vyint, vzint, u0int
  real*8, allocatable, dimension(:)         :: rhosint,hint,Pint
  real*8, allocatable, dimension(:)         :: lapm1int,shiftxint,shiftyint,shiftzint
  real*8, allocatable, dimension(:)         :: gxxint,gxyint,gxzint,gyyint,gyzint,gzzint,phiint
  real*8, allocatable, dimension(:)         :: sbxint,sbyint,sbzint,sbtint
  !
  foundflag = HorizonWasFound(horizon_number)
  ! Horizon not found, set the flux to 0
  if (foundflag .ne. 1) then
     F_E_fluid = 0.d0
     F_E_em = 0.d0
     f_M0 = 0.d0
     F_J_fluid = 0.d0
     F_J_em = 0.d0
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
  allocate(u0int(ntot))
  allocate(rhosint(ntot))
  allocate(hint(ntot))
  allocate(Pint(ntot))
  allocate(lapm1int(ntot))
  allocate(shiftxint(ntot))
  allocate(shiftyint(ntot))
  allocate(shiftzint(ntot))
  allocate(gxxint(ntot))
  allocate(gxyint(ntot))
  allocate(gxzint(ntot))
  allocate(gyyint(ntot))
  allocate(gyzint(ntot))
  allocate(gzzint(ntot))
  allocate(phiint(ntot))
  allocate(sbxint(ntot))
  allocate(sbyint(ntot))
  allocate(sbzint(ntot))
  allocate(sbtint(ntot))
  PI = 3.14159265358979323846D0
  dphi = 2.d0 * PI / N_phi
  if (Symmetry==EQUATORIAL) then
     dcostheta = 1.d0 / N_theta
     sym_factor = 2.d0
  else
     dcostheta = 2.d0 / N_theta
     sym_factor = 1.d0
  end if
  write(*,*) "***************"
  write(*,*) "horizon_number: ",horizon_number
! Get the origin of BH!  foundflag = HorizonLocalCoordinateOrigin(horizon_number,xh,yh,zh)
  foundflag = HorizonLocalCoordinateOrigin(horizon_number,xh_local,yh_local,zh_local)
  foundflag = HorizonCentroid(horizon_number,xh_centroid,yh_centroid,zh_centroid);
  write(*,*) "xh_last: ",xh_last
  write(*,*) "yh_last: ",yh_last
  write(*,*) "zh_last: ",zh_last
  write(*,*) "t_last: ",t_last
  write(*,*) "xh_local: ",xh_local
  write(*,*) "yh_local: ",yh_local
  write(*,*) "zh_local: ",zh_local
  write(*,*) "xh_centroid: ",xh_centroid
  write(*,*) "yh_centroid: ",yh_centroid
  write(*,*) "zh_centroid: ",zh_centroid
  xh=xh_local
  yh=yh_local
  zh=zh_local
  xh_dot = (xh - xh_last)/dt
  yh_dot = (yh - yh_last)/dt
  zh_dot = (zh - zh_last)/dt
  write(*,*) "xh_dot: ",xh_dot
  write(*,*) "yh_dot: ",yh_dot
  write(*,*) "zh_dot: ",zh_dot
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
  write(*,*) "ah_raddi(1): ",ah_radii(1)
  write(*,*) "ah_raddi(100): ",ah_radii(100)
  write(*,*) "ah_raddi(200): ",ah_radii(200)
  write(*,*) "***************"
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
  ! These are now updated outside this subroutine.
  !!t_last = t
  !!xh_last = xh
  !!yh_last = yh
  !!zh_last = zh
! Interpolate the grid functions onto the AH
  call CCTK_VarIndex(vindex,"mhd_evolve::vx")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,vxint)
  call CCTK_VarIndex(vindex,"mhd_evolve::vy")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,vyint)
  call CCTK_VarIndex(vindex,"mhd_evolve::vz")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,vzint)
  call CCTK_VarIndex(vindex,"mhd_evolve::u0")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,u0int)
  call CCTK_VarIndex(vindex,"mhd_evolve::rho_star")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,rhosint)
  call CCTK_VarIndex(vindex,"mhd_evolve::h")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,hint)
  call CCTK_VarIndex(vindex,"mhd_evolve::P")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,Pint)
  call CCTK_VarIndex(vindex,"lapse::lapm1")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,lapm1int)
  call CCTK_VarIndex(vindex,"shift::shiftx")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,shiftxint)
  call CCTK_VarIndex(vindex,"shift::shifty")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,shiftyint)
  call CCTK_VarIndex(vindex,"shift::shiftz")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,shiftzint)
  call CCTK_VarIndex(vindex,"bssn::gxx")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gxxint)
  call CCTK_VarIndex(vindex,"bssn::gxy")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gxyint)
  call CCTK_VarIndex(vindex,"bssn::gxz")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gxzint)
  call CCTK_VarIndex(vindex,"bssn::gyy")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gyyint)
  call CCTK_VarIndex(vindex,"bssn::gyz")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gyzint)
  call CCTK_VarIndex(vindex,"bssn::gzz")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gzzint)
  call CCTK_VarIndex(vindex,"bssn::phi")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,phiint)
  call CCTK_VarIndex(vindex,"mhd_evolve::sbx")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,sbxint)
  call CCTK_VarIndex(vindex,"mhd_evolve::sby")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,sbyint)
  call CCTK_VarIndex(vindex,"mhd_evolve::sbz")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,sbzint)
  call CCTK_VarIndex(vindex,"mhd_evolve::sbt")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,sbtint)
  !
  ! integrate on the horizon
  !
  F_E_fluid = 0.d0
  F_E_em = 0.d0
  F_M0 = 0.d0
  F_J_fluid = 0.d0
  F_J_em = 0.d0
  do i=1,ntot
     ! read in data
     vx_i  = vxint(i)
     vy_i  = vyint(i)
     vz_i  = vzint(i)
     u0_i  = u0int(i)
     lapse_i = lapm1int(i)+1.d0
     shiftx_i = shiftxint(i)
     shifty_i = shiftyint(i)
     shiftz_i = shiftzint(i)
     gxx_i = gxxint(i)
     gxy_i = gxyint(i)
     gxz_i = gxzint(i)
     gyy_i = gyyint(i)
     gyz_i = gyzint(i)
     gzz_i = gzzint(i)
     psi4_i = exp(4.d0*phiint(i))
     psi6_i = exp(6.d0*phiint(i))
     rhos_i = rhosint(i)
     h_i = hint(i)
     P_i = Pint(i)
     bx_i = sbxint(i)
     by_i = sbyint(i)
     bz_i = sbzint(i)
     bt_i = sbtint(i)
     shift_x_i = psi4_i*(gxx_i*shiftx_i+gxy_i*shifty_i+gxz_i*shiftz_i)
     shift_y_i = psi4_i*(gxy_i*shiftx_i+gyy_i*shifty_i+gyz_i*shiftz_i)
     shift_z_i = psi4_i*(gxz_i*shiftx_i+gyz_i*shifty_i+gzz_i*shiftz_i)
     shift2_i = shift_x_i*shiftx_i + shift_y_i*shifty_i + shift_z_i*shiftz_i
     shiftv_i = shift_x_i*vx_i + shift_y_i*vy_i + shift_z_i*vz_i
     shiftb_i = shift_x_i*bx_i + shift_y_i*by_i + shift_z_i*bz_i
     u_0_i = u0_i*(-lapse_i*lapse_i + shift2_i + shiftv_i)
     u_x_i = u0_i*psi4_i*(gxx_i*(shiftx_i+vx_i) +&
          gxy_i*(shifty_i+vy_i) +&
          gxz_i*(shiftz_i+vz_i))
     u_y_i = u0_i*psi4_i*(gxy_i*(shiftx_i+vx_i) +&
          gyy_i*(shifty_i+vy_i) +&
          gyz_i*(shiftz_i+vz_i))
     u_phi_i = -pointcoords(i,2)*u_x_i + pointcoords(i,1)*u_y_i
     b_t_i = bt_i*(-lapse_i*lapse_i + shift2_i) + shiftb_i
     b_x_i = shift_x_i*bt_i + psi4_i * (gxx_i*bx_i+gxy_i*by_i+gxz_i*bz_i)
     b_y_i = shift_y_i*bt_i + psi4_i * (gxy_i*bx_i+gyy_i*by_i+gyz_i*bz_i)
     b_z_i = shift_z_i*bt_i + psi4_i * (gxz_i*bx_i+gyz_i*by_i+gzz_i*bz_i)
     b_phi_i = -pointcoords(i,2)*b_x_i + pointcoords(i,1)*b_y_i
     b2 = bt_i*b_t_i+bx_i*b_x_i+by_i*b_y_i+bz_i*b_z_i
     F_E_fluidt = rhos_i*h_i*u_0_i + lapse_i*psi6_i*P_i
     F_E_fluidx = rhos_i*h_i*vx_i*u_0_i
     F_E_fluidy = rhos_i*h_i*vy_i*u_0_i
     F_E_fluidz = rhos_i*h_i*vz_i*u_0_i
     F_E_emt = lapse_i*psi6_i*(0.5d0*b2 + b2*u0_i*u_0_i - bt_i*b_t_i)
     F_E_emx = lapse_i*psi6_i*(b2*u0_i*vx_i*u_0_i - bx_i*b_t_i)
     F_E_emy = lapse_i*psi6_i*(b2*u0_i*vy_i*u_0_i - by_i*b_t_i)
     F_E_emz = lapse_i*psi6_i*(b2*u0_i*vz_i*u_0_i - bz_i*b_t_i)
     F_J_fluidt = rhos_i*h_i*u_phi_i
     F_J_fluidx = rhos_i*h_i*vx_i*u_phi_i-lapse_i*psi6_i*pointcoords(i,2)*P_i
     F_J_fluidy = rhos_i*h_i*vy_i*u_phi_i+lapse_i*psi6_i*pointcoords(i,1)*P_i
     F_J_fluidz = rhos_i*h_i*vz_i*u_phi_i
     F_J_emt = lapse_i*psi6_i*(u0_i*u_phi_i*b2 - bt_i*b_phi_i)
     F_J_emx = lapse_i*psi6_i*(-0.5d0*b2*pointcoords(i,2)+&
          u0_i*vx_i*u_phi_i*b2 - bx_i*b_phi_i)
     F_J_emy = lapse_i*psi6_i*(0.5d0*b2*pointcoords(i,1)+&
          u0_i*vy_i*u_phi_i*b2 - by_i*b_phi_i)
     F_J_emz =lapse_i*psi6_i*(u0_i*vz_i*u_phi_i*b2 - bz_i*b_phi_i)
     F_M0_ijk = (nt_d(i) + nx_d(i)*vx_i + ny_d(i)*vy_i + nz_d(i)*vz_i)*rhos_i
     F_E_fluid_ijk = -(nt_d(i)*F_E_fluidt+nx_d(i)*F_E_fluidx+ny_d(i)*F_E_fluidy+nz_d(i)*F_E_fluidz)
     F_E_em_ijk = -(nt_d(i)*F_E_emt + nx_d(i)*F_E_emx + ny_d(i)*F_E_emy + nz_d(i)*F_E_emz)
     F_J_fluid_ijk = -(nt_d(i)*F_J_fluidt+nx_d(i)*F_J_fluidx+ny_d(i)*F_J_fluidy+nz_d(i)*F_J_fluidz)
     F_J_em_ijk = -(nt_d(i)*F_J_emt + nx_d(i)*F_J_emx + ny_d(i)*F_J_emy + nz_d(i)*F_J_emz)
     dS = dphi * dcostheta * sym_factor * ah_radii(i)**2
     F_M0  = F_M0 + F_M0_ijk*dS
     F_E_fluid = F_E_fluid + F_E_fluid_ijk*dS
     F_E_em = F_E_em + F_E_em_ijk*dS
     F_J_fluid = F_J_fluid + F_J_fluid_ijk*dS
     F_J_em = F_J_em + F_J_em_ijk*dS
  end do
  ! deallocate memory
  deallocate(pointcoords)
  deallocate(nt_d, nx_d, ny_d, nz_d, ah_radii, x_ah1, y_ah1, z_ah1)
  deallocate(vxint, vyint, vzint, rhosint,hint,Pint)
  deallocate(lapm1int,shiftxint,shiftyint,shiftzint)
  deallocate(gxxint,gxyint,gxzint,gyyint,gyzint,gzzint,phiint)
  deallocate(sbxint,sbyint,sbzint,sbtint)
end subroutine flux_through_nonstatic_BH
