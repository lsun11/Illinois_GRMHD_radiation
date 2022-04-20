#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine BH_local_diagnostics(cctkGH,Mirr,M_BH,J_BH, xi_err, & 
                       horizon_number,N_theta,N_phi,Symmetry,found_horizon)
  implicit none
  DECLARE_CCTK_FUNCTIONS
  CCTK_POINTER :: cctkGH
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  real*8 :: Mirr,M_BH,J_BH
  real*8 :: PI, xh,yh,zh, sym_factor 
  real*8 :: sintheta,costheta,phiangle,dphi,dcostheta,xi_err
  integer :: i,j, horizon_number, foundflag, n,ntot
  integer :: N_theta,N_phi, Symmetry, found_horizon
  real*8, allocatable, dimension(:)         :: xinterp,yinterp,zinterp
  real*8, allocatable, dimension(:)         :: ah_radii, x_ah1,y_ah1,z_ah1
  real*8, allocatable, dimension(:)         :: Psi4int,gxxint,gxyint,gxzint
  real*8, allocatable, dimension(:)         :: gyyint,gyzint,gzzint
  real*8, allocatable, dimension(:)         :: Kint,Kxxint,Kxyint,Kxzint
  real*8, allocatable, dimension(:)         :: Kyyint,Kyzint,Kzzint
  real*8, allocatable, dimension(:)         :: RKt,RKp,dArea
  real*8, allocatable, dimension(:)         :: qtt,qtp,qpp,detq
  real*8, allocatable, dimension(:)         :: quptt,quptp,quppp,R_2d
  real*8, allocatable, dimension(:)         :: Gam_ttt,Gam_ttp,Gam_tpp
  real*8, allocatable, dimension(:)         :: Gam_ptt,Gam_ptp,Gam_ppp
  real*8, allocatable, dimension(:)         :: xi_theta, xi_phi
!
  foundflag = HorizonWasFound(horizon_number)
  ! Horizon not found, quit
  if (foundflag .ne. 1) then
     found_horizon = 0
     Mirr = 0.d0
     M_BH = 0.d0
     J_BH = 0.d0
     xi_err = 0.d0
     return
  end if
  
  found_horizon = 1

!  if (Symmetry .ne. EQUATORIAL .and. Symmetry .ne. NO_SYMM) then 
  if (Symmetry .ne. EQUATORIAL) then 
     write(*,*) 'Symmetry not supported in BH_local_diagnostics'
     stop
  end if 
 
  sym_factor = 2.d0
!  if (Symmetry==NO_SYMM) sym_factor = 1.d0
  ntot = N_theta*N_phi

  ! allocate memory
  allocate(xinterp(ntot))
  allocate(yinterp(ntot))
  allocate(zinterp(ntot))
  allocate(ah_radii(ntot))
  allocate(x_ah1(ntot))
  allocate(y_ah1(ntot))
  allocate(z_ah1(ntot))
  allocate(Psi4int(ntot))
  allocate(gxxint(ntot))
  allocate(gxyint(ntot))
  allocate(gxzint(ntot))
  allocate(gyyint(ntot))
  allocate(gyzint(ntot))
  allocate(gzzint(ntot))
  allocate(Kint(ntot))
  allocate(Kxxint(ntot))
  allocate(Kxyint(ntot))
  allocate(Kxzint(ntot))
  allocate(Kyyint(ntot))
  allocate(Kyzint(ntot))
  allocate(Kzzint(ntot))
  allocate(RKt(ntot))
  allocate(RKp(ntot))
  allocate(qtt(ntot))
  allocate(qtp(ntot))
  allocate(qpp(ntot))
  allocate(detq(ntot))
  allocate(quptt(ntot))
  allocate(quptp(ntot))
  allocate(quppp(ntot))
  allocate(Gam_ttt(ntot))
  allocate(Gam_ttp(ntot))
  allocate(Gam_tpp(ntot))
  allocate(Gam_ptt(ntot))
  allocate(Gam_ptp(ntot))
  allocate(Gam_ppp(ntot))
  allocate(R_2d(ntot))
  allocate(xi_theta(ntot))
  allocate(xi_phi(ntot))
  allocate(dArea(ntot))

  PI = 3.14159265358979323846D0
  dphi = 2.d0 * PI / N_phi
  dcostheta = 1.d0 / N_theta
!  if (Symmetry==NO_SYMM) dcostheta = 2.d0/N_theta

! Get the origin of BH
  foundflag = HorizonLocalCoordinateOrigin(horizon_number,xh,yh,zh)

  do i=1,N_theta
     costheta = 1.d0 - (i - 0.5d0)*dcostheta
     sintheta = sqrt(1.d0 - costheta*costheta)
     do j=1,N_phi
	n = j + (i-1)*N_phi
        phiangle = (j - 0.5d0)*dphi
        x_ah1(n) = xh + sintheta*cos(phiangle)
        y_ah1(n) = yh + sintheta*sin(phiangle)
        z_ah1(n) = zh + costheta
     end do
  end do

  ! Find horizon radii
  foundflag = HorizonRadiusInDirection(horizon_number,ntot,x_ah1,y_ah1,z_ah1,ah_radii)

  ! Now set the points on the horizon surface 
  do i=1,N_theta
     costheta = 1.d0 - (i - 0.5d0)*dcostheta
     sintheta = sqrt(1.d0 - costheta*costheta)
     do j=1,N_phi
    	n = j + (i-1)*N_phi
        phiangle = (j - 0.5d0)*dphi
        xinterp(n) = xh + ah_radii(n)*sintheta*cos(phiangle)
        yinterp(n) = yh + ah_radii(n)*sintheta*sin(phiangle)
        zinterp(n) = zh + ah_radii(n)*costheta
     end do
  end do

  ! Perform interpolation 
  call interpolation_for_BH_local_diagnostics(cctkGH,xinterp,yinterp,zinterp, &
        Psi4int,gxxint,gxyint,gxzint,gyyint,gyzint,gzzint, &
        Kint,Kxxint,Kxyint,Kxzint,Kyyint,Kyzint,Kzzint, ntot)
  

  ! Compute the induced metric q_ij = gamma_ij - r_i r_j and 
  ! transform to spherical coordinates q_ab, where a,b = theta, phi.
  ! Calculate n_i = psi^6 \partial_i f / (df/dr), where f=0 on the AH surface.
  ! Also compute: 
  !   det(q_ab) &  q^{ab};
  !   dArea = gamma^{ij} r_i n_j so that
  !   AH area = Integrate[ dArea  dtheta dphi];
  !   RK^a = K^{ja} n_j so that 
  !   J = Integrate[ xi_a RK^a  dcostheta dphi],
  !   where xi^a is the (approximate) axial Killing vector.
  ! 
  call BH_local_diagnostics_aux(N_theta,N_phi,ntot, &
        dcostheta,dphi, ah_radii, Psi4int,&
        gxxint,gxyint,gxzint,gyyint,gyzint,gzzint, &
        Kxxint,Kxyint,Kxzint,Kyyint,Kyzint,Kzzint, &
        qtt,qtp,qpp,detq,quptt,quptp,quppp, dArea, &
        RKt,RKp, Symmetry)

  ! Compute Christoffel symbols associated with q_ab
  call compute_Gamma_abc(Gam_ttt,Gam_ttp,Gam_tpp,Gam_ptt,Gam_ptp,Gam_ppp, &
        qtt,qtp,qpp,quptt,quptp,quppp, dcostheta, dphi, N_theta, N_phi, ntot, & 
	Symmetry)
  
  ! Compute Ricci scalar associated with q_ab
  call compute_Ricci_2D(R_2d,Gam_ttt,Gam_ttp,Gam_tpp,Gam_ptt,Gam_ptp,Gam_ppp, &
        qtt,qtp,qpp,detq,quptt,quptp,quppp, & 
        dcostheta, dphi, N_theta, N_phi, ntot, Symmetry)

  ! Find the approximate axial Killing vector xi_a 
  call find_axial_Killing_vector(xi_theta,xi_phi,xi_err,Gam_ttt,Gam_ttp, & 
  	Gam_tpp,Gam_ptt,Gam_ptp,Gam_ppp, detq,quptt,quptp,quppp, R_2d, & 
	dcostheta, dphi, N_theta, N_phi, ntot, Symmetry)

  ! Perform surface integrations to compute  Mirr, J_BH and M_BH
  call BH_diagnostics_surf_integration(Mirr, M_BH, J_BH, dcostheta, dphi, &
                xi_theta,xi_phi,RKt,RKp,dArea, ntot, sym_factor)

! Release memory
  deallocate(xinterp,yinterp,zinterp,ah_radii,x_ah1,y_ah1,z_ah1)
  deallocate(Psi4int,gxxint,gxyint,gxzint,gyyint,gyzint,gzzint)
  deallocate(Kint,Kxxint,Kxyint,Kxzint,Kyyint,Kyzint,Kzzint)
  deallocate(RKt,RKp,qtt,qtp,qpp,detq,quptt,quptp,quppp)
  deallocate(Gam_ttt,Gam_ttp,Gam_tpp,Gam_ptt,Gam_ptp,Gam_ppp,R_2d)
  deallocate(xi_theta,xi_phi,dArea)

end subroutine BH_local_diagnostics


subroutine interpolation_for_BH_local_diagnostics(cctkGH, xinterp,yinterp,zinterp, & 
	Psi4int,gxxint,gxyint,gxzint,gyyint,gyzint,gzzint, & 
	Kint,Kxxint,Kxyint,Kxzint,Kyyint,Kyzint,Kzzint, ntot)
  implicit none
  DECLARE_CCTK_FUNCTIONS
  CCTK_POINTER                             :: cctkGH
  integer :: ntot,i
  real*8, dimension(ntot) :: xinterp,yinterp,zinterp
  real*8, dimension(ntot) :: Psi4int,gxxint,gxyint,gxzint,gyyint,gyzint,gzzint
  real*8, dimension(ntot) :: Kint,Kxxint,Kxyint,Kxzint,Kyyint,Kyzint,Kzzint
  real*8, parameter :: f1o3 = 1.d0/3.d0
  integer :: N_dims,interpolation_order, ierr,interp_handle
  integer :: param_table_handle,coord_system_handle
  integer, parameter :: N_input_arrays = 14, N_output_arrays = 14
  integer,dimension(N_input_arrays)  :: input_array_type_codes,input_array_varindices
  integer,dimension(N_output_arrays) :: output_array_type_codes
  character(60)                             :: options_string
  CCTK_POINTER,dimension(N_output_arrays)   :: output_array_pointers
  CCTK_POINTER, dimension(3)                :: interp_coords
!
  N_dims = 3
  interpolation_order = 2

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
  call CCTK_VarIndex (input_array_varindices(2), "bssn::gxx")
  call CCTK_VarIndex (input_array_varindices(3), "bssn::gxy")
  call CCTK_VarIndex (input_array_varindices(4), "bssn::gxz")
  call CCTK_VarIndex (input_array_varindices(5), "bssn::gyy")
  call CCTK_VarIndex (input_array_varindices(6), "bssn::gyz")
  call CCTK_VarIndex (input_array_varindices(7), "bssn::gzz")
  call CCTK_VarIndex (input_array_varindices(8), "bssn::trK")
  call CCTK_VarIndex (input_array_varindices(9), "bssn::Axx")
  call CCTK_VarIndex (input_array_varindices(10), "bssn::Axy")
  call CCTK_VarIndex (input_array_varindices(11), "bssn::Axz")
  call CCTK_VarIndex (input_array_varindices(12), "bssn::Ayy")
  call CCTK_VarIndex (input_array_varindices(13), "bssn::Ayz")
  call CCTK_VarIndex (input_array_varindices(14), "bssn::Azz")

  output_array_pointers(1) = CCTK_PointerTo(Psi4int)
  output_array_pointers(2) = CCTK_PointerTo(gxxint)
  output_array_pointers(3) = CCTK_PointerTo(gxyint)
  output_array_pointers(4) = CCTK_PointerTo(gxzint)
  output_array_pointers(5) = CCTK_PointerTo(gyyint)
  output_array_pointers(6) = CCTK_PointerTo(gyzint)
  output_array_pointers(7) = CCTK_PointerTo(gzzint)
  output_array_pointers(8) = CCTK_PointerTo(Kint)
  output_array_pointers(9) = CCTK_PointerTo(Kxxint)
  output_array_pointers(10) = CCTK_PointerTo(Kxyint)
  output_array_pointers(11) = CCTK_PointerTo(Kxzint)
  output_array_pointers(12) = CCTK_PointerTo(Kyyint)
  output_array_pointers(13) = CCTK_PointerTo(Kyzint)
  output_array_pointers(14) = CCTK_PointerTo(Kzzint)

  interp_coords(1) = CCTK_PointerTo(xinterp)
  interp_coords(2) = CCTK_PointerTo(yinterp)
  interp_coords(3) = CCTK_PointerTo(zinterp)

  ! Perform interpolation:
  call CCTK_InterpGridArrays(ierr,cctkGH,N_dims,interp_handle, &
       param_table_handle,coord_system_handle, &
       ntot,input_array_type_codes,interp_coords, &
       N_input_arrays,input_array_varindices, &
       N_output_arrays, output_array_type_codes, output_array_pointers)

  ! Convert conformal metric to physical metric
  do i=1,ntot
     Psi4int(i) = exp(4.d0*Psi4int(i))
     gxxint(i) = Psi4int(i)*gxxint(i)
     gxyint(i) = Psi4int(i)*gxyint(i)
     gxzint(i) = Psi4int(i)*gxzint(i)
     gyyint(i) = Psi4int(i)*gyyint(i)
     gyzint(i) = Psi4int(i)*gyzint(i)
     gzzint(i) = Psi4int(i)*gzzint(i)
     Kxxint(i) = Psi4int(i)*Kxxint(i) + f1o3*gxxint(i)*Kint(i)
     Kxyint(i) = Psi4int(i)*Kxyint(i) + f1o3*gxyint(i)*Kint(i)
     Kxzint(i) = Psi4int(i)*Kxzint(i) + f1o3*gxzint(i)*Kint(i)
     Kyyint(i) = Psi4int(i)*Kyyint(i) + f1o3*gyyint(i)*Kint(i)
     Kyzint(i) = Psi4int(i)*Kyzint(i) + f1o3*gyzint(i)*Kint(i)
     Kzzint(i) = Psi4int(i)*Kzzint(i) + f1o3*gzzint(i)*Kint(i)
  end do

end subroutine interpolation_for_BH_local_diagnostics


subroutine BH_local_diagnostics_aux(N_theta,N_phi,ntot, & 
        dcostheta,dphi, ah_radii, Psi4int,& 
	gxxint,gxyint,gxzint,gyyint,gyzint,gzzint, & 
        Kxxint,Kxyint,Kxzint,Kyyint,Kyzint,Kzzint, &
	qtt,qtp,qpp,detq,quptt,quptp,quppp, dArea, & 
	RKt,RKp,Symmetry)
  implicit none
  integer :: N_theta,N_phi,ntot,n,i,j, Symmetry
  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1
  real*8, dimension(ntot) :: ah_radii,RKt,RKp,Psi4int,dArea
  real*8, dimension(ntot) :: gxxint,gxyint,gxzint,gyyint,gyzint,gzzint
  real*8, dimension(ntot) :: Kxxint,Kxyint,Kxzint,Kyyint,Kyzint,Kzzint
  real*8, dimension(ntot) :: qtt,qtp,qpp,detq,quptt,quptp,quppp
  real*8 :: costheta,sintheta,dcostheta,dphi,phiangle
  real*8 :: cosphi,sinphi,dmudx,dmudy,dmudz,dphidx,dphidy,dphidz
  real*8 :: rL,r2L,Lxt,Lxp,Lyt,Lyp,Lzt,Lzp
  real*8 :: Psi6,Psi6r2_dfdr,dRdmu,dRdphi
  real*8 :: rxL,ryL,rzL, qxxL,qxyL,qxzL,qyyL,qyzL,qzzL
  real*8 :: detqL,qttL,qtpL,qppL,qupttL,quptpL,qupppL
  real*8 :: gxxL,gxyL,gxzL,gyyL,gyzL,gzzL, nxL,nyL,nzL
  real*8 :: gupxxL,gupxyL,gupxzL,gupyyL,gupyzL,gupzzL,detg
  real*8 :: KxxL,KxyL,KxzL,KyyL,KyzL,KzzL
  real*8 :: Rupx,Rupy,Rupz,nupx,nupy,nupz,RKx,RKy,RKz,RK_t,RK_p
  real*8 :: rnorm,dRdtheta
  real*8, parameter :: pi = 3.1415926535897932d0
  integer :: ind0,ind1,ind2,indm1,indm2
!
  do i=1,N_theta
     costheta = 1.d0 - (i - 0.5d0)*dcostheta
     sintheta = sqrt(1.d0 - costheta*costheta)
     do j=1,N_phi
	n = j + (i-1)*N_phi
        phiangle = (j - 0.5d0)*dphi
        cosphi = cos(phiangle)
        sinphi = sin(phiangle)
	rL = ah_radii(n)
	r2L = rL*rL
	Psi6 = ( sqrt(Psi4int(n)) )**3

        ! dR/dmu
        if (i==1) then 
   	   ind0 = j
	   ind1 = j + N_phi
	   ind2 = j + 2*N_phi
	   dRdmu = (1.5d0*ah_radii(ind0) - 2.d0*ah_radii(ind1) + 0.5d0*ah_radii(ind2)) / dcostheta
        else if (i==N_theta) then 
!	   if (Symmetry==EQUATORIAL) then 
	      ind1 = j + (i-1)*N_phi  ! Equatorial sym --> ah_radii(ind1) = ah_radii(ind0)
	      indm1 = j + (i-2)*N_phi
	      dRdmu = 0.5d0*(ah_radii(indm1)-ah_radii(ind1))/dcostheta
	   !!!else 
	   !!!   ind0 = j + (i-1)*N_phi
	   !!!   indm1 = j + (i-2)*N_phi
	   !!!   indm2 = j + (i-3)*N_phi
	   !!!   dRdmu = -(1.5d0*ah_radii(ind0) - 2.d0*ah_radii(indm1) + 0.5d0*ah_radii(indm2)) / dcostheta
	   !!!end if
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

        dmudx = -costheta*sintheta*cosphi/rL
	dmudy = -costheta*sintheta*sinphi/rL
	dmudz = sintheta*sintheta/rL
        dphidx = -sinphi/(rL*sintheta)
  	dphidy = cosphi/(rL*sintheta)

        ! Let f = sqrt[(x-xh)^2+(y-yh)^2+(z-zh)^2] - R(mu,phi), where mu=cos(theta).
        ! Compute r_i = partial f / partial x^i
	rxL = sintheta*cosphi - dRdmu*dmudx - dRdphi*dphidx
        ryL = sintheta*sinphi - dRdmu*dmudy - dRdphi*dphidy
	rzL = costheta - dRdmu*dmudz

	! Compute n_i = psi^6 \partial_i f / (df/dr)
	Psi6r2_dfdr = r2L*Psi6/((rxL*cosphi + ryL*sinphi)*sintheta + rzL*costheta)
	nxL = rxL*Psi6r2_dfdr
	nyL = ryL*Psi6r2_dfdr
	nzL = rzL*Psi6r2_dfdr

        ! Compute the inverse of the spatial metric
	gxxL = gxxint(n)
	gxyL = gxyint(n)
        gxzL = gxzint(n)
        gyyL = gyyint(n)
        gyzL = gyzint(n)
        gzzL = gzzint(n)
        detg = gxxL * gyyL * gzzL + gxyL * gyzL * gxzL + & 
		gxzL * gxyL * gyzL - gxzL * gyyL * gxzL & 
		- gxyL * gxyL * gzzL - gxxL * gyzL * gyzL
        gupxxL = ( gyyL * gzzL - gyzL * gyzL )/detg
	gupxyL = - ( gxyL * gzzL - gyzL * gxzL )/detg
	gupxzL =   ( gxyL * gyzL - gyyL * gxzL )/detg
	gupyyL =   ( gxxL * gzzL - gxzL * gxzL )/detg
	gupyzL = - ( gxxL * gyzL - gxyL * gxzL )/detg
	gupzzL =   ( gxxL * gyyL - gxyL * gxyL )/detg

        !Now normalize r_i by the metric
	rnorm = sqrt(gupxxL*rxL*rxL + 2.d0*gupxyL*rxL*ryL + & 
		2.d0*gupxzL*rxL*rzL + gupyyL*ryL*ryL + & 
		2.d0*gupyzL*ryL*rzL + gupzzL*rzL*rzL)
        rxL = rxL/rnorm
	ryL = ryL/rnorm
	rzL = rzL/rnorm

  	! Now compute the induced metric q_ij 
	qxxL = gxxL - rxL*rxL
	qxyL = gxyL - rxL*ryL
	qxzL = gxzL - rxL*rzL
	qyyL = gyyL - ryL*ryL
	qyzL = gyzL - ryL*rzL
	qzzL = gzzL - rzL*rzL

        ! Compute q_{theta theta}, q_{theta phi} and q_{phi phi}
        ! First calculate the transformation matrix 
	! L^i_a = \partial x^i/\partial x^a, [x^i=(x,y,z), x^a=(theta,phi)]
	dRdtheta = -dRdmu*sintheta
	Lxt = dRdtheta*sintheta*cosphi + rL*costheta*cosphi
	Lxp = dRdphi*sintheta*cosphi - rL*sintheta*sinphi
	Lyt = dRdtheta*sintheta*sinphi + rL*costheta*sinphi
	Lyp = dRdphi*sintheta*sinphi + rL*sintheta*cosphi
	Lzt = dRdtheta*costheta  -rL*sintheta
	Lzp = dRdphi*costheta
	qttL = Lxt*Lxt*qxxL + 2.d0*Lxt*Lyt*qxyL + 2.d0*Lxt*Lzt*qxzL + &
                Lyt*Lyt*qyyL + 2.d0*Lyt*Lzt*qyzL + Lzt*Lzt*qzzL
	qtpL = Lxt*Lxp*qxxL + (Lxt*Lyp+Lxp*Lyt)*qxyL + &
                (Lxt*Lzp+Lxp*Lzt)*qxzL + Lyt*Lyp*qyyL + &
                (Lyt*Lzp+Lyp*Lzt)*qyzL + Lzt*Lzp*qzzL
	qppL = Lxp*Lxp*qxxL + 2.d0*Lxp*Lyp*qxyL + 2.d0*Lxp*Lzp*qxzL + &
                 Lyp*Lyp*qyyL + 2.d0*Lyp*Lzp*qyzL + Lzp*Lzp*qzzL
	
	! Compute det(q) and q^{theta theta}, q^{theta phi}, q^{phi phi}
	detqL = qttL*qppL - qtpL*qtpL 
	qtt(n) = qttL
	qtp(n) = qtpL
	qpp(n) = qppL
	detq(n) = detqL
	qupttL = qppL/detqL
	quptpL = -qtpL/detqL
	qupppL = qttL/detqL
	quptt(n) = qupttL
	quptp(n) = quptpL
	quppp(n) = qupppL

	! Compute dArea = gamma^{ij} r_i n_j
        ! Note that AH area = Integrate[ dArea  dtheta dphi]
	Rupx = gupxxL*rxL + gupxyL*ryL + gupxzL*rzL
        Rupy = gupxyL*rxL + gupyyL*ryL + gupyzL*rzL
        Rupz = gupxzL*rxL + gupyzL*ryL + gupzzL*rzL
	dArea(n) = Rupx*nxL + Rupy*nyL + Rupz*nzL

	! Compute RK^a = K^{aj} n_j
	! Note that J = Integrate[ xi_a RK^a  dcostheta dphi], 
	! where xi^a is the (approximate) axial Killing vector.
	KxxL = Kxxint(n)
	KxyL = Kxyint(n)
	KxzL = Kxzint(n)
	KyyL = Kyyint(n)
	KyzL = Kyzint(n)
	KzzL = Kzzint(n)
	nupx = gupxxL*nxL + gupxyL*nyL + gupxzL*nzL
        nupy = gupxyL*nxL + gupyyL*nyL + gupyzL*nzL
        nupz = gupxzL*nxL + gupyzL*nyL + gupzzL*nzL
	RKx = nupx*KxxL + nupy*KxyL + nupz*KxzL
        RKy = nupx*KxyL + nupy*KyyL + nupz*KyzL
        RKz = nupx*KxzL + nupy*KyzL + nupz*KzzL
	RK_t = Lxt*RKx + Lyt*RKy + Lzt*RKz
	RK_p = Lxp*RKx + Lyp*RKy + Lzp*RKz
        RKt(n) = qupttL*RK_t + quptpL*RK_p
	RKp(n) = quptpL*Rk_t + qupppL*RK_p

     end do
  end do

end subroutine BH_local_diagnostics_aux


subroutine compute_Gamma_abc(Gam_ttt,Gam_ttp,Gam_tpp,Gam_ptt,Gam_ptp,Gam_ppp, & 
        qtt,qtp,qpp,quptt,quptp,quppp, &
        dcostheta, dphi, N_theta, N_phi, ntot, Symmetry)
  implicit none
  integer :: N_theta, N_phi, ntot, Symmetry
  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1
  real*8, dimension(ntot) :: Gam_ttt,Gam_ttp,Gam_tpp,Gam_ptt,Gam_ptp,Gam_ppp
  real*8, dimension(ntot) :: qtt,qtp,qpp,quptt,quptp,quppp
  integer :: n,i,j,ind0,ind1,ind2,indm1,indm2
  real*8 :: dcostheta, dphi,costheta,sintheta,phiangle,cosphi,sinphi
  real*8 :: qttt,qttp,qtpt,qtpp,qppt,qppp
  real*8 :: Gam_down_ttt,Gam_down_ttp,Gam_down_tpp
  real*8 :: Gam_down_ptt,Gam_down_ptp,Gam_down_ppp
  real*8 :: qupttL,quptpL,qupppL
!
  do i=1,N_theta
     costheta = 1.d0 - (i - 0.5d0)*dcostheta
     sintheta = sqrt(1.d0 - costheta*costheta)
     do j=1,N_phi
        n = j + (i-1)*N_phi
        phiangle = (j - 0.5d0)*dphi
        cosphi = cos(phiangle)
        sinphi = sin(phiangle)

        ! d q_{ij} / dtheta
        if (i==1) then
           ind0 = j
           ind1 = j + N_phi
           ind2 = j + 2*N_phi
           qttt = (1.5d0*qtt(ind0) - 2.d0*qtt(ind1) + 0.5d0*qtt(ind2)) / dcostheta * (-sintheta)
           qtpt = (1.5d0*qtp(ind0) - 2.d0*qtp(ind1) + 0.5d0*qtp(ind2)) / dcostheta * (-sintheta)
           qppt = (1.5d0*qpp(ind0) - 2.d0*qpp(ind1) + 0.5d0*qpp(ind2)) / dcostheta * (-sintheta)
        else if (i==N_theta) then
!	   if (Symmetry==EQUATORIAL) then 
              ind1 = j + (i-1)*N_phi  
              indm1 = j + (i-2)*N_phi
              qttt = 0.5d0*(qtt(indm1)-qtt(ind1))/dcostheta * (-sintheta)
              ! Note that q_{theta phi} is anti-symmetric across z=0 plane
	      qtpt = 0.5d0*(qtp(indm1)+qtp(ind1))/dcostheta * (-sintheta)
	      qppt = 0.5d0*(qpp(indm1)-qpp(ind1))/dcostheta * (-sintheta)
	   !!!else
	   !!!   ind0 = j + (i-1)*N_phi
	   !!!   indm1 = j + (i-2)*N_phi
	   !!!   indm2 = j + (i-3)*N_phi
	   !!!   qttt = -(1.5d0*qtt(ind0) - 2.d0*qtt(indm1) + 0.5d0*qtt(indm2)) / dcostheta * (-sintheta)
           !!!   qtpt = -(1.5d0*qtp(ind0) - 2.d0*qtp(indm1) + 0.5d0*qtp(indm2)) / dcostheta * (-sintheta)
           !!!   qppt = -(1.5d0*qpp(ind0) - 2.d0*qpp(indm1) + 0.5d0*qpp(indm2)) / dcostheta * (-sintheta)
	   !!!end if
        else
           ind1 = j + i*N_phi
           indm1 = j + (i-2)*N_phi
           qttt = 0.5d0*(qtt(indm1)-qtt(ind1))/dcostheta * (-sintheta)
           qtpt = 0.5d0*(qtp(indm1)-qtp(ind1))/dcostheta * (-sintheta)
           qppt = 0.5d0*(qpp(indm1)-qpp(ind1))/dcostheta * (-sintheta)
        end if

	! d q_{ij} / dphi
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
        qttp = 0.5d0*(qtt(ind1)-qtt(indm1))/dphi
	qtpp = 0.5d0*(qtp(ind1)-qtp(indm1))/dphi
        qppp = 0.5d0*(qpp(ind1)-qpp(indm1))/dphi

        ! Compute Gamma_{ijk} 
	Gam_down_ttt = 0.5d0*qttt
	Gam_down_ttp = 0.5d0*qttp
	Gam_down_tpp = qtpp - 0.5d0*qppt
	Gam_down_ptt = qtpt - 0.5d0*qttp
	Gam_down_ptp = 0.5d0*qppt
	Gam_down_ppp = 0.5d0*qppp

	! Compute Gamma^i_{jk}
	qupttL = quptt(n)
	quptpL = quptp(n)
	qupppL = quppp(n)
	Gam_ttt(n) = qupttL*Gam_down_ttt + quptpL*Gam_down_ptt
	Gam_ttp(n) = qupttL*Gam_down_ttp + quptpL*Gam_down_ptp
	Gam_tpp(n) = qupttL*Gam_down_tpp + quptpL*Gam_down_ppp
	Gam_ptt(n) = quptpL*Gam_down_ttt + qupppL*Gam_down_ptt
        Gam_ptp(n) = quptpL*Gam_down_ttp + qupppL*Gam_down_ptp
	Gam_ppp(n) = quptpL*Gam_down_tpp + qupppL*Gam_down_ppp

     end do
  end do

end subroutine compute_Gamma_abc

subroutine compute_Ricci_2D(R_2d,Gam_ttt,Gam_ttp,Gam_tpp,Gam_ptt,Gam_ptp,Gam_ppp, &
        qtt,qtp,qpp,detq,quptt,quptp,quppp, dcostheta, dphi, N_theta, N_phi, ntot, Symmetry)
  implicit none
  integer :: N_theta, N_phi, ntot,Symmetry
  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1
  real*8, dimension(ntot) :: R_2d,Gam_ttt,Gam_ttp,Gam_tpp,Gam_ptt,Gam_ptp,Gam_ppp
  real*8, dimension(ntot) :: qtt,qtp,qpp,quptt,quptp,quppp,detq
  integer :: n,i,j,ind0,ind1,ind2,indm1,indm2
  real*8 :: dcostheta, dphi,costheta,sintheta,phiangle,cosphi,sinphi
  real*8 :: Gam_tppt, Gam_ttpp, Gam_pppt, Gam_ptpp
  real*8 :: detqL,qttL,qtpL,qppL,qupttL,quptpL,qupppL
  real*8 :: Gam_tttL,Gam_ttpL,Gam_tppL
  real*8 :: Gam_pttL,Gam_ptpL,Gam_pppL
  real*8 :: Rt_ptp, Rp_ptp,Rtptp
!
  do i=1,N_theta
     costheta = 1.d0 - (i - 0.5d0)*dcostheta
     sintheta = sqrt(1.d0 - costheta*costheta)
     do j=1,N_phi
        n = j + (i-1)*N_phi
        phiangle = (j - 0.5d0)*dphi
        cosphi = cos(phiangle)
        sinphi = sin(phiangle)

        ! Compute: partial_theta Gamma^{theta}_{phi phi},
	!          partial_theta Gamma^{phi}_{phi phi}.
        if (i==1) then
           ind0 = j
           ind1 = j + N_phi
           ind2 = j + 2*N_phi
           Gam_tppt = (1.5d0*Gam_tpp(ind0) - 2.d0*Gam_tpp(ind1) + & 
		0.5d0*Gam_tpp(ind2)) / dcostheta * (-sintheta)
           Gam_pppt = (1.5d0*Gam_ppp(ind0) - 2.d0*Gam_ppp(ind1) + & 
		0.5d0*Gam_ppp(ind2)) / dcostheta * (-sintheta)

        else if (i==N_theta) then
!	   if (Symmetry==EQUATORIAL) then 
              ind1 = j + (i-1)*N_phi
              indm1 = j + (i-2)*N_phi
	      ! Note that Gamma^{theta}_{phi phi} is anti-symmetric across z=0 plane
              Gam_tppt = 0.5d0*(Gam_tpp(indm1)+Gam_tpp(ind1))/dcostheta * (-sintheta)
              Gam_pppt = 0.5d0*(Gam_ppp(indm1)-Gam_ppp(ind1))/dcostheta * (-sintheta)
	   !!!else
	   !!!   ind0 = j + (i-1)*N_phi
	   !!!   indm1 = j + (i-2)*N_phi
	   !!!   indm2 = j + (i-3)*N_phi
	   !!!   Gam_ttpt = -(1.5d0*Gam_ttp(ind0) - 2.d0*Gam_ttp(indm1) + &
           !!!     0.5d0*Gam_ttp(indm2)) / dcostheta * (-sintheta)
           !!!   Gam_ptpt = -(1.5d0*Gam_ptp(ind0) - 2.d0*Gam_ptp(indm1) + &
           !!!     0.5d0*Gam_ptp(indm2)) / dcostheta * (-sintheta)
	   !!!end if
        else
           ind1 = j + i*N_phi
           indm1 = j + (i-2)*N_phi
           Gam_tppt = 0.5d0*(Gam_tpp(indm1)-Gam_tpp(ind1))/dcostheta * (-sintheta)
           Gam_pppt = 0.5d0*(Gam_ppp(indm1)-Gam_ppp(ind1))/dcostheta * (-sintheta)
        end if

	! Compute: partial_phi Gamma^{theta}_{theta phi},
        !          partial_phi Gamma^{phi}_{theta phi}.
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
        Gam_ttpp = 0.5d0*(Gam_ttp(ind1)-Gam_ttp(indm1))/dphi
        Gam_ptpp = 0.5d0*(Gam_ptp(ind1)-Gam_ptp(indm1))/dphi

	! Compute R^{theta}_{phi theta phi} and R^{phi}_{phi theta phi}
	Gam_tttL = Gam_ttt(n)
	Gam_ttpL = Gam_ttp(n)
     	Gam_tppL = Gam_tpp(n)
	Gam_pttL = Gam_ptt(n)
        Gam_ptpL = Gam_ptp(n)
        Gam_pppL = Gam_ppp(n)
	qttL = qtt(n)
	qtpL = qtp(n)
	detqL = detq(n)
	Rt_ptp = Gam_tppt - Gam_ttpp + Gam_tttL*Gam_tppL + Gam_ttpL*Gam_pppL & 
		- Gam_ttpL*Gam_ttpL - Gam_tppL*Gam_ptpL
	Rp_ptp = Gam_pppt - Gam_ptpp + Gam_pttL*Gam_tppL - Gam_ptpL*Gam_ttpL

	! Compute R_{theta phi theta phi} and Ricci scalar (R_2d)
	Rtptp = qttL*Rt_ptp + qtpL*Rp_ptp
	R_2d(n) = 2.d0*Rtptp/detqL

     end do 
  end do

end subroutine compute_Ricci_2D

subroutine find_axial_Killing_vector(xi_theta,xi_phi,xi_err,Gam_ttt,Gam_ttp, &
        Gam_tpp,Gam_ptt,Gam_ptp,Gam_ppp, detq,quptt,quptp,quppp, R_2d, &
        dcostheta, dphi, N_theta, N_phi, ntot, Symmetry)
  implicit none
  integer :: N_theta, N_phi, ntot, Symmetry,i,j, ind, int_mu
  real*8 :: dcostheta, dphi, mu, h, L_mid, xi_err
  real*8, dimension(ntot) :: xi_theta,xi_phi
  real*8, dimension(ntot) :: Gam_ttt,Gam_ttp,Gam_tpp,Gam_ptt,Gam_ptp,Gam_ppp
  real*8, dimension(ntot) :: detq,quptt,quptp,quppp, R_2d
  real*8, dimension(3) :: xi,xiout
  integer :: ind_mid
  real*8 :: xit0,xip0,qupttL,quptpL,qupppL,xi_norm2,dxit,dxip,dxi_norm2,phi
  external dxi_L_dmu, dxi_L_dphi
! 
! Find the approximate axial Killing vector at the equator, assuming equatorial 
! symmetry
! 
  call compute_axial_Killing_vector_eq(xi_phi,quppp,dphi,N_theta,N_phi, ntot)

! Now Killing transport xi_a out of the equator along lines of constant phi
  int_mu = 1
  ind_mid = 1 + (N_theta/2 - 1)*N_phi

  do j=1,N_phi

    ! Killing transport from mu =0 to mu = dcostheta/2
     i = N_theta
     ind = j + (i-1)*N_phi
    ! xi(1) = xi_theta, xi(2) = xi_phi, xi(3) = L
     xi(1) = 0.d0
     xi(2) = xi_phi(ind)
     xi(3) = 0.d0
     mu = 0.d0
     h = 0.5d0*dcostheta
     call BH_local_diagnostics_rk2(xi,mu,h,xiout,dxi_L_dmu, &
                Gam_ttt,Gam_ttp,Gam_tpp,Gam_ptt,Gam_ptp,Gam_ppp, &
                detq,quptt,quptp,quppp, R_2d, dcostheta, dphi, &
                N_theta, N_phi, ntot, i,j, int_mu, Symmetry)    

     mu = h 
     xi_theta(ind) = xiout(1)
     xi_phi(ind) = xiout(2)
     h = dcostheta

    ! Killing transport to other mu's
     do i=N_theta-1,1,-1
        xi(1) = xiout(1)
	xi(2) = xiout(2)
	xi(3) = xiout(3)
	call BH_local_diagnostics_rk2(xi,mu,h,xiout,dxi_L_dmu, &
                Gam_ttt,Gam_ttp,Gam_tpp,Gam_ptt,Gam_ptp,Gam_ppp, &
                detq,quptt,quptp,quppp, R_2d, dcostheta, dphi, &
                N_theta, N_phi, ntot, i,j, int_mu, Symmetry)
	 
	mu = mu + h
	ind = j + (i-1)*N_phi
	xi_theta(ind) = xiout(1)
	xi_phi(ind) = xiout(2)
        ! Used for checking accuracy later
        if (ind==ind_mid) L_mid = xiout(3)
     end do
     
  end do

! Finally, estimate the accuracy of the approximate Killing vector 
! by integrating the Killing transport equation along the closed 
! curve mu = 1 - (N_theta/2 - 0.5)*dcostheta
!
  int_mu = 0
  i = N_theta/2
  xit0 = xi_theta(ind_mid)
  xip0 = xi_phi(ind_mid)
  xi(1) = xit0 
  xi(2) = xip0
  xi(3) = L_mid
  h = dphi
  phi = 0.5d0*dphi
  do j=1,N_phi
     call BH_local_diagnostics_rk2(xi,phi,h,xiout,dxi_L_dphi, &
                Gam_ttt,Gam_ttp,Gam_tpp,Gam_ptt,Gam_ptp,Gam_ppp, &
                detq,quptt,quptp,quppp, R_2d, dcostheta, dphi, &
                N_theta, N_phi, ntot, i,j, int_mu, Symmetry)
     xi(1) = xiout(1)
     xi(2) = xiout(2)
     xi(3) = xiout(3)
     phi = phi + dphi
  end do
  dxit = xi(1) - xit0
  dxip = xi(2) - xip0
  qupttL = quptt(ind_mid)
  quptpL = quptp(ind_mid)
  qupppL = quppp(ind_mid)
  xi_norm2 = qupttL*xit0*xit0 + 2.d0*quptpL*xit0*xip0 + qupppL*xip0*xip0
  dxi_norm2 = qupttL*dxit*dxit + 2.d0*quptpL*dxit*dxip + qupppL*dxip*dxip
  xi_err = sqrt(dxi_norm2/xi_norm2)

end subroutine find_axial_Killing_vector

! Compute the approximate axial Killing vvector at the equator in 
! equatorial symmetry: xi_theta = 0, 
!		       xi_phi = C sqrt(q_{phi phi}), with 
!		       C = Integrate[ sqrt(q_{phi phi}) dphi]/(2 Pi)
!
subroutine compute_axial_Killing_vector_eq(xi_phi,quppp,dphi,N_theta,N_phi,ntot)
  implicit none
  integer :: N_theta,N_phi, ntot, j,ind0,ind1
  real*8, dimension(ntot) :: xi_phi,quppp
  real*8 :: sqppL,C,dphi
  real*8, parameter :: f1o2pi=0.15915494309189533577d0   ! 1/(2 Pi)
!
  C = 0.d0

  do j=1,N_phi
     ! Interpolate q_{phi phi} at the equator (mu = 0)
     !  to 4th order and then take the square root. 
     ! Note that q_{phi phi} = 1/q^{phi phi} at the equator 
     !  in equatorial symmetry
     ind0 = j + (N_theta-1)*N_phi
     ind1 = j + (N_theta-2)*N_phi
     sqppL = sqrt(0.125d0*(9.d0/quppp(ind0) - 1.d0/quppp(ind1)))

     xi_phi(ind0) = sqppL 
     C = C + sqppL
  end do

  C = C*f1o2pi*dphi

  ! Normalize xi_phi
  do j=1,N_phi
     ind0 = j + (N_theta-1)*N_phi
     xi_phi(ind0) = C*xi_phi(ind0)
  end do
  
end subroutine compute_axial_Killing_vector_eq


! Provide: partial xi_theta / partial phi, 
!          partial xi_phi / partial phi, 
! 	   partial L / partial phi
!  In the following xi(1) = xi_theta, xi(2) = xi_phi, xi(3) = L; 
!		    dxi(j) = partial xi(j) / partial phi . 
!  It is assumed that theta is fixed (fixed i0) during integration.
!
subroutine dxi_L_dphi(xi,dxi,Gam_ttt,Gam_ttp, &
        Gam_tpp,Gam_ptt,Gam_ptp,Gam_ppp, detq,quptt,quptp,quppp, R_2d, &
        dcostheta, dphi, mu, phi, N_theta, N_phi, ntot, i0,j0, Symmetry)
  implicit none
  integer :: N_theta, N_phi, ntot, i0,j0
  real*8 :: dcostheta, dphi, mu, phi
  real*8, dimension(3) :: xi,dxi
  real*8, dimension(ntot) :: Gam_ttt,Gam_ttp,Gam_tpp,Gam_ptt,Gam_ptp,Gam_ppp
  real*8, dimension(ntot) :: detq,quptt,quptp,quppp, R_2d
  integer, parameter :: m = 4 ! interpolation order
  real*8, dimension(m) :: phin,ttpn,ptpn,tppn,pppn,detqn,qupttn,quptpn,R_2dn
  integer :: jm2,n,ind_phi,ind,Symmetry
  real*8 :: phi0,djm2,tmp
  real*8 :: Gam_ttpL,Gam_ptpL,Gam_tppL,Gam_pppL,sdetqL,qupttL,quptpL,R_2dL
!
  djm2 = phi/dphi + 0.5d0 - dint(m/2.d0+0.01d0)
  jm2 = floor(djm2+1.d-10)
  phi0 = phi - (dint(m/2.d0+0.01d0) + (djm2-jm2) )*dphi

  ! Interpolation
  do n=1,m
     phi0 = phi0 + dphi
     ind_phi = jm2+n
     if (ind_phi .lt. 1) ind_phi = ind_phi + N_phi
     if (ind_phi .gt. N_phi) ind_phi = ind_phi - N_phi
     ind = ind_phi + (i0-1)*N_phi
     phin(n) = phi0 
     ttpn(n) = Gam_ttp(ind)
     ptpn(n) = Gam_ptp(ind)
     tppn(n) = Gam_tpp(ind)
     pppn(n) = Gam_ppp(ind)
     detqn(n) = detq(ind)
     qupttn(n) = quptt(ind)
     quptpn(n) = quptp(ind)
     R_2dn(n) = R_2d(ind)
  end do
  call polint(phin,ttpn,m,phi,Gam_ttpL,tmp)
  call polint(phin,ptpn,m,phi,Gam_ptpL,tmp)
  call polint(phin,tppn,m,phi,Gam_tppL,tmp)
  call polint(phin,pppn,m,phi,Gam_pppL,tmp)
  call polint(phin,detqn,m,phi,sdetqL,tmp)
  sdetqL = sqrt(sdetqL)
  call polint(phin,qupttn,m,phi,qupttL,tmp)
  call polint(phin,quptpn,m,phi,quptpL,tmp)
  call polint(phin,R_2dn,m,phi,R_2dL,tmp)

  ! Now compute partial derivatives
  dxi(1) = Gam_ttpL*xi(1) + Gam_ptpL*xi(2) - sdetqL*xi(3)
  dxi(2) = Gam_tppL*xi(1) + Gam_pppL*xi(2) 
  dxi(3) = 0.5d0*R_2dL*sdetqL*(qupttL*xi(1) + quptpL*xi(2))

end subroutine dxi_L_dphi

! Provide: partial xi_theta / partial theta, 
!          partial xi_phi / partial theta, 
! 	   partial L / partial theta
!  In the following xi(1) = xi_theta, xi(2) = xi_phi, xi(3) = L; 
!		    dxi(j) = partial xi(j) / partial theta . 
!  It is assumed that phi is fixed (fixed j0) during integration.
!
subroutine dxi_L_dmu(xi,dxi,Gam_ttt,Gam_ttp, &
        Gam_tpp,Gam_ptt,Gam_ptp,Gam_ppp, detq,quptt,quptp,quppp, R_2d, &
        dcostheta, dphi, mu, phi, N_theta, N_phi, ntot, i0,j0, Symmetry)
  implicit none
  integer :: N_theta, N_phi, ntot, i0,j0, Symmetry
  real*8 :: dcostheta, dphi, mu, phi, one_minus_mu
  real*8, dimension(3) :: xi,dxi
  real*8, dimension(ntot) :: Gam_ttt,Gam_ttp,Gam_tpp,Gam_ptt,Gam_ptp,Gam_ppp
  real*8, dimension(ntot) :: detq,quptt,quptp,quppp, R_2d
  integer, parameter :: m = 4 ! interpolation order
  real*8, dimension(m) :: one_minus_mu_n,tttn,pttn,ttpn,ptpn,detqn,quptpn,qupppn,R_2dn
  integer :: jm2,n,ind,jm2pn
  real*8 :: djm2,tmp,fac,one_minus_mu0
  real*8 :: Gam_tttL,Gam_pttL,Gam_ttpL,Gam_ptpL,sdetqL,quptpL,qupppL,R_2dL
!
  one_minus_mu = 1.d0-mu
  djm2 = one_minus_mu/dcostheta + 0.5d0 - dint(m/2.d0+0.01d0)
  jm2 = int(djm2+1.d-10)
  jm2 = max(jm2,0)
  one_minus_mu0 = (jm2-0.5d0)*dcostheta

  ! Interpolation
  do n=1,m
     one_minus_mu0 = one_minus_mu0 + dcostheta
     one_minus_mu_n(n) = one_minus_mu0
     jm2pn = jm2+n
     ! Use equatorial symmetry 
     if (jm2pn .le. N_theta) then 
        ind = j0 + (jm2pn-1)*N_phi
        tttn(n) = Gam_ttt(ind)
        pttn(n) = Gam_ptt(ind)
        ttpn(n) = Gam_ttp(ind)
        ptpn(n) = Gam_ptp(ind)
        detqn(n) = detq(ind)
        quptpn(n) = quptp(ind)
        qupppn(n) = quppp(ind)
        R_2dn(n) = R_2d(ind)
     else
	ind = j0 + (2*N_theta - jm2pn)*N_phi
	tttn(n) = -Gam_ttt(ind)
	pttn(n) = Gam_ptt(ind)
	ttpn(n) = Gam_ttp(ind)
	ptpn(n) = -Gam_ptp(ind)
	detqn(n) = detq(ind)
	quptpn(n) = -quptp(ind)
	qupppn(n) = quppp(ind)
	R_2dn(n) = R_2d(ind)
     end if
  end do
  call polint(one_minus_mu_n,tttn,m,one_minus_mu,Gam_tttL,tmp)
  call polint(one_minus_mu_n,pttn,m,one_minus_mu,Gam_pttL,tmp)
  call polint(one_minus_mu_n,ttpn,m,one_minus_mu,Gam_ttpL,tmp)
  call polint(one_minus_mu_n,ptpn,m,one_minus_mu,Gam_ptpL,tmp)
  call polint(one_minus_mu_n,detqn,m,one_minus_mu,sdetqL,tmp)
  sdetqL = sqrt(sdetqL)
  call polint(one_minus_mu_n,quptpn,m,one_minus_mu,quptpL,tmp)
  call polint(one_minus_mu_n,qupppn,m,one_minus_mu,qupppL,tmp)
  call polint(one_minus_mu_n,R_2dn,m,one_minus_mu,R_2dL,tmp)

  ! Now compute partial derivatives
  fac = -1.d0/sqrt(1.d0-mu*mu)
  dxi(1) = (Gam_tttL*xi(1) + Gam_pttL*xi(2) )*fac
  dxi(2) = (Gam_ttpL*xi(1) + Gam_ptpL*xi(2)  + sdetqL*xi(3))*fac
  dxi(3) = (-0.5d0*R_2dL*sdetqL*(quptpL*xi(1) + qupppL*xi(2)))*fac

end subroutine dxi_L_dmu

! This is the standard 2th order Runge-Kutta algorithm.
      subroutine BH_local_diagnostics_rk2(xi,x,h,xiout,derivs, & 
                Gam_ttt,Gam_ttp,Gam_tpp,Gam_ptt,Gam_ptp,Gam_ppp, & 
		detq,quptt,quptp,quppp, R_2d, dcostheta, dphi, & 
                N_theta, N_phi, ntot, i0,j0, int_mu, Symmetry)
      implicit none
      integer :: int_mu,i0,j0,ntot,N_theta,N_phi, Symmetry
      REAL*8, dimension(3) :: xi,dxi,k1,k2, xi_mid, xiout
      real*8, dimension(ntot) :: Gam_ttt,Gam_ttp,Gam_tpp,Gam_ptt,Gam_ptp,Gam_ppp
      real*8, dimension(ntot) :: detq,quptt,quptp,quppp, R_2d
      real*8 :: dcostheta, dphi, mu, phi, x, h, x_mid
      EXTERNAL derivs
!
      ! If int_mu = 1, integrate in the mu direction (mu = x, h = dcostheta). 
      ! If int_mu = 0, integrate in the phi direction (phi = x, h = dphi).
      mu = int_mu*x
      phi = (1-int_mu)*x
      call derivs(xi,dxi,Gam_ttt,Gam_ttp, &
        Gam_tpp,Gam_ptt,Gam_ptp,Gam_ppp, detq,quptt,quptp,quppp, R_2d, &
        dcostheta, dphi, mu, phi, N_theta, N_phi, ntot, i0,j0, Symmetry)

      x_mid = x + 0.5d0*h
      xi_mid = xi + 0.5d0*h*dxi 
      mu = int_mu*x_mid
      phi = (1-int_mu)*x_mid
      call derivs(xi_mid,dxi,Gam_ttt,Gam_ttp, &
        Gam_tpp,Gam_ptt,Gam_ptp,Gam_ppp, detq,quptt,quptp,quppp, R_2d, &
        dcostheta, dphi, mu, phi, N_theta, N_phi, ntot, i0,j0, Symmetry)

      xiout = xi + h*dxi

      end subroutine BH_local_diagnostics_rk2

subroutine BH_diagnostics_surf_integration(Mirr, M_BH, J_BH, dcostheta, dphi, &
                xi_theta,xi_phi,RKt,RKp,dArea, &
                ntot, sym_factor)
  implicit none
  integer :: ntot,n
  real*8 :: Mirr, M_BH, J_BH, dcostheta, dphi, sym_factor
  real*8, dimension(ntot) :: xi_theta,xi_phi,RKt,RKp,dArea
  real*8, parameter :: pi = 3.1415926535897932d0
  real*8 :: fac
!
  fac = sym_factor*dcostheta*dphi
  Mirr = 0.d0
  J_BH = 0.d0
  do n=1,ntot
     J_BH = J_BH + RKt(n)*xi_theta(n) + RKp(n)*xi_phi(n)
     Mirr = Mirr + dArea(n)
  end do

  Mirr = 0.25d0*sqrt(Mirr*fac/pi)
  J_BH = 0.125d0*J_BH*fac/pi
  M_BH = sqrt( Mirr**2 + (0.5d0*J_BH/Mirr)**2 )

end subroutine BH_diagnostics_surf_integration
