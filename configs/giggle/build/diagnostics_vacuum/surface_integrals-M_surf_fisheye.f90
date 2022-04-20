!-----------------------------------------------------------------------------
! Calculate ADM mass via surface integral
!-----------------------------------------------------------------------------
subroutine M_surf_integral(cctkGH,output_integral)
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
  
  REAL*8  excision_radius
  REAL*8  run_time
  INTEGER*4 Symmetry
  INTEGER*4 bssn_enable
  INTEGER*4 cowling_enable
  INTEGER*4 excision_enable
  INTEGER*4 fisheye_enable
  INTEGER*4 iter_count
  INTEGER*4 number_of_mol_ministeps
  INTEGER*4 rot_metric
  INTEGER*4 trA_detg_enforce
  COMMON /cctk_params_global/excision_radius,run_time,Symmetry,bssn_enable,cowling_enable,excision_enable,fisheye_enable,iter_count&
  &,number_of_mol_ministeps,rot_metric,trA_detg_enforce
  REAL*8  BH_Vol_Excise_Radius
  REAL*8  dcostheta
  REAL*8  ddrbddr
  REAL*8  dphi
  REAL*8  drbdr
  REAL*8  inner_volInt_radius
  REAL*8  rbr
  REAL*8  rsurf2
  REAL*8  surf_radius
  INTEGER*4 Compute_VolIntegrands_Every
  INTEGER*4 N_phi
  INTEGER*4 N_theta
  INTEGER*4 WhichIntegral
  INTEGER*4 arithsurf
  INTEGER*4 enable_Jz_constraint
  INTEGER*4 enable_M_constraint
  INTEGER*4 enable_P_constraint
  INTEGER*4 nsurf
  INTEGER*4 ntot
  INTEGER*4 num_BHs
  INTEGER*4 numphi
  INTEGER*4 numtheta
  INTEGER*4 scaledsurf
  INTEGER*4 sym_factor
  COMMON /diagnostics_vacuumrest/BH_Vol_Excise_Radius,dcostheta,ddrbddr,dphi,drbdr,inner_volInt_radius,rbr,rsurf2,surf_radius,Compu&
  &te_VolIntegrands_Every,N_phi,N_theta,WhichIntegral,arithsurf,enable_Jz_constraint,enable_M_constraint,enable_P_constraint,nsurf,&
  &ntot,num_BHs,numphi,numtheta,scaledsurf,sym_factor
  REAL*8  CCTKH0
  REAL*8  CCTKH1
  REAL*8  CCTKH2
  REAL*8  CCTKH3
  REAL*8  CCTKH13
  REAL*8  CCTKH14
  REAL*8  CCTKH15
  REAL*8  CCTKH16
  REAL*8  CCTKH17
  REAL*8  CCTKH18
  REAL*8  CCTKH19
  REAL*8  CCTKH20
  integer*8  bitant_plane
  integer*8  domain
  integer*8  CCTKH4
  integer*8  CCTKH5
  integer*8  CCTKH12
  INTEGER*4 CCTKH6
  INTEGER*4 CCTKH7
  INTEGER*4 CCTKH8
  INTEGER*4 CCTKH9
  INTEGER*4 CCTKH10
  INTEGER*4 CCTKH11
  COMMON /GRIDrest/CCTKH0,CCTKH1,CCTKH2,CCTKH3,CCTKH13,CCTKH14,CCTKH15,CCTKH16,CCTKH17,CCTKH18,CCTKH19,CCTKH20,bitant_plane,domain,&
  &CCTKH4,CCTKH5,CCTKH12,CCTKH6,CCTKH7,CCTKH8,CCTKH9,CCTKH10,CCTKH11
  REAL*8  CCTKH23
  REAL*8  CCTKH24
  REAL*8  CCTKH25
  integer*8  slicing_type
  INTEGER*4 CCTKH21
  INTEGER*4 CCTKH22
  INTEGER*4 CCTKH26
  INTEGER*4 CCTKH27
  INTEGER*4 CCTKH28
  INTEGER*4 CCTKH29
  INTEGER*4 CCTKH30
  COMMON /LAPSErest/CCTKH23,CCTKH24,CCTKH25,slicing_type,CCTKH21,CCTKH22,CCTKH26,CCTKH27,CCTKH28,CCTKH29,CCTKH30
  INTEGER*4 Spatial_Gauge
  INTEGER*4 CCTKH31
  INTEGER*4 CCTKH32
  COMMON /SHIFTrest/Spatial_Gauge,CCTKH31,CCTKH32
  
  ! Variables in the function call
  integer*8                             :: cctkGH
  real*8                                   :: xcenter,ycenter,zcenter
  real*8                                   :: output_integral
  ! Variables needed for interpolation
  integer*8, dimension(3)               :: interp_coords
  character(60)                            :: options_string
  integer                                  :: interpolation_order,nchars
  integer                                  :: ierr,N_dims,interp_handle,param_table_handle,coord_system_handle,N_interp_points
  integer                                  :: N_input_arrays,N_output_arrays
  real*8, dimension(N_theta*N_phi)         :: xinterp,yinterp,zinterp
  integer,dimension(13)                    :: input_array_type_codes,input_array_varindices,output_array_type_codes
  integer*8,dimension(13)               :: output_array_pointers
  ! Output arrays, dummy indices, parameters
  real*8,dimension(N_theta*N_phi)          :: Psiint,gupxxint,gupxyint,gupxzint,gupyyint,gupyzint,gupzzint
  real*8,dimension(N_theta*N_phi)          :: phixint,phiyint,phizint,Gamxint,Gamyint,Gamzint
  real*8                                   :: phiangle,costheta,sintheta,PI,nxp,nyp,nzp
  real*8                                   :: rphys,fac1,fac2,f4,df4,df4x,df4y,df4z,nx,ny,nz
  real*8                                   :: jxx,jxy,jxz,jyy,jyz,jzz,jxxx,jyyy,jzzz,jxyz
  real*8                                   :: jxxy,jxyy,jxxz,jxzz,jyyz,jyzz,gamx,gamy,gamz
  integer                                  :: i,j,n
  PI = 3.14159265358979323846D0
  N_dims = 3
  interpolation_order = 2
  N_interp_points = N_theta*N_phi
  interp_handle = -1
!  call CCTK_InterpHandle (interp_handle, uniform cartesian)
  call CCTK_InterpHandle (interp_handle, "Lagrange polynomial interpolation")
  if (interp_handle .lt. 0) then
     call CCTK_Warn(0,46,"surface_integrals-M_surf_fisheye.F90","diagnostics_vacuum","Cannot get handle for interpolation ! Forgot &
  &to activate an implementation providing interpolation operators ??")
  endif
  param_table_handle = -1
  options_string = "order = " // char(ichar('0') + interpolation_order)
  call Util_TableCreateFromString (param_table_handle, options_string)
  if (param_table_handle .lt. 0) then
     call CCTK_Warn(0,53,"surface_integrals-M_surf_fisheye.F90","diagnostics_vacuum","Cannot create parameter table for interpolato&
  &r")
  endif
  coord_system_handle = -1
  call CCTK_CoordSystemHandle (coord_system_handle, "cart3d")
  if (coord_system_handle .lt. 0) then
     call CCTK_Warn(0,59,"surface_integrals-M_surf_fisheye.F90","diagnostics_vacuum","Cannot get handle for cart3d coordinate syste&
  &m ! Forgot to activate an implementation providing coordinates ??")
  endif
  input_array_type_codes = 107
  output_array_type_codes = 107
  ! Specify interpolation input arrays, output arrays:
  N_input_arrays = 13
  N_output_arrays = 13
  call CCTK_VarIndex (input_array_varindices(1), "bssn::phi")
  call CCTK_VarIndex (input_array_varindices(2), "bssn::gupxx")
  call CCTK_VarIndex (input_array_varindices(3), "bssn::gupxy")
  call CCTK_VarIndex (input_array_varindices(4), "bssn::gupxz")
  call CCTK_VarIndex (input_array_varindices(5), "bssn::gupyy")
  call CCTK_VarIndex (input_array_varindices(6), "bssn::gupyz")
  call CCTK_VarIndex (input_array_varindices(7), "bssn::gupzz")
  call CCTK_VarIndex (input_array_varindices(8), "bssn::phix")
  call CCTK_VarIndex (input_array_varindices(9), "bssn::phiy")
  call CCTK_VarIndex (input_array_varindices(10), "bssn::phiz")
  call CCTK_VarIndex (input_array_varindices(11), "bssn::Gammax")
  call CCTK_VarIndex (input_array_varindices(12), "bssn::Gammay")
  call CCTK_VarIndex (input_array_varindices(13), "bssn::Gammaz")
  output_array_pointers(1) = CCTK_PointerTo(Psiint)
  output_array_pointers(2) = CCTK_PointerTo(gupxxint)
  output_array_pointers(3) = CCTK_PointerTo(gupxyint)
  output_array_pointers(4) = CCTK_PointerTo(gupxzint)
  output_array_pointers(5) = CCTK_PointerTo(gupyyint)
  output_array_pointers(6) = CCTK_PointerTo(gupyzint)
  output_array_pointers(7) = CCTK_PointerTo(gupzzint)
  output_array_pointers(8) = CCTK_PointerTo(phixint)
  output_array_pointers(9) = CCTK_PointerTo(phiyint)
  output_array_pointers(10) = CCTK_PointerTo(phizint)
  output_array_pointers(11) = CCTK_PointerTo(Gamxint)
  output_array_pointers(12) = CCTK_PointerTo(Gamyint)
  output_array_pointers(13) = CCTK_PointerTo(Gamzint)
  ! Set up interpolation coordinate arrays:
  n = 1
  do i=1,N_theta
     costheta = 1.D0 - (i - 0.5D0)*dcostheta
     sintheta = sqrt(1.D0 - costheta*costheta)
     do j=1,N_phi
        phiangle = (j - 0.5D0)*dphi
        if(N_phi==1) phiangle = 0.D0
        xinterp(n) = surf_radius*sintheta*cos(phiangle)
        yinterp(n) = surf_radius*sintheta*sin(phiangle)
        zinterp(n) = surf_radius*costheta
        n = n + 1
     end do
  end do
  interp_coords(1) = CCTK_PointerTo(xinterp)
  interp_coords(2) = CCTK_PointerTo(yinterp)
  interp_coords(3) = CCTK_PointerTo(zinterp)
  ! Perform interpolation:
  call CCTK_InterpGridArrays(ierr,cctkGH,N_dims,interp_handle, &
       param_table_handle,coord_system_handle, &
       N_interp_points,input_array_type_codes,interp_coords, &
       N_input_arrays,input_array_varindices, &
       N_output_arrays, output_array_type_codes, output_array_pointers)
  Psiint = exp(Psiint)
  !Were looking at an r=constant surface, and fisheye transfs. are spherically symmetric, so rphys,dR,and ddR are fixed.
  rphys = surf_radius*rbr
  fac1 = (rbr)**(1.D0/3.D0) * drbdr**(1.D0/6.D0)
  fac2 = (drbdr/rphys-1.D0/surf_radius)/3.D0 + ddrbddr/(drbdr*6.D0)
  ! we need to calculate \tilde{Gamma}^i from \bar{\tilde{Gamma}}^i, via jacobian and derivatives
  ! f4 is F in my notes, and df4*nx gives you dF/dx
  f4=fac1**(-4)
  df4=f4/3.0d0*(-4.0d0/rphys*(drbdr-rbr)-2.0d0*ddrbddr/drbdr)
!  n=ntot/50
  if(ierr.ge.0) then
     output_integral = 0.D0
     !|~~~~~> Integrate the dot product of Dphi with the surface normal
     do i=1,ntot
        !Convert phix, phiy, and phiz to physical coordinates:
        nx= xinterp(i)/surf_radius
        ny= yinterp(i)/surf_radius
        nz= zinterp(i)/surf_radius
        df4x=df4*nx
        df4y=df4*ny
        df4z=df4*nz
        phixint(i) = phixint(i) - nx*fac2
        phiyint(i) = phiyint(i) - ny*fac2
        phizint(i) = phizint(i) - nz*fac2
        !
        Jxx=rbr+nx*nx*(drbdr-rbr)
        Jxy=nx*ny*(drbdr-rbr)
        Jxz=nx*nz*(drbdr-rbr)
        Jyy=rbr+ny*ny*(drbdr-rbr)
        Jyz=ny*nz*(drbdr-rbr)
        Jzz=rbr+nz*nz*(drbdr-rbr)
        jxxx=nx**3*ddrbddr+(drbdr-rbr)/surf_radius*(3.0d0*nx-3.0d0*nx**3)
        jyyy=ny**3*ddrbddr+(drbdr-rbr)/surf_radius*(3.0d0*ny-3.0d0*ny**3)
        jzzz=nz**3*ddrbddr+(drbdr-rbr)/surf_radius*(3.0d0*nz-3.0d0*nz**3)
        jxxy=nx**2*ny*ddrbddr+(drbdr-rbr)/surf_radius*(ny-3.0d0*nx**2*ny)
        jxyy=nx*ny**2*ddrbddr+(drbdr-rbr)/surf_radius*(nx-3.0d0*nx*ny**2)
        jxxz=nx**2*nz*ddrbddr+(drbdr-rbr)/surf_radius*(nz-3.0d0*nx**2*nz)
        jxzz=nx*nz**2*ddrbddr+(drbdr-rbr)/surf_radius*(nx-3.0d0*nx*nz**2)
        jyyz=ny**2*nz*ddrbddr+(drbdr-rbr)/surf_radius*(nz-3.0d0*ny**2*nz)
        jyzz=ny*nz**2*ddrbddr+(drbdr-rbr)/surf_radius*(ny-3.0d0*ny*nz**2)
        jxyz=nx*ny*nz*(ddrbddr-3.0d0*(drbdr-rbr)/surf_radius)
        Gamx=f4*(jxx*gamxint(i)+jxy*gamyint(i)+jxz*gamzint(i))+ &
             gupxxint(i)*(0.5d0*df4x*jxx-f4*jxxx)+ &
             gupyyint(i)*(0.5d0*df4y*jxy-f4*jxyy)+ &
             gupzzint(i)*(0.5d0*df4z*jxz-f4*jxzz)+ &
             gupxyint(i)*(0.5d0*(df4x*jxy+df4y*jxx)-2.0d0*f4*jxxy)+ &
             gupxzint(i)*(0.5d0*(df4x*jxz+df4z*jxx)-2.0d0*f4*jxxz)+ &
             gupyzint(i)*(0.5d0*(df4y*jxz+df4z*jxy)-2.0d0*f4*jxyz)
        Gamy=f4*(jxy*gamxint(i)+jyy*gamyint(i)+jyz*gamzint(i))+ &
             gupxxint(i)*(0.5d0*df4x*jxy-f4*jxxy)+ &
             gupyyint(i)*(0.5d0*df4y*jyy-f4*jyyy)+ &
             gupzzint(i)*(0.5d0*df4z*jyz-f4*jyzz)+ &
             gupxyint(i)*(0.5d0*(df4x*jyy+df4y*jxy)-2.0d0*f4*jxyy)+ &
             gupxzint(i)*(0.5d0*(df4x*jyz+df4z*jxy)-2.0d0*f4*jxyz)+ &
             gupyzint(i)*(0.5d0*(df4y*jyz+df4z*jyy)-2.0d0*f4*jyyz)
        Gamz=f4*(jxz*gamxint(i)+jyz*gamyint(i)+jzz*gamzint(i))+ &
             gupxxint(i)*(0.5d0*df4x*jxz-f4*jxxz)+ &
             gupyyint(i)*(0.5d0*df4y*jyz-f4*jyyz)+ &
             gupzzint(i)*(0.5d0*df4z*jzz-f4*jzzz)+ &
             gupxyint(i)*(0.5d0*(df4x*jyz+df4y*jxz)-2.0d0*f4*jxyz)+ &
             gupxzint(i)*(0.5d0*(df4x*jzz+df4z*jxz)-2.0d0*f4*jxzz)+ &
             gupyzint(i)*(0.5d0*(df4y*jzz+df4z*jyz)-2.0d0*f4*jyzz)
        output_integral = output_integral + Psiint(i) * fac1 * &
             (nx*(gupxxint(i)*phixint(i)+gupxyint(i)*phiyint(i)+gupxzint(i)*phizint(i))+ &
             ny*(gupxyint(i)*phixint(i)+gupyyint(i)*phiyint(i)+gupyzint(i)*phizint(i))+ &
             nz*(gupxzint(i)*phixint(i)+gupyzint(i)*phiyint(i)+gupzzint(i)*phizint(i)))- &
             rbr**2/8.0d0*(nx*gamx+ny*Gamy+nz*Gamz)
     end do
     !|  Multiply by -1/2\pi * surface element.
     output_integral = output_integral * -0.5D0/PI * surf_radius * surf_radius * dphi * dcostheta*sym_factor
  else
     output_integral = -9999.D0
  end if
end subroutine M_surf_integral
