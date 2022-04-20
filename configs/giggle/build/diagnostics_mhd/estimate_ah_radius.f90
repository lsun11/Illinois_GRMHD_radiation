subroutine estimate_ah_radius(horizon_number,Symmetry,found_horizon,xh,yh,zh,rad_est)
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
  
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  real*8 :: xh,yh,zh
  integer :: horizon_number, foundflag,ntot,i
  integer :: Symmetry, found_horizon
  real*8, allocatable, dimension(:)         :: ah_radii, x_ah,y_ah,z_ah
  real*8  :: rad_est
  !
  foundflag = HorizonWasFound(horizon_number)
  ! Horizon not found, 
  if (foundflag .ne. 1) then
     found_horizon = 0
     return
  end if
  found_horizon = 1
  if (Symmetry .ne. EQUATORIAL .and. Symmetry .ne. NO_SYMM) then
     write(*,*) 'Symmetry not supported in estimate_ah_radius'
     stop
  end if
  ntot = 6
  ! allocate memory
  allocate(ah_radii(ntot))
  allocate(x_ah(ntot))
  allocate(y_ah(ntot))
  allocate(z_ah(ntot))
  ! Get the origin of BH
  foundflag = HorizonLocalCoordinateOrigin(horizon_number,xh,yh,zh)
  !+X direction
  x_ah(1) = xh+1.d0
  y_ah(1) = yh
  z_ah(1) = zh
  !-X direction
  x_ah(2) = xh-1.d0
  y_ah(2) = yh
  z_ah(2) = zh
  !+Y direction
  x_ah(3) = xh
  y_ah(3) = yh+1.d0
  z_ah(3) = zh
  !-Y direction
  x_ah(4) = xh
  y_ah(4) = yh-1.d0
  z_ah(4) = zh
  !+Z direction
  x_ah(5) = xh
  y_ah(5) = yh
  z_ah(5) = zh+1.d0
  !-Z direction
  x_ah(6) = xh
  y_ah(6) = yh
  z_ah(6) = zh-1.d0
  ! Find horizon radii
  foundflag = HorizonRadiusInDirection(horizon_number,ntot,x_ah,y_ah,z_ah,ah_radii)
  rad_est = 0.d0
  do i=1,6
     rad_est=rad_est+ah_radii(i)/6.d0
  end do
  ! deallocate memory
  deallocate(ah_radii, x_ah, y_ah, z_ah)
end subroutine Estimate_ah_radius
