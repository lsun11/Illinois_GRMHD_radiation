subroutine surface_dens_2d(cctkGH,time,N_X,N_Z,&
     X_max,Z_max,Symmetry,myproc)
  implicit none
!  INTEGER cctk_dim&&                           INTEGER cctk_gsh(cctk_dim),cctk_lsh(cctk_dim)&&                           INTEGER cc
!  REAL*8  excision_radius&&REAL*8  run_time&&INTEGER*4 Symmetry&&INTEGER*4 bssn_enable&&INTEGER*4 cowling_enable&&INTEGER*4 excisio
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
  character                                :: filename_2d*50
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  real*8 :: time
  real*8 :: Sigma_avgL,SigmaL
  real*8 :: rr_i, dS, PI, f1ospi
  real*8 :: cosphi,sinphi,sym_factor
  real*8 :: sintheta,costheta,nn,xn,yn,phiangle,dphi,dcostheta
  integer :: i,j,k,interp_order,n,ntot
  integer :: ind0,ind1,ind2,indm1,indm2,vindex,myproc
  integer :: N_X,N_Y,N_Z, Symmetry
  real*8 :: X_max,Y_max,Z_max
  real*8 :: XL,YL,ZL,dX,dY,dZ
  real*8, allocatable, dimension(:,:)        :: pointcoords
  real*8, allocatable, dimension(:)          :: rhosint,lapm1int,phiint
  if (Symmetry .ne. EQUATORIAL) then
     write(*,*) 'Symmetry not supported in surface_density_profile'
     stop
  end if
  write(*,*) "inside 2d"
  N_Y=N_X
  Y_max=X_max
  sym_factor = 2.d0
  ntot = N_Z*N_Y*N_X
  ! allocate memory  
  allocate(pointcoords(ntot,3))
  allocate(rhosint(ntot))
  allocate(lapm1int(ntot))
  allocate(phiint(ntot))
  PI = 3.14159265358979323844D0
  dX = 2.d0*X_max / N_X
  dY = 2.d0*Y_max / N_Y
  dZ = Z_max / N_Z
  n = 1
!!!!  !$omp parallel do
  do i=1,N_X
     XL = -X_max+(i-0.5)*dX
     do j=1,N_Y
        YL = -Y_max+(j-0.5)*dY
        do k=1,N_Z
           ZL = dZ*(k-0.5)
           pointcoords(n,1) = XL
           pointcoords(n,2) = YL
           pointcoords(n,3) = ZL
           n = n + 1
        end do
     end do
  end do
 !!!!! !$omp end parallel do
  ! Interpolate the grid function
  call CCTK_VarIndex(vindex,"mhd_evolve::rho_star")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,rhosint)
  call CCTK_VarIndex(vindex,"lapse::lapm1")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,lapm1int)
  call CCTK_VarIndex(vindex,"bssn::phi")
  call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,phiint)
  !
  ! integrate along z
  !
  if (myproc.eq.0) then
  write( filename_2d, '(a14,E10.4E2,a4)' )  'surf_dens_2d_t',time,'.dat'
  open(UNIT=16,FILE=filename_2d,STATUS="REPLACE")
  write(16,*) ""
  write(16,*) '#2d surface density'
  n=1
  write(16,*) ""
  write(16,*) ""
  write(16,*) '#Time =',time
  write(16,'(t2,a6,t31,a5)') '#varpi','Sigma'
  do i=1,N_X
     XL = -X_max+(i-0.5)*dX
     do j=1,N_Y
        YL = -Y_max+(j-0.5)*dY
        SigmaL = 0.d0
        do k=1,N_Z
           ! read in data
           SigmaL = SigmaL + rhosint(n)*dZ*(lapm1int(n)+1.d0)*exp(6*phiint(n))
           n=n+1
        end do
        write(16,'(3e18.10)') XL,YL,SigmaL
     end do
  write(16,*) ""
  end do
  close(16)
 endif
 ! deallocate memory
 deallocate(pointcoords)
 deallocate(rhosint,lapm1int,phiint)
end subroutine surface_dens_2d
