!-------------------------------------------
! Excise regions interior to TWO punctures.
!-------------------------------------------
subroutine vol_integrand_excise_two_punctures(ex, X, Y, Z, &
     PsiRes, Symmetry,xc1,yc1,zc1,ah_radius1,xc2,yc2,zc2,ah_radius2)
  implicit none
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: PsiRes
  real*8                                      :: xc1,yc1,zc1,ah_radius1,xc2,yc2,zc2,ah_radius2
  integer                                     :: Symmetry
!
  integer         :: imin, jmin, kmin, imax, jmax, kmax
  integer         :: i,j,k
  real*8          :: F1o8, ONE, TWO, ZERO, Psil, Psi5, h3, EIGHT, PI
  integer         :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  real*8 :: rp1,rp2
  parameter ( ZERO = 0.D0, ONE = 1.D0, F1o8 = 1.D0 / 8.D0, TWO = 2.D0 )
  parameter ( EIGHT = 8.D0 )
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3)
  parameter(AXISYM = 4)
  PI = acos(-ONE)
!
! Input translation
!
  imin = lbound(PsiRes,1)
  jmin = lbound(PsiRes,2)
  kmin = lbound(PsiRes,3)
  imax = ubound(PsiRes,1)
  jmax = ubound(PsiRes,2)
  kmax = ubound(PsiRes,3)

  if (Symmetry==AXISYM) then
     jmin=2
     jmax=2
  end if
!
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
           rp1=sqrt((x(i,j,k)-xc1)**2+(y(i,j,k)-yc1)**2+(z(i,j,k)-zc1)**2)
           rp2=sqrt((x(i,j,k)-xc2)**2+(y(i,j,k)-yc2)**2+(z(i,j,k)-zc2)**2)
           if(rp1.lt.(ah_radius1) .or. rp2.lt.(ah_radius2)) PsiRes(i,j,k)=0.D0
       end do
     end do
  end do
  
  return
end subroutine vol_integrand_excise_two_punctures

!------------------------------------------------------------------
! Excise all but two thick spherical shells surrounding punctures.
!------------------------------------------------------------------
subroutine vol_integrand_two_spherical_shells(ex, X, Y, Z, &
     PsiRes, Symmetry,xc1,yc1,zc1,radius1a,radius1b,xc2,yc2,zc2,radius2a,radius2b)
  implicit none
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: PsiRes
  real*8                                      :: xc1,yc1,zc1,radius1a,radius1b,xc2,yc2,zc2,radius2a,radius2b
  integer                                     :: Symmetry
!
  integer         :: imin, jmin, kmin, imax, jmax, kmax
  integer         :: i,j,k
  real*8          :: F1o8, ONE, TWO, ZERO, Psil, Psi5, h3, EIGHT, PI
  integer         :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  real*8 :: rp1,rp2
  parameter ( ZERO = 0.D0, ONE = 1.D0, F1o8 = 1.D0 / 8.D0, TWO = 2.D0 )
  parameter ( EIGHT = 8.D0 )
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3)
  parameter(AXISYM = 4)
  PI = acos(-ONE)
!
! Input translation
!
  imin = lbound(PsiRes,1)
  jmin = lbound(PsiRes,2)
  kmin = lbound(PsiRes,3)
  imax = ubound(PsiRes,1)
  jmax = ubound(PsiRes,2)
  kmax = ubound(PsiRes,3)

  if (Symmetry==AXISYM) then
     jmin=2
     jmax=2
  end if
!
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
           rp1=sqrt((x(i,j,k)-xc1)**2+(y(i,j,k)-yc1)**2+(z(i,j,k)-zc1)**2)
           rp2=sqrt((x(i,j,k)-xc2)**2+(y(i,j,k)-yc2)**2+(z(i,j,k)-zc2)**2)
           if((rp1.lt.radius1b .and. rp1.gt.radius1a) .or. (rp2.lt.radius2b .and. rp2.gt.radius2a)) then
              ! When we are in one of the spherical_shells, we leave psires alone!
           else
              PsiRes(i,j,k)=0.D0
           end if
       end do
     end do
  end do
  
  return
end subroutine vol_integrand_two_spherical_shells

!-----------------------------------------------------------
! Excise all but thick spherical shell around ONE puncture.
!-----------------------------------------------------------
subroutine vol_integrand_one_spherical_shell(ex, X, Y, Z, &
     PsiRes, Symmetry,xc,yc,zc,shell_radius1,shell_radius2)
  implicit none
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: PsiRes
  real*8                                      :: xc,yc,zc,shell_radius1,shell_radius2
  integer                                     :: Symmetry
!
  integer         :: imin, jmin, kmin, imax, jmax, kmax
  integer         :: i,j,k
  real*8          :: F1o8, ONE, TWO, ZERO, Psil, Psi5, h3, EIGHT, PI
  integer         :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  real*8 :: radialpoint
  parameter ( ZERO = 0.D0, ONE = 1.D0, F1o8 = 1.D0 / 8.D0, TWO = 2.D0 )
  parameter ( EIGHT = 8.D0 )
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3)
  parameter(AXISYM = 4)
  PI = acos(-ONE)
!
! Input translation
!
  imin = lbound(PsiRes,1)
  jmin = lbound(PsiRes,2)
  kmin = lbound(PsiRes,3)
  imax = ubound(PsiRes,1)
  jmax = ubound(PsiRes,2)
  kmax = ubound(PsiRes,3)

  if (Symmetry==AXISYM) then
     jmin=2
     jmax=2
  end if

  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
           radialpoint=sqrt((x(i,j,k)-xc)**2+(y(i,j,k)-yc)**2+(z(i,j,k)-zc)**2)
           if(radialpoint .lt.(shell_radius1) .or. radialpoint.gt.(shell_radius2)) PsiRes(i,j,k)=0.D0
       end do
     end do
  end do
  
  return
end subroutine vol_integrand_one_spherical_shell
