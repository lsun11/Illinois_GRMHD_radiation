!
!-----------------------------------------------------------------------------
!
!  Kreiss-Oliger dissipation 
!
!-----------------------------------------------------------------------------
!
subroutine kreiss_oliger(ex,X,Z,f,dlap,dT,C_ko)
   implicit none
   integer, dimension(3)                        :: ex
   real*8, dimension(ex(1),ex(2),ex(3))                :: X,Z,f,dlap
   real*8                                            :: h4,dT,fac,dX,dZ,C_ko
   integer                                        :: i,j,k,imin,imax,kmin,kmax
!
  imin = lbound(f,1)
  kmin = lbound(f,3)
  imax = ubound(f,1)
  kmax = ubound(f,3)
  dX = X(imin+1,1,1)-X(imin,1,1)
  dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
  fac=-C_ko*dX*dX*dZ*dZ/dT/16.d0
!
  j=2
  do k=kmin,kmax
     do i=imin,imax
        f(i,j,k)=f(i,j,k)+fac*dlap(i,j,k)
     end do
  end do
end subroutine kreiss_oliger
!
!-----------------------------------------------------------------------------
!
!  Kreiss-Oliger dissipation (3D)
!
!-----------------------------------------------------------------------------
!
subroutine kreiss_oliger3d(ex,X,Y,Z,f,dlap,dT,C_ko)
   implicit none
   integer, dimension(3)                        :: ex
   real*8, dimension(ex(1),ex(2),ex(3))         :: X,Y,Z,f,dlap
   real*8                                       :: h4,dT,fac,dX,dY,dZ,C_ko
   integer                                      :: i,j,k,imin,imax,kmin,kmax
   integer                                      :: jmin, jmax
!
  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)
  imax = ubound(f,1)
  jmax = ubound(f,2)
  kmax = ubound(f,3)
  dX = X(imin+1,1,1)-X(imin,1,1)
  dY = Y(1,jmin+1,1)-Y(1,jmin,1)
  dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
  fac=-C_ko/dT/16.d0*(dX*dY*dZ)**(4.d0/3.d0)
!
  f=f+fac*dlap
end subroutine kreiss_oliger3d
!-----------------------------------------------------------------------------
!
! Compute the flat space Laplacian (3D)
!
!-----------------------------------------------------------------------------
!
subroutine flat_lap3d(ex,X,Y,Z,f,lapf,SYMx,SYMy,SYMz,Symmetry)
  implicit none
  integer, dimension(3)                         :: ex
  real*8, dimension(ex(1),ex(2),ex(3))          :: X,Y,Z,f,lapf
  integer                                       :: i,j,k,imin,imax,kmin,kmax
  integer                                       :: jmin,jmax
  real*8                                        :: dX,dY,dZ,ddx,ddy,ddz
  real*8                                        :: fxx,fyy,fzz,r1,r2
  real*8                                        :: SYMx,SYMz,SYMy
  integer                                       :: Symmetry,NO_SYMM,EQUATORIAL,OCTANT,PI_SYMM,AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
!
  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)
  imax = ubound(f,1)
  jmax = ubound(f,2)
  kmax = ubound(f,3)
  dX = X(imin+1,1,1)-X(imin,1,1)
  dY = Y(1,jmin+1,1)-Y(1,jmin,1)
  dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
  if(Symmetry == OCTANT) then
     if(X(1,1,1).lt.0.D0) then
        imin = imin + 1
     end if
     if(Y(1,1,1).lt.0.D0) then
        jmin = jmin + 1
     end if
     if(Z(1,1,1).lt.0.D0) then
        kmin = kmin + 1
     end if
  else if(Symmetry == AXISYM) then
     if(X(1,1,1).lt.0.D0) then
        imin = imin + 1
     end if
     if(Z(1,1,1).lt.0.D0) then
        kmin = kmin + 1
     end if
  else if(Symmetry == EQUATORIAL) then
     if(Z(1,1,1).lt.0.D0) then
        kmin = kmin + 1
     end if
  end if
  ddx=1.d0/dX/dX
  ddy=1.d0/dY/dY
  ddz=1.d0/dZ/dZ
!
     do k=kmin,kmax
        do j=jmin,jmax
          do i=imin,imax
             !Is there a bug here in EQUATORIAL symmetry?!:
             if (i==imin) then
                fxx=(f(i+1,j,k)-2.d0*f(i,j,k)+SYMx*f(i,j,k))*ddx
             else if (i==imax) then
              r1 = sqrt(X(i,1,1)**2 + Y(1,j,1)**2 + Z(1,1,k)**2)
              r2 = sqrt((X(i,1,1)+dX)**2 + Y(1,j,1)**2 + Z(1,1,k)**2)
              fxx=(f(i-1,j,k)+(r1/r2 - 2.d0)*f(i,j,k))*ddx
             else
              fxx=(f(i+1,j,k)-2.d0*f(i,j,k)+f(i-1,j,k))*ddx
             end if
             if (j==jmin) then
                fyy=(f(i,j+1,k)-2.d0*f(i,j,k)+SYMy*f(i,j,k))*ddy
             else if (j==jmax) then
              r1 = sqrt(X(i,1,1)**2 + Y(1,j,1)**2 + Z(1,1,k)**2)
              r2 = sqrt(X(i,1,1)**2 + (Y(1,j,1)+dY)**2 + Z(1,1,k)**2)
              fyy=(f(i,j-1,k)+(r1/r2 - 2.d0)*f(i,j,k))*ddy
             else
              fyy=(f(i,j+1,k)-2.d0*f(i,j,k)+f(i,j-1,k))*ddx
             end if
           if (k==kmin) then
              fzz=(f(i,j,k+1)-2.d0*f(i,j,k)+SYMz*f(i,j,k))*ddz
           else if (k==kmax) then
              r1 = sqrt(X(i,1,1)**2 + Y(1,j,1)**2 + Z(1,1,k)**2)
              r2 = sqrt(X(i,1,1)**2 + Y(1,j,1)**2 + (Z(1,1,k)+dZ)**2)
              fzz=(f(i,j,k-1)+(r1/r2 - 2.d0)*f(i,j,k))*ddz
           else
              fzz=(f(i,j,k+1)-2.d0*f(i,j,k)+f(i,j,k-1))*ddz
           end if
           lapf(i,j,k)=fxx+fyy+fzz
        end do
      end do
     end do
end subroutine flat_lap3d
!
!-----------------------------------------------------------------------------
!
! Compute the flat space Laplacian
!
!-----------------------------------------------------------------------------
!
subroutine flat_lap(ex,X,Z,f,lapf,SYMx,SYMz,Symmetry)
  implicit none
  integer, dimension(3)                         :: ex
  real*8, dimension(ex(1),ex(2),ex(3))          :: X,Z,f,lapf
  integer                                       :: i,j,k,imin,imax,kmin,kmax
  real*8                                        :: dX,dZ,ddx,ddz
  real*8                                        :: fx,fxx,fz,fzz,r1,r2
  real*8                                        :: SYMx,SYMz
  integer                                       :: Symmetry,NO_SYMM,EQUATORIAL,OCTANT,PI_SYMM,AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
!
  imin = lbound(f,1)
  kmin = lbound(f,3)
  imax = ubound(f,1)
  kmax = ubound(f,3)
  dX = X(imin+1,1,1)-X(imin,1,1)
  dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
  ddx=1.d0/dX/dX
  ddz=1.d0/dZ/dZ
  if(Symmetry == AXISYM) then
     if(X(1,1,1).lt.0.D0) then
        imin = imin + 1
     end if
     if(Z(1,1,1).lt.0.D0) then
        kmin = kmin + 1
     end if
  end if
  j=2
!
     do k=kmin,kmax
        do i=imin,imax
           if (i==imin) then
              fx=(f(i+1,j,k)-SYMx*f(i,j,k))/2.d0/dX
              fxx=(f(i+1,j,k)-2.d0*f(i,j,k)+SYMx*f(i,j,k))*ddx
           else if (i==imax) then
              fx=(f(i,j,k)-f(i-1,j,k))/dX
!!$              fxx=(f(i-1,j,k)-f(i,j,k))*ddx
              r1 = sqrt(X(i,1,1)**2 + Z(1,1,k)**2)
              r2 = sqrt((X(i,1,1)+dX)**2 + Z(1,1,k)**2)
              fxx=(f(i-1,j,k)+(r1/r2 - 2.d0)*f(i,j,k))*ddx
           else
              fx=(f(i+1,j,k)-f(i-1,j,k))/2.d0/dX
              fxx=(f(i+1,j,k)-2.d0*f(i,j,k)+f(i-1,j,k))*ddx
           end if
           if (k==kmin) then
              fzz=(f(i,j,k+1)-2.d0*f(i,j,k)+SYMz*f(i,j,k))*ddz
           else if (k==kmax) then
!!$              fzz=(f(i,j,k-1)-f(i,j,k))*ddz
              r1 = sqrt(X(i,1,1)**2 + Z(1,1,k)**2)
              r2 = sqrt(X(i,1,1)**2 + (Z(1,1,k)+dZ)**2)
              fzz=(f(i,j,k-1)+(r1/r2 - 2.d0)*f(i,j,k))*ddz
           else
              fzz=(f(i,j,k+1)-2.d0*f(i,j,k)+f(i,j,k-1))*ddz
           end if
           lapf(i,j,k)=fxx+fx/X(i,1,1)+fzz
        end do
     end do
end subroutine flat_lap
