!---------------------------------------------------------------
! A subroutine to correct the scalar f
! on the y/= 0 planes in the axisymmetric case
!---------------------------------------------------------------
subroutine ezaxibc_scalar(ex,X,Y,Z,f,dm,Symmetry)
  implicit none
  integer, dimension(3)                 :: ex
  real*8, dimension(ex(1),ex(2),ex(3))  :: X,Y,Z,f
  integer :: dm(1:ex(1),1:ex(2),1:ex(3))
  integer, intent(in)                   :: Symmetry
  integer                      :: i, imin, imax
  integer                      :: j, jmin, jmax
  integer                      :: k, kmin, kmax
  real*8                       :: dX, dY, dZ
  real*8, dimension(4)         :: xa, ya
  real*8                       :: ftemp
  real*8                       :: r
  real*8                       :: err
  integer                      :: im, l
  integer :: ksym,jsym,isym,diax,diay,diaz
  real*8, parameter    :: ZERO  = 0.D0, ONE  = 1.D0, TWO = 2.D0, FOUR = 4.D0
  real*8  :: sf(1:3),rr
  integer, parameter :: inEZ=0, onbd=1, outEZ=2
!
! Input translation
!
  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)
  imax = ubound(f,1)
  jmax = ubound(f,2)
  kmax = ubound(f,3)
  dX = X(imin + 1,1,1) - X(imin,1,1)
  dY = Y(1,jmin + 1,1) - Y(1,jmin,1)
  dZ = Z(1,1,kmin + 1) - Z(1,1,kmin)
  if(X(1,1,1).lt.0.D0) then
     f(1,:,:) = f(2,:,:)
!     imin = imin + 1
  end if
  if(Z(1,1,1).lt.0.D0) then
     f(:,:,1) = f(:,:,2)
!     kmin = kmin + 1
  end if
  j = jmin + 1
  do k = kmin, kmax
    do i = imin, imax
!!$       if(dm(i,j,k)==outEZ) then
      r = sqrt(X(i,1,1)*X(i,1,1)+dY*dY)
      !interpolate f to the necessary x value
      if (i-imin < 1) then
        im = imin
      elseif (imax-i < 3) then
        im = imax - 3
      else
        im = i-1
      endif
!      if (i-imin < 2) then 
!        im = imin
!      elseif (imax-i < 3) then
!        im = imax - 3
!      else
!        im = i-1
!      endif
      if (dm(im,j,k)==inEZ) im = im+1
!      if (i==imin) then
!        ! Apply axisymmetry across the inner boundary
!        xa(1) = X(im,1,1)-dX
!        ya(1) = f(im,j,k)
!        do l=2,4
!          xa(l) = X(im+l-2,1,1)
!          ya(l) = f(im+l-2,j,k)
!        end do
!      else
        do l=1,4
          xa(l) = X(im+l-1,1,1)
          ya(l) = f(im+l-1,j,k)
        end do
!      end if 
      call polint(xa,ya,4,r,ftemp,err)
      !rotate in positive direction and assign f
      f(i,j+1,k) = ftemp
      !rotate in negative direction and assign f 
      f(i,j-1,k) =  ftemp
!!$      end if
    enddo
  enddo
  !now apply symmetry BCs:
!  if(X(1,1,1).lt.0.D0) then
!     f(1,:,:) = f(2,:,:)
!  end if
end subroutine ezaxibc_scalar
!---------------------------------------------------------------
! A subroutine to correct the vector f^a
! on the y/= 0 planes in the axisymmetric case
!---------------------------------------------------------------
subroutine ezaxibc_vector(ex,X,Y,Z,fx,fy,fz,dm,Symmetry)
  implicit none
  integer, dimension(3)                 :: ex
  real*8, dimension(ex(1),ex(2),ex(3))  :: X,Y,Z,fx, fy, fz
  integer :: dm(1:ex(1),1:ex(2),1:ex(3))
  integer, intent(in)                   :: Symmetry
  integer                      :: i, imin, imax
  integer                      :: j, jmin, jmax
  integer                      :: k, kmin, kmax
  real*8                       :: dX, dY, dZ
  real*8, dimension(4)         :: xa, ya
  real*8                       :: fxtemp, fytemp, fztemp
  real*8                       :: r
  real*8                       :: err
  integer                      :: im, l
  real*8                       :: cosine, sine
  integer :: ksym,jsym,isym,diax,diay,diaz
  real*8, parameter    :: ZERO  = 0.D0, ONE  = 1.D0, TWO = 2.D0, FOUR = 4.D0
  real*8  :: sf(1:3),rr
  integer, parameter :: inEZ=0, onbd=1, outEZ=2
!
! Input translation
!
  imin = lbound(fx,1)
  jmin = lbound(fx,2)
  kmin = lbound(fx,3)
  imax = ubound(fx,1)
  jmax = ubound(fx,2)
  kmax = ubound(fx,3)
  dX = X(imin + 1,1,1) - X(imin,1,1)
  dY = Y(1,jmin + 1,1) - Y(1,jmin,1)
  dZ = Z(1,1,kmin + 1) - Z(1,1,kmin)
  if(X(1,1,1).lt.0.D0) then
     fx(1,:,:) = -fx(2,:,:)
     fy(1,:,:) = -fy(2,:,:)
     fz(1,:,:) = fz(2,:,:)
!     imin = imin + 1
  end if
  if(Z(1,1,1).lt.0.D0) then
     fx(:,:,1) = fx(:,:,2)
     fy(:,:,1) = fy(:,:,2)
     fz(:,:,1) = -fz(:,:,2)
!     kmin = kmin + 1
  end if
  !set j = jmin+1
  j = jmin + 1
  do k = kmin, kmax
    do i = imin, imax
!!$       if(dm(i,j,k)==outEZ) then
      r = sqrt(X(i,1,1)*X(i,1,1)+dY*dY)
      !interpolate fx, fy, and fz to the necessary x value
      if (i-imin < 1) then
        im = imin
      elseif (imax-i < 3) then
        im = imax - 3
      else
        im = i-1
      endif
      if (dm(im,j,k)==inEZ) im = im+1
!      if (i==imin) then
!        ! Apply axisymmetry across the inner boundary
!        xa(1) = X(im,1,1)-dX
!        ya(1) = -fx(im,j,k)
!        do l=2,4
!          xa(l) = X(im+l-2,1,1)
!          ya(l) = fx(im+l-2,j,k)
!        end do
!      else
        do l=1,4
          xa(l) = X(im+l-1,1,1)
          ya(l) = fx(im+l-1,j,k)
        end do
!      end if 
      call polint(xa,ya,4,r,fxtemp,err)
!      if (i==imin) then
!        xa(1) = X(im,1,1)-dX
!        ya(1) = -fy(im,j,k)
!        do l=2,4
!          xa(l) = X(im+l-2,1,1)
!          ya(l) = fy(im+l-2,j,k)
!        end do
!      else
        do l=1,4
          ya(l) = fy(im+l-1,j,k)
        end do
!      end if
      call polint(xa,ya,4,r,fytemp,err)
!     if (i==imin) then
!        ya(1) = fz(im,j,k)
!        do l=2,4
!          ya(l) = fz(im+l-2,j,k)
!        end do
!      else
        do l=1,4
          ya(l) = fz(im+l-1,j,k)
        end do
!      end if
      call polint(xa,ya,4,r,fztemp,err)
      cosine = X(i,1,1)/r
      sine = dY/r
      !rotate in positive direction and assign f^a
      fx(i,j+1,k) = cosine*fxtemp - sine*fytemp
      fy(i,j+1,k) = sine*fxtemp + cosine*fytemp
      fz(i,j+1,k) = fztemp
      !rotate in negative direction and assign f^a 
      fx(i,j-1,k) =  cosine*fxtemp + sine*fytemp
      fy(i,j-1,k) = -sine*fxtemp + cosine*fytemp
      fz(i,j-1,k) = fztemp
!!$      end if
    enddo
  enddo
end subroutine ezaxibc_vector
!---------------------------------------------------------------
! A subroutine to correct the symmetric tensor fij
! on the y/= 0 planes in the axisymmetric case
!---------------------------------------------------------------
subroutine ezaxibc_tensor(ex,x,y,z,fxx,fxy,fxz,fyy,fyz,fzz,dm,Symmetry)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: x,y,z,fxx,fxy,fxz,fyy,fyz,fzz
  integer :: dm(1:ex(1),1:ex(2),1:ex(3))
  integer, intent(in)                      :: Symmetry
  real*8,dimension(1,3)                   :: pointcoords
!
  integer                    :: i, j, k, l
  integer                    :: imin, jmin, kmin, imax, jmax, kmax, jm, im, jj
  real*8                     :: fxx_int, fxy_int, fxz_int, fyy_int, fyz_int, fzz_int
  real*8                     :: pm, df, theta, cs, sn, dx, dy, dz, r, two
  real*8, dimension(4)             :: xa, fxxa, fxya, fxza, fyya, fyza, fzza
  real*8                     :: rxx,ryx,rxy,ryy
  real*8                     :: sxx,syx,sxy,syy
  integer :: ksym,jsym,isym,diax,diay,diaz
  real*8  :: sf(1:3),rr
  integer, parameter :: inEZ=0, onbd=1, outEZ=2
  two = 2.d0
!
  imin = lbound(fxx,1)
  jmin = lbound(fxx,2)
  kmin = lbound(fxx,3)
  imax = ubound(fxx,1)
  jmax = ubound(fxx,2)
  kmax = ubound(fxx,3)
  dX = X(imin + 1,1,1) - X(imin,1,1)
  dY = Y(1,jmin + 1,1) - Y(1,jmin,1)
  dZ = Z(1,1,kmin + 1) - Z(1,1,kmin)
  if(X(1,1,1).lt.0.D0) then
     fxx(1,2,:) = fxx(2,2,:)
     fxy(1,2,:) = fxy(2,2,:)
     fxz(1,2,:) = -fxz(2,2,:)
     fyy(1,2,:) = fyy(2,2,:)
     fyz(1,2,:) = -fyz(2,2,:)
     fzz(1,2,:) = fzz(2,2,:)
!     imin = imin + 1
  end if
  if(Z(1,1,1).lt.0.D0) then
     fxx(:,2,1) = fxx(:,2,2)
     fxy(:,2,1) = fxy(:,2,2)
     fxz(:,2,1) = -fxz(:,2,2)
     fyy(:,2,1) = fyy(:,2,2)
     fyz(:,2,1) = -fyz(:,2,2)
     fzz(:,2,1) = fzz(:,2,2)
!     kmin = kmin + 1
  end if
  !set j = jmin+1
  j = jmin + 1
  do k = kmin, kmax
     do i = imin, imax
!!$       if(dm(i,j,k)==outEZ) then
        r = sqrt(X(i,1,1)*X(i,1,1)+dY*dY)
        !interpolate fx, fy, and fz to the necessary x value
        if (i-imin < 1) then
           im = imin
        elseif (imax-i < 3) then
           im = imax - 3
        else
           im = i-1
        endif
      if (dm(im,j,k)==inEZ) im = im+1
      if(im+4-1 .gt. imax) then
         !Following function uses return_funcvals_at_points() inside gw_extraction thorn to get interpolation stencil at points
         !  This is a global grid operation, but shouldnt happen very often.
         if(im+4-1-1 .gt. imax) then
            write(*,*) "interpolation stencil WAYYYY out of bounds!",im+4-1,j,k,imax
            stop
         end if
         pointcoords(1,1) = X(im+4-1-1,1,1)+dX
         pointcoords(1,2) = Y(1,j,1)
         pointcoords(1,3) = Z(1,1,k)
         do l=1,3
            xa(l) = X(im+l-1,1,1)
            fxxa(l) = fxx(im+l-1,j,k)
            fxya(l) = fxy(im+l-1,j,k)
            fxza(l) = fxz(im+l-1,j,k)
            fyya(l) = fyy(im+l-1,j,k)
            fyza(l) = fyz(im+l-1,j,k)
            fzza(l) = fzz(im+l-1,j,k)
         end do
         xa(4) = pointcoords(1,1)
         !TODO FIXME: add cctkGH,cctk_gsh,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal
!         call return_funcvals_at_points(cctkGH,1,ext,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoor
!         call return_funcvals_at_points(cctkGH,1,ext,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoor
!         call return_funcvals_at_points(cctkGH,1,ext,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoor
!         call return_funcvals_at_points(cctkGH,1,ext,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoor
!         call return_funcvals_at_points(cctkGH,1,ext,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoor
!         call return_funcvals_at_points(cctkGH,1,ext,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoor
         !write(*,*) Applied kludge for interpolation stencil, since it was out of bounds!,im,j,k,imax
         write(*,*) "TODO FIXME: add cctkGH,cctk_gsh,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal"
         write(*,*) "TODO.... Need to apply kludge for interpolation stencil, since it was out of bounds!",im,j,k,imax
         stop
      else
         do l=1,4
            xa(l) = X(im+l-1,1,1)
            fxxa(l) = fxx(im+l-1,j,k)
            fxya(l) = fxy(im+l-1,j,k)
            fxza(l) = fxz(im+l-1,j,k)
            fyya(l) = fyy(im+l-1,j,k)
            fyza(l) = fyz(im+l-1,j,k)
            fzza(l) = fzz(im+l-1,j,k)
         end do
      end if
        call polint(xa,fxxa,4,r,fxx_int,df)
        call polint(xa,fxya,4,r,fxy_int,df)
        call polint(xa,fxza,4,r,fxz_int,df)
        call polint(xa,fyya,4,r,fyy_int,df)
        call polint(xa,fyza,4,r,fyz_int,df)
        call polint(xa,fzza,4,r,fzz_int,df)
        rxx = X(i,1,1)/r
        ryx = dY/r
        rxy = -ryx
        ryy = rxx
        sxx = ryy
        syx = -ryx
        sxy = -rxy
        syy = rxx
        fxx(i,j+1,k) = fxx_int*sxx*sxx + 2*fxy_int*sxx*syx + fyy_int*syx*syx
        fyy(i,j+1,k) = fxx_int*sxy*sxy + 2*fxy_int*sxy*syy + fyy_int*syy*syy
        fxy(i,j+1,k) = fxx_int*sxx*sxy + fxy_int*(sxy*syx+sxx*syy) + fyy_int*syx*syy
        fxz(i,j+1,k) = fxz_int*sxx + fyz_int*syx
        fyz(i,j+1,k) = fxz_int*sxy + fyz_int*syy
        fzz(i,j+1,k) = fzz_int
        fxx(i,j-1,k) = fxx_int*rxx*rxx + 2*fxy_int*rxx*ryx + fyy_int*ryx*ryx
        fyy(i,j-1,k) = fxx_int*rxy*rxy + 2*fxy_int*rxy*ryy + fyy_int*ryy*ryy
        fxy(i,j-1,k) = fxx_int*rxx*rxy + fxy_int*(rxy*ryx+rxx*ryy) + fyy_int*ryx*ryy
        fxz(i,j-1,k) = fxz_int*rxx + fyz_int*ryx
        fyz(i,j-1,k) = fxz_int*rxy + fyz_int*ryy
        fzz(i,j-1,k) = fzz_int
        if(i==2 .and. j==2 .and. k==15) then
!           write(*,*) blah,dm(im,j,k),inEZ
           write(*,*) "hello axisym-excision.F90",xa,fxxa
           write(*,*) "hello axisym-excision.F90",fxx_int,fxy_int,fxz_int,fyy_int,fyz_int,fzz_int
           write(*,*) "hello axisym-excision.F90",fxx(i,j-1,k),fxy(i,j-1,k),fxz(i,j-1,k),fyy(i,j-1,k),fyz(i,j-1,k),fzz(i,j-1,k)
           write(*,*) "hello axisym-excision.F90",fxx(i,j+1,k),fxy(i,j+1,k),fxz(i,j+1,k),fyy(i,j+1,k),fyz(i,j+1,k),fzz(i,j+1,k)
        end if
!!$     end if
     end do
  end do
  if(1==0 .and. X(1,1,1).lt.0.D0) then
     fxx(1,3,:) = fxx(2,3,:)
     fxy(1,3,:) = fxy(2,3,:)
     fxz(1,3,:) = -fxz(2,3,:)
     fyy(1,3,:) = fyy(2,3,:)
     fyz(1,3,:) = -fyz(2,3,:)
     fzz(1,3,:) = fzz(2,3,:)
     fxx(1,1,:) = fxx(2,1,:)
     fxy(1,1,:) = fxy(2,1,:)
     fxz(1,1,:) = -fxz(2,1,:)
     fyy(1,1,:) = fyy(2,1,:)
     fyz(1,1,:) = -fyz(2,1,:)
     fzz(1,1,:) = fzz(2,1,:)
  end if
  if(1==0 .and. Z(1,1,1).lt.0.D0) then
     fxx(:,3,1) = fxx(:,3,2)
     fxy(:,3,1) = fxy(:,3,2)
     fxz(:,3,1) = -fxz(:,3,2)
     fyy(:,3,1) = fyy(:,3,2)
     fyz(:,3,1) = -fyz(:,3,2)
     fzz(:,3,1) = fzz(:,3,2)
     fxx(:,1,1) = fxx(:,1,2)
     fxy(:,1,1) = fxy(:,1,2)
     fxz(:,1,1) = -fxz(:,1,2)
     fyy(:,1,1) = fyy(:,1,2)
     fyz(:,1,1) = -fyz(:,1,2)
     fzz(:,1,1) = fzz(:,1,2)
  end if
end subroutine ezaxibc_tensor
