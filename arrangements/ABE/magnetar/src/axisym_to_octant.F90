!---------------------------------------------------------------
! A subroutine to correct the scalar f
! on the y/= 0 planes in the axisymmetric case
!---------------------------------------------------------------

subroutine axi_to_oct_scalar(ex,X,Y,Z,fin,f)
  implicit none

  integer, dimension(3)                 :: ex
  real*8, dimension(ex(1),ex(2),ex(3))  :: X,Y,Z,f
  real*8, dimension(ex(1)*ex(3))        :: fin

  integer                      :: i, imin, imax
  integer                      :: j, jmin, jmax
  integer                      :: k, kmin, kmax
  real*8                       :: dX, dY, dZ
  real*8, dimension(4)         :: xa, ya
  real*8                       :: ftemp
  real*8                       :: r    !the radius
  real*8                       :: err  !error from the interpolation
  integer                      :: im, l
  integer :: ksym,jsym,isym,diax,diay,diaz
  real*8, parameter    :: ZERO  = 0.D0, ONE  = 1.D0, TWO = 2.D0, FOUR = 4.D0

  !
  ! Input translation
  !
  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)
  imax = ubound(f,1)
  jmax = ubound(f,2)
  kmax = ubound(f,3)

!  imin = 3
!  jmin = 3
!  kmin = 3

  dX = X(imin + 1,1,1) - X(imin,1,1)
  dY = Y(1,jmin + 1,1) - Y(1,jmin,1)
  dZ = Z(1,1,kmin + 1) - Z(1,1,kmin)

  if(X(1,1,1).lt.0.D0) then
!     f(1,:,:) = fin(4,:,:)
!     f(2,:,:) = fin(3,:,:)
     !     imin = imin + 1
  end if
  if(Z(1,1,1).lt.0.D0) then
!     f(:,:,1) = fin(:,:,4)
!     f(:,:,2) = fin(:,:,3)
     !     kmin = kmin + 1
  end if

  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax

           r = sqrt(X(i,j,k)*X(i,j,k)+Y(i,j,k)*Y(i,j,k))

           im = r/dX
           !interpolate f to the necessary x value
           if (i-imin < 1) then 
              im = imin
           elseif (imax-i < 3) then
              im = imax - 3
           else
              im = i-1
           endif

           do l=1,4
              xa(l) = X(im+l,1,1)  !X(im+l-1,1,1)
              if(im.lt.1) then 
!                 write(*,*) "hello!",im
              end if
              ya(l) = fin((im+l)+ex(1)*k) !fin(im+l-1,2,k)
           end do

           call polint(xa,ya,4,r,ftemp,err)

           if(j==3 .and. ftemp==0.D0) then 
              !write(*,*) "hi",i,j,xa,ya,r,ftemp
           end if

           !rotate in positive direction and assign f
           f(i,j,k) = ftemp

           !rotate in negative direction and assign f 
!           f(i,j-1,k) =  ftemp

        enddo
     enddo
  enddo

  !now apply symmetry BC's:
  !  if(X(1,1,1).lt.0.D0) then
  !     f(1,:,:) = f(2,:,:)
  !  end if
end subroutine axi_to_oct_scalar


!---------------------------------------------------------------
! A subroutine to correct the vector f^a
! on the y/= 0 planes in the axisymmetric case
!---------------------------------------------------------------

subroutine axi_to_oct_vector(ex,X,Y,Z,fxin,fyin,fzin,fx,fy,fz)
  implicit none

  integer, dimension(3)                 :: ex
  real*8, dimension(ex(1),ex(2),ex(3))  :: X,Y,Z,fx, fy, fz
  real*8, dimension(ex(1)*ex(3))        :: fxin, fyin, fzin

  integer                      :: i, imin, imax
  integer                      :: j, jmin, jmax
  integer                      :: k, kmin, kmax
  real*8                       :: dX, dY, dZ
  real*8, dimension(4)         :: xa, ya
  real*8                       :: fxtemp, fytemp, fztemp
  real*8                       :: r    !the radius
  real*8                       :: err  !error from the interpolation
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

  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax

           r = sqrt(X(i,j,k)*X(i,j,k)+Y(i,j,k)*Y(i,j,k))

           im = r/dX

           !interpolate fx, fy, and fz to the necessary x value
           if (i-imin < 1) then 
              im = imin
           elseif (imax-i < 3) then
              im = imax - 3
           else
              im = i-1
           endif

           do l=1,4
              xa(l) = X((im+l),1,1)
              ya(l) = fxin((im+l)+ex(1)*k)
           end do

           call polint(xa,ya,4,r,fxtemp,err)

           do l=1,4
              ya(l) = fyin((im+l)+ex(1)*k)
           end do

           call polint(xa,ya,4,r,fytemp,err)

           do l=1,4
              ya(l) = fzin((im+l)+ex(1)*k)
           end do

           call polint(xa,ya,4,r,fztemp,err)

           cosine = X(i,1,1)/r
           sine = Y(i,j,k)/r

           !rotate in positive direction and assign f^a
           fx(i,j+1,k) = cosine*fxtemp - sine*fytemp
           fy(i,j+1,k) = sine*fxtemp + cosine*fytemp
           fz(i,j+1,k) = fztemp

           !rotate in negative direction and assign f^a 
           fx(i,j-1,k) =  cosine*fxtemp + sine*fytemp
           fy(i,j-1,k) = -sine*fxtemp + cosine*fytemp
           fz(i,j-1,k) = fztemp

        enddo
     enddo
  enddo
end subroutine axi_to_oct_vector
!---------------------------------------------------------------
! A subroutine to correct the symmetric tensor fij
! on the y/= 0 planes in the axisymmetric case
!---------------------------------------------------------------
subroutine axi_to_oct_tensor(ex,x,y,z,fxxin,fxyin,fxzin,fyyin,fyzin,fzzin,fxx,fxy,fxz,fyy,fyz,fzz)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: x,y,z,fxx,fxy,fxz,fyy,fyz,fzz
  real*8, dimension(ex(1)*ex(3))           :: fxxin,fxyin,fxzin,fyyin,fyzin,fzzin
  !
  integer                    :: i, j, k, l
  integer                    :: imin, jmin, kmin, imax, jmax, kmax, jm, im, jj
  real*8                     :: fxx_int, fxy_int, fxz_int, fyy_int, fyz_int, fzz_int
  real*8                     :: pm, df, theta, cs, sn, dX, dY, dZ, r, two
  real*8, dimension(4)	     :: xa, fxxa, fxya, fxza, fyya, fyza, fzza

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

!  if(X(1,1,1).lt.0.D0) then
!     fxx(1,2,:) = fxxin(2,2,:)
!     fxy(1,2,:) = fxyin(2,2,:)
!     fxz(1,2,:) = -fxzin(2,2,:)
!     fyy(1,2,:) = fyyin(2,2,:)
!     fyz(1,2,:) = -fyzin(2,2,:)
!     fzz(1,2,:) = fzzin(2,2,:)
!  end if
!  if(Z(1,1,1).lt.0.D0) then
!     fxx(:,2,1) = fxxin(:,2,2)
!     fxy(:,2,1) = fxyin(:,2,2)
!     fxz(:,2,1) = -fxzin(:,2,2)
!     fyy(:,2,1) = fyyin(:,2,2)
!     fyz(:,2,1) = -fyzin(:,2,2)
!     fzz(:,2,1) = fzzin(:,2,2)
!  end if

  imin = 3
  jmin = 3
  kmin = 3
  
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax

           r = sqrt(X(i,j,k)*X(i,j,k)+Y(i,j,k)*Y(i,j,k))

           im = r/dX
           !interpolate fx, fy, and fz to the necessary x value
           if (i-imin < 1) then 
              im = imin
           elseif (imax-i < 3) then
              im = imax - 3
           else
              im = i-1
           endif

           do l=1,4
              xa(l) = X((im+l),1,1) !im+l-1,1,1)
              fxxa(l) = fxxin((im+l)+ex(1)*k)
              fxya(l) = fxyin((im+l)+ex(1)*k)
              fxza(l) = fxzin((im+l)+ex(1)*k)
              fyya(l) = fyyin((im+l)+ex(1)*k)
              fyza(l) = fyzin((im+l)+ex(1)*k)
              fzza(l) = fzzin((im+l)+ex(1)*k)
           end do

           call polint(xa,fxxa,4,r,fxx_int,df)
           call polint(xa,fxya,4,r,fxy_int,df)
           call polint(xa,fxza,4,r,fxz_int,df)
           call polint(xa,fyya,4,r,fyy_int,df)
           call polint(xa,fyza,4,r,fyz_int,df)
           call polint(xa,fzza,4,r,fzz_int,df)

           rxx = X(i,1,1)/r
           ryx = Y(i,j,k)/r
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
              write(*,*) "hello axisym-excision.F90",xa,fxxa
              write(*,*) "hello axisym-excision.F90",fxx_int,fxy_int,fxz_int,fyy_int,fyz_int,fzz_int
              write(*,*) "hello axisym-excision.F90",fxx(i,j-1,k),fxy(i,j-1,k),fxz(i,j-1,k),fyy(i,j-1,k),fyz(i,j-1,k),fzz(i,j-1,k)
              write(*,*) "hello axisym-excision.F90",fxx(i,j+1,k),fxy(i,j+1,k),fxz(i,j+1,k),fyy(i,j+1,k),fyz(i,j+1,k),fzz(i,j+1,k)
           end if


        end do
     end do
  end do
end subroutine axi_to_oct_tensor
