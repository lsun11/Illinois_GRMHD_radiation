!---------------------------------------------------------------
! A subroutine to correct the scalar f
! on the y/= 0 planes in the axisymmetric case
!---------------------------------------------------------------

subroutine axi_to_eq_scalar(ex,X,Y,Z,fin,f)
  implicit none

  integer, dimension(3)                          :: ex
  real*8, dimension(ex(1),ex(2),ex(3))           :: X,Y,Z,f
  real*8, dimension(ex(1),ex(2),ex(3))           :: fin

  integer                                        :: i, imin, imax
  integer                                        :: j, jmin, jmax
  integer                                        :: k, kmin, kmax
  real*8                                         :: dX, dY, dZ
  real*8, dimension(4)                           :: xa, ya
  real*8                                         :: ftemp
  real*8                                         :: r    !the radius
  real*8                                         :: err  !error from the interpolation
  integer                                        :: im, l
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

  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax

           r = sqrt(X(i,1,1)*X(i,1,1)+Y(1,j,1)*Y(1,j,1))      

           !interpolate f to the necessary x value
           if (i-imin < 1) then 
              im = imin
           elseif (imax-i < 3) then
              im = imax - 3
           else
              im = i-1
           endif

           do l=1,4
              xa(l) = X(im+l-1,2,1)
              ya(l) = fin(im+l-1,2,k)
           end do

           call polint(xa,ya,4,r,ftemp,err)

           !rotate and assign f
           f(i,j,k) = ftemp

        enddo
     enddo
  enddo
end subroutine axi_to_eq_scalar


!---------------------------------------------------------------
! A subroutine to correct the vector f^a
! on the y/= 0 planes in the axisymmetric case
!---------------------------------------------------------------

subroutine axi_to_eq_vector(ex,X,Y,Z,fxin,fyin,fzin,fx,fy,fz)
  implicit none

  integer, dimension(3)                 :: ex
  real*8, dimension(ex(1),ex(2),ex(3))           :: X,Y,Z,fx, fy, fz
  real*8, dimension(ex(1),ex(2),ex(3))           :: fxin,fyin,fzin

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

           r = sqrt(X(i,1,1)*X(i,1,1)+Y(1,j,1)*Y(1,j,1))      

           !interpolate fx, fy, and fz to the necessary x value
           if (i-imin < 1) then 
              im = imin
           elseif (imax-i < 3) then
              im = imax - 3
           else
              im = i-1
           endif

           do l=1,4
              xa(l) = X(im+l-1,1,1)
              ya(l) = fxin(im+l-1,2,k)
           end do

           call polint(xa,ya,4,r,fxtemp,err)

           do l=1,4
              ya(l) = fyin(im+l-1,2,k)
           end do

           call polint(xa,ya,4,r,fytemp,err)
           do l=1,4
              ya(l) = fzin(im+l-1,2,k)
           end do

           call polint(xa,ya,4,r,fztemp,err)

           cosine = X(i,1,1)/r
           sine = Y(1,j,1)/r

           !rotate and assign f^a
           fx(i,j,k) = cosine*fxtemp - sine*fytemp
           fy(i,j,k) = sine*fxtemp + cosine*fytemp
           fz(i,j,k) = fztemp
        enddo
     enddo
  enddo
end subroutine axi_to_eq_vector

!---------------------------------------------------------------
! A subroutine to correct the symmetric tensor fij
! on the y/= 0 planes in the axisymmetric case
!---------------------------------------------------------------
subroutine axi_to_eq_tensor(ex,X,Y,Z,fxxin,fxyin,fxzin,fyyin,fyzin,fzzin,fxx,fxy,fxz,fyy,fyz,fzz)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: x,y,z,fxx,fxy,fxz,fyy,fyz,fzz
  real*8, dimension(100,3,100)             :: Xin,Yin,Zin,fxxin,fxyin,fxzin,fyyin,fyzin,fzzin
  real*8,dimension(1,3)                   :: pointcoords
  !
  integer                    :: i, j, k, l
  integer                    :: imin, jmin, kmin, imax, jmax, kmax, jm, im, jj
  real*8                     :: fxx_int, fxy_int, fxz_int, fyy_int, fyz_int, fzz_int
  real*8                     :: pm, df, theta, cs, sn, dx, dy, dz, r, two
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

  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax

           r = sqrt(X(i,1,1)*X(i,1,1)+Y(1,j,1)*Y(1,j,1))      

           !interpolate fx, fy, and fz to the necessary x value
           if (i-imin < 1) then 
              im = imin
           elseif (imax-i < 3) then
              im = imax - 3
           else
              im = i-1
           endif

           do l=1,4
              xa(l) = Xin(im+l-1,1,1)
              fxxa(l) = fxxin(im+l-1,2,k)
              fxya(l) = fxyin(im+l-1,2,k)
              fxza(l) = fxzin(im+l-1,2,k)
              fyya(l) = fyyin(im+l-1,2,k)
              fyza(l) = fyzin(im+l-1,2,k)
              fzza(l) = fzzin(im+l-1,2,k)
           end do

           call polint(xa,fxxa,4,r,fxx_int,df)
           call polint(xa,fxya,4,r,fxy_int,df)
           call polint(xa,fxza,4,r,fxz_int,df)
           call polint(xa,fyya,4,r,fyy_int,df)
           call polint(xa,fyza,4,r,fyz_int,df)
           call polint(xa,fzza,4,r,fzz_int,df)

           rxx = X(i,1,1)/r
           ryx = Y(1,j,1)/r
           rxy = -ryx 
           ryy = rxx

           sxx = ryy
           syx = -ryx 
           sxy = -rxy 
           syy = rxx

           fxx(i,j,k) = fxx_int*sxx*sxx + 2*fxy_int*sxx*syx + fyy_int*syx*syx
           fyy(i,j,k) = fxx_int*sxy*sxy + 2*fxy_int*sxy*syy + fyy_int*syy*syy
           fxy(i,j,k) = fxx_int*sxx*sxy + fxy_int*(sxy*syx+sxx*syy) + fyy_int*syx*syy
           fxz(i,j,k) = fxz_int*sxx + fyz_int*syx
           fyz(i,j,k) = fxz_int*sxy + fyz_int*syy
           fzz(i,j,k) = fzz_int
        end do
     end do
  end do

end subroutine axi_to_eq_tensor
