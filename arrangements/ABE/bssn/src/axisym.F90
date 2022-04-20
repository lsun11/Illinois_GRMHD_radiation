!---------------------------------------------------------------
! A subroutine to correct the scalar f
! on the y/= 0 planes in the axisymmetric case
!---------------------------------------------------------------

subroutine axibc_scalar(ex,X,Y,Z,f)
  implicit none

  integer, dimension(3)                 :: ex
  real*8, dimension(ex(1),ex(2),ex(3))  :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))  :: f

  integer                      :: i, imin, imax
  integer                      :: j, jmin, jmax
  integer                      :: k, kmin, kmax
  real*8                       :: dX, dY, dZ
  real*8, dimension(4)         :: xa, ya
  real*8                       :: ftemp
  real*8                       :: r    !the radius
  real*8                       :: err  !error from the interpolation
  integer                      :: im, l

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

  !added below to fix for ghost zones across axis...  we compute zones across axis at the end of this routine.
  do i = 0, imin
     if(X(imin,1,1) .lt. 0.D0) then
        imin = imin + 1     
     end if
  end do

!  write(*,*) "axibc: IMIN = ",imin

  !set j = jmin+1
  j = jmin + 1

  do k = kmin, kmax
    do i = imin, imax

      r = sqrt(X(i,1,1)*X(i,1,1)+dY*dY)      

      !interpolate f to the necessary x value
      if (i-imin < 2) then 
        im = imin
      elseif (imax-i < 3) then
        im = imax - 3
      else
        im = i-1
      endif

      if (i==imin) then
       ! Apply axisymmetry across the inner boundary
        xa(1) = X(im,1,1)-dX
        ya(1) = f(im,j,k)
        do l=2,4
          xa(l) = X(im+l-2,1,1)
          ya(l) = f(im+l-2,j,k)
        end do
      else
        do l=1,4
          xa(l) = X(im+l-1,1,1)
          ya(l) = f(im+l-1,j,k)
        end do
      end if 
         
      call polint_nr_axi(xa,ya,4,r,ftemp,err)

      !rotate in positive direction and assign f
      f(i,j+1,k) = ftemp

      !rotate in negative direction and assign f 
      f(i,j-1,k) =  ftemp

    enddo
  enddo

  if(imin == 2) then
     i = imin - 1
     f(i,:,:) = f(i+1,:,:)
  end if
  if(imin == 3) then
     i = imin - 2
     f(i,:,:) = f(i+3,:,:)
     
     i = imin - 1 
     f(i,:,:) = f(i+1,:,:)
  end if

end subroutine axibc_scalar


!---------------------------------------------------------------
! A subroutine to correct the vector f^a
! on the y/= 0 planes in the axisymmetric case
!---------------------------------------------------------------

subroutine axibc_vector(ex,X,Y,Z,fx,fy,fz)
  implicit none

  integer, dimension(3)                 :: ex
  real*8, dimension(ex(1),ex(2),ex(3))  :: X,Y,Z,fx, fy, fz

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

  !added below to fix for ghost zones across axis...  we compute zones across axis at the end of this routine..
  do i = 0, imin
     if(X(imin,1,1) .lt. 0.D0) then
        imin = imin + 1     
     end if
  end do

  !set j = jmin+1
  j = jmin + 1

  do k = kmin, kmax
    do i = imin, imax

      r = sqrt(X(i,1,1)*X(i,1,1)+dY*dY)      

      !interpolate fx, fy, and fz to the necessary x value
      if (i-imin < 2) then 
        im = imin
      elseif (imax-i < 3) then
        im = imax - 3
      else
        im = i-1
      endif

      if (i==imin) then
        ! Apply axisymmetry across the inner boundary
        xa(1) = X(im,1,1)-dX
        ya(1) = -fx(im,j,k)
        do l=2,4
          xa(l) = X(im+l-2,1,1)
          ya(l) = fx(im+l-2,j,k)
        end do
      else
        do l=1,4
          xa(l) = X(im+l-1,1,1)
          ya(l) = fx(im+l-1,j,k)
        end do
      end if 
         
      call polint_nr_axi(xa,ya,4,r,fxtemp,err)

      if (i==imin) then
	xa(1) = X(im,1,1)-dX
        ya(1) = -fy(im,j,k)
        do l=2,4
	  xa(l) = X(im+l-2,1,1)
          ya(l) = fy(im+l-2,j,k)
        end do
      else
        do l=1,4
          ya(l) = fy(im+l-1,j,k)
        end do
      end if
    
      call polint_nr_axi(xa,ya,4,r,fytemp,err)
 
      if (i==imin) then
        ya(1) = fz(im,j,k)
        do l=2,4
          ya(l) = fz(im+l-2,j,k)
        end do
      else
        do l=1,4
          ya(l) = fz(im+l-1,j,k)
        end do
      end if
 
      call polint_nr_axi(xa,ya,4,r,fztemp,err)

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

    enddo
  enddo
  if(imin == 2) then
     
     !    fx(i,jmin+2,:) = -fy(i+1,jmin+2,:)
     !    fy(i,jmin+2,:) = fx(i+1,jmin+2,:)
     !    fz(i,jmin+2,:) = fz(i+1,jmin+2,:)
     !
     !    fx(i,jmin,:) = fy(i+1,jmin,:)
     !    fy(i,jmin,:) = -fx(i+1,jmin,:)
     !    fz(i,jmin,:) = fz(i+1,jmin,:)
     
     i = imin - 1
     j = jmin + 1
     do k=kmin,kmax
!        cosine = 0.D0
!        sine = 1.D0
        cosine = 0.6D0
        sine = 0.8D0
        
        fxtemp = fx(imin,j,k)
        fytemp = fy(imin,j,k)
        fztemp = fz(imin,j,k)
        
        !rotate in positive direction and assign f^a
        fx(i,j+1,k) = cosine*fxtemp - sine*fytemp
        fy(i,j+1,k) = sine*fxtemp + cosine*fytemp
        fz(i,j+1,k) = fztemp

        fx(i,j,k) = -fx(i+1,j,k)
        fy(i,j,k) = -fy(i+1,j,k)
        fz(i,j,k) = fz(i+1,j,k)
        
        !rotate in negative direction and assign f^a 
        fx(i,j-1,k) =  cosine*fxtemp + sine*fytemp
        fy(i,j-1,k) = -sine*fxtemp + cosine*fytemp
        fz(i,j-1,k) = fztemp
     end do


     !        fx(i,:,:) = -fy(i+1,:,:)
     !        fy(i,:,:) = fx(i+1,:,:)
     !        fz(i,:,:) = fz(i+1,:,:)
     !     write(*,*) "AXISYM: ",fx(1,1,3),fy(2,1,3)
  end if
  if(imin == 3) then
     write(*,*) "2 GHOST ZONES NOT SUPPORTED (YET)!"
     i = imin - 2
     fx(i,:,:) = -fx(i+3,:,:)
     fy(i,:,:) = -fy(i+3,:,:)
     fz(i,:,:) = fz(i+3,:,:)
     
     i = imin - 1 
     fx(i,:,:) = -fx(i+1,:,:)
     fy(i,:,:) = -fy(i+1,:,:)
     fz(i,:,:) = fz(i+1,:,:)
  end if

end subroutine axibc_vector

!---------------------------------------------------------------
! A subroutine to correct the symmetric tensor fij
! on the y/= 0 planes in the axisymmetric case
!---------------------------------------------------------------
subroutine axibc_tensor(ex,x,y,z,fxx,fxy,fxz,fyy,fyz,fzz)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: x,y,z,fxx,fxy,fxz,fyy,fyz,fzz
!
  integer                    :: i, j, k, l
  integer                    :: imin, jmin, kmin, imax, jmax, kmax, jm, im, jj
  real*8                     :: fxx_int, fxy_int, fxz_int, fyy_int, fyz_int, fzz_int
  real*8                     :: pm, df, theta, cs, sn, dx, dy, dz, r, two
  real*8, dimension(4)	     :: xa, fxxa, fxya, fxza, fyya, fyza, fzza
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

  !added below to fix for ghost zones across axis...  we compute zones across axis at the end of this routine..
  do i = 0, imin
     if(X(imin,1,1) .lt. 0.D0) then
        imin = imin + 1     
     end if
  end do
  !  write(*,*) "IMIN = ",imin

  !set j = jmin+1
  j = jmin + 1

  do k = kmin, kmax
     do i = imin, imax

        r = sqrt(X(i,1,1)*X(i,1,1)+dY*dY)      
        
        !interpolate fx, fy, and fz to the necessary x value
        if (i-imin < 2) then 
           im = imin
        elseif (imax-i < 3) then
           im = imax - 3
        else
           im = i-1
        endif
        
        if (i==imin) then
           !write(*,*) im,i
           ! Apply axisymmetry across the inner boundary
           xa(1) = X(im,1,1)-dX
           fxxa(1) = fxx(im,j,k)
           fxya(1) = fxy(im,j,k)
           fxza(1) = -fxz(im,j,k)
           fyya(1) = fyy(im,j,k)
           fyza(1) = -fyz(im,j,k)
           fzza(1) = fzz(im,j,k)
           do l=2,4
              xa(l) = X(im+l-2,1,1)
              fxxa(l) = fxx(im+l-2,j,k)
              fxya(l) = fxy(im+l-2,j,k)
              fxza(l) = fxz(im+l-2,j,k)
              fyya(l) = fyy(im+l-2,j,k)
              fyza(l) = fyz(im+l-2,j,k)
              fzza(l) = fzz(im+l-2,j,k)
           end do
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
         
        call polint_nr_axi(xa,fxxa,4,r,fxx_int,df)
        call polint_nr_axi(xa,fxya,4,r,fxy_int,df)
        call polint_nr_axi(xa,fxza,4,r,fxz_int,df)
        call polint_nr_axi(xa,fyya,4,r,fyy_int,df)
        call polint_nr_axi(xa,fyza,4,r,fyz_int,df)
        call polint_nr_axi(xa,fzza,4,r,fzz_int,df)
        
        cs = X(i,1,1)/r
        sn = dY/r
        
        fxx(i,j+1,k) = cs*cs*fxx_int - two*cs*sn*fxy_int &
             + sn*sn*fyy_int 
        fxy(i,j+1,k) = cs*sn*(fxx_int - fyy_int) &
             + (cs*cs - sn*sn)*fxy_int
        fxz(i,j+1,k) = cs*fxz_int - sn*fyz_int
        fyy(i,j+1,k) = sn*sn*fxx_int + two*cs*sn*fxy_int &
             + cs*cs*fyy_int
        fyz(i,j+1,k) = sn*fxz_int + cs*fyz_int
        fzz(i,j+1,k) = fzz_int
        
        sn = -dY/r
        
        fxx(i,j-1,k) = cs*cs*fxx_int - two*cs*sn*fxy_int &
             + sn*sn*fyy_int 
        fxy(i,j-1,k) = cs*sn*(fxx_int - fyy_int) &
             + (cs*cs - sn*sn)*fxy_int
        fxz(i,j-1,k) = cs*fxz_int - sn*fyz_int
        fyy(i,j-1,k) = sn*sn*fxx_int + two*cs*sn*fxy_int &
             + cs*cs*fyy_int
        fyz(i,j-1,k) = sn*fxz_int + cs*fyz_int
        fzz(i,j-1,k) = fzz_int
        
     end do
  end do

  if(imin == 2) then
     i = imin - 1
     j = jmin + 1

     do k=kmin,kmax        
        cs = 0.6D0
        sn = 0.8D0
!        cs = 0.0D0
!        sn = 1.0D0
        fxx_int = fxx(imin,j,k)
        fxy_int = fxy(imin,j,k)
        fxz_int = fxz(imin,j,k)
        fyy_int = fyy(imin,j,k)
        fyz_int = fyz(imin,j,k)
        fzz_int = fzz(imin,j,k)
        
        fxx(i,j+1,k) = cs*cs*fxx_int - two*cs*sn*fxy_int &
             + sn*sn*fyy_int 
        fxy(i,j+1,k) = cs*sn*(fxx_int - fyy_int) &
             + (cs*cs - sn*sn)*fxy_int
        fxz(i,j+1,k) = cs*fxz_int - sn*fyz_int
        fyy(i,j+1,k) = sn*sn*fxx_int + two*cs*sn*fxy_int &
             + cs*cs*fyy_int
        fyz(i,j+1,k) = sn*fxz_int + cs*fyz_int
        fzz(i,j+1,k) = fzz_int

        fxx(i,j,k) = fxx(i+1,j,k)
        fxy(i,j,k) = fxy(i+1,j,k)
        fxz(i,j,k) = -fxz(i+1,j,k)
        fyy(i,j,k) = fyy(i+1,j,k)
        fyz(i,j,k) = -fyz(i+1,j,k)
        fzz(i,j,k) = fzz(i+1,j,k)     
        

        sn = -0.8D0
!        sn = -1.0D0
        fxx(i,j-1,k) = cs*cs*fxx_int - two*cs*sn*fxy_int &
             + sn*sn*fyy_int 
        fxy(i,j-1,k) = cs*sn*(fxx_int - fyy_int) &
             + (cs*cs - sn*sn)*fxy_int
        fxz(i,j-1,k) = cs*fxz_int - sn*fyz_int
        fyy(i,j-1,k) = sn*sn*fxx_int + two*cs*sn*fxy_int &
             + cs*cs*fyy_int
        fyz(i,j-1,k) = sn*fxz_int + cs*fyz_int
        fzz(i,j-1,k) = fzz_int
     end do
!     OLD JUNK:
!     fxx(i,:,:) = fxx(i+1,:,:)
!     fxy(i,:,:) = -fxy(i+1,:,:)
!     fxz(i,:,:) = -fxz(i+1,:,:)
!     fyy(i,:,:) = fyy(i+1,:,:)
!     fyz(i,:,:) = -fyz(i+1,:,:)
!     fzz(i,:,:) = fzz(i+1,:,:)

  end if


  if(imin == 3) then
     do i=imin,imax
        write(*,*) "2 GHOST ZONES NOT SUPPORTED (YET)!"
     end do
     i = imin - 2
     fxx(i,:,:) = fxx(i+3,:,:)
     fxy(i,:,:) = fxy(i+3,:,:)
     fxz(i,:,:) = -fxz(i+3,:,:)
     fyy(i,:,:) = fyy(i+3,:,:)
     fyz(i,:,:) = -fyz(i+3,:,:)
     fzz(i,:,:) = fzz(i+3,:,:)     
     
     i = imin - 1 
     fxx(i,:,:) = fxx(i+1,:,:)
     fxy(i,:,:) = fxy(i+1,:,:)
     fxz(i,:,:) = -fxz(i+1,:,:)
     fyy(i,:,:) = fyy(i+1,:,:)
     fyz(i,:,:) = -fyz(i+1,:,:)
     fzz(i,:,:) = fzz(i+1,:,:)  
  end if

end subroutine axibc_tensor

!---------------------------------------------------------------
! A subroutine to correct the general tensor fij
! on the y/= 0 planes in the axisymmetric case
!---------------------------------------------------------------
subroutine axibc_gen_tensor(ex,X,Y,Z, &
     fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z,fxx,fxy,fxz
  real*8, dimension(ex(1),ex(2),ex(3))     :: fyx,fyy,fyz
  real*8, dimension(ex(1),ex(2),ex(3))     :: fzx,fzy,fzz
!
  integer                    :: i, j, k, n, im2
  integer                    :: imin, jmin, kmin, imax, jmax, kmax, jm, im, jj
  real*8                     :: fxx_int, fxy_int, fxz_int
  real*8                     :: fyx_int, fyy_int, fyz_int
  real*8                     :: fzx_int, fzy_int, fzz_int
  real*8                     :: pm, df, theta, dX, dY, dZ, cs, sn
  real*8, dimension(4)	     :: xa, fxxa, fxya, fxza
  real*8, dimension(4)	     :: fyxa, fyya, fyza
  real*8, dimension(4)	     :: fzxa, fzya, fzza
  real*8                     :: HALF, ONE, TWO, ZERO
  parameter(HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0)
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
  
  !set j = jmin+1
  j = jmin + 1

  do k = kmin, kmax
     do i = imin, imax
        
        pm=sqrt(X(i,1,1)*X(i,1,1)+dY*dY)
        if (i<imin+2) then 
           im=imin  
        else if (i>imax-3) then
           im=imax-3
        else
           im=i-1
        end if
        if (i==imin) then
           ! Apply axisymmetry across the inner boundary
           xa(1) = X(im,1,1)-dX
           fxxa(1) = fxx(im,j,k)
           fxya(1) = fxy(im,j,k)
           fxza(1) = -fxz(im,j,k)
           fyxa(1) = fyx(im,j,k)
           fyya(1) = fyy(im,j,k)
           fyza(1) = -fyz(im,j,k)
           fzxa(1) = -fzx(im,j,k)
           fzya(1) = -fzy(im,j,k)
           fzza(1) = fzz(im,j,k)
           do jj=2,4
              xa(jj) = X(im+jj-2,1,1)
              fxxa(jj)=fxx(im+jj-2,j,k)
              fxya(jj)=fxy(im+jj-2,j,k)
              fxza(jj)=fxz(im+jj-2,j,k)
              fyxa(jj)=fyx(im+jj-2,j,k)
              fyya(jj)=fyy(im+jj-2,j,k)
              fyza(jj)=fyz(im+jj-2,j,k)
              fzxa(jj)=fzx(im+jj-2,j,k)
              fzya(jj)=fzy(im+jj-2,j,k)
              fzza(jj)=fzz(im+jj-2,j,k)
           end do
        else
           do jj=0,3
              xa(jj+1)=X(im+jj,1,1)
              fxxa(jj+1)=fxx(im+jj,j,k)
              fxya(jj+1)=fxy(im+jj,j,k)
              fxza(jj+1)=fxz(im+jj,j,k)
              fyxa(jj+1)=fyx(im+jj,j,k)
              fyya(jj+1)=fyy(im+jj,j,k)
              fyza(jj+1)=fyz(im+jj,j,k)
              fzxa(jj+1)=fzx(im+jj,j,k)
              fzya(jj+1)=fzy(im+jj,j,k)
              fzza(jj+1)=fzz(im+jj,j,k)
           end do
        end if
        ! Interpolate, 1D 
        call polint_nr_axi(xa,fxxa,4,pm,fxx_int,df)
        call polint_nr_axi(xa,fxya,4,pm,fxy_int,df)
        call polint_nr_axi(xa,fxza,4,pm,fxz_int,df)
        call polint_nr_axi(xa,fyxa,4,pm,fyx_int,df)
        call polint_nr_axi(xa,fyya,4,pm,fyy_int,df)
        call polint_nr_axi(xa,fyza,4,pm,fyz_int,df)
        call polint_nr_axi(xa,fzxa,4,pm,fzx_int,df)
        call polint_nr_axi(xa,fzya,4,pm,fzy_int,df)
        call polint_nr_axi(xa,fzza,4,pm,fzz_int,df)

        cs = X(i,1,1)/pm
        sn = dY/pm

        fxx(i,j+1,k) = cs*cs*fxx_int - cs*sn*(fxy_int + fyx_int) &
             + sn*sn*fyy_int
        fxy(i,j+1,k) = cs*sn*(fxx_int - fyy_int) + cs*cs*fxy_int &
             - sn*sn*fyx_int
        fxz(i,j+1,k) = cs*fxz_int - sn*fyz_int
        fyx(i,j+1,k) = cs*sn*(fxx_int - fyy_int) + cs*cs*fyx_int &
             - sn*sn*fxy_int
        fyy(i,j+1,k) = sn*sn*fxx_int + cs*sn*(fxy_int + fyx_int) &
             + cs*cs*fyy_int
        fyz(i,j+1,k) = sn*fxz_int + cs*fyz_int
        fzx(i,j+1,k) = cs*fzx_int - sn*fzy_int
        fzy(i,j+1,k) = sn*fzx_int + cs*fzy_int
        fzz(i,j+1,k) = fzz_int
     
        sn = -sn
        
        fxx(i,j-1,k) = cs*cs*fxx_int - cs*sn*(fxy_int + fyx_int) &
             + sn*sn*fyy_int
        fxy(i,j-1,k) = cs*sn*(fxx_int - fyy_int) + cs*cs*fxy_int &
             - sn*sn*fyx_int
        fxz(i,j-1,k) = cs*fxz_int - sn*fyz_int
        fyx(i,j-1,k) = cs*sn*(fxx_int - fyy_int) + cs*cs*fyx_int &
             - sn*sn*fxy_int
        fyy(i,j-1,k) = sn*sn*fxx_int + cs*sn*(fxy_int + fyx_int) &
             + cs*cs*fyy_int
        fyz(i,j-1,k) = sn*fxz_int + cs*fyz_int
        fzx(i,j-1,k) = cs*fzx_int - sn*fzy_int
        fzy(i,j-1,k) = sn*fzx_int + cs*fzy_int
        fzz(i,j-1,k) = fzz_int
        
     enddo
  enddo
!
end subroutine axibc_gen_tensor
!
!-----------------------------------------------------------------------------
!
!  Numerical Recipe interpolation subroutine 
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE polint_nr_axi(xa,ya,n,x,y,dy)
      IMPLICIT NONE
      integer                    :: n,NMAX,i,m,ns
      real*8			 :: dy,x,y
      real*8, dimension(n)	 :: xa,ya
      PARAMETER (NMAX=10)
      real*8, dimension(NMAX)    :: c,d
      real*8                     :: den,dif,dift,ho,hp,w
!
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) then 
            write(*,*) 'failure in polint_nr_axi'
	    stop
	  end if
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END SUBROUTINE polint_nr_axi
! 
!-----------------------------------------------------------------------------
!
!  f -> X*f
!
!-----------------------------------------------------------------------------
!
  subroutine mult_x(ex,X,f)
  implicit none
  integer, dimension(3)			:: ex
  real*8, dimension(ex(1),ex(2),ex(3))  :: X,f
  integer                               :: i,j,k
  integer                               :: imin, jmin, kmin, imax, jmax, kmax
!
  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)
  imax = ubound(f,1)
  jmax = ubound(f,2)
  kmax = ubound(f,3)
!
     do k=kmin,kmax
        do j=jmin,jmax
           do i=imin,imax
              f(i,j,k)=X(i,1,1)*f(i,j,k)
           end do
        end do
     end do
end subroutine mult_x

! 
!-----------------------------------------------------------------------------
!
!  f -> X^2*f
!
!-----------------------------------------------------------------------------
!
  subroutine mult_x2(ex,X,f)
  implicit none
  integer, dimension(3)			:: ex
  real*8, dimension(ex(1),ex(2),ex(3))  :: X,f
  integer                               :: i,j,k
  integer                               :: imin, jmin, kmin, imax, jmax, kmax
!
  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)
  imax = ubound(f,1)
  jmax = ubound(f,2)
  kmax = ubound(f,3)
!
     do k=kmin,kmax
        do j=jmin,jmax
           do i=imin,imax
              f(i,j,k)=X(i,1,1)*X(i,1,1)*f(i,j,k)
           end do
        end do
     end do
end subroutine mult_x2

!
!-----------------------------------------------------------------------------
!
!  f1 -> f1/X, f2 -> f2/X
!
!-----------------------------------------------------------------------------
!
subroutine div_x(ex,X,f1,f2)
  implicit none
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin, imax, jmax, kmax
  integer, dimension(3)              :: ex
  real*8, dimension(ex(1),ex(2),ex(3)) :: X,f1,f2
!
  imin = lbound(f1,1)
  jmin = lbound(f1,2)
  kmin = lbound(f1,3)
  imax = ubound(f1,1)
  jmax = ubound(f1,2)
  kmax = ubound(f1,3)
!
     do k=kmin,kmax
        do j=jmin,jmax
           do i=imin,imax
              f1(i,j,k)=f1(i,j,k)/X(i,1,1)
              f2(i,j,k)=f2(i,j,k)/X(i,1,1)
           end do
        end do
     end do
end subroutine div_x

!
!-----------------------------------------------------------------------------
!
!  f1 -> f1/X^2, f2 -> f2/X^2
!
!-----------------------------------------------------------------------------
!
subroutine div_x2(ex,X,f1,f2)
  implicit none
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin, imax, jmax, kmax
  integer, dimension(3)              :: ex
  real*8, dimension(ex(1),ex(2),ex(3)) :: X,f1,f2
!
  imin = lbound(f1,1)
  jmin = lbound(f1,2)
  kmin = lbound(f1,3)
  imax = ubound(f1,1)
  jmax = ubound(f1,2)
  kmax = ubound(f1,3)
!
     do k=kmin,kmax
        do j=jmin,jmax
           do i=imin,imax
              f1(i,j,k)=f1(i,j,k)/X(i,1,1)/X(i,1,1)
              f2(i,j,k)=f2(i,j,k)/X(i,1,1)/X(i,1,1)
           end do
        end do
     end do
end subroutine div_x2

!
!-----------------------------------------------------------------------------
!
!  S_rhs -> S_rhs + sign * Si * vj / X
!
!-----------------------------------------------------------------------------
!
subroutine addsvox(ex,X,Si,vj,S_rhs,sign)
  implicit none
  integer , dimension(3)              :: ex
  real*8, dimension(ex(1),ex(2),ex(3)) ::X,Si,vj,S_rhs
  real*8 			      :: sign
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin, imax, jmax, kmax
!
  imin = lbound(Si,1)
  jmin = lbound(Si,2)
  kmin = lbound(Si,3)
  imax = ubound(Si,1)
  jmax = ubound(Si,2)
  kmax = ubound(Si,3)
!
     do k=kmin,kmax
        do j=jmin,jmax
           do i=imin,imax
              S_rhs(i,j,k)=S_rhs(i,j,k)+sign*Si(i,j,k)*vj(i,j,k)/X(i,1,1)
           end do
        end do
     end do
end subroutine addSvox
