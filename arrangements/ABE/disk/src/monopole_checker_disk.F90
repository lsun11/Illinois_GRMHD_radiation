!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Look for magnetic monopoles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_absdivb_disk(ext,dX,dY,dZ,X,Z,Bx,By,Bz,phi,absdivB,Sym_Bz,Symmetry,x_shell)
  implicit none
  integer, dimension(3)				:: ext
  real*8, dimension(ext(1),ext(2),ext(3))        :: X,Z
  real*8, dimension(ext(1),ext(2),ext(3))        :: Bx,By,Bz,phi, absdivB
  integer					:: i,imin,imax,j,jmin,jmax
  integer                                        :: k,kmin,kmax,Symmetry
  real*8						:: dX,dY,dZ,f1odX,f1odY,f1odZ
  real*8						:: Sym_Bz, Zmin
  real*8, dimension(ext(1),ext(2),ext(3))        :: Bxx,Byy,Bzz,tmp,e6phi
  integer                                        :: i_bd
  real*8                                         :: x_shell
  logical                                        :: have_inner_boundary
  !
  imin = lbound(Bx,1)+1
  imax = ubound(Bx,1)
  jmin = lbound(Bx,2)
  jmax = ubound(Bx,2)
  kmin = lbound(Bx,3)+1
  kmax = ubound(Bx,3)
  f1odX = 1.d0/(dX)
  f1odY = 1.d0/(dY)
  f1odZ = 1.d0/(dZ)
  Zmin = Z(imin,jmin,kmin+1)

  e6phi = exp(6.d0*phi)
  Bx = Bx*e6phi
  By = By*e6phi
  Bz = Bz*e6phi

  if(Symmetry==4) then

     !
     ! Bxx
     tmp(imin,:,:) = 0.d0
     do i=imin+1,imax
        tmp(i,:,:) = f1odX*(X(i,:,:)*Bx(i,:,:) - X(i-1,:,:)*Bx(i-1,:,:))
     end do
     ! average in Z
     if (Zmin .gt. 0.d0) then
        Bxx(:,:,kmin) = tmp(:,:,kmin)*(1.d0 - Sym_Bz)/2.d0
     else
        Bxx(:,:,kmin) = 0.d0
     end if
     do k=kmin+1,kmax
        Bxx(:,:,k) = (tmp(:,:,k) + tmp(:,:,k-1))/2.d0
     end do
     !
     ! Bzz
     do i=imin,imax
        if (Zmin .gt. 0.d0) then
           tmp(i,:,kmin) = X(i,:,kmin)*f1odZ*(1.d0-Sym_Bz)*Bz(i,:,kmin)
        else
           tmp(i,:,kmin) = 0.d0
        end if
     end do
     do k=kmin+1,kmax
        do i=imin,imax
           tmp(i,:,k) = X(i,:,k)*f1odZ*(Bz(i,:,k) - Bz(i,:,k-1))
        end do
     end do
     ! average in X
     Bzz(imin,:,:) = 0.d0
     do i=imin+1,imax
        Bzz(i,:,:) = (tmp(i,:,:) + tmp(i-1,:,:))/2.d0
     end do
     absdivB = abs(Bxx + Bzz)

  else if(Symmetry==1) then

     !
     ! Bxx
     Bxx(imin,:,:) = 0.d0
     do i=imin+1,imax
        Bxx(i,:,:) = f1odX*(Bx(i,:,:) - Bx(i-1,:,:))
     end do
     ! average in Y
     tmp(:,jmin,:) = 0.d0
     do j=jmin+1,jmax
        tmp(:,j,:) = (Bxx(:,j,:) + Bxx(:,j-1,:))/2.d0
     end do
     ! average in Z
     Bxx(:,:,kmin) = tmp(:,:,kmin)*(1.d0 - Sym_Bz)/2.d0
     do k=kmin+1,kmax
        Bxx(:,:,k) = (tmp(:,:,k) + tmp(:,:,k-1))/2.d0
     end do

     !
     ! Byy
     Byy(:,jmin,:) = 0.d0
     do j=jmin+1,jmax
        Byy(:,j,:) = f1odY*(By(:,j,:) - By(:,j-1,:))
     end do
     ! average in X
     tmp(imin,:,:) = 0.d0
     do i=imin+1,imax
        tmp(i,:,:) = (Byy(i,:,:) + Byy(i-1,:,:))/2.d0
     end do
     ! average in Z
     Byy(:,:,kmin) = tmp(:,:,kmin)*(1.d0 - Sym_Bz)/2.d0
     do k=kmin+1,kmax
        Byy(:,:,k) = (tmp(:,:,k) + tmp(:,:,k-1))/2.d0
     end do

     !
     ! Bzz
     Bzz(:,:,kmin) = f1odZ*(1.d0 - Sym_Bz)*Bz(:,:,kmin)
     do k=kmin,kmax
        Bzz(:,:,k) = f1odZ*(Bz(:,:,k) - Bz(:,:,k-1))
     end do
     ! average in X
     tmp(imin,:,:) = 0.d0
     do i=imin+1,imax
        tmp(i,:,:) = (Bzz(i,:,:) + Bzz(i-1,:,:))/2.d0
     end do
     ! average in Y
     Bzz(:,jmin,:) = 0.d0
     do j=jmin+1,jmax
        Bzz(:,j,:) = (tmp(:,j,:) + tmp(:,j-1,:))/2.d0
     end do

     absdivB = abs(Bxx + Byy + Bzz)

  else if (Symmetry==0) then

     !
     ! Bxx
     Bxx(imin,:,:) = 0.d0
     do i=imin+1,imax
        Bxx(i,:,:) = f1odX*(Bx(i,:,:) - Bx(i-1,:,:))
     end do
     ! average in Y
     tmp(:,jmin,:) = 0.d0
     do j=jmin+1,jmax
        tmp(:,j,:) = (Bxx(:,j,:) + Bxx(:,j-1,:))/2.d0
     end do
     ! average in Z
     Bxx(:,:,kmin) = 0.d0
     do k=kmin+1,kmax
        Bxx(:,:,k) = (tmp(:,:,k) + tmp(:,:,k-1))/2.d0
     end do

     !
     ! Byy
     Byy(:,jmin,:) = 0.d0
     do j=jmin+1,jmax
        Byy(:,j,:) = f1odY*(By(:,j,:) - By(:,j-1,:))
     end do
     ! average in X
     tmp(imin,:,:) = 0.d0
     do i=imin+1,imax
        tmp(i,:,:) = (Byy(i,:,:) + Byy(i-1,:,:))/2.d0
     end do
     ! average in Z
     Byy(:,:,kmin) = 0.d0
     do k=kmin+1,kmax
        Byy(:,:,k) = (tmp(:,:,k) + tmp(:,:,k-1))/2.d0
     end do

     !
     ! Bzz
     Bzz(:,:,kmin) = 0.d0
     do k=kmin,kmax
        Bzz(:,:,k) = f1odZ*(Bz(:,:,k) - Bz(:,:,k-1))
     end do
     ! average in X
     tmp(imin,:,:) = 0.d0
     do i=imin+1,imax
        tmp(i,:,:) = (Bzz(i,:,:) + Bzz(i-1,:,:))/2.d0
     end do
     ! average in Y
     Bzz(:,jmin,:) = 0.d0
     do j=jmin+1,jmax
        Bzz(:,j,:) = (tmp(:,j,:) + tmp(:,j-1,:))/2.d0
     end do

     absdivB = abs(Bxx + Byy + Bzz)

  else
     write(*,*) "Warning:  monopole locator not implemented for this symmetry"

  end if

  Bx = Bx/e6phi
  By = By/e6phi
  Bz = Bz/e6phi

end subroutine compute_absdivb_disk
