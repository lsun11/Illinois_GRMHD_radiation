
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Look for magnetic monopoles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_absdivb(ext,dX,dY,dZ,X,Z,Bx,By,Bz,phi,absdivB,Sym_Bz,Symmetry, &
     Bxx,Byy,Bzz,tmp,e6phi)
 implicit none
 integer, dimension(3)				:: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Z
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bx,By,Bz,phi, absdivB
 integer					:: i,imin,imax,j,jmin,jmax
 integer                                        :: k,kmin,kmax,Symmetry
 real*8						:: dX,dY,dZ,f1odX,f1odY,f1odZ
 real*8						:: Sym_Bz, Zmin
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bxx,Byy,Bzz,tmp,e6phi
!
 imin = lbound(Bx,1)
 imax = ubound(Bx,1)
 jmin = lbound(Bx,2)
 jmax = ubound(Bx,2)
 kmin = lbound(Bx,3)
 kmax = ubound(Bx,3)
 f1odX = 1.d0/(dX)
 f1odY = 1.d0/(dY)
 f1odZ = 1.d0/(dZ)

 if(Symmetry==1 .and. Z(1,1,1).lt.0.D0) kmin = kmin + 1
 if(Symmetry==4) then
    if(X(1,1,1).lt.0.D0) imin = imin + 1
    if(Z(1,1,1).lt.0.D0) kmin = kmin + 1
 end if
 Zmin = Z(kmin,1,1)
 
 e6phi = exp(6.d0*phi)
 Bx = Bx*e6phi
 By = By*e6phi
 Bz = Bz*e6phi

 if(Symmetry==4) then

    !
    ! Bxx
    tmp(imin,:,:) = 0.d0
    do i=imin+1,imax
       tmp(i,:,:) = f1odX*(X(i,1,1)*Bx(i,:,:) - X(i-1,1,1)*Bx(i-1,:,:))
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
          tmp(i,:,kmin) = X(i,1,1)*f1odZ*(1.d0-Sym_Bz)*Bz(i,:,kmin)
       else
	  tmp(i,:,kmin) = 0.d0
       end if
    end do
    do k=kmin+1,kmax
       do i=imin,imax
          tmp(i,:,k) = X(i,1,1)*f1odZ*(Bz(i,:,k) - Bz(i,:,k-1))
       end do
    end do
    ! average in X
    Bzz(imin,:,:) = 0.d0
    do i=imin+1,imax
       Bzz(i,:,:) = (tmp(i,:,:) + tmp(i-1,:,:))/2.d0
    end do
    absdivB = abs(Bxx + Bzz)

! *** TEST ***
!    absdivB(imin,:,:) = 0.d0
    absdivB(imax,:,:) = 0.d0
    absdivB(:,jmin,:) = 0.d0
    absdivB(:,jmax,:) = 0.d0
!    absdivB(:,:,kmin) = 0.d0
    absdivB(:,:,kmax) = 0.d0
! ************

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
    do k=kmin+1,kmax
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

! *** TEST ***
!    absdivB(imin,:,:) = 0.d0
    absdivB(imax,:,:) = 0.d0
!    absdivB(:,jmin,:) = 0.d0
    absdivB(:,jmax,:) = 0.d0
    absdivB(:,:,kmin) = 0.d0
    absdivB(:,:,kmax) = 0.d0
! ************

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
    do k=kmin+1,kmax
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

end subroutine compute_absdivb


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Look for magnetic monopoles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_absdivb2(ext,dX,dY,dZ,X,Z,Bx,By,Bz,phi,absdivB,Sym_Bz,Symmetry, &
     Bxx,Byy,Bzz,tmp,e6phi)
 implicit none
 integer, dimension(3)				:: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Z
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bx,By,Bz,phi, absdivB
 integer					:: i,imin,imax,j,jmin,jmax
 integer                                        :: k,kmin,kmax,Symmetry
 real*8						:: dX,dY,dZ,f1odX,f1odY,f1odZ
 real*8						:: Sym_Bz, Zmin
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bxx,Byy,Bzz,tmp,e6phi
!
 imin = lbound(Bx,1)
 imax = ubound(Bx,1)
 jmin = lbound(Bx,2)
 jmax = ubound(Bx,2)
 kmin = lbound(Bx,3)
 kmax = ubound(Bx,3)
 f1odX = 1.d0/(dX)
 f1odY = 1.d0/(dY)
 f1odZ = 1.d0/(dZ)

 if(Symmetry==1 .and. Z(1,1,1).lt.0.D0) kmin = kmin + 1
 if(Symmetry==4) then
    if(X(1,1,1).lt.0.D0) imin = imin + 1
    if(Z(1,1,1).lt.0.D0) kmin = kmin + 1
 end if
 Zmin = Z(kmin,1,1)
 
 e6phi = exp(6.d0*phi)
 Bx = Bx*e6phi
 By = By*e6phi
 Bz = Bz*e6phi


    !
    ! Bxx
    Bxx(imin,:,:) = 0.d0
    do i=imin+1,imax-1
       Bxx(i,:,:) = f1odX*(Bx(i+1,:,:) - Bx(i-1,:,:))*0.5D0
    end do
    Bxx(imin,:,:) = f1odX*(Bx(imin+1,:,:) - Bx(imin,:,:))*0.5D0

    !
    ! Byy
    Byy(:,jmin,:) = 0.d0
    do j=jmin+1,jmax-1
       Byy(:,j,:) = f1odY*(By(:,j+1,:) - By(:,j-1,:))*0.5D0
    end do
    Byy(:,jmin,:) = f1odY*(By(:,jmin+1,:) - By(:,jmin,:))*0.5D0

    !
    ! Bzz
    Bzz(:,:,kmin) = 0.d0
    do k=kmin+1,kmax-1
       Bzz(:,:,k) = f1odZ*(Bz(:,:,k+1) - Bz(:,:,k-1))*0.5D0
    end do
    Bzz(:,:,kmin) = f1odZ*(Bz(:,:,kmin+1) - Sym_Bz*Bz(:,:,kmin))*0.5D0

    absdivB = abs(Bxx + Byy + Bzz)

! *** TEST ***
    absdivB(imin,:,:) = 0.d0
    absdivB(imax,:,:) = 0.d0
    absdivB(:,jmin,:) = 0.d0
    absdivB(:,jmax,:) = 0.d0
!    absdivB(:,:,kmin) = 0.d0
    absdivB(:,:,kmax) = 0.d0
! ************ 

 Bx = Bx/e6phi
 By = By/e6phi
 Bz = Bz/e6phi

end subroutine compute_absdivb2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Look for magnetic monopoles (for CT scheme of Del Zanna)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_absdivb_ct3_old(ext,dX,dY,dZ,Bxtilde,Bytilde,Bztilde,absdivB)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bxtilde,Bytilde,Bztilde,absdivB
 integer                                        :: i,j,k
 real*8                                         :: dX,dY,dZ
!
 !$omp parallel do 
 do k=2,ext(3)-1
    do j=2,ext(2)-1
       do i=2,ext(1)-1
	  absdivb(i,j,k) = abs( ( Bxtilde(i,j,k) - Bxtilde(i-1,j,k) ) / dX + &
			   ( Bytilde(i,j,k) - Bytilde(i,j-1,k) ) / dY + & 
			   ( Bztilde(i,j,k) - Bztilde(i,j,k-1) ) / dZ )
       end do
    end do
 end do
 !$omp end parallel do
 
 !$omp parallel do
 do k=1,ext(3)
    do j=1,ext(2)
       absdivb(1,j,k) = 0.d0
       absdivb(ext(1),j,k) = 0.d0
    end do
 end do
 !$omp end parallel do

 !$omp parallel do
 do k=1,ext(3)
    do i=1,ext(1)
       absdivb(i,1,k) = 0.d0
       absdivb(i,ext(2),k) = 0.d0
    end do
 end do
 !$omp end parallel do

 !$omp parallel do
 do j=1,ext(2)
    do i=1,ext(1)
       absdivb(i,j,1) = 0.d0
       absdivb(i,j,ext(3)) = 0.d0
    end do
 end do
 !$omp end parallel do

end subroutine compute_absdivb_ct3_old

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Look for magnetic monopoles (for CT scheme of Del Zanna)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_absdivb_ct3(ext,dX,dY,dZ,phi,Bx,By,Bz, & 
                               Bxtilde,Bytilde,Bztilde,absdivB)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: phi,Bx,By,Bz,absdivB
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bxtilde,Bytilde,Bztilde
 integer                                        :: i,j,k
 real*8                                         :: dX,dY,dZ
!
 !$omp parallel do
 do k=2,ext(3)-1
    do j=2,ext(2)-1
       do i=2,ext(1)-1
	  Bxtilde(i,j,k) = Bx(i,j,k) * exp(3.d0*(phi(i,j,k)+phi(i+1,j,k)))
          Bytilde(i,j,k) = By(i,j,k) * exp(3.d0*(phi(i,j,k)+phi(i,j+1,k)))
          Bztilde(i,j,k) = Bz(i,j,k) * exp(3.d0*(phi(i,j,k)+phi(i,j,k+1)))
       end do
    end do
 end do
 !$omp end parallel do

 !$omp parallel do 
 do k=3,ext(3)-1
    do j=3,ext(2)-1
       do i=3,ext(1)-1
	  absdivb(i,j,k) = abs( ( Bxtilde(i,j,k) - Bxtilde(i-1,j,k) ) / dX + &
			   ( Bytilde(i,j,k) - Bytilde(i,j-1,k) ) / dY + & 
			   ( Bztilde(i,j,k) - Bztilde(i,j,k-1) ) / dZ )
       end do
    end do
 end do
 !$omp end parallel do

 !$omp parallel do
 do k=1,ext(3)
    do j=1,ext(2)
       absdivb(1,j,k) = 0.d0
       absdivb(2,j,k) = 0.d0
       absdivb(ext(1),j,k) = 0.d0
    end do
 end do
 !$omp end parallel do

 !$omp parallel do
 do k=1,ext(3)
    do i=1,ext(1)
       absdivb(i,1,k) = 0.d0
       absdivb(i,2,k) = 0.d0
       absdivb(i,ext(2),k) = 0.d0
    end do
 end do
 !$omp end parallel do

 !$omp parallel do
 do j=1,ext(2)
    do i=1,ext(1)
       absdivb(i,j,1) = 0.d0
       absdivb(i,j,2) = 0.d0
       absdivb(i,j,ext(3)) = 0.d0
    end do
 end do
 !$omp end parallel do

end subroutine compute_absdivb_ct3
