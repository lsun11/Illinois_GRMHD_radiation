!-----------------------------------------------------------!
!-----------------------------------------------------------!
subroutine find_ur_ul(ext,X,Y,Z,U,Ur,Ul,nabla_U,m,Sym_arg,Reconstruction,Symmetry)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8, dimension(ext(1),ext(2),ext(3))        :: U,Ur,Ul,nabla_U
 integer                                        :: m, i,j,k,mmax,mmin
 integer                                        :: Reconstruction,SPPM
 real*8                                         :: Sym_arg
 integer                            :: Symmetry, NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
 parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
 parameter (SPPM = 5)
!
 if (Reconstruction == SPPM) then
    call find_ur_ul_sppm(ext,X,Y,Z,U,Ur,Ul,nabla_U,m)
    return
 end if
 mmin = lbound(U,m)
 mmax = ubound(U,m)
 if(Symmetry == OCTANT) then
    if(m==1 .and. X(1,1,1) .lt. 0.D0) then
       mmin = lbound(U,1)+1
    else if(m==2 .and. Y(1,1,1) .lt. 0.D0) then
       mmin = lbound(U,2)+1
    else if(m==3 .and. Z(1,1,1) .lt. 0.D0) then
       mmin = lbound(U,3)+1
    end if
 else if(Symmetry == AXISYM) then
    if(m==1 .and. X(1,1,1) .lt. 0.D0) then
       mmin = lbound(U,1)+1
    else if(m==3 .and. Z(1,1,1) .lt. 0.D0) then
       mmin = lbound(U,3)+1
    end if
 else if(Symmetry == EQUATORIAL) then
    if(m==3 .and. Z(1,1,1) .lt. 0.D0) then
       mmin = lbound(U,3)+1
    end if
 end if
 Ur = U - 0.5d0*nabla_U
 if (m==1) then
    Ul(mmin,:,:) = Sym_arg*Ur(mmin,:,:)
!    i = mmin
!    if (Sym_arg > 0.d0) then
!       Ul(i,:,:) = Ur(i,:,:)
!    else
!       Ul(i,:,:) = 0.d0
!       Ur(i,:,:) = 0.d0
!    end if
    do i = mmin+1,mmax
       Ul(i,:,:) = U(i-1,:,:) + 0.5d0*nabla_U(i-1,:,:)
    end do
    return
 end if
 if (m==2) then
    Ul(:,mmin,:) = Sym_arg*Ur(:,mmin,:)
!    j = mmin
!    if (Sym_arg > 0.d0) then
!       Ul(:,j,:) = Ur(:,j,:)
!    else 
!       Ul(:,j,:) = 0.d0
!       Ur(:,j,:) = 0.d0
!    end if
    do j = mmin+1,mmax
       Ul(:,j,:) = U(:,j-1,:) + 0.5d0*nabla_U(:,j-1,:)
    end do
    return
 end if
 if (m==3) then
    Ul(:,:,mmin) = Sym_arg*Ur(:,:,mmin)
!    k = mmin
!    if (Sym_arg > 0.d0) then
!       Ul(:,:,k) = Ur(:,:,k)
!    else
!       Ul(:,:,k) = 0.d0
!       Ur(:,:,k) = 0.d0
!    end if
    do k = mmin+1,mmax
       Ul(:,:,k) = U(:,:,k-1) + 0.5d0*nabla_U(:,:,k-1)
    end do
    return
 end if
end subroutine find_ur_ul
!--------------------------------------------------------------
! Note that in this subroutine, Ur is evaluated at m-1/2, 
! whereas Ul is evaluated at m+1/2
!--------------------------------------------------------------
subroutine find_ur_ul_sppm(ext,X,Y,Z,U,Ur,Ul,dU,m,Symmetry)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8, dimension(ext(1),ext(2),ext(3))        :: U,Ur,Ul,dU
 real*8, dimension(ext(1),ext(2),ext(3))        :: dpU, dmU
 integer                                        :: m, i,j,k,mmax,mmin
 real*8, parameter        :: third = 1.d0/3.d0, sixth = 1.d0/6.d0, zero = 0.d0
 real*8, parameter      :: b = 2.d0
 integer                            :: Symmetry, NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
 parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
!
 mmin = lbound(U,m)
 mmax = ubound(U,m)
!!$ if(Symmetry == OCTANT) then 
!!$    if(m==1 .and. X(1,1,1) .lt. 0.D0) then
!!$       mmin = lbound(U,1)+1
!!$    else if(m==2 .and. Y(1,1,1) .lt. 0.D0) then
!!$       mmin = lbound(U,2)+1
!!$    else if(m==3 .and. Z(1,1,1) .lt. 0.D0) then
!!$       mmin = lbound(U,3)+1     
!!$    end if
!!$ else if(Symmetry == AXISYM) then
!!$    if(m==1 .and. X(1,1,1) .lt. 0.D0) then
!!$       mmin = lbound(U,1)+1
!!$    else if(m==3 .and. Z(1,1,1) .lt. 0.D0) then
!!$       mmin = lbound(U,3)+1     
!!$    end if
!!$ else if(Symmetry == EQUATORIAL) then
!!$    if(m==3 .and. Z(1,1,1) .lt. 0.D0) then
!!$       mmin = lbound(U,3)+1     
!!$    end if
!!$ end if
 if (m==1) then
    do i=mmin,mmax-1
       where (dU(i,:,:)*dU(i+1,:,:) .le. zero)
          dpU(i,:,:) = zero
          dmU(i,:,:) = zero
       elsewhere
          dpU(i,:,:) = min( b*abs(dU(i+1,:,:)), abs(dU(i,:,:)) )
          dmU(i,:,:) = min( b*abs(dU(i,:,:)), abs(dU(i+1,:,:)) )
          where (dU(i,:,:) .lt. zero)
            dpU(i,:,:) = -dpU(i,:,:)
            dmU(i,:,:) = -dmU(i,:,:)
          end where
       end where
    end do
    dpU(mmax,:,:) = dpU(mmax-1,:,:)
    dmU(mmax,:,:) = dmU(mmax-1,:,:)
    Ur = U - third*dpU - sixth*dmU
    Ul = U + sixth*dpU + third*dmU
    return
 end if
 if (m==2) then
    do j=mmin,mmax-1
       where (dU(:,j,:)*dU(:,j+1,:) .le. zero)
          dpU(:,j,:) = zero
          dmU(:,j,:) = zero
       elsewhere
          dpU(:,j,:) = min( b*abs(dU(:,j+1,:)), abs(dU(:,j,:)) )
          dmU(:,j,:) = min( b*abs(dU(:,j,:)), abs(dU(:,j+1,:)) )
          where (dU(:,j,:) .lt. zero)
            dpU(:,j,:) = -dpU(:,j,:)
            dmU(:,j,:) = -dmU(:,j,:)
          end where
       end where
    end do
    dpU(:,mmax,:) = dpU(:,mmax-1,:)
    dmU(:,mmax,:) = dmU(:,mmax-1,:)
    Ur = U - third*dpU - sixth*dmU
    Ul = U + sixth*dpU + third*dmU
    return
 end if
 if (m==3) then
    do k=mmin,mmax-1
       where (dU(:,:,k)*dU(:,:,k+1) .le. zero)
          dpU(:,:,k) = zero
          dmU(:,:,k) = zero
       elsewhere
          dpU(:,:,k) = min( b*abs(dU(:,:,k+1)), abs(dU(:,:,k)) )
          dmU(:,:,k) = min( b*abs(dU(:,:,k)), abs(dU(:,:,k+1)) )
          where (dU(:,:,k) .lt. zero)
            dpU(:,:,k) = -dpU(:,:,k)
            dmU(:,:,k) = -dmU(:,:,k)
          end where
       end where
    end do
    dpU(:,:,mmax) = dpU(:,:,mmax-1)
    dmU(:,:,mmax) = dmU(:,:,mmax-1)
    Ur = U - third*dpU - sixth*dmU
    Ul = U + sixth*dpU + third*dmU
    return
 end if
end subroutine find_ur_ul_sppm
!--------------------------------------------------!
! 
!--------------------------------------------------!
!
subroutine compute_nabla(ext,dX,dY,dZ,U,nabla_U,m,Sym_arg,X,Y,Z,LIMITER,Symmetry)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: U,nabla_U
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 integer                                        :: m
 integer                                        :: LIMITER
 integer                                        :: ZERO,MC,VAN_LEER,MINMOD,dU
 real*8                                         :: Sym_arg
 real*8                                         :: dX, dY, dZ
 integer                            :: Symmetry, NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
 parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
 parameter(ZERO=0, MC=1, VAN_LEER=2, MINMOD=3, dU = 4)
!
 if (LIMITER==ZERO) then
    nabla_U = 0.d0
    return
 elseif (LIMITER==MC) then
    call compute_nabla_mc(ext,dX,dY,dZ,U,nabla_U,m,Sym_arg,X,Y,Z,Symmetry)
    return
 elseif (LIMITER==VAN_LEER) then
    call compute_nabla_van_leer(ext,U,nabla_U,m,Sym_arg,X,Y,Z,Symmetry)
    return
 elseif (LIMITER==MINMOD) then
    call compute_nabla_minmod(ext,U,nabla_U,m,Sym_arg,X,Y,Z,Symmetry)
    return
 elseif (LIMITER==dU) then
    call compute_du(ext,U,nabla_U,m,Sym_arg,X,Y,Z,Symmetry)
    return
 end if
end subroutine compute_nabla
subroutine compute_du(ext,U,dU,m,Sym_arg,X,Y,Z,Symmetry)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: U,dU
 integer                                        :: m,i,j,k,mmax,mmin
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8                                         :: Sym_arg
 integer                            :: Symmetry, NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
 parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
!
 mmax = ubound(U,m)
 mmin = lbound(U,m)
 if(Symmetry == OCTANT) then
    if(m==1 .and. X(1,1,1).lt.0.D0) then
       mmin = lbound(U,1)+1
    else if(m==2 .and. Y(1,1,1).lt.0.D0) then
       mmin = lbound(U,2)+1
    else if(m==3 .and. Z(1,1,1).lt.0.D0) then
       mmin = lbound(U,3)+1
    end if
 else if(Symmetry == AXISYM) then
    if(m==1 .and. X(1,1,1).lt.0.D0) then
       mmin = lbound(U,1)+1
    else if(m==3 .and. Z(1,1,1).lt.0.D0) then
       mmin = lbound(U,3)+1
    end if
 else if(Symmetry == EQUATORIAL) then
    if(m==3 .and. Z(1,1,1).lt.0.D0) then
       mmin = lbound(U,3)+1
    end if
 end if
 if (m==1) then
    do i = mmin+1,mmax
       dU(i,:,:) = U(i,:,:) - U(i-1,:,:)
    end do
    dU(mmin,:,:) = (1.d0-Sym_arg)*U(mmin,:,:)
    return
 end if
 if (m==2) then
    do j = mmin+1,mmax
       dU(:,j,:) = U(:,j,:) - U(:,j-1,:)
    end do
    dU(:,mmin,:) = (1.d0-Sym_arg)*U(:,mmin,:)
    return
 end if
 if (m==3) then
    do k = mmin+1,mmax
       dU(:,:,k) = U(:,:,k) - U(:,:,k-1)
    end do
    if (Z(1,1,mmin) .ge. 0.d0) then
       dU(:,:,mmin) = (1.d0-Sym_arg)*U(:,:,mmin)
    else
       dU(:,:,mmin) = dU(:,:,mmin+1)
    end if
    return
 end if
end subroutine compute_du
!--------------------------------------------------!
! Minmod limiter                                   !
!--------------------------------------------------!
!
subroutine compute_nabla_minmod(ext,U,nabla_U,m,Sym_arg,X,Y,Z,Symmetry)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: U,nabla_U
 integer                                        :: m,i,j,k,mmax,mmin
 real*8, dimension(ext(1),ext(2),ext(3))        :: dU
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8                                         :: Sym_arg
 integer                            :: Symmetry, NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
 parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
!
 mmax = ubound(U,m)
 mmin = lbound(U,m)
 if(Symmetry == OCTANT) then
    if(m==1 .and. X(1,1,1).lt.0.D0) then
       mmin = lbound(U,1)+1
    else if(m==2 .and. Y(1,1,1).lt.0.D0) then
       mmin = lbound(U,2)+1
    else if(m==3 .and. Z(1,1,1).lt.0.D0) then
       mmin = lbound(U,3)+1
    end if
 else if(Symmetry == AXISYM) then
    if(m==1 .and. X(1,1,1).lt.0.D0) then
       mmin = lbound(U,1)+1
    else if(m==3 .and. Z(1,1,1).lt.0.D0) then
       mmin = lbound(U,3)+1
    end if
 else if(Symmetry == EQUATORIAL) then
    if(m==3 .and. Z(1,1,1).lt.0.D0) then
       mmin = lbound(U,3)+1
    end if
 end if
 if (m==1) then
    do i = mmin+1,mmax
       dU(i,:,:) = U(i,:,:) - U(i-1,:,:)
    end do
    if (Sym_arg < 0.d0) then
       dU(mmin,:,:) = 2.d0 * U(mmin,:,:)
    else
       dU(mmin,:,:) = 0.d0
    end if
    do i = mmin,mmax-1
       where (dU(i,:,:)*dU(i+1,:,:) .gt. 0.d0)
          where (abs(dU(i,:,:)) .lt. abs(dU(i+1,:,:)))
            nabla_U(i,:,:) = dU(i,:,:)
          elsewhere
            nabla_U(i,:,:) = dU(i+1,:,:)
          end where
       elsewhere
          nabla_U(i,:,:) = 0.d0
       end where
    end do
    nabla_U(mmax,:,:) = dU(mmax,:,:)
    return
 end if
 if (m==2) then
    do j = mmin+1,mmax
       dU(:,j,:) = U(:,j,:) - U(:,j-1,:)
    end do
    if (Sym_arg < 0.d0) then
       dU(:,mmin,:) = 2.d0 * U(:,mmin,:)
    else
       dU(:,mmin,:) = 0.d0
    end if
    do j = mmin,mmax-1
       where (dU(:,j,:)*dU(:,j+1,:) .gt. 0.d0)
          where (abs(dU(:,j,:)) .lt. abs(dU(:,j+1,:)))
            nabla_U(:,j,:) = dU(:,j,:)
          elsewhere
            nabla_U(:,j,:) = dU(:,j+1,:)
          end where
       elsewhere
          nabla_U(:,j,:) = 0.d0
       end where
    end do
    nabla_U(:,mmax,:) = dU(:,mmax,:)
    return
 end if
 if (m==3) then
    do k = mmin+1,mmax
       dU(:,:,k) = U(:,:,k) - U(:,:,k-1)
    end do
    if (Sym_arg < 0.d0) then
       dU(:,:,mmin) = 2.d0 * U(:,:,mmin)
    else
       dU(:,:,mmin) = 0.d0
    end if
    do k = mmin,mmax-1
       where (dU(:,:,k)*dU(:,:,k+1) .gt. 0.d0)
          where (abs(dU(:,:,k)) .lt. abs(dU(:,:,k+1)))
            nabla_U(:,:,k) = dU(:,:,k)
          elsewhere
            nabla_U(:,:,k) = dU(:,:,k+1)
          end where
       elsewhere
          nabla_U(:,:,k) = 0.d0
       end where
    end do
    nabla_U(:,:,mmax) = dU(:,:,mmax)
    if (Z(1,1,mmin) .lt. 0.d0 .or. Symmetry==NO_SYMM) nabla_U(:,:,mmin) = dU(:,:,mmin+1)
    return
 end if
end subroutine compute_nabla_minmod
!--------------------------------------------------!
! van Leer limiter                                 !
!--------------------------------------------------!
!
subroutine compute_nabla_van_leer(ext,U,nabla_U,m,Sym_arg,X,Y,Z,Symmetry)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: U,nabla_U
 integer                                        :: m,i,j,k,mmax,mmin
 real*8, dimension(ext(1),ext(2),ext(3))        :: dU
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8                                         :: Sym_arg
 integer                            :: Symmetry, NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
 parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
!
 mmax = ubound(U,m)
 mmin = lbound(U,m)
 if(Symmetry == OCTANT) then
    if(m==1 .and. X(1,1,1).lt.0.D0) then
       mmin = lbound(U,1)+1
    else if(m==2 .and. Y(1,1,1).lt.0.D0) then
       mmin = lbound(U,2)+1
    else if(m==3 .and. Z(1,1,1).lt.0.D0) then
       mmin = lbound(U,3)+1
    end if
 else if(Symmetry == AXISYM) then
    if(m==1 .and. X(1,1,1).lt.0.D0) then
       mmin = lbound(U,1)+1
    else if(m==3 .and. Z(1,1,1).lt.0.D0) then
       mmin = lbound(U,3)+1
    end if
 else if(Symmetry == EQUATORIAL) then
    if(m==3 .and. Z(1,1,1).lt.0.D0) then
       mmin = lbound(U,3)+1
    end if
 end if
 if (m==1) then
    do i = mmin+1,mmax
       dU(i,:,:) = U(i,:,:) - U(i-1,:,:)
    end do
    if (Sym_arg < 0.d0) then
       dU(mmin,:,:) = 2.d0 * U(mmin,:,:)
    else
       dU(mmin,:,:) = 0.d0
    end if
    do i = mmin,mmax-1
       where (dU(i,:,:)*dU(i+1,:,:) .gt. 0.d0)
         nabla_U(i,:,:) = 2.d0*dU(i,:,:)/(dU(i,:,:)+dU(i+1,:,:))*dU(i+1,:,:)
       elsewhere
         nabla_U(i,:,:) = 0.d0
       end where
    end do
    nabla_U(mmax,:,:) = dU(mmax,:,:)
 return
 end if
 if (m==2) then
    do j = mmin+1,mmax
       dU(:,j,:) = U(:,j,:) - U(:,j-1,:)
    end do
    if (Sym_arg < 0.d0) then
       dU(:,mmin,:) = 2.d0 * U(:,mmin,:)
    else
       dU(:,mmin,:) = 0.d0
    end if
    do j = mmin,mmax-1
       where (dU(:,j,:)*dU(:,j+1,:) .gt. 0.d0)
         nabla_U(:,j,:) = 2.d0*dU(:,j,:)/(dU(:,j,:)+dU(:,j+1,:))*dU(:,j+1,:)
       elsewhere
         nabla_U(:,j,:) = 0.d0
       end where
    end do
    nabla_U(:,mmax,:) = dU(:,mmax,:)
 return
 end if
 if (m==3) then
    do k = mmin+1,mmax
       dU(:,:,k) = U(:,:,k) - U(:,:,k-1)
    end do
    if (Sym_arg < 0.d0) then
       dU(:,:,mmin) = 2.d0 * U(:,:,mmin)
    else
       dU(:,:,mmin) = 0.d0
    end if
    do k = mmin,mmax-1
       where (dU(:,:,k)*dU(:,:,k+1) .gt. 0.d0)
         nabla_U(:,:,k) = 2.d0*dU(:,:,k)/(dU(:,:,k)+dU(:,:,k+1))*dU(:,:,k+1)
       elsewhere
         nabla_U(:,:,k) = 0.d0
       end where
    end do
    nabla_U(:,:,mmax) = dU(:,:,mmax)
    if (Z(1,1,mmin) .lt. 0.d0 .or. Symmetry==NO_SYMM)  nabla_U(:,:,mmin) = dU(:,:,mmin+1)
 return
 end if
end subroutine compute_nabla_van_leer
!--------------------------------------------------!
! MC limiter                                           !
!--------------------------------------------------!
!
subroutine compute_nabla_mc(ext,dX,dY,dZ,U,nabla_U,m,Sym_arg,X,Y,Z,Symmetry)
 implicit none
  interface
    SUBROUTINE polint(xa,ya,n,x,y,dy)
      IMPLICIT NONE
      integer                    :: n
      real*8                     :: dy,x,y
      real*8, dimension(n)       :: xa,ya
    end SUBROUTINE polint
  end interface
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: U, nabla_U
 integer                                        :: m,i,j,k,mmax,mmin
 real*8, dimension(ext(1),ext(2),ext(3))        :: mdU,dU
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8                                         :: Sym_arg
 integer                                        :: ii, l = 3
 real*8, dimension(3)                           :: xa, ya
 real*8                                         :: dnabla, dX, dY, dZ
 integer                                        :: imin, imax, jmin, jmax, kmin, kmax
 integer                            :: Symmetry, NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
 parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
!
 mmax = ubound(U,m)
 mmin = lbound(U,m)
 imax = ubound(U,1)
 imin = lbound(U,1)
 jmax = ubound(U,2)
 jmin = lbound(U,2)
 kmax = ubound(U,3)
 kmin = lbound(U,3)
 if(Symmetry == OCTANT) then
    if(m==1 .and. X(1,1,1).lt.0.D0) then
       mmin = lbound(U,1)+1
    else if(m==2 .and. Y(1,1,1).lt.0.D0) then
       mmin = lbound(U,2)+1
    else if(m==3 .and. Z(1,1,1).lt.0.D0) then
       mmin = lbound(U,3)+1
    end if
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
    if(m==1 .and. X(1,1,1).lt.0.D0) then
       mmin = lbound(U,1)+1
    else if(m==3 .and. Z(1,1,1).lt.0.D0) then
       mmin = lbound(U,3)+1
    end if
    if(X(1,1,1).lt.0.D0) then
       imin = imin + 1
    end if
    if(Z(1,1,1).lt.0.D0) then
       kmin = kmin + 1
    end if
 else if(Symmetry == EQUATORIAL) then
    if(m==3 .and. Z(1,1,1).lt.0.D0) then
       mmin = lbound(U,3)+1
    end if
    if(Z(1,1,1).lt.0.D0) then
       kmin = kmin + 1
    end if
 end if
 if (m==1) then
    do i = mmin+1,mmax
       dU(i,:,:) = U(i,:,:) - U(i-1,:,:)
    end do
    if (Sym_arg < 0.d0) then
       dU(mmin,:,:) = 2.d0 * U(imin,:,:)
    else
       dU(mmin,:,:) = 0.d0
    end if
    do k = kmin, kmax
       do j=jmin,jmax
          do i = mmin,mmax-1
             if (dU(i,j,k)*dU(i+1,j,k) .gt. 0.d0) then
                mdU(i,j,k) = 0.5d0 * ( dU(i,j,k) + dU(i+1,j,k) )
!                if(i==50 .and. j==2 .and. k==30) then 
!                   write(*,*) HIII2a,mdU(50,2,30),0.5d0 * ( dU(50,2,30) + dU(51,2,30) )
!                end if
                if ( (abs(dU(i,j,k)) .lt. abs(dU(i+1,j,k))) .and. &
                     (2.d0*abs(dU(i,j,k)) .lt. abs(mdU(i,j,k))) ) then
                   nabla_U(i,j,k) = 2.d0*dU(i,j,k)
                elseif ( (abs(dU(i+1,j,k)) .lt. abs(dU(i,j,k))) .and. &
                     (2.d0*abs(dU(i+1,j,k)) .lt. abs(mdU(i,j,k))) ) then
                   nabla_U(i,j,k) = 2.d0*dU(i+1,j,k)
                else
                   nabla_U(i,j,k) = mdU(i,j,k)
                end if
             else
                nabla_U(i,j,k) = 0.d0
             end if
          end do
          nabla_U(mmax,j,k) = dU(mmax,j,k)
       end do
    end do
!    write(*,*) HIII1,m,dU(50,2,30)
!    write(*,*) HIII2b,mdU(50,2,30),0.5d0 * ( dU(50,2,30) + dU(51,2,30) )
!    write(*,*) HIII2,m,nabla_U(50,2,30),2.d0*dU(50,2,30),2.d0*dU(51,2,30),mdU(50,2,30)!,0.5d0 * ( dU(50,2,30) + dU(51,2,30) )
    return
 end if
 if (m==2) then
    do j = mmin+1,mmax
       dU(:,j,:) = U(:,j,:) - U(:,j-1,:)
    end do
    if (Sym_arg < 0) then
       dU(:,mmin,:) = 2.d0*U(:,mmin,:)
    else
       dU(:,mmin,:) = 0.d0
    end if
    do k = kmin,kmax
       do i = imin,imax
          do j = mmin,mmax-1
             if (dU(i,j,k)*dU(i,j+1,k) .gt. 0.d0) then
                mdU(i,j,k) = 0.5d0 * ( dU(i,j,k) + dU(i,j+1,k) )
                if ( (abs(dU(i,j,k)) .lt. abs(dU(i,j+1,k))) .and. &
                     (2.d0*abs(dU(i,j,k)) .lt. abs(mdU(i,j,k))) ) then
                   nabla_U(i,j,k) = 2.d0*dU(i,j,k)
                else if ( (abs(dU(i,j+1,k)) .lt. abs(dU(i,j,k))) .and. &
                     (2.d0*abs(dU(i,j+1,k)) .lt. abs(mdU(i,j,k))) ) then
                   nabla_U(i,j,k) = 2.d0*dU(i,j+1,k)
                else
                   nabla_U(i,j,k) = mdU(i,j,k)
                end if
             else
                nabla_U(i,j,k) = 0.d0
             end if
          end do
          nabla_U(i,mmax,k) = dU(i,mmax,k)
       end do
    end do
    return
 end if
 if (m==3) then
    do k = mmin+1,mmax
       dU(:,:,k) = U(:,:,k) - U(:,:,k-1)
    end do
    if (Sym_arg < 0) then
       dU(:,:,mmin) = 2.d0*U(:,:,mmin)
    else
       dU(:,:,mmin) = 0.d0
    end if
    do j=jmin,jmax
       do i=imin, imax
          do k = mmin,mmax-1
             if (dU(i,j,k)*dU(i,j,k+1) .gt. 0.d0) then
                mdU(i,j,k) = 0.5d0 * ( dU(i,j,k) + dU(i,j,k+1) )
                if ( (abs(dU(i,j,k)) .lt. abs(dU(i,j,k+1))) .and. &
                     (2.d0*abs(dU(i,j,k)) .lt. abs(mdU(i,j,k))) ) then
                   nabla_U(i,j,k) = 2.d0*dU(i,j,k)
                elseif ( (abs(dU(i,j,k+1)) .lt. abs(dU(i,j,k))) .and. &
                     (2.d0*abs(dU(i,j,k+1)) .lt. abs(mdU(i,j,k))) ) then
                   nabla_U(i,j,k) = 2.d0*dU(i,j,k+1)
                else
                   nabla_U(i,j,k) = mdU(i,j,k)
                end if
             else
                nabla_U(i,j,k) = 0.d0
             end if
          end do
       end do
       nabla_U(:,:,mmax) = dU(:,:,mmax)
       if (Z(1,1,mmin) .lt. 0.d0 .or. Symmetry==NO_SYMM) nabla_U(:,:,mmin) = dU(:,:,mmin+1)
    end do
    return
 end if
end subroutine compute_nabla_mc
subroutine compute_nabla_du(ext,U,nabla_U,m,Sym_arg,X,Y,Z,Symmetry)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: U, nabla_U
 integer                                        :: m,i,j,k,mmax,mmin
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8                                         :: Sym_arg
 real*8                                         :: dnabla, dX, dY, dZ
 integer                                        :: imin,imax,jmin,jmax,kmin,kmax
 real*8, parameter                                :: b = 2.d0
 integer                            :: Symmetry, NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
 parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
!
 mmax = ubound(U,m)
 mmin = lbound(U,m)
 if(Symmetry == OCTANT) then
    if(m==1 .and. X(1,1,1).lt.0.D0) then
       mmin = lbound(U,1)+1
    else if(m==2 .and. Y(1,1,1).lt.0.D0) then
       mmin = lbound(U,2)+1
    else if(m==3 .and. Z(1,1,1).lt.0.D0) then
       mmin = lbound(U,3)+1
    end if
 else if(Symmetry == AXISYM) then
    if(m==1 .and. X(1,1,1).lt.0.D0) then
       mmin = lbound(U,1)+1
    else if(m==3 .and. Z(1,1,1).lt.0.D0) then
       mmin = lbound(U,3)+1
    end if
 else if(Symmetry == EQUATORIAL) then
    if(m==3 .and. Z(1,1,1).lt.0.D0) then
       mmin = lbound(U,3)+1
    end if
 end if
 if (m==1) then
    nabla_U(mmin,:,:) = (1.d0-Sym_arg)*U(mmin,:,:)
    do i = mmin+1,mmax
       nabla_U(i,:,:) = U(i,:,:) - U(i-1,:,:)
    end do
 elseif (m==2) then
    nabla_U(:,mmin,:) = (1.d0-Sym_arg)*U(:,mmin,:)
    do j = mmin+1,mmax
       nabla_U(:,j,:) = U(:,j,:) - U(:,j-1,:)
    end do
 else
    do k = mmin+1,mmax
       nabla_U(:,:,k) = U(:,:,k) - U(:,:,k-1)
    end do
    if (Z(1,1,mmin) .lt. 0.d0) then
       nabla_U(:,:,mmin) = nabla_U(:,:,mmin+1)
    else
       nabla_U(:,:,mmin) = (1.d0-Sym_arg)*U(:,:,mmin)
    end if
 end if
end subroutine compute_nabla_du
subroutine compute_u0(ext,u0,rho_star,w,lapse)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: u0,rho_star,w,lapse
! 
 where (rho_star .gt. 0.d0)
   u0 = w/rho_star/(1.d0+lapse)
 elsewhere
   u0 = 1.d0
 end where
end subroutine compute_u0
!-----------------------------------------------------------------!
! Compute 1 + u_0                                                  !
!-----------------------------------------------------------------!
!
subroutine compute_one_plus_u_0(ext,one_plus_u_0,one_plus_h_u_0, &
                                P,h,rho_b,n,u0,vx,vy,vz, &
                                lapse,shiftx,shifty,shiftz, &
                                phi,gxx,gxy,gxz,gyy,gyz,gzz)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: u0,vx,vy,vz,P,rho_b,h
 real*8, dimension(ext(1),ext(2),ext(3))        :: lapse,shiftx,shifty,shiftz
 real*8, dimension(ext(1),ext(2),ext(3))        :: phi,gxx,gxy,gxz,gyy,gyz,gzz
 real*8, dimension(ext(1),ext(2),ext(3))        :: one_plus_u_0,one_plus_h_u_0
 real*8                                                :: n
!
 real*8, dimension(ext(1),ext(2),ext(3))        :: u0_1
!
 integer                                        :: i, imin,imax
!
! Compute gamma_{ij} (v^i+beta^i)(v^j+beta^j) and store in u0_1
 where (rho_b .gt. 0.d0)
   u0_1 = exp(6.d0*phi)*( gxx*(vx+shiftx)**2 + &
                2.d0*gxy*(vx+shiftx)*(vy+shifty) + &
                2.d0*gxz*(vx+shiftx)*(vz+shiftz) + &
                gyy*(vy+shifty)**2 + &
                2.d0*gyz*(vy+shifty)*(vz+shiftz) + &
                gzz*(vz+shiftz)**2 )
! Compute sqrt( alpha^2 - gamma_{ij} (v^i+beta^1)(v^j+beta^j) ) and store in 
! one_plus_u_0
   one_plus_u_0 = sqrt( (1.d0+lapse)**2 - u0_1 )
! Compute u^0-1 and store in u0_1
   u0_1 = (u0_1 - lapse*(2.d0+lapse) )/one_plus_u_0/(1.d0+one_plus_u_0)
! Compute u^i u_i and store in one_plus_u_0
   one_plus_u_0 = u0**2 * exp(6.d0*phi) * ( gxx*vx*(shiftx+vx) + &
        gxy*(vx*(shifty+vy)+vy*(shiftx+vx)) + &
        gxz*(vx*(shiftz+vz)+vz*(shiftx+vx)) + &
        gyy*vy*(shifty+vy) + gyz*(vy*(shiftz+vz)+vz*(shifty+vy)) + &
        gzz*vz*(shiftz+vz) )
! Finally, compute 1 + u_0 and 1 + h u_0
   one_plus_u_0 = (u0_1 - one_plus_u_0)/u0
   one_plus_h_u_0 = one_plus_u_0 + (1.d0+n)*P/rho_b*(one_plus_u_0-1.d0)
 elsewhere
   one_plus_u_0 = 0.d0
   one_plus_h_u_0 = 0.d0
 end where
  one_plus_u_0 = 1.d0-u0
  one_plus_h_u_0 = 1.d0-h*u0
end subroutine compute_one_plus_u_0
!--------------------------------------------------!
! cmin and cmax                                           !
!--------------------------------------------------!
!
subroutine find_cmin_cmax(ext,cplus,cminus,cmin,cmax,m)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: cplus,cminus,cmin,cmax
 integer                                         :: m,i,j,k,mmin,mmax
!
 mmin = lbound(cplus,m)
 mmax = ubound(cplus,m)
 if (m==1) then
    cmax(mmin,:,:) = max( 0.d0, cplus(mmin,:,:))
    cmin(mmin,:,:) = -min( 0.d0, cminus(mmin,:,:))
    do i = mmin+1,mmax
       cmax(i,:,:) = max( 0.d0, cplus(i-1,:,:), cplus(i,:,:))
       cmin(i,:,:) = -min( 0.d0, cminus(i-1,:,:), cminus(i,:,:))
    end do
    return
 end if
 if (m==2) then
    cmax(:,mmin,:) = max( 0.d0, cplus(:,mmin,:))
    cmin(:,mmin,:) = -min( 0.d0, cminus(:,mmin,:))
    do j = mmin+1,mmax
       cmax(:,j,:) = max( 0.d0, cplus(:,j-1,:), cplus(:,j,:))
       cmin(:,j,:) = -min( 0.d0, cminus(:,j-1,:), cminus(:,j,:))
    end do
    return
 end if
 if (m==3) then
    cmax(:,:,mmin) = max( 0.d0, cplus(:,:,mmin))
    cmin(:,:,mmin) = -min( 0.d0, cminus(:,:,mmin))
    do k = mmin+1,mmax
       cmax(:,:,k) = max( 0.d0, cplus(:,:,k-1), cplus(:,:,k))
       cmin(:,:,k) = -min( 0.d0, cminus(:,:,k-1), cminus(:,:,k))
    end do
    return
 end if
!!$ cmin = 1.d0
!!$ cmax = 1.d0
end subroutine find_cmin_cmax
!
!
subroutine compute_cmax_cmin_hybrid(ext,v02r,v02l,cmax,cmin,rho_br,rho_bl, &
                Pr,Pl, vxr,vxl,vyr,vyl,vzr,vzl, Bxr,Bxl,Byr,Byl,Bzr,Bzl, &
                alpha,betax,betay,betaz,phi,gxx,gxy,gxz,gyy,gyz,gzz, &
                gupxx,gupyy,gupzz,m,neos,ergo_star,ergo_sigma,rho_tab, P_tab, eps_tab, k_tab, &
                gamma_tab, gamma_th, &
                temp1,temp2,temp3,temp4,temp5, &
                cplusr,cminusr, u0,cplusl,cminusl)
!
 implicit none
 interface
   subroutine find_v02_hybrid(ext,v02,rho_b,P,u0,vx,vy,vz, &
        neos,ergo_star,ergo_sigma,rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, &
        sqBx,sqBy,sqBz, &
        lapse,shiftx,shifty,shiftz, phi,gxx,gxy,gxz,gyy,gyz,gzz, &
        temp1,temp2,temp3,temp4,temp5)
     implicit none
     integer, dimension(3)                          :: ext
     real*8, dimension(ext(1),ext(2),ext(3))        :: vx,vy,vz,u0
     real*8, dimension(ext(1),ext(2),ext(3))        :: sqBx,sqBy,sqBz
     real*8, dimension(ext(1),ext(2),ext(3))        :: rho_b,P,v02
     real*8, dimension(ext(1),ext(2),ext(3))        :: lapse,phi
     real*8, dimension(ext(1),ext(2),ext(3))        :: shiftx,shifty,shiftz
     real*8, dimension(ext(1),ext(2),ext(3))        :: gxx,gxy,gxz,gyy,gyz,gzz
     real*8, dimension(ext(1),ext(2),ext(3))        :: temp1,temp2,temp3,temp4,temp5
     integer                                        :: neos,ergo_star
     real*8, dimension(neos)                        :: rho_tab, P_tab, eps_tab
     real*8, dimension(neos+1)                      :: k_tab, gamma_tab
     real*8                                         :: gamma_th,ergo_sigma
   end subroutine find_v02_hybrid
 end interface
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: v02r,v02l,cmax,cmin
 real*8, dimension(ext(1),ext(2),ext(3))        :: rho_br,rho_bl,Pr,Pl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vxr,vxl,vyr,vyl,vzr,vzl
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bxr,Bxl,Byr,Byl,Bzr,Bzl
 real*8, dimension(ext(1),ext(2),ext(3))        :: alpha,betax,betay,betaz
 real*8, dimension(ext(1),ext(2),ext(3))        :: phi,gxx,gxy,gxz,gyy,gyz,gzz
 real*8, dimension(ext(1),ext(2),ext(3))        :: temp1,temp2,temp3,temp4,temp5
 real*8, dimension(ext(1),ext(2),ext(3))        :: gupxx,gupyy,gupzz
 real*8, dimension(ext(1),ext(2),ext(3))        :: cplusr,cminusr, u0
 real*8, dimension(ext(1),ext(2),ext(3))        :: cplusl,cminusl
 integer                                        :: neos,ergo_star
 real*8, dimension(neos)                        :: rho_tab, P_tab, eps_tab
 real*8, dimension(neos+1)                      :: k_tab, gamma_tab
 real*8                                         :: gamma_th,ergo_sigma
 integer                                        :: m
!
 call find_v02_hybrid(ext,v02r,rho_br,Pr,u0,vxr,vyr,vzr, &
        neos,ergo_star,ergo_sigma,rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, &
        Bxr,Byr,Bzr, alpha,betax,betay,betaz, phi,gxx,gxy,gxz,gyy,gyz,gzz, &
        temp1,temp2,temp3,temp4,temp5)
 if (m==1) then
    call find_cplus_cminus(ext,cplusr,cminusr,v02r,u0, &
                vxr,alpha,betax,phi,gupxx)
 else if (m==2) then
    call find_cplus_cminus(ext,cplusr,cminusr,v02r,u0, &
                vyr,alpha,betay,phi,gupyy)
 else
    call find_cplus_cminus(ext,cplusr,cminusr,v02r,u0, &
                vzr,alpha,betaz,phi,gupzz)
 end if
 call find_v02_hybrid(ext,v02l,rho_bl,Pl,u0,vxl,vyl,vzl, &
        neos,ergo_star,ergo_sigma,rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, &
        Bxl,Byl,Bzl, alpha,betax,betay,betaz, phi,gxx,gxy,gxz,gyy,gyz,gzz, &
        temp1,temp2,temp3,temp4,temp5)
 if (m==1) then
    call find_cplus_cminus(ext,cplusl,cminusl,v02l,u0, &
                vxl,alpha,betax,phi,gupxx)
 else if (m==2) then
    call find_cplus_cminus(ext,cplusl,cminusl,v02l,u0, &
                vyl,alpha,betay,phi,gupyy)
 else
    call find_cplus_cminus(ext,cplusl,cminusl,v02l,u0, &
                vzl,alpha,betaz,phi,gupzz)
 end if
 cmax = max(0.d0,cplusl,cplusr)
 cmin = -min(0.d0,cminusl,cminusr)
end subroutine compute_cmax_cmin_hybrid
!--------------------------------------------------!
! HLL flux                                            !
!--------------------------------------------------!
!
subroutine flux_hll(ext, Ur,Ul, Fr,Fl, F, cmax,cmin)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: Ur,Ul,F,cmin,cmax,Fl,Fr
!
 F = (cmin*Fr + cmax*Fl - cmin*cmax*(Ur-Ul) )/(cmax + cmin)
end subroutine flux_hll
!
!
subroutine compute_cmax_cmin(ext,v02r,v02l,cmax,cmin,rho_br,rho_bl, &
                Pr,Pl, vxr,vxl,vyr,vyl,vzr,vzl, Bxr,Bxl,Byr,Byl,Bzr,Bzl, &
                alpha,betax,betay,betaz,phi,gxx,gxy,gxz,gyy,gyz,gzz, &
                gupxx,gupyy,gupzz,n,m)
!
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: v02r,v02l,cmax,cmin
 real*8, dimension(ext(1),ext(2),ext(3))        :: rho_br,rho_bl,Pr,Pl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vxr,vxl,vyr,vyl,vzr,vzl
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bxr,Bxl,Byr,Byl,Bzr,Bzl
 real*8, dimension(ext(1),ext(2),ext(3))        :: alpha,betax,betay,betaz
 real*8, dimension(ext(1),ext(2),ext(3))        :: phi,gxx,gxy,gxz,gyy,gyz,gzz
 real*8, dimension(ext(1),ext(2),ext(3))        :: gupxx,gupyy,gupzz
 real*8, dimension(ext(1),ext(2),ext(3))        :: cplusr,cminusr, u0
 real*8, dimension(ext(1),ext(2),ext(3))        :: cplusl,cminusl
 real*8                                                :: n
 integer                                        :: m
!
 call find_v02(ext,v02r,rho_br,Pr,n,u0,vxr,vyr,vzr, Bxr,Byr,Bzr, &
        alpha,betax,betay,betaz, phi,gxx,gxy,gxz,gyy,gyz,gzz)
 if (m==1) then
    call find_cplus_cminus(ext,cplusr,cminusr,v02r,u0, &
                vxr,alpha,betax,phi,gupxx)
 else if (m==2) then
    call find_cplus_cminus(ext,cplusr,cminusr,v02r,u0, &
                vyr,alpha,betay,phi,gupyy)
 else
    call find_cplus_cminus(ext,cplusr,cminusr,v02r,u0, &
                vzr,alpha,betaz,phi,gupzz)
 end if
 call find_v02(ext,v02l,rho_bl,Pl,n,u0,vxl,vyl,vzl, Bxl,Byl,Bzl, &
        alpha,betax,betay,betaz, phi,gxx,gxy,gxz,gyy,gyz,gzz)
 if (m==1) then
    call find_cplus_cminus(ext,cplusl,cminusl,v02l,u0, &
                vxl,alpha,betax,phi,gupxx)
 else if (m==2) then
    call find_cplus_cminus(ext,cplusl,cminusl,v02l,u0, &
                vyl,alpha,betay,phi,gupyy)
 else
    call find_cplus_cminus(ext,cplusl,cminusl,v02l,u0, &
                vzl,alpha,betaz,phi,gupzz)
 end if
 cmax = max(0.d0,cplusl,cplusr)
 cmin = -min(0.d0,cminusl,cminusr)
end subroutine compute_cmax_cmin
!----------------------------------------------------!
! Find maximum and minimum velocities along x^i-axis !
!----------------------------------------------------!
!
subroutine find_cplus_cminus(ext,cplus,cminus,v02,u0, &
                 vi,lapse,shifti,phi,gupii)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: v02,cplus,cminus,u0,vi
 real*8, dimension(ext(1),ext(2),ext(3))        :: lapse,shifti,gupii,phi
 real*8                                         :: cp
 integer                                        :: imin,imax,jmin,jmax,kmin,kmax
 integer                                        :: i,j,k
!
! Auxiliary arrays
!
 real*8, dimension(ext(1),ext(2),ext(3))        :: a,b,c,detm
!
 imin = lbound(u0,1)
 jmin = lbound(u0,2)
 kmin = lbound(u0,3)
 imax = ubound(u0,1)
 jmax = ubound(u0,2)
 kmax = ubound(u0,3)
 a = u0**2 * (1.d0-v02) + v02/(1.d0+lapse)**2
 b = 2.d0* ( shifti/(1.d0+lapse)**2 * v02 - u0**2 * vi * (1.d0-v02) )
 c = (u0*vi)**2 * (1.d0-v02) - v02 * ( exp(-4.d0*phi)*gupii - &
        (shifti/(1.d0+lapse))**2 )
! write(*,*) cplusli3,shifti(10,2,10),lapse(10,2,10)
 detm = sqrt( max(b*b - 4.d0*a*c,0.d0) )
 cplus = (detm-b)/2.d0/a
 cminus = -(detm+b)/2.d0/a
!!$ write(*,*) Inside find_cplus_cminus:
!!$ write(*,*) cplus(1,1,1)  = , cplus(1,1,1)
!!$ write(*,*) cminus(1,1,1) = , cminus(1,1,1)
!!$ write(*,*) detm(1,1,1)   = , detm(1,1,1)
!!$ write(*,*) b(1,1,1)      = , b(1,1,1)
!!$ write(*,*) a(1,1,1)      = , a(1,1,1)
!!$ write(*,*) c(1,1,1)      = , c(1,1,1)
!!$ write(*,*) phi(1,1,1)    = , phi(1,1,1)
!!$ write(*,*) gupii(1,1,1)  = , gupii(1,1,1)
 do k=kmin,kmax
    do j=jmin,jmax
       do i=imin,imax
          if (cplus(i,j,k) .lt. cminus(i,j,k)) then
             cp = cminus(i,j,k)
             cminus(i,j,k) = cplus(i,j,k)
             cplus(i,j,k) = cp
          end if
       end do
    end do
 end do
end subroutine find_cplus_cminus
