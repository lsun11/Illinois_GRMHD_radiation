!-----------------------------------------------------------------!
! compute \Delta_0f and \Delta_0^2f  (f = {rho_b,P,v^i})
!-----------------------------------------------------------------!
!
subroutine find_centderivs(ext,X,Y,Z,rho,d1rho,d2rho,P,d1P,d2P, &
     vx,d1vx,d2vx,vy,d1vy,d2vy,vz,d1vz,d2vz,m,Symmetry)
 implicit none
 interface
    subroutine ddx12(ext,X,f,d1f,d2f,dX,sym,Symmetry)
      implicit none
      integer, dimension(3)                    :: ext
      real*8, dimension(ext(1),ext(2),ext(3))  :: X
      real*8, dimension(ext(1),ext(2),ext(3))  :: f,d1f,d2f
      real*8                                   :: sym,dX
      integer                                  :: Symmetry
    end subroutine ddx12
    subroutine ddy12(ext,Y,f,d1f,d2f,dY,sym,Symmetry)
      implicit none
      integer, dimension(3)                    :: ext
      real*8, dimension(ext(1),ext(2),ext(3))  :: Y
      real*8, dimension(ext(1),ext(2),ext(3))  :: f,d1f,d2f
      real*8                                   :: sym,dY
      integer                                  :: Symmetry
    end subroutine ddy12
    subroutine ddz12(ext,Z,f,d1f,d2f,dZ,sym,Symmetry)
      implicit none
      integer, dimension(3)                    :: ext
      real*8, dimension(ext(1),ext(2),ext(3))  :: Z
      real*8, dimension(ext(1),ext(2),ext(3))  :: f,d1f,d2f
      real*8                                   :: sym,dZ
      integer                                  :: Symmetry
    end subroutine ddz12
 end interface
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8, dimension(ext(1),ext(2),ext(3))        :: rho,d1rho,d2rho
 real*8, dimension(ext(1),ext(2),ext(3))        :: P,d1P,d2P
 real*8, dimension(ext(1),ext(2),ext(3))        :: vx,d1vx,d2vx
 real*8, dimension(ext(1),ext(2),ext(3))        :: vy,d1vy,d2vy
 real*8, dimension(ext(1),ext(2),ext(3))        :: vz,d1vz,d2vz
 integer, intent(in)                            :: m,Symmetry
 integer                                        :: i,j,k
 real*8                     :: ONE, TWO, SYM, ANTI
 parameter(ONE = 1.0D0, TWO = 2.D0, SYM = 1.D0, ANTI = -1.D0)
!
! X derivatives
!
 if(m==1) then
    call ddx12(ext,X,rho,d1rho,d2rho,ONE,SYM,Symmetry)
    call ddx12(ext,X,P,  d1P,  d2P,  ONE,SYM,Symmetry)
    call ddx12(ext,X,vz, d1vz, d2vz, ONE,SYM,Symmetry)
    if(Symmetry==2) then
       call ddx12(ext,X,vx,d1vx,d2vx,ONE,ANTI,Symmetry)
       call ddx12(ext,X,vy,d1vy,d2vy,ONE,SYM,Symmetry)
    else if(Symmetry==4) then
       call ddx12(ext,X,vx,d1vx,d2vx,ONE,ANTI,Symmetry)
       call ddx12(ext,X,vy,d1vy,d2vy,ONE,ANTI,Symmetry)
    else
       call ddx12(ext,X,vx,d1vx,d2vx,ONE,SYM,Symmetry)
       call ddx12(ext,X,vy,d1vy,d2vy,ONE,SYM,Symmetry)
    end if

!*********Testing***********
!!$    do i=1,32
!!$       write(*,*) d2rho(i,2,3)
!!$    end do

 else if(m==2) then
    if(Symmetry.ne.4) then
       call ddy12(ext,Y,rho,d1rho,d2rho,ONE,SYM,Symmetry)
       call ddy12(ext,Y,P,  d1P,  d2P,  ONE,SYM,Symmetry)
       call ddy12(ext,Y,vz, d1vz, d2vz, ONE,SYM,Symmetry)
       call ddy12(ext,Y,vx, d1vx, d2vx, ONE,SYM,Symmetry)
       if(Symmetry==2) then
          call ddy12(ext,Y,vy,d1vy,d2vy,ONE,ANTI,Symmetry)
       else
          call ddy12(ext,Y,vy,d1vy,d2vy,ONE,SYM,Symmetry)
       end if
    end if
 else if(m==3) then
    call ddz12(ext,Z,rho,d1rho,d2rho,ONE,SYM,Symmetry)
    call ddz12(ext,Z,P,  d1P,  d2P,  ONE,SYM,Symmetry)
    call ddz12(ext,Z,vx, d1vx, d2vx, ONE,SYM,Symmetry)
    call ddz12(ext,Z,vy, d1vy, d2vy, ONE,SYM,Symmetry)
    if(Symmetry.ne.0) then
       call ddz12(ext,Z,vz,d1vz,d2vz,ONE,ANTI,Symmetry)
    else
       call ddz12(ext,Z,vz,d1vz,d2vz,ONE,SYM,Symmetry)
    end if
 else
    write(*,*) "m not in proper range"
 end if

! write(*,*) "inside centderivs:",m,d1vx(15,2,2),d1vy(15,2,2),d1vz(15,2,2),d2vx(15,2,2),d2vy(15,2,2),d2vz(15,2,2),vx(2,2,3),vx(15,2,2),vx(2,2,4),vx(15,2,2)+vx(2,2,3),vx(15,2,2)-vx(2,2,3)

!
end subroutine find_centderivs

!-----------------------------------------------------------------------------
!
! Computes second partial derivative of function f with respect to x
!
!-----------------------------------------------------------------------------
subroutine ddx12(ext,X,f,d1f,d2f,dX,sym,Symmetry)
!
! Input parameters:
!
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: X
  real*8, dimension(ext(1),ext(2),ext(3))  :: f,d1f,d2f
  real*8                                   :: sym,dX
! sym  :    symmetry imposed at boundaries ( 1: symmetric, -1: antisymmetric)
  integer                    :: i,Symmetry
  integer                    :: imin, imax
  real*8                     :: dxm1,dxm2
  real*8                     :: HALF, ONE, TWO
  integer                            :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(HALF = 0.5D0, ONE = 1.0D0, TWO = 2.D0)
  imin = lbound(f,1)
  imax = ubound(f,1)
  if(Symmetry == OCTANT .or. Symmetry == AXISYM) then
     if(X(1,1,1) .lt. 0.D0) then
        imin = imin + 1
     end if
  end if
! Derivative with respect to x
  dxm1 = HALF/dX
  dxm2 = ONE/(dX*dX)
  d1f(imin,:,:) = dxm1 * (f(imin+1,:,:) - sym*f(imin,:,:))
  d2f(imin,:,:) = dxm2 * (f(imin+1,:,:) + (sym - TWO)*f(imin,:,:))
  do i = imin + 1, imax - 1
     d1f(i,:,:) = dxm1 * (f(i+1,:,:) - f(i-1,:,:))
     d2f(i,:,:) = dxm2 * (f(i+1,:,:) + f(i-1,:,:) - TWO*f(i,:,:))
  enddo
  d1f(imax,:,:) = TWO * dxm1 * (f(imax,:,:) - f(imax-1,:,:))
  d2f(imax,:,:) = d2f(imax-1,:,:)
end subroutine ddx12

!-----------------------------------------------------------------------------
!
! Computes partial derivative of function f with respect to y
!
!-----------------------------------------------------------------------------
subroutine ddy12(ext,Y,f,d1f,d2f,dY,sym,Symmetry)
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: Y
  real*8, dimension(ext(1),ext(2),ext(3))  :: f,d1f,d2f
  real*8                                   :: sym,dY
! sym  :    symmetry imposed at boundaries ( 1: symmetric, -1: antisymmetric)
  integer                    :: j,Symmetry
  integer                    :: jmin, jmax
  real*8                     :: dym1,dym2
  real*8                     :: HALF, ONE, TWO
  integer                            :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(HALF = 0.5D0, ONE = 1.0D0, TWO = 2.D0)
  jmin = lbound(f,2)
  jmax = ubound(f,2)
  if(Symmetry == OCTANT) then
     if(Y(1,1,1) .lt. 0.D0) then
        jmin = jmin + 1
     end if
  end if
! Derivative with respect to y
  dym1 = HALF/dY
  dym2 = ONE/(dY*dY)
  d1f(:,jmin,:) = dym1 * (f(:,jmin+1,:) - sym*f(:,jmin,:))
  d2f(:,jmin,:) = dym2 * (f(:,jmin+1,:) + (sym - TWO)*f(:,jmin,:))
  do j = jmin + 1, jmax - 1
     d1f(:,j,:) = dym1 * (f(:,j+1,:) - f(:,j-1,:))
     d2f(:,j,:) = dym2 * (f(:,j+1,:) + f(:,j-1,:) - TWO*f(:,j,:))
  enddo
  d1f(:,jmax,:) = TWO * dym1 * (f(:,jmax,:) - f(:,jmax-1,:))
  d2f(:,jmax,:) = d2f(:,jmax-1,:)
end subroutine ddy12

!-----------------------------------------------------------------------------
!
! Computes partial derivative of function f with respect to z
!
!-----------------------------------------------------------------------------
subroutine ddz12(ext,Z,f,d1f,d2f,dZ,sym,Symmetry)
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: Z
  real*8, dimension(ext(1),ext(2),ext(3))  :: f,d1f,d2f
  real*8                                   :: sym,dZ
! sym  :    symmetry imposed at boundaries ( 1: symmetric, -1: antisymmetric)
  integer                    :: k,Symmetry
  integer                    :: kmin, kmax  
  real*8                     :: dzm1, dzm2
  real*8                     :: HALF, ONE, TWO
  integer                            :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(HALF = 0.5D0, ONE = 1.0D0, TWO = 2.D0)
  kmin = lbound(f,3)
  kmax = ubound(f,3)
  if(Symmetry == OCTANT .or. Symmetry == AXISYM .or. Symmetry == EQUATORIAL) then
     if(Z(1,1,1) .lt. 0.D0) then
        kmin = kmin + 1
     end if
  end if

!
! Derivative with respect to z
!
  dzm1 = HALF/dZ
  dzm2 = ONE/(dZ*dZ)
  d1f(:,:,kmin) = dzm1 * (f(:,:,kmin+1) - sym*f(:,:,kmin))
  d2f(:,:,kmin) = dzm2 * (f(:,:,kmin+1) + (sym - TWO)*f(:,:,kmin))
  do k = kmin + 1, kmax - 1
     d1f(:,:,k) = dzm1 * (f(:,:,k+1) - f(:,:,k-1))
     d2f(:,:,k) = dzm2 * (f(:,:,k+1) + f(:,:,k-1) - TWO*f(:,:,k))
  enddo
  d1f(:,:,kmax) = TWO * dzm1 * (f(:,:,kmax) - f(:,:,kmax-1))
  d2f(:,:,kmax) = d2f(:,:,kmax-1)
end subroutine ddz12


!-----------------------------------------------------------------!
! compute \tilde Q_i(x_{i\pm 1/2})
!-----------------------------------------------------------------!
!
subroutine find_ur_ul_ceno(ext,rho,d1rho,d2rho,rhor,rhol,P,d1P,d2P,Pr,Pl, &
     vx,d1vx,d2vx,vxr,vxl,vy,d1vy,d2vy,vyr,vyl,vz,d1vz,d2vz,vzr,vzl,m,Symmetry, &
     alpha_f, betax_f, betay_f, betaz_f, phi_f, &
     gxx_f, gxy_f, gxz_f, gyy_f, gyz_f, gzz_f)
 implicit none
 interface
    subroutine ceno_ur_ul(ext,f,d1f,d2f,fr,fl,m,sym)
      integer, dimension(3)                     :: ext
      real*8, dimension(ext(1),ext(2),ext(3))   :: f,d1f,d2f,fr,fl
      integer, intent(in)                       :: m
      real*8                                    :: sym
    end subroutine ceno_ur_ul
 end interface
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: rho,d1rho,d2rho,rhor,rhol
 real*8, dimension(ext(1),ext(2),ext(3))        :: P,d1P,d2P,Pr,Pl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vx,d1vx,d2vx,vxr,vxl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vy,d1vy,d2vy,vyr,vyl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vz,d1vz,d2vz,vzr,vzl
 real*8, dimension(ext(1),ext(2),ext(3))        :: alpha_f, phi_f
 real*8, dimension(ext(1),ext(2),ext(3))        :: betax_f, betay_f, betaz_f
 real*8, dimension(ext(1),ext(2),ext(3))        :: gxx_f, gxy_f, gxz_f
 real*8, dimension(ext(1),ext(2),ext(3))        :: gyy_f, gyz_f, gzz_f
 integer                                        :: m,Symmetry
 real*8                                         :: ONE, TWO, SYM, ANTI
 real*8, parameter                              :: fac = 0.99d0
 real*8, dimension(ext(1),ext(2),ext(3))        :: temp
 parameter(ONE = 1.0D0, TWO = 2.D0, SYM = 1.D0, ANTI = -1.D0)
  
 call ceno_ur_ul(ext,rho,d1rho,d2rho,rhor,rhol,m,SYM)
 call ceno_ur_ul(ext,P,d1P,d2P,Pr,Pl,m,SYM)

 if(m==1) then
    call ceno_ur_ul(ext,vz,d1vz,d2vz,vzr,vzl,m,SYM)
    if(Symmetry==2) then
       call ceno_ur_ul(ext,vx,d1vx,d2vx,vxr,vxl,m,ANTI)
       call ceno_ur_ul(ext,vy,d1vy,d2vy,vyr,vyl,m,SYM)
    else if(Symmetry==4) then
       call ceno_ur_ul(ext,vx,d1vx,d2vx,vxr,vxl,m,ANTI)
       call ceno_ur_ul(ext,vy,d1vy,d2vy,vyr,vyl,m,ANTI)
    else
       call ceno_ur_ul(ext,vx,d1vx,d2vx,vxr,vxl,m,SYM)
       call ceno_ur_ul(ext,vy,d1vy,d2vy,vyr,vyl,m,SYM)
    end if
 else if(m==2) then
    if(Symmetry.ne.4) then
       call ceno_ur_ul(ext,vz,d1vz,d2vz,vzr,vzl,m,SYM)
       call ceno_ur_ul(ext,vx,d1vx,d2vx,vxr,vxl,m,SYM)
       if(Symmetry==2) then
          call ceno_ur_ul(ext,vy,d1vy,d2vy,vyr,vyl,m,ANTI)
       else
          call ceno_ur_ul(ext,vy,d1vy,d2vy,vyr,vyl,m,SYM)
       end if
    end if
 else if(m==3) then
    call ceno_ur_ul(ext,vx,d1vx,d2vx,vxr,vxl,m,SYM)
    call ceno_ur_ul(ext,vy,d1vy,d2vy,vyr,vyl,m,SYM)
    if(Symmetry==0) then
        call ceno_ur_ul(ext,vz,d1vz,d2vz,vzr,vzl,m,SYM)
     else
        call ceno_ur_ul(ext,vz,d1vz,d2vz,vzr,vzl,m,ANTI)
     end if
  else
     write(*,*) "m not in proper range"
  end if

  temp = exp(4.d0*phi_f)*(gxx_f*(vxl + betax_f)**2 + &
       2.d0*gxy_f*(vxl + betax_f)*(vyl + betay_f) +         &
       2.d0*gxz_f*(vxl + betax_f)*(vzl + betaz_f) +         &
       gyy_f*(vyl + betay_f)**2 +                           &
       2.d0*gyz_f*(vyl + betay_f)*(vzl + betaz_f) +         &
       gzz_f*(vzl + betaz_f)**2 )/(alpha_f+1.d0)**2
! *** Check for superluminal velocity ***
  where (temp .gt. 1.d0)
    vxl = (vxl + betax_f)*sqrt(fac/temp)-betax_f
    vyl = (vyl + betay_f)*sqrt(fac/temp)-betay_f
    vzl = (vzl + betaz_f)*sqrt(fac/temp)-betaz_f
!!$    vxl = vxl * sqrt(fac/temp)
!!$    vyl = vyl * sqrt(fac/temp)
!!$    vzl = vzl * sqrt(fac/temp)
  end where

  temp = exp(4.d0*phi_f)*(gxx_f*(vxr + betax_f)**2 + &
       2.d0*gxy_f*(vxr + betax_f)*(vyr + betay_f) +         &
       2.d0*gxz_f*(vxr + betax_f)*(vzr + betaz_f) +         &
       gyy_f*(vyr + betay_f)**2 +                           &
       2.d0*gyz_f*(vyr + betay_f)*(vzr + betaz_f) +         &
       gzz_f*(vzr + betaz_f)**2 )/(alpha_f+1.d0)**2
! *** Check for superluminal velocity ***
  where (temp .gt. 1.d0)
    vxr = (vxr + betax_f)*sqrt(fac/temp)-betax_f
    vyr = (vyr + betay_f)*sqrt(fac/temp)-betay_f
    vzr = (vzr + betaz_f)*sqrt(fac/temp)-betaz_f
!!$    vxr = vxr * sqrt(fac/temp)
!!$    vyr = vyr * sqrt(fac/temp)
!!$    vzr = vzr * sqrt(fac/temp)
  end where
end subroutine find_ur_ul_ceno

subroutine ceno_ur_ul(ext,f,d1f,d2f,fr,fl,m,sym)
  implicit none
  integer, dimension(3)                     :: ext
  real*8, dimension(ext(1),ext(2),ext(3))   :: f,d1f,d2f,fr,fl
  integer, intent(in)                       :: m
  real*8                                    :: sym
  integer              :: imin,jmin,kmin,imax,jmax,kmax
  integer              :: i, j, k, k0
  real*8               :: HALF,ONE,FRONT,BACK,f9o8,f3o2,f1o8
  real*8               :: Q(-1:1,-1:1)
  real*8               :: d(-1:1),alpha(-1:1)
  real*8               :: dk0
  parameter(HALF = 0.5D0, ONE = 1.D0, FRONT = 1.D0, BACK = 1.D0)
  parameter(f9o8 = 1.125d0, f3o2 = 1.5d0, f1o8 = 0.125d0) 

 imax = ubound(f,1)
 imin = lbound(f,1)
 jmax = ubound(f,2)
 jmin = lbound(f,2)
 kmax = ubound(f,3)
 kmin = lbound(f,3)
 alpha(-1) = ONE
 alpha(0)  = 0.7D0
 alpha(1)  = ONE

! Note:  Q(k,n) = Q^k_i(x_{i+m/2})
 do k=kmin,kmax
    do j=jmin,jmax
       do i=imin,imax
          
          if(m==1) then
             if(i==imin) then
                Q(-1,-1) =  sym * (f(i,j,k) - HALF * d1f(i,j,k) + f1o8 * d2f(i,j,k))
                Q(-1,1)  =  sym * (f(i,j,k) - f3o2 * d1f(i,j,k) + f9o8 * d2f(i,j,k))
             else
                Q(-1,-1) =  f(i-1,j,k) + HALF * d1f(i-1,j,k) + f1o8 * d2f(i-1,j,k)
                Q(-1,1)  =  f(i-1,j,k) + f3o2 * d1f(i-1,j,k) + f9o8 * d2f(i-1,j,k)
             end if
             Q(0,-1)  =  f(i,j,k) - HALF * d1f(i,j,k) + f1o8 * d2f(i,j,k)
             Q(0,1)   =  f(i,j,k) + HALF * d1f(i,j,k) + f1o8 * d2f(i,j,k)
             if(i==imax) then
                Q(1,-1)  =  Q(0,-1)
                Q(1,1)   =  Q(0,1)
             else
                Q(1,-1)  =  f(i+1,j,k) - f3o2 * d1f(i+1,j,k) + f9o8 * d2f(i+1,j,k)
                Q(1,1)   =  f(i+1,j,k) - HALF * d1f(i+1,j,k) + f1o8 * d2f(i+1,j,k)
             end if
             
          else if(m==2) then
             if(j==jmin) then
                Q(-1,-1) =  sym * (f(i,j,k) - HALF * d1f(i,j,k) + f1o8 * d2f(i,j,k))
                Q(-1,1)  =  sym * (f(i,j,k) - f3o2 * d1f(i,j,k) + f9o8 * d2f(i,j,k))
             else
                Q(-1,-1) =  f(i,j-1,k) + HALF * d1f(i,j-1,k) + f1o8 * d2f(i,j-1,k)
                Q(-1,1)  =  f(i,j-1,k) + f3o2 * d1f(i,j-1,k) + f9o8 * d2f(i,j-1,k)
             end if
             Q(0,-1)  =  f(i,j,k) - HALF * d1f(i,j,k) + f1o8 * d2f(i,j,k)
             Q(0,1)   =  f(i,j,k) + HALF * d1f(i,j,k) + f1o8 * d2f(i,j,k)
             if(j==jmax) then
                Q(1,-1)  =  Q(0,-1)
                Q(1,1)   =  Q(0,1)
             else
                Q(1,-1)  =  f(i,j+1,k) - f3o2 * d1f(i,j+1,k) + f9o8 * d2f(i,j+1,k)
                Q(1,1)   =  f(i,j+1,k) - HALF * d1f(i,j+1,k) + f1o8 * d2f(i,j+1,k)
             end if
             
          else
             if(k==kmin) then
                Q(-1,-1) =  sym * (f(i,j,k) - HALF * d1f(i,j,k) + f1o8 * d2f(i,j,k))
                Q(-1,1)  =  sym * (f(i,j,k) - f3o2 * d1f(i,j,k) + f9o8 * d2f(i,j,k))
             else
                Q(-1,-1) =  f(i,j,k-1) + HALF * d1f(i,j,k-1) + f1o8 * d2f(i,j,k-1)
                Q(-1,1)  =  f(i,j,k-1) + f3o2 * d1f(i,j,k-1) + f9o8 * d2f(i,j,k-1)
             end if
             Q(0,-1)  =  f(i,j,k) - HALF * d1f(i,j,k) + f1o8 * d2f(i,j,k)
             Q(0,1)   =  f(i,j,k) + HALF * d1f(i,j,k) + f1o8 * d2f(i,j,k)
             if(k==kmax) then
                Q(1,-1)  =  Q(0,-1)
                Q(1,1)   =  Q(0,1)
             else
                Q(1,-1)  =  f(i,j,k+1) - f3o2 * d1f(i,j,k+1) + f9o8 * d2f(i,j,k+1)
                Q(1,1)   =  f(i,j,k+1) - HALF * d1f(i,j,k+1) + f1o8 * d2f(i,j,k+1)
             end if

          end if

          if( (Q(-1,-1) - fr(i,j,k))*(Q(0,-1)  - fr(i,j,k)) > 0 .and. &
              (Q(0,-1)  - fr(i,j,k))*(Q(1,-1)  - fr(i,j,k)) > 0 ) then
             if (m==1) then
                if (i==imax) then
                   d(-1)    =  alpha(-1)*abs(Q(-1,-1) - fr(i,j,k))
                   d(0)     =  alpha(0) *abs(Q(0,-1)  - fr(i,j,k))
                   d(1)     =  alpha(1) *abs(Q(1,-1)  - fr(i,j,k))
                else
                   d(-1)    =  alpha(-1)*( abs(Q(-1,-1) - fr(i,j,k)) + abs(Q(-1,1) - fl(i+1,j,k)) )
                   d(0)     =  alpha(0) *( abs(Q(0,-1)  - fr(i,j,k)) + abs(Q(0,1)  - fl(i+1,j,k)) )
                   d(1)     =  alpha(1) *( abs(Q(1,-1)  - fr(i,j,k)) + abs(Q(1,1)  - fl(i+1,j,k)) )
                end if
             else if (m==2) then
                if (j==jmax) then
                   d(-1)    =  alpha(-1)*abs(Q(-1,-1) - fr(i,j,k))
                   d(0)     =  alpha(0) *abs(Q(0,-1)  - fr(i,j,k))
                   d(1)     =  alpha(1) *abs(Q(1,-1)  - fr(i,j,k))
                else
                   d(-1)    =  alpha(-1)*( abs(Q(-1,-1) - fr(i,j,k)) + abs(Q(-1,1) - fl(i,j+1,k)) )
                   d(0)     =  alpha(0) *( abs(Q(0,-1)  - fr(i,j,k)) + abs(Q(0,1)  - fl(i,j+1,k)) )
                   d(1)     =  alpha(1) *( abs(Q(1,-1)  - fr(i,j,k)) + abs(Q(1,1)  - fl(i,j+1,k)) )
                end if
             else
                if (k==kmax) then
                   d(-1)    =  alpha(-1)*abs(Q(-1,-1) - fr(i,j,k))
                   d(0)     =  alpha(0) *abs(Q(0,-1)  - fr(i,j,k))
                   d(1)     =  alpha(1) *abs(Q(1,-1)  - fr(i,j,k))
                else
                   d(-1)    =  alpha(-1)*( abs(Q(-1,-1) - fr(i,j,k)) + abs(Q(-1,1) - fl(i,j,k+1)) )
                   d(0)     =  alpha(0) *( abs(Q(0,-1)  - fr(i,j,k)) + abs(Q(0,1)  - fl(i,j,k+1)) )
                   d(1)     =  alpha(1) *( abs(Q(1,-1)  - fr(i,j,k)) + abs(Q(1,1)  - fl(i,j,k+1)) )
                end if
             end if

             k0 = 0
             dk0 = d(0)
             if(d(-1).lt.dk0) then
                k0 = -1
                dk0 = d(-1)
             end if
             if(d(1).lt.dk0) then
                k0 = 1
                dk0 = d(1)
             end if

             if (m==1) then
                if (i==imax) then
                   fr(i,j,k) = Q(k0,-1)
                else if (i==imin) then
                   fl(i+1,j,k) = Q(k0,1)
                   fr(i,j,k) = Q(k0,-1)
                   if (sym > 0.d0) then
                      fl(i,j,k) = fr(i,j,k)
                   else
                      fr(i,j,k) = 0.d0
                      fl(i,j,k) = 0.d0
                   end if
                else
                   fl(i+1,j,k) = Q(k0,1)
                   fr(i,j,k) = Q(k0,-1)
                end if
             else if (m==2) then
                if (j==jmax) then
                   fr(i,j,k) = Q(k0,-1)
                else if (j==jmin) then
                   fl(i,j+1,k) = Q(k0,1)
                   fr(i,j,k) = Q(k0,-1)
                   if (sym > 0.d0) then
                      fl(i,j,k) = fr(i,j,k)
                   else
                      fr(i,j,k) = 0.d0
                      fl(i,j,k) = 0.d0
                   end if
                else
                   fl(i,j+1,k) = Q(k0,1)
                   fr(i,j,k) = Q(k0,-1)
                end if                
             else 
                if (k==kmax) then
                   fr(i,j,k) = Q(k0,-1)
                else if (k==kmin) then
                   fl(i,j,k+1) = Q(k0,1)
                   fr(i,j,k) = Q(k0,-1)
                   if (sym > 0.d0) then
                      fl(i,j,k) = fr(i,j,k)
                   else
                      fr(i,j,k) = 0.d0
                      fl(i,j,k) = 0.d0
                   end if
                else
                   fl(i,j,k+1) = Q(k0,1)
                   fr(i,j,k) = Q(k0,-1)
                end if                
             end if

          end if
       end do
    end do
 end do

end subroutine ceno_ur_ul

subroutine second_dir(ext,f,d2f,m,sym)
  implicit none
  integer, dimension(3)                     :: ext
  real*8, dimension(ext(1),ext(2),ext(3))   :: f,d2f
  integer, intent(in)                       :: m
  real*8                                    :: sym
 real*8                     :: ONE, TWO
 parameter(ONE = 1.0D0, TWO = 2.D0)

  if(m==1) then
    call ddx2b(ext,f,d2f,ONE,sym)     
  else if(m==2) then
    call ddy2b(ext,f,d2f,ONE,sym)     
  else if(m==3) then
    call ddz2b(ext,f,d2f,ONE,sym)
  else
     write(*,*) "m not in proper range"
  end if

end subroutine second_dir

subroutine hat_flux(ext,f,d2f,m,sym)
  implicit none
  integer, dimension(3)                     :: ext
  real*8, dimension(ext(1),ext(2),ext(3))   :: f,d2f,tild2f
  integer, intent(in)                       :: m
  real*8                                    :: sym,d2min
  real*8                     :: ZERO, ONE, TWO, GAM1
  integer              :: imin,jmin,kmin,imax,jmax,kmax
  integer              :: i, j, k, k0
  real*8               :: alpha(-1:1),mm3
  EXTERNAL mm3
  parameter(ZERO = 0.D0, ONE = 1.0D0, TWO = 2.D0, GAM1 = 1.D0/24.D0)
  imax = ubound(f,1)
  imin = lbound(f,1)
  jmax = ubound(f,2)
  jmin = lbound(f,2)
  kmax = ubound(f,3)
  kmin = lbound(f,3)
  alpha(-1) = ONE
  alpha(0)  = 0.7D0
  alpha(1)  = ONE
  
  do k=kmin,kmax
     do j=jmin,jmax
        do i=imin,imax
           
           if(m==1) then
              if(i==imin) then
                 d2min = sym*d2f(i+1,j,k)
                 tild2f(i,j,k) = mm3(d2min,d2f(i,j,k),d2f(i+1,j,k),alpha)
              else if(i==imax) then
                 tild2f(i,j,k) = d2f(i-1,j,k)
              else
                 tild2f(i,j,k) = mm3(d2f(i-1,j,k),d2f(i,j,k),d2f(i+1,j,k),alpha)
              end if
           else if(m==2) then
              if(j==jmin) then
                 d2min = sym*d2f(i,j+1,k)
                 tild2f(i,j,k) = mm3(d2min,d2f(i,j,k),d2f(i,j+1,k),alpha)
              else if(j==jmax) then
                 tild2f(i,j,k) = d2f(i,j-1,k)                 
              else
                tild2f(i,j,k) = mm3(d2f(i,j-1,k),d2f(i,j,k),d2f(i,j+1,k),alpha)
              end if
           else
              if(k==kmin) then
                 d2min = sym*d2f(i,j,k+1)
                 tild2f(i,j,k) = mm3(d2min,d2f(i,j,k),d2f(i,j,k+1),alpha)
              else if(k==kmax) then
                 tild2f(i,j,k) = d2f(i,j,k-1)
              else
                tild2f(i,j,k) = mm3(d2f(i,j,k-1),d2f(i,j,k),d2f(i,j,k+1),alpha)
              end if
           end if
           
        end do
     end do
  end do

  f = f - GAM1 * tild2f

end subroutine hat_flux


function mm3(a,b,c,alpha)
  implicit none
  real*8  :: mm3,a,b,c,alpha(-1:1),minQ

  if(a*b.gt.0.d0 .and. b*c.gt.0.d0) then
     mm3 = b
     minQ = alpha(0)*b
     if(abs(alpha(-1)*a).lt.abs(minQ)) then
        mm3 = a
        minQ = alpha(-1)*a
     end if
     if(abs(alpha(1)*c).lt.abs(minQ)) then
        mm3 = c
     end if
  else
     mm3 = 0.d0
  end if
  return
end function mm3

subroutine ceno_floor(ext,n,rhob_floor,P_floor,rho_br,rho_bl,Pr,Pl)
  implicit none
  integer, dimension(3)                     :: ext
  real*8, dimension(ext(1),ext(2),ext(3))   :: rhob_floor,P_floor
  real*8, dimension(ext(1),ext(2),ext(3))   :: rho_br, rho_bl, Pr, Pl
  integer                                   :: imin,jmin,kmin,imax,jmax,kmax
  integer                                   :: i, j, k
  real*8, parameter                         :: TINY = 1.D-6, ONE = 1.d0
  real*8, parameter                         :: fact = 0.5d0
  real*8                                    :: n, gam
  imax = ubound(Pr,1)
  imin = lbound(Pr,1)
  jmax = ubound(Pr,2)
  jmin = lbound(Pr,2)
  kmax = ubound(Pr,3)
  kmin = lbound(Pr,3)

  gam = ONE + ONE/n

  do k=kmin,kmax
     do j=jmin,jmax
        do i=imin,imax

!!$           if ((rho_bl(i,j,k) < rhob_floor(i,j,k)) .and. (Pl(i,j,k) < P_floor(i,j,k))) then
!!$              rho_bl(i,j,k) = rhob_floor(i,j,k)
!!$              Pl(i,j,k) = P_floor(i,j,k)
!!$           end if
!!$           
!!$           if ((rho_bl(i,j,k) > rhob_floor(i,j,k)) .and. (Pl(i,j,k) < ((ONE-TINY)*(rho_bl(i,j,k))**gam))) then
!!$              Pl(i,j,k) = (rho_bl(i,j,k))**gam
!!$           end if
!!$           
!!$           if ((rho_bl(i,j,k) < rhob_floor(i,j,k)) .and. (Pl(i,j,k) > P_floor(i,j,k))) then
!!$              rho_bl(i,j,k) = rhob_floor(i,j,k)
!!$              if (Pl(i,j,k) < (ONE-TINY)*(rho_bl(i,j,k))**gam) then
!!$                 Pl(i,j,k) = (rho_bl(i,j,k))**gam
!!$              end if
!!$           end if
!!$           
!!$           if ((rho_br(i,j,k) < rhob_floor(i,j,k)) .and. (Pr(i,j,k) < P_floor(i,j,k))) then
!!$              rho_br(i,j,k) = rhob_floor(i,j,k)
!!$              Pr(i,j,k) = P_floor(i,j,k)
!!$           end if
!!$           
!!$           if ((rho_br(i,j,k) > rhob_floor(i,j,k)) .and. (Pr(i,j,k) < ((ONE-TINY)*(rho_br(i,j,k))**gam))) then
!!$              Pr(i,j,k) = (rho_br(i,j,k))**gam
!!$           end if
!!$           
!!$           if ((rho_br(i,j,k) < rhob_floor(i,j,k)) .and. (Pr(i,j,k) > P_floor(i,j,k))) then
!!$              rho_br(i,j,k) = rhob_floor(i,j,k)
!!$              if (Pr(i,j,k) < (ONE-TINY)*(rho_br(i,j,k))**gam) then
!!$                 Pr(i,j,k) = (rho_br(i,j,k))**gam
!!$              end if
!!$           end if
!!$
!!$           if (Pr(i,j,k) .lt. 0.d0) then
!!$              Pr(i,j,k) = rho_br(i,j,k)**gam
!!$           end if
!!$
!!$           if (Pl(i,j,k) .lt. 0.d0) then
!!$              Pl(i,j,k) = rho_bl(i,j,k)**gam
!!$           end if

           if (Pl(i,j,k) .lt. fact*rho_bl(i,j,k)**gam) then
              Pl(i,j,k) = fact*rho_bl(i,j,k)**gam
           end if

           if (Pr(i,j,k) .lt. fact*rho_br(i,j,k)**gam) then
              Pr(i,j,k) = fact*rho_br(i,j,k)**gam
           end if

        end do
     end do
  end do

end subroutine ceno_floor
