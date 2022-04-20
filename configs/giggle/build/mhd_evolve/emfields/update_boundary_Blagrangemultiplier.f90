!----------------------------------------------------------------------------------------
! Update outer boundaries for hyperbolic divergence cleaning variable Blagrangemultiplier
!----------------------------------------------------------------------------------------
subroutine update_boundary_Blagrangemultiplier(ex,fake_ex,X,Y,Z, PhysicalR, RadiusDerivative, &
     dT,dX,dY,dZ, alpha, &
     phi_old,gxx_old,gxy_old,gxz_old,gyy_old,gyz_old,gzz_old, &
     Blagrangemultiplier_old,Blagrangemultiplier_new, &
     Symmetry, &
     have_bdry_min,have_bdry_max)
  implicit none
  !
  ! Input parameters:
  !
  integer, dimension(3)                    :: ex,fake_ex,have_bdry_min,have_bdry_max
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: PhysicalR,RadiusDerivative
  real*8                                   :: dT,dX,dY,dZ
  real*8                                   :: Xcenter,Ycenter,Zcenter
  real*8, dimension(ex(1),ex(2),ex(3))     :: alpha
  real*8, dimension(ex(1),ex(2),ex(3))     :: Blagrangemultiplier_old,Blagrangemultiplier_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: phi_old,gxx_old,gxy_old,gxz_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: gyy_old,gyz_old,gzz_old
  integer                                  :: Symmetry
  !
  ! Other variables:
  !
  integer                    :: i, j, k
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  integer                    :: il,jl,kl,iu,ju,ku
  real*8                     :: Dlx,Dly,Dlz,r,r1,Rp,Rp1,r2
  real*8                     :: facr,facRp,fac,fac1,fac2,rp2
  real*8                     :: f_int_Blagrangemultiplier
  real*8                     :: HALF, ONE, TWO, THREE, ZERO
  real*8, parameter          :: f1o3=1.d0/3.d0,f1o6=1.d0/6.d0
  parameter(HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0, THREE = 3.D0, ZERO = 0.D0)
  real*8, parameter          :: SYM = 1.D0, ANT = -1.D0
  ! The following were added for the purpose of axisymmetrization:
  integer                    :: jlb, jub
  logical                      :: fish_to_phys
  imin = ex(1)-fake_ex(1)+1
  jmin = ex(2)-fake_ex(2)+1
  kmin = ex(3)-fake_ex(3)+1
  imax = fake_ex(1)
  jmax = fake_ex(2)
  kmax = fake_ex(3)
  !Following is USEFUL when nghostzone > 1 !!!!
  !write(*,*) Calling BBBBupd_bdry with imin,imax =,imin,imax
  ! Change:  Set bounds on y loops
  !WARNING: FOLLOWING WONT WORK IN AXISYMMETRY.  Have modified below code to ignore jlb,jub, which is necessary in axisymmetry
  jlb = jmin
  jub = jmax - 1
  Xcenter = 0.D0
  Ycenter = 0.D0
  Zcenter = 0.D0
!!$
!!$  if (have_bdry_max(1)==1) then
!!$     Blagrangemultiplier_new(imax,:,:) = Blagrangemultiplier_new(imax-1,:,:)
!!$  end if
!!$
!!$  if (have_bdry_min(1)==1) then
!!$     Blagrangemultiplier_new(imin,:,:) = Blagrangemultiplier_new(imin+1,:,:)
!!$  end if
!!$
!!$  if (have_bdry_max(2)==1) then
!!$     Blagrangemultiplier_new(:,jmax,:) = Blagrangemultiplier_new(:,jmax-1,:)
!!$  end if
!!$
!!$  if (have_bdry_min(2)==1) then
!!$     Blagrangemultiplier_new(:,jmin,:) = Blagrangemultiplier_new(:,jmin+1,:)
!!$  end if
!!$
!!$  if (have_bdry_max(3)==1) then
!!$     Blagrangemultiplier_new(:,:,kmax) = Blagrangemultiplier_new(:,:,kmax-1)
!!$  end if
!!$
!!$  if (have_bdry_min(3)==1) then
!!$     Blagrangemultiplier_new(:,:,kmin) = Blagrangemultiplier_new(:,:,kmin+1)
!!$  end if
!!$
  !
  ! x = max
  if(have_bdry_max(1)==1) then
     i = imax
!!!$omp parallel do
     do k = kmin, kmax - 1
        do j = jmin, jmax - 1
           call weights2(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),RadiusDerivative(i,j,k), &
                & alpha(i,j,k),phi_old(i,j,k),&
                & dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il&
                &,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, Xcenter,&
                & Ycenter, Zcenter,gxx_old(i,j,k),gxy_old(i,j,k),gxz_old(i,j,k),gyy_old(i,j,k),gyz_old(i,j,k),gzz_old(i,j,k))
           !
           ! Interpolate 
           !
           call interp(ex,Blagrangemultiplier_old  ,f_int_Blagrangemultiplier  ,il,iu,jl,ju,kl,ku,Dlx&
                &,Dly,Dlz,SYM,SYM,SYM)
           r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
           r1 = facr*r
           Rp = PhysicalR(i,j,k)
           Rp1 = facRp*Rp
           Blagrangemultiplier_new(i,j,k)   = f_int_Blagrangemultiplier * facRp
           !if(k==kmin .and. j==jmin+1) write(*,*) max,facRp,Dlx,Dly,Dlz,ONE - (alpha(i,j,k) + ONE)*exp(-TWO*phi_old(i,j,k))*dT/(sqr
        enddo
     enddo
!!!$omp end parallel do
  end if
  !
  ! y = ymax.  Do not do this for axisymmetric cases.
  !
  if(have_bdry_max(2)==1) then
     j = jmax
!!!$omp parallel do
     do k = kmin, kmax - 1
        do i = imin, imax
           call weights2(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),RadiusDerivative(i,j,k),alpha(i,j,k),phi_old(i,j,k),&
                & dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il&
                &,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, Xcenter,&
                & Ycenter, Zcenter,gxx_old(i,j,k),gxy_old(i,j,k),gxz_old(i,j,k),gyy_old(i,j,k),gyz_old(i,j,k),gzz_old(i,j,k))
           !
           ! Interpolate 
           !
           call interp(ex,Blagrangemultiplier_old  ,f_int_Blagrangemultiplier  ,il,iu,jl,ju,kl,ku,Dlx&
                &,Dly,Dlz,SYM,SYM,SYM)
           r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
           r1 = facr*r
           Rp = PhysicalR(i,j,k)
           Rp1 = facRp*Rp
           Blagrangemultiplier_new(i,j,k)   = f_int_Blagrangemultiplier * facRp
        enddo
     enddo
!!!$omp end parallel do
  endif
  !
  ! z = zmax
  !
  !WARNING: FOLLOWING WONT WORK IN AXISYMMETRY
  !jub = jub + 1
  if(have_bdry_max(3)==1) then
     k = kmax
!!!$omp parallel do
     do j = jmin, jmax
        do i = imin, imax
           call weights2(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),RadiusDerivative(i,j,k),alpha(i,j,k),phi_old(i,j,k),&
                & dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il&
                &,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, Xcenter,&
                & Ycenter, Zcenter,gxx_old(i,j,k),gxy_old(i,j,k),gxz_old(i,j,k),gyy_old(i,j,k),gyz_old(i,j,k),gzz_old(i,j,k))
           !
           ! Interpolate 
           !
           call interp(ex,Blagrangemultiplier_old  ,f_int_Blagrangemultiplier  ,il,iu,jl,ju,kl,ku,Dlx&
                &,Dly,Dlz,SYM,SYM,SYM)
           r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
           r1 = facr*r
           Rp = PhysicalR(i,j,k)
           Rp1 = facRp*Rp
           Blagrangemultiplier_new(i,j,k)   = f_int_Blagrangemultiplier * facRp
        enddo
     enddo
!!!$omp end parallel do
  end if
  if(have_bdry_min(1)==1) then
     !
     ! x = xmin
     !
     i = imin
!!!$omp parallel do
     do k = kmin, kmax-1
        do j = jmin, jmax-1
           call weights2(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),RadiusDerivative(i,j,k),alpha(i,j,k),phi_old(i,j,k)&
                &, dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin&
                &,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry,&
                & Xcenter, Ycenter, Zcenter,gxx_old(i,j,k),gxy_old(i,j,k),gxz_old(i,j,k),gyy_old(i,j,k),gyz_old(i,j,k),gzz_old(i,j,&
  &k))
           !
           ! Interpolate 
           !
           call interp(ex,Blagrangemultiplier_old  ,f_int_Blagrangemultiplier  ,il,iu,jl,ju,kl,ku&
                &,Dlx,Dly,Dlz,SYM,SYM,SYM)
           r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
           r1 = facr*r
           Rp = PhysicalR(i,j,k)
           Rp1 = facRp*Rp
           Blagrangemultiplier_new(i,j,k)   = f_int_Blagrangemultiplier * facRp
           !if(k==kmin .and. j==jmin+1) write(*,*) min,facRp,Dlx,Dly,Dlz,ONE - (alpha(i,j,k) + ONE)*exp(-TWO*phi_old(i,j,k))*dT/(sqr
        enddo
     enddo
!!!$omp end parallel do
  endif
  if(have_bdry_min(2)==1) then
     !
     ! y = ymin
     !
     j = jmin
!!!$omp parallel do
     do k = kmin, kmax-1
        do i = imin +1, imax-1
           call weights2(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),RadiusDerivative(i,j,k),alpha(i,j,k),phi_old(i,j,k)&
                &, dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin&
                &,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry,&
                & Xcenter, Ycenter, Zcenter,gxx_old(i,j,k),gxy_old(i,j,k),gxz_old(i,j,k),gyy_old(i,j,k),gyz_old(i,j,k),gzz_old(i,j,&
  &k))
           !
           ! Interpolate 
           !
           call interp(ex,Blagrangemultiplier_old  ,f_int_Blagrangemultiplier  ,il,iu,jl,ju,kl,ku,Dlx&
                &,Dly,Dlz,SYM,SYM,SYM)
           r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
           r1 = facr*r
           Rp = PhysicalR(i,j,k)
           Rp1 = facRp*Rp
           Blagrangemultiplier_new(i,j,k)   = f_int_Blagrangemultiplier * facRp
        enddo
     enddo
!!!$omp end parallel do
  endif
  if(have_bdry_min(3)==1) then
     !
     ! z = zmin
     !
     k = kmin
     do j = jmin+1, jmax-1
        do i = imin+1, imax-1
           call weights2(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),RadiusDerivative(i,j,k),alpha(i,j,k),phi_old(i,j,k)&
                &, dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin&
                &,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry,&
                & Xcenter, Ycenter, Zcenter,gxx_old(i,j,k),gxy_old(i,j,k),gxz_old(i,j,k),gyy_old(i,j,k),gyz_old(i,j,k),gzz_old(i,j,&
  &k))
           !
           ! Interpolate 
           !
           call interp(ex,Blagrangemultiplier_old  ,f_int_Blagrangemultiplier  ,il,iu,jl,ju,kl,ku,Dlx&
                &,Dly,Dlz,SYM,SYM,SYM)
           r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
           r1 = facr*r
           Rp = PhysicalR(i,j,k)
           Rp1 = facRp*Rp
           Blagrangemultiplier_new(i,j,k)   = f_int_Blagrangemultiplier * facRp
        enddo
     enddo
  endif
end subroutine update_boundary_Blagrangemultiplier
