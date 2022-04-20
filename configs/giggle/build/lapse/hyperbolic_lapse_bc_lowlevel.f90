!-----------------------------------------------------------------------------
!
! Apply boundary conditions
!
!-----------------------------------------------------------------------------
subroutine hyper_lapse_bc(ex, fake_ex, dT, x, y, z, PhysicalRadius, RadiusDerivative, &
     lapse_new, dtlapse_new, &
     lapse, dtlapse, phi, &
     Xcenter,Ycenter,Zcenter,Symmetry, &
     gxx, gxy, gxz, gyy, gyz, gzz, &
     have_bdry_min,have_bdry_max)
  implicit none
  !
  ! Input parameters:
  !
  integer, dimension(3)                    :: ex, fake_ex, have_bdry_min, have_bdry_max
  real*8                                   :: dT
  real*8, dimension(ex(1),ex(2),ex(3))     :: PhysicalRadius, RadiusDerivative
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: phi,lapse,dtlapse
  real*8, dimension(ex(1),ex(2),ex(3))     :: lapse_new,dtlapse_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: gxx,gxy,gxz,gyy,gyz,gzz
  integer                                  :: Symmetry
  integer                    :: i, j, k
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  integer                    :: il,jl,kl,iu,ju,ku
  integer                    :: NO_SYMM, EQUATORIAL, OCTANT
  integer                    :: PI_SYMM, AXISYM
  real*8                     :: Courant,dX,dY,dZ,Dlx,Dly,Dlz,f_int,facr,facRp
  real*8                     :: HALF, ONE, TWO, ZERO
  integer                    :: bc
  parameter(HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0)
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2)
  parameter(PI_SYMM = 3, AXISYM = 4)
  real*8, parameter          :: SYM = 1.D0, ANT = -1.D0
  real*8                     :: Xcenter, Ycenter, Zcenter
  Xcenter = ZERO
  Ycenter = ZERO
  Zcenter = ZERO
  !
  imin = ex(1)-fake_ex(1)+1
  jmin = ex(2)-fake_ex(2)+1
  kmin = ex(3)-fake_ex(3)+1
  imax = fake_ex(1)
  jmax = fake_ex(2)
  kmax = fake_ex(3)
  ! Following probably needs to be fixed for Carpet compatibility:
  dX = x(imin+1,1,1) - x(imin,1,1)
  dY = y(1,jmin+1,1) - y(1,jmin,1)
  dZ = z(1,1,kmin+1) - z(1,1,kmin)
  !
  ! x = xmax
  !
  if(have_bdry_max(1)==1) then
     i = imax
     do k = kmin, kmax - 1
        do j = jmin, jmax - 1
           call weights2(x(i,j,k),y(i,j,k),z(i,j,k),PhysicalRadius(i,j,k),RadiusDerivative(i,j,k), &
                lapse(i,j,k),phi(i,j,k), &
                dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                kmin,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, &
                Xcenter, Ycenter, Zcenter, &
                gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
           call interp(ex,lapse,f_int,il,iu,jl,ju,kl,ku, &
                Dlx,Dly,Dlz,SYM,SYM,SYM)
           lapse_new(i,j,k) = f_int * facRp
           call interp(ex,dtlapse,f_int,il,iu,jl,ju,kl,ku, &
                Dlx,Dly,Dlz,SYM,SYM,SYM)
           dtlapse_new(i,j,k) = f_int * facRp
        end do
     end do
  end if
  !
  ! y = ymax
  !
  if(have_bdry_max(2)==1) then
     j = jmax
     do k = kmin, kmax - 1
        do i = imin, imax
           call weights2(x(i,j,k),y(i,j,k),z(i,j,k),PhysicalRadius(i,j,k),RadiusDerivative(i,j,k), &
                lapse(i,j,k),phi(i,j,k), &
                dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                kmin,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, &
                Xcenter, Ycenter, Zcenter, &
                gxx(i,j,k),gxy(i,j,k),gxz(i,j,k), &
                gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
           call interp(ex,lapse,f_int,il,iu,jl,ju,kl,ku, &
                Dlx,Dly,Dlz,SYM,SYM,SYM)
           lapse_new(i,j,k) = f_int * facRp
           call interp(ex,dtlapse,f_int,il,iu,jl,ju,kl,ku, &
                Dlx,Dly,Dlz,SYM,SYM,SYM)
           dtlapse_new(i,j,k) = f_int * facRp
        end do
     end do
  end if
  !
  ! z = zmax
  !
  if(have_bdry_max(3)==1) then
     k = kmax
     do j = jmin, jmax
        do i = imin, imax
           call weights2(x(i,j,k),y(i,j,k),z(i,j,k),PhysicalRadius(i,j,k),RadiusDerivative(i,j,k),&
                lapse(i,j,k),phi(i,j,k), &
                dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                kmin,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, &
                Xcenter, Ycenter, Zcenter, &
                gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
           call interp(ex,lapse,f_int,il,iu,jl,ju,kl,ku, &
                Dlx,Dly,Dlz,SYM,SYM,SYM)
           lapse_new(i,j,k) = f_int * facRp
           call interp(ex,dtlapse,f_int,il,iu,jl,ju,kl,ku, &
                Dlx,Dly,Dlz,SYM,SYM,SYM)
           dtlapse_new(i,j,k) = f_int * facRp
        end do
     end do
  end if
  !
  ! x = xmin
  !
  if(have_bdry_min(1)==1) then
     i = imin
     do k = kmin, kmax - 1
        do j = jmin, jmax - 1
           call weights2(x(i,j,k),y(i,j,k),z(i,j,k),PhysicalRadius(i,j,k),RadiusDerivative(i,j,k),&
                lapse(i,j,k),phi(i,j,k), &
                dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                kmin,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, &
                Xcenter, Ycenter, Zcenter, &
                gxx(i,j,k),gxy(i,j,k),gxz(i,j,k), &
                gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
           call interp(ex,lapse,f_int,il,iu,jl,ju,kl,ku, &
                Dlx,Dly,Dlz,SYM,SYM,SYM)
           lapse_new(i,j,k) = f_int * facRp
           call interp(ex,dtlapse,f_int,il,iu,jl,ju,kl,ku, &
                Dlx,Dly,Dlz,SYM,SYM,SYM)
           dtlapse_new(i,j,k) = f_int * facRp
        end do
     end do
  end if
  !
  ! y = ymin
  !
  if(have_bdry_min(2)==1) then
     j = jmin
     do k = kmin, kmax - 1
        do i = imin, imax - 1
           call weights2(x(i,j,k),y(i,j,k),z(i,j,k),PhysicalRadius(i,j,k),RadiusDerivative(i,j,k),&
                lapse(i,j,k),phi(i,j,k), &
                dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                kmin,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, &
                Xcenter, Ycenter, Zcenter, &
                gxx(i,j,k),gxy(i,j,k),gxz(i,j,k), &
                gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
           call interp(ex,lapse,f_int,il,iu,jl,ju,kl,ku, &
                Dlx,Dly,Dlz,SYM,SYM,SYM)
           lapse_new(i,j,k) = f_int * facRp
           call interp(ex,dtlapse,f_int,il,iu,jl,ju,kl,ku, &
                Dlx,Dly,Dlz,SYM,SYM,SYM)
           dtlapse_new(i,j,k) = f_int * facRp
        end do
     end do
  end if
  !
  ! z = zmin
  !
  if(have_bdry_min(3)==1) then
     k = kmin
     do j = jmin, jmax
        do i = imin, imax
           call weights2(x(i,j,k),y(i,j,k),z(i,j,k),PhysicalRadius(i,j,k),RadiusDerivative(i,j,k),&
                lapse(i,j,k),phi(i,j,k), &
                dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                kmin,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, &
                Xcenter, Ycenter, Zcenter, &
                gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
           call interp(ex,lapse,f_int,il,iu,jl,ju,kl,ku, &
                Dlx,Dly,Dlz,SYM,SYM,SYM)
           lapse_new(i,j,k) = f_int * facRp
           call interp(ex,dtlapse,f_int,il,iu,jl,ju,kl,ku, &
                Dlx,Dly,Dlz,SYM,SYM,SYM)
           dtlapse_new(i,j,k) = f_int * facRp
        end do
     end do
  endif
end subroutine hyper_lapse_bc
