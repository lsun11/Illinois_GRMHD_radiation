!-----------------------------------------------------------------------------
! Apply boundary conditions
!-----------------------------------------------------------------------------
subroutine hbpuncture_shift_bc(ex, fake_ex, dT, dx,dy,dz, &
     x, y, z, r_fish, PhysicalRadius, RadiusDerivative, &
     shiftx_new, shifty_new, shiftz_new, &
     shiftx_old, shifty_old, shiftz_old, &
     dtshiftx_new, dtshifty_new, dtshiftz_new, &
     dtshiftx_old, dtshifty_old, dtshiftz_old, &
     lapse, phi, Symmetry,Xcenter,Ycenter,Zcenter, &
     gxx,gxy,gxz,gyy,gyz,gzz, &
     shift_bc_radial_falloff_power,shift_bc_type, &
     have_bdry_min,have_bdry_max,use_trans_fish_phys_new,fisheye_enable)
  implicit none
  !
  ! Input parameters:
  !
  integer, dimension(3)                    :: ex,fake_ex,have_bdry_min,have_bdry_max
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z,r_fish
  real*8, dimension(ex(1),ex(2),ex(3))     :: PhysicalRadius,RadiusDerivative
  real*8                                   :: dT,dx,dy,dz
  real*8, dimension(ex(1),ex(2),ex(3))     :: phi,lapse
  real*8, dimension(ex(1),ex(2),ex(3))     :: shiftx_old,shifty_old,shiftz_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: shiftx_new,shifty_new,shiftz_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: dtshiftx_old,dtshifty_old,dtshiftz_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: dtshiftx_new,dtshifty_new,dtshiftz_new
  integer                                  :: shift_bc_radial_falloff_power
  integer                                  :: Symmetry,shift_bc_type,use_trans_fish_phys_new,fisheye_enable
  real*8, dimension(ex(1),ex(2),ex(3))     :: gxx,gxy,gxz,gyy,gyz,gzz
  integer                    :: i, j, k
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  integer                    :: il,jl,kl,iu,ju,ku
  integer                    :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  integer                    :: WAVELIKE, FREEZE
  real*8                     :: Dlx,Dly,Dlz,f_int,facr,facRp
  real*8                     :: HALF, ONE, TWO, ZERO
  parameter(HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0)
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(WAVELIKE = 0, FREEZE = 1)
  real*8, parameter          :: SYM = 1.D0, ANT = -1.D0
  real*8                     :: Xcenter, Ycenter, Zcenter
  logical                    :: fish_to_phys
  !
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

  if(shift_bc_type==WAVELIKE) then
     if(fisheye_enable==1) then
        fish_to_phys = .TRUE.
        call trans_fish_phys_vector(ex,fake_ex, X, Y, Z, r_fish, PhysicalRadius, RadiusDerivative, &
             shiftx_old, shifty_old, shiftz_old, &
             fish_to_phys,2,have_bdry_min,have_bdry_max,use_trans_fish_phys_new)
     end if
     ! Transform components of the shift from fisheye to "physical" coordinates
     !     fish_to_phys = .TRUE.
     !     call trans_fish_phys_vector(ex, X, Y, Z, PhysicalR, RDerivative, &
     !                        shiftx_old, shifty_old, shiftz_old, &
     !                        Symmetry,fish_to_phys,2)
     !
     ! x = xmax
     !
     if(have_bdry_max(1)==1) then
        i = imax
        do k = kmin, kmax - 1
           do j = jmin, jmax - 1
              call weights2(x(i,j,k),y(1,j,1),z(i,j,k), &
                   PhysicalRadius(i,j,k),RadiusDerivative(i,j,k), &
                   lapse(i,j,k),phi(i,j,k), &
                   dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                   kmin,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, &
                   Xcenter, Ycenter, Zcenter, &
                   gxx(i,j,k),gxy(i,j,k),gxz(i,j,k), &
                   gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
              call interp(ex,shiftx_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,ANT,SYM,SYM)
              shiftx_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,shifty_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,SYM,ANT,SYM)
              shifty_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,shiftz_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,SYM,SYM,ANT)
              shiftz_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,dtshiftx_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,ANT,SYM,SYM)
              dtshiftx_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,dtshifty_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,SYM,ANT,SYM)
              dtshifty_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,dtshiftz_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,SYM,SYM,ANT)
              dtshiftz_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
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
              call weights2(x(i,j,k),y(1,j,1),z(i,j,k), &
                   PhysicalRadius(i,j,k),RadiusDerivative(i,j,k), &
                   lapse(i,j,k),phi(i,j,k), &
                   dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                   kmin,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, &
                   Xcenter, Ycenter, Zcenter, &
                   gxx(i,j,k),gxy(i,j,k),gxz(i,j,k), &
                   gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
              call interp(ex,shiftx_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,ANT,SYM,SYM)
              shiftx_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,shifty_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,SYM,ANT,SYM)
              shifty_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,shiftz_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,SYM,SYM,ANT)
              shiftz_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,dtshiftx_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,ANT,SYM,SYM)
              dtshiftx_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,dtshifty_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,SYM,ANT,SYM)
              dtshifty_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,dtshiftz_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,SYM,SYM,ANT)
              dtshiftz_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
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
              call weights2(x(i,j,k),y(1,j,1),z(i,j,k), &
                   PhysicalRadius(i,j,k),RadiusDerivative(i,j,k), &
                   lapse(i,j,k),phi(i,j,k), &
                   dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                   kmin,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, &
                   Xcenter, Ycenter, Zcenter, &
                   gxx(i,j,k),gxy(i,j,k),gxz(i,j,k), &
                   gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
              call interp(ex,shiftx_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,ANT,SYM,SYM)
              shiftx_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              if(Symmetry==AXISYM) then
                 call interp(ex,shifty_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,ANT,ANT,SYM)
              else
                 call interp(ex,shifty_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,SYM,ANT,SYM)
              end if
              shifty_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,shiftz_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,SYM,SYM,ANT)
              shiftz_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,dtshiftx_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,ANT,SYM,SYM)
              dtshiftx_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              if(Symmetry==AXISYM) then
                 call interp(ex,dtshifty_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,ANT,ANT,SYM)
              else
                 call interp(ex,dtshifty_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,SYM,ANT,SYM)
              end if
              dtshifty_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,dtshiftz_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,SYM,SYM,ANT)
              dtshiftz_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
           end do
        end do
     end if
     !
     ! z=zmin
     !
     if(have_bdry_min(3)==1) then
        k = kmin
        do j = jmin, jmax
           do i = imin, imax
              call weights2(x(i,j,k),y(1,j,1),z(i,j,k), &
                   PhysicalRadius(i,j,k),RadiusDerivative(i,j,k), &
                   lapse(i,j,k),phi(i,j,k), &
                   dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                   kmin,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, &
                   Xcenter, Ycenter, Zcenter, &
                   gxx(i,j,k),gxy(i,j,k),gxz(i,j,k), &
                   gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
              call interp(ex,shiftx_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,ANT,SYM,SYM)
              shiftx_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              if(Symmetry==AXISYM) then
                 call interp(ex,shifty_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,ANT,ANT,SYM)
              else
                 call interp(ex,shifty_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,SYM,ANT,SYM)
              end if
              shifty_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,shiftz_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,SYM,SYM,ANT)
              shiftz_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,dtshiftx_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,ANT,SYM,SYM)
              dtshiftx_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              if(Symmetry==AXISYM) then
                 call interp(ex,dtshifty_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,ANT,ANT,SYM)
              else
                 call interp(ex,dtshifty_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,SYM,ANT,SYM)
              end if
              dtshifty_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,dtshiftz_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,SYM,SYM,ANT)
              dtshiftz_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
           end do
        end do
     endif


     !
     ! x = xmin
     !
     if(have_bdry_min(1)==1) then
        i = imin
        do k = kmin, kmax - 1
           do j = jmin, jmax - 1
              call weights2(x(i,j,k),y(1,j,1),z(i,j,k), &
                   PhysicalRadius(i,j,k),RadiusDerivative(i,j,k), &
                   lapse(i,j,k),phi(i,j,k), &
                   dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                   kmin,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, &
                   Xcenter, Ycenter, Zcenter, &
                   gxx(i,j,k),gxy(i,j,k),gxz(i,j,k), &
                   gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
              call interp(ex,shiftx_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,ANT,SYM,SYM)
              shiftx_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,shifty_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,SYM,ANT,SYM)
              shifty_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,shiftz_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,SYM,SYM,ANT)
              shiftz_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,dtshiftx_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,ANT,SYM,SYM)
              dtshiftx_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,dtshifty_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,SYM,ANT,SYM)
              dtshifty_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,dtshiftz_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,SYM,SYM,ANT)
              dtshiftz_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
           end do
        end do
     end if

     !
     ! y = ymin
     !
     if(have_bdry_min(2)==1) then
        j = jmin
        do k = kmin, kmax - 1
           do i = imin+1, imax - 1
              !           Dude!  We've already updated i=imin!
              !           do i = imin, imax - 1
              call weights2(x(i,j,k),y(1,j,1),z(i,j,k), &
                   PhysicalRadius(i,j,k),RadiusDerivative(i,j,k), &
                   lapse(i,j,k),phi(i,j,k), &
                   dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                   kmin,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, &
                   Xcenter, Ycenter, Zcenter, &
                   gxx(i,j,k),gxy(i,j,k),gxz(i,j,k), &
                   gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
              call interp(ex,shiftx_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,ANT,SYM,SYM)
              shiftx_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,shifty_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,SYM,ANT,SYM)
              shifty_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,shiftz_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,SYM,SYM,ANT)
              shiftz_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,dtshiftx_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,ANT,SYM,SYM)
              dtshiftx_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,dtshifty_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,SYM,ANT,SYM)
              dtshifty_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
              call interp(ex,dtshiftz_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,SYM,SYM,ANT)
              dtshiftz_new(i,j,k) = f_int * facRp**(shift_bc_radial_falloff_power)
           end do
        end do
     end if

     ! Transform the components of the shift back to fisheye coordinates
     if(fisheye_enable==1) then
        fish_to_phys = .FALSE.
        call trans_fish_phys_vector(ex, fake_ex, X, Y, Z, r_fish, PhysicalRadius, RadiusDerivative, &
             shiftx_old, shifty_old, shiftz_old, &
             fish_to_phys,2,have_bdry_min,have_bdry_max,use_trans_fish_phys_new)
        
        call trans_fish_phys_vector(ex, fake_ex,X, Y, Z, r_fish, PhysicalRadius, RadiusDerivative, &
             shiftx_new, shifty_new, shiftz_new, &
             fish_to_phys,0,have_bdry_min,have_bdry_max,use_trans_fish_phys_new)
     end if

  else if(shift_bc_type==FREEZE) then
     !
     ! x = xmax
     !
     if(have_bdry_max(1)==1) then
        shiftx_new(imax,:,:) = shiftx_old(imax,:,:)
        shifty_new(imax,:,:) = shifty_old(imax,:,:)
        shiftz_new(imax,:,:) = shiftz_old(imax,:,:)
        dtshiftx_new(imax,:,:) = dtshiftx_old(imax,:,:)
        dtshifty_new(imax,:,:) = dtshifty_old(imax,:,:)
        dtshiftz_new(imax,:,:) = dtshiftz_old(imax,:,:)
     end if
     !
     ! y = ymax
     !
     if(have_bdry_max(2)==1) then
        shiftx_new(:,jmax,:) = shiftx_old(:,jmax,:)
        shifty_new(:,jmax,:) = shifty_old(:,jmax,:)
        shiftz_new(:,jmax,:) = shiftz_old(:,jmax,:)
        dtshiftx_new(:,jmax,:) = dtshiftx_old(:,jmax,:)
        dtshifty_new(:,jmax,:) = dtshifty_old(:,jmax,:)
        dtshiftz_new(:,jmax,:) = dtshiftz_old(:,jmax,:)
     end if
     !
     ! z = zmax
     !
     if(have_bdry_max(3)==1) then
        shiftx_new(:,:,kmax) = shiftx_old(:,:,kmax)
        shifty_new(:,:,kmax) = shifty_old(:,:,kmax)
        shiftz_new(:,:,kmax) = shiftz_old(:,:,kmax)
        dtshiftx_new(:,:,kmax) = dtshiftx_old(:,:,kmax)
        dtshifty_new(:,:,kmax) = dtshifty_old(:,:,kmax)
        dtshiftz_new(:,:,kmax) = dtshiftz_old(:,:,kmax)
     end if
     !
     ! x = xmin
     !
     if(have_bdry_min(1)==1) then
        shiftx_new(imin,:,:) = shiftx_old(imin,:,:)
        shifty_new(imin,:,:) = shifty_old(imin,:,:)
        shiftz_new(imin,:,:) = shiftz_old(imin,:,:)
        dtshiftx_new(imin,:,:) = dtshiftx_old(imin,:,:)
        dtshifty_new(imin,:,:) = dtshifty_old(imin,:,:)
        dtshiftz_new(imin,:,:) = dtshiftz_old(imin,:,:)
     end if

     !
     ! y = ymin
     !
     if(have_bdry_min(2)==1) then
        shiftx_new(:,jmin,:) = shiftx_old(:,jmin,:)
        shifty_new(:,jmin,:) = shifty_old(:,jmin,:)
        shiftz_new(:,jmin,:) = shiftz_old(:,jmin,:)
        dtshiftx_new(:,jmin,:) = dtshiftx_old(:,jmin,:)
        dtshifty_new(:,jmin,:) = dtshifty_old(:,jmin,:)
        dtshiftz_new(:,jmin,:) = dtshiftz_old(:,jmin,:)
     end if

     !
     ! z = zmin
     !
     if(have_bdry_min(3)==1) then
        shiftx_new(:,:,kmin) = shiftx_old(:,:,kmin)
        shifty_new(:,:,kmin) = shifty_old(:,:,kmin)
        shiftz_new(:,:,kmin) = shiftz_old(:,:,kmin)
        dtshiftx_new(:,:,kmin) = dtshiftx_old(:,:,kmin)
        dtshifty_new(:,:,kmin) = dtshifty_old(:,:,kmin)
        dtshiftz_new(:,:,kmin) = dtshiftz_old(:,:,kmin)
     end if

  end if
end subroutine hbpuncture_shift_bc
