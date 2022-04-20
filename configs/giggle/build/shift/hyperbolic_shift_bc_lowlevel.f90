!------------------------------------------
! Boundary conditions for hyperbolic shift
!------------------------------------------
subroutine hyper_shift_bc(ex,fake_ex, dT,dx,dy,dz, &
     x, y, z, r_fish, PhysicalRadius, RadiusDerivative, &
     shiftx_new, shifty_new, shiftz_new, &
     shiftx_old, shifty_old, shiftz_old, &
     lapse, phi, Symmetry,Xcenter,Ycenter,Zcenter, &
     gxx,gxy,gxz,gyy,gyz,gzz,shift_bc_type, &
     have_bdry_min,have_bdry_max,use_trans_fish_phys_new,fisheye_enable)
  implicit none
  !
  ! Input parameters:
  !
  integer, dimension(3)                    :: ex,fake_ex,have_bdry_min,have_bdry_max
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z,r_fish
  real*8, dimension(ex(1),ex(2),ex(3))     :: PhysicalRadius, RadiusDerivative
  real*8                                   :: dT,dx,dy,dz
  real*8, dimension(ex(1),ex(2),ex(3))     :: phi,lapse
  real*8, dimension(ex(1),ex(2),ex(3))     :: shiftx_old,shifty_old,shiftz_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: shiftx_new,shifty_new,shiftz_new
  integer                                  :: Symmetry,use_trans_fish_phys_new,fisheye_enable
  real*8, dimension(ex(1),ex(2),ex(3))     :: gxx,gxy,gxz,gyy,gyz,gzz
  integer                    :: i, j, k
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  integer                    :: il,jl,kl,iu,ju,ku
  integer                    :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  integer                    :: WAVELIKE, FREEZE
  real*8                     :: Dlx,Dly,Dlz,f_int,facr,facRp
  real*8                     :: HALF, ONE, TWO, ZERO
  integer                    :: shift_bc_type
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
     ! Transform components of the shift from fisheye to physical coordinates
     if(fisheye_enable==1) then
        fish_to_phys = .TRUE.
        call trans_fish_phys_vector(ex,fake_ex, X, Y, Z, r_fish,PhysicalRadius, RadiusDerivative, &
             shiftx_old, shifty_old, shiftz_old, &
             fish_to_phys,2,have_bdry_min,have_bdry_max,use_trans_fish_phys_new)
     end if
     !
     ! x = xmax
     !
     if(have_bdry_max(1)==1) then
        i = imax
        do k = kmin, kmax - 1
           do j = jmin, jmax - 1
              call weights2(x(i,j,k),y(1,j,1),z(i,j,k),PhysicalRadius(i,j,k),RadiusDerivative(i,j,k), &
                   lapse(i,j,k),phi(i,j,k), &
                   dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                   kmin,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, &
                   Xcenter, Ycenter, Zcenter, &
                   gxx(i,j,k),gxy(i,j,k),gxz(i,j,k), &
                   gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
              if (Symmetry.ne.PI_SYMM) then
                 call interp(ex,shiftx_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,ANT,SYM,SYM)
                 shiftx_new(i,j,k) = f_int * facRp * facRp
                 call interp(ex,shifty_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,SYM,ANT,SYM)
                 shifty_new(i,j,k) = f_int * facRp * facRp
                 call interp(ex,shiftz_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,SYM,SYM,ANT)
                 shiftz_new(i,j,k) = f_int * facRp * facRp
              end if
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
              call weights2(x(i,j,k),y(1,j,1),z(i,j,k), PhysicalRadius(i,j,k),RadiusDerivative(i,j,k),&
                   lapse(i,j,k),phi(i,j,k), &
                   dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                   kmin,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, &
                   Xcenter, Ycenter, Zcenter, &
                   gxx(i,j,k),gxy(i,j,k),gxz(i,j,k), &
                   gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
              if (Symmetry.ne.PI_SYMM) then
                 call interp(ex,shiftx_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,ANT,SYM,SYM)
                 shiftx_new(i,j,k) = f_int * facRp * facRp
                 call interp(ex,shifty_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,SYM,ANT,SYM)
                 shifty_new(i,j,k) = f_int * facRp * facRp
                 call interp(ex,shiftz_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,SYM,SYM,ANT)
                 shiftz_new(i,j,k) = f_int * facRp * facRp
              end if
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
              call weights2(x(i,j,k),y(1,j,1),z(i,j,k), PhysicalRadius(i,j,k),RadiusDerivative(i,j,k),&
                   lapse(i,j,k),phi(i,j,k), &
                   dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                   kmin,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, &
                   Xcenter, Ycenter, Zcenter, &
                   gxx(i,j,k),gxy(i,j,k),gxz(i,j,k), &
                   gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
              if (Symmetry.ne.PI_SYMM) then
                 call interp(ex,shiftx_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,ANT,SYM,SYM)
                 shiftx_new(i,j,k) = f_int * facRp * facRp
                 if(Symmetry==AXISYM) then
                    call interp(ex,shifty_old,f_int,il,iu,jl,ju,kl,ku, &
                         Dlx,Dly,Dlz,ANT,ANT,SYM)
                 else
                    call interp(ex,shifty_old,f_int,il,iu,jl,ju,kl,ku, &
                         Dlx,Dly,Dlz,SYM,ANT,SYM)
                 end if
                 shifty_new(i,j,k) = f_int * facRp * facRp
                 call interp(ex,shiftz_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,SYM,SYM,ANT)
                 shiftz_new(i,j,k) = f_int * facRp * facRp
              end if
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
              call weights2(x(i,j,k),y(1,j,1),z(i,j,k),PhysicalRadius(i,j,k),RadiusDerivative(i,j,k), &
                   lapse(i,j,k),phi(i,j,k), &
                   dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                   kmin,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, &
                   Xcenter, Ycenter, Zcenter, &
                   gxx(i,j,k),gxy(i,j,k),gxz(i,j,k), &
                   gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
              if (Symmetry.ne.PI_SYMM) then
                 call interp(ex,shiftx_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,ANT,SYM,SYM)
                 shiftx_new(i,j,k) = f_int * facRp * facRp
                 call interp(ex,shifty_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,SYM,ANT,SYM)
                 shifty_new(i,j,k) = f_int * facRp * facRp
                 call interp(ex,shiftz_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,SYM,SYM,ANT)
                 shiftz_new(i,j,k) = f_int * facRp * facRp
              end if
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
              call weights2(x(i,j,k),y(1,j,1),z(i,j,k), PhysicalRadius(i,j,k),RadiusDerivative(i,j,k),&
                   lapse(i,j,k),phi(i,j,k), &
                   dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                   kmin,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, &
                   Xcenter, Ycenter, Zcenter, &
                   gxx(i,j,k),gxy(i,j,k),gxz(i,j,k), &
                   gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
              if (Symmetry.ne.PI_SYMM) then
                 call interp(ex,shiftx_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,ANT,SYM,SYM)
                 shiftx_new(i,j,k) = f_int * facRp * facRp
                 call interp(ex,shifty_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,SYM,ANT,SYM)
                 shifty_new(i,j,k) = f_int * facRp * facRp
                 call interp(ex,shiftz_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,SYM,SYM,ANT)
                 shiftz_new(i,j,k) = f_int * facRp * facRp
              end if
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
              call weights2(x(i,j,k),y(1,j,1),z(i,j,k), PhysicalRadius(i,j,k),RadiusDerivative(i,j,k),&
                   lapse(i,j,k),phi(i,j,k), &
                   dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, &
                   kmin,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, &
                   Xcenter, Ycenter, Zcenter, &
                   gxx(i,j,k),gxy(i,j,k),gxz(i,j,k), &
                   gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
              call interp(ex,shiftx_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,ANT,SYM,SYM)
              shiftx_new(i,j,k) = f_int * facRp * facRp
              if(Symmetry==AXISYM) then
                 call interp(ex,shifty_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,ANT,ANT,SYM)
              else
                 call interp(ex,shifty_old,f_int,il,iu,jl,ju,kl,ku, &
                      Dlx,Dly,Dlz,SYM,ANT,SYM)
              end if
              shifty_new(i,j,k) = f_int * facRp * facRp
              call interp(ex,shiftz_old,f_int,il,iu,jl,ju,kl,ku, &
                   Dlx,Dly,Dlz,SYM,SYM,ANT)
              shiftz_new(i,j,k) = f_int * facRp * facRp
           end do
        end do
     endif
     ! Transform the components of the shift back to fisheye coordinates
     if(fisheye_enable==1) then
        fish_to_phys = .FALSE.
        call trans_fish_phys_vector(ex,fake_ex, X, Y, Z,r_fish, PhysicalRadius, RadiusDerivative, &
             shiftx_old, shifty_old, shiftz_old, &
             fish_to_phys,2,have_bdry_min,have_bdry_max,use_trans_fish_phys_new)
        call trans_fish_phys_vector(ex,fake_ex, X, Y, Z,r_fish, PhysicalRadius, RadiusDerivative, &
             shiftx_new, shifty_new, shiftz_new, &
             fish_to_phys,0,have_bdry_min,have_bdry_max,use_trans_fish_phys_new)
     end if
     ! For axisymmetric runs, rotate the shift to the y .ne. 0 planes
     !  if (Symmetry==AXISYM) then
     !     call axibc_vector(ex,X,Y,Z,shiftx_new,shifty_new,shiftz_new)
     !  end if
  else if(shift_bc_type==FREEZE) then
     !
     ! x = xmax
     !
     if(have_bdry_max(1)==1) then
        i = imax
        do k = kmin, kmax - 1
           do j = jmin, jmax - 1
              shiftx_new(i,j,k) = shiftx_old(i,j,k)
              shifty_new(i,j,k) = shifty_old(i,j,k)
              shiftz_new(i,j,k) = shiftz_old(i,j,k)
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
              shiftx_new(i,j,k) = shiftx_old(i,j,k)
              shifty_new(i,j,k) = shifty_old(i,j,k)
              shiftz_new(i,j,k) = shiftz_old(i,j,k)
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
              shiftx_new(i,j,k) = shiftx_old(i,j,k)
              shifty_new(i,j,k) = shifty_old(i,j,k)
              shiftz_new(i,j,k) = shiftz_old(i,j,k)
           end do
        end do
     end if
     !
     ! x = xmin
     !
     if(have_bdry_min(1)==1) then
        shiftx_new(imin,:,:) = shiftx_old(imin,:,:)
        shifty_new(imin,:,:) = shifty_old(imin,:,:)
        shiftz_new(imin,:,:) = shiftz_old(imin,:,:)
     end if
     !
     ! y = ymin
     !
     if(have_bdry_min(2)==1) then
        shiftx_new(:,jmin,:) = shiftx_old(:,jmin,:)
        shifty_new(:,jmin,:) = shifty_old(:,jmin,:)
        shiftz_new(:,jmin,:) = shiftz_old(:,jmin,:)
     end if
     !
     ! z = zmin
     !
     if(have_bdry_min(3)==1) then
        shiftx_new(:,:,kmin) = shiftx_old(:,:,kmin)
        shifty_new(:,:,kmin) = shifty_old(:,:,kmin)
        shiftz_new(:,:,kmin) = shiftz_old(:,:,kmin)
     end if
  end if
end subroutine hyper_shift_bc
!-----------------------------------------------------------------------------
!
! Subtract off the radial component of the shift.
!
!-----------------------------------------------------------------------------
!
subroutine subtract_radial2(ex,X, Y, Z, PhysicalR,  &
     shiftx, shifty, shiftz, psi6_cent, psi6_init, Rscale)
  implicit none
  integer, dimension(3),intent(in)                :: ex
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: X,Y,Z,PhysicalR
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: shiftx, shifty, shiftz
  real*8                                          :: psi6_cent, psi6_init,Rscale
  !
  ! other variables
  !
  integer     :: imin, jmin, kmin, imax, jmax, kmax
  integer     :: i,j,k
  real*8      :: R,shiftR,M,fac
  !
  imin = lbound(shiftx,1)
  jmin = lbound(shiftx,2)
  kmin = lbound(shiftx,3)
  imax = ubound(shiftx,1)
  jmax = ubound(shiftx,2)
  kmax = ubound(shiftx,3)
  do k=kmin, kmax
     do j=jmin, jmax
        do i=imin, imax
           R = sqrt(X(i,1,1)*X(i,1,1) + Y(1,j,1)*Y(1,j,1) + Z(1,1,k)*Z(1,1,k))
           shiftR = (X(i,1,1)*shiftx(i,j,k) + Y(1,j,1)*shifty(i,j,k) &
                + Z(1,1,k)*shiftz(i,j,k))/R
           if(psi6_cent<2.d0*psi6_init) then
              fac = 0.D0
           else if(psi6_cent<4.d0*psi6_init) then
              fac = -1.D0 + 0.5D0*psi6_cent/psi6_init
           else
              fac = 1.D0
           end if
           fac = fac / (1.D0 + (PhysicalR(i,j,k)/Rscale)**4)
           shiftx(i,j,k) = shiftx(i,j,k) - shiftR * (X(i,1,1)/R) * fac
           shifty(i,j,k) = shifty(i,j,k) - shiftR * (Y(1,j,1)/R) * fac
           shiftz(i,j,k) = shiftz(i,j,k) - shiftR * (Z(1,1,k)/R) * fac
        end do
     end do
  end do
end subroutine subtract_radial2
!
! Transform the components of a vector
! between fisheye coordinates and physical coordinates
! near the boundary of the grid
!
subroutine trans_fish_phys_vector(ex,fake_ex, X, Y, Z, r_fish,PhysicalR, RadiusDerivative, &
     vx,vy,vz, &
     fish_to_phys,lb,have_bdry_min,have_bdry_max,use_trans_fish_phys_new)
  implicit none
  integer, dimension(3)                    :: ex,fake_ex,have_bdry_min,have_bdry_max
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z,r_fish,PhysicalR,RadiusDerivative
  real*8, dimension(ex(1),ex(2),ex(3))     :: vx,vy,vz
  integer                                  :: lb,use_trans_fish_phys_new
  !
  integer                                  :: i,j,k,imax,jmax,kmax
  integer                                  :: imin,jmin,kmin
  logical                                  :: fish_to_phys
  logical, dimension(ex(1),ex(2),ex(3))    :: tree
  !
  imin = ex(1)-fake_ex(1)+1
  jmin = ex(2)-fake_ex(2)+1
  kmin = ex(3)-fake_ex(3)+1
  imax = fake_ex(1)
  jmax = fake_ex(2)
  kmax = fake_ex(3)
  ! We dont want to do the transformation more than once at each grid point.
  ! Thats why its necessary to mark the grid point where the
  ! transformation has been carried out. Tree is the logical variable
  ! to make the mark.
  !
  tree = .FALSE.
  ! xmax
  if(have_bdry_max(1)==1) then
     do k=kmin,kmax
        do j=jmin,jmax
           do i=imax-lb,imax
              ! Recognize the tree?
              if (tree(i,j,k) .eqv. .FALSE.) then
                 if(use_trans_fish_phys_new==0) then
                    call trans_fish_phys_vj(ex,X,Y,Z,r_fish,PhysicalR,RadiusDerivative, &
                         vx,vy,vz,i,j,k,fish_to_phys)
                 else
                    call trans_fish_phys_vj_new(ex,X,Y,Z,r_fish,PhysicalR,RadiusDerivative, &
                         vx,vy,vz,i,j,k,fish_to_phys)
                 end if
                 tree(i,j,k) = .TRUE.
              end if
           end do
        end do
     end do
  end if
  ! ymax
  if(have_bdry_max(2)==1) then
     do k=kmin,kmax
        do j=jmax-lb,jmax
           do i=imin,imax
              ! Recognize the tree?
              if (tree(i,j,k) .eqv. .FALSE.) then
                 if(use_trans_fish_phys_new==0) then
                    call trans_fish_phys_vj(ex,X,Y,Z,r_fish,PhysicalR,RadiusDerivative, &
                         vx,vy,vz,i,j,k,fish_to_phys)
                 else
                    call trans_fish_phys_vj_new(ex,X,Y,Z,r_fish,PhysicalR,RadiusDerivative, &
                         vx,vy,vz,i,j,k,fish_to_phys)
                 end if
                 tree(i,j,k) = .TRUE.
              end if
           end do
        end do
     end do
  end if
  ! zmax
  if(have_bdry_max(3)==1) then
     do k=kmax-lb,kmax
        do j=jmin,jmax
           do i=imin,imax
              ! Recognize the tree?
              if (tree(i,j,k) .eqv. .FALSE.) then
                 if(use_trans_fish_phys_new==0) then
                    call trans_fish_phys_vj(ex,X,Y,Z,r_fish,PhysicalR,RadiusDerivative, &
                         vx,vy,vz,i,j,k,fish_to_phys)
                 else
                    call trans_fish_phys_vj_new(ex,X,Y,Z,r_fish,PhysicalR,RadiusDerivative, &
                         vx,vy,vz,i,j,k,fish_to_phys)
                 end if
                 tree(i,j,k) = .TRUE.
              end if
           end do
        end do
     end do
  end if
  ! xmin
  if(have_bdry_min(1)==1) then
     do k=kmin,kmax
        do j=jmin,jmax
           do i=imin,imin+lb
              ! Recognize the tree?
              if (tree(i,j,k) .eqv. .FALSE.) then
                 if(use_trans_fish_phys_new==0) then
                    call trans_fish_phys_vj(ex,X,Y,Z,r_fish,PhysicalR,RadiusDerivative, &
                         vx,vy,vz,i,j,k,fish_to_phys)
                 else
                    call trans_fish_phys_vj_new(ex,X,Y,Z,r_fish,PhysicalR,RadiusDerivative, &
                         vx,vy,vz,i,j,k,fish_to_phys)
                 end if
                 tree(i,j,k) = .TRUE.
              end if
           end do
        end do
     end do
  end if
  ! ymin 
  if(have_bdry_min(2)==1) then
     do k=kmin,kmax
        do j=jmin,jmin+lb
           do i=imin,imax
              ! Recognize the tree?
              if (tree(i,j,k) .eqv. .FALSE.) then
                 if(use_trans_fish_phys_new==0) then
                    call trans_fish_phys_vj(ex,X,Y,Z,r_fish,PhysicalR,RadiusDerivative, &
                         vx,vy,vz,i,j,k,fish_to_phys)
                 else
                    call trans_fish_phys_vj_new(ex,X,Y,Z,r_fish,PhysicalR,RadiusDerivative, &
                         vx,vy,vz,i,j,k,fish_to_phys)
                 end if
                 tree(i,j,k) = .TRUE.
              end if
           end do
        end do
     end do
  end if
  ! zmin
  if(have_bdry_min(3)==1) then
     do k=kmin,kmin+lb
        do j=jmin,jmax
           do i=imin,imax
              ! Recognize the tree?
              if (tree(i,j,k) .eqv. .FALSE.) then
                 if(use_trans_fish_phys_new==0) then
                    call trans_fish_phys_vj(ex,X,Y,Z,r_fish,PhysicalR,RadiusDerivative, &
                         vx,vy,vz,i,j,k,fish_to_phys)
                 else
                    call trans_fish_phys_vj_new(ex,X,Y,Z,r_fish,PhysicalR,RadiusDerivative, &
                         vx,vy,vz,i,j,k,fish_to_phys)
                 end if
                 tree(i,j,k) = .TRUE.
              end if
           end do
        end do
     end do
  end if
end subroutine trans_fish_phys_vector
subroutine trans_fish_phys_vj(ex, X, Y, Z, r_fish, PhysicalR, RadiusDerivative, &
     vx,vy,vz,i,j,k,fish_to_phys)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z,r_fish,PhysicalR,RadiusDerivative
  real*8, dimension(ex(1),ex(2),ex(3))     :: vx,vy,vz
  real*8                                   :: Rp,dr,r,xi,yi,zi,fac
  real*8                                   :: Jxx,Jxy,Jxz,Jyy,Jyz,Jzz
  real*8                                   :: vxi,vyi,vzi
  integer                                  :: i,j,k
  logical                                  :: fish_to_phys
  !
  xi = X(i,1,1)
  yi = Y(1,j,1)
  zi = Z(1,1,k)
  vxi = vx(i,j,k)
  vyi = vy(i,j,k)
  vzi = vz(i,j,k)
  r = r_fish(i,j,k)
  Rp = PhysicalR(i,j,k)
  dR = RadiusDerivative(i,j,k)
  ! transformation matrix
  if (fish_to_phys) then
     fac = dR - Rp/r
     Jxx = Rp/r + xi/r * xi/r * fac
     Jxy = xi/r * yi/r * fac
     Jxz = xi/r * zi/r * fac;
     Jyy = Rp/r + yi/r * yi/r * fac
     Jyz = yi/r * zi/r * fac
     Jzz = Rp/r + zi/r * zi/r * fac
  else
     fac = 1.d0/dR - r/Rp
     Jxx = r/Rp + xi/r * xi/r * fac
     Jxy = xi/r * yi/r * fac
     Jxz = xi/r * zi/r * fac;
     Jyy = r/Rp + yi/r * yi/r * fac
     Jyz = yi/r * zi/r * fac
     Jzz = r/Rp + zi/r * zi/r * fac
  end if
  vx(i,j,k) = Jxx*vxi + Jxy*vyi + Jxz*vzi
  vy(i,j,k) = Jxy*vxi + Jyy*vyi + Jyz*vzi
  vz(i,j,k) = Jxz*vxi + Jyz*vyi + Jzz*vzi
end subroutine trans_fish_phys_vj
subroutine trans_fish_phys_vj_new(ex, X, Y, Z, r_fish, PhysicalR, RadiusDerivative, &
     vx,vy,vz,i,j,k,fish_to_phys)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z,r_fish,PhysicalR,RadiusDerivative
  real*8, dimension(ex(1),ex(2),ex(3))     :: vx,vy,vz
  real*8                                   :: Rp,dr,r,xi,yi,zi,fac
  real*8                                   :: Jxx,Jxy,Jxz,Jyy,Jyz,Jzz
  real*8                                   :: vxi,vyi,vzi
  real*8                                   :: rpr,rinv,rpinv,drinv,nx,ny,nz
  integer                                  :: i,j,k
  logical                                  :: fish_to_phys
  !
  xi = X(i,1,1)
  yi = Y(1,j,1)
  zi = Z(1,1,k)
  vxi = vx(i,j,k)
  vyi = vy(i,j,k)
  vzi = vz(i,j,k)
  r = r_fish(i,j,k)
  Rp = PhysicalR(i,j,k)
  dR = RadiusDerivative(i,j,k)
  rinv=1.0d0/r
  nx=xi*rinv
  ny=yi*rinv
  nz=zi*rinv
  ! transformation matrix
  if (fish_to_phys) then
     rpr=Rp/r
     fac = dR - rpr
     Jxx = rpr + nx * nx * fac
     Jxy = nx * ny * fac
     Jxz = nx * nz * fac;
     Jyy = rpr + ny * ny * fac
     Jyz = ny * nz * fac
     Jzz = rpr + nz * nz * fac
  else
     rpr=r/Rp
     fac = 1.d0/dR - rpr
     Jxx = rpr + nx * nx * fac
     Jxy = nx * ny * fac
     Jxz = nx * nz * fac;
     Jyy = rpr + ny * ny * fac
     Jyz = ny * nz * fac
     Jzz = rpr + nz * nz * fac
  end if
  vx(i,j,k) = Jxx*vxi + Jxy*vyi + Jxz*vzi
  vy(i,j,k) = Jxy*vxi + Jyy*vyi + Jyz*vzi
  vz(i,j,k) = Jxz*vxi + Jyz*vyi + Jzz*vzi
end subroutine trans_fish_phys_vj_new
