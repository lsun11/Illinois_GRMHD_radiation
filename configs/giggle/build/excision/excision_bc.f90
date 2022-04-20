!-----------------------------------------------------------------------------
!
! $Id: excision_bc.f90,v 1.1.1.1 2006/02/23 17:48:41 zetienne Exp $
!
!-----------------------------------------------------------------------------
!-------------------------------------------------------------------------+
!
! Use Alcubierres simple excision method
!
! the differences are that the excision zone in the codes is a sphere,
! not a cube. Therefore an interpolation is needed.
!
!-------------------------------------------------------------------------+
subroutine metric_excision_bc(ex,x,y,z,phi_rhs,trK_rhs,&
                     gxx_rhs,gxy_rhs,gxz_rhs,gyy_rhs,gyz_rhs,gzz_rhs,&
                     Axx_rhs,Axy_rhs,Axz_rhs,Ayy_rhs,Ayz_rhs,Azz_rhs,&
                     Gmx_rhs,Gmy_rhs,Gmz_rhs,&
                     symmetry,dm)
  implicit none
!~~~~~> Input arguments
  integer, intent(in) :: ex(1:3),symmetry
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::x,y,z
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::phi_rhs,trK_rhs
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::gxx_rhs,gxy_rhs,gxz_rhs
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::gyy_rhs,gyz_rhs,gzz_rhs
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::Axx_rhs,Axy_rhs,Axz_rhs
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::Ayy_rhs,Ayz_rhs,Azz_rhs
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::Gmx_rhs,Gmy_rhs,Gmz_rhs
  integer :: dm(1:ex(1),1:ex(2),1:ex(3))
!~~~~~> local variable
  real*8,dimension(ex(1),ex(2),ex(3))   ::jx,jy,jz
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k,ni,flax,flay,flaz
  integer :: ksym,jsym,isym
  real*8  :: sf(1:3),rr
  integer, parameter :: inEZ=0, onbd=1, outEZ=2
  real*8,  parameter :: zero=0.d0, one=1.d0, TWO=2.D0
  integer, parameter :: nosymmetry=0, equatorial=1, octant=2, pi_symm=3
  integer, parameter :: axisym = 4
  real*8             :: dx, dy, dz, PI, fac
!~~~~~~> interface
!  interface
!   subroutine find_excision_zone(ex,x,y,z,dm,innerbd,Symmetry)
!    implicit none
!    integer, intent(in)  :: ex(1:3)
!    real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::x,y,z
!    integer, intent(out) :: dm(1:ex(1),1:ex(2),1:ex(3))
!    real*8               :: innerbd
!    integer, intent(in)  :: Symmetry
!   end subroutine find_excision_zone
!
!  end interface
!~~~~~~> the computational domain
  PI = acos(-one)
  imin = lbound(phi_rhs,1)
  jmin = lbound(phi_rhs,2)
  kmin = lbound(phi_rhs,3)
  imax = ubound(phi_rhs,1)
  jmax = ubound(phi_rhs,2)
  kmax = ubound(phi_rhs,3)
  dx = x(imin+1,1,1) - x(imin,1,1)
  dy = y(1,jmin+1,1) - y(1,jmin,1)
  dz = z(1,1,kmin+1) - z(1,1,kmin)
  if(Symmetry == OCTANT) then
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
     if(X(1,1,1).lt.0.D0) then
        imin = imin + 1
     end if
     if(Z(1,1,1).lt.0.D0) then
        kmin = kmin + 1
     end if
  else if(Symmetry == EQUATORIAL) then
     if(Z(1,1,1).lt.0.D0) then
        kmin = kmin + 1
     end if
  end if
!~~~~~~> Find the grid points on the boundary
!  call find_excision_zone(ex,x,y,z,dm,innerbd,symmetry)
  write(*,*) "excision_bc.f90:",Axx_rhs(10,2,2)
!~~~~~~> inner boundary condition
  if( symmetry == nosymmetry ) then
     isym = imin + 1
     jsym = jmin + 1
     ksym = kmin + 1
  else if( symmetry == equatorial ) then
     isym = imin + 1
     jsym = jmin + 1
     ksym = kmin
  else if( symmetry == octant ) then
     isym = imin
     jsym = jmin
     ksym = kmin
  else if( symmetry == pi_symm ) then
     isym = imin + 1
     jsym = jmin
     ksym = kmin
  endif
  if(Symmetry==axisym) then
     do k = kmin, kmax-1, 1
        do j = jmin, jmax, 1
           do i = imin, imax-1, 1
              if( dm(i,j,k) == onbd ) then
                 ni = 2
                 flax = 1
                 flaz = 1
                 if(x(i,1,1) < zero) flax = -1
                 if(z(1,1,k) < zero) flaz = -1
                 sf(1) = x(i,1,1) * flax
                 sf(3) = z(1,1,k) * flaz
                 if( dm(i+flax,j     ,k     ) == onbd ) sf(1) = zero
                 if( dm(i     ,j     ,k+flaz) == onbd ) sf(3) = zero
                 if(dm(i+flax,j,k) == outEZ) sf(3) = zero
                 if((dm(i+flax,j,k) == onbd).and.(dm(i,j,k+flaz) == onbd)) then
                    ni = 1
                    if( dm(i+flax,j,k+flaz) /= outEZ ) then
                         write(*,*) 'something wrong in find_excision_zone'
                         stop
                          end if
                 endif
                 if( ni == 2 ) then
                    rr = sf(1) + sf(3)
                    sf(1) = sf(1) / rr
                    sf(3) = sf(3) / rr
                    phi_rhs(i,j,k) = &
                         sf(1) * phi_rhs(i+flax,j,k) + &
                         sf(3) * phi_rhs(i,j,k+flaz)
                    trK_rhs(i,j,k) = &
                         sf(1) * trK_rhs(i+flax,j,k) + &
                         sf(3) * trK_rhs(i,j,k+flaz)
                    gxx_rhs(i,j,k) = &
                         sf(1) * gxx_rhs(i+flax,j,k) + &
                         sf(3) * gxx_rhs(i,j,k+flaz)
                    gxy_rhs(i,j,k) = &
                         sf(1) * gxy_rhs(i+flax,j,k) + &
                         sf(3) * gxy_rhs(i,j,k+flaz)
                    gxz_rhs(i,j,k) = &
                         sf(1) * gxz_rhs(i+flax,j,k) + &
                         sf(3) * gxz_rhs(i,j,k+flaz)
                    gyy_rhs(i,j,k) = &
                         sf(1) * gyy_rhs(i+flax,j,k) + &
                         sf(3) * gyy_rhs(i,j,k+flaz)
                    gyz_rhs(i,j,k) = &
                         sf(1) * gyz_rhs(i+flax,j,k) + &
                         sf(3) * gyz_rhs(i,j,k+flaz)
                    gzz_rhs(i,j,k) = &
                         sf(1) * gzz_rhs(i+flax,j,k) + &
                         sf(3) * gzz_rhs(i,j,k+flaz)
                    Axx_rhs(i,j,k) = &
                         sf(1) * Axx_rhs(i+flax,j,k) + &
                         sf(3) * Axx_rhs(i,j,k+flaz)
                    Axy_rhs(i,j,k) = &
                         sf(1) * Axy_rhs(i+flax,j,k) + &
                         sf(3) * Axy_rhs(i,j,k+flaz)
                    Axz_rhs(i,j,k) = &
                         sf(1) * Axz_rhs(i+flax,j,k) + &
                         sf(3) * Axz_rhs(i,j,k+flaz)
                    Ayy_rhs(i,j,k) = &
                         sf(1) * Ayy_rhs(i+flax,j,k) + &
                         sf(3) * Ayy_rhs(i,j,k+flaz)
                    Ayz_rhs(i,j,k) = &
                         sf(1) * Ayz_rhs(i+flax,j,k) + &
                         sf(3) * Ayz_rhs(i,j,k+flaz)
                    Azz_rhs(i,j,k) = &
                         sf(1) * Azz_rhs(i+flax,j,k) + &
                         sf(3) * Azz_rhs(i,j,k+flaz)
                    Gmx_rhs(i,j,k) = &
                         sf(1) * Gmx_rhs(i+flax,j,k) + &
                         sf(3) * Gmx_rhs(i,j,k+flaz)
                    Gmy_rhs(i,j,k) = &
                         sf(1) * Gmy_rhs(i+flax,j,k) + &
                         sf(3) * Gmy_rhs(i,j,k+flaz)
                    Gmz_rhs(i,j,k) = &
                         sf(1) * Gmz_rhs(i+flax,j,k) + &
                         sf(3) * Gmz_rhs(i,j,k+flaz)
                 else if( ni == 1 ) then
                    phi_rhs(i,j,k) = phi_rhs(i+flax,j,k+flaz)
                    trK_rhs(i,j,k) = trK_rhs(i+flax,j,k+flaz)
                    gxx_rhs(i,j,k) = gxx_rhs(i+flax,j,k+flaz)
                    gxy_rhs(i,j,k) = gxy_rhs(i+flax,j,k+flaz)
                    gxz_rhs(i,j,k) = gxz_rhs(i+flax,j,k+flaz)
                    gyy_rhs(i,j,k) = gyy_rhs(i+flax,j,k+flaz)
                    gyz_rhs(i,j,k) = gyz_rhs(i+flax,j,k+flaz)
                    gzz_rhs(i,j,k) = gzz_rhs(i+flax,j,k+flaz)
                    Axx_rhs(i,j,k) = Axx_rhs(i+flax,j,k+flaz)
                    Axy_rhs(i,j,k) = Axy_rhs(i+flax,j,k+flaz)
                    Axz_rhs(i,j,k) = Axz_rhs(i+flax,j,k+flaz)
                    Ayy_rhs(i,j,k) = Ayy_rhs(i+flax,j,k+flaz)
                    Ayz_rhs(i,j,k) = Ayz_rhs(i+flax,j,k+flaz)
                    Azz_rhs(i,j,k) = Azz_rhs(i+flax,j,k+flaz)
                    Gmx_rhs(i,j,k) = Gmx_rhs(i+flax,j,k+flaz)
                    Gmy_rhs(i,j,k) = Gmy_rhs(i+flax,j,k+flaz)
                    Gmz_rhs(i,j,k) = Gmz_rhs(i+flax,j,k+flaz)
                 endif
              else if( dm(i,j,k) == inEZ ) then
                 phi_rhs(i,j,k) = zero
                 trK_rhs(i,j,k) = zero
                 gxx_rhs(i,j,k) = zero
                 gxy_rhs(i,j,k) = zero
                 gxz_rhs(i,j,k) = zero
                 gyy_rhs(i,j,k) = zero
                 gyz_rhs(i,j,k) = zero
                 gzz_rhs(i,j,k) = zero
                 Axx_rhs(i,j,k) = zero
                 Axy_rhs(i,j,k) = zero
                 Axz_rhs(i,j,k) = zero
                 Ayy_rhs(i,j,k) = zero
                 Ayz_rhs(i,j,k) = zero
                 Azz_rhs(i,j,k) = zero
                 Gmx_rhs(i,j,k) = zero
                 Gmy_rhs(i,j,k) = zero
                 Gmz_rhs(i,j,k) = zero
              endif
           enddo
        enddo
     end do
  else
     do k = ksym, kmax-1, 1
        do j = jsym, jmax-1, 1
           do i = isym, imax-1, 1
              if( dm(i,j,k) == onbd ) then
                 ni = 2
                 flax = 1
                 flay = 1
                 flaz = 1
                 if(x(i,1,1) < zero) flax = -1
                 if(y(1,j,1) < zero) flay = -1
                 if(z(1,1,k) < zero) flaz = -1
                 sf(1) = x(i,1,1) * flax
                 sf(2) = y(1,j,1) * flay
                 sf(3) = z(1,1,k) * flaz
                 if( dm(i+flax,j     ,k     ) == onbd ) sf(1) = zero
                 if( dm(i     ,j+flay,k     ) == onbd ) sf(2) = zero
                 if( dm(i     ,j     ,k+flaz) == onbd ) sf(3) = zero
                 if((dm(i+flax,j,k) == outEZ) .or. (dm(i,j+flay,k) == outEZ)) sf(3) = zero
                 if(( dm(i+flax,j     ,k     ) == onbd ) .and. &
                      ( dm(i     ,j+flay,k     ) == onbd ) .and. &
                      ( dm(i     ,j     ,k+flaz) == onbd )) then
                    ni = 1
                    if( dm(i     ,j+flay,k+flaz) == outEZ ) sf(2) = one
                    if( dm(i+flax,j     ,k+flaz) == outEZ ) sf(3) = one
!!$       if( dm(i+flax,j+flay,k     ) == outEZ ) sf(1) = one
                    if( dm(i+flax,j+flay,k     ) == outEZ ) then
                       sf(1) = one
                       sf(2) = zero
                       sf(3) = zero
                    endif
                    if(( dm(i+flax,j+flay,k) /= outEZ ) .and. &
                         ( dm(i,j+flay,k+flaz) /= outEZ ) .and. &
                         ( dm(i+flax,j,k+flaz) /= outEZ )) then
                         write(*,*) 'something wrong in find_excision_zone'
                         stop
                    end if
                 endif
                 if( ni == 2 ) then
                    rr = sf(1) + sf(2) + sf(3)
                    sf(1) = sf(1) / rr
                    sf(2) = sf(2) / rr
                    sf(3) = sf(3) / rr
                    phi_rhs(i,j,k) = sf(1) * phi_rhs(i+flax,j,k) + &
                         sf(2) * phi_rhs(i,j+flay,k) + &
                         sf(3) * phi_rhs(i,j,k+flaz)
                    trK_rhs(i,j,k) = sf(1) * trK_rhs(i+flax,j,k) + &
                         sf(2) * trK_rhs(i,j+flay,k) + &
                         sf(3) * trK_rhs(i,j,k+flaz)
                    gxx_rhs(i,j,k) = sf(1) * gxx_rhs(i+flax,j,k) + &
                         sf(2) * gxx_rhs(i,j+flay,k) + &
                         sf(3) * gxx_rhs(i,j,k+flaz)
                    gxy_rhs(i,j,k) = sf(1) * gxy_rhs(i+flax,j,k) + &
                         sf(2) * gxy_rhs(i,j+flay,k) + &
                         sf(3) * gxy_rhs(i,j,k+flaz)
                    gxz_rhs(i,j,k) = sf(1) * gxz_rhs(i+flax,j,k) + &
                         sf(2) * gxz_rhs(i,j+flay,k) + &
                         sf(3) * gxz_rhs(i,j,k+flaz)
                    gyy_rhs(i,j,k) = sf(1) * gyy_rhs(i+flax,j,k) + &
                         sf(2) * gyy_rhs(i,j+flay,k) + &
                         sf(3) * gyy_rhs(i,j,k+flaz)
                    gyz_rhs(i,j,k) = sf(1) * gyz_rhs(i+flax,j,k) + &
                         sf(2) * gyz_rhs(i,j+flay,k) + &
                         sf(3) * gyz_rhs(i,j,k+flaz)
                    gzz_rhs(i,j,k) = sf(1) * gzz_rhs(i+flax,j,k) + &
                         sf(2) * gzz_rhs(i,j+flay,k) + &
                         sf(3) * gzz_rhs(i,j,k+flaz)
                    Axx_rhs(i,j,k) = sf(1) * Axx_rhs(i+flax,j,k) + &
                         sf(2) * Axx_rhs(i,j+flay,k) + &
                         sf(3) * Axx_rhs(i,j,k+flaz)
                    Axy_rhs(i,j,k) = sf(1) * Axy_rhs(i+flax,j,k) + &
                         sf(2) * Axy_rhs(i,j+flay,k) + &
                         sf(3) * Axy_rhs(i,j,k+flaz)
                    Axz_rhs(i,j,k) = sf(1) * Axz_rhs(i+flax,j,k) + &
                         sf(2) * Axz_rhs(i,j+flay,k) + &
                         sf(3) * Axz_rhs(i,j,k+flaz)
                    Ayy_rhs(i,j,k) = sf(1) * Ayy_rhs(i+flax,j,k) + &
                         sf(2) * Ayy_rhs(i,j+flay,k) + &
                         sf(3) * Ayy_rhs(i,j,k+flaz)
                    Ayz_rhs(i,j,k) = sf(1) * Ayz_rhs(i+flax,j,k) + &
                         sf(2) * Ayz_rhs(i,j+flay,k) + &
                         sf(3) * Ayz_rhs(i,j,k+flaz)
                    Azz_rhs(i,j,k) = sf(1) * Azz_rhs(i+flax,j,k) + &
                         sf(2) * Azz_rhs(i,j+flay,k) + &
                         sf(3) * Azz_rhs(i,j,k+flaz)
                    Gmx_rhs(i,j,k) = sf(1) * Gmx_rhs(i+flax,j,k) + &
                         sf(2) * Gmx_rhs(i,j+flay,k) + &
                         sf(3) * Gmx_rhs(i,j,k+flaz)
                    Gmy_rhs(i,j,k) = sf(1) * Gmy_rhs(i+flax,j,k) + &
                         sf(2) * Gmy_rhs(i,j+flay,k) + &
                         sf(3) * Gmy_rhs(i,j,k+flaz)
                    Gmz_rhs(i,j,k) = sf(1) * Gmz_rhs(i+flax,j,k) + &
                         sf(2) * Gmz_rhs(i,j+flay,k) + &
                         sf(3) * Gmz_rhs(i,j,k+flaz)
                 else if( ni == 1 ) then
                    rr = sf(1) + sf(2) +sf(3)
                    sf(1) = sf(1) / rr
                    sf(2) = sf(2) / rr
                    sf(3) = sf(3) / rr
                    phi_rhs(i,j,k) = sf(1) * phi_rhs(i+flax,j+flay,k) + &
                         sf(2) * phi_rhs(i,j+flay,k+flaz) + &
                         sf(3) * phi_rhs(i+flax,j,k+flaz)
                    trK_rhs(i,j,k) = sf(1) * trK_rhs(i+flax,j+flay,k) + &
                         sf(2) * trK_rhs(i,j+flay,k+flaz) + &
                         sf(3) * trK_rhs(i+flax,j,k+flaz)
                    gxx_rhs(i,j,k) = sf(1) * gxx_rhs(i+flax,j+flay,k) + &
                         sf(2) * gxx_rhs(i,j+flay,k+flaz) + &
                         sf(3) * gxx_rhs(i+flax,j,k+flaz)
                    gxy_rhs(i,j,k) = sf(1) * gxy_rhs(i+flax,j+flay,k) + &
                         sf(2) * gxy_rhs(i,j+flay,k+flaz) + &
                         sf(3) * gxy_rhs(i+flax,j,k+flaz)
                    gxz_rhs(i,j,k) = sf(1) * gxz_rhs(i+flax,j+flay,k) + &
                         sf(2) * gxz_rhs(i,j+flay,k+flaz) + &
                         sf(3) * gxz_rhs(i+flax,j,k+flaz)
                    gyy_rhs(i,j,k) = sf(1) * gyy_rhs(i+flax,j+flay,k) + &
                         sf(2) * gyy_rhs(i,j+flay,k+flaz) + &
                         sf(3) * gyy_rhs(i+flax,j,k+flaz)
                    gyz_rhs(i,j,k) = sf(1) * gyz_rhs(i+flax,j+flay,k) + &
                         sf(2) * gyz_rhs(i,j+flay,k+flaz) + &
                         sf(3) * gyz_rhs(i+flax,j,k+flaz)
                    gzz_rhs(i,j,k) = sf(1) * gzz_rhs(i+flax,j+flay,k) + &
                         sf(2) * gzz_rhs(i,j+flay,k+flaz) + &
                         sf(3) * gzz_rhs(i+flax,j,k+flaz)
                    Axx_rhs(i,j,k) = sf(1) * Axx_rhs(i+flax,j+flay,k) + &
                         sf(2) * Axx_rhs(i,j+flay,k+flaz) + &
                         sf(3) * Axx_rhs(i+flax,j,k+flaz)
                    Axy_rhs(i,j,k) = sf(1) * Axy_rhs(i+flax,j+flay,k) + &
                         sf(2) * Axy_rhs(i,j+flay,k+flaz) + &
                         sf(3) * Axy_rhs(i+flax,j,k+flaz)
                    Axz_rhs(i,j,k) = sf(1) * Axz_rhs(i+flax,j+flay,k) + &
                         sf(2) * Axz_rhs(i,j+flay,k+flaz) + &
                         sf(3) * Axz_rhs(i+flax,j,k+flaz)
                    Ayy_rhs(i,j,k) = sf(1) * Ayy_rhs(i+flax,j+flay,k) + &
                         sf(2) * Ayy_rhs(i,j+flay,k+flaz) + &
                         sf(3) * Ayy_rhs(i+flax,j,k+flaz)
                    Ayz_rhs(i,j,k) = sf(1) * Ayz_rhs(i+flax,j+flay,k) + &
                         sf(2) * Ayz_rhs(i,j+flay,k+flaz) + &
                         sf(3) * Ayz_rhs(i+flax,j,k+flaz)
                    Azz_rhs(i,j,k) = sf(1) * Azz_rhs(i+flax,j+flay,k) + &
                         sf(2) * Azz_rhs(i,j+flay,k+flaz) + &
                         sf(3) * Azz_rhs(i+flax,j,k+flaz)
                    Gmx_rhs(i,j,k) = sf(1) * Gmx_rhs(i+flax,j+flay,k) + &
                         sf(2) * Gmx_rhs(i,j+flay,k+flaz) + &
                         sf(3) * Gmx_rhs(i+flax,j,k+flaz)
                    Gmy_rhs(i,j,k) = sf(1) * Gmy_rhs(i+flax,j+flay,k) + &
                         sf(2) * Gmy_rhs(i,j+flay,k+flaz) + &
                         sf(3) * Gmy_rhs(i+flax,j,k+flaz)
                    Gmz_rhs(i,j,k) = sf(1) * Gmz_rhs(i+flax,j+flay,k) + &
                         sf(2) * Gmz_rhs(i,j+flay,k+flaz) + &
                         sf(3) * Gmz_rhs(i+flax,j,k+flaz)
                 endif
              else if( dm(i,j,k) == inEZ ) then
                 phi_rhs(i,j,k) = zero
                 trK_rhs(i,j,k) = zero
                 gxx_rhs(i,j,k) = zero
                 gxy_rhs(i,j,k) = zero
                 gxz_rhs(i,j,k) = zero
                 gyy_rhs(i,j,k) = zero
                 gyz_rhs(i,j,k) = zero
                 gzz_rhs(i,j,k) = zero
                 Axx_rhs(i,j,k) = zero
                 Axy_rhs(i,j,k) = zero
                 Axz_rhs(i,j,k) = zero
                 Ayy_rhs(i,j,k) = zero
                 Ayz_rhs(i,j,k) = zero
                 Azz_rhs(i,j,k) = zero
                 Gmx_rhs(i,j,k) = zero
                 Gmy_rhs(i,j,k) = zero
                 Gmz_rhs(i,j,k) = zero
              endif
           enddo
        enddo
     enddo
  end if
  write(*,*) "excision_bc.f90:",Axx_rhs(10,2,2)
  return
end subroutine metric_excision_bc
!-------------------------------------------------------------------------+
!
! Use Alcubierres simple excision method for scalar
!
! excision region : sphere
!
!-------------------------------------------------------------------------+
subroutine scalar_excision_bc(ex,x,y,z,f,symmetry,dm)
  implicit none
!~~~~~~> Input arguments
  integer, intent(in) :: ex(1:3),symmetry
  integer :: dm(1:ex(1),1:ex(2),1:ex(3))
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::x,y,z
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::f
  real*8    :: innerbd
!~~~~~~> local variable
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k,flax,flay,flaz,ni
  integer :: ksym,jsym,isym,diax,diay,diaz
  real*8  :: sf(1:3),rr
  integer, parameter :: inEZ=0, onbd=1, outEZ=2
  real*8,  parameter :: zeo=0.d0, one=1.d0, two=2.d0
  integer, parameter :: nosymmetry=0, equatorial=1, octant=2
  integer, parameter :: pi_symm=3, axisym=4
!~~~~~~> interface
  interface
!     subroutine find_excision_zone(ex,x,y,z,dm,innerbd,Symmetry)
!       implicit none
!       integer, intent(in)  :: ex(1:3)
!       real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::x,y,z
!       integer, intent(out) :: dm(1:ex(1),1:ex(2),1:ex(3))
!       real*8               :: innerbd
!       integer, intent(in)  :: Symmetry
!     end subroutine find_excision_zone
  end interface
  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)
  imax = ubound(f,1)
  jmax = ubound(f,2)
  kmax = ubound(f,3)
!  if(Symmetry == OCTANT) then 
!     if(X(1,1,1).lt.0.D0) then
!        imin = imin + 1
!     end if
!     if(Y(1,1,1).lt.0.D0) then
!        jmin = jmin + 1
!     end if
!     if(Z(1,1,1).lt.0.D0) then
!        kmin = kmin + 1
!     end if
!  else if(Symmetry == AXISYM) then
!     if(X(1,1,1).lt.0.D0) then
!        imin = imin + 1
!     end if
!     if(Z(1,1,1).lt.0.D0) then
!        kmin = kmin + 1
!     end if
!  else if(Symmetry == EQUATORIAL) then
!     if(Z(1,1,1).lt.0.D0) then
!        kmin = kmin + 1
!     end if
!  end if
!~~~~~~> Find the grid points on the boundary
!  call find_excision_zone(ex,x,y,z,dm,innerbd,symmetry)
!~~~~~~> inner boundary condition
  if( symmetry == nosymmetry ) then
     isym = imin + 1
     jsym = jmin + 1
     ksym = kmin + 1
  else if( symmetry == equatorial ) then
     isym = imin + 1
     jsym = jmin + 1
     ksym = kmin
  else if( symmetry == octant ) then
     isym = imin
     jsym = jmin
     ksym = kmin
  else if( symmetry == pi_symm ) then
     isym = imin + 1
     jsym = jmin
     ksym = kmin
  endif
  if(symmetry .ne. axisym) then
     do k = ksym, kmax-1, 1
        do j = jsym, jmax-1, 1
           do i = isym, imax-1, 1
              if( dm(i,j,k) == onbd ) then
                 ni = 2
                 flax = 1
                 flay = 1
                 flaz = 1
                 if(x(i,1,1) < zeo) flax = -1
                 if(y(1,j,1) < zeo) flay = -1
                 if(z(1,1,k) < zeo) flaz = -1
                 sf(1) = x(i,1,1) * flax
                 sf(2) = y(1,j,1) * flay
                 sf(3) = z(1,1,k) * flaz
                 if( dm(i+flax,j     ,k     ) == onbd ) sf(1) = zeo
                 if( dm(i     ,j+flay,k     ) == onbd ) sf(2) = zeo
                 if( dm(i     ,j     ,k+flaz) == onbd ) sf(3) = zeo
                 if((dm(i+flax,j,k) == outEZ) .or. (dm(i,j+flay,k) == outEZ)) &
                      sf(3) = zeo
                 if(( dm(i+flax,j     ,k     ) == onbd ) .and. &
                      ( dm(i     ,j+flay,k     ) == onbd ) .and. &
                      ( dm(i     ,j     ,k+flaz) == onbd )) then
                    ni = 1
                    if( dm(i     ,j+flay,k+flaz) == outEZ ) sf(2) = one
                    if( dm(i+flax,j     ,k+flaz) == outEZ ) sf(3) = one
                    if( dm(i+flax,j+flay,k     ) == outEZ ) then
                       sf(1) = one
                       sf(2) = zeo
                       sf(3) = zeo
                    endif
                    if(( dm(i+flax,j+flay,k) /= outEZ ) .and. &
                         ( dm(i,j+flay,k+flaz) /= outEZ ) .and. &
                         ( dm(i+flax,j,k+flaz) /= outEZ )) then
                         write(*,*) 'something wrong in find_excision_zone'
                         stop
                    end if
                 endif
                 if( ni == 2 ) then
                    rr = sf(1) + sf(2) + sf(3)
                    sf(1) = sf(1) / rr
                    sf(2) = sf(2) / rr
                    sf(3) = sf(3) / rr
                    f(i,j,k) = sf(1) * f(i+flax,j,k) + &
                         sf(2) * f(i,j+flay,k) + &
                         sf(3) * f(i,j,k+flaz)
                 else if( ni == 1 ) then
                    rr = sf(1) + sf(2) +sf(3)
                    sf(1) = sf(1) / rr
                    sf(2) = sf(2) / rr
                    sf(3) = sf(3) / rr
                    f(i,j,k) = sf(1) * f(i+flax,j+flay,k) + &
                         sf(2) * f(i,j+flay,k+flaz) + &
                         sf(3) * f(i+flax,j,k+flaz)
                 endif
              else if( dm(i,j,k) == inEZ ) then
                 f(i,j,k) = zeo
              endif
           enddo
        enddo
     enddo
  else
     do k = kmin, kmax-1, 1
        do j = jmin, jmax, 1
           do i = imin, imax-1, 1
              if( dm(i,j,k) == onbd ) then
                 ni = 2
                 flax = 1
                 flaz = 1
                 if(x(i,1,1) < zeo) flax = -1
                 if(z(1,1,k) < zeo) flaz = -1
                   if (k==kmin) flaz = max(flaz,0)
                 sf(1) = x(i,1,1) * flax
                 sf(3) = z(1,1,k) * flaz
                 if( dm(i+flax,j     ,k     ) == onbd ) sf(1) = zeo
                 if( dm(i     ,j     ,k+flaz) == onbd ) sf(3) = zeo
                 if(dm(i+flax,j,k) == outEZ) sf(3) = zeo
                 if((dm(i+flax,j,k) == onbd).and.(dm(i,j,k+flaz) == onbd)) then
                    ni = 1
                    if( dm(i+flax,j,k+flaz) /= outEZ .and. flaz .ne. 0 ) then
                         write(*,*) 'something wrong in find_excision_zone'
                         stop
                    end if
                 endif
                 if( ni == 2 ) then
                    rr = sf(1) + sf(3)
                    sf(1) = sf(1) / rr
                    sf(3) = sf(3) / rr
                    f(i,j,k) = sf(1) * f(i+flax,j,k) + &
                         sf(3) * f(i,j,k+flaz)
                 else if( ni == 1 ) then
                    f(i,j,k) = f(i+flax,j,k+flaz)
                 endif
              else if( dm(i,j,k) == inEZ ) then
                 f(i,j,k) = zeo
              end if
           end do
        end do
     end do
  end if
  return
end subroutine scalar_excision_bc
!-------------------------------------------------------------------------+
!
! Use Alcubierres simple excision method for vector components
!                                             (or for 3 scalars)
! excision region : sphere
!
!-------------------------------------------------------------------------+
subroutine vector_excision_bc(ex,x,y,z,fx,fy,fz,symmetry,dm)
  implicit none
!~~~~~~> Input arguments
  integer, intent(in) :: ex(1:3),symmetry
  integer :: dm(1:ex(1),1:ex(2),1:ex(3))
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::x,y,z
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::fx,fy,fz
!~~~~~~> local variable
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k,flax,flay,flaz,ni
  integer :: ksym,jsym,isym,diax,diay,diaz
  real*8  :: sf(1:3),rr
  integer, parameter :: inEZ=0, onbd=1, outEZ=2
  real*8,  parameter :: zeo=0.d0, one=1.d0, two=2.d0
  integer, parameter :: nosymmetry=0, equatorial=1, octant=2
  integer, parameter :: pi_symm=3, axisym=4
!~~~~~~> interface
  interface
    ! subroutine find_excision_zone(ex,x,y,z,dm,innerbd,Symmetry)
    !   implicit none
    !   integer, intent(in)  :: ex(1:3)
    !   real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::x,y,z
    !   integer, intent(out) :: dm(1:ex(1),1:ex(2),1:ex(3))
    !   real*8               :: innerbd
    !   integer, intent(in)  :: Symmetry
    ! end subroutine find_excision_zone
  end interface
  imin = lbound(fx,1)
  jmin = lbound(fx,2)
  kmin = lbound(fx,3)
  imax = ubound(fx,1)
  jmax = ubound(fx,2)
  kmax = ubound(fx,3)
!~~~~~~> Find the grid points on the boundary
!  call find_excision_zone(ex,x,y,z,dm,innerbd,symmetry)
!  if(Symmetry == OCTANT) then 
!     if(X(1,1,1).lt.0.D0) then
!        imin = imin + 1
!     end if
!     if(Y(1,1,1).lt.0.D0) then
!        jmin = jmin + 1
!     end if
!     if(Z(1,1,1).lt.0.D0) then
!        kmin = kmin + 1
!     end if
!  else if(Symmetry == AXISYM) then
!     if(X(1,1,1).lt.0.D0) then
!        imin = imin + 1
!     end if
!     if(Z(1,1,1).lt.0.D0) then
!        kmin = kmin + 1
!     end if
!  else if(Symmetry == EQUATORIAL) then
!     if(Z(1,1,1).lt.0.D0) then
!        kmin = kmin + 1
!     end if
!  end if
!~~~~~~> inner boundary condition
  if( symmetry == nosymmetry ) then
     isym = imin + 1
     jsym = jmin + 1
     ksym = kmin + 1
  else if( symmetry == equatorial ) then
     isym = imin + 1
     jsym = jmin + 1
     ksym = kmin
  else if( symmetry == octant ) then
     isym = imin
     jsym = jmin
     ksym = kmin
  else if( symmetry == pi_symm ) then
     isym = imin + 1
     jsym = jmin
     ksym = kmin
  endif
  if(symmetry .ne. axisym) then
     do k = ksym, kmax-1, 1
        do j = jsym, jmax-1, 1
           do i = isym, imax-1, 1
              if( dm(i,j,k) == onbd ) then
                 ni = 2
                 flax = 1
                 flay = 1
                 flaz = 1
                 if(x(i,1,1) < zeo) flax = -1
                 if(y(1,j,1) < zeo) flay = -1
                 if(z(1,1,k) < zeo) flaz = -1
                 sf(1) = x(i,1,1) * flax
                 sf(2) = y(1,j,1) * flay
                 sf(3) = z(1,1,k) * flaz
                 if( dm(i+flax,j     ,k     ) == onbd ) sf(1) = zeo
                 if( dm(i     ,j+flay,k     ) == onbd ) sf(2) = zeo
                 if( dm(i     ,j     ,k+flaz) == onbd ) sf(3) = zeo
                 if((dm(i+flax,j,k) == outEZ) .or. (dm(i,j+flay,k) == outEZ)) &
                      sf(3) = zeo
                 if(( dm(i+flax,j     ,k     ) == onbd ) .and. &
                      ( dm(i     ,j+flay,k     ) == onbd ) .and. &
                      ( dm(i     ,j     ,k+flaz) == onbd )) then
                    ni = 1
                    if( dm(i     ,j+flay,k+flaz) == outEZ ) sf(2) = one
                    if( dm(i+flax,j     ,k+flaz) == outEZ ) sf(3) = one
                    if( dm(i+flax,j+flay,k     ) == outEZ ) then
                       sf(1) = one
                       sf(2) = zeo
                       sf(3) = zeo
                    endif
                    if(( dm(i+flax,j+flay,k) /= outEZ ) .and. &
                         ( dm(i,j+flay,k+flaz) /= outEZ ) .and. &
                         ( dm(i+flax,j,k+flaz) /= outEZ ) ) then
                         write(*,*) 'something wrong in find_excision_zone'
                         stop
                    end if
                 endif
                 if( ni == 2 ) then
                    rr = sf(1) + sf(2) + sf(3)
                    sf(1) = sf(1) / rr
                    sf(2) = sf(2) / rr
                    sf(3) = sf(3) / rr
                    fx(i,j,k) = sf(1) * fx(i+flax,j,k) + &
                         sf(2) * fx(i,j+flay,k) + &
                         sf(3) * fx(i,j,k+flaz)
                    if(i==2 .and. j==2 .and. k==11) then
                       write(*,*) "hi1", sf(1)
                       write(*,*) "hi2", fx(i+flax,j,k)
                       write(*,*) "hi3", sf(2)
                       write(*,*) "hi4", fx(i,j+flay,k)
                       write(*,*) "hi5", sf(3)
                       write(*,*) "hi6", fx(i,j,k+flaz)
                       write(*,*) "hi7", fx(i,j,k)
                    end if
                    fy(i,j,k) = sf(1) * fy(i+flax,j,k) + &
                         sf(2) * fy(i,j+flay,k) + &
                         sf(3) * fy(i,j,k+flaz)
                    fz(i,j,k) = sf(1) * fz(i+flax,j,k) + &
                         sf(2) * fz(i,j+flay,k) + &
                         sf(3) * fz(i,j,k+flaz)
                 else if( ni == 1 ) then
                    rr = sf(1) + sf(2) +sf(3)
                    sf(1) = sf(1) / rr
                    sf(2) = sf(2) / rr
                    sf(3) = sf(3) / rr
                    fx(i,j,k) = sf(1) * fx(i+flax,j+flay,k) + &
                         sf(2) * fx(i,j+flay,k+flaz) + &
                         sf(3) * fx(i+flax,j,k+flaz)
                    fy(i,j,k) = sf(1) * fy(i+flax,j+flay,k) + &
                         sf(2) * fy(i,j+flay,k+flaz) + &
                         sf(3) * fy(i+flax,j,k+flaz)
                    fz(i,j,k) = sf(1) * fz(i+flax,j+flay,k) + &
                         sf(2) * fz(i,j+flay,k+flaz) + &
                         sf(3) * fz(i+flax,j,k+flaz)
                 endif
              else if( dm(i,j,k) == inEZ ) then
                 fx(i,j,k) = zeo
                 fy(i,j,k) = zeo
                 fz(i,j,k) = zeo
              endif
           enddo
        enddo
     enddo
  else
     do k = kmin, kmax-1, 1
        do j = jmin, jmax, 1
           do i = imin, imax-1, 1
              if( dm(i,j,k) == onbd ) then
                 ni = 2
                 flax = 1
                 flaz = 1
                 if(x(i,1,1) < zeo) flax = -1
                 if(z(1,1,k) < zeo) flaz = -1
                 if (k==kmin) flaz = max(flaz,0)
                 sf(1) = x(i,1,1) * flax
                 sf(3) = z(1,1,k) * flaz
                 if( dm(i+flax,j     ,k     ) == onbd ) sf(1) = zeo
                 if( dm(i     ,j     ,k+flaz) == onbd ) sf(3) = zeo
                 if(dm(i+flax,j,k) == outEZ) sf(3) = zeo
                 if((dm(i+flax,j,k) == onbd).and.(dm(i,j,k+flaz) == onbd)) then
                    ni = 1
                    if( dm(i+flax,j,k+flaz) /= outEZ .and. flaz .ne. 0 ) then
                         write(*,*) 'something wrong in find_excision_zone'
                         stop
                    end if
                 endif
                 if( ni == 2 ) then
                    rr = sf(1) + sf(3)
                    sf(1) = sf(1) / rr
                    sf(3) = sf(3) / rr
                    fx(i,j,k) = sf(1) * fx(i+flax,j,k) + &
                         sf(3) * fx(i,j,k+flaz)
                    if(i==2 .and. j==2 .and. k==11) then
                       write(*,*) "hia1", sf(1)
                       write(*,*) "hia2", fx(i+flax,j,k)
                       write(*,*) "hia5", sf(3)
                       write(*,*) "hia6", fx(i,j,k+flaz)
                       write(*,*) "hia7", fx(i,j,k)
                    end if
                    fy(i,j,k) = sf(1) * fy(i+flax,j,k) + &
                         sf(3) * fy(i,j,k+flaz)
                    fz(i,j,k) = sf(1) * fz(i+flax,j,k) + &
                         sf(3) * fz(i,j,k+flaz)
                 else if( ni == 1 ) then
                    fx(i,j,k) = fx(i+flax,j,k+flaz)
                    if(i==2 .and. j==2 .and. k==11) then
                       write(*,*) "hib2", fx(i+flax,j,k+flaz)
                       write(*,*) "hib7", fx(i,j,k)
                    end if
                    fy(i,j,k) = fy(i+flax,j,k+flaz)
                    fz(i,j,k) = fz(i+flax,j,k+flaz)
                 endif
              else if( dm(i,j,k) == inEZ ) then
                 fx(i,j,k) = zeo
                 fy(i,j,k) = zeo
                 fz(i,j,k) = zeo
              end if
           end do
        end do
     end do
  end if
  return
end subroutine vector_excision_bc
!-----------------------------------------------------------------------------
!
! Remove matter inside and on the excision zone
!
!-----------------------------------------------------------------------------
subroutine remove_interior(ex,x,y,z,rho_star,dm,symmetry)
  implicit none
  integer, dimension(3)                               :: ex
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::x,y,z
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: rho_star
  integer :: dm(1:ex(1),1:ex(2),1:ex(3))
  integer :: excise_flag
  integer :: symmetry
! other variables
  integer              :: imin,jmin,kmin,imax,jmax,kmax
  integer              :: i, j, k
  integer :: ksym,jsym,isym,diax,diay,diaz
  real*8, parameter    :: ZERO  = 0.D0, ONE  = 1.D0, TWO = 2.D0, FOUR = 4.D0
  real*8  :: sf(1:3),rr
  integer, parameter :: inEZ=0, onbd=1, outEZ=2
  integer, parameter :: nosymmetry=0, equatorial=1, octant=2, pi_symm=3
  interface
!     subroutine find_excision_zone(ex,x,y,z,dm,innerbd,Symmetry)
!       implicit none
!       integer, intent(in)  :: ex(1:3)
!       real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::x,y,z
!       integer, intent(out) :: dm(1:ex(1),1:ex(2),1:ex(3))
!       real*8               :: innerbd
!       integer, intent(in)  :: Symmetry
!     end subroutine find_excision_zone
  end interface
  imin = lbound(rho_star,1)
  jmin = lbound(rho_star,2)
  kmin = lbound(rho_star,3)
  imax = ubound(rho_star,1)
  jmax = ubound(rho_star,2)
  kmax = ubound(rho_star,3)
!~~~~~~> Find the grid points on the boundary
!  call find_excision_zone(ex,x,y,z,dm,innerbd,symmetry)
  do k = kmin, kmax, 1
     do j = jmin, jmax, 1
        do i = imin, imax, 1
           if( dm(i,j,k)==inEZ.or.dm(i,j,k)==onbd ) then
              rho_star(i,j,k) = 0.d0
           end if
        end do
     end do
  end do
  return
end subroutine remove_interior
!-----------------------------------------------------------------------------
!
! Remove matter inside the excision zone
!
!-----------------------------------------------------------------------------
subroutine remove_interior2(ex,x,y,z,rho_star,dm,symmetry)
  implicit none
  integer, dimension(3)                               :: ex
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::x,y,z
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: rho_star
  integer :: excise_flag
  integer, intent(in) :: symmetry
  integer :: dm(1:ex(1),1:ex(2),1:ex(3))
! other variables
  integer              :: imin,jmin,kmin,imax,jmax,kmax
  integer              :: i, j, k
  real*8, parameter    :: ZERO  = 0.D0, ONE  = 1.D0, TWO = 2.D0, FOUR = 4.D0
  real*8  :: sf(1:3),rr
  integer, parameter :: inEZ=0, onbd=1, outEZ=2
  integer, parameter :: nosymmetry=0, equatorial=1, octant=2, pi_symm=3
  interface
     subroutine find_excision_zone(ex,x,y,z,dm,innerbd,Symmetry)
       implicit none
       integer, intent(in)  :: ex(1:3)
       real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::x,y,z
       integer, intent(out) :: dm(1:ex(1),1:ex(2),1:ex(3))
       real*8               :: innerbd
       integer, intent(in)  :: Symmetry
     end subroutine find_excision_zone
  end interface
  imin = lbound(rho_star,1)
  jmin = lbound(rho_star,2)
  kmin = lbound(rho_star,3)
  imax = ubound(rho_star,1)
  jmax = ubound(rho_star,2)
  kmax = ubound(rho_star,3)
!~~~~~~> Find the grid points on the boundary
!  call find_excision_zone(ex,x,y,z,dm,innerbd,symmetry)
  do k = kmin, kmax, 1
     do j = jmin, jmax, 1
        do i = imin, imax, 1
           if( dm(i,j,k)==inEZ ) then
              rho_star(i,j,k) = 0.d0
           end if
        end do
     end do
  end do
  return
end subroutine remove_interior2
!-------------------------------------------------------------------------+
!
! Find out the grid points on the boundary
!
!-------------------------------------------------------------------------+
  subroutine find_excision_zone(ex,x,y,z,dm,innerbd,Symmetry)
  implicit none
!~~~~~~> Input arguments
  integer, intent(in)  :: ex(1:3)
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::x,y,z
  integer, intent(out) :: dm(1:ex(1),1:ex(2),1:ex(3))
  real*8,  intent(in)  :: innerbd
  integer, intent(in)  :: Symmetry
!~~~~~~> Local variable
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k,m,n,zeq,xeq,yeq
  real*8  :: rr,rin
  integer, parameter :: onbd=1, inEZ=0, outEZ=2
  real*8   :: inbd
  real*8,  parameter :: outcut=300.d0, zeo=0.d0,cubd=1.001d0
!~~~~~~> interface
  interface
   real*8 function radius(x,y,z,Symmetry)
    implicit none
    real*8, intent(in) :: x,y,z
    integer, intent(in) :: Symmetry
   end function radius
  end interface
  inbd = innerbd * 0.8D0
!~~~~~~> Computational domain
  imin = lbound(x,1)
  jmin = lbound(x,2)
  kmin = lbound(x,3)
  imax = ubound(x,1)
  jmax = ubound(x,2)
  kmax = ubound(x,3)
!~~~~~~> Initialize dm
  dm = outEZ
!~~~~~~> Start finding....
  do k = kmin, kmax, 1
   do j = jmin, jmax, 1
    do i = imin, imax, 1
     rr = radius(x(i,1,1),y(1,j,1),z(1,1,k),Symmetry)
     if( rr < outcut ) then
      if( rr < inbd ) then
       dm(i,j,k) = inEZ
      else
       do m = -1,1,1
        do n = -1,1,2
         xeq = i + n
         yeq = j + m
         if( xeq < imin ) xeq = imin
         if( xeq > imax ) xeq = imax
         if( yeq < jmin ) yeq = jmin
         if( yeq > jmax ) yeq = jmax
         rin = radius(x(xeq,1,1),y(1,yeq,1),z(1,1,k),Symmetry)
         if( rin < inbd ) dm(i,j,k) = onbd
         yeq = j + n
         zeq = k + m
         if( yeq < jmin ) yeq = jmin
         if( yeq > jmax ) yeq = jmax
         if( zeq < kmin ) zeq = kmin
         if( zeq > kmax ) zeq = kmax
         rin = radius(x(i,1,1),y(1,yeq,1),z(1,1,zeq),Symmetry)
         if( rin < inbd ) dm(i,j,k) = onbd
         zeq = k + n
         xeq = i + m
         if( zeq < kmin ) zeq = kmin
         if( zeq > kmax ) zeq = kmax
         if( xeq < imin ) xeq = imin
         if( xeq > imax ) xeq = imax
         rin = radius(x(xeq,1,1),y(1,j,1),z(1,1,zeq),Symmetry)
         if( rin < inbd ) dm(i,j,k) = onbd
        enddo
       enddo
      endif
     endif
    enddo
   enddo
  enddo
  return
  end subroutine find_excision_zone
