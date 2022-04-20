!-----------------------------------------------------------------------------
!
! $Id: translate_matter.f90,v 1.11 2001/07/03 19:59:53 mduez Exp $
!
!-----------------------------------------------------------------------------
!
! Translate matter source terms
!
!-----------------------------------------------------------------------------
subroutine translate_matter(ex, Nx, Ny, Nz, ntot, &
     X, Y, Z, Xdest, Ydest, Zdest, &
     rho_old,S_old,rho_new,S_new, &
     Sxx_old, Sxx_new, &
     rho_b_old,u0_old,rho_b_new,u0_new, &
     lapse_old,phi_old,lapse_new,phi_new)
  implicit none
!
! Input parameters:
!
  integer, dimension(3)                    :: ex
  integer                                  :: Nx,Ny,Nz,ntot
  real*8, dimension(Nx)                    :: X
  real*8, dimension(Ny)                    :: Y
  real*8, dimension(Nz)                    :: Z
  real*8, dimension(ex(1)*ex(2)*ex(3))     :: Xdest,Ydest,Zdest
  real*8, dimension(ntot/4)                :: rho_old,S_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: rho_new,S_new
  real*8, dimension(ntot/4)                :: Sxx_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: Sxx_new
  real*8, dimension(ntot/4)                :: rho_b_old,u0_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: rho_b_new,u0_new
  real*8, dimension(ntot/4)                :: lapse_old,phi_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: lapse_new,phi_new
!
! Other variables:
!
  integer                    :: i, j, k, index
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  integer                    :: il,jl,kl,iu,ju,ku
  real*8                     :: Courant,dX,dY,dZ,Delx,Dely,Delz,f_int,fac,temp
  real*8                     :: HALF, ONE, TWO
  integer                    :: nx_adm, ny_adm
  parameter(HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0)
!
! where the action is...
!
  imin = lbound(rho_new,1)
  jmin = lbound(rho_new,2)
  kmin = lbound(rho_new,3)
  imax = ubound(rho_new,1)
  jmax = ubound(rho_new,2)
  kmax = ubound(rho_new,3)
  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)
  nx_adm = Nx
  ny_adm = Ny
!
!
  do i = imin, imax
     do j = jmin, jmax
        do k = kmin, kmax
           index = ex(1)*ex(2)*(k-1) + ex(1)*(j-1) + i
           call weight(Xdest(index),Ydest(index),Zdest(index), &
                X(1),Y(1),Z(1), &
                dX,dY,dZ,i,j,k,imin,imax,jmin,jmax,Nx/2,Ny/2, &
                il,iu,jl,ju,Delx,Dely)
!
! Interpolate 
!
           call interp4(ntot,Nx/2,Ny/2,nx_adm,ny_adm,lapse_old,f_int,il,iu,jl,ju,k,Delx,Dely,ONE,ONE)
           lapse_new(i,j,k) = f_int
           call interp4(ntot,Nx/2,Ny/2,nx_adm,ny_adm,phi_old,f_int,il,iu,jl,ju,k,Delx,Dely,ONE,ONE)
           phi_new(i,j,k) = f_int
           call interp4(ntot,Nx/2,Ny/2,nx_adm,ny_adm,rho_old,f_int,il,iu,jl,ju,k,Delx,Dely,ONE,ONE)
           rho_new(i,j,k) = f_int
           call interp4(ntot,Nx/2,Ny/2,nx_adm,ny_adm,S_old,f_int,il,iu,jl,ju,k,Delx,Dely,ONE,ONE)
           S_new(i,j,k) = f_int
           call interp4(ntot,Nx/2,Ny/2,nx_adm,ny_adm,Sxx_old,f_int,il,iu,jl,ju,k,Delx,Dely,ONE,ONE)
           Sxx_new(i,j,k) = f_int
           call interp4(ntot,Nx/2,Ny/2,nx_adm,ny_adm,rho_b_old,f_int,il,iu,jl,ju,k,Delx,Dely,ONE,ONE)
           rho_b_new(i,j,k) = f_int
           call interp4(ntot,Nx/2,Ny/2,nx_adm,ny_adm,u0_old,f_int,il,iu,jl,ju,k,Delx,Dely,ONE,ONE)
           u0_new(i,j,k) = f_int
        enddo
     enddo
  enddo
end subroutine translate_matter
!
!-----------------------------------------------------------------------------
!
! Interpolation for point i,j,k
!
!-----------------------------------------------------------------------------
!
subroutine weight(Xdest,Ydest,Zdest, &
     Xgmin,Ygmin,Zgmin, &
     dX,dY,dZ,i,j,k, &
     imin,imax,jmin,jmax,Nx,Ny, &
     il,iu,jl,ju,Delx,Dely)
  implicit none
  real*8, intent(in)   :: Xdest,Ydest,Zdest
  real*8, intent(in)   :: Xgmin,Ygmin,Zgmin,dX,dY,dZ
  integer, intent(in)  :: i,j,k
  integer, intent(in)  :: imin,imax,jmin,jmax
  integer, intent(in)  :: Nx,Ny
  integer, intent(out) :: il,iu,jl,ju
  real*8, intent(out)  :: Delx,Dely
!
  real*8               :: xl,xu,yl,yu,ONE,TWO
  parameter ( ONE = 1.D0, TWO = 2.D0 )
!
  il = floor((Xdest-Xgmin)/dX) + 1
  jl = floor((Ydest-Ygmin)/dY) + 1
  if (il < 2) then
     il = 2
  end if
  if (jl < 2) then
     jl = 2
  end if
  if (il > 2*Nx-1) then
     il = 2*Nx-1
  end if
  if (jl > 2*Ny-1) then
     jl = 2*Ny-1
  end if
!
  iu = il + 1
  ju = jl + 1
  xl = Xgmin + dX*(il-1)
  xu = xl + dX
  yl = Ygmin + dY*(jl - 1)
  yu = yl + dY
!
! get upper values
!
! find coefficients for interpolation
!
  Delx = (Xdest - xl)/(xu - xl)
  Dely = (Ydest - yl)/(yu - yl)
  return
end subroutine weight
!
!-----------------------------------------------------------------------------
!
! Interpolate function f using weights Delx and Dely
!
!-----------------------------------------------------------------------------
!!
subroutine interp4(ntot,Nx,Ny,nx_adm,ny_adm,f,f_int,il,iu,jl,ju,kl,Dx,Dy,SymX,SymY)
  implicit none
!
! Input:
!
  integer                                  :: ntot,Nx,Ny,nx_adm,ny_adm
  real*8, dimension(ntot)                  :: f
  real*8                                   :: f_int
  integer                                  :: il,iu,jl,ju,kl
  real*8                                   :: Dx,Dy,SymX,SymY
!
! Other variables:
!
  integer     :: iuc, ilc, juc, jlc, klc
  integer     :: iuuc, illc, juuc, jllc
  integer     :: ll_ll, ll_l, ll_u, ll_uu
  integer     :: l_ll, l_l, l_u, l_uu
  integer     :: u_ll, u_l, u_u, u_uu
  integer     :: uu_ll, uu_l, uu_u, uu_uu
  integer     :: Lx, Ly
  real*8      :: f_1x1y, f_2x1y, f_3x1y, f_4x1y
  real*8      :: f_1x2y, f_2x2y, f_3x2y, f_4x2y
  real*8      :: f_1x3y, f_2x3y, f_3x3y, f_4x3y
  real*8      :: f_1x4y, f_2x4y, f_3x4y, f_4x4y
  real*8      :: f_1x, f_2x, f_3x, f_4x
  real*8      :: f_int1, f_int2
  real*8      :: l,two,six,F3o4,F1o4
  real*8      :: sxu,syu,sxuu,syuu,sxl,syl,sxll,syll
  parameter ( l = 1.D0, two = 2.D0, six = 6.D0 )
  parameter ( F3o4 = 3.D0/4.D0, F1o4 = 1.D0/4.D0 )
!
! First find the appropriate indicies as they were assigned in C.
!
  iuc  = iu - 1
  iuuc = iuc + 1
  ilc  = il - 1
  illc = ilc - 1
  juc  = ju - 1
  juuc = juc + 1
  jlc  = jl - 1
  jllc = jlc - 1
  klc  = kl - 1
  Lx = Nx
  Ly = Ny
!
! Now we reflect onto the appropriate quadrant
!  
  if (iuc<Lx) then
     iuc = Lx - iuc - 1
     sxu = SymX
  else
     iuc = iuc - Lx
     sxu = l
  end if
  if (ilc<Lx) then
     ilc = Lx - ilc - 1
     sxl = SymX
  else
     ilc = ilc - Lx
     sxl = l
  end if
  if (illc<Lx) then
     illc = Lx - illc - 1
     sxll = SymX
  else
     illc = illc - Lx
     sxll = l
  end if
  if (iuuc<Lx) then
     iuuc = Lx - iuuc - 1
     sxuu = Symx
  else
     iuuc = iuuc - Lx
     sxuu = l
  end if
  if (juc<Ly) then
     juc = Ly - juc - 1
     syu = Symy
  else
     juc = juc - Ly
     syu = l
  end if
  if (jlc<Ly) then
     jlc = Ly - jlc - 1
     syl = Symy
  else
     jlc = jlc - Ly
     syl = l
  end if
  if (jllc<Ly) then
     jllc = Ly - jllc - 1
     syll = Symy
  else
     jllc = jllc - Ly
     syll = l
  end if
  if (juuc<Ly) then
     juuc = Ly - juuc - 1
     syuu = Symy
  else
     juuc = juuc - Ly
     syuu = l
  end if
!
! No point in interpolating in the vacuum
!
  if (f(Ny*Nx*klc + Nx*jlc + ilc + 1)==0.0 &
     .and.f(Ny*Nx*klc + Nx*jlc + iuc + 1)==0.0 &
     .and.f(Ny*Nx*klc + Nx*juc + ilc + 1)==0.0 &
     .and.f(Ny*Nx*klc + Nx*juc + iuc + 1)==0.0) then
     f_int = 0.D0
  else
!
!! Now for some 2D cubic interpolation
!
     ll_ll  = Ny*Nx*klc + Nx*(jllc) + (illc) + 1
     ll_l   = Ny*Nx*klc + Nx*jlc    + (illc) + 1
     ll_u   = Ny*Nx*klc + Nx*juc    + (illc) + 1
     ll_uu  = Ny*Nx*klc + Nx*(juuc) + (illc) + 1
     l_ll   = Ny*Nx*klc + Nx*(jllc) + (ilc)  + 1
     l_l    = Ny*Nx*klc + Nx*jlc    + (ilc)  + 1
     l_u    = Ny*Nx*klc + Nx*juc    + (ilc)  + 1
     l_uu   = Ny*Nx*klc + Nx*(juuc) + (ilc)  + 1
     u_ll   = Ny*Nx*klc + Nx*(jllc) + (iuc)  + 1
     u_l    = Ny*Nx*klc + Nx*jlc    + (iuc)  + 1
     u_u    = Ny*Nx*klc + Nx*juc    + (iuc)  + 1
     u_uu   = Ny*Nx*klc + Nx*(juuc) + (iuc)  + 1
     uu_ll  = Ny*Nx*klc + Nx*(jllc) + (iuuc) + 1
     uu_l   = Ny*Nx*klc + Nx*jlc    + (iuuc) + 1
     uu_u   = Ny*Nx*klc + Nx*juc    + (iuuc) + 1
     uu_uu  = Ny*Nx*klc + Nx*(juuc) + (iuuc) + 1
     f_1x = (Dy*(Dy-l)*(two-Dy)/six)*f(ll_ll)*sxll*syll + &
          ((l+Dy)*(l-Dy)*(two-Dy)/two)*f(ll_l)*sxll*syl + &
          (Dy*(l+Dy)*(two-Dy)/two)*f(ll_u)*sxll*syu + &
          (Dy*(l+Dy)*(Dy-l)/six)*f(ll_uu)*sxll*syuu
     f_2x = (Dy*(Dy-l)*(two-Dy)/six)*f(l_ll)*sxl*syll + &
          ((l+Dy)*(l-Dy)*(two-Dy)/two)*f(l_l)*sxl*syl + &
          (Dy*(l+Dy)*(two-Dy)/two)*f(l_u)*sxl*syu + &
          (Dy*(l+Dy)*(Dy-l)/six)*f(l_uu)*sxl*syuu
     f_3x = (Dy*(Dy-l)*(two-Dy)/six)*f(u_ll)*sxu*syll + &
          ((l+Dy)*(l-Dy)*(two-Dy)/two)*f(u_l)*sxu*syl + &
          (Dy*(l+Dy)*(two-Dy)/two)*f(u_u)*sxu*syu + &
          (Dy*(l+Dy)*(Dy-l)/six)*f(u_uu)*sxu*syuu
     f_4x = (Dy*(Dy-l)*(two-Dy)/six)*f(uu_ll)*sxuu*syll + &
          ((l+Dy)*(l-Dy)*(two-Dy)/two)*f(uu_l)*sxuu*syl + &
          (Dy*(l+Dy)*(two-Dy)/two)*f(uu_u)*sxuu*syu + &
          (Dy*(l+Dy)*(Dy-l)/six)*f(uu_uu)*sxuu*syuu
     f_int1 = (Dx*(Dx-l)*(two-Dx)/six)*f_1x + &
          ((l+Dx)*(l-Dx)*(two-Dx)/two)*f_2x + &
          (Dx*(l+Dx)*(two-Dx)/two)*f_3x + &
          (Dx*(l+Dx)*(Dx-l)/six)*f_4x
     f_int = f_int1
  end if
!
  return
end subroutine interp4
!
!
!-----------------------------------------------------------------------------
!
! Interpolate function f using weights Delx and Dely
!
!-----------------------------------------------------------------------------
!!
subroutine interp1(ntot,Nx,Ny,nx_adm,ny_adm,f,f_int,il,iu,jl,ju,kl,Dx,Dy,SymX,SymY)
  implicit none
!
! Input:
!
  integer                                  :: ntot,Nx,Ny,nx_adm,ny_adm
  real*8, dimension(ntot)                  :: f
  real*8                                   :: f_int
  integer                                  :: il,iu,jl,ju,kl
  real*8                                   :: Dx,Dy,SymX,SymY
!
! Other variables:
!
  integer     :: iuc, ilc, juc, jlc, klc
  integer     :: iuuc, illc, juuc, jllc
  integer     :: ll_ll, ll_l, ll_u, ll_uu
  integer     :: l_ll, l_l, l_u, l_uu
  integer     :: u_ll, u_l, u_u, u_uu
  integer     :: uu_ll, uu_l, uu_u, uu_uu
  integer     :: Lx, Ly
  real*8      :: f_1x1y, f_2x1y, f_3x1y, f_4x1y
  real*8      :: f_1x2y, f_2x2y, f_3x2y, f_4x2y
  real*8      :: f_1x3y, f_2x3y, f_3x3y, f_4x3y
  real*8      :: f_1x4y, f_2x4y, f_3x4y, f_4x4y
  real*8      :: f_1x, f_2x, f_3x, f_4x
  real*8      :: f_int1, f_int2
  real*8      :: l,two,six,F3o4,F1o4
  real*8      :: sxu,syu,sxuu,syuu,sxl,syl,sxll,syll
  parameter ( l = 1.D0, two = 2.D0, six = 6.D0 )
  parameter ( F3o4 = 3.D0/4.D0, F1o4 = 1.D0/4.D0 )
!
! First find the appropriate indicies as they were assigned in C.
!
  iuc  = iu - 1
  ilc  = il - 1
  juc  = ju - 1
  jlc  = jl - 1
  klc  = kl - 1
  Lx = Nx
  Ly = Ny
!
! Now we reflect onto the appropriate quadrant
!  
  if (iuc<Lx) then
     iuc = Lx - iuc - 1
     sxu = SymX
  else
     iuc = iuc - Lx
     sxu = l
  end if
  if (ilc<Lx) then
     ilc = Lx - ilc - 1
     sxl = SymX
  else
     ilc = ilc - Lx
     sxl = l
  end if
  if (juc<Ly) then
     juc = Ly - juc - 1
     syu = Symy
  else
     juc = juc - Ly
     syu = l
  end if
  if (jlc<Ly) then
     jlc = Ly - jlc - 1
     syl = Symy
  else
     jlc = jlc - Ly
     syl = l
  end if
!
! No point in interpolating in the vacuum
!
  if (f(Ny*Nx*klc + Nx*jlc + ilc + 1)==0.0 &
     .and.f(Ny*Nx*klc + Nx*jlc + iuc + 1)==0.0 &
     .and.f(Ny*Nx*klc + Nx*juc + ilc + 1)==0.0 &
     .and.f(Ny*Nx*klc + Nx*juc + iuc + 1)==0.0) then
     f_int = 0.D0
  else
!
!! Now for some 2D linear interpolation
!
     l_l    = Ny*Nx*klc + Nx*jlc + (ilc) + 1
     l_u    = Ny*Nx*klc + Nx*juc + (ilc) + 1
     u_l    = Ny*Nx*klc + Nx*jlc + (iuc) + 1
     u_u    = Ny*Nx*klc + Nx*juc + (iuc) + 1
     f_int  = (l - Dx) * ( (l - Dy)*f(l_l) + Dy * f(l_u) ) &
              +    Dx  * ( (1 - Dy)*f(u_l) + Dy * f(u_u) )
  end if
!
  return
end subroutine interp1
!
!-----------------------------------------------------------------------------
!
! Translate matter source terms for the binary neutron star problem
!
!-----------------------------------------------------------------------------
! Notes:  Nx,Ny,Nz are number of grid points on reference octant.
! nx_adm,ny_adm refer to the grid of the grid functions.  ntot_adm = nx_adm*ny_adm*nz_adm.
! X,Y,Z is the reference grid--expanded so that it takes up all 4 octants.
! Xdest,Ydest,Zdest are destination arrays---the gridfunction grid rotated in the xy plane.
! shiftx_old is a 1D array describing the shiftx data on the reference octant.  
! shiftx_new is a 3D array, the piece of the gridfunction accessible to this processor.
!-----------------------------------------------------------------------------
subroutine translate_bns(ex, Nx, Ny, Nz, nx_adm, ny_adm, ntot, X, Y, Z, Xdest, &
  Ydest, Zdest,shiftx_old,shifty_old,shiftz_old,shiftx_new, &
  shifty_new,shiftz_new,v2_old,v2_new,lapse_old,lapse_new, &
  phi_old,phi_new,lbx,lby,lbz)
  implicit none
!
! Input parameters:
!
  integer, dimension(3)                    :: ex
  integer                                  :: Nx,Ny,Nz,nx_adm,ny_adm,ntot
  real*8, dimension(2*Nx)                  :: X
  real*8, dimension(2*Ny)                  :: Y
  real*8, dimension(2*Nz)                  :: Z
  real*8, dimension(ex(1)*ex(2)*ex(3))     :: Xdest,Ydest,Zdest
  real*8, dimension(Nx*Ny*Nz)              :: shiftx_old,shifty_old,shiftz_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: shiftx_new,shifty_new,shiftz_new
  real*8, dimension(Nx*Ny*Nz)              :: v2_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: v2_new
  real*8, dimension(Nx*Ny*Nz)              :: lapse_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: lapse_new
  real*8, dimension(Nx*Ny*Nz)              :: phi_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: phi_new
  integer                                  :: lbx,lby,lbz
!
! Other variables:
!
  integer                    :: i, j, k, index
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  integer                    :: il,jl,kl,iu,ju,ku
  real*8                     :: Courant,dX,dY,dZ,Delx,Dely,Delz,f_int,fac
  real*8                     :: HALF, ONE, TWO
  real*8                     :: symx, symy
  parameter(HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0)
!
! where the action is...
!
  imin = lbound(shiftx_new,1)
  jmin = lbound(shiftx_new,2)
  kmin = lbound(shiftx_new,3)
  imax = ubound(shiftx_new,1)
  jmax = ubound(shiftx_new,2)
  kmax = ubound(shiftx_new,3)
  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)
!
!
  do i = imin, imax
     do j = jmin, jmax
        do k = kmin, kmax
           index = ex(1)*ex(2)*(k-1) + ex(1)*(j-1) + i
           call weight(Xdest(index),Ydest(index),Zdest(index), &
                X(1),Y(1),Z(1), &
                dX,dY,dZ,i,j,k,imin,imax,jmin,jmax,Nx,Ny, &
                il,iu,jl,ju,Delx,Dely)
           kl = k + lbz
!
! Interpolate 
!
           symx = ONE
           symy = -ONE
           call interp4(ntot,Nx,Ny,nx_adm,ny_adm,shiftx_old,f_int,il,iu,jl,ju,kl,Delx,Dely,symx,symy)
           shiftx_new(i,j,k) = f_int
           symx = -ONE
           symy = ONE
           call interp4(ntot,Nx,Ny,nx_adm,ny_adm,shifty_old,f_int,il,iu,jl,ju,kl,Delx,Dely,symx,symy)
           shifty_new(i,j,k) = f_int
           symx = -ONE
           symy = -ONE
           call interp4(ntot,Nx,Ny,nx_adm,ny_adm,shiftz_old,f_int,il,iu,jl,ju,kl,Delx,Dely,symx,symy)
           shiftz_new(i,j,k) = f_int
           call interp4(ntot,Nx,Ny,nx_adm,ny_adm,lapse_old,f_int,il,iu,jl,ju,kl,Delx,Dely,ONE,ONE)
           lapse_new(i,j,k) = f_int
           call interp4(ntot,Nx,Ny,nx_adm,ny_adm,phi_old,f_int,il,iu,jl,ju,kl,Delx,Dely,ONE,ONE)
           phi_new(i,j,k) = f_int
           call interp4(ntot,Nx,Ny,nx_adm,ny_adm,v2_old,f_int,il,iu,jl,ju,kl,Delx,Dely,ONE,ONE)
           v2_new(i,j,k) = f_int
        enddo
     enddo
  enddo
end subroutine translate_bns
!
!-----------------------------------------------------------------------------
!
! Rotate the Stress Tensor
!
!-----------------------------------------------------------------------------
subroutine rotate_tensor(ex,R,Sx,Sy,Sxx,Sxy,Sxz,Syy,Syz)
  implicit none
!
! Input parameters:
!
  integer, dimension(3)                    :: ex
  real*8, dimension(2,2)                   :: R
  real*8, dimension(ex(1),ex(2),ex(3))     :: Sx,Sy
  real*8, dimension(ex(1),ex(2),ex(3))     :: Sxx,Sxy,Sxz
  real*8, dimension(ex(1),ex(2),ex(3))     :: Syy,Syz
!
! Other variables:
!
  integer                    :: i, j, k, index
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  real*8                     :: xnew,ynew,xxnew,xynew,xznew,yynew,yznew
!
  imin = lbound(Sx,1)
  jmin = lbound(Sx,2)
  kmin = lbound(Sx,3)
  imax = ubound(Sx,1)
  jmax = ubound(Sx,2)
  kmax = ubound(Sx,3)
  do i = imin, imax
     do j = jmin, jmax
        do k = kmin, kmax
           xnew = R(1,1)*Sx(i,j,k) + R(2,1)*Sy(i,j,k)
           ynew = R(1,2)*Sx(i,j,k) + R(2,2)*Sy(i,j,k)
           Sx(i,j,k) = xnew
           Sy(i,j,k) = ynew
           xxnew = &
                R(1,1)*R(1,1)*Sxx(i,j,k) + R(1,1)*R(2,1)*Sxy(i,j,k) + &
                R(2,1)*R(1,1)*Sxy(i,j,k) + R(2,1)*R(2,1)*Syy(i,j,k)
           xynew = &
                R(1,1)*R(1,2)*Sxx(i,j,k) + R(1,1)*R(2,2)*Sxy(i,j,k) + &
                R(2,1)*R(1,2)*Sxy(i,j,k) + R(2,1)*R(2,2)*Syy(i,j,k)
           xznew = &
                R(1,1)*Sxz(i,j,k) + R(2,1)*Syz(i,j,k)
           yynew = &
                R(1,2)*R(1,2)*Sxx(i,j,k) + R(1,2)*R(2,2)*Sxy(i,j,k) + &
                R(2,2)*R(1,2)*Sxy(i,j,k) + R(2,2)*R(2,2)*Syy(i,j,k)
           yznew = &
                R(1,2)*Sxz(i,j,k) + R(2,2)*Syz(i,j,k)
           Sxx(i,j,k) = xxnew
           Sxy(i,j,k) = xynew
           Sxz(i,j,k) = xznew
           Syy(i,j,k) = yynew
           Syz(i,j,k) = yznew
        enddo
     enddo
  enddo
end subroutine rotate_tensor
!-----------------------------------------------------------------------------
!
! Rotate the Shift
!
!-----------------------------------------------------------------------------
subroutine rotate_shift(ex,R,shiftx,shifty)
  implicit none
!
! Input parameters:
!
  integer, dimension(3)                    :: ex
  real*8, dimension(4)                     :: R
  real*8, dimension(ex(1),ex(2),ex(3))     :: shiftx,shifty
!
! Other variables:
!
  integer                    :: i, j, k, index
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  real*8                     :: R11,R12,R21,R22,xnew,ynew
!
  R11 = R(1)
  R12 = R(2)
  R21 = R(3)
  R22 = R(4)
  imin = lbound(shiftx,1)
  jmin = lbound(shiftx,2)
  kmin = lbound(shiftx,3)
  imax = ubound(shiftx,1)
  jmax = ubound(shiftx,2)
  kmax = ubound(shiftx,3)
  do k = kmin, kmax
    do j = jmin, jmax
      do i = imin, imax
           xnew = R11*shiftx(i,j,k) + R12*shifty(i,j,k)
           ynew = R21*shiftx(i,j,k) + R22*shifty(i,j,k)
           shiftx(i,j,k) = xnew
           shifty(i,j,k) = ynew
        enddo
     enddo
  enddo
end subroutine rotate_shift
!-----------------------------------------------------------------------------
! Translate irrotational binary longitudinal terms.
!-----------------------------------------------------------------------------
subroutine Test()
  implicit none
!  integer                    :: ex
end subroutine Test
subroutine Translate_irr(ex, Nx, Ny, Nz, nx_adm, ny_adm, ntot, X, Y, Z, Xdest, &
  Ydest, Zdest,shiftx_old,shifty_old,shiftz_old,shiftx_new, &
  shifty_new,shiftz_new,v2_old,v2_new,lapse_old,lapse_new, &
  phi_old,phi_new,ux_old,uy_old,uz_old,ux_new,uy_new,uz_new,lbz)
  implicit none
!
! Input parameters:
!
  integer, dimension(3)                    :: ex
  integer                                  :: Nx,Ny,Nz,nx_adm,ny_adm,ntot,lbz
  real*8, dimension(2*Nx)                  :: X
  real*8, dimension(2*Ny)                  :: Y
  real*8, dimension(2*Nz)                  :: Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: Xdest,Ydest,Zdest
  !real*8, dimension(ex(1)*ex(2)*ex(3))     :: Xdest,Ydest,Zdest
  real*8, dimension(Nx*Ny*Nz)              :: shiftx_old,shifty_old,shiftz_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: shiftx_new,shifty_new,shiftz_new
  real*8, dimension(Nx*Ny*Nz)              :: v2_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: v2_new
  real*8, dimension(Nx*Ny*Nz)              :: lapse_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: lapse_new
  real*8, dimension(Nx*Ny*Nz)              :: phi_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: phi_new
  real*8, dimension(Nx*Ny*Nz)              :: ux_old,uy_old,uz_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: ux_new,uy_new,uz_new
!
! Other variables:
!
  integer                    :: i, j, k, index
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  integer                    :: il,jl,kl,iu,ju,ku
  real*8                     :: Courant,dX,dY,dZ,Delx,Dely,Delz,f_int,fac
  real*8                     :: HALF, ONE, TWO
  real*8                     :: symx, symy
  parameter(HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0)
!
! where the action is...
!
  imin = lbound(shiftx_new,1)
  jmin = lbound(shiftx_new,2)
  !this is just a guess, for EQUATORIAL symmetry!
  if(lbz .eq. -1) then
     kmin = lbound(shiftx_new,3)+1
  else
     kmin = lbound(shiftx_new,3)
  end if
  imax = ubound(shiftx_new,1)
  jmax = ubound(shiftx_new,2)
  kmax = ubound(shiftx_new,3)
  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)
!
!
  !this is just a guess, for EQUATORIAL symmetry!
  ! lbz = -1
!subroutine weight(Xdest,Ydest,Zdest, &
!     Xgmin,Ygmin,Zgmin, &
!     dX,dY,dZ,i,j,k, &
!     imin,imax,jmin,jmax,Nx,Ny, &
!     il,iu,jl,ju,Delx,Dely)
  v2_new = 0.D0
  shiftx_new = 0.D0
  shifty_new = 0.D0
  shiftz_new = 0.D0
  lapse_new = 0.D0
  phi_new = 0.D0
  ux_new = 0.D0
  uy_new = 0.D0
  uz_new = 0.D0
  do i = imin, imax
     do j = jmin, jmax
        do k = kmin, kmax
           !index = ex(1)*ex(2)*(k-1) + ex(1)*(j-1) + i
           call weight(Xdest(i,j,k),Ydest(i,j,k),Zdest(i,j,k), &
                X(1),Y(1),Z(1), &
                dX,dY,dZ,i,j,k,imin,imax,jmin,jmax,Nx,Ny, &
                il,iu,jl,ju,Delx,Dely)
           kl = k + lbz
!
! Interpolate 
!
           symx = ONE
           symy = -ONE
           call interp4(ntot,Nx,Ny,nx_adm,ny_adm,shiftx_old,f_int,il,iu,jl,ju,kl,Delx,Dely,symx,symy)
           shiftx_new(i,j,k) = f_int
           symx = -ONE
           symy = ONE
           call interp4(ntot,Nx,Ny,nx_adm,ny_adm,shifty_old,f_int,il,iu,jl,ju,kl,Delx,Dely,symx,symy)
           shifty_new(i,j,k) = f_int
           symx = -ONE
           symy = -ONE
           call interp4(ntot,Nx,Ny,nx_adm,ny_adm,shiftz_old,f_int,il,iu,jl,ju,kl,Delx,Dely,symx,symy)
           shiftz_new(i,j,k) = f_int
           call interp4(ntot,Nx,Ny,nx_adm,ny_adm,lapse_old,f_int,il,iu,jl,ju,kl,Delx,Dely,ONE,ONE)
           lapse_new(i,j,k) = f_int
           call interp4(ntot,Nx,Ny,nx_adm,ny_adm,phi_old,f_int,il,iu,jl,ju,kl,Delx,Dely,ONE,ONE)
           phi_new(i,j,k) = f_int
           call interp4(ntot,Nx,Ny,nx_adm,ny_adm,v2_old,f_int,il,iu,jl,ju,kl,Delx,Dely,ONE,ONE)
           v2_new(i,j,k) = f_int
           symx = ONE
           symy = -ONE
           call interp4(ntot,Nx,Ny,nx_adm,ny_adm,ux_old,f_int,il,iu,jl,ju,kl,Delx,Dely,symx,symy)
           ux_new(i,j,k) = f_int
           symx = -ONE
           symy = ONE
           call interp4(ntot,Nx,Ny,nx_adm,ny_adm,uy_old,f_int,il,iu,jl,ju,kl,Delx,Dely,symx,symy)
           uy_new(i,j,k) = f_int
           symx = -ONE
           symy = -ONE
           call interp4(ntot,Nx,Ny,nx_adm,ny_adm,uz_old,f_int,il,iu,jl,ju,kl,Delx,Dely,symx,symy)
           uz_new(i,j,k) = f_int
        enddo
     enddo
  enddo
end subroutine translate_irr
subroutine new_translate_irr(ext,nx,ny,nz,xmax,ymax,zmax, &
                dx,dy,dz, q_s, q, X,Y,Z,PhysicalR, symx,symy,symz, &
                m1,m2,m3,n,lbz)
  implicit none
  integer, dimension(3)                         :: ext
  integer                                       :: nx,ny,nz,m1,m2,m3,n
  real*8, dimension(nx*ny*nz)                   :: q_s
  real*8, dimension(ext(1),ext(2),ext(3))       :: q,PhysicalR,X,Y,Z
  real*8                                        :: xmax,ymax,zmax,dx,dy,dz
  integer                                       :: i,j,k, imin,imax, ia,ja,ka,lbz
  integer                                       :: jmin,jmax,kmin,kmax
  integer                                       :: ii,jj,kk,ind,ind1,ind2,ind3
  real*8                                        :: symx,symy,symz,xi,yi,zi,ri
  real*8                                        :: axi,ayi,azi,lamx,lamy,lamz
  real*8                                        :: xa,ya,za,ra
  integer, parameter                                :: m=4
  real*8, dimension(m)                                :: x1(m),x2(m),x3(m)
  real*8, dimension(m,m,m)                        :: yt
  real*8                                        :: symx1,symy1,symz1
!
  imin = lbound(q,1)
  jmin = lbound(q,2)
  kmin = lbound(q,3)
  imax = ubound(q,1)
  jmax = ubound(q,2)
  kmax = ubound(q,3)
  !this is just a guess, for EQUATORIAL symmetry!
  if(lbz .eq. -1) then
     kmin = lbound(q,3)+1
  else
     kmin = lbound(q,3)
  end if
  do k = kmin,kmax
     do j = jmin,jmax
        do i = imin,imax
           ri = sqrt(X(i,1,1)*X(i,1,1) + Y(1,j,1)*Y(1,j,1) + Z(1,1,k)*Z(1,1,k))
           !TEMPORARY, UNTIL FISHEYE IS SET UP!
           PhysicalR(i,j,k) = ri
           xi = X(i,1,1)/ri * PhysicalR(i,j,k)
           yi = Y(1,j,1)/ri * PhysicalR(i,j,k)
           zi = Z(1,1,k)/ri * PhysicalR(i,j,k)
           axi = abs(xi)
           ayi = abs(yi)
           azi = abs(zi)
           if (axi .ge. xmax .or. ayi .ge. ymax .or. azi .gt. zmax) then
              ! Assign the value according to the asymptotic behavior 
              ! q ~ x^m1 y^m2 z^m3 / r^n
              lamx = xmax/(X(i,1,1)/ri)
              lamy = ymax/(Y(1,j,1)/ri)
              lamz = zmax/(Z(1,1,k)/ri)
              if (lamx .le. lamy .and. lamx .le. lamz) then
                 ia = nx
                 ja = max(int(lamx*(Y(1,j,1)/ri)/dy) + 1, 2)
                 ka = max(int(lamx*(Z(1,1,k)/ri)/dz) + 1, 2)
              end if
              if (lamy .le. lamx .and. lamy .le. lamz) then
                 ia = max(int(lamy*(X(i,1,1)/ri)/dx) + 1, 2)
                 ja = ny
                 ka = max(int(lamy*(Z(1,1,k)/ri)/dz) + 1, 2)
              end if
              if (lamz .le. lamx .and. lamz .le. lamy) then
                    ia = max(int(lamz*(X(i,1,1)/ri)/dx) + 1, 2)
                 ja = max(int(lamz*(Y(1,j,1)/ri)/dy) + 1, 2)
                 ka = nz
              end if
              xa = dble(ia-1)*dx
              ya = dble(ja-1)*dy
              za = dble(ka-1)*dz
              ra = sqrt(xa*xa + ya*ya + za*za)
              ind = ia-1 + nx*(ja-1) + nx*ny*(ka-1) + 1
              q(i,j,k) = q_s(ind) * (axi/xa)**m1 * (ayi/ya)**m2 * &
                         (azi/za)**m3 * (ra/PhysicalR(i,j,k))**n
           else
              ! Point is inside the initial data grid, do 3D interpolation
              ia = int(axi/dx)-1
              ja = int(ayi/dy)-1
              ka = int(azi/dz)-1
              if (ia+4 .gt. nx) ia = nx-4
              if (ja+4 .gt. ny) ja = ny-4
              if (ka+4 .gt. nz) ka = nz-4
              symx1 = 1.d0
              symy1 = 1.d0
              symz1 = 1.d0
              if (ia .eq. -1) symx1 = symx
              if (ja .eq. -1) symy1 = symy
              if (ka .eq. -1) symz1 = symz
              do ii=1,4
                 x1(ii) = dble(ia+ii-1)*dx
                 x2(ii) = dble(ja+ii-1)*dy
                 x3(ii) = dble(ka+ii-1)*dz
              end do
              do kk=1,m
                 ind3 = iabs(ka+kk-1)
                 do jj=1,m
                    ind2 = iabs(ja+jj-1)
                    do ii=1,m
                       ind1 = iabs(ia+ii-1)
                       ind = ind1 + nx*ind2 + nx*ny*ind3 + 1
                       yt(ii,jj,kk) = q_s(ind)
                    end do
                 end do
              end do
              yt(1,:,:) = symx1*yt(1,:,:)
              yt(:,1,:) = symy1*yt(:,1,:)
              yt(:,:,1) = symz1*yt(:,:,1)
              call polin3(x1,x2,x3,yt,m,m,m,axi,ayi,azi,q(i,j,k))
           end if
           if (xi .lt. 0.d0) q(i,j,k) = symx*q(i,j,k)
           if (yi .lt. 0.d0) q(i,j,k) = symy*q(i,j,k)
           if (zi .lt. 0.d0) q(i,j,k) = symz*q(i,j,k)
       end do
     end do
  end do
end subroutine new_translate_irr
      SUBROUTINE polin3(x1a,x2a,x3a,ya,m1,m2,m3,x1,x2,x3,y)
      implicit none
      interface
        subroutine polint(xa,ya,n,x,y,dy)
         IMPLICIT NONE
         integer                    :: n
         real*8                     :: dy,x,y
         real*8, dimension(n)       :: xa,ya
        end subroutine polint
      end interface
      INTEGER :: m1,m2,m3
      REAL*8 :: dy,x1,x2,x3,y
      real*8 :: x1a(m1),x2a(m2),x3a(m3),ya(m1,m2,m3)
!     USES polint
      INTEGER :: i,j,k
      REAL*8 :: y1tmp(m2,m3),y2tmp(m1),y3tmp(m2),y4tmp(m3)
!
      do k=1,m3
         do j=1,m2
            do i=1,m1
               y2tmp(i) = ya(i,j,k)
            end do
            call polint(x1a,y2tmp,m1,x1,y1tmp(j,k),dy)
         end do
      end do
      do k=1,m3
         do j=1,m2
            y3tmp(j) = y1tmp(j,k)
         end do
         call polint(x2a,y3tmp,m2,x2,y4tmp(k),dy)
      end do
      call polint(x3a,y4tmp,m3,x3,y,dy)
      end SUBROUTINE polin3

!***************************************************************!
! 							        !
! Use Keisuke's initial data. Read data from files              !
!                                                               !
!***************************************************************!
subroutine read_KTid(ext, nx,ny,nz, dX,dY,dZ, xmin,ymin,zmin, X,Y,Z, &
		     alpham1, betax,betay,betaz, &
		     phi, Axx,Axy,Axz,Ayy,Ayz,Azz, &
	   	     rho_b, u_x,u_y,u_z, Symmetry)
  implicit none
  integer, dimension(3)                         :: ext
  real*8, dimension(ext(1),ext(2),ext(3))       :: X,Y,Z
  integer                                       :: nx,ny,nz, i,j,k, ii,jj
  integer                                       :: imin,imax,jmin,jmax,kmin,kmax
  integer                                       :: gi,gj,gk,ind,ind_read
  integer                                       :: Symmetry
  integer, parameter 			        :: EQUATORIAL = 1
  real*8					:: dX,dY,dZ, xmin,ymin,zmin,temp
  real*8                                        :: psim4
  real*8, dimension(ext(1),ext(2),ext(3))	:: alpham1, betax,betay,betaz
  real*8, dimension(ext(1),ext(2),ext(3))       :: phi, rho_b, u_x,u_y,u_z
  real*8, dimension(ext(1),ext(2),ext(3))       :: Axx,Axy,Axz,Ayy,Ayz,Azz
  logical					:: exit_do
!
  imin = lbound(phi,1)
  imax = ubound(phi,1)
  jmin = lbound(phi,2)
  jmax = ubound(phi,2)
  kmin = lbound(phi,3)
  kmax = ubound(phi,3)

! Open data files
  open(11,file='Psi.id',status='old')
  open(12,file='lapse.id',status='old')
  open(13,file='shiftx.id',status='old')
  open(14,file='shifty.id',status='old')
  open(15,file='shiftz.id',status='old')
  open(16,file='Kxx.id',status='old')
  open(17,file='Kxy.id',status='old')
  open(18,file='Kxz.id',status='old')
  open(19,file='Kyy.id',status='old')
  open(20,file='Kyz.id',status='old')
  open(21,file='Kzz.id',status='old')
  open(22,file='rho_b.id',status='old')
  open(23,file='u_x.id',status='old')
  open(24,file='u_y.id',status='old')
  open(25,file='u_z.id',status='old')

!  if (Symmetry==EQUATORIAL .and. Z(1,1,kmin) .lt. 0.d0) kmin = kmin + 1

  if (Symmetry==EQUATORIAL) then
     exit_do = .FALSE.
     do 
        if (Z(1,1,kmin) .lt. 0.d0) then 
	   kmin = kmin +1
	else 
	   exit_do = .TRUE.
	end if
	if (exit_do) exit
     end do
  else
    write(*,*) 'Symmetry not supported in read_KTid!'
  end if

  ind_read = 0
  do k=kmin,kmax
     gk = int((Z(1,1,k)-zmin)/dZ + 1.d-5) 
     do j=jmin,jmax
  	gj = int((Y(1,j,1)-ymin)/dY + 1.d-5)
	do i=imin,imax
	   gi = int((X(i,1,1)-xmin)/dX + 1.d-5) 
	   ind = gi + gj*nx + gk*nx*ny
	   if (ind .lt. ind_read) then
	      write(*,*) "ind < ind_read in read_ktid..."
	      write(*,*) "This is not supposed to happen. Something's wrong!"
	      stop
	   end if
	   if (ind .gt. ind_read) then
	      do ii = ind_read, ind-1
	         do jj=11,25
		    read(jj,*) 
		 end do
	      end do
	      ind_read = ind
	   end if 
	   read(11,*) temp
	   phi(i,j,k) = log(temp)
           psim4 = 1.d0/temp**4
           read(12,*) temp
     	   alpham1(i,j,k) = temp - 1.d0
	   read(13,*) temp
	   betax(i,j,k) = temp
	   read(14,*) temp
	   betay(i,j,k) = temp
	   read(15,*) temp
	   betaz(i,j,k) = temp
	   read(16,*) temp
	   Axx(i,j,k) = temp*psim4
	   read(17,*) temp
	   Axy(i,j,k) = temp*psim4
	   read(18,*) temp
	   Axz(i,j,k) = temp*psim4
	   read(19,*) temp
	   Ayy(i,j,k) = temp*psim4
	   read(20,*) temp
	   Ayz(i,j,k) = temp*psim4
	   read(21,*) temp
	   Azz(i,j,k) = temp*psim4
	   read(22,*) temp
	   rho_b(i,j,k) = temp
	   read(23,*) temp
	   u_x(i,j,k) = temp
           read(24,*) temp
           u_y(i,j,k) = temp
           read(25,*) temp
           u_z(i,j,k) = temp
	   ind_read = ind_read + 1
	end do
     end do
  end do

! Close data files
  do ii=11,25
     close(ii)
  end do

!baaad:  Use cartsymgn!
! Set data in the lower Z ghost zone
!  kmin = lbound(phi,3)
!  if (Symmetry==EQUATORIAL .and. Z(1,1,kmin) .lt. 0.d0) then 
!     alpham1(:,:,kmin) = alpham1(:,:,kmin+1)
!       betax(:,:,kmin) = betax(:,:,kmin+1)
!       betay(:,:,kmin) = betay(:,:,kmin+1)
!       betaz(:,:,kmin) = -betaz(:,:,kmin+1)
!         phi(:,:,kmin) = phi(:,:,kmin+1)
!	 Axx(:,:,kmin) = Axx(:,:,kmin+1)
!         Axy(:,:,kmin) = Axy(:,:,kmin+1)
!         Axz(:,:,kmin) = -Axz(:,:,kmin+1)
!         Ayy(:,:,kmin) = Ayy(:,:,kmin+1)
!         Ayz(:,:,kmin) = -Ayz(:,:,kmin+1)
!         Azz(:,:,kmin) = Azz(:,:,kmin+1)
!       rho_b(:,:,kmin) = rho_b(:,:,kmin+1)
!     	 u_x(:,:,kmin) = u_x(:,:,kmin)
!         u_y(:,:,kmin) = u_y(:,:,kmin)
!         u_z(:,:,kmin) = -u_z(:,:,kmin)
!  end if
end subroutine read_KTid
