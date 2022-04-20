!-----------------------------------------------------------------------------
!
! $Id: ks_initial_metric.f90,v 1.1.1.1 2006/02/23 17:48:41 zetienne Exp $
!
!-----------------------------------------------------------------------------
!---------------------------------------------------+
!
! Set up the mask function for mass and angular momentum integrations.
!
!---------------------------------------------------+
subroutine set_masks(ex,x,y,z,mskf,hcmskf,excision_radius,Symmetry)
  implicit none
!~~~~~> input variables
  integer, intent(in)  :: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: x,y,z
  real*8,  intent(in)  :: excision_radius
  integer, intent(in)  :: Symmetry
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: mskf,hcmskf
!~~~~~> local variables
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: curr
  real*8, parameter :: zeo=0.d0,one=1.d0
!~~~~~> Interface
  interface
    real*8 function get_radius(x,y,z,Symmetry)
      implicit none
      real*8,  intent(in) :: x,y,z
      integer, intent(in) :: Symmetry
    end function get_radius
  end interface
!~~~~~> Input translation
  imin = lbound(mskf,1)
  jmin = lbound(mskf,2)
  kmin = lbound(mskf,3)
  imax = ubound(mskf,1)
  jmax = ubound(mskf,2)
  kmax = ubound(mskf,3)
!~~~~~> assign the value to the mask function.
    mskf = one
  hcmskf = one
  do k = kmin,kmax,1
    do j = jmin,jmax,1
      do i = imin,imax,1
        curr = get_radius(x(i,1,1),y(1,j,1),z(1,1,k),Symmetry)
!        if(curr < 1.1*excision_radius)   mskf(i,j,k) = zeo
        if(curr < 2.d0*excision_radius)   mskf(i,j,k) = zeo
        if(curr < excision_radius)   hcmskf(i,j,k) = zeo
      end do
    end do
  end do
  return
end subroutine set_masks
!---------------------------------------------------+
!
! Mask the variables for mass and angular momentum volume integrations.
!
!---------------------------------------------------+
subroutine mask_integrands(ex,mskf,axx,axy,axz,ayy,ayz,azz,k,kx,ky,kz, &
     rho,rho_star,sx,sy,sz,trr, &
     gaxxx,gaxxy,gaxxz,gaxyy,gaxyz,gaxzz, &
     gayxx,gayxy,gayxz,gayyy,gayyz,gayzz, &
     gazxx,gazxy,gazxz,gazyy,gazyz,gazzz,myres  )
  implicit none
!~~~~~> input variables
  integer, intent(in)  :: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in)    :: mskf
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: axx,axy,axz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: ayy,ayz,azz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: k,  kx,ky,kz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: rho,sx,sy,sz,trr
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: gaxxx,gaxxy,gaxxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: gaxyy,gaxyz,gaxzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: gayxx,gayxy,gayxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: gayyy,gayyz,gayzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: gazxx,gazxy,gazxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: gazyy,gazyz,gazzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: rho_star,myres
!~~~~~> mask the variables
  axx   = axx   * mskf
  axy   = axy   * mskf
  axz   = axz   * mskf
  ayy   = ayy   * mskf
  ayz   = ayz   * mskf
  azz   = azz   * mskf
  k     = k     * mskf
  kx    = kx    * mskf
  ky    = ky    * mskf
  kz    = kz    * mskf
  rho   = rho   * mskf
  rho_star   = rho_star   * mskf
  sx    = sx    * mskf
  sy    = sy    * mskf
  sz    = sz    * mskf
  trr   = trr   * mskf
  gaxxx = gaxxx * mskf
  gaxxy = gaxxy * mskf
  gaxxz = gaxxz * mskf
  gaxyy = gaxyy * mskf
  gaxyz = gaxyz * mskf
  gaxzz = gaxzz * mskf
  gayxx = gayxx * mskf
  gayxy = gayxy * mskf
  gayxz = gayxz * mskf
  gayyy = gayyy * mskf
  gayyz = gayyz * mskf
  gayzz = gayzz * mskf
  gazxx = gazxx * mskf
  gazxy = gazxy * mskf
  gazxz = gazxz * mskf
  gazyy = gazyy * mskf
  gazyz = gazyz * mskf
  gazzz = gazzz * mskf
  return
end subroutine mask_integrands
!---------------------------------------------------+
!
! Mask the variables for mass and angular momentum volume integrations.
!
!---------------------------------------------------+
subroutine mask_integrands2(ext,mskf,ahmskf,axx,axy,axz,ayy,ayz,azz,k,kx,ky,kz, &
     rho,rho_star,sx,sy,sz,s,trr, &
     gaxxx,gaxxy,gaxxz,gaxyy,gaxyz,gaxzz, &
     gayxx,gayxy,gayxz,gayyy,gayyz,gayzz, &
     gazxx,gazxy,gazxz,gazyy,gazyz,gazzz,myres, &
     st_y, Bx, Bz, Ex, Ez)
  implicit none
!~~~~~> input variables
  integer, intent(in)  :: ext(1:3)
  real*8, dimension(ext(1),ext(2),ext(3)), intent(in)    :: mskf,ahmskf
  real*8, dimension(ext(1),ext(2),ext(3)), intent(inout) :: axx,axy,axz
  real*8, dimension(ext(1),ext(2),ext(3)), intent(inout) :: ayy,ayz,azz
  real*8, dimension(ext(1),ext(2),ext(3)), intent(inout) :: k,  kx,ky,kz
  real*8, dimension(ext(1),ext(2),ext(3)), intent(inout) :: rho,sx,sy,sz,trr,s
  real*8, dimension(ext(1),ext(2),ext(3)), intent(inout) :: gaxxx,gaxxy,gaxxz
  real*8, dimension(ext(1),ext(2),ext(3)), intent(inout) :: gaxyy,gaxyz,gaxzz
  real*8, dimension(ext(1),ext(2),ext(3)), intent(inout) :: gayxx,gayxy,gayxz
  real*8, dimension(ext(1),ext(2),ext(3)), intent(inout) :: gayyy,gayyz,gayzz
  real*8, dimension(ext(1),ext(2),ext(3)), intent(inout) :: gazxx,gazxy,gazxz
  real*8, dimension(ext(1),ext(2),ext(3)), intent(inout) :: gazyy,gazyz,gazzz
  real*8, dimension(ext(1),ext(2),ext(3)), intent(inout) :: rho_star,myres
  real*8, dimension(ext(1),ext(2),ext(3)), intent(inout) :: st_y,Bx,Bz,Ex,Ez
!~~~~~> mask the variables
  axx   = axx   * mskf
  axy   = axy   * mskf
  axz   = axz   * mskf
  ayy   = ayy   * mskf
  ayz   = ayz   * mskf
  azz   = azz   * mskf
  k     = k     * mskf
  kx    = kx    * mskf
  ky    = ky    * mskf
  kz    = kz    * mskf
  rho   = rho   * mskf
!  rho_star   = rho_star   * mskf
  rho_star   = rho_star   * ahmskf
  sx    = sx    * mskf
  sy    = sy    * mskf
  sz    = sz    * mskf
  s     = s     * mskf
  trr   = trr   * mskf
  gaxxx = gaxxx * mskf
  gaxxy = gaxxy * mskf
  gaxxz = gaxxz * mskf
  gaxyy = gaxyy * mskf
  gaxyz = gaxyz * mskf
  gaxzz = gaxzz * mskf
  gayxx = gayxx * mskf
  gayxy = gayxy * mskf
  gayxz = gayxz * mskf
  gayyy = gayyy * mskf
  gayyz = gayyz * mskf
  gayzz = gayzz * mskf
  gazxx = gazxx * mskf
  gazxy = gazxy * mskf
  gazxz = gazxz * mskf
  gazyy = gazyy * mskf
  gazyz = gazyz * mskf
  gazzz = gazzz * mskf
  st_y  = st_y  * mskf
  Bx    = Bx    * mskf
  Bz    = Bz    * mskf
  Ex    = Ex    * mskf
  Ez    = Ez    * mskf
  return
end subroutine mask_integrands2
!---------------------------------------------------+
!
! Mask the variables for mass and angular momentum volume integrations.
!
!---------------------------------------------------+
subroutine mask_integrands3(ext,ahmskf,st_y,Bx,By,Bz,Ex,Ey,Ez)
  implicit none
!~~~~~> input variables
  integer, intent(in)                                    :: ext(1:3)
  real*8, dimension(ext(1),ext(2),ext(3)), intent(in)    :: ahmskf
  real*8, dimension(ext(1),ext(2),ext(3)), intent(inout) :: st_y,Bx,By,Bz,Ex,Ey,Ez
!~~~~~> mask the variables
  st_y  = st_y  * ahmskf
  Bx    = Bx    * ahmskf
  By    = By    * ahmskf
  Bz    = Bz    * ahmskf
  Ex    = Ex    * ahmskf
  Ey    = Ey    * ahmskf
  Ez    = Ez    * ahmskf
  return
end subroutine mask_integrands3
!---------------------------------------------------+
!
! Mask a variable in the Hamiltonian constraint
!
!---------------------------------------------------+
subroutine hcmskact(ex,mskf,f)
  implicit none
!~~~~~> input variables
  integer, intent(in)  :: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in)    :: mskf
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: f
!~~~~~> mask the variable
  f = f * mskf
  return
end subroutine hcmskact
!---------------------------------------------------------------------
!
! get_radius
!
!---------------------------------------------------------------------
real*8 function get_radius(x,y,z,Symmetry)
  implicit none
  real*8,  intent(in) :: x,y,z
  integer, intent(in) :: Symmetry
  integer, parameter  :: AXISYM = 4
  if (Symmetry == AXISYM) then
     get_radius = sqrt(x*x + z*z)
  else
     get_radius = sqrt(x*x + y*y + z*z)
  end if
end function get_radius
