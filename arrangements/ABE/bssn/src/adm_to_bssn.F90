!-----------------------------------------------------------------------------
!
! $Id: adm_to_bssn.f90,v 1.1.1.1 2006/02/23 17:48:41 zetienne Exp $
!
!-----------------------------------------------------------------------------
!
! convert physical to tilde metric, compute phi, and invert metric
!
!-----------------------------------------------------------------------------
subroutine convert(ext,phi,gxx,gxy,gxz,gyy,gyz,gzz,gupxx,gupxy,gupxz, &
     gupyy,gupyz,gupzz,det,detmin,detmax)
  implicit none
!
! Input parameters:
!
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: phi,gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ext(1),ext(2),ext(3))  :: gupxx,gupxy,gupxz
  real*8, dimension(ext(1),ext(2),ext(3))  :: gupyy,gupyz,gupzz
  real*8, dimension(ext(1),ext(2),ext(3))  :: det
  real*8                                   :: detmin, detmax
!
! Other variables: 
!
  real*8                                   :: F1o3, F1o12
  parameter( F1o3 = 1.D0 / 3.D0, F1o12 = 1.D0/12.D0  )

!
! compute and check determinant
!
!  write(*,*) "gijs: ", minval(gxx), minval(gxy), minval(gxz), minval(gyy), minval(gyz), minval(gzz)
  det =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz &
       - gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  detmin = minval(det)
  detmax = maxval(det)

  write(*,*) "maxmin determinants: ", detmin, detmax
! determine phi
!

  phi = F1o12 * log(det)

!
! divide metric by determinant
!
  gxx = gxx / det**F1o3
  gxy = gxy / det**F1o3
  gxz = gxz / det**F1o3
  gyy = gyy / det**F1o3
  gyz = gyz / det**F1o3
  gzz = gzz / det**F1o3

! Invert metric...
!
  gupxx =   ( gyy * gzz - gyz * gyz )! / det
  gupxy = - ( gxy * gzz - gyz * gxz )! / det
  gupxz =   ( gxy * gyz - gyy * gxz )! / det
  gupyy =   ( gxx * gzz - gxz * gxz )! / det
  gupyz = - ( gxx * gyz - gxy * gxz )! / det
  gupzz =   ( gxx * gyy - gxy * gxy )! / det

end subroutine convert

!-----------------------------------
! Compute Gamma^i using GenericFD.h
!-----------------------------------
#include "GenericFD.h"

subroutine setgamma_v2(ex, nghostzones, dX, dY, dZ, &
     phi, gxx, gxy, gxz, gyy, gyz, gzz, & 
     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Gamx, Gamy, Gamz, temp, &
     Symmetry)
  implicit none
  interface
     subroutine ddx(ext,f,df,dX,sym)
       implicit none
       integer, dimension(3)                    :: ext
       real*8, dimension(ext(1),ext(2),ext(3))  :: f,df
       real*8                                   :: sym,dX
     end subroutine ddx
     subroutine ddy(ext,f,df,dY,sym)
       implicit none
       integer, dimension(3)                    :: ext
       real*8, dimension(ext(1),ext(2),ext(3))  :: f,df
       real*8                                   :: sym,dY
     end subroutine ddy
     subroutine ddz(ext,f,df,dZ,sym)
       implicit none
       integer, dimension(3)                    :: ext
       real*8, dimension(ext(1),ext(2),ext(3))  :: f,df
       real*8                                   :: sym,dZ
     end subroutine ddz
  end interface
  integer, dimension(3)                 :: ex,nghostzones
  real*8                                :: dX,dY,dZ
  real*8, dimension(ex(1),ex(2),ex(3))  :: phi,gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ex(1),ex(2),ex(3))  :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3))  :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3))  :: Gamx, Gamy, Gamz
  real*8, dimension(ex(1),ex(2),ex(3))  :: temp
  integer                               :: Symmetry

!  real*8, dimension(ex(1),ex(2),ex(3))  :: temp
  integer                               :: istart, jstart, kstart, iend, jend, kend
  integer                               :: i,j,k
! 1st of 2 needed #includes for GenericFD.h:
#include "../../GenFD_decl_varF90.h"
  real*8                                :: F1o3, F1o12, SYM, ANTI
  integer                               :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter ( F1o3 = 1.D0/3.D0, F1o12 = 1.D0/12.D0 )
  parameter ( SYM = 1.D0, ANTI = - 1.D0 )
  parameter (NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

! Determine Gamx TO SECOND ORDER, first, since it fills the ghostzones
!
! Gamx:
!
  call ddx(ex,gupxx,temp,dX,SYM)
  Gamx = - temp
  call ddy(ex,gupxy,temp,dY,ANTI)
  Gamx = Gamx - temp
  call ddz(ex,gupxz,temp,dZ,ANTI)
  Gamx = Gamx - temp

!
! Gamy:
!
  if(Symmetry==AXISYM) then
     call ddx(ex,gupxy,temp,dX,SYM)
  else
     call ddx(ex,gupxy,temp,dX,ANTI)
  end if
  Gamy = - temp
  call ddy(ex,gupyy,temp,dY,SYM)
  Gamy = Gamy - temp
  call ddz(ex,gupyz,temp,dZ,ANTI)
  Gamy = Gamy - temp
!
! Gamz:
!
  call ddx(ex,gupxz,temp,dX,ANTI)
  Gamz = - temp
  call ddy(ex,gupyz,temp,dY,ANTI)
  Gamz = Gamz - temp
  call ddz(ex,gupzz,temp,dZ,SYM)
  Gamz = Gamz - temp

  istart = nghostzones(1)+1
  jstart = nghostzones(2)+1
  kstart = nghostzones(3)+1
  
  iend = ex(1) - nghostzones(1)
  jend = ex(2) - nghostzones(2)
  kend = ex(3) - nghostzones(3)

  ! Now fill in the interior of the grid with appropriate order derivatives.

  ! 2nd of 2 needed #includes for GenericFD.h:
#include "../../GenFD_set_varF90.h"


  !Following lines needed since nghostzones[0] = ORDER, and
  !   not ORDER-1 in axisymmetry
  !   (so that rotation can be done on multiprocessor runs)
  if(Symmetry==4) then
     istart = istart - 1
     iend = iend + 1
  end if

  write(*,*) "HELLO INSIDE setgamma_v2. nghostzones, nghostzones = ",nghostzones,nghostzones

  do k=kstart,kend
     do j=jstart,jend
        do i=istart,iend
           Gamx(i,j,k) = -(D1gf(gupxx,i,j,k) + D2gf(gupxy,i,j,k) + D3gf(gupxz,i,j,k))
           Gamy(i,j,k) = -(D1gf(gupxy,i,j,k) + D2gf(gupyy,i,j,k) + D3gf(gupyz,i,j,k))
           Gamz(i,j,k) = -(D1gf(gupxz,i,j,k) + D2gf(gupyz,i,j,k) + D3gf(gupzz,i,j,k))
        end do
     end do
  end do

end subroutine setgamma_v2


!-------------------------------------------------------------------
! Compute Gamma^i using old, second-order derivative functions only
!-------------------------------------------------------------------
subroutine setgamma(ex, X, Y, Z, temp, &
     phi, gxx, gxy, gxz, gyy, gyz, gzz, & 
     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Gamx, Gamy, Gamz, &
     Symmetry) !, glob_imin,glob_jmin,glob_kmin,Nx,Nz,  &
  !     PI_SYMM_gupxy, PI_SYMM_gupyy, PI_SYMM_gupyz)
  implicit none
  interface
     subroutine ddx(ext,f,df,dX,sym)
       implicit none
       integer, dimension(3)                    :: ext
       real*8, dimension(ext(1),ext(2),ext(3))  :: f,df
       real*8                                   :: sym,dX
     end subroutine ddx
     subroutine ddy(ext,f,df,dY,sym)
       implicit none
       integer, dimension(3)                    :: ext
       real*8, dimension(ext(1),ext(2),ext(3))  :: f,df
       real*8                                   :: sym,dY
     end subroutine ddy
!     subroutine ddy_pi(ext,f,df,dY,SYM,   &
!		glob_imin,glob_jmin,glob_kmin,Nx,Nz,  &
!		PI_SYMM_f)
!       implicit none
!       integer, dimension(3)                    :: ext
!       real*8, dimension(ext(1),ext(2),ext(3))  :: f,df
!       real*8                                   :: dY,SYM
!       integer                                  :: Nx,Nz
!       integer                                  :: glob_imin,glob_jmin,glob_kmin
!       real*8, dimension(Nx+1,Nz+1)             :: PI_SYMM_f
!     end subroutine ddy_pi
     subroutine ddz(ext,f,df,dZ,sym)
       implicit none
       integer, dimension(3)                    :: ext
       real*8, dimension(ext(1),ext(2),ext(3))  :: f,df
       real*8                                   :: sym,dZ
     end subroutine ddz
  end interface
  integer, dimension(3)                 :: ex
  real*8, dimension(ex(1),ex(2),ex(3))  :: X,Y,Z,temp
  real*8, dimension(ex(1),ex(2),ex(3))  :: phi,gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ex(1),ex(2),ex(3))  :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3))  :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3))  :: Gamx, Gamy, Gamz
  integer                               :: Symmetry
!  integer                               :: Nx,Nz
!  integer                               :: glob_imin,glob_jmin,glob_kmin
!  real*8, dimension(Nx+1,Nz+1)          :: PI_SYMM_gupxy, PI_SYMM_gupyy, PI_SYMM_gupyz
!
!  real*8, dimension(ex(1),ex(2),ex(3))  :: temp
  real*8                                :: dX, dY, dZ
  real*8                                :: gupxxx,gupxyy,gupxzz
  real*8                                :: gupyxx,gupyyy,gupyzz
  real*8                                :: gupzxx,gupzyy,gupzzz
  integer                               :: imin, jmin, kmin, imax, jmax, kmax
  integer                               :: i,j,k
  real*8                                :: F1o3, F1o12, SYM, ANTI
  integer                               :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter ( F1o3 = 1.D0/3.D0, F1o12 = 1.D0/12.D0 )
  parameter ( SYM = 1.D0, ANTI = - 1.D0 )
  parameter (NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
!
! Input translation
!
  imin = lbound(gxx,1)
  jmin = lbound(gxx,2)
  kmin = lbound(gxx,3)
  imax = ubound(gxx,1)
  jmax = ubound(gxx,2)
  kmax = ubound(gxx,3)
  dX = X(imin + 1,1,1) - X(imin,1,1)
  dY = Y(1,jmin + 1,1) - Y(1,jmin,1)
  dZ = Z(1,1,kmin + 1) - Z(1,1,kmin)

!
! Determine Gamx...
!
! Gamx:
!
  call ddx(ex,gupxx,temp,dX,SYM)
  Gamx = - temp
  if (Symmetry.ne.PI_SYMM) then
    call ddy(ex,gupxy,temp,dY,ANTI)
!  elseif(Symmetry==PI_SYMM) then
!    call ddy_pi(ex,gupxy,temp,dY,ANTI,  &
!		glob_imin,glob_jmin,glob_kmin,Nx,Nz,  &
!		PI_SYMM_gupxy)
  endif
  Gamx = Gamx - temp
  call ddz(ex,gupxz,temp,dZ,ANTI)
  Gamx = Gamx - temp

!11  format(6f19.15)

!
! Gamy:
!
  if(Symmetry==AXISYM) then
     call ddx(ex,gupxy,temp,dX,SYM)
  else
     call ddx(ex,gupxy,temp,dX,ANTI)
  end if
  Gamy = - temp
  if (Symmetry.ne.PI_SYMM) then
    call ddy(ex,gupyy,temp,dY,SYM)
!  elseif(Symmetry==PI_SYMM) then
!    call ddy_pi(ex,gupyy,temp,dY,SYM,  &
!		glob_imin,glob_jmin,glob_kmin,Nx,Nz,  &
!		PI_SYMM_gupyy)
  endif
  Gamy = Gamy - temp
  call ddz(ex,gupyz,temp,dZ,ANTI)
  Gamy = Gamy - temp
!
! Gamz:
!
  call ddx(ex,gupxz,temp,dX,ANTI)
  Gamz = - temp
  if (Symmetry.ne.PI_SYMM) then
    call ddy(ex,gupyz,temp,dY,ANTI)
!  elseif(Symmetry==PI_SYMM) then
!    call ddy_pi(ex,gupyz,temp,dY,ANTI,  &
!		glob_imin,glob_jmin,glob_kmin,Nx,Nz,  &
!		PI_SYMM_gupyz)
  endif
  Gamz = Gamz - temp
  call ddz(ex,gupzz,temp,dZ,SYM)
  Gamz = Gamz - temp

! In axisymmetry, the Gamma^a values on the y/=0 planes will
! need to be corrected. 

!  if (Symmetry == AXISYM) then
!    call axibc_vector(ex,X,Y,Z,Gamx,Gamy,Gamz)
!  endif

end subroutine setgamma
