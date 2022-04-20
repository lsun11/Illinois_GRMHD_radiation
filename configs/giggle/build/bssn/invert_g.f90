!-----------------------------------------------------------------------------
!
! $Id: invert_g.F90  Exp $
!
!-----------------------------------------------------------------------------
!
! Invert matrix
!
!-----------------------------------------------------------------------------
subroutine invert_g(ext,gxx,gxy,gxz,gyy,gyz,gzz,gupxx,gupxy,gupxz, &
     gupyy,gupyz,gupzz,detmin,detmax)
  implicit none
!
! Input parameters:
!
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ext(1),ext(2),ext(3))  :: gupxx,gupxy,gupxz
  real*8, dimension(ext(1),ext(2),ext(3))  :: gupyy,gupyz,gupzz
  real*8                                   :: detmin, detmax
!
! Other parameters:
!
  real*8, dimension(ext(1),ext(2),ext(3))  :: detinv
!
! compute and check determinant
!
  detinv =  1.D0/(gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz &
       - gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz)
!  detmin = minval(det)
!  detmax = maxval(det)
! Invert metric...
!
  gupxx =   ( gyy * gzz - gyz * gyz ) * detinv
  gupxy = - ( gxy * gzz - gyz * gxz ) * detinv
  gupxz =   ( gxy * gyz - gyy * gxz ) * detinv
  gupyy =   ( gxx * gzz - gxz * gxz ) * detinv
  gupyz = - ( gxx * gyz - gxy * gxz ) * detinv
  gupzz =   ( gxx * gyy - gxy * gxy ) * detinv
!  write(*,*) invert g: ,gupxy(2,2,15),det(2,2,15),gxx(2,2,15),gxy(2,2,15),gxz(2,2,15),gyy(2,2,15),gyz(2,2,15),gzz(2,2,15)
end subroutine invert_g
