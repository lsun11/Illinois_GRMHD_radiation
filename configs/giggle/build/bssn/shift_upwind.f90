!----------------------------------------------------------------------
! Compute advection terms in right hand sides of field equations
! See Shibata in Prog. Theor. Phys., 101, 1199 (1999) (also gr-qc/9905058)
! Note:  centered advective term has already been included in ()_rhs by 
! compute_rhs.f90.
!
!-----------------------------------------------------------------------------
subroutine shift_upwind1(ex,dT,X,Y,Z, &
     gxx,gxy,gxz,gyy,gyz,gzz, &
     gxx_rhs,gxy_rhs,gxz_rhs,gyy_rhs,gyz_rhs,gzz_rhs, &
     Axx,Axy,Axz,Ayy,Ayz,Azz, &
     Axx_rhs,Axy_rhs,Axz_rhs,Ayy_rhs,Ayz_rhs,Azz_rhs, &
     phi, phix, phiy, phiz, phi_rhs, &
     trK, trK_rhs, &
     Gamx, Gamy, Gamz, &
     Gamx_rhs, Gamy_rhs, Gamz_rhs, &
     betax, betay, betaz, temp1, temp2, temp3, Symmetry)
  implicit none
  interface
     subroutine ddx3(ext,f,ddf,s,dX,sym)
       implicit none
       integer, dimension(3)                    :: ext
       real*8, dimension(ext(1),ext(2),ext(3))  :: f,ddf,s
       real*8                                   :: sym,dX
     end subroutine ddx3
     subroutine ddy3(ext,f,ddf,s,dY,sym)
       implicit none
       integer, dimension(3)                    :: ext
       real*8, dimension(ext(1),ext(2),ext(3))  :: f,ddf,s
       real*8                                   :: sym,dY
     end subroutine ddy3
     subroutine ddz3(ext,f,ddf,s,dZ,sym)
       implicit none
       integer, dimension(3)                    :: ext
       real*8, dimension(ext(1),ext(2),ext(3))  :: f,ddf,s
       real*8                                   :: sym,dZ
     end subroutine ddz3
   end interface
!
! Input parameters:
!
  integer, dimension(3)                :: ex
  real*8                               :: dT
  real*8, dimension(ex(1),ex(2),ex(3)) :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gxy,gxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyy,gyz,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx_rhs,gxy_rhs,gxz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyy_rhs,gyz_rhs,gzz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)) :: Axx,Axy,Axz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Axx_rhs,Axy_rhs,Axz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)) :: Ayy_rhs,Ayz_rhs,Azz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)) :: phi,phi_rhs
  real*8, dimension(ex(1),ex(2),ex(3)) :: phix,phiy,phiz
  real*8, dimension(ex(1),ex(2),ex(3)) :: trK, trK_rhs
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamx,Gamy,Gamz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamx_rhs,Gamy_rhs,Gamz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)) :: betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3)) :: temp1, temp2, temp3
  integer                              :: Symmetry
!
! Other variables:
! 
  integer                            :: imin, jmin, kmin, imax, jmax, kmax, i, j, k
  real*8                             :: dX, dY, dZ, s, v
  real*8                             :: F1o3, F1o6, ONE, TWO, FOUR, ZERO
  real*8                             :: F2o3, SIX, EIGHT, HALF, PI
  real*8                             :: SYM, ANTI
  integer                            :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM
  parameter (NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3)
  parameter ( ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0, F1o6 = 1.D0/6.D0 )
  parameter ( ZERO = 0.D0, F1o3 = 1.D0/3.D0, F2o3 = 2.D0/3.D0, SIX = 6.D0 )
  parameter ( SYM = 1.D0, ANTI = - 1.D0, EIGHT = 8.D0, HALF = 0.5D0 )
  PI = acos(-ONE)
!  write(*,*)b_p:,trK(10,10,11),trK_rhs(10,10,11),phi_rhs(10,10,11)
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
  if(Symmetry.eq.PI_SYMM) then
     write(6,*)'PI_SYMM upwind routines not programmed'
     stop
endif
!
! advect gxx
!
  call ddx3(ex,gxx,temp2,temp3,dX,SYM)
  gxx_rhs = gxx_rhs  &
          + ( (ONE - temp3)*abs(betax)*dX + temp3*dT*betax**2 )*HALF*temp2
    call ddy3(ex,gxx,temp2,temp3,dY,SYM)
  gxx_rhs = gxx_rhs  &
          + ( (ONE - temp3)*abs(betay)*dY + temp3*dT*betay**2 )*HALF*temp2
  call ddz3(ex,gxx,temp2,temp3,dZ,SYM)
  gxx_rhs = gxx_rhs  &
          + ( (ONE - temp3)*abs(betaz)*dZ + temp3*dT*betaz**2 )*HALF*temp2
!
! advect gxy
!
  call ddx3(ex,gxy,temp2,temp3,dX,ANTI)
  gxy_rhs = gxy_rhs  &
          + ( (ONE - temp3)*abs(betax)*dX + temp3*dT*betax**2 )*HALF*temp2
    call ddy3(ex,gxy,temp2,temp3,dY,ANTI)
  gxy_rhs = gxy_rhs  &
          + ( (ONE - temp3)*abs(betay)*dY + temp3*dT*betay**2 )*HALF*temp2
  call ddz3(ex,gxy,temp2,temp3,dZ,SYM)
  gxy_rhs = gxy_rhs  &
          + ( (ONE - temp3)*abs(betaz)*dZ + temp3*dT*betaz**2 )*HALF*temp2
!
! advect gxz
!
  call ddx3(ex,gxz,temp2,temp3,dX,ANTI)
  gxz_rhs = gxz_rhs  &
          + ( (ONE - temp3)*abs(betax)*dX + temp3*dT*betax**2 )*HALF*temp2
    call ddy3(ex,gxz,temp2,temp3,dY,SYM)
  gxz_rhs = gxz_rhs  &
          + ( (ONE - temp3)*abs(betay)*dY + temp3*dT*betay**2 )*HALF*temp2
  call ddz3(ex,gxz,temp2,temp3,dZ,ANTI)
  gxz_rhs = gxz_rhs  &
          + ( (ONE - temp3)*abs(betaz)*dZ + temp3*dT*betaz**2 )*HALF*temp2
!
! advect gyy
!
  call ddx3(ex,gyy,temp2,temp3,dX,SYM)
  gyy_rhs = gyy_rhs  &
          + ( (ONE - temp3)*abs(betax)*dX + temp3*dT*betax**2 )*HALF*temp2
    call ddy3(ex,gyy,temp2,temp3,dY,SYM)
  gyy_rhs = gyy_rhs  &
          + ( (ONE - temp3)*abs(betay)*dY + temp3*dT*betay**2 )*HALF*temp2
  call ddz3(ex,gyy,temp2,temp3,dZ,SYM)
  gyy_rhs = gyy_rhs  &
          + ( (ONE - temp3)*abs(betaz)*dZ + temp3*dT*betaz**2 )*HALF*temp2
!
! advect gyz
!
  call ddx3(ex,gyz,temp2,temp3,dX,SYM)
  gyz_rhs = gyz_rhs  &
          + ( (ONE - temp3)*abs(betax)*dX + temp3*dT*betax**2 )*HALF*temp2
    call ddy3(ex,gyz,temp2,temp3,dY,ANTI)
  gyz_rhs = gyz_rhs  &
          + ( (ONE - temp3)*abs(betay)*dY + temp3*dT*betay**2 )*HALF*temp2
  call ddz3(ex,gyz,temp2,temp3,dZ,ANTI)
  gyz_rhs = gyz_rhs  &
          + ( (ONE - temp3)*abs(betaz)*dZ + temp3*dT*betaz**2 )*HALF*temp2
!
! advect gzz
!
  call ddx3(ex,gzz,temp2,temp3,dX,SYM)
  gzz_rhs = gzz_rhs  &
          + ( (ONE - temp3)*abs(betax)*dX + temp3*dT*betax**2 )*HALF*temp2
    call ddy3(ex,gzz,temp2,temp3,dY,SYM)
  gzz_rhs = gzz_rhs  &
          + ( (ONE - temp3)*abs(betay)*dY + temp3*dT*betay**2 )*HALF*temp2
  call ddz3(ex,gzz,temp2,temp3,dZ,SYM)
  gzz_rhs = gzz_rhs  &
          + ( (ONE - temp3)*abs(betaz)*dZ + temp3*dT*betaz**2 )*HALF*temp2
!
! advect Axx
!
  call ddx3(ex,Axx,temp2,temp3,dX,SYM)
  Axx_rhs = Axx_rhs  &
          + ( (ONE - temp3)*abs(betax)*dX + temp3*dT*betax**2 )*HALF*temp2
    call ddy3(ex,Axx,temp2,temp3,dY,SYM)
  Axx_rhs = Axx_rhs  &
          + ( (ONE - temp3)*abs(betay)*dY + temp3*dT*betay**2 )*HALF*temp2
  call ddz3(ex,Axx,temp2,temp3,dZ,SYM)
  Axx_rhs = Axx_rhs  &
          + ( (ONE - temp3)*abs(betaz)*dZ + temp3*dT*betaz**2 )*HALF*temp2
!
! advect Axy
!
  call ddx3(ex,Axy,temp2,temp3,dX,ANTI)
  Axy_rhs = Axy_rhs  &
          + ( (ONE - temp3)*abs(betax)*dX + temp3*dT*betax**2 )*HALF*temp2
    call ddy3(ex,Axy,temp2,temp3,dY,ANTI)
  Axy_rhs = Axy_rhs  &
          + ( (ONE - temp3)*abs(betay)*dY + temp3*dT*betay**2 )*HALF*temp2
  call ddz3(ex,Axy,temp2,temp3,dZ,SYM)
  Axy_rhs = Axy_rhs  &
          + ( (ONE - temp3)*abs(betaz)*dZ + temp3*dT*betaz**2 )*HALF*temp2
!
! advect Axz
!
  call ddx3(ex,Axz,temp2,temp3,dX,ANTI)
  Axz_rhs = Axz_rhs  &
          + ( (ONE - temp3)*abs(betax)*dX + temp3*dT*betax**2 )*HALF*temp2
    call ddy3(ex,Axz,temp2,temp3,dY,SYM)
  Axz_rhs = Axz_rhs  &
          + ( (ONE - temp3)*abs(betay)*dY + temp3*dT*betay**2 )*HALF*temp2
  call ddz3(ex,Axz,temp2,temp3,dZ,ANTI)
  Axz_rhs = Axz_rhs  &
          + ( (ONE - temp3)*abs(betaz)*dZ + temp3*dT*betaz**2 )*HALF*temp2
!
! advect Ayy
!
  call ddx3(ex,Ayy,temp2,temp3,dX,SYM)
  Ayy_rhs = Ayy_rhs  &
          + ( (ONE - temp3)*abs(betax)*dX + temp3*dT*betax**2 )*HALF*temp2
    call ddy3(ex,Ayy,temp2,temp3,dY,SYM)
  Ayy_rhs = Ayy_rhs  &
          + ( (ONE - temp3)*abs(betay)*dY + temp3*dT*betay**2 )*HALF*temp2
  call ddz3(ex,Ayy,temp2,temp3,dZ,SYM)
  Ayy_rhs = Ayy_rhs  &
          + ( (ONE - temp3)*abs(betaz)*dZ + temp3*dT*betaz**2 )*HALF*temp2
!
! advect Ayz
!
  call ddx3(ex,Ayz,temp2,temp3,dX,SYM)
  Ayz_rhs = Ayz_rhs  &
          + ( (ONE - temp3)*abs(betax)*dX + temp3*dT*betax**2 )*HALF*temp2
    call ddy3(ex,Ayz,temp2,temp3,dY,ANTI)
  Ayz_rhs = Ayz_rhs  &
          + ( (ONE - temp3)*abs(betay)*dY + temp3*dT*betay**2 )*HALF*temp2
  call ddz3(ex,Ayz,temp2,temp3,dZ,ANTI)
  Ayz_rhs = Ayz_rhs  &
          + ( (ONE - temp3)*abs(betaz)*dZ + temp3*dT*betaz**2 )*HALF*temp2
!
! advect Azz
!
  call ddx3(ex,Azz,temp2,temp3,dX,SYM)
  Azz_rhs = Azz_rhs  &
          + ( (ONE - temp3)*abs(betax)*dX + temp3*dT*betax**2 )*HALF*temp2
    call ddy3(ex,Azz,temp2,temp3,dY,SYM)
  Azz_rhs = Azz_rhs  &
          + ( (ONE - temp3)*abs(betay)*dY + temp3*dT*betay**2 )*HALF*temp2
  call ddz3(ex,Azz,temp2,temp3,dZ,SYM)
  Azz_rhs = Azz_rhs  &
          + ( (ONE - temp3)*abs(betaz)*dZ + temp3*dT*betaz**2 )*HALF*temp2
!
! advect phi
!
  call ddx3(ex,phi,temp2,temp3,dX,SYM)
  phi_rhs = phi_rhs  &
          + ( (ONE - temp3)*abs(betax)*dX + temp3*dT*betax**2 )*HALF*temp2
    call ddy3(ex,phi,temp2,temp3,dY,SYM)
  phi_rhs = phi_rhs  &
          + ( (ONE - temp3)*abs(betay)*dY + temp3*dT*betay**2 )*HALF*temp2
  call ddz3(ex,phi,temp2,temp3,dZ,SYM)
  phi_rhs = phi_rhs  &
          + ( (ONE - temp3)*abs(betaz)*dZ + temp3*dT*betaz**2 )*HALF*temp2
!
! advect trK
!
  call ddx3(ex,trK,temp2,temp3,dX,SYM)
  trK_rhs = trK_rhs  &
          + ( (ONE - temp3)*abs(betax)*dX + temp3*dT*betax**2 )*HALF*temp2
!  write(*,*)trk1:,trK_rhs(10,10,11),betax(10,10,11),temp2(10,10,11),temp3(10,10,11),dX,dY,dZ
   call ddy3(ex,trK,temp2,temp3,dY,SYM)
  trK_rhs = trK_rhs  &
          + ( (ONE - temp3)*abs(betay)*dY + temp3*dT*betay**2 )*HALF*temp2
!  write(6,*)trk2:,trK_rhs(10,10,11),betay(10,10,11),temp2(10,10,11),temp3(10,10,11)
  call ddz3(ex,trK,temp2,temp3,dZ,SYM)
  trK_rhs = trK_rhs  &
          + ( (ONE - temp3)*abs(betaz)*dZ + temp3*dT*betaz**2 )*HALF*temp2
!  write(*,*)trk3:,trK_rhs(10,10,11)
!
! advect Gamx
!
  call ddx3(ex,Gamx,temp2,temp3,dX,ANTI)
  Gamx_rhs = Gamx_rhs  &
          + ( (ONE - temp3)*abs(betax)*dX + temp3*dT*betax**2 )*HALF*temp2
    call ddy3(ex,Gamx,temp2,temp3,dY,SYM)
  Gamx_rhs = Gamx_rhs  &
          + ( (ONE - temp3)*abs(betay)*dY + temp3*dT*betay**2 )*HALF*temp2
  call ddz3(ex,Gamx,temp2,temp3,dZ,SYM)
  Gamx_rhs = Gamx_rhs  &
          + ( (ONE - temp3)*abs(betaz)*dZ + temp3*dT*betaz**2 )*HALF*temp2
!
! advect Gamy
!
  call ddx3(ex,Gamy,temp2,temp3,dX,SYM)
  Gamy_rhs = Gamy_rhs  &
          + ( (ONE - temp3)*abs(betax)*dX + temp3*dT*betax**2 )*HALF*temp2
    call ddy3(ex,Gamy,temp2,temp3,dY,ANTI)
  Gamy_rhs = Gamy_rhs  &
          + ( (ONE - temp3)*abs(betay)*dY + temp3*dT*betay**2 )*HALF*temp2
  call ddz3(ex,Gamy,temp2,temp3,dZ,SYM)
  Gamy_rhs = Gamy_rhs  &
          + ( (ONE - temp3)*abs(betaz)*dZ + temp3*dT*betaz**2 )*HALF*temp2
!
! advect Gamz
!
  call ddx3(ex,Gamz,temp2,temp3,dX,SYM)
  Gamz_rhs = Gamz_rhs  &
          + ( (ONE - temp3)*abs(betax)*dX + temp3*dT*betax**2 )*HALF*temp2
    call ddy3(ex,Gamz,temp2,temp3,dY,SYM)
  Gamz_rhs = Gamz_rhs  &
          + ( (ONE - temp3)*abs(betay)*dY + temp3*dT*betay**2 )*HALF*temp2
  call ddz3(ex,Gamz,temp2,temp3,dZ,ANTI)
  Gamz_rhs = Gamz_rhs  &
          + ( (ONE - temp3)*abs(betaz)*dZ + temp3*dT*betaz**2 )*HALF*temp2
end subroutine shift_upwind1
