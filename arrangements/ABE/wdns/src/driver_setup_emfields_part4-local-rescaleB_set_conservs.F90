#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine driver_setup_emfields_part4_local_rescaleB_set_conservs(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext
  real*8                                   :: dX,dY,dZ,P_max,Ab
  integer                                  :: imin,imax,jmin,jmax,kmin,kmax
  integer                                  :: i,j,k
  real*8                                   :: psim6,pomega2,al,sqrtg,sqrtg4,B2s
  real*8                                   :: fs4pi,B_xs,B_ys,B_zs,sb0,sb2,xn,yn
  real*8                                   :: sb_x,sb_y,sb_z,psi4,u_x,u_y,u_z
  real*8                                   :: fac,psin
  integer                                  :: AXISYM,EQUATORIAL
  integer                                  :: OCTANT
  real*8, parameter                        :: SYM = 1.d0, ANTI = -1.d0
  CCTK_REAL reduction_value
  parameter(EQUATORIAL = 1, OCTANT = 2, AXISYM = 4)
  
  if(CCTK_ITERATION .eq. ITERATION_TO_RESET_MAGNETIC_FIELDS .or. CCTK_ITERATION .eq. 0) then
     !
     fs4pi = sqrt(4.d0*acos(-1.d0))
     ext = cctk_lsh

     dX = X(2,1,1) - X(1,1,1)
     dY = Y(1,2,1) - Y(1,1,1)
     dZ = Z(1,1,2) - Z(1,1,1)

     imin = 1
     jmin = 1
     kmin = 1
     imax = ext(1)
     jmax = ext(2)
     kmax = ext(3)

     ! This was necessary in an older version.  In this version, B is rescaled, and conservatives are set in part2.  Part3 is now only a diagnostic.
!!$  !Rescale Bi's:
!!$  fac = sqrt(betam1/wdns_avg_betam1)
!!$  write(*,*) "setup_emfields_part4FAC=",fac
!!$  Bx = Bx*fac
!!$  By = By*fac
!!$  Bz = Bz*fac
!!$
!!$  !     !For emfields, we assume that you've set Bx, By, Bz (the UN-tilded B^i's)
!!$  !     Bxtilde = Bx
!!$  !     Bytilde = By
!!$  !     Bztilde = Bz
!!$  !
!!$  !     ! Here, we convert B^i to tilde B^i
!!$  !     b2bt = 1.D0
!!$  !     call convert_b(ext,Bxtilde,Bytilde,Bztilde,phi,b2bt)
!!$
!!$  ! Finally compute mhd_st_i and tau
!!$  do k = kmin,kmax
!!$     do j = jmin,jmax
!!$        do i=imin,imax
!!$           psim6 = exp(-6.d0*phi(i,j,k))
!!$           psi4 = exp(4.d0*phi(i,j,k))
!!$           ! Compute b^0 and b_i
!!$           al = 1.d0 + lapm1(i,j,k)
!!$           sqrtg = 1.d0/psim6
!!$           sqrtg4 = al * sqrtg
!!$           B2s = psi4*(gxx(i,j,k)*Bx(i,j,k)**2 + &
!!$                2.d0*gxy(i,j,k)*Bx(i,j,k)*By(i,j,k) + &
!!$                2.d0*gxz(i,j,k)*Bx(i,j,k)*Bz(i,j,k) + &
!!$                gyy(i,j,k)*By(i,j,k)**2 + & 
!!$                2.d0*gyz(i,j,k)*By(i,j,k)*Bz(i,j,k) + & 
!!$                gzz(i,j,k)*Bz(i,j,k)**2)/(fs4pi*al)**2
!!$           psin = psi4/(al*fs4pi)
!!$           B_xs  = psin * (gxx(i,j,k) * Bx(i,j,k) + gxy(i,j,k) * By(i,j,k) + &
!!$                gxz(i,j,k) * Bz(i,j,k))
!!$           B_ys  = psin * (gxy(i,j,k) * Bx(i,j,k) + gyy(i,j,k) * By(i,j,k) + &
!!$                gyz(i,j,k) * Bz(i,j,k))
!!$           B_zs  = psin * (gxz(i,j,k) * Bx(i,j,k) + gyz(i,j,k) * By(i,j,k) + &
!!$                gzz(i,j,k) * Bz(i,j,k))
!!$           psin = psi4*u0(i,j,k)
!!$           u_x = ( gxx(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) +  &
!!$                gxy(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) +  &
!!$                gxz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) )*psin
!!$           u_y = ( gxy(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) +  &
!!$                gyy(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) +  &
!!$                gyz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) )*psin
!!$           u_z = ( gxz(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) +  &
!!$                gyz(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) +  &
!!$                gzz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) )*psin
!!$           sb0 = (u_x*Bx(i,j,k) + u_y*By(i,j,k) + u_z*Bz(i,j,k))/(fs4pi*al)
!!$           sb2 = (B2s + sb0**2)/u0(i,j,k)**2
!!$           sb_x = (B_xs + u_x*sb0)/u0(i,j,k)
!!$           sb_y = (B_ys + u_y*sb0)/u0(i,j,k)
!!$           sb_z = (B_zs + u_z*sb0)/u0(i,j,k)
!!$           mhd_st_x(i,j,k) = st_x(i,j,k)+sqrtg4*(u0(i,j,k)*sb2*u_x-sb0*sb_x)
!!$           mhd_st_y(i,j,k) = st_y(i,j,k)+sqrtg4*(u0(i,j,k)*sb2*u_y-sb0*sb_y)
!!$           mhd_st_z(i,j,k) = st_z(i,j,k)+sqrtg4*(u0(i,j,k)*sb2*u_z-sb0*sb_z)
!!$           if(i==1 .and. j==1 .and. k==1) then
!!$              write(*,*) "TAU111: ",tau(i,j,k)
!!$           end if
!!$           tau(i,j,k) = tau(i,j,k) + sqrtg*( sb2*(al*u0(i,j,k))**2 &
!!$                - sb2*0.5d0 - (al*sb0)**2 )
!!$           if(i==1 .and. j==1 .and. k==1) then
!!$              write(*,*) "ATAU111: ",sqrtg,sb2,al,u0(i,j,k),sb0,tau(i,j,k)
!!$           end if
!!$
!!$        end do
!!$     end do
!!$  end do

!!$  !---------------------------------------
!!$  ! CHECK FOR NANS
!!$  do k=1,cctk_lsh(3)
!!$     do j=1,cctk_lsh(2)
!!$        do i=1,cctk_lsh(1)
!!$           if(isnan(Bz(i,j,k))) then
!!$              write(*,*) "id: found a NAN in Bz at",i,j,k
!!$              stop
!!$           end if
!!$           if(isnan(By(i,j,k))) then
!!$              write(*,*) "id: found a NAN in By at",i,j,k
!!$              stop
!!$           end if
!!$           if(isnan(Bx(i,j,k))) then
!!$              write(*,*) "id: found a NAN in Bx at",i,j,k,cctk_lsh
!!$              stop
!!$           end if
!!$        end do
!!$     end do
!!$  end do
!!$  !---------------------------------------
  end if

end subroutine driver_setup_emfields_part4_local_rescaleB_set_conservs
