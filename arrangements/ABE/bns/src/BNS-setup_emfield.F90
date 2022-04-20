#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine BNS_setup_emfield(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  interface
     subroutine gderivs_oct(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
     end subroutine gderivs_oct
     subroutine gderivs_eq(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
     end subroutine gderivs_eq
     subroutine gderivs_axi(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
     end subroutine gderivs_axi
  end interface

  integer, dimension(3)                    :: ext
  real*8                                   :: dX,dY,dZ,P_max,Ab
  integer 				   :: handle,index,ierr
  integer                                  :: imin,imax,jmin,jmax,kmin,kmax
  integer                                  :: i,j,k,dummy
  real*8,allocatable,dimension(:,:,:)      :: A_phi,A_phix,A_phiy,A_phiz
  real*8                                   :: psim6,pomega2,al,sqrtg,sqrtg4,B2s
  real*8                                   :: fs4pi,B_xs,B_ys,B_zs,sb0,sb2,xn,yn
  real*8                                   :: sb_x,sb_y,sb_z,psi4,u_x,u_y,u_z
  integer                                  :: adjimin, adjjmin, adjkmin
  integer                                  :: adjimax, adjjmax, adjkmax
  real*8                                   :: multfactor,xmax,dIntegral,pmV,pV
  real*8                                   :: avg_betam1,fac,psin
  integer                                  :: AXISYM,EQUATORIAL
  integer                                  :: OCTANT
  real*8, parameter                        :: SYM = 1.d0, ANTI = -1.d0
  CCTK_REAL reduction_value
  parameter(EQUATORIAL = 1, OCTANT = 2, AXISYM = 4)
!
  fs4pi = sqrt(4.d0*acos(-1.d0))
  ext(1) = cctk_lsh(1)
  ext(2) = cctk_lsh(2)
  ext(3) = cctk_lsh(3)
  dX = X(2,1,1) - X(1,1,1)
  dY = Y(1,2,1) - Y(1,1,1)
  dZ = Z(1,1,2) - Z(1,1,1)

  allocate(A_phi(ext(1),ext(2),ext(3)))
  allocate(A_phix(ext(1),ext(2),ext(3)))
  allocate(A_phiy(ext(1),ext(2),ext(3)))
  allocate(A_phiz(ext(1),ext(2),ext(3)))

  imin = lbound(P,1)
  imax = ubound(P,1)
  jmin = lbound(P,2)
  jmax = ubound(P,2)
  kmin = lbound(P,3)
  kmax = ubound(P,3)

! Find P_max
  call CCTK_VarIndex(index,"mhd_evolve::P")
  call CCTK_ReductionHandle(handle,"maximum")
  if (handle .gt. 0) then
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,reduction_value,1,index)
  else
     call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
  end if
  P_max = reduction_value

! Compute magnetic vector potential A_{\phi} with Ab=1
! Note that A_{\phi} is calculated from the center of each NS
  call BNS_compute_Aphi(A_phi,ext,X,Y,Z,r,PhysicalRadius,P,P_max,p_c,xpc1,xpc2,Sym_Bz)
!
! Compute the derivatives of A_phi
!
 if (Symmetry==EQUATORIAL) then
    call gderivs_eq(ext,A_phi,A_phix,A_phiy,A_phiz,dX,dY,dZ,SYM,SYM ,Sym_Bz )
 else
    write(*,*) 'Symmetry type not supported in BNS_setup_emfield'
    stop
 end if
!
! Now compute B^i for each NS according to 
!  B^x = (-xn/pomega^2) e^(-6 phi) * A_{phi,z};
!  B^y = (-y/pomega^2) e^(-6 phi) * A_{phi,z};
!  B^z = e^(-6 phi) * (xn A_{phi,x} + y A_{phi,y})/pomega^2;
!  pomega^2 = xn^2 + y^2, xn = x - x_c, x_c = center of the NS
!
 do k = kmin,kmax
    do j = jmin,jmax
       do i=imin,imax
          psim6 = exp(-6.d0*phi(i,j,k))
	  psi4 = exp(4.d0*phi(i,j,k))
	  if (X(i,j,k) .gt. 0.d0) then
	     xn = X(i,j,k) - xpc1_fish
	  else
	     xn = X(i,j,k) - xpc2_fish
	  end if
   	  yn = Y(i,j,k)
          pomega2 = xn**2 + yn**2
          Bx(i,j,k) = -xn/pomega2 * psim6 * A_phiz(i,j,k)
          By(i,j,k) = -yn/pomega2 * psim6 * A_phiz(i,j,k)
          Bz(i,j,k) = psim6/pomega2 * (xn*A_phix(i,j,k) +  &
                        yn*A_phiy(i,j,k))
! Compute b^0 and b_i
          al = 1.d0 + lapm1(i,j,k)
          sqrtg = 1.d0/psim6
          sqrtg4 = al * sqrtg
          B2s = psi4*(gxx(i,j,k)*Bx(i,j,k)**2 + &
                2.d0*gxy(i,j,k)*Bx(i,j,k)*By(i,j,k) + &
                2.d0*gxz(i,j,k)*Bx(i,j,k)*Bz(i,j,k) + &
                gyy(i,j,k)*By(i,j,k)**2 + & 
                2.d0*gyz(i,j,k)*By(i,j,k)*Bz(i,j,k) + &
                gzz(i,j,k)*Bz(i,j,k)**2)/(fs4pi*al)**2
	  psin = psi4*u0(i,j,k)
          u_x = ( gxx(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) +  &
                  gxy(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) +  & 
                  gxz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) )*psin
 	  u_y = ( gxy(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) +  & 
                  gyy(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) +  & 
		  gyz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) )*psin
	  u_z = ( gxz(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) +  &
		  gyz(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) +  & 
		  gzz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) )*psin
          sb0 = (u_x*Bx(i,j,k) + u_y*By(i,j,k) + u_z*Bz(i,j,k))/(fs4pi*al)
          sb2 = (B2s + sb0**2)/u0(i,j,k)**2
          ! Temporarily store P_mag=b^2/2 to Ex, and e^(6 phi) to Ey
	  Ex(i,j,k) = sb2*0.5d0
          Ey(i,j,k) = sqrtg
       end do
    end do
  end do
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vars')

  deallocate(A_phi, A_phix, A_phiy, A_phiz)

! Renormalize the B-field so that the mass-averaged P_mag/P = betam1
  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"grid::X")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xmax,1,index)
  call compute_adj_minmax(ext,X,Y,Z,Symmetry, &
       adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax, &
       multfactor,cctk_nghostzones,xmax)
  call wavg(ext,dIntegral,X,Y,Z,Ex,Ey,Symmetry, &
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax)
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,pmV,CCTK_VARIABLE_REAL)
  call wavg(ext,dIntegral,X,Y,Z,P,Ey,Symmetry, &
          adjimin,adjimax,adjjmin,adjjmax,adjkmin,adjkmax)
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dIntegral,pV,CCTK_VARIABLE_REAL)
  avg_betam1 = pmV/pV
  fac = sqrt(betam1/avg_betam1)
  Bx = Bx*fac
  By = By*fac
  Bz = Bz*fac

! Finally compute mhd_st_i and tau
 do k = kmin,kmax
    do j = jmin,jmax
       do i=imin,imax
          psim6 = exp(-6.d0*phi(i,j,k))
	  psi4 = exp(4.d0*phi(i,j,k))
          ! Compute b^0 and b_i
          al = 1.d0 + lapm1(i,j,k)
          sqrtg = 1.d0/psim6
          sqrtg4 = al * sqrtg
          B2s = psi4*(gxx(i,j,k)*Bx(i,j,k)**2 + &
                2.d0*gxy(i,j,k)*Bx(i,j,k)*By(i,j,k) + &
                2.d0*gxz(i,j,k)*Bx(i,j,k)*Bz(i,j,k) + &
                gyy(i,j,k)*By(i,j,k)**2 + & 
		2.d0*gyz(i,j,k)*By(i,j,k)*Bz(i,j,k) + & 
                gzz(i,j,k)*Bz(i,j,k)**2)/(fs4pi*al)**2
          psin = psi4/(al*fs4pi)
          B_xs  = psin * (gxx(i,j,k) * Bx(i,j,k) + gxy(i,j,k) * By(i,j,k) + &
                         gxz(i,j,k) * Bz(i,j,k))
          B_ys  = psin * (gxy(i,j,k) * Bx(i,j,k) + gyy(i,j,k) * By(i,j,k) + &
                         gyz(i,j,k) * Bz(i,j,k))
          B_zs  = psin * (gxz(i,j,k) * Bx(i,j,k) + gyz(i,j,k) * By(i,j,k) + &
                         gzz(i,j,k) * Bz(i,j,k))
          psin = psi4*u0(i,j,k)
          u_x = ( gxx(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) +  &
                  gxy(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) +  &
                  gxz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) )*psin
          u_y = ( gxy(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) +  &
                  gyy(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) +  &
                  gyz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) )*psin
          u_z = ( gxz(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) +  &
                  gyz(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) +  &
                  gzz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) )*psin
          sb0 = (u_x*Bx(i,j,k) + u_y*By(i,j,k) + u_z*Bz(i,j,k))/(fs4pi*al)
          sb2 = (B2s + sb0**2)/u0(i,j,k)**2
          sb_x = (B_xs + u_x*sb0)/u0(i,j,k)
          sb_y = (B_ys + u_y*sb0)/u0(i,j,k)
          sb_z = (B_zs + u_z*sb0)/u0(i,j,k)
          mhd_st_x(i,j,k) = st_x(i,j,k)+sqrtg4*(u0(i,j,k)*sb2*u_x-sb0*sb_x)
          mhd_st_y(i,j,k) = st_y(i,j,k)+sqrtg4*(u0(i,j,k)*sb2*u_y-sb0*sb_y)
          mhd_st_z(i,j,k) = st_z(i,j,k)+sqrtg4*(u0(i,j,k)*sb2*u_z-sb0*sb_z)
          tau(i,j,k) = tau(i,j,k) + sqrtg*( sb2*(al*u0(i,j,k))**2 &
                           - sb2*0.5d0 - (al*sb0)**2 )
       end do
    end do
 end do

end subroutine BNS_setup_emfield

!------------------------------------------------------------------------
! Compute A_phi = pm^2 max( P/P_max - p_c, 0) for Sym_Bz = 1  and 
!         A_phi = pm^2 * z/r * max( P/P_max - p_c, 0) for Sym_Bz = -1.
! where pm^2 = (x-x_c1)^2 + y^2 for x>0, 
!       pm^2 = (x-x_c2)^2 + y^2 for x<0. 
! x_c1 and x_c2 are the position of the center of the right and left 
!       neutron stars respectively.
!------------------------------------------------------------------------
subroutine BNS_compute_Aphi(A_phi,ext,X,Y,Z,r,PhysR,P,P_max,p_c,xpc1,xpc2,Sym_Bz)
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: A_phi,X,Y,Z,P,r,PhysR
  real*8                                   :: P_max,p_c,Sym_Bz,pm2,xpc1,xpc2
  real*8                                   :: xp,yp,xpc1p,xpc2p
  integer                                  :: imin,imax,jmin,jmax,kmin,kmax
  integer                                  :: i,j,k
!
  imin = lbound(P,1)
  imax = ubound(P,1)
  jmin = lbound(P,2)
  jmax = ubound(P,2)
  kmin = lbound(P,3)
  kmax = ubound(P,3)
  do k=kmin,kmax
     do j=jmin,jmax
	do i=imin,imax
	   xp = X(i,j,k)/r(i,j,k)*PhysR(i,j,k)
	   yp = Y(i,j,k)/r(i,j,k)*PhysR(i,j,k)
	   if (xp .gt. 0.d0) then
	      pm2 = (xp-xpc1)**2 + yp**2
	   else
	      pm2 = (xp-xpc2)**2 + yp**2
	   end if   
	   A_phi(i,j,k) = pm2*max(P(i,j,k)/P_max - p_c, 0.d0)
           if (Sym_Bz .lt. 0.d0) A_phi(i,j,k) = A_phi(i,j,k)*Z(i,j,k)/r(i,j,k)
	end do
     end do
  end do
   
end subroutine BNS_compute_Aphi
