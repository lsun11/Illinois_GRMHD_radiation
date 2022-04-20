#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine WDNS_setup_emfield_part2_local_setBi(CCTK_ARGUMENTS)
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
  real*8                                   :: dX,dY,dZ,Ab
  integer 				   :: handle,index,ierr
  integer                                  :: imin,imax,jmin,jmax,kmin,kmax
  integer                                  :: i,j,k,dummy
  integer                                  :: im1,jm1,km1,ip1,jp1,kp1
  real*8,allocatable,dimension(:,:,:)      :: A_phi,A_phix,A_phiy,A_phiz
  real*8                                   :: psim6,psim6_s,pomega2,al,sqrtg,sqrtg4,B2s
  real*8                                   :: fs4pi,B_xs,B_ys,B_zs,sb0,sb2,xn,yn
  real*8                                   :: sb_x,sb_y,sb_z,psi4,u_x,u_y,u_z
  real*8                                   :: psin,x_NS_CoM_coord,y_NS_CoM_coord
  real*8                                   :: Yijk,Yijkp1,Yijp1k,Yijp1kp1,Yip1jk,Yip1jkp1,Yip1jp1k,Yip1jp1kp1
  real*8                                   :: Xijk,Xijkp1,Xijp1k,Xijp1kp1,Xip1jk,Xip1jkp1,Xip1jp1k,Xip1jp1kp1
  real*8				   :: Ax000,Ax001,Ax010,Ax011,Ax100,Ax101,Ax110,Ax111
  real*8                                   :: Ay000,Ay001,Ay010,Ay011,Ay100,Ay101,Ay110,Ay111
  integer                                  :: AXISYM,EQUATORIAL
  integer                                  :: OCTANT
  real*8, parameter                        :: SYM = 1.d0, ANTI = -1.d0
  CCTK_REAL reduction_value
  parameter(EQUATORIAL = 1, OCTANT = 2, AXISYM = 4)
  !


  write(*,*) 'Inside WDNS_setup_emfield_part2_local_setBi'

  fs4pi = sqrt(4.d0*acos(-1.d0))
  ext = cctk_lsh

  dX = X(2,1,1) - X(1,1,1)
  dY = Y(1,2,1) - Y(1,1,1)
  dZ = Z(1,1,2) - Z(1,1,1)

  allocate(A_phi(ext(1),ext(2),ext(3)))
  allocate(A_phix(ext(1),ext(2),ext(3)))
  allocate(A_phiy(ext(1),ext(2),ext(3)))
  allocate(A_phiz(ext(1),ext(2),ext(3)))

  imin = 1
  jmin = 1
  kmin = 1
  imax = ext(1)
  jmax = ext(2)
  kmax = ext(3)

  ! Compute magnetic vector potential A_{\phi} with Ab=1
  ! Note that A_{\phi} is calculated from the center of each NS

  !     x_NS_CoM_coord = CoMx_VolInt/CoM_VolInt_denominator
  !     y_NS_CoM_coord = CoMy_VolInt/CoM_VolInt_denominator

  x_NS_CoM_coord = initial_ns_coord_x
  y_NS_CoM_coord = initial_ns_coord_y

  write(*,*) "INSIDE EMFIELDS SETUP: x,y of ns:",x_NS_CoM_coord,y_NS_CoM_coord
  write(*,*) "INSIDE EMFIELDS SETUP: Pmax:",wdns_P_max,p_c

  if(em_field_type==1 .and. constrained_transport_scheme .ne. 3) then
     write(*,*) "Sorry, em_field_type==1 (toroidal fields) not yet supported for constrained_transport_scheme != 3."
     stop
  end if
  
  if (em_field_type==1 .and. Sym_Bz .gt. 0.d0) then 
     write(*,*) "Sorry, em_field_type==1 (toroidal fields) not yet supported with Sym_Bz==1.  You can set Sym_Bz=-1 !"
     stop
  end if

  if (enable_trace_field_line==1 .and. em_field_type==1) then 
     write(*,*) "Sorry, the initial data for the field line tracer variables currently only supports the poloidal initial data. You have to write your own initial data for the toroidal configuration."
  end if

  !
  ! Compute the vector potential A_phi 
  !

  if (em_field_type .ne. 1) then 
     call WDNS_compute_Aphi(ext,X,Y,Z,PhysicalRadius,P,A_phi,Ax,Ay,Az, & 
       mhd_psi_line, mhd_u_psi, mhd_chi_line, mhd_u_chi, &
       betam1,p_c,x_NS_CoM_coord,y_NS_CoM_coord,wdns_P_max,Sym_Bz, & 
       enable_trace_field_line,constrained_transport_scheme,em_field_type)
  end if

  !
  ! Compute the vector potential Ax and Ay for near toroidal configuration
  !
  if (em_field_type==1) then 
     call WDNS_compute_A_toroidal(ext,X,Y,Z,P,Ax,Ay,Az,phi,betam1, &
        p_c,x_NS_CoM_coord,y_NS_CoM_coord,r0,wdns_P_max)
  end if

  !
  ! Compute the derivatives of A_phi, if we are not using constrained_transport_scheme==3
  !
  if (constrained_transport_scheme .ne. 3) then 
     if (Symmetry==EQUATORIAL) then
        call gderivs_eq(ext,A_phi,A_phix,A_phiy,A_phiz,dX,dY,dZ,SYM,SYM ,Sym_Bz )
     else
        write(*,*) 'Symmetry type not supported in WDNS_setup_emfield'
        stop
     end if
  end if

  ! Compute B^i on staggered grid and temporarily store them in A_phix,A_phiy,A_phiz
  !
  if (constrained_transport_scheme==3) then 
     do k=1,ext(3)
        do j=1,ext(2)
           do i=1,ext(1)
              im1 = max(i-1,1)
              jm1 = max(j-1,1)
              km1 = max(k-1,1)
              ip1 = min(i+1,ext(1))
              jp1 = min(j+1,ext(2))
              kp1 = min(k+1,ext(3))

              psim6_s = exp(-3.d0 * (phi(i,j,k) + phi(ip1,j,k)) )
              A_phix(i,j,k) = ( (Az(i,j,k)-Az(i,jm1,k))/dY   & 
                   - (Ay(i,j,k)-Ay(i,j,km1))/dZ ) * psim6_s

              psim6_s = exp(-3.d0 * (phi(i,j,k) + phi(i,jp1,k)) )
              A_phiy(i,j,k) = ( (Ax(i,j,k)-Ax(i,j,km1))/dZ & 
                   - (Az(i,j,k)-Az(im1,j,k))/dX ) * psim6_s

              psim6_s = exp(-3.d0 * (phi(i,j,k) + phi(i,j,kp1)) )
              A_phiz(i,j,k) = ( (Ay(i,j,k)-Ay(im1,j,k))/dX &
                   - (Ax(i,j,k)-Ax(i,jm1,k))/dY ) * psim6_s

           end do
        end do
     end do



     ! Now compute B^i on unstaggered grid by simple averge
     do k=1,ext(3)
        do j=1,ext(2)
           do i=1,ext(1)
              im1 = max(i-1,1)
              jm1 = max(j-1,1)
              km1 = max(k-1,1)
              Bx(i,j,k) = 0.5d0* (A_phix(i,j,k) + A_phix(im1,j,k))
              By(i,j,k) = 0.5d0* (A_phiy(i,j,k) + A_phiy(i,jm1,k))
              Bz(i,j,k) = 0.5d0* (A_phiz(i,j,k) + A_phiz(i,j,km1))
           end do
        end do
     end do
  end if



  !
  ! Now compute B^i according to (exercise for the readers)
  !  B^x = (-x/pomega^2) e^(-6 phi) * A_{phi,z}; 
  !  B^y = (-y/pomega^2) e^(-6 phi) * A_{phi,z};
  !  B^z = e^(-6 phi) * (x A_{phi,x} + y A_{phi,y})/pomega^2; 
  !  pomega^2 = x^2 + y^2
  !
  ! and then calculate mhd_st_i and tau
  !
  do k = kmin,kmax
     do j = jmin,jmax
        do i = imin,imax
           xn = X(i,j,k) - x_NS_CoM_coord
           yn = Y(i,j,k) - y_NS_CoM_coord

           psim6 = exp(-6.d0*phi(i,j,k))
           pomega2 = xn**2 + yn**2
	   if (constrained_transport_scheme==3) then 
	      ! do nothing since B^i has been computed 
   	   elseif (em_field_type==0 .or. i==imax .or. j==jmax .or. k==kmax) then 
              Bx(i,j,k) = -xn/pomega2 * psim6 * A_phiz(i,j,k)
              By(i,j,k) = -yn/pomega2 * psim6 * A_phiz(i,j,k)
              Bz(i,j,k) = psim6/pomega2 * (xn*A_phix(i,j,k) +  & 
                   yn*A_phiy(i,j,k))
	   else
	      ! Compute Bx, By, Bz that satisfy div(B)=0 to machine precision.
  	      ! Here is the recipe, which has been verified by Mathematica: 
	      !
	      ! \tilde{B}^x(i,j,k)=-(Ay(i,j,k+1)-Ay(i,j,k) + Ay(i+1,j,k+1)-Ay(i+1,j,k) + Ay(i,j+1,k+1)-Ay(i,j+1,k) + Ay(i+1,j+1,k+1)-Ay(i+1,j+1,k))/(4*dz) + (Az(i,j+1,k)-Az(i,j,k) + Az(i+1,j+1,k)-Az(i+1,j,k) + Az(i,j+1,k+1)-Az(i,j,k+1) + Az(i+1,j+1,k+1)-Az(i+1,j,k+1))/(4*dy)
	      ! \tilde{B}^y(i,j,k)=(Ax(i,j,k+1)-Ax(i,j,k) + Ax(i+1,j,k+1)-Ax(i+1,j,k) + Ax(i,j+1,k+1)-Ax(i,j+1,k) + Ax(i+1,j+1,k+1)-Ax(i+1,j+1,k))/(4*dz) - (Az(i+1,j,k)-Az(i,j,k) + Az(i+1,j+1,k)-Az(i,j+1,k) + Az(i+1,j,k+1)-Az(i,j,k+1) + Az(i+1,j+1,k+1)-Az(i,j+1,k+1))/(4*dx)
	      ! \tilde{B}^z(i,j,k):=(Ay(i+1,j,k)-Ay(i,j,k) + Ay(i+1,j+1,k)-Ay(i,j+1,k) + Ay(i+1,j,k+1)-Ay(i,j,k+1) + Ay(i+1,j+1,k+1)-Ay(i,j+1,k+1))/(4*dx) - (Ax(i,j+1,k)-Ax(i,j,k) + Ax(i+1,j+1,k)-Ax(i+1,j,k) + Ax(i,j+1,k+1)-Ax(i,j,k+1) + Ax(i+1,j+1,k+1)-Ax(i+1,j,k+1))/(4*dy)
              ! Here \tilde{B}^i = psi^6 B^i, and 
	      ! Ax(i,j,k), Ay(i,j,k), and Az(i,j,k) are the 
              ! 3 (covariant) components of the vector potential at the corner 
	      ! point (x_i-dx/2, y_j-dy/2, z_k-dz/2).
    	      ! In the present case, we set Ax = -A_phi * y/(x^2+y^2), 
	      !  Ay = A_phi * x/(x^2+y^2), and Az = 0.
	      ! Note that B^i can't be set this way at the boundary points
	      ! becuase of the sentcil structure. Presumably, this 
	      ! is not a problem when we have enough ghostzones.
	      !
	      ! In the following, Axabc denotes Ax(i+a,j+b,k+c) and 
	      ! similarly for Ayabc.
	      ! 

              ! Here we set the coordinates, relative to the origin, which in this case is the center of the NS.
              Yijk      = yn
              Yijkp1    = yn
              Yijp1k    = yn + dY
              Yijp1kp1  = yn + dY
              Yip1jk    = yn
              Yip1jkp1  = yn
              Yip1jp1k  = yn + dY
              Yip1jp1kp1= yn + dY

              Xijk      = xn
              Xijkp1    = xn
              Xijp1k    = xn
              Xijp1kp1  = xn
              Xip1jk    = xn + dX
              Xip1jkp1  = xn + dX
              Xip1jp1k  = xn + dX
              Xip1jp1kp1= xn + dX

              Ax000 = -A_phi(i,j,k)* (Yijk-0.5d0*dY)/ & 
                   ((Xijk-0.5d0*dX)**2 + (Yijk-0.5d0*dY)**2)
              Ax001 = -A_phi(i,j,k+1)*(Yijkp1-0.5d0*dY)/ & 
                   ((Xijkp1-0.5d0*dX)**2 + (Yijkp1-0.5d0*dY)**2)
              Ax010 = -A_phi(i,j+1,k)* (Yijp1k-0.5d0*dY)/ & 
                   ((Xijp1k-0.5d0*dX)**2 + (Yijp1k-0.5d0*dY)**2)
              Ax011 = -A_phi(i,j+1,k+1)*(Yijp1kp1-0.5d0*dY)/ & 
                   ((Xijp1kp1-0.5d0*dX)**2 + (Yijp1kp1-0.5d0*dY)**2)
              Ax100 = -A_phi(i+1,j,k)*(Yip1jk-0.5d0*dY)/ & 
                   ((Xip1jk-0.5d0*dX)**2 + (Yip1jk-0.5d0*dY)**2)
              Ax101 = -A_phi(i+1,j,k+1)*(Yip1jkp1-0.5d0*dY)/ & 
                   ((Xip1jkp1-0.5d0*dX)**2 + (Yip1jkp1-0.5d0*dY)**2)
              Ax110 = -A_phi(i+1,j+1,k)*(Yip1jp1k-0.5d0*dY)/ & 
                   ((Xip1jp1k-0.5d0*dX)**2 + (Yip1jp1k-0.5d0*dY)**2)
              Ax111 = -A_phi(i+1,j+1,k+1)*(Yip1jp1kp1-0.5d0*dY)/ & 
                   ((Xip1jp1kp1-0.5d0*dX)**2 + (Yip1jp1kp1-0.5d0*dY)**2)
              Ay000 = A_phi(i,j,k)* (Xijk-0.5d0*dX)/ &
                   ((Xijk-0.5d0*dX)**2 + (Yijk-0.5d0*dY)**2)
              Ay001 = A_phi(i,j,k+1)*(Xijkp1-0.5d0*dX)/ &
                   ((Xijkp1-0.5d0*dX)**2 + (Yijkp1-0.5d0*dY)**2)
              Ay010 = A_phi(i,j+1,k)* (Xijp1k-0.5d0*dX)/ &
                   ((Xijp1k-0.5d0*dX)**2 + (Yijp1k-0.5d0*dY)**2)
              Ay011 = A_phi(i,j+1,k+1)*(Xijp1kp1-0.5d0*dX)/ &
                   ((Xijp1kp1-0.5d0*dX)**2 + (Yijp1kp1-0.5d0*dY)**2)
              Ay100 = A_phi(i+1,j,k)*(Xip1jk-0.5d0*dX)/ &
                   ((Xip1jk-0.5d0*dX)**2 + (Yip1jk-0.5d0*dY)**2)
              Ay101 = A_phi(i+1,j,k+1)*(Xip1jkp1-0.5d0*dX)/ &
                   ((Xip1jkp1-0.5d0*dX)**2 + (Yip1jkp1-0.5d0*dY)**2)
              Ay110 = A_phi(i+1,j+1,k)*(Xip1jp1k-0.5d0*dX)/ &
                   ((Xip1jp1k-0.5d0*dX)**2 + (Yip1jp1k-0.5d0*dY)**2)
              Ay111 = A_phi(i+1,j+1,k+1)*(Xip1jp1kp1-0.5d0*dX)/ &
                   ((Xip1jp1kp1-0.5d0*dX)**2 + (Yip1jp1kp1-0.5d0*dY)**2)

	      Bx(i,j,k) = -( (Ay001-Ay000) + (Ay101-Ay100) + (Ay011-Ay010) + (Ay111-Ay110) ) * 0.25d0/dZ * psim6
	      By(i,j,k) = ( (Ax001-Ax000) + (Ax101-Ax100) + (Ax011-Ax010) + (Ax111-Ax110) ) * 0.25d0/dZ * psim6
	      Bz(i,j,k) = ( (Ay100-Ay000) + (Ay110-Ay010) + (Ay101-Ay001) + (Ay111-Ay011) ) * 0.25d0/dX * psim6 - ( (Ax010-Ax000) + (Ax110-Ax100) + (Ax011-Ax001) + (Ax111-Ax101) ) * 0.25d0/dY * psim6
	   end if
           ! Compute b^0 and b_i
           al = 1.d0 + lapm1(i,j,k)
           sqrtg = 1.d0/psim6
           sqrtg4 = al * sqrtg
           B2s = exp(4.d0*phi(i,j,k))*(gxx(i,j,k)*Bx(i,j,k)**2 + & 
		2.d0*gxy(i,j,k)*Bx(i,j,k)*By(i,j,k) + & 
		2.d0*gxz(i,j,k)*Bx(i,j,k)*Bz(i,j,k) + &
        	gyy(i,j,k)*By(i,j,k)**2 + 2.d0*gyz(i,j,k)*By(i,j,k)*Bz(i,j,k) + & 
		gzz(i,j,k)*Bz(i,j,k)**2)/(fs4pi*al)**2
           psin = exp(4.d0*phi(i,j,k))/al/fs4pi
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


           sb0 = (u_x*Bx(i,j,k) + u_y*By(i,j,k) + &
                u_z*Bz(i,j,k))/fs4pi/al
           sb2 = (B2s + sb0**2)/u0(i,j,k)**2
           !Branson's way of ensuring b^2/P is what we want: 
           !  Call this routine once to calibrate sb2/P, then average it.
           !   Scale sb2 accordingly, then call this function again to obtain the desired sb2/P
           ! Here's the necessary line of code.  You'll need to add the input parameters as well.
           !bsq(i,j,k) = sb2/P(i,j,k)

           sb_x = (B_xs + u_x*sb0)/u0(i,j,k)
           sb_y = (B_ys + u_y*sb0)/u0(i,j,k)
           sb_z = (B_zs + u_z*sb0)/u0(i,j,k)
           ! Now compute mhd_st_i and tau
           mhd_st_x(i,j,k) = st_x(i,j,k)+sqrtg4*(u0(i,j,k)*sb2*u_x-sb0*sb_x)
           mhd_st_y(i,j,k) = st_y(i,j,k)+sqrtg4*(u0(i,j,k)*sb2*u_y-sb0*sb_y)
           mhd_st_z(i,j,k) = st_z(i,j,k)+sqrtg4*(u0(i,j,k)*sb2*u_z-sb0*sb_z)
           tau(i,j,k) = tau(i,j,k) + sqrtg*( sb2*(al*u0(i,j,k))**2 &
                - sb2*0.5d0 - (al*sb0)**2 )

           ! CHECK NEXT FUNCTION IN SCHEDULE.CCL for reason why temp1 and temp2 are defined.
           ! Temporarily store e^(6 phi)*b^2/2 (b^2/2 = P_mag) to temp1, and e^(6 phi)*P(i,j,k) to temp2, e^(6 phi) to temp3 
           temp1(i,j,k) = sb2*0.5d0*sqrtg
           !if(sb2.gt.0.D0) write(*,*) sb2
           temp2(i,j,k) = P(i,j,k)*sqrtg
           if (rho_b(i,j,k) .gt. rho_b_atm*1.d5) then
              temp3(i,j,k) = 1.d0
           else
              temp3(i,j,k) = 0.d0
           end if
        end do
     end do
  end do


	if (em_gauge.ne.0) psi6phi = 0.d0


	write(*,*) "Checking for nans in Bx, By, Bz as initial data"
		 call check_for_nans(cctkGH,cctk_lsh,Bx,By,Bz)

	write(*,*) "Checking for nans in Ax, Ay, Az as initial data"
		 call check_for_nans(cctkGH,cctk_lsh,Ax,Ay,Az)




!!$  !
!!$  ! Now compute UNSCALED B^i for each NS according to 
!!$  !  B^x = (-xn/pomega^2) e^(-6 phi) * A_{phi,z};
!!$  !  B^y = (-y/pomega^2) e^(-6 phi) * A_{phi,z};
!!$  !  B^z = e^(-6 phi) * (xn A_{phi,x} + y A_{phi,y})/pomega^2;
!!$  !  pomega^2 = xn^2 + y^2, xn = x - x_c, x_c = center of the NS
!!$  !  
!!$  !  Note that we rescale B^i to the desired amplitude in driver_setup_emfields_part4...
!!$  !
!!$  do k = kmin,kmax
!!$     do j = jmin,jmax
!!$        do i=imin,imax
!!$           psim6 = exp(-6.d0*phi(i,j,k))
!!$           psi4 = exp(4.d0*phi(i,j,k))
!!$
!!$           xn = X(i,j,k) - x_NS_CoM_coord
!!$           yn = Y(i,j,k) - y_NS_CoM_coord
!!$
!!$           pomega2 = xn**2 + yn**2
!!$
!!$           if(pomega2==0.D0) then
!!$              write(*,*) "BAD POMEGA2",i,j,k
!!$           end if
!!$
!!$           Bx(i,j,k) = -xn/pomega2 * psim6 * A_phiz(i,j,k)
!!$           By(i,j,k) = -yn/pomega2 * psim6 * A_phiz(i,j,k)
!!$           Bz(i,j,k) = psim6/pomega2 * (xn*A_phix(i,j,k) +  &
!!$                yn*A_phiy(i,j,k))
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
!!$           ! Temporarily store e^(6 phi)*b^2/2 (b^2/2 = P_mag) to temp1, and e^(6 phi)*P(i,j,k) to temp2
!!$           temp1(i,j,k) = sb2*0.5d0*sqrtg
!!$           !if(sb2.gt.0.D0) write(*,*) sb2
!!$           temp2(i,j,k) = P(i,j,k)*sqrtg
!!$
!!$        end do
!!$     end do
!!$  end do

  deallocate(A_phi, A_phix, A_phiy, A_phiz)

end subroutine WDNS_setup_emfield_part2_local_setBi


! OLD CODE, ONLY WORKS WITH constrained_transport_scheme!=3
!!$!------------------------------------------------------------------------
!!$! Compute A_phi = pm^2 max( P/wdns_P_max - p_c, 0) for Sym_Bz = 1  and 
!!$!         A_phi = pm^2 * z/r * max( P/wdns_P_max - p_c, 0) for Sym_Bz = -1.
!!$! where pm^2 = (x-x_ns_com)^2 + (y-y_ns_com)^2 
!!$! x_c1 and x_c2 are the position of the center of the right and left 
!!$!       neutron stars respectively.
!!$!------------------------------------------------------------------------
!!$subroutine WDNS_compute_Aphi(A_phi,ext,X,Y,Z,r,PhysR,P,wdns_P_max,p_c,x_NS_CoM_coord,y_NS_CoM_coord,Sym_Bz)
!!$  implicit none
!!$  integer, dimension(3)                    :: ext
!!$  real*8, dimension(ext(1),ext(2),ext(3))  :: A_phi,X,Y,Z,P,r,PhysR
!!$  real*8                                   :: wdns_P_max,p_c,Sym_Bz,pm2,x_NS_CoM_coord,y_NS_CoM_coord
!!$  real*8                                   :: xp,yp,x_NS_CoM_coordp,y_NS_CoM_coordp
!!$  integer                                  :: imin,imax,jmin,jmax,kmin,kmax
!!$  integer                                  :: i,j,k
!!$  !
!!$  imin = lbound(P,1)
!!$  imax = ubound(P,1)
!!$  jmin = lbound(P,2)
!!$  jmax = ubound(P,2)
!!$  kmin = lbound(P,3)
!!$  kmax = ubound(P,3)
!!$  do k=kmin,kmax
!!$     do j=jmin,jmax
!!$	do i=imin,imax
!!$           xp = X(i,j,k)
!!$           yp = Y(i,j,k)
!!$
!!$           pm2 = (xp-x_NS_CoM_coord)**2 + (yp-y_NS_CoM_coord)**2
!!$
!!$	   A_phi(i,j,k) = pm2*max(P(i,j,k)/wdns_P_max - p_c, 0.d0)
!!$           if (Sym_Bz .lt. 0.d0) A_phi(i,j,k) = A_phi(i,j,k)*Z(i,j,k)/r(i,j,k)
!!$	end do
!!$     end do
!!$  end do
!!$
!!$end subroutine WDNS_compute_Aphi


!------------------------------------------------------------------------
! Compute A_phi = pm^2 * sqrt(8 pi wdns_P_max * betam1) * max( P/wdns_P_max - p_c, 0) for Sym_Bz = 1  and 
!         A_phi = pm^2 * z/r * sqrt(8 pi wdns_P_max * betam1) * max( P/wdns_P_max - p_c, 0) for Sym_Bz = -1 

!   A_phi = (x^2+y^2) * sqrt(8 pi P_max * betam1) * 
!		max( P/P_max - p_c, 0) 

! where pm^2 = (x-x_ns_com)^2 + (y-y_ns_com)^2 
! x_ns_com is the position of the center of the neutron star.
!------------------------------------------------------------------------
subroutine WDNS_compute_Aphi(ext,X,Y,Z,PhysR,P,A_phi,Ax,Ay,Az, & 
     mhd_psi_line, mhd_u_psi, mhd_chi_line, mhd_u_chi, &
     betam1,p_c,x_NS_CoM_coord,y_NS_CoM_coord,wdns_P_max,Sym_Bz, & 
     enable_trace_field_line,constrained_transport_scheme,em_field_type)
  implicit none
  integer, dimension(3)                          :: ext
  real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z,PhysR
  real*8, dimension(ext(1),ext(2),ext(3))        :: A_phi, P, Ax,Ay,Az
  real*8, dimension(ext(1),ext(2),ext(3))        :: mhd_psi_line, mhd_u_psi, mhd_chi_line, mhd_u_chi
  real*8                                         :: x_NS_CoM_coord,y_NS_CoM_coord
  real*8					 :: betam1,p_c,fac, pomega2
  real*8                                         :: wdns_P_max,xp,yp,zp,r
  real*8 					 :: xp_s,yp_s,zp_s,Aphi_s,P_s,phi
  integer					 :: i,j,k,em_field_type
  integer                                        :: fisheye_enable,ip1,jp1,kp1
  integer                                        :: constrained_transport_scheme
  integer                                        :: enable_trace_field_line
  real*8                                         :: Sym_Bz,P_corner,hdX,hdY,hdZ
  !
  ! fac = sqrt(8.d0*acos(-1.d0)*wdns_P_max*betam1)/Riso**2
  fac = sqrt(8.d0*acos(-1.d0)*wdns_P_max*betam1)
  !n_b = 2.d0/3.d0

  hdX = 0.5d0*(X(2,1,1) - X(1,1,1))
  hdY = 0.5d0*(Y(1,2,1) - Y(1,1,1))
  hdZ = 0.5d0*(Z(1,1,2) - Z(1,1,1))

  if ( abs(sqrt(X(1,1,1)**2 + Y(1,1,1)**2 + Z(1,1,1)**2)-PhysR(1,1,1)) .gt. hdX*1.d-3) then 
     fisheye_enable = 1
  else
     fisheye_enable = 0
  end if
  if (fisheye_enable==1 .and. (em_field_type==2 .or. constrained_transport_scheme==3)) then 
     write(*,*) 'em_field_type=2 and constrained_transport_scheme=3 does not support fisheye (yet).'
     stop
  end if

  do k = 1,ext(3)
     do j = 1,ext(2)
        do i = 1,ext(1)
           r = sqrt(X(i,1,1)**2 + Y(1,j,1)**2 + Z(1,1,k)**2)
           xp = X(i,1,1)/r * PhysR(i,j,k)
           yp = Y(1,j,1)/r * PhysR(i,j,k)
           !!pomega2 = xp**2 + yp**2 + hdX*1.d-13
           pomega2 = (xp-x_NS_CoM_coord)**2 + (yp-y_NS_CoM_coord)**2 + hdX*1.d-13
	   if (constrained_transport_scheme==3) then 
	      ! setup A_i on staggered grid
	      xp_s = X(i,1,1)+hdX 
	      yp_s = Y(1,j,1)+hdY
	      zp_s = Z(1,1,k)+hdZ

	      ip1 = min(i+1,ext(1))
	      jp1 = min(j+1,ext(2))
	      kp1 = min(k+1,ext(3))

	      xp = X(i,1,1)
	      yp = yp_s
	      zp = zp_s
              !!pomega2 = xp**2 + yp**2 + hdX*1.d-13	
              pomega2 = (xp-x_NS_CoM_coord)**2 + (yp-y_NS_CoM_coord)**2 + hdX*1.d-13
	      P_s = 0.25d0*( P(i,j,k) + P(i,jp1,k) + P(i,j,kp1) + P(i,jp1,kp1) )
              !Aphi_s = pomega2 * fac * max(P_s/wdns_P_max - p_c, 0.d0)**1
              Aphi_s = pomega2 * fac * max(P_s/wdns_P_max - p_c, 0.d0)**2
              if (Sym_Bz .lt. 0.d0) Aphi_s = Aphi_s*zp/sqrt(pomega2)
	      !!Ax(i,j,k) = -yp/pomega2*Aphi_s
	      Ax(i,j,k) = -(yp-y_NS_CoM_coord)/pomega2*Aphi_s

	      xp = xp_s
	      yp = Y(1,j,1)
	      zp = zp_s
	      !!pomega2 = xp**2 + yp**2 + hdX*1.d-13
              pomega2 = (xp-x_NS_CoM_coord)**2 + (yp-y_NS_CoM_coord)**2 + hdX*1.d-13
	      P_s = 0.25d0*( P(i,j,k) + P(ip1,j,k) + P(i,j,kp1) + P(ip1,j,kp1) )
              !Aphi_s = pomega2 * fac * max(P_s/wdns_P_max - p_c, 0.d0)**1
              Aphi_s = pomega2 * fac * max(P_s/wdns_P_max - p_c, 0.d0)**2
              if (Sym_Bz .lt. 0.d0) Aphi_s = Aphi_s*zp/sqrt(pomega2)
              !!Ay(i,j,k) = xp/pomega2*Aphi_s
              Ay(i,j,k) = (xp-x_NS_CoM_coord)/pomega2*Aphi_s

	      Az(i,j,k) = 0.d0

	      A_phi(i,j,k) = 0.d0



 	      if (enable_trace_field_line==1) then 
		 pomega2 = (X(i,j,k)-x_NS_CoM_coord)**2 + (Y(i,j,k)-y_NS_CoM_coord)**2 + hdX*1.d-13
		 !Aphi_s = pomega2 * fac * max((P(i,j,k)/wdns_P_max - p_c), 0.d0)**1
                 Aphi_s = pomega2 * fac * max((P(i,j,k)/wdns_P_max - p_c), 0.d0)**2
	         phi = atan2(Y(i,j,k)-y_NS_CoM_coord,X(i,j,k)-x_NS_CoM_coord)
	         mhd_psi_line(i,j,k) = Aphi_s*cos(phi)
		 mhd_chi_line(i,j,k) = Aphi_s*sin(phi)
		 mhd_u_psi(i,j,k) = 0.d0
		 mhd_u_chi(i,j,k) = 0.d0
	      end if

           elseif (em_field_type==0 .or. i==1 .or. j==1 .or. k==1) then 
              !A_phi(i,j,k) = pomega2 * fac * max(P(i,j,k)/wdns_P_max - p_c, 0.d0)**1
              A_phi(i,j,k) = pomega2 * fac * max(P(i,j,k)/wdns_P_max - p_c, 0.d0)**2
              if (enable_trace_field_line==1) then
                 phi = atan2(Y(i,j,k)-y_NS_CoM_coord,X(i,j,k)-x_NS_CoM_coord)
                 mhd_psi_line(i,j,k) = A_phi(i,j,k)*cos(phi)
                 mhd_chi_line(i,j,k) = A_phi(i,j,k)*sin(phi)
                 mhd_u_psi(i,j,k) = 0.d0
                 mhd_u_chi(i,j,k) = 0.d0
              end if
              ! Branson's alternate ways of setting Aphi.  n_b sets the power-law dependence on input parameter P.
              ! Alt way 1: A_phi(i,j,k) = pomega2 * fac * max(P(i,j,k)**(1.d0/n_b) - (wdns_P_max*p_c)**(1.d0/n_b), 0.d0)
              ! Alt way 2: A_phi(i,j,k) = pomega2 * fac * max((P(i,j,k)/(p_c*wdns_P_max))**(1.d0/n_b) - 1.d0, 0.d0)
              if (Sym_Bz .lt. 0.d0) A_phi(i,j,k) = A_phi(i,j,k)*Z(1,1,k)/r*PhysR(i,j,k)/sqrt(pomega2)
           else
	      ! Note that for this initial data thorn, A_phi(i,j,k) 
              ! stores A_phi at the corner point
	      ! x_i - dx/2, y_j-dy/2, z_j-dz/2
	      ! We can only do that for i>1, j>1 and k>1. I assume 
	      ! this will be fine as long as we have enough ghostzones.
              P_corner = 0.125d0*(P(i,j,k) + P(i-1,j,k) +  &
                   P(i,j-1,k) + P(i-1,j-1,k) + &
                   P(i,j,k-1) + P(i-1,j,k-1) + &
                   P(i,j-1,k-1) + P(i-1,j-1,k-1) )
	      xp = X(i,j,k) - hdX
	      yp = Y(i,j,k) - hdY
	      zp = Z(i,j,k) - hdZ
              !!pomega2 = xp**2 + yp**2 + hdX*1.d-13
              pomega2 = (xp-x_NS_CoM_coord)**2 + (yp-y_NS_CoM_coord)**2 + hdX*1.d-13
	      !A_phi(i,j,k) = pomega2 * fac * max(P_corner/wdns_P_max - p_c, 0.d0)**1
              A_phi(i,j,k) = pomega2 * fac * max(P_corner/wdns_P_max - p_c, 0.d0)**2
	      if (Sym_Bz .lt. 0.d0) A_phi(i,j,k) = A_phi(i,j,k)*zp/sqrt(pomega2)
              if (enable_trace_field_line==1) then
                 pomega2 = (X(i,j,k)-x_NS_CoM_coord)**2 + (Y(i,j,k)-y_NS_CoM_coord)**2 + hdX*1.d-13
                 !Aphi_s = pomega2 * fac * max((P(i,j,k)/wdns_P_max - p_c), 0.d0)**1
                 Aphi_s = pomega2 * fac * max((P(i,j,k)/wdns_P_max - p_c), 0.d0)**2
                 phi = atan2(Y(i,j,k)-y_NS_CoM_coord,X(i,j,k)-x_NS_CoM_coord)
                 mhd_psi_line(i,j,k) = Aphi_s*cos(phi)
                 mhd_chi_line(i,j,k) = Aphi_s*sin(phi)
                 mhd_u_psi(i,j,k) = 0.d0
                 mhd_u_chi(i,j,k) = 0.d0
              end if
           end if
        end do
     end do
  end do

end subroutine WDNS_compute_Aphi

! Set up vector potential for toroidal B field according to 
! Ax = ((x-x_ns)z/r0) sqrt[ 8 pi betam1 psi^4 max(P - p_c Pmax , 0) ]
! Ay = ((y-y_ns)z/r0) sqrt[ 8 pi betam1 psi^4 max(P - p_c Pmax , 0) ]
! Az = 0
subroutine WDNS_compute_A_toroidal(ext,X,Y,Z,P,Ax,Ay,Az,phi,betam1, &
     p_c,x_NS_CoM_coord,y_NS_CoM_coord,r0,wdns_P_max) 
   implicit none
   integer, dimension(3) :: ext
   real*8, dimension(ext(1),ext(2),ext(3)) :: X,Y,Z,P,Ax,Ay,Az,phi
   real*8 :: betam1,p_c,x_NS_CoM_coord,y_NS_CoM_coord,wdns_P_max,r0
   real*8 :: zs,xi,yj,hdx,hdy,hdz,psi4s,Ps,Ps_Pcut
   integer :: i,j,k,ip1,jp1,kp1
   real*8, parameter :: f8pi = 25.13274122871834590768d0 ! 8 pi
!
   hdX = 0.5d0*(X(2,1,1) - X(1,1,1))
   hdY = 0.5d0*(Y(1,2,1) - Y(1,1,1))
   hdZ = 0.5d0*(Z(1,1,2) - Z(1,1,1))

   Az = 0.d0
   do k=1,ext(3)
      do j=1,ext(2)
	 do i=1,ext(1)
	    xi = X(i,j,k) - x_NS_CoM_coord
	    yj = Y(i,j,k) - y_NS_CoM_coord
	    zs = Z(i,j,k) + hdz

 	    ip1 = min(i+1,ext(1))
            jp1 = min(j+1,ext(2))
	    kp1 = min(k+1,ext(3))
    	    
            Ps = 0.25d0*(P(i,j,k) + P(i,jp1,k) + P(i,j,kp1) + P(i,jp1,kp1) )
            Ps_Pcut = max(Ps - p_c*wdns_P_max, 0.d0)
	    psi4s = exp( phi(i,j,k) + phi(i,jp1,k) + phi(i,j,kp1) + phi(i,jp1,kp1) )
	    Ax(i,j,k) = xi*zs/r0*sqrt(f8pi*betam1*psi4s*Ps_Pcut)

            Ps = 0.25d0*(P(i,j,k) + P(ip1,j,k) + P(i,j,kp1) + P(ip1,j,kp1) )
            Ps_Pcut = max(Ps - p_c*wdns_P_max, 0.d0)
            psi4s = exp( phi(i,j,k) + phi(ip1,j,k) + phi(i,j,kp1) + phi(ip1,j,kp1) )
	    Ay(i,j,k) = yj*zs/r0*sqrt(f8pi*betam1*psi4s*Ps_Pcut)
	 end do
      end do
   end do
end subroutine WDNS_compute_A_toroidal
