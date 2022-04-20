!------------------------------------------------------------
! Setup the initial EM fields from a vector potential 
!------------------------------------------------------------
!
subroutine setup_poloidal_emfields(ext,Riso,P_max,p_c,betam1,X,Y,Z,PhysR,phi, &
     alpham1,shiftx,shifty,shiftz, & 
     gxx,gxy,gxz,gyy,gyz,gzz,rho_b, P,Bx,By,Bz, Ax,Ay,Az, &
     Bx_stagger, By_stagger, Bz_stagger, &
     st_x,st_y,st_z,mhd_st_x,mhd_st_y,mhd_st_z,tau,u0, &
     vx,vy,vz,Symmetry, Sym_Bz, &
     A_phi,A_phix,A_phiy,A_phiz, Aphi_scaling_factor, & 
     constrained_transport_scheme,em_field_type, u_x,u_y,u_z,psi_n)
  implicit none
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
  integer, dimension(3)				:: ext
  real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z,PhysR
  real*8, dimension(ext(1),ext(2),ext(3))	:: rho_b,P,Bx,By,Bz
  real*8, dimension(ext(1),ext(2),ext(3))        :: mhd_st_x,mhd_st_y,mhd_st_z,tau
  real*8, dimension(ext(1),ext(2),ext(3))        :: st_x,st_y,st_z
  real*8, dimension(ext(1),ext(2),ext(3))        :: u0,vx,vy,vz
  real*8, dimension(ext(1),ext(2),ext(3))        :: u_x,u_y,u_z
  real*8					        :: B_xs, B_ys, B_zs, B2s
  real*8, dimension(ext(1),ext(2),ext(3))        :: phi,psi_n,alpham1
  real*8, dimension(ext(1),ext(2),ext(3))        :: Ax,Ay,Az
  real*8, dimension(ext(1),ext(2),ext(3))        :: Bx_stagger, By_stagger, Bz_stagger
  real*8, dimension(ext(1),ext(2),ext(3))        :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ext(1),ext(2),ext(3))        :: shiftx,shifty,shiftz
  real*8						:: Riso,p_c,betam1,P_max
  real*8                                         :: dX,dY,dZ,al,sqrtg4,sqrtg
  real*8                                         :: Aphi_scaling_factor,psim6_s
  real*8                                         :: sb0,sb2,sb_x,sb_y,sb_z
  integer					:: Symmetry,constrained_transport_scheme

  real*8, dimension(ext(1),ext(2),ext(3))        :: A_phi,A_phix,A_phiy,A_phiz
  !
  integer                                        :: AXISYM,EQUATORIAL
  integer                                        :: OCTANT
  integer                                        :: i,j,k,imin,imax
  integer                                        :: im1,jm1,km1,ip1,jp1,kp1
  integer                                        :: jmin,jmax,kmin,kmax
  integer                                        :: em_field_type
  real*8, parameter 				:: SYM = 1.d0, ANTI = -1.d0
  real*8						:: pomega2,psim6, fs4pi,psin
  real*8                                         :: Sym_Bz
  real*8					 :: Ax000,Ax001,Ax010,Ax011,Ax100,Ax101,Ax110,Ax111
  real*8                                         :: Ay000,Ay001,Ay010,Ay011,Ay100,Ay101,Ay110,Ay111
  parameter(EQUATORIAL = 1, OCTANT = 2, AXISYM = 4)
  !
  fs4pi = sqrt(4.d0*acos(-1.d0))
  imin = 1
  jmin = 1
  kmin = 1
  imax = ext(1)
  jmax = ext(2)
  kmax = ext(3)

  dX = X(imin+1,1,1)-X(imin,1,1)
  dY = Y(1,jmin+1,1)-Y(1,jmin,1)
  dZ = Z(1,1,kmin+1)-Z(1,1,kmin)

  psi_n = exp(4.d0*phi)*u0
  u_x = (gxx*(shiftx+vx) + gxy*(shifty+vy) + gxz*(shiftz+vz))*psi_n
  u_y = (gxy*(shiftx+vx) + gyy*(shifty+vy) + gyz*(shiftz+vz))*psi_n
  u_z = (gxz*(shiftx+vx) + gyz*(shifty+vy) + gzz*(shiftz+vz))*psi_n

  !
  ! Compute the vector potential A_phi 
  !
  call compute_Aphi(ext,X,Y,Z,PhysR,P,A_phi,Ax,Ay,Az,Riso,betam1,p_c, & 
                    P_max,Sym_Bz, & 
                    constrained_transport_scheme,em_field_type)
  ! call compute_Aphi2(ext,X,Y,Z,PhysR,rho_b,A_phi,Riso,betam1,p_c,rho_b_max,Sym_Bz)
  A_phi = A_phi*Aphi_scaling_factor
  ! 
  ! Compute the derivatives of A_phi
  !
  if (constrained_transport_scheme .ne. 3) then 
     if (Symmetry==OCTANT) then
        call gderivs_oct(ext,A_phi,A_phix,A_phiy,A_phiz,dX,dY,dZ,SYM,SYM,Sym_Bz)
     elseif (Symmetry==EQUATORIAL) then
        call gderivs_eq(ext,A_phi,A_phix,A_phiy,A_phiz,dX,dY,dZ,SYM,SYM ,Sym_Bz )
     elseif (Symmetry==AXISYM) then
        call gderivs_axi(ext,A_phi,A_phix,A_phiy,A_phiz,dX,dY,dZ,SYM,SYM,Sym_Bz)
     end if
  end if

  ! Compute B^i on staggered grid 
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
	      Bx_stagger(i,j,k) = ( (Az(i,j,k)-Az(i,jm1,k))/dY   & 
			       - (Ay(i,j,k)-Ay(i,j,km1))/dZ ) * psim6_s

              psim6_s = exp(-3.d0 * (phi(i,j,k) + phi(i,jp1,k)) )
              By_stagger(i,j,k) = ( (Ax(i,j,k)-Ax(i,j,km1))/dZ & 
                               - (Az(i,j,k)-Az(im1,j,k))/dX ) * psim6_s

              psim6_s = exp(-3.d0 * (phi(i,j,k) + phi(i,j,kp1)) )
              Bz_stagger(i,j,k) = ( (Ay(i,j,k)-Ay(im1,j,k))/dX &
                               - (Ax(i,j,k)-Ax(i,jm1,k))/dY ) * psim6_s
	      
	   end do	
	end do
     end do

     ! Now compute B^i on unstaggered grid by simple average
     do k=1,ext(3)
	do j=1,ext(2)
	   do i=1,ext(1)
	      im1 = max(i-1,1)
	      jm1 = max(j-1,1)
	      km1 = max(k-1,1)
	      Bx(i,j,k) = 0.5d0* (Bx_stagger(i,j,k) + Bx_stagger(im1,j,k))
              By(i,j,k) = 0.5d0* (By_stagger(i,j,k) + By_stagger(i,jm1,k))
              Bz(i,j,k) = 0.5d0* (Bz_stagger(i,j,k) + Bz_stagger(i,j,km1))
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
           psim6 = exp(-6.d0*phi(i,j,k))
           pomega2 = X(i,1,1)**2 + Y(1,j,1)**2
	   if (constrained_transport_scheme==3) then 
	      ! do nothing since B^i has been computed 
   	   elseif (em_field_type==0 .or. i==imax .or. j==jmax .or. k==kmax) then 
              Bx(i,j,k) = -X(i,1,1)/pomega2 * psim6 * A_phiz(i,j,k)
              By(i,j,k) = -Y(1,j,1)/pomega2 * psim6 * A_phiz(i,j,k)
              Bz(i,j,k) = psim6/pomega2 * (X(i,1,1)*A_phix(i,j,k) +  & 
                   Y(1,j,1)*A_phiy(i,j,k))
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
	      Ax000 = -A_phi(i,j,k)* (Y(i,j,k)-0.5d0*dY)/ & 
                 ((X(i,j,k)-0.5d0*dX)**2 + (Y(i,j,k)-0.5d0*dY)**2)
	      Ax001 = -A_phi(i,j,k+1)*(Y(i,j,k+1)-0.5d0*dY)/ & 
 		 ((X(i,j,k+1)-0.5d0*dX)**2 + (Y(i,j,k+1)-0.5d0*dY)**2)
              Ax010 = -A_phi(i,j+1,k)* (Y(i,j+1,k)-0.5d0*dY)/ & 
		 ((X(i,j+1,k)-0.5d0*dX)**2 + (Y(i,j+1,k)-0.5d0*dY)**2)
	      Ax011 = -A_phi(i,j+1,k+1)*(Y(i,j+1,k+1)-0.5d0*dY)/ & 
		 ((X(i,j+1,k+1)-0.5d0*dX)**2 + (Y(i,j+1,k+1)-0.5d0*dY)**2)
	      Ax100 = -A_phi(i+1,j,k)*(Y(i+1,j,k)-0.5d0*dY)/ & 
     		 ((X(i+1,j,k)-0.5d0*dX)**2 + (Y(i+1,j,k)-0.5d0*dY)**2)
	      Ax101 = -A_phi(i+1,j,k+1)*(Y(i+1,j,k+1)-0.5d0*dY)/ & 
		 ((X(i+1,j,k+1)-0.5d0*dX)**2 + (Y(i+1,j,k+1)-0.5d0*dY)**2)
	      Ax110 = -A_phi(i+1,j+1,k)*(Y(i+1,j+1,k)-0.5d0*dY)/ & 
		 ((X(i+1,j+1,k)-0.5d0*dX)**2 + (Y(i+1,j+1,k)-0.5d0*dY)**2)
	      Ax111 = -A_phi(i+1,j+1,k+1)*(Y(i+1,j+1,k+1)-0.5d0*dY)/ & 
		 ((X(i+1,j+1,k+1)-0.5d0*dX)**2 + (Y(i+1,j+1,k+1)-0.5d0*dY)**2)
              Ay000 = A_phi(i,j,k)* (X(i,j,k)-0.5d0*dX)/ &
                  ((X(i,j,k)-0.5d0*dX)**2 + (Y(i,j,k)-0.5d0*dY)**2)
              Ay001 = A_phi(i,j,k+1)*(X(i,j,k+1)-0.5d0*dX)/ &
                  ((X(i,j,k+1)-0.5d0*dX)**2 + (Y(i,j,k+1)-0.5d0*dY)**2)
              Ay010 = A_phi(i,j+1,k)* (X(i,j+1,k)-0.5d0*dX)/ &
                  ((X(i,j+1,k)-0.5d0*dX)**2 + (Y(i,j+1,k)-0.5d0*dY)**2)
              Ay011 = A_phi(i,j+1,k+1)*(X(i,j+1,k+1)-0.5d0*dX)/ &
                  ((X(i,j+1,k+1)-0.5d0*dX)**2 + (Y(i,j+1,k+1)-0.5d0*dY)**2)
              Ay100 = A_phi(i+1,j,k)*(X(i+1,j,k)-0.5d0*dX)/ &
                  ((X(i+1,j,k)-0.5d0*dX)**2 + (Y(i+1,j,k)-0.5d0*dY)**2)
              Ay101 = A_phi(i+1,j,k+1)*(X(i+1,j,k+1)-0.5d0*dX)/ &
                  ((X(i+1,j,k+1)-0.5d0*dX)**2 + (Y(i+1,j,k+1)-0.5d0*dY)**2)
              Ay110 = A_phi(i+1,j+1,k)*(X(i+1,j+1,k)-0.5d0*dX)/ &
                  ((X(i+1,j+1,k)-0.5d0*dX)**2 + (Y(i+1,j+1,k)-0.5d0*dY)**2)
              Ay111 = A_phi(i+1,j+1,k+1)*(X(i+1,j+1,k+1)-0.5d0*dX)/ &
                  ((X(i+1,j+1,k+1)-0.5d0*dX)**2 + (Y(i+1,j+1,k+1)-0.5d0*dY)**2)

	      Bx(i,j,k) = -( (Ay001-Ay000) + (Ay101-Ay100) + (Ay011-Ay010) + (Ay111-Ay110) ) * 0.25d0/dZ * psim6
	      By(i,j,k) = ( (Ax001-Ax000) + (Ax101-Ax100) + (Ax011-Ax010) + (Ax111-Ax110) ) * 0.25d0/dZ * psim6
	      Bz(i,j,k) = ( (Ay100-Ay000) + (Ay110-Ay010) + (Ay101-Ay001) + (Ay111-Ay011) ) * 0.25d0/dX * psim6 - ( (Ax010-Ax000) + (Ax110-Ax100) + (Ax011-Ax001) + (Ax111-Ax101) ) * 0.25d0/dY * psim6
	   end if
           ! Compute b^0 and b_i
           al = 1.d0 + alpham1(i,j,k)
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
           sb0 = (u_x(i,j,k)*Bx(i,j,k) + u_y(i,j,k)*By(i,j,k) + &
                u_z(i,j,k)*Bz(i,j,k))/fs4pi/al
           sb2 = (B2s + sb0**2)/u0(i,j,k)**2
           !Branson's way of ensuring b^2/P is what we want: 
           !  Call this routine once to calibrate sb2/P, then average it.
           !   Scale sb2 accordingly, then call this function again to obtain the desired sb2/P
           ! Here's the necessary line of code.  You'll need to add the input parameters as well.
           !bsq(i,j,k) = sb2/P(i,j,k)
           sb_x = (B_xs + u_x(i,j,k)*sb0)/u0(i,j,k)
           sb_y = (B_ys + u_y(i,j,k)*sb0)/u0(i,j,k)
           sb_z = (B_zs + u_z(i,j,k)*sb0)/u0(i,j,k)
           ! Now compute mhd_st_i and tau
           mhd_st_x(i,j,k) = st_x(i,j,k)+sqrtg4*(u0(i,j,k)*sb2*u_x(i,j,k)-sb0*sb_x)
           mhd_st_y(i,j,k) = st_y(i,j,k)+sqrtg4*(u0(i,j,k)*sb2*u_y(i,j,k)-sb0*sb_y)
           mhd_st_z(i,j,k) = st_z(i,j,k)+sqrtg4*(u0(i,j,k)*sb2*u_z(i,j,k)-sb0*sb_z)
           tau(i,j,k) = tau(i,j,k) + sqrtg*( sb2*(al*u0(i,j,k))**2 &
                - sb2*0.5d0 - (al*sb0)**2 )

        end do
     end do
  end do
end subroutine setup_poloidal_emfields

!------------------------------------------------------------
! Compute the vector potential A_phi: 
!   A_phi = (x^2+y^2) * sqrt(8 pi P_max * betam1) * 
!		max( P/P_max - p_c, 0) 
!------------------------------------------------------------
subroutine compute_Aphi(ext,X,Y,Z,PhysR,P,A_phi,Ax,Ay,Az,Riso,betam1, & 
                         p_c,P_max,Sym_Bz, & 
                         constrained_transport_scheme,em_field_type)
  implicit none
  integer, dimension(3)                          :: ext
  real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z,PhysR
  real*8, dimension(ext(1),ext(2),ext(3))        :: A_phi, P, Ax,Ay,Az
  real*8					:: betam1,p_c,fac, pomega2
  real*8                                         :: Riso,P_max,xp,yp,zp,r
  real*8 					:: xp_s,yp_s,zp_s,Aphi_s,P_s
  integer					:: i,j,k,em_field_type
  integer                                       :: fisheye_enable,ip1,jp1,kp1
  integer                                       :: constrained_transport_scheme
  real*8                                         :: Sym_Bz,P_corner,hdX,hdY,hdZ
  !
  ! fac = sqrt(8.d0*acos(-1.d0)*P_max*betam1)/Riso**2
  fac = sqrt(8.d0*acos(-1.d0)*P_max*betam1)
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
           pomega2 = xp**2 + yp**2 + hdX*1.d-13
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
	      pomega2 = xp**2 + yp**2 + hdX*1.d-13	
	      P_s = 0.25d0*( P(i,j,k) + P(i,jp1,k) + P(i,j,kp1) + P(i,jp1,kp1) )
	      Aphi_s = pomega2 * fac * max(P_s/P_max - p_c, 0.d0)
              if (Sym_Bz .lt. 0.d0) Aphi_s = Aphi_s*zp/sqrt(pomega2)
	      Ax(i,j,k) = -yp/pomega2*Aphi_s

	      xp = xp_s
	      yp = Y(1,j,1)
	      zp = zp_s
	      pomega2 = xp**2 + yp**2 + hdX*1.d-13
	      P_s = 0.25d0*( P(i,j,k) + P(ip1,j,k) + P(i,j,kp1) + P(ip1,j,kp1) )
              Aphi_s = pomega2 * fac * max(P_s/P_max - p_c, 0.d0)
              if (Sym_Bz .lt. 0.d0) Aphi_s = Aphi_s*zp/sqrt(pomega2)
              Ay(i,j,k) = xp/pomega2*Aphi_s

	      Az(i,j,k) = 0.d0

	      A_phi(i,j,k) = 0.d0
 
           elseif (em_field_type==0 .or. i==1 .or. j==1 .or. k==1) then 
              A_phi(i,j,k) = pomega2 * fac * max(P(i,j,k)/P_max - p_c, 0.d0)
              ! Branson's alternate ways of setting Aphi.  n_b sets the power-law dependence on input parameter P.
              ! Alt way 1: A_phi(i,j,k) = pomega2 * fac * max(P(i,j,k)**(1.d0/n_b) - (P_max*p_c)**(1.d0/n_b), 0.d0)
              ! Alt way 2: A_phi(i,j,k) = pomega2 * fac * max((P(i,j,k)/(p_c*P_max))**(1.d0/n_b) - 1.d0, 0.d0)
              if (Sym_Bz .lt. 0.d0) A_phi(i,j,k) = A_phi(i,j,k)*Z(1,1,k)/r*PhysR(i,j,k)/sqrt(pomega2)
           else
	      ! Note that for this initial data, A_phi(i,j,k) 
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
              pomega2 = xp**2 + yp**2 + hdX*1.d-13
	      A_phi(i,j,k) = pomega2 * fac * max(P_corner/P_max - p_c, 0.d0)
	      if (Sym_Bz .lt. 0.d0) A_phi(i,j,k) = A_phi(i,j,k)*zp/sqrt(pomega2)
           end if
        end do
     end do
  end do

  ! write(*,*) "hello. inside computeaphi: ",A_phi(30,30,2),P(30,30,2),P_max,p_c,betam1,fac,1.D0,Sym_Bz

end subroutine compute_Aphi
