!------------------------------------------------------------
! Setup an initial toroidal EM field 
!   B^{\phi} = A_b sqrt( max( P/Pmax - p_c, 0) )
!------------------------------------------------------------
!
subroutine initial_toroidal_emfields(ext,P_max,p_c,Ab,X,Y,Z,phi, &
		alpham1,shiftx,shifty,shiftz, & 
		gxx,gxy,gxz,gyy,gyz,gzz,rho_b, P,Bx,By,Bz, &
		st_x,st_y,st_z,mhd_st_x,mhd_st_y,mhd_st_z,tau,u0, &
                vx,vy,vz,Symmetry, Sym_Bz)
 implicit none
 integer, dimension(3)				:: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8, dimension(ext(1),ext(2),ext(3))	:: rho_b,P,Bx,By,Bz
 real*8, dimension(ext(1),ext(2),ext(3))        :: mhd_st_x,mhd_st_y,mhd_st_z,tau
 real*8, dimension(ext(1),ext(2),ext(3))        :: st_x,st_y,st_z
 real*8, dimension(ext(1),ext(2),ext(3))        :: u0,vx,vy,vz
 real*8, dimension(ext(1),ext(2),ext(3))        :: u_x,u_y,u_z
 real*8					        :: B_xs, B_ys, B_zs, B2s
 real*8, dimension(ext(1),ext(2),ext(3))        :: phi,psi_n,alpham1
 real*8, dimension(ext(1),ext(2),ext(3))        :: gxx,gxy,gxz,gyy,gyz,gzz
 real*8, dimension(ext(1),ext(2),ext(3))        :: shiftx,shifty,shiftz
 real*8						:: Riso,p_c,Ab,P_max
 real*8                                         :: dX,dY,dZ,al,sqrtg4,sqrtg
 real*8                                         :: sb0,sb2,sb_x,sb_y,sb_z
 integer					:: Symmetry
!
 integer                                        :: AXISYM,EQUATORIAL
 integer                                        :: OCTANT
 integer                                        :: i,j,k,imin,imax
 integer                                        :: jmin,jmax,kmin,kmax
 real*8, parameter 				:: SYM = 1.d0, ANTI = -1.d0
 real*8						:: psim6, fs4pi,psin
 real*8                                         :: Sym_Bz, Bphi
 parameter(EQUATORIAL = 1, OCTANT = 2, AXISYM = 4)
!
! if (Sym_Bz .gt. 0.d0) then
!    write(*,*) 'Sym_Bz = ',Sym_Bz
!    write(*,*) 'Sym_Bz must be negative here' 
!    stop
! end if


 fs4pi = sqrt(4.d0*acos(-1.d0))
 imin = lbound(P,1)
 imax = ubound(P,1)
 jmin = lbound(P,2)
 jmax = ubound(P,2)
 kmin = lbound(P,3)
 kmax = ubound(P,3)

 dX = X(imin+1,1,1)-X(imin,1,1)
 dY = Y(1,jmin+1,1)-Y(1,jmin,1)
 dZ = Z(1,1,kmin+1)-Z(1,1,kmin)

 psi_n = exp(4.d0*phi)*u0
 u_x = (gxx*(shiftx+vx) + gxy*(shifty+vy) + gxz*(shiftz+vz))*psi_n
 u_y = (gxy*(shiftx+vx) + gyy*(shifty+vy) + gyz*(shiftz+vz))*psi_n
 u_z = (gxz*(shiftx+vx) + gyz*(shifty+vy) + gzz*(shiftz+vz))*psi_n

!
! Compute B^i according to
!  B^x = -y B^{\phi}
!  B^y =  x B^{\phi}
!  B^z = 0
!
! and then calculate mhd_st_i and tau
!
 do i = imin,imax
    do j = jmin,jmax
       do k=kmin,kmax
	  psim6 = exp(-6.d0*phi(i,j,k))
	  Bphi = Ab * sqrt( max(P(i,j,k)/P_max - p_c ,0.d0) )
	  if (Sym_Bz .gt. 0.d0) Bphi = Z(1,1,k)*Bphi
          Bx(i,j,k) = -Y(1,j,1) * Bphi
	  By(i,j,k) = X(i,1,1) * Bphi
	  Bz(i,j,k) = 0.d0
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
          sb_x = (B_xs + u_x(i,j,k)*sb0)/u0(i,j,k)
          sb_y = (B_ys + u_y(i,j,k)*sb0)/u0(i,j,k)
          sb_z = (B_zs + u_z(i,j,k)*sb0)/u0(i,j,k)
! if (j==2) then
!    E_em = E_em + 4.d0*acos(-1.d0)*X(i,1,1)*sqrtg* & 
!		( ( (al*u0(i,j,k))**2 - 0.5d0)*sb2 - (al*sb0)**2 ) * dX*dZ
!    if (k==kmax) write(11,*) X(i,1,1),E_em
! end if 
! ************
! Now compute mhd_st_i and tau
	  mhd_st_x(i,j,k) = st_x(i,j,k)+sqrtg4*(u0(i,j,k)*sb2*u_x(i,j,k)-sb0*sb_x)
          mhd_st_y(i,j,k) = st_y(i,j,k)+sqrtg4*(u0(i,j,k)*sb2*u_y(i,j,k)-sb0*sb_y)
          mhd_st_z(i,j,k) = st_z(i,j,k)+sqrtg4*(u0(i,j,k)*sb2*u_z(i,j,k)-sb0*sb_z)
          tau(i,j,k) = tau(i,j,k) + sqrtg*( sb2*(al*u0(i,j,k))**2 &
                           - sb2*0.5d0 - (al*sb0)**2 )
       end do
    end do
 end do
end subroutine initial_toroidal_emfields
