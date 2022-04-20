!-----------------------------------------------------------------------------
!
! $Id
!
!-----------------------------------------------------------------------------
!
! Integrate Omega S_phi to get T
!
!-----------------------------------------------------------------------------
subroutine tint(ex, dT, &
     X, Y, Z, &
     phi, v_x, v_y, st_x, st_y, Symmetry, &
     adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax)
  implicit none
! Input parameters:
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: phi,v_x,v_y,st_x,st_y
  integer                                     :: adjimin, adjjmin, adjkmin
  integer                                     :: adjimax, adjjmax, adjkmax
  integer                                     :: Symmetry
! output:
  real*8                                      :: dT
! Other variables:
  real*8                             :: dX, dY, dZ
  real*8                             :: dV,R,cs,sn,pi,FOUR,ONE
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin
  integer                            :: imax, jmax, kmax
  integer                            :: AXISYM
  parameter(AXISYM = 4)
  parameter(ONE = 1.D0, FOUR = 4.D0)
  pi = acos( - one )
!
! Coordinate grid size, putting adjustments on such that
! it excludes the ghost zones.
!
  imin = lbound(phi,1) - adjimin
  jmin = lbound(phi,2) - adjjmin
  kmin = lbound(phi,3) - adjkmin
  imax = ubound(phi,1) - adjimax
  jmax = ubound(phi,2) - adjjmax
  kmax = ubound(phi,3) - adjkmax
!
  dX = X(imin+1,1,1)-X(imin,1,1)
  dY = Y(1,jmin+1,1)-Y(1,jmin,1)
  dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
  dV = dX * dY * dZ
  if (Symmetry==AXISYM) then
     jmin=2
     jmax=2
     dV=dX*dZ*FOUR*PI
  end if
!
!-----------------------------------------------------------------------------
! integrate
!-----------------------------------------------------------------------------
  dT = 0.D0
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
           if (Symmetry == AXISYM) then
              dT = dT + 0.5D0 * v_y(i,j,k)  &
                * st_y(i,j,k) * abs(X(i,1,1))
           else
              R = sqrt(X(i,1,1)*X(i,1,1) + Y(1,j,1)*Y(1,j,1))
              cs = X(i,1,1)/R
              sn = Y(1,j,1)/R
              dT = dT + 0.5D0 * (cs*v_y(i,j,k) - sn*v_x(i,j,k)) &
                * (cs*st_y(i,j,k) - sn*st_x(i,j,k)) 
           end if
        end do ! i-loop
     end do ! j-loop
  end do ! k-loop
  dT = dT * dV
  return
end subroutine tint
!-----------------------------------------------------------------------------
!
! Integrate \rho_star \epsilon to get the total internal energy
!
!-----------------------------------------------------------------------------
subroutine minternal(ex, dT, &
     X, Y, Z, &
     rho_star, h, n_poly, Symmetry, &
     adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax)
  implicit none
! Input parameters:
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: rho_star,h
  integer                                     :: adjimin, adjjmin, adjkmin
  integer                                     :: adjimax, adjjmax, adjkmax
  integer                                     :: Symmetry
  real*8                                      :: n_poly
! output:
  real*8                                      :: dT
! Other variables:
  real*8                             :: dX, dY, dZ, Gamma
  real*8                             :: dV,R,cs,sn,pi,FOUR,ONE
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin
  integer                            :: imax, jmax, kmax
  integer                            :: AXISYM
  parameter(AXISYM = 4)
  parameter(ONE = 1.D0, FOUR = 4.D0)
  pi = acos( - one )
  Gamma   = ONE + ONE/n_poly
!
! Coordinate grid size, putting adjustments on such that
! it excludes the ghost zones.
!
  imin = lbound(rho_star,1) - adjimin
  jmin = lbound(rho_star,2) - adjjmin
  kmin = lbound(rho_star,3) - adjkmin
  imax = ubound(rho_star,1) - adjimax
  jmax = ubound(rho_star,2) - adjjmax
  kmax = ubound(rho_star,3) - adjkmax
!
  dX = X(imin+1,1,1)-X(imin,1,1)
  dY = Y(1,jmin+1,1)-Y(1,jmin,1)
  dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
  dV = dX * dY * dZ
  if (Symmetry==AXISYM) then
     jmin=2
     jmax=2
     dV=dX*dZ*FOUR*PI
  end if
!
!-----------------------------------------------------------------------------
! integrate
!-----------------------------------------------------------------------------
  dT = 0.D0
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
           if (Symmetry == AXISYM) then
              dT = dT + rho_star(i,j,k)*((h(i,j,k)-ONE)/Gamma) &
                   * abs(X(i,1,1))
           else
              dT = dT + rho_star(i,j,k)*((h(i,j,k)-ONE)/Gamma)
           end if
        end do ! i-loop
     end do ! j-loop
  end do ! k-loop
  dT = dT * dV
  return
end subroutine minternal

!-----------------------------------------------------------------------------
!
! Integrate \rho_star \epsilon to get the total internal energy (for the
!   hybrid EOS)
!
!-----------------------------------------------------------------------------
subroutine minternal_hybrid(ex, dT, X, Y, Z, rho_tiny, &
     rho_star,rho_b,P,neos,ergo_star, ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th, &
     Symmetry, adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax)
  implicit none
  interface
    subroutine compute_pcold_epscold(rhob, P_cold, eps_cold, &
                neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
     implicit none
     integer :: neos,ergo_star
     real*8  :: rhob, P_cold, eps_cold,ergo_sigma
     real*8, dimension(neos) :: rho_tab,P_tab,eps_tab
     real*8, dimension(neos+1) :: k_tab, gamma_tab
    end subroutine compute_pcold_epscold
  end interface
! Input parameters:
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: rho_star,P,rho_b
  integer                                     :: adjimin, adjjmin, adjkmin
  integer                                     :: adjimax, adjjmax, adjkmax
  integer                                     :: Symmetry,neos,ergo_star
  real*8, dimension(neos)                     :: rho_tab,P_tab,eps_tab
  real*8, dimension(neos+1)                   :: k_tab, gamma_tab
! output:
  real*8                                      :: dT
! Other variables:
  real*8                             :: dX, dY, dZ
  real*8                             :: dV,R,cs,sn,pi,FOUR,ONE, rho_tiny
  real*8                             :: rhob, P_cold, eps_cold, eps,gamma_th,ergo_sigma
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin
  integer                            :: imax, jmax, kmax
  integer                            :: AXISYM
  parameter(AXISYM = 4)
  parameter(ONE = 1.D0, FOUR = 4.D0)
  pi = acos( - one )
!
! Coordinate grid size, putting adjustments on such that
! it excludes the ghost zones.
!
  imin = lbound(rho_star,1) - adjimin
  jmin = lbound(rho_star,2) - adjjmin
  kmin = lbound(rho_star,3) - adjkmin
  imax = ubound(rho_star,1) - adjimax
  jmax = ubound(rho_star,2) - adjjmax
  kmax = ubound(rho_star,3) - adjkmax
!
  dX = X(imin+1,1,1)-X(imin,1,1)
  dY = Y(1,jmin+1,1)-Y(1,jmin,1)
  dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
  dV = dX * dY * dZ
  if (Symmetry==AXISYM) then
     jmin=2
     jmax=2
     dV=dX*dZ*FOUR*PI
  end if
!
!-----------------------------------------------------------------------------
! integrate
!-----------------------------------------------------------------------------
  dT = 0.D0
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax

           rhob = rho_b(i,j,k)
           if (rhob .gt. rho_tiny) then
              call compute_pcold_epscold(rhob, P_cold, eps_cold, &
                     neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
              eps = eps_cold + (P(i,j,k)-P_cold)/(gamma_th-1.d0)/rhob
           else
              eps = 0.d0
           end if

           if (Symmetry == AXISYM) then
              dT = dT + rho_star(i,j,k)*eps* abs(X(i,1,1))
           else
              dT = dT + rho_star(i,j,k)*eps
           end if
        end do ! i-loop
     end do ! j-loop
  end do ! k-loop
  dT = dT * dV
  return
end subroutine minternal_hybrid

!-----------------------------------------------------------------------------
!
! Integrate \rho_0\epsilon_ad to get the adiabatic part of the internal energy
!
!-----------------------------------------------------------------------------
subroutine minternal_ad(ex, dT, &
     X, Y, Z, &
     rho_star, rho_b, rho_tiny, n_poly, Symmetry, &
     adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax)
  implicit none
! Input parameters:
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: rho_star,rho_b
  integer                                     :: adjimin, adjjmin, adjkmin
  integer                                     :: adjimax, adjjmax, adjkmax
  integer                                     :: Symmetry
  real*8                                      :: n_poly, rho_tiny
! output:
  real*8                                      :: dT
! Other variables:
  real*8                             :: dX, dY, dZ, gam1, dV
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin
  integer                            :: imax, jmax, kmax
  integer                            :: AXISYM
  parameter(AXISYM = 4)
  gam1   = 1.d0/n_poly
!
! Coordinate grid size, putting adjustments on such that
! it excludes the ghost zones.
!
  imin = lbound(rho_star,1) - adjimin
  jmin = lbound(rho_star,2) - adjjmin
  kmin = lbound(rho_star,3) - adjkmin
  imax = ubound(rho_star,1) - adjimax
  jmax = ubound(rho_star,2) - adjjmax
  kmax = ubound(rho_star,3) - adjkmax
!
  dX = X(imin+1,1,1)-X(imin,1,1)
  dY = Y(1,jmin+1,1)-Y(1,jmin,1)
  dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
  dV = dX * dY * dZ
  if (Symmetry==AXISYM) then
     jmin=2
     jmax=2
     dV=dX*dZ*4.d0*acos(-1.d0)
  end if
!
!-----------------------------------------------------------------------------
! integrate
!-----------------------------------------------------------------------------
  dT = 0.D0
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
	   if (rho_b(i,j,k) .gt. 0.d0) then 
              if (Symmetry == AXISYM) then
                 dT = dT + rho_star(i,j,k)/gam1* &
			X(i,1,1) * rho_b(i,j,k)**gam1 
              else
                 dT = dT + rho_star(i,j,k)/gam1* &
			rho_b(i,j,k)**gam1
              end if
	    end if
        end do ! i-loop
     end do ! j-loop
  end do ! k-loop
  dT = dT * dV
end subroutine minternal_ad

!-----------------------------------------------------------------------------
!
! Integrate \rho_*\epsilon_ad to get the adiabatic part of the internal energy
!   (for hybrid EOS)
!
!-----------------------------------------------------------------------------
subroutine minternal_ad_hybrid(ex, dT, rho_tiny, X, Y, Z, &
     rho_star,rho_b,neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
     Symmetry,adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax)
  implicit none
  interface
    subroutine compute_pcold_epscold(rhob, P_cold, eps_cold, &
                neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
     implicit none
     integer :: neos,ergo_star
     real*8  :: rhob, P_cold, eps_cold,ergo_sigma
     real*8, dimension(neos) :: rho_tab,P_tab,eps_tab
     real*8, dimension(neos+1) :: k_tab, gamma_tab
    end subroutine compute_pcold_epscold
  end interface
! Input parameters:
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: rho_star,rho_b
  integer                                     :: adjimin, adjjmin, adjkmin
  integer                                     :: adjimax, adjjmax, adjkmax
  integer                                     :: Symmetry,neos,ergo_star
  real*8, dimension(neos)                     :: rho_tab,P_tab,eps_tab
  real*8, dimension(neos+1)                   :: k_tab, gamma_tab
! output:
  real*8                                      :: dT
! Other variables:
  real*8                             :: dX, dY, dZ, dV
  real*8                             :: P_cold, eps_cold, rho_tiny,ergo_sigma
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin
  integer                            :: imax, jmax, kmax
  integer                            :: AXISYM
  parameter(AXISYM = 4)
!
! Coordinate grid size, putting adjustments on such that
! it excludes the ghost zones.
!
  imin = lbound(rho_star,1) - adjimin
  jmin = lbound(rho_star,2) - adjjmin
  kmin = lbound(rho_star,3) - adjkmin
  imax = ubound(rho_star,1) - adjimax
  jmax = ubound(rho_star,2) - adjjmax
  kmax = ubound(rho_star,3) - adjkmax
!
  dX = X(imin+1,1,1)-X(imin,1,1)
  dY = Y(1,jmin+1,1)-Y(1,jmin,1)
  dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
  dV = dX * dY * dZ
  if (Symmetry==AXISYM) then
     jmin=2
     jmax=2
     dV=dX*dZ*4.d0*acos(-1.d0)
  end if
!
!-----------------------------------------------------------------------------
! integrate
!-----------------------------------------------------------------------------
  dT = 0.D0
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
           if (rho_b(i,j,k) .gt. rho_tiny) then
              call compute_pcold_epscold(rho_b(i,j,k), P_cold, eps_cold, &
                     neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
              if (Symmetry == AXISYM) then
                 dT = dT + rho_star(i,j,k)* X(i,1,1) * eps_cold
              else
                 dT = dT + rho_star(i,j,k) * eps_cold
              end if
            end if
        end do ! i-loop
     end do ! j-loop
  end do ! k-loop
  dT = dT * dV
end subroutine minternal_ad_hybrid

!-----------------------------------------------------------------------------
!
! Integrate R*rho_star to get <R>
!
!-----------------------------------------------------------------------------
subroutine mrint(ex, dR, &
     X, Y, Z, rho_b, Symmetry, &
     adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax)
  implicit none
! Input parameters:
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: rho_b
  integer                                     :: adjimin, adjjmin, adjkmin
  integer                                     :: adjimax, adjjmax, adjkmax
  integer                                     :: Symmetry
! output:
  real*8                                      :: dR
! Other variables:
  real*8                             :: dX, dY, dZ
  real*8                             :: dV,R
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin
  integer                            :: imax, jmax, kmax
  integer, parameter		     :: AXISYM=4
!
! Coordinate grid size, putting adjustments on such that
! it excludes the ghost zones.
!
  imin = lbound(rho_b,1) - adjimin
  jmin = lbound(rho_b,2) - adjjmin
  kmin = lbound(rho_b,3) - adjkmin
  imax = ubound(rho_b,1) - adjimax
  jmax = ubound(rho_b,2) - adjjmax
  kmax = ubound(rho_b,3) - adjkmax
!
  dX = X(imin+1,1,1)-X(imin,1,1)
  dY = Y(1,jmin+1,1)-Y(1,jmin,1)
  dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
  dV = dX * dY * dZ
  if (Symmetry==AXISYM) then
     dV=dX*dZ*4.d0*acos(-1.d0)
     jmin=2
     jmax=2
  end if
!
!-----------------------------------------------------------------------------
! integrate
!-----------------------------------------------------------------------------
  dR = 0.D0
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
           R = sqrt(X(i,1,1)*X(i,1,1) + Y(1,j,1)*Y(1,j,1) + Z(1,1,k)*Z(1,1,k))
	   if (Symmetry==AXISYM) then
	      dR = dR + dV  * R * rho_b(i,j,k) * abs(X(i,1,1))
	   else
              dR = dR + dV  * R * rho_b(i,j,k)
	   end if
        end do ! i-loop
     end do ! j-loop
  end do ! k-loop
  return
end subroutine mrint

!-----------------------------------------------------------------------------
!
! Compute Ixx
!
!-----------------------------------------------------------------------------
subroutine ixx(ex, dT, &
     X, Y, Z, rho_star, &
     adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax)
  implicit none
! Input parameters:
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: rho_star
  integer                                     :: adjimin, adjjmin, adjkmin
  integer                                     :: adjimax, adjjmax, adjkmax
! output:
  real*8                                      :: dT
! Other variables:
  real*8                             :: dX, dY, dZ
  real*8                             :: dV,R,cs,sn
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin
  integer                            :: imax, jmax, kmax
  imin = lbound(rho_star,1) - adjimin
  jmin = lbound(rho_star,2) - adjjmin
  kmin = lbound(rho_star,3) - adjkmin
  imax = ubound(rho_star,1) - adjimax
  jmax = ubound(rho_star,2) - adjjmax
  kmax = ubound(rho_star,3) - adjkmax
  dX = X(imin+1,1,1)-X(imin,1,1)
  dY = Y(1,jmin+1,1)-Y(1,jmin,1)
  dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
  dV = dX * dY * dZ
  dT = 0.D0
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
           dT = dT + dV  * rho_star(i,j,k) * X(i,1,1) * X(i,1,1)
        end do ! i-loop
     end do ! j-loop
  end do ! k-loop
  return
end subroutine ixx
!-----------------------------------------------------------------------------
!
! Compute Iyy
!
!-----------------------------------------------------------------------------
subroutine iyy(ex, dT, &
     X, Y, Z, rho_star, &
     adjimin, adjimax, adjjmin, adjjmax, adjkmin, adjkmax)
  implicit none
! Input parameters:
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: rho_star
  integer                                     :: adjimin, adjjmin, adjkmin
  integer                                     :: adjimax, adjjmax, adjkmax
! output:
  real*8                                      :: dT
! Other variables:
  real*8                             :: dX, dY, dZ
  real*8                             :: dV,R,cs,sn
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin
  integer                            :: imax, jmax, kmax
  imin = lbound(rho_star,1) - adjimin
  jmin = lbound(rho_star,2) - adjjmin
  kmin = lbound(rho_star,3) - adjkmin
  imax = ubound(rho_star,1) - adjimax
  jmax = ubound(rho_star,2) - adjjmax
  kmax = ubound(rho_star,3) - adjkmax
  dX = X(imin+1,1,1)-X(imin,1,1)
  dY = Y(1,jmin+1,1)-Y(1,jmin,1)
  dZ = Z(1,1,kmin+1)-Z(1,1,kmin)
  dV = dX * dY * dZ
  dT = 0.D0
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
           dT = dT + dV  * rho_star(i,j,k) * Y(1,j,1) * Y(1,j,1)
        end do ! i-loop
     end do ! j-loop
  end do ! k-loop
  return
end subroutine iyy
!-----------------------------------------------------------------------------
!
! Compute T_{0 r}
!
!-----------------------------------------------------------------------------
subroutine find_tr0(ex, X, Y, Z, rho_b, h, P, rho_star, &
     st_x, st_y, st_z, &
     lapse, shiftx, shifty, shiftz, phi, &
     gxx, gxy, gxz, gyy, gyz, gzz, &
     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Tr0)
  implicit none
! Input parameters:
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: rho_b,h,P,phi,Tr0
  real*8, dimension(ex(1),ex(2),ex(3))        :: st_x,st_y,st_z,rho_star
  real*8, dimension(ex(1),ex(2),ex(3))        :: lapse,shiftx,shifty,shiftz
  real*8, dimension(ex(1),ex(2),ex(3))        :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ex(1),ex(2),ex(3))        :: gupxx, gupxy, gupxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: gupyy, gupyz, gupzz
! Other variables:
  real*8, dimension(ex(1),ex(2),ex(3))        :: ut,ux,uy,uz,al,psim4
  real*8, dimension(ex(1),ex(2),ex(3))        :: bx,by,bz
  real*8                             :: r, Tiny
  integer                            :: i, j, k
  integer                            :: imin, jmin, kmin
  integer                            :: imax, jmax, kmax
  real*8, parameter                  :: ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0
!
  Tiny = 1.d-8
  imin = lbound(phi,1)
  jmin = lbound(phi,2)
  kmin = lbound(phi,3)
  imax = ubound(phi,1)
  jmax = ubound(phi,2)
  kmax = ubound(phi,3)

  al = lapse + ONE
  psim4 = exp(4.D0*phi)
  psim4 = ONE / psim4
  where (rho_b > Tiny)
     ! Calculate h with Gamma-law EOS

     ux = st_x / (rho_star*h)
     uy = st_y / (rho_star*h)
     uz = st_z / (rho_star*h)
  
     ! Calculate u_t using u_{\alpha}u^{\alpha} = -1
     ut = -al*sqrt(ONE + psim4*(gupxx*ux*ux + TWO*gupxy*ux*uy + &
          TWO*gupxz*ux*uz + gupyy*uy*uy + TWO*gupyz*uy*uz + gupzz*uz*uz))
     ut = ut + shiftx*ux + shifty*uy + shiftz*uz
  elsewhere 
     ux = ZERO
     uy = ZERO
     uz = ZERO
     ut = -al
  end where

! calculate \beta_i
  psim4 = ONE / psim4
  bx = psim4 * (shiftx*gxx + shifty*gxy + shiftz*gxz)
  by = psim4 * (shiftx*gxy + shifty*gyy + shiftz*gyz)
  bz = psim4 * (shiftx*gxz + shifty*gyz + shiftz*gzz)

! compute u_r and g_{0r}
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
           r = sqrt(X(i,1,1)**2 + Y(1,j,1)**2 + Z(1,1,k)**2)
! store u_r in ux
           ux(i,j,k) = (ux(i,j,k)*X(i,1,1) + uy(i,j,k)*Y(1,j,1) &
                + uz(i,j,k)*Z(1,1,k))/r
! store g_{0r} in uy
           uy(i,j,k) = (bx(i,j,k)*X(i,1,1) + by(i,j,k)*Y(1,j,1) &
                + bz(i,j,k)*Z(1,1,k))/r
        end do
     end do
  end do
! compute T_{r0}
  Tr0 = rho_b*h*ux*ut + P*uy
!
  return
end subroutine find_tr0
!-----------------------------------------------------------------------------
!
! phi -> phi + const * PsiRes
!
!-----------------------------------------------------------------------------
subroutine addham(ex, phi, PsiRes)
  implicit none
! Input parameters:
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: phi,PsiRes
!
  real*8 :: const
!
  const = 2.d-6
  phi = phi + const * PsiRes
!
  return
end subroutine addham

!-----------------------------------------------!
! Compute Omega and Omega^2		        !
!-----------------------------------------------!
subroutine compute_omega(ex, X,Y, v_x,v_y, omega, omega_square)
  implicit none
!
! Input parameters:
!
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y
  real*8, dimension(ex(1),ex(2),ex(3))     :: v_x,v_y
!
! Other variables:
!
  real*8, dimension(ex(1),ex(2),ex(3))     :: omega, omega_square
  integer				   :: imin,jmin,imax,jmax,i,j,k
  integer                                  :: kmin,kmax
  real*8				   :: R_cy2
!
  imin = lbound(v_x,1)
  imax = ubound(v_x,1)
  jmin = lbound(v_x,2)
  jmax = ubound(v_x,2)
  kmin = lbound(v_x,3)
  kmax = ubound(v_x,3)  
  do k=kmin,kmax
     do j=jmin,jmax
	do i=imin,imax
	   R_cy2 = X(i,1,1)*X(i,1,1) + Y(1,j,1)*Y(1,j,1)
	   omega(i,j,k) = (X(i,1,1)*v_y(i,j,k)-Y(1,j,1)*v_x(i,j,k))/R_cy2
	   omega_square(i,j,k) = omega(i,j,k)*omega(i,j,k)
	end do 
     end do 
  end do 
end subroutine compute_omega

subroutine compute_vol_intg(ex,rho_b,phi,e6phi)
  implicit none
  integer, dimension(3)				:: ex
  real*8, dimension(ex(1),ex(2),ex(3))		:: rho_b,e6phi,phi
  integer 					:: imin,imax,jmin,jmax,kmin,kmax
  integer					:: i,j,k
!
  imin = lbound(rho_b,1)
  imax = ubound(rho_b,1)
  jmin = lbound(rho_b,2)
  jmax = ubound(rho_b,2)
  kmin = lbound(rho_b,3)
  kmax = ubound(rho_b,3)
  e6phi = exp(6.d0*phi)
  do k=kmin,kmax
     do j=jmin,jmax
	do i=imin,imax
	   if (rho_b(i,j,k)==0.d0) e6phi(i,j,k) = 0.d0
	end do
     end do
  end do
end subroutine compute_vol_intg

subroutine compute_inverse_length2(ex,X,Y,Z,alpha,phi,gupxx,gupxy,gupxz, &
			gupyy,gupyz,gupzz,v_x,v_y,il2,Symmetry,rho_star,rho_tiny)
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
     subroutine gderivs_pi(ex,f,fx,fy,fz,dX,dY,dZ,SYM3,   &
                     glob_imin,glob_jmin,glob_kmin,Nx,Nz,  &
                     PI_SYMM_f)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8                                   :: dX,dY,dZ,SYM3
       integer                                  :: Nx,Nz
       integer                                  :: glob_imin,glob_jmin,glob_kmin
       real*8, dimension(Nx+1,Nz+1)             :: PI_SYMM_f
     end subroutine gderivs_pi
  end interface

  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: alpha,phi,il2,rho_star
  real*8, dimension(ex(1),ex(2),ex(3))        :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3))        :: omega,omegax,omegay,omegaz,v_y,v_x
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8				      :: dX,dY,dZ,rho_tiny
  integer				      :: i,j,k,imax,imin,jmax,jmin,kmin,kmax
  integer                                     :: Symmetry,AXISYM,EQUATORIAL
  real*8, parameter                           :: SYM = 1.d0, ANTI = -1.d0
  logical kill
 PARAMETER (EQUATORIAL = 1,AXISYM = 4)
!
  imin = lbound(phi,1)
  jmin = lbound(phi,2)
  kmin = lbound(phi,3)
  imax = ubound(phi,1)
  jmax = ubound(phi,2)
  kmax = ubound(phi,3)

  dX = X(imin+1,1,1)-X(imin,1,1)
  dY = Y(1,jmin+1,1)-Y(1,jmin,1)
  dZ = Z(1,1,kmin+1)-Z(1,1,kmin)

  do k=kmin,kmax
     do j=jmin,jmax
	do i=imin,imax
           omega(i,j,k) = (X(i,1,1)*v_y(i,j,k)-Y(1,j,1)*v_x(i,j,k))/(X(i,1,1)*X(i,1,1)+Y(1,j,1)*Y(1,j,1))
	end do
     end do
  end do

  if (Symmetry==AXISYM) then
     call gderivs_axi(ex,omega,omegax,omegay,omegaz,dX,dY,dZ,SYM,SYM,SYM)
  else
     call gderivs_eq(ex,omega,omegax,omegay,omegaz,dX,dY,dZ,SYM,SYM,SYM)
  end if

  omegay=0.d0

  do k=kmin+1,kmax-1
     do j=jmin+1,jmax-1
	do i=imin+1,imax-1
	   kill=.FALSE.
	   if (omega(i,j,k) .eq. 0.d0) kill=.TRUE.
	   if (omega(i+1,j,k) .eq. 0.d0) kill=.TRUE.
	   if (omega(i-1,j,k) .eq. 0.d0) kill=.TRUE.
	   if (omega(i,j,k+1) .eq. 0.d0) kill=.TRUE.
           if (omega(i,j,k-1) .eq. 0.d0) kill=.TRUE.
	   if (kill) then
	      omegax(i,j,k)=0.d0
	      omegaz(i,j,k)=0.d0
	   end if
	end do
     end do
  end do

  where (rho_star>rho_tiny) 
    il2 = (1.d0+alpha)*exp(6.d0*phi)* ( gupxx*gupxx*omegax*omegax + &
	   2.d0*gupxx*gupxy*omegax*omegay + 2.d0*gupxx*gupxz*omegax*omegaz + &
	   gupyy*gupyy*omegay*omegay + 2.d0*gupyy*gupyz*omegay*omegaz + &
	   gupzz*gupzz*omegaz*omegaz)/omega/omega
  elsewhere 
    il2 = 0.d0
  end where

end subroutine compute_inverse_length2
