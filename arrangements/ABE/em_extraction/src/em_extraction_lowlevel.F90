#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine em_extraction_lowlevel(cctkGH,cctk_lsh,cctk_nghostzones, &
     dX,dY,dZ,r,x,y,z, &
     gxx,gxy,gxz,gyy,gyz,gzz, &
     lapm1,shiftx,shifty,shiftz, &
     phi,psi, &
     Bx,By,Bz,Ex,Ey,Ez, &
     phi0re,phi0im,phi2re,phi2im)
  implicit none
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !Input variables:
  CCTK_POINTER                                           :: cctkGH
  integer,dimension(3)                                   :: cctk_lsh,cctk_nghostzones
  real*8                                                 :: dX,dY,dZ
  real*8, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: r,x,y,z
  real*8, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: gxx,gxy,gxz,gyy,gyz,gzz,Bx,By,Bz,Ex,Ey,Ez
  real*8, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: lapm1,shiftx,shifty,shiftz
  real*8, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: phi,psi

  !Output gridfunctions:
  real*8, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: phi0re,phi0im,phi2re,phi2im


  integer nx,ny,nz
  CCTK_REAL one,two,fac
  integer i,j,k,local_spatial_order
  integer istart,iend,jstart,jend,kstart,kend

  real*8 phixL,phiyL,phizL,func
  real*8 psi_4thpower
  real*8 gxxL,gxyL,gxzL,gyyL,gyzL,gzzL
  real*8 psixL,psiyL,psizL
  real*8 shiftxL,shiftyL,shiftzL,lapse
  real*8 g4tt, g4tx, g4ty, g4tz, g4xx, g4xy, g4xz, g4yy, g4yz, g4zz
  real*8 Ftx, Fty, Ftz, Fxy, Fxz, Fyz, F_tx, F_ty, F_tz, F_xy, F_xz, F_yz
  real*8 BxL, ByL, BzL, B_xL, B_yL, B_zL, ExL, EyL, EzL, xL, yL, zL
  real*8 r2, dthx, dthy,dthz, dphx, dphy, dphz, rhoL
  real*8 drdotdth, drx, dry, drz, thx,thy, thz, th2, phx, phy,phz, ph2,drdotdph, thdotdph
  real*8 hatt_t, hatt_x, hatt_y, hatt_z
  real*8 hatr_x, hatr_y, hatr_z
  real*8 hatph_x, hatph_y, hatph_z
  real*8 hatth_x, hatth_y, hatth_z
  real*8 norm_hat_r, norm_hat_th, norm_hat_ph
  real*8 hatr_dot_hatth, hatr_dot_hatph, hatth_dot_hatph
  real*8 detgam, shift_xL, shift_yL, shift_zL
  real*8 F_hatt_hatth,F_hatr_hatth,F_hatt_hatph,F_hatr_hatph, Rephi0,Imphi0 ,Rephi2,Imphi2

  write(*,*) "em_extraction begins..."

  !     Set up variables used in the grid loop for the physical grid points
  istart = cctk_nghostzones(1)+1;
  jstart = cctk_nghostzones(2)+1;
  kstart = cctk_nghostzones(3)+1;
  iend = cctk_lsh(1) - cctk_nghostzones(1);
  jend = cctk_lsh(2) - cctk_nghostzones(2);
  kend = cctk_lsh(3) - cctk_nghostzones(3);

  !     define nx,ny,nz from cctk_lsh(DIM)
  nx=cctk_lsh(1)
  ny=cctk_lsh(2)
  nz=cctk_lsh(3)

  two = 2.0D0
  one = 1.0D0


!!!$omp parallel do
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
           if(r(i,j,k).gt.compute_phi2_min_radius .and. r(i,j,k).lt.compute_phi2_max_radius) then

              psi_4thpower = exp(phi(i,j,k))
              psi_4thpower =psi_4thpower*psi_4thpower*psi_4thpower*psi_4thpower
! Vasilis says: We shouldn't need the following 6 lines
!              local_spatial_order=4
!              if(i==1.or.j==1.or.k==1.or.i==cctk_lsh(1).or.j==cctk_lsh(2).or.k==cctk_lsh(3)) then
!                 local_spatial_order=0
!              else if(i==2.or.j==2.or.k==2.or.i==cctk_lsh(1)-1.or.j==cctk_lsh(2)-1.or.k==cctk_lsh(3)-1) then
!                 local_spatial_order=2
!              end if

              !     Define local variables to calculate the orthonormal tetrad

	      ! The 3 metric
              gxxL=gxx(i,j,k)*psi_4thpower
              gxyL=gxy(i,j,k)*psi_4thpower
              gxzL=gxz(i,j,k)*psi_4thpower
              gyyL=gyy(i,j,k)*psi_4thpower
              gyzL=gyz(i,j,k)*psi_4thpower
              gzzL=gzz(i,j,k)*psi_4thpower
	      ! The shift vector
	      shiftxL = shiftx(i,j,k)
	      shiftyL = shifty(i,j,k)
	      shiftzL = shiftz(i,j,k)
	      ! The lapse function
	      lapse = lapm1(i,j,k) + 1.d0

	      ! The coordinates
              xL=x(i,j,k)
              yL=y(i,j,k)
              zL=z(i,j,k)

!===============================================================================
              ! FIRST: Set up the tetrad vectors assuming spherical coordinates
	      ! So that in flat spacetime, \hat r = e_{\hat r},
              ! \hat \phi = e_{\hat \phi}, \hat \theta = e_{\hat \theta}  
	      ! centered at the origin of the coordinate system
!===============================================================================

	! Components of the radial vector
	drx = xL;
	dry = yL;
	drz = zL;

	! Magnitude of (radial) position vector squared

	r2 = gxxL*xL*xL + gyyL*yL*yL + gzzL*zL*zL + 2.d0*(gxyL*xL*yL + gxzL*xL*zL + gyzL*yL*zL)  

	! Flat-space distance from Z-axis
	rhoL = sqrt(xL*xL+yL*yL)
	! cartesian components of partial_theta
	dthx = xL*zL/rhoL
	dthy = yL*zL/rhoL
	dthz = -rhoL

	! Start Gram-Schmidt orhonormalization process
	! The timelike unit vector normal to t=const hypersurfaces
	! will serve as the orthonormal timelike \hat t^a vector
        ! The radial vector \partial_r = x^i\partial_i is chosen next
        ! and since it is purely spatial it is already normal to \hat t^a
	! so we only need to normalize it. 
	! Next we orthonormalize \partial_\theta
	
	! \partial_r \dot \partial_\theta

	drdotdth = drx*dthx*gxxL + dry*dthx*gxyL + drx*dthy*gxyL + drz*dthx*gxzL &
	         + drx*dthz*gxzL + dry*dthy*gyyL + drz*dthy*gyzL + dry*dthz*gyzL + drz*dthz*gzzL

	! Cartesian Components of theta vector normal to \hat r, which is also normal to r
	thx = dthx - (drdotdth/r2)*drx
	thy = dthy - (drdotdth/r2)*dry
	thz = dthz - (drdotdth/r2)*drz

	! (Magnitude of theta)^2 
	th2 = gxxL*thx*thx + gyyL*thy*thy + gzzL*thz*thz &
	      + 2.d0*(gxyL*thx*thy + gxzL*thx*thz + gyzL*thy*thz)  

	! Cartesian Components of vector \partial_\phi
	dphx = -yL
	dphy = xL
	dphz = 0.d0


	! \partial_r \dot \partial_\phi
	drdotdph = dphx*drx*gxxL + dphy*drx*gxyL + dphx*dry*gxyL + dphz*drx*gxzL &
	         + dphx*drz*gxzL + dphy*dry*gyyL + dphz*dry*gyzL + dphy*drz*gyzL + dphz*drz*gzzL

	! \partial_\phi \dot \partial_\theta
	thdotdph = dphx*dthx*gxxL + dphy*dthx*gxyL + dphx*dthy*gxyL + dphz*dthx*gxzL & 
	         + dphx*dthz*gxzL + dphy*dthy*gyyL + dphz*dthy*gyzL + dphy*dthz*gyzL + dphz*dthz*gzzL

	! Components of phi vector which is normal to both \hat\theta and \hat r
	phx = dphx - (drdotdph/r2)*drx - (thdotdph/th2)*thx
	phy = dphy - (drdotdph/r2)*dry - (thdotdph/th2)*thy
	phx = dphx - (drdotdph/r2)*drz - (thdotdph/th2)*thz

	! (magnitude of \phi)^2
	ph2 = gxxL*phx*phx + gyyL*phy*phy + gzzL*phz*phz & 
	     + 2.d0*(gxyL*phx*phy + gxzL*phx*phz + gyzL*phy*phz)  


	! Orhonormal tetrad

	! \hat t^a = n^a
	hatt_t = 1.d0/lapse
	hatt_x = -shiftxL/lapse
	hatt_y = -shiftyL/lapse
	hatt_z = -shiftzL/lapse

	! \hat r
	hatr_x = drx/sqrt(r2)
	hatr_y = dry/sqrt(r2)
	hatr_z = drz/sqrt(r2)

	! \hat \theta
	hatth_x = thx/sqrt(th2) 
	hatth_y = thy/sqrt(th2) 
	hatth_z = thz/sqrt(th2)

	! \hat \phi
	hatph_x = phx/sqrt(ph2) 
	hatph_y = phy/sqrt(ph2) 
	hatph_z = phz/sqrt(ph2)


	! Check orthonormality of the tetrad

	norm_hat_r = gxxL*hatr_x*hatr_x + 2.d0*gxyL*hatr_x*hatr_y + gyyL*hatr_y*hatr_y & 
	             + 2.d0*gxzL*hatr_x*hatr_z + 2.d0*gyzL*hatr_y*hatr_z + gzzL*hatr_z*hatr_z

	norm_hat_th = gxxL*hatth_x*hatth_x + 2.d0*gxyL*hatth_x*hatth_y + gyyL*hatth_y*hatth_y &
	              + 2.d0*gxzL*hatth_x*hatth_z + 2.d0*gyzL*hatth_y*hatth_z + gzzL*hatth_z*hatth_z

	norm_hat_ph = gxxL*hatph_x*hatph_x + 2.d0*gxyL*hatph_x*hatph_y + gyyL*hatph_y*hatph_y &
	              + 2.d0*gxzL*hatph_x*hatph_z + 2.d0*gyzL*hatph_y*hatph_z + gzzL*hatph_z*hatph_z

	if ( abs(norm_hat_r-1.d0) > 1.e-13 .or. abs(norm_hat_th-1.d0) > 1.e-13 .or. abs(norm_hat_ph-1.d0) > 1.e-13) then
	  write(*,*) "Stopping: tetrad is not normalized properly in em_extraction_lowlevel"
	  write(*,*) "|hat r|=", norm_hat_r, "|hat th|=", norm_hat_th, "|hat ph|=", norm_hat_ph 
   	  stop
	endif

	hatr_dot_hatth = gxxL*hatr_x*hatth_x + gxyL*hatr_y*hatth_x + gxzL*hatr_z*hatth_x & 
	               + gxyL*hatr_x*hatth_y + gyyL*hatr_y*hatth_y + gyzL*hatr_z*hatth_y & 
	     	       + gxzL*hatr_x*hatth_z + gyzL*hatr_y*hatth_z + gzzL*hatr_z*hatth_z

	hatr_dot_hatph = gxxL*hatr_x*hatph_x + gxyL*hatr_y*hatph_x + gxzL*hatr_z*hatph_x & 
	               + gxyL*hatr_x*hatph_y + gyyL*hatr_y*hatph_y + gyzL*hatr_z*hatph_y & 
		       + gxzL*hatr_x*hatph_z + gyzL*hatr_y*hatph_z + gzzL*hatr_z*hatph_z

	hatth_dot_hatph = gxxL*hatth_x*hatph_x + gxyL*hatth_y*hatph_x + gxzL*hatth_z*hatph_x & 
			+ gxyL*hatth_x*hatph_y + gyyL*hatth_y*hatph_y + gyzL*hatth_z*hatph_y &
		        + gxzL*hatth_x*hatph_z + gyzL*hatth_y*hatph_z + gzzL*hatth_z*hatph_z

	if ( abs(hatr_dot_hatth) > 1.e-13 .or. abs(hatr_dot_hatph) > 1.e-13 .or. abs(hatth_dot_hatph) > 1.e-13) then
	  write(*,*) "Stopping: tetrad is not orthogonal in em_extraction_lowlevel"
	  write(*,*) "at grid point x, y, z = ", xL, yL, zL
	  write(*,*) "|hat r dot hat th|=", hatr_dot_hatth, "|hat r dot hat ph|=", hatr_dot_hatph, "|hat th dot hat ph|=", hatth_dot_hatph
          write(*,*) "|hat thx|=", hatth_x,  "|hat thy|=", hatth_y,  "|hat thz|=", hatth_z 
          write(*,*) "|dthx|=", dthx,  "|dthy|=", dthy,  "|dthz|=", dthz 
          write(*,*) "|thx|=", thx,  "|thy|=", thy,  "|thz|=", thz 
          write(*,*) "|gxx|=", gxxL, "|gxy|=", gxyL, "|gxz|=", gxzL
          write(*,*) "|gyy|=", gyyL, "|gyz|=", gyzL, "|gzz|=", gzzL
          write(*,*) "psi4=", psi_4thpower, "psi4_b=", exp(4.d0*phi(i,j,k))
          write(*,*) "|gxx|=", gxx(i,j,k), "|gxy|=", gxy(i,j,k), "|gxz|=", gxz(i,j,k)
          write(*,*) "|gyy|=", gyy(i,j,k), "|gyz|=", gyz(i,j,k), "|gzz|=", gzz(i,j,k)

	  stop
	endif
	! =======================================================================================
         !          End of orthonormal tetrad calculation
	! =======================================================================================

	! We have the components of the orthonormal tetrad in the cartesian coordinate basis
	! The following are (contravariant) vector components
	! hatt = (hatt_t, hatt_x hatt_y, hatt_z) = n^a ,i.e., the time like normal vector	
	! hatr = (0, hatr_x, hatr_y, hatr_z)
	! hatth = (0, hatth_x, hatth_y, hatth_z)
	! hatph = (0, hatph_x, hatph_y, hatph_z)

	if( xL*yL.ne.0.d0) then
        !     Calculate only if not on the z-axis for the "radial" tetrad case


		! Now, we turn to the calculation of the Faraday tensor F_\mu\nu

		BxL = Bx(i,j,k)   ! These are the B^i
		ByL = By(i,j,k)
		BzL = Bz(i,j,k)

		ExL = Ex(i,j,k)   ! These are the E^i
		EyL = Ey(i,j,k)
		EzL = Ez(i,j,k)

		! First we need the covariant components of the B-field
		! B_i = \gamma_ij B^j

		B_xL = gxxL*BxL + gxyL*ByL + gxzL*BzL
		B_yL = gxyL*BxL + gyyL*ByL + gyzL*BzL
		B_zL = gxzL*BxL + gyzL*ByL + gzzL*BzL

		! We will also need the determinant of the 3-metric

		detgam = -gxzL*gxzL*gyyL + 2.d0*gxyL*gxzL*gyzL & 
		         -gxxL*gyzL*gyzL - gxyL*gxyL*gzzL + gxxL*gyyL*gzzL

		! Now compute the Maxwell tensor, F^{munu} (contravariant components) 

		Ftx = hatt_t*ExL
		Fty = hatt_t*EyL
		Ftz = hatt_t*EzL
		Fxy = hatt_x*EyL - hatt_y*ExL + B_zL/sqrt(detgam)
		Fxz = hatt_x*EzL - hatt_z*ExL - B_yL/sqrt(detgam)
		Fyz = hatt_y*EzL - hatt_z*EyL + B_xL/sqrt(detgam)


		! Compute 4-metric
		shift_xL = gxxL*shiftxL + gxyL*shiftyL + gxzL*shiftzL
		shift_yL = gxyL*shiftxL + gyyL*shiftyL + gyzL*shiftzL
		shift_zL = gxzL*shiftxL + gyzL*shiftyL + gzzL*shiftzL


		g4tt = -lapse*lapse + shift_xL*shiftxL + shift_yL*shiftyL + shift_zL*shiftzL
		g4tx = shift_xL
		g4ty = shift_yL
		g4tz = shift_zL
		g4xx = gxxL
		g4xy = gxyL
		g4xz = gxzL
		g4yy = gyyL
		g4yz = gyzL
		g4zz = gzzL


		! Compute F_munu (covariant components)

		F_tx = -Ftx*g4tx*g4tx - Fty*g4tx*g4ty - Ftz*g4tx*g4tz + Ftx*g4tt*g4xx & 
		       -Fxy*g4ty*g4xx - Fxz*g4tz*g4xx + Fty*g4tt*g4xy + Fxy*g4tx*g4xy & 	
		       -Fyz*g4tz*g4xy + Ftz*g4tt*g4xz + Fxz*g4tx*g4xz + Fyz*g4ty*g4xz
		F_ty = -Ftx*g4tx*g4ty - Fty*g4ty*g4ty - Ftz*g4ty*g4tz + Ftx*g4tt*g4xy & 
		       -Fxy*g4ty*g4xy - Fxz*g4tz*g4xy + Fty*g4tt*g4yy + Fxy*g4tx*g4yy & 
		       -Fyz*g4tz*g4yy + Ftz*g4tt*g4yz + Fxz*g4tx*g4yz + Fyz*g4ty*g4yz
		F_tz = -Ftx*g4tx*g4tz - Fty*g4ty*g4tz - Ftz*g4tz*g4tz + Ftx*g4tt*g4xz & 
		       -Fxy*g4ty*g4xz - Fxz*g4tz*g4xz + Fty*g4tt*g4yz + Fxy*g4tx*g4yz & 
		       -Fyz*g4tz*g4yz + Ftz*g4tt*g4zz + Fxz*g4tx*g4zz + Fyz*g4ty*g4zz
		F_xy = -Ftx*g4ty*g4xx + Ftx*g4tx*g4xy - Fty*g4ty*g4xy - Fxy*g4xy*g4xy & 
		       -Ftz*g4ty*g4xz - Fxz*g4xy*g4xz + Fty*g4tx*g4yy + Fxy*g4xx*g4yy & 
		       -Fyz*g4xz*g4yy + Ftz*g4tx*g4yz + Fxz*g4xx*g4yz + Fyz*g4xy*g4yz
		F_xz = -Ftx*g4tz*g4xx - Fty*g4tz*g4xy + Ftx*g4tx*g4xz - Ftz*g4tz*g4xz & 
		       -Fxy*g4xy*g4xz - Fxz*g4xz*g4xz + Fty*g4tx*g4yz + Fxy*g4xx*g4yz & 
		       -Fyz*g4xz*g4yz + Ftz*g4tx*g4zz + Fxz*g4xx*g4zz + Fyz*g4xy*g4zz
		F_yz = -Ftx*g4tz*g4xy + Ftx*g4ty*g4xz - Fty*g4tz*g4yy - Fxy*g4xz*g4yy & 
		       +Fty*g4ty*g4yz - Ftz*g4tz*g4yz + Fxy*g4xy*g4yz - Fxz*g4xz*g4yz & 
		       -Fyz*g4yz*g4yz + Ftz*g4ty*g4zz + Fxz*g4xy*g4zz + Fyz*g4yy*g4zz


		! Calculate the Faraday tensor F dotted into the orhonormal tetrad 
		! In slot notation of two-forms F(\hat t,\hat \theta)
		F_hatt_hatth = F_tx*hatth_x*hatt_t + F_ty*hatth_y*hatt_t + F_tz*hatth_z*hatt_t & 
		             + F_xy*hatth_y*hatt_x + F_xz*hatth_z*hatt_x - F_xy*hatth_x*hatt_y & 
		             + F_yz*hatth_z*hatt_y - F_xz*hatth_x*hatt_z - F_yz*hatth_y*hatt_z
	
		! In slot notation of two-forms F(\hat r,\hat \theta)
		F_hatr_hatth = -F_xy*hatr_y*hatth_x - F_xz*hatr_z*hatth_x + F_xy*hatr_x*hatth_y & 
		              - F_yz*hatr_z*hatth_y + F_xz*hatr_x*hatth_z + F_yz*hatr_y*hatth_z

		! In slot notation of two-forms F(\hat t,\hat \phi)
		F_hatt_hatph = F_tx*hatph_x*hatt_t + F_ty*hatph_y*hatt_t + F_tz*hatph_z*hatt_t & 
		             + F_xy*hatph_y*hatt_x + F_xz*hatph_z*hatt_x - F_xy*hatph_x*hatt_y & 
			     + F_yz*hatph_z*hatt_y - F_xz*hatph_x*hatt_z - F_yz*hatph_y*hatt_z
	
		! In slot notation of two-forms F(\hat r,\hat \phi)
		F_hatr_hatph = -F_xy*hatr_y*hatph_x - F_xz*hatr_z*hatph_x + F_xy*hatr_x*hatph_y & 
		              - F_yz*hatr_z*hatph_y + F_xz*hatr_x*hatph_z + F_yz*hatr_y*hatph_z


		Rephi0 = 0.5d0*(F_hatt_hatth + F_hatr_hatth)
		Imphi0 = 0.5d0*(F_hatt_hatph + F_hatr_hatph)

		Rephi2 = 0.5d0*(-F_hatt_hatth + F_hatr_hatth)
		Imphi2 = 0.d05*(F_hatt_hatph - F_hatr_hatph)



	        !     If required scale by radial coordinate
	        if (scale_with_radius_em == 1) then
	            fac = r(i,j,k)**radius_power_em
	         else
	            fac = 1.0
	         end if

	         !     Phis we do pairwise
	          phi0re(i,j,k) = fac*Rephi0
	          phi0im(i,j,k) = fac*Imphi0
	          phi2re(i,j,k) = fac*Rephi2
	          phi2im(i,j,k) = fac*Imphi2
            else	 
                 !     Set values to zero without calculating on the z-axis
	          phi0re(i,j,k) = 0.d0
        	  phi0im(i,j,k) = 0.d0
         	  phi2re(i,j,k) = 0.d0
  	          phi2im(i,j,k) = 0.d0
            endif
       else ! (r(i,j,k).lt.compute_phi2_min_radius .and. r(i,j,k).gt.compute_phi2_max_radius)
	   phi0re(i,j,k) = -100000.d0
           phi0im(i,j,k) = -100000.d0
           phi2re(i,j,k) = -100000.d0
  	   phi2im(i,j,k) = -100000.d0
       end if


        enddo
     enddo
  enddo
!!!$omp end parallel do

  write(*,*) "outside em_extraction_lowlevel()"
  return
end subroutine em_extraction_lowlevel
