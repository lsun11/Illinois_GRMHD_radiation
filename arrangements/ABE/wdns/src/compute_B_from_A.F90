!------------------------------------------------------------------------------- 
! When this function is called, Ai's are defined everywhere from prolongation.
! For safety, we apply symmetry BC's on the Ai's before computing Bi's from the
! Ai's.
!-------------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine wdns_compute_B_from_A(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer               :: stagger_x,stagger_y,stagger_z,dummy,i,j,k
  integer               :: ip1,jp1,kp1
  real*8  		:: dX,dY,dZ, r_plus_z1, r_plus_z2
  real*8                :: v2,alp_u02,u02,u0L,u_x,u_y,u_z,bt_x,bt_y,bt_z,B2L,udotB,b2,vxi,vyi,vzi,gxxi,gxyi,gxzi,gyyi,gyzi,gzzi,p4,horizdirn_x,horizdirn_y,horizdirn_z,AH_radius_minimum
  integer, dimension(3) :: Ax_syms,Ay_syms,Az_syms,Bx_syms,By_syms,Bz_syms,psi6phi_syms
  integer, parameter    :: NO_SYMM = 0, EQUATORIAL = 1
  real*8, parameter     :: F1o4pi = 0.07957747154594766788
  integer, parameter    :: COPY = 1, FREEZE = 2, EXTRAP = 3, QUAD = 4, PLANAR = 5

! Storage variables:
  integer, parameter                 :: mmm = 17
  real*8, dimension(mmm)             :: AUX
  real*8                             :: rho_b_l, P_l, vx_l, vy_l, vz_l, Bx_l, By_l, Bz_l
  real*8                             :: rho_s_l, tau_l, mhd_st_x_l, mhd_st_y_l, mhd_st_z_l
  real*8                             :: st_x_l, st_y_l, st_z_l, w_l, h_l, u0_l, Ex_l, Ey_l, Ez_l
  real*8                             :: sbt_l, sbx_l, sby_l, sbz_l, rho_l, Sx_l, Sy_l, Sz_l
  real*8                             :: Sxx_l, Sxy_l, Sxz_l, Syy_l, Syz_l, Szz_l


!
! First apply symmetry BC on Ai
!
  !!call CartSymGN(dummy,cctkGH,'mhd_evolve::em_Ax')
  !!call CartSymGN(dummy,cctkGH,'mhd_evolve::em_Ay')
  !!call CartSymGN(dummy,cctkGH,'mhd_evolve::em_Az')

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  if (Symmetry==NO_SYMM) then 
     Ax_syms = 1
     Ay_syms = 1
     Az_syms = 1
     Bx_syms = 1
     By_syms = 1
     Bz_syms = 1
  elseif (Symmetry==EQUATORIAL) then
     Ax_syms = 1
     Ay_syms = 1
     Az_syms = 1
     Bx_syms = 1
     By_syms = 1
     Bz_syms = 1
     if (Sym_Bz .gt. 0.d0) then 
	Az_syms(3) = -1
	Bx_syms(3) = -1
        By_syms(3) = -1
     else
	Ax_syms(3) = -1
	Ay_syms(3) = -1
        Bz_syms(3) = -1
     end if

     psi6phi_syms(1) = 1 
     psi6phi_syms(2) = 1
     psi6phi_syms(3) = 1

  else
     write(*,*) 'Symmetry no supported yet for A evolution'
     stop
  end if

  !
  ! First apply symmetry BC on Ai and psi^6 Phi
  !
  ! Ax is defined on the semi-staggered grid (i,j+1/2,k+1/2)
  stagger_x = 0
  stagger_y = 1
  stagger_z = 1
  call mhd_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,Ax,Ax_syms,stagger_x, &
                        stagger_y,stagger_z)

  ! Ay is defined on the semi-staggered grid (i+1/2,j,k+1/2)
  stagger_x = 1
  stagger_y = 0
  stagger_z = 1
  call mhd_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,Ay,Ay_syms,stagger_x, &
                        stagger_y,stagger_z)

  ! Az is defined on the semi-staggered grid (i+1/2,j+1/2,k)
  stagger_x = 1
  stagger_y = 1
  stagger_z = 0
  call mhd_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,Az,Az_syms,stagger_x, &
                        stagger_y,stagger_z)

  ! Psi^6 Phi is defined on staggered grid (i+1/2,j+1/2,k+1/2)
  stagger_x = 1
  stagger_y = 1
  stagger_z = 1
  call mhd_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,psi6phi,psi6phi_syms,stagger_x, &
                        stagger_y,stagger_z)

  ! Compute the \tilde{B}^i on staggered grid (from 
  ! curl Ai) and temporarily store them in Bir
  !$omp parallel do
  do k=2,cctk_lsh(3)
     do j=2,cctk_lsh(2)
        do i=2,cctk_lsh(1)
           Bxr(i,j,k) = (Az(i,j,k)-Az(i,j-1,k))/dY &
                            - (Ay(i,j,k)-Ay(i,j,k-1))/dZ 
           Byr(i,j,k) = (Ax(i,j,k)-Ax(i,j,k-1))/dZ &
                            - (Az(i,j,k)-Az(i-1,j,k))/dX 
           Bzr(i,j,k) = (Ay(i,j,k)-Ay(i-1,j,k))/dX &
                            - (Ax(i,j,k)-Ax(i,j-1,k))/dY 

        end do
     end do
  end do
  !$omp end parallel do

  ! Apply BC on staggered \tilde{B}^i
  call apply_bfield_bc_ijkmin(cctk_lsh,Bxr,Byr,Bzr,EM_BC)

  !  Bx_stagger is defined at (i+1/2,j,k)
  stagger_x = 1
  stagger_y = 0
  stagger_z = 0
  call mhd_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,Bxr,Bx_syms,stagger_x, &
                        stagger_y,stagger_z)

  !  By_stagger is defined at (i,j+1/2,k)
  stagger_x = 0
  stagger_y = 1
  stagger_z = 0
  call mhd_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,Byr,By_syms,stagger_x, &
                        stagger_y,stagger_z)

  !  Bz_stagger is defined at (i,j,k+1/2)
  stagger_x = 0
  stagger_y = 0
  stagger_z = 1
  call mhd_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,Bzr,Bz_syms,stagger_x, &
                        stagger_y,stagger_z)

  ! Compute B^i on staggered grid 
  ! Note that at the outer boundary, phi->0, it may be ok to set 
  ! phi(imax,:,:)=phi(imax-1,:,:) etc.
  !$omp parallel do private (ip1,jp1,kp1)
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
	   ip1 = min(i+1,cctk_lsh(1))
           jp1 = min(j+1,cctk_lsh(2))
           kp1 = min(k+1,cctk_lsh(3))
           Bx_stagger(i,j,k) = Bxr(i,j,k)*exp(-3.d0*(phi(i,j,k)+phi(ip1,j,k)))
           By_stagger(i,j,k) = Byr(i,j,k)*exp(-3.d0*(phi(i,j,k)+phi(i,jp1,k)))
           Bz_stagger(i,j,k) = Bzr(i,j,k)*exp(-3.d0*(phi(i,j,k)+phi(i,j,kp1)))
        end do
     end do
  end do
  !$omp end parallel do

  ! Apply symmetry BC's on the staggered B^i's
  stagger_x = 1
  stagger_y = 0
  stagger_z = 0
  call mhd_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,Bx_stagger,Bx_syms,stagger_x, &
                        stagger_y,stagger_z)

  stagger_x = 0
  stagger_y = 1
  stagger_z = 0
  call mhd_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,By_stagger,By_syms,stagger_x, &
                        stagger_y,stagger_z)

  stagger_x = 0
  stagger_y = 0
  stagger_z = 1
  call mhd_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,Bz_stagger,Bz_syms,stagger_x, &
                        stagger_y,stagger_z)

  ! Now compute B^i on unstaggered grid by simple average.
  ! Note that B^i is continuous along x^i [see Del Zanna, 
  ! Bucciantini & Londrillo A&A 400, 397 (2003)], so averging 
  ! should be pretty good.
  !$omp parallel do
  do k=2,cctk_lsh(3)
     do j=2,cctk_lsh(2)
        do i=2,cctk_lsh(1)
           Bx(i,j,k) = 0.5d0*(Bx_stagger(i,j,k)+Bx_stagger(i-1,j,k))
           By(i,j,k) = 0.5d0*(By_stagger(i,j,k)+By_stagger(i,j-1,k))
           Bz(i,j,k) = 0.5d0*(Bz_stagger(i,j,k)+Bz_stagger(i,j,k-1))
        end do
     end do
  end do
  !$omp end parallel do

  ! Above loop filled in Bi's everywhere, except the imin, jmin, and kmin planes.
  !   We apply boundary conditions to fill those in.  In general, this will result
  !   in incorrect values set inside symmetry ghostzones.  However, we fill in the
  !   symmetry ghostzones next, overwriting incorrect values.
  call apply_bfield_bc_ijkmin(cctk_lsh,Bx,By,Bz,EM_BC)
  stagger_x = 0
  stagger_y = 0
  stagger_z = 0
  call mhd_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,Bx,Bx_syms,stagger_x, &
                        stagger_y,stagger_z)
  call mhd_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,By,By_syms,stagger_x, &
                        stagger_y,stagger_z)
  call mhd_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,Bz,Bz_syms,stagger_x, &
                        stagger_y,stagger_z)


! Finally we compute b^2, for diagnostic purposes
  call mhd_evolve_compute_b2_cpp(cctkGH,cctk_lsh, phi, lapm1, &
       shiftx,shifty,shiftz,vx,vy,vz,Bx,By,Bz, & 
       gxx, gxy, gxz, gyy, gyz, gzz, temp1)

  !$omp parallel do
  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
           if(i.gt.3 .and. j.gt.3 .and. k.gt.3 .and. &
                i.lt.cctk_lsh(1)-2 .and. j.lt.cctk_lsh(2)-2 .and. k.lt.cctk_lsh(3)-2) then
              smallb2(i,j,k) = temp1(i,j,k)
           else 
              smallb2(i,j,k) = -2.D0
           end if
        end do
     end do
  end do
  !$omp end parallel do
 
end subroutine wdns_compute_B_from_A
