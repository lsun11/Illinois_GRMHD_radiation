!------------------------------------------------------------------------------- 
! When this function is called, Ai's are defined everywhere from prolongation.
! For safety, we apply symmetry BC's on the Ai's before computing Bi's from the
! Ai's.
!-------------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine compute_B_from_A(CCTK_ARGUMENTS)
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

  ! Ai's NaN fixers:
  real*8  :: Ax_avg,Ay_avg,Az_avg
  integer :: numpts

!
! First apply symmetry BC on Ai
!
  !!call CartSymGN(dummy,cctkGH,'mhd_evolve::em_Ax')
  !!call CartSymGN(dummy,cctkGH,'mhd_evolve::em_Ay')
  !!call CartSymGN(dummy,cctkGH,'mhd_evolve::em_Az')

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  !$omp parallel do private(Ax_avg,Ay_avg,Az_avg)
  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
           ! ZACH SAYS: This should happen RARELY.
           if(isnan(Ax(i,j,k)).eqv. .true.) then
              write(*,*) "FIXING Ax!",i,j,k
              Ax_avg=0; numpts=0
              if(i.gt.1 .and. isnan(Ax(i-1,j,k)).eqv. .false.) then
                 Ax_avg=Ax_avg + Ax(i-1,j,k)
                 numpts=numpts + 1
              end if
              if(j.gt.1 .and. isnan(Ax(i,j-1,k)).eqv. .false.) then
                 Ax_avg=Ax_avg + Ax(i,j-1,k)
                 numpts=numpts + 1
              end if
              if(k.gt.1 .and. isnan(Ax(i,j,k-1)).eqv. .false.) then
                 Ax_avg=Ax_avg + Ax(i,j,k-1)
                 numpts=numpts + 1
              end if
              if(i.lt.cctk_lsh(1) .and. isnan(Ax(i+1,j,k)).eqv. .false.) then
                 Ax_avg=Ax_avg + Ax(i+1,j,k)
                 numpts=numpts + 1
              end if
              if(j.lt.cctk_lsh(2) .and. isnan(Ax(i,j+1,k)).eqv. .false.) then
                 Ax_avg=Ax_avg + Ax(i,j+1,k)
                 numpts=numpts + 1
              end if
              if(k.lt.cctk_lsh(3) .and. isnan(Ax(i,j,k+1)).eqv. .false.) then
                 Ax_avg=Ax_avg + Ax(i,j,k+1)
                 numpts=numpts + 1
              end if
              if(numpts==0) then
                 Ax(i,j,k) = 0.D0
              else
                 Ax(i,j,k) = Ax_avg/numpts
              end if
	      write(*,*) "Ax(i,j,k), numpts=", Ax(i,j,k),numpts
           end if

           ! ZACH SAYS: This should happen RARELY.
           if(isnan(Ay(i,j,k)).eqv. .true.) then
              write(*,*) "FIXING Ay!",i,j,k
              Ay_avg=0; numpts=0
              if(i.gt.1 .and. isnan(Ay(i-1,j,k)).eqv. .false.) then
                 Ay_avg=Ay_avg + Ay(i-1,j,k)
                 numpts=numpts + 1
              end if
              if(j.gt.1 .and. isnan(Ay(i,j-1,k)).eqv. .false.) then
                 Ay_avg=Ay_avg + Ay(i,j-1,k)
                 numpts=numpts + 1
              end if
              if(k.gt.1 .and. isnan(Ay(i,j,k-1)).eqv. .false.) then
                 Ay_avg=Ay_avg + Ay(i,j,k-1)
                 numpts=numpts + 1
              end if
              if(i.lt.cctk_lsh(1) .and. isnan(Ay(i+1,j,k)).eqv. .false.) then
                 Ay_avg=Ay_avg + Ay(i+1,j,k)
                 numpts=numpts + 1
              end if
              if(j.lt.cctk_lsh(2) .and. isnan(Ay(i,j+1,k)).eqv. .false.) then
                 Ay_avg=Ay_avg + Ay(i,j+1,k)
                 numpts=numpts + 1
              end if
              if(k.lt.cctk_lsh(3) .and. isnan(Ay(i,j,k+1)).eqv. .false.) then
                 Ay_avg=Ay_avg + Ay(i,j,k+1)
                 numpts=numpts + 1
              end if
              if(numpts==0) then
                 Ay(i,j,k) = 0.D0
              else
                 Ay(i,j,k) = Ay_avg/numpts
              end if
           end if

           ! ZACH SAYS: This should happen RARELY.
           if(isnan(Az(i,j,k)).eqv. .true.) then
              write(*,*) "FIXING Az!",i,j,k
              Az_avg=0; numpts=0
              if(i.gt.1 .and. isnan(Az(i-1,j,k)).eqv. .false.) then
                 Az_avg=Az_avg + Az(i-1,j,k)
                 numpts=numpts + 1
              end if
              if(j.gt.1 .and. isnan(Az(i,j-1,k)).eqv. .false.) then
                 Az_avg=Az_avg + Az(i,j-1,k)
                 numpts=numpts + 1
              end if
              if(k.gt.1 .and. isnan(Az(i,j,k-1)).eqv. .false.) then
                 Az_avg=Az_avg + Az(i,j,k-1)
                 numpts=numpts + 1
              end if
              if(i.lt.cctk_lsh(1) .and. isnan(Az(i+1,j,k)).eqv. .false.) then
                 Az_avg=Az_avg + Az(i+1,j,k)
                 numpts=numpts + 1
              end if
              if(j.lt.cctk_lsh(2) .and. isnan(Az(i,j+1,k)).eqv. .false.) then
                 Az_avg=Az_avg + Az(i,j+1,k)
                 numpts=numpts + 1
              end if
              if(k.lt.cctk_lsh(3) .and. isnan(Az(i,j,k+1)).eqv. .false.) then
                 Az_avg=Az_avg + Az(i,j,k+1)
                 numpts=numpts + 1
              end if
              if(numpts==0) then
                 Az(i,j,k) = 0.D0
              else
                 Az(i,j,k) = Az_avg/numpts
              end if
           end if
        end do
     end do
  end do
  !$omp end parallel do


  if (Symmetry==NO_SYMM) then 
     Ax_syms = 1
     Ay_syms = 1
     Az_syms = 1
     Bx_syms = 1
     By_syms = 1
     Bz_syms = 1
     psi6phi_syms(1) = 1
     psi6phi_syms(2) = 1
     psi6phi_syms(3) = 1
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

     ! Special BC for radial magnetized Bondi test: impose (r+z) Ax 
     ! and (r+z) Ay to be symmetry across the z-axis
     ! 
     if (A_BC_rad_mag_bondi_enable==1) then 
	Ax_syms = 1
	Ay_syms = 1
	Az_syms = 1

        ! Multiply Ax and Ay by (r+z)
	!$omp parallel do 
	do k=1,cctk_lsh(3)
	   do j=1,cctk_lsh(2)
	      do i=1,cctk_lsh(1)
	         Ax(i,j,k) = Ax(i,j,k) * ( sqrt(X(i,j,k)**2 + & 
						(Y(i,j,k)+0.5d0*dY)**2 + &
						(Z(i,j,k)+0.5d0*dZ)**2) + & 
					   Z(i,j,k)+0.5d0*dZ )
	         Ay(i,j,k) = Ay(i,j,k) * ( sqrt((X(i,j,k)+0.5d0*dX)**2 + &
                                                Y(i,j,k)**2 + &
                                                (Z(i,j,k)+0.5d0*dZ)**2) + &
                                           Z(i,j,k)+0.5d0*dZ ) 
	      end do
	   end do
	end do
	!$omp end parallel do
     end if
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

  if (Symmetry==EQUATORIAL .and. A_BC_rad_mag_bondi_enable==1) then 
     ! Divide Ax and Ay by (r+z)
     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              Ax(i,j,k) = Ax(i,j,k) / ( sqrt(X(i,j,k)**2 + &
                                                (Y(i,j,k)+0.5d0*dY)**2 + &
                                                (Z(i,j,k)+0.5d0*dZ)**2) + &
                                           Z(i,j,k)+0.5d0*dZ )  
              Ay(i,j,k) = Ay(i,j,k) / ( sqrt((X(i,j,k)+0.5d0*dX)**2 + &
                                                Y(i,j,k)**2 + &
                                                (Z(i,j,k)+0.5d0*dZ)**2) + &
                                           Z(i,j,k)+0.5d0*dZ ) 
           end do
        end do
     end do
     !$omp end parallel do
  end if

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
 


  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           if(isnan(Ax(i,j,k)).eqv. .true.) then
              write(*,*) "At the end of compute_B_from_A.F90, Ax is NAN",i,j,k
           end if    

           if(isnan(Ay(i,j,k)).eqv. .true.) then
              write(*,*) "At the end of compute_B_from_A.F90, Ay is NAN",i,j,k
           end if

           if(isnan(Az(i,j,k)).eqv. .true.) then
              write(*,*) "At the end of compute_B_from_A.F90, Az is NAN",i,j,k
           end if

        end do
     end do
  end do


end subroutine compute_B_from_A






subroutine mhd_do_symmetry_gz_one_gridfunction_stagger(ext,X,Y,Z,Symmetry,gridfunc,gridfunc_syms, & 
                              stagger_x,stagger_y,stagger_z)
  implicit none

  integer, dimension(3)                    :: ext,gridfunc_syms
  real*8,dimension(ext(1),ext(2),ext(3))   :: X,Y,Z,gridfunc
  real*8                                   :: dX,dY,dZ
  integer				   :: stagger_x,stagger_y,stagger_z
  real*8                                   :: x_offset,y_offset,z_offset
  integer                                  :: Symmetry
  integer :: i

  dX = X(2,1,1) - X(1,1,1)
  dY = Y(1,2,1) - Y(1,1,1)
  dZ = Z(1,1,2) - Z(1,1,1)

  if(Symmetry.eq.1) then
     z_offset = dZ*0.5d0*stagger_z
     if(abs(z(1,1,2)+z_offset).lt.dz*0.1) then
        gridfunc(:,:,1) = gridfunc_syms(3)*gridfunc(:,:,3)
     else if(abs(z(1,1,3)+z_offset).lt.dz*0.1) then
        do i=1,2
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,6-i)
        end do
     else if(abs(z(1,1,4)+z_offset).lt.dz*0.1) then
        do i=1,3
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,8-i)
        end do
     else if(abs(z(1,1,5)+z_offset).lt.dz*0.1) then
        do i=1,4
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,10-i)
        end do
     else if(abs(z(1,1,6)+z_offset).lt.dz*0.1) then
        do i=1,5
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,12-i)
        end do
        !=============================================
     else if(abs(z(1,1,2)+z_offset-dZ*0.5D0).lt.dz*0.1) then
        gridfunc(:,:,1) = gridfunc_syms(3)*gridfunc(:,:,2)
     else if(abs(z(1,1,3)+z_offset-dZ*0.5D0).lt.dz*0.1) then
        do i=1,2
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,5-i)
        end do
     else if(abs(z(1,1,4)+z_offset-dZ*0.5D0).lt.dz*0.1) then
        do i=1,3
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,7-i)
        end do
     else if(abs(z(1,1,5)+z_offset-dZ*0.5D0).lt.dz*0.1) then
        do i=1,4
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,9-i)
        end do
     else if(abs(z(1,1,6)+z_offset-dZ*0.5D0).lt.dz*0.1) then
        do i=1,5
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,11-i)
        end do
     else
        ! Uncomment the following only on 1 processor:
        !write(*,*) "This grid type is not supported!  Z-axis coordinates = ",ext,gridfunc_syms,z(1,1,:)
        !stop
     end if
  else if (Symmetry .ne. 0) then
     write(*,*) "Sorry, Symmetry!=1 not currently supported!"
     stop
  end if
end subroutine mhd_do_symmetry_gz_one_gridfunction_stagger

subroutine apply_bfield_bc_ijkmin(ext,Bx,By,Bz,EM_BC)
  implicit none
  integer, dimension(3) :: ext
  integer :: i,j,k, EM_BC
  integer, parameter :: COPY = 1, FREEZE = 2, EXTRAP = 3, QUAD = 4, PLANAR = 5
  real*8, dimension(ext(1),ext(2),ext(3)) :: Bx,By,Bz

  if (EM_BC==FREEZE) then 
     write(*,*) 'Frozen EM BC not supported for constrained_transport_scheme=3'
     stop
  end if

  if (EM_BC==COPY) then 
     !$omp parallel do
     do k=2,ext(3)
	do j=2,ext(2)
	   Bx(1,j,k) = Bx(2,j,k)
           By(1,j,k) = By(2,j,k)
           Bz(1,j,k) = Bz(2,j,k)
	end do
     end do
     !$omp end parallel do
     !$omp parallel do
     do k=2,ext(3)
        do i=1,ext(1)
           Bx(i,1,k) = Bx(i,2,k)
           By(i,1,k) = By(i,2,k)
           Bz(i,1,k) = Bz(i,2,k)
        end do
     end do
     !$omp end parallel do
     !$omp parallel do
     do j=1,ext(2)
        do i=1,ext(1)
           Bx(i,j,1) = Bx(i,j,2)
           By(i,j,1) = By(i,j,2)
           Bz(i,j,1) = Bz(i,j,2)
        end do
     end do
     !$omp end parallel do
  end if

  if (EM_BC==EXTRAP) then 
     !$omp parallel do
     do k=2,ext(3)
        do j=2,ext(2)
           Bx(1,j,k) = 2.d0*Bx(2,j,k) - Bx(3,j,k)
           By(1,j,k) = 2.d0*By(2,j,k) - By(3,j,k)
           Bz(1,j,k) = 2.d0*Bz(2,j,k) - Bz(3,j,k)
        end do
     end do
     !$omp end parallel do
     !$omp parallel do
     do k=2,ext(3)
        do i=1,ext(1)
           Bx(i,1,k) = 2.d0*Bx(i,2,k) - Bx(i,3,k)
           By(i,1,k) = 2.d0*By(i,2,k) - By(i,3,k)
           Bz(i,1,k) = 2.d0*Bz(i,2,k) - Bz(i,3,k)
        end do
     end do
     !$omp end parallel do
     !$omp parallel do
     do j=1,ext(2)
        do i=1,ext(1)
           Bx(i,j,1) = 2.d0*Bx(i,j,2) - Bx(i,j,3)
           By(i,j,1) = 2.d0*By(i,j,2) - By(i,j,3)
           Bz(i,j,1) = 2.d0*Bz(i,j,2) - Bz(i,j,3)
        end do
     end do
     !$omp end parallel do
  end if

  if (EM_BC==PLANAR) then
     write(*,*) "EM_BC=PLANAR not yet supported!"
     stop
  end if

  if (EM_BC==QUAD) then
     !$omp parallel do
     do k=2,ext(3)
        do j=2,ext(2)
           Bx(1,j,k) = 3.d0*Bx(2,j,k) - 3.d0*Bx(3,j,k) + Bx(4,j,k)
           By(1,j,k) = 3.d0*By(2,j,k) - 3.d0*By(3,j,k) + By(4,j,k)
           Bz(1,j,k) = 3.d0*Bz(2,j,k) - 3.d0*Bz(3,j,k) + Bz(4,j,k)
        end do
     end do
     !$omp end parallel do
     !$omp parallel do
     do k=2,ext(3)
        do i=1,ext(1)
           Bx(i,1,k) = 3.d0*Bx(i,2,k) - 3.d0*Bx(i,3,k) + Bx(i,4,k)
           By(i,1,k) = 3.d0*By(i,2,k) - 3.d0*By(i,3,k) + By(i,4,k)
           Bz(i,1,k) = 3.d0*Bz(i,2,k) - 3.d0*Bz(i,3,k) + Bz(i,4,k)
        end do
     end do
     !$omp end parallel do
     !$omp parallel do
     do j=1,ext(2)
        do i=1,ext(1)
           Bx(i,j,1) = 3.d0*Bx(i,j,2) - 3.d0*Bx(i,j,3) + Bx(i,j,4)
           By(i,j,1) = 3.d0*By(i,j,2) - 3.d0*By(i,j,3) + By(i,j,4)
           Bz(i,j,1) = 3.d0*Bz(i,j,2) - 3.d0*Bz(i,j,3) + Bz(i,j,4)
        end do
     end do
     !$omp end parallel do
  end if

end subroutine apply_bfield_bc_ijkmin
