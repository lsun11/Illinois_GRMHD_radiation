!------------------------------------------------------------------------------- 
! When this function is called, Ai's are defined everywhere from prolongation.
! For safety, we apply symmetry BC's on the Ai's before bhns_computing Bi's from the
! Ai's.
!-------------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bhns_compute_B_from_A(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer               :: stagger_x,stagger_y,stagger_z,dummy,i,j,k
  integer               :: ip1,jp1,kp1
  real*8  		:: dX,dY,dZ, r_plus_z1, r_plus_z2
  integer, dimension(3) :: Ax_syms,Ay_syms,Az_syms,Bx_syms,By_syms,Bz_syms
  integer, parameter    :: NO_SYMM = 0, EQUATORIAL = 1
  integer, parameter    :: COPY = 1, FREEZE = 2, EXTRAP = 3, QUAD = 4, PLANAR = 5
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
  else
     write(*,*) 'Symmetry no supported yet for A evolution'
     stop
  end if

  !
  ! First apply symmetry BC on Ai
  !
  ! Ax is defined on the semi-staggered grid (i,j+1/2,k+1/2)
  stagger_x = 0
  stagger_y = 1
  stagger_z = 1
  call bhns_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,Ax,Ax_syms,stagger_x, &
                        stagger_y,stagger_z)

  ! Ay is defined on the semi-staggered grid (i+1/2,j,k+1/2)
  stagger_x = 1
  stagger_y = 0
  stagger_z = 1
  call bhns_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,Ay,Ay_syms,stagger_x, &
                        stagger_y,stagger_z)

  ! Az is defined on the semi-staggered grid (i+1/2,j+1/2,k)
  stagger_x = 1
  stagger_y = 1
  stagger_z = 0
  call bhns_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,Az,Az_syms,stagger_x, &
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
  call bhns_apply_bfield_bc_ijkmin(cctk_lsh,Bxr,Byr,Bzr,EM_BC)

  !  Bx_stagger is defined at (i+1/2,j,k)
  stagger_x = 1
  stagger_y = 0
  stagger_z = 0
  call bhns_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,Bxr,Bx_syms,stagger_x, &
                        stagger_y,stagger_z)

  !  By_stagger is defined at (i,j+1/2,k)
  stagger_x = 0
  stagger_y = 1
  stagger_z = 0
  call bhns_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,Byr,By_syms,stagger_x, &
                        stagger_y,stagger_z)

  !  Bz_stagger is defined at (i,j,k+1/2)
  stagger_x = 0
  stagger_y = 0
  stagger_z = 1
  call bhns_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
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
  call bhns_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,Bx_stagger,Bx_syms,stagger_x, &
                        stagger_y,stagger_z)

  stagger_x = 0
  stagger_y = 1
  stagger_z = 0
  call bhns_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,By_stagger,By_syms,stagger_x, &
                        stagger_y,stagger_z)

  stagger_x = 0
  stagger_y = 0
  stagger_z = 1
  call bhns_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
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
  call bhns_apply_bfield_bc_ijkmin(cctk_lsh,Bx,By,Bz,EM_BC)
  stagger_x = 0
  stagger_y = 0
  stagger_z = 0
  call bhns_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,Bx,Bx_syms,stagger_x, &
                        stagger_y,stagger_z)
  call bhns_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,By,By_syms,stagger_x, &
                        stagger_y,stagger_z)
  call bhns_do_symmetry_gz_one_gridfunction_stagger(cctk_lsh,X,Y,Z, &
                        Symmetry,Bz,Bz_syms,stagger_x, &
                        stagger_y,stagger_z)

end subroutine bhns_compute_B_from_A

subroutine bhns_do_symmetry_gz_one_gridfunction_stagger(ext,X,Y,Z,Symmetry,gridfunc,gridfunc_syms, & 
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
end subroutine bhns_do_symmetry_gz_one_gridfunction_stagger

subroutine bhns_apply_bfield_bc_ijkmin(ext,Bx,By,Bz,EM_BC)
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

end subroutine bhns_apply_bfield_bc_ijkmin
