!--------------------------------------------------------------------------
! Reconstruction driver (slope limiter, metric & primitve facevals) : v2.0
!--------------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine reconstruction_v2(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !OTHER PARAMETERS: 
  real*8                :: dX,dY,dZ
  integer, dimension(3) :: ext
  integer, dimension(3) :: U_syms1,U_syms2,U_syms3,U_syms4
  integer, dimension(3) :: U_syms5,U_syms6,U_syms7,U_syms8
  integer, dimension(3) :: U_syms

  integer   :: i,j,k,im1,jm1,km1

  CCTK_REAL :: reduction_value
  integer   :: dummy,handle,index,ierr,cell_centering_enabled
  integer   :: PPM_PLUS, PPM, CENO,MC,SPPM
  integer   :: SCALAR,VECTORX,VECTORY,VECTORZ,BVECTORX,BVECTORY,BVECTORZ
  real*8    :: P_atm,eps_atm
  parameter(PPM_PLUS = 1, PPM = 2,CENO = 3,MC = 4, SPPM = 5)
  parameter(SCALAR = 1,VECTORX = 2,VECTORY = 3,VECTORZ = 4,BVECTORX = 5,BVECTORY = 6,BVECTORZ = 7)

  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ! Set the U_syms symmetry vectors so that all necessary ghostzone values may be filled
  call set_U_syms(Symmetry,VECTORX,Sym_Bz,U_syms1)
  call set_U_syms(Symmetry,VECTORY,Sym_Bz,U_syms2)
  call set_U_syms(Symmetry,VECTORZ,Sym_Bz,U_syms3)
  call set_U_syms(Symmetry,SCALAR,Sym_Bz,U_syms4)
  call set_U_syms(Symmetry,SCALAR,Sym_Bz,U_syms5)
  call set_U_syms(Symmetry,BVECTORX,Sym_Bz,U_syms6)
  call set_U_syms(Symmetry,BVECTORY,Sym_Bz,U_syms7)
  call set_U_syms(Symmetry,BVECTORZ,Sym_Bz,U_syms8)

  if(Symmetry.ne.4 .and. Symmetry.ne.0) then
     call cell_centered_or_not(cctk_lsh,X,Y,Z,cell_centering_enabled)
  else
     cell_centering_enabled=1
  end if

  !if(m==1 .and. iter_count==1) write(*,*) "CELL CENTERING ENABLED? (1=yes)",cell_centering_enabled

  ! compute_dU_everywhere depends on U being set EVERYWHERE, even ghostzones *(except on y!=0 plane in axisym),
  ! so here we fill in symmetry ghostzones for all the U's
  if(cell_centering_enabled==1) then
     call apply_U_syms(cctk_lsh,cctk_nghostzones,Symmetry,X,Y,Z,U_syms1,vx)
     call apply_U_syms(cctk_lsh,cctk_nghostzones,Symmetry,X,Y,Z,U_syms2,vy)
     call apply_U_syms(cctk_lsh,cctk_nghostzones,Symmetry,X,Y,Z,U_syms3,vz)

     call apply_U_syms(cctk_lsh,cctk_nghostzones,Symmetry,X,Y,Z,U_syms4,rho_b)
     call apply_U_syms(cctk_lsh,cctk_nghostzones,Symmetry,X,Y,Z,U_syms5,P)

     if(em_evolve_enable==1) then
        call apply_U_syms(cctk_lsh,cctk_nghostzones,Symmetry,X,Y,Z,U_syms6,Bx)
        call apply_U_syms(cctk_lsh,cctk_nghostzones,Symmetry,X,Y,Z,U_syms7,By)
        call apply_U_syms(cctk_lsh,cctk_nghostzones,Symmetry,X,Y,Z,U_syms8,Bz)
     end if
  else if(cell_centering_enabled==0) then
     call do_symmetry_gz_one_gridfunction(cctk_lsh,X,Y,Z,Symmetry,vx,U_syms1)
     call do_symmetry_gz_one_gridfunction(cctk_lsh,X,Y,Z,Symmetry,vy,U_syms2)
     call do_symmetry_gz_one_gridfunction(cctk_lsh,X,Y,Z,Symmetry,vz,U_syms3)

     call do_symmetry_gz_one_gridfunction(cctk_lsh,X,Y,Z,Symmetry,rho_b,U_syms4)
     call do_symmetry_gz_one_gridfunction(cctk_lsh,X,Y,Z,Symmetry,P,    U_syms5)

     if(em_evolve_enable==1) then
        call do_symmetry_gz_one_gridfunction(cctk_lsh,X,Y,Z,Symmetry,Bx,U_syms6)
        call do_symmetry_gz_one_gridfunction(cctk_lsh,X,Y,Z,Symmetry,By,U_syms7)
        call do_symmetry_gz_one_gridfunction(cctk_lsh,X,Y,Z,Symmetry,Bz,U_syms8)
     end if
  end if

  !Compute dU = U(i) - U(i-1) everywhere, where dir'n of i given by m (storing dU's in temp's)
  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms1,vx,temp1)
  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms2,vy,temp2)
  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms3,vz,temp3)

  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms4,rho_b,temp4)

  if (reconstruct_Pthermal_instead_of_P==0) then 
     call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms5,P,temp5)
  else
     ! Convert P to P_thermal and store in P_thermal
     call convert_P_to_Pthermal_cpp(cctkGH, cctk_lsh, rho_b, P, P_thermal, &
                                neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, & 
			        k_tab, gamma_tab)
     call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms5,P_thermal,temp5)
  end if

  if(em_evolve_enable==1) then
     ! Note: We don't need to reconstruct B^m for the Del Zanna CT scheme
     if (m==1 .and. constrained_transport_scheme==3) then 
        !$omp parallel do
	do k=1,ext(3)
	   do j=1,ext(2)
	      do i=1,ext(1)
		 temp6(i,j,k) = 0.d0
	      end do
	   end do
        end do
        !$omp end parallel do
     else
        call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms6,Bx,temp6)
     end if

     if (m==2 .and. constrained_transport_scheme==3) then
        !$omp parallel do
        do k=1,ext(3)
           do j=1,ext(2)
              do i=1,ext(1)
                 temp7(i,j,k) = 0.d0
              end do
           end do
        end do
        !$omp end parallel do
     else
        call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms7,By,temp7)
     end if

     if (m==3 .and. constrained_transport_scheme==3) then
        !$omp parallel do
        do k=1,ext(3)
           do j=1,ext(2)
              do i=1,ext(1)
                 temp8(i,j,k) = 0.d0
              end do
           end do
        end do
        !$omp end parallel do
     else
        call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms8,Bz,temp8)
     end if
  end if

  if(Reconstruction.ne.SPPM) then
     !Apply slope limiter to dU's
     call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp1,dvx_m)
     call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp2,dvy_m)
     call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp3,dvz_m)

     call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp4,drho_b_m)
     call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp5,dP_m)

     if(em_evolve_enable==1) then
	! Note: We don't need to reconstruct B^m for the Del Zanna CT scheme
	if (m==1 .and. constrained_transport_scheme==3) then 
	   !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
                    dBx_m(i,j,k) = 0.D0
                 end do
              end do
           end do
           !$omp end parallel do
	else
           call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp6,dBx_m)
        end if

        if (m==2 .and. constrained_transport_scheme==3) then
           !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
                    dBy_m(i,j,k) = 0.D0
                 end do
              end do
           end do
           !$omp end parallel do
        else
           call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp7,dBy_m)
        end if

        if (m==3 .and. constrained_transport_scheme==3) then
           !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
                    dBz_m(i,j,k) = 0.D0
                 end do
              end do
           end do
           !$omp end parallel do
        else
           call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp8,dBz_m)
        end if
     else
        !$omp parallel do
        do k=1,cctk_lsh(3)
           do j=1,cctk_lsh(2)
              do i=1,cctk_lsh(1)
                 dBx_m(i,j,k) = 0.D0
                 dBy_m(i,j,k) = 0.D0
                 dBz_m(i,j,k) = 0.D0
              end do
           end do
        end do
        !$omp end parallel do
     end if
  else
     !SPPM takes dU as input, applies midmod slope limiter to it, and then directly computes Ur and Ul.
     dvx_m = temp1
     dvy_m = temp2
     dvz_m = temp3
     drho_b_m = temp4
     dP_m = temp5

     if(em_evolve_enable==1) then
        dBx_m = temp6
        dBy_m = temp7
        dBz_m = temp8
     else
        !$omp parallel do
        do k=1,cctk_lsh(3)
           do j=1,cctk_lsh(2)
              do i=1,cctk_lsh(1)
                 dBx_m(i,j,k) = 0.D0
                 dBy_m(i,j,k) = 0.D0
                 dBz_m(i,j,k) = 0.D0
              end do
           end do
        end do
        !$omp end parallel do
     end if
  end if

!!$  do k=1,ext(3)
!!$     do j=1,ext(2)
!!$        do i=1,ext(1)
!!$           if(abs(Y(i,j,k)+(dY*0.5D0)).lt.dY*0.001 .and. abs(X(i,j,k)+(dX*0.5D0)).lt.dX*0.001 .and. abs(Z(i,j,k)+(dX*9.5D0)).lt.dZ*0.001) then
!!$              write(*,*) "vy:",i,j,k,vy(i,j,k),temp2(i,j,k)
!!$           end if
!!$        end do
!!$     end do
!!$  end do
  ! We must sync all the nabla_U's since Ur, Ul depend on ghostzone values:
  !  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_nablas')
  ! Symmetry ghostzones are already filled in the above routines!
  !  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_nablas')

  !!if(excision_enable == 1) then
  !!   call scalar_excision_bc(ext,X,Y,Z,drho_b_m,Symmetry,excision_zone_gf)
  !!   call scalar_excision_bc(ext,X,Y,Z,dP_m,Symmetry,excision_zone_gf)
  !!   call vector_excision_bc(ext,X,Y,Z,dvx_m,dvy_m,dvz_m,Symmetry,excision_zone_gf)
  !!   if (em_evolve_enable==1) &
  !!      call vector_excision_bc(ext,X,Y,Z,dBx_m,dBy_m,dBz_m,Symmetry,excision_zone_gf)
  !!end if

  if(Reconstruction.ne.SPPM) then
     ! Finally, compute Ur and Ul:
     call compute_Ur_Ul_general(m,cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z,Symmetry,cell_centering_enabled,U_syms1,vx,dvx_m,vxr,vxl)
     call compute_Ur_Ul_general(m,cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z,Symmetry,cell_centering_enabled,U_syms2,vy,dvy_m,vyr,vyl)
     call compute_Ur_Ul_general(m,cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z,Symmetry,cell_centering_enabled,U_syms3,vz,dvz_m,vzr,vzl)

     call compute_Ur_Ul_general(m,cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z,Symmetry,cell_centering_enabled,U_syms4,rho_b,drho_b_m,rho_br,rho_bl)
     if (reconstruct_Pthermal_instead_of_P==0) then 
        call compute_Ur_Ul_general(m,cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z,Symmetry,cell_centering_enabled,U_syms5,P,dP_m,Pr,Pl)
     else
        call compute_Ur_Ul_general(m,cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z,Symmetry,cell_centering_enabled,U_syms5,P_thermal,dP_m,Pr,Pl)
     end if
  else
     ! Finally, compute Ur and Ul:
     call compute_Ur_Ul_SPPM(m,cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z,Symmetry,cell_centering_enabled,U_syms1,vx,dvx_m,vxr,vxl)
     call compute_Ur_Ul_SPPM(m,cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z,Symmetry,cell_centering_enabled,U_syms2,vy,dvy_m,vyr,vyl)
     call compute_Ur_Ul_SPPM(m,cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z,Symmetry,cell_centering_enabled,U_syms3,vz,dvz_m,vzr,vzl)

     call compute_Ur_Ul_SPPM(m,cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z,Symmetry,cell_centering_enabled,U_syms4,rho_b,drho_b_m,rho_br,rho_bl)
     if (reconstruct_Pthermal_instead_of_P==0) then
        call compute_Ur_Ul_SPPM(m,cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z,Symmetry,cell_centering_enabled,U_syms5,P,dP_m,Pr,Pl)
     else
        call compute_Ur_Ul_SPPM(m,cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z,Symmetry,cell_centering_enabled,U_syms5,P_thermal,dP_m,Pr,Pl)
     end if
  end if

  if (reconstruct_Pthermal_instead_of_P==1) then 
     call convert_Pthermal_to_P_cpp(cctkGH, cctk_lsh, &
                                rho_br, rho_bl, Pr, Pl, neos, &
                                rho_tab, P_tab, eps_tab, k_tab, gamma_tab)
  end if

  if(em_evolve_enable==1) then
     if(Reconstruction.ne.SPPM) then
	! Note: We don't need to reconstruct B^m for the Del Zanna CT scheme
	if (m==1 .and. constrained_transport_scheme==3) then 
	   !$omp parallel do private (im1)
	   do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
	            im1 = max(i-1,1)
	            Bxr(i,j,k) = Bx_stagger(im1,j,k)
		    Bxl(i,j,k) = Bx_stagger(im1,j,k)
                 end do
              end do
           end do
	   !$omp end parallel do
	else
           call compute_Ur_Ul_general(m,cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z,Symmetry,cell_centering_enabled,U_syms6,Bx,dBx_m,Bxr,Bxl)
        end if

        if (m==2 .and. constrained_transport_scheme==3) then
           !$omp parallel do private (jm1)
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
	            jm1 = max(j-1,1)
                    Byr(i,j,k) = By_stagger(i,jm1,k)
                    Byl(i,j,k) = By_stagger(i,jm1,k)
                 end do
              end do
           end do
           !$omp end parallel do
        else
           call compute_Ur_Ul_general(m,cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z,Symmetry,cell_centering_enabled,U_syms7,By,dBy_m,Byr,Byl)
	end if

        if (m==3 .and. constrained_transport_scheme==3) then
           !$omp parallel do private (km1)
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
	            km1 = max(k-1,1)
                    Bzr(i,j,k) = Bz_stagger(i,j,km1)
                    Bzl(i,j,k) = Bz_stagger(i,j,km1)
                 end do
              end do
           end do
           !$omp end parallel do
        else
           call compute_Ur_Ul_general(m,cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z,Symmetry,cell_centering_enabled,U_syms8,Bz,dBz_m,Bzr,Bzl)
	end if
     else
        ! Note: We don't need to reconstruct B^m for the Del Zanna CT scheme
        if (m==1 .and. constrained_transport_scheme==3) then
           !$omp parallel do private (im1)
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
	     	    im1 = max(i-1,1)
                    Bxr(i,j,k) = Bx_stagger(im1,j,k)
                    Bxl(i,j,k) = Bx_stagger(im1,j,k)
                 end do
              end do
           end do
           !$omp end parallel do
        else
           call compute_Ur_Ul_SPPM(m,cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z,Symmetry,cell_centering_enabled,U_syms6,Bx,dBx_m,Bxr,Bxl)
        end if

        if (m==2 .and. constrained_transport_scheme==3) then
           !$omp parallel do private (jm1)
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
		    jm1 = max(j-1,1)
                    Byr(i,j,k) = By_stagger(i,jm1,k)
                    Byl(i,j,k) = By_stagger(i,jm1,k)
                 end do
              end do
           end do
           !$omp end parallel do
        else
           call compute_Ur_Ul_SPPM(m,cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z,Symmetry,cell_centering_enabled,U_syms7,By,dBy_m,Byr,Byl)
        end if

        if (m==3 .and. constrained_transport_scheme==3) then
           !$omp parallel do private (km1)
           do k=2,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
	 	    km1 = max(k-1,1)
                    Bzr(i,j,k) = Bz_stagger(i,j,km1)
                    Bzl(i,j,k) = Bz_stagger(i,j,km1)
                 end do
              end do
           end do
           !$omp end parallel do
        else
           call compute_Ur_Ul_SPPM(m,cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z,Symmetry,cell_centering_enabled,U_syms8,Bz,dBz_m,Bzr,Bzl)
	end if
     end if

  else
     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              Bxr(i,j,k) = 0.D0
              Bxl(i,j,k) = 0.D0
              Byr(i,j,k) = 0.D0
              Byl(i,j,k) = 0.D0
              Bzr(i,j,k) = 0.D0
              Bzl(i,j,k) = 0.D0
           end do
        end do
     end do
     !$omp end parallel do
  end if

  if(Reconstruction==SPPM) then
     if(m==1) then
        !$omp parallel do
        do k=1,cctk_lsh(3)
           do j=1,cctk_lsh(2)
              do i=cctk_lsh(1),2,-1
                 rho_bl(i,j,k) = rho_bl(i-1,j,k)
                 Pl(i,j,k) = Pl(i-1,j,k)
                 vxl(i,j,k) = vxl(i-1,j,k)
                 vyl(i,j,k) = vyl(i-1,j,k)
                 vzl(i,j,k) = vzl(i-1,j,k)
              end do
           end do
        end do
        !$omp end parallel do

        rho_bl(1,:,:) = rho_bl(2,:,:)
        Pl(1,:,:) = Pl(2,:,:)
        vxl(1,:,:) = vxl(2,:,:)
        vyl(1,:,:) = vyl(2,:,:)
        vzl(1,:,:) = vzl(2,:,:)
     else if(m==2) then
        !$omp parallel do
        do k=1,cctk_lsh(3)
           do j=cctk_lsh(2),2,-1
              do i=1,cctk_lsh(1)
                 rho_bl(i,j,k) = rho_bl(i,j-1,k)
                 Pl(i,j,k) = Pl(i,j-1,k)
                 vxl(i,j,k) = vxl(i,j-1,k)
                 vyl(i,j,k) = vyl(i,j-1,k)
                 vzl(i,j,k) = vzl(i,j-1,k)
              end do
           end do
        end do
        !$omp end parallel do

        rho_bl(:,1,:) = rho_bl(:,2,:)
        Pl(:,1,:) = Pl(:,2,:)
        vxl(:,1,:) = vxl(:,2,:)
        vyl(:,1,:) = vyl(:,2,:)
        vzl(:,1,:) = vzl(:,2,:)
     else if(m==3) then
        !$omp parallel do
        do k=cctk_lsh(3),2,-1
           do j=1,cctk_lsh(2)
              do i=1,cctk_lsh(1)
                 rho_bl(i,j,k) = rho_bl(i,j,k-1)
                 Pl(i,j,k) = Pl(i,j,k-1)
                 vxl(i,j,k) = vxl(i,j,k-1)
                 vyl(i,j,k) = vyl(i,j,k-1)
                 vzl(i,j,k) = vzl(i,j,k-1)
              end do
           end do
        end do
        !$omp end parallel do

        rho_bl(:,:,1) = rho_bl(:,:,2)
        Pl(:,:,1) = Pl(:,:,2)
        vxl(:,:,1) = vxl(:,:,2)
        vyl(:,:,1) = vyl(:,:,2)
        vzl(:,:,1) = vzl(:,:,2)
     end if

     !     call sppm_shift_hydro_ul(ext,Z,rho_br,rho_bl,Pr,Pl, &
     !          vxr,vxl,vyr,vyl,vzr,vzl,m,Symmetry)
     !     call sppm_shift_emfields_ul(ext,Z,Bxr,Bxl,Byr,Byl, &
     !          Bzr,Bzl,m,Symmetry,Sym_Bz)
     !call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')
     !call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')
     !     call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_B_quantities')
     !call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_B_quantities')
  end if

  !!call compute_pcold_epscold(rho_b_atm,P_atm,eps_atm, & 
  !!                  neos,ergo_star, ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
  !!!$omp parallel do
  !!do k=1,cctk_lsh(3)
  !!   do j=1,cctk_lsh(2)
  !!      do i=1,cctk_lsh(1)
  !!         Pl(i,j,k) = max( Pl(i,j,k), P_atm)
  !!         rho_br(i,j,k) = max( rho_br(i,j,k), rho_b_atm)
  !!         Pr(i,j,k) = max( Pr(i,j,k), P_atm)
  !!         rho_bl(i,j,k) = max( rho_bl(i,j,k), rho_b_atm)
  !!      end do
  !!   end do
  !!end do
  !!!$omp end parallel do
 
!  where (Pl .lt. 0.d0 .or. rho_bl .lt. 0.d0)
!     Pl = 0.d0
!     rho_bl = 0.d0
!  end where
!
!  where (Pr .lt. 0.d0 .or. rho_br .lt. 0.d0)
!     Pr = 0.d0
!     rho_br = 0.d0
!  end where

  !call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_rho_br_rho_bl')

  ! Symmetry ghostzones are already filled in the above routines!
  !call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_rho_br_rho_bl')

  ! Average to find metric quantities at the interfaces

  ! IF PPM RECONSTRUCTION IS ENABLED, THIS FUNCTION (driver_reconstruction_v2()) IS NOT CALLED.  THEREFORE WE COMMENT THE FOLLOWING SECTION OUT.  Soon we will simply delete these lines of code.
!!$  if (Reconstruction == PPM) then 
!!$     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms4,lapm1,lapm1_f)
!!$     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms1,shiftx,shiftx_f)
!!$     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms2,shifty,shifty_f)
!!$     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms3,shiftz,shiftz_f)
!!$!!!!!$     call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms4,lapm1,lapm1_f)
!!$!!!!!$     call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms1,shiftx,shiftx_f)
!!$!!!!!$     call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms2,shifty,shifty_f)
!!$!!!!!$     call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms3,shiftz,shiftz_f)
!!$
!!$     U_syms = U_syms1*U_syms1
!!$     !call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gxx,gxx_f)
!!$     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gxx,gxx_f)
!!$     U_syms = U_syms1*U_syms2
!!$     !call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gxy,gxy_f)
!!$     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gxy,gxy_f)
!!$     U_syms = U_syms1*U_syms3
!!$     !call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gxz,gxz_f)
!!$     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gxz,gxz_f)
!!$     U_syms = U_syms2*U_syms2
!!$     !call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gyy,gyy_f)
!!$     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gyy,gyy_f)
!!$     U_syms = U_syms2*U_syms3
!!$     !call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gyz,gyz_f) 
!!$     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gyz,gyz_f) 
!!$     U_syms = U_syms3*U_syms3
!!$     !call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gzz,gzz_f)
!!$     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gzz,gzz_f)
!!$
!!$     !call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms4,phi,phi_f)
!!$     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms4,phi,phi_f)
!!$
!!$     !call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gupxx,gupxx_f)
!!$     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gupxx,gupxx_f)
!!$     !call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gupyy,gupyy_f)
!!$     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gupyy,gupyy_f)
!!$     !call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gupzz,gupzz_f)
!!$     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gupzz,gupzz_f)
!!$  else 
     call compute_face_avg(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,lapm1,lapm1_f)
     call compute_face_avg(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,shiftx,shiftx_f)
     call compute_face_avg(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,shifty,shifty_f)
     call compute_face_avg(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,shiftz,shiftz_f)

     call compute_face_avg(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,gxx,gxx_f)
     call compute_face_avg(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,gxy,gxy_f)
     call compute_face_avg(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,gxz,gxz_f)
     call compute_face_avg(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,gyy,gyy_f)
     call compute_face_avg(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,gyz,gyz_f)
     call compute_face_avg(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,gzz,gzz_f)

     call compute_face_avg(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,phi,phi_f)

     call compute_face_avg(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,gupxx,gupxx_f)
     call compute_face_avg(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,gupyy,gupyy_f)
     call compute_face_avg(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,gupzz,gupzz_f)
!!$  end if

  !  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_metric_facevals')
  ! Symmetry ghostzones are already filled in the above routines!
  !call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_metric_facevals')

  !!if(excision_enable == 1) then
  !!   call scalar_excision_bc(ext,X,Y,Z, &
  !!        lapm1_f,Symmetry,excision_zone_gf)
  !!   call scalar_excision_bc(ext,X,Y,Z, &
  !!        phi_f,Symmetry,excision_zone_gf)
  !!   call vector_excision_bc(ext,X,Y,Z, &
  !!        shiftx_f,shifty_f,shiftz_f, &
  !!        Symmetry,excision_zone_gf)
  !!   call vector_excision_bc(ext,X,Y,Z, &
  !!        gxx_f,gxy_f,gxz_f, &
  !!        Symmetry,excision_zone_gf)
  !!   call vector_excision_bc(ext,X,Y,Z, &
  !!        gyy_f,gyz_f,gzz_f, &
  !!        Symmetry,excision_zone_gf)
  !!   call vector_excision_bc(ext,X,Y,Z, &
  !!        gupxx_f,gupyy_f,gupzz_f, &
  !!        Symmetry,excision_zone_gf)
  !!end if

  ! Implement PPM reconstruction.  IF PPM RECONSTRUCTION IS ENABLED, THIS FUNCTION (driver_reconstruction_v2()) IS NOT CALLED.  THEREFORE WE COMMENT THIS SECTION OUT.  Soon we will simply delete these lines of code.
!!$  if (Reconstruction.eq.PPM) then 
!!$
!!$     ! Compute center differences d rho and d^2 rho everywhere 
!!$     ! and store them in temp1 and temp2
!!$     call find_centderivs_f_cpp(cctkGH,ext,cctk_nghostzones,rho_b,temp1,temp2,m,U_syms4, & 
!!$          Symmetry,cell_centering_enabled)
!!$     !call find_centderivs_f(ext,cctk_nghostzones,rho_b,temp1,temp2,m,U_syms4, & 
!!$     !           Symmetry,cell_centering_enabled)
!!$
!!$     ! Copy rho_br and rho_bl computed from MC to temp3 and temp4 
!!$     !$omp parallel do
!!$     do k=1,cctk_lsh(3)
!!$        do j=1,cctk_lsh(2)
!!$           do i=1,cctk_lsh(1)
!!$              temp3(i,j,k) = rho_br(i,j,k)
!!$	      temp4(i,j,k) = rho_bl(i,j,k)
!!$           end do
!!$        end do
!!$     end do
!!$     !$omp end parallel do     
!!$
!!$     call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,rho_b,drho_b_m, & 
!!$          rho_br,rho_bl,m,U_syms4, &
!!$          Symmetry,cell_centering_enabled)
!!$     if (reconstruct_Pthermal_instead_of_P==0) then 
!!$        call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,P,dP_m, &
!!$             Pr,Pl,m,U_syms5, &
!!$             Symmetry,cell_centering_enabled)
!!$     else
!!$        call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,P_thermal,dP_m, &
!!$             Pr,Pl,m,U_syms5, &
!!$             Symmetry,cell_centering_enabled)
!!$     end if
!!$     call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,vx,dvx_m, &
!!$          vxr,vxl,m,U_syms1, &
!!$          Symmetry,cell_centering_enabled)
!!$     call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,vy,dvy_m, &
!!$          vyr,vyl,m,U_syms2, &
!!$          Symmetry,cell_centering_enabled)
!!$     call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,vz,dvz_m, &
!!$          vzr,vzl,m,U_syms3, &
!!$          Symmetry,cell_centering_enabled)
!!$
!!$     if (em_evolve_enable==1) then
!!$        ! Note: We don't need to reconstruct B^m for the Del Zanna CT scheme
!!$        if (m==1 .and. constrained_transport_scheme==3) then
!!$           !$omp parallel do private (im1)
!!$           do k=1,cctk_lsh(3)
!!$              do j=1,cctk_lsh(2)
!!$                 do i=1,cctk_lsh(1)
!!$		    im1 = max(i-1,1)
!!$                    Bxr(i,j,k) = Bx_stagger(im1,j,k)
!!$                    Bxl(i,j,k) = Bx_stagger(im1,j,k)
!!$                 end do
!!$              end do
!!$           end do
!!$           !$omp end parallel do
!!$        else
!!$           call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,Bx,dBx_m, &
!!$             Bxr,Bxl,m,U_syms6, &
!!$             Symmetry,cell_centering_enabled)
!!$	end if
!!$
!!$        if (m==2 .and. constrained_transport_scheme==3) then
!!$           !$omp parallel do private (jm1)
!!$           do k=1,cctk_lsh(3)
!!$              do j=1,cctk_lsh(2)
!!$                 do i=1,cctk_lsh(1)
!!$		    jm1 = max(j-1,1)
!!$                    Byr(i,j,k) = By_stagger(i,jm1,k)
!!$                    Byl(i,j,k) = By_stagger(i,jm1,k)
!!$                 end do
!!$              end do
!!$           end do
!!$           !$omp end parallel do
!!$        else
!!$           call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,By,dBy_m, &
!!$             Byr,Byl,m,U_syms7, &
!!$             Symmetry,cell_centering_enabled)
!!$	end if
!!$
!!$        if (m==3 .and. constrained_transport_scheme==3) then
!!$           !$omp parallel do private (km1)
!!$           do k=1,cctk_lsh(3)
!!$              do j=1,cctk_lsh(2)
!!$                 do i=1,cctk_lsh(1)
!!$		    km1 = max(k-1,1)
!!$                    Bzr(i,j,k) = Bz_stagger(i,j,km1)
!!$                    Bzl(i,j,k) = Bz_stagger(i,j,km1)
!!$                 end do
!!$              end do
!!$           end do
!!$           !$omp end parallel do
!!$        else
!!$           call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,Bz,dBz_m, &
!!$             Bzr,Bzl,m,U_syms8, &
!!$             Symmetry,cell_centering_enabled)
!!$	end if
!!$
!!$     end if
!!$!!!!!$     call ppm_find_face_vals_v2(ext,cctk_nghostzones,rho_b,drho_b_m, & 
!!$!!!!!$          rho_br,rho_bl,m,U_syms4, &
!!$!!!!!$          Symmetry,cell_centering_enabled)
!!$!!!!!$     call ppm_find_face_vals_v2(ext,cctk_nghostzones,P,dP_m, &
!!$!!!!!$          Pr,Pl,m,U_syms5, &
!!$!!!!!$          Symmetry,cell_centering_enabled)
!!$!!!!!$     call ppm_find_face_vals_v2(ext,cctk_nghostzones,vx,dvx_m, &
!!$!!!!!$          vxr,vxl,m,U_syms1, &
!!$!!!!!$          Symmetry,cell_centering_enabled)
!!$!!!!!$     call ppm_find_face_vals_v2(ext,cctk_nghostzones,vy,dvy_m, &
!!$!!!!!$          vyr,vyl,m,U_syms2, &
!!$!!!!!$          Symmetry,cell_centering_enabled)
!!$!!!!!$     call ppm_find_face_vals_v2(ext,cctk_nghostzones,vz,dvz_m, &
!!$!!!!!$          vzr,vzl,m,U_syms3, &
!!$!!!!!$          Symmetry,cell_centering_enabled)
!!$!!!!!$
!!$!!!!!$     if (em_evolve_enable==1) then
!!$!!!!!$        call ppm_find_face_vals_v2(ext,cctk_nghostzones,Bx,dBx_m, &
!!$!!!!!$             Bxr,Bxl,m,U_syms6, &
!!$!!!!!$             Symmetry,cell_centering_enabled)
!!$!!!!!$        call ppm_find_face_vals_v2(ext,cctk_nghostzones,By,dBy_m, &
!!$!!!!!$             Byr,Byl,m,U_syms7, &
!!$!!!!!$             Symmetry,cell_centering_enabled)
!!$!!!!!$        call ppm_find_face_vals_v2(ext,cctk_nghostzones,Bz,dBz_m, &
!!$!!!!!$             Bzr,Bzl,m,U_syms8, &
!!$!!!!!$             Symmetry,cell_centering_enabled)
!!$!!!!!$     end if
!!$
!!$     ! Steppen rho_b
!!$     call ppm_steepen_rho_v2_cpp(cctkGH,ext,rho_b,temp1,temp2,rho_br,rho_bl, & 
!!$          temp3,temp4,P,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th, & 
!!$          neos,m,Symmetry)
!!$!!!!!$     call ppm_steepen_rho_v2(ext,rho_b,temp1,temp2,rho_br,rho_bl, & 
!!$!!!!!$          temp3,temp4,P,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th, & 
!!$!!!!!$          neos,m,Symmetry)
!!$
!!$     ! Impose monotonicity 
!!$     call ppm_monotonize_v2_cpp(cctkGH,ext,rho_b,rho_br,rho_bl,Symmetry)
!!$     if (reconstruct_Pthermal_instead_of_P==0) then 
!!$        call ppm_monotonize_v2_cpp(cctkGH,ext,P,Pr,Pl,Symmetry)
!!$     else
!!$        call ppm_monotonize_v2_cpp(cctkGH,ext,P_thermal,Pr,Pl,Symmetry)
!!$     end if
!!$     call ppm_monotonize_v2_cpp(cctkGH,ext,vx,vxr,vxl,Symmetry)
!!$     call ppm_monotonize_v2_cpp(cctkGH,ext,vy,vyr,vyl,Symmetry)
!!$     call ppm_monotonize_v2_cpp(cctkGH,ext,vz,vzr,vzl,Symmetry)
!!$
!!$     if (em_evolve_enable==1) then
!!$	if (m .ne. 1 .or. constrained_transport_scheme .ne. 3) &
!!$           call ppm_monotonize_v2_cpp(cctkGH,ext,Bx,Bxr,Bxl,Symmetry)
!!$
!!$        if (m .ne. 2 .or. constrained_transport_scheme .ne. 3) &
!!$           call ppm_monotonize_v2_cpp(cctkGH,ext,By,Byr,Byl,Symmetry)
!!$
!!$        if (m .ne. 3 .or. constrained_transport_scheme .ne. 3) &
!!$           call ppm_monotonize_v2_cpp(cctkGH,ext,Bz,Bzr,Bzl,Symmetry)
!!$     end if
!!$!!!!!$     call ppm_monotonize_v2(ext,rho_b,rho_br,rho_bl,Symmetry)
!!$!!!!!$     call ppm_monotonize_v2(ext,P,Pr,Pl,Symmetry)
!!$!!!!!$     call ppm_monotonize_v2(ext,vx,vxr,vxl,Symmetry)
!!$!!!!!$     call ppm_monotonize_v2(ext,vy,vyr,vyl,Symmetry)
!!$!!!!!$     call ppm_monotonize_v2(ext,vz,vzr,vzl,Symmetry)
!!$!!!!!$
!!$!!!!!$     if (em_evolve_enable==1) then
!!$!!!!!$        call ppm_monotonize_v2(ext,Bx,Bxr,Bxl,Symmetry)
!!$!!!!!$        call ppm_monotonize_v2(ext,By,Byr,Byl,Symmetry)
!!$!!!!!$        call ppm_monotonize_v2(ext,Bz,Bzr,Bzl,Symmetry)
!!$!!!!!$     end if
!!$
!!$     ! Compute ftilde and store it to temp5
!!$     !BUG:
!!$     !     call ppm_ftilde_v2_cpp(cctkGH,ext,cctk_nghostzones,temp5,P,vx,vy,vz,U_syms5, & 
!!$     !          U_syms1,U_syms2,U_syms3, &
!!$     !          m,Symmetry,cell_centering_enabled)
!!$     call ppm_ftilde_v2(ext,cctk_nghostzones,temp5,P,vx,vy,vz,U_syms5, & 
!!$          U_syms1,U_syms2,U_syms3, &
!!$          m,Symmetry,cell_centering_enabled)
!!$
!!$     !     write(*,*) "ftilde:",temp5(10,10,:)
!!$     ! Apply flattening
!!$     call ppm_flatten_v2_cpp(cctkGH,ext,temp5,rho_b,rho_br,rho_bl,P,P_thermal,Pr,Pl, &
!!$          vx,vxr,vxl,vy,vyr,vyl,vz,vzr,vzl,Bx,Bxr,Bxl,By,Byr,Byl,Bz,Bzr,Bzl, &
!!$          em_evolve_enable,m,Symmetry,constrained_transport_scheme, & 
!!$          reconstruct_Pthermal_instead_of_P)
!!$!!!!!$     call ppm_flatten_v2(ext,temp5,rho_b,rho_br,rho_bl,P,Pr,Pl, &
!!$!!!!!$          vx,vxr,vxl,vy,vyr,vyl,vz,vzr,vzl,Bx,Bxr,Bxl,By,Byr,Byl,Bz,Bzr,Bzl, &
!!$!!!!!$          em_evolve_enable,m,Symmetry)
!!$
!!$  if (reconstruct_Pthermal_instead_of_P==1) then
!!$     call convert_Pthermal_to_P_cpp(cctkGH, cctk_lsh, &
!!$                                rho_br, rho_bl, Pr, Pl, neos, &
!!$                                rho_tab, P_tab, eps_tab, k_tab, gamma_tab)
!!$  end if
!!$
!!$     ! Shift the stencil
!!$     call ppm_shift_v2_cpp(cctkGH,ext,rho_br,rho_bl,temp6,m,Symmetry)
!!$     call ppm_shift_v2_cpp(cctkGH,ext,Pr,Pl,temp6,m,Symmetry)
!!$     call ppm_shift_v2_cpp(cctkGH,ext,vxr,vxl,temp6,m,Symmetry)
!!$     call ppm_shift_v2_cpp(cctkGH,ext,vyr,vyl,temp6,m,Symmetry)
!!$     call ppm_shift_v2_cpp(cctkGH,ext,vzr,vzl,temp6,m,Symmetry)
!!$
!!$     if (em_evolve_enable==1) then
!!$	if (m .ne. 1 .or. constrained_transport_scheme .ne. 3) & 
!!$           call ppm_shift_v2_cpp(cctkGH,ext,Bxr,Bxl,temp6,m,Symmetry)
!!$
!!$        if (m .ne. 2 .or. constrained_transport_scheme .ne. 3) &
!!$           call ppm_shift_v2_cpp(cctkGH,ext,Byr,Byl,temp6,m,Symmetry)
!!$
!!$        if (m .ne. 3 .or. constrained_transport_scheme .ne. 3) &
!!$          call ppm_shift_v2_cpp(cctkGH,ext,Bzr,Bzl,temp6,m,Symmetry)
!!$     end if
!!$!!!!!$     call ppm_shift_v2(ext,rho_br,rho_bl,temp6,m,Symmetry)
!!$!!!!!$     call ppm_shift_v2(ext,Pr,Pl,temp6,m,Symmetry)
!!$!!!!!$     call ppm_shift_v2(ext,vxr,vxl,temp6,m,Symmetry)
!!$!!!!!$     call ppm_shift_v2(ext,vyr,vyl,temp6,m,Symmetry)
!!$!!!!!$     call ppm_shift_v2(ext,vzr,vzl,temp6,m,Symmetry)
!!$!!!!!$
!!$!!!!!$     if (em_evolve_enable==1) then
!!$!!!!!$        call ppm_shift_v2(ext,Bxr,Bxl,temp6,m,Symmetry)
!!$!!!!!$        call ppm_shift_v2(ext,Byr,Byl,temp6,m,Symmetry)
!!$!!!!!$        call ppm_shift_v2(ext,Bzr,Bzl,temp6,m,Symmetry)
!!$!!!!!$     end if
!!$  end if

  if (Reconstruction.eq.PPM_PLUS .or. Reconstruction.eq.CENO) then 
     write(*,*) 'Reconstruction method not supported yet.' 
     stop
  end if

  !  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')
  ! Symmetry ghostzones are already filled in the above routines!
  !call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')
  !  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_B_quantities')
  ! Symmetry ghostzones are already filled in the above routines!
  !call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_B_quantities')

end subroutine reconstruction_v2
