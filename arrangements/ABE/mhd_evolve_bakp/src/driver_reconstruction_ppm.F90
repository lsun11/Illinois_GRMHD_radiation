!--------------------------------------------------------------------------
! Reconstruction driver (slope limiter, metric & primitve facevals) : v2.0
!--------------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine reconstruction_ppm_part1(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !OTHER PARAMETERS: 
  real*8                :: dX,dY,dZ
  integer, dimension(3) :: ext
  integer, dimension(3) :: U_syms1,U_syms2,U_syms3,U_syms4
  integer, dimension(3) :: U_syms5,U_syms6,U_syms7,U_syms8
  integer, dimension(3) :: U_syms9,U_syms10,U_syms11,U_syms12, U_syms13
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


if(1==0) then

  do k=1, ext(3)
     do j=1, ext(2)
        do i=1, ext(1)

if (i==27.and.j==24.and.k==16) then
   write(*,*) "1.In driver_reconstruction_ppm_part1, at the start, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==25.and.j==14.and.k==19) then
   write(*,*) "2.In driver_reconstruction_ppm_part1, at the start, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==49.and.j==2.and.k==2) then
    write(*,*) " In driver_reconstruction_ppm_part1, at the start, tau(i,j,k) is ", tau(i,j,k)
end if

if (i==49.and.j==2.and.k==2) then
    write(*,*) " In driver_reconstruction_ppm_part1, at the start, tau_rhs(i,j,k) is ", tau_rhs(i,j,k)
end if


if (i==49.and.j==2.and.k==2) then
    write(*,*) " In driver_reconstruction_ppm_part1, at the start, tau_rad(i,j,k) is ", tau_rad(i,j,k)
end if

if (i==49.and.j==2.and.k==2) then
    write(*,*) " In driver_reconstruction_ppm_part1, at the start, tau_rad_rhs(i,j,k) is ", tau_rad_rhs(i,j,k)
end if


         end do
    end do
  end do

end if


  ! Set the U_syms symmetry vectors so that all necessary ghostzone values may be filled
  call set_U_syms(Symmetry,VECTORX,Sym_Bz,U_syms1)
  call set_U_syms(Symmetry,VECTORY,Sym_Bz,U_syms2)
  call set_U_syms(Symmetry,VECTORZ,Sym_Bz,U_syms3)
  call set_U_syms(Symmetry,SCALAR,Sym_Bz,U_syms4)
  call set_U_syms(Symmetry,SCALAR,Sym_Bz,U_syms5)
  call set_U_syms(Symmetry,BVECTORX,Sym_Bz,U_syms6)
  call set_U_syms(Symmetry,BVECTORY,Sym_Bz,U_syms7)
  call set_U_syms(Symmetry,BVECTORZ,Sym_Bz,U_syms8)
  call set_U_syms(Symmetry,SCALAR,Sym_Bz,U_syms9)
  call set_U_syms(Symmetry,VECTORX,Sym_Bz,U_syms10)
  call set_U_syms(Symmetry,VECTORY,Sym_Bz,U_syms11)
  call set_U_syms(Symmetry,VECTORZ,Sym_Bz,U_syms12)

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
     
     if(rad_evolve_enable==1) then
        call apply_U_syms(cctk_lsh,cctk_nghostzones,Symmetry,X,Y,Z,U_syms9,E_rad)
	
	call apply_U_syms(cctk_lsh,cctk_nghostzones,Symmetry,X,Y,Z,U_syms10,F_radx)
     	call apply_U_syms(cctk_lsh,cctk_nghostzones,Symmetry,X,Y,Z,U_syms11,F_rady)
     	call apply_U_syms(cctk_lsh,cctk_nghostzones,Symmetry,X,Y,Z,U_syms12,F_radz)

	call apply_U_syms(cctk_lsh,cctk_nghostzones,Symmetry,X,Y,Z,U_syms13,P_rad)
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

     if(rad_evolve_enable==1) then
        call do_symmetry_gz_one_gridfunction(cctk_lsh,X,Y,Z,Symmetry,E_rad,U_syms9)

        call do_symmetry_gz_one_gridfunction(cctk_lsh,X,Y,Z,Symmetry,F_radx,U_syms10)
        call do_symmetry_gz_one_gridfunction(cctk_lsh,X,Y,Z,Symmetry,F_rady,U_syms11)
        call do_symmetry_gz_one_gridfunction(cctk_lsh,X,Y,Z,Symmetry,F_radz,U_syms12)
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

   if(rad_evolve_enable==1) then
     ! Note: temp9 is taken. use temp13 for E_rad.
     call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms9,E_rad,temp13)
     call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms10,F_radx,temp10)
     call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms11,F_rady,temp11)
     call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms12,F_radz,temp12)
     call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms12,P_rad,temp14)
   end if

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

  if (rad_evolve_enable == 1) then
     call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp1,tau_rad_flux)
     call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp1,S_radx_flux)
     call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp2,S_rady_flux)
     call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp3,S_radz_flux)
  end if 

  ! Compute rho_br_mc and rho_bl_mc:
  call compute_Ur_Ul_general(m,cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z,Symmetry,cell_centering_enabled,U_syms4,rho_b,drho_b_m,rho_br,rho_bl)

  if(1==0) then
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
  end if

  if(1==1) then
     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms4,lapm1,lapm1_f)
     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms1,shiftx,shiftx_f)
     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms2,shifty,shifty_f)
     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms3,shiftz,shiftz_f)

     U_syms = U_syms1*U_syms1
     !call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gxx,gxx_f)
     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gxx,gxx_f)
     U_syms = U_syms1*U_syms2
     !call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gxy,gxy_f)
     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gxy,gxy_f)
     U_syms = U_syms1*U_syms3
     !call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gxz,gxz_f)
     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gxz,gxz_f)
     U_syms = U_syms2*U_syms2
     !call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gyy,gyy_f)
     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gyy,gyy_f)
     U_syms = U_syms2*U_syms3
     !call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gyz,gyz_f) 
     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gyz,gyz_f) 
     U_syms = U_syms3*U_syms3
     !call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gzz,gzz_f)
     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gzz,gzz_f)

     !call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms4,phi,phi_f)
     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms4,phi,phi_f)

     !call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gupxx,gupxx_f)
     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gupxx,gupxx_f)
     !call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gupyy,gupyy_f)
     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gupyy,gupyy_f)
     !call compute_face_avg_ppm(m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gupzz,gupzz_f)
     call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms,gupzz,gupzz_f)
  end if

  ! Compute center differences d rho and d^2 rho everywhere 
  ! and store them in temp1 and temp2, respectively
  call find_centderivs_f_cpp(cctkGH,ext,cctk_nghostzones,rho_b,temp1,temp2,m,U_syms4, & 
       Symmetry,cell_centering_enabled)

  ! Copy rho_br and rho_bl computed from MC to temp3 and temp4 
  !$omp parallel do
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
           temp3(i,j,k) = rho_br(i,j,k)
           temp4(i,j,k) = rho_bl(i,j,k)
        end do
     end do
  end do
  !$omp end parallel do     

  call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,rho_b,drho_b_m, & 
       rho_br,rho_bl,m,U_syms4, &
       Symmetry,cell_centering_enabled)
  if (reconstruct_Pthermal_instead_of_P==0) then 
     call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,P,dP_m, &
          Pr,Pl,m,U_syms5, &
          Symmetry,cell_centering_enabled)
  else
     call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,P_thermal,dP_m, &
          Pr,Pl,m,U_syms5, &
          Symmetry,cell_centering_enabled)
  end if
  call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,vx,dvx_m, &
       vxr,vxl,m,U_syms1, &
       Symmetry,cell_centering_enabled)
  call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,vy,dvy_m, &
       vyr,vyl,m,U_syms2, &
       Symmetry,cell_centering_enabled)
  call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,vz,dvz_m, &
       vzr,vzl,m,U_syms3, &
       Symmetry,cell_centering_enabled)

  if (em_evolve_enable==1) then
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
        call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,Bx,dBx_m, &
          Bxr,Bxl,m,U_syms6, &
          Symmetry,cell_centering_enabled)
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
        call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,By,dBy_m, &
          Byr,Byl,m,U_syms7, &
          Symmetry,cell_centering_enabled)
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
        call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,Bz,dBz_m, &
          Bzr,Bzl,m,U_syms8, &
          Symmetry,cell_centering_enabled)
     end if

  else
     !$omp parallel do
      do k=1,cctk_lsh(3)
         do j=1,cctk_lsh(2)
            do i=1,cctk_lsh(1)
	       Bxl(i,j,k) = 0.d0
               Bxr(i,j,k) = 0.d0
               Byl(i,j,k) = 0.d0
               Byr(i,j,k) = 0.d0
               Bzl(i,j,k) = 0.d0
               Bzr(i,j,k) = 0.d0
            end do
         end do
      end do
     !$omp end parallel do
  end if
  
  if (rad_evolve_enable == 1) then
     call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,E_rad,tau_rad_flux, &
          E_radr,E_radl,m,U_syms9, &
          Symmetry,cell_centering_enabled)

     call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,F_radx,S_radx_flux, &
          F_radxr,F_radxl,m,U_syms10, &
          Symmetry,cell_centering_enabled)
     call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,F_rady,S_rady_flux, &
          F_radyr,F_radyl,m,U_syms11, &
          Symmetry,cell_centering_enabled) 
     call ppm_find_face_vals_v2_cpp(cctkGH,ext,cctk_nghostzones,F_radz,S_radx_flux, &
          F_radzr,F_radzl,m,U_syms12, &
          Symmetry,cell_centering_enabled)
  end if




  ! Steppen rho_b
  call ppm_steepen_rho_v2_cpp(cctkGH,ext,rho_b,temp1,temp2,rho_br,rho_bl, & 
       temp3,temp4,P,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th, & 
       neos,m,Symmetry)

  ! Move the monotonicity procedure to the last step
  !! Impose monotonicity 
  !!call ppm_monotonize_v2_cpp(cctkGH,ext,rho_b,rho_br,rho_bl,Symmetry)
  !!if (reconstruct_Pthermal_instead_of_P==0) then 
  !!   call ppm_monotonize_v2_cpp(cctkGH,ext,P,Pr,Pl,Symmetry)
  !!else
  !!   call ppm_monotonize_v2_cpp(cctkGH,ext,P_thermal,Pr,Pl,Symmetry)
  !!end if
  !!call ppm_monotonize_v2_cpp(cctkGH,ext,vx,vxr,vxl,Symmetry)
  !!call ppm_monotonize_v2_cpp(cctkGH,ext,vy,vyr,vyl,Symmetry)
  !!call ppm_monotonize_v2_cpp(cctkGH,ext,vz,vzr,vzl,Symmetry)

  !!if (em_evolve_enable==1) then
  !!   if (m .ne. 1 .or. constrained_transport_scheme .ne. 3) &
  !!      call ppm_monotonize_v2_cpp(cctkGH,ext,Bx,Bxr,Bxl,Symmetry)

  !!   if (m .ne. 2 .or. constrained_transport_scheme .ne. 3) &
  !!      call ppm_monotonize_v2_cpp(cctkGH,ext,By,Byr,Byl,Symmetry)

  !!   if (m .ne. 3 .or. constrained_transport_scheme .ne. 3) &
  !!      call ppm_monotonize_v2_cpp(cctkGH,ext,Bz,Bzr,Bzl,Symmetry)
  !!end if

  ! Compute ftilde and store it to temp6
  !BUG:
  !     call ppm_ftilde_v2_cpp(cctkGH,ext,cctk_nghostzones,temp6,P,vx,vy,vz,U_syms5, &
  !          U_syms1,U_syms2,U_syms3, &
  !          m,Symmetry,cell_centering_enabled)
  call ppm_ftilde_v2(ext,cctk_nghostzones,temp6,P,vx,vy,vz,U_syms5, &
       U_syms1,U_syms2,U_syms3, &
       m,Symmetry,cell_centering_enabled)

 ! Need to do a sync before going to the rest of the reconstruction

end subroutine reconstruction_ppm_part1

subroutine reconstruction_ppm_part2(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !OTHER PARAMETERS:
  integer, dimension(3) :: ext
  integer   :: i,j,k,im1,jm1,km1

  real*8    :: P_atm,eps_atm

  real*8                :: horizdirn_x,horizdirn_y,horizdirn_z,AH_radius_minimum


  ext = cctk_lsh


if(1==0) then

  do k=1, ext(3)
     do j=1, ext(2)
        do i=1, ext(1)

if (i==27.and.j==24.and.k==16) then
   write(*,*) "1.In driver_reconstruction_ppm_part2, at the very end, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==25.and.j==14.and.k==19) then
   write(*,*) "2.In driver_reconstruction_ppm_part2, at the very end, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==49.and.j==2.and.k==2) then
    write(*,*) " In driver_reconstruction_ppm_part2, at the very end, tau(i,j,k) is ", tau(i,j,k)
end if

if (i==49.and.j==2.and.k==2) then
    write(*,*) " In driver_reconstruction_ppm_part2, at the very end, tau_rhs(i,j,k) is ", tau_rhs(i,j,k)
end if



if (i==49.and.j==2.and.k==2) then
    write(*,*) " In driver_reconstruction_ppm_part2, at the very end, tau_rad(i,j,k) is ", tau_rad(i,j,k)
end if

if (i==49.and.j==2.and.k==2) then
    write(*,*) " In driver_reconstruction_ppm_part2, at the very end, tau_rad_rhs(i,j,k) is ", tau_rad_rhs(i,j,k)
end if



         end do
    end do
  end do

end if



  ! Compute f and store it to temp9
  if(use_std_PPM_instead_of_efficient_3_gz_version.eq.0) then
     !$omp parallel do
     do k=1,ext(3)
        do j=1,ext(2)
           do i=1,ext(1)
              temp9(i,j,k) = temp6(i,j,k)
           end do
        end do
     end do
     !$omp end parallel do
  else
     call ppm_f_cpp(cctkGH,ext,m,temp6,temp9,P)
  end if
  ! Apply flattening
  !!call ppm_flatten_v2_cpp(cctkGH,ext,temp6,rho_b,rho_br,rho_bl,P,P_thermal,Pr,Pl, &
  !!     vx,vxr,vxl,vy,vyr,vyl,vz,vzr,vzl,Bx,Bxr,Bxl,By,Byr,Byl,Bz,Bzr,Bzl, &
  !!     em_evolve_enable,m,Symmetry,constrained_transport_scheme, &
  !!     reconstruct_Pthermal_instead_of_P)

  call ppm_flatten_1var_cpp(cctkGH,ext,m,temp9,rho_b,rho_br,rho_bl)
  if (reconstruct_Pthermal_instead_of_P==1) then 
     call ppm_flatten_1var_cpp(cctkGH,ext,m,temp9,P_thermal,Pr,Pl)
  else
     call ppm_flatten_1var_cpp(cctkGH,ext,m,temp9,P,Pr,Pl)
  end if
  call ppm_flatten_1var_cpp(cctkGH,ext,m,temp9,vx,vxr,vxl)
  call ppm_flatten_1var_cpp(cctkGH,ext,m,temp9,vy,vyr,vyl)
  call ppm_flatten_1var_cpp(cctkGH,ext,m,temp9,vz,vzr,vzl)

  if (em_evolve_enable==1) then
     if (m .ne. 1 .or. constrained_transport_scheme .ne. 3) &
	call ppm_flatten_1var_cpp(cctkGH,ext,m,temp9,Bx,Bxr,Bxl)

     if (m .ne. 2 .or. constrained_transport_scheme .ne. 3) &
        call ppm_flatten_1var_cpp(cctkGH,ext,m,temp9,By,Byr,Byl)

     if (m .ne. 3 .or. constrained_transport_scheme .ne. 3) &
        call ppm_flatten_1var_cpp(cctkGH,ext,m,temp9,Bz,Bzr,Bzl)
  end if

  if (rad_evolve_enable==1) then
      call ppm_flatten_1var_cpp(cctkGH,ext,m,temp9,E_rad,E_radr,E_radl)
      call ppm_flatten_1var_cpp(cctkGH,ext,m,temp9,F_radx,F_radxr,F_radxl)
      call ppm_flatten_1var_cpp(cctkGH,ext,m,temp9,F_rady,F_radyr,F_radyl)
      call ppm_flatten_1var_cpp(cctkGH,ext,m,temp9,F_radz,F_radzr,F_radzl)
      call ppm_flatten_1var_cpp(cctkGH,ext,m,temp9,P_rad,P_radr,P_radl)
  end if


  ! Impose monotonicity
  call ppm_monotonize_v2_cpp(cctkGH,ext,rho_b,rho_br,rho_bl,Symmetry)
  if (reconstruct_Pthermal_instead_of_P==0) then
     call ppm_monotonize_v2_cpp(cctkGH,ext,P,Pr,Pl,Symmetry)
  else
     call ppm_monotonize_v2_cpp(cctkGH,ext,P_thermal,Pr,Pl,Symmetry)
  end if
  call ppm_monotonize_v2_cpp(cctkGH,ext,vx,vxr,vxl,Symmetry)
  call ppm_monotonize_v2_cpp(cctkGH,ext,vy,vyr,vyl,Symmetry)
  call ppm_monotonize_v2_cpp(cctkGH,ext,vz,vzr,vzl,Symmetry)

  if (em_evolve_enable==1) then
     if (m .ne. 1 .or. constrained_transport_scheme .ne. 3) &
        call ppm_monotonize_v2_cpp(cctkGH,ext,Bx,Bxr,Bxl,Symmetry)

     if (m .ne. 2 .or. constrained_transport_scheme .ne. 3) &
        call ppm_monotonize_v2_cpp(cctkGH,ext,By,Byr,Byl,Symmetry)

     if (m .ne. 3 .or. constrained_transport_scheme .ne. 3) &
        call ppm_monotonize_v2_cpp(cctkGH,ext,Bz,Bzr,Bzl,Symmetry)
  end if



if(1==0) then

  do k=1, ext(3)
     do j=1, ext(2)
        do i=1, ext(1)

if (i==27.and.j==24.and.k==16) then
   write(*,*) "1.In driver_reconstruction_ppm_part2, at the middle, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==25.and.j==14.and.k==19) then
   write(*,*) "2.In driver_reconstruction_ppm_part2, at the middle, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==49.and.j==2.and.k==2) then
    write(*,*) " In driver_reconstruction_ppm_part2, at the middle, tau(i,j,k) is ", tau(i,j,k)
end if

if (i==49.and.j==2.and.k==2) then
    write(*,*) " In driver_reconstruction_ppm_part2, at the middle, tau_rhs(i,j,k) is ", tau_rhs(i,j,k)
end if


if (i==49.and.j==2.and.k==2) then
    write(*,*) " In driver_reconstruction_ppm_part2, at the middle, tau_rad(i,j,k) is ", tau_rad(i,j,k)
end if

if (i==49.and.j==2.and.k==2) then
    write(*,*) " In driver_reconstruction_ppm_part2, at the middle, tau_rad_rhs(i,j,k) is ", tau_rad_rhs(i,j,k)
end if




         end do
    end do
  end do

end if



 if (rad_evolve_enable==1) then
    call  ppm_monotonize_v2_cpp(cctkGH,ext,E_rad,E_radr,E_radl,Symmetry)

    call  ppm_monotonize_v2_cpp(cctkGH,ext,F_radx,F_radxr,F_radxl,Symmetry)
    call  ppm_monotonize_v2_cpp(cctkGH,ext,F_rady,F_radyr,F_radyl,Symmetry)
    call  ppm_monotonize_v2_cpp(cctkGH,ext,F_radz,F_radzr,F_radzl,Symmetry)
 end if


  if (reconstruct_Pthermal_instead_of_P==1) then
     call convert_Pthermal_to_P_cpp(cctkGH, cctk_lsh, &
                                rho_br, rho_bl, Pr, Pl, neos, &
                                rho_tab, P_tab, eps_tab, k_tab, gamma_tab)
  end if

  ! Shift the stencil
  call ppm_shift_v2_cpp(cctkGH,ext,rho_br,rho_bl,temp6,m,Symmetry)
  call ppm_shift_v2_cpp(cctkGH,ext,Pr,Pl,temp6,m,Symmetry)
  call ppm_shift_v2_cpp(cctkGH,ext,vxr,vxl,temp6,m,Symmetry)
  call ppm_shift_v2_cpp(cctkGH,ext,vyr,vyl,temp6,m,Symmetry)
  call ppm_shift_v2_cpp(cctkGH,ext,vzr,vzl,temp6,m,Symmetry)

  if (em_evolve_enable==1) then
     if (m .ne. 1 .or. constrained_transport_scheme .ne. 3) &
        call ppm_shift_v2_cpp(cctkGH,ext,Bxr,Bxl,temp6,m,Symmetry)

     if (m .ne. 2 .or. constrained_transport_scheme .ne. 3) &
      call ppm_shift_v2_cpp(cctkGH,ext,Byr,Byl,temp6,m,Symmetry)

     if (m .ne. 3 .or. constrained_transport_scheme .ne. 3) &
       call ppm_shift_v2_cpp(cctkGH,ext,Bzr,Bzl,temp6,m,Symmetry)
  end if

 if (rad_evolve_enable==1) then
    call ppm_shift_v2_cpp(cctkGH,ext,E_radr,E_radl,temp6,m,Symmetry)

    call ppm_shift_v2_cpp(cctkGH,ext,F_radxr,F_radxl,temp6,m,Symmetry)
    call ppm_shift_v2_cpp(cctkGH,ext,F_radyr,F_radyl,temp6,m,Symmetry)
    call ppm_shift_v2_cpp(cctkGH,ext,F_radzr,F_radzl,temp6,m,Symmetry)

    call ppm_shift_v2_cpp(cctkGH,ext,P_radr,P_radl,temp6,m,Symmetry)
 end if



if(1==0) then

  do k=1, ext(3)
     do j=1, ext(2)
        do i=1, ext(1)

if (i==27.and.j==24.and.k==16) then
   write(*,*) "1.In driver_reconstruction_ppm_part2, at the very end, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==25.and.j==14.and.k==19) then
   write(*,*) "2.In driver_reconstruction_ppm_part2, at the very end, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==49.and.j==2.and.k==2) then
    write(*,*) " In driver_reconstruction_ppm_part2, at the very end, tau(i,j,k) is ", tau(i,j,k)
end if

if (i==49.and.j==2.and.k==2) then
    write(*,*) " In driver_reconstruction_ppm_part2, at the very end, tau_rhs(i,j,k) is ", tau_rhs(i,j,k)
end if


if (i==49.and.j==2.and.k==2) then
    write(*,*) " In driver_reconstruction_ppm_part2, at the very end, tau_rad(i,j,k) is ", tau_rad(i,j,k)
end if

if (i==49.and.j==2.and.k==2) then
    write(*,*) " In driver_reconstruction_ppm_part2, at the very end, tau_rad_rhs(i,j,k) is ", tau_rad_rhs(i,j,k)
end if


         end do
    end do
  end do

end if




end subroutine reconstruction_ppm_part2
