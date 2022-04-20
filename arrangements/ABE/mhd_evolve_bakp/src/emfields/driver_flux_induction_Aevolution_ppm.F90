!---------------------------------------------------------------
! Compute Ai_rhs (fluxes for the induction equation) using HRSC 
!---------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!=====================================================================================
! In this step, we reconstruct the flux F3 = vx Bytilde - vy Bxtilde 
! at (i+1/2,j+1/2,k). The fluxes F3rr, F3rl, F3lr and F3ll (see 
! Del Zanna, Bucciantini & Londrillo A&A 400, 397 (2003) Eq 44 and 
! earlier for details, F3=E_z) are stored in fxy, fxz, fyx and fyz, 
! respectively. The values of Bxtilde^r-Bxtilde^l and 
! Bytilde^r-Bytilde^l are stored in fzx and fzy, respectively. 
! These variables will be used at beginning of step 2. Values of 
! the x-direction of cmax and cmin are stored in temp10 and MONOPOLE, 
! to be used in steps 2 and 3.
!=====================================================================================
subroutine flux_induction_Aevolution_ppm_step1(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  integer  ::  dummy,i,j,k,cell_centering_enabled,ip1,jp1
  integer, dimension(3) :: U_syms1,U_syms2,U_syms3,U_syms4,U_syms5,U_syms6,U_syms7
  integer, parameter :: VECTORX = 2,VECTORY = 3, VECTORZ = 4
  integer, parameter :: BVECTORX = 5,BVECTORY = 6, SCALAR = 1

  ! Set the U_syms symmetry vectors so that all necessary ghostzone values may be filled
  call set_U_syms(Symmetry,VECTORX,Sym_Bz,U_syms1)
  call set_U_syms(Symmetry,VECTORY,Sym_Bz,U_syms2)
  call set_U_syms(Symmetry,VECTORZ,Sym_Bz,U_syms3)
  call set_U_syms(Symmetry,SCALAR,Sym_Bz,U_syms4)
  call set_U_syms(Symmetry,SCALAR,Sym_Bz,U_syms5)
  call set_U_syms(Symmetry,BVECTORX,Sym_Bz,U_syms6)
  call set_U_syms(Symmetry,BVECTORY,Sym_Bz,U_syms7)

  if(Symmetry.ne.4 .and. Symmetry.ne.0) then
     call cell_centered_or_not(cctk_lsh,X,Y,Z,cell_centering_enabled)
  else
     cell_centering_enabled=1
  end if

  !Zach says: I don't think Carpet respects CCTK_SyncGroup()'s.  At best, these are inert.
  !call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')
  !!call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')

  ! Reconstruct By at (i+1/2,j+1/2,k)

  ! By_stagger stores values of By on a semi-staggered grid, at points
  ! (i,j+1/2,k).  Here, we compute dU for By, at points (i-1/2,j+1/2,k):
  m = 1
  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms7,By_stagger,temp1)
  !Apply slope limiter
  call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp1,dBy_m)

  ! Compute Byr (stored in ftzx) and Byl (stored in ftzy)
  call ppm_find_face_vals_v2_cpp(cctkGH,cctk_lsh,cctk_nghostzones,By_stagger,dBy_m, &
       ftzx,ftzy,m,U_syms7,Symmetry,cell_centering_enabled)
  call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,By_stagger,ftzx,ftzy,Symmetry)
  ! Apply flattening: 
  !Note that f has been computed in the hydro step and stored in temp9
  if (em_evolve_enable .ne. 2) then 
     call ppm_flatten_1var_cpp(cctkGH,cctk_lsh,m,temp9,By_stagger,ftzx,ftzy)
  end if

  ! Compute e^(6 phi) at (i+1/2,j+1/2,k) and stored in temp11
  m = 2
  call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms4,phi_f,temp1)
  !omp parallel do private (ip1,jp1)
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
           ip1 = min(i+1,cctk_lsh(1))
           jp1 = min(j+1,cctk_lsh(2))
           temp11(i,j,k) = exp(6.d0*temp1(ip1,jp1,k))
        end do
     end do
  end do
  !omp end parallel do

  ! Reconstruct Bx, vx and vy at (i+1/2,j+1/2,k).

  ! Note that we have already reconstructed vx and vy along the x-direction, 
  ! at (i-1/2,j,k). That result is stored in v{x,y}{r,l}.  Bx_stagger data
  ! are on a semi-staggered grid, at (i+1/2,j,k).  Here, we compute dU for 
  ! vx and vy at (i-1/2,j-1/2,k), and Bx at (i+1/2,j-1/2,k), and store the 
  ! result in temp{1,2,6,7,5}:
  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms1,vxr,temp1)
  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms1,vxl,temp2)
  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms2,vyr,temp6)
  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms2,vyl,temp7)
  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms6,Bx_stagger,temp5)

  !Next, we apply the MC slope limiter, storing the result in temp{3,4,8,9},dBx_m
  call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp1,temp3)
  call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp2,temp4)
  call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp6,temp8)
  call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp7,temp9)
  call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp5,dBx_m)

  ! Finally, from the slope-limited dU data, we compute Ul and Ur,
  ! using a host of fij's and ftij's as temporary storage:
  !========================
  ! variable | where stored
  !==========|=============
  ! vxrr     | fxy
  ! vxrl     | fxz
  ! vxlr     | fyx
  ! vxll     | fyz
  ! Bxr      | fzx
  ! Bxl      | fzy
  ! vyrr     | ftxy
  ! vyrl     | ftxz
  ! vylr     | ftyx
  ! vyll     | ftyz
  call ppm_find_face_vals_v2_cpp(cctkGH,cctk_lsh,cctk_nghostzones,vxr,temp3, &
       fxy,fxz,m,U_syms1,Symmetry,cell_centering_enabled)  
  call ppm_find_face_vals_v2_cpp(cctkGH,cctk_lsh,cctk_nghostzones,vxl,temp4, &
       fyx,fyz,m,U_syms1,Symmetry,cell_centering_enabled)
  call ppm_find_face_vals_v2_cpp(cctkGH,cctk_lsh,cctk_nghostzones,vyr,temp8, &
       ftxy,ftxz,m,U_syms2,Symmetry,cell_centering_enabled)
  call ppm_find_face_vals_v2_cpp(cctkGH,cctk_lsh,cctk_nghostzones,vyl,temp9, &
       ftyx,ftyz,m,U_syms2,Symmetry,cell_centering_enabled)
  call ppm_find_face_vals_v2_cpp(cctkGH,cctk_lsh,cctk_nghostzones,Bx_stagger,dBx_m, &
       fzx,fzy,m,U_syms6,Symmetry,cell_centering_enabled)

  ! Monotonicity will be imposed last
  !! Impose monotonicity
  !!call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vxr,fxy,fxz,Symmetry)
  !!call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vxl,fyx,fyz,Symmetry)
  !!call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vyr,ftxy,ftxz,Symmetry)
  !!call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vyl,ftyx,ftyz,Symmetry)  
  !!call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,Bx_stagger,fzx,fzy,Symmetry)

  ! Compute ftilde and store it to temp1
  if (em_evolve_enable .ne. 2) then 
     call ppm_ftilde_v2(cctk_lsh,cctk_nghostzones,temp1,P,vx,vy,vz,U_syms5, &
          U_syms1,U_syms2,U_syms3, &
          m,Symmetry,cell_centering_enabled)
  end if

  ! Need to do a sync before going to the rest of the reconstruction

end subroutine flux_induction_Aevolution_ppm_step1

subroutine flux_induction_Aevolution_ppm_step1p5(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !OTHER PARAMETERS:
  integer   :: i,j,k,ip1,jp1,kp1

  real*8                :: horizdirn_x,horizdirn_y,horizdirn_z,AH_radius_minimum


  !========================
  ! variable | where stored
  !==========|=============
  ! vxrr     | fxy
  ! vxrl     | fxz
  ! vxlr     | fyx
  ! vxll     | fyz
  ! Bxr      | fzx
  ! Bxl      | fzy
  ! vyrr     | ftxy
  ! vyrl     | ftxz
  ! vylr     | ftyx
  ! vyll     | ftyz


  ! Compute f and store it to temp2
  if (em_evolve_enable .ne. 2) then 
     ! Compute f and store it to temp2
     if(use_std_PPM_instead_of_efficient_3_gz_version.eq.0) then
        !$omp parallel do
        do k=1,cctk_lsh(3)
           do j=1,cctk_lsh(2)
              do i=1,cctk_lsh(1)
                 temp2(i,j,k) = temp1(i,j,k)
              end do
           end do
        end do
        !$omp end parallel do
     else
        call ppm_f_cpp(cctkGH,cctk_lsh,m,temp1,temp2,P)
     end if

     ! Apply flattening
     call ppm_flatten_1var_cpp(cctkGH,cctk_lsh,m,temp2,vxr,fxy,fxz)
     call ppm_flatten_1var_cpp(cctkGH,cctk_lsh,m,temp2,vxl,fyx,fyz)
     call ppm_flatten_1var_cpp(cctkGH,cctk_lsh,m,temp2,Bx_stagger,fzx,fzy)
     call ppm_flatten_1var_cpp(cctkGH,cctk_lsh,m,temp2,vyr,ftxy,ftxz)
     call ppm_flatten_1var_cpp(cctkGH,cctk_lsh,m,temp2,vyl,ftyx,ftyz)
  end if

  ! Impose monotonicity
  call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vxr,fxy,fxz,Symmetry)
  call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vxl,fyx,fyz,Symmetry)
  call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vyr,ftxy,ftxz,Symmetry)
  call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vyl,ftyx,ftyz,Symmetry)
  call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,Bx_stagger,fzx,fzy,Symmetry)

  !! Shift the stencil, use temp3 as a temporary storage
  !!call ppm_shift_v2_cpp(cctkGH,cctk_lsh,fxy,fxz,temp3,m,Symmetry)
  !!call ppm_shift_v2_cpp(cctkGH,cctk_lsh,fyx,fyz,temp3,m,Symmetry)
  !!call ppm_shift_v2_cpp(cctkGH,cctk_lsh,fzx,fzy,temp3,m,Symmetry)
  !!call ppm_shift_v2_cpp(cctkGH,cctk_lsh,ftxy,ftxz,temp3,m,Symmetry)
  !!call ppm_shift_v2_cpp(cctkGH,cctk_lsh,ftyx,ftyz,temp3,m,Symmetry)

!!! Shift the stencil for Byr (stored in ftzx) and Byl (stored in ftzy)
  !! m = 1
  !!call ppm_shift_v2_cpp(cctkGH,cctk_lsh,ftzx,ftzy,temp3,m,Symmetry)

  ! Note: I am skipping the ppm_shift step since we are going to shift the 
  ! stencil again below. I am combining the two shifts into one. 
  ! By skipping the ppm_shift step, the variables change its meaning. 
  ! Here is the new map: 
  !=======================================
  ! gf    |   variable  |   staggering    |
  ! ===== | =========== | =============== |
  ! fxy   |    vxrl     | (i-1/2,j+1/2,k) |
  ! fxz   |    vxrr     | (i-1/2,j-1/2,k) |
  ! fyx   |    vxll     | (i-1/2,j+1/2,k) |
  ! fyz   |    vxlr     | (i-1/2,j-1/2,k) |
  ! fzx   |    Bxl      | (i+1/2,j+1/2,k) |
  ! fzy   |    Bxr      | (i+1/2,j-1/2,k) |
  ! ftxy  |    vyrl     | (i-1/2,j+1/2,k) |
  ! ftxz  |    vyrr     | (i-1/2,j-1/2,k) |
  ! ftyx  |    vyll     | (i-1/2,j+1/2,k) |
  ! ftyz  |    vylr     | (i-1/2,j-1/2,k) |
  ! ftzx  |    Byl      | (i+1/2,j+1/2,k) |
  ! ftzy  |    Byr      | (i-1/2,j+1/2,k) |


  ! Shift the stencil and change B^i_rl to tildeB^i_rl
  ! (Note: We want the reconstructed variables to be at 
  !   (i+1/2,j+1/2,k), but all the reconstructed v's 
  !   and B's are at the position listed above. 
  !   Hence we need to shift the stencil.)
  !========================
  ! variable | where stored
  !==========|=============
  ! vxrr     | temp1
  ! vxrl     | temp2
  ! vxlr     | temp3
  ! vxll     | temp4
  ! vyrr     | temp5
  ! vyrl     | temp6
  ! vylr     | temp7
  ! vyll     | temp8
  !omp parallel do private(ip1,jp1)
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
	   ip1 = min(i+1,cctk_lsh(1))
	   jp1 = min(j+1,cctk_lsh(2))
           temp2(i,j,k) = fxy(ip1,j,k)
           temp1(i,j,k) = fxz(ip1,jp1,k)
           temp4(i,j,k) = fyx(ip1,j,k)
           temp3(i,j,k) = fyz(ip1,jp1,k)
           temp6(i,j,k) = ftxy(ip1,j,k)
           temp5(i,j,k) = ftxz(ip1,jp1,k)
           temp8(i,j,k) = ftyx(ip1,j,k)
           temp7(i,j,k) = ftyz(ip1,jp1,k)
           Bxl(i,j,k) = fzx(i,j,k)*temp11(i,j,k)
           Bxr(i,j,k) = fzy(i,jp1,k)*temp11(i,j,k)
           Byl(i,j,k) = ftzx(i,j,k)*temp11(i,j,k)
           Byr(i,j,k) = ftzy(ip1,j,k)*temp11(i,j,k)
        end do
     end do
  end do
  !omp end parallel do
  ! Switch to first-order scheme inside the AH
  !!if(num_BHs.gt.0.5D0) then
  !!   !     Get horizon radius in direction where it is likely to be minimized:
  !!   horizdirn_x = 0.D0
  !!   horizdirn_y = 0.D0
  !!   horizdirn_z = 100000.D0
  !!   call get_ah_radius_in_dirn(cctkGH,horizdirn_x,horizdirn_y,horizdirn_z,AH_radius_minimum);
  !!   if(AH_radius_minimum.lt.0.D0) then
  !!      AH_radius_minimum = bh_radius_z(1)
  !!      write(*,*) "WARNING:  Couldn't find horizon, so using last known BH radius & position to limit B^2 inside the horizon.  Radius=",AH_radius_minimum,num_BHs,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1)
  !!   end if

  !!   !        write(*,*) "Radius=",AH_radius_minimum,num_BHs,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1)


  !!   ! Hard-coded 0.8 safety factor here:
  !!   AH_radius_minimum = 0.8D0 * AH_radius_minimum

  !!   ! Switch to first-order scheme inside AH
  !!   !omp parallel do private(ip1,jp1)
  !!   do k=1,cctk_lsh(3)
  !!      do j=1,cctk_lsh(2)
  !!         do i=1,cctk_lsh(1)
  !!            if(sqrt((x(i,j,k)-bh_posn_x(1))**2 + (y(i,j,k)-bh_posn_y(1))**2 + (z(i,j,k)-bh_posn_z(1))**2) .lt. AH_radius_minimum) then
  !!               ip1 = min(i+1,cctk_lsh(1))
  !!               jp1 = min(j+1,cctk_lsh(2))
  !!               temp1(i,j,k) = vx(ip1,jp1,k)
  !!               temp2(i,j,k) = vx(ip1,j,k)
  !!               temp3(i,j,k) = vx(i,jp1,k)
  !!               temp4(i,j,k) = vx(i,j,k)
  !!               temp5(i,j,k) = vy(ip1,jp1,k)
  !!               temp6(i,j,k) = vy(ip1,j,k)
  !!               temp7(i,j,k) = vy(i,jp1,k)
  !!               temp8(i,j,k) = vy(i,j,k)
  !!               Bxl(i,j,k) = Bx_stagger(i,j,k)*temp11(i,j,k)
  !!               Bxr(i,j,k) = Bx_stagger(i,jp1,k)*temp11(i,j,k)
  !!               Byl(i,j,k) = By_stagger(i,j,k)*temp11(i,j,k)
  !!               Byr(i,j,k) = By_stagger(ip1,j,k)*temp11(i,j,k)
  !!            end if
  !!         end do
  !!      end do
  !!   end do
  !!   !omp end parallel do
  !!end if
  ! **************

  ! Compute the 4 reconstructed F3=(v cross Btilde)^z values at (i+1/2,j+1/2,k).
  ! Also compute Bxr - Bxl (stored in fzx) and Byr - Byl (stored in fzy). 
  ! (Note: the reconstruction order is first x and then y)
  ! Variables are stored as follows::
  !========================
  ! variable  | where stored
  !===========|=============
  ! F3rr      | fxy
  ! F3rl      | fxz
  ! F3lr      | fyx
  ! F3ll      | fyz
  ! F3ll      | fyz
  ! Bxr - Bxl | fzx
  ! Byr - Byl | fzy
  !omp parallel do
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
	   fxy(i,j,k) = temp1(i,j,k)*Byr(i,j,k) - temp5(i,j,k)*Bxr(i,j,k)
           fxz(i,j,k) = temp2(i,j,k)*Byr(i,j,k) - temp6(i,j,k)*Bxl(i,j,k)
	   fyx(i,j,k) = temp3(i,j,k)*Byl(i,j,k) - temp7(i,j,k)*Bxr(i,j,k)
	   fyz(i,j,k) = temp4(i,j,k)*Byl(i,j,k) - temp8(i,j,k)*Bxl(i,j,k)
	   fzx(i,j,k) = Bxr(i,j,k) - Bxl(i,j,k)
	   fzy(i,j,k) = Byr(i,j,k) - Byl(i,j,k)
        end do
     end do
  end do
  !omp end parallel do

  !Zach says: I don't think Carpet respects CCTK_SyncGroup()'s.  At best, these are inert.
  !!call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::em_fijs')
  !!call CartSymGN(dummy,cctkGH,'mhd_evolve::em_fijs') 

  ! In the previous hydro routine, cmax and cmin were computed 
  ! in the x-direction, at (i-1/2,j,k).  Here we shift the 
  ! stencil so that cmax/cmin at (i+1/2,j,k) is stored at temp10/MONOPOLE,
  ! respectively.  This way we can later compute Az_rhs (and later
  ! Ay_rhs) using the 2D HLL flux (Del Zanna, Bucciantini & Londrillo
  ! A&A 400, 397 (2003), Eq. (44).  Note that cmax/cmin (\alpha^{\pm}
  ! as defined in the paper) is at a slightly DIFFERENT point than that 
  ! described in the Del Zanna et al paper (e.g., (i+1/2,j,k) instead of
  ! (i+1/2,j+1/2,k) for F3).  Yuk Tung discussed this point with Shibata,
  ! who found that the effect is negligible.
  ! Finally, store cmax and cmin in temp10 and MONOPOLE (to be used in steps 2 and 3)
  !omp parallel do private (ip1)
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
	   ip1 = min(i+1,cctk_lsh(1))
	   temp10(i,j,k) = cmax(ip1,j,k)
	   MONOPOLE(i,j,k) = cmin(ip1,j,k)
        end do
     end do
  end do
  !omp end parallel do

end subroutine flux_induction_Aevolution_ppm_step1p5

!=====================================================================================
! In this step, first Az_rhs is computed.  Next, we reconstruct 
! the flux F1 = vy Bztilde - vz Bytilde at (i,j+1/2,k+1/2).
! The fluxes F1rr, F1rl, F1lr and F1ll (see
! Del Zanna, Bucciantini & Londrillo A&A 400, 397 (2003) for details)
! are stored in fxy, fxz, fyx and fyz, respectively. The values of
! Bytilde^r-Bytilde^l and Bztilde^r-Bztilde^l are stored in
! fzx and fzy, respectively. These variables will be used at
! beginning of step 3. Values of the y-direction of cmax
! and cmin are stored in ftxy and ftxz, to be used in
! step 3.
!=====================================================================================
subroutine flux_induction_Aevolution_ppm_step2(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer  ::  dummy,i,j,k, cell_centering_enabled,jp1,kp1
  integer, dimension(3) :: U_syms1,U_syms2,U_syms3,U_syms4,U_syms5
  integer, dimension(3) :: U_syms6,U_syms7,U_syms8
  integer, parameter :: VECTORX = 2,VECTORY = 3,VECTORZ = 4
  integer, parameter :: BVECTORX = 5, BVECTORY = 6,BVECTORZ = 7, SCALAR = 1

  ! In the previous hydro routine, cmax and cmin were computed 
  ! in the y-direction, at (i,j-1/2,k).  Here we shift the 
  ! stencil so that cmax/cmin at (i,j+1/2,k) is stored at ftxy/ftxz,
  ! respectively.  This way we can compute Az_rhs (and later Ax_rhs)
  ! using the 2D HLL flux (Del Zanna, Bucciantini & Londrillo A&A 
  ! 400, 397 (2003), Eq. (44).  Note that cmax/cmin (\alpha^{\pm}
  ! as defined in the paper) is at a slightly DIFFERENT point than that 
  ! described in the Del Zanna et al paper (e.g., (i+1/2,j,k) instead of
  ! (i+1/2,j+1/2,k) for F3).  Yuk Tung discussed this point with Shibata,
  ! who found that the effect is negligible.
  !omp parallel do private (jp1)
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
	   jp1 = min(j+1,cctk_lsh(2))
           ftxy(i,j,k) = cmax(i,jp1,k)
           ftxz(i,j,k) = cmin(i,jp1,k)
        end do
     end do
  end do
  !omp end parallel do

  ! Compute Az_rhs=-Ez=-F3 using the 2D HLL flux.
  ! See Del Zanna, Bucciantini & Londrillo A&A 400, 397 (2003), Eq. (44)
  ! 
  ! In the previous step, 
  ! the fluxes F3=(v cross Btilde)^z at (i+1/2,j+1/2,k)
  ! have been stored in fxy (F3rr), fxz (F3rl), fyx (F3lr) and fyz (F3ll). 
  ! Also, Bxr - Bxl is stored in fzx, Byr - Byl is stored in fzy; 
  ! cmax, cmin in the x-direction are stored in temp10 and MONOPOLE.
  ! cmax and cmin in the y-direction have just been computed and stored 
  ! in ftxy and ftxz, albeit at slightly different points than prescribed
  ! in Del Zanna et al. (see comment above).  In summary:
  !========================
  ! variable  | where stored
  !===========|=============
  ! F3rr      | fxy
  ! F3rl      | fxz
  ! F3lr      | fyx
  ! F3ll      | fyz
  ! Bxr - Bxl | fzx
  ! Byr - Byl | fzy
  call flux_hll2d_cpp(cctkGH,cctk_lsh,temp10, &
       MONOPOLE, ftxy, ftxz, &
       fxy, fxz, fyx, fyz, fzx, fzy, & 
       Bztilde_or_Az_rhs,use_central_scheme_instead_of_hll)

  !Zach says: I don't think Carpet respects CCTK_SyncGroup()'s.  At best, these are inert.
  !call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')
  !!call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')

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

  ! Reconstruct Bz at (i,j+1/2,k+1/2)
  ! Bz_stagger stores values of Bz on a semi-staggered grid, at points
  ! (i,j,k+1/2).  Here, we compute dU for Bz, at points (i,j+1/2,k-1/2):
  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms8,Bz_stagger,temp1)
  !Apply slope limiter
  call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp1,dBz_m)

  ! Compute Bzr (stored in ftzx) and Bzl (stored in ftzy)
  call ppm_find_face_vals_v2_cpp(cctkGH,cctk_lsh,cctk_nghostzones,Bz_stagger,dBz_m, &
       ftzx,ftzy,m,U_syms8,Symmetry,cell_centering_enabled)
  ! Apply flattening:
  !Note that f has been computed in the hydro step and stored in temp9
  if (em_evolve_enable .ne. 2) then 
     call ppm_flatten_1var_cpp(cctkGH,cctk_lsh,m,temp9,Bz_stagger,ftzx,ftzy)
  end if
  call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,Bz_stagger,ftzx,ftzy,Symmetry)

  ! Compute e^(6 phi) at (i,j+1/2,k+1/2) and stored in temp11
  m = 3
  call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms4,phi_f,temp1)
  !omp parallel do private (jp1,kp1)
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
	   jp1 = min(j+1,cctk_lsh(2))
	   kp1 = min(k+1,cctk_lsh(3))
	   temp11(i,j,k) = exp(6.d0*temp1(i,jp1,kp1) )
        end do
     end do
  end do
  !omp end parallel do

  ! Reconstruct By, vy and vz at (i,j+1/2,k+1/2).
  ! Note that we have already reconstructed vy and vz along the y-direction, 
  ! at (i,j-1/2,k). That result is stored in v{y,z}{r,l}.  By_stagger data
  ! are on a semi-staggered grid, at (i,j+1/2,k).  Here, we compute dU for 
  ! vy and vz at (i,j-1/2,k-1/2), and By at (i,j-1/2,k+1/2), and store the 
  ! result in temp{1,2,6,7,5}:
  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms2,vyr,temp1)
  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms2,vyl,temp2)
  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms3,vzr,temp6)
  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms3,vzl,temp7)
  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms7,By_stagger,temp5)

  !Next, we apply the MC slope limiter, storing the result in temp{3,4,8,9},dBy_m
  call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp1,temp3)
  call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp2,temp4)
  call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp6,temp8)
  call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp7,temp9)
  call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp5,dBy_m)

  ! Finally, from the slope-limited dU data, we compute Ul and Ur,
  ! using a host of fij's and ftij's as temporary storage:
  !========================
  ! variable | where stored
  !==========|=============
  ! vyrr     | fxy
  ! vyrl     | fxz
  ! vylr     | fyx
  ! vyll     | fyz
  ! Byr      | fzx
  ! Byl      | fzy
  ! vzrr     | ftxy
  ! vzrl     | ftxz
  ! vzlr     | ftyx
  ! vzll     | ftyz
  call ppm_find_face_vals_v2_cpp(cctkGH,cctk_lsh,cctk_nghostzones,vyr,temp3, &
       fxy,fxz,m,U_syms2,Symmetry,cell_centering_enabled)
  call ppm_find_face_vals_v2_cpp(cctkGH,cctk_lsh,cctk_nghostzones,vyl,temp4, &
       fyx,fyz,m,U_syms2,Symmetry,cell_centering_enabled)
  call ppm_find_face_vals_v2_cpp(cctkGH,cctk_lsh,cctk_nghostzones,vzr,temp8, &
       ftxy,ftxz,m,U_syms3,Symmetry,cell_centering_enabled)
  call ppm_find_face_vals_v2_cpp(cctkGH,cctk_lsh,cctk_nghostzones,vzl,temp9, &
       ftyx,ftyz,m,U_syms3,Symmetry,cell_centering_enabled)
  call ppm_find_face_vals_v2_cpp(cctkGH,cctk_lsh,cctk_nghostzones,By_stagger,dBy_m, &
       fzx,fzy,m,U_syms7,Symmetry,cell_centering_enabled)

  ! Impose monotonicity last
!!! Impose monotonicity
  !!call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vyr,fxy,fxz,Symmetry)
  !!call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vyl,fyx,fyz,Symmetry)
  !!call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vzr,ftxy,ftxz,Symmetry)
  !!call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vzl,ftyx,ftyz,Symmetry)
  !!call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,By_stagger,fzx,fzy,Symmetry)

  ! Compute ftilde and store it to temp1
  if (em_evolve_enable .ne. 2) then
     call ppm_ftilde_v2(cctk_lsh,cctk_nghostzones,temp1,P,vx,vy,vz,U_syms5, &
          U_syms1,U_syms2,U_syms3, &
          m,Symmetry,cell_centering_enabled)
  end if

  ! Need to do a sync before going to the rest of the reconstruction

end subroutine flux_induction_Aevolution_ppm_step2

subroutine flux_induction_Aevolution_ppm_step2p5(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !OTHER PARAMETERS:
  integer   :: i,j,k,ip1,jp1,kp1

  real*8                :: horizdirn_x,horizdirn_y,horizdirn_z,AH_radius_minimum

  !========================
  ! variable | where stored
  !==========|=============
  ! vyrr     | fxy
  ! vyrl     | fxz
  ! vylr     | fyx
  ! vyll     | fyz
  ! Byr      | fzx
  ! Byl      | fzy
  ! vzrr     | ftxy
  ! vzrl     | ftxz
  ! vzlr     | ftyx
  ! vzll     | ftyz

  ! Compute f and store it to temp2
  if (em_evolve_enable .ne. 2) then 
     ! Compute f and store it to temp2
     if(use_std_PPM_instead_of_efficient_3_gz_version.eq.0) then
        !$omp parallel do
        do k=1,cctk_lsh(3)
           do j=1,cctk_lsh(2)
              do i=1,cctk_lsh(1)
                 temp2(i,j,k) = temp1(i,j,k)
              end do
           end do
        end do
        !$omp end parallel do
     else
        call ppm_f_cpp(cctkGH,cctk_lsh,m,temp1,temp2,P)
     end if

     ! Apply flattening
     call ppm_flatten_1var_cpp(cctkGH,cctk_lsh,m,temp2,vyr,fxy,fxz)
     call ppm_flatten_1var_cpp(cctkGH,cctk_lsh,m,temp2,vyl,fyx,fyz)
     call ppm_flatten_1var_cpp(cctkGH,cctk_lsh,m,temp2,By_stagger,fzx,fzy)
     call ppm_flatten_1var_cpp(cctkGH,cctk_lsh,m,temp2,vzr,ftxy,ftxz)
     call ppm_flatten_1var_cpp(cctkGH,cctk_lsh,m,temp2,vzl,ftyx,ftyz)
  end if

  ! Impose monotonicity
  call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vyr,fxy,fxz,Symmetry)
  call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vyl,fyx,fyz,Symmetry)
  call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vzr,ftxy,ftxz,Symmetry)
  call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vzl,ftyx,ftyz,Symmetry)
  call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,By_stagger,fzx,fzy,Symmetry)

  ! Note: I am skipping the ppm_shift step since we are going to shift the
  ! stencil again below. I am combining the two shifts into one.
  ! By skipping the ppm_shift step, the variables change its meaning.
  ! Here is the new map:
  !=======================================
  ! gf    |   variable  |   staggering    |
  ! ===== | =========== | =============== |
  ! fxy   |    vyrl     | (i,j-1/2,k+1/2) |
  ! fxz   |    vyrr     | (i,j-1/2,k-1/2) |
  ! fyx   |    vyll     | (i,j-1/2,k+1/2) |
  ! fyz   |    vylr     | (i,j-1/2,k-1/2) |
  ! fzx   |    Byl      | (i,j+1/2,k+1/2) |
  ! fzy   |    Byr      | (i,j+1/2,k-1/2) |
  ! ftxy  |    vzrl     | (i,j-1/2,k+1/2) |
  ! ftxz  |    vzrr     | (i,j-1/2,k-1/2) |
  ! ftyx  |    vzll     | (i,j-1/2,k+1/2) |
  ! ftyz  |    vzlr     | (i,j-1/2,k-1/2) |
  ! ftzx  |    Bzl      | (i,j+1/2,k+1/2) |
  ! ftzy  |    Bzr      | (i,j-1/2,k+1/2) |

  ! Shift the stencil and change B^i_rl to tildeB^i_rl
  ! (Note: We want the reconstructed variables to be at 
  !   (i,j+1/2,k+1/2), but all the reconstructed v's 
  !   and B's are located as listed above. 
  !   Hence we need to shift the stencil.)
  !========================
  ! variable | where stored
  !==========|=============
  ! vyrr     | temp1
  ! vyrl     | temp2
  ! vylr     | temp3
  ! vyll     | temp4
  ! vzrr     | temp5
  ! vzrl     | temp6
  ! vzlr     | temp7
  ! vzll     | temp8
  !omp parallel do private (jp1,kp1)
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
	   jp1 = min(j+1,cctk_lsh(2))
	   kp1 = min(k+1,cctk_lsh(3))
           temp2(i,j,k) = fxy(i,jp1,k)
           temp1(i,j,k) = fxz(i,jp1,kp1)
           temp4(i,j,k) = fyx(i,jp1,k)
           temp3(i,j,k) = fyz(i,jp1,kp1)
           temp6(i,j,k) = ftxy(i,jp1,k)
           temp5(i,j,k) = ftxz(i,jp1,kp1)
           temp8(i,j,k) = ftyx(i,jp1,k)
           temp7(i,j,k) = ftyz(i,jp1,kp1)
           Byl(i,j,k) = fzx(i,j,k)*temp11(i,j,k)
           Byr(i,j,k) = fzy(i,j,kp1)*temp11(i,j,k)
           Bzl(i,j,k) = ftzx(i,j,k)*temp11(i,j,k)
           Bzr(i,j,k) = ftzy(i,jp1,k)*temp11(i,j,k)
        end do
     end do
  end do
  !omp end parallel do

  ! Switch to first-order scheme inside the AH
  !!if(num_BHs.gt.0.5D0) then
  !!   !     Get horizon radius in direction where it is likely to be minimized:
  !!   horizdirn_x = 0.D0
  !!   horizdirn_y = 0.D0
  !!   horizdirn_z = 100000.D0
  !!   call get_ah_radius_in_dirn(cctkGH,horizdirn_x,horizdirn_y,horizdirn_z,AH_radius_minimum);
  !!   if(AH_radius_minimum.lt.0.D0) then
  !!      AH_radius_minimum = bh_radius_z(1)
  !!      write(*,*) "WARNING:  Couldn't find horizon, so using last known BH radius & position to limit B^2 inside the horizon.  Radius=",AH_radius_minimum,num_BHs,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1)
  !!   end if

  !!   ! Hard-coded 0.8 safety factor here:
  !!   AH_radius_minimum = 0.8D0 * AH_radius_minimum

  !!   ! Switch to first-order scheme inside AH
  !!   !omp parallel do private (jp1,kp1)
  !!   do k=1,cctk_lsh(3)
  !!      do j=1,cctk_lsh(2)
  !!         do i=1,cctk_lsh(1)
  !!            if(sqrt((x(i,j,k)-bh_posn_x(1))**2 + (y(i,j,k)-bh_posn_y(1))**2 + (z(i,j,k)-bh_posn_z(1))**2) .lt. AH_radius_minimum) then
  !!               jp1 = min(j+1,cctk_lsh(2))
  !!               kp1 = min(k+1,cctk_lsh(3))
  !!               temp1(i,j,k) = vy(i,jp1,kp1)
  !!               temp2(i,j,k) = vy(i,jp1,k)
  !!               temp3(i,j,k) = vy(i,j,kp1)
  !!               temp4(i,j,k) = vy(i,j,k)
  !!               temp5(i,j,k) = vz(i,jp1,kp1)
  !!               temp6(i,j,k) = vz(i,jp1,k)
  !!               temp7(i,j,k) = vz(i,j,kp1)
  !!               temp8(i,j,k) = vz(i,j,k)
  !!               Byl(i,j,k) = By_stagger(i,j,k)*temp11(i,j,k)
  !!               Byr(i,j,k) = By_stagger(i,j,kp1)*temp11(i,j,k)
  !!               Bzl(i,j,k) = Bz_stagger(i,j,k)*temp11(i,j,k)
  !!               Bzr(i,j,k) = Bz_stagger(i,jp1,k)*temp11(i,j,k)
  !!            end if
  !!         end do
  !!      end do
  !!   end do
  !!   !omp end parallel do
  !!end if
  ! ************

  ! Compute the 4 reconstructed F1=(v cross Btilde)^x values at (i,j+1/2,k+1/2).
  ! Also compute Byr - Byl (stored in fzx) and Bzr - Bzl (stored in fzy). 
  ! (Note: the reconstruction order is first y and then z)
  ! Variables are stored as follows::
  !========================
  ! variable  | where stored
  !===========|=============
  ! F1rr      | fxy
  ! F1rl      | fxz
  ! F1lr      | fyx
  ! F1ll      | fyz
  ! Byr - Byl | fzx
  ! Bzr - Bzl | fzy
  !omp parallel do
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
	   fxy(i,j,k) = temp1(i,j,k)*Bzr(i,j,k) - temp5(i,j,k)*Byr(i,j,k)
           fxz(i,j,k) = temp2(i,j,k)*Bzr(i,j,k) - temp6(i,j,k)*Byl(i,j,k)
	   fyx(i,j,k) = temp3(i,j,k)*Bzl(i,j,k) - temp7(i,j,k)*Byr(i,j,k)
	   fyz(i,j,k) = temp4(i,j,k)*Bzl(i,j,k) - temp8(i,j,k)*Byl(i,j,k)
	   fzx(i,j,k) = Byr(i,j,k) - Byl(i,j,k)
	   fzy(i,j,k) = Bzr(i,j,k) - Bzl(i,j,k)
        end do
     end do
  end do
  !omp end parallel do

  !Zach says: I don't think Carpet respects CCTK_SyncGroup()'s.  At best, these are inert.
  !call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::em_fijs')
  !!call CartSymGN(dummy,cctkGH,'mhd_evolve::em_fijs') 

  ! In the previous hydro routine, cmax and cmin were computed 
  ! in the x-direction, at (i-1/2,j,k).  Here we shift the 
  ! stencil so that cmax/cmin at (i+1/2,j,k) is stored at temp10/MONOPOLE,
  ! respectively.  This way we can later compute Az_rhs (and later
  ! Ay_rhs) using the 2D HLL flux (Del Zanna, Bucciantini & Londrillo
  ! A&A 400, 397 (2003), Eq. (44).  Note that cmax/cmin (\alpha^{\pm}
  ! as defined in the paper) is at a slightly DIFFERENT point than that 
  ! described in the Del Zanna et al paper (e.g., (i+1/2,j,k) instead of
  ! (i+1/2,j+1/2,k) for F3.  Yuk Tung discussed this point with Shibata,
  ! who found that the effect is negligible.
  ! Finally, store cmax/cmin in ftxy/ftxz, respectively (to be used in step 3)
  !omp parallel do private (jp1)
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
	   jp1 = min(j+1,cctk_lsh(2))
	   ftxy(i,j,k) = cmax(i,jp1,k)
	   ftxz(i,j,k) = cmin(i,jp1,k)
        end do
     end do
  end do
  !omp end parallel do

end subroutine flux_induction_Aevolution_ppm_step2p5


!=====================================================================================
! In this step, first Ax_rhs is computed.  Next, we reconstruct 
! the flux F2 = vz Bxtilde - vx Bztilde at (i+1/2,j,k+1/2).
! Finally, we compute Ay_rhs
!=====================================================================================
subroutine flux_induction_Aevolution_ppm_step3(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer  ::  dummy,i,j,k,cell_centering_enabled,ip1,jp1,kp1
  integer, dimension(3) :: U_syms1,U_syms2,U_syms3,U_syms4,U_syms5
  integer, dimension(3) :: U_syms6,U_syms7,U_syms8
  integer, parameter :: VECTORX = 2,VECTORY = 3,VECTORZ = 4
  integer, parameter :: BVECTORX = 5, BVECTORY = 6,BVECTORZ = 7, SCALAR = 1

  ! In the previous hydro routine, cmax and cmin were computed 
  ! in the z-direction, at (i,j,k-1/2).  Here we shift the 
  ! stencil so that cmax/cmin is at (i,j,k+1/2).  This way we 
  ! can compute Ax_rhs (and later Ay_rhs) using the 2D HLL flux
  ! (Del Zanna, Bucciantini & Londrillo A&A 400, 397 (2003), Eq. (44).  
  ! Note that cmax/cmin (\alpha^{\pm} as defined in the paper) is at 
  ! a slightly DIFFERENT point than that described in the Del Zanna et 
  ! al paper (e.g., (i+1/2,j,k) instead of (i+1/2,j+1/2,k) for F3).  
  ! Yuk Tung discussed this point with Shibata, who found that the
  ! effect is negligible.
  !omp parallel do private (kp1)
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
           kp1 = min(k+1,cctk_lsh(3))
           temp1(i,j,k) = cmax(i,j,kp1)
           temp2(i,j,k) = cmin(i,j,kp1)
        end do
     end do
  end do
  !omp end parallel do

  !omp parallel do
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
           cmax(i,j,k) = temp1(i,j,k)
           cmin(i,j,k) = temp2(i,j,k)
        end do
     end do
  end do
  !omp end parallel do

  ! Compute Ax_rhs=-Ex=-F1 using the 2D HLL flux.
  ! See Del Zanna, Bucciantini & Londrillo A&A 400, 397 (2003), Eq. (44)

  ! In the previous step, 
  ! the fluxes F1=(v cross Btilde)^x at (i,j+1/2,k+1/2)
  ! have been stored in fxy (F1rr), fxz (F1rl), fyx (F1lr) and fyz (F1ll). 
  ! Also, Byr - Byl is stored in fzx, Bzr - Bzl is stored in fzy; 
  ! cmax, cmin in the y-direction are stored in ftxy and ftxz, respectively.
  ! cmax and cmin in the z-direction have just been computed and stored 
  ! in cmax and cmin, albeit at slightly different points than prescribed
  ! in Del Zanna et al. (see comment above).  In summary:
  !========================
  ! variable  | where stored
  !===========|=============
  ! F1rr      | fxy
  ! F1rl      | fxz
  ! F1lr      | fyx
  ! F1ll      | fyz
  ! Byr - Byl | fzx
  ! Bzr - Bzl | fzy
  call flux_hll2d_cpp(cctkGH,cctk_lsh,ftxy, ftxz, cmax,cmin, &
       fxy, fxz, fyx, fyz, fzx, fzy, & 
       Bxtilde_or_Ax_rhs,use_central_scheme_instead_of_hll)

  !Zach says: I don't think Carpet respects CCTK_SyncGroup()'s.  At best, these are inert.
  !call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')
  !!call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')

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

  ! Reconstruct Bx at (i+1/2,j,k+1/2)
  ! Bx_stagger stores values of Bx on a semi-staggered grid, at points
  ! (i+1/2,j,k).  Here, we compute dU for Bx, at points (i+1/2,j,k-1/2):
  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms6,Bx_stagger,temp1)
  !Apply slope limiter
  call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp1,dBx_m)

  ! Compute Bxr (stored in ftzx) and Bxl (stored in ftzy)
  call ppm_find_face_vals_v2_cpp(cctkGH,cctk_lsh,cctk_nghostzones,Bx_stagger,dBx_m, &
       ftzx,ftzy,m,U_syms6,Symmetry,cell_centering_enabled)
  ! Apply flattening:
  !Note that f has been computed in the hydro step and stored in temp9
  if (em_evolve_enable .ne. 2) then 
     call ppm_flatten_1var_cpp(cctkGH,cctk_lsh,m,temp9,Bx_stagger,ftzx,ftzy)
  end if
  call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,Bx_stagger,ftzx,ftzy,Symmetry)

  ! Compute e^(6 phi) at (i+1/2,j,k+1/2) and stored in temp11
  m = 1
  call compute_face_avg_ppm_cpp(cctkGH,m,cctk_lsh,cctk_nghostzones,cell_centering_enabled,Symmetry,U_syms4,phi_f,temp1)
  !omp parallel do private (ip1,kp1)
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
	   kp1 = min(k+1,cctk_lsh(3))
	   ip1 = min(i+1,cctk_lsh(1))
	   temp11(i,j,k) = exp( 6.d0*temp1(ip1,j,kp1) )
        end do
     end do
  end do
  !omp end parallel do

  ! Reconstruct Bz, vz and vx at (i+1/2,j,k+1/2).
  ! Note that we have already reconstructed vz and vx along the z-direction, 
  ! at (i,j,k-1/2). That result is stored in v{z,x}{r,l}.  Bz_stagger data
  ! are on a semi-staggered grid, at (i,j,k+1/2).  Here, we compute dU for 
  ! vz and vx at (i-1/2,j,k-1/2), and Bz at (i-1/2,j,k+1/2), and store the 
  ! result in temp{1,2,6,7,5}:
  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms3,vzr,temp1)
  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms3,vzl,temp2)
  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms1,vxr,temp6)
  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms1,vxl,temp7)
  call compute_dU_everywhere(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,cell_centering_enabled,U_syms8,Bz_stagger,temp5)

  !Next, we apply the MC slope limiter, storing the result in temp{3,4,8,9},dBz_m
  call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp1,temp3)
  call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp2,temp4)
  call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp6,temp8)
  call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp7,temp9)
  call mc_slope_limiter(m,cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,temp5,dBz_m)

  ! Finally, from the slope-limited dU data, we compute Ul and Ur,
  ! using a host of fij's and ftij's as temporary storage:
  !========================
  ! variable | where stored
  !==========|=============
  ! vzrr     | fxy
  ! vzrl     | fxz
  ! vzlr     | fyx
  ! vzll     | fyz
  ! Bzr      | fzx
  ! Bzl      | fzy
  ! vxrr     | ftxy
  ! vxrl     | ftxz
  ! vxlr     | ftyx
  ! vxll     | ftyz
  call ppm_find_face_vals_v2_cpp(cctkGH,cctk_lsh,cctk_nghostzones,vzr,temp3, &
       fxy,fxz,m,U_syms3,Symmetry,cell_centering_enabled)
  call ppm_find_face_vals_v2_cpp(cctkGH,cctk_lsh,cctk_nghostzones,vzl,temp4, &
       fyx,fyz,m,U_syms3,Symmetry,cell_centering_enabled)
  call ppm_find_face_vals_v2_cpp(cctkGH,cctk_lsh,cctk_nghostzones,vxr,temp8, &
       ftxy,ftxz,m,U_syms1,Symmetry,cell_centering_enabled)
  call ppm_find_face_vals_v2_cpp(cctkGH,cctk_lsh,cctk_nghostzones,vxl,temp9, &
       ftyx,ftyz,m,U_syms1,Symmetry,cell_centering_enabled)
  call ppm_find_face_vals_v2_cpp(cctkGH,cctk_lsh,cctk_nghostzones,Bz_stagger,dBz_m, &
       fzx,fzy,m,U_syms8,Symmetry,cell_centering_enabled)

  ! Impose monotonicity last
!!! Impose monotonicity
  !!call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vzr,fxy,fxz,Symmetry)
  !!call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vzl,fyx,fyz,Symmetry)
  !!call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vxr,ftxy,ftxz,Symmetry)
  !!call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vxl,ftyx,ftyz,Symmetry)
  !!call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,Bz_stagger,fzx,fzy,Symmetry)

  ! Compute ftilde and store it to temp1
  if (em_evolve_enable .ne. 2) then 
     call ppm_ftilde_v2(cctk_lsh,cctk_nghostzones,temp1,P,vx,vy,vz,U_syms5, &
          U_syms1,U_syms2,U_syms3, &
          m,Symmetry,cell_centering_enabled)
  end if

  ! Need to do a sync before going to the rest of the reconstruction

end subroutine flux_induction_Aevolution_ppm_step3

subroutine flux_induction_Aevolution_ppm_step3p5(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !OTHER PARAMETERS:
  integer   :: i,j,k,ip1,jp1,kp1

  real*8                :: horizdirn_x,horizdirn_y,horizdirn_z,AH_radius_minimum

  !========================
  ! variable | where stored
  !==========|=============
  ! vzrr     | fxy
  ! vzrl     | fxz
  ! vzlr     | fyx
  ! vzll     | fyz
  ! Bzr      | fzx
  ! Bzl      | fzy
  ! vxrr     | ftxy
  ! vxrl     | ftxz
  ! vxlr     | ftyx
  ! vxll     | ftyz

  ! Compute f and store it to temp2
  if (em_evolve_enable .ne. 2) then 
     ! Compute f and store it to temp2
     if(use_std_PPM_instead_of_efficient_3_gz_version.eq.0) then
        !$omp parallel do
        do k=1,cctk_lsh(3)
           do j=1,cctk_lsh(2)
              do i=1,cctk_lsh(1)
                 temp2(i,j,k) = temp1(i,j,k)
              end do
           end do
        end do
        !$omp end parallel do
     else
        call ppm_f_cpp(cctkGH,cctk_lsh,m,temp1,temp2,P)
     end if

     ! Apply flattening
     call ppm_flatten_1var_cpp(cctkGH,cctk_lsh,m,temp2,vzr,fxy,fxz)
     call ppm_flatten_1var_cpp(cctkGH,cctk_lsh,m,temp2,vzl,fyx,fyz)
     call ppm_flatten_1var_cpp(cctkGH,cctk_lsh,m,temp2,Bz_stagger,fzx,fzy)
     call ppm_flatten_1var_cpp(cctkGH,cctk_lsh,m,temp2,vxr,ftxy,ftxz)
     call ppm_flatten_1var_cpp(cctkGH,cctk_lsh,m,temp2,vxl,ftyx,ftyz)
  end if

  ! Impose monotonicity
  call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vzr,fxy,fxz,Symmetry)
  call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vzl,fyx,fyz,Symmetry)
  call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vxr,ftxy,ftxz,Symmetry)
  call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,vxl,ftyx,ftyz,Symmetry)
  call ppm_monotonize_v2_cpp(cctkGH,cctk_lsh,Bz_stagger,fzx,fzy,Symmetry)

  ! Note: I am skipping the ppm_shift step since we are going to shift the
  ! stencil again below. I am combining the two shifts into one.
  ! By skipping the ppm_shift step, the variables change its meaning.
  ! Here is the new map:
  !=======================================
  ! gf    |   variable  |   staggering    |
  ! ===== | =========== | =============== |
  ! fxy   |    vzrl     | (i+1/2,j,k-1/2) |
  ! fxz   |    vzrr     | (i-1/2,j,k-1/2) |
  ! fyx   |    vzll     | (i+1/2,j,k-1/2) |
  ! fyz   |    vzlr     | (i-1/2,j,k-1/2) |
  ! fzx   |    Bzl      | (i+1/2,j,k+1/2) |
  ! fzy   |    Bzr      | (i-1/2,j,k+1/2) |
  ! ftxy  |    vxrl     | (i+1/2,j,k-1/2) |
  ! ftxz  |    vxrr     | (i-1/2,j,k-1/2) |
  ! ftyx  |    vxll     | (i+1/2,j,k-1/2) |
  ! ftyz  |    vxlr     | (i-1/2,j,k-1/2) |
  ! ftzx  |    Bxl      | (i+1/2,j,k+1/2) |
  ! ftzy  |    Bxr      | (i+1/2,j,k-1/2) |


  ! Shift the stencil and change B^i_rl to tildeB^i_rl
  ! (Note: We want the reconstructed variables to be at
  !   (i+1/2,j,k+1/2), but the reconstructed v's
  !   and B's are located at positions listed above.
  !   Hence we need to shift the stencil)
  !========================
  ! variable | where stored
  !==========|=============
  ! vzrr     | temp1
  ! vzrl     | temp2
  ! vzlr     | temp3
  ! vzll     | temp4
  ! vxrr     | temp5
  ! vxrl     | temp6
  ! vxlr     | temp7
  ! vxll     | temp8
  !omp parallel do private (ip1,kp1)
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
	   ip1 = min(i+1,cctk_lsh(1))
	   kp1 = min(k+1,cctk_lsh(3))
           temp2(i,j,k) = fxy(i,j,kp1)
           temp1(i,j,k) = fxz(ip1,j,kp1)
           temp4(i,j,k) = fyx(i,j,kp1)
           temp3(i,j,k) = fyz(ip1,j,kp1)
           temp6(i,j,k) = ftxy(i,j,kp1)
           temp5(i,j,k) = ftxz(ip1,j,kp1)
           temp8(i,j,k) = ftyx(i,j,kp1)
           temp7(i,j,k) = ftyz(ip1,j,kp1)
           Bzl(i,j,k) = fzx(i,j,k)*temp11(i,j,k)
           Bzr(i,j,k) = fzy(ip1,j,k)*temp11(i,j,k)
           Bxl(i,j,k) = ftzx(i,j,k)*temp11(i,j,k)
           Bxr(i,j,k) = ftzy(i,j,kp1)*temp11(i,j,k)
        end do
     end do
  end do
  !omp end parallel do

  ! Switch to first-order scheme inside the AH
  !!if(num_BHs.gt.0.5D0) then
  !!   !     Get horizon radius in direction where it is likely to be minimized:
  !!   horizdirn_x = 0.D0
  !!   horizdirn_y = 0.D0
  !!   horizdirn_z = 100000.D0
  !!   call get_ah_radius_in_dirn(cctkGH,horizdirn_x,horizdirn_y,horizdirn_z,AH_radius_minimum);
  !!   if(AH_radius_minimum.lt.0.D0) then
  !!      AH_radius_minimum = bh_radius_z(1)
  !!      write(*,*) "WARNING:  Couldn't find horizon, so using last known BH radius & position to limit B^2 inside the horizon.  Radius=",AH_radius_minimum,num_BHs,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1)
  !!   end if

  !!   !        write(*,*) "Radius=",AH_radius_minimum,num_BHs,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1)


  !!   ! Hard-coded 0.8 safety factor here:
  !!   AH_radius_minimum = 0.8D0 * AH_radius_minimum

  !!   ! Switch to first-order scheme inside AH
  !!   !omp parallel do private (ip1,kp1)
  !!   do k=1,cctk_lsh(3)
  !!      do j=1,cctk_lsh(2)
  !!         do i=1,cctk_lsh(1)
  !!            if(sqrt((x(i,j,k)-bh_posn_x(1))**2 + (y(i,j,k)-bh_posn_y(1))**2 + (z(i,j,k)-bh_posn_z(1))**2) .lt. AH_radius_minimum) then
  !!               ip1 = min(i+1,cctk_lsh(1))
  !!               kp1 = min(k+1,cctk_lsh(3))
  !!               temp1(i,j,k) = vz(ip1,j,kp1)
  !!               temp2(i,j,k) = vz(i,j,kp1)
  !!               temp3(i,j,k) = vz(ip1,j,k)
  !!               temp4(i,j,k) = vz(i,j,k)
  !!               temp5(i,j,k) = vx(ip1,j,kp1)
  !!               temp6(i,j,k) = vx(i,j,kp1)
  !!               temp7(i,j,k) = vx(ip1,j,k)
  !!               temp8(i,j,k) = vx(i,j,k)
  !!               Bzl(i,j,k) = Bz_stagger(i,j,k)*temp11(i,j,k)
  !!               Bzr(i,j,k) = Bz_stagger(ip1,j,k)*temp11(i,j,k)
  !!               Bxl(i,j,k) = Bx_stagger(i,j,k)*temp11(i,j,k)
  !!               Bxr(i,j,k) = Bx_stagger(i,j,kp1)*temp11(i,j,k)
  !!            end if
  !!         end do
  !!      end do
  !!   end do
  !!   !omp end parallel do
  !!end if

  ! Compute the 4 reconstructed F2=(v cross Btilde)^y values at (i+1/2,j,k+1/2).
  ! Also compute Bzr - Bzl (stored in fzx) and Bxr - Bxl (stored in fzy). 
  ! (Note: the reconstruction order is first z and then x)
  ! Variables are stored as follows::
  !========================
  ! variable  | where stored
  !===========|=============
  ! F2rr      | fxy
  ! F2rl      | fxz
  ! F2lr      | fyx
  ! F2ll      | fyz
  ! Bzr - Bzl | fzx
  ! Bxr - Bxl | fzy
  !omp parallel do
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
	   fxy(i,j,k) = temp1(i,j,k)*Bxr(i,j,k) - temp5(i,j,k)*Bzr(i,j,k)
           fxz(i,j,k) = temp2(i,j,k)*Bxr(i,j,k) - temp6(i,j,k)*Bzl(i,j,k)
	   fyx(i,j,k) = temp3(i,j,k)*Bxl(i,j,k) - temp7(i,j,k)*Bzr(i,j,k)
	   fyz(i,j,k) = temp4(i,j,k)*Bxl(i,j,k) - temp8(i,j,k)*Bzl(i,j,k)
	   fzx(i,j,k) = Bzr(i,j,k) - Bzl(i,j,k)
	   fzy(i,j,k) = Bxr(i,j,k) - Bxl(i,j,k)
        end do
     end do
  end do
  !omp end parallel do

  !Zach says: I don't think Carpet respects CCTK_SyncGroup()'s.  At best, these are inert.
  !call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::em_fijs')
  !call CartSymGN(dummy,cctkGH,'mhd_evolve::em_fijs') 

  ! Compute Ay_rhs=-Ey=-F2 using the 2D HLL flux.
  ! See Del Zanna, Bucciantini & Londrillo A&A 400, 397 (2003), Eq. (44)   
  !
  ! The fluxes F2=(v cross Btilde)^y at (i+1/2,j,k+1/2)
  ! have been stored in fxy (F2rr), fxz (F2rl), fyx (F2lr) and fyz (F2ll). 
  ! Also, Bzr - Bzl is stored in fzx, Bxr - Bxl is stored in fzy; 
  ! cmax, cmin in the x-direction are stored in temp11 and MONOPOLE, respectively.
  ! cmax and cmin in the z-direction have just been computed and stored 
  ! in cmax and cmin, albeit at slightly different points than prescribed
  ! in Del Zanna et al. (see comment above).  In summary:
  !========================
  ! variable  | where stored
  !===========|=============
  ! F2rr      | fxy
  ! F2rl      | fxz
  ! F2lr      | fyx
  ! F2ll      | fyz
  ! Bzr - Bzl | fzx
  ! Bxr - Bxl | fzy
  call flux_hll2d_cpp(cctkGH,cctk_lsh,cmax, cmin, temp10,MONOPOLE, &
       fxy, fxz, fyx, fyz, fzx, fzy, &
       Bytilde_or_Ay_rhs,use_central_scheme_instead_of_hll)

end subroutine flux_induction_Aevolution_ppm_step3p5
