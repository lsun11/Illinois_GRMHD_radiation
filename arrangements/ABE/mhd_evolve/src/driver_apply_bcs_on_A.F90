#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
subroutine driver_apply_bcs_on_A(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext,fake_ext,thisproc_have_global_bdry_min,thisproc_have_global_bdry_max
  real*8                :: dT,dX,dY,dZ,levelnumber,b2bt
  real*8                :: rho_fail_max_step,M_fail_step,randomnumber
  integer               :: proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax
  integer               :: glob_imax,glob_jmax,glob_kmax
  integer               :: i,j,k,vindex,NFont,NFont_L
  integer               :: handle,dummy,ierr
  integer               :: repairs_needed,ignore_ghostzones
  CCTK_REAL             :: reduction_value
  integer               :: AXISYM
  parameter(AXISYM = 4)

  ext = cctk_lsh

  levelnumber = cctk_levfac(1)
  levelnumber = log(levelnumber)/log(2.D0)+1.D0

  glob_imax = ext(1)
  glob_jmax = ext(2)
  glob_kmax = ext(3)
  proc_imin = -100
  proc_jmin = -100
  proc_kmin = -100
  proc_imax = -1
  proc_jmax = -1
  proc_kmax = -1

  if(cctk_iteration.gt.0) then
     dT = CCTK_DELTA_TIME
     dX = CCTK_DELTA_SPACE(1)
     dY = CCTK_DELTA_SPACE(2)
     dZ = CCTK_DELTA_SPACE(3)

     thisproc_have_global_bdry_min = have_global_bdry_min(int(levelnumber),:)
     thisproc_have_global_bdry_max = have_global_bdry_max(int(levelnumber),:)

     if(thisproc_have_global_bdry_min(1)==1) proc_imin = 0
     if(thisproc_have_global_bdry_min(2)==1) proc_jmin = 0
     if(thisproc_have_global_bdry_min(3)==1) proc_kmin = 0
     if(thisproc_have_global_bdry_max(1)==1) proc_imax = glob_imax
     if(thisproc_have_global_bdry_max(2)==1) proc_jmax = glob_jmax
     if(thisproc_have_global_bdry_max(3)==1) proc_kmax = glob_kmax


     if(levelnumber==1.D0) then
      ! This is a quick and dirty fix: set thisproc_have_global_bdry_min/max so that
      !    ALL boundaries are updated, EXCEPT in symmetry ghostzones!
        call global_edges(dX,dY,dZ, &
             X(1,1,1),Y(1,1,1),Z(1,1,1), &
             X(1,1,1),Y(1,1,1),Z(1,1,1), &
             X(cctk_lsh(1),1,1),Y(1,cctk_lsh(2),1),Z(1,1,cctk_lsh(3)), &
             X(cctk_lsh(1),1,1),Y(1,cctk_lsh(2),1),Z(1,1,cctk_lsh(3)), &
             thisproc_have_global_bdry_min,thisproc_have_global_bdry_max,Symmetry)
      end if
        
      
      do k=1,cctk_lsh(3)
         do j=1,cctk_lsh(2)
            do i=1,cctk_lsh(1)
               if (isnan(Ax(i,j,k))) then
                  write(*,*) "1. Inside driver_apply_bcs_on_A.F90, Ax is nan !!!!!", X(i,j,k),Y(i,j,k),Z(i,j,k)
               end if
            end do
         end do
      end do
      


      do i=1,cctk_nghostzones(1)
         fake_ext = cctk_lsh - cctk_nghostzones + i
         if(Symmetry==AXISYM) fake_ext(1) = fake_ext(1) + 1
         call emfields_bc_Ai_psi6phi(ext,fake_ext,X,Y,Z, &
              psi6phi_p,Ax_p,Ay_p,Az_p,psi6phi,Ax,Ay,Az,Symmetry, &
              thisproc_have_global_bdry_min,thisproc_have_global_bdry_max,EM_BC)
      end do


      do k=1,cctk_lsh(3)
         do j=1,cctk_lsh(2)
            do i=1,cctk_lsh(1)
               if (isnan(Ax(i,j,k))) then
                  write(*,*) "2. Inside driver_apply_bcs_on_A.F90, Ax is nan !!!!!", X(i,j,k),Y(i,j,k),Z(i,j,k)
               end if
            end do
         end do
      end do


      
  end if

end subroutine driver_apply_bcs_on_A
