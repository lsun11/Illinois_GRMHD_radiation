!--------------------------------------------------------------
! Driver for primitives and boundary conditions routines, v2.0
!--------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
subroutine mhd_emfields_bcs(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext,fake_ext,thisproc_have_global_bdry_min,thisproc_have_global_bdry_max
  real*8                :: dX,dY,dZ,levelnumber,b2bt
  real*8                :: rho_fail_max_step,M_fail_step,randomnumber
  integer               :: proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax
  integer               :: glob_imax,glob_jmax,glob_kmax
  integer               :: i,j,k,vindex,NFont,NFont_L
  integer               :: handle,dummy,ierr
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

     !     if(levelnumber==1 .or. primitives_solver.ne.1) then
     ! Should check if there are NaNs first...

     if(em_evolve_enable==1) then
        if(levelnumber==1.D0) then
           ! First apply BC's to Bitilde's
           ! This is a quick and dirty fix: set thisproc_have_global_bdry_min/max so that
           !    ALL boundaries are updated, EXCEPT in symmetry ghostzones!
           call global_edges(dX,dY,dZ, &
                X(1,1,1),Y(1,1,1),Z(1,1,1), &
                X(1,1,1),Y(1,1,1),Z(1,1,1), &
                X(cctk_lsh(1),1,1),Y(1,cctk_lsh(2),1),Z(1,1,cctk_lsh(3)), &
                X(cctk_lsh(1),1,1),Y(1,cctk_lsh(2),1),Z(1,1,cctk_lsh(3)), &
                thisproc_have_global_bdry_min,thisproc_have_global_bdry_max,Symmetry)
           do i=1,cctk_nghostzones(1)
              fake_ext = cctk_lsh - cctk_nghostzones + i
              if(Symmetry==AXISYM) fake_ext(1) = fake_ext(1) + 1
              call emfields_bc_newv2(ext,fake_ext,X,Y,Z, &
                   Bxtilde_p,Bytilde_p,Bztilde_p,Bxtilde,Bytilde,Bztilde,Symmetry, &
                   thisproc_have_global_bdry_min,thisproc_have_global_bdry_max,EM_BC)
           end do
        end if

     end if
  end if
end subroutine mhd_emfields_bcs
