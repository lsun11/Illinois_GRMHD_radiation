!------------------------------------------------------------
! Master update boundary condition (bc) driver for the lapse
!------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine lapse_update_bc(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer                                       :: dummy
  integer, dimension(3)                         :: ext
  real*8                                        :: dX,dY,dZ
!
  dX=CCTK_DELTA_SPACE(1)
  dY=CCTK_DELTA_SPACE(2)
  dZ=CCTK_DELTA_SPACE(3)

  ext = cctk_lsh

  call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')

  if(cctk_iteration.gt.0) then
!     call global_edges(dX,dY,dZ, &
!          Xglobal(1),Yglobal(1),Zglobal(1), &
!          Xlocal1d(1),Ylocal1d(1),Zlocal1d(1), &
!          Xglobal(cctk_gsh(1)),Yglobal(cctk_gsh(2)),Zglobal(cctk_gsh(3)), &
!          Xlocal1d(cctk_lsh(1)),Ylocal1d(cctk_lsh(2)),Zlocal1d(cctk_lsh(3)), &
!          have_bdry_min,have_bdry_max,Symmetry)
     
     if (CCTK_Equals(slicing_type,"frozen").eq.1 .or. cowling_enable==1) then
        ! do nothing; already updated
     else if (CCTK_Equals(slicing_type,"geodesic").eq.1) then
        ! do nothing; already updated 
     else if (CCTK_Equals(slicing_type,"harmonic").eq.1) then
        ! do nothing; already updated 
     else if (CCTK_Equals(slicing_type,"hyperbolic").eq.1) then
        call hyperbolic_lapse_bc_driver(CCTK_PASS_FTOF)
     else if (CCTK_Equals(slicing_type,"opl").eq.1) then
        call opl_lapse_bc_driver(CCTK_PASS_FTOF)
     else if (CCTK_Equals(slicing_type,"opl_loglapse").eq.1) then
        call opl_lapse_bc_driver(CCTK_PASS_FTOF)
     end if

!     call CCTK_SyncGroup(dummy,cctkGH,'lapse::lapse_vars')
  end if

  call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')

end subroutine lapse_update_bc
