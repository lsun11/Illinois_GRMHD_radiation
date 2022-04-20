#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine BHNS_setup_emfield_part1_Pmax(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer                                  :: handle,vindex,ierr
  real*8                                   :: reduction_value

  bhnsinsertBNOW=0

if(CCTK_ITERATION .eq. ITERATION_TO_INSERT_MAGNETIC_FIELDS .or. (ITERATION_TO_INSERT_MAGNETIC_FIELDS.lt.0 .and. CCTK_ITERATION .eq. 0) ) then


   print *, "CCTK_ITERATION=",CCTK_ITERATION, "ITERATION_TO_INSERT_MAGNETIC_FIELDS=",ITERATION_TO_INSERT_MAGNETIC_FIELDS
     bhnsinsertBNOW=1

     ! Find P_max
     call CCTK_VarIndex(vindex,"mhd_evolve::P")
     call CCTK_ReductionHandle(handle,"maximum")
     if (handle .gt. 0) then
        call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,reduction_value,1,vindex)
     else
        call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
     end if
     bhns_P_max = reduction_value

     write(*,*) "INSIDE BHNS_setup_emfield_part1_Pmax.  PMAX = ",bhns_P_max

     if(unequalmass.eq.1) print *,"bhns_P_max2=",bhns_P_max2
     ! Find rhob_max
     call CCTK_VarIndex(vindex,"mhd_evolve::rho_b")
     call CCTK_ReductionHandle(handle,"maximum")
     if (handle .gt. 0) then
        call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,reduction_value,1,vindex)
     else
        call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
     end if
     bhns_rhob_max = reduction_value
     
     write(*,*) " driver_setup_emfields_part1_Pmax.F90 (x_NS,y_NS)=",CoMx_VolInt,CoM_VolInt_denominator,CoMy_VolInt,CoM_VolInt_denominator    
     write(*,*) "INSIDE BHNS_setup_emfield_part1_Pmax.  PMAX = ",bhns_P_max, "rhobMAX = ", bhns_rhob_max
  end if

  write(*,*) "bhnsinsertBNOW",bhnsinsertBNOW
  write(*,*) "EXIT BHNS_setup_emfield_part1_Pmax",bhnsinsertBNOW
end subroutine BHNS_setup_emfield_part1_Pmax
