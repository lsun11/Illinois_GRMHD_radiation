!--------------------------------------------------------
! Okay, we've read in the initial data from files.
! Now we set up all other required variables, including:
!  emfields, BSSN variables, primitives, etc.
!--------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine magnetar_compute_P_max(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer               :: ierr,index,handle
  integer               :: Nfont, Nfont_l
  CCTK_REAL             :: reduction_value

  integer               :: i,j,k

  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  integer :: ONE,ZERO
  parameter(ONE = 1.D0, ZERO = 0.D0)

  write(*,*) "part3: compute p_max..."

  call CCTK_VarIndex(index,"mhd_evolve::P")
  call CCTK_ReductionHandle(handle,"maximum")
  if (handle .gt. 0) then
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,reduction_value,1,index)
     if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH).eq.0) then
        print *,"Maximum value of P is ",reduction_value
     end if
  else
     call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
  end if
  magnetar_P_max = reduction_value 
  pfloor = P_fact*magnetar_P_max

end subroutine magnetar_compute_P_max
