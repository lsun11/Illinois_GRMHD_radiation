#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bhns_initialdata_global2(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext

  real*8 			           :: dT,dX,dY,dZ
  real*8 			           :: n1,n2,n3,mf

  integer :: vindex,handle,ierr
  real*8  :: rho_max,tau_max
  real*8  :: ONE,ZERO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(ONE = 1.D0, ZERO = 0.D0)

  INTEGER :: i,j,k
  character :: filename*30,c2*2,c1
  CCTK_REAL :: reduction_value

  real*8, dimension(1,3)                   :: pointcoords
  real*8 :: tauinterp
  pointcoords(1,1) = 2.305D0
  pointcoords(1,2) =-4.25D-2
  pointcoords(1,3) = 4.25D-2
  call CCTK_VarIndex(vindex,"mhd_evolve::rho_star")
  call interp_driver_carp(cctkGH,1,pointcoords,vindex,tauinterp)
  write(*,*) "initdataglobal rs tau:",tauinterp

  call CCTK_VarIndex(vindex,"mhd_evolve::tau")
  call interp_driver_carp(cctkGH,1,pointcoords,vindex,tauinterp)
  write(*,*) "initdataglobal tau:",tauinterp


  ! Set tau_atm for primitives solver.
  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(vindex,"mhd_evolve::tau")
  if (handle .gt. 0) then
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,reduction_value,1,vindex)
     if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH).eq.0) then
        print *,"Maximum value of tau is ",reduction_value
     end if
  else
     call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
  end if
  tau_max = reduction_value
  tau_atm = tau_fact*tau_max

  write(*,*) "TAU_ATM: ",tau_atm

end subroutine bhns_initialdata_global2
