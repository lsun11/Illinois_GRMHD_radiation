!------------------------------------------
! Timestepping driver for hyperbolic lapse
!------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine hyperbolic_lapse_timestepping(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS   

  interface
     subroutine hyperbolic_lapse_rhs(cctkGH,cctk_nghostzones,ext, &
          hyper_a1,hyper_a2,hyper_a3, &
          r,PhysicalRadius,RadiusDerivative, &
          phi,lapm1,lapm1_rhs, &
          lapset,lapset_rhs, &
          trK_rhs,trK,trK_init)
       implicit none
       CCTK_POINTER                            :: cctkGH
       integer, dimension(3)                   :: ext,cctk_nghostzones
       real*8                                  :: hyper_a1,hyper_a2,hyper_a3
       real*8, dimension(ext(1),ext(2),ext(3)) :: r, PhysicalRadius, RadiusDerivative
       real*8, dimension(ext(1),ext(2),ext(3)) :: phi,lapm1,lapm1_rhs
       real*8, dimension(ext(1),ext(2),ext(3)) :: lapset,lapset_rhs
       real*8, dimension(ext(1),ext(2),ext(3)) :: trK_rhs,trK,trK_init
     end subroutine hyperbolic_lapse_rhs
 end interface
  integer                                   :: n1,n2,n3,m,i,dummy
  integer, dimension(3)                     :: ext,fake_ext
  real*8                                    :: dT,dX
  real*8                                    :: HALF,ld_eps,ld_c
  integer                    :: PI_SYMM, AXISYM
  parameter(PI_SYMM = 3, AXISYM = 4)
  parameter ( HALF = 0.5D0 )
!
  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dT = CCTK_DELTA_TIME

!  call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')
  if(excision_enable == 1) then
     call scalar_excision_bc(ext,X,Y,Z,lapset, &
          Symmetry,excision_zone_gf)
  end if

!  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')
!  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_aux_restrict')

  call hyperbolic_lapse_rhs(cctkGH,cctk_nghostzones,cctk_lsh, &
       hyper_a1,hyper_a2,hyper_a3, &
       r,PhysicalRadius,RadiusDerivative, &
       phi,lapm1,lapm1_rhs, &
       lapset,lapset_rhs, &
       trK_rhs,trK,trK_init)

end subroutine hyperbolic_lapse_timestepping
