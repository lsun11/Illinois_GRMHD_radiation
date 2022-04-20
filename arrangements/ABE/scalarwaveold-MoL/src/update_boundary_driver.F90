!-----------------------------------------------------------------------------
! Update outer boundaries
!-----------------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine scalarwaveMoL_update_boundary(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,fake_ext
  real*8                                   :: dT,dx,dy,dz
  integer                                  :: i,AXISYM,temp

  AXISYM = 4

  ext = cctk_lsh

  dT = CCTK_DELTA_TIME
  dx = CCTK_DELTA_SPACE(1)
  dy = CCTK_DELTA_SPACE(2)
  dz = CCTK_DELTA_SPACE(3)

  call CartSymVN(temp,cctkGH,'scalarwaveMoL::phi')
  call CartSymVN(temp,cctkGH,'scalarwaveMoL::phidot')

  write(*,*) "UPDATING BOUNDARIES"
  do i=1,cctk_nghostzones(2)
     fake_ext = cctk_lsh - cctk_nghostzones + i
     if(scalarwave_Symmetry==AXISYM) fake_ext(1) = fake_ext(1) + 1
     call scalarwaveMoL_update_boundary_lowlevel(ext, fake_ext,X, Y, Z, &
          dX,dY,dZ,dT, &
          phi_p,phi,phidot_p,phidot,scalarwave_Symmetry)
     !write(*,*) "hi1.",i,phi_p(cctk_lsh(1)-cctk_nghostzones(1)+i,11,11),phi(cctk_lsh(1)-cctk_nghostzones(1)+i,11,11)
  end do

!  timestep_iteration = timestep_iteration + 1

end subroutine scalarwaveMoL_update_boundary
