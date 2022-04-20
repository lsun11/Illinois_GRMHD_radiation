!----------------------------------------
! Set up rest mass integrand for 1D shock
!----------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine mhd_shock_rest_mass_integrand1d(CCTK_ARGUMENTS)
  implicit none

  !~~~~~> Input parameters:
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !~~~~~> Other variables:
  real*8                             :: dX, dY, dZ, u0_analytic, fudgefactor
  integer                            :: i,j,k
  integer                            :: AXISYM
  parameter(AXISYM = 4)

  !-----------------------------------------------------------------------------
  ! Compute integrand
  !-----------------------------------------------------------------------------
  shock_rest_mass_integrand = 0.D0
  k = cctk_nghostzones(3)+1
  j = cctk_nghostzones(2)+1

  ! We must multiply by cctk_levfac's here because the Carpet weight function assumes 3D integrals.
  !    You see, when a CCTK_Reduce() sum is called, it actually sums the integrand*weight function.
  !    The weight function is set to 1/(cctk_levfac(1)*cctk_levfac(2)*cctk_levfac(3)) in non-
  !    ghostzone regions of the grid, where e.g., cctk_levfac(1)=1 on the coarsest level and
  !    it increases from there by factors of 2 as the grids get finer.
  !    Thus for a 1D integral along the x direction, we need to multiply by cctk_levfac(2)*cctk_levfac(3).

  ! In addition, we must multiply by 4.D0 on all levels because we set up the grid as follows:
  !
  !CoordBase::boundary_size_z_lower = 3
  !CoordBase::boundary_shiftout_z_lower = 1
  !CoordBase::boundary_size_y_lower = 3
  !CoordBase::boundary_shiftout_y_lower = 1
  !
  !    Any time there is a boundary_shiftout_z_lower, Carpet assumes that this a symmetry boundary, and will 
  !    multiply the weight function by a factor of 0.5 for each "symmetry".  Here, Carpet assumes that there
  !    are 2 symmetries here: one in y and one in z.  Thus the weight function on the coarsest level will be 
  !    set to 0.25, and to 0.25*0.25 on the finer level (due to the additional cctk_levfac factors of 2.  Note
  !    that these symmetry factors work well in the case of, e.g., equatorial symmetry, where the z=0 plane is
  !    included on the grid, since we multiply the final result by symmfactor=2, and we don't want to double-count
  !    the z=0 points.
  fudgefactor=4.D0*cctk_levfac(2)*cctk_levfac(3)

  do i = 1,cctk_lsh(1)
     shock_rest_mass_integrand(i,j,k) = rho_star(i,j,k)*fudgefactor
  end do ! i-loop

end subroutine mhd_shock_rest_mass_integrand1d
