#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


subroutine stationary_puncture_id_nonasymptoticisotropic(CCTK_ARGUMENTS)

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext
  real*8 			           :: dT,dX,dY,dZ

  real*8                                   :: psib,lapm1b,rb
  integer :: index,dummy
  integer :: ierr
  real*8  :: ONE,ZERO
  parameter(ONE = 1.D0, ZERO = 0.D0)

  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  phi = log(1.D0 + 1.D0/(2.D0*PhysicalRadius))
  psi = exp(phi)
  lapm1 = (1.D0 + 1.D0/(2.D0*PhysicalRadius))**(-2) - 1.D0
  shiftx = 0.D0
  shifty = 0.D0
  shiftz = 0.D0
  !======================================
  ! Set everything else to flat data
  !======================================
  gxx = ONE
  gxy = ZERO
  gxz = ZERO
  gyy = ONE
  gyz = ZERO
  gzz = ONE
  trK = ZERO
  gupxx = ONE
  gupxy = ZERO
  gupxz = ZERO
  gupyy = ONE
  gupyz = ZERO
  gupzz = ONE

  Gammax = ZERO
  Gammay = ZERO
  Gammaz = ZERO

  !======================================
  ! Set everything else to Zero!
  !======================================
  Gammaxxx = ZERO
  Gammaxxy = ZERO
  Gammaxxz = ZERO
  Gammaxyy = ZERO
  Gammaxyz = ZERO
  Gammaxzz = ZERO
  Gammayxx = ZERO
  Gammayxy = ZERO
  Gammayxz = ZERO
  Gammayyy = ZERO
  Gammayyz = ZERO
  Gammayzz = ZERO
  Gammazxx = ZERO
  Gammazxy = ZERO
  Gammazxz = ZERO
  Gammazyy = ZERO
  Gammazyz = ZERO
  Gammazzz = ZERO
  rho = ZERO

  Axx = ZERO
  Axy = ZERO
  Axz = ZERO
  Ayy = ZERO
  Ayz = ZERO
  Azz = ZERO
  trK = ZERO    

  !We don't need sync's or sym_gz fillers for this case!

end subroutine stationary_puncture_id_nonasymptoticisotropic
