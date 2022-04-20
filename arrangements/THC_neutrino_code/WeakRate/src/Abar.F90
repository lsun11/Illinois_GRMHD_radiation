#include <weaktableindex.h>
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

CCTK_INT FUNCTION NucleiAbarImpl(rho, temp, ye, abar)
    use table3d_mod
    use lk_interpolations
    use weak_constants

    implicit none

    CCTK_REAL, INTENT(IN)  :: rho, temp, ye !density, temperature, Ye
    CCTK_REAL, INTENT(OUT) :: abar

    CCTK_REAL :: rho_cgs, lrho0, temp0, ltemp0, ye0
    INTEGER :: enforceTableBounds
    INTEGER :: boundsErr

    NucleiAbarImpl = 0

    !Convert the density in cgs
    rho_cgs  = rho * cactus2cgsRho
    temp0 = temp
    ye0   = ye

    boundsErr = enforceTableBounds(rho_cgs,temp0,ye0)

    !TODO add a return error message when out of bounds
    if (boundsErr.eq.-1) then
        NucleiAbarImpl = -1
    end if

    !Converto to cgs units and take the log10
    lrho0  = log10(rho_cgs)
    ltemp0 = log10(temp0)

    abar = linearInterpolation3d(lrho0,ltemp0,ye0,ABAR)

    return
END FUNCTION NucleiAbarImpl

