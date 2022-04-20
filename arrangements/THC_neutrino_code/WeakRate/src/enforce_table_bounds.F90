#include "cctk.h"
#include "cctk_Parameters.h"

INTEGER FUNCTION enforceTableBounds(rho, temp, ye)
    use table3d_mod
    implicit none
    real*8, intent(inout) :: rho, temp, ye
    DECLARE_CCTK_PARAMETERS

    enforceTableBounds = 0

    if(rho.lt.eos_rhomin) then
       rho = eos_rhomin
       if(use_rho_min_ext.eq.0) then
          enforceTableBounds = -1
       endif
    endif
    if(rho.gt.eos_rhomax) then
       rho = eos_rhomax
       if(use_rho_max_ext.eq.0) then
          enforceTableBounds = -1
       endif
    endif

    if(temp.lt.eos_tempmin) then
       temp = eos_tempmin
       if(use_temp_min_ext.eq.0) then
          enforceTableBounds = -1
       endif
    endif
    if(temp.gt.eos_tempmax) then
       temp = eos_tempmax
       if(use_temp_max_ext.eq.0) then
          enforceTableBounds = -1
       endif
    endif

    if(ye.lt.eos_yemin) then
       ye = eos_yemin
       if(use_ye_min_ext.eq.0) then
          enforceTableBounds = -1
       endif
    endif
    if(ye.gt.eos_yemax) then
       ye = eos_yemax
       if(use_ye_max_ext.eq.0) then
          enforceTableBounds = -1
       endif
    endif

    return
end function enforceTableBounds
