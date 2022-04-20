//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2015, David Radice <dradice@caltech.edu>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include <assert.h>
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils_macro.h"

void THC_LK_CalcOpacity(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_LK_CalcOpacity");
    }

#pragma omp parallel
    {
        UTILS_LOOP3(thc_lk_calc_opacity,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            int ierr = NeutrinoOpacity(
                    rho[ijk], temperature[ijk], Y_e[ijk],
                    &kappa_0_nue[ijk],
                    &kappa_0_nua[ijk],
                    &kappa_0_nux[ijk],
                    &kappa_1_nue[ijk],
                    &kappa_1_nua[ijk],
                    &kappa_1_nux[ijk]);
            assert(!ierr);
            assert(isfinite(kappa_0_nue[ijk]));
            assert(isfinite(kappa_0_nua[ijk]));
            assert(isfinite(kappa_0_nux[ijk]));
            assert(isfinite(kappa_1_nue[ijk]));
            assert(isfinite(kappa_1_nua[ijk]));
            assert(isfinite(kappa_1_nux[ijk]));
        } UTILS_ENDLOOP3(thc_lk_calc_opacity);
    }
}
