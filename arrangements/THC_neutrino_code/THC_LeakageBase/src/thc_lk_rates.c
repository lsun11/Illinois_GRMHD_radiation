//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2014, David Radice <dradice@caltech.edu>
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
#include <float.h>
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils_macro.h"
#include "thc_lk_rates.h"

#define POW2(X) ((X)*(X))

CCTK_REAL thc_lk_calc_eff_rate(
        CCTK_REAL const R_free,
        CCTK_REAL const dens,
        CCTK_REAL const kappa,
        CCTK_REAL const tau) {
    DECLARE_CCTK_PARAMETERS

    if(dens <= DBL_EPSILON*R_free || kappa < DBL_EPSILON) {
        return 0;
    }
    CCTK_REAL const itloss = R_free/dens;
    CCTK_REAL const lambda = 1.0/kappa;
    CCTK_REAL const tdiff  = DiffFact * lambda * POW2(tau);

    return R_free/(1 + itloss * tdiff);
}

void THC_LK_Rates(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_LK_Rates");
    }

#pragma omp parallel
    {
        UTILS_LOOP3(thc_lk_rates,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            /* Compute free rates */
            CCTK_REAL R_nue, R_nua, R_nux;
            CCTK_REAL Q_nue, Q_nua, Q_nux;
            int ierr = NeutrinoEmission(rho[ijk], temperature[ijk], Y_e[ijk],
                    &R_nue, &R_nua, &R_nux, &Q_nue, &Q_nua, &Q_nux);
            assert(!ierr);
            assert(isfinite(R_nue));
            assert(isfinite(R_nua));
            assert(isfinite(R_nux));
            assert(isfinite(Q_nue));
            assert(isfinite(Q_nua));
            assert(isfinite(Q_nux));

            /* Store the free rates if so requested */
            if(store_free_rates) {
                R_free_nue[ijk] = R_nue;
                R_free_nua[ijk] = R_nua;
                R_free_nux[ijk] = R_nux;
                Q_free_nue[ijk] = Q_nue;
                Q_free_nua[ijk] = Q_nua;
                Q_free_nux[ijk] = Q_nux;
            }

            /* Computes effective rates */
            CCTK_REAL num_nue, num_nua, num_nux;
            CCTK_REAL ene_nue, ene_nua, ene_nux;
            ierr = NeutrinoDensity(rho[ijk], temperature[ijk], Y_e[ijk],
                    &num_nue, &num_nua, &num_nux,
                    &ene_nue, &ene_nua, &ene_nux);
            assert(!ierr);
            assert(isfinite(num_nue));
            assert(isfinite(num_nua));
            assert(isfinite(num_nux));
            assert(isfinite(ene_nue));
            assert(isfinite(ene_nua));
            assert(isfinite(ene_nux));

            R_eff_nue[ijk] = thc_lk_calc_eff_rate(R_nue, num_nue,
                    kappa_0_nue[ijk], optd_0_nue[ijk]);
            R_eff_nua[ijk] = thc_lk_calc_eff_rate(R_nua, num_nua,
                    kappa_0_nua[ijk], optd_0_nua[ijk]);
            R_eff_nux[ijk] = thc_lk_calc_eff_rate(R_nux, num_nux,
                    kappa_0_nux[ijk], optd_0_nux[ijk]);
            Q_eff_nue[ijk] = thc_lk_calc_eff_rate(Q_nue, ene_nue,
                    kappa_1_nue[ijk], optd_1_nue[ijk]);
            Q_eff_nua[ijk] = thc_lk_calc_eff_rate(Q_nua, ene_nua,
                    kappa_1_nua[ijk], optd_1_nua[ijk]);
            Q_eff_nux[ijk] = thc_lk_calc_eff_rate(Q_nux, ene_nux,
                    kappa_1_nux[ijk], optd_1_nux[ijk]);
            assert(isfinite(R_eff_nue[ijk]));
            assert(isfinite(R_eff_nua[ijk]));
            assert(isfinite(R_eff_nux[ijk]));
            assert(isfinite(Q_eff_nue[ijk]));
            assert(isfinite(Q_eff_nua[ijk]));
            assert(isfinite(Q_eff_nux[ijk]));
        } UTILS_ENDLOOP3(thc_lk_rates);
    }
}
