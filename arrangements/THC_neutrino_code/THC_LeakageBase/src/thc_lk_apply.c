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
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils_macro.h"

#define POW2(X) ((X)*(X))

void THC_LK_NoAbsorption(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_LK_NoAbsorption");
    }

    int const siz = UTILS_GFSIZE(cctkGH);
    memset(abs_number, 0, siz*sizeof(CCTK_REAL));
    memset(abs_energy, 0, siz*sizeof(CCTK_REAL));
}

void THC_LK_Apply(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_LK_Apply");
    }

    int const siz = UTILS_GFSIZE(cctkGH);
    CCTK_REAL const * velx = &vel[0*siz];
    CCTK_REAL const * vely = &vel[1*siz];
    CCTK_REAL const * velz = &vel[2*siz];

    CCTK_REAL * sconx = &scon[0*siz];
    CCTK_REAL * scony = &scon[1*siz];
    CCTK_REAL * sconz = &scon[2*siz];

    CCTK_REAL const mb = AverageBaryonMass();

    CCTK_REAL const dt = CCTK_DELTA_TIME;

#pragma omp parallel
    {
        UTILS_LOOP3(thc_lk_apply,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            CCTK_REAL const v_x = gxx[ijk] * velx[ijk] + gxy[ijk] * vely[ijk] +
              gxz[ijk] * velz[ijk];
            CCTK_REAL const v_y = gxy[ijk] * velx[ijk] + gyy[ijk] * vely[ijk] +
              gyz[ijk] * velz[ijk];
            CCTK_REAL const v_z = gxz[ijk] * velx[ijk] + gyz[ijk] * vely[ijk] +
              gzz[ijk] * velz[ijk];

            if(rho[ijk] > density_cutoff) {
                /* Sanity checks */
                assert(isfinite(abs_number[ijk]));
                assert(isfinite(abs_energy[ijk]));

                /* Total lepton number / energy absorption (emission) */
                CCTK_REAL const R = abs_number[ijk] -
                    mb*(R_eff_nue[ijk] - R_eff_nua[ijk]);
                assert(isfinite(R));
                CCTK_REAL const Q = abs_energy[ijk] -
                    (Q_eff_nue[ijk] + Q_eff_nua[ijk] + Q_eff_nux[ijk]);
                assert(isfinite(Q));

                /* Compute RHS */
                CCTK_REAL const densxp_dot = alp[ijk] * volform[ijk] * R;
                CCTK_REAL const sconx_dot = alp[ijk] * volform[ijk] *
                    w_lorentz[ijk] * v_x * Q;
                CCTK_REAL const scony_dot = alp[ijk] * volform[ijk] *
                    w_lorentz[ijk] * v_y * Q;
                CCTK_REAL const sconz_dot = alp[ijk] * volform[ijk] *
                    w_lorentz[ijk] * v_z * Q;
                CCTK_REAL const tau_dot = alp[ijk] * volform[ijk] *
                    w_lorentz[ijk] * Q;

                /* Update composition */
                CCTK_REAL const dens_tot = densxp[ijk] + densxn[ijk];
                densxp[ijk] = densxp[ijk] / (1.0 - densxp_dot/densxp[ijk] * dt);
                densxn[ijk] = dens_tot - densxp[ijk];

                /* Update momentum */
                sconx[ijk] = sconx[ijk] + dt * sconx_dot;
                scony[ijk] = scony[ijk] + dt * scony_dot;
                sconz[ijk] = sconz[ijk] + dt * sconz_dot;

                /* Update the energy */
                tau[ijk] = tau[ijk] / (1.0 - tau_dot/tau[ijk] * dt);

                /* Store luminosity on the volume */
                if(store_neu_luminosity) {
                    luminosity_nue[ijk] = alp[ijk] *
                      alp[ijk] * volform[ijk] * w_lorentz[ijk] * Q_eff_nue[ijk];
                    luminosity_nua[ijk] = alp[ijk] *
                      alp[ijk] * volform[ijk] * w_lorentz[ijk] * Q_eff_nua[ijk];
                    luminosity_nux[ijk] = alp[ijk] *
                      alp[ijk] * volform[ijk] * w_lorentz[ijk] * Q_eff_nux[ijk];
                }
            }
        } UTILS_ENDLOOP3(thc_lk_apply);
    }
}
