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
#include <float.h>
#include <math.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils_macro.h"

#define SQ(X) ((X)*(X))

void THC_LK_UpdateOpticalDepth(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_LK_UpdateOpticalDepth");
    }

    int const oi = CCTK_GFINDEX3D(cctkGH,1,0,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
    int const oj = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
    int const ok = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);

    CCTK_REAL const dx = CCTK_DELTA_SPACE(0);
    CCTK_REAL const dy = CCTK_DELTA_SPACE(1);
    CCTK_REAL const dz = CCTK_DELTA_SPACE(2);

    CCTK_REAL * tau_vec[] = {
        optd_0_nue, optd_0_nua, optd_0_nux,
        optd_1_nue, optd_1_nua, optd_1_nux
    };
    CCTK_REAL * tau_vec_p[] = {
        optd_0_nue_p, optd_0_nua_p, optd_0_nux_p,
        optd_1_nue_p, optd_1_nua_p, optd_1_nux_p
    };
    CCTK_REAL const * kappa_vec[] = {
        kappa_0_nue, kappa_0_nua, kappa_0_nux,
        kappa_1_nue, kappa_1_nua, kappa_1_nux
    };

    /* Updates the optical detph as in Nielsen et al. 2014 */
#pragma omp parallel
    {
        UTILS_LOOP3(thc_lk_update_optd,
                k, cctk_nghostzones[2], cctk_lsh[2] - cctk_nghostzones[2],
                j, cctk_nghostzones[1], cctk_lsh[1] - cctk_nghostzones[1],
                i, cctk_nghostzones[0], cctk_lsh[0] - cctk_nghostzones[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            for(int v = 0; v < 6; ++v) {
                tau_vec[v][ijk] = DBL_MAX;
                for(int nk = -1; nk <= 1; ++nk)
                for(int nj = -1; nj <= 1; ++nj)
                for(int ni = -1; ni <= 1; ++ni) {
                    if(ni == 0 && nj == 0 && nk == 0) {
                        continue;
                    }
                    CCTK_REAL const delta[3] = {ni*dx, nj*dy, nk*dz};
                    CCTK_REAL const metric[9] = {
                        0.5*(gxx[ijk] + gxx[ijk + ni*oi + nj*oj + nk*ok]),
                        0.5*(gxy[ijk] + gxy[ijk + ni*oi + nj*oj + nk*ok]),
                        0.5*(gxz[ijk] + gxz[ijk + ni*oi + nj*oj + nk*ok]),
                        0.5*(gxy[ijk] + gxy[ijk + ni*oi + nj*oj + nk*ok]),
                        0.5*(gyy[ijk] + gyy[ijk + ni*oi + nj*oj + nk*ok]),
                        0.5*(gyz[ijk] + gyz[ijk + ni*oi + nj*oj + nk*ok]),
                        0.5*(gxz[ijk] + gxz[ijk + ni*oi + nj*oj + nk*ok]),
                        0.5*(gyz[ijk] + gyz[ijk + ni*oi + nj*oj + nk*ok]),
                        0.5*(gzz[ijk] + gzz[ijk + ni*oi + nj*oj + nk*ok])
                    };
                    CCTK_REAL dl2 = 0;
                    for(int i = 0; i < 3; ++i)
                    for(int j = 0; j < 3; ++j) {
                        dl2 += metric[3*i + j]*delta[i]*delta[j];
                    }
                    assert(dl2 > 0);
                    CCTK_REAL const dl = sqrt(dl2);

                    CCTK_REAL const kappa_avg = 0.5*(kappa_vec[v][ijk] +
                            kappa_vec[v][ijk + ni*oi + nj*oj + nk*ok]);

                    CCTK_REAL const path_optd = kappa_avg*dl +
                        tau_vec_p[v][ijk + ni*oi + nj*oj + nk*ok];

                    if(tau_vec[v][ijk] > path_optd) {
                        tau_vec[v][ijk] = path_optd;
                    }
                }
                assert(tau_vec[v][ijk] < DBL_MAX);
            }
        } UTILS_ENDLOOP3(thc_lk_update_optd);
    }
}
