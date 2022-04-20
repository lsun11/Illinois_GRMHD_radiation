//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2012, David Radice <david.radice@aei.mpg.de>
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


#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils.hh"
#include "thc_macro.hh"

// Adapted from GRHydro_Tmunu.F90
extern "C" void THC_AddToTmunu(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_AddToTmunu");
    }

    int const siz = UTILS_GFSIZE(cctkGH);

    CCTK_REAL const * velx = &vel[0*siz];
    CCTK_REAL const * vely = &vel[1*siz];
    CCTK_REAL const * velz = &vel[2*siz];

#pragma omp parallel
    {
        UTILS_LOOP3(thc_tmunu,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            CCTK_REAL const atmof = UTILS_BITMASK_CHECK_FLAG(bitmask[ijk],
                    THC_FLAG_ATMOSPHERE) ? 0.0 : 1.0;

            CCTK_REAL const rho_h = rho[ijk]*(1 + eps[ijk]) + press[ijk];

            CCTK_REAL const v_x = gxx[ijk]*velx[ijk] + gxy[ijk]*vely[ijk] +
                gxz[ijk]*velz[ijk];
            CCTK_REAL const v_y = gxy[ijk]*velx[ijk] + gyy[ijk]*vely[ijk] +
                gyz[ijk]*velz[ijk];
            CCTK_REAL const v_z = gxz[ijk]*velx[ijk] + gyz[ijk]*vely[ijk] +
                gzz[ijk]*velz[ijk];

            CCTK_REAL const beta_x = gxx[ijk]*betax[ijk] + gxy[ijk]*betay[ijk] +
                gxz[ijk]*betaz[ijk];
            CCTK_REAL const beta_y = gxy[ijk]*betax[ijk] + gyy[ijk]*betay[ijk] +
                gyz[ijk]*betaz[ijk];
            CCTK_REAL const beta_z = gxz[ijk]*betax[ijk] + gyz[ijk]*betay[ijk] +
                gzz[ijk]*betaz[ijk];
            CCTK_REAL const beta2 = beta_x*betax[ijk] + beta_y*betay[ijk] +
                beta_z*betaz[ijk];

            CCTK_REAL const u_t = w_lorentz[ijk]*(-alp[ijk] + beta_x*velx[ijk] +
                    beta_y*vely[ijk] + beta_z*velz[ijk]);
            CCTK_REAL const u_x = w_lorentz[ijk]*v_x;
            CCTK_REAL const u_y = w_lorentz[ijk]*v_y;
            CCTK_REAL const u_z = w_lorentz[ijk]*v_z;

            eTtt[ijk] += atmof * (rho_h * u_t * u_t + press[ijk]*(beta2 -
                    SQ(alp[ijk])));

            eTtx[ijk] += atmof * (rho_h * u_t * u_x + press[ijk]*beta_x);
            eTty[ijk] += atmof * (rho_h * u_t * u_y + press[ijk]*beta_y);
            eTtz[ijk] += atmof * (rho_h * u_t * u_z + press[ijk]*beta_z);

            eTxx[ijk] += atmof * (rho_h * u_x * u_x + press[ijk]*gxx[ijk]);
            eTxy[ijk] += atmof * (rho_h * u_x * u_y + press[ijk]*gxy[ijk]);
            eTxz[ijk] += atmof * (rho_h * u_x * u_z + press[ijk]*gxz[ijk]);
            eTyy[ijk] += atmof * (rho_h * u_y * u_y + press[ijk]*gyy[ijk]);
            eTyz[ijk] += atmof * (rho_h * u_y * u_z + press[ijk]*gyz[ijk]);
            eTzz[ijk] += atmof * (rho_h * u_z * u_z + press[ijk]*gzz[ijk]);
        } UTILS_ENDLOOP3(thc_tmunu);
    }
}
