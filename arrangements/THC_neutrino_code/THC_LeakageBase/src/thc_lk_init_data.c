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


#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils_macro.h"

void THC_LK_InitData(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

      /*
    if(verbose) {
        CCTK_INFO("THC_LK_InitData");
    }
      */
    int gfsiz = UTILS_GFSIZE(cctkGH);

    /* Initialize everything to zero apart from the optical depth that is
     * initialized in a more complex way */
    memset(kappa_0_nue,    0, gfsiz*sizeof(CCTK_REAL));
    memset(kappa_0_nua,    0, gfsiz*sizeof(CCTK_REAL));
    memset(kappa_0_nux,    0, gfsiz*sizeof(CCTK_REAL));
    memset(kappa_1_nue,    0, gfsiz*sizeof(CCTK_REAL));
    memset(kappa_1_nua,    0, gfsiz*sizeof(CCTK_REAL));
    memset(kappa_1_nux,    0, gfsiz*sizeof(CCTK_REAL));
    if(store_free_rates) {
        memset(R_free_nue,     0, gfsiz*sizeof(CCTK_REAL));
        memset(R_free_nua,     0, gfsiz*sizeof(CCTK_REAL));
        memset(R_free_nux,     0, gfsiz*sizeof(CCTK_REAL));
        memset(Q_free_nue,     0, gfsiz*sizeof(CCTK_REAL));
        memset(Q_free_nua,     0, gfsiz*sizeof(CCTK_REAL));
        memset(Q_free_nux,     0, gfsiz*sizeof(CCTK_REAL));
    }
    memset(R_eff_nue,      0, gfsiz*sizeof(CCTK_REAL));
    memset(R_eff_nua,      0, gfsiz*sizeof(CCTK_REAL));
    memset(R_eff_nux,      0, gfsiz*sizeof(CCTK_REAL));
    memset(Q_eff_nue,      0, gfsiz*sizeof(CCTK_REAL));
    memset(Q_eff_nua,      0, gfsiz*sizeof(CCTK_REAL));
    memset(Q_eff_nux,      0, gfsiz*sizeof(CCTK_REAL));
    memset(abs_number,     0, gfsiz*sizeof(CCTK_REAL));
    memset(abs_energy,     0, gfsiz*sizeof(CCTK_REAL));
    if(store_neu_luminosity) {
        memset(luminosity_nue, 0, gfsiz*sizeof(CCTK_REAL));
        memset(luminosity_nua, 0, gfsiz*sizeof(CCTK_REAL));
        memset(luminosity_nux, 0, gfsiz*sizeof(CCTK_REAL));
    }
}
