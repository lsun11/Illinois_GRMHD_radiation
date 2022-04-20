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


#ifndef THC_LEAKAGEBASE_RATES_H
#define THC_LEAKAGEBASE_RATES_H

#include "cctk.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Computes the effective emission rate as in Ruffert & Janka 1996
 */
CCTK_REAL thc_lk_calc_eff_rate(
        CCTK_REAL const R_free,         /* free rate */
        CCTK_REAL const dens,           /* number or energy density */
        CCTK_REAL const kappa,          /* local opacity */
        CCTK_REAL const tau);           /* optical depth */

#ifdef __cplusplus
}
#endif

#endif
