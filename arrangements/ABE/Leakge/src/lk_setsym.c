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


#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "Symmetry.h"

void LK_SetSym(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    int const sym_scal[3] = {1, 1, 1};
    SetCartSymGN(cctkGH, sym_scal, "Leakage::leakage_optd");
    SetCartSymGN(cctkGH, sym_scal, "Leakage::leakage_eff_rates");
}
