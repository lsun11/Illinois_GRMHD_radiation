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


#include <algorithm>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils_macro.h"

extern "C" void THC_InitSource(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_InitSource");
    }

    int const siz = UTILS_GFSIZE(cctkGH);

    if(CCTK_Equals(eos_type, "barotropic")) {
        std::fill(&rhs_dens[0], &rhs_dens[siz], 0.0);
    }
    if(CCTK_Equals(eos_type, "ideal")) {
        std::fill(&rhs_dens[0], &rhs_dens[siz], 0.0);
    }
    if(CCTK_Equals(eos_type, "nuclear")) {
        std::fill(&rhs_densxn[0], &rhs_densxn[siz], 0.0);
        std::fill(&rhs_densxp[0], &rhs_densxp[siz], 0.0);
    }

    std::fill(&rhs_scon[0], &rhs_scon[3*siz], 0.0);

    if(!CCTK_Equals(eos_type, "barotropic")) {
        std::fill(&rhs_tau [0], &rhs_tau[siz], 0.0);
    }
}
