//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2016, David Radice <dradice@caltech.edu>
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

#include "utils.hh"

extern "C"
void THC_ResetDensgain(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    int const siz = UTILS_GFSIZE(cctkGH);
    if(timelevels > 0) {
        std::fill(&densgain[0], &densgain[siz], 0.0);
    }
    if(timelevels > 1) {
        std::fill(&densgain_p[0], &densgain_p[siz], 0.0);
    }
    if(timelevels > 2) {
        std::fill(&densgain_p_p[0], &densgain_p_p[siz], 0.0);
    }
}
