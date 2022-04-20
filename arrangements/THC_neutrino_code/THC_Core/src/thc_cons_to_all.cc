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
#include <cassert>
#include <cstring>
#include <limits>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "hrscc_config.hh"
#include "utils.hh"

#include "thc_claws.hh"
#include "thc_macro.hh"
#include "thc_printer.hh"

namespace {

void cons_to_all(cGH const * const cctkGH) {
    DECLARE_CCTK_PARAMETERS

    if(CCTK_Equals(eos_type, "barotropic")) {
        hrscc::cons_to_all<thc::barotropic::Euler>(cctkGH);
    }
    else if(CCTK_Equals(eos_type, "ideal")) {
        hrscc::cons_to_all<thc::ideal::Euler>(cctkGH);
    }
    else if(CCTK_Equals(eos_type, "nuclear")) {
        hrscc::cons_to_all<thc::nuclear::Euler>(cctkGH);
    }
    else if(CCTK_Equals(eos_type, "ultrarelativistic")) {
        hrscc::cons_to_all<thc::ultrarelativistic::Euler>(cctkGH);
    }
    else {
        abort();
    }
}

}

extern "C" void THC_ConsToAll(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    static int substep_old = -1;

    if(verbose) {
        CCTK_INFO("THC_ConsToAll");
    }

    int const siz = UTILS_GFSIZE(cctkGH);
    if(hydro_excision) {
        for(int ijk = 0; ijk < siz; ++ijk) {
            bitmask[ijk] =
                hydro_excision_mask[ijk] ?  THC_FLAG_IS_INSIDE_BH : 0;
        }
    }
    else {
        std::fill(&bitmask[0], &bitmask[siz], 0);
    }

    thc::Printer::start(
            "[INFO|THC|ConsToAll]: ",
            "[WARN|THC|ConsToAll]: ",
            "[ERR|THC|ConsToAll]: ",
            c2a_max_num_msg, c2a_max_num_msg);
    cons_to_all(cctkGH);
    thc::Printer::stop();

    if(!CCTK_Equals(eos_type, "ultrarelativistic")) {
        int const substep = *static_cast<CCTK_INT *>(
                CCTK_VarDataPtr(cctkGH, 0, "MoL::MoL_Intermediate_Step"));
        assert(substep < max_num_substeps);
        // If substep == substep_old, we are not in a new RK substep, but
        // something else is happening (e.g., prolongation or restriction) and
        // we do not need to update c2a_densgain_db.
        if(substep >= 0 && substep != substep_old) {
            std::memcpy(&c2a_densgain_db[substep*siz],
                    &c2a_densgain[0], siz*sizeof(CCTK_REAL));
            substep_old = substep;
        }
    }
}
