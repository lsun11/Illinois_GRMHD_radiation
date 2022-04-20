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


#include <cassert>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "hrscc_config.hh"
#include "utils.hh"

#include "thc_claws.hh"
#include "thc_printer.hh"

namespace {

void rhs(cGH const * const cctkGH) {
    DECLARE_CCTK_PARAMETERS

    if(CCTK_Equals(eos_type, "barotropic")) {
        hrscc::compute_rhs<thc::barotropic::Euler>(cctkGH);
    }
    else if(CCTK_Equals(eos_type, "ideal")) {
        hrscc::compute_rhs<thc::ideal::Euler>(cctkGH);
    }
    else if(CCTK_Equals(eos_type, "nuclear")) {
        hrscc::compute_rhs<thc::nuclear::Euler>(cctkGH);
    }
    else if(CCTK_Equals(eos_type, "ultrarelativistic")) {
        hrscc::compute_rhs<thc::ultrarelativistic::Euler>(cctkGH);
    }
    else {
        abort();
    }
}

}

extern "C" void THC_RHS(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_RHS");
    }

    thc::Printer::start(
            "[INFO|THC|RHS]: ",
            "[WARN|THC|RHS]: ",
            "[ERR|THC|RHS]: ");
    rhs(cctkGH);
    thc::Printer::stop();

    if(!CCTK_Equals(eos_type, "ultrarelativistic")) {
        int const gsiz = UTILS_GFSIZE(cctkGH);
        CCTK_REAL const idt = 1.0 / CCTK_DELTA_TIME;

        CCTK_INT const * MoL_Intermediate_Step = static_cast<CCTK_INT *>(
                CCTK_VarDataPtr(cctkGH, 0, "MoL::MoL_Intermediate_Step"));
        assert(*MoL_Intermediate_Step > 0);
        int const offset = (*MoL_Intermediate_Step - 1)*gsiz;

        for(int ijk = 0; ijk < gsiz; ++ijk) {
            rhs_densgain[ijk] = idt * c2a_densgain_db[offset + ijk];
        }
    }
}
