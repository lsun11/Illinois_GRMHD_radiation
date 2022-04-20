//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2013, David Radice <dradice@caltech.edu>
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


#include <limits>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "thc_atmosphere.hh"
#include "thc_claws.hh"

using namespace thc;

extern "C" void THC_HRSCCRegister(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_HRSCCRegister");
    }

    if(CCTK_Equals(eos_type, "barotropic")) {
        EulerBase<barotropic::Euler>::get_var_idx();
        barotropic::Euler::conserved_lbound[0] = atmosphere::rho;
    }
    else if(CCTK_Equals(eos_type, "ideal")) {
        EulerBase<ideal::Euler>::get_var_idx();
        ideal::Euler::conserved_lbound[0] = atmosphere::rho;
    }
    else if(CCTK_Equals(eos_type, "nuclear")) {
        EulerBase<nuclear::Euler>::get_var_idx();
        nuclear::Euler::conserved_lbound[0] = atmosphere::rho;
        nuclear::Euler::conserved_lbound[1] = atmosphere::rho;
    }
    else if(CCTK_Equals(eos_type, "ultrarelativistic")) {
        EulerBase<ultrarelativistic::Euler>::get_var_idx();
        ultrarelativistic::Euler::conserved_lbound[0] = atmosphere::rho;
        ultrarelativistic::Euler::gamma = ultrarelativistic_gamma;
    }
    else {
        abort();
    }
}
