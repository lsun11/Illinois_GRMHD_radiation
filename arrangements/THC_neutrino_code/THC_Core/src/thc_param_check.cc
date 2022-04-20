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


#include <cstdio>
#include <cstring>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define MIN(A, B) ((A) < (B) ? (A) : (B))
#define MINVALUE(X) (MIN(MIN(X[0], X[1]), X[2]))

extern "C" void THC_ParamCheck(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    // Check that the user knows what (s)he is doing
    if(!give_wrong_results) {
        CCTK_PARAMWARN("Correct results are not supported!");
    }

    // Check that the options are consistent with GRHD (if we are doing GRHD)
    if(CCTK_Equals(physics, "GRHD")) {
        if(CCTK_Equals(eos_type, "ultrarelativistic")) {
            CCTK_PARAMWARN("GRHD has not been implemented for "
                    "ultrarelativistic hydrodynamics!");
        }
        if(CCTK_Equals(scheme, "DG")) {
            CCTK_PARAMWARN("GRHD with DG has not been implemented yet!");
        }
    }
    // GRMHD is not implemented yet
    if(CCTK_Equals(physics, "GRMHD")) {
        CCTK_PARAMWARN("GRMHD has not been implemented yet!");
    }

    // Barotropic EoS evolution not been implemented yet
    if(CCTK_Equals(eos_type, "barotropic")) {
        CCTK_PARAMWARN("Barotropic EoS evolution has not been "
                "implemented yet!");
    }

    // For nuclear EoS we need entropy, temperature and Y_e storage to be active
    if(CCTK_Equals(eos_type, "nuclear")) {
        if(CCTK_Equals(initial_entropy, "None")) {
            CCTK_PARAMWARN("THC requires the entropy storage to be "
                "activated!");
        }
        if(CCTK_Equals(initial_temperature, "None")) {
            CCTK_PARAMWARN("THC requires the temperature storage to be "
                "activated!");
        }
        if(CCTK_Equals(initial_Y_e, "None")) {
            CCTK_PARAMWARN("THC requires the Y_e storage to be activated!");
        }
    }

    // All the parameters needed *only* for the evolution should be checked
    // after this point
    if(!CCTK_Equals(evolution_method, "THCode")) {
        return;
    }

    // For nuclear EoS we evolve entropy, temperature and Y_e
    if(CCTK_Equals(eos_type, "nuclear")) {
        if(!CCTK_Equals(entropy_evolution_method, "THCode")) {
            CCTK_PARAMWARN("entropy_evolution_method must be set to "
                    "\"THCode\"");
        }
        if(!CCTK_Equals(temperature_evolution_method, "THCode")) {
            CCTK_PARAMWARN("temperature_evolution_method must be set to "
                    "\"THCode\"");
        }
        if(!CCTK_Equals(Y_e_evolution_method, "THCode")) {
            CCTK_PARAMWARN("Y_e_evolution_method must be set to "
                    "\"THCode\"");
        }
    }
    else {
        if(CCTK_Equals(entropy_evolution_method, "THCode")) {
            CCTK_PARAMWARN("entropy_evolution_method is set to "
                    "\"THCode\", but entropy can be evolved only with "
                    "\"nuclear\" eos_type");
        }
        if(CCTK_Equals(temperature_evolution_method, "THCode")) {
            CCTK_PARAMWARN("temperature_evolution_method is set to "
                    "\"THCode\", but temperture can be evolved only with "
                    "\"nuclear\" eos_type");
        }
        if(CCTK_Equals(Y_e_evolution_method, "THCode")) {
            CCTK_PARAMWARN("Y_e_evolution_method is set to \"THCode\", "
                    "but Y_e can be evolved only with \"nuclear\" eos_type");
        }
    }

    // Check that we have enough ghostpoints for GR RHS
    if(CCTK_Equals(physics, "GRHD")) {
        int nghostzones = MINVALUE(cctk_nghostzones);
        if(fd_order > nghostzones) {
            char buf[BUFSIZ];
            std::snprintf(buf, BUFSIZ, "You have set fd_order = %d, this "
                    "requires at least %d ghostzones, but "
                    "cctk_nghostzones = %d", fd_order, 2*fd_order, nghostzones);
            CCTK_PARAMWARN(buf);
        }
    }

    // Check if we are using positivity preserving
    if(pplim) {
        if(atmo_tolerance > 0 && !atmo_I_know_what_I_am_doing) {
            CCTK_PARAMWARN("The positivity preserving limiter only works "
                    "if atmo_tolerance = 0");
        }
    }

    // Check that if refluxing is requested everything needed is in place
    if(refluxing) {
        if(!CCTK_IsThornActive("Refluxing") ||
                !CCTK_IsThornActive("THC_Refluxing")) {
            CCTK_PARAMWARN("Refluxing requres the thorns \"Refluxing\" and "
                    "\"THC_Refluxing\" to be active");
        }
    }

    // Check that we have enough space to store the c2a_densgain history
    if(!CCTK_Equals(eos_type, "ultrarelativistic")) {
        CCTK_INT const * MoL_Intermediate_Steps = static_cast<CCTK_INT const *>(
                CCTK_ParameterGet("MoL_Intermediate_Steps", "MoL", NULL));
        if(*MoL_Intermediate_Steps > max_num_substeps) {
            CCTK_PARAMWARN("THC_Core::max_num_substeps should be larger than "
                    "MoL::MoL_Intermediate_Steps");
        }
    }
}
