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

extern "C" void THC_MoLRegister(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_MoLRegister");
    }

    int ierr = 0;

    if(!CCTK_Equals(evolution_method, "THCode")) {
        return;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Register evolved variables
    ///////////////////////////////////////////////////////////////////////////
    if(CCTK_Equals(eos_type, "barotropic")) {
        ierr |= MoLRegisterEvolved(CCTK_VarIndex("THC_Core::dens"),
                CCTK_VarIndex("THC_Core::rhs_dens"));
    }
    if(CCTK_Equals(eos_type, "ideal")) {
        ierr |= MoLRegisterEvolved(CCTK_VarIndex("THC_Core::dens"),
                CCTK_VarIndex("THC_Core::rhs_dens"));
    }
    if(CCTK_Equals(eos_type, "nuclear")) {
        ierr |= MoLRegisterEvolved(CCTK_VarIndex("THC_Core::densxn"),
                CCTK_VarIndex("THC_Core::rhs_densxn"));
        ierr |= MoLRegisterEvolved(CCTK_VarIndex("THC_Core::densxp"),
                CCTK_VarIndex("THC_Core::rhs_densxp"));
    }

    ierr |= MoLRegisterEvolved(CCTK_VarIndex("THC_Core::scon[0]"),
            CCTK_VarIndex("THC_Core::rhs_scon[0]"));
    ierr |= MoLRegisterEvolved(CCTK_VarIndex("THC_Core::scon[1]"),
            CCTK_VarIndex("THC_Core::rhs_scon[1]"));
    ierr |= MoLRegisterEvolved(CCTK_VarIndex("THC_Core::scon[2]"),
            CCTK_VarIndex("THC_Core::rhs_scon[2]"));

    if(!CCTK_Equals(eos_type, "barotropic")) {
        ierr |= MoLRegisterEvolved(CCTK_VarIndex("THC_Core::tau"),
                CCTK_VarIndex("THC_Core::rhs_tau"));
    }

    if(!CCTK_Equals(eos_type, "ultrarelativistic")) {
        ierr |= MoLRegisterEvolved(CCTK_VarIndex("THC_Core::densgain"),
                CCTK_VarIndex("THC_Core::rhs_densgain"));
    }
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // Register constrained variables
    ///////////////////////////////////////////////////////////////////////////
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("HydroBase::rho"));
    if(CCTK_Equals(eos_type, "nuclear")) {
        ierr |= MoLRegisterConstrained(CCTK_VarIndex("THC_Core::dens"));
    }

    ierr |= MoLRegisterConstrained(CCTK_VarIndex("HydroBase::vel[0]"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("HydroBase::vel[1]"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("HydroBase::vel[2]"));
    if(!CCTK_Equals(eos_type, "ultrarelativistic")) {
        ierr |= MoLRegisterConstrained(CCTK_VarIndex("THC_Core::zvec[0]"));
        ierr |= MoLRegisterConstrained(CCTK_VarIndex("THC_Core::zvec[1]"));
        ierr |= MoLRegisterConstrained(CCTK_VarIndex("THC_Core::zvec[2]"));
    }

    ierr |= MoLRegisterConstrained(CCTK_VarIndex("HydroBase::eps"));

    ierr |= MoLRegisterConstrained(CCTK_VarIndex("HydroBase::press"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("HydroBase::w_lorentz"));
    if(CCTK_Equals(eos_type, "nuclear")) {
        ierr |= MoLRegisterConstrained(CCTK_VarIndex("HydroBase::Y_e"));
        ierr |= MoLRegisterConstrained(CCTK_VarIndex("HydroBase::entropy"));
        ierr |= MoLRegisterConstrained(CCTK_VarIndex("HydroBase::temperature"));
    }
    ///////////////////////////////////////////////////////////////////////////

    if(ierr) {
        CCTK_ERROR("Could not register with MoL");
    }
}
