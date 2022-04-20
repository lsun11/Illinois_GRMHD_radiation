//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2013, David Radice <david.radice@aei.mpg.de>
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

#include "global_eos.h"

#include "thc_atmosphere.hh"
#include "thc_claws.hh"

CCTK_REAL thc::atmosphere::rho     = 0;
CCTK_REAL thc::atmosphere::eps     = 0;
CCTK_REAL thc::atmosphere::Ye      = 0;

CCTK_REAL thc::atmosphere::csound  = 0;
CCTK_REAL thc::atmosphere::entropy = 0;
CCTK_REAL thc::atmosphere::temp    = 0;
CCTK_REAL thc::atmosphere::press   = 0;

CCTK_REAL thc::atmosphere::tol     = 0;
bool thc::atmosphere::dynamic      = false;

extern "C" void THC_SetupAtmosphere(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_SetupAtmosphere");
    }

    thc::atmosphere::rho      = atmo_rho;
    thc::atmosphere::eps      = atmo_eps;
    thc::atmosphere::Ye       = atmo_Y_e;

    thc::atmosphere::csound   = 0;
    thc::atmosphere::entropy  = 0;
    thc::atmosphere::press    = 0;
    thc::atmosphere::temp     = atmo_temperature;

    thc::atmosphere::tol      = atmo_tolerance;
    thc::atmosphere::dynamic  = atmo_is_dynamic;

    if(CCTK_Equals(eos_type, "barotropic")) {
        CCTK_ERROR("Barotropic atmosphere not implemented yet!");
    }
    else if(CCTK_Equals(eos_type, "ideal")) {
        whizza::eos_thermal_impl const & eos =
            whizza::global_eos_thermal::get_eos().implementation();
        whizza::eos_thermal_impl::status status;
        eos.check_rho_eps_ye(thc::atmosphere::rho, thc::atmosphere::eps,
                0.5, status);
        if(status.failed) {
            CCTK_ERROR(status.err_msg.c_str());
        }
        thc::atmosphere::press = eos.press_from_valid_rho_eps_ye(
                thc::atmosphere::rho, thc::atmosphere::eps, 0.5);
        thc::atmosphere::csound = eos.csnd_from_valid_rho_eps_ye(
                thc::atmosphere::rho, thc::atmosphere::eps, 0.5);
    }
    else if(CCTK_Equals(eos_type, "nuclear")) {
        whizza::eos_thermal_impl const & eos =
            whizza::global_eos_thermal::get_eos().implementation();
        whizza::eos_thermal_impl::status status;
        eos.check_rho_temp_ye(thc::atmosphere::rho, thc::atmosphere::temp,
                thc::atmosphere::Ye, status);
        if(status.failed) {
            CCTK_ERROR(status.err_msg.c_str());
        }
        thc::atmosphere::eps = eos.eps_from_valid_rho_temp_ye(
                thc::atmosphere::rho, thc::atmosphere::temp,
                thc::atmosphere::Ye);
        thc::atmosphere::csound = eos.csnd_from_valid_rho_eps_ye(
                thc::atmosphere::rho, thc::atmosphere::eps,
                thc::atmosphere::Ye);
        thc::atmosphere::entropy = eos.entropy_from_valid_rho_eps_ye(
                thc::atmosphere::rho, thc::atmosphere::eps,
                thc::atmosphere::Ye);
        thc::atmosphere::press = eos.press_from_valid_rho_eps_ye(
                thc::atmosphere::rho, thc::atmosphere::eps,
                thc::atmosphere::Ye);
    }
    else if(CCTK_Equals(eos_type, "ultrarelativistic")) {
        thc::atmosphere::press = (ultrarelativistic_gamma - 1.0)*
            thc::atmosphere::rho;
    }
    else {
        CCTK_ERROR("Unknown eos_type!");
    }
}
