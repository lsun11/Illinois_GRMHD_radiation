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


#include <algorithm>
#include <cmath>
#include <cstring>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "hrscc_config.hh"
#include "utils.hh"

#include "thc_atmosphere.hh"
#include "thc_claws.hh"
#include "thc_printer.hh"

namespace {

void prim_to_all(cGH const * const cctkGH) {
    DECLARE_CCTK_PARAMETERS

    if(CCTK_Equals(eos_type, "barotropic")) {
        hrscc::prim_to_all<thc::barotropic::Euler>(cctkGH);
    }
    else if(CCTK_Equals(eos_type, "ideal")) {
        hrscc::prim_to_all<thc::ideal::Euler>(cctkGH);
    }
    else if(CCTK_Equals(eos_type, "nuclear")) {
        hrscc::prim_to_all<thc::nuclear::Euler>(cctkGH);
    }
    else if(CCTK_Equals(eos_type, "ultrarelativistic")) {
        hrscc::prim_to_all<thc::ultrarelativistic::Euler>(cctkGH);
    }
    else {
        abort();
    }
}

}

extern "C" void THC_PrimToAll(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_PrimToAll");
    }

    // Setup bitmask
    int const siz = UTILS_GFSIZE(cctkGH);
    std::fill(&bitmask[0], &bitmask[siz], 0);

#ifdef THC_DEBUG
    // Initialize conservatives
    if(!CCTK_Equals(eos_type, "ultrarelativistic")) {
        std::fill(&dens[0], &dens[siz], 0.0);
    }
    if(CCTK_Equals(eos_type, "nuclear")) {
        std::fill(&densxn[0], &densxn[siz], 0.0);
        std::fill(&densxp[0], &densxp[siz], 0.0);
    }
    std::fill(&scon[0], &scon[3*siz], 0.0);
    if(!CCTK_Equals(eos_type, "barotropic")) {
        std::fill(&tau[0], &tau[siz], 0.0);
    }

    // Initialize extra variables
    std::fill(&press[0], &press[siz], 0.0);
    std::fill(&w_lorentz[0], &w_lorentz[siz], 1.0);
    std::fill(&csound[0], &csound[siz], 0.0);
    if(!CCTK_Equals(eos_type, "ultrarelativistic")) {
        std::fill(&c2a_nbiter[0], &c2a_nbiter[siz], 0.0);
        std::fill(&volform[0], &volform[siz], 0.0);
    }
    if(CCTK_Equals(eos_type, "nuclear")) {
        std::fill(&entropy[0], &entropy[siz], 0.0);
        std::fill(&temperature[0], &temperature[siz], 0.0);
    }
#endif

    thc::Printer::start(
            "[INFO|THC|PrimToAll]: ",
            "[WARN|THC|PrimToAll]: ",
            "[ERR|THC|PrimToAll]: ");
    prim_to_all(cctkGH);
    thc::Printer::stop();
}

extern "C" void THC_PrimToAllInitial(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_PrimToAllInitial");
    }

    int const siz = UTILS_GFSIZE(cctkGH);

    bool const nuclear = CCTK_Equals(eos_type, "nuclear");
    bool const ultra = CCTK_Equals(eos_type, "ultrarelativistic");

    // Setup bitmask
    std::fill(&bitmask[0], &bitmask[siz], 0);

    // Setup atmosphere
    for(int ijk = 0; ijk < siz; ++ijk) {
        if(rho[ijk] < thc::atmosphere::rho*(1 + thc::atmosphere::tol)) {
            UTILS_BITMASK_SET_FLAG(bitmask[ijk], THC_FLAG_ATMOSPHERE);
            rho[ijk]         = thc::atmosphere::rho;
            vel[ijk + 0*siz] = 0;
            vel[ijk + 1*siz] = 0;
            vel[ijk + 2*siz] = 0;
            eps[ijk]         = thc::atmosphere::eps;
            if(nuclear) {
                Y_e[ijk]     = thc::atmosphere::Ye;
            }
        }
    }

    // Compute zvec if necessary
    if(!ultra) {
        for(int ijk = 0; ijk < siz; ++ijk) {
            CCTK_REAL const velx = vel[ijk + 0*siz];
            CCTK_REAL const vely = vel[ijk + 1*siz];
            CCTK_REAL const velz = vel[ijk + 2*siz];

            CCTK_REAL const v_x = gxx[ijk]*velx + gxy[ijk]*vely
                + gxz[ijk]*velz;
            CCTK_REAL const v_y = gxy[ijk]*velx + gyy[ijk]*vely
                + gyz[ijk]*velz;
            CCTK_REAL const v_z = gxz[ijk]*velx + gyz[ijk]*vely
                + gzz[ijk]*velz;

            CCTK_REAL const v2 = v_x*velx + v_y*vely + v_z*velz;
            CCTK_REAL const W  = std::sqrt(1.0/(1.0 - v2));

            zvec[ijk + 0*siz] = W*velx;
            zvec[ijk + 1*siz] = W*vely;
            zvec[ijk + 2*siz] = W*velz;
        }
    }

    // Initialize dens gain
    if(!ultra) {
        std::fill(&densgain[0], &densgain[siz], 0.0);
        std::fill(&c2a_densgain_db[0],
                &c2a_densgain_db[siz*max_num_substeps], 0.0);
    }

#ifdef THC_DEBUG
    // Initialize conservatives
    if(!CCTK_Equals(eos_type, "ultrarelativistic")) {
        std::fill(&dens[0], &dens[siz], 0.0);
    }
    if(CCTK_Equals(eos_type, "nuclear")) {
        std::fill(&densxn[0], &densxn[siz], 0.0);
        std::fill(&densxp[0], &densxp[siz], 0.0);
    }
    std::fill(&scon[0], &scon[3*siz], 0.0);
    if(!CCTK_Equals(eos_type, "barotropic")) {
        std::fill(&tau[0], &tau[siz], 0.0);
    }

    // Initialize extra variables
    std::fill(&press[0], &press[siz], 0.0);
    std::fill(&w_lorentz[0], &w_lorentz[siz], 1.0);
    std::fill(&csound[0], &csound[siz], 0.0);
    if(!CCTK_Equals(eos_type, "ultrarelativistic")) {
        std::fill(&c2a_nbiter[0], &c2a_nbiter[siz], 0.0);
        std::fill(&volform[0], &volform[siz], 0.0);
    }
    if(CCTK_Equals(eos_type, "nuclear")) {
        std::fill(&entropy[0], &entropy[siz], 0.0);
        std::fill(&temperature[0], &temperature[siz], 0.0);
    }
#endif

    thc::Printer::start(
            "[INFO|THC|PrimToAll]: ",
            "[WARN|THC|PrimToAll]: ",
            "[ERR|THC|PrimToAll]: ");
    prim_to_all(cctkGH);
    thc::Printer::stop();
}
