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

#include "Symmetry.h"

extern "C" void THC_SetSym(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    int const sym_scal[3] = { 1,  1,  1};
    int const sym_vecx[3] = {-1,  1,  1};
    int const sym_vecy[3] = { 1, -1,  1};
    int const sym_vecz[3] = { 1,  1, -1};

    if(CCTK_Equals(eos_type, "barotropic")) {
        SetCartSymVN(cctkGH, sym_scal, "THC_Core::dens");
    }
    if(CCTK_Equals(eos_type, "ideal")) {
        SetCartSymVN(cctkGH, sym_scal, "THC_Core::dens");
    }
    if(CCTK_Equals(eos_type, "nuclear")) {
        SetCartSymVN(cctkGH, sym_scal, "THC_Core::densxn");
        SetCartSymVN(cctkGH, sym_scal, "THC_Core::densxp");
    }

    SetCartSymVN(cctkGH, sym_vecx, "THC_Core::scon[0]");
    SetCartSymVN(cctkGH, sym_vecy, "THC_Core::scon[1]");
    SetCartSymVN(cctkGH, sym_vecz, "THC_Core::scon[2]");

    if(!CCTK_Equals(eos_type, "barotropic")) {
        SetCartSymVN(cctkGH, sym_scal, "THC_Core::tau");
    }

    if(!CCTK_Equals(eos_type, "ultrarelativistic")) {
        SetCartSymVN(cctkGH, sym_scal, "THC_Core::densgain");
    }
}
