//  CPPUtils: C++ utilities for Cactus
//  Copyright (C) 2014, David Radice <dradice@caltech.edu>
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


#ifndef UTILS_VALENCIA_HH
#define UTILS_VALENCIA_HH

#include <cctk.h>

namespace utils {
namespace valencia {

//! Computes the four-velocity
void uvel(
        CCTK_REAL alp,
        CCTK_REAL betax,
        CCTK_REAL betay,
        CCTK_REAL betaz,
        CCTK_REAL w_lorentz,
        CCTK_REAL velx,
        CCTK_REAL vely,
        CCTK_REAL velz,
        CCTK_REAL u[4]);


} // namespace valencia
} // namespace utils

#endif
