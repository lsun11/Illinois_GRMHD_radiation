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


#include "utils_valencia.hh"

namespace utils {
namespace valencia {

void uvel(
        CCTK_REAL alp,
        CCTK_REAL betax,
        CCTK_REAL betay,
        CCTK_REAL betaz,
        CCTK_REAL w_lorentz,
        CCTK_REAL velx,
        CCTK_REAL vely,
        CCTK_REAL velz,
        CCTK_REAL u[4]) {
    CCTK_REAL const ialp = 1.0/alp;
    u[0] = w_lorentz*ialp;
    u[1] = w_lorentz*(velx - betax*ialp);
    u[2] = w_lorentz*(vely - betay*ialp);
    u[3] = w_lorentz*(velz - betaz*ialp);
}

} // namespace utils
} // namespace valencia
