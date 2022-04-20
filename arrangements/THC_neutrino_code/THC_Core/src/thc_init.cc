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

extern "C" int THC_Init(void) {
    CCTK_RegisterBanner("Templated Hydrodynamics Code: a modular "
            "(relativistic) CFD code for Cactus.");
    return 0;
}
