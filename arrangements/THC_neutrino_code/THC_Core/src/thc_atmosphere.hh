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


#ifndef THC_ATMOSPHERE_HH
#define THC_ATMOSPHERE_HH

#include "cctk.h"

namespace thc {

class atmosphere {
    public:
        static CCTK_REAL rho;
        static CCTK_REAL eps;
        static CCTK_REAL Ye;

        static CCTK_REAL csound;
        static CCTK_REAL entropy;
        static CCTK_REAL press;
        static CCTK_REAL temp;

        static CCTK_REAL tol;
        static bool dynamic;
};

}

#endif
