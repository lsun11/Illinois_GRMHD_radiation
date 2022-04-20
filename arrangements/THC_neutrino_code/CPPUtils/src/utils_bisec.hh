//  CPPUtils: C++ utilities for Cactus
//  Copyright (C) 2011, David Radice <david.radice@aei.mpg.de>
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


#ifndef UTILS_BISEC_HH
#define UTILS_BISEC_HH

#include <cctk.h>

#include <utils_functional.hh>

namespace utils {

//! bisection routine namespace
namespace bisec {

enum status {
    //! the final residual is within the required accuracy
    ok,
    //! \brief the bisection has failed (because the target function
    //! doesn't change sign at the boundary of the interval)
    failed
};

//! Finds the root of a function in a given interval
status rootfinder(
        //! [in] target function
        unary_function<CCTK_REAL, CCTK_REAL> const & f,
        //! [in] left side of the target interval
        CCTK_REAL point_min,
        //! [in] right side of the target interval
        CCTK_REAL point_max,
        //! [in] number of iterations
        int iterations,
        //! [out] approximated root
        CCTK_REAL & point);

} // namespace bisec

} // namespace utils

#endif
