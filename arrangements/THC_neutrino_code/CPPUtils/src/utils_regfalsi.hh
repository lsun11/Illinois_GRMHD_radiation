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


#ifndef UTILS_REGULA_FALSI_HH
#define UTILS_REGULA_FALSI_HH

#include <cctk.h>

#include <utils_functional.hh>

namespace utils {

//! Regula falsi root-finder namespace
namespace regfalsi {

enum status {
    //! the final residual is within the required accuracy
    ok,
    //! the root-finder exceeded the maximum number of iterations
    out_of_iterations,
    //! \brief the root-finder has failed (because the target function
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
        //! [in] absolute value of the minimum residual
        CCTK_REAL abs_min_residual,
        //! [in] number of iterations
        int max_iterations,
        //! [in,out] initial guess in input, final result in output
        CCTK_REAL & point,
        //! [out] final residual
        CCTK_REAL & residual,
        //! [out] number of iterations
        int & iterations
        );

} // namespace regula falsi

} // namespace utils

#endif
