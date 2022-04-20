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


#ifndef UTILS_NEWTON_HH
#define UTILS_NEWTON_HH

#include <limits>

#include <cctk.h>

#include <utils_functional.hh>

namespace utils {

//! Newton-Raphson root-finder namespace
namespace newton {

//! exit-status of the Newton-Raphson root finder
enum status {
    //! the final residual is within the required accuracy
    ok,
    //! the root finder exceeded the maximum number of iterations
    out_of_iterations,
    //! the root finder is stuck somewhere
    stuck,
    //! the root finder diverged away from the search interval
    diverged,
    //! the root finder exited in an unexpected way
    bug
};

//! Newton-Raphson root finder
newton::status rootfinder(
        //! [in] target function
        unary_function<CCTK_REAL, CCTK_REAL> const & f,
        //! [in] derivative of the target function
        unary_function<CCTK_REAL, CCTK_REAL> const & df,
        //! [in] minimum step
        CCTK_REAL min_step,
        //! [in] absolute value of the minimum residual
        CCTK_REAL abs_min_residual,
        //! [in] maximum number of iterations
        int max_iterations,
        //! [in,out] initial guess point in input, final point in output
        CCTK_REAL & point,
        //! [out] final step size (with sign)
        CCTK_REAL & step,
        //! [out] final residual (with sign)
        CCTK_REAL & residual,
        //! [out] number of iterations
        int & iteration,
        //! [in] minimum value for the point
        CCTK_REAL point_min = - std::numeric_limits<double>::infinity(),
        //! [in] maximum value for the point
        CCTK_REAL point_max = std::numeric_limits<double>::infinity()
        );

} // namespace newton
} // namespace utils

#endif

