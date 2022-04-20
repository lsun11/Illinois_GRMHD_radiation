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


#include <cmath>

#include <utils_newton.hh>

namespace utils {
namespace newton {

newton::status rootfinder(
        unary_function<CCTK_REAL, CCTK_REAL> const & f,
        unary_function<CCTK_REAL, CCTK_REAL> const & df,
        CCTK_REAL min_step,
        CCTK_REAL abs_min_residual,
        int max_iterations,
        CCTK_REAL & point,
        CCTK_REAL & step,
        CCTK_REAL & residual,
        int & iteration,
        CCTK_REAL point_min,
        CCTK_REAL point_max) {
    CCTK_REAL old_point;

    residual = f(point);
    iteration = 0;

    do {
        ++iteration;
        old_point = point;
        point = old_point - residual/df(old_point);
        residual = f(point);
        step = point - old_point;
        if(point < point_min || point > point_max) {
            if(std::abs(residual) < abs_min_residual) {
                return newton::ok;
            }
            else {
                return newton::diverged;
            }
        }
    } while(std::abs(residual) > abs_min_residual
                && std::abs(step) > min_step
                && iteration < max_iterations);

    if(std::abs(residual) < abs_min_residual) {
        return newton::ok;
    }
    else if(std::abs(step) < min_step) {
        return newton::stuck;
    }
    else if(iteration == max_iterations) {
        return newton::out_of_iterations;
    }

    return newton::bug;
}

} // namespace newton
} // namespace utils
