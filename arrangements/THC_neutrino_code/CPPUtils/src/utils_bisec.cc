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

#include <utils_bisec.hh>

namespace utils {
namespace bisec {

status rootfinder(
        unary_function<CCTK_REAL, CCTK_REAL> const & f,
        CCTK_REAL point_min,
        CCTK_REAL point_max,
        int iterations,
        CCTK_REAL & point) {
    point = 0.5*(point_min + point_max);

    CCTK_REAL fa = f(point_min);
    CCTK_REAL fb = f(point_max);
    CCTK_REAL fx = f(point);

    CCTK_REAL rf = (point_max - point_min) * std::pow(0.5, iterations);

    if(fa*fb >= 0) {
        if(std::abs(fa) < rf) {
            point = point_min;
            return bisec::ok;
        }
        else if(std::abs(fb) < rf) {
            point = point_max;
            return bisec::ok;
        }
        else {
            return bisec::failed;
        }
    }

    CCTK_INT k = 0;
    do {
        ++k;

        if(fa*fx < 0) {
            point_max = point;
            fb = fx;
        }
        else if(fa*fx > 0) {
            point_min = point;
            fa = fx;
        }
        else {
            if(fa == 0) {
                point = point_min;
            }
            break;
        }

        point = 0.5*(point_min + point_max);
        fx = f(point);
    } while(k <= iterations);

    return bisec::ok;
}

} // namespace bisec
} // namespace utils
