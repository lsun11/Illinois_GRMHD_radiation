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

#include <utils_regfalsi.hh>

namespace utils {
namespace regfalsi {

status rootfinder(
        unary_function<CCTK_REAL, CCTK_REAL> const & f,
        CCTK_REAL point_min,
        CCTK_REAL point_max,
        CCTK_REAL abs_min_residual,
        int max_iterations,
        CCTK_REAL & point,
        CCTK_REAL & residual,
        int & iterations) {
    // x, f(x) s.t. f(x) > 0
    CCTK_REAL xp;
    CCTK_REAL fp;

    // x, f(x) s.t. f(x) < 0
    CCTK_REAL xn;
    CCTK_REAL fn;

    // xs, f(xs) s.t. f(xs) f(x) < 0
    CCTK_REAL xs;
    CCTK_REAL fs;

    // current point
    CCTK_REAL & xx = point;
    CCTK_REAL & fx = residual;

    fx = f(point_min);
    if(fx < 0) {
        xn = point_min;
        fn = fx;

        xp = point_max;
        fp = f(xp);
    }
    else {
        xn = point_max;
        fn = f(xn);

        xp = point_min;
        fp = fx;
    }

    if(fp*fn >= 0) {
        if(std::abs(fn) < abs_min_residual) {
            xx = xn;
            return regfalsi::ok;
        }
        else if(std::abs(fp) < abs_min_residual) {
            xx = xp;
            return regfalsi::ok;
        }
        else {
            return regfalsi::failed;
        }
    }

    // Use the given point as initial guess
    if(xx < point_min || xx > point_max) {
        xx = 0.5*(point_min + point_max);
    }
    fx = f(xx);
    if(fx < 0) {
        xn = xx;
        fn = fx;
    }
    else {
        xp = xx;
        fp = fx;
    }

    if(fx*fn < 0) {
        xs = xn;
        fs = fn;
    }
    else {
        xs = xp;
        fs = fp;
    }

    iterations = 1;
    while(iterations < max_iterations && std::abs(fx) > abs_min_residual) {
        ++iterations;

        xx = xx - (xx - xs)/(fx - fs)*fx;
        fx = f(xx);

        if(fx < 0) {
            xn = xx;
            fn = fx;
        }
        else {
            xp = xx;
            fp = fx;
        }

        if(fx*fn < 0) {
            xs = xn;
            fs = fn;
        }
        else {
            xs = xp;
            fs = fp;
        }
    }

    if(std::abs(fx) < abs_min_residual) {
        return regfalsi::ok;
    }

    return regfalsi::out_of_iterations;
}

} // namespace regfalsi
} // namespace utils
