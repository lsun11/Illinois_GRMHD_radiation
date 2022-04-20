//  CPPUtils: C++ utilities for Cactus
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


#ifndef UTILS_QFORM_HH
#define UTILS_QFORM_HH

namespace utils {

//! Computes a generalized dot product, optimized for small arrays of known size
/*!
 *  \tparam T type of the arrays
 *  \tparam N size of the arrays
 *
 *  This computes
 *  \f[
        g_{ab} x^a y^b
    \f]
 */
template<typename T, int N>
class gdot {
    public:
        //! Evaluates the quadratic form
        static inline T eval(
                T const * g,            //!< [in] metric
                int const ldg,          //!< [in] leading dimension of g
                T const * x,            //!< [in] x vector
                int const incx,         //!< [in] increment of the x vector
                T const * y,            //!< [in] y vector
                int const incy          //!< [in] increment of the y vector
                ) {
            T sum(0);
            for(int a = 0; a < N; ++a) {
                sum += dot<T,N>::eval(&g[a*ldg], 1, y, incy)*x[a*incx];
            }
            return sum;
        }
        //! Evaluates the quadratic form: simple version
        static inline T simple(
                T const * g,            //!< [in] metric
                T const * x,            //!< [in] x vector
                T const * y             //!< [in] y vector
                ) {
            return eval(g, N, x, 1, y, 1);
        }
};

} // namespace

#endif
