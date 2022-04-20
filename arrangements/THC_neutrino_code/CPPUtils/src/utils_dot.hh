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


#include <utils_blas_small.hh>

#ifndef UTILS_DOT_HH
#define UTILS_DOT_HH

extern "C" float cblas_sdot(int const N, float const * X, int const incX,
        float const * Y, int const incY);

extern "C" double cblas_ddot(int const N, double const * X, int const incX,
        double const * Y, int const incY);

namespace utils {
namespace details {

//! BLAS dot optimized for arrays of known size
/*!
 *  \tparam T type of the arrays
 *  \tparam N size of the arrays
 *  \tparam small optimize for small arrays
 */
template<typename T, int N, bool small>
class _dot;

//! BLAS dot optimized for large arrays of known size
/*!
 *  \tparam N size of the arrays
 */
template<int N>
class _dot<float, N, false> {
    public:
        //! Evaluate the dot product between x and y
        /*!
         *  This has the same interface of the *dot from the BLAS specification
         */
        static inline float eval(
                float const * x,        //!< [in] first array
                int incx,               //!< [in] increment for the first array
                float const * y,        //!< [in] second array
                int incy) {             //!< [in] increment for the second array
            return cblas_sdot(N, x, incx, y, incy);
        }
        //! Simplified interface for the dot product
        static inline float simple(
                float const * x,        //!< [in] first array
                float const * y) {      //!< [in] second array
            return cblas_sdot(N, x, 1, y, 1);
        }
};

//! BLAS dot optimized for large arrays of known size
/*!
 *  \tparam N size of the arrays
 */
template<int N>
class _dot<double, N, false> {
    public:
        //! Evaluate the dot product between x and y
        /*!
         *  This has the same interface of the *dot from the BLAS specification
         */
        static inline double eval(
                double const * x,       //!< [in] first array
                int incx,               //!< [in] increment for the first array
                double const * y,       //!< [in] second array
                int incy) {             //!< [in] increment for the second array
            return cblas_ddot(N, x, incx, y, incy);
        }
        //! Simplified interface for the dot product
        static inline double simple(
                double const * x,       //!< [in] first array
                double const * y) {     //!< [in] second array
            return cblas_ddot(N, x, 1, y, 1);
        }
};

//! BLAS dot optimized for *small* arrays of known size
/*!
 *  \tparam T type of the arrays
 *  \tparam N size of the arrays
 */
template<typename T, int N>
class _dot<T, N, true> {
    public:
        //! Evaluate the dot product between x and y
        /*!
         *  This has the same interface of the *dot from the BLAS specification
         */
        static inline T eval(
                T const * x,          //!< [in] first array
                int incx,             //!< [in] increment for the first array
                T const * y,          //!< [in] second array
                int incy) {           //!< [in] increment for the second array
            return x[0]*y[0] + _dot<T, N-1, true>::eval(&x[incx], incx,
                    &y[incy], incy);
        }
        //! Simplified interface for the dot product
        static inline T simple(
                T const * x,          //!< [in] first array
                T const * y) {        //!< [in] second array
            return x[0]*y[0] + _dot<T, N-1, true>::simple(&x[1], &y[1]);
        }
};

//! Specialization of dot<T, N, true> for N == 1
template<typename T>
class _dot<T, 1, true> {
    public:
        static inline T eval(
                T const * x,          //!< [in] first array
                int,                  //!< [in] increment for the first array
                T const * y,          //!< [in] second array
                int) {                //!< [in] increment for the second array
            return x[0]*y[0];
        }
        static inline T simple(
                T const * x,          //!< [in] first array
                T const * y) {        //!< [in] second array
            return x[0]*y[0];
        }
};

} // namespace details

//! BLAS dot product interface
template<typename T, int N>
class dot {
    public:
        //! Evaluate the dot product between x and y
        /*!
         *  This has the same interface of the *dot from the BLAS specification
         */
        static inline T eval(
                T const * x,          //!< [in] first array
                int incx,             //!< [in] increment for the first array
                T const * y,          //!< [in] second array
                int incy) {           //!< [in] increment for the second array
            return details::_dot<T, N, details::blas_small<N>::value>::eval(x,
                    incx, y, incy);
        }
        //! Simplified interface for the dot product
        static inline T simple(
                T const * x,          //!< [in] first array
                T const * y) {        //!< [in] second array
            return details::_dot<T, N, details::blas_small<N>::value>::simple(
                    x, y);
        }
};

} // namespace utils

#endif
