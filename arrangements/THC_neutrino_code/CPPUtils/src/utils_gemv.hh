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


#ifndef UTILS_GEMV_HH
#define UTILS_GEMV_HH

#include <utils_blas_small.hh>
#include <utils_dot.hh>

#ifndef CBLAS_ENUM_DEFINED_H
   #define CBLAS_ENUM_DEFINED_H
   enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102 };
   enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113,
                         AtlasConj=114};
#endif

extern "C" void cblas_sgemv(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 const float *X, const int incX, const float beta,
                 float *Y, const int incY);

extern "C" void cblas_dgemv(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY);

namespace utils {
namespace details {

//! BLAS *gemv for arrays of known size
/*!
 *  \tparam T type of the arrays
 *  \tparam tranpose treat the matrix as being transposed
 *  \tparam M number of rows of the matrix
 *  \tparam N number of colums of the matrix
 *  \tparam small optimize for small vectors
 */
template<typename T, bool transpose, int M, int N, bool small>
class _gemv;

//! Specialization of gemv for large matrices
template<bool transpose, int M, int N>
class _gemv<float, transpose, M, N, false> {
    public:
        //! BLAS dgemv
        static inline void eval(
                float alpha,            //!< [in] alpha constant
                float const * A,        //!< [in] A matrix
                int lda,                //!< [in] length of the first dim of A
                float const * x,        //!< [in] x vector
                int incx,               //!< [in] increment of the x vector
                float beta,             //!< [in] beta constant
                float * y,              //!< [in, out] y vector
                int incy                //!< [in] increment of the y vector
                ) {
            cblas_sgemv(CblasRowMajor, (transpose ? CblasTrans : CblasNoTrans),
                M, N, alpha, A, lda, x, incx, beta, y, incy);
        }
        //! Simplified interface for vector-matrix multiplication
        static inline void simple(
                float const * A,        //!< [in] A matrix
                float const * x,        //!< [in] x vector
                float * y               //!< [in, out] y vector
                ) {
            cblas_sgemv(CblasRowMajor, (transpose ? CblasTrans : CblasNoTrans),
                M, N, 1, A, N, x, 1, 0, y, 1);
        }
};

//! Specialization of gemv for large matrices
template<bool transpose, int M, int N>
class _gemv<double, transpose, M, N, false> {
    public:
        //! BLAS dgemv
        static inline void eval(
                double alpha,           //!< [in] alpha constant
                double const * A,       //!< [in] A matrix
                int lda,                //!< [in] length of the first dim of A
                double const * x,       //!< [in] x vector
                int incx,               //!< [in] increment of the x vector
                double beta,            //!< [in] beta constant
                double * y,             //!< [in, out] y vector
                int incy                //!< [in] increment of the y vector
                ) {
            cblas_dgemv(CblasRowMajor, (transpose ? CblasTrans : CblasNoTrans),
                M, N, alpha, A, lda, x, incx, beta, y, incy);
        }
        //! Simplified interface for vector-matrix multiplication
        static inline void simple(
                double const * A,       //!< [in] A matrix
                double const * x,       //!< [in] x vector
                double * y              //!< [in, out] y vector
                ) {
            cblas_dgemv(CblasRowMajor, (transpose ? CblasTrans : CblasNoTrans),
                M, N, 1, A, N, x, 1, 0, y, 1);
        }
};

//! Specialization of gemv for the case of non-transposed *small* matrices
/*!
 *  \tparam T type of the arrays
 *  \tparam M number of rows of the matrix
 *  \tparam N number of colums of the matrix
 */
template<typename T, int M, int N>
class _gemv<T, false, M, N, true> {
    public:
        //! BLAS Tgemv
        static inline void eval(
                T alpha,                //!< [in] alpha constant
                T const * A,            //!< [in] A matrix
                int lda,                //!< [in] length of the first dim of A
                T const * x,            //!< [in] x vector
                int incx,               //!< [in] increment of the x vector
                T beta,                 //!< [in] beta constant
                T * y,                  //!< [in, out] y vector
                int incy                //!< [in] increment of the y vector
                ) {
            for(int i = 0; i < M; ++i) {
                y[i*incy] = alpha * dot<T, N>::eval(&A[i*lda], 1, x, incx)
                    + beta * y[i*incy];
            }
        }
        //! Simplified interface for vector-matrix multiplication
        static inline void simple(
                T const * A,            //!< [in] A matrix
                T const * x,            //!< [in] x vector
                T * y                   //!< [in, out] y vector
                ) {
            for(int i = 0; i < M; ++i) {
                y[i] = dot<T, N>::simple(&A[i*N], x);
            }
        }
};

//! specialization of gemv for the case of transposed *small* matrices
/*!
 *  \tparam T type of the arrays
 *  \tparam M number of rows of the matrix
 *  \tparam N number of colums of the matrix
 */
template<typename T, int M, int N>
class _gemv<T, true, M, N, true> {
    public:
        //! BLAS Tgemv
        static inline void eval(
                T alpha,                //!< [in] alpha constant
                T const * A,            //!< [in] A matrix
                int lda,                //!< [in] length of the first dim of A
                T const * x,            //!< [in] x vector
                int incx,               //!< [in] increment of the x vector
                T beta,                 //!< [in] beta constant
                T * y,                  //!< [in, out] y vector
                int incy                //!< [in] increment of the y vector
                ) {
            for(int i = 0; i < N; ++i) {
                y[i*incy] = alpha * dot<T, M>::eval(&A[i], lda, x, incx) +
                    beta * y[i*incy];
            }
        }
        //! Simplified interface for vector-matrix multiplication
        static inline void simple(
                T const * A,            //!< [in] A matrix
                T const * x,            //!< [in] x vector
                T * y                   //!< [out] y vector
                ) {
            for(int i = 0; i < N; ++i) {
                y[i] = dot<T, M>::eval(&A[i], N, x, 1);
            }
        }
};

} // namespace details

//! BLAS *gemv for arrays of known size
/*!
 *  \tparam T type of the arrays
 *  \tparam tranpose treat the matrix as being transposed
 *  \tparam M number of rows of the matrix
 *  \tparam N number of colums of the matrix
 */
template<typename T, bool transpose, int M, int N>
class gemv {
    public:
        //! BLAS Tgemv
        static inline void eval(
                T alpha,                //!< [in] alpha constant
                T const * A,            //!< [in] A matrix
                int lda,                //!< [in] length of the first dim of A
                T const * x,            //!< [in] x vector
                int incx,               //!< [in] increment of the x vector
                T beta,                 //!< [in] beta constant
                T * y,                  //!< [in, out] y vector
                int incy                //!< [in] increment of the y vector
                ) {
            details::_gemv<T, transpose, M, N, details::blas_small<M>::value>::
                eval(alpha, A, lda, x, incx, beta, y, incy);
        }
        //! Simplified interface for vector-matrix multiplication
        static inline void simple(
                T const * A,            //!< [in] A matrix
                T const * x,            //!< [in] x vector
                T * y                   //!< [in, out] y vector
                ) {
            details::_gemv<T, transpose, M, N, details::blas_small<M>::value>
                ::simple(A,x,y);
        }
};

} // namespace utils

#endif
