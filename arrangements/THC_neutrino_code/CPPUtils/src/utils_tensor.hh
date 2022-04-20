//  CPPUtils: C++ utilities for Cactus
//  Copyright (C) 2015, 2016, David Radice <dradice@caltech.edu>
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


#ifndef UTILS_TENSOR_HH
#define UTILS_TENSOR_HH

#include <algorithm>

#include <cctk.h>

#include <utils_metric.hh>
#include <utils_pow.hh>
#include <utils_valencia.hh>

namespace utils {
//! Classes for dealing with tensors and tensor fields
namespace tensor {

//! Kronecher delta symbol
inline CCTK_REAL delta(int a, int b) {
    return a == b ? 1.0 : 0.0;
}

///////////////////////////////////////////////////////////////////////////////
//! Generic tensor
///////////////////////////////////////////////////////////////////////////////
/*!
 *  \tparam T data type
 *  \tparam ndim number of dimensions
 *  \tparam rank tensor rank
 */
template<typename T, int ndim, int rank>
class generic;

template<typename T, int ndim_>
class generic<T, ndim_, 1> {
    public:
        //! Data type
        typedef T data_t;
        //! Rank of the tensor
        static int const rank = 1;
        //! Number of dimensions
        static int const ndim = ndim_;
        //! Number of degrees of freedom
        static int const ndof = ndim;

        //! Computes the compressed index associated with the given indices
        /*!
         *  The data is stored in a row major format
         */
        static inline int multiindex(int const a) {
            return a;
        }

        //! Access the components of a tensor using a compressed index
        /*!
         *  The data is stored in a row major format
         */
        inline T & operator[](int const a) {
#ifdef CPPUTILS_DEBUG
            assert(a >= 0 && a < ndof);
#endif
            return m_data[a];
        }
        //! Access the components of a tensor using a compressed index
        /*!
         *  The data is stored in a row major format
         */
        inline T const & operator[](int const a) const {
#ifdef CPPUTILS_DEBUG
            assert(a >= 0 && a <ndof);
#endif
            return m_data[a];
        }

        //! Access the components of a tensor using the natural indices
        inline T & operator()(int const a) {
            return this->operator[](a);
        }
        //! Access the components of a tensor using the natural indices
        inline T const & operator()(int const a) const {
            return this->operator[](a);
        }
    private:
        T m_data[ndim];
};

template<typename T, int ndim_>
class generic<T, ndim_, 2> {
    public:
        //! Data type
        typedef T data_t;
        //! Rank of the tensor
        static int const rank = 2;
        //! Number of dimensions
        static int const ndim = ndim_;
        //! Number of degrees of freedom
        static int const ndof = ndim*ndim;

        //! Computes the compressed index associated with the given indices
        /*!
         *  The data is stored in a row major format
         */
        static inline int multiindex(int const a, int const b) {
            return a*ndim + b;
        }

        //! Access the components of a tensor using a compressed index
        /*!
         *  The data is stored in a row major format
         */
        inline T & operator[](int const a) {
#ifdef CPPUTILS_DEBUG
            assert(a >= 0 && a < ndof);
#endif
            return m_data[a];
        }
        //! Access the components of a tensor using a compressed index
        /*!
         *  The data is stored in a row major format
         */
        inline T const & operator[](int const a) const {
#ifdef CPPUTILS_DEBUG
            assert(a >= 0 && a < ndof);
#endif
            return m_data[a];
        }

        //! Access the components of a tensor using the natural indices
        inline T & operator()(int const a, int const b) {
            return this->operator[](multiindex(a,b));
        }
        //! Access the components of a tensor using the natural indices
        inline T const & operator()(int const a, int const b) const {
            return this->operator[](multiindex(a,b));
        }
    private:
        T m_data[ndof];
};

template<typename T, int ndim_>
class generic<T, ndim_, 3> {
    public:
        //! Data type
        typedef T data_t;
        //! Rank of the tensor
        static int const rank = 3;
        //! Number of dimensions
        static int const ndim = ndim_;
        //! Number of degrees of freedom
        static int const ndof = ndim*ndim*ndim;

        //! Computes the compressed index associated with the given indices
        /*!
         *  The data is stored in a row major format
         */
        static inline int multiindex(int const a, int const b, int const c) {
            return a*ndim*ndim + b*ndim + c;
        }

        //! Access the components of a tensor using a compressed index
        /*!
         *  The data is stored in a row major format
         */
        inline T & operator[](int const a) {
#ifdef CPPUTILS_DEBUG
            assert(a >= 0 && a < ndof);
#endif
            return m_data[a];
        }
        //! Access the components of a tensor using a compressed index
        /*!
         *  The data is stored in a row major format
         */
        inline T const & operator[](int const a) const {
#ifdef CPPUTILS_DEBUG
            assert(a >= 0 && a < ndof);
#endif
            return m_data[a];
        }

        //! Access the components of a tensor using the natural indices
        inline T & operator()(int const a, int const b, int const c) {
            return this->operator[](multiindex(a,b,c));
        }
        //! Access the components of a tensor using the natural indices
        inline T const & operator()(int const a, int const b,
                int const c) const {
            return this->operator[](multiindex(a,b,c));
        }
    private:
        T m_data[ndof];
};

///////////////////////////////////////////////////////////////////////////////
//! Symmetric tensor with respect to the last 2 indices
///////////////////////////////////////////////////////////////////////////////
/*!
 *  \tparam T data type
 *  \tparam ndim number of dimensions
 *  \tparam rank tensor rank
 */
template<typename T, int ndim, int rank>
class symmetric2;

template<typename T, int ndim_>
class symmetric2<T, ndim_, 2> {
    public:
        //! Data type
        typedef T data_t;
        //! Rank of the tensor
        static int const rank = 2;
        //! Number of dimensions
        static int const ndim = ndim_;
        //! Number of degrees of freedom
        static int const ndof = (ndim*(ndim + 1))/2;

        //! Computes the compressed index associated with the given indices
        /*!
         *  The data is stored in a row major format
         *  Symmetric components g(i,j) with i > j are not stored so that, for
         *  instance a symmetric2<3, 2> tensor is stored as [Txx, Txy, Txz,
         *  Tyy, Tyz, Tzz].
         */
        static inline int multiindex(int const a, int const b) {
            if(b < a) {
                return multiindex(b, a);
            }
            int const offset = ndim*a - (a*(a - 1))/2;
            return offset + b - a;
        }

        //! Access the components of a tensor using a compressed index
        /*!
         *  The data is stored in a row major format.
         *  Symmetric components g(i,j) with i > j are not stored so that, for
         *  instance a symmetric2<3, 2> tensor is stored as [Txx, Txy, Txz,
         *  Tyy, Tyz, Tzz].
         */
        inline T & operator[](int const a) {
#ifdef CPPUTILS_DEBUG
            assert(a >= 0 && a < ndof);
#endif
            return m_data[a];
        }
        //! Access the components of a tensor using a compressed index
        /*!
         *  The data is stored in a row major format.
         *  Symmetric components g(i,j) with i > j are not stored so that, for
         *  instance a symmetric2<3, 2> tensor is stored as [Txx, Txy, Txz,
         *  Tyy, Tyz, Tzz].
         */
        inline T const & operator[](int const a) const {
#ifdef CPPUTILS_DEBUG
            assert(a >= 0 && a < ndof);
#endif
            return m_data[a];
        }

        //! Access the components of a tensor using the natural indices
        inline T & operator()(int const a, int const b) {
            return this->operator[](multiindex(a,b));
        }
        //! Access the components of a tensor using the natural indices
        inline T const & operator()(int const a, int const b) const {
            return this->operator[](multiindex(a,b));
        }
    private:
        T m_data[(ndim*(ndim + 1))/2];
};

template<typename T, int ndim_>
class symmetric2<T, ndim_, 3> {
    public:
        //! Data type
        typedef T data_t;
        //! Rank of the tensor
        static int const rank = 3;
        //! Number of dimensions
        static int const ndim = ndim_;
        //! Number of degrees of freedom
        static int const ndof = ndim*(ndim*(ndim + 1))/2;

        //! Computes the compressed index associated with the given indices
        /*!
         *  The data is stored in a row major format.
         *  Symmetric components g(i,j) with i > j are not stored so that, for
         *  instance a symmetric2<2, 3> tensor is stored as  [Txxx, Txxy, Txyy,
         *  Tyxx, Tyxy, Tyyy]
         */
        static inline int multiindex(int const a, int const b, int const c) {
            if(c < b) {
                return multiindex(a, c, b);
            }
            int const offset = ndim*b - (b*(b - 1))/2;
            return a*(ndim*(ndim + 1))/2 + offset + c - b;
        }

        //! Access the components of a tensor using a compressed index
        /*!
         *  The data is stored in a row major format.
         *  Symmetric components g(i,j) with i > j are not stored so that, for
         *  instance a symmetric2<2, 3> tensor is stored as  [Txxx, Txxy, Txyy,
         *  Tyxx, Tyxy, Tyyy]
         */
        inline T & operator[](int const a) {
#ifdef CPPUTILS_DEBUG
            assert(a >= 0 && a < ndof);
#endif
            return m_data[a];
        }
        //! Access the components of a tensor using a compressed index
        /*!
         *  The data is stored in a row major format.
         *  Symmetric components g(i,j) with i > j are not stored so that, for
         *  instance a symmetric2<2, 3> tensor is stored as  [Txxx, Txxy, Txyy,
         *  Tyxx, Tyxy, Tyyy]
         */
        inline T const & operator[](int const a) const {
#ifdef CPPUTILS_DEBUG
            assert(a >= 0 && a < ndof);
#endif
            return m_data[a];
        }

        //! Access the components of a tensor using the natural indices
        inline T & operator()(int const a, int const b, int const c) {
            return this->operator[](multiindex(a,b,c));
        }
        //! Access the components of a tensor using the natural indices
        inline T const & operator()(int const a, int const b, int const c) const {
            return this->operator[](multiindex(a,b,c));
        }
    private:
        T m_data[ndof];
};

///////////////////////////////////////////////////////////////////////////////
//! Metric tensor
///////////////////////////////////////////////////////////////////////////////
template<int ndim>
class metric;

template<int ndim>
class inv_metric;

//! Spatial metric
template<>
class metric<3>: public symmetric2<CCTK_REAL, 3, 2> {
    public:
        //! Computes the metric determinant at the given point
        inline CCTK_REAL det() const {
            return utils::metric::spatial_det(
                    this->operator[](0),
                    this->operator[](1),
                    this->operator[](2),
                    this->operator[](3),
                    this->operator[](4),
                    this->operator[](5));
        }
};

//! Inverse of the spatial metric
template<>
class inv_metric<3>: public symmetric2<CCTK_REAL, 3, 2> {
    public:
        //! Constructs the inverse metric from the metric
        inline void from_metric(
                CCTK_REAL const gxx,
                CCTK_REAL const gxy,
                CCTK_REAL const gxz,
                CCTK_REAL const gyy,
                CCTK_REAL const gyz,
                CCTK_REAL const gzz) {
            CCTK_REAL const det = utils::metric::spatial_det(
                    gxx, gxy, gxz, gyy, gyz, gzz);
            CCTK_REAL uxx, uxy, uxz, uyy, uyz, uzz;
            utils::metric::spatial_inv(det,
                    gxx, gxy, gxz, gyy, gyz, gzz,
                    &uxx, &uxy, &uxz, &uyy, &uyz, &uzz);
            this->operator[](0) = uxx;
            this->operator[](1) = uxy;
            this->operator[](2) = uxz;
            this->operator[](3) = uyy;
            this->operator[](4) = uyz;
            this->operator[](5) = uzz;
        }
        //! Constructs the inverse metric from the metric and the spatial det
        inline void from_metric(
                CCTK_REAL const gxx,
                CCTK_REAL const gxy,
                CCTK_REAL const gxz,
                CCTK_REAL const gyy,
                CCTK_REAL const gyz,
                CCTK_REAL const gzz,
                CCTK_REAL const det) {
            CCTK_REAL uxx, uxy, uxz, uyy, uyz, uzz;
            utils::metric::spatial_inv(det,
                    gxx, gxy, gxz, gyy, gyz, gzz,
                    &uxx, &uxy, &uxz, &uyy, &uyz, &uzz);
            this->operator[](0) = uxx;
            this->operator[](1) = uxy;
            this->operator[](2) = uxz;
            this->operator[](3) = uyy;
            this->operator[](4) = uyz;
            this->operator[](5) = uzz;
        }
        //! Constructs the inverse metric from the metric
        inline void from_metric(metric<3> const & g) {
            CCTK_REAL const det = g.det();
            CCTK_REAL uxx, uxy, uxz, uyy, uyz, uzz;
            utils::metric::spatial_inv(det,
                    g[0], g[1], g[2], g[3], g[4], g[5],
                    &uxx, &uxy, &uxz, &uyy, &uyz, &uzz);
            this->operator[](0) = uxx;
            this->operator[](1) = uxy;
            this->operator[](2) = uxz;
            this->operator[](3) = uyy;
            this->operator[](4) = uyz;
            this->operator[](5) = uzz;
        }
        //! Constructs the inverse metric from the metric and the spatial det
        inline void from_metric_det(metric<3> const & g, CCTK_REAL const det) {
            CCTK_REAL uxx, uxy, uxz, uyy, uyz, uzz;
            utils::metric::spatial_inv(det,
                    g[0], g[1], g[2], g[3], g[4], g[5],
                    &uxx, &uxy, &uxz, &uyy, &uyz, &uzz);
            this->operator[](0) = uxx;
            this->operator[](1) = uxy;
            this->operator[](2) = uxz;
            this->operator[](3) = uyy;
            this->operator[](4) = uyz;
            this->operator[](5) = uzz;
        }
};

//! Spacetime metric
template<>
class metric<4>: public symmetric2<CCTK_REAL, 4, 2> {
    public:
        //! Construct the spacetime metric from the ADM quantities
        inline void from_adm(
                CCTK_REAL const alp,
                CCTK_REAL const betax,
                CCTK_REAL const betay,
                CCTK_REAL const betaz,
                CCTK_REAL const gxx,
                CCTK_REAL const gxy,
                CCTK_REAL const gxz,
                CCTK_REAL const gyy,
                CCTK_REAL const gyz,
                CCTK_REAL const gzz) {
            CCTK_REAL g[16];
            utils::metric::spacetime(alp, betax, betay, betaz, gxx, gxy, gxz,
                    gyy, gyz, gzz, &g[0]);
            for(int a = 0; a < 4; ++a)
            for(int b = a; b < 4; ++b) {
                this->operator()(a, b) = g[4*a + b];
            }
        }
};

//! Spacetime inverse metric
template<>
class inv_metric<4>: public symmetric2<CCTK_REAL, 4, 2> {
    public:
        //! Construct the spacetime metric from the ADM quantities
        inline void from_adm(
                CCTK_REAL const alp,
                CCTK_REAL const betax,
                CCTK_REAL const betay,
                CCTK_REAL const betaz,
                CCTK_REAL const gxx,
                CCTK_REAL const gxy,
                CCTK_REAL const gxz,
                CCTK_REAL const gyy,
                CCTK_REAL const gyz,
                CCTK_REAL const gzz) {
            CCTK_REAL u[16];
            utils::metric::spacetime_upper(alp, betax, betay,
                    betaz, gxx, gxy, gxz, gyy, gyz, gzz, &u[0]);
            for(int a = 0; a < 4; ++a)
            for(int b = a; b < 4; ++b) {
                this->operator()(a, b) = u[4*a + b];
            }
        }
};

///////////////////////////////////////////////////////////////////////////////
// Special tensor fields
///////////////////////////////////////////////////////////////////////////////
//! Fluid four velocity (as a vector)
class fluid_velocity_field_const {
    public:
        //! Initialize the fluid 4 velocity
        fluid_velocity_field_const(
                CCTK_REAL const * alp,
                CCTK_REAL const * betax,
                CCTK_REAL const * betay,
                CCTK_REAL const * betaz,
                CCTK_REAL const * w_lorentz,
                CCTK_REAL const * velx,
                CCTK_REAL const * vely,
                CCTK_REAL const * velz) {
            m_data[0] = alp;
            m_data[1] = betax;
            m_data[2] = betay;
            m_data[3] = betaz;
            m_data[4] = w_lorentz;
            m_data[5] = velx;
            m_data[6] = vely;
            m_data[7] = velz;
#ifdef CPPUTILS_DEBUG
            for(int i = 0; i < 8; ++i) {
                assert(m_data[i]);
            }
#endif
        }

        //! Evaluate the fluid four velocity at a given location
        inline void get(ptrdiff_t const ijk, generic<CCTK_REAL,4,1> * u) const {
            utils::valencia::uvel(m_data[0][ijk], m_data[1][ijk],
                    m_data[2][ijk], m_data[3][ijk], m_data[4][ijk],
                    m_data[5][ijk], m_data[6][ijk], m_data[7][ijk],
                    &u->operator[](0));
        }
    private:
        CCTK_REAL const * m_data[8];
};

//! Class describing the geometry of the ADM slicing
class slicing_geometry_const {
    public:
        //! Initialize the slicing geometry from the ADM quantities
        slicing_geometry_const(
                CCTK_REAL const * alp,
                CCTK_REAL const * betax,
                CCTK_REAL const * betay,
                CCTK_REAL const * betaz,
                CCTK_REAL const * gxx,
                CCTK_REAL const * gxy,
                CCTK_REAL const * gxz,
                CCTK_REAL const * gyy,
                CCTK_REAL const * gyz,
                CCTK_REAL const * gzz,
                CCTK_REAL const * kxx,
                CCTK_REAL const * kxy,
                CCTK_REAL const * kxz,
                CCTK_REAL const * kyy,
                CCTK_REAL const * kyz,
                CCTK_REAL const * kzz,
                CCTK_REAL const * volform) {
            m_data[0]  = alp;
            m_data[1]  = betax;
            m_data[2]  = betay;
            m_data[3]  = betaz;
            m_data[4]  = gxx;
            m_data[5]  = gxy;
            m_data[6]  = gxz;
            m_data[7]  = gyy;
            m_data[8]  = gyz;
            m_data[9]  = gzz;
            m_data[10] = kxx;
            m_data[11] = kxy;
            m_data[12] = kxz;
            m_data[13] = kyy;
            m_data[14] = kyz;
            m_data[15] = kzz;
            m_data[16] = volform;
#ifdef CPPUTILS_DEBUG
            for(int i = 0; i < 16; ++i) {
                assert(m_data[i]);
            }
#endif
        }

        //! \brief Get the normal vector to the spacelike hyper-surface of the
        //! ADM decomposition
        inline void get_normal(ptrdiff_t const ijk,
                generic<CCTK_REAL, 4, 1> * n) const {
            n->operator[](0) = 1.0/m_data[0][ijk];
            n->operator[](1) = - m_data[1][ijk]/m_data[0][ijk];
            n->operator[](2) = - m_data[2][ijk]/m_data[0][ijk];
            n->operator[](3) = - m_data[3][ijk]/m_data[0][ijk];
        }

        //! Get a component of the spatial metric on the entire grid
        CCTK_REAL const * get_space_metric_comp(int const a, int const b) const {
            return m_data[4 + metric<3>::multiindex(a,b)];
        }
        //! Get a component of the extrinsic curvature on the entire grid
        CCTK_REAL const * get_extr_curv_comp(int const a, int const b) const {
            return m_data[10 + symmetric2<CCTK_REAL,3,2>::multiindex(a,b)];

        }

        //! Get the spatial metric at a given location
        inline void get_metric(ptrdiff_t const ijk, metric<3> * g) const {
            g->operator[](0) = m_data[4][ijk];
            g->operator[](1) = m_data[5][ijk];
            g->operator[](2) = m_data[6][ijk];
            g->operator[](3) = m_data[7][ijk];
            g->operator[](4) = m_data[8][ijk];
            g->operator[](5) = m_data[9][ijk];
        }
        //! Get inverse spatial metric at a given location
        inline void get_inv_metric(ptrdiff_t const ijk,
                inv_metric<3> * u) const {
            u->from_metric(m_data[4][ijk], m_data[5][ijk], m_data[6][ijk],
                    m_data[7][ijk], m_data[8][ijk], m_data[9][ijk],
                    m_data[16][ijk]*m_data[16][ijk]);
        }

        //! Get the spacetime metric at a given location
        inline void get_metric(ptrdiff_t const ijk, metric<4> * g) const {
            g->from_adm(m_data[0][ijk], m_data[1][ijk], m_data[2][ijk],
                    m_data[3][ijk], m_data[4][ijk], m_data[5][ijk],
                    m_data[6][ijk], m_data[7][ijk], m_data[8][ijk],
                    m_data[9][ijk]);
        }

        //! Get the spacetime metric at a given location
        inline void get_inv_metric(ptrdiff_t const ijk,
                inv_metric<4> * u) const {
            u->from_adm(m_data[0][ijk], m_data[1][ijk], m_data[2][ijk],
                    m_data[3][ijk], m_data[4][ijk], m_data[5][ijk],
                    m_data[6][ijk], m_data[7][ijk], m_data[8][ijk],
                    m_data[9][ijk]);
        }

        //! Get the extrinsic curvature at a given location
        inline void get_extr_curv(ptrdiff_t const ijk,
                symmetric2<CCTK_REAL, 3, 2> * k) const {
            k->operator[](0) = m_data[10][ijk];
            k->operator[](1) = m_data[11][ijk];
            k->operator[](2) = m_data[12][ijk];
            k->operator[](3) = m_data[13][ijk];
            k->operator[](4) = m_data[14][ijk];
            k->operator[](5) = m_data[15][ijk];
        }
    private:
        CCTK_REAL const * m_data[17];
};

} // namespace tensor
} // namespace utils

#endif
