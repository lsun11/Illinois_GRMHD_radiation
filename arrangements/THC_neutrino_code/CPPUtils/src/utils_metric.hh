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


#ifndef UTILS_METRIC_HH
#define UTILS_METRIC_HH

#include <cctk.h>

namespace utils {
//! Misc utilities for manipulating metric tensors
/*!
 *  NOTE: these are low level routines, see utils_tensor for a higher-level API
 */
namespace metric {

//! Construct a spatial metric
void space(
        CCTK_REAL const gxx,
        CCTK_REAL const gxy,
        CCTK_REAL const gxz,
        CCTK_REAL const gyy,
        CCTK_REAL const gyz,
        CCTK_REAL const gzz,
        CCTK_REAL gamma[9]);

//! Computes the determinant of the metric
CCTK_REAL spatial_det(
        CCTK_REAL const gxx,
        CCTK_REAL const gxy,
        CCTK_REAL const gxz,
        CCTK_REAL const gyy,
        CCTK_REAL const gyz,
        CCTK_REAL const gzz);
//! Computes the determinant of the metric
CCTK_REAL spatial_det(
        CCTK_REAL const gamma[9]);

//! Computes the inverse of the spatial metric
void spatial_inv(
        CCTK_REAL const det,
        CCTK_REAL const gxx,
        CCTK_REAL const gxy,
        CCTK_REAL const gxz,
        CCTK_REAL const gyy,
        CCTK_REAL const gyz,
        CCTK_REAL const gzz,
        CCTK_REAL * uxx,
        CCTK_REAL * uxy,
        CCTK_REAL * uxz,
        CCTK_REAL * uyy,
        CCTK_REAL * uyz,
        CCTK_REAL * uzz);
//! Computes the inverse of the spatial metric
void spatial_inv(
        CCTK_REAL const det,
        CCTK_REAL const gamma[9],
        CCTK_REAL ugamma[9]);

//! Construct a spacetime metric
void spacetime(
        CCTK_REAL const alp,
        CCTK_REAL const betax,
        CCTK_REAL const betay,
        CCTK_REAL const betaz,
        CCTK_REAL const gxx,
        CCTK_REAL const gxy,
        CCTK_REAL const gxz,
        CCTK_REAL const gyy,
        CCTK_REAL const gyz,
        CCTK_REAL const gzz,
        CCTK_REAL g[16]);

//! Computes the inverse of the spacetime metric
void spacetime_upper(
        CCTK_REAL const alp,
        CCTK_REAL const betax,
        CCTK_REAL const betay,
        CCTK_REAL const betaz,
        CCTK_REAL const gxx,
        CCTK_REAL const gxy,
        CCTK_REAL const gxz,
        CCTK_REAL const gyy,
        CCTK_REAL const gyz,
        CCTK_REAL const gzz,
        CCTK_REAL u[16]);

} // namespace metric
} // namespace utils

#endif
