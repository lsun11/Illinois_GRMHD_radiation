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


#include "utils_metric.hh"

#include "utils_gemv.hh"
#include "utils_pow.hh"

#define SQ(X) ((X)*(X))

namespace utils {
namespace metric {

void space(
        CCTK_REAL const gxx,
        CCTK_REAL const gxy,
        CCTK_REAL const gxz,
        CCTK_REAL const gyy,
        CCTK_REAL const gyz,
        CCTK_REAL const gzz,
        CCTK_REAL gamma[9]) {
    gamma[0] = gxx;
    gamma[1] = gxy;
    gamma[2] = gxz;
    gamma[3] = gxy;
    gamma[4] = gyy;
    gamma[5] = gyz;
    gamma[6] = gxz;
    gamma[7] = gyz;
    gamma[8] = gzz;
}

CCTK_REAL spatial_det(
        CCTK_REAL const gxx,
        CCTK_REAL const gxy,
        CCTK_REAL const gxz,
        CCTK_REAL const gyy,
        CCTK_REAL const gyz,
        CCTK_REAL const gzz) {
    return - SQ(gxz)*gyy + 2*gxy*gxz*gyz - gxx*SQ(gyz) - SQ(gxy)*gzz +
        gxx*gyy*gzz;
}
CCTK_REAL spatial_det(
        CCTK_REAL const gamma[9]) {
    return spatial_det(gamma[0], gamma[1], gamma[2],
            gamma[4], gamma[5], gamma[8]);
}

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
        CCTK_REAL * uzz) {
    *uxx = (-SQ(gyz) + gyy*gzz)/det;
    *uxy = (gxz*gyz  - gxy*gzz)/det;
    *uyy = (-SQ(gxz) + gxx*gzz)/det;
    *uxz = (-gxz*gyy + gxy*gyz)/det;
    *uyz = (gxy*gxz  - gxx*gyz)/det;
    *uzz = (-SQ(gxy) + gxx*gyy)/det;
}
void spatial_inv(
        CCTK_REAL const det,
        CCTK_REAL const gamma[9],
        CCTK_REAL ugamma[9]) {
    CCTK_REAL uxx, uxy, uxz, uyy, uyz, uzz;
    spatial_inv(det, gamma[0], gamma[1], gamma[2],
            gamma[4], gamma[5], gamma[8], &uxx, &uxy,
            &uxz, &uyy, &uyz, &uzz);
    space(uxx, uxy, uxz, uyy, uyz, uzz, ugamma);
}

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
        CCTK_REAL g[16]) {
    g[5]  = gxx;
    g[6]  = gxy;
    g[7]  = gxz;
    g[9]  = gxy;
    g[10] = gyy;
    g[11] = gyz;
    g[13] = gxz;
    g[14] = gyz;
    g[15] = gzz;

    CCTK_REAL betaup[3] = {betax, betay, betaz};
    CCTK_REAL betadw[3] = {0, 0, 0};
    gemv<CCTK_REAL, false, 3, 3>::eval(1.0, &g[5], 4, &betaup[0],
            1, 0.0, &betadw[0], 1);

    g[0] = - pow<2>(alp) + betadw[0]*betaup[0] + betadw[1]*betaup[1] +
            betadw[2]*betaup[2];

    g[1] = betadw[0];
    g[2] = betadw[1];
    g[3] = betadw[2];

    g[4]  = betadw[0];
    g[8]  = betadw[1];
    g[12] = betadw[2];
}

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
        CCTK_REAL u[16]) {
    u[0] = - 1.0/pow<2>(alp);

    CCTK_REAL const det = spatial_det(gxx, gxy, gxz, gyy, gyz, gzz);

    CCTK_REAL uxx, uxy, uxz, uyy, uyz, uzz;
    spatial_inv(det, gxx, gxy, gxz, gyy, gyz, gzz,
            &uxx, &uxy, &uxz, &uyy, &uyz, &uzz);
    u[5]  = uxx + betax * betax * u[0];
    u[6]  = uxy + betax * betay * u[0];
    u[7]  = uxz + betax * betaz * u[0];
    u[9]  = uxy + betax * betay * u[0];
    u[10] = uyy + betay * betay * u[0];
    u[11] = uyz + betay * betaz * u[0];
    u[13] = uxz + betax * betaz * u[0];
    u[14] = uyz + betay * betaz * u[0];
    u[15] = uzz + betaz * betaz * u[0];

    u[1] = betax * (-u[0]);
    u[2] = betay * (-u[0]);
    u[3] = betaz * (-u[0]);

    u[4]  = u[1];
    u[8]  = u[2];
    u[12] = u[3];
}

} // namespace metric
} // namespace utils
