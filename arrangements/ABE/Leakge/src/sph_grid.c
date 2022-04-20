//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2015, David Radice <dradice@caltech.edu>
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


#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "cctk_Arguments.h"

#include "sph_grid.h"

#ifndef M_PI
#define M_PI   3.14159265358979323846264338327950288
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923132169163975144
#endif

void sph_grid_init(
        SphericalGrid ** out,
        CCTK_REAL const origin[3],
        CCTK_REAL const rmax,
        int const nrad,
        int const ntheta,
        int const nphi,
        bool endpoint) {
    assert(rmax > 0);
    assert(nrad > 1);
    assert(ntheta > 2);
    assert(nphi > 5);

    SphericalGrid * grid = malloc(sizeof(SphericalGrid));

    memcpy(grid->origin, origin, 3*sizeof(CCTK_REAL));
    grid->rmax = rmax;
    grid->nrad = nrad;
    grid->nray = ntheta*nphi;
    grid->ntheta = ntheta;
    grid->nphi = nphi;
    grid->endpoint = endpoint;
    grid->dr = rmax/(nrad - 1);
    grid->dtheta = M_PI/(ntheta - endpoint);
    grid->dphi = 2*M_PI/(nphi - endpoint);

    *out = grid;
}
void sph_grid_free(SphericalGrid * grid) {
    free(grid);
}

void sph_grid_get_origin(SphericalGrid const * grid, CCTK_REAL origin[3]) {
    memcpy(origin, grid->origin, 3*sizeof(CCTK_REAL));
}
CCTK_REAL sph_grid_get_rmax(SphericalGrid const * grid) {
    return grid->rmax;
}

int sph_grid_get_nrad(SphericalGrid const * grid) {
    assert(grid);
    return grid->nrad;
}
int sph_grid_get_ntheta(SphericalGrid const * grid) {
    assert(grid);
    return grid->ntheta;
}
int sph_grid_get_nphi(SphericalGrid const * grid) {
    assert(grid);
    return grid->nphi;
}

CCTK_REAL sph_grid_get_dr(SphericalGrid const * grid) {
    assert(grid);
    return grid->dr;
}
CCTK_REAL sph_grid_get_dtheta(SphericalGrid const * grid) {
    assert(grid);
    return grid->dtheta;
}
CCTK_REAL sph_grid_get_dphi(SphericalGrid const * grid) {
    assert(grid);
    return grid->dphi;
}
void sph_grid_get_delta(
        SphericalGrid const * grid,
        CCTK_REAL * dr,
        CCTK_REAL * dtheta,
        CCTK_REAL * dphi) {
    *dr = sph_grid_get_dr(grid);
    *dtheta = sph_grid_get_dtheta(grid);
    *dphi = sph_grid_get_dphi(grid);
}

int sph_grid_get_iray(
        SphericalGrid const * grid,
        int const itheta,
        int const iphi) {
    return itheta + iphi * grid->ntheta;
}
int sph_grid_get_itheta(
        SphericalGrid const * grid,
        int const iray) {
    return iray % grid->ntheta;
}
int sph_grid_get_iphi(
        SphericalGrid const * grid,
        int const iray) {
    return iray / grid->ntheta;
}

CCTK_REAL sph_grid_get_r(SphericalGrid const * grid, int const i) {
    assert(grid);
    assert(i >= 0 && i < grid->nrad);
    return grid->dr*i;
}
CCTK_REAL sph_grid_get_theta(SphericalGrid const * grid, int const i) {
    assert(grid);
    assert(i >= 0 && i < grid->nray);
    return grid->dtheta * (1.0 + 0.5*(grid->endpoint ? 0 : 1)) *
        sph_grid_get_itheta(grid, i);
}
CCTK_REAL sph_grid_get_phi(SphericalGrid const * grid, int const i) {
    assert(grid);
    assert(i >= 0 && i < grid->nray);
    return grid->dphi * sph_grid_get_iphi(grid, i);
}
void sph_grid_get_r_theta_phi(
        SphericalGrid const * grid,
        int const irad,
        int const iray,
        CCTK_REAL * r,
        CCTK_REAL * theta,
        CCTK_REAL * phi) {
    *r     = sph_grid_get_r(grid, irad);
    *theta = sph_grid_get_theta(grid, iray);
    *phi   = sph_grid_get_phi(grid, iray);
}
void sph_grid_get_x_y_z(
        SphericalGrid const * grid,
        int const irad,
        int const iray,
        CCTK_REAL * x,
        CCTK_REAL * y,
        CCTK_REAL * z) {
    CCTK_REAL r, theta, phi;
    sph_grid_get_r_theta_phi(grid, irad, iray, &r, &theta, &phi);
    coord_sph_to_cart(r, theta, phi, x, y, z);
    *x += grid->origin[0];
    *y += grid->origin[1];
    *z += grid->origin[2];
}

void coord_sph_to_cart(
        CCTK_REAL const r,
        CCTK_REAL const theta,
        CCTK_REAL const phi,
        CCTK_REAL * x,
        CCTK_REAL * y,
        CCTK_REAL * z) {
    assert(r >= 0);
    assert(theta >= 0 && theta <= M_PI);
    assert(phi >= 0 && phi <= 2*M_PI);
    *x = r*sin(theta)*cos(phi);
    *y = r*sin(theta)*sin(phi);
    *z = r*cos(theta);
}
void coord_cart_to_sph(
        CCTK_REAL const x,
        CCTK_REAL const y,
        CCTK_REAL const z,
        CCTK_REAL * r,
        CCTK_REAL * theta,
        CCTK_REAL * phi) {
    *r = sqrt(x*x + y*y + z*z);
    if(*r > DBL_EPSILON) {
        *theta = acos(z/(*r));
        *phi = atan2(y, x);
        if(*phi < 0) {
            *phi += 2*M_PI;
        }
    }
    else {
        *theta = M_PI_2;
        *phi = 0;
    }
}

void coord_sph_to_cart_J(
        CCTK_REAL const r,
        CCTK_REAL const theta,
        CCTK_REAL const phi,
        CCTK_REAL * J) {
    CCTK_REAL const sin_theta = sin(theta);
    CCTK_REAL const cos_theta = cos(theta);
    CCTK_REAL const sin_phi   = sin(phi);
    CCTK_REAL const cos_phi   = cos(phi);

    J[0] =       sin_theta * cos_phi;
    J[1] =   r * cos_theta * cos_phi;
    J[2] = - r * sin_theta * sin_phi;
    J[3] =       sin_theta * sin_phi;
    J[4] =   r * cos_theta * sin_phi;
    J[5] =   r * sin_theta * cos_phi;
    J[6] =       cos_theta;
    J[7] = - r * sin_theta;
    J[8] = 0;
}
void coord_cart_to_sph_J(
        CCTK_REAL const r,
        CCTK_REAL const theta,
        CCTK_REAL const phi,
        CCTK_REAL * J) {
    CCTK_REAL const ir        = 1.0/r;
    CCTK_REAL const sin_theta = sin(theta);
    CCTK_REAL const cos_theta = cos(theta);
    CCTK_REAL const csc_theta = 1.0/sin_theta;
    CCTK_REAL const sin_phi   = sin(phi);
    CCTK_REAL const cos_phi   = cos(phi);

    J[0] =        sin_theta * cos_phi;
    J[1] =        sin_theta * sin_phi;
    J[2] =        cos_theta;
    J[3] =   ir * cos_theta * cos_phi;
    J[4] =   ir * cos_theta * sin_phi;
    J[5] = - ir * sin_theta;
    J[6] = - ir * csc_theta * sin_phi;
    J[7] =   ir * csc_theta * cos_phi;
    J[8] = 0;
}


