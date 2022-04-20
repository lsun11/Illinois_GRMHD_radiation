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


#ifndef THC_SPH_GRID_H
#define THC_SPH_GRID_H

#include <stdbool.h>

#include "cctk.h"

#define THC_SPH_INDEX(grid, irad, iray)                                      \
    ((iray)*(grid->nrad) + (irad))

#ifdef __cplusplus
extern "C" {
#endif

//! Structure representing a spherical grid
typedef struct {
    //! Center of the spherical grid
    CCTK_REAL origin[3];
    //! Maximum radial extent of the spherical grid
    CCTK_REAL rmax;

    //! Number of radial points on the spherical grid
    int nrad;
    //! Number of angular points (both theta and phi) on the spherical grid
    int nray;
    //! Number of angular points in the inclination direction
    int ntheta;
    //! Number of angular points in the azimutal direction
    int nphi;
    //! If true the grid also includes the poles and the points at phi = 2 pi
    bool endpoint;

    //! Grid spacing in the radial direction
    CCTK_REAL dr;
    //! Grid spacing in the inclination direction
    CCTK_REAL dtheta;
    //! Grid spacing in the azimutal direction
    CCTK_REAL dphi;
} SphericalGrid;

//! Initialize a spherical grid
void thc_sph_grid_init(
        //! [out] Newly allocated SphericalGrid
        SphericalGrid ** grid,
        //! [in] Origin of the spherical grid to create
        CCTK_REAL const origin[3],
        //! [in] Maximum radius of the spherical grid to create
        CCTK_REAL const rmax,
        //! [in] Number of radial points in the spherical grid to create
        int const nrad,
        //! [in] Number of angular points in the inclination direction
        int const ntheta,
        //! [in] Number of points in the azimutal direction
        int const nphi,
        //! [in] \brief If true the grid will include the points at the poles
        //! and at phi = 2 pi
        bool endpoint);
//! De-allocates a spherical grid
void thc_sph_grid_free(
        SphericalGrid * grid);

//! Gets the origin of the spherical grid
void thc_sph_grid_get_origin(
        //! [in] SphericalGrid object
        SphericalGrid const * grid,
        //! [out] Center of the spherical grid
        CCTK_REAL origin[3]);
//! Gets the maximum radius of the spherical grid
CCTK_REAL thc_sph_grid_get_rmax(SphericalGrid const * grid);

//! Gets the number of points in the radial direction of the spherical grid
int thc_sph_grid_get_nrad(SphericalGrid const * grid);
//! Gets the number of points in the inclination direction of the spherical grid
int thc_sph_grid_get_ntheta(SphericalGrid const * grid);
//! Gets the number of points in the azimutal direction of the spherical grid
int thc_sph_grid_get_nphi(SphericalGrid const * grid);

//! Gets the spacing in the radial direction
CCTK_REAL thc_sph_grid_get_dr(SphericalGrid const * grid);
//! Gets the spacing in the inclination direction
CCTK_REAL thc_sph_grid_get_dtheta(SphericalGrid const * grid);
//! Gets the spacing in the azimutal direction
CCTK_REAL thc_sph_grid_get_dphi(SphericalGrid const * grid);
//! Gets all of the spacings
void thc_sph_grid_get_delta(
        SphericalGrid const * grid,
        CCTK_REAL * dr,
        CCTK_REAL * dtheta,
        CCTK_REAL * dphi);

//! Convert an angular index from 2D (theta, phi) to 1D (ray)
int thc_sph_grid_get_iray(
        SphericalGrid const * grid,
        int const itheta,
        int const iphi);
//! Gets the theta index associated with a given ray index
int thc_sph_grid_get_itheta(
        SphericalGrid const * grid,
        int const iray);
//! Gets the phi index associated with a given ray index
int thc_sph_grid_get_iphi(
        SphericalGrid const * grid,
        int const iray);

//! Gets the radial coordinate at the given grid point
CCTK_REAL thc_sph_grid_get_r(
        SphericalGrid const * grid,
        int const irad);
//! Gets the inclination of the given ray
CCTK_REAL thc_sph_grid_get_theta(
        SphericalGrid const * grid,
        int const iray);
//! Gets the azimuth of the given ray
CCTK_REAL thc_sph_grid_get_phi(
        SphericalGrid const * grid,
        int const iray);
//! Gets the coordinates of the given grid point
void thc_sph_grid_get_r_theta_phi(
        SphericalGrid const * grid,
        int const irad,
        int const iray,
        CCTK_REAL * r,
        CCTK_REAL * theta,
        CCTK_REAL * phi);
//! Gets the location of a given grid point in cartesian coordinates
void thc_sph_grid_get_x_y_z(
        SphericalGrid const * grid,
        int const irad,
        int const iray,
        CCTK_REAL * x,
        CCTK_REAL * y,
        CCTK_REAL * z);

//! Converts from spherical to cartesian coordinates
/*!
 *  Note: this is a generic coordinate transformation and does not take into
 *  account the possible shift of the spherical grid
 */
void thc_coord_sph_to_cart(
        CCTK_REAL const r,
        CCTK_REAL const theta,
        CCTK_REAL const phi,
        CCTK_REAL * x,
        CCTK_REAL * y,
        CCTK_REAL * z);
//! Converts from cartesian to spherical coordinates
/*!
 * Note: this is a generic coordinate transformation and does not take into
 * account the possible shift of the spherical grid
 */
void thc_coord_cart_to_sph(
        CCTK_REAL const x,
        CCTK_REAL const y,
        CCTK_REAL const z,
        CCTK_REAL * r,
        CCTK_REAL * theta,
        CCTK_REAL * phi);

//! \brief Computes the Jacobian of the coordinate transformation from
//! spherical to cartesian
/*!
 *  J_{ij} = d (x, y, z)_i / d (r, theta, phi)^j
 *
 *  Note: this is a generic coordinate transformation and does not take into
 *  account the possible shift of the spherical grid
 */
void thc_coord_sph_to_cart_J(
        CCTK_REAL const r,
        CCTK_REAL const theta,
        CCTK_REAL const phi,
        CCTK_REAL J[9]);
//! \brief Computes the Jacobian of the coordinate transformation from
//! cartesian to spherical
/*!
 *  J_{ij} = d (r, theta, phi)_i / d (x, y, z)^j
 *
 *  Note: this is a generic coordinate transformation and does not take into
 *  account the possible shift of the spherical grid
 */
void thc_coord_cart_to_sph_J(
        CCTK_REAL const r,
        CCTK_REAL const theta,
        CCTK_REAL const phi,
        CCTK_REAL J[9]);

#ifdef __cplusplus
}
#endif

#endif
