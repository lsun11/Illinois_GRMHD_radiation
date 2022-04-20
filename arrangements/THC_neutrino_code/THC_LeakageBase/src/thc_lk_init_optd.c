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


// Almost all the grid functions can be found in THCCORE/THC_CORE/src/thc_sph_grid.h and THCCORE/THC_CORE/src/thc_sph_grid.c


#include <assert.h>
#include <math.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include <util_Table.h>

#include "thc_sph_grid.h"
#include "utils_macro.h"

#define length(X) (sizeof((X))/sizeof(*(X)))
#define MIN(X,Y) ((X)<(Y)?(X):(Y))
#define SQ(X) ((X)*(X))

static struct {
    SphericalGrid * grid;
    CCTK_REAL * metric[6];
    CCTK_REAL * grr;
    CCTK_REAL * xcor[3];
    CCTK_REAL * rho;
    CCTK_REAL * temp;
    CCTK_REAL * Y_e;
    CCTK_REAL * kappa[6];
    CCTK_REAL * tau[6];
} * __optd_data[2] = {NULL, NULL};

void THC_LK_CalcOpticalDepth(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

      /*
    if(verbose) {
        CCTK_INFO("THC_LK_CalcOpticalDepth");
    }
      */
      
    int const gsiz = nrad*ntheta*nphi;

    CCTK_REAL const * centers[2] = {&center_grid1[0], &center_grid2[0]};
    for(int gidx = 0; gidx < num_sph_grids; ++gidx) {
        __optd_data[gidx] = malloc(sizeof(*__optd_data[gidx]));

        /* Initialize the grid */
        thc_sph_grid_init(&__optd_data[gidx]->grid, centers[gidx], rmax,
                nrad, ntheta, nphi, true);
        CCTK_REAL const drad = thc_sph_grid_get_dr(__optd_data[gidx]->grid);
        for(int v = 0; v < 6; ++v) {
            __optd_data[gidx]->metric[v] = malloc(gsiz*sizeof(CCTK_REAL));
        }
        __optd_data[gidx]->grr = malloc(gsiz*sizeof(CCTK_REAL));
        for(int d = 0; d < 3; ++d) {
            __optd_data[gidx]->xcor[d] = malloc(gsiz*sizeof(CCTK_REAL));
        }
        __optd_data[gidx]->rho = malloc(gsiz*sizeof(CCTK_REAL));
        __optd_data[gidx]->temp = malloc(gsiz*sizeof(CCTK_REAL));
        __optd_data[gidx]->Y_e = malloc(gsiz*sizeof(CCTK_REAL));
        for(int v = 0; v < 6; ++v) {
            __optd_data[gidx]->kappa[v] = malloc(gsiz*sizeof(CCTK_REAL));
            __optd_data[gidx]->tau[v] = malloc(gsiz*sizeof(CCTK_REAL));
        }

        /* Compute the coordinates */
        for(int iray = 0; iray < ntheta*nphi; ++iray)
        for(int irad = 0; irad < nrad; ++irad) {
            int const idx = THC_SPH_INDEX(__optd_data[gidx]->grid, irad, iray);
            thc_sph_grid_get_x_y_z(__optd_data[gidx]->grid, irad, iray,
                    &__optd_data[gidx]->xcor[0][idx],
                    &__optd_data[gidx]->xcor[1][idx],
                    &__optd_data[gidx]->xcor[2][idx]);
        }

        /* Interpolate data onto the spherical grid */
        int const interp_handle = CCTK_InterpHandle(interpolator);
        assert(interp_handle >= 0);
        int const options_handle =
            Util_TableCreateFromString(interpolator_options);
        assert(options_handle >= 0);
        int const coords_handle = CCTK_CoordSystemHandle("cart3d");
        assert(coords_handle >= 0);

        void const ** interp_coords =
            (void const **)&__optd_data[gidx]->xcor[0];
        int const npoints = gsiz;
        CCTK_INT const input_array_indices[] = {
            CCTK_VarIndex("ADMBase::gxx"),
            CCTK_VarIndex("ADMBase::gxy"),
            CCTK_VarIndex("ADMBase::gxz"),
            CCTK_VarIndex("ADMBase::gyy"),
            CCTK_VarIndex("ADMBase::gyz"),
            CCTK_VarIndex("ADMBase::gzz"),
            CCTK_VarIndex("HydroBase::rho"),
            CCTK_VarIndex("HydroBase::temperature"),
            CCTK_VarIndex("HydroBase::Y_e")
        };
        int const ninputs = length(input_array_indices);

        CCTK_INT const output_array_types[] = {
            CCTK_VARIABLE_REAL,
            CCTK_VARIABLE_REAL,
            CCTK_VARIABLE_REAL,
            CCTK_VARIABLE_REAL,
            CCTK_VARIABLE_REAL,
            CCTK_VARIABLE_REAL,
            CCTK_VARIABLE_REAL,
            CCTK_VARIABLE_REAL,
            CCTK_VARIABLE_REAL,
        };
        assert(ninputs == length(output_array_types));

        void * output_arrays[] = {
            (void *)__optd_data[gidx]->metric[0],
            (void *)__optd_data[gidx]->metric[1],
            (void *)__optd_data[gidx]->metric[2],
            (void *)__optd_data[gidx]->metric[3],
            (void *)__optd_data[gidx]->metric[4],
            (void *)__optd_data[gidx]->metric[5],
            (void *)__optd_data[gidx]->rho,
            (void *)__optd_data[gidx]->temp,
            (void *)__optd_data[gidx]->Y_e,
        };
        assert(ninputs == length(output_arrays));

        int ierr = CCTK_InterpGridArrays(cctkGH, 3, interp_handle,
                options_handle, coords_handle, npoints, CCTK_VARIABLE_REAL,
                interp_coords, ninputs, input_array_indices,
                ninputs, output_array_types, output_arrays);
        assert(!ierr);
        Util_TableDestroy(options_handle);

        /* Compute opacities */
        for(int iray = 0; iray < ntheta*nphi; ++iray)
        for(int irad = 0; irad < nrad; ++irad) {
            int const idx = THC_SPH_INDEX(__optd_data[gidx]->grid, irad, iray);
            int ierr = NeutrinoOpacity(
                __optd_data[gidx]->rho[idx],
                __optd_data[gidx]->temp[idx],
                __optd_data[gidx]->Y_e[idx],
                &__optd_data[gidx]->kappa[0][idx],
                &__optd_data[gidx]->kappa[1][idx],
                &__optd_data[gidx]->kappa[2][idx],
                &__optd_data[gidx]->kappa[3][idx],
                &__optd_data[gidx]->kappa[4][idx],
                &__optd_data[gidx]->kappa[5][idx]);
            assert(!ierr);
            assert(isfinite(__optd_data[gidx]->kappa[0][idx]));
            assert(isfinite(__optd_data[gidx]->kappa[1][idx]));
            assert(isfinite(__optd_data[gidx]->kappa[2][idx]));
            assert(isfinite(__optd_data[gidx]->kappa[3][idx]));
            assert(isfinite(__optd_data[gidx]->kappa[4][idx]));
            assert(isfinite(__optd_data[gidx]->kappa[5][idx]));
        }

        /* Compute grr */
        for(int iray = 0; iray < ntheta*nphi; ++iray)
        for(int irad = 0; irad < nrad; ++irad) {
            int const idx = THC_SPH_INDEX(__optd_data[gidx]->grid, irad, iray);
            CCTK_REAL * grr = &__optd_data[gidx]->grr[idx];

            CCTK_REAL r, theta, phi;
            thc_sph_grid_get_r_theta_phi(__optd_data[gidx]->grid, irad, iray,
                    &r, &theta, &phi);

            /* Three metric in cartesian coordinates  */
            CCTK_REAL gamma_c[9];
            gamma_c[0] = __optd_data[gidx]->metric[0][idx]; /* gxx */
            gamma_c[1] = __optd_data[gidx]->metric[1][idx]; /* gxy */
            gamma_c[2] = __optd_data[gidx]->metric[2][idx]; /* gxz */
            gamma_c[3] = __optd_data[gidx]->metric[1][idx]; /* gxy */
            gamma_c[4] = __optd_data[gidx]->metric[3][idx]; /* gyy */
            gamma_c[5] = __optd_data[gidx]->metric[4][idx]; /* gyz */
            gamma_c[6] = __optd_data[gidx]->metric[2][idx]; /* gxz */
            gamma_c[7] = __optd_data[gidx]->metric[4][idx]; /* gyz */
            gamma_c[8] = __optd_data[gidx]->metric[5][idx]; /* gzz */

            CCTK_REAL Jac[9];
            thc_coord_sph_to_cart_J(r, theta, phi, &Jac[0]);

            *grr = 0;
            for(int i = 0; i < 3; ++i)
            for(int j = 0; j < 3; ++j) {
                *grr +=  Jac[3*i] * Jac[3*j] * gamma_c[3*i + j];
            }
        }

        /* Compute optical depth */
        for(int iray = 0; iray < ntheta*nphi; ++iray) {
            for(int v = 0; v < 6; ++v) {
                int const idx = THC_SPH_INDEX(__optd_data[gidx]->grid, 0, iray);
                CCTK_REAL * grr = &__optd_data[gidx]->grr[idx];
                CCTK_REAL * tau = &__optd_data[gidx]->tau[v][idx];
                CCTK_REAL * kappa = &__optd_data[gidx]->kappa[v][idx];

                tau[nrad-1] = 0;
                for(int irad = nrad - 2; irad >= 0; --irad) {
                    tau[irad] = tau[irad + 1] +
                        drad * (0.5*(kappa[irad] + kappa[irad+1])) *
                            sqrt(0.5*(grr[irad] + grr[irad+1]));
                }
            }
        }
    }
}

void THC_LK_InterpOpticalDepth(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

      /*
    if(verbose) {
        CCTK_INFO("THC_LK_InterpOpticalDepth");
    }
      */
      
    CCTK_REAL const * centers[2] = {&center_grid1[0], &center_grid2[0]};

    /* Allocate interpolation points */
    assert(cctk_lsh[0] == cctk_ash[0]);
    assert(cctk_lsh[1] == cctk_ash[1]);
    assert(cctk_lsh[2] == cctk_ash[2]);
    int const gfsiz     = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];
    CCTK_REAL * p_r     = malloc(gfsiz*sizeof(CCTK_REAL));
    CCTK_REAL * p_theta = malloc(gfsiz*sizeof(CCTK_REAL));
    CCTK_REAL * p_phi   = malloc(gfsiz*sizeof(CCTK_REAL));

    /* We interpolate from each of the spherical grids and then select the
     * pointwise value of the optical depth in each point by looking at which
     * of the two centers is closer */
    CCTK_REAL * my_optd[2][6];
    for(int v = 0; v < 6; ++v) {
        my_optd[0][v] = malloc(gfsiz*sizeof(CCTK_REAL));
        my_optd[1][v] = malloc(gfsiz*sizeof(CCTK_REAL));
    }

    for(int gidx = 0; gidx < num_sph_grids; ++gidx) {
        /* Compute interpolation points */
#pragma omp parallel
        {
            UTILS_LOOP3(thc_lk_calc_interp_points,
                    k, 0, cctk_lsh[2],
                    j, 0, cctk_lsh[1],
                    i, 0, cctk_lsh[0]) {
                int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
                thc_coord_cart_to_sph(
                        x[ijk] - centers[gidx][0],
                        y[ijk] - centers[gidx][1],
                        z[ijk] - centers[gidx][2],
                        &p_r[ijk], &p_theta[ijk], &p_phi[ijk]);
                /* This effectively ensures all of the points out of the spherical
                 * grid are filled with a zeroth order extrapolation */
                p_r[ijk] = MIN(p_r[ijk], rmax);
            } UTILS_ENDLOOP3(thc_lk_calc_interp_points);
        }

        /* Do the actual interpolation */
        int const interp_handle = CCTK_InterpHandle(interpolator);
        assert(interp_handle >= 0);
        int const options_handle =
            Util_TableCreateFromString(interpolator_options);
        assert(options_handle >= 0);

        CCTK_REAL coord_origin[3] = {0, 0, 0};
        CCTK_REAL coord_delta[3];
        thc_sph_grid_get_delta(__optd_data[gidx]->grid, &coord_delta[0],
                &coord_delta[1], &coord_delta[2]);
        void const * interp_coords[] = {p_r, p_theta, p_phi};

        void const ** input_arrays = (void const **)&__optd_data[gidx]->tau[0];
        int const ninput = 6;

        CCTK_INT const input_array_codes[] = {
            CCTK_VARIABLE_REAL,
            CCTK_VARIABLE_REAL,
            CCTK_VARIABLE_REAL,
            CCTK_VARIABLE_REAL,
            CCTK_VARIABLE_REAL,
            CCTK_VARIABLE_REAL
        };
        assert(length(input_array_codes) == ninput);

        CCTK_INT const input_array_dims[] = {nrad, ntheta, nphi};

        CCTK_INT const output_array_codes[] = {
            CCTK_VARIABLE_REAL,
            CCTK_VARIABLE_REAL,
            CCTK_VARIABLE_REAL,
            CCTK_VARIABLE_REAL,
            CCTK_VARIABLE_REAL,
            CCTK_VARIABLE_REAL
        };
        assert(length(output_array_codes) == ninput);

        void ** output_arrays = (void **)&my_optd[gidx][0];

        int ierr = CCTK_InterpLocalUniform(3, interp_handle, options_handle,
                coord_origin, coord_delta, gfsiz, CCTK_VARIABLE_REAL,
                interp_coords, ninput, input_array_dims, input_array_codes,
                input_arrays, ninput, output_array_codes, output_arrays);
        assert(!ierr);
        Util_TableDestroy(options_handle);
    }

    CCTK_REAL * optd_vec[] = {
        optd_0_nue, optd_0_nua, optd_0_nux,
        optd_1_nue, optd_1_nua, optd_1_nux
    };
    CCTK_REAL * optd_vec_p[] = {
        optd_0_nue_p, optd_0_nua_p, optd_0_nux_p,
        optd_1_nue_p, optd_1_nua_p, optd_1_nux_p
    };
    CCTK_REAL * optd_vec_p_p[] = {
        optd_0_nue_p_p, optd_0_nua_p_p, optd_0_nux_p_p,
        optd_1_nue_p_p, optd_1_nua_p_p, optd_1_nux_p_p
    };

    /* Choose which value of the optical depth to use in each point */
    if(num_sph_grids == 1) {
        for(int v = 0; v < 6; ++v) {
            memcpy(optd_vec[v], my_optd[0][v], gfsiz*sizeof(CCTK_REAL));
        }
    }
    else {
        for(int ijk = 0; ijk < gfsiz; ++ijk) {
            CCTK_REAL const d0 = SQ(x[ijk] - center_grid1[0]) +
                SQ(y[ijk] - center_grid1[1]) + SQ(z[ijk] - center_grid1[2]);
            CCTK_REAL const d1 = SQ(x[ijk] - center_grid2[0]) +
                SQ(y[ijk] - center_grid2[1]) + SQ(z[ijk] - center_grid2[2]);
            if(d0 < d1) {
                for(int v = 0; v < 6; ++v) {
                    optd_vec[v][ijk] = my_optd[0][v][ijk];
                }
            }
            else {
                for(int v = 0; v < 6; ++v) {
                    optd_vec[v][ijk] = my_optd[1][v][ijk];
                }
            }
        }
    }

    /* Check optical depth */
    for(int ijk = 0; ijk < gfsiz; ++ijk) {
        assert(isfinite(optd_0_nue[ijk]));
        assert(isfinite(optd_0_nua[ijk]));
        assert(isfinite(optd_0_nux[ijk]));
        assert(isfinite(optd_1_nue[ijk]));
        assert(isfinite(optd_1_nua[ijk]));
        assert(isfinite(optd_1_nux[ijk]));
    }

    /* Copy the optical depth to the previous timelevels */
    if(timelevels > 1) {
        for(int v = 0; v < 6; ++v) {
            memcpy(optd_vec_p[v], optd_vec[v], gfsiz*sizeof(CCTK_REAL));
        }
        if(timelevels > 2) {
            for(int v = 0; v < 6; ++v) {
                memcpy(optd_vec_p_p[v], optd_vec[v], gfsiz*sizeof(CCTK_REAL));
            }
        }
    }

    /* Cleanup */
    free(p_r);
    free(p_theta);
    free(p_phi);
    for(int v = 0; v < 6; ++v) {
        free(my_optd[0][v]);
        free(my_optd[1][v]);
    }
}

void THC_LK_OpticalDepthCleanup(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

      /*
    if(verbose) {
        CCTK_INFO("THC_LK_OpticalDepthCleanup");
    }
      */
      
    for(int gidx = 0; gidx < num_sph_grids; ++gidx) {
        assert(__optd_data[gidx]);
        thc_sph_grid_free(__optd_data[gidx]->grid);
        for(int v = 0; v < 6; ++v) {
            free(__optd_data[gidx]->metric[v]);
        }
        free(__optd_data[gidx]->grr);
        free(__optd_data[gidx]->rho);
        free(__optd_data[gidx]->temp);
        free(__optd_data[gidx]->Y_e);
        for(int v = 0; v < 6; ++v) {
            free(__optd_data[gidx]->kappa[v]);
            free(__optd_data[gidx]->tau[v]);
        }
        for(int d = 0; d < 3; ++d) {
            free(__optd_data[gidx]->xcor[d]);
        }
        free(__optd_data[gidx]);
        __optd_data[gidx] = NULL;
    }
}
