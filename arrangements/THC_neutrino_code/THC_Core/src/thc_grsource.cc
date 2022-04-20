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


#include <cmath>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "finite_difference.h"
#include "utils.hh"

#define SQ(X) ((X)*(X))

extern "C" void THC_GRSource(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_GRSource");
    }

    int const gsiz = UTILS_GFSIZE(cctkGH);
    CCTK_REAL const * velx = &vel[0*gsiz];
    CCTK_REAL const * vely = &vel[1*gsiz];
    CCTK_REAL const * velz = &vel[2*gsiz];
    CCTK_REAL * rhs_sconx  = &rhs_scon[0*gsiz];
    CCTK_REAL * rhs_scony  = &rhs_scon[1*gsiz];
    CCTK_REAL * rhs_sconz  = &rhs_scon[2*gsiz];

    CCTK_REAL const dx  = CCTK_DELTA_SPACE(0);
    CCTK_REAL const dy  = CCTK_DELTA_SPACE(1);
    CCTK_REAL const dz  = CCTK_DELTA_SPACE(2);
    CCTK_REAL const idx = 1.0/dx;
    CCTK_REAL const idy = 1.0/dy;
    CCTK_REAL const idz = 1.0/dz;

    bool const barotropic = CCTK_Equals(eos_type, "barotropic");
    bool const ideal      = CCTK_Equals(eos_type, "ideal");
    bool const nuclear    = CCTK_Equals(eos_type, "nuclear");

    if(!(barotropic || ideal || nuclear)) {
        CCTK_ERROR("This is a bug in thc_grsource.cc");
    }

#pragma omp parallel
    {
        UTILS_LOOP3(thc_grsource,
                k, cctk_nghostzones[2], cctk_lsh[2]-cctk_nghostzones[2],
                j, cctk_nghostzones[1], cctk_lsh[1]-cctk_nghostzones[1],
                i, cctk_nghostzones[0], cctk_lsh[0]-cctk_nghostzones[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            CCTK_REAL const det = utils::metric::spatial_det(gxx[ijk],
                    gxy[ijk], gxz[ijk], gyy[ijk], gyz[ijk], gzz[ijk]);
            CCTK_REAL const sqrt_det = std::sqrt(det);
            CCTK_REAL uxx, uxy, uxz, uyy, uyz, uzz;
            utils::metric::spatial_inv(det, gxx[ijk], gxy[ijk], gxz[ijk],
                    gyy[ijk], gyz[ijk], gzz[ijk], &uxx, &uxy, &uxz, &uyy,
                    &uyz, &uzz);

            CCTK_REAL const rho_h_W2 = (rho[ijk]*(1. + eps[ijk]) + press[ijk])*
                SQ(w_lorentz[ijk]);

            CCTK_REAL const velxshift = velx[ijk] - betax[ijk]/alp[ijk];
            CCTK_REAL const velyshift = vely[ijk] - betay[ijk]/alp[ijk];
            CCTK_REAL const velzshift = velz[ijk] - betaz[ijk]/alp[ijk];

            CCTK_REAL const vlowx = gxx[ijk]*velx[ijk] + gxy[ijk]*vely[ijk] +
                gxz[ijk]*velz[ijk];
            CCTK_REAL const vlowy = gxy[ijk]*velx[ijk] + gyy[ijk]*vely[ijk] +
                gyz[ijk]*velz[ijk];
            CCTK_REAL const vlowz = gxz[ijk]*velx[ijk] + gyz[ijk]*vely[ijk] +
                gzz[ijk]*velz[ijk];

            // Contravariant components of the stress-energy tensor
            CCTK_REAL const T00 = (rho_h_W2 - press[ijk])/SQ(alp[ijk]);

            CCTK_REAL const T0x = rho_h_W2*velxshift/alp[ijk] +
                press[ijk]*betax[ijk]/SQ(alp[ijk]);
            CCTK_REAL const T0y = rho_h_W2*velyshift/alp[ijk] +
                press[ijk]*betay[ijk]/SQ(alp[ijk]);
            CCTK_REAL const T0z = rho_h_W2*velzshift/alp[ijk] +
                press[ijk]*betaz[ijk]/SQ(alp[ijk]);

            CCTK_REAL const Txx = rho_h_W2*velxshift*velxshift +
                press[ijk]*(uxx - betax[ijk]*betax[ijk]/SQ(alp[ijk]));
            CCTK_REAL const Txy = rho_h_W2*velxshift*velyshift +
                press[ijk]*(uxy - betax[ijk]*betay[ijk]/SQ(alp[ijk]));
            CCTK_REAL const Txz = rho_h_W2*velxshift*velzshift +
                press[ijk]*(uxz - betax[ijk]*betaz[ijk]/SQ(alp[ijk]));
            CCTK_REAL const Tyy = rho_h_W2*velyshift*velyshift +
                press[ijk]*(uyy - betay[ijk]*betay[ijk]/SQ(alp[ijk]));
            CCTK_REAL const Tyz = rho_h_W2*velyshift*velzshift +
                press[ijk]*(uyz - betay[ijk]*betaz[ijk]/SQ(alp[ijk]));
            CCTK_REAL const Tzz = rho_h_W2*velzshift*velzshift +
                press[ijk]*(uzz - betaz[ijk]*betaz[ijk]/SQ(alp[ijk]));

            // Derivatives of the lapse, metric and shift
            CCTK_REAL const dx_alp = idx*cdiff_x(cctkGH, alp, i, j, k, fd_order);
            CCTK_REAL const dy_alp = idy*cdiff_y(cctkGH, alp, i, j, k, fd_order);
            CCTK_REAL const dz_alp = idz*cdiff_z(cctkGH, alp, i, j, k, fd_order);

            CCTK_REAL const dx_betax = idx*cdiff_x(cctkGH, betax, i, j, k, fd_order);
            CCTK_REAL const dx_betay = idx*cdiff_x(cctkGH, betay, i, j, k, fd_order);
            CCTK_REAL const dx_betaz = idx*cdiff_x(cctkGH, betaz, i, j, k, fd_order);
            CCTK_REAL const dy_betax = idy*cdiff_y(cctkGH, betax, i, j, k, fd_order);
            CCTK_REAL const dy_betay = idy*cdiff_y(cctkGH, betay, i, j, k, fd_order);
            CCTK_REAL const dy_betaz = idy*cdiff_y(cctkGH, betaz, i, j, k, fd_order);
            CCTK_REAL const dz_betax = idz*cdiff_z(cctkGH, betax, i, j, k, fd_order);
            CCTK_REAL const dz_betay = idz*cdiff_z(cctkGH, betay, i, j, k, fd_order);
            CCTK_REAL const dz_betaz = idz*cdiff_z(cctkGH, betaz, i, j, k, fd_order);

            CCTK_REAL const dx_gxx = idx*cdiff_x(cctkGH, gxx, i, j, k, fd_order);
            CCTK_REAL const dx_gxy = idx*cdiff_x(cctkGH, gxy, i, j, k, fd_order);
            CCTK_REAL const dx_gxz = idx*cdiff_x(cctkGH, gxz, i, j, k, fd_order);
            CCTK_REAL const dx_gyy = idx*cdiff_x(cctkGH, gyy, i, j, k, fd_order);
            CCTK_REAL const dx_gyz = idx*cdiff_x(cctkGH, gyz, i, j, k, fd_order);
            CCTK_REAL const dx_gzz = idx*cdiff_x(cctkGH, gzz, i, j, k, fd_order);
            CCTK_REAL const dy_gxx = idy*cdiff_y(cctkGH, gxx, i, j, k, fd_order);
            CCTK_REAL const dy_gxy = idy*cdiff_y(cctkGH, gxy, i, j, k, fd_order);
            CCTK_REAL const dy_gxz = idy*cdiff_y(cctkGH, gxz, i, j, k, fd_order);
            CCTK_REAL const dy_gyy = idy*cdiff_y(cctkGH, gyy, i, j, k, fd_order);
            CCTK_REAL const dy_gyz = idy*cdiff_y(cctkGH, gyz, i, j, k, fd_order);
            CCTK_REAL const dy_gzz = idy*cdiff_y(cctkGH, gzz, i, j, k, fd_order);
            CCTK_REAL const dz_gxx = idz*cdiff_z(cctkGH, gxx, i, j, k, fd_order);
            CCTK_REAL const dz_gxy = idz*cdiff_z(cctkGH, gxy, i, j, k, fd_order);
            CCTK_REAL const dz_gxz = idz*cdiff_z(cctkGH, gxz, i, j, k, fd_order);
            CCTK_REAL const dz_gyy = idz*cdiff_z(cctkGH, gyy, i, j, k, fd_order);
            CCTK_REAL const dz_gyz = idz*cdiff_z(cctkGH, gyz, i, j, k, fd_order);
            CCTK_REAL const dz_gzz = idz*cdiff_z(cctkGH, gzz, i, j, k, fd_order);

            // Contract the shift with the extrinsic curvature
            CCTK_REAL const shiftshiftk =
                   betax[ijk]*betax[ijk]*kxx[ijk] +
                   betay[ijk]*betay[ijk]*kyy[ijk] +
                   betaz[ijk]*betaz[ijk]*kzz[ijk] +
                2*(betax[ijk]*betay[ijk]*kxy[ijk] +
                   betax[ijk]*betaz[ijk]*kxz[ijk] +
                   betay[ijk]*betaz[ijk]*kyz[ijk]);

            CCTK_REAL const shiftkx = kxx[ijk]*betax[ijk] +
                kxy[ijk]*betay[ijk] + kxz[ijk]*betaz[ijk];
            CCTK_REAL const shiftky = kxy[ijk]*betax[ijk] +
                kyy[ijk]*betay[ijk] + kyz[ijk]*betaz[ijk];
            CCTK_REAL const shiftkz = kxz[ijk]*betax[ijk] +
                kyz[ijk]*betay[ijk] + kzz[ijk]*betaz[ijk];

            // Contract the matter terms with the extrinsic curvature
            CCTK_REAL const sumTK = Txx*kxx[ijk] + Tyy*kyy[ijk] + Tzz*kzz[ijk] +
                2*(Txy*kxy[ijk] + Txz*kxz[ijk] + Tyz*kyz[ijk]);

            // Source term for tau
            CCTK_REAL const tau_source = T00 *
               (shiftshiftk - (betax[ijk]*dx_alp + betay[ijk]*dy_alp +
                                betaz[ijk]*dz_alp)) +
               T0x * (-dx_alp + 2*shiftkx) +
               T0y * (-dy_alp + 2*shiftky) +
               T0z * (-dz_alp + 2*shiftkz) +
               sumTK;

            // Contract the shift with the derivatives of the metric
            CCTK_REAL const halfshiftdgx = 0.5*(betax[ijk]*betax[ijk]*dx_gxx +
                    betay[ijk]*betay[ijk]*dx_gyy + betaz[ijk]*betaz[ijk]*dx_gzz) +
                betax[ijk]*betay[ijk]*dx_gxy + betax[ijk]*betaz[ijk]*dx_gxz +
                betay[ijk]*betaz[ijk]*dx_gyz;
            CCTK_REAL const halfshiftdgy = 0.5*(betax[ijk]*betax[ijk]*dy_gxx +
                    betay[ijk]*betay[ijk]*dy_gyy + betaz[ijk]*betaz[ijk]*dy_gzz) +
                betax[ijk]*betay[ijk]*dy_gxy + betax[ijk]*betaz[ijk]*dy_gxz +
                betay[ijk]*betaz[ijk]*dy_gyz;
            CCTK_REAL const halfshiftdgz = 0.5*(betax[ijk]*betax[ijk]*dz_gxx +
                    betay[ijk]*betay[ijk]*dz_gyy + betaz[ijk]*betaz[ijk]*dz_gzz) +
                betax[ijk]*betay[ijk]*dz_gxy + betax[ijk]*betaz[ijk]*dz_gxz +
                betay[ijk]*betaz[ijk]*dz_gyz;

            // Contract the matter with derivatives of the metric
            CCTK_REAL const halfTdgx = 0.5*(Txx*dx_gxx + Tyy*dx_gyy + Tzz*dx_gzz) +
                Txy*dx_gxy + Txz*dx_gxz + Tyz*dx_gyz;
            CCTK_REAL const halfTdgy = 0.5*(Txx*dy_gxx + Tyy*dy_gyy + Tzz*dy_gzz) +
                Txy*dy_gxy + Txz*dy_gxz + Tyz*dy_gyz;
            CCTK_REAL const halfTdgz = 0.5*(Txx*dz_gxx + Tyy*dz_gyy + Tzz*dz_gzz) +
                Txy*dz_gxy + Txz*dz_gxz + Tyz*dz_gyz;

            // Momentum sources
            CCTK_REAL const sx_source = T00*
                (halfshiftdgx - alp[ijk]*dx_alp) +
                T0x*(betax[ijk]*dx_gxx + betay[ijk]*dx_gxy + betaz[ijk]*dx_gxz) +
                T0y*(betax[ijk]*dx_gxy + betay[ijk]*dx_gyy + betaz[ijk]*dx_gyz) +
                T0z*(betax[ijk]*dx_gxz + betay[ijk]*dx_gyz + betaz[ijk]*dx_gzz) +
                halfTdgx + rho_h_W2 *
                (vlowx*dx_betax + vlowy*dx_betay + vlowz*dx_betaz)/alp[ijk];
            CCTK_REAL const sy_source = T00*
                (halfshiftdgy - alp[ijk]*dy_alp) +
                T0x*(betax[ijk]*dy_gxx + betay[ijk]*dy_gxy + betaz[ijk]*dy_gxz) +
                T0y*(betax[ijk]*dy_gxy + betay[ijk]*dy_gyy + betaz[ijk]*dy_gyz) +
                T0z*(betax[ijk]*dy_gxz + betay[ijk]*dy_gyz + betaz[ijk]*dy_gzz) +
                halfTdgy + rho_h_W2 *
                (vlowx*dy_betax + vlowy*dy_betay + vlowz*dy_betaz)/alp[ijk];
            CCTK_REAL const sz_source = T00*
                (halfshiftdgz - alp[ijk]*dz_alp) +
                T0x*(betax[ijk]*dz_gxx + betay[ijk]*dz_gxy + betaz[ijk]*dz_gxz) +
                T0y*(betax[ijk]*dz_gxy + betay[ijk]*dz_gyy + betaz[ijk]*dz_gyz) +
                T0z*(betax[ijk]*dz_gxz + betay[ijk]*dz_gyz + betaz[ijk]*dz_gzz) +
                halfTdgz + rho_h_W2 *
                (vlowx*dz_betax + vlowy*dz_betay + vlowz*dz_betaz)/alp[ijk];

            rhs_sconx[ijk] += alp[ijk] * sqrt_det * sx_source;
            rhs_scony[ijk] += alp[ijk] * sqrt_det * sy_source;
            rhs_sconz[ijk] += alp[ijk] * sqrt_det * sz_source;

            if(!barotropic) {
                rhs_tau[ijk] += alp[ijk] * sqrt_det * tau_source;
            }
        } UTILS_ENDLOOP3(thc_grsource);
    }
}
