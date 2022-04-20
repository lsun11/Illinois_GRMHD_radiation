//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
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


#ifndef THC_MACRO_HH
#define THC_MACRO_HH

#include "hrscc.hh"
#include "utils.hh"

#ifdef THC_DEBUG
#warning "THC is compiled in debugging mode"
#endif

#define SQ(X) (utils::pow<2>((X)))
#define CB(X) (utils::pow<3>((X)))

#define UNUSED(X) ((void)X)

#define ISFINITE(X) (!std::isnan(X) && !std::isinf(X))

#define THC_FLAG_ATMOSPHERE                 1
#define THC_FLAG_C2A_FAILURE                2
#define THC_FLAG_RHO_LT_RHOMIN              4
#define THC_FLAG_RHO_GT_RHOMAX              8
#define THC_FLAG_EPS_LT_EPSMIN             16
#define THC_FLAG_EPS_GT_EPSMAX             32
#define THC_FLAG_YE_LT_YEMIN               64
#define THC_FLAG_YE_GT_YEMAX              128
#define THC_FLAG_CSOUND_IS_CRAP           256
#define THC_FLAG_IS_INSIDE_BH             512
#define THC_FLAG_NOT_FINITE              1024

namespace thc {

// this template is used to help when rotating velocity/fluxes
template<hrscc::policy::direction_t dir>
class index;

template<>
class index<hrscc::policy::x> {
    public:
        enum {x = 0};
        enum {y = 1};
        enum {z = 2};

        enum {xx = 0};
        enum {xy = 1};
        enum {xz = 2};
        enum {yy = 3};
        enum {yz = 4};
        enum {zz = 5};
};

template<>
class index<hrscc::policy::y> {
    public:
        enum {x = 1};
        enum {y = 2};
        enum {z = 0};

        enum {xx = 3};
        enum {xy = 4};
        enum {xz = 1};
        enum {yy = 5};
        enum {yz = 2};
        enum {zz = 0};
};

template<>
class index<hrscc::policy::z> {
    public:
        enum {x = 2};
        enum {y = 0};
        enum {z = 1};

        enum {xx = 5};
        enum {xy = 2};
        enum {xz = 4};
        enum {yy = 0};
        enum {yz = 1};
        enum {zz = 3};
};

} // namespace

// handy alias for the ultrarelativistic hydrodynamics
#define ULR_MAKE_ALIAS(dir,observer)                                           \
    cGH const * cctkGH     = observer.cctkGH;                                  \
    CCTK_REAL & tau        = observer.conserved[0];                            \
    CCTK_REAL & sconx      = observer.conserved[1+thc::index<dir>::x];         \
    CCTK_REAL & scony      = observer.conserved[1+thc::index<dir>::y];         \
    CCTK_REAL & sconz      = observer.conserved[1+thc::index<dir>::z];         \
    CCTK_REAL & rho        = observer.primitive[0];                            \
    CCTK_REAL & velx       = observer.primitive[1+thc::index<dir>::x];         \
    CCTK_REAL & vely       = observer.primitive[1+thc::index<dir>::y];         \
    CCTK_REAL & velz       = observer.primitive[1+thc::index<dir>::z];         \
    CCTK_REAL & flux_tau   = observer.flux[dir][0];                            \
    CCTK_REAL & flux_sconx = observer.flux[dir][1+thc::index<dir>::x];         \
    CCTK_REAL & flux_scony = observer.flux[dir][1+thc::index<dir>::y];         \
    CCTK_REAL & flux_sconz = observer.flux[dir][1+thc::index<dir>::z];         \
    CCTK_REAL & eps        = observer.field[0];                                \
    CCTK_REAL & press      = observer.field[1];                                \
    CCTK_REAL & w_lorentz  = observer.field[2];                                \
    CCTK_REAL & csound     = observer.field[3];                                \
    CCTK_REAL * const right_eigenvector[4] = {                                 \
        &observer.right_eigenvector[0][0],                                     \
        &observer.right_eigenvector[1+thc::index<dir>::x][0],                  \
        &observer.right_eigenvector[1+thc::index<dir>::y][0],                  \
        &observer.right_eigenvector[1+thc::index<dir>::z][0]};                 \
    CCTK_REAL * const left_eigenvector[4][4] = {                               \
        {                                                                      \
            &observer.left_eigenvector[0][0],                                  \
            &observer.left_eigenvector[0][1+thc::index<dir>::x],               \
            &observer.left_eigenvector[0][1+thc::index<dir>::y],               \
            &observer.left_eigenvector[0][1+thc::index<dir>::z]                \
        },                                                                     \
        {                                                                      \
            &observer.left_eigenvector[1][0],                                  \
            &observer.left_eigenvector[1][1+thc::index<dir>::x],               \
            &observer.left_eigenvector[1][1+thc::index<dir>::y],               \
            &observer.left_eigenvector[1][1+thc::index<dir>::z]                \
        },                                                                     \
        {                                                                      \
            &observer.left_eigenvector[2][0],                                  \
            &observer.left_eigenvector[2][1+thc::index<dir>::x],               \
            &observer.left_eigenvector[2][1+thc::index<dir>::y],               \
            &observer.left_eigenvector[2][1+thc::index<dir>::z]                \
        },                                                                     \
        {                                                                      \
            &observer.left_eigenvector[3][0],                                  \
            &observer.left_eigenvector[3][1+thc::index<dir>::x],               \
            &observer.left_eigenvector[3][1+thc::index<dir>::y],               \
            &observer.left_eigenvector[3][1+thc::index<dir>::z]                \
        }                                                                      \
        };                                                                     \
    CCTK_INT & bitmask = observer.bitmask[0];                                  \
                                                                               \
    UNUSED(cctkGH);                                                            \
    UNUSED(sconx);                                                             \
    UNUSED(scony);                                                             \
    UNUSED(sconz);                                                             \
    UNUSED(tau);                                                               \
    UNUSED(rho);                                                               \
    UNUSED(velx);                                                              \
    UNUSED(vely);                                                              \
    UNUSED(velz);                                                              \
    UNUSED(flux_sconx);                                                        \
    UNUSED(flux_scony);                                                        \
    UNUSED(flux_sconz);                                                        \
    UNUSED(flux_tau);                                                          \
    UNUSED(eps);                                                               \
    UNUSED(press);                                                             \
    UNUSED(w_lorentz);                                                         \
    UNUSED(csound);                                                            \
    UNUSED(right_eigenvector);                                                 \
    UNUSED(left_eigenvector);                                                  \
    UNUSED(bitmask)

// handy alias for the barotropic EoS
#define BRT_MAKE_ALIAS(dir,observer)                                           \
    cGH const * cctkGH      = observer.cctkGH;                                 \
    CCTK_REAL & dens        = observer.conserved[0];                           \
    CCTK_REAL & sconx       = observer.conserved[1+thc::index<dir>::x];        \
    CCTK_REAL & scony       = observer.conserved[1+thc::index<dir>::y];        \
    CCTK_REAL & sconz       = observer.conserved[1+thc::index<dir>::z];        \
    CCTK_REAL & rho         = observer.primitive[0];                           \
    CCTK_REAL & zvecx       = observer.primitive[1+thc::index<dir>::x];        \
    CCTK_REAL & zvecy       = observer.primitive[1+thc::index<dir>::y];        \
    CCTK_REAL & zvecz       = observer.primitive[1+thc::index<dir>::z];        \
    CCTK_REAL & flux_dens   = observer.flux[dir][0];                           \
    CCTK_REAL & flux_sconx  = observer.flux[dir][1+thc::index<dir>::x];        \
    CCTK_REAL & flux_scony  = observer.flux[dir][1+thc::index<dir>::y];        \
    CCTK_REAL & flux_sconz  = observer.flux[dir][1+thc::index<dir>::z];        \
    CCTK_REAL & alp         = observer.field[0];                               \
    CCTK_REAL & betax       = observer.field[1+thc::index<dir>::x];            \
    CCTK_REAL & betay       = observer.field[1+thc::index<dir>::y];            \
    CCTK_REAL & betaz       = observer.field[1+thc::index<dir>::z];            \
    CCTK_REAL & gxx         = observer.field[4+thc::index<dir>::xx];           \
    CCTK_REAL & gxy         = observer.field[4+thc::index<dir>::xy];           \
    CCTK_REAL & gxz         = observer.field[4+thc::index<dir>::xz];           \
    CCTK_REAL & gyy         = observer.field[4+thc::index<dir>::yy];           \
    CCTK_REAL & gyz         = observer.field[4+thc::index<dir>::yz];           \
    CCTK_REAL & gzz         = observer.field[4+thc::index<dir>::zz];           \
    CCTK_REAL & eps         = observer.field[10];                              \
    CCTK_REAL & press       = observer.field[11];                              \
    CCTK_REAL & velx        = observer.field[12+thc::index<dir>::x];           \
    CCTK_REAL & vely        = observer.field[12+thc::index<dir>::y];           \
    CCTK_REAL & velz        = observer.field[12+thc::index<dir>::z];           \
    CCTK_REAL & w_lorentz   = observer.field[15];                              \
    CCTK_REAL & c2a_nbiter  = observer.field[16];                              \
    CCTK_REAL & csound      = observer.field[17];                              \
    CCTK_REAL & densgain    = observer.field[18];                              \
    CCTK_REAL & volform     = observer.field[19];                              \
    CCTK_REAL * const right_eigenvector[4] = {                                 \
        &observer.right_eigenvector[0][0],                                     \
        &observer.right_eigenvector[1+thc::index<dir>::x][0],                  \
        &observer.right_eigenvector[1+thc::index<dir>::y][0],                  \
        &observer.right_eigenvector[1+thc::index<dir>::z][0]};                 \
    CCTK_REAL * const left_eigenvector[4][4] = {                               \
        {                                                                      \
            &observer.left_eigenvector[0][0],                                  \
            &observer.left_eigenvector[0][1+thc::index<dir>::x],               \
            &observer.left_eigenvector[0][1+thc::index<dir>::y],               \
            &observer.left_eigenvector[0][1+thc::index<dir>::z]                \
        },                                                                     \
        {                                                                      \
            &observer.left_eigenvector[1][0],                                  \
            &observer.left_eigenvector[1][1+thc::index<dir>::x],               \
            &observer.left_eigenvector[1][1+thc::index<dir>::y],               \
            &observer.left_eigenvector[1][1+thc::index<dir>::z]                \
        },                                                                     \
        {                                                                      \
            &observer.left_eigenvector[2][0],                                  \
            &observer.left_eigenvector[2][1+thc::index<dir>::x],               \
            &observer.left_eigenvector[2][1+thc::index<dir>::y],               \
            &observer.left_eigenvector[2][1+thc::index<dir>::z]                \
        },                                                                     \
        {                                                                      \
            &observer.left_eigenvector[3][0],                                  \
            &observer.left_eigenvector[3][1+thc::index<dir>::x],               \
            &observer.left_eigenvector[3][1+thc::index<dir>::y],               \
            &observer.left_eigenvector[3][1+thc::index<dir>::z]                \
        }                                                                      \
        };                                                                     \
    CCTK_INT & bitmask = observer.bitmask[0];                                  \
                                                                               \
    UNUSED(cctkGH);                                                            \
    UNUSED(dens);                                                              \
    UNUSED(sconx);                                                             \
    UNUSED(scony);                                                             \
    UNUSED(sconz);                                                             \
    UNUSED(rho);                                                               \
    UNUSED(zvecx);                                                             \
    UNUSED(zvecy);                                                             \
    UNUSED(zvecz);                                                             \
    UNUSED(flux_dens);                                                         \
    UNUSED(flux_sconx);                                                        \
    UNUSED(flux_scony);                                                        \
    UNUSED(flux_sconz);                                                        \
    UNUSED(alp);                                                               \
    UNUSED(betax);                                                             \
    UNUSED(betay);                                                             \
    UNUSED(betaz);                                                             \
    UNUSED(gxx);                                                               \
    UNUSED(gxy);                                                               \
    UNUSED(gxz);                                                               \
    UNUSED(gyy);                                                               \
    UNUSED(gyz);                                                               \
    UNUSED(gzz);                                                               \
    UNUSED(eps);                                                               \
    UNUSED(press);                                                             \
    UNUSED(velx);                                                              \
    UNUSED(vely);                                                              \
    UNUSED(velz);                                                              \
    UNUSED(w_lorentz);                                                         \
    UNUSED(c2a_nbiter);                                                        \
    UNUSED(csound);                                                            \
    UNUSED(densgain);                                                          \
    UNUSED(volform);                                                           \
    UNUSED(right_eigenvector);                                                 \
    UNUSED(left_eigenvector);                                                  \
    UNUSED(bitmask)

// handy alias for the ideal EoS
#define IDL_MAKE_ALIAS(dir,observer)                                           \
    cGH const * cctkGH      = observer.cctkGH;                                 \
    CCTK_REAL & dens        = observer.conserved[0];                           \
    CCTK_REAL & sconx       = observer.conserved[1+thc::index<dir>::x];        \
    CCTK_REAL & scony       = observer.conserved[1+thc::index<dir>::y];        \
    CCTK_REAL & sconz       = observer.conserved[1+thc::index<dir>::z];        \
    CCTK_REAL & tau         = observer.conserved[4];                           \
    CCTK_REAL & rho         = observer.primitive[0];                           \
    CCTK_REAL & zvecx       = observer.primitive[1+thc::index<dir>::x];        \
    CCTK_REAL & zvecy       = observer.primitive[1+thc::index<dir>::y];        \
    CCTK_REAL & zvecz       = observer.primitive[1+thc::index<dir>::z];        \
    CCTK_REAL & eps         = observer.primitive[4];                           \
    CCTK_REAL & flux_dens   = observer.flux[dir][0];                           \
    CCTK_REAL & flux_sconx  = observer.flux[dir][1+thc::index<dir>::x];        \
    CCTK_REAL & flux_scony  = observer.flux[dir][1+thc::index<dir>::y];        \
    CCTK_REAL & flux_sconz  = observer.flux[dir][1+thc::index<dir>::z];        \
    CCTK_REAL & flux_tau    = observer.flux[dir][4];                           \
    CCTK_REAL & alp         = observer.field[0];                               \
    CCTK_REAL & betax       = observer.field[1+thc::index<dir>::x];            \
    CCTK_REAL & betay       = observer.field[1+thc::index<dir>::y];            \
    CCTK_REAL & betaz       = observer.field[1+thc::index<dir>::z];            \
    CCTK_REAL & gxx         = observer.field[4+thc::index<dir>::xx];           \
    CCTK_REAL & gxy         = observer.field[4+thc::index<dir>::xy];           \
    CCTK_REAL & gxz         = observer.field[4+thc::index<dir>::xz];           \
    CCTK_REAL & gyy         = observer.field[4+thc::index<dir>::yy];           \
    CCTK_REAL & gyz         = observer.field[4+thc::index<dir>::yz];           \
    CCTK_REAL & gzz         = observer.field[4+thc::index<dir>::zz];           \
    CCTK_REAL & press       = observer.field[10];                              \
    CCTK_REAL & velx        = observer.field[11+thc::index<dir>::x];           \
    CCTK_REAL & vely        = observer.field[11+thc::index<dir>::y];           \
    CCTK_REAL & velz        = observer.field[11+thc::index<dir>::z];           \
    CCTK_REAL & w_lorentz   = observer.field[14];                              \
    CCTK_REAL & c2a_nbiter  = observer.field[15];                              \
    CCTK_REAL & csound      = observer.field[16];                              \
    CCTK_REAL & densgain    = observer.field[17];                              \
    CCTK_REAL & volform     = observer.field[18];                              \
    CCTK_REAL * const right_eigenvector[5] = {                                 \
        &observer.right_eigenvector[0][0],                                     \
        &observer.right_eigenvector[1+thc::index<dir>::x][0],                  \
        &observer.right_eigenvector[1+thc::index<dir>::y][0],                  \
        &observer.right_eigenvector[1+thc::index<dir>::z][0],                  \
        &observer.right_eigenvector[4][0]                                      \
    };                                                                         \
    CCTK_REAL * const left_eigenvector[5][5] = {                               \
        {                                                                      \
            &observer.left_eigenvector[0][0],                                  \
            &observer.left_eigenvector[0][1+thc::index<dir>::x],               \
            &observer.left_eigenvector[0][1+thc::index<dir>::y],               \
            &observer.left_eigenvector[0][1+thc::index<dir>::z],               \
            &observer.left_eigenvector[0][4]                                   \
        },                                                                     \
        {                                                                      \
            &observer.left_eigenvector[1][0],                                  \
            &observer.left_eigenvector[1][1+thc::index<dir>::x],               \
            &observer.left_eigenvector[1][1+thc::index<dir>::y],               \
            &observer.left_eigenvector[1][1+thc::index<dir>::z],               \
            &observer.left_eigenvector[1][4]                                   \
        },                                                                     \
        {                                                                      \
            &observer.left_eigenvector[2][0],                                  \
            &observer.left_eigenvector[2][1+thc::index<dir>::x],               \
            &observer.left_eigenvector[2][1+thc::index<dir>::y],               \
            &observer.left_eigenvector[2][1+thc::index<dir>::z],               \
            &observer.left_eigenvector[2][4]                                   \
        },                                                                     \
        {                                                                      \
            &observer.left_eigenvector[3][0],                                  \
            &observer.left_eigenvector[3][1+thc::index<dir>::x],               \
            &observer.left_eigenvector[3][1+thc::index<dir>::y],               \
            &observer.left_eigenvector[3][1+thc::index<dir>::z],               \
            &observer.left_eigenvector[3][4]                                   \
        },                                                                     \
        {                                                                      \
            &observer.left_eigenvector[4][0],                                  \
            &observer.left_eigenvector[4][1+thc::index<dir>::x],               \
            &observer.left_eigenvector[4][1+thc::index<dir>::y],               \
            &observer.left_eigenvector[4][1+thc::index<dir>::z],               \
            &observer.left_eigenvector[4][4]                                   \
        }                                                                      \
        };                                                                     \
    CCTK_INT & bitmask = observer.bitmask[0];                                  \
                                                                               \
    UNUSED(cctkGH);                                                            \
    UNUSED(dens);                                                              \
    UNUSED(sconx);                                                             \
    UNUSED(scony);                                                             \
    UNUSED(sconz);                                                             \
    UNUSED(tau);                                                               \
    UNUSED(rho);                                                               \
    UNUSED(zvecx);                                                             \
    UNUSED(zvecy);                                                             \
    UNUSED(zvecz);                                                             \
    UNUSED(eps);                                                               \
    UNUSED(flux_dens);                                                         \
    UNUSED(flux_sconx);                                                        \
    UNUSED(flux_scony);                                                        \
    UNUSED(flux_sconz);                                                        \
    UNUSED(flux_tau);                                                          \
    UNUSED(alp);                                                               \
    UNUSED(betax);                                                             \
    UNUSED(betay);                                                             \
    UNUSED(betaz);                                                             \
    UNUSED(gxx);                                                               \
    UNUSED(gxy);                                                               \
    UNUSED(gxz);                                                               \
    UNUSED(gyy);                                                               \
    UNUSED(gyz);                                                               \
    UNUSED(gzz);                                                               \
    UNUSED(press);                                                             \
    UNUSED(velx);                                                              \
    UNUSED(vely);                                                              \
    UNUSED(velz);                                                              \
    UNUSED(w_lorentz);                                                         \
    UNUSED(c2a_nbiter);                                                        \
    UNUSED(csound);                                                            \
    UNUSED(densgain);                                                          \
    UNUSED(volform);                                                           \
    UNUSED(right_eigenvector);                                                 \
    UNUSED(left_eigenvector);                                                  \
    UNUSED(bitmask)

// handy alias for the nuclear EoS
#define NUC_MAKE_ALIAS(dir,observer)                                           \
    cGH const * cctkGH      = observer.cctkGH;                                 \
    CCTK_REAL & densxn      = observer.conserved[0];                           \
    CCTK_REAL & densxp      = observer.conserved[1];                           \
    CCTK_REAL & sconx       = observer.conserved[2+thc::index<dir>::x];        \
    CCTK_REAL & scony       = observer.conserved[2+thc::index<dir>::y];        \
    CCTK_REAL & sconz       = observer.conserved[2+thc::index<dir>::z];        \
    CCTK_REAL & tau         = observer.conserved[5];                           \
    CCTK_REAL & rho         = observer.primitive[0];                           \
    CCTK_REAL & zvecx       = observer.primitive[1+thc::index<dir>::x];        \
    CCTK_REAL & zvecy       = observer.primitive[1+thc::index<dir>::y];        \
    CCTK_REAL & zvecz       = observer.primitive[1+thc::index<dir>::z];        \
    CCTK_REAL & eps         = observer.primitive[4];                           \
    CCTK_REAL & Y_e         = observer.primitive[5];                           \
    CCTK_REAL & flux_densxn = observer.flux[dir][0];                           \
    CCTK_REAL & flux_densxp = observer.flux[dir][1];                           \
    CCTK_REAL & flux_sconx  = observer.flux[dir][2+thc::index<dir>::x];        \
    CCTK_REAL & flux_scony  = observer.flux[dir][2+thc::index<dir>::y];        \
    CCTK_REAL & flux_sconz  = observer.flux[dir][2+thc::index<dir>::z];        \
    CCTK_REAL & flux_tau    = observer.flux[dir][5];                           \
    CCTK_REAL & alp         = observer.field[0];                               \
    CCTK_REAL & betax       = observer.field[1+thc::index<dir>::x];            \
    CCTK_REAL & betay       = observer.field[1+thc::index<dir>::y];            \
    CCTK_REAL & betaz       = observer.field[1+thc::index<dir>::z];            \
    CCTK_REAL & gxx         = observer.field[4+thc::index<dir>::xx];           \
    CCTK_REAL & gxy         = observer.field[4+thc::index<dir>::xy];           \
    CCTK_REAL & gxz         = observer.field[4+thc::index<dir>::xz];           \
    CCTK_REAL & gyy         = observer.field[4+thc::index<dir>::yy];           \
    CCTK_REAL & gyz         = observer.field[4+thc::index<dir>::yz];           \
    CCTK_REAL & gzz         = observer.field[4+thc::index<dir>::zz];           \
    CCTK_REAL & entropy     = observer.field[10];                              \
    CCTK_REAL & press       = observer.field[11];                              \
    CCTK_REAL & temperature = observer.field[12];                              \
    CCTK_REAL & velx        = observer.field[13+thc::index<dir>::x];           \
    CCTK_REAL & vely        = observer.field[13+thc::index<dir>::y];           \
    CCTK_REAL & velz        = observer.field[13+thc::index<dir>::z];           \
    CCTK_REAL & w_lorentz   = observer.field[16];                              \
    CCTK_REAL & c2a_nbiter  = observer.field[17];                              \
    CCTK_REAL & csound      = observer.field[18];                              \
    CCTK_REAL & dens        = observer.field[19];                              \
    CCTK_REAL & densgain    = observer.field[20];                              \
    CCTK_REAL & volform     = observer.field[21];                              \
    CCTK_REAL * const right_eigenvector[6] = {                                 \
        &observer.right_eigenvector[0][0],                                     \
        &observer.right_eigenvector[1][0],                                     \
        &observer.right_eigenvector[2+thc::index<dir>::x][0],                  \
        &observer.right_eigenvector[2+thc::index<dir>::y][0],                  \
        &observer.right_eigenvector[2+thc::index<dir>::z][0],                  \
        &observer.right_eigenvector[5][0]                                      \
    };                                                                         \
    CCTK_REAL * const left_eigenvector[6][6] = {                               \
        {                                                                      \
            &observer.left_eigenvector[0][0],                                  \
            &observer.left_eigenvector[0][1],                                  \
            &observer.left_eigenvector[0][2+thc::index<dir>::x],               \
            &observer.left_eigenvector[0][2+thc::index<dir>::y],               \
            &observer.left_eigenvector[0][2+thc::index<dir>::z],               \
            &observer.left_eigenvector[0][5]                                   \
        },                                                                     \
        {                                                                      \
            &observer.left_eigenvector[1][0],                                  \
            &observer.left_eigenvector[1][1],                                  \
            &observer.left_eigenvector[1][2+thc::index<dir>::x],               \
            &observer.left_eigenvector[1][2+thc::index<dir>::y],               \
            &observer.left_eigenvector[1][2+thc::index<dir>::z],               \
            &observer.left_eigenvector[1][5]                                   \
        },                                                                     \
        {                                                                      \
            &observer.left_eigenvector[2][0],                                  \
            &observer.left_eigenvector[2][1],                                  \
            &observer.left_eigenvector[2][2+thc::index<dir>::x],               \
            &observer.left_eigenvector[2][2+thc::index<dir>::y],               \
            &observer.left_eigenvector[2][2+thc::index<dir>::z],               \
            &observer.left_eigenvector[2][5]                                   \
        },                                                                     \
        {                                                                      \
            &observer.left_eigenvector[3][0],                                  \
            &observer.left_eigenvector[3][1],                                  \
            &observer.left_eigenvector[3][2+thc::index<dir>::x],               \
            &observer.left_eigenvector[3][2+thc::index<dir>::y],               \
            &observer.left_eigenvector[3][2+thc::index<dir>::z],               \
            &observer.left_eigenvector[3][5]                                   \
        },                                                                     \
        {                                                                      \
            &observer.left_eigenvector[4][0],                                  \
            &observer.left_eigenvector[4][1],                                  \
            &observer.left_eigenvector[4][2+thc::index<dir>::x],               \
            &observer.left_eigenvector[4][2+thc::index<dir>::y],               \
            &observer.left_eigenvector[4][2+thc::index<dir>::z],               \
            &observer.left_eigenvector[4][5]                                   \
        },                                                                     \
        {                                                                      \
            &observer.left_eigenvector[5][0],                                  \
            &observer.left_eigenvector[5][1],                                  \
            &observer.left_eigenvector[5][2+thc::index<dir>::x],               \
            &observer.left_eigenvector[5][2+thc::index<dir>::y],               \
            &observer.left_eigenvector[5][2+thc::index<dir>::z],               \
            &observer.left_eigenvector[5][5]                                   \
        }                                                                      \
        };                                                                     \
    CCTK_INT & bitmask = observer.bitmask[0];                                  \
                                                                               \
    UNUSED(cctkGH);                                                            \
    UNUSED(densxn);                                                            \
    UNUSED(densxp);                                                            \
    UNUSED(sconx);                                                             \
    UNUSED(scony);                                                             \
    UNUSED(sconz);                                                             \
    UNUSED(tau);                                                               \
    UNUSED(rho);                                                               \
    UNUSED(zvecx);                                                             \
    UNUSED(zvecy);                                                             \
    UNUSED(zvecz);                                                             \
    UNUSED(eps);                                                               \
    UNUSED(Y_e);                                                               \
    UNUSED(flux_densxn);                                                       \
    UNUSED(flux_densxp);                                                       \
    UNUSED(flux_sconx);                                                        \
    UNUSED(flux_scony);                                                        \
    UNUSED(flux_sconz);                                                        \
    UNUSED(flux_tau);                                                          \
    UNUSED(alp);                                                               \
    UNUSED(betax);                                                             \
    UNUSED(betay);                                                             \
    UNUSED(betaz);                                                             \
    UNUSED(gxx);                                                               \
    UNUSED(gxy);                                                               \
    UNUSED(gxz);                                                               \
    UNUSED(gyy);                                                               \
    UNUSED(gyz);                                                               \
    UNUSED(gzz);                                                               \
    UNUSED(entropy);                                                           \
    UNUSED(press);                                                             \
    UNUSED(temperature);                                                       \
    UNUSED(velx);                                                              \
    UNUSED(vely);                                                              \
    UNUSED(velz);                                                              \
    UNUSED(w_lorentz);                                                         \
    UNUSED(c2a_nbiter);                                                        \
    UNUSED(csound);                                                            \
    UNUSED(dens);                                                              \
    UNUSED(densgain);                                                          \
    UNUSED(volform);                                                           \
    UNUSED(right_eigenvector);                                                 \
    UNUSED(left_eigenvector);                                                  \
    UNUSED(bitmask)

#define L_eigenvector *left_eigenvector
#define R_eigenvector right_eigenvector
#define p_eigenvalue observer.eigenvalue

#endif
