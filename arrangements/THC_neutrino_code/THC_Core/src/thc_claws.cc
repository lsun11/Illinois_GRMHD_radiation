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


#include <limits>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "thc_claws.hh"

//////////////////////////////////////////////////////////////////////////////
// Ultrarelativistic EoS
//////////////////////////////////////////////////////////////////////////////
namespace thc {

template<>
char const * EulerBase<ultrarelativistic::Euler>::conserved_name[4] = {
    "THC_Core::tau",
    "THC_Core::scon[0]",
    "THC_Core::scon[1]",
    "THC_Core::scon[2]"
};

template<>
char const * EulerBase<ultrarelativistic::Euler>::primitive_name[4] = {
    "HydroBase::rho",
    "HydroBase::vel[0]",
    "HydroBase::vel[1]",
    "HydroBase::vel[2]"
};
template<>
char const * EulerBase<ultrarelativistic::Euler>::field_name[4] = {
    "HydroBase::eps",
    "HydroBase::press",
    "HydroBase::w_lorentz",
    "THC_Core::csound",
};
template<>
char const * EulerBase<ultrarelativistic::Euler>::rhs_name[4] = {
    "THC_Core::rhs_tau",
    "THC_Core::rhs_scon[0]",
    "THC_Core::rhs_scon[1]",
    "THC_Core::rhs_scon[2]"
};
template<>
char const * EulerBase<ultrarelativistic::Euler>::bitmask_name[1] = {
    "THC_Core::bitmask"
};

} // namespace thc

namespace hrscc {

template<>
int CLaw<thc::ultrarelativistic::Euler>::conserved_idx[4] = {0, 0, 0, 0};
template<>
int CLaw<thc::ultrarelativistic::Euler>::primitive_idx[4] = {0, 0, 0, 0};
template<>
int CLaw<thc::ultrarelativistic::Euler>::rhs_idx[4] = {0, 0, 0, 0};
template<>
int CLaw<thc::ultrarelativistic::Euler>::field_idx[4] = {0, 0, 0, 0};
template<>
int CLaw<thc::ultrarelativistic::Euler>::bitmask_idx[1] = {0};
template<>
int CLaw<thc::ultrarelativistic::Euler>::num_flux_idx[4*3] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};
template<>
CCTK_REAL CLaw<thc::ultrarelativistic::Euler>::conserved_lbound[4] = {
    std::numeric_limits<CCTK_REAL>::quiet_NaN(),
    std::numeric_limits<CCTK_REAL>::quiet_NaN(),
    std::numeric_limits<CCTK_REAL>::quiet_NaN(),
    std::numeric_limits<CCTK_REAL>::quiet_NaN()
};

} // namespace hrscc
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Barotropic EoS
//////////////////////////////////////////////////////////////////////////////
namespace thc {

template<>
char const * EulerBase<barotropic::Euler>::conserved_name[4] = {
    "THC_Core::dens",
    "THC_Core::scon[0]",
    "THC_Core::scon[1]",
    "THC_Core::scon[2]"
};

template<>
char const * EulerBase<barotropic::Euler>::primitive_name[4] = {
    "HydroBase::rho",
    "THC_Core::zvec[0]",
    "THC_Core::zvec[1]",
    "THC_Core::zvec[2]"
};
template<>
char const * EulerBase<barotropic::Euler>::field_name[20] = {
    "ADMBase::alp",
    "ADMBase::betax",
    "ADMBase::betay",
    "ADMBase::betaz",
    "ADMBase::gxx",
    "ADMBase::gxy",
    "ADMBase::gxz",
    "ADMBase::gyy",
    "ADMBase::gyz",
    "ADMBase::gzz",
    "HydroBase::eps",
    "HydroBase::press",
    "HydroBase::vel[0]",
    "HydroBase::vel[1]",
    "HydroBase::vel[2]",
    "HydroBase::w_lorentz",
    "THC_Core::c2a_nbiter",
    "THC_Core::csound",
    "THC_Core::c2a_densgain",
    "THC_Core::volform"
};
template<>
char const * EulerBase<barotropic::Euler>::rhs_name[4] = {
    "THC_Core::rhs_dens",
    "THC_Core::rhs_scon[0]",
    "THC_Core::rhs_scon[1]",
    "THC_Core::rhs_scon[2]"
};
template<>
char const * EulerBase<barotropic::Euler>::bitmask_name[1] = {
    "THC_Core::bitmask"
};

} // namespace thc

namespace hrscc {

template<>
int CLaw<thc::barotropic::Euler>::conserved_idx[4] = {0, 0, 0, 0};
template<>
int CLaw<thc::barotropic::Euler>::primitive_idx[4] = {0, 0, 0, 0};
template<>
int CLaw<thc::barotropic::Euler>::rhs_idx[4] = {0, 0, 0, 0};
template<>
int CLaw<thc::barotropic::Euler>::field_idx[20] =
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
template<>
int CLaw<thc::barotropic::Euler>::bitmask_idx[1] = {0};
template<>
int CLaw<thc::barotropic::Euler>::num_flux_idx[4*3] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};
template<>
CCTK_REAL CLaw<thc::barotropic::Euler>::conserved_lbound[4] = {
    std::numeric_limits<CCTK_REAL>::quiet_NaN(),
    std::numeric_limits<CCTK_REAL>::quiet_NaN(),
    std::numeric_limits<CCTK_REAL>::quiet_NaN(),
    std::numeric_limits<CCTK_REAL>::quiet_NaN()
};

} // namespace hrscc

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Ideal EoS
//////////////////////////////////////////////////////////////////////////////
namespace thc {

template<>
char const * EulerBase<ideal::Euler>::conserved_name[5] = {
    "THC_Core::dens",
    "THC_Core::scon[0]",
    "THC_Core::scon[1]",
    "THC_Core::scon[2]",
    "THC_Core::tau"
};

template<>
char const * EulerBase<ideal::Euler>::primitive_name[5] = {
    "HydroBase::rho",
    "THC_Core::zvec[0]",
    "THC_Core::zvec[1]",
    "THC_Core::zvec[2]",
    "HydroBase::eps"
};
template<>
char const * EulerBase<ideal::Euler>::field_name[19] = {
    "ADMBase::alp",
    "ADMBase::betax",
    "ADMBase::betay",
    "ADMBase::betaz",
    "ADMBase::gxx",
    "ADMBase::gxy",
    "ADMBase::gxz",
    "ADMBase::gyy",
    "ADMBase::gyz",
    "ADMBase::gzz",
    "HydroBase::press",
    "HydroBase::vel[0]",
    "HydroBase::vel[1]",
    "HydroBase::vel[2]",
    "HydroBase::w_lorentz",
    "THC_Core::c2a_nbiter",
    "THC_Core::csound",
    "THC_Core::c2a_densgain",
    "THC_Core::volform"
};
template<>
char const * EulerBase<ideal::Euler>::rhs_name[5] = {
    "THC_Core::rhs_dens",
    "THC_Core::rhs_scon[0]",
    "THC_Core::rhs_scon[1]",
    "THC_Core::rhs_scon[2]",
    "THC_Core::rhs_tau",
};
template<>
char const * EulerBase<ideal::Euler>::bitmask_name[1] = {
    "THC_Core::bitmask"
};

} // namespace thc

namespace hrscc {

template<>
int CLaw<thc::ideal::Euler>::conserved_idx[5] = {0, 0, 0, 0, 0};
template<>
int CLaw<thc::ideal::Euler>::primitive_idx[5] = {0, 0, 0, 0, 0};
template<>
int CLaw<thc::ideal::Euler>::rhs_idx[5] = {0, 0, 0, 0, 0};
template<>
int CLaw<thc::ideal::Euler>::field_idx[19] =
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
template<>
int CLaw<thc::ideal::Euler>::bitmask_idx[1] = {0};
template<>
int CLaw<thc::ideal::Euler>::num_flux_idx[5*3] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};
template<>
CCTK_REAL CLaw<thc::ideal::Euler>::conserved_lbound[5] = {
    std::numeric_limits<CCTK_REAL>::quiet_NaN(),
    std::numeric_limits<CCTK_REAL>::quiet_NaN(),
    std::numeric_limits<CCTK_REAL>::quiet_NaN(),
    std::numeric_limits<CCTK_REAL>::quiet_NaN(),
    std::numeric_limits<CCTK_REAL>::quiet_NaN()
};

} // namespace hrscc

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Nuclear EoS
//////////////////////////////////////////////////////////////////////////////
namespace thc {

template<>
char const * EulerBase<nuclear::Euler>::conserved_name[6] = {
    "THC_Core::densxn",
    "THC_Core::densxp",
    "THC_Core::scon[0]",
    "THC_Core::scon[1]",
    "THC_Core::scon[2]",
    "THC_Core::tau"
};

template<>
char const * EulerBase<nuclear::Euler>::primitive_name[6] = {
    "HydroBase::rho",
    "THC_Core::zvec[0]",
    "THC_Core::zvec[1]",
    "THC_Core::zvec[2]",
    "HydroBase::eps",
    "HydroBase::Y_e"
};
template<>
char const * EulerBase<nuclear::Euler>::field_name[22] = {
    "ADMBase::alp",
    "ADMBase::betax",
    "ADMBase::betay",
    "ADMBase::betaz",
    "ADMBase::gxx",
    "ADMBase::gxy",
    "ADMBase::gxz",
    "ADMBase::gyy",
    "ADMBase::gyz",
    "ADMBase::gzz",
    "HydroBase::entropy",
    "HydroBase::press",
    "HydroBase::temperature",
    "HydroBase::vel[0]",
    "HydroBase::vel[1]",
    "HydroBase::vel[2]",
    "HydroBase::w_lorentz",
    "THC_Core::c2a_nbiter",
    "THC_Core::csound",
    "THC_Core::dens",
    "THC_Core::c2a_densgain",
    "THC_Core::volform"
};
template<>
char const * EulerBase<nuclear::Euler>::rhs_name[6] = {
    "THC_Core::rhs_densxn",
    "THC_Core::rhs_densxp",
    "THC_Core::rhs_scon[0]",
    "THC_Core::rhs_scon[1]",
    "THC_Core::rhs_scon[2]",
    "THC_Core::rhs_tau",
};
template<>
char const * EulerBase<nuclear::Euler>::bitmask_name[1] = {
    "THC_Core::bitmask"
};

} // namespace thc

namespace hrscc {

template<>
int CLaw<thc::nuclear::Euler>::conserved_idx[6] = {0, 0, 0, 0, 0, 0};
template<>
int CLaw<thc::nuclear::Euler>::primitive_idx[6] = {0, 0, 0, 0, 0, 0};
template<>
int CLaw<thc::nuclear::Euler>::rhs_idx[6] = {0, 0, 0, 0, 0, 0};
template<>
int CLaw<thc::nuclear::Euler>::field_idx[22] =
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
template<>
int CLaw<thc::nuclear::Euler>::bitmask_idx[1] = {0};
template<>
int CLaw<thc::nuclear::Euler>::num_flux_idx[6*3] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};
template<>
CCTK_REAL CLaw<thc::nuclear::Euler>::conserved_lbound[6] = {
    std::numeric_limits<CCTK_REAL>::quiet_NaN(),
    std::numeric_limits<CCTK_REAL>::quiet_NaN(),
    std::numeric_limits<CCTK_REAL>::quiet_NaN(),
    std::numeric_limits<CCTK_REAL>::quiet_NaN(),
    std::numeric_limits<CCTK_REAL>::quiet_NaN(),
    std::numeric_limits<CCTK_REAL>::quiet_NaN()
};

} // namespace hrscc

//////////////////////////////////////////////////////////////////////////////

CCTK_REAL thc::ultrarelativistic::Euler::gamma = 1.0/3.0;
