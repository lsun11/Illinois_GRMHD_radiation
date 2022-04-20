//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2013, David Radice <dradice@caltech.edu>
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


#ifndef THC_CLAWS_BAROTROPIC_HH
#define THC_CLAWS_BAROTROPIC_HH

namespace thc {
namespace barotropic {

class Euler: public EulerBase<Euler> {
    public:
        typedef hrscc::CLaw<Euler> claw;

        inline void prim_to_all(
                hrscc::Observer<claw> & observer
                ) const {
#pragma omp critical
            CCTK_ERROR("barotropic::Euler::prim_to_all not implemented!");
        }

        inline void cons_to_all(
                hrscc::Observer<claw> & observer
                ) const {
#pragma omp critical
            CCTK_ERROR("barotropic::Euler::cons_to_all not implemented!");
        }

        template<hrscc::policy::direction_t dir>
        inline void fluxes(
                hrscc::Observer<claw> & observer
                ) const {
#pragma omp critical
            CCTK_ERROR("barotropic::Euler::fluxes not implemented!");
        }

        template<hrscc::policy::direction_t dir>
        inline void eigenvalues(
                hrscc::Observer<claw> & observer
                ) const {
#pragma omp critical
            CCTK_ERROR("barotropic::Euler::eigenvalues not implemented!");
        }

        template<hrscc::policy::direction_t dir>
        inline void eig(
                hrscc::Observer<claw> & observer
                ) const {
#pragma omp critical
            CCTK_ERROR("barotropic::Euler::eig not implemented!");
        }
};

} // namespace barotropic
} // namespace thc

#endif
