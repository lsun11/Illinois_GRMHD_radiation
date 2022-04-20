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


#ifndef THC_CLAWS_HH
#define THC_CLAWS_HH

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "con2prim.h"
#include "eos_thermal.h"
#include "global_eos_thermal.h"
#include "hrscc.hh"
#include "utils.hh"

#include "thc_atmosphere.hh"
#include "thc_macro.hh"
#include "thc_printer.hh"
#include "thc_traits.hh"

#define THC_CLAWS_ASSERT(X)                                                   \
    if(!(X)) {                                                                \
        std::ostringstream ss;                                                \
        this->print_location(observer, ss);                                   \
        Printer::print_err(ss.str());                                         \
    }                                                                         \
    assert((X))

#define THC_CLAWS_CHECK(X)                                                    \
    THC_CLAWS_ASSERT(ISFINITE(X))

#define THC_CLAWS_CHECK_ALL                                                   \
    THC_CLAWS_ASSERT(this->check_all(observer))

namespace thc {

template<typename D>
class EulerBase: public hrscc::CLaw<D> {
    public:
        typedef hrscc::Observer<hrscc::CLaw<D> > Observer;
        enum {nequations = hrscc::CLaw<D>::nequations};
        enum {nexternal = hrscc::CLaw<D>::nexternal};
        enum {nbitmasks = hrscc::CLaw<D>::nbitmasks};

        static char const * conserved_name[nequations];
        static char const * primitive_name[nequations];
        static char const * rhs_name[nequations];
        static char const * field_name[nexternal];
        static char const * bitmask_name[nbitmasks];

        static void get_var_idx() {
            DECLARE_CCTK_PARAMETERS
            for(int i = 0; i < nequations; ++i) {
                hrscc::CLaw<D>::conserved_idx[i] =
                    utils::cctk::var_index(conserved_name[i]);
            }
            for(int i = 0; i < nequations; ++i) {
                hrscc::CLaw<D>::primitive_idx[i] =
                    utils::cctk::var_index(primitive_name[i]);
            }
            for(int i = 0; i < nequations; ++i) {
                hrscc::CLaw<D>::rhs_idx[i] =
                    utils::cctk::var_index(rhs_name[i]);
            }
            for(int i = 0; i < nexternal; ++i) {
                hrscc::CLaw<D>::field_idx[i] =
                    utils::cctk::var_index(field_name[i]);
            }
            for(int i = 0; i < nbitmasks; ++i) {
                hrscc::CLaw<D>::bitmask_idx[i] =
                    utils::cctk::var_index(bitmask_name[i]);
            }
            if(refluxing) {
                char vname[512];
                for(int i = 0; i < 3*nequations; ++i) {
                    std::snprintf(vname, 512, "Refluxing::flux[%d]", i);
                    hrscc::CLaw<D>::num_flux_idx[i] =
                        utils::cctk::var_index(vname);
                }
            }
        }

        bool check_cons(Observer const & observer) const {
            for(int i = 0; i < nequations; ++i) {
                if(!ISFINITE(observer.conserved[i])) {
                    return false;
                }
            }
            return true;
        }
        bool check_prim(Observer const & observer) const {
            for(int i = 0; i < nequations; ++i) {
                if(!ISFINITE(observer.primitive[i])) {
                    return false;
                }
            }
            return true;
        }
        bool check_fields(Observer const & observer) const {
            for(int i = 0; i < nexternal; ++i) {
                if(!ISFINITE(observer.field[i])) {
                    return false;
                }
            }
            return true;
        }
        template<hrscc::policy::direction_t dir>
        bool check_fluxes(Observer const & observer) const {
            for(int i = 0; i < nequations; ++i) {
                if(!ISFINITE(observer.flux[dir][i])) {
                    return false;
                }
            }
            return true;
        }
        bool check_eigenvalues(Observer const & observer) const {
            for(int i = 0; i < nequations; ++i) {
                if(!ISFINITE(observer.eigenvalue[i])) {
                    return false;
                }
            }
            return true;
        }
        bool check_eigenvectors(Observer const & observer) const {
            for(int i = 0; i < nequations; ++i)
            for(int j = 0; j < nequations; ++j) {
                if(!ISFINITE(observer.left_eigenvector[i][j])) {
                    return false;
                }
                if(!ISFINITE(observer.right_eigenvector[i][j])) {
                    return false;
                }
            }
            return true;
        }
        bool check_all(Observer const & observer) const {
            if(!this->check_cons(observer)) {
                return false;
            }
            if(!this->check_prim(observer)) {
                return false;
            }
            if(!this->check_fields(observer)) {
                return false;
            }
            if(!this->check_fluxes<hrscc::policy::x>(observer)) {
                return false;
            }
            if(!this->check_fluxes<hrscc::policy::y>(observer)) {
                return false;
            }
            if(!this->check_fluxes<hrscc::policy::z>(observer)) {
                return false;
            }
            if(!this->check_eigenvalues(observer)) {
                return false;
            }
            if(!this->check_eigenvectors(observer)) {
                return false;
            }
            return true;
        }

        void print_location(
                hrscc::Observer<hrscc::CLaw<D> > & observer,
                std::ostream & ss) const {
            ss << "Iteration = "
               << observer.cctkGH->cctk_iteration
               << std::endl;
            ss << "Reflevel = "
               << ilogb(observer.cctkGH->cctk_levfac[0])
               << std::endl;
            ss << "(i, j, k) = ("
               << observer.i
               << ", "
               << observer.j
               << ", "
               << observer.k
               << ")\n";
            ss << "shift = ("
               << observer.shift[0]
               << " "
               << observer.shift[1]
               << " "
               << observer.shift[2]
               << ")\n";
            ss << "(x, y, z) = ("
               << observer.x
               << ", "
               << observer.y
               << ", "
               << observer.z
               << ")\n";
            for(int i = 0; i < nequations; ++i) {
                ss << EulerBase<D>::conserved_name[i]
                   << " = "
                   << observer.conserved[i]
                   << std::endl;
            }
            ss << std::endl;
            for(int i = 0; i < nequations; ++i) {
                ss << EulerBase<D>::primitive_name[i]
                   << " = "
                   << observer.primitive[i]
                   << std::endl;
            }
            ss << std::endl;
            for(int i = 0; i < nexternal; ++i) {
                ss << EulerBase<D>::field_name[i]
                   << " = "
                   << observer.field[i]
                   << std::endl;
            }
            ss << "THC_Core::bitmask = "
               << observer.bitmask[0]
               << std::endl;
        }
};

} // namespace thc

#include "thc_claws/barotropic.hh"
#include "thc_claws/ideal.hh"
#include "thc_claws/nuclear.hh"
#include "thc_claws/ultrarelativistic.hh"

#endif
