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


#ifndef UTILS_FILTER_HH
#define UTILS_FILTER_HH

//! Machine precision: \f$ \log \epsilon \f$
#define UTILS_FILTER_LOGEPS     35

#include <cmath>
#include <limits>

#include "cctk.h"

#include "utils_pow.hh"

namespace utils {
//! Linear filters namespace
namespace filter {

class generic {
    public:
        virtual CCTK_REAL operator()(CCTK_REAL) const = 0;
};

//! Dummy, do nothing filter
class dummy: public generic {
    public:
        //! Returns 1
        inline CCTK_REAL operator()(CCTK_REAL) const {
            return 1.0;
        }
};

//! First-order linear filter
class linear: public generic {
    public:
        inline CCTK_REAL operator()(CCTK_REAL eta) const {
            return 1.0 - eta;
        }
};

//! Lanczos filter: \f$ \sigma(\eta) = \frac{\sin\pi\eta}{\pi\eta} \f$
class lanczos: public generic {
    public:
        //! Evaluate the filter
        inline CCTK_REAL operator()(CCTK_REAL eta) const {
            if(eta < 1e-16) {
                return 1.0;
            }
            return std::sin(M_PI*eta)/(M_PI*eta);
        }
};

//! Modified Lanczos filter: \f$ \sigma(\eta) = \frac{\sin \eta}{\eta} \f$
class mlanczos: public generic {
    public:
        //! Evaluate the filter
        inline CCTK_REAL operator()(CCTK_REAL eta) const {
            if(eta < 1e-16) {
                return 1.0;
            }
            return std::sin(eta)/eta;
        }
};

//! Raised cosine filter: \f$ \sigma(\eta) = 1/2 (1+\cos \pi \eta) \f$
class raised_cosine: public generic {
    public:
        //! Evaluate the filter
        inline CCTK_REAL operator()(CCTK_REAL eta) const {
            return 0.5*(1.0+std::cos(M_PI*eta));
        }
};

//! Sharped cosine filter
class sharped_cosine: public generic {
    public:
        //! Evaluate the filter
        inline CCTK_REAL operator()(CCTK_REAL eta) const {
            CCTK_REAL const s = 0.5*(1.0 + std::cos(M_PI*eta));
            return s*s*s*s*(35 - 84*s + 70*s*s - 20*s*s*s);
        }
};

//! logerfc filter
/*!
 *  It is an approximate Vandeven filter, see:
 *
 *  J. B. Boyd, The Erfc-Log Filter and the Asymptotics of the Euler
 *  and Vandeven Sequence Accelerations,In Proceedings
 *  of the Third International Conference on Spectral and
 *  High Order Methods (1996) 267-276
 *
 *  equation (15).
 */
class erfclog: public generic {
    public:
        //! Constructor
        erfclog(
                unsigned p = 8                //! [in] filter order
                );
        //! Evaluate the filter
        CCTK_REAL operator()(CCTK_REAL eta) const;
        //! Filter order
        unsigned p;
};


//! Exponential filter: \f$ \sigma(\eta) = \exp(-\eta^p) \f$
/*!
 *  The strength of this filter should be in the form
 *  \f$ \frac{C_p (N+1) \Delta t}{h} \f$, see Meister et al. (2009),
 *  On Spectral Filtering for Discontinuous Galerkin
 *  Methods on Unstructured Triangular Grids.
 */
class exp: public generic {
    public:
        //! Constructor
        exp(
                unsigned p = 4
                );
        //! Filter evaluation
        inline CCTK_REAL operator()(CCTK_REAL eta) const {
                return std::exp(-UTILS_FILTER_LOGEPS*std::pow(eta,
                            static_cast<int>(p)));
        }
        //! Filter order
        unsigned p;
};

//! Cut-off filter
class cut_off: public generic {
    public:
        //! Constructor
        cut_off(CCTK_REAL etac = 0.1);
        //! Filter evaluation
        inline CCTK_REAL operator()(CCTK_REAL eta) const {
            return eta < etac ? 1.0 : 0.0;
        }

        //! Cut-off "frequency"
        CCTK_REAL etac;
};

//! Exponential cut-off filter
/*!
 *  The strength of this filter should be in the form
 *  \f$ \frac{C_p (N+1) \Delta t}{h} \f$, see Meister et al. (2009),
 *  On Spectral Filtering for Discontinuous Galerkin
 *  Methods on Unstructured Triangular Grids.
 */
class exp_cut_off: public generic {
    public:
        //! Constructor
        exp_cut_off(
                unsigned p = 4,
                CCTK_REAL etac = 0.5
                );
        //! Filter evaluation
        inline CCTK_REAL operator()(CCTK_REAL eta) const {
            return (eta>etac)?
                (std::exp(-UTILS_FILTER_LOGEPS*std::pow(eta-etac,
                        static_cast<int>(p)))):
                (1.0);
        }
        //! Filter order
        unsigned p;
        //! Cut-off "frequency"
        CCTK_REAL etac;
};

//! Spherical harmonics spline filter
/*!
 *  This filter only makes sense for a spherical harmonics expansion
 */
class sspline: public generic {
    public:
        //! Filter evaluation
        inline CCTK_REAL operator()(CCTK_REAL eta) const {
            return 1.0 / (1.0 + utils::pow<4>(eta));
        }
};

}
}

#endif
