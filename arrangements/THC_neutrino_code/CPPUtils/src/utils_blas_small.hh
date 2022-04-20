//  CPPUtils: C++ utilities for Cactus
//  Copyright (C) 2011, David Radice <david.radice@aei.mpg.de>
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


#ifndef UTILS_BITS_BLAS_SMALL_HH
#define UTILS_BITS_BLAS_SMALL_HH

#include <boost/preprocessor/repetition/repeat.hpp>

#include <utils_config.hh>

namespace utils {
namespace details {

//! Choose between small-size-optimized BLAS and standard BLAS implementations
template<int N>
class blas_small {
    public:
        static bool const value = false;
};

#define UTILS_DECLARE_AS_SMALL(z, N, unused)                                   \
    template<>                                                                 \
    class blas_small<N> {                                                      \
        public:                                                                \
            static bool const value = true;                                    \
    };

BOOST_PP_REPEAT(UTILS_CONFIG_BLAS_SMALL, UTILS_DECLARE_AS_SMALL, ~)

#undef UTILS_DECLARE_AS_SMALL

} // namespace details
} // namespace utils

#endif
