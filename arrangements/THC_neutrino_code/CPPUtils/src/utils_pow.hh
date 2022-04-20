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


#ifndef UTILS_POW_HH
#define UTILS_POW_HH

namespace utils {

#define UTILS_POW(T)                                                           \
/*! Computes x**N using template meta-programming */                           \
template<int N>                                                                \
inline T pow(T x) {                                                            \
    return T(x*pow<N-1>(x));                                                   \
}                                                                              \
                                                                               \
/*! Specialized version for x**1 */                                            \
template<>                                                                     \
inline T pow<1>(T x) {                                                         \
    return T(x);                                                               \
}                                                                              \
                                                                               \
/*! Specialized version for x**0 */                                            \
template<>                                                                     \
inline T pow<0>(T) {                                                           \
    return T(1);                                                               \
}

UTILS_POW(int)
UTILS_POW(unsigned)
UTILS_POW(float)
UTILS_POW(double)

#undef UTILS_POW

} // namespace

#endif
