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


#ifndef UTILS_MACRO_H
#define UTILS_MACRO_H

#include <stdlib.h>

#include <cctk_Loop.h>

//! custom loop macro
#define UTILS_LOOP3(NAME,I,SI,EI,J,SJ,EJ,K,SK,EK)                              \
    _Pragma("omp for collapse(3)")                                             \
    for(int I = SI; I < EI; ++I)                                               \
    for(int J = SJ; J < EJ; ++J)                                               \
    for(int K = SK; K < EK; ++K)

//! custom loop macro
#define UTILS_ENDLOOP3(NAME)

//! layout memory size
#define UTILS_GFSIZE(CGH)                                                      \
    (CGH->cctk_ash[0]*CGH->cctk_ash[1]*CGH->cctk_ash[2])

//! find a unified index using the given stride
#define UTILS_GFINDEX3D(STRIDE,I,J,K)                                          \
    ((I)*(STRIDE[0]) + (J)*(STRIDE[1]) + (K)*(STRIDE[2]))

//! make a flag
#define UTILS_BITMASK_MAKE_FLAG(INDEX)                                         \
    (1 << INDEX)

//! set a flag in a bitmask
#define UTILS_BITMASK_SET_FLAG(BITMASK,FLAG)                                   \
    (BITMASK) |= (FLAG)

//! set all the flags in a bitmask
#define UTILS_BITMASK_SET_ALL_FLAGS(BITMASK)                                   \
    (BITMASK = ~(0))

//! unset a flag in a bitmask
#define UTILS_BITMASK_UNSET_FLAG(BITMASK,FLAG)                                 \
    (BITMASK) &= ~(FLAG)

//! unset all flags in a bitmask
#define UTILS_BITMASK_UNSET_ALL_FLAGS(BITMASK)                                 \
    (BITMASK) = 0

//! check if a flag is set in a bitmask
#define UTILS_BITMASK_CHECK_FLAG(BITMASK,FLAG)                                 \
    (((BITMASK) & (FLAG)) == (FLAG))

//! extract a random float in the given interval
#define UTILS_RANDOM(MIN,MAX)                                                  \
    ((MIN) + (CCTK_REAL)rand()/((CCTK_REAL)RAND_MAX) * ((MAX) - (MIN)))

//! sign function
#define UTILS_SIGN(X)                                                          \
    ((X < 0) ? (-1) : (1))

#endif
