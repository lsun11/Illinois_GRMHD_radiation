//  CPPUtils: C++ utilities for Cactus
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


#ifndef UTILS_CACTUS_HH
#define UTILS_CACTUS_HH

#include "cctk.h"

namespace utils {
namespace cctk {

//! Wrapper around CCTK_VarIndex
int var_index(char const * name);

//! Wrapper around CCTK_VarDataPtr
void * var_data_ptr(
        cGH const * cctkGH,
        int const timelevel,
        char const * name);

//! Wrapper around CCTK_VarDataPtrI
void * var_data_ptr_i(
        cGH const * cctkGH,
        int const timelevel,
        int const vidx);

} // namespace cctk
} // namespace utils

#endif
