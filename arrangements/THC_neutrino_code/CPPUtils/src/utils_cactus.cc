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


#include <cassert>

#include "cctk.h"
#include "utils_cactus.hh"

namespace utils {
namespace cctk {

int var_index(char const * name) {
    int i = CCTK_VarIndex(name);
    if(i < 0) {
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Could not find the index of \"%s\"", name);
    }
    return i;
}

void * var_data_ptr(
        cGH const * cctkGH,
        int const timelevel,
        char const * name) {
    return var_data_ptr_i(cctkGH, timelevel, var_index(name));
}

void * var_data_ptr_i(
        cGH const * cctkGH,
        int const timelevel,
        int const vidx) {
    assert(vidx >= 0);
    void * ptr = CCTK_VarDataPtrI(cctkGH, timelevel, vidx);
    if(NULL == ptr) {
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Could not get data for variable \"%d\"", vidx);
    }
    return ptr;
}

} // namespace cctk
} // namespace utils
