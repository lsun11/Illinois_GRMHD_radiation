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


#ifndef UTILS_FUNCTIONAL_HH
#define UTILS_FUNCTIONAL_HH

#include <functional>

namespace utils {

//! An abstract function
template<typename arg_t, typename res_t>
class unary_function: std::unary_function<arg_t, res_t> {
    public:
        virtual res_t operator()(arg_t) const = 0;
};

}

#endif
