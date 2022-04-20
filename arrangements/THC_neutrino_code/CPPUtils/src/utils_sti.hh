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


#ifndef UTILS_STI_HH
#define UTILS_STI_HH

namespace utils {

//! construct for static type identification
template<typename A, typename B>
class is_same_type {
    public:
        static bool const value = false;
};

//! specialization of \e is_same_type<A,B> to the case in which A=B
template<typename A>
class is_same_type<A,A> {
    public:
        static bool const value = true;
};

//! returns true if \e b is of type \e A
template<typename A, typename B>
inline bool is_a(B const & b) {
    return is_same_type<A,B>::value;
}

} // namespace

#endif
