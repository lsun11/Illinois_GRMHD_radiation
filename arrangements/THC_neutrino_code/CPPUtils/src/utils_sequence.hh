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


#ifndef UTILS_SEQUENCE_HH
#define UTILS_SEQUENCE_HH

namespace utils {

//! predicates on small fixed size sequences
template<typename T, int N>
class sequence {
    public:
        //! returns true if the sequence is strictly increasing
        static inline bool is_strictly_increasing(T const x[N]) {
            return x[0] < x[1] &&
                sequence<T, N-1>::is_strictly_increasing(&x[1]);
        }

        //! returns true if the sequence is increasing
        static inline bool is_increasing(T const x[N]) {
            return x[0] <= x[1] && sequence<T, N-1>::is_increasing(&x[1]);
        }

        //! returns true if the sequence is strictly increasing
        static inline bool is_strictly_decreasing(T const x[N]) {
            return x[0] > x[1] &&
                sequence<T, N-1>::is_strictly_decreasing(&x[1]);
        }

        //! returns true if the sequence is increasing
        static inline bool is_decreasing(T const x[N]) {
            return x[0] >= x[1] && sequence<T, N-1>::is_decreasing(&x[1]);
        }

        //! returns true if the sequence is strictly monotonic
        static inline bool is_strictly_monotonic(T const x[N]) {
            return sequence<T, N>::is_strictly_increasing(x) ||
                sequence<T, N>::is_strictly_decreasing(x);
        }

        //! returns true if the sequence is monotonic
        static inline bool is_monotonic(T const x[N]) {
            return sequence<T, N>::is_increasing(x) ||
                sequence<T, N>::is_decreasing(x);
        }
};

//! specialization of \e sequence<T,N> to N=2
template<typename T>
class sequence<T, 2> {
    public:
        //! returns true if the sequence is strictly increasing
        static inline bool is_strictly_increasing(T const x[2]) {
            return x[0] < x[1];
        }

        //! returns true if the sequence is increasing
        static inline bool is_increasing(T const x[2]) {
            return x[0] <= x[1];
        }

        //! returns true if the sequence is strictly increasing
        static inline bool is_strictly_decreasing(T const x[2]) {
            return x[0] > x[1];
        }

        //! returns true if the sequence is increasing
        static inline bool is_decreasing(T const x[2]) {
            return x[0] >= x[1];
        }

        //! returns true if the sequence is strictly monotonic
        static inline bool is_strictly_monotonic(T const x[2]) {
            return x[0] != x[1];
        }

        //! returns true if the sequence is monotonic
        static inline bool is_monotonic(T const x[2]) {
            return true;
        }
};

} // namespace

#endif
