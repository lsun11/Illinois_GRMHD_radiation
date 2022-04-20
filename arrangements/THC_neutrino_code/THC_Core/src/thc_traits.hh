//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
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


#ifndef THC_TRAITS_HH
#define THC_TRAITS_HH

namespace thc {

namespace barotropic {
class Euler;
}

namespace ideal {
class Euler;
}

namespace nuclear {
class Euler;
}

namespace ultrarelativistic {
class Euler;
}

}

namespace hrscc {

template<>
class traits<thc::barotropic::Euler> {
    public:
        enum {nequations = 4};
        enum {nexternal = 20};
        enum {nbitmasks = 1};
        static bool const pure = false;
};

template<>
class traits<thc::ideal::Euler> {
    public:
        enum {nequations = 5};
        enum {nexternal = 19};
        enum {nbitmasks = 1};
        static bool const pure = false;
};

template<>
class traits<thc::nuclear::Euler> {
    public:
        enum {nequations = 6};
        enum {nexternal = 22};
        enum {nbitmasks = 1};
        static bool const pure = false;
};

template<>
class traits<thc::ultrarelativistic::Euler> {
    public:
        enum {nequations = 4};
        enum {nexternal = 4};
        enum {nbitmasks = 1};
        static bool const pure = false;
};

}

#endif
