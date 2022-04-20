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


#include <algorithm>
#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>

#include <cctk.h>

#include "utils_filter.hh"

#include "utils_pow.hh"

namespace utils {
namespace filter {

erfclog::erfclog(unsigned p): p(p) {}

CCTK_REAL erfclog::operator()(CCTK_REAL eta) const {
	if(eta < std::numeric_limits<CCTK_REAL>::epsilon()) {
		return 1.0;
	}
	if(eta > 1 - std::numeric_limits<CCTK_REAL>::epsilon()) {
		return 0.0;
	}

	CCTK_REAL tmp;
	CCTK_REAL bar_eta = eta - 0.5;

    CCTK_REAL const den = 4.0*utils::pow<2>(bar_eta);
    if(den < std::numeric_limits<CCTK_REAL>::epsilon()) {
        return 0.5;
    }

	tmp = -std::log(1.0-4.0*utils::pow<2>(bar_eta));
	tmp /= den;
	tmp = std::sqrt(tmp);
	tmp *= 2.0*std::sqrt(CCTK_REAL(p))*bar_eta;

	return 0.5*erfc(tmp);
}

exp::exp(unsigned p): p(p) {}

cut_off::cut_off(CCTK_REAL etac): etac(etac) {}

exp_cut_off::exp_cut_off(unsigned p,CCTK_REAL etac): p(p), etac(etac) {}

}
}
