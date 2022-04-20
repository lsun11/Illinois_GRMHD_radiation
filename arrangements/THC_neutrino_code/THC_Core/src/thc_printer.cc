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


#include "thc_printer.hh"

#include <algorithm>
#include <cstring>
#include <sstream>

namespace thc {

std::string Printer::_M_prefix_info;
std::string Printer::_M_prefix_warn;
std::string Printer::_M_prefix_err;

int Printer::_M_n_info;
int Printer::_M_n_warn;
int Printer::_M_max_n_info;
int Printer::_M_max_n_warn;

#ifdef _OPENMP
omp_lock_t Printer::_M_lock;
#endif

char * Printer::_M_scratch;
int Printer::_M_scratch_siz;

void Printer::start(
        std::string const & prefix_info,
        std::string const & prefix_warn,
        std::string const & prefix_err,
        int max_n_info,
        int max_n_warn) {
#ifdef _OPENMP
    omp_init_lock(&_M_lock);
#endif
    _M_n_info     = 0;
    _M_n_warn     = 0;

    _M_prefix_info = prefix_info;
    _M_prefix_warn = prefix_warn;
    _M_prefix_err  = prefix_err;

    _M_max_n_info = max_n_info;
    _M_max_n_warn = max_n_warn;

    _M_scratch     = new char[BUFSIZ];
    _M_scratch_siz = BUFSIZ;
}

void Printer::print_info(std::string const & msg) {
    _M_n_info++;
#ifdef _OPENMP
    // This seems redundant, but avoids unnecessary calls to omp_set_lock
    if(_M_max_n_info > 0 && _M_n_info > _M_max_n_info) {
        return;
    }
    omp_set_lock(&_M_lock);
#endif
    if(_M_max_n_info < 0 || _M_n_info <= _M_max_n_info) {
        _M_print(stdout, _M_prefix_info, msg);
    }
#ifdef _OPENMP
    omp_unset_lock(&_M_lock);
#endif
}

void Printer::print_warn(std::string const & msg) {
    _M_n_warn++;
#ifdef _OPENMP
    // This seems redundant, but avoids unnecessary calls to omp_set_lock
    if(_M_max_n_warn > 0 && _M_n_warn > _M_max_n_warn) {
        return;
    }
    omp_set_lock(&_M_lock);
#endif
    if(_M_max_n_warn < 0 || _M_n_warn <= _M_max_n_warn) {
        _M_print(stderr, _M_prefix_warn, msg);
    }
#ifdef _OPENMP
    omp_unset_lock(&_M_lock);
#endif
}

void Printer::print_err(std::string const & msg) {
#ifdef _OPENMP
    omp_set_lock(&_M_lock);
#endif
    _M_print(stdout, _M_prefix_err, msg);
    _M_print(stderr, _M_prefix_err, msg);
#ifdef _OPENMP
    omp_unset_lock(&_M_lock);
#endif
}

void Printer::stop() {
#ifdef _OPENMP
    omp_destroy_lock(&_M_lock);
#endif
    if(_M_max_n_info >= 0 && _M_n_info > _M_max_n_info) {
        std::ostringstream ss;
        ss << "Additional " << _M_n_info - _M_max_n_info << " "
           << "messagges not shown";
        _M_print(stdout, _M_prefix_info, ss.str());
    }
    if(_M_max_n_warn >= 0 && _M_n_warn > _M_max_n_warn) {
        std::ostringstream ss;
        ss << "Additional " << _M_n_warn - _M_max_n_warn << " "
           << "messages not shown";
        _M_print(stderr, _M_prefix_warn, ss.str());
    }
    delete[] _M_scratch;
    _M_scratch_siz = 0;
}

void Printer::_M_print(
        FILE * outfile,
        std::string const & prefix,
        std::string const & msg) {
    // Allocates space in the scratch buffer if needed
    int psiz = prefix.size();
    int msiz = msg.size();
    int cnt  = std::count(msg.begin(), msg.end(), '\n');
    // We need to store:
    // 1. A copy of msg + '\0'
    // 2. A new version of msg where we add to all the newlines a
    //    copy of the prefix + '\n' + '\n' + '\0'
    int rsiz = 2*msiz + (cnt + 1)*psiz + 4;
    if(rsiz > _M_scratch_siz) {
        delete[] _M_scratch;
        _M_scratch = new char[rsiz];
        _M_scratch_siz = rsiz;
    }

    // Copy msg into the scratch area
    char * omsg = _M_scratch;
    std::strcpy(omsg, msg.c_str());

    // Assemble message
    char * nmsg = omsg + msiz + 1;
    *nmsg = '\0';
    char * p = std::strtok(omsg, "\n");
    while(p) {
        nmsg = std::strcat(nmsg, prefix.c_str());
        nmsg = std::strcat(nmsg, p);
        nmsg = std::strcat(nmsg, "\n");
        p = std::strtok(NULL, "\n");
    }
    nmsg = std::strcat(nmsg, "\n");

    fputs(nmsg, outfile);
    fflush(outfile);
}

} // namespace
