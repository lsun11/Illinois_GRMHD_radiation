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


#ifndef THC_PRINTER_HH
#define THC_PRINTER_HH

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <string>

namespace thc {

//! Custom printing facility
/*!
 *  Thread-safe simple output
 */
class Printer {
    public:
        //! Initialize the printer
        /*!
         *  \param prefix_info Prefix for info messages
         *  \param prefix_warn Prefix for warn messages
         *  \param prefix_err  Prefix for error messagges
         *  \param max_n_info  Maximum number of info messages to print
         *                     (ignored if negative)
         *  \param max_n_warn  Maximum number of warn messages to print
         *                     (ignored if negative)
         *
         *  This should be called from outside a parallel block!
         *
         *  NOTE: There is no limit on error messagges
         */
        static void start(
                std::string const & prefix_info = "",
                std::string const & prefix_warn = "",
                std::string const & prefix_err  = "",
                int max_n_info = -1,
                int max_n_warn = -1);
        //! Prints an info message
        static void print_info(std::string const & msg);
        //! Prints a warning
        static void print_warn(std::string const & msg);
        //! Prints an error
        static void print_err(std::string const & msg);
        //! Stops the printer and print number of discarded messages
        /*!
         *  This should be called from outside a parallel block!
         */
        static void stop();
    private:
        static void _M_print(
                FILE * outfile,
                std::string const & prefix,
                std::string const & msg);
    private:
        static std::string _M_prefix_info;
        static std::string _M_prefix_warn;
        static std::string _M_prefix_err;

        static int _M_n_info;
        static int _M_n_warn;
        static int _M_max_n_info;
        static int _M_max_n_warn;

#ifdef _OPENMP
        static omp_lock_t _M_lock;
#endif

        static char * _M_scratch;
        static int _M_scratch_siz;
};

} // namespace

#endif
