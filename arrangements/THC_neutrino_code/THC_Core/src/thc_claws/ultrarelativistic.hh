//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2013, David Radice <david.radice@aei.mpg.de>
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


#ifndef THC_CLAWS_ULTRARELATIVISTIC_HH
#define THC_CLAWS_ULTRARELATIVISTIC_HH

namespace thc {

namespace ultrarelativistic {

class Euler: public EulerBase<Euler> {
    public:
        typedef hrscc::CLaw<Euler> claw;

        inline void prim_to_all(
                hrscc::Observer<claw> & observer
                ) const {
            ULR_MAKE_ALIAS(hrscc::policy::x, observer);

            // By setting \eps = 0 we can use the classical expression for the
            // enthalpy, h = 1 + eps + p/rho_b, even in the ultrarelativistic
            // case using rho instead of rho_b
            eps = 0;
            csound = std::sqrt(gamma - 1);

            w_lorentz = 1.0 / std::sqrt(1 - SQ(velx) - SQ(vely) - SQ(velz));
            press = (gamma - 1)*rho;

            tau   = (rho + press)*SQ(w_lorentz) - press;
            sconx = (rho + press)*SQ(w_lorentz)*velx;
            scony = (rho + press)*SQ(w_lorentz)*vely;
            sconz = (rho + press)*SQ(w_lorentz)*velz;
        }

        inline void set_to_atmosphere(
                hrscc::Observer<claw> & observer
                ) const {
            ULR_MAKE_ALIAS(hrscc::policy::x, observer);
            DECLARE_CCTK_PARAMETERS

            if(atmosphere::dynamic) {
                rho = atmosphere::rho;

                CCTK_REAL const chix = sconx / (gamma*rho);
                CCTK_REAL const chiy = scony / (gamma*rho);
                CCTK_REAL const chiz = sconz / (gamma*rho);
                CCTK_REAL const chi2 = SQ(chix) + SQ(chiy) + SQ(chiz);

                if(chi2 > std::numeric_limits<CCTK_REAL>::epsilon()) {
                    velx = chix/(2*chi2)*(std::sqrt(1 + 4*chi2) - 1);
                    vely = chiy/(2*chi2)*(std::sqrt(1 + 4*chi2) - 1);
                    velz = chiz/(2*chi2)*(std::sqrt(1 + 4*chi2) - 1);
                }
                else {
                    velx = 0;
                    vely = 0;
                    velz = 0;
                }

                w_lorentz = 1/std::sqrt(1 - SQ(velx) - SQ(vely) - SQ(velz));
            }
            else {
                rho  = atmosphere::rho;
                velx = 0;
                vely = 0;
                velz = 0;
            }

            this->prim_to_all(observer);
        }

        // From Hawke 2001
        inline void cons_to_all(
                hrscc::Observer<claw> & observer
                ) const {
            ULR_MAKE_ALIAS(hrscc::policy::x, observer);
            DECLARE_CCTK_PARAMETERS
            std::ostringstream ss;

            if(UTILS_BITMASK_CHECK_FLAG(bitmask, THC_FLAG_IS_INSIDE_BH)) {
                UTILS_BITMASK_SET_FLAG(bitmask, THC_FLAG_ATMOSPHERE);
                this->set_to_atmosphere(observer);
                return;
            }

            if(c2a_fix_nans) {
                if(!ISFINITE(sconx) || !ISFINITE(scony) || !ISFINITE(sconz) ||
                   !ISFINITE(tau)) {
                    UTILS_BITMASK_SET_FLAG(bitmask, THC_FLAG_NOT_FINITE);
                    this->set_to_atmosphere(observer);
                    return;
                }
                else {
                    UTILS_BITMASK_UNSET_FLAG(bitmask, THC_FLAG_NOT_FINITE);
                }
            }

            // By setting \eps = 0 we can use the classical expression for the
            // enthalpy, h = 1 + eps + p/rho_b, even in the ultrarelativistic
            // case using rho instead of rho_b
            eps = 0;
            csound = std::sqrt(gamma - 1);

            if(tau < atmosphere::rho*(1 + atmosphere::tol)) {
                UTILS_BITMASK_SET_FLAG(bitmask, THC_FLAG_ATMOSPHERE);
                this->set_to_atmosphere(observer);
                return;
            }

            CCTK_REAL const beta = 0.5 - 0.25*gamma;
            CCTK_REAL const S2 = SQ(sconx) + SQ(scony) + SQ(sconz);

            CCTK_REAL const Delta = 4*SQ(beta)*SQ(tau) +
                        (gamma - 1)*(SQ(tau) - S2);

            if(Delta > 0) {
                rho = 1/(gamma - 1)*(-2*beta*tau + std::sqrt(Delta));
            }
            else {
                ss << "ConsToAll failure\n";
                this->print_location(observer, ss);
                if(c2a_kill_on_failure) {
                    Printer::print_err(ss.str());
#pragma omp critical
                    CCTK_ERROR("Aborting the run");
                }
                else {
                    ss << "Setting the point to atmosphere\n";
                    UTILS_BITMASK_SET_FLAG(bitmask, THC_FLAG_ATMOSPHERE);
                    this->set_to_atmosphere(observer);
                    Printer::print_err(ss.str());
                    return;
                }
            }

            CCTK_REAL const chix = sconx / (gamma*rho);
            CCTK_REAL const chiy = scony / (gamma*rho);
            CCTK_REAL const chiz = sconz / (gamma*rho);
            CCTK_REAL const chi2 = SQ(chix) + SQ(chiy) + SQ(chiz);

            if(chi2 > std::numeric_limits<CCTK_REAL>::epsilon()) {
                velx = chix/(2*chi2)*(std::sqrt(1 + 4*chi2) - 1);
                vely = chiy/(2*chi2)*(std::sqrt(1 + 4*chi2) - 1);
                velz = chiz/(2*chi2)*(std::sqrt(1 + 4*chi2) - 1);
            }
            else {
                velx = 0;
                vely = 0;
                velz = 0;
            }

            w_lorentz = 1/std::sqrt(1 - SQ(velx) - SQ(vely) - SQ(velz));
            press = (gamma - 1)*rho;

            if(rho < atmosphere::rho*(1 + atmosphere::tol)) {
                UTILS_BITMASK_SET_FLAG(bitmask, THC_FLAG_ATMOSPHERE);
                this->set_to_atmosphere(observer);
            }
        }

        template<hrscc::policy::direction_t dir>
        inline void fluxes(
                hrscc::Observer<claw> & observer
                ) const {
            ULR_MAKE_ALIAS(dir, observer);

            flux_tau   = sconx;
            flux_sconx = sconx*velx + press;
            flux_scony = scony*velx;
            flux_sconz = sconz*velx;
        }

        // From Hawke 2001
        template<hrscc::policy::direction_t dir>
        inline void eigenvalues(
                hrscc::Observer<claw> & observer
                ) const {
            ULR_MAKE_ALIAS(dir, observer);

            CCTK_REAL const v2 = SQ(velx) + SQ(vely) + SQ(velz);

            p_eigenvalue[0] = 1/(1 - SQ(velx)*SQ(csound))*(velx*
                    (1 - SQ(csound)) - csound*std::sqrt((1 - v2)*
                        (1 - SQ(velx) - SQ(csound)*(v2 - SQ(velx)))));
            p_eigenvalue[1] = velx;
            p_eigenvalue[2] = velx;
            p_eigenvalue[3] = 1/(1 - SQ(velx)*SQ(csound))*(velx*
                    (1 - SQ(csound)) + csound*std::sqrt((1 - v2)*
                        (1 - SQ(velx) - SQ(csound)*(v2 - SQ(velx)))));
        }

        // From Hawke 2001
        // The second component of l_+ and l_- has been recomputed with
        // Mathematica since it appears to be wrong in the reference above
        template<hrscc::policy::direction_t dir>
        inline void eig(
                hrscc::Observer<claw> & observer
                ) const {
            ULR_MAKE_ALIAS(dir, observer);
            this->eigenvalues<dir>(observer);

            ///////////////////////////////////////////////////////////////////
            // Some usefull constants
            ///////////////////////////////////////////////////////////////////
            CCTK_REAL const v2 = SQ(velx) + SQ(vely) + SQ(velz);
            CCTK_REAL const lambda_m = p_eigenvalue[0];
            CCTK_REAL const lambda_p = p_eigenvalue[3];
            CCTK_REAL const Delta = CB(1 - SQ(velx))*(gamma - 1)*(lambda_m -
                    lambda_p)*(1 - v2)*(1 - (gamma - 1)*v2);

            ///////////////////////////////////////////////////////////////////
            // Left eigenvectors
            ///////////////////////////////////////////////////////////////////
            L_eigenvector[0][0] = - 1/Delta*((1 - gamma)*lambda_p*CB(velx)*v2 +
                    3*(gamma - 1)*velx*(SQ(vely) + SQ(velz))*(velx - lambda_p) +
                    gamma*(1 - SQ(velx)) + gamma*(SQ(vely) + SQ(velz)) +
                    gamma*lambda_p*velx*(3*SQ(velx) - 2) +
                    lambda_p*velx*(3 - 4*SQ(velx)) + SQ(velx)*SQ(velx) -
                    (SQ(vely) + SQ(velz)) - 1);
            L_eigenvector[0][1] = - 1/Delta*(gamma*(SQ(velx)*(1 - SQ(velx))*
                    (lambda_p - velx) + (SQ(vely) + SQ(velz))*(3*lambda_p -
                    velx)*SQ(velx) + (SQ(vely) + SQ(velz))*(lambda_p -
                    3*velx)) + (SQ(velx)*(SQ(velx) - 3*(SQ(vely) + SQ(velz))) -
                    (SQ(vely) + SQ(velz)) - 1)*(lambda_p - velx) +
                    2*velx*(1 - SQ(velx))*(SQ(vely) + SQ(velz)));
            //L_eigenvector[0][1] = - 1/Delta*1/(1 - v2)*(velx*((1 + v2)*
            //        (1 - v2*(gamma - 1)) + 2*(1 - v2*(gamma - 1))*
            //        SQ(velx)*SQ(velx) + 2*(2*v2*(gamma - 1) - gamma)*
            //        (SQ(vely) + SQ(velz)) + SQ(velx)*((3 + v2)*(-1 +
            //        v2*(gamma - 1)) - 2*(gamma - 2)*SQ(vely) - 2*(gamma - 2)*
            //        SQ(velz))) + lambda_p*(2*(-1 + v2*(gamma - 1))*
            //        SQ(velx)*SQ(velx) + (-1 + v2*(gamma - 1))*
            //        (1 + v2 - 2*SQ(vely) - 2*SQ(velz)) + SQ(velx)*(-
            //        (3 + v2)*(-1 + v2*(gamma - 1)) + 2*(-3 + v2 + 2*gamma -
            //        v2*gamma)*SQ(vely) + 2*(-3 + v2 + 2*gamma - v2*gamma)*
            //        SQ(velz))));
            L_eigenvector[0][2] = - 1/Delta*(-2*(gamma - 1)*vely*(1 - SQ(velx))*
                    (1 - 2*lambda_p*velx + SQ(velx)));
            L_eigenvector[0][3] = - 1/Delta*(-2*(gamma - 1)*velz*(1 - SQ(velx))*
                    (1 - 2*lambda_p*velx + SQ(velx)));

            L_eigenvector[1][0] = 1/(1 - v2)*((-vely*(1 + v2 - 2*SQ(velx)))/
                    (1 - SQ(velx)));
            L_eigenvector[1][1] = 1/(1 - v2)*((velx*vely*(1 + v2 - 2*SQ(velx)))/
                    (1 - SQ(velx)));
            L_eigenvector[1][2] = 1/(1 - v2)*(1 - v2 + 2*SQ(vely));
            L_eigenvector[1][3] = 2/(1 - v2)*vely*velz;

            L_eigenvector[2][0] = 1/(1 - v2)*((-velz*(1 + v2 - 2*SQ(velx)))/
                    (1 - SQ(velx)));
            L_eigenvector[2][1] = 1/(1 - v2)*((velx*velz*(1 + v2 - 2*SQ(velx)))/
                    (1 - SQ(velx)));
            L_eigenvector[2][2] = 2/(1 - v2)*vely*velz;
            L_eigenvector[2][3] = 1/(1 - v2)*(1 - v2 + 2*SQ(velz));

            L_eigenvector[3][0] = + 1/Delta*((1 - gamma)*lambda_m*CB(velx)*v2 +
                    3*(gamma - 1)*velx*(SQ(vely) + SQ(velz))*(velx - lambda_m) +
                    gamma*(1 - SQ(velx)) + gamma*(SQ(vely) + SQ(velz)) +
                    gamma*lambda_m*velx*(3*SQ(velx) - 2) +
                    lambda_m*velx*(3 - 4*SQ(velx)) + SQ(velx)*SQ(velx) -
                    (SQ(vely) + SQ(velz)) - 1);
            L_eigenvector[3][1] = + 1/Delta*(gamma*(SQ(velx)*(1 - SQ(velx))*
                    (lambda_m - velx) + (SQ(vely) + SQ(velz))*(3*lambda_m -
                    velx)*SQ(velx) + (SQ(vely) + SQ(velz))*(lambda_m -
                    3*velx)) + (SQ(velx)*(SQ(velx) - 3*(SQ(vely) + SQ(velz))) -
                    (SQ(vely) + SQ(velz)) - 1)*(lambda_m - velx) +
                    2*velx*(1 - SQ(velx))*(SQ(vely) + SQ(velz)));
            //L_eigenvector[3][1] = + 1/Delta*1/(1 - v2)*(velx*((1 + v2)*
            //        (1 - v2*(gamma - 1)) + 2*(1 - v2*(gamma - 1))*
            //        SQ(velx)*SQ(velx) + 2*(2*v2*(gamma - 1) - gamma)*
            //        (SQ(vely) + SQ(velz)) + SQ(velx)*((3 + v2)*(-1 +
            //        v2*(gamma - 1)) - 2*(gamma - 2)*SQ(vely) - 2*(gamma - 2)*
            //        SQ(velz))) + lambda_m*(2*(-1 + v2*(gamma - 1))*
            //        SQ(velx)*SQ(velx) + (-1 + v2*(gamma - 1))*
            //        (1 + v2 - 2*SQ(vely) - 2*SQ(velz)) + SQ(velx)*(-
            //        (3 + v2)*(-1 + v2*(gamma - 1)) + 2*(-3 + v2 + 2*gamma -
            //        v2*gamma)*SQ(vely) + 2*(-3 + v2 + 2*gamma - v2*gamma)*
            //        SQ(velz))));
            L_eigenvector[3][2] = + 1/Delta*(-2*(gamma - 1)*vely*(1 - SQ(velx))*
                    (1 - 2*lambda_m*velx + SQ(velx)));
            L_eigenvector[3][3] = + 1/Delta*(-2*(gamma - 1)*velz*(1 - SQ(velx))*
                    (1 - 2*lambda_m*velx + SQ(velx)));

            ///////////////////////////////////////////////////////////////////
            // Right eigenvectors
            ///////////////////////////////////////////////////////////////////
            R_eigenvector[0][0] = (1 - SQ(velx))*(velx - lambda_m)*
                    (1 - (gamma - 1)*v2);
            R_eigenvector[1][0] = (1 - SQ(velx))*((lambda_m*velx - 1)*
                    (gamma - 1)*(1 - v2) + velx*(2 - gamma)*(velx - lambda_m));
            R_eigenvector[2][0] = vely*((gamma - 1)*(lambda_m*(1 + SQ(velx))*
                    v2 - lambda_m*SQ(velx) - 2*velx*v2) + (2 - gamma)*SQ(velx)*
                    (lambda_m - velx) + gamma*velx - lambda_m);
            R_eigenvector[3][0] = velz*((gamma - 1)*(lambda_m*(1 + SQ(velx))*
                    v2 - lambda_m*SQ(velx) - 2*velx*v2) + (2 - gamma)*SQ(velx)*
                    (lambda_m - velx) + gamma*velx - lambda_m);

            R_eigenvector[0][1] = 2*vely/(1 + v2 - 2*SQ(velx));
            R_eigenvector[1][1] = 2*velx*vely/(1 + v2 - 2*SQ(velx));
            R_eigenvector[2][1] = 1;
            R_eigenvector[3][1] = 0;

            R_eigenvector[0][2] = 2*velz/(1 + v2 - 2*SQ(velx));
            R_eigenvector[1][2] = 2*velx*velz/(1 + v2 - 2*SQ(velx));
            R_eigenvector[2][2] = 0;
            R_eigenvector[3][2] = 1;

            R_eigenvector[0][3] = (1 - SQ(velx))*(velx - lambda_p)*
                    (1 - (gamma - 1)*v2);
            R_eigenvector[1][3] = (1 - SQ(velx))*((lambda_p*velx - 1)*
                    (gamma - 1)*(1 - v2) + velx*(2 - gamma)*(velx - lambda_p));
            R_eigenvector[2][3] = vely*((gamma - 1)*(lambda_p*(1 + SQ(velx))*
                    v2 - lambda_p*SQ(velx) - 2*velx*v2) + (2 - gamma)*SQ(velx)*
                    (lambda_p - velx) + gamma*velx - lambda_p);
            R_eigenvector[3][3] = velz*((gamma - 1)*(lambda_p*(1 + SQ(velx))*
                    v2 - lambda_p*SQ(velx) - 2*velx*v2) + (2 - gamma)*SQ(velx)*
                    (lambda_p - velx) + gamma*velx - lambda_p);
        }

        static CCTK_REAL gamma;
};

} // namespace ultrarelativistic

} // namespace thc

#endif
