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


#ifndef THC_CLAWS_NUCLEAR_HH
#define THC_CLAWS_NUCLEAR_HH

namespace thc {
namespace nuclear {

class Euler: public EulerBase<Euler> {
    public:
        typedef hrscc::CLaw<Euler> claw;

        inline void prim_to_all(
                hrscc::Observer<claw> & observer
                ) const {
            NUC_MAKE_ALIAS(hrscc::policy::x, observer);
            DECLARE_CCTK_PARAMETERS

#ifdef THC_DEBUG
            THC_CLAWS_CHECK_ALL;
#endif

            // Compute volume form
            whizza::pz_metric3 g;
            whizza::collect(g.lo, gxx, gxy, gxz, gyy, gyz, gzz);
            volform = std::sqrt(whizza::determinant(g.lo));

            // If we are in the atmosphere we do not need to compute anything
            if(UTILS_BITMASK_CHECK_FLAG(bitmask, THC_FLAG_ATMOSPHERE)) {
                this->set_to_atmosphere(observer);
                return;
            }

            // Compute velocity and Lorentz factor
            CCTK_REAL const z_x = gxx*zvecx + gxy*zvecy + gxz*zvecz;
            CCTK_REAL const z_y = gxy*zvecx + gyy*zvecy + gyz*zvecz;
            CCTK_REAL const z_z = gxz*zvecx + gyz*zvecy + gzz*zvecz;

            CCTK_REAL const W2  = z_x*zvecx + z_y*zvecy + z_z*zvecz + 1;
            w_lorentz = std::sqrt(W2);

            velx = zvecx/w_lorentz;
            vely = zvecy/w_lorentz;
            velz = zvecz/w_lorentz;

            CCTK_REAL const v_x = z_x/w_lorentz;
            CCTK_REAL const v_y = z_y/w_lorentz;
            CCTK_REAL const v_z = z_z/w_lorentz;

            // Get EoS object
            whizza::eos_thermal_impl const & eos =
                whizza::global_eos_thermal::get_eos().implementation();

            // Enforce allowed range for rho
            whizza::eos_thermal::range rgrho = eos.range_rho();
            if(rho < rgrho.min) {
                UTILS_BITMASK_SET_FLAG(bitmask, THC_FLAG_RHO_LT_RHOMIN);
                rho = rgrho.min;
            }
            else if(rho > rgrho.max) {
                UTILS_BITMASK_SET_FLAG(bitmask, THC_FLAG_RHO_GT_RHOMAX);
                rho = rgrho.max;
            }

            // Enforce allowed range for Y_e
            whizza::eos_thermal::range rgY_e = eos.range_ye();
            if(Y_e < rgY_e.min) {
                UTILS_BITMASK_SET_FLAG(bitmask, THC_FLAG_YE_LT_YEMIN);
                Y_e = rgY_e.min;
            }
            else if(Y_e > rgY_e.max) {
                UTILS_BITMASK_SET_FLAG(bitmask, THC_FLAG_YE_GT_YEMAX);
                Y_e = rgY_e.max;
            }

            // Enforce allowed range for eps
            whizza::eos_thermal::range rgeps = eos.range_eps(rho, Y_e);
            if(eps < rgeps.min) {
                UTILS_BITMASK_SET_FLAG(bitmask, THC_FLAG_EPS_LT_EPSMIN);
                eps = rgeps.min;
            }
            else if(eps > rgeps.max) {
                UTILS_BITMASK_SET_FLAG(bitmask, THC_FLAG_EPS_GT_EPSMAX);
                eps = rgeps.max;
            }

            // Compute temperature, entropy, press and csound
            temperature = eos.temp_from_valid_rho_eps_ye(rho, eps, Y_e);
            entropy = eos.entropy_from_valid_rho_temp_ye(rho, temperature, Y_e);
            press = eos.press_from_valid_rho_temp_ye(rho, temperature, Y_e);
            csound = eos.csnd_from_valid_rho_temp_ye(rho, temperature, Y_e);
            if(std::isnan(csound) || csound > 1) {
                UTILS_BITMASK_SET_FLAG(bitmask,THC_FLAG_CSOUND_IS_CRAP);
                csound = 0;
            }

            // Compute conserved variables
            CCTK_REAL const h = 1 + eps + press/rho;
            densxn = volform*(1 - Y_e)*rho*w_lorentz;
            densxp = volform*Y_e*rho*w_lorentz;
            sconx  = volform*rho*h*W2*v_x;
            scony  = volform*rho*h*W2*v_y;
            sconz  = volform*rho*h*W2*v_z;
            tau    = volform*rho*h*W2 - volform*press - (densxn + densxp);

            dens   = densxn + densxp;

#ifdef THC_DEBUG
            THC_CLAWS_CHECK_ALL;
#endif
        }

        inline void set_to_atmosphere(
                hrscc::Observer<claw> & observer
                ) const {
            NUC_MAKE_ALIAS(hrscc::policy::x, observer);
            DECLARE_CCTK_PARAMETERS

            rho   = atmosphere::rho;
            zvecx = 0;
            zvecy = 0;
            zvecz = 0;
            eps   = atmosphere::eps;
            Y_e   = atmosphere::Ye;

            velx = 0;
            vely = 0;
            velz = 0;
            w_lorentz = 1;

            csound      = atmosphere::csound;
            entropy     = atmosphere::entropy;
            press       = atmosphere::press;
            temperature = atmosphere::temp;

            densgain = volform*rho - densxp - densxn;

            densxn = volform*(1-Y_e)*rho;
            densxp = volform*Y_e*rho;
            sconx  = 0;
            scony  = 0;
            sconz  = 0;
            tau    = volform*rho*(1 + eps) - (densxn + densxp);

            dens   = densxn + densxp;

#ifdef THC_DEBUG
            THC_CLAWS_CHECK_ALL;
#endif
        }

        inline void cons_to_all(
                hrscc::Observer<claw> & observer
                ) const {
            NUC_MAKE_ALIAS(hrscc::policy::x, observer);
            DECLARE_CCTK_PARAMETERS
            std::ostringstream ss;

            if(UTILS_BITMASK_CHECK_FLAG(bitmask, THC_FLAG_IS_INSIDE_BH)) {
                UTILS_BITMASK_SET_FLAG(bitmask, THC_FLAG_ATMOSPHERE);
                this->set_to_atmosphere(observer);
                return;
            }

            if(c2a_fix_nans) {
                if(!ISFINITE(densxn) || !ISFINITE(densxp) ||
                   !ISFINITE(sconx) || !ISFINITE(scony) || !ISFINITE(sconz) ||
                   !ISFINITE(tau)) {
                    UTILS_BITMASK_SET_FLAG(bitmask, THC_FLAG_NOT_FINITE);
                    this->set_to_atmosphere(observer);
                    return;
                }
                else {
                    UTILS_BITMASK_UNSET_FLAG(bitmask, THC_FLAG_NOT_FINITE);
                }
            }

#ifdef THC_DEBUG
            THC_CLAWS_CHECK_ALL;
#endif

            densgain = 0;

            try {
                // Basic atmosphere checks on the total density
                dens = densxn + densxp;

                if(dens < atmosphere::rho*(1.0 + atmosphere::tol)) {
                    UTILS_BITMASK_SET_FLAG(bitmask, THC_FLAG_ATMOSPHERE);
                    this->set_to_atmosphere(observer);
                    return;
                }

                // Setup the c2p object
                whizza::atmosphere atmo(
                        atmosphere::rho,
                        atmosphere::eps,
                        atmosphere::Ye,
                        atmosphere::press,
                        atmosphere::rho*(1.0 + atmosphere::tol));

                whizza::eos_thermal const eos =
                    whizza::global_eos_thermal::get_eos();
                whizza::eos_thermal_impl::status status;

                whizza::con2prim const c2p_strict(eos, c2a_rho_strict,
                        c2a_zmax, atmo, c2a_accuracy, c2a_max_iter,
                        c2a_fail_if_eps_lt_min, c2a_fail_if_eps_gt_max,
                        c2a_fail_if_ye_lt_min, c2a_fail_if_ye_gt_max);
                whizza::con2prim const c2p_safe(eos, eos.range_rho().max,
                        c2a_zmax, atmo, c2a_accuracy, c2a_max_iter,
                        c2a_fail_if_eps_lt_min, c2a_fail_if_eps_gt_max,
                        c2a_fail_if_ye_lt_min, c2a_fail_if_ye_gt_max);

                whizza::pz_metric3 g;
                whizza::collect(g.lo, gxx, gxy, gxz, gyy, gyz, gzz);
                g.complete_from_lower();
                volform = std::sqrt(g.det);

                whizza::cons_vars cv;
                cv.collect(dens, tau, densxp, sconx, scony, sconz);

                whizza::prim_vars pv;

                whizza::con2prim::report report;
                if(alp < c2a_BH_alp) {
                    c2p_safe(pv, cv, g, report);
                }
                else {
                    c2p_strict(pv, cv, g, report);
                }

                if(report.failed) {
                    UTILS_BITMASK_SET_FLAG(bitmask, THC_FLAG_C2A_FAILURE);
                    ss << "ConsToAll failure!\n" << report.err_msg << std::endl;
                    this->print_location(observer, ss);

                    if(c2a_kill_on_failure) {
                        Printer::print_err(ss.str());
#pragma omp critical
                        CCTK_ERROR("Aborting the run");
                    }

                    if(c2a_set_to_nan_on_failure) {
                        CCTK_REAL const dens_old = dens;
                        cv.scatter(dens, tau, densxp, sconx, scony, sconz);
                        densgain = dens - dens_old;
                        densxn = dens - densxp;
                        pv.scatter(rho, eps, Y_e, press, velx, vely, velz,
                                   w_lorentz);
                        zvecx = w_lorentz*velx;
                        zvecy = w_lorentz*vely;
                        zvecz = w_lorentz*velz;
                        Printer::print_warn(ss.str());
                        return;
                    }
                    else {
                        ss << "Adjusting conservative variables\n";
                        cv.collect(dens, tau, densxp, sconx, scony, sconz);
                        c2p_safe(pv, cv, g, report);
                        if(report.failed) {
                            ss << "Adjusting failed\n";
                            ss << "Setting point to atmosphere\n";
                            UTILS_BITMASK_SET_FLAG(bitmask,THC_FLAG_ATMOSPHERE);
                            this->set_to_atmosphere(observer);
                            Printer::print_warn(ss.str());
                            return;
                        }
                    }
                }

                c2a_nbiter = static_cast<CCTK_REAL>(report.iters);

                if(report.set_atmo) {
                    UTILS_BITMASK_SET_FLAG(bitmask, THC_FLAG_ATMOSPHERE);
                    this->set_to_atmosphere(observer);
                    if(ss.tellp() > 0) {
                        Printer::print_warn(ss.str());
                    }
                    return;
                }

                // Compute the primitives
                pv.scatter(rho, eps, Y_e, press, velx, vely, velz,
                       w_lorentz);

                zvecx = w_lorentz*velx;
                zvecy = w_lorentz*vely;
                zvecz = w_lorentz*velz;

                // Check that the variables are in their allowed ranges
                whizza::eos_thermal_impl const & ieos = eos.implementation();
                ieos.check_rho_eps_ye(rho, eps, Y_e, status);
#if THC_DEBUG
                if(status.failed) {
                    if(c2a_fix_nans) {
                        this->set_to_atmosphere(observer);
                        return;
                    }
                    this->print_location(observer, ss);
                    Printer::print_err(ss.str());
#pragma omp critical
                    CCTK_ERROR(status.err_msg.c_str());
                }
#endif

                // Compute other thermodynamical quantities
                temperature = ieos.temp_from_valid_rho_eps_ye(rho, eps, Y_e);
                entropy = ieos.entropy_from_valid_rho_temp_ye(rho,
                        temperature, Y_e);
                csound = ieos.csnd_from_valid_rho_temp_ye(rho,
                        temperature, Y_e);
                if(std::isnan(csound) || csound > 1) {
                    UTILS_BITMASK_SET_FLAG(bitmask,THC_FLAG_CSOUND_IS_CRAP);
                    csound = 0;
                }

                if(c2a_fix_conservatives && report.adjust_cons) {
                    CCTK_REAL const dens_old = dens;
                    cv.scatter(dens, tau, densxp, sconx, scony, sconz);
                    densgain = dens - dens_old;
                    densxn = dens - densxp;
                }
            }
            catch(std::exception & e) {
                if(c2a_fix_nans) {
                    this->set_to_atmosphere(observer);
                    return;
                }
                this->print_location(observer, ss);
                Printer::print_err(ss.str());
#pragma omp critical
                CCTK_ERROR(e.what());
            }

#ifdef THC_DEBUG
            THC_CLAWS_CHECK_ALL;
#endif
            if(ss.tellp() > 0) {
                Printer::print_warn(ss.str());
            }
        }

        template<hrscc::policy::direction_t dir>
        inline void fluxes(
                hrscc::Observer<claw> & observer
                ) const {
            NUC_MAKE_ALIAS(dir, observer);

#ifdef THC_DEBUG
            THC_CLAWS_CHECK_ALL;
#endif

            flux_densxn = alp*(densxn*(velx - betax/alp));
            flux_densxp = alp*(densxp*(velx - betax/alp));
            flux_sconx  = alp*(sconx* (velx - betax/alp) + volform*press);
            flux_scony  = alp*(scony* (velx - betax/alp));
            flux_sconz  = alp*(sconz* (velx - betax/alp));
            flux_tau    = alp*(tau*   (velx - betax/alp) + volform*press*velx);

#ifdef THC_DEBUG
            THC_CLAWS_CHECK_ALL;
#endif
        }

        template<hrscc::policy::direction_t dir>
        inline void eigenvalues(
                hrscc::Observer<claw> & observer
                ) const {
            NUC_MAKE_ALIAS(dir, observer);

#ifdef THC_DEBUG
            THC_CLAWS_CHECK_ALL;
#endif

            if(UTILS_BITMASK_CHECK_FLAG(bitmask, THC_FLAG_CSOUND_IS_CRAP)) {
                csound = 0;
                observer.reset_eigenvalues();
                return;
            }

#define THC_GENERAL_EULER_COMPUTE_EIGENVALUES                                  \
            whizza::pz_metric3 g;                                              \
            whizza::collect(g.lo, gxx, gxy, gxz, gyy, gyz, gzz);               \
            g.complete_from_lower();                                           \
                                                                               \
            CCTK_REAL const v_x = gxx*velx + gxy*vely + gxz*velz;              \
            CCTK_REAL const v_y = gxy*velx + gyy*vely + gyz*velz;              \
            CCTK_REAL const v_z = gxz*velx + gyz*vely + gzz*velz;              \
            CCTK_REAL const v2  = v_x*velx + v_y*vely + v_z*velz;              \
                                                                               \
            p_eigenvalue[0] = alp/(1-v2*SQ(csound))*(velx*(1-SQ(csound)) -     \
                    csound*std::sqrt((1-v2)*(g.up(0,0)*(1-v2*SQ(csound)) -     \
                    velx*velx*(1-SQ(csound))))) - betax;                       \
            p_eigenvalue[1] = alp*velx - betax;                                \
            p_eigenvalue[2] = alp*velx - betax;                                \
            p_eigenvalue[3] = alp*velx - betax;                                \
            p_eigenvalue[4] = alp/(1-v2*SQ(csound))*(velx*(1-SQ(csound)) +     \
                    csound*std::sqrt((1-v2)*(g.up(0,0)*(1-v2*SQ(csound)) -     \
                    velx*velx*(1-SQ(csound))))) - betax

            THC_GENERAL_EULER_COMPUTE_EIGENVALUES;

#ifdef THC_DEBUG
            THC_CLAWS_CHECK_ALL;
#endif
        }

        template<hrscc::policy::direction_t dir>
        inline void eig(
                hrscc::Observer<claw> & observer
                ) const {
#pragma omp critical
            CCTK_ERROR("nuclear::Euler::eig not implemented!");
        }
};

} // namespace nuclear
} // namespace thc


#endif
