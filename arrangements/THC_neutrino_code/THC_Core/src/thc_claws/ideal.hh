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


#ifndef THC_CLAWS_IDEAL_HH
#define THC_CLAWS_IDEAL_HH

namespace thc {

namespace ideal {

class Euler: public EulerBase<Euler> {
    public:
        typedef hrscc::CLaw<Euler> claw;

        void prim_to_all(
                hrscc::Observer<claw> & observer
                ) const {
            IDL_MAKE_ALIAS(hrscc::policy::x, observer);
            DECLARE_CCTK_PARAMETERS

#ifdef THC_DEBUG
            THC_CLAWS_CHECK_ALL;
#endif

            // Compute volume form
            whizza::pz_metric3 g;
            whizza::collect(g.lo, gxx, gxy, gxz, gyy, gyz, gzz);
            volform = std::sqrt(whizza::determinant(g.lo));

            // If we are in the atmosphere we do not need to compute anything
            if(UTILS_BITMASK_CHECK_FLAG(bitmask, THC_FLAG_ATMOSPHERE) ||
               rho <  atmosphere::rho*(1.0 + atmosphere::tol)) {
                this->set_to_atmosphere(observer);
                return;
            }

            CCTK_REAL const fake_Y_e = 0.5;

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

            whizza::eos_thermal const eos =
                whizza::global_eos_thermal::get_eos();
            whizza::eos_thermal_impl::status status;

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

            // Enforce allowed range for eps
            whizza::eos_thermal::range rgeps = eos.range_eps(rho, fake_Y_e);
            if(eps < rgeps.min) {
                UTILS_BITMASK_SET_FLAG(bitmask, THC_FLAG_EPS_LT_EPSMIN);
                eps = rgeps.min;
            }
            else if(eps > rgeps.max) {
                UTILS_BITMASK_SET_FLAG(bitmask, THC_FLAG_EPS_GT_EPSMAX);
                eps = rgeps.max;
            }

            CCTK_REAL chi, kappa;
            eos.press_derivs_from_rho_eps_ye(press, chi, kappa, rho,
                    eps, fake_Y_e, status);
#ifdef THC_DEBUG
            if(status.failed) {
                std::ostringstream ss;
                this->print_location(observer, ss);
                Printer::print_err(ss.str());
#pragma omp critical
                CCTK_ERROR(status.err_msg.c_str());
            }
#endif

            CCTK_REAL const h = 1 + eps + press/rho;
            csound = std::sqrt((chi + press/SQ(rho)*kappa) / h);
            if(std::isnan(csound) || csound > 1) {
                UTILS_BITMASK_SET_FLAG(bitmask, THC_FLAG_CSOUND_IS_CRAP);
            }

            // Compute conserved variables
            dens  = volform*rho*w_lorentz;
            sconx = volform*rho*h*W2*v_x;
            scony = volform*rho*h*W2*v_y;
            sconz = volform*rho*h*W2*v_z;
            tau   = volform*rho*h*W2 - volform*press - dens;
            //CCTK_REAL const v2  = v_x*velx + v_y*vely + v_z*velz;
            //tau = volform*W2*(v2*(rho*w_lorentz/(1.0 + w_lorentz) + press) +
            //        rho*eps);

#ifdef THC_DEBUG
            THC_CLAWS_CHECK_ALL;
#endif
        }

        inline void set_to_atmosphere(
                hrscc::Observer<claw> & observer
                ) const {
            IDL_MAKE_ALIAS(hrscc::policy::x, observer);
            DECLARE_CCTK_PARAMETERS

            rho   = atmosphere::rho;
            zvecx = 0;
            zvecy = 0;
            zvecz = 0;
            eps   = atmosphere::eps;

            velx = 0;
            vely = 0;
            velz = 0;
            w_lorentz = 1;

            csound = atmosphere::csound;
            press  = atmosphere::press;

            densgain = volform*rho - dens;

            dens  = volform*rho;
            sconx = 0;
            scony = 0;
            sconz = 0;
            tau   = volform*rho*(1 + eps) - dens;

#ifdef THC_DEBUG
            THC_CLAWS_CHECK_ALL;
#endif
        }

        void cons_to_all(
                hrscc::Observer<claw> & observer
                ) const {
            IDL_MAKE_ALIAS(hrscc::policy::x, observer);
            DECLARE_CCTK_PARAMETERS
            std::ostringstream ss;

#ifdef THC_DEBUG
            THC_CLAWS_CHECK_ALL;
#endif

            densgain = 0;

            if(UTILS_BITMASK_CHECK_FLAG(bitmask, THC_FLAG_IS_INSIDE_BH)) {
                UTILS_BITMASK_SET_FLAG(bitmask, THC_FLAG_ATMOSPHERE);
                this->set_to_atmosphere(observer);
                return;
            }

            if(c2a_fix_nans) {
                if(!ISFINITE(dens) ||
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

            CCTK_REAL fake_Y_e = 0.5;
            CCTK_REAL fake_dens_Y_e = dens*fake_Y_e;

            try {
                if(dens < atmosphere::rho*(1.0 + atmosphere::tol)) {
                    UTILS_BITMASK_SET_FLAG(bitmask, THC_FLAG_ATMOSPHERE);
                    this->set_to_atmosphere(observer);
                    return;
                }

                whizza::atmosphere atmo(
                        atmosphere::rho,
                        atmosphere::eps,
                        fake_Y_e,
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
                cv.collect(dens, tau, fake_dens_Y_e, sconx, scony, sconz);

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
                        cv.scatter(dens, tau, fake_dens_Y_e,
                                sconx, scony, sconz);
                        pv.scatter(rho, eps, fake_Y_e, press, velx, vely, velz,
                                   w_lorentz);
                        zvecx = w_lorentz*velx;
                        zvecy = w_lorentz*vely;
                        zvecz = w_lorentz*velz;
                        densgain = dens - dens_old;
                        Printer::print_warn(ss.str());
                        return;
                    }
                    else {
                        ss << "Adjusting conservative variables\n";
                        cv.collect(dens, tau, fake_dens_Y_e, sconx, scony,
                                sconz);
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

                pv.scatter(rho, eps, fake_Y_e, press, velx, vely, velz,
                       w_lorentz);
                zvecx = w_lorentz*velx;
                zvecy = w_lorentz*vely;
                zvecz = w_lorentz*velz;

                CCTK_REAL chi, kappa;
                eos.press_derivs_from_rho_eps_ye(press, chi, kappa, rho,
                        eps, fake_Y_e, status);
#ifdef THC_DEBUG
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

                CCTK_REAL const h = 1 + eps + press/rho;
                csound = std::sqrt((chi + press/SQ(rho)*kappa) / h);
                if(std::isnan(csound) || csound > 1) {
                    UTILS_BITMASK_SET_FLAG(bitmask,THC_FLAG_CSOUND_IS_CRAP);
                }

                if(c2a_fix_conservatives && report.adjust_cons) {
                    CCTK_REAL const dens_old = dens;
                    cv.scatter(dens, tau, fake_dens_Y_e, sconx, scony, sconz);
                    densgain = dens - dens_old;
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
            IDL_MAKE_ALIAS(dir, observer);

#ifdef THC_DEBUG
            THC_CLAWS_CHECK_ALL;
#endif

            flux_dens  = alp*(dens* (velx - betax/alp));
            flux_sconx = alp*(sconx*(velx - betax/alp) + volform*press);
            flux_scony = alp*(scony*(velx - betax/alp));
            flux_sconz = alp*(sconz*(velx - betax/alp));
            flux_tau   = alp*(tau * (velx - betax/alp) + volform*press*velx);

#ifdef THC_DEBUG
            THC_CLAWS_CHECK_ALL;
#endif
        }

        template<hrscc::policy::direction_t dir>
        void eigenvalues(
                hrscc::Observer<claw> & observer
                ) const {
            IDL_MAKE_ALIAS(dir, observer);

#ifdef THC_DEBUG
            THC_CLAWS_CHECK_ALL;
#endif

            if(UTILS_BITMASK_CHECK_FLAG(bitmask, THC_FLAG_CSOUND_IS_CRAP)) {
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

        // From L.Rezzolla and O.Zanotti 2013
        template<hrscc::policy::direction_t dir>
        void eig(
                hrscc::Observer<claw> & observer
                ) const {
            IDL_MAKE_ALIAS(dir, observer);
            DECLARE_CCTK_PARAMETERS

#ifdef THC_DEBUG
            THC_CLAWS_CHECK_ALL;
#endif

            CCTK_REAL fake_Y_e = 0.5;

            if(UTILS_BITMASK_CHECK_FLAG(bitmask, THC_FLAG_CSOUND_IS_CRAP) ||
                    rho < eig_min_rho || csound < eig_min_csound ||
                    alp < eig_min_alp) {
                observer.reset_eigenvalues();
                observer.reset_eigenvectors();
                return;
            }

            ///////////////////////////////////////////////////////////////////
            // Call EOS
            ///////////////////////////////////////////////////////////////////
            whizza::eos_thermal const eos =
                whizza::global_eos_thermal::get_eos();
            whizza::eos_thermal_impl::status status;

            CCTK_REAL chi, kappa;
            eos.press_derivs_from_rho_eps_ye(press, chi, kappa, rho,
                    eps, fake_Y_e, status);
#ifdef THC_DEBUG
            if(status.failed) {
                std::ostringstream ss;
                this->print_location(observer, ss);
                Printer::print_err(ss.str());
#pragma omp critical
                CCTK_ERROR(status.err_msg.c_str());
            }
#endif

            CCTK_REAL const h = 1 + eps + press/rho;

            ///////////////////////////////////////////////////////////////////
            // Compute eigenvalues
            ///////////////////////////////////////////////////////////////////
            THC_GENERAL_EULER_COMPUTE_EIGENVALUES;

#undef THC_GENERAL_EULER_COMPUTE_EIGENVALUES

            CCTK_REAL const det = g.det;
            ///////////////////////////////////////////////////////////////////
            // Some usefull constants
            ///////////////////////////////////////////////////////////////////
            CCTK_REAL const W2 = SQ(w_lorentz);

            CCTK_REAL const tilde_kappa = kappa/rho;
            CCTK_REAL const Kappa = tilde_kappa/(tilde_kappa - SQ(csound));

            CCTK_REAL const lambda_m = p_eigenvalue[0];
            CCTK_REAL const Lambda_m = (lambda_m + betax)/alp;
            CCTK_REAL const A_m = (g.up(0,0) - velx*velx)/
                    (g.up(0,0) - velx*Lambda_m);
            CCTK_REAL const V_m = (velx - Lambda_m)/
                    (g.up(0,0) - velx*Lambda_m);
            CCTK_REAL const C_m = v_x - V_m;

            CCTK_REAL const lambda_p = p_eigenvalue[4];
            CCTK_REAL const Lambda_p = (lambda_p + betax)/alp;
            CCTK_REAL const A_p = (g.up(0,0) - velx*velx)/
                    (g.up(0,0) - velx*Lambda_p);
            CCTK_REAL const V_p = (velx - Lambda_p)/
                    (g.up(0,0) - velx*Lambda_p);
            CCTK_REAL const C_p = v_x - V_p;

            CCTK_REAL const Gamma_xx = gyy*gzz - SQ(gyz);
            CCTK_REAL const Gamma_xy = -(gxy*gzz - gxz*gyz);
            CCTK_REAL const Gamma_xz = gxy*gyz - gyy*gxz;
            CCTK_REAL const xi = Gamma_xx - det*SQ(velx);
            CCTK_REAL const Delta = CB(h)*w_lorentz*(Kappa - 1)*(C_p - C_m)*xi;

            CCTK_REAL const N_m = (1 - Kappa)*(-det*velx + V_p*(W2*xi -
                        Gamma_xx)) - Kappa*W2*V_p*xi;
            CCTK_REAL const N_p = (1 - Kappa)*(-det*velx + V_m*(W2*xi -
                        Gamma_xx)) - Kappa*W2*V_m*xi;

            CCTK_REAL const c_m = h*h/Delta;
            CCTK_REAL const c_1 = w_lorentz/(Kappa - 1);
            CCTK_REAL const c_2 = 1/(h*xi);
            CCTK_REAL const c_3 = 1/(h*xi);
            CCTK_REAL const c_p = - h*h/Delta;

#if THC_DEBUG
            THC_CLAWS_CHECK(c_m);
            THC_CLAWS_CHECK(c_1);
            THC_CLAWS_CHECK(c_2);
            THC_CLAWS_CHECK(c_3);
            THC_CLAWS_CHECK(c_p);
#endif

            ///////////////////////////////////////////////////////////////////
            // Left eigenvectors
            ///////////////////////////////////////////////////////////////////
            L_eigenvector[0][0] = c_m*(h*w_lorentz*V_p*xi + N_m);
            L_eigenvector[0][1] = c_m*(Gamma_xx*(1 - Kappa*A_p) + (2*Kappa - 1)
                    *V_p*(W2*velx*xi - Gamma_xx*velx));
            L_eigenvector[0][2] = c_m*(Gamma_xy*(1 - Kappa*A_p) + (2*Kappa - 1)
                    *V_p*(W2*vely*xi - Gamma_xy*velx));
            L_eigenvector[0][3] = c_m*(Gamma_xz*(1 - Kappa*A_p) + (2*Kappa - 1)
                    *V_p*(W2*velz*xi - Gamma_xz*velx));
            L_eigenvector[0][4] = c_m*N_m;

            L_eigenvector[1][0] = c_1*(h - w_lorentz);
            L_eigenvector[1][1] = c_1*w_lorentz*velx;
            L_eigenvector[1][2] = c_1*w_lorentz*vely;
            L_eigenvector[1][3] = c_1*w_lorentz*velz;
            L_eigenvector[1][4] = - c_1*w_lorentz;

            L_eigenvector[2][0] = c_2*(-gzz*v_y + gyz*v_z);
            L_eigenvector[2][1] = c_2*(velx*(gzz*v_y - gyz*v_z));
            L_eigenvector[2][2] = c_2*( gzz*(1 - v_x*velx) + gxz*v_z*velx);
            L_eigenvector[2][3] = c_2*(-gyz*(1 - v_x*velx) - gxz*v_y*velx);
            L_eigenvector[2][4] = c_2*(-gzz*v_y + gyz*v_z);

            L_eigenvector[3][0] = c_3*(-gyy*v_z + gyz*v_y);
            L_eigenvector[3][1] = c_3*(velx*(gyy*v_z - gyz*v_y));
            L_eigenvector[3][2] = c_3*(-gyz*(1 - v_x*velx) - gxy*v_z*velx);
            L_eigenvector[3][3] = c_3*( gyy*(1 - v_x*velx) + gxy*v_y*velx);
            L_eigenvector[3][4] = c_3*(-gyy*v_z + gyz*v_y);

            L_eigenvector[4][0] = c_p*(h*w_lorentz*V_m*xi + N_p);
            L_eigenvector[4][1] = c_p*(Gamma_xx*(1 - Kappa*A_m) + (2*Kappa - 1)
                    *V_m*(W2*velx*xi - Gamma_xx*velx));
            L_eigenvector[4][2] = c_p*(Gamma_xy*(1 - Kappa*A_m) + (2*Kappa - 1)
                    *V_m*(W2*vely*xi - Gamma_xy*velx));
            L_eigenvector[4][3] = c_p*(Gamma_xz*(1 - Kappa*A_m) + (2*Kappa - 1)
                    *V_m*(W2*velz*xi - Gamma_xz*velx));
            L_eigenvector[4][4] = c_p*N_p;

            ///////////////////////////////////////////////////////////////////
            // Right eigenvectors
            ///////////////////////////////////////////////////////////////////
            R_eigenvector[0][0] = 1;
            R_eigenvector[1][0] = h*w_lorentz*(v_x - V_m);
            R_eigenvector[2][0] = h*w_lorentz*v_y;
            R_eigenvector[3][0] = h*w_lorentz*v_z;
            R_eigenvector[4][0] = h*w_lorentz*A_m - 1;

            R_eigenvector[0][1] = Kappa/(h*w_lorentz);
            R_eigenvector[1][1] = v_x;
            R_eigenvector[2][1] = v_y;
            R_eigenvector[3][1] = v_z;
            R_eigenvector[4][1] = 1 - Kappa/(h*w_lorentz);

            R_eigenvector[0][2] = w_lorentz*v_y;
            R_eigenvector[1][2] = h*(gxy + 2*W2*v_x*v_y);
            R_eigenvector[2][2] = h*(gyy + 2*W2*v_y*v_y);
            R_eigenvector[3][2] = h*(gyz + 2*W2*v_z*v_y);
            R_eigenvector[4][2] = w_lorentz*v_y*(2*h*w_lorentz - 1);

            R_eigenvector[0][3] = w_lorentz*v_z;
            R_eigenvector[1][3] = h*(gxz + 2*W2*v_x*v_z);
            R_eigenvector[2][3] = h*(gyz + 2*W2*v_y*v_z);
            R_eigenvector[3][3] = h*(gzz + 2*W2*v_z*v_z);
            R_eigenvector[4][3] = w_lorentz*v_z*(2*h*w_lorentz - 1);

            R_eigenvector[0][4] = 1;
            R_eigenvector[1][4] = h*w_lorentz*(v_x - V_p);
            R_eigenvector[2][4] = h*w_lorentz*v_y;
            R_eigenvector[3][4] = h*w_lorentz*v_z;
            R_eigenvector[4][4] = h*w_lorentz*A_p - 1;

#ifdef THC_DEBUG
            THC_CLAWS_CHECK_ALL;
#endif
        }
};

} // namespace ideal

} // namespace thc

#undef THC_CLAWS_GENERAL_ASSERT

#endif
