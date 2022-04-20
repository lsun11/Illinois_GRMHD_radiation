#include <weaktableindex.h>
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

INTEGER FUNCTION Emissions_cgs(rho, temp, ye,&
                             R_nue, R_nua, R_nux, Q_nue, Q_nua, Q_nux)
    use Inf_NaN_Detection
    use table3d_mod
    use lk_interpolations
    use weak_constants
    implicit none

    REAL*8, INTENT(IN)  :: rho, temp, ye
    REAL*8, INTENT(OUT) :: R_nue, R_nua, R_nux, Q_nue, Q_nua, Q_nux

    REAL*8 :: mass_fact_cgs

    REAL*8 :: beta, gamma, gamma_const, &
              block_factor_e, block_factor_a, block_factor_x, &
              enr_p, enr_m

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !All the emission rates are expressed per unit of volume
    !The number rates (R) are in  1.0 / (sec cm^3)
    !The energy rates (Q) are in  MeV/ (sec cm^3)

    REAL*8 :: Rbeta_nue, Rbeta_nua, &
              Qbeta_nue, Qbeta_nua, &
              Rpair, Qpair_Factor, pair_const,&
              Rpair_nue, Rpair_nua, Rpair_nux, &
              Qpair_nue, Qpair_nua, Qpair_nux, &
              Rgamma, Qgamma_Factor,&
              Rplasm_nue, Rplasm_nua, Rplasm_nux, &
              Qplasm_nue, Qplasm_nua, Qplasm_nux, &
              Rbrem, Qbrem

    REAL*8 :: fermi3, fermi4, fermi4o3, fermi5o4

    REAL*8 :: nb, lrho, ltemp, &
              mu_n, mu_p, mu_e, &
              abar, zbar, xp, xn, xh

    REAL*8 :: eta_nue, eta_nua, eta_nux, &
              eta_n, eta_p, eta_hat, eta_e, &
              eta_np, eta_pn, enr_tilde_p, &
              enr_tilde_m

    Emissions_cgs = 0

#define WEAK_RATES_ITS_ME
#include "inc/weak_rates_guts.inc"
#undef WEAK_RATES_ITS_ME

    !B5 B6 B7 Energy moments of electron and positrons
    enr_p = 8.0d0 * pi / hc_mevcm**3 * temp**4 * FERMI3(eta_e)
    enr_m = 8.0d0 * pi / hc_mevcm**3 * temp**4 * FERMI3(-eta_e)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !Emission of nue by the beta-process (electron and positron capture
    ! see eq.(27) (B1) (B2) see also Tubbs&Schramm, 1975; Bruenn, 1985

    !Coefficient
    beta = pi * clite * (1.0d0 + 3.0d0 * alpha**2) *&
                sigma_0 / (hc_mevcm**3 * me_mev**2)

    !Blocking factors (removed them to match Ott's leakage)
    block_factor_e = 1.0d0 + exp(eta_nue - FERMI5O4(eta_e))
    block_factor_a = 1.0d0 + exp(eta_nua - FERMI5O4(-eta_e))

    !neu electron capture rate
    Rbeta_nue = beta * eta_pn * temp**5 *  FERMI4(eta_e)/block_factor_e

    !neu electron capture energy rate
    Qbeta_nue = Rbeta_nue * temp * FERMI5O4(eta_e)

    !neu positron capture rate
    Rbeta_nua = beta * eta_np * temp**5 * FERMI4(-eta_e)/block_factor_a

    !neu positron capture energy rate
    Qbeta_nua = Rbeta_nua * temp * FERMI5O4(-eta_e)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !e-e+ pair processes from Ruffert et al. Ruffert 1998 formula (B9)
    !nu of the electron

    block_factor_e = 1.0d0 + exp(eta_nue - 0.5d0 *&
                ( FERMI4O3(eta_e) + FERMI4O3(-eta_e) ))

    !anti-nu of the electron
    block_factor_a = 1.0d0 + exp(eta_nua - 0.5d0 *&
                ( FERMI4O3(eta_e) + FERMI4O3(-eta_e) ))

    !nu of tau and mu
    block_factor_x = 1.0d0+exp(eta_nux-0.5d0 *&
                ( FERMI4O3(eta_e) + FERMI4O3(-eta_e) ))

    !B8 electron-positron pair annihilation
    pair_const = ((sigma_0 * clite) / me_mev**2) * enr_m * enr_p

    enr_tilde_m = 8.0d0 * pi / hc_mevcm**3 *temp**5 *&
              FERMI4(eta_e)
    enr_tilde_p = 8.0d0 * pi / hc_mevcm**3 *temp**5 *&
              FERMI4(-eta_e)
    !B8
    Rpair =  pair_const / (36.0d0 * block_factor_e * block_factor_a) *&
                ((Cv-Ca)**2 + (Cv+Ca)**2)

    Rpair_nue = Rpair
    Rpair_nua = Rpair

    Qpair_Factor = 0.5d0 * (temp * (FERMI4O3(-eta_e) + FERMI4O3(eta_e)))

    !Matching the factor in Ott's leakage
    !Qpair_Factor = 0.5d0*(enr_tilde_m*enr_p+enr_m*enr_tilde_p)/(enr_m*enr_p)

    Qpair_nue = Rpair * Qpair_Factor
    Qpair_nua = Rpair * Qpair_Factor

    !B10
    Rpair =  pair_const/(9.0d0 * block_factor_x**2) *&
             ((Cv-Ca)**2+(Cv + Ca - 2.0d0)**2)

    Rpair_nux = Rpair
    Qpair_nux = Rpair * Qpair_Factor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !plasmon decay from Ruffert et al.
    !Definition of gamma
    gamma = gamma_0 * sqrt((pi**2 + 3.0d0 * eta_e**2)/3.0d0)

    !Blocking factor for plasmon decay
    block_factor_e = 1.0d0 + &
                  exp(eta_nue -(1.0d0 + 0.5d0 * gamma**2/(1.0d0 + gamma)))
    block_factor_a = 1.0d0 + &
                  exp(eta_nua -(1.0d0 + 0.5d0 * gamma**2/(1.0d0 + gamma)))
    block_factor_x = 1.0d0 + &
                  exp(eta_nux -(1.0d0 + 0.5d0 * gamma**2/(1.0d0 + gamma)))

     gamma_const = pi**3 * sigma_0 * clite * temp**8 /&
                  (me_mev**2 * 3.0d0 * fsc * hc_mevcm**6)*&
                  gamma**6*exp(-gamma)*(1.0d0+gamma)

    !B11
    Rgamma = Cv**2 * gamma_const / (block_factor_e * block_factor_a)
    Qgamma_Factor = 0.5d0 * temp * (2.0d0 + gamma**2 / (1.0d0+gamma))

    Rplasm_nue = Rgamma
    Qplasm_nue = Rgamma * Qgamma_Factor

    Rplasm_nua = Rgamma
    Qplasm_nua = Rgamma * Qgamma_Factor

    !B12
    Rgamma = (Cv-1.0d0)**2 * 4.0d0 * gamma_const / block_factor_x**2
    Rplasm_nux = Rgamma
    Qplasm_nux = Rgamma * Qgamma_Factor

    ! Bremsstrahlung fitting formula described in
    ! A. Burrows et al. Nuclear Physics A 777 (2006) 356-394
    !TODO check the units MeV or erg
    Rbrem = 0.231d0 * (2.0778d2/mev_to_erg) * 0.5d0 * &
             (xn**2+ xp**2 + 28.0d0/3.0d0 * xn * xp) * &
                    rho**2 * temp**(4.5d0)
    Qbrem = Rbrem * temp / 0.231d0 * 0.504d0

    !Remeber the rates are in  MeV/ sec / cm^3
    R_nue = Rbeta_nue + Rpair_nue + Rplasm_nue !+ Rbrem
    R_nua = Rbeta_nua + Rpair_nua + Rplasm_nua !+ Rbrem
    R_nux = Rpair_nux + Rplasm_nux !+ 4.0d0 * Rbrem

    Q_nue = Qbeta_nue + Qpair_nue + Qplasm_nue !+ Qbrem
    Q_nua = Qbeta_nua + Qpair_nua + Qplasm_nua !+ Qbrem
    Q_nux = Qpair_nux + Qplasm_nux !+ 4.0d0 * Qbrem

    if (isnan(R_nue).or.isinf(R_nue)) then
       write(*,*) "Emissions_cgs: NaN/Inf in R_nue", rho, temp, ye
       Emissions_cgs = -1
    endif

    if (isnan(R_nua).or.isinf(R_nua)) then
       write(*,*) "Emissions_cgs: NaN/Inf in R_nua", rho, temp, ye
       Emissions_cgs = -1
    endif

    if (isnan(R_nux).or.isinf(R_nux)) then
       write(*,*) "Emissions_cgs: NaN/Inf in R_nux", rho, temp, ye
       Emissions_cgs = -1
    endif

    if (isnan(Q_nue).or.isinf(Q_nue)) then
       write(*,*) "Emissions_cgs: NaN/Inf in Q_nue", rho, temp, ye
       Emissions_cgs = -1
    endif

    if (isnan(Q_nua).or.isinf(Q_nua)) then
       write(*,*) "Emissions_cgs: NaN/Inf in Q_nua", rho, temp, ye
       Emissions_cgs = -1
    endif

    if (isnan(Q_nux).or.isinf(Q_nux)) then
       write(*,*) "Emissions_cgs: NaN/Inf in Q_nux", rho, temp, ye
       Emissions_cgs = -1
    endif

    return
END FUNCTION Emissions_cgs
