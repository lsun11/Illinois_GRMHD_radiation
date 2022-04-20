#include <weaktableindex.h>
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

INTEGER FUNCTION Opacities_cgs(rho, temp, ye, &
                                kappa_0_nue, kappa_0_nua, kappa_0_nux, &
                                kappa_1_nue, kappa_1_nua, kappa_1_nux)

    !----------------------------------------------------------------------
    !OPACITIES_CGS(RHO, TEMP, YE)
    ! Function to compute the total opacity coefficients in cgs (cm^-1).
    ! Input::     rho (g cm^-3), temp (MeV), ye
    ! Output::    kappa_0_nue electron nu opct. for number transport.
    !             kappa_0_nua electron anti-nu opct. for number transport.
    !             kappa_0_nux heavy nu opct. for number transport.
    !             kappa_1_nue electron nu opct. for energy transport.
    !             kappa_1_nua electron anti-nu opct. for energy transport.
    !             kappa_1_nux heavy nu opct. for energy tranport.
    !----------------------------------------------------------------------

    use Inf_NaN_Detection
    use table3d_mod
    use weak_constants

    IMPLICIT NONE

    CCTK_REAL, INTENT(IN)  :: rho, temp, ye
    CCTK_REAL, INTENT(OUT) :: kappa_0_nue, kappa_0_nua, kappa_0_nux,&
                              kappa_1_nue, kappa_1_nua, kappa_1_nux

    CCTK_REAL :: kappa_0_nue_abs, kappa_0_nua_abs, kappa_0_nux_abs,&
                 kappa_1_nue_abs, kappa_1_nua_abs, kappa_1_nux_abs,&
                 kappa_0_nue_sct, kappa_0_nua_sct, kappa_0_nux_sct,&
                 kappa_1_nue_sct, kappa_1_nua_sct, kappa_1_nux_sct

    INTEGER :: Absorption_cgs, Scattering_cgs
    INTEGER :: err = 0

    character(len=100) warnline

    Opacities_cgs = 0

    kappa_0_nue_abs = 0.0d0
    kappa_0_nua_abs = 0.0d0
    kappa_0_nux_abs = 0.0d0
    kappa_1_nue_abs = 0.0d0
    kappa_1_nua_abs = 0.0d0
    kappa_1_nux_abs = 0.0d0

    kappa_0_nue_sct = 0.0d0
    kappa_0_nua_sct = 0.0d0
    kappa_0_nux_sct = 0.0d0
    kappa_1_nue_sct = 0.0d0
    kappa_1_nua_sct = 0.0d0
    kappa_1_nux_sct = 0.0d0


    err = Absorption_cgs(rho, temp, ye, &
                         kappa_0_nue_abs, kappa_0_nua_abs, kappa_0_nux_abs,&
                         kappa_1_nue_abs, kappa_1_nua_abs, kappa_1_nux_abs)

    if (err.ne.0) then
       write(*,*) "Opacities_cgs: Problem in Absorption_cgs"
       Opacities_cgs = -1
    endif

     err = Scattering_cgs(rho, temp, ye, &
                          kappa_0_nue_sct, kappa_0_nua_sct, kappa_0_nux_sct,&
                          kappa_1_nue_sct, kappa_1_nua_sct, kappa_1_nux_sct)

    if (err.ne.0) then
       write(*,*) "Opacities_cgs: Problem in Scattering_cgs"
       Opacities_cgs = -1
    endif

    if (isnan(kappa_0_nue_abs).or.isinf(kappa_0_nue_abs)) then
       write(*,*) "Opacities_cgs: NaN/Inf in kappa_0_nue_abs"
       Opacities_cgs = -1
    endif

    if (isnan(kappa_0_nua_abs).or.isinf(kappa_0_nua_abs)) then
       write(*,*) "Opacities_cgs: NaN/Inf in kappa_0_nua_abs"
       Opacities_cgs = -1
    endif

    if (isnan(kappa_0_nux_abs).or.isinf(kappa_0_nux_abs)) then
       write(*,*) "Opacities_cgs: NaN/Inf in kappa_0_nux_abs"
       Opacities_cgs = -1
    endif

    if (isnan(kappa_0_nue_sct).or.isinf(kappa_0_nue_sct)) then
       write(*,*) "Opacities_cgs: NaN/Inf in kappa_0_nue_sct"
       Opacities_cgs = -1
    endif

    if (isnan(kappa_0_nua_sct).or.isinf(kappa_0_nua_sct)) then
       write(*,*) "Opacities_cgs: NaN/Inf in kappa_0_nua_sct"
       Opacities_cgs = -1
    endif

    if (isnan(kappa_0_nux_sct).or.isinf(kappa_0_nux_sct)) then
       write(*,*) "Opacities_cgs: NaN/Inf in kappa_0_nux_sct"
       Opacities_cgs = -1
    endif

    kappa_0_nue = kappa_0_nue_abs + kappa_0_nue_sct
    kappa_0_nua = kappa_0_nua_abs + kappa_0_nua_sct
    kappa_0_nux = kappa_0_nux_abs + kappa_0_nux_sct

    kappa_1_nue = kappa_1_nue_abs + kappa_1_nue_sct
    kappa_1_nua = kappa_1_nua_abs + kappa_1_nua_sct
    kappa_1_nux = kappa_1_nux_abs + kappa_1_nux_sct

    return
END FUNCTION Opacities_cgs

INTEGER FUNCTION Absorption_cgs(rho, temp, ye,&
                                 kappa_0_nue, kappa_0_nua, kappa_0_nux, &
                                 kappa_1_nue, kappa_1_nua, kappa_1_nux)

    use Inf_NaN_Detection
    use table3d_mod
    use weak_constants
    use lk_interpolations

    IMPLICIT NONE

    CCTK_REAL, INTENT(IN) :: rho, temp, ye
    CCTK_REAL, INTENT(OUT) :: kappa_0_nue, kappa_0_nua, kappa_0_nux, &
                              kappa_1_nue, kappa_1_nua, kappa_1_nux

    CCTK_REAL :: zeta_nue_abs_n, zeta_nue_abs_p, zeta_nue_abs_h, &
                 zeta_nua_abs_n, zeta_nua_abs_p, zeta_nua_abs_h, &
                 zeta_nux_abs_n, zeta_nux_abs_p, zeta_nux_abs_h, &
                 zeta_nue, zeta_nua, zeta_nux

    CCTK_REAL :: mass_fact_cgs
    CCTK_REAL :: block_factor, abs_zeta

    CCTK_REAL :: FERMI5O4, FERMI5O3, FERMI4O2

    CCTK_REAL :: nb, lrho, ltemp, &
                 mu_n, mu_p, mu_e, &
                 abar, zbar, xp, xn, xh

    CCTK_REAL :: eta_nue, eta_nua, eta_nux, &
                 eta_n, eta_p, eta_hat, eta_e, &
                 eta_np, eta_pn

    Absorption_cgs = 0

#define WEAK_RATES_ITS_ME
#include "inc/weak_rates_guts.inc"
#undef WEAK_RATES_ITS_ME


    !---------------------------------------------------------------------
    ! ABSOPTION
    !---------------------------------------------------------------------

    abs_zeta = (1.0d0 + 3.0d0 * alpha**2) * 0.25d0 * sigma_0/me_erg**2

    block_factor = 1.0d0 + exp(eta_e - FERMI5O4(eta_nue))
    zeta_nue_abs_n = eta_np * abs_zeta / block_factor
    zeta_nua_abs_n = 0.0d0 !no absorption of e-anti-nu on neutrons
    zeta_nux_abs_n = 0.0d0 !no absorption of heavy-nu on neutrons

    block_factor = 1.0d0 + exp(-eta_e - FERMI5O4(eta_nua))
    zeta_nue_abs_p = 0.0d0 !no absorption of e-nu on protons
    zeta_nua_abs_p = eta_pn * abs_zeta/block_factor
    zeta_nux_abs_p = 0.0d0 !no absorption of heavy-nu on neutrons

    zeta_nue_abs_h = 0.0d0  !no absorption on nuclei
    zeta_nua_abs_h = 0.0d0  !no absorption on nuclei
    zeta_nux_abs_h = 0.0d0  !no absorption on nuclei

    zeta_nue = zeta_nue_abs_n + zeta_nue_abs_p + zeta_nue_abs_h
    zeta_nua = zeta_nua_abs_n + zeta_nua_abs_p + zeta_nua_abs_h
    zeta_nux = zeta_nux_abs_n + zeta_nux_abs_p + zeta_nux_abs_h

    kappa_0_nue = zeta_nue *&
                        (temp * mev_to_erg)**2 * FERMI4O2(eta_nue)
    kappa_0_nua = zeta_nua *&
                        (temp * mev_to_erg)**2 * FERMI4O2(eta_nua)
    kappa_0_nux = zeta_nux *&
                        (temp * mev_to_erg)**2 * FERMI4O2(eta_nux)

    kappa_1_nue = zeta_nue *&
                        (temp * mev_to_erg)**2 * FERMI5O3(eta_nue)
    kappa_1_nua = zeta_nua *&
                        (temp * mev_to_erg)**2 * FERMI5O3(eta_nua)
    kappa_1_nux = zeta_nux *&
                        (temp * mev_to_erg)**2 * FERMI5O3(eta_nux)

    if (isnan(kappa_0_nue).or.isinf(kappa_0_nue)) then
       write(*,*) "Absorption_cgs: NaN/Inf in kappa_0_nue", rho, temp, ye
       Absorption_cgs = -1
    endif

    if (isnan(kappa_0_nua).or.isinf(kappa_0_nua)) then
       write(*,*) "Absorption_cgs: NaN/Inf in kappa_0_nua", rho, temp, ye
       Absorption_cgs = -1
    endif

    if (isnan(kappa_0_nux).or.isinf(kappa_0_nux)) then
       write(*,*) "Absorption_cgs: NaN/Inf in kappa_0_nux", rho, temp, ye
       Absorption_cgs = -1
    endif

    if (isnan(kappa_1_nue).or.isinf(kappa_1_nue)) then
       write(*,*) "Absorption_cgs: NaN/Inf in kappa_1_nue", rho, temp, ye
       Absorption_cgs = -1
    endif

    if (isnan(kappa_1_nua).or.isinf(kappa_1_nua)) then
       write(*,*) "Absorption_cgs: NaN/Inf in kappa_1_nua", rho, temp, ye
       Absorption_cgs = -1
    endif

    if (isnan(kappa_1_nux).or.isinf(kappa_1_nux)) then
       write(*,*) "Absorption_cgs: NaN/Inf in kappa_1_nux", rho, temp, ye
       Absorption_cgs = -1
    endif

    return
END FUNCTION Absorption_cgs


INTEGER FUNCTION Scattering_cgs(rho, temp, ye,&
                                 kappa_0_nue, kappa_0_nua, kappa_0_nux, &
                                 kappa_1_nue, kappa_1_nua, kappa_1_nux)

    use Inf_NaN_Detection
    use table3d_mod
    use weak_constants
    use lk_interpolations

    IMPLICIT NONE

    CCTK_REAL, INTENT(IN) :: rho, temp, ye
    CCTK_REAL, INTENT(OUT) :: kappa_0_nue, kappa_0_nua, kappa_0_nux, &
                              kappa_1_nue, kappa_1_nua, kappa_1_nux

    CCTK_REAL :: zeta_nue_sct_n, zeta_nue_sct_p, zeta_nue_sct_h, &
                 zeta_nua_sct_n, zeta_nua_sct_p, zeta_nua_sct_h, &
                 zeta_nux_sct_n, zeta_nux_sct_p, zeta_nux_sct_h, &
                 zeta_nue, zeta_nua, zeta_nux

    CCTK_REAL :: mass_fact_cgs

    CCTK_REAL :: scttr_cff, scttr_cff_n, scttr_cff_p

    CCTK_REAL :: FERMI5O3, FERMI4O2

    CCTK_REAL :: nb, lrho, ltemp, &
                 mu_n, mu_p, mu_e, &
                 abar, zbar, xp, xn, xh

    CCTK_REAL :: eta_nue, eta_nua, eta_nux, &
                 eta_n, eta_p, eta_hat, eta_e, &
                 eta_np, eta_pn

    Scattering_cgs = 0

#define WEAK_RATES_ITS_ME
#include "inc/weak_rates_guts.inc"
#undef WEAK_RATES_ITS_ME

    !------------------------------------------------------------------------
    !SCATTERING
    !------------------------------------------------------------------------
    !Eqs (A17) with different species plus the A21 which multiplies
    !the mass fraction
    !The formula reported in Rosswog paper is wrong. A term is missing
    !(1+3g_A^2)/4 that is roughly ~1.4.
    !Corrected version from Ruffert et al. paper (TODO update using Barrows 2006)

    scttr_cff_p = nb * ((1.0d0 + 5.0d0 * alpha**2) / 24.0d0) * sigma_0/me_erg**2
    scttr_cff_n = nb * ((4.0d0 * (Cv-1.0d0)**2 + 5.0d0 * alpha**2) / 24.0d0) * &
                  sigma_0/me_erg**2

    !Neutrino nucleon scattering coefficients
    !electron neutrinos
    zeta_nue_sct_n = xn*scttr_cff_n
    zeta_nue_sct_p = xp*scttr_cff_p
    !electron antineutrinos
    zeta_nua_sct_n = xn*scttr_cff_n
    zeta_nua_sct_p = xp*scttr_cff_p
    !tau and mu neutrinos
    zeta_nux_sct_n = xn*scttr_cff_n
    zeta_nux_sct_p = xp*scttr_cff_p


    !Coherent neutrinos nucleus scattering (Shapiro & Teukolsky 1983,
    !sin^2Theta_w has been approx by 0.25)
    !XXX it has only 1 factor of A because zeta multiples the number fraction,
    !not the mass fractions
    scttr_cff = nb * 0.0625d0 * sigma_0 / me_erg**2 * abar * &
                      (1.0d0 - zbar / abar)**2

    !On heavy nuclei
    zeta_nue_sct_h = xh*scttr_cff
    zeta_nua_sct_h = xh*scttr_cff
    zeta_nux_sct_h = xh*scttr_cff


    zeta_nue = zeta_nue_sct_n + zeta_nue_sct_p + zeta_nue_sct_h
    zeta_nua = zeta_nua_sct_n + zeta_nua_sct_p + zeta_nua_sct_h
    zeta_nux = zeta_nux_sct_n + zeta_nux_sct_p + zeta_nux_sct_h

    kappa_0_nue = zeta_nue *&
                        (temp * mev_to_erg)**2 * FERMI4O2(eta_nue)
    kappa_0_nua = zeta_nua *&
                        (temp * mev_to_erg)**2 * FERMI4O2(eta_nua)
    kappa_0_nux = zeta_nux *&
                        (temp * mev_to_erg)**2 * FERMI4O2(eta_nux)

    kappa_1_nue = zeta_nue *&
                        (temp * mev_to_erg)**2 * FERMI5O3(eta_nue)
    kappa_1_nua = zeta_nua *&
                        (temp * mev_to_erg)**2 * FERMI5O3(eta_nua)
    kappa_1_nux = zeta_nux *&
                        (temp * mev_to_erg)**2 * FERMI5O3(eta_nux)

    if (isnan(kappa_0_nue).or.isinf(kappa_0_nue)) then
       write(*,*) "Scattering_cgs: NaN/Inf in kappa_0_nue", rho, temp, ye
       Scattering_cgs = -1
    endif

    if (isnan(kappa_0_nua).or.isinf(kappa_0_nua)) then
       write(*,*) "Scattering_cgs: NaN/Inf in kappa_0_nua", rho, temp, ye
       Scattering_cgs = -1
    endif

    if (isnan(kappa_0_nux).or.isinf(kappa_0_nux)) then
       write(*,*) "Scattering_cgs: NaN/Inf in kappa_0_nux", rho, temp, ye
       Scattering_cgs = -1
    endif

    if (isnan(kappa_1_nue).or.isinf(kappa_1_nue)) then
       write(*,*) "Scattering_cgs: NaN/Inf in kappa_1_nue", rho, temp, ye
       Scattering_cgs = -1
    endif

    if (isnan(kappa_1_nua).or.isinf(kappa_1_nua)) then
       write(*,*) "Scattering_cgs: NaN/Inf in kappa_1_nua", rho, temp, ye
       Scattering_cgs = -1
    endif

    if (isnan(kappa_1_nux).or.isinf(kappa_1_nux)) then
       write(*,*) "Scattering_cgs: NaN/Inf in kappa_1_nux", rho, temp, ye
       Scattering_cgs = -1
    endif

    return
END FUNCTION Scattering_cgs
