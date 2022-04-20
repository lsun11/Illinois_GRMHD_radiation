#ifndef WEAK_RATES_ITS_ME
#error "This file should not be included by the end user!"
#endif

      !Density is assumed to be in cgs units and
      !the temperature in MeV

      lrho  = log10(rho)
      ltemp = log10(temp)

      !Compute the baryon number density (mass_fact is given in MeV)
       mass_fact_cgs = mass_fact * mev_to_erg / (clite*clite)
       nb = rho / mass_fact_cgs

      !Interpolate the chemical potentials (stored in MeV in the table)
       mu_n = linearInterpolation3d(lrho, ltemp, ye, MU_N)
       mu_p = linearInterpolation3d(lrho, ltemp, ye, MU_P)
       mu_e = linearInterpolation3d(lrho, ltemp, ye, MU_E)

      !Interpolate the fractions
       abar = linearInterpolation3d(lrho, ltemp, ye, ABAR)
       zbar = linearInterpolation3d(lrho, ltemp, ye, ZBAR)
       xp   = linearInterpolation3d(lrho, ltemp, ye, XP)
       xn   = linearInterpolation3d(lrho, ltemp, ye, XN)
       xh   = linearInterpolation3d(lrho, ltemp, ye, XH)

      !Compute the neutrino degeneracy assuming that neutrons and
      !protons chemical potentials DO NOT include the rest mass density
      ! eta_nue = (mu_p + mu_e - mu_n - Qnp) / temp
      !Compute the neutrino degeneracy assuming that neutrons and
      !protons chemical potentials includes the rest mass density
      !This is the correct formula for stellarcollapse.org tables
       eta_nue = (mu_p + mu_e - mu_n) / temp
       eta_nua = -eta_nue
       eta_nux = 0.0d0
       eta_e   = mu_e / temp

      !Neutron and proton degeneracy
       eta_n   = mu_n / temp
       eta_p   = mu_p / temp
      !Difference in the degeneracy parameters without
      !neutron-proton rest mass difference
       eta_hat = eta_n - eta_p - Qnp / temp

      !Janka takes into account the Pauli blocking effect for
      !degenerate nucleons as in Bruenn (1985). Ruffert et al. Eq. (A8)
      !xp = xp / (1.0d0 + 2.0d0 / 3.0d0 * (max(eta_p, 0.0d0)))
      !xn = xn / (1.0d0 + 2.0d0 / 3.0d0 * (max(eta_n, 0.0d0)))

      !Consistency check on the fractions
       xp = max(0.0d0, xp)
       xn = max(0.0d0, xn)
       xh = max(0.0d0, xh)
       abar = max(0.0d0, abar)
       zbar = max(0.0d0, zbar)

      !eta takes into account the nucleon final state blocking
      !(at high density)
       eta_np = nb * (xp-xn) / (exp(-eta_hat) - 1.0d0)
       eta_pn = nb * (xn-xp) / (exp(eta_hat) - 1.0d0)

      !There is no significant defferences between Rosswog (prev. formula)
      !and Janka's prescriptions
      ! eta_np = nb * ((2.0d0 * ye-1.0d0) / (exp(eta_p-eta_n) - 1.0d0))
      ! eta_pn = eta_np * exp(eta_p-eta_n)

     !See Bruenn (ApJSS 58 1985) formula 3.1, non degenerate matter limit.
      if  (rho < 2.0d11) then
        eta_pn = nb * xp
        eta_np = nb * xn
      endif

     !Consistency Eqs (A9) (Rosswog's paper) they should be positive
      eta_pn = max(eta_pn, 0.0d0)
      eta_np = max(eta_np, 0.0d0)

