/*! \file con2prim.h
\brief Class definitions for conservative to primitive framework.
\author Wolfgang Kastaun

This header defines a function object con2prim used  to convert conservative
to primitive variables and enforce an artificial atmosphere.
This header also defines a class prim_vars representing primitive variables,
a class cons_vars representing conserved variables, and a class atmosphere
describing an artificial atmosphere.
*/

#ifndef CON2PRIM_H
#define CON2PRIM_H

#include "smtensor.h"
#include "eos_thermal.h"
#include <string>

namespace whizza {

typedef sm_tensor1<pz_real, 3, true> pz_vec3u;
typedef sm_tensor1<pz_real, 3, false> pz_vec3l;
typedef sm_tensor2_sym<pz_real, 3, false, false> pz_symt3l;
typedef sm_metric<pz_real, 3> pz_metric3;


/// Structure to represent the primitive variables needed for evolution
struct prim_vars {
  pz_real rho;      ///< Rest mass density
  pz_real eps;      ///< Specific internal energy
  pz_real ye;       ///< Electron fraction
  pz_real press;    ///< Pressure
  pz_vec3u vel;     ///< 3-velocity \f$ v^i \f$
  pz_real w_lor;    ///< Lorentz factor

  ///Convenience method to set all members at once
  void collect(pz_real rho_, pz_real eps_, pz_real ye_, pz_real press_,
               pz_real velx_, pz_real vely_, pz_real velz_, pz_real w_lor_);

  ///Convenience method to copy all members into single variables
  void scatter(pz_real& rho_, pz_real& eps_, pz_real& ye_, pz_real& press_,
               pz_real& velx_, pz_real& vely_, pz_real& velz_, pz_real& w_lor_) const;
};

/// Structure to represent the evolved conserved variables
struct cons_vars {
  pz_real dens;       ///< Conserved density \f$ D \f$
  pz_real tau;        ///< Conserved energy \f$ \tau \f$
  pz_real tracer_ye;  ///< Conserved tracer \f$ Y_e^T \f$ for the electron fraction
  pz_vec3l scon;      ///< Conserved momentum \f$ S_i \f$

  /// Compute conserved variables from primitives and 3-metric
  void from_prim(const prim_vars& pv, const pz_metric3& g);

  /// Compute conserved variables from partially initialized primitives.
  void from_prim_part(prim_vars& pv, const pz_symt3l& glo, const double vol_elem);

  ///Convenience method to set all members at once
  void collect(pz_real dens_, pz_real tau_, pz_real tracer_ye_,
               pz_real sconx_, pz_real scony_, pz_real sconz_);

  ///Convenience method to copy all members into single variables
  void scatter(pz_real& dens_, pz_real& tau_, pz_real& tracer_ye_,
               pz_real& sconx_, pz_real& scony_, pz_real& sconz_) const;

};

/// Class representing an artificial atmosphere.
struct atmosphere {
  const pz_real rho;      ///< Rest mass density of the atmosphere
  const pz_real eps;      ///< Specific internal energy
  const pz_real ye;       ///< Electron fraction
  const pz_real press;    ///< Pressure
  const pz_real rho_cut;  ///< Cutoff density below which we set to atmosphere

  /// Constructor.
  atmosphere(pz_real rho_, pz_real eps_, pz_real ye_, pz_real press_, pz_real rho_cut_);

  /// Set primitive vars to atmosphere.
  void set(prim_vars& pv) const;

  /// Set conserved vars to atmosphere.
  void set(
    cons_vars& cv,          ///< Conserved vars
    const pz_real vol_elem  ///< Volume element of the 3-metric
  ) const;

  /// Set primitives and conserved variables to atmosphere.
  void set(prim_vars& pv, cons_vars& cv, const pz_real vol_elem) const;
};


/// Class representing conservative to primitive conversion
/**
Besides the actual conversion, this also defines a policy how to handle unphysical
conserved variables, i.e. conserved variables which do not correspond to any sane
primitive variables. Sane means that the velocity is below some speed limit below
the speed of light, that density, specific internal energy, as well as electron
fraction are in the validity range of the EOS.
We distinguish completely valid input, unphysical input that may be adjusted, and
unphysical input that is considered an error in the evolution.
The first two lead to fully consistent conserved and primitive variables.
In the last case, conservative and primitives are always explicitly set to NAN.
Finally, this class supports the enforcing an artificial atmosphere below a certain
density.
**/
class con2prim {
  public:
	typedef eos_thermal::range range;

  /// Information about adjustments and errors during con2prim
  struct report {
    pz_real dtau;         ///< Change in conserved energy due to adjustments (including atmosphere enforcement).
    int iters;            ///< Number of calls to the EOS needed for the root finding.
    bool failed;          ///< If the input was invalid according to the error policy.
    bool adjust_cons;     ///< Whether the conserved variables were adjusted.
    bool set_atmo;        ///< Whether the artificial atmosphere was enforced.
    std::string err_msg;  ///< In case failed==true, this contains an error message (else undefined).

    /// Default constructor
    report() : dtau(0), iters(0), failed(false), adjust_cons(false), set_atmo(false) {}
  };

  /// Function object representing the root function.
  /** This contains all the fixed parameters defining the function.
      It also remembers intermediate results from the last evaluation,
      which we do not want to recompute for performance reasons.
      In particular, it remembers if specific internal energy and density
      were in the validity region of the EOS, or if they needed to be limited.
  **/
  struct froot {
    enum validity {VALID, TOO_LARGE, TOO_SMALL};

    const eos_thermal_impl& eos;      ///< The EOS.
    pz_real ye;                       ///< Electron fraction, limited to valid range (fixed parameter).
    const range rho_range;            ///< Valid density interval of the EOS.
    const pz_real d;                  ///< Fixed parameter \f$ d = \frac{D}{\sqrt{\det(g_{ij})}} \f$
    const pz_real q;                  ///< Fixed parameter \f$ q = \frac{\tau}{D} \f$
    const pz_real r;                  ///< Fixed parameter \f$ r = \frac{ \sqrt{S_i S^i} }{D} \f$
    pz_real w_lor;                    ///< Last evaluation: Lorentz factor.
    pz_real rho;                      ///< Last evaluation: density.
    pz_real eps;                      ///< Last evaluation: specific internal	energy.
    pz_real press;                    ///< Last evaluation: pressure.
    pz_real h;                        ///< Last evaluation: \f$ \bar{h}(z) \f$
    validity eps_valid;               ///< Last evaluation: if \f$ \epsilon \f$ was in the EOS validity range.
    validity rho_valid;               ///< Last evaluation: if \f$ \rho \f$ was in the EOS validity range.
    validity ye_valid;                ///< If initial \f$ Y_e \f$ was in the EOS validity range.
    eos_thermal::status eos_errs;     ///< Last evaluation: EOS error status

    /// Constructor
    froot(
      const eos_thermal& eos_,        ///< The EOS
      pz_real ye_,                    ///< Electron fraction
      pz_real d_,                     ///< \f$ d = \frac{D}{\sqrt{\det(g_{ij})}} \f$
      pz_real q_,                     ///< \f$ q = \frac{\tau}{D} \f$
      pz_real r_                      ///< \f$ r = \frac{ \sqrt{S_i S^i} }{D} \f$
    );

    /// The root function
    pz_real operator()(pz_real z);

    /// The convergence criterion for root finding.
    static bool converged(pz_real z, pz_real res, pz_real acc);
  };

  /// Constructor. The parameters define the fixed parameters of the conversion (EOS, error handling, atmosphere).
  con2prim(
    const eos_thermal& eos_,   ///< The EOS.
    pz_real rho_strict_,       ///< Allow certain adjustments only below this density
    pz_real z_lim_,            ///< Maximum allowed velocity in terms of \f$ z = W v \f$
    const atmosphere& atmo_,   ///< Artificial atmosphere settings
    pz_real acc_,              ///< Accuracy for root finding.
    int max_iter_,             ///< Maximum allowed iterations for root finding.
    bool fail_if_eps_lt_min = false,  ///< Fail if \f$\rho > \rho_{\mathrm{strict}}\f$
                                      ///  and \f$\eps < \eps_{\min}\f$
    bool fail_if_eps_gt_max = true,   ///< Fail if \f$\rho > \rho_{\mathrm{strict}}\f$
                                      ///  and \f$\eps > \eps_{\max}\f$
    bool fail_if_ye_lt_min = true,    ///< Fail if \f$\rho > \rho_{\mathrm{strict}}\f$
                                      ///  and \f$Y_e < (Y_e)_{\min}\f$
    bool fail_if_ye_gt_max = true     ///< Fail if \f$\rho > \rho_{\mathrm{strict}}\f$
                                      ///  and \f$Ye > (Y_e)_{\max}\f$
  );

  /// Convert from conserved to primitive variables
  void operator()(prim_vars& pv, cons_vars& cv, const pz_metric3& gu, report& rep) const;

  private:

  const eos_thermal eos;          ///< The EOS.
  const pz_real rho_strict;       ///< Allow certain adjustments only below this density
  const pz_real z_lim;            ///< Maximum allowed velocity in terms of \f$ z = W v \f$
  pz_real k_lim;                  ///< \f$ k \f$ corresponding to maximum allowed velocity.
  const atmosphere atmo;          ///< Artificial atmosphere settings
  const pz_real acc;              ///< Accuracy for root finding.
  const int max_iter;             ///< Maximum allowed iterations for root finding.

  bool fail_if_eps_lt_min;
  bool fail_if_eps_gt_max;
  bool fail_if_ye_lt_min;
  bool fail_if_ye_gt_max;

  /// Set primitives and conserved to atmosphere and consider this in the report.
  void set_to_atmo(prim_vars& pv, cons_vars& cv, const pz_real vol_elem, report& rep) const;

  /// Set primitives and conserved to NaN and note failure in the report.
  void set_to_nan(prim_vars& pv, cons_vars& cv, report& rep) const;

  /// Handle failure: set everything to NaN, note failure and set error message in report
  void set_fail(prim_vars& pv, cons_vars& cv, report& rep, std::string msg) const;

  /// Find the root in the given interval.
  pz_real findroot(
    froot& f,           ///< Function object for which to find the root
    const pz_real zmin, ///< Left initial bracket.
    const pz_real zmax, ///< Right initial bracket
    report& rep         ///< Store information on number of iterations or errors.
  ) const;
};

}


#endif

