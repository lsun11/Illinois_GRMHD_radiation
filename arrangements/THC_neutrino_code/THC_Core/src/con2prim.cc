/*! \file con2prim.cc
\brief Implementation of the classes con2prim, cons_vars, prim_vars, and atmosphere.
\author Wolfgang Kastaun

This contains all the code needed to convert from conserved to primitive variables.
The main method is con2prim::operator(),
which uses con2prim::findroot() to find the
root of a function object con2prim::froot in order to find the solution.
\see  <a HREF="../../doc/con2prim/con2prim.pdf">con2prim.pdf</a>
*/

#include "con2prim.h"
#include <limits>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <gsl/gsl_math.h>

using namespace whizza;
using namespace std;


ostream& operator<<(ostream& os, const con2prim::froot& f)
{
  os << "d = " << f.d << ", q = " << f.q << ", r = " << f.r;
  return os;
}

/** Find the root of the function object in the given interval.
We use the Illinois variant of the Regula Falsi.
This method keeps a bracket \f$ [x_l, x_r] \f$ around the root which is getting smaller each iteration.
Each step, the position of the root is estimated by linear interpolation
\f[ x_b = \frac{ x_l f_r - x_r f_l}{f_r - f_l},\qquad f_l = f(x_l), f_r=f(x_r) \f]
Since the root was bracketed, the new point is inside the interval.
We then evaluate the function at the new point, \f$ f_b = f(x_b) \f$.
If the root is bracketed by the left and the new point, the right
bracket and function value is set to the new point, else the left one.
To speed up convergence, we check if the same side \f$ s \f$ remained unchanged two times in row.
In that case, we pretend the function value there is only half as big,
i.e. \f$ f_s \to \frac{1}{2} f_s \f$.
This has the effect of moving the next best guess closer to the "lazy" side (the Illinois trick),
until it crosses the actual root and provides a better bracketing.

The convergence criterium is outsourced in the froot::converged() method.
The maximum number of iterations is limited by the con2prim::max_iter member.
Note the last evaluated point is always the best guess. We rely on that because
froot remembers intermediate expressions we re-use later.

\returns The location of the root
**/
pz_real con2prim::findroot(froot& f, const pz_real zmin, const pz_real zmax, report& rep) const
{
  enum lazy_side_t {LEFT, RIGHT, NONE};

  lazy_side_t lazy_side = NONE;
  pz_real x_left        = zmin;
  pz_real x_right       = zmax;
  pz_real f_left        = f(x_left);
  rep.iters             = 1;
  if (froot::converged(x_left, f_left, acc)) {
    return x_left;
  }
  pz_real f_right       = f(x_right);
  rep.iters++;
  pz_real x_best        = x_right;
  pz_real f_best        = f_right;
  if (f_left * f_right > 0) {
    stringstream ss;
    ss << "Con2Prim: root not bracketed (zmin = " << x_left << ", zmax = " << x_right << "f(zmin) = "
       << f_left << "f(zmax) = " << f_right << ", " << f << ")";
    rep.err_msg = ss.str();
    rep.failed  = true;
    return -1;
  }
  while (!froot::converged(x_best, f_best, acc))
  {
    x_best = (x_left * f_right - x_right * f_left) / (f_right - f_left);
    f_best = f(x_best);

    if (f_best * f_right > 0)
    {
      x_right   = x_best;
      f_right   = f_best;
      if (lazy_side == LEFT) f_left *= 0.5;
      lazy_side = LEFT;
    }
    else
    {
      x_left    = x_best;
      f_left    = f_best;
      if (lazy_side == RIGHT) f_right *= 0.5;
      lazy_side = RIGHT;
    }

    if (++rep.iters > max_iter) {
      stringstream ss;
      ss << "Con2Prim: root finding did not converge." << "(steps = "
         << rep.iters << ", z_best = " << x_best << ", " << f;
      rep.err_msg = ss.str();
      rep.failed  = true;
      return -1;
    }
  };
  return x_best;
}

/**
Besides trivial bookkeeping, we compute the quantity
\f[ k_L = \frac{2 v_L}{1 + v_L^2},   \qquad
    v_L = \frac{z_L}{\sqrt{1 + z_L^2}} \f]
which is later used in a first check based on conservatives alone
if the momentum is too high.
**/
con2prim::con2prim(const eos_thermal& eos_,  pz_real rho_strict_,
                   pz_real z_lim_, const atmosphere& atmo_, pz_real acc_, int max_iter_,
                   bool _fail_if_eps_lt_min, bool _fail_if_eps_gt_max,
                   bool _fail_if_ye_lt_min,  bool _fail_if_ye_gt_max
                   )
: eos(eos_), rho_strict(rho_strict_), z_lim(z_lim_),
  atmo(atmo_), acc(acc_), max_iter(max_iter_),
  fail_if_eps_lt_min(_fail_if_eps_lt_min),
  fail_if_eps_gt_max(_fail_if_eps_gt_max),
  fail_if_ye_lt_min(_fail_if_ye_lt_min),
  fail_if_ye_gt_max(_fail_if_ye_gt_max)
{
  pz_real v_lim = z_lim / sqrt(1.0 + z_lim*z_lim);
  k_lim         = 2.0 * v_lim / (1.0 + v_lim*v_lim);
}

/**
This sets the parameters and EOS defining the root function.
We also compute the electron fraction limited to the allowed range of the EOS.
**/
con2prim::froot::froot(const eos_thermal& eos_, pz_real ye_,
                       pz_real d_, pz_real q_, pz_real r_)
: eos(eos_.implementation()), rho_range(eos_.range_rho()), d(d_), q(q_), r(r_)
{
  const range ye_range = eos_.range_ye();
  if (ye_ > ye_range.max) {
    ye        = ye_range.max;
    ye_valid  = TOO_LARGE;
  }
  else if (ye_ < ye_range.min) {
    ye        = ye_range.min;
    ye_valid  = TOO_SMALL;
  }
  else {
    ye        = ye_;
    ye_valid  = VALID;
  }
}

/**
From the input \f$ z = W v \f$, we compute
\f{eqnarray*}{
  \tilde{W}(z)     &=& \sqrt{1+z^2} \\
  \tilde{\rho}(z)  &=& \frac{d}{\tilde{W}(z)} \\
  \tilde{\epsilon}(z) &=& \tilde{W}(z) q - z r + \frac{z^2}{1+\tilde{W}(z)}
\f}
We then compute values limited to the allowed EOS range
\f{eqnarray*}{
  \hat{\rho}(z)       &=& \max(\rho_0, \min(\rho_1, \tilde{\rho}(z))) \\
  \hat{\epsilon}(z)   &=& \max(\epsilon_0(\hat{\rho}(z)), \min(\epsilon_1(\hat{\rho}(z)), \tilde{\epsilon}(z)))
\f}
Those we insert into the EOS to compute a pressure
\f$ \hat{P} = P(\hat{\rho}, \hat{\epsilon}, \hat{Y}_e) \f$
and the quantity
\f$ \hat{a} = \frac{\hat{P}}{\hat{\rho}\left( 1 + \hat{\epsilon} \right) } \f$.
We then compute a "half-limited" enthalpy
\f[ \bar{h}(z)    = (1+\tilde{\epsilon}(z))(1 + \hat{a}(z)) \f]
The root function is finally given by
\f[ f(z) = z - \frac{r}{\bar{h}(z)} \f]
For later use, we store density, specific energy, if we had to limit them, and \f$ \bar{h} \f$
in data members.
\returns The root function \f$ f(z) \f$.
**/
pz_real con2prim::froot::operator()(pz_real z)
{
  w_lor           = sqrt(1.0 + z*z);
  rho             = d / w_lor;
  if (rho < rho_range.min) {
    rho = rho_range.min;
    rho_valid = TOO_SMALL;
  }
	else if (rho > rho_range.max) {
    rho = rho_range.max;
    rho_valid = TOO_LARGE;
  }
  else rho_valid = VALID;

  pz_real eps_raw = q * w_lor + z*(z / (w_lor + 1.0) - r);

  const eos_thermal::range reps = eos.range_eps(rho, ye);

  if (eps_raw <= reps.min) {
    eps_valid = TOO_SMALL;
    eps       = reps.min;
  }
  else if (eps_raw >= reps.max) {
    eps_valid = TOO_LARGE;
    eps       = reps.max;
  } else {
    eps_valid = VALID;
    eps       = eps_raw;
  }

  press           = eos.press_from_rho_eps_ye(rho, eps, ye, eos_errs);
  if (eos_errs.failed) { //Should never happen here, ranges are already checked
    string msg= string("con2prim::froot: unexpected EOS failure (") + eos_errs.err_msg + ")";
    throw logic_error(msg);
  }
  pz_real a       = press / (rho * (1.0 + eps));
  h               = (1.0 + eps_raw)*(1.0 + a);

  return z - r / h;
}

/**
Currently, we use the \f$ f(z) \f$ as estimate for the error in \f$ z \f$.
For a wide range of parameters, the slope of the root function is of order one.
**/
bool con2prim::froot::converged(pz_real z, pz_real res, pz_real acc)
{
  return fabs(res) < acc;
  //return fabs(res*z)  <= (1.0+z*z)*acc;
}

void con2prim::set_to_atmo(prim_vars& pv, cons_vars& cv, const pz_real vol_elem, report& rep) const
{
  pz_real tau0    = cv.tau;
  atmo.set(pv, cv, vol_elem);
  rep.dtau        = cv.tau - tau0;
  rep.adjust_cons = true;
  rep.set_atmo    = true;
}

void con2prim::set_to_nan(prim_vars& pv, cons_vars& cv, report& rep) const
{
  rep.adjust_cons = true;
  pz_real nan     = numeric_limits<float>::quiet_NaN();
  pv.rho = pv.eps = pv.ye = pv.press = pv.vel(0) = pv.vel(1) = pv.vel(2) = pv.w_lor = nan;
  cv.dens =cv.tau = cv.tracer_ye = cv.scon(0) = cv.scon(1) = cv.scon(2) = nan;
}


void con2prim::set_fail(prim_vars& pv, cons_vars& cv, report& rep, string msg) const
{
  rep.failed      = true;
  rep.err_msg     = msg;
  set_to_nan(pv, cv, rep);
}

/**
This converts conserved variables to primitives and enforces an artificial atmosphere.
If the input data is invalid or unphysical in a way we consider a fatal error,
conserved and primitives are set to NAN, the errs.failed flag is set, and an error
message is stored in errs.err_msg. In the following this is denoted as "failing".
For unphysical input we decide to repair, the conserved variables are adjusted
and errs.adjust_const flag is set.
If the rest mass density is below the atmo.rho_cut, conserved and primitives are set to
atmosphere and the flag errs.set_atmo is set.
For any adjustment, including atmosphere enforcement, errs.dtau is set to the
resulting change in conserved energy.

We first check if the metric volume element is a finite positive number.
If not, fail.
Then we check if the conserved density (divided by the volume element) is below
the atmosphere cutoff. If yes, the rest mass density must be smaller as well,
and we set to atmosphere.
Next, we compute
\f[ r_i = \frac{1}{D} S_i, \qquad r^i = g^{ij} r_j \f]
\f[ r   = \sqrt{r_i r^i}, \qquad
    q   = \frac{\tau}{D}, \qquad
    Y_e = \frac{Y_e^T}{D} \f]
where \f$ Y_e^T \f$ is the conserved tracer variable for the electron fraction.
We then check if those are finite numbers, else we fail.
Now we enforce bounds energy and momentum have to satisfy for physical input.
Then we set up the root function defined by the above scalars (see froot)
and compute the initial bracketing of the root from
\f[
  z_{-} = \frac{k/2}{\sqrt{1 - k^2/4}} ,\qquad
  z_{+} = \frac{k}{\sqrt{1 - k^2}}     ,\qquad
  k \equiv \frac{r}{1+q}
\f]
We then search for the root using con2prim::findroot().
There is always a solution in the above interval, but the root finding algorithm could still fail
by converging too slowly. In that case, we fail.
Otherwise we check if the density is below the atmosphere cut,
in which case we set to the artificial atmosphere.
Next, we check if the solution is physical.
First, if the density is above the EOS range, we fail.
If the specific internal energy had to be limited because it was too small,
we always allow it, since this happens frequently when evolving at zero temperature.
If it was too big, we fail if the density is above rho_strict, else we use the limited value.
If the electron fraction is out of range, we fail if the density is above rho_strict,
else we use the limited value.
If the velocity is above our speed limit, we fail for densities above rho_strict.
Else we rescale the velocity, compute a new Lorentz factor, and recompute the rest mass density
keeping the conserved density fixed. We also limit the specific energy to the range
corresponding to the changed rest mass density.

Finally, if any adjustments were applied before, we recompute the conserved variables
from the adjusted primitives.
The primitive variables are set using the values stored in the froot function object,
which remembers intermediate results from the last evaluation.
The 3-velocity is computed from
\f[ v^i = \frac{r^i}{\tilde{W} \bar{h}} \f]
see froot::operator() for definitions of \f$ \tilde{W}, \bar{h} \f$.

\see <a HREF="../../doc/con2prim/con2prim.pdf">con2prim.pdf</a>
**/
void con2prim::operator()(prim_vars& pv, cons_vars& cv, const pz_metric3& g, report& errs) const
{
  errs.dtau         = 0.0;
  errs.iters        = 0;
  errs.adjust_cons  = false;
  errs.set_atmo     = false;
  errs.failed       = false;

  if ((0==gsl_finite(g.vol_elem)) || (g.vol_elem <= 0)) {
    stringstream ss;
    ss << "Con2Prim: problem with metric (det(gij)= " << g.det << ", dV= " << g.vol_elem  << ")";
    set_fail(pv, cv, errs, ss.str());
    return;
  }

  const pz_real d = cv.dens / g.vol_elem;

  if (d <= atmo.rho_cut) {
    set_to_atmo(pv, cv, g.vol_elem, errs);
    return;
  }

  const pz_vec3l rl   = cv.scon / cv.dens;
  pz_vec3u ru         = g.raise(rl);
  pz_real r           = sqrt(ru*rl);
  pz_real q           = cv.tau / cv.dens;
  const pz_real ye0   = cv.tracer_ye / cv.dens;

  if ((0==gsl_finite(d)) || (0==gsl_finite(r))|| (0==gsl_finite(q)) || (0==gsl_finite(ye0))) {
    stringstream ss;
    ss << "Con2Prim: input contains NaNs/INFs (d =  " << d
       << ", r = " << r << ", q = " << q << ", ye = " << ye0 << endl;
    set_fail(pv, cv, errs, ss.str());
    return;
  }


  if (q < 0) {
    //cerr << "adj q=" << q <<endl;
    errs.adjust_cons    = true;
    q = 0;
  }
  const pz_real r_lim = k_lim * (1.0 + q);

  if (r > r_lim) {
    //cerr << "adj r=" <<r << ", rlim = " << r_lim << endl;
    errs.adjust_cons    = true;
    ru                 *= (r_lim/r);
    r                   = r_lim;
  }

  const pz_real k     = r / (1.0 + q);
  const pz_real z_min = (k / 2.0) / sqrt(1.0 - k*k / 4.0);
  // analytic bracket plus safety margin needed because of roundoff errors.
  const pz_real z_max = 1e-6 + k / sqrt(1.0 - k*k);

  froot f(eos, ye0, d, q, r);
  const pz_real z     = findroot(f, z_min, z_max, errs);
  if (errs.failed) {
    set_to_nan(pv, cv, errs);
    return;
  }

  if (f.rho < atmo.rho_cut) {
    set_to_atmo(pv, cv, g.vol_elem, errs);
    return;
  }
  if (f.rho_valid == froot::TOO_LARGE) {
    stringstream ss;
    ss << "Con2Prim: rho > rho_max  (" << f << ")";
    set_fail(pv, cv, errs, ss.str());
    return;
  }

  pv.rho    = f.rho;
  pv.eps    = f.eps;
  pv.ye     = f.ye;
  pv.press  = f.press;
  pv.vel    = ru / (f.w_lor * f.h);
  pv.w_lor  = f.w_lor;

  if (f.eps_valid != froot::VALID) {
    //cerr << "eps_valid=" << f.eps_valid << "rho="<<f.rho <<", rhostrict= "<< rho_strict <<endl;
    if((f.rho >= rho_strict) &&
        (
          (fail_if_eps_lt_min && f.eps_valid == froot::TOO_SMALL) ||
          (fail_if_eps_gt_max && f.eps_valid == froot::TOO_LARGE)
        )) {
      stringstream ss;
      ss << "Con2Prim: eps out of range (rho = "<< f.rho << ", " << f << ")";
      set_fail(pv, cv, errs, ss.str());
      return;
    }
    errs.adjust_cons = true;
  }
  if (f.ye_valid != froot::VALID) {
    if((f.rho >= rho_strict) &&
        (
          (fail_if_ye_lt_min && f.ye_valid == froot::TOO_SMALL) ||
          (fail_if_ye_gt_max && f.ye_valid == froot::TOO_LARGE)
        )) {
      stringstream ss;
      ss << "Con2Prim: ye out of range (rho = " << f.rho << ", " << f <<")";
      set_fail(pv, cv, errs, ss.str());
      return;
    }
    errs.adjust_cons = true;
  }

  if (z >= z_lim) {
    //cerr << "adj z=" << z <<endl;
    pv.w_lor      = sqrt(1.0 + z_lim*z_lim);
    pv.rho        = d / pv.w_lor;
    pv.vel        = ru * (z_lim / (pv.w_lor * r));
    if (pv.rho >= rho_strict) {
      stringstream ss;
      ss << "Con2Prim: speed limit exceeded (rho = " << pv.rho << ", z=" << z << ", " << f << ")";
      set_fail(pv, cv, errs, ss.str());
      return;
    }
    range rga     = eos.range_eps(pv.rho, pv.ye);
    pv.eps        = (pv.eps <= rga.min) ? rga.min : ((pv.eps>=rga.max) ? rga.max : pv.eps);
    eos_thermal::status eos_errs;
    pv.press      = eos.press_from_rho_eps_ye(pv.rho, pv.eps, pv.ye, eos_errs);
    assert(!eos_errs.failed);
    errs.adjust_cons = true;
  }
  if (errs.adjust_cons) {
    pz_real tau0 = cv.tau;
    cv.from_prim(pv, g);
    errs.dtau    = cv.tau - tau0;
  }
}

void prim_vars::collect(pz_real rho_, pz_real eps_, pz_real ye_, pz_real press_,
               pz_real velx_, pz_real vely_, pz_real velz_, pz_real w_lor_)
{
  rho     = rho_;
  eps     = eps_;
  ye      = ye_;
  press   = press_;
  vel(0)  = velx_;
  vel(1)  = vely_;
  vel(2)  = velz_;
  w_lor   = w_lor_;
}

void prim_vars::scatter(pz_real& rho_, pz_real& eps_, pz_real& ye_, pz_real& press_,
               pz_real& velx_, pz_real& vely_, pz_real& velz_, pz_real& w_lor_) const
{
  rho_    = rho;
  eps_    = eps;
  ye_     = ye;
  press_  = press;
  velx_   = vel(0);
  vely_   = vel(1);
  velz_   = vel(2);
  w_lor_  = w_lor;
}

/**
The conserved variables are computed from
\f{eqnarray*}{
  D     &=& V_c W \rho \\
  S_i   &=& V_c W^2 \rho h v_i \\
  \tau  &=& V_c W^2 \left( v^2 \left( \rho \frac{W}{1 + W} + P \right) + \rho \epsilon \right) \\
  Y_e^T &=& D Y_e
\f}
where \f$ V_c = \sqrt{\det(g_{ij})} \f$ is the 3-metric volume element.
Note the formula for \f$ \tau \f$ is written in a way which is accurate also in
the Newtonian limit, i.e. for \f$ v^2 \ll 1, \epsilon \ll 1 \f$.
**/
void cons_vars::from_prim(const prim_vars& pv, const pz_metric3& g)
{
  const pz_vec3l vl   = g.lower(pv.vel);
  const pz_real hrho  = pv.rho * (1.0 + pv.eps) + pv.press;
  const pz_real v2    = vl * pv.vel;
	const pz_real wl2   = pv.w_lor * pv.w_lor;
  dens                = g.vol_elem * pv.w_lor * pv.rho;
  scon                = (g.vol_elem * wl2 * hrho) * vl;
  tau                 = g.vol_elem * wl2 * (v2*(pv.rho * pv.w_lor / (1.0 + pv.w_lor)
                          + pv.press) + pv.rho * pv.eps);
  tracer_ye           = dens * pv.ye;
}

/**
Like cons_vars::from_prim(), but we do not use the Lorentz factor from the primitives.
Instead, we recompute it from the velocity and set it in the primitives.
**/
void cons_vars::from_prim_part(prim_vars& pv, const pz_symt3l& glo, const double vol_elem)
{
  const pz_vec3l vl   = glo * pv.vel;
  const pz_real hrho  = pv.rho * (1.0 + pv.eps) + pv.press;
  const pz_real v2    = vl * pv.vel;
  pv.w_lor            = 1.0 / sqrt(1.0 - v2);
	const pz_real wl2   = pv.w_lor * pv.w_lor;
  dens                = vol_elem * pv.w_lor * pv.rho;
  scon                = (vol_elem * wl2 * hrho) * vl;
  tau                 = vol_elem * wl2 * (v2*(pv.rho * pv.w_lor / (1.0 + pv.w_lor)
                          + pv.press) + pv.rho * pv.eps);
  tracer_ye           = dens * pv.ye;
}


void cons_vars::collect(pz_real dens_, pz_real tau_, pz_real tracer_ye_,
               pz_real sconx_, pz_real scony_, pz_real sconz_)
{
  dens      = dens_;
  tau       = tau_;
  tracer_ye = tracer_ye_;
  scon(0)   = sconx_;
  scon(1)   = scony_;
  scon(2)   = sconz_;
}

void cons_vars::scatter(pz_real& dens_, pz_real& tau_, pz_real& tracer_ye_,
               pz_real& sconx_, pz_real& scony_, pz_real& sconz_) const
{
  dens_       = dens;
  tau_        = tau;
  tracer_ye_  = tracer_ye;
  sconx_      = scon(0);
  scony_      = scon(1);
  sconz_      = scon(2);
}


atmosphere::atmosphere(pz_real rho_, pz_real eps_, pz_real ye_, pz_real press_, pz_real rho_cut_)
: rho(rho_), eps(eps_), ye(ye_), press(press_), rho_cut(rho_cut_) {}

void atmosphere::set(prim_vars& pv) const
{
  pv.rho    = rho;
  pv.eps    = eps;
  pv.ye     = ye;
  pv.press  = press;
  pv.vel    = ZERO;
  pv.w_lor  = 1.0;
}

void atmosphere::set(cons_vars& cv, pz_real vol_elem) const
{
  cv.dens       = vol_elem * rho;
  cv.tau        = cv.dens * eps;
  cv.tracer_ye  = cv.dens * ye;
  cv.scon       = ZERO;
}

void atmosphere::set(prim_vars& pv, cons_vars& cv, const pz_real vol_elem) const
{
  set(pv);
  set(cv, vol_elem);
}


