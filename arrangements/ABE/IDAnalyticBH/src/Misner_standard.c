 /*@@
   @file      Misner_standard.c
   @date      March 1997
   @author    Joan Masso
   @desc
              Set up initial data for two Misner black holes
   @enddesc
   @history
   @hdate Sun Oct 17 11:05:48 1999 @hauthor Tom Goodale
   @hdesc Converted to C
   @endhistory
   @version   $Header: /cactus/CactusEinstein/IDAnalyticBH/src/Misner_standard.c,v 1.24 2004/05/17 12:28:57 goodale Exp $
 @@*/

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "IDAnalyticBH.h"

static const char *rcsid = "$Header: /cactus/CactusEinstein/IDAnalyticBH/src/Misner_standard.c,v 1.24 2004/05/17 12:28:57 goodale Exp $";

CCTK_FILEVERSION(CactusEinstein_IDAnalyticBH_Misner_standard_c)

 /*@@
   @routine    Misner_standard
   @date
   @author     Joan Masso, Ed Seidel
   @desc
            Initialize the metric with a time symmetrical
            black hole spacetime containing
            two axially symmetric misner black holes with a
            mass/length parameter mu. The mass is computed.
            The spacetime line element has the form:
                 $$ ds^2 = -dt^2 + \Psi^4 (dx^2+dy^2+dz^2) $$
            and only $\Psi$ differs.
            (Conformal factor from Karen Camarda)
   @enddesc

   @par      mu
   @pdesc    Misner parameter.
   @ptype    real
   @pcomment Values less than 1.8 do not really correspond to two
            black holes, as there is an initial single event horizon
            surrounding the throats. So, with low values of mu we also
            have distorted single black holes.
   @endpar

   @par      nmax
   @pdesc    Summation limit for the misner series in the 'twobh' case.
   @ptype    integer
   @endpar
   @history
   @hdate Fri Apr 26 10:04:05 2002 @hauthor Tom Goodale
   @hdesc Changed to use new StaticConformal stuff
   @endhistory 

@@*/
void Misner_standard(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int i, n;
  int npoints;
  CCTK_REAL *csch, *coth, inv_r1, inv_r2;
  CCTK_REAL x_squared, y_squared, xy_squared;
  CCTK_REAL inv_r1_cubed, inv_r2_cubed;
  CCTK_REAL inv_r1_5, inv_r2_5;
  CCTK_REAL inv_psi;
  CCTK_INT powfac;
  CCTK_REAL adm_mass;
  const CCTK_REAL zero = 0.0, one = 1.0, three = 3.0;
  int make_conformal_derivs = 0;

  CCTK_VInfo(CCTK_THORNSTRING,
             "setting up Misner initial data");

  /* Check if we should create and store conformal factor stuff */
  if(CCTK_EQUALS(metric_type, "static conformal"))
  {
    if      (CCTK_EQUALS(conformal_storage,"factor"))
    {
      *conformal_state = 1;
      make_conformal_derivs = 0;
    }
    else if (CCTK_EQUALS(conformal_storage,"factor+derivs"))
    {
      *conformal_state = 2;
      make_conformal_derivs = 1;
    }
    else if (CCTK_EQUALS(conformal_storage,"factor+derivs+2nd derivs"))
    {
      *conformal_state = 3;
      make_conformal_derivs = 1;
    }
    else
    {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
"Misner_Standard(): impossible value for conformal_storage=\"%s\"!",
                 conformal_storage);                            /*NOTREACHED*/
    }
  }      

  /* total number of points on this processor */
  npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];


  /*     Initialize so we can accumulate
   *     -------------------------------
   */
  if (make_conformal_derivs)
  {
    IDAnalyticBH_zero_CCTK_REAL_array(npoints, psix);
    IDAnalyticBH_zero_CCTK_REAL_array(npoints, psiy);
    IDAnalyticBH_zero_CCTK_REAL_array(npoints, psiz);

    if(*conformal_state > 2)
    {
      IDAnalyticBH_zero_CCTK_REAL_array(npoints, psixx);
      IDAnalyticBH_zero_CCTK_REAL_array(npoints, psixy);
      IDAnalyticBH_zero_CCTK_REAL_array(npoints, psixz);
      IDAnalyticBH_zero_CCTK_REAL_array(npoints, psiyy);
      IDAnalyticBH_zero_CCTK_REAL_array(npoints, psiyz);
      IDAnalyticBH_zero_CCTK_REAL_array(npoints, psizz);
    }
  }

  csch = (CCTK_REAL *) malloc (2 * (nmax + 1) * sizeof (CCTK_REAL));
  coth = csch + nmax + 1;

  /*     compute the ADM mass
   *     --------------------
   */
  adm_mass = zero;
  for(n = 1; n <= nmax; n++)
  {
    csch[n] = one / sinh(mu*n);
    coth[n] = one / tanh(mu*n);
    adm_mass   += 4.0 * csch[n];
  }
  CCTK_VInfo(CCTK_THORNSTRING, "ADM mass is %f", (double) adm_mass);


  for(i = 0; i < npoints; i++)
  {
    psi [i] = one;

    x_squared  = SQR(x[i]);
    y_squared  = SQR(y[i]);
    xy_squared = x_squared + y_squared;

    for(n = nmax; n >= 1; n--)
    {
      inv_r1 = one / sqrt(xy_squared+SQR(z[i]+coth[n]));
      inv_r2 = one / sqrt(xy_squared+SQR(z[i]-coth[n]));

      psi[i] += csch[n]*(inv_r1 + inv_r2);

      if (make_conformal_derivs)
      {
        inv_r1_cubed = inv_r1 * inv_r1 * inv_r1;
        inv_r2_cubed = inv_r2 * inv_r2 * inv_r2;
        inv_r1_5     = inv_r1 * inv_r1 * inv_r1_cubed;
        inv_r2_5     = inv_r2 * inv_r2 * inv_r2_cubed;
        psix[i]  +=  -x[i] * (inv_r2_cubed + inv_r1_cubed) * csch[n];
        psiy[i]  +=  -y[i] * (inv_r2_cubed + inv_r1_cubed) * csch[n];
        psiz[i]  +=  (-(z[i]-coth[n])*inv_r2_cubed - (z[i]+coth[n])*inv_r1_cubed) * csch[n];


        if(*conformal_state > 2)
        {
          psixx[i] +=  (three*x_squared*(inv_r1_5 + inv_r2_5)
                        - inv_r1_cubed - inv_r2_cubed) * csch[n];
          psixy[i] +=  three*x[i]*y[i]*(inv_r1_5 + inv_r2_5) * csch[n];
          psixz[i] +=  (three*x[i]*(z[i] - coth[n])*inv_r2_5
                        + three*x[i]*(z[i] + coth[n])*inv_r1_5) * csch[n];
          psiyy[i] +=  (three*y_squared*(inv_r1_5 + inv_r2_5)
                        - inv_r1_cubed - inv_r2_cubed) * csch[n];
          psiyz[i] +=  (three*y[i]*(z[i] - coth[n])*inv_r2_5
                      + three*y[i]*(z[i] + coth[n])*inv_r1_5) * csch[n];
          psizz[i] += (-inv_r2_cubed+three*SQR(z[i] - coth[n])*inv_r2_5
                       + three*SQR(z[i] + coth[n])*inv_r1_5 - inv_r1_cubed) * csch[n];
        }
      }
    }

    /*     Cactus convention
     *     -----------------
     */
    if (make_conformal_derivs)
    {
      inv_psi = one / psi[i];

      psix[i]  *= inv_psi;
      psiy[i]  *= inv_psi;
      psiz[i]  *= inv_psi;

      if(*conformal_state > 2)
      {
        psixx[i] *= inv_psi;
        psixy[i] *= inv_psi;
        psixz[i] *= inv_psi;
        psiyy[i] *= inv_psi;
        psiyz[i] *= inv_psi;
        psizz[i] *= inv_psi;
      }
    }
  }

  /*     Should initialize lapse to Cadez value if possible
   *     --------------------------------------------------
   */

  if (CCTK_Equals(initial_lapse,"cadez"))
  {
    CCTK_INFO("Initialise with cadez lapse");

    for(i = 0; i < npoints; i++)
    {
      xy_squared = SQR(x[i]) + SQR(y[i]);

      lapm1[i] = one-1.0;

      powfac = 1;

      for(n = 1; n <= nmax; n++)
      {
        inv_r1 = one / sqrt(xy_squared+SQR(z[i]+coth[n]));
        inv_r2 = one / sqrt(xy_squared+SQR(z[i]-coth[n]));
        powfac = -powfac;

        lapm1[i] += powfac * csch[n] * (inv_r1 + inv_r2);
      }

      lapm1[i] /= psi[i];
    }
  }

  /*     Metric depends on conformal state
   *     ---------------------------------
   */

  if (CCTK_EQUALS(metric_type, "static conformal"))
  {
    for(i = 0; i < npoints; i++)
    {
      gxx[i] = one;
      gyy[i] = one;
      gzz[i] = one;
    }
  }
  else
  {
    for(i = 0; i < npoints; i++)
    {
      gxx[i] = psi[i] * psi[i] * psi[i] * psi[i];
      gyy[i] = gxx[i];
      gzz[i] = gxx[i];
    }
  }
  IDAnalyticBH_zero_CCTK_REAL_array(npoints, gxy);
  IDAnalyticBH_zero_CCTK_REAL_array(npoints, gxz);
  IDAnalyticBH_zero_CCTK_REAL_array(npoints, gyz);

  /*     Time-symmetric data
   *     -------------------
   */
  IDAnalyticBH_zero_CCTK_REAL_array(npoints, kxx);
  IDAnalyticBH_zero_CCTK_REAL_array(npoints, kxy);
  IDAnalyticBH_zero_CCTK_REAL_array(npoints, kxz);
  IDAnalyticBH_zero_CCTK_REAL_array(npoints, kyy);
  IDAnalyticBH_zero_CCTK_REAL_array(npoints, kyz);
  IDAnalyticBH_zero_CCTK_REAL_array(npoints, kzz);

  if (csch)
  {
    free (csch);
  }
}
