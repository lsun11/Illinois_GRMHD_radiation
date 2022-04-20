
 /*@@
   @file      Schwarzschild.c
   @date      Sun Oct 17 10:35:41 1999
   @author    Tom Goodale
   @desc
              C version of Scwhwarzschild lapse routine
   @enddesc
   @version   $Id: Schwarzschild.c,v 1.22 2004/07/06 15:53:55 schnetter Exp $
 @@*/

#include <math.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "IDAnalyticBH.h"

static const char *rcsid = "$Header: /cactus/CactusEinstein/IDAnalyticBH/src/Schwarzschild.c,v 1.22 2004/07/06 15:53:55 schnetter Exp $";
CCTK_FILEVERSION(CactusEinstein_IDAnalyticBH_Schwarzschild_c)


void Schwarzschild(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  const CCTK_REAL zero = 0.0, one = 1.0, two = 2.0, three = 3.0;
  CCTK_REAL tmp, r_squared, r_cubed;
  int make_conformal_derivs = 0;
  int i, npoints;

  CCTK_VInfo(CCTK_THORNSTRING,
             "setting up Schwarzschild initial");

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
"Schwarzschild(): impossible value for conformal_storage=\"%s\"!",
                 conformal_storage);                            /*NOTREACHED*/
    }
  }      


  npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

  if(CCTK_EQUALS(metric_type, "static conformal"))
  {
    for (i = 0; i < npoints; i++)
    {
      /*        Compute conformal factor */
      CCTK_REAL const rr = pow(pow(r[i], 4) + pow(epsilon, 4), 0.25);

      psi[i] = (one + mass / (two*rr));

      if(make_conformal_derivs)
      {
        /*        derivatives of psi / psi */
        r_squared = rr*rr;
        r_cubed   = rr*r_squared;
        tmp = mass / (two*r_cubed*psi[i]);


        psix[i] = -x[i]*tmp;
        psiy[i] = -y[i]*tmp;
        psiz[i] = -z[i]*tmp;

        if(*conformal_state > 2)
        {
          tmp = mass/(two*r_squared*r_cubed*psi[i]);
          psixy[i] = three*x[i]*y[i]*tmp;
          psixz[i] = three*x[i]*z[i]*tmp;
          psiyz[i] = three*y[i]*z[i]*tmp;
      
          psixx[i]  = (three*x[i]*x[i] - r_squared)*tmp;
          psiyy[i]  = (three*y[i]*y[i] - r_squared)*tmp;
          psizz[i]  = (three*z[i]*z[i] - r_squared)*tmp;
        }
      }
      gxx[i] = one;
      gyy[i] = one;
      gzz[i] = one;
      gxy[i] = zero;
      gxz[i] = zero;
      gyz[i] = zero;
    }
  }
  else
  {
    for (i = 0; i < npoints; i++)
    {
      CCTK_REAL const rr = pow(pow(r[i], 4) + pow(epsilon, 4), 0.25);

      tmp = one + mass / (two*rr);
      gxx[i] = tmp*tmp*tmp*tmp;
      gyy[i] = gxx[i];
      gzz[i] = gxx[i];
      gxy[i] = zero;
      gxz[i] = zero;
      gyz[i] = zero;
    }
  }
  
  /*     If the initial lapse is not one ... */
  if (CCTK_Equals(initial_lapse,"schwarzschild")
      || CCTK_Equals(initial_lapse,"schwarz"))
  {
    CCTK_INFO("Initialise with Schwarzschild lapse");

    for (i = 0; i < npoints; i++)
    {
      lapm1[i] = (two*r[i] - mass) / (two*r[i] + mass) - 1.0;
    }
  }

  /*     time symmetric initial slice */
  IDAnalyticBH_zero_CCTK_REAL_array(npoints, kxx);
  IDAnalyticBH_zero_CCTK_REAL_array(npoints, kxy);
  IDAnalyticBH_zero_CCTK_REAL_array(npoints, kxz);
  IDAnalyticBH_zero_CCTK_REAL_array(npoints, kyy);
  IDAnalyticBH_zero_CCTK_REAL_array(npoints, kyz);
  IDAnalyticBH_zero_CCTK_REAL_array(npoints, kzz);
}
