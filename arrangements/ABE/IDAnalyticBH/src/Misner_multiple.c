 /*@@
   @file      Misner_multiple.F
   @date      
   @author    Carsten Gundlach
   @desc 
              Set up initial data for multiple Misner black holes
   @enddesc 
   @version   $Header: /cactus/CactusEinstein/IDAnalyticBH/src/Misner_multiple.c,v 1.21 2004/05/17 12:28:56 goodale Exp $
 @@*/

#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "IDAnalyticBH.h"

static const char *rcsid = "$Header: /cactus/CactusEinstein/IDAnalyticBH/src/Misner_multiple.c,v 1.21 2004/05/17 12:28:56 goodale Exp $";

CCTK_FILEVERSION(CactusEinstein_IDAnalyticBH_Misner_multiple_c)

 /*@@
   @routine    Misner_multiple
   @date       
   @author     Carsten Gundlach
   @desc 
               Set up initial data for multiple Misner black holes
   @enddesc 
   @calls      MisnerEvalPsi
   @history
   @hdate Fri Apr 26 10:04:05 2002 @hauthor Tom Goodale
   @hdesc Changed to use new StaticConformal stuff
   @hdate 6.May.2003   @hauthor Jonathan Thornburg
   @hdesc code cleanups
   @endhistory 
@@*/
void Misner_multiple(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int i, npoints;
  CCTK_REAL xval, yval, zval;
  CCTK_REAL inv_psi, tmp0, tmp1, tmp2, tmp3, tmp4;
  const CCTK_REAL nm_eps             = 1e-6;   /* finite differencing step*/
  const CCTK_REAL halved_inv_nm_eps  = 0.5 / nm_eps;
  const CCTK_REAL inv_nm_eps_squared = 1.0 / SQR(nm_eps);
  int make_conformal_derivs = 0;

  CCTK_VInfo(CCTK_THORNSTRING,
             "setting up Misner initial data for %d black holes",
             (int)misner_nbh);

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
"Misner_multiple(): impossible value for conformal_storage=\"%s\"!",
                 conformal_storage);                            /*NOTREACHED*/
    }
  }


  npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

  /*     Initialize C global variables
   *     -----------------------------
   */

  Misner_init(misner_nbh, mu, nmax);

  /*     Get value of psi pointwise from a C function
   *     --------------------------------------------
   */

  for(i = 0; i < npoints; i++)
  {
    xval = x[i];
    yval = y[i];
    zval = z[i];

    tmp0 = MisnerEvalPsi(xval, yval, zval);
    psi[i] = tmp0;
           
    /*              Only calculate derivatives of psi if required
     *              ---------------------------------------------
     */
    if (make_conformal_derivs)
    {
      tmp1 = MisnerEvalPsi(xval+nm_eps,yval,zval);
      tmp2 = MisnerEvalPsi(xval-nm_eps,yval,zval);
      psix[i] = (tmp1-tmp2) * halved_inv_nm_eps;

      if(*conformal_state > 2)
      {
        psixx[i] = (tmp1+tmp2-2.0*tmp0) * inv_nm_eps_squared;
      }

      tmp1 = MisnerEvalPsi(xval,yval+nm_eps,zval);
      tmp2 = MisnerEvalPsi(xval,yval-nm_eps,zval);
      psiy[i]  = (tmp1-tmp2) * halved_inv_nm_eps;

      if(*conformal_state > 2)
      {
        psiyy[i] = (tmp1+tmp2-2.0*tmp0) * inv_nm_eps_squared;
      }

      tmp1 = MisnerEvalPsi(xval,yval,zval+nm_eps);
      tmp2 = MisnerEvalPsi(xval,yval,zval-nm_eps);
      psiz[i] = (tmp1-tmp2) * halved_inv_nm_eps;

      if(*conformal_state > 2)
      {
        psizz[i] = (tmp1+tmp2-2.0*tmp0) * inv_nm_eps_squared;
      }

      if(*conformal_state > 2)
      {
        tmp1 = MisnerEvalPsi(xval+nm_eps,yval+nm_eps,zval);
        tmp2 = MisnerEvalPsi(xval+nm_eps,yval-nm_eps,zval);
        tmp3 = MisnerEvalPsi(xval-nm_eps,yval+nm_eps,zval);
        tmp4 = MisnerEvalPsi(xval-nm_eps,yval-nm_eps,zval);
        psixy[i] = 0.25*(tmp1-tmp2-tmp3+tmp4) * inv_nm_eps_squared;

        tmp1 = MisnerEvalPsi(xval,yval+nm_eps,zval+nm_eps);
        tmp2 = MisnerEvalPsi(xval,yval-nm_eps,zval+nm_eps);
        tmp3 = MisnerEvalPsi(xval,yval+nm_eps,zval-nm_eps);
        tmp4 = MisnerEvalPsi(xval,yval-nm_eps,zval-nm_eps);
        psiyz[i] = 0.25*(tmp1-tmp2-tmp3+tmp4) * inv_nm_eps_squared;

        tmp1 = MisnerEvalPsi(xval+nm_eps,yval,zval+nm_eps);
        tmp2 = MisnerEvalPsi(xval+nm_eps,yval,zval-nm_eps);
        tmp3 = MisnerEvalPsi(xval-nm_eps,yval,zval+nm_eps);
        tmp4 = MisnerEvalPsi(xval-nm_eps,yval,zval-nm_eps);
        psixz[i] = 0.25*(tmp1-tmp2-tmp3+tmp4) * inv_nm_eps_squared;
      }
    }
  }

  /*     Cactus conventions
   *     ------------------
   */

  if (make_conformal_derivs)
  {
    for(i = 0; i < npoints; i++)
    {
      inv_psi = 1.0 / psi[i];

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

  /*     Metric depends on conformal state
   *     ---------------------------------
   */

  if (CCTK_EQUALS(metric_type, "static conformal"))
  {
    for(i = 0; i < npoints; i++)
    {
      gxx[i] = 1.0;
      gyy[i] = 1.0;
      gzz[i] = 1.0;
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

}
