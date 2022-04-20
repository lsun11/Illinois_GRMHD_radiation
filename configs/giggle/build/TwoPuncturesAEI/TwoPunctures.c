/* TwoPunctures:  File  "TwoPunctures.c"*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "TP_utilities.h"
#include "TwoPunctures.h"



/* Swap two variables */
static inline
void swap (CCTK_REAL * restrict const a, CCTK_REAL * restrict const b)
{
  CCTK_REAL const t = *a; *a=*b; *b=t;
}
#undef SWAP
#define SWAP(a,b) (swap(&(a),&(b)))

/* -------------------------------------------------------------------*/
void
TwoPuncturesAEI (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  * mp = par_m_plus;
  * mm = par_m_minus;

  enum GRID_SETUP_METHOD { GSM_Taylor_expansion, GSM_evaluation };
  enum GRID_SETUP_METHOD gsm;


  int antisymmetric_lapse, averaged_lapse, pmn_lapse, brownsville_lapse;

  int const nvar = 1, n1 = npoints_A, n2 = npoints_B, n3 = npoints_phi;

  int imin[3], imax[3];
  int const ntotal = n1 * n2 * n3 * nvar;
#if 0
  int percent10 = 0;
#endif
  static CCTK_REAL *F = NULL;
  static derivs u, v, cf_v;
  //  double *cf_v = (double *)malloc(sizeof(double)*ntotal);
  CCTK_REAL admMass, aux;

  if (! F) {
    /* Solve only when called for the first time */
    F = dvector (0, ntotal - 1);
    allocate_derivs (&u, ntotal);
    allocate_derivs (&v, ntotal);
    allocate_derivs (&cf_v, ntotal);

    if (use_sources) {
      CCTK_INFO ("Solving puncture equation for BH-NS/NS-NS system");
    } else {
      CCTK_INFO ("Solving puncture equation for BH-BH system");
    }
    CCTK_VInfo (CCTK_THORNSTRING, "b = %g", par_b);
    
    /* initialise to 0 */
    for (int j = 0; j < ntotal; j++)
    {
      cf_v.d0[j] = 0.0;
      cf_v.d1[j] = 0.0;
      cf_v.d2[j] = 0.0;
      cf_v.d3[j] = 0.0;
      cf_v.d11[j] = 0.0;
      cf_v.d12[j] = 0.0;
      cf_v.d13[j] = 0.0;
      cf_v.d22[j] = 0.0;
      cf_v.d23[j] = 0.0;
      cf_v.d33[j] = 0.0;
      v.d0[j] = 0.0;
      v.d1[j] = 0.0;
      v.d2[j] = 0.0;
      v.d3[j] = 0.0;
      v.d11[j] = 0.0;
      v.d12[j] = 0.0;
      v.d13[j] = 0.0;
      v.d22[j] = 0.0;
      v.d23[j] = 0.0;
      v.d33[j] = 0.0;
    }
    /* call for external initial guess */
    /*
    if (use_external_initial_guess)
    {
      set_initial_guess(cctkGH, v);
    }
    */

    /* If bare masses are not given, iteratively solve for them given the 
       target ADM masses target_M_plus and target_M_minus and with initial 
       guesses given by par_m_plus and par_m_minus. */
    if(!(give_bare_mass)) {
      CCTK_REAL tmp, Mp_adm, Mm_adm, Mp_adm_err, Mm_adm_err, up, um;
      char valbuf[100];

      CCTK_REAL M_p = target_M_plus;
      CCTK_REAL M_m = target_M_minus;

      CCTK_VInfo (CCTK_THORNSTRING, "Attempting to find bare masses.");
      CCTK_VInfo (CCTK_THORNSTRING, "Target ADM masses: M_p=%g and M_m=%g",
                  (double) M_p, (double) M_m);
      CCTK_VInfo (CCTK_THORNSTRING, "ADM mass tolerance: %g", (double) adm_tol);

      /* Loop until both ADM masses are within adm_tol of their target */
      do {
        CCTK_VInfo (CCTK_THORNSTRING, "Bare masses: mp=%.15g, mm=%.15g",
                    (double)*mp, (double)*mm);
        Newton (cctkGH, nvar, n1, n2, n3, v, Newton_tol, 1);
	printf("hi %e\n",v.d0[0]);
        F_of_v (cctkGH, nvar, n1, n2, n3, v, F, u);
	//SpecCoef(n1, n2, n3, 0, v.d0, cf_v.d0);

	up = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v, par_b, 0., 0.);
	um = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v,-par_b, 0., 0.);
	
	//up = PunctIntPolAtArbitPositionFast(0, nvar, n1, n2, n3, cf_v, par_b, 0., 0.);       
	//um = PunctIntPolAtArbitPositionFast(0, nvar, n1, n2, n3, cf_v, -par_b, 0., 0.);
	
        /* Calculate the ADM masses from the current bare mass guess */
        Mp_adm = (1 + up) * *mp + *mp * *mm / (4. * par_b);
        Mm_adm = (1 + um) * *mm + *mp * *mm / (4. * par_b);

        /* Check how far the current ADM masses are from the target */
        Mp_adm_err = fabs(M_p-Mp_adm);
        Mm_adm_err = fabs(M_m-Mm_adm);
        CCTK_VInfo (CCTK_THORNSTRING, "ADM mass error: M_p_err=%.15g, M_m_err=%.15g",
                    (double)Mp_adm_err, (double)Mm_adm_err);
        
        /* Invert the ADM mass equation and update the bare mass guess so that
           it gives the correct target ADM masses */
        tmp = -4*par_b*( 1 + um + up + um*up ) + 
                sqrt(16*par_b*M_m*(1 + um)*(1 + up) + 
                  pow(-M_m + M_p + 4*par_b*(1 + um)*(1 + up),2));
        *mp = (tmp + M_p - M_m)/(2.*(1 + up));
        *mm = (tmp - M_p + M_m)/(2.*(1 + um));
        
        /* Set the par_m_plus and par_m_minus parameters */
        sprintf (valbuf,"%.17g", (double) *mp);
        CCTK_ParameterSet ("par_m_plus", "TwoPuncturesAEI", valbuf);
        
        sprintf (valbuf,"%.17g", (double) *mm);
        CCTK_ParameterSet ("par_m_minus", "TwoPuncturesAEI", valbuf);

	printf("hi %e\n",v.d2[0]);

        
      } while ( (Mp_adm_err > adm_tol) ||
                (Mm_adm_err > adm_tol) );
                
      CCTK_VInfo (CCTK_THORNSTRING, "Found bare masses.");
    }

    Newton (cctkGH, nvar, n1, n2, n3, v, Newton_tol, Newton_maxit);
 
    F_of_v (cctkGH, nvar, n1, n2, n3, v, F, u);

    SpecCoef(n1, n2, n3, 0, v.d0, cf_v.d0);

    CCTK_VInfo (CCTK_THORNSTRING,
		  "The two puncture masses are mp=%.17g and mm=%.17g",
                (double) *mp, (double) *mm);

    /* print out ADM mass, eq.: \Delta M_ADM=2*r*u=4*b*V for A=1,B=0,phi=0 */
    //    admMass = (*mp + *mm
    //         - 4*par_b*PunctEvalAtArbitPosition(v.d0, 0, 1, 0, 0, nvar, n1, n2, n3));

    //    aux = PunctEvalAtArbitPosition(v.d0, 0, 1, 0, 0, nvar, n1, n2, n3);
    aux = PunctEvalAtArbitPositionFast(cf_v.d0, 0, 1.0, 0.0, 0.0, nvar, n1, n2, n3);

    admMass = (*mp + *mm
               - 4*par_b*aux);
   
    CCTK_VInfo (CCTK_THORNSTRING, "The total ADM mass is %g", (double) admMass);
    
    char filename[100];
    sprintf(filename,"psi_bbh-proc%d.dat",CCTK_MyProc(cctkGH));
    FILE *outfile = fopen(filename,"w");
    fwrite(mm,sizeof(double),1,outfile);
    fwrite(mp,sizeof(double),1,outfile);
    fclose(outfile);

  } 

  printf("BEFORE FILE READ-IN: %e %e\n",*mm,*mp);
  char filenamein[100];
  sprintf(filenamein,"psi_bbh-proc%d.dat",CCTK_MyProc(cctkGH));
  FILE *infile = fopen(filenamein,"r");
  fread(mm,sizeof(double),1,infile);
  fread(mp,sizeof(double),1,infile);
  fclose(infile);
  printf("AFTER FILE %s READ-IN: %e %e\n",filenamein,*mm,*mp);

  char valbuf2[100];
  /* Set the par_m_plus and par_m_minus parameters */
  sprintf (valbuf2,"%.17g", (double) *mp);
  CCTK_ParameterSet ("par_m_plus", "TwoPuncturesAEI", valbuf2);
  
  sprintf (valbuf2,"%.17g", (double) *mm);
  CCTK_ParameterSet ("par_m_minus", "TwoPuncturesAEI", valbuf2);
  
  if (CCTK_EQUALS(grid_setup_method, "Taylor expansion"))
  {
    gsm = GSM_Taylor_expansion;
  }
  else if (CCTK_EQUALS(grid_setup_method, "evaluation"))
  {
    gsm = GSM_evaluation;
  }
  else
  {
    CCTK_WARN (0, "internal error");
  }

  antisymmetric_lapse = CCTK_EQUALS(initial_lapse, "twopunctures-antisymmetric");
  averaged_lapse = CCTK_EQUALS(initial_lapse, "twopunctures-averaged");
	pmn_lapse = CCTK_EQUALS(initial_lapse, "psi^n");
  if (pmn_lapse)
		CCTK_VInfo(CCTK_THORNSTRING, "Setting initial lapse to psi^%f profile.",
               (double)initial_lapse_psi_exponent);
  brownsville_lapse = CCTK_EQUALS(initial_lapse, "brownsville");
  if (brownsville_lapse)
    CCTK_VInfo(CCTK_THORNSTRING, 
               "Setting initial lapse to a Brownsville-style profile "
               "with exp %f.",
               (double)initial_lapse_psi_exponent);

  int conformal_state_TPAEI = 1;

  CCTK_INFO ("Interpolating result");
  /*
  if (CCTK_EQUALS(metric_type, "static conformal")) {
    if (CCTK_EQUALS(conformal_storage, "factor")) {
      *conformal_state_TPAEI = 1;
    } else if (CCTK_EQUALS(conformal_storage, "factor+derivs")) {
      *conformal_state_TPAEI = 2;
    } else if (CCTK_EQUALS(conformal_storage, "factor+derivs+2nd derivs")) {
      *conformal_state_TPAEI = 3;
    }
  } else {
    *conformal_state_TPAEI = 0;
  }
  */

  for (int d = 0; d < 3; ++ d)
  {
    /*
    imin[d] = 0           + (cctk_bbox[2*d  ] ? 0 : cctk_nghostzones[d]);
    imax[d] = cctk_lsh[d] - (cctk_bbox[2*d+1] ? 0 : cctk_nghostzones[d]);
    */
    imin[d] = 0;
    imax[d] = cctk_lsh[d];
  }


  // Calculate Spectral coefficients once and for all
  // SpecCoef(int n1, int n2, int n3, int ivar, double *v, double *cf)
  SpecCoef(n1, n2, n3, 0, v.d0, cf_v.d0);

#pragma omp parallel for
  for (int k = imin[2]; k < imax[2]; ++k)
  {
    for (int j = imin[1]; j < imax[1]; ++j)
    {
      for (int i = imin[0]; i < imax[0]; ++i)
      {
#if 0
        /* We can't output this when running in parallel */
        if (percent10 != 10*(i+j*cctk_lsh[0]+k*cctk_lsh[0]*cctk_lsh[1]) /
                            (cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2]))
        {
            percent10 = 10*(i+j*cctk_lsh[0]+k*cctk_lsh[0]*cctk_lsh[1]) /
                           (cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2]);
            CCTK_VInfo(CCTK_THORNSTRING, "%3d%% done", percent10*10);
        }
#endif

        const int ind = CCTK_GFINDEX3D (cctkGH, i, j, k);

        CCTK_REAL x1, y1, z1;
        x1 = x[ind] - center_offset[0];
        y1 = y[ind] - center_offset[1];
        z1 = z[ind] - center_offset[2];

        /* We implement swapping the x and z coordinates as follows.
           The bulk of the code that performs the actual calculations
           is unchanged.  This code looks only at local variables.
           Before the bulk --i.e., here-- we swap all x and z tensor
           components, and after the code --i.e., at the end of this
           main loop-- we swap everything back.  */
        if (swap_xz) {
          /* Swap the x and z coordinates */
          SWAP (x1, z1);
        }

        CCTK_REAL r_plus
          = sqrt(pow(x1 - par_b, 2) + pow(y1, 2) + pow(z1, 2));
        CCTK_REAL r_minus
          = sqrt(pow(x1 + par_b, 2) + pow(y1, 2) + pow(z1, 2));

        CCTK_REAL U;
        switch (gsm)
        {
        case GSM_Taylor_expansion:
          U = PunctTaylorExpandAtArbitPosition
            (0, nvar, n1, n2, n3, v, x1, y1, z1);
          break;
        case GSM_evaluation:
	  //          U = PunctIntPolAtArbitPosition
	  //            (0, nvar, n1, n2, n3, v, x1, y1, z1);
          U = PunctIntPolAtArbitPositionFast
            (0, nvar, n1, n2, n3, cf_v, x1, y1, z1);
          break;
        default:
          assert (0);
        }
        r_plus = pow (pow (r_plus, 4) + pow (TP_epsilon, 4), 0.25);
        r_minus = pow (pow (r_minus, 4) + pow (TP_epsilon, 4), 0.25);
        if (r_plus < TP_Tiny)
            r_plus = TP_Tiny;
        if (r_minus < TP_Tiny)
            r_minus = TP_Tiny;
        CCTK_REAL psi1 = 1
          + 0.5 * *mp / r_plus
          + 0.5 * *mm / r_minus + U;
#define EXTEND(M,r) \
          ( M * (3./8 * pow(r, 4) / pow(TP_Extend_Radius, 5) - \
                 5./4 * pow(r, 2) / pow(TP_Extend_Radius, 3) + \
                 15./8 / TP_Extend_Radius))
        if (r_plus < TP_Extend_Radius) {
          psi1 = 1
             + 0.5 * EXTEND(*mp,r_plus)
             + 0.5 * *mm / r_minus + U;
        }
        if (r_minus < TP_Extend_Radius) {
          psi1 = 1
             + 0.5 * EXTEND(*mm,r_minus)
             + 0.5 * *mp / r_plus + U;
        }
        CCTK_REAL static_psi = 1;
        
        CCTK_REAL Aij[3][3];
        BY_Aijofxyz (x1, y1, z1, Aij);

        CCTK_REAL old_alp;
        if (multiply_old_lapse)
            old_alp = lapm1[ind];

        if ((conformal_state_TPAEI > 0) || (pmn_lapse) || (brownsville_lapse)) {

          CCTK_REAL xp, yp, zp, rp, ir;
          CCTK_REAL s1, s3, s5;
          CCTK_REAL p, px, py, pz, pxx, pxy, pxz, pyy, pyz, pzz;
          p = 1.0;
          px = py = pz = 0.0;
          pxx = pxy = pxz = 0.0;
          pyy = pyz = pzz = 0.0;

          /* first puncture */
          xp = x1 - par_b;
          yp = y1;
          zp = z1;
          rp = sqrt (xp*xp + yp*yp + zp*zp);
          rp = pow (pow (rp, 4) + pow (TP_epsilon, 4), 0.25);
          if (rp < TP_Tiny)
              rp = TP_Tiny;
          ir = 1.0/rp;

          if (rp < TP_Extend_Radius) {
            ir = EXTEND(1., rp);
          }

          s1 = 0.5* *mp *ir;
          s3 = -s1*ir*ir;
          s5 = -3.0*s3*ir*ir;

          p += s1;

          px += xp*s3;
          py += yp*s3;
          pz += zp*s3;

          pxx += xp*xp*s5 + s3;
          pxy += xp*yp*s5;
          pxz += xp*zp*s5;
          pyy += yp*yp*s5 + s3;
          pyz += yp*zp*s5;
          pzz += zp*zp*s5 + s3;

          /* second puncture */
          xp = x1 + par_b;
          yp = y1;
          zp = z1;
          rp = sqrt (xp*xp + yp*yp + zp*zp);
          rp = pow (pow (rp, 4) + pow (TP_epsilon, 4), 0.25);
          if (rp < TP_Tiny)
              rp = TP_Tiny;
          ir = 1.0/rp;

          if (rp < TP_Extend_Radius) {
            ir = EXTEND(1., rp);
          }

          s1 = 0.5* *mm *ir;
          s3 = -s1*ir*ir;
          s5 = -3.0*s3*ir*ir;

          p += s1;

          px += xp*s3;
          py += yp*s3;
          pz += zp*s3;

          pxx += xp*xp*s5 + s3;
          pxy += xp*yp*s5;
          pxz += xp*zp*s5;
          pyy += yp*yp*s5 + s3;
          pyz += yp*zp*s5;
          pzz += zp*zp*s5 + s3;

          if (conformal_state_TPAEI >= 1) {
            static_psi = p;
            psi[ind] = static_psi + U;
          }

          if (pmn_lapse)
            lapm1[ind] = pow(psi[ind], initial_lapse_psi_exponent) - 1.0;
	  //lapm1[ind] = pow(p, initial_lapse_psi_exponent) - 1.0;
          if (brownsville_lapse)
            lapm1[ind] = 2.0/(1.0+pow(p, initial_lapse_psi_exponent)) - 1.0;

        } /* if conformal-state > 0 */
          
        puncture_u[ind] = U;

        gxx[ind] = 1.0; //pow (psi1 / static_psi, 4);
        gxy[ind] = 0;
        gxz[ind] = 0;
        gyy[ind] = 1.0; //pow (psi1 / static_psi, 4);
        gyz[ind] = 0;
        gzz[ind] = 1.0; //pow (psi1 / static_psi, 4);

        kxx[ind] = Aij[0][0] / pow(psi1, 2);
        kxy[ind] = Aij[0][1] / pow(psi1, 2);
        kxz[ind] = Aij[0][2] / pow(psi1, 2);
        kyy[ind] = Aij[1][1] / pow(psi1, 2);
        kyz[ind] = Aij[1][2] / pow(psi1, 2);
        kzz[ind] = Aij[2][2] / pow(psi1, 2);

        if (antisymmetric_lapse || averaged_lapse) {
          lapm1[ind] =
            ((1.0 -0.5* *mp /r_plus -0.5* *mm/r_minus)
            /(1.0 +0.5* *mp /r_plus +0.5* *mm/r_minus)) - 1.0;

          if (r_plus < TP_Extend_Radius) {
            lapm1[ind] =
              ((1.0 -0.5*EXTEND(*mp, r_plus) -0.5* *mm/r_minus)
              /(1.0 +0.5*EXTEND(*mp, r_plus) +0.5* *mm/r_minus)) - 1.0;
          }
          if (r_minus < TP_Extend_Radius) {
            lapm1[ind] =
              ((1.0 -0.5*EXTEND(*mm, r_minus) -0.5* *mp/r_plus)
              /(1.0 +0.5*EXTEND(*mp, r_minus) +0.5* *mp/r_plus)) - 1.0;
          }
          
          if (averaged_lapse) {
            lapm1[ind] = 0.5 * (1.0 + (lapm1[ind]+1.0)) - 1.0;
          }
        }
        if (multiply_old_lapse)
          lapm1[ind] = (lapm1[ind]+1.0)*(old_alp+1.0) - 1.0;

        if (swap_xz) {
          /* Swap the x and z components of all tensors */
          SWAP (gxx[ind], gzz[ind]);
          SWAP (gxy[ind], gyz[ind]);
          SWAP (kxx[ind], kzz[ind]);
          SWAP (kxy[ind], kyz[ind]);
        } /* if swap_xz */

      } /* for i */
    }   /* for j */
  }     /* for k */

  /*
  if (use_sources && rescale_sources)
  {
    Rescale_Sources(cctkGH,
                    cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2],
                    x, y, z,
                    (conformal_state_TPAEI > 0) ? psi : NULL,
                    gxx, gyy, gzz,
                    gxy, gxz, gyz);
  }
  */

  if (0) {
    /* Keep the result around for the next time */
    free_dvector (F, 0, ntotal - 1);
    free_derivs (&u, ntotal);
    free_derivs (&v, ntotal);
  }
}
