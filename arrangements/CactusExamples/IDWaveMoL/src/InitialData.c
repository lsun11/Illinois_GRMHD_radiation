 /*@@
   @file      InitialData.c
   @date      
   @author    Werner Benger   
   @desc 
              Initial data for the 3D Wave Equation
	      Derived from Tom Goodale
   @enddesc 
   @version $Header: /cactusdevcvs/CactusExamples/IDWaveMoL/src/InitialData.c,v 1.3 2005/08/23 19:53:12 schnetter Exp $
 @@*/

#include <math.h>

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

/*
#ifndef DEBUG_MOL
#define DEBUG_MOL
#endif
*/

static const char *rcsid = "$Header: /cactusdevcvs/CactusExamples/IDWaveMoL/src/InitialData.c,v 1.3 2005/08/23 19:53:12 schnetter Exp $";

CCTK_FILEVERSION(CactusExamples_IDWaveMoL_InitialData_c)

static inline CCTK_REAL sqr(const CCTK_REAL val)
{
  return val*val;
}

void IDWaveMoL_InitialData(CCTK_ARGUMENTS);


 /*@@
   @routine    IDWaveMoL_InitialData
   @date       
   @author     Tom Goodale
   @desc 
               Set up initial data for the wave equation
   @enddesc 
   @calls      
   @calledby   
   @history 
   @hdate Mon Oct 11 11:48:03 1999 @hauthor Werner Benger
   @hdesc  Converted to C++
   @hdate Mon Oct 11 11:48:20 1999 @hauthor Tom Goodale
   @hdesc Added the rest of the initial data.
   @hdate Thu Feb 17 09:22:01 2000 @hauthor Tom Goodale
   @hdesc Converted to C
   @endhistory 

@@*/

void IDWaveMoL_InitialData(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  int i,j,k;

  CCTK_REAL omega;
  CCTK_REAL arg;
  int index;
  CCTK_REAL T, X, Y, Z, R;
  const CCTK_REAL pi = acos(-1.0);

#ifdef DEBUG_MOL
  printf("Initial data says sigma = %g, radius = %g\n", sigma, radius);
#endif

  if (CCTK_Equals(initial_data, "gaussian"))
  {
    for(k=0; k<cctk_lsh[2]; k++)
    {
      for(j=0; j<cctk_lsh[1]; j++)
      {
        for(i=0; i<cctk_lsh[0]; i++)
        {
          index = CCTK_GFINDEX3D(cctkGH,i,j,k);

          X = x[index] - centrex;
          Y = y[index] - centrey;
          Z = z[index] - centrez;

          R = sqrt(X*X + Y*Y + Z*Z);

          phi[index] = amplitude*exp( - sqr( (R - radius) / sigma ) );
          phit[index] = 0.0;
          if (R > 1.e-10) {
            phix[index] = -2.0*(R-radius)/sqr(sigma)*X/R*phi[index];
            phiy[index] = -2.0*(R-radius)/sqr(sigma)*Y/R*phi[index];
            phiz[index] = -2.0*(R-radius)/sqr(sigma)*Z/R*phi[index];
          } 
          else {
            phix[index] = 0.0;
            phiy[index] = 0.0;
            phiz[index] = 0.0;
          } 
        }
      }
    }
  }
  else if (CCTK_Equals(initial_data, "plane"))
  {
    omega = sqrt(sqr(kx) + sqr(ky) + sqr(kz));
    T = cctk_time - offsett;

    for(k=0; k<cctk_lsh[2]; k++)
    {
      for(j=0; j<cctk_lsh[1]; j++)
      {
        for(i=0; i<cctk_lsh[0]; i++)
        {
          index = CCTK_GFINDEX3D(cctkGH,i,j,k);

          X = x[index] - offsetx;
          Y = y[index] - offsety;
          Z = z[index] - offsetz;

          arg = 2*pi * (kx*X + ky*Y + kz*Z - omega*T);

          phi[index]  = slopet*T + slopex*X + slopey*Y + slopez*Z
                               +        amplitude *         cos (arg);
          phit[index] = slopet + 2*pi * amplitude * omega * sin (arg);
          phix[index] = slopex - 2*pi * amplitude * kx    * sin (arg);
          phiy[index] = slopey - 2*pi * amplitude * ky    * sin (arg);
          phiz[index] = slopez - 2*pi * amplitude * kz    * sin (arg);
        }
      }
    }
  }
  else 
  {
    CCTK_WARN(0, "Value of parameter \"initial_data\" not recognized.");
  }

}

