 /*@@
   @file      Simple.c
   @date      September 4 1999
   @author    Gabrielle Allen
   @desc 
   Standard specification of timestep
   @enddesc 
 @@*/

#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header: /cactusdevcvs/CactusBase/Time/src/Simple.c,v 1.17 2004/10/05 20:01:05 yye00 Exp $";

CCTK_FILEVERSION(CactusBase_Time_Simple_c)

void Time_Simple(CCTK_ARGUMENTS);

void Time_Simple(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL min_spacing;
  int d;

  /* Calculate the minimum grid spacing */
  min_spacing = cctk_delta_space[0];
  for (d=1; d<cctk_dim; ++d)
  {
    min_spacing = (min_spacing<cctk_delta_space[d]
                   ? min_spacing : cctk_delta_space[d]);
  }

  /* Calculate the timestep */
  cctkGH->cctk_delta_time = dtfac*min_spacing;

  if (verbose)
  {
    CCTK_VInfo(CCTK_THORNSTRING,
               "Using a simple Courant condition to set then timestep");
    CCTK_VInfo(CCTK_THORNSTRING,
               "  ... using a dtfac of %g", (double)dtfac);
    CCTK_VInfo(CCTK_THORNSTRING,
               "  ... using a minimum spacing of %g", (double)min_spacing);
  }

  CCTK_VInfo(CCTK_THORNSTRING,
             "Timestep set to %g (courant_static)",
             (double)(cctkGH->cctk_delta_time/cctkGH->cctk_timefac));
}
