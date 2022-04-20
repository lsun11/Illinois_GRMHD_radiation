 /*@@
   @file      Initialise.c
   @date      May 12 2001
   @author    Gabrielle Allen
   @desc 
   Initialise grid variables
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"

static const char *rcsid = "$Header: /cactusdevcvs/CactusBase/Time/src/Initialise.c,v 1.3 2004/01/15 15:50:25 schnetter Exp $";

CCTK_FILEVERSION(CactusBase_Time_Initialise_c)

void Time_Initialise(CCTK_ARGUMENTS);

void Time_Initialise(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS

  *courant_wave_speed = 0;  
  *courant_min_time = 0;
  *courant_dt = 0;
}
