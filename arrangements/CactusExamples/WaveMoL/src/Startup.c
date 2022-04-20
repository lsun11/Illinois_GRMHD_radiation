 /*@@
   @file      Startup.c
   @date      Fri Nov  9 13:27:48 2001
   @author    Ian Hawke
   @desc 
   Register banner - a straight copy of WaveToy.
   @enddesc 
 @@*/


#include "cctk.h"

static const char *rcsid = "$Header: /cactusdevcvs/CactusExamples/WaveMoL/src/Startup.c,v 1.1.1.1 2004/06/26 14:49:02 hawke Exp $";

CCTK_FILEVERSION(CactusExamples_WaveMoL_Startup_c)

int WaveMoL_Startup(void);

 /*@@
   @routine    WaveMoL_Startup
   @date       Fri Nov  9 13:28:58 2001
   @author     Ian Hawke
   @desc 
   Again a straight copy of WaveToyC
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

int WaveMoL_Startup(void)
{

   const char *banner = "WaveMoL: Evolutions of a Scalar Field";

   CCTK_RegisterBanner(banner);

   return 0;
}
