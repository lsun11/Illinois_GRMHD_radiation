 /*@@
   @file      Startup.c
   @date      Sat 10th November 2001
   @author    Gabrielle Allen
   @desc
              Startup routines for AHFinder.
   @enddesc
   @version   $Id: Startup.c,v 1.2 2004/05/17 12:28:56 goodale Exp $
 @@*/

#include "cctk.h"

static const char *rcsid = "$Header: /cactus/CactusEinstein/AHFinder/src/Startup.c,v 1.2 2004/05/17 12:28:56 goodale Exp $";

CCTK_FILEVERSION(CactusEinstein_AHFinder_Startup_c)

int AHFinder_TimeForOutput (const cGH *GH, int vindex);

/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
void AHFinder_Startup (void);

 /*@@
   @routine   AHfinder_Startup
   @date      Sat 10th November 2001
   @author    Gabrielle Allen
   @desc
              The startup registration routine for AHFinder.
              Registers AHFinder as an IO Method and provide the
              only necessary method TimeForOutput
   @enddesc
   @calls     CCTK_RegisterGHExtensionSetupGH
@@*/
void AHFinder_Startup (void)
{
  int handle;

  handle = CCTK_RegisterIOMethod ("AHFinder");
  CCTK_RegisterIOMethodTimeToOutput (handle, AHFinder_TimeForOutput);
}




