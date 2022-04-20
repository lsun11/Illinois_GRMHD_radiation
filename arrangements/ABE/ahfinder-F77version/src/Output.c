 /*@@
   @file      Output.c
   @date      November 10 2001
   @author    Gabrielle Allen
   @desc
              Functions to deal with making AHFinder an IO method
   @enddesc
   @version   $Id: Output.c,v 1.2 2002/01/04 11:39:32 tradke Exp $
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header: /cactus/CactusEinstein/AHFinder/src/Output.c,v 1.2 2002/01/04 11:39:32 tradke Exp $";

CCTK_FILEVERSION(CactusEinstein_AHFinder_Output_c)

static int ahf_ncall=0;

int AHFinder_TimeForOutput (const cGH *GH, int vindex);

 /*@@
   @routine    AHFinder_TimeForOutput
   @date       June 31 1999
   @author     Gabrielle Allen
   @desc
               Decides if it is time to output a variable using info output
   @enddesc
   @calls      CheckSteerableParameters

   @var        GH
   @vdesc      Pointer to CCTK GH
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        vindex
   @vdesc      index of variable to check for output
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               true/false (1 or 0) if analysis should be called
   @endreturndesc
@@*/
int AHFinder_TimeForOutput (const cGH *GH, int vindex)
{
  DECLARE_CCTK_PARAMETERS

  int retval=1;
  int ahfafter;

  if (vindex==CCTK_VarIndex("ahfinder::triggervar"))
  {

    if (ahf_findaftertime <= 0) 
    {
      ahfafter = ahf_findafter;
      if (GH->cctk_iteration < ahfafter) {retval=0;}
    }
    else
    {
      if (GH->cctk_time < ahf_findaftertime) {retval=0;}
      if (ahf_ncall == 0) {ahfafter = GH->cctk_iteration;}
    }
  
    if ( ((GH->cctk_iteration-ahfafter) % ahf_findevery) ) 
    {
      retval = 0;
    }

    if (retval==1) {ahf_ncall = 0;}

  }
  else
  {
    retval = 0;
  }
      
  return retval;
}





