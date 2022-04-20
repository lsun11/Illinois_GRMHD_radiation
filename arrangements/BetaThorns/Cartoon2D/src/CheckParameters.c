/*@@
   @file      CheckParameters.c
   @date      Wed Nov  3 10:17:46 MET 1999
   @author    Sai Iyer
   @desc 
              Check Cartoon2D parameters
   @enddesc 
 @@*/

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

static const char *rcsid = "$Id: CheckParameters.c,v 1.8 2004/05/14 10:52:19 schnetter Exp $";

CCTK_FILEVERSION(Development_Cartoon2D_CheckParameters_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void Cartoon2D_CheckParameters(CCTK_ARGUMENTS);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    Cartoon2D_CheckParameters
   @date       Wed Nov  3 10:17:46 MET 1999
   @author     Sai Iyer
   @desc 
               Check Cartoon2D parameters
   @enddesc 
   @calls      
   @calledby   
   @history 

   @endhistory 

@@*/

void Cartoon2D_CheckParameters(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if(order > 2*cctk_nghostzones[0])
  {
    CCTK_PARAMWARN("Ghostzone width in x direction too small.");
  }

  return;
}

