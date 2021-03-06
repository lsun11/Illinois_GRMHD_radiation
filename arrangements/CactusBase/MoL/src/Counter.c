 /*@@
   @file      Counter.c
   @date      Mon May 20 09:52:33 2002
   @author    Ian Hawke
   @desc 
   Routines setting and altering the loop counter.
   @enddesc 
   @version   $Header: /cactusdevcvs/CactusBase/MoL/src/Counter.c,v 1.11 2004/07/01 10:57:18 hawke Exp $
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

/* #include "carpet.h" */

static const char *rcsid = "$Header: /cactusdevcvs/CactusBase/MoL/src/Counter.c,v 1.11 2004/07/01 10:57:18 hawke Exp $";

CCTK_FILEVERSION(CactusBase_MoL_Counter_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

int MoL_SetCounter(CCTK_ARGUMENTS);

int MoL_DecrementCounter(CCTK_ARGUMENTS);

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
   @routine    MoL_SetCounter
   @date       Mon May 20 09:54:39 2002
   @author     Ian Hawke
   @desc 
   Initially set the counter to the number of intermediate steps.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

int MoL_SetCounter(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  printf("Start MoL_SetCounter!!!!! \n");


  *MoL_Intermediate_Step = MoL_Intermediate_Steps;

   
/* #ifdef HAVE_CARPET */
  if ((*MoL_Intermediate_Step))
  {
    /* Disable prolongating during the iterations */
/*     CarpetEnableProlongating (0); */
    if (disable_prolongation)
    {
      if (CCTK_IsFunctionAliased("EnableProlongating"))
      {
        EnableProlongating(0);
      }
      else
      {
        CCTK_WARN(4, "Cannot disable prolongation as function"
                  " \"EnableProlongating\" is not provided by any thorn!");
      }
    }
  }
/* #endif */



  return 0;
}

 /*@@
   @routine    MoL_DecrementCounter
   @date       Mon May 20 09:55:13 2002
   @author     Ian Hawke
   @desc 
   During the loop decrement the counter
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

int MoL_DecrementCounter(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  (*MoL_Intermediate_Step) --;

/* #ifdef HAVE_CARPET */
  if (! (*MoL_Intermediate_Step))
  {
    /* Re-enable prolongating before the final PostStep */
/*     CarpetEnableProlongating (1); */
    if (disable_prolongation)
    {
      if (CCTK_IsFunctionAliased("EnableProlongating"))
      {
        EnableProlongating(1);
      }
      else
      {
        CCTK_WARN(4, "Cannot enable prolongation as function"
                  " \"EnableProlongating\" is not provided by any thorn!");
      }
    }
  }
/* #endif */

  return 0;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
