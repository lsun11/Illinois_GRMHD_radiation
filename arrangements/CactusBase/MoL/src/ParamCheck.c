 /*@@
   @file      ParamCheck.c
   @date      Mon May 20 09:50:55 2002
   @author    Ian Hawke
   @desc 
   Basic parameter checking for thorn MoL.
   @enddesc 
   @version   $Header: /cactusdevcvs/CactusBase/MoL/src/ParamCheck.c,v 1.14 2006/12/05 17:44:41 schnetter Exp $
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

static const char *rcsid = "$Header: /cactusdevcvs/CactusBase/MoL/src/ParamCheck.c,v 1.14 2006/12/05 17:44:41 schnetter Exp $";

CCTK_FILEVERSION(CactusBase_MoL_ParamCheck_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

int MoL_ParamCheck(CCTK_ARGUMENTS);

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
   @routine    MoL_ParamCheck
   @date       Mon May 20 09:56:05 2002
   @author     Ian Hawke
   @desc 
   Basic parameter checking.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

int MoL_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT options_table, ierr, GenericIntermediateSteps;

  if (CCTK_Equals(ODE_Method, "Generic"))
  {
    if (MoL_Num_Scratch_Levels < MoL_Intermediate_Steps - 1)
    {
      CCTK_PARAMWARN("When using a generic solver the number "
                     "of scratch levels must be at least the "
                     "number of intermediate steps - 1");
    }
    if ( (CCTK_Equals(Generic_Type, "Classic RK3"))&&
         ((!(MoL_Intermediate_Steps == 3))||(!(MoL_Num_Scratch_Levels > 1))) )
    {
      CCTK_PARAMWARN("When using the classic RK3 evolver the "
                     "number of intermediate steps must be 3 "
                     "and the number of scratch levels at least 2");
    }
    if (CCTK_Equals(Generic_Type, "Table"))
    {
      options_table =
        Util_TableCreateFromString(Generic_Method_Descriptor);
      if (options_table < 0)
      {
        CCTK_WARN(0, "Failed to create table from "
                  "Generic_Method_Descriptor!");
      }
      ierr = Util_TableGetInt(options_table,
                              &GenericIntermediateSteps,
                              "GenericIntermediateSteps");
      if (ierr < 1)
      {
        if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY)
        {
          CCTK_WARN(0, "When using the generic table options "
                    "you must set \"GenericIntermediateSteps\" in "
                    "the options table");
        }
        else
        {
          CCTK_WARN(0, "Table error - check with Ian.");
        }
      }
      if (MoL_Intermediate_Steps != GenericIntermediateSteps)
      {
        CCTK_PARAMWARN("The number of intermediate steps must "
                       "equal the number specified in the table options!");
      }
      ierr = Util_TableDestroy(options_table);
    }
  }

  if ( (CCTK_Equals(ODE_Method, "RK2"))&&(!(MoL_Intermediate_Steps == 2)) )
  {
    CCTK_PARAMWARN("When using the efficient RK2 evolver the "
                   "number of intermediate steps must be 2");
  }

  if ( (CCTK_Equals(ODE_Method, "RK3"))&&(!(MoL_Intermediate_Steps == 3)) )
  {
    CCTK_PARAMWARN("When using the efficient RK3 evolver the "
                   "number of intermediate steps must be 3");
  }

  if ( (CCTK_Equals(ODE_Method, "RK4")) && ( (!(MoL_Intermediate_Steps == 4))
       || (!(MoL_Num_Scratch_Levels > 0)) ) )
  {
    CCTK_PARAMWARN("When using the efficient RK4 evolver the "
                   "number of intermediate steps must be 4, and"
                   " the number of scratch levels at least 1");
  }

  if ( (CCTK_Equals(ODE_Method, "RK45") || CCTK_Equals(ODE_Method, "RK45CK")) &&
       ( !((MoL_Intermediate_Steps == 6)||(MoL_Num_Scratch_Levels > 5)) ) )
  {
    CCTK_PARAMWARN("When using the RK45 or RK45CK evolver, the "
                   "number of intermediate steps must be 6 " 
                   "and the number of scratch levels at least 6.");
  }

  if ( (CCTK_Equals(ODE_Method, "RK65")) &&
       ( !((MoL_Intermediate_Steps == 8)||(MoL_Num_Scratch_Levels > 7)) ) )
  {
    CCTK_PARAMWARN("When using the RK65 evolver the "
                   "number of intermediate steps must be 8 " 
                   "and the number of scratch levels at least 8.");
  }

  if ( (CCTK_Equals(ODE_Method, "RK87")) &&
       ( !((MoL_Intermediate_Steps == 13)||(MoL_Num_Scratch_Levels > 12)) ) )
  {
    CCTK_PARAMWARN("When using the RK87 evolver the "
                   "number of intermediate steps must be 13 " 
                   "and the number of scratch levels at least 13.");
  }

  if (adaptive_stepsize)
  {
    if (CCTK_Equals(ODE_Method, "RK45")||CCTK_Equals(ODE_Method, "RK45CK")||CCTK_Equals(ODE_Method, "RK65")||CCTK_Equals(ODE_Method, "RK87"))
    {
      /* everything is fine, do nothing */
    }
    else
    {
      CCTK_PARAMWARN("Adaptive time step sizes are only possible with the RK45, RK45CK, RK65, and RK87 solvers");
    }
  }

  return 0;
  
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
