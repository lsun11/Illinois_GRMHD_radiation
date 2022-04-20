 /*@@
   @file      StepSize.c
   @date      Tue Sep 07 2004
   @author    Erik Schnetter
   @desc 
   Control the time step size.
   @enddesc 
   @version   $Header: /cactusdevcvs/CactusBase/MoL/src/StepSize.c,v 1.3 2006/01/23 10:39:58 hawke Exp $
 @@*/

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "ExternalVariables.h"

static const char *rcsid = "$Header: /cactusdevcvs/CactusBase/MoL/src/StepSize.c,v 1.3 2006/01/23 10:39:58 hawke Exp $";

CCTK_FILEVERSION(CactusBase_MoL_StepSize_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MoL_StartLoop(CCTK_ARGUMENTS);

void MoL_InitAdaptiveError(CCTK_ARGUMENTS);
void MoL_FindAdaptiveError(CCTK_ARGUMENTS);
void MoL_ReduceAdaptiveError(CCTK_ARGUMENTS);

void MoL_SetEstimatedDt(CCTK_ARGUMENTS);

void MoL_FinishLoop(CCTK_ARGUMENTS);

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
   @routine    MoL_StartLoop
   @date       Tue Sep 07 2004
   @author     Erik Schnetter
   @desc 
   Start the step size control loop, so that at least one iteration is done.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void
MoL_StartLoop(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  printf("Start MoL_StartLoop!!!!! \n");

  *MoL_Stepsize_Bad = 1;

  if (adaptive_stepsize)
  {
    *EstimatedDt = cctkGH->cctk_delta_time;
  }
  
}

 /*@@
   @routine    MoL_InitAdaptiveError
   @date       Tue Sep 07 2004
   @author     Erik Schnetter
   @desc 
   Initialize error counters for adaptive stepsize control
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

static inline CCTK_REAL
square (CCTK_REAL const x)
{
  return x * x;
}

void MoL_InitAdaptiveError(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* Initialise global error */
  *Error = 0;
  *Count = 0;
}

 /*@@
   @routine    MoL_FindAdaptiveError
   @date       Thu Jan 27 10:22:26 2005
   @author     Erik Schnetter
   @desc 
   Compute the error local to this component/patch/...
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_FindAdaptiveError(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL const * restrict UpdateVar;
  CCTK_REAL const * restrict RHSVar;
  CCTK_REAL const * restrict ErrorVar;

  int imin[3], imax[3];

  int var;
  int i, j, k;
  int d;

  assert (cctk_dim <= 3);
  for (d = 0; d < cctk_dim; d++)
  {
    imin[d] = cctk_bbox[2*d] ? 0 : cctk_nghostzones[d];
    imax[d] = cctk_lsh[d] - (cctk_bbox[2*d+1] ? 0 : cctk_nghostzones[d]);
  }

  /* Calculate absolute error */
  for (var = 0; var < MoLNumEvolvedVariables; var++)
  {

    UpdateVar = CCTK_VarDataPtrI(cctkGH, 0, EvolvedVariableIndex[var]);
    RHSVar = CCTK_VarDataPtrI(cctkGH, 0, RHSVariableIndex[var]);
    ErrorVar
      = CCTK_VarDataPtrI(cctkGH, 0, 
                         CCTK_FirstVarIndex("MOL::ERRORESTIMATE") + var);

    assert (cctk_dim == 3);
    for (k = imin[2]; k < imax[2]; k++)
    {
      for (j = imin[1]; j < imax[1]; j++)
      {
        for (i = imin[0]; i < imax[0]; i++)
        {
          int const index = CCTK_GFINDEX3D(cctkGH, i, j, k);
          CCTK_REAL const scale
            = (square(maximum_absolute_error)
               + square(maximum_relative_error * UpdateVar[index])
               + square(maximum_relative_error * RHS_error_weight
                        * (*Original_Delta_Time) * RHSVar[index]));
          *Error += square(ErrorVar[index]) / scale;
        }
      }
    }

  } /* for var */

  *Count
    += (MoLNumEvolvedVariables
        * (imax[0] - imin[0]) * (imax[1] - imin[1]) * (imax[2] - imin[2]));
}

 /*@@
   @routine    MoL_ReduceAdaptiveError
   @date       Thu Jan 27 10:23:14 2005
   @author     Erik Schnetter
   @desc 
   Find the global error estimate. 
   Change the timestep based on the error estimate.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_ReduceAdaptiveError(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int redop;

  CCTK_REAL red_local[2], red_global[2];
  CCTK_REAL p1, p2;
  int ierr;

  /* Get global result over all processors */
  redop = CCTK_ReductionHandle ("sum");
  assert (redop >= 0);

  red_local[0] = *Error;
  red_local[1] = *Count;
  ierr = CCTK_ReduceLocArrayToArray1D
    (cctkGH, -1, redop, red_local, red_global, 2, CCTK_VARIABLE_REAL);
  assert (ierr == 0);
  *Error = red_global[0];
  *Count = red_global[1];

  /* Calculate L2-norm */
  *Error = sqrt(*Error / *Count);
  if (! CCTK_EQUALS(verbose, "none"))
  {
    CCTK_VInfo (CCTK_THORNSTRING, "Integration accuracy quotient is %g", (double)*Error);
  }

  if ( CCTK_EQUALS(ODE_Method,"RK45") )
  {
     p1 = 0.2;
     p2 = 0.25;
  }

  if ( CCTK_EQUALS(ODE_Method,"RK65") )
  {
     p1 = 1.0/6.0;
     p2 = 0.2;
  }

  if ( CCTK_EQUALS(ODE_Method,"RK87") )
  {
     p1 = 0.125;
     p2 = 1.0/7.0;
  }

  /* Decide whether to accept this step */
  *MoL_Stepsize_Bad = *Error > 1;

  if (*MoL_Stepsize_Bad)
  {
    /* The error is too large; reject the time step and reduce the
       step size */
    cctkGH->cctk_time -= cctkGH->cctk_delta_time;
    if (! CCTK_EQUALS(verbose, "none"))
    {
      CCTK_VInfo (CCTK_THORNSTRING, "*** REJECTING TIME STEP ***");
    }

    cctkGH->cctk_delta_time
      = safety_factor * (*Original_Delta_Time) / pow(*Error, p1);
    if (! CCTK_EQUALS(verbose, "none"))
    {
      CCTK_VInfo (CCTK_THORNSTRING, "Setting time step to %g", (double)cctkGH->cctk_delta_time);
    }

    if (cctkGH->cctk_delta_time < (*Original_Delta_Time) / maximum_decrease)
    {
      /* No more than a factor of 10 decrease */
      cctkGH->cctk_delta_time = (*Original_Delta_Time) / maximum_decrease;
      if (! CCTK_EQUALS(verbose, "none"))
      {
        CCTK_VInfo (CCTK_THORNSTRING, "   Time step reduction too large; clamping time step to %g", (double)cctkGH->cctk_delta_time);
      }
    }

    cctkGH->cctk_time += cctkGH->cctk_delta_time;
  }
  else
  {
    /* The error is acceptable; estimate the next step size */
    *EstimatedDt = 
      safety_factor * (*Original_Delta_Time) / pow(*Error, p2);

    if (*EstimatedDt > (*Original_Delta_Time) * maximum_increase)
    {
      /* No more than a factor of 5 increase */
      *EstimatedDt = (*Original_Delta_Time) * maximum_increase;
      if (! CCTK_EQUALS(verbose, "none"))
      {
        CCTK_VInfo (CCTK_THORNSTRING, "   Time step increase too large; clamping time step to %g", (double)(*EstimatedDt));
      }
    }
  }
}

 /*@@
   @routine    MoL_SetEstimatedDt
   @date       Thu Jan 27 14:05:08 2005
   @author     Ian Hawke
   @desc 
   Actually set the timestep in PostStep. 
   This avoids problems when the timestep is changed in the middle of the
   evolution loop.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_SetEstimatedDt(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  cctkGH->cctk_delta_time = *EstimatedDt;

  if (! CCTK_EQUALS(verbose, "none"))
  {
    CCTK_VInfo (CCTK_THORNSTRING, "Setting time step to %g", 
                (double)cctkGH->cctk_delta_time);
  }

}

 /*@@
   @routine    MoL_FinishLoop
   @date       Thu Jan 27 10:24:19 2005
   @author     Erik Schnetter
   @desc 
   Loop control if adaptive timestepping is not used.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_FinishLoop(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* Keep time step size unchanged */
  *MoL_Stepsize_Bad = 0;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
