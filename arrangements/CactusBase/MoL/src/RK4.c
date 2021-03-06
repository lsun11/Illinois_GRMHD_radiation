 /*@@
   @file      RK4.c
   @date      Fri July 14, 2006
   @author    Yosef Zlochower
   @desc 
   A routine to perform RK4 evolution. Mostly copied from
   genericRK.c
   @enddesc 
   @version   $Header: /cactusdevcvs/CactusBase/MoL/src/RK4.c,v 1.3 2008/09/23 15:46:30 schnetter Exp $
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <stdio.h>
#include "ExternalVariables.h"

static const char *rcsid = "$Header: /cactusdevcvs/CactusBase/MoL/src/RK4.c,v 1.3 2008/09/23 15:46:30 schnetter Exp $";

CCTK_FILEVERSION(CactusBase_MoL_RK4_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MoL_RK4Add(CCTK_ARGUMENTS);

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
   @routine    MoL_RK4Add
   @date       
   @author     
   @desc 
   Performs a single step of a RK4 type time
   integration.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_RK4Add(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
    
  cGroupDynamicData arraydata;
  CCTK_INT groupindex, ierr;
  CCTK_INT arraytotalsize, arraydim;

  /* FIXME */

#ifdef MOLDOESCOMPLEX

CCTK_WARN(0, "not implemented");
  CCTK_COMPLEX Complex_alpha, Complex_beta, Complex_Delta_Time;
  CCTK_COMPLEX * restrict OldComplexVar;
  CCTK_COMPLEX * restrict UpdateComplexVar;
  CCTK_COMPLEX const * restrict RHSComplexVar;
  CCTK_COMPLEX * restrict ScratchComplexVar;

#endif

  static CCTK_INT scratchspace_firstindex = -99;
  CCTK_INT index, var, scratchstep, alphaindex, scratchindex;
  CCTK_INT totalsize, singlearraysize;
  CCTK_REAL alpha, beta;
  CCTK_REAL * restrict UpdateVar;
  CCTK_REAL * restrict OldVar;
  CCTK_REAL const * restrict RHSVar;
  CCTK_REAL * restrict ScratchVar;

  CCTK_INT arrayscratchlocation;

  totalsize = 1;
  for (arraydim = 0; arraydim < cctk_dim; arraydim++)
  {
    totalsize *= cctk_lsh[arraydim];
  }

  if (scratchspace_firstindex == -99)
  {
    scratchspace_firstindex = CCTK_FirstVarIndex("MOL::SCRATCHSPACE");
  }

  switch (MoL_Intermediate_Steps - (*MoL_Intermediate_Step))
  {
    case 0:
      alpha = 1.0 / 3.0;
      beta  = 0.5;
      break;
    case 1:
      alpha = 2.0 / 3.0;
      beta  = 0.5;
      break;
    case 2:
      alpha = 1.0 / 3.0;
      beta  = 1.0;
      break;
    case 3:
      alpha = 1.0;
      beta  = 1.0 / 6.0;
  }

  /* FIXME */

#ifdef MOLDOESCOMPLEX
  
  Complex_Delta_Time = CCTK_Cmplx((*Original_Delta_Time) / cctkGH->cctk_timefac, 0);
  Complex_beta = CCTK_Cmplx(beta, 0);

#endif

  /* Real GFs */

  for (var = 0; var < MoLNumEvolvedVariables; var++)
  {
    
    UpdateVar = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, 
                                              EvolvedVariableIndex[var]);
    OldVar = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 1, 
                                              EvolvedVariableIndex[var]);
    RHSVar = (CCTK_REAL const *)CCTK_VarDataPtrI(cctkGH, 0, 
                                                 RHSVariableIndex[var]);
/* #define MOLDEBUG 1 */
#ifdef MOLDEBUG
    printf("In generic RK. Variable %d (%s). RHS %d (%s). beta %g.\n",
           EvolvedVariableIndex[var],
           CCTK_VarName(EvolvedVariableIndex[var]),
           RHSVariableIndex[var],
           CCTK_VarName(RHSVariableIndex[var]),
           beta);
#endif

#pragma omp parallel for
    for (index = 0; index < totalsize; index++)
    {
      UpdateVar[index] = OldVar[index] + 
        (*Original_Delta_Time) / cctkGH->cctk_timefac * beta * RHSVar[index];
#ifdef MOLDEBUG
      if (CCTK_EQUALS(verbose,"extreme"))
      {
        printf("Variable: %d. Index: %d. dt: %f. beta %f. RHS: %f. q: %f.\n",
               var, index, (*Original_Delta_Time) / cctkGH->cctk_timefac, beta, RHSVar[index], 
               UpdateVar[index]);
      }
#endif
    }
    ScratchVar = CCTK_VarDataPtrI(cctkGH, 0, 
                                      scratchspace_firstindex
                                      + var );
    
    /* scratch storage */
    if ((*MoL_Intermediate_Step) == MoL_Intermediate_Steps)
    {
#pragma omp parallel for
      for (index = 0; index < totalsize; index++)
      {
        ScratchVar[index] = 0;
      }
    }

    if ((*MoL_Intermediate_Step)>1)
    {
#pragma omp parallel for
      for (index = 0; index < totalsize; index++)
      {
        ScratchVar[index] += alpha * UpdateVar[index];
      }
    }
    else
    {
#pragma omp parallel for
      for (index = 0; index < totalsize; index++)
      {
        UpdateVar[index] += ScratchVar[index] - 4.0 / 3.0 * OldVar[index];
      }
    }

  }


  /* Complex GFs */

  /* FIXME */

#ifdef MOLDOESCOMPLEX
CCTK_WARN(0, "not done");
#endif

  /* Real arrays */

  arrayscratchlocation = 0;

  for (var = 0; var < MoLNumEvolvedArrayVariables; var++)
  {
    
    UpdateVar = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, 
                                              EvolvedArrayVariableIndex[var]);
    OldVar = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 1, 
                                              EvolvedArrayVariableIndex[var]);
    RHSVar = (CCTK_REAL const *)CCTK_VarDataPtrI(cctkGH, 0, 
                                                 RHSArrayVariableIndex[var]);
    
    groupindex = CCTK_GroupIndexFromVarI(EvolvedArrayVariableIndex[var]);
    ierr = CCTK_GroupDynamicData(cctkGH, groupindex,
                                 &arraydata);
    if (ierr)
    {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, 
                 "The driver does not return group information "
                 "for group '%s'.", 
                 CCTK_GroupName(groupindex));
    }
    arraytotalsize = 1;
    for (arraydim = 0; arraydim < arraydata.dim; arraydim++)
    {
      arraytotalsize *= arraydata.lsh[arraydim];
    }

    ScratchVar = &ArrayScratchSpace[arrayscratchlocation];

#pragma omp parallel for
    for (index = 0; index < arraytotalsize; index++)
    {
      UpdateVar[index] = OldVar[index] +
        (*Original_Delta_Time) / cctkGH->cctk_timefac * beta * RHSVar[index];
    }

    if ((*MoL_Intermediate_Step) == MoL_Intermediate_Steps)
    {
#pragma omp parallel for
      for (index = 0; index < arraytotalsize; index++)
      {
        ScratchVar[index] = 0;
      }
    }
    
    if ((*MoL_Intermediate_Step)>1)
    {
#pragma omp parallel for
      for (index = 0; index < arraytotalsize; index++)
      {
        ScratchVar[index] += alpha * UpdateVar[index];
      }
    }
    else
    {
#pragma omp parallel for
      for (index = 0; index < arraytotalsize; index++)
      {
        UpdateVar[index] += ScratchVar[index] - 4.0 / 3.0 * OldVar[index];
      }
    }
    arrayscratchlocation +=  arraytotalsize;
  }

  return;
}
