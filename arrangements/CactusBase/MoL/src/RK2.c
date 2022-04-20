 /*@@
   @file      RK2.c
   @date      Sun May 26 04:13:45 2002
   @author    Ian Hawke
   @desc 
   A specialized second order Runge-Kutta time integrator. This is
   the integrator that Shu refers to as the optimal TVD second 
   order method (see reference in documentation). It is equivalent
   to Heun's predictor-corrector method, or the MacCormack method.
   @enddesc 
   @version   $Header: /cactusdevcvs/CactusBase/MoL/src/RK2.c,v 1.12 2008/09/23 15:46:30 schnetter Exp $
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "ExternalVariables.h"

static const char *rcsid = "$Header: /cactusdevcvs/CactusBase/MoL/src/RK2.c,v 1.12 2008/09/23 15:46:30 schnetter Exp $";

CCTK_FILEVERSION(CactusBase_MoL_RK2_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MoL_RK2Add(CCTK_ARGUMENTS);

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
   @routine    MoL_RK2Add
   @date       Sun May 26 04:17:23 2002
   @author     Ian Hawke
   @desc 
   Performs second order Runge-Kutta time integration.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_RK2Add(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  cGroupDynamicData arraydata;
  CCTK_INT groupindex, ierr;
  CCTK_INT arraytotalsize, arraydim;
  
  CCTK_INT index, var;
  CCTK_INT totalsize;
  CCTK_REAL const * restrict OldVar;
  CCTK_REAL       * restrict UpdateVar;
  CCTK_REAL const * restrict RHSVar;
  
  /* FIXME */

#ifdef MOLDOESCOMPLEX

  CCTK_COMPLEX const * restrict OldComplexVar;
  CCTK_COMPLEX       * restrict UpdateComplexVar;
  CCTK_COMPLEX const * restrict RHSComplexVar;
  CCTK_COMPLEX Complex_Delta_Time = CCTK_Cmplx(CCTK_DELTA_TIME, 0);
  CCTK_COMPLEX Complex_Half = CCTK_Cmplx(0.5, 0);

#endif

#ifdef MOLDEBUG
  printf("Inside RK2.\nStep %d.\nRefinement %d.\nTimestep %g.\n"
         "Spacestep %g.\nTime %g\n",
         MoL_Intermediate_Steps - *MoL_Intermediate_Step + 1,
         *cctk_levfac,
         CCTK_DELTA_TIME,
         CCTK_DELTA_SPACE(0),
         cctk_time);
#endif  

  totalsize = 1;
  for (arraydim = 0; arraydim < cctk_dim; arraydim++)
  {
    totalsize *= cctk_lsh[arraydim];
  }

  switch (*MoL_Intermediate_Step)
  {
  
    case 2:
      {
        for (var = 0; var < MoLNumEvolvedVariables; var++)
        {
          UpdateVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 0,
                                                   EvolvedVariableIndex[var]);
          RHSVar = (CCTK_REAL const*)CCTK_VarDataPtrI(cctkGH, 0, 
                                                RHSVariableIndex[var]);
          
#pragma omp parallel for
          for (index = 0; index < totalsize; index++)
          {
            UpdateVar[index] += CCTK_DELTA_TIME * RHSVar[index];
          }
        }

        for (var = 0; var < MoLNumEvolvedArrayVariables; var++)
        {
          UpdateVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 0,
                                                   EvolvedArrayVariableIndex[var]);
          RHSVar = (CCTK_REAL const*)CCTK_VarDataPtrI(cctkGH, 0, 
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

#pragma omp parallel for
          for (index = 0; index < arraytotalsize; index++)
          {
            UpdateVar[index] += CCTK_DELTA_TIME * RHSVar[index];
          }
        }

  /* FIXME */

#ifdef MOLDOESCOMPLEX

        for (var = 0; var < MoLNumEvolvedComplexVariables; var++)
        {
          UpdateComplexVar = (CCTK_COMPLEX*)CCTK_VarDataPtrI(cctkGH, 0,
                                                             EvolvedComplexVariableIndex[var]);
          RHSComplexVar = (CCTK_COMPLEX const*)CCTK_VarDataPtrI(cctkGH, 0, 
                                                          RHSComplexVariableIndex[var]);
          
#pragma omp parallel for
          for (index = 0; index < totalsize; index++)
          {
            UpdateComplexVar[index] = CCTK_CmplxAdd(UpdateComplexVar[index],
                                                    CCTK_CmplxMul(Complex_Delta_Time,
                                                                  RHSComplexVar[index]));
          }
        }

#endif

        break;
      }
    case 1:
      {
        for (var = 0; var < MoLNumEvolvedVariables; var++)
        {
          OldVar = (CCTK_REAL const*)CCTK_VarDataPtrI(cctkGH, 1,
                                                EvolvedVariableIndex[var]);
          UpdateVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 0,
                                                   EvolvedVariableIndex[var]);
          RHSVar = (CCTK_REAL const*)CCTK_VarDataPtrI(cctkGH, 0, 
                                                RHSVariableIndex[var]);
          
#pragma omp parallel for
          for (index = 0; index < totalsize; index++)
          {
            UpdateVar[index] = 0.5 * (OldVar[index] + UpdateVar[index]) + 
                                      CCTK_DELTA_TIME * RHSVar[index];
          }
        }

        for (var = 0; var < MoLNumEvolvedArrayVariables; var++)
        {
          OldVar = (CCTK_REAL const*)CCTK_VarDataPtrI(cctkGH, 1,
                                                   EvolvedArrayVariableIndex[var]);
          UpdateVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 0,
                                                   EvolvedArrayVariableIndex[var]);
          RHSVar = (CCTK_REAL const*)CCTK_VarDataPtrI(cctkGH, 0, 
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
          
#pragma omp parallel for
          for (index = 0; index < arraytotalsize; index++)
          {
            UpdateVar[index] = 0.5 * (OldVar[index] + UpdateVar[index]) + 
                                      CCTK_DELTA_TIME * RHSVar[index];
          }
        }

  /* FIXME */

#ifdef MOLDOESCOMPLEX

        for (var = 0; var < MoLNumEvolvedComplexVariables; var++)
        {
          OldComplexVar = (CCTK_COMPLEX const*)CCTK_VarDataPtrI(cctkGH, 1,
                                                       EvolvedComplexVariableIndex[var]);
          UpdateComplexVar = (CCTK_COMPLEX*)CCTK_VarDataPtrI(cctkGH, 0,
                                                             EvolvedComplexVariableIndex[var]);
          RHSComplexVar = (CCTK_COMPLEX const*)CCTK_VarDataPtrI(cctkGH, 0, 
                                                          RHSComplexVariableIndex[var]);
          
#pragma omp parallel for
          for (index = 0; index < totalsize; index++)
          {
            UpdateComplexVar[index] =
              CCTK_CmplxAdd(CCTK_CmplxMul(Complex_Half,
                                          (CCTK_CmplxAdd(OldComplexVar[index],
                                                         UpdateComplexVar[index]))), 
                            CCTK_CmplxMul(Complex_Delta_Time, 
                                          RHSComplexVar[index])); 
          }
        }

#endif

        break;
      }
    default:
      {
        CCTK_WARN(0, "RK2 expects MoL_Intermediate_Step to be "
                  "in [1,2]. This should be caught at ParamCheck - "
                  "bug Ian!");
        break;
      }
      
  }

  return;

}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
