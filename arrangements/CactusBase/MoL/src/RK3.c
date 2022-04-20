 /*@@
   @file      RK3.c
   @date      Tue Jul 22 00:38:47 2003
   @author    Ian Hawke
   @desc 
   A specialized third order Runge-Kutta time integrator. This is
   the integrator that Shu refers to as the optimal TVD third 
   order method (see reference in documentation).
   @enddesc 
   @version   $Header: /cactusdevcvs/CactusBase/MoL/src/RK3.c,v 1.8 2008/09/23 15:46:30 schnetter Exp $
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "ExternalVariables.h"

static const char *rcsid = "$Header: /cactusdevcvs/CactusBase/MoL/src/RK3.c,v 1.8 2008/09/23 15:46:30 schnetter Exp $";

CCTK_FILEVERSION(CactusBase_MoL_RK3_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MoL_RK3Add(CCTK_ARGUMENTS);

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
   @routine    MoL_RK3Add
   @date       Tue Jul 22 00:39:55 2003
   @author     Ian Hawke
   @desc 
   Performs third order Runge-Kutta time integration.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

 @@*/

void MoL_RK3Add(CCTK_ARGUMENTS)
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
  CCTK_COMPLEX Complex_Third = CCTK_Cmplx(1.0/3.0, 0);
  CCTK_COMPLEX Complex_TwoThird = CCTK_Cmplx(2.0/3.0, 0);
  CCTK_COMPLEX Complex_Quarter = CCTK_Cmplx(0.25, 0);
  CCTK_COMPLEX Complex_ThreeQuarter = CCTK_Cmplx(0.75, 0);

#endif

#ifdef MOLDEBUG
  printf("Inside RK3.\nStep %d.\nRefinement %d.\nTimestep %g.\n"
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

    case 3:
      {
        for (var = 0; var < MoLNumEvolvedVariables; var++)
        {
          UpdateVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 0,
                                                   EvolvedVariableIndex[var]);
          RHSVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 0, 
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
          RHSVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 0, 
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
          RHSComplexVar = (CCTK_COMPLEX*)CCTK_VarDataPtrI(cctkGH, 0, 
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
  
    case 2:
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
            UpdateVar[index] = 0.25 * (3*OldVar[index] +
                                       UpdateVar[index]) + 
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
            UpdateVar[index] = 0.25*(3*OldVar[index] + UpdateVar[index]) +
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
              CCTK_CmplxAdd(CCTK_CmplxAdd(CCTK_CmplxMul(Complex_ThreeQuarter,
                                                        OldComplexVar[index]),
                                          CCTK_CmplxMul(ComplexQuarter,
                                                        UpdateComplexVar[index])), 
                            CCTK_CmplxMul(CCTK_CmplxMul(Complex_Delta_Time, 
                                                        RHSComplexVar[index]))); 
          }
        }

#endif

        break;
      }
    case 1:
      {
        const CCTK_REAL one_third = 1.0 / 3.0;

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
            UpdateVar[index] = (OldVar[index] + 2*UpdateVar[index]) * one_third
              + CCTK_DELTA_TIME * RHSVar[index];
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
            UpdateVar[index] = (OldVar[index] + 2*UpdateVar[index]) * one_third
              + CCTK_DELTA_TIME * RHSVar[index];
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
              CCTK_CmplxAdd(CCTK_CmplxAdd(CCTK_CmplxMul(Complex_Third,
                                                        OldComplexVar[index]),
                                          CCTK_CmplxMul(ComplexTwoThird,
                                                        UpdateComplexVar[index])), 
                            CCTK_CmplxMul(Complex_Delta_Time, 
                                          RHSComplexVar[index])); 
          }
        }

#endif

        break;
      }
    default:
      {
        CCTK_WARN(0, "RK3 expects MoL_Intermediate_Step to be "
                  "in [1,3]. This should be caught at ParamCheck - bug Ian!");
        break;
      }
      
  }

  return;

}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
