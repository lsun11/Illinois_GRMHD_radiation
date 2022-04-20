/*@@
   @file       ScheduleIDScalarWaveMoL.c
   @author     Automatically generated by CreateScheduleBindings.pl
   @desc
               Creates the schedule and parameter recovery bindings 
               for thorn IDScalarWaveMoL
   @enddesc
@@*/

#define THORN_IS_IDScalarWaveMoL

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_ScheduleBindings.h"

/* prototypes for schedule bindings functions to be registered */
/* Note that this is a cheat, we just need a function pointer. */
extern int CCTK_FNAME(IDScalarWaveMoL_InitialData)(void);
extern int CCTK_FNAME(IDScalarWaveMoL_Errors)(void);


void CCTKi_BindingsSchedule_IDScalarWaveMoL(void);
void CCTKi_BindingsSchedule_IDScalarWaveMoL(void)
{
  DECLARE_CCTK_PARAMETERS
  {
    int cctkschedulei_tlevelarray[] = {1,0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(IDScalarWaveMoL_InitialData),
                           "IDScalarWaveMoL_InitialData",
                           "IDScalarWaveMoL",
                           "IDScalarWaveMoL",
                           "Initial data for the scalar field",
                           "CCTK_INITIAL",
                           "Fortran",
                           1,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           0,  /* Number of SYNC     groups   */
                           0,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           0,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "WaveToyMoL::scalarevolve");
  }
  {
    int cctkschedulei_tlevelarray[] = {1,1,0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(IDScalarWaveMoL_Errors),
                           "IDScalarWaveMoL_Errors",
                           "IDScalarWaveMoL",
                           "IDScalarWaveMoL",
                           "Calculate errors of the scalar field",
                           "CCTK_ANALYSIS",
                           "Fortran",
                           2,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           1,  /* Number of TRIGGERS groups   */
                           0,  /* Number of SYNC     groups   */
                           0,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           0,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "IDScalarWaveMoL::scalarevolveerror",
                           "WaveToyMoL::scalarevolve",
                           "IDScalarWaveMoL::scalarevolveerror");
  }
}

/*@@
  @routine    CCTKi_BindingsParameterRecovery_IDScalarWaveMoL
  @author     Automatically generated by CreateScheduleBindings.pl
  @desc
              Creates the parameter recovery bindings for thorn IDScalarWaveMoL
  @enddesc
@@*/

int CCTKi_BindingsParameterRecovery_IDScalarWaveMoL(void);
int CCTKi_BindingsParameterRecovery_IDScalarWaveMoL(void)
{
  /* this thorn doesn't define any parameter recovery routines */
  return (0);
}

