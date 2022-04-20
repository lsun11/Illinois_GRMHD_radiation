/*@@
   @file       ScheduleWaveToyMoL.c
   @author     Automatically generated by CreateScheduleBindings.pl
   @desc
               Creates the schedule and parameter recovery bindings 
               for thorn WaveToyMoL
   @enddesc
@@*/

#define THORN_IS_WaveToyMoL

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_ScheduleBindings.h"

/* prototypes for schedule bindings functions to be registered */
/* Note that this is a cheat, we just need a function pointer. */
extern int wavetoymol_startup_(void);
extern int wavetoymol_initsymbound_(void);
extern int wavetoymol_registervars_(void);
extern int wavetoymol_calcrhs_(void);
extern int wavetoymol_boundaries_(void);
extern int wavetoymol_estimateerror_(void);
extern int wavetoymol_estimateerrorboundaries_(void);


void CCTKi_BindingsSchedule_WaveToyMoL(void);
void CCTKi_BindingsSchedule_WaveToyMoL(void)
{
  DECLARE_CCTK_PARAMETERS
  CCTKi_ScheduleGroupStorage("WaveToyMoL::scalarevolve",3);
  CCTKi_ScheduleGroupStorage("WaveToyMoL::scalarevolvedot",1);
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)wavetoymol_startup_,
                           "WaveToyMol_Startup",
                           "WaveToyMoL",
                           "WaveToyMoL",
                           "Register banner",
                           "CCTK_STARTUP",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           0,  /* Number of SYNC     groups   */
                           1,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           0,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "meta");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)wavetoymol_initsymbound_,
                           "WaveToyMoL_InitSymBound",
                           "WaveToyMoL",
                           "WaveToyMoL",
                           "Schedule symmetries",
                           "CCTK_BASEGRID",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           0,  /* Number of SYNC     groups   */
                           1,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           0,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "meta");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)wavetoymol_registervars_,
                           "WaveToyMoL_RegisterVars",
                           "WaveToyMoL",
                           "WaveToyMoL",
                           "Register variables for MoL",
                           "MoL_Register",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           0,  /* Number of SYNC     groups   */
                           1,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           0,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "meta");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)wavetoymol_calcrhs_,
                           "WaveToyMoL_CalcRHS",
                           "WaveToyMoL",
                           "WaveToyMoL",
                           "Calculate RHS for MoL",
                           "MoL_CalcRHS",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           0,  /* Number of SYNC     groups   */
                           0,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           0,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */);
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)wavetoymol_boundaries_,
                           "WaveToyMoL_Boundaries",
                           "WaveToyMoL",
                           "WaveToyMoL",
                           "Select boundary conditions in MoL",
                           "MoL_PostStep",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           1,  /* Number of SYNC     groups   */
                           1,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           0,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "WaveToyMoL::scalarevolve",
                           "level");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleGroup("ApplyBCs",
"ApplyBCs",
                        "WaveToyMoL",
                        "WaveToyMoL",
                        "Apply boundary conditions in MoL",
                        "MoL_PostStep",
                        0,  /* Number of STORAGE  groups   */
                        0,  /* Number of COMM     groups   */
                        0,  /* Number of TRIGGERS groups   */
                        0,  /* Number of SYNC     groups   */
                        0,  /* Number of Options           */
                        0,  /* Number of BEFORE  routines  */
                        1,  /* Number of AFTER   routines  */
                        0,  /* Number of WHILE   variables */
                        0,  /* Number of IF   variables */
                        cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                        "WaveToyMoL_Boundaries");
  }
if (estimate_error)
{
  CCTKi_ScheduleGroupStorage("WaveToyMoL::scalarevolveerrorestimate",1);
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)wavetoymol_estimateerror_,
                           "WaveToyMoL_EstimateError",
                           "WaveToyMoL",
                           "WaveToyMoL",
                           "Estimate the truncation error",
                           "CCTK_POSTSTEP",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           0,  /* Number of SYNC     groups   */
                           0,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           0,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */);
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)wavetoymol_estimateerrorboundaries_,
                           "WaveToyMoL_EstimateErrorBoundaries",
                           "WaveToyMoL",
                           "WaveToyMoL",
                           "Select boundary conditions for the truncation error",
                           "CCTK_POSTSTEP",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           1,  /* Number of SYNC     groups   */
                           1,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           1,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "WaveToyMoL::scalarevolveerrorestimate",
                           "level",
                           "WaveToyMoL_EstimateError");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleGroup("ApplyBCs",
"ApplyBCs",
                        "WaveToyMoL",
                        "WaveToyMoL",
                        "Apply boundary conditions",
                        "CCTK_POSTSTEP",
                        0,  /* Number of STORAGE  groups   */
                        0,  /* Number of COMM     groups   */
                        0,  /* Number of TRIGGERS groups   */
                        0,  /* Number of SYNC     groups   */
                        0,  /* Number of Options           */
                        0,  /* Number of BEFORE  routines  */
                        1,  /* Number of AFTER   routines  */
                        0,  /* Number of WHILE   variables */
                        0,  /* Number of IF   variables */
                        cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                        "WaveToyMoL_EstimateError");
  }
}
}

/*@@
  @routine    CCTKi_BindingsParameterRecovery_WaveToyMoL
  @author     Automatically generated by CreateScheduleBindings.pl
  @desc
              Creates the parameter recovery bindings for thorn WaveToyMoL
  @enddesc
@@*/

int CCTKi_BindingsParameterRecovery_WaveToyMoL(void);
int CCTKi_BindingsParameterRecovery_WaveToyMoL(void)
{
  /* this thorn doesn't define any parameter recovery routines */
  return (0);
}
