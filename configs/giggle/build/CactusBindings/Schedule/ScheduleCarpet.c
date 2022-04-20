/*@@
   @file       ScheduleCarpet.c
   @author     Automatically generated by CreateScheduleBindings.pl
   @desc
               Creates the schedule and parameter recovery bindings 
               for thorn Carpet
   @enddesc
@@*/

#define THORN_IS_Carpet

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_ScheduleBindings.h"

/* prototypes for schedule bindings functions to be registered */
/* Note that this is a cheat, we just need a function pointer. */
extern int CarpetMultiModelStartup(void);
extern int CarpetStartup(void);
extern int CarpetParamCheck(void);
extern int CarpetRefineTimeStep(void);


void CCTKi_BindingsSchedule_Carpet(void);
void CCTKi_BindingsSchedule_Carpet(void)
{
  DECLARE_CCTK_PARAMETERS
  CCTKi_ScheduleGroupStorage("Carpet::timing",1);
  CCTKi_ScheduleGroupStorage("Carpet::timing2",1);
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CarpetMultiModelStartup,
                           "MultiModel_Startup",
                           "Carpet",
                           "Driver",
                           "Multi-model Startup routine",
                           "CCTK_STARTUP",
                           "C",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           0,  /* Number of SYNC     groups   */
                           0,  /* Number of Options           */
                           1,  /* Number of BEFORE  routines  */
                           0,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "Driver_Startup");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CarpetStartup,
                           "Driver_Startup",
                           "Carpet",
                           "Driver",
                           "Startup routine",
                           "CCTK_STARTUP",
                           "C",
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
    CCTKi_ScheduleFunction((void *)CarpetParamCheck,
                           "CarpetParamCheck",
                           "Carpet",
                           "Driver",
                           "Parameter checking routine",
                           "CCTK_PARAMCHECK",
                           "C",
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
if (refine_timestep)
{
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CarpetRefineTimeStep,
                           "CarpetRefineTimeStep",
                           "Carpet",
                           "Driver",
                           "Correct time step size for spacing on finer grids",
                           "CCTK_BASEGRID",
                           "C",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           0,  /* Number of SYNC     groups   */
                           1,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           1,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "singlemap",
                           "Time_Simple");
  }
}
}

/*@@
  @routine    CCTKi_BindingsParameterRecovery_Carpet
  @author     Automatically generated by CreateScheduleBindings.pl
  @desc
              Creates the parameter recovery bindings for thorn Carpet
  @enddesc
@@*/

int CCTKi_BindingsParameterRecovery_Carpet(void);
int CCTKi_BindingsParameterRecovery_Carpet(void)
{
  /* this thorn doesn't define any parameter recovery routines */
  return (0);
}
