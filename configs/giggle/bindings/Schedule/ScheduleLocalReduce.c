/*@@
   @file       ScheduleLocalReduce.c
   @author     Automatically generated by CreateScheduleBindings.pl
   @desc
               Creates the schedule and parameter recovery bindings 
               for thorn LocalReduce
   @enddesc
@@*/

#define THORN_IS_LocalReduce

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_ScheduleBindings.h"

/* prototypes for schedule bindings functions to be registered */
/* Note that this is a cheat, we just need a function pointer. */
extern int LocalReduce_Startup(void);


void CCTKi_BindingsSchedule_LocalReduce(void);
void CCTKi_BindingsSchedule_LocalReduce(void)
{
  DECLARE_CCTK_PARAMETERS
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)LocalReduce_Startup,
                           "LocalReduce_Startup",
                           "LocalReduce",
                           "LocalReduce",
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
}

/*@@
  @routine    CCTKi_BindingsParameterRecovery_LocalReduce
  @author     Automatically generated by CreateScheduleBindings.pl
  @desc
              Creates the parameter recovery bindings for thorn LocalReduce
  @enddesc
@@*/

int CCTKi_BindingsParameterRecovery_LocalReduce(void);
int CCTKi_BindingsParameterRecovery_LocalReduce(void)
{
  /* this thorn doesn't define any parameter recovery routines */
  return (0);
}
