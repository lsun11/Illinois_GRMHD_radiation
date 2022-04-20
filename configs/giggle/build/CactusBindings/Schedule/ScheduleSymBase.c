/*@@
   @file       ScheduleSymBase.c
   @author     Automatically generated by CreateScheduleBindings.pl
   @desc
               Creates the schedule and parameter recovery bindings 
               for thorn SymBase
   @enddesc
@@*/

#define THORN_IS_SymBase

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_ScheduleBindings.h"

/* prototypes for schedule bindings functions to be registered */
/* Note that this is a cheat, we just need a function pointer. */
extern int SymBase_Startup(void);
extern int SymBase_Statistics(void);
extern int SymBase_Check(void);


void CCTKi_BindingsSchedule_SymBase(void);
void CCTKi_BindingsSchedule_SymBase(void)
{
  DECLARE_CCTK_PARAMETERS
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)SymBase_Startup,
                           "SymBase_Startup",
                           "SymBase",
                           "SymBase",
                           "Register GH Extension for SymBase",
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
    CCTKi_ScheduleGroup("SymBase_Wrapper",
"SymBase_Wrapper",
                        "SymBase",
                        "SymBase",
                        "Wrapper group for SymBase",
                        "CCTK_WRAGH",
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
    CCTKi_ScheduleGroup("SymmetryRegister",
"SymmetryRegister",
                        "SymBase",
                        "SymBase",
                        "Register your symmetries here",
                        "SymBase_Wrapper",
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
    CCTKi_ScheduleFunction((void *)SymBase_Statistics,
                           "SymBase_Statistics",
                           "SymBase",
                           "SymBase",
                           "Print symmetry boundary face descriptions",
                           "SymBase_Wrapper",
                           "C",
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
                           "SymmetryRegister");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)SymBase_Check,
                           "SymBase_Check",
                           "SymBase",
                           "SymBase",
                           "Check whether the driver set up the grid consistently",
                           "CCTK_BASEGRID",
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
  @routine    CCTKi_BindingsParameterRecovery_SymBase
  @author     Automatically generated by CreateScheduleBindings.pl
  @desc
              Creates the parameter recovery bindings for thorn SymBase
  @enddesc
@@*/

int CCTKi_BindingsParameterRecovery_SymBase(void);
int CCTKi_BindingsParameterRecovery_SymBase(void)
{
  /* this thorn doesn't define any parameter recovery routines */
  return (0);
}
