/*@@
   @file       ScheduleCarpetEvolutionMask.c
   @author     Automatically generated by CreateScheduleBindings.pl
   @desc
               Creates the schedule and parameter recovery bindings 
               for thorn CarpetEvolutionMask
   @enddesc
@@*/

#define THORN_IS_CarpetEvolutionMask

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_ScheduleBindings.h"

/* prototypes for schedule bindings functions to be registered */
/* Note that this is a cheat, we just need a function pointer. */
extern int EvolutionMaskBase_InitEvolutionMask(void);
extern int CarpetEvolutionMaskSetup(void);
extern int enforce_evolution_mask(void);


void CCTKi_BindingsSchedule_CarpetEvolutionMask(void);
void CCTKi_BindingsSchedule_CarpetEvolutionMask(void)
{
  DECLARE_CCTK_PARAMETERS
  CCTKi_ScheduleGroupStorage("CarpetEvolutionMask::evolution_mask",1);
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleGroup("EvolutionMaskBase_SetupEvolutionMask",
"EvolutionMaskBase_SetupEvolutionMask",
                        "CarpetEvolutionMask",
                        "CarpetEvolutionMask",
                        "Set up the mask function",
                        "CCTK_BASEGRID",
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
    CCTKi_ScheduleGroup("EvolutionMaskBase_SetupEvolutionMask",
"EvolutionMaskBase_SetupEvolutionMask",
                        "CarpetEvolutionMask",
                        "CarpetEvolutionMask",
                        "Set up the mask function",
                        "CCTK_POSTREGRIDINITIAL",
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
    CCTKi_ScheduleGroup("EvolutionMaskBase_SetupEvolutionMask",
"EvolutionMaskBase_SetupEvolutionMask",
                        "CarpetEvolutionMask",
                        "CarpetEvolutionMask",
                        "Set up the mask function",
                        "CCTK_POSTREGRID",
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
    CCTKi_ScheduleFunction((void *)EvolutionMaskBase_InitEvolutionMask,
                           "EvolutionMaskBase_InitEvolutionMask",
                           "CarpetEvolutionMask",
                           "CarpetEvolutionMask",
                           "Initialise the mask function",
                           "EvolutionMaskBase_SetupEvolutionMask",
                           "C",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           0,  /* Number of SYNC     groups   */
                           2,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           0,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "global",
                           "loop-local");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleGroup("SetupEvolutionMask",
"SetupEvolutionMask",
                        "CarpetEvolutionMask",
                        "CarpetEvolutionMask",
                        "Set up the weight function (schedule other routines in here)",
                        "EvolutionMaskBase_SetupEvolutionMask",
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
                        "EvolutionMaskBase_InitEvolutionMask");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CarpetEvolutionMaskSetup,
                           "CarpetEvolutionMaskSetup",
                           "CarpetEvolutionMask",
                           "CarpetEvolutionMask",
                           "Set up the mask function for the restriction regions",
                           "SetupEvolutionMask",
                           "C",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           0,  /* Number of SYNC     groups   */
                           2,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           0,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "global",
                           "loop-singlemap");
  }
if (enforce_mask)
{
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)enforce_evolution_mask,
                           "enforce_evolution_mask",
                           "CarpetEvolutionMask",
                           "CarpetEvolutionMask",
                           "Enforce Evolution Mask",
                           "MoL_PostRHS",
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
}

/*@@
  @routine    CCTKi_BindingsParameterRecovery_CarpetEvolutionMask
  @author     Automatically generated by CreateScheduleBindings.pl
  @desc
              Creates the parameter recovery bindings for thorn CarpetEvolutionMask
  @enddesc
@@*/

int CCTKi_BindingsParameterRecovery_CarpetEvolutionMask(void);
int CCTKi_BindingsParameterRecovery_CarpetEvolutionMask(void)
{
  /* this thorn doesn't define any parameter recovery routines */
  return (0);
}
