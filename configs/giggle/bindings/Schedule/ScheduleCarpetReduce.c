/*@@
   @file       ScheduleCarpetReduce.c
   @author     Automatically generated by CreateScheduleBindings.pl
   @desc
               Creates the schedule and parameter recovery bindings 
               for thorn CarpetReduce
   @enddesc
@@*/

#define THORN_IS_CarpetReduce

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_ScheduleBindings.h"

/* prototypes for schedule bindings functions to be registered */
/* Note that this is a cheat, we just need a function pointer. */
extern int CarpetReduceStartup(void);
extern int MaskBase_InitMask(void);
extern int CoordBase_SetupMask(void);
extern int CarpetMaskSetup(void);


void CCTKi_BindingsSchedule_CarpetReduce(void);
void CCTKi_BindingsSchedule_CarpetReduce(void)
{
  DECLARE_CCTK_PARAMETERS
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CarpetReduceStartup,
                           "CarpetReduceStartup",
                           "CarpetReduce",
                           "reduce",
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
  CCTKi_ScheduleGroupStorage("CarpetReduce::weight",1);
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleGroup("MaskBase_SetupMask",
"MaskBase_SetupMask",
                        "CarpetReduce",
                        "reduce",
                        "Set up the weight function",
                        "CCTK_BASEGRID",
                        0,  /* Number of STORAGE  groups   */
                        0,  /* Number of COMM     groups   */
                        0,  /* Number of TRIGGERS groups   */
                        0,  /* Number of SYNC     groups   */
                        0,  /* Number of Options           */
                        0,  /* Number of BEFORE  routines  */
                        2,  /* Number of AFTER   routines  */
                        0,  /* Number of WHILE   variables */
                        0,  /* Number of IF   variables */
                        cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                        "SpatialCoordinates",
                        "SphericalSurface_Setup");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleGroup("MaskBase_SetupMask",
"MaskBase_SetupMask",
                        "CarpetReduce",
                        "reduce",
                        "Set up the weight function",
                        "CCTK_POSTREGRIDINITIAL",
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
                        "SpatialCoordinates");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleGroup("MaskBase_SetupMask",
"MaskBase_SetupMask",
                        "CarpetReduce",
                        "reduce",
                        "Set up the weight function",
                        "CCTK_POSTREGRID",
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
                        "SpatialCoordinates");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleGroup("MaskBase_SetupMask",
"MaskBase_SetupMask",
                        "CarpetReduce",
                        "reduce",
                        "Set up the weight function",
                        "CCTK_POST_RECOVER_VARIABLES",
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
    CCTKi_ScheduleFunction((void *)MaskBase_InitMask,
                           "MaskBase_InitMask",
                           "CarpetReduce",
                           "reduce",
                           "Initialise the weight function",
                           "MaskBase_SetupMask",
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
    CCTKi_ScheduleGroup("SetupMask",
"SetupMask",
                        "CarpetReduce",
                        "reduce",
                        "Set up the weight function (schedule other routines in here)",
                        "MaskBase_SetupMask",
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
                        "MaskBase_InitMask");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CoordBase_SetupMask,
                           "CoordBase_SetupMask",
                           "CarpetReduce",
                           "reduce",
                           "Set up the outer boundaries of the weight function",
                           "SetupMask",
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
    CCTKi_ScheduleFunction((void *)CarpetMaskSetup,
                           "CarpetMaskSetup",
                           "CarpetReduce",
                           "reduce",
                           "Set up the weight function for the restriction regions",
                           "SetupMask",
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
}

/*@@
  @routine    CCTKi_BindingsParameterRecovery_CarpetReduce
  @author     Automatically generated by CreateScheduleBindings.pl
  @desc
              Creates the parameter recovery bindings for thorn CarpetReduce
  @enddesc
@@*/

int CCTKi_BindingsParameterRecovery_CarpetReduce(void);
int CCTKi_BindingsParameterRecovery_CarpetReduce(void)
{
  /* this thorn doesn't define any parameter recovery routines */
  return (0);
}
