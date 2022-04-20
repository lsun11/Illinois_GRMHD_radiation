/*@@
   @file       ScheduleCarpetRegrid2.c
   @author     Automatically generated by CreateScheduleBindings.pl
   @desc
               Creates the schedule and parameter recovery bindings 
               for thorn CarpetRegrid2
   @enddesc
@@*/

#define THORN_IS_CarpetRegrid2

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_ScheduleBindings.h"

/* prototypes for schedule bindings functions to be registered */
/* Note that this is a cheat, we just need a function pointer. */
extern int CarpetRegrid2_ParamCheck(void);
extern int CarpetRegrid2_Initialise(void);


void CCTKi_BindingsSchedule_CarpetRegrid2(void);
void CCTKi_BindingsSchedule_CarpetRegrid2(void)
{
  DECLARE_CCTK_PARAMETERS
  CCTKi_ScheduleGroupStorage("CarpetRegrid2::last_iteration",1);
  CCTKi_ScheduleGroupStorage("CarpetRegrid2::last_map",1);
  CCTKi_ScheduleGroupStorage("CarpetRegrid2::active",1);
  CCTKi_ScheduleGroupStorage("CarpetRegrid2::num_levels",1);
  CCTKi_ScheduleGroupStorage("CarpetRegrid2::positions",1);
  CCTKi_ScheduleGroupStorage("CarpetRegrid2::radii",1);
  CCTKi_ScheduleGroupStorage("CarpetRegrid2::old_active",1);
  CCTKi_ScheduleGroupStorage("CarpetRegrid2::old_positions",1);
  CCTKi_ScheduleGroupStorage("CarpetRegrid2::old_num_levels",1);
  CCTKi_ScheduleGroupStorage("CarpetRegrid2::old_radii",1);
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CarpetRegrid2_ParamCheck,
                           "CarpetRegrid2_ParamCheck",
                           "CarpetRegrid2",
                           "CarpetRegrid2",
                           "Check parameters",
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
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CarpetRegrid2_Initialise,
                           "CarpetRegrid2_Initialise",
                           "CarpetRegrid2",
                           "CarpetRegrid2",
                           "Initialise locations of refined regions",
                           "CCTK_WRAGH",
                           "C",
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
                           "global");
  }
}

/*@@
  @routine    CCTKi_BindingsParameterRecovery_CarpetRegrid2
  @author     Automatically generated by CreateScheduleBindings.pl
  @desc
              Creates the parameter recovery bindings for thorn CarpetRegrid2
  @enddesc
@@*/

int CCTKi_BindingsParameterRecovery_CarpetRegrid2(void);
int CCTKi_BindingsParameterRecovery_CarpetRegrid2(void)
{
  /* this thorn doesn't define any parameter recovery routines */
  return (0);
}