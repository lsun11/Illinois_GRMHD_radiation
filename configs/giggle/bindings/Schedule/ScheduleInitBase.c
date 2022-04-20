/*@@
   @file       ScheduleInitBase.c
   @author     Automatically generated by CreateScheduleBindings.pl
   @desc
               Creates the schedule and parameter recovery bindings 
               for thorn InitBase
   @enddesc
@@*/

#define THORN_IS_InitBase

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_ScheduleBindings.h"

/* prototypes for schedule bindings functions to be registered */
/* Note that this is a cheat, we just need a function pointer. */


void CCTKi_BindingsSchedule_InitBase(void);
void CCTKi_BindingsSchedule_InitBase(void)
{
  DECLARE_CCTK_PARAMETERS
}

/*@@
  @routine    CCTKi_BindingsParameterRecovery_InitBase
  @author     Automatically generated by CreateScheduleBindings.pl
  @desc
              Creates the parameter recovery bindings for thorn InitBase
  @enddesc
@@*/

int CCTKi_BindingsParameterRecovery_InitBase(void);
int CCTKi_BindingsParameterRecovery_InitBase(void)
{
  /* this thorn doesn't define any parameter recovery routines */
  return (0);
}
