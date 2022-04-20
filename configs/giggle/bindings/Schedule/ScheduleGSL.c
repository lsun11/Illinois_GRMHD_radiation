/*@@
   @file       ScheduleGSL.c
   @author     Automatically generated by CreateScheduleBindings.pl
   @desc
               Creates the schedule and parameter recovery bindings 
               for thorn GSL
   @enddesc
@@*/

#define THORN_IS_GSL

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_ScheduleBindings.h"

/* prototypes for schedule bindings functions to be registered */
/* Note that this is a cheat, we just need a function pointer. */


void CCTKi_BindingsSchedule_GSL(void);
void CCTKi_BindingsSchedule_GSL(void)
{
  DECLARE_CCTK_PARAMETERS
}

/*@@
  @routine    CCTKi_BindingsParameterRecovery_GSL
  @author     Automatically generated by CreateScheduleBindings.pl
  @desc
              Creates the parameter recovery bindings for thorn GSL
  @enddesc
@@*/

int CCTKi_BindingsParameterRecovery_GSL(void);
int CCTKi_BindingsParameterRecovery_GSL(void)
{
  /* this thorn doesn't define any parameter recovery routines */
  return (0);
}

