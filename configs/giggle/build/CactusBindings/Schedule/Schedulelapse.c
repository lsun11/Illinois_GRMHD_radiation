/*@@
   @file       Schedulelapse.c
   @author     Automatically generated by CreateScheduleBindings.pl
   @desc
               Creates the schedule and parameter recovery bindings 
               for thorn lapse
   @enddesc
@@*/

#define THORN_IS_lapse

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_ScheduleBindings.h"

/* prototypes for schedule bindings functions to be registered */
/* Note that this is a cheat, we just need a function pointer. */
extern int Lapse_Startup(void);
extern int Lapse_CheckParameters(void);
extern int Lapse_InitSymBound(void);
extern int setup_initial_lapse_(void);
extern int lapse_postinitialdata_(void);
extern int lapse_RegisterVars(void);
extern int lapse_timestepping_(void);
extern int lapse_update_bc_(void);
extern int lapse_postbc_(void);


void CCTKi_BindingsSchedule_lapse(void);
void CCTKi_BindingsSchedule_lapse(void)
{
  DECLARE_CCTK_PARAMETERS
  CCTKi_ScheduleGroupStorage("lapse::lapse_vars",3);
  CCTKi_ScheduleGroupStorage("lapse::lapse_derivatives",1);
  CCTKi_ScheduleGroupStorage("lapse::lapse_vars_temp",1);
  CCTKi_ScheduleGroupStorage("lapse::lapse_vars_aux",1);
  CCTKi_ScheduleGroupStorage("lapse::lapse_rhs",1);
  CCTKi_ScheduleGroupStorage("BSSN::BSSN_vars",3);
  CCTKi_ScheduleGroupStorage("BSSN::BSSN_rhs",1);
  CCTKi_ScheduleGroupStorage("BSSN::BSSN_gupij",1);
  CCTKi_ScheduleGroupStorage("BSSN::have_global_bdry",1);
  CCTKi_ScheduleGroupStorage("fisheye::fisheye_vars",1);
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)Lapse_Startup,
                           "Lapse_Startup",
                           "lapse",
                           "lapse",
                           "Register banner",
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
    CCTKi_ScheduleFunction((void *)Lapse_CheckParameters,
                           "Lapse_CheckParameters",
                           "lapse",
                           "lapse",
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
    CCTKi_ScheduleFunction((void *)Lapse_InitSymBound,
                           "Lapse_InitSymBound",
                           "lapse",
                           "lapse",
                           "Schedule symmetries",
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
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)setup_initial_lapse_,
                           "lapse_initialdata",
                           "lapse",
                           "lapse",
                           "Initial data for lapse",
                           "CCTK_INITIAL",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           1,  /* Number of SYNC     groups   */
                           0,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           0,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "lapse::lapse_vars");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)lapse_postinitialdata_,
                           "lapsepostid",
                           "lapse",
                           "lapse",
                           "Compute post-initialdata quantities",
                           "ABE_PostInitial",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           0,  /* Number of SYNC     groups   */
                           2,  /* Number of Options           */
                           1,  /* Number of BEFORE  routines  */
                           1,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "GLOBAL",
                           "loop-local",
                           "empostid",
                           "shiftpostid");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)lapse_RegisterVars,
                           "lapse_RegisterVars",
                           "lapse",
                           "lapse",
                           "Register variables for MoL",
                           "MoL_Register",
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
                           "META");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)lapse_timestepping_,
                           "lapse_rhs",
                           "lapse",
                           "lapse",
                           "Evaluate RHS for lapse",
                           "MoL_CalcRHS",
                           "Fortran",
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
                           "bssn_rhs");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)lapse_update_bc_,
                           "lapse_update_bc",
                           "lapse",
                           "lapse",
                           "Update lapse bc's",
                           "ABE_PostStep",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           1,  /* Number of SYNC     groups   */
                           0,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           0,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "lapse::lapse_vars");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)lapse_postbc_,
                           "lapse_postbc",
                           "lapse",
                           "lapse",
                           "Update lapse derivatives",
                           "ABE_PostStep",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           1,  /* Number of SYNC     groups   */
                           0,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           1,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "lapse::lapse_derivatives",
                           "lapse_update_bc");
  }
}

/*@@
  @routine    CCTKi_BindingsParameterRecovery_lapse
  @author     Automatically generated by CreateScheduleBindings.pl
  @desc
              Creates the parameter recovery bindings for thorn lapse
  @enddesc
@@*/

int CCTKi_BindingsParameterRecovery_lapse(void);
int CCTKi_BindingsParameterRecovery_lapse(void)
{
  /* this thorn doesn't define any parameter recovery routines */
  return (0);
}