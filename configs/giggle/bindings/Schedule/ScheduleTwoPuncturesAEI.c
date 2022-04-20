/*@@
   @file       ScheduleTwoPuncturesAEI.c
   @author     Automatically generated by CreateScheduleBindings.pl
   @desc
               Creates the schedule and parameter recovery bindings 
               for thorn TwoPuncturesAEI
   @enddesc
@@*/

#define THORN_IS_TwoPuncturesAEI

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_ScheduleBindings.h"

/* prototypes for schedule bindings functions to be registered */
/* Note that this is a cheat, we just need a function pointer. */
extern int TwoPuncturesAEI(void);
extern int CCTK_FNAME(TwoPuncturesAEI_initialdata)(void);
extern int CCTK_FNAME(TwoPuncturesAEI_diagnostics)(void);


void CCTKi_BindingsSchedule_TwoPuncturesAEI(void);
void CCTKi_BindingsSchedule_TwoPuncturesAEI(void)
{
  DECLARE_CCTK_PARAMETERS
  CCTKi_ScheduleGroupStorage("BSSN::BSSN_vars",3);
  CCTKi_ScheduleGroupStorage("BSSN::BSSN_gupij",1);
  CCTKi_ScheduleGroupStorage("BSSN::BSSN_matter",1);
  CCTKi_ScheduleGroupStorage("BSSN::BSSN_AH",1);
  CCTKi_ScheduleGroupStorage("BSSN::BSSN_aux_restrict2",1);
  CCTKi_ScheduleGroupStorage("BSSN::phi_derivs",1);
  CCTKi_ScheduleGroupStorage("BSSN::BSSN_diag_restrict",1);
  CCTKi_ScheduleGroupStorage("lapse::lapse_vars",3);
  CCTKi_ScheduleGroupStorage("lapse::lapse_derivatives",1);
  CCTKi_ScheduleGroupStorage("shift::shift_vars",3);
  CCTKi_ScheduleGroupStorage("gw_extraction::gw_moment_arrays",1);
  CCTKi_ScheduleGroupStorage("fisheye::fisheye_vars",1);
  CCTKi_ScheduleGroupStorage("diagnostics_vacuum::surf_params",1);
  CCTKi_ScheduleGroupStorage("diagnostics_vacuum::bh_posns",1);
  CCTKi_ScheduleGroupStorage("diagnostics_vacuum::volIntegrals",1);
  CCTKi_ScheduleGroupStorage("TwoPuncturesAEI::puncture_u",1);
  CCTKi_ScheduleGroupStorage("TwoPuncturesAEI::bare_mass",1);
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)TwoPuncturesAEI,
                           "TwoPuncturesAEI_readinfiles",
                           "TwoPuncturesAEI",
                           "TwoPuncturesAEI",
                           "Create puncture black hole initial data",
                           "CCTK_INITIAL",
                           "C",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           1,  /* Number of SYNC     groups   */
                           0,  /* Number of Options           */
                           1,  /* Number of BEFORE  routines  */
                           0,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "BSSN::BSSN_vars",
                           "lapse_initialdata");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(TwoPuncturesAEI_initialdata),
                           "TwoPuncturesAEI_initialdata",
                           "TwoPuncturesAEI",
                           "TwoPuncturesAEI",
                           "Set up TwoPuncturesAEI initial data on individual, local grids - part 1",
                           "CCTK_INITIAL",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           1,  /* Number of SYNC     groups   */
                           0,  /* Number of Options           */
                           1,  /* Number of BEFORE  routines  */
                           1,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "BSSN::BSSN_vars",
                           "lapse_initialdata",
                           "TwoPuncturesAEI_readinfiles");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(TwoPuncturesAEI_diagnostics),
                           "poststep",
                           "TwoPuncturesAEI",
                           "TwoPuncturesAEI",
                           "Evaluate diagnostic integrals",
                           "CCTK_ANALYSIS",
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
                           "GLOBAL");
  }
}

/*@@
  @routine    CCTKi_BindingsParameterRecovery_TwoPuncturesAEI
  @author     Automatically generated by CreateScheduleBindings.pl
  @desc
              Creates the parameter recovery bindings for thorn TwoPuncturesAEI
  @enddesc
@@*/

int CCTKi_BindingsParameterRecovery_TwoPuncturesAEI(void);
int CCTKi_BindingsParameterRecovery_TwoPuncturesAEI(void)
{
  /* this thorn doesn't define any parameter recovery routines */
  return (0);
}
