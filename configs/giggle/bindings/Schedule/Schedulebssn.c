/*@@
   @file       Schedulebssn.c
   @author     Automatically generated by CreateScheduleBindings.pl
   @desc
               Creates the schedule and parameter recovery bindings 
               for thorn bssn
   @enddesc
@@*/

#define THORN_IS_bssn

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_ScheduleBindings.h"

/* prototypes for schedule bindings functions to be registered */
/* Note that this is a cheat, we just need a function pointer. */
extern int BSSN_Startup(void);
extern int BSSN_CheckParameters(void);
extern int BSSN_RegisterVars(void);
extern int BSSN_InitSymBound(void);
extern int CCTK_FNAME(bssn_set_have_global_bdry_minmax)(void);
extern int CCTK_FNAME(bssn_set_have_global_bdry_minmax)(void);
extern int CCTK_FNAME(BSSN_PostInitialData)(void);
extern int CCTK_FNAME(print_time)(void);
extern int CCTK_FNAME(BSSN_timestepping_Cowling)(void);
extern int CCTK_FNAME(BSSN_timestepping)(void);
extern int CCTK_FNAME(driver_bssn_update_boundary)(void);
extern int CCTK_FNAME(driver_ricci_constraints)(void);
extern int CCTK_FNAME(driver_bssn_post_regrid)(void);
extern int CCTK_FNAME(driver_bssn_post_regrid)(void);
extern int CCTK_FNAME(driver_bssn_post_regrid)(void);
extern int CCTK_FNAME(driver_bssn_post_regrid)(void);
extern int CCTK_FNAME(set_metric_rotation)(void);
extern int CCTK_FNAME(cowling_poststep)(void);
extern int CCTK_FNAME(driver_post_update_boundary)(void);
extern int CCTK_FNAME(BSSN_Gauge_Derivs)(void);
extern int CCTK_FNAME(increment_iter_count)(void);


void CCTKi_BindingsSchedule_bssn(void);
void CCTKi_BindingsSchedule_bssn(void)
{
  DECLARE_CCTK_PARAMETERS
  CCTKi_ScheduleGroupStorage("BSSN::BSSN_vars",3);
  CCTKi_ScheduleGroupStorage("BSSN::BSSN_rhs",1);
  CCTKi_ScheduleGroupStorage("BSSN::have_global_bdry",1);
  CCTKi_ScheduleGroupStorage("BSSN::BSSN_diag_restrict",1);
  CCTKi_ScheduleGroupStorage("BSSN::BSSN_gupij",1);
  CCTKi_ScheduleGroupStorage("BSSN::BSSN_matter",1);
  CCTKi_ScheduleGroupStorage("BSSN::BSSN_aux_restrict2",1);
  CCTKi_ScheduleGroupStorage("BSSN::phi_derivs",1);
  CCTKi_ScheduleGroupStorage("BSSN::BSSN_aux_private",1);
  CCTKi_ScheduleGroupStorage("bssn::BSSN_adm_mass",1);
  CCTKi_ScheduleGroupStorage("BSSN::BSSN_AH",1);
  CCTKi_ScheduleGroupStorage("BSSN::BSSN_refbd",3);
  CCTKi_ScheduleGroupStorage("lapse::lapse_vars",3);
  CCTKi_ScheduleGroupStorage("lapse::lapse_derivatives",1);
  CCTKi_ScheduleGroupStorage("lapse::lapse_vars_aux",1);
  CCTKi_ScheduleGroupStorage("shift::shift_vars",3);
  CCTKi_ScheduleGroupStorage("shift::shift_vars_temp",1);
  CCTKi_ScheduleGroupStorage("fisheye::fisheye_vars",1);
  CCTKi_ScheduleGroupStorage("excision::excision_int_gfs",1);
  CCTKi_ScheduleGroupStorage("BSSN::metric_spher_pol_1",1);
  CCTKi_ScheduleGroupStorage("BSSN::metric_spher_pol_2",1);
  CCTKi_ScheduleGroupStorage("BSSN::metric_spher_pol_3",1);
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)BSSN_Startup,
                           "BSSN_Startup",
                           "bssn",
                           "BSSN",
                           "Print a banner onscreen indicating that bssn thorn is active",
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
    CCTKi_ScheduleFunction((void *)BSSN_CheckParameters,
                           "BSSN_CheckParameters",
                           "bssn",
                           "BSSN",
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
    CCTKi_ScheduleFunction((void *)BSSN_RegisterVars,
                           "BSSN_RegisterVars",
                           "bssn",
                           "BSSN",
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
    CCTKi_ScheduleFunction((void *)BSSN_InitSymBound,
                           "BSSN_InitSymBound",
                           "bssn",
                           "BSSN",
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
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(bssn_set_have_global_bdry_minmax),
                           "bssn_set_have_global_bdry_minmax",
                           "bssn",
                           "BSSN",
                           "Set have_global_bdry_max/min (used for update_boundary)",
                           "CCTK_POSTINITIAL",
                           "Fortran",
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
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(bssn_set_have_global_bdry_minmax),
                           "bssn_set_have_global_bdry_minmax",
                           "bssn",
                           "BSSN",
                           "Set have_global_bdry_max/min (used for update_boundary) when restarting from checkpoint, since have_global_bdry_max/min is different on each processor, and checkpoint will only save proc. zero!",
                           "CCTK_POST_RECOVER_VARIABLES",
                           "Fortran",
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
    CCTKi_ScheduleGroup("ABE_PostInitial",
"ABE_PostInitial",
                        "bssn",
                        "BSSN",
                        "ABE Post-initial data routines",
                        "CCTK_POSTPOSTINITIAL",
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
                        "MoL_PostStep");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(BSSN_PostInitialData),
                           "postid",
                           "bssn",
                           "BSSN",
                           "Compute post-initialdata quantities",
                           "ABE_PostInitial",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           0,  /* Number of SYNC     groups   */
                           2,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           1,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "GLOBAL",
                           "loop-local",
                           "shift_initialdata");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(print_time),
                           "print_time",
                           "bssn",
                           "BSSN",
                           "Print the iteration number and time.",
                           "MoL_PreStep",
                           "Fortran",
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
if(cowling_enable == 1 || bssn_enable == 0) {
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(BSSN_timestepping_Cowling),
                           "bssn_rhs",
                           "bssn",
                           "BSSN",
                           "Set BSSN rhs's to zero for Cowling",
                           "MoL_CalcRHS",
                           "Fortran",
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
else {
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(BSSN_timestepping),
                           "bssn_rhs",
                           "bssn",
                           "BSSN",
                           "Set BSSN rhs's",
                           "MoL_CalcRHS",
                           "Fortran",
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
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleGroup("ABE_PostStep",
"ABE_PostStep",
                        "bssn",
                        "BSSN",
                        "ABE post-step (i.e., after gridfunction update)",
                        "MoL_Step",
                        0,  /* Number of STORAGE  groups   */
                        0,  /* Number of COMM     groups   */
                        0,  /* Number of TRIGGERS groups   */
                        0,  /* Number of SYNC     groups   */
                        0,  /* Number of Options           */
                        1,  /* Number of BEFORE  routines  */
                        1,  /* Number of AFTER   routines  */
                        0,  /* Number of WHILE   variables */
                        0,  /* Number of IF   variables */
                        cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                        "MoL_ResetDeltaTime",
                        "MoL_PostStep");
  }
if(cowling_enable == 0) {
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(driver_bssn_update_boundary),
                           "bssn_update_bc",
                           "bssn",
                           "BSSN",
                           "Update outer boundary.",
                           "ABE_PostStep",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           2,  /* Number of SYNC     groups   */
                           0,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           0,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "BSSN::BSSN_vars",
                           "BSSN::BSSN_refbd");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(driver_ricci_constraints),
                           "bssn_ricci_const",
                           "bssn",
                           "BSSN",
                           "Compute Ricci, constraints.",
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
                           "BSSN::phi_derivs",
                           "bssn_update_bc");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(driver_bssn_post_regrid),
                           "bssn_postregrid",
                           "bssn",
                           "BSSN",
                           "postregridinitial: Set auxiliary BSSN quantities over all grids after Carpet moves any grid.",
                           "CCTK_POSTREGRIDINITIAL",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           9,  /* Number of SYNC     groups   */
                           2,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           1,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "BSSN::phi_derivs",
                           "BSSN::BSSN_diag_restrict",
                           "BSSN::BSSN_gupij",
                           "BSSN::BSSN_vars",
                           "BSSN::BSSN_matter",
                           "BSSN::BSSN_aux_restrict2",
                           "BSSN::BSSN_aux_private",
                           "BSSN::BSSN_AH",
                           "BSSN::BSSN_refbd",
                           "LEVEL",
                           "LOOP-LOCAL",
                           "fish_postregrid_update0");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(driver_bssn_post_regrid),
                           "bssn_postregrid",
                           "bssn",
                           "BSSN",
                           "postrestrictinitial:  Set auxiliary BSSN quantities over all grids after Carpet moves any grid.",
                           "CCTK_POSTRESTRICTINITIAL",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           9,  /* Number of SYNC     groups   */
                           2,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           1,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "BSSN::phi_derivs",
                           "BSSN::BSSN_diag_restrict",
                           "BSSN::BSSN_gupij",
                           "BSSN::BSSN_vars",
                           "BSSN::BSSN_matter",
                           "BSSN::BSSN_aux_restrict2",
                           "BSSN::BSSN_aux_private",
                           "BSSN::BSSN_AH",
                           "BSSN::BSSN_refbd",
                           "LEVEL",
                           "LOOP-LOCAL",
                           "fish_postregrid_update0");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(driver_bssn_post_regrid),
                           "bssn_postregrid",
                           "bssn",
                           "BSSN",
                           "Set auxiliary BSSN quantities over all grids after Carpet moves any grid.",
                           "CCTK_POSTRESTRICT",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           9,  /* Number of SYNC     groups   */
                           2,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           1,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "BSSN::phi_derivs",
                           "BSSN::BSSN_diag_restrict",
                           "BSSN::BSSN_gupij",
                           "BSSN::BSSN_vars",
                           "BSSN::BSSN_matter",
                           "BSSN::BSSN_aux_restrict2",
                           "BSSN::BSSN_aux_private",
                           "BSSN::BSSN_AH",
                           "BSSN::BSSN_refbd",
                           "LEVEL",
                           "LOOP-LOCAL",
                           "fish_postregrid_update0");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(driver_bssn_post_regrid),
                           "bssn_postregrid",
                           "bssn",
                           "BSSN",
                           "Set auxiliary BSSN quantities over all grids after Carpet moves any grid.",
                           "CCTK_POSTREGRID",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           9,  /* Number of SYNC     groups   */
                           2,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           1,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "BSSN::phi_derivs",
                           "BSSN::BSSN_diag_restrict",
                           "BSSN::BSSN_gupij",
                           "BSSN::BSSN_vars",
                           "BSSN::BSSN_matter",
                           "BSSN::BSSN_aux_restrict2",
                           "BSSN::BSSN_aux_private",
                           "BSSN::BSSN_AH",
                           "BSSN::BSSN_refbd",
                           "LEVEL",
                           "LOOP-LOCAL",
                           "fish_postregrid_update");
  }
} else if (rot_metric==1) {
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(set_metric_rotation),
                           "metric_rotation",
                           "bssn",
                           "BSSN",
                           "rotate the metric analytically",
                           "ABE_PostStep",
                           "Fortran",
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
} else if (cowling_enable==1) {
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(cowling_poststep),
                           "cowling_poststep",
                           "bssn",
                           "BSSN",
                           "May need this to do prolongation for BSSN_refbd",
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
                           "BSSN::BSSN_refbd");
  }
}
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(driver_post_update_boundary),
                           "bssn_post_bc",
                           "bssn",
                           "BSSN",
                           "Post bc update: If bc_type==7, compute ADM mass.",
                           "ABE_PostStep",
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
                           "bssn_ricci_const");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(BSSN_Gauge_Derivs),
                           "gaugederivs",
                           "bssn",
                           "BSSN",
                           "Compute time derivative of lapse and shift (might change it later)",
                           "ABE_PostStep",
                           "Fortran",
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
                           "shift_update_bc",
                           "lapse_update_bc");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(increment_iter_count),
                           "iter_count_update",
                           "bssn",
                           "BSSN",
                           "Increment iter_count, which gives the current MoL step",
                           "ABE_PostStep",
                           "Fortran",
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
                           "bssn_update_bc",
                           "gaugederivs");
  }
}

/*@@
  @routine    CCTKi_BindingsParameterRecovery_bssn
  @author     Automatically generated by CreateScheduleBindings.pl
  @desc
              Creates the parameter recovery bindings for thorn bssn
  @enddesc
@@*/

int CCTKi_BindingsParameterRecovery_bssn(void);
int CCTKi_BindingsParameterRecovery_bssn(void)
{
  /* this thorn doesn't define any parameter recovery routines */
  return (0);
}

