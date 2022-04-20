/*@@
   @file       ScheduleOS_collapse.c
   @author     Automatically generated by CreateScheduleBindings.pl
   @desc
               Creates the schedule and parameter recovery bindings 
               for thorn OS_collapse
   @enddesc
@@*/

#define THORN_IS_OS_collapse

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_ScheduleBindings.h"

/* prototypes for schedule bindings functions to be registered */
/* Note that this is a cheat, we just need a function pointer. */
extern int CCTK_FNAME(read_OS_metric_inputfile_driver)(void);
extern int CCTK_FNAME(OS_metric_initialdata)(void);
extern int CCTK_FNAME(OS_initialdata_driver)(void);
extern int CCTK_FNAME(OS_initialdata_local)(void);
extern int CCTK_FNAME(OS_initialdata_global)(void);
extern int CCTK_FNAME(OS_initialdata_local2)(void);
extern int CCTK_FNAME(OS_initialdata_global2)(void);
extern int CCTK_FNAME(driver_h17)(void);
extern int CCTK_FNAME(OS_diagnostics)(void);


void CCTKi_BindingsSchedule_OS_collapse(void);
void CCTKi_BindingsSchedule_OS_collapse(void)
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
  CCTKi_ScheduleGroupStorage("shift::shift_vars",3);
  CCTKi_ScheduleGroupStorage("mhd_evolve::mhd_conservatives",3);
  CCTKi_ScheduleGroupStorage("mhd_evolve::mhd_rhs",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::mhd_primitives",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::mhd_vs",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::disk_atmosphere",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::eos_params1",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::eos_params2",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::mhdscalar",1);
  CCTKi_ScheduleGroupStorage("gw_extraction::gw_moment_arrays",1);
  CCTKi_ScheduleGroupStorage("fisheye::fisheye_vars",1);
  CCTKi_ScheduleGroupStorage("diagnostics_vacuum::surf_params",1);
  CCTKi_ScheduleGroupStorage("diagnostics_vacuum::volIntegrals",1);
  CCTKi_ScheduleGroupStorage("diagnostics_mhd::volIntegrals_mhd",1);
  CCTKi_ScheduleGroupStorage("excision::excision_int_gfs",1);
  CCTKi_ScheduleGroupStorage("OS_collapse::particle_tracer",1);
  CCTKi_ScheduleGroupStorage("OS_collapse::pcle_stuff",1);
  CCTKi_ScheduleGroupStorage("OS_collapse::more_pcle_stuff",1);
  CCTKi_ScheduleGroupStorage("OS_collapse::v_previous",1);
  CCTKi_ScheduleGroupStorage("OS_collapse::OS_center_diagnostics",1);
  CCTKi_ScheduleGroupStorage("OS_collapse::h17",1);
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(read_OS_metric_inputfile_driver),
                           "metric_readinfiles",
                           "OS_collapse",
                           "OS_collapse",
                           "Read initial Psi and lapse from input file",
                           "CCTK_INITIAL",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           0,  /* Number of SYNC     groups   */
                           1,  /* Number of Options           */
                           1,  /* Number of BEFORE  routines  */
                           0,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "LOCAL",
                           "lapse_initialdata");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(OS_metric_initialdata),
                           "OS_metric_initialdata",
                           "OS_collapse",
                           "OS_collapse",
                           "Set up TwoPunctures initial data on individual, local grids - part 1",
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
                           "metric_readinfiles");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(OS_initialdata_driver),
                           "bondi_id",
                           "OS_collapse",
                           "OS_collapse",
                           "Set up the bondi solution as initial data",
                           "CCTK_INITIAL",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           0,  /* Number of SYNC     groups   */
                           1,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           1,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "LOCAL",
                           "id_part2");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(OS_initialdata_local),
                           "id_local",
                           "OS_collapse",
                           "OS_collapse",
                           "Fill in metric quantities",
                           "CCTK_INITIAL",
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
                           "BSSN::BSSN_vars",
                           "bondi_id");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleGroup("OS_postid",
"OS_postid",
                        "OS_collapse",
                        "OS_collapse",
                        "Finish up BHNS initial data.  Need to schedule several GLOBAL function calls that don't work in CCTK_INITIAL. :(",
                        "CCTK_POSTPOSTINITIAL",
                        0,  /* Number of STORAGE  groups   */
                        0,  /* Number of COMM     groups   */
                        0,  /* Number of TRIGGERS groups   */
                        0,  /* Number of SYNC     groups   */
                        0,  /* Number of Options           */
                        2,  /* Number of BEFORE  routines  */
                        0,  /* Number of AFTER   routines  */
                        0,  /* Number of WHILE   variables */
                        0,  /* Number of IF   variables */
                        cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                        "ABE_PostInitial",
                        "MoL_PostStep");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(OS_initialdata_global),
                           "third_initialdata",
                           "OS_collapse",
                           "OS_collapse",
                           "Set rho_b_atm, etc.",
                           "OS_postid",
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
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(OS_initialdata_local2),
                           "fourth_initialdata",
                           "OS_collapse",
                           "OS_collapse",
                           "Fill in matter quantities.  Note that rho_b depends on rho_b_atm, and most hydro vars depend on rho_b",
                           "OS_postid",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           1,  /* Number of SYNC     groups   */
                           2,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           1,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "BSSN::BSSN_vars",
                           "GLOBAL",
                           "LOOP-LOCAL",
                           "third_initialdata");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(OS_initialdata_global2),
                           "fifth_initialdata",
                           "OS_collapse",
                           "OS_collapse",
                           "Set tau_atm, which depends on tau, which depends on rho_b, etc.",
                           "OS_postid",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           0,  /* Number of SYNC     groups   */
                           1,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           1,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "GLOBAL",
                           "fourth_initialdata");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(driver_h17),
                           "poststep",
                           "OS_collapse",
                           "OS_collapse",
                           "RHS of Equation H.17 in 'numerical relativity', by baumgarte and shapiro",
                           "CCTK_ANALYSIS",
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
                           "OS_collapse::h17");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(OS_diagnostics),
                           "OS_diagnostics",
                           "OS_collapse",
                           "OS_collapse",
                           "Evaluate diagnostic integrals",
                           "CCTK_ANALYSIS",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           0,  /* Number of SYNC     groups   */
                           1,  /* Number of Options           */
                           0,  /* Number of BEFORE  routines  */
                           1,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "GLOBAL",
                           "poststep");
  }
}

/*@@
  @routine    CCTKi_BindingsParameterRecovery_OS_collapse
  @author     Automatically generated by CreateScheduleBindings.pl
  @desc
              Creates the parameter recovery bindings for thorn OS_collapse
  @enddesc
@@*/

int CCTKi_BindingsParameterRecovery_OS_collapse(void);
int CCTKi_BindingsParameterRecovery_OS_collapse(void)
{
  /* this thorn doesn't define any parameter recovery routines */
  return (0);
}
