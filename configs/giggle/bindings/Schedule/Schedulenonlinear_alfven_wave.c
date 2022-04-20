/*@@
   @file       Schedulenonlinear_alfven_wave.c
   @author     Automatically generated by CreateScheduleBindings.pl
   @desc
               Creates the schedule and parameter recovery bindings 
               for thorn nonlinear_alfven_wave
   @enddesc
@@*/

#define THORN_IS_nonlinear_alfven_wave

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_ScheduleBindings.h"

/* prototypes for schedule bindings functions to be registered */
/* Note that this is a cheat, we just need a function pointer. */
extern int CCTK_FNAME(nonlinear_alfven_wave_initialdata_local)(void);
extern int CCTK_FNAME(nonlinear_alfven_wave_initialdata_global)(void);
extern int CCTK_FNAME(nonlinear_alfven_wave_analytic_time)(void);


void CCTKi_BindingsSchedule_nonlinear_alfven_wave(void);
void CCTKi_BindingsSchedule_nonlinear_alfven_wave(void)
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
  CCTKi_ScheduleGroupStorage("mhd_evolve::em_conservativex",3);
  CCTKi_ScheduleGroupStorage("mhd_evolve::em_conservativey",3);
  CCTKi_ScheduleGroupStorage("mhd_evolve::em_conservativez",3);
  CCTKi_ScheduleGroupStorage("mhd_evolve::em_Ax",3);
  CCTKi_ScheduleGroupStorage("mhd_evolve::em_Ay",3);
  CCTKi_ScheduleGroupStorage("mhd_evolve::em_Az",3);
  CCTKi_ScheduleGroupStorage("mhd_evolve::em_rhsx",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::em_rhsy",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::em_rhsz",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::disk_atmosphere",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::eos_params1",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::eos_params2",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::mhdscalar",1);
  CCTKi_ScheduleGroupStorage("fisheye::fisheye_vars",1);
  CCTKi_ScheduleGroupStorage("excision::excision_int_gfs",1);
  CCTKi_ScheduleGroupStorage("gw_extraction::gw_moment_arrays",1);
  CCTKi_ScheduleGroupStorage("diagnostics_vacuum::volIntegrals",1);
  CCTKi_ScheduleGroupStorage("diagnostics_mhd::volIntegrals_mhd",1);
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(nonlinear_alfven_wave_initialdata_local),
                           "initialdata",
                           "nonlinear_alfven_wave",
                           "nonlinear_alfven_wave",
                           "Set up initial data for nonlinear_alfven_wave star evolver",
                           "CCTK_INITIAL",
                           "Fortran",
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
                           "lapse_initialdata");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(nonlinear_alfven_wave_initialdata_global),
                           "third_initialdata",
                           "nonlinear_alfven_wave",
                           "nonlinear_alfven_wave",
                           "Set rho_b_atm, etc.",
                           "ABE_PostInitial",
                           "Fortran",
                           0,  /* Number of STORAGE  groups   */
                           0,  /* Number of COMM     groups   */
                           0,  /* Number of TRIGGERS groups   */
                           0,  /* Number of SYNC     groups   */
                           1,  /* Number of Options           */
                           1,  /* Number of BEFORE  routines  */
                           1,  /* Number of AFTER   routines  */
                           0,  /* Number of WHILE   variables */
                           0,  /* Number of IF   variables */
                           cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */,
                           "GLOBAL",
                           "postid",
                           "second_initialdata");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(nonlinear_alfven_wave_analytic_time),
                           "compute_alfven_wave_analytic",
                           "nonlinear_alfven_wave",
                           "nonlinear_alfven_wave",
                           "Compute analytic solution of the nonlinear Alfven wave",
                           "CCTK_ANALYSIS",
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
                           "sanitycheck_restore_Aij");
  }
}

/*@@
  @routine    CCTKi_BindingsParameterRecovery_nonlinear_alfven_wave
  @author     Automatically generated by CreateScheduleBindings.pl
  @desc
              Creates the parameter recovery bindings for thorn nonlinear_alfven_wave
  @enddesc
@@*/

int CCTKi_BindingsParameterRecovery_nonlinear_alfven_wave(void);
int CCTKi_BindingsParameterRecovery_nonlinear_alfven_wave(void)
{
  /* this thorn doesn't define any parameter recovery routines */
  return (0);
}

