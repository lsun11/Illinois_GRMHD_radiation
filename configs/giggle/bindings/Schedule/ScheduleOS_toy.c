/*@@
   @file       ScheduleOS_toy.c
   @author     Automatically generated by CreateScheduleBindings.pl
   @desc
               Creates the schedule and parameter recovery bindings 
               for thorn OS_toy
   @enddesc
@@*/

#define THORN_IS_OS_toy

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_ScheduleBindings.h"

/* prototypes for schedule bindings functions to be registered */
/* Note that this is a cheat, we just need a function pointer. */
extern int CCTK_FNAME(mhd_OS_initialdata_local)(void);
extern int CCTK_FNAME(OS_rad_copy_to_prev_timelevel)(void);
extern int CCTK_FNAME(mhd_OS_initialdata_global)(void);
extern int CCTK_FNAME(particle_tracer_toy)(void);
extern int CCTK_FNAME(OS_rad_diag_integrals)(void);


void CCTKi_BindingsSchedule_OS_toy(void);
void CCTKi_BindingsSchedule_OS_toy(void)
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
  CCTKi_ScheduleGroupStorage("mhd_evolve::em_Blagrangemultiplier",3);
  CCTKi_ScheduleGroupStorage("mhd_evolve::em_Blagrangemultiplier_rhs",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::disk_atmosphere",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::eos_params1",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::eos_params2",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::mhdscalar",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::rad_conservatives",3);
  CCTKi_ScheduleGroupStorage("mhd_evolve::rad_conservatives_rhs",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::rad_primitives",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::radscalar",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::rad_pressure",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::mhd_nosync",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::micphys_conservatives",3);
  CCTKi_ScheduleGroupStorage("mhd_evolve::micphys_conservatives_rhs",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::microphys_primitives",1);
  CCTKi_ScheduleGroupStorage("fisheye::fisheye_vars",1);
  CCTKi_ScheduleGroupStorage("excision::excision_int_gfs",1);
  CCTKi_ScheduleGroupStorage("gw_extraction::gw_moment_arrays",1);
  CCTKi_ScheduleGroupStorage("diagnostics_vacuum::surf_params",1);
  CCTKi_ScheduleGroupStorage("diagnostics_vacuum::volIntegrals",1);
  CCTKi_ScheduleGroupStorage("diagnostics_mhd::volIntegrals_mhd",1);
  CCTKi_ScheduleGroupStorage("OS_toy::mhd_OS_private",1);
  CCTKi_ScheduleGroupStorage("OS_toy::mhd_OS_VolInt",1);
  CCTKi_ScheduleGroupStorage("OS_toy::particle_tracer_coord",1);
  CCTKi_ScheduleGroupStorage("OS_toy::pcle_stuff",1);
  CCTKi_ScheduleGroupStorage("OS_toy::more_pcle_stuff",1);
  CCTKi_ScheduleGroupStorage("OS_toy::v_previous",1);
  CCTKi_ScheduleGroupStorage("OS_toy::OS_center_diagnostics",1);
  CCTKi_ScheduleGroupStorage("OS_toy::analytic",1);
  CCTKi_ScheduleGroupStorage("OS_toy::Tracer",1);
  CCTKi_ScheduleGroupStorage("mhd_evolve::OS_stellar_surface",1);
  CCTKi_ScheduleGroupStorage("movingbox::volIntegrals_movingbox",1);
if (enable_OS_collapse==1){
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(mhd_OS_initialdata_local),
                           "OS_initialdata",
                           "OS_toy",
                           "OS_toy",
                           "Set up initial data for mhd_shock star evolver",
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
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(OS_rad_copy_to_prev_timelevel),
                           "OS_second_initialdata",
                           "OS_toy",
                           "OS_toy",
                           "Copy to previous timelevel",
                           "CCTK_INITIAL",
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
                           "OS_initialdata");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(mhd_OS_initialdata_global),
                           "OS_third_initialdata",
                           "OS_toy",
                           "OS_toy",
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
                           "OS_second_initialdata");
  }
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(particle_tracer_toy),
                           "tracer_toy",
                           "OS_toy",
                           "OS_toy",
                           "do particle tracer",
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
  {
    int cctkschedulei_tlevelarray[] = {0};
    CCTKi_ScheduleFunction((void *)CCTK_FNAME(OS_rad_diag_integrals),
                           "OS_poststep",
                           "OS_toy",
                           "OS_toy",
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
                           "tracer_toy");
  }
}
}

/*@@
  @routine    CCTKi_BindingsParameterRecovery_OS_toy
  @author     Automatically generated by CreateScheduleBindings.pl
  @desc
              Creates the parameter recovery bindings for thorn OS_toy
  @enddesc
@@*/

int CCTKi_BindingsParameterRecovery_OS_toy(void);
int CCTKi_BindingsParameterRecovery_OS_toy(void)
{
  /* this thorn doesn't define any parameter recovery routines */
  return (0);
}

