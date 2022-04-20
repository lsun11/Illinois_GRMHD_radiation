#include "cctk.h"
#include "cctk_Arguments.h"

static CCTK_REAL mol_dt;


void MoL_OldBdry_SetDt (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS

  CCTK_REAL dt;
  CCTK_REAL old_time;

  old_time = *Original_Time - *Original_Delta_Time / cctkGH->cctk_timefac;

  mol_dt = cctkGH->cctk_delta_time;

  dt = cctkGH->cctk_time - old_time;

  cctkGH->cctk_delta_time = dt * cctkGH->cctk_timefac;
}

void MoL_OldBdry_ResetDt (CCTK_ARGUMENTS)
{
  cctkGH->cctk_delta_time = mol_dt;
}
