#include "cctk_Config.h"
#include "CParameterStructNames.h"

struct
{
  CCTK_REAL excision_radius;
  CCTK_REAL run_time;
  CCTK_INT Symmetry;
  CCTK_INT bssn_enable;
  CCTK_INT cowling_enable;
  CCTK_INT excision_enable;
  CCTK_INT fisheye_enable;
  CCTK_INT iter_count;
  CCTK_INT number_of_mol_ministeps;
  CCTK_INT rot_metric;
  CCTK_INT trA_detg_enforce;
} GLOBAL_PARAMETER_STRUCT;
