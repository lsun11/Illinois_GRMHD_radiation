/*@@
   @header  CarpetInterp_CParameters.h
   @author  Automatically generated by CreateParameterBindings.pl
   @desc
            Declares parameters of thorn CarpetInterp
   @enddesc
 @@*/


#ifndef _CARPETINTERP_PARAMETERS_H_
#define _CARPETINTERP_PARAMETERS_H_ 1

#include "ParameterCGlobal.h"
#include "ParameterCPrivateCARPETINTERP.h"
#include "ParameterCRestrictedCACTUS.h"

#define DECLARE_CCTK_PARAMETERS \
  DECLARE_GLOBAL_PARAMETER_STRUCT_PARAMS \
  DECLARE_PRIVATE_CARPETINTERP_STRUCT_PARAMS \
  CCTK_REAL const cctk_initial_time = RESTRICTED_CACTUS_STRUCT.cctk_initial_time; \
  enum { \
      dummy_CACTUS_cctk_initial_time = sizeof( cctk_initial_time ) \
  };

#endif  /* _CARPETINTERP_PARAMETERS_H_ */

