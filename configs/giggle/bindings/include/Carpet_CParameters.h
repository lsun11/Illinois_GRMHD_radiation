/*@@
   @header  Carpet_CParameters.h
   @author  Automatically generated by CreateParameterBindings.pl
   @desc
            Declares parameters of thorn Carpet
   @enddesc
 @@*/


#ifndef _CARPET_PARAMETERS_H_
#define _CARPET_PARAMETERS_H_ 1

#include "ParameterCGlobal.h"
#include "ParameterCRestrictedDRIVER.h"
#include "ParameterCPrivateCARPET.h"
#include "ParameterCRestrictedCACTUS.h"
#include "ParameterCRestrictedINITBASE.h"
#include "ParameterCRestrictedIO.h"

#define DECLARE_CCTK_PARAMETERS \
  DECLARE_GLOBAL_PARAMETER_STRUCT_PARAMS \
  DECLARE_RESTRICTED_DRIVER_STRUCT_PARAMS \
  DECLARE_PRIVATE_CARPET_STRUCT_PARAMS \
  CCTK_INT const terminate_next = RESTRICTED_CACTUS_STRUCT.terminate_next; \
  const char * const terminate = RESTRICTED_CACTUS_STRUCT.terminate; \
  CCTK_INT const cctk_itlast = RESTRICTED_CACTUS_STRUCT.cctk_itlast; \
  CCTK_REAL const cctk_initial_time = RESTRICTED_CACTUS_STRUCT.cctk_initial_time; \
  CCTK_REAL const cctk_final_time = RESTRICTED_CACTUS_STRUCT.cctk_final_time; \
  CCTK_REAL const max_runtime = RESTRICTED_CACTUS_STRUCT.max_runtime; \
  const char * const initial_data_setup_method = RESTRICTED_INITBASE_STRUCT.initial_data_setup_method; \
  const char * const out_dir = RESTRICTED_IO_STRUCT.out_dir; \
  enum { \
      dummy_CACTUS_terminate_next = sizeof( terminate_next ) \
    , dummy_CACTUS_terminate = sizeof( terminate ) \
    , dummy_CACTUS_cctk_itlast = sizeof( cctk_itlast ) \
    , dummy_CACTUS_cctk_initial_time = sizeof( cctk_initial_time ) \
    , dummy_CACTUS_cctk_final_time = sizeof( cctk_final_time ) \
    , dummy_CACTUS_max_runtime = sizeof( max_runtime ) \
    , dummy_INITBASE_initial_data_setup_method = sizeof( initial_data_setup_method ) \
    , dummy_IO_out_dir = sizeof( out_dir ) \
  };

#endif  /* _CARPET_PARAMETERS_H_ */
