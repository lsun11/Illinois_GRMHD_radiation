/*@@
   @header  CarpetLib_CParameters.h
   @author  Automatically generated by CreateParameterBindings.pl
   @desc
            Declares parameters of thorn CarpetLib
   @enddesc
 @@*/


#ifndef _CARPETLIB_PARAMETERS_H_
#define _CARPETLIB_PARAMETERS_H_ 1

#include "ParameterCGlobal.h"
#include "ParameterCPrivateCARPETLIB.h"
#include "ParameterCRestrictedIO.h"

#define DECLARE_CCTK_PARAMETERS \
  DECLARE_GLOBAL_PARAMETER_STRUCT_PARAMS \
  DECLARE_PRIVATE_CARPETLIB_STRUCT_PARAMS \
  const char * const out_dir = RESTRICTED_IO_STRUCT.out_dir; \
  enum { \
      dummy_IO_out_dir = sizeof( out_dir ) \
  };

#endif  /* _CARPETLIB_PARAMETERS_H_ */

