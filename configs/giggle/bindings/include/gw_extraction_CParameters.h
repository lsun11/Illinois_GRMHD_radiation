/*@@
   @header  gw_extraction_CParameters.h
   @author  Automatically generated by CreateParameterBindings.pl
   @desc
            Declares parameters of thorn gw_extraction
   @enddesc
 @@*/


#ifndef _GW_EXTRACTION_PARAMETERS_H_
#define _GW_EXTRACTION_PARAMETERS_H_ 1

#include "ParameterCGlobal.h"
#include "ParameterCRestrictedGW_EXTRACTION.h"
#include "ParameterCPrivateGW_EXTRACTION.h"
#include "ParameterCRestrictedDIAGNOSTICS_VACUUM.h"
#include "ParameterCRestrictedIO.h"

#define DECLARE_CCTK_PARAMETERS \
  DECLARE_GLOBAL_PARAMETER_STRUCT_PARAMS \
  DECLARE_RESTRICTED_GW_EXTRACTION_STRUCT_PARAMS \
  DECLARE_PRIVATE_GW_EXTRACTION_STRUCT_PARAMS \
  CCTK_REAL const surf_radius = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.surf_radius; \
  CCTK_INT const numtheta = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.numtheta; \
  CCTK_INT const numphi = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.numphi; \
  CCTK_REAL const rbr = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.rbr; \
  CCTK_INT const out_every = RESTRICTED_IO_STRUCT.out_every; \
  enum { \
      dummy_DIAGNOSTICS_VACUUM_surf_radius = sizeof( surf_radius ) \
    , dummy_DIAGNOSTICS_VACUUM_numtheta = sizeof( numtheta ) \
    , dummy_DIAGNOSTICS_VACUUM_numphi = sizeof( numphi ) \
    , dummy_DIAGNOSTICS_VACUUM_rbr = sizeof( rbr ) \
    , dummy_IO_out_every = sizeof( out_every ) \
  };

#endif  /* _GW_EXTRACTION_PARAMETERS_H_ */

