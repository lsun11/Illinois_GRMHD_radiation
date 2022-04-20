/*@@
   @header  bssn_CParameters.h
   @author  Automatically generated by CreateParameterBindings.pl
   @desc
            Declares parameters of thorn bssn
   @enddesc
 @@*/


#ifndef _BSSN_PARAMETERS_H_
#define _BSSN_PARAMETERS_H_ 1

#include "ParameterCGlobal.h"
#include "ParameterCRestrictedBSSN.h"
#include "ParameterCPrivateBSSN.h"
#include "ParameterCRestrictedDIAGNOSTICS_MHD.h"
#include "ParameterCRestrictedDIAGNOSTICS_VACUUM.h"
#include "ParameterCRestrictedEXCISION.h"
#include "ParameterCRestrictedGRID.h"
#include "ParameterCRestrictedIO.h"
#include "ParameterCRestrictedLAPSE.h"
#include "ParameterCRestrictedMETHODOFLINES.h"
#include "ParameterCRestrictedSHIFT.h"

#define DECLARE_CCTK_PARAMETERS \
  DECLARE_GLOBAL_PARAMETER_STRUCT_PARAMS \
  DECLARE_RESTRICTED_BSSN_STRUCT_PARAMS \
  DECLARE_PRIVATE_BSSN_STRUCT_PARAMS \
  CCTK_INT const Nphi_points = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.Nphi_points; \
  CCTK_INT const ntot = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.ntot; \
  CCTK_INT const N_phi = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.N_phi; \
  CCTK_INT const N_theta = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.N_theta; \
  CCTK_INT const sym_factor = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.sym_factor; \
  CCTK_REAL const dcostheta = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.dcostheta; \
  CCTK_REAL const dphi = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.dphi; \
  CCTK_REAL const surf_radius = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.surf_radius; \
  CCTK_REAL const rbr = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.rbr; \
  CCTK_REAL const drbdr = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.drbdr; \
  CCTK_REAL const ddrbddr = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.ddrbddr; \
  CCTK_REAL const C_ko = RESTRICTED_EXCISION_STRUCT.C_ko; \
  const char * const domain = RESTRICTED_GRID_STRUCT.domain; \
  const char * const bitant_plane = RESTRICTED_GRID_STRUCT.bitant_plane; \
  CCTK_INT const out_every = RESTRICTED_IO_STRUCT.out_every; \
  const char * const slicing_type = RESTRICTED_LAPSE_STRUCT.slicing_type; \
  CCTK_INT const MoL_Num_Evolved_Vars = RESTRICTED_METHODOFLINES_STRUCT.MoL_Num_Evolved_Vars; \
  CCTK_INT const MoL_Num_Constrained_Vars = RESTRICTED_METHODOFLINES_STRUCT.MoL_Num_Constrained_Vars; \
  CCTK_INT const MoL_Num_SaveAndRestore_Vars = RESTRICTED_METHODOFLINES_STRUCT.MoL_Num_SaveAndRestore_Vars; \
  CCTK_INT const Spatial_Gauge = RESTRICTED_SHIFT_STRUCT.Spatial_Gauge; \
  enum { \
      dummy_DIAGNOSTICS_MHD_Nphi_points = sizeof( Nphi_points ) \
    , dummy_DIAGNOSTICS_VACUUM_ntot = sizeof( ntot ) \
    , dummy_DIAGNOSTICS_VACUUM_N_phi = sizeof( N_phi ) \
    , dummy_DIAGNOSTICS_VACUUM_N_theta = sizeof( N_theta ) \
    , dummy_DIAGNOSTICS_VACUUM_sym_factor = sizeof( sym_factor ) \
    , dummy_DIAGNOSTICS_VACUUM_dcostheta = sizeof( dcostheta ) \
    , dummy_DIAGNOSTICS_VACUUM_dphi = sizeof( dphi ) \
    , dummy_DIAGNOSTICS_VACUUM_surf_radius = sizeof( surf_radius ) \
    , dummy_DIAGNOSTICS_VACUUM_rbr = sizeof( rbr ) \
    , dummy_DIAGNOSTICS_VACUUM_drbdr = sizeof( drbdr ) \
    , dummy_DIAGNOSTICS_VACUUM_ddrbddr = sizeof( ddrbddr ) \
    , dummy_EXCISION_C_ko = sizeof( C_ko ) \
    , dummy_GRID_domain = sizeof( domain ) \
    , dummy_GRID_bitant_plane = sizeof( bitant_plane ) \
    , dummy_IO_out_every = sizeof( out_every ) \
    , dummy_LAPSE_slicing_type = sizeof( slicing_type ) \
    , dummy_METHODOFLINES_MoL_Num_Evolved_Vars = sizeof( MoL_Num_Evolved_Vars ) \
    , dummy_METHODOFLINES_MoL_Num_Constrained_Vars = sizeof( MoL_Num_Constrained_Vars ) \
    , dummy_METHODOFLINES_MoL_Num_SaveAndRestore_Vars = sizeof( MoL_Num_SaveAndRestore_Vars ) \
    , dummy_SHIFT_Spatial_Gauge = sizeof( Spatial_Gauge ) \
  };

#endif  /* _BSSN_PARAMETERS_H_ */

