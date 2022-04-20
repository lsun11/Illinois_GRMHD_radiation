/*@@
   @header  disk_powerlaw_CParameters.h
   @author  Automatically generated by CreateParameterBindings.pl
   @desc
            Declares parameters of thorn disk_powerlaw
   @enddesc
 @@*/


#ifndef _DISK_POWERLAW_PARAMETERS_H_
#define _DISK_POWERLAW_PARAMETERS_H_ 1

#include "ParameterCGlobal.h"
#include "ParameterCPrivateDISK_POWERLAW.h"
#include "ParameterCRestrictedBSSN.h"
#include "ParameterCRestrictedDIAGNOSTICS_MHD.h"
#include "ParameterCRestrictedDIAGNOSTICS_VACUUM.h"
#include "ParameterCRestrictedGRID.h"
#include "ParameterCRestrictedGW_EXTRACTION.h"
#include "ParameterCRestrictedIO.h"
#include "ParameterCRestrictedMHD_EVOLVE.h"

#define DECLARE_CCTK_PARAMETERS \
  DECLARE_GLOBAL_PARAMETER_STRUCT_PARAMS \
  DECLARE_PRIVATE_DISK_POWERLAW_STRUCT_PARAMS \
  CCTK_REAL const xbh1_initial = RESTRICTED_BSSN_STRUCT.xbh1_initial; \
  CCTK_REAL const xbh2_initial = RESTRICTED_BSSN_STRUCT.xbh2_initial; \
  CCTK_REAL const rah1_initial = RESTRICTED_BSSN_STRUCT.rah1_initial; \
  CCTK_REAL const rah2_initial = RESTRICTED_BSSN_STRUCT.rah2_initial; \
  CCTK_REAL const binary_orb_freq = RESTRICTED_BSSN_STRUCT.binary_orb_freq; \
  CCTK_INT const enable_lower_order_at_boundaries = RESTRICTED_BSSN_STRUCT.enable_lower_order_at_boundaries; \
  CCTK_INT const N_rad_ray = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.N_rad_ray; \
  CCTK_REAL const rhosurf_rmin = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.rhosurf_rmin; \
  CCTK_REAL const rhosurf_rmax = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.rhosurf_rmax; \
  CCTK_REAL const inner_lum_rad_ratio = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.inner_lum_rad_ratio; \
  CCTK_REAL const lum_outer_rad = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.lum_outer_rad; \
  CCTK_REAL const surf_radius = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.surf_radius; \
  CCTK_REAL const rbr = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.rbr; \
  CCTK_REAL const drbdr = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.drbdr; \
  CCTK_REAL const ddrbddr = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.ddrbddr; \
  CCTK_INT const scaledsurf = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.scaledsurf; \
  CCTK_INT const nsurf = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.nsurf; \
  CCTK_INT const arithsurf = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.arithsurf; \
  CCTK_REAL const rsurf2 = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.rsurf2; \
  CCTK_INT const N_theta = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.N_theta; \
  CCTK_INT const N_phi = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.N_phi; \
  CCTK_INT const sym_factor = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.sym_factor; \
  const char * const domain = RESTRICTED_GRID_STRUCT.domain; \
  const char * const bitant_plane = RESTRICTED_GRID_STRUCT.bitant_plane; \
  CCTK_REAL const radius_GW = RESTRICTED_GW_EXTRACTION_STRUCT.radius_GW; \
  CCTK_REAL const theta_GW = RESTRICTED_GW_EXTRACTION_STRUCT.theta_GW; \
  CCTK_REAL const phi_GW = RESTRICTED_GW_EXTRACTION_STRUCT.phi_GW; \
  CCTK_INT const out_every = RESTRICTED_IO_STRUCT.out_every; \
  CCTK_INT const primitives_solver = RESTRICTED_MHD_EVOLVE_STRUCT.primitives_solver; \
  CCTK_INT const enable_shocktest_primitive_mode = RESTRICTED_MHD_EVOLVE_STRUCT.enable_shocktest_primitive_mode; \
  CCTK_INT const use_HARM_primitives = RESTRICTED_MHD_EVOLVE_STRUCT.use_HARM_primitives; \
  CCTK_INT const tau_stildefix_enable = RESTRICTED_MHD_EVOLVE_STRUCT.tau_stildefix_enable; \
  CCTK_REAL const sdots_o_rhot = RESTRICTED_MHD_EVOLVE_STRUCT.sdots_o_rhot; \
  CCTK_INT const constrained_transport_scheme = RESTRICTED_MHD_EVOLVE_STRUCT.constrained_transport_scheme; \
  CCTK_INT const force_font_fix_fail = RESTRICTED_MHD_EVOLVE_STRUCT.force_font_fix_fail; \
  CCTK_INT const enable_primitives_disk = RESTRICTED_MHD_EVOLVE_STRUCT.enable_primitives_disk; \
  CCTK_INT const enable_HARM_energyvariable = RESTRICTED_MHD_EVOLVE_STRUCT.enable_HARM_energyvariable; \
  CCTK_REAL const rho_b_atm = RESTRICTED_MHD_EVOLVE_STRUCT.rho_b_atm; \
  CCTK_REAL const rho_fact = RESTRICTED_MHD_EVOLVE_STRUCT.rho_fact; \
  CCTK_REAL const tau_atm = RESTRICTED_MHD_EVOLVE_STRUCT.tau_atm; \
  CCTK_REAL const tau_fact = RESTRICTED_MHD_EVOLVE_STRUCT.tau_fact; \
  CCTK_REAL const P_fact = RESTRICTED_MHD_EVOLVE_STRUCT.P_fact; \
  CCTK_INT const neos = RESTRICTED_MHD_EVOLVE_STRUCT.neos; \
  CCTK_REAL const pfloor = RESTRICTED_MHD_EVOLVE_STRUCT.pfloor; \
  CCTK_REAL const gamma_th = RESTRICTED_MHD_EVOLVE_STRUCT.gamma_th; \
  CCTK_INT const eps_flag = RESTRICTED_MHD_EVOLVE_STRUCT.eps_flag; \
  CCTK_REAL const Sym_Bz = RESTRICTED_MHD_EVOLVE_STRUCT.Sym_Bz; \
  CCTK_REAL const K_poly = RESTRICTED_MHD_EVOLVE_STRUCT.K_poly; \
  CCTK_INT const em_evolve_enable = RESTRICTED_MHD_EVOLVE_STRUCT.em_evolve_enable; \
  CCTK_REAL const rho_b_max = RESTRICTED_MHD_EVOLVE_STRUCT.rho_b_max; \
  enum { \
      dummy_BSSN_xbh1_initial = sizeof( xbh1_initial ) \
    , dummy_BSSN_xbh2_initial = sizeof( xbh2_initial ) \
    , dummy_BSSN_rah1_initial = sizeof( rah1_initial ) \
    , dummy_BSSN_rah2_initial = sizeof( rah2_initial ) \
    , dummy_BSSN_binary_orb_freq = sizeof( binary_orb_freq ) \
    , dummy_BSSN_enable_lower_order_at_boundaries = sizeof( enable_lower_order_at_boundaries ) \
    , dummy_DIAGNOSTICS_MHD_N_rad_ray = sizeof( N_rad_ray ) \
    , dummy_DIAGNOSTICS_MHD_rhosurf_rmin = sizeof( rhosurf_rmin ) \
    , dummy_DIAGNOSTICS_MHD_rhosurf_rmax = sizeof( rhosurf_rmax ) \
    , dummy_DIAGNOSTICS_MHD_inner_lum_rad_ratio = sizeof( inner_lum_rad_ratio ) \
    , dummy_DIAGNOSTICS_MHD_lum_outer_rad = sizeof( lum_outer_rad ) \
    , dummy_DIAGNOSTICS_VACUUM_surf_radius = sizeof( surf_radius ) \
    , dummy_DIAGNOSTICS_VACUUM_rbr = sizeof( rbr ) \
    , dummy_DIAGNOSTICS_VACUUM_drbdr = sizeof( drbdr ) \
    , dummy_DIAGNOSTICS_VACUUM_ddrbddr = sizeof( ddrbddr ) \
    , dummy_DIAGNOSTICS_VACUUM_scaledsurf = sizeof( scaledsurf ) \
    , dummy_DIAGNOSTICS_VACUUM_nsurf = sizeof( nsurf ) \
    , dummy_DIAGNOSTICS_VACUUM_arithsurf = sizeof( arithsurf ) \
    , dummy_DIAGNOSTICS_VACUUM_rsurf2 = sizeof( rsurf2 ) \
    , dummy_DIAGNOSTICS_VACUUM_N_theta = sizeof( N_theta ) \
    , dummy_DIAGNOSTICS_VACUUM_N_phi = sizeof( N_phi ) \
    , dummy_DIAGNOSTICS_VACUUM_sym_factor = sizeof( sym_factor ) \
    , dummy_GRID_domain = sizeof( domain ) \
    , dummy_GRID_bitant_plane = sizeof( bitant_plane ) \
    , dummy_GW_EXTRACTION_radius_GW = sizeof( radius_GW ) \
    , dummy_GW_EXTRACTION_theta_GW = sizeof( theta_GW ) \
    , dummy_GW_EXTRACTION_phi_GW = sizeof( phi_GW ) \
    , dummy_IO_out_every = sizeof( out_every ) \
    , dummy_MHD_EVOLVE_primitives_solver = sizeof( primitives_solver ) \
    , dummy_MHD_EVOLVE_enable_shocktest_primitive_mode = sizeof( enable_shocktest_primitive_mode ) \
    , dummy_MHD_EVOLVE_use_HARM_primitives = sizeof( use_HARM_primitives ) \
    , dummy_MHD_EVOLVE_tau_stildefix_enable = sizeof( tau_stildefix_enable ) \
    , dummy_MHD_EVOLVE_sdots_o_rhot = sizeof( sdots_o_rhot ) \
    , dummy_MHD_EVOLVE_constrained_transport_scheme = sizeof( constrained_transport_scheme ) \
    , dummy_MHD_EVOLVE_force_font_fix_fail = sizeof( force_font_fix_fail ) \
    , dummy_MHD_EVOLVE_enable_primitives_disk = sizeof( enable_primitives_disk ) \
    , dummy_MHD_EVOLVE_enable_HARM_energyvariable = sizeof( enable_HARM_energyvariable ) \
    , dummy_MHD_EVOLVE_rho_b_atm = sizeof( rho_b_atm ) \
    , dummy_MHD_EVOLVE_rho_fact = sizeof( rho_fact ) \
    , dummy_MHD_EVOLVE_tau_atm = sizeof( tau_atm ) \
    , dummy_MHD_EVOLVE_tau_fact = sizeof( tau_fact ) \
    , dummy_MHD_EVOLVE_P_fact = sizeof( P_fact ) \
    , dummy_MHD_EVOLVE_neos = sizeof( neos ) \
    , dummy_MHD_EVOLVE_pfloor = sizeof( pfloor ) \
    , dummy_MHD_EVOLVE_gamma_th = sizeof( gamma_th ) \
    , dummy_MHD_EVOLVE_eps_flag = sizeof( eps_flag ) \
    , dummy_MHD_EVOLVE_Sym_Bz = sizeof( Sym_Bz ) \
    , dummy_MHD_EVOLVE_K_poly = sizeof( K_poly ) \
    , dummy_MHD_EVOLVE_em_evolve_enable = sizeof( em_evolve_enable ) \
    , dummy_MHD_EVOLVE_rho_b_max = sizeof( rho_b_max ) \
  };

#endif  /* _DISK_POWERLAW_PARAMETERS_H_ */

