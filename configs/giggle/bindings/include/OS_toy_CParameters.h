/*@@
   @header  OS_toy_CParameters.h
   @author  Automatically generated by CreateParameterBindings.pl
   @desc
            Declares parameters of thorn OS_toy
   @enddesc
 @@*/


#ifndef _OS_TOY_PARAMETERS_H_
#define _OS_TOY_PARAMETERS_H_ 1

#include "ParameterCGlobal.h"
#include "ParameterCPrivateOS_TOY.h"
#include "ParameterCRestrictedCACTUS.h"
#include "ParameterCRestrictedDIAGNOSTICS_MHD.h"
#include "ParameterCRestrictedDIAGNOSTICS_VACUUM.h"
#include "ParameterCRestrictedGRID.h"
#include "ParameterCRestrictedIO.h"
#include "ParameterCRestrictedMHD_EVOLVE.h"

#define DECLARE_CCTK_PARAMETERS \
  DECLARE_GLOBAL_PARAMETER_STRUCT_PARAMS \
  DECLARE_PRIVATE_OS_TOY_STRUCT_PARAMS \
  CCTK_INT const cctk_itlast = RESTRICTED_CACTUS_STRUCT.cctk_itlast; \
  CCTK_INT const N_rad_ray = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.N_rad_ray; \
  CCTK_REAL const rhosurf_rmin = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.rhosurf_rmin; \
  CCTK_REAL const rhosurf_rmax = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.rhosurf_rmax; \
  CCTK_REAL const surf_radius = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.surf_radius; \
  CCTK_INT const scaledsurf = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.scaledsurf; \
  CCTK_REAL const rbr = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.rbr; \
  CCTK_REAL const drbdr = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.drbdr; \
  CCTK_REAL const ddrbddr = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.ddrbddr; \
  CCTK_INT const nsurf = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.nsurf; \
  CCTK_INT const arithsurf = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.arithsurf; \
  CCTK_REAL const rsurf2 = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.rsurf2; \
  CCTK_INT const N_theta = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.N_theta; \
  CCTK_INT const N_phi = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.N_phi; \
  CCTK_INT const sym_factor = RESTRICTED_DIAGNOSTICS_VACUUM_STRUCT.sym_factor; \
  CCTK_INT const out_every = RESTRICTED_IO_STRUCT.out_every; \
  CCTK_INT const tau_stildefix_enable = RESTRICTED_MHD_EVOLVE_STRUCT.tau_stildefix_enable; \
  CCTK_INT const constrained_transport_scheme = RESTRICTED_MHD_EVOLVE_STRUCT.constrained_transport_scheme; \
  CCTK_INT const force_font_fix_fail = RESTRICTED_MHD_EVOLVE_STRUCT.force_font_fix_fail; \
  CCTK_INT const enable_primitives_disk = RESTRICTED_MHD_EVOLVE_STRUCT.enable_primitives_disk; \
  CCTK_INT const primitives_solver = RESTRICTED_MHD_EVOLVE_STRUCT.primitives_solver; \
  CCTK_INT const enable_shocktest_primitive_mode = RESTRICTED_MHD_EVOLVE_STRUCT.enable_shocktest_primitive_mode; \
  CCTK_INT const enable_HARM_energyvariable = RESTRICTED_MHD_EVOLVE_STRUCT.enable_HARM_energyvariable; \
  CCTK_REAL const rho_b_atm = RESTRICTED_MHD_EVOLVE_STRUCT.rho_b_atm; \
  CCTK_REAL const rho_fact = RESTRICTED_MHD_EVOLVE_STRUCT.rho_fact; \
  CCTK_REAL const tau_atm = RESTRICTED_MHD_EVOLVE_STRUCT.tau_atm; \
  CCTK_REAL const tau_fact = RESTRICTED_MHD_EVOLVE_STRUCT.tau_fact; \
  CCTK_REAL const pfloor = RESTRICTED_MHD_EVOLVE_STRUCT.pfloor; \
  CCTK_REAL const P_fact = RESTRICTED_MHD_EVOLVE_STRUCT.P_fact; \
  CCTK_INT const neos = RESTRICTED_MHD_EVOLVE_STRUCT.neos; \
  CCTK_REAL const gamma_th = RESTRICTED_MHD_EVOLVE_STRUCT.gamma_th; \
  CCTK_INT const eps_flag = RESTRICTED_MHD_EVOLVE_STRUCT.eps_flag; \
  CCTK_REAL const Sym_Bz = RESTRICTED_MHD_EVOLVE_STRUCT.Sym_Bz; \
  CCTK_REAL const K_poly = RESTRICTED_MHD_EVOLVE_STRUCT.K_poly; \
  CCTK_REAL const M_B = RESTRICTED_MHD_EVOLVE_STRUCT.M_B; \
  CCTK_INT const rad_evolve_enable = RESTRICTED_MHD_EVOLVE_STRUCT.rad_evolve_enable; \
  CCTK_REAL const rad_opacity_abs = RESTRICTED_MHD_EVOLVE_STRUCT.rad_opacity_abs; \
  CCTK_REAL const rad_opacity_sct = RESTRICTED_MHD_EVOLVE_STRUCT.rad_opacity_sct; \
  CCTK_INT const enable_OS_collapse = RESTRICTED_MHD_EVOLVE_STRUCT.enable_OS_collapse; \
  CCTK_INT const rad_closure_scheme = RESTRICTED_MHD_EVOLVE_STRUCT.rad_closure_scheme; \
  CCTK_REAL const Erad_cut = RESTRICTED_MHD_EVOLVE_STRUCT.Erad_cut; \
  CCTK_REAL const Erad_atm_cut = RESTRICTED_MHD_EVOLVE_STRUCT.Erad_atm_cut; \
  CCTK_REAL const gamma_OS = RESTRICTED_MHD_EVOLVE_STRUCT.gamma_OS; \
  enum { \
      dummy_CACTUS_cctk_itlast = sizeof( cctk_itlast ) \
    , dummy_DIAGNOSTICS_MHD_N_rad_ray = sizeof( N_rad_ray ) \
    , dummy_DIAGNOSTICS_MHD_rhosurf_rmin = sizeof( rhosurf_rmin ) \
    , dummy_DIAGNOSTICS_MHD_rhosurf_rmax = sizeof( rhosurf_rmax ) \
    , dummy_DIAGNOSTICS_VACUUM_surf_radius = sizeof( surf_radius ) \
    , dummy_DIAGNOSTICS_VACUUM_scaledsurf = sizeof( scaledsurf ) \
    , dummy_DIAGNOSTICS_VACUUM_rbr = sizeof( rbr ) \
    , dummy_DIAGNOSTICS_VACUUM_drbdr = sizeof( drbdr ) \
    , dummy_DIAGNOSTICS_VACUUM_ddrbddr = sizeof( ddrbddr ) \
    , dummy_DIAGNOSTICS_VACUUM_nsurf = sizeof( nsurf ) \
    , dummy_DIAGNOSTICS_VACUUM_arithsurf = sizeof( arithsurf ) \
    , dummy_DIAGNOSTICS_VACUUM_rsurf2 = sizeof( rsurf2 ) \
    , dummy_DIAGNOSTICS_VACUUM_N_theta = sizeof( N_theta ) \
    , dummy_DIAGNOSTICS_VACUUM_N_phi = sizeof( N_phi ) \
    , dummy_DIAGNOSTICS_VACUUM_sym_factor = sizeof( sym_factor ) \
    , dummy_IO_out_every = sizeof( out_every ) \
    , dummy_MHD_EVOLVE_tau_stildefix_enable = sizeof( tau_stildefix_enable ) \
    , dummy_MHD_EVOLVE_constrained_transport_scheme = sizeof( constrained_transport_scheme ) \
    , dummy_MHD_EVOLVE_force_font_fix_fail = sizeof( force_font_fix_fail ) \
    , dummy_MHD_EVOLVE_enable_primitives_disk = sizeof( enable_primitives_disk ) \
    , dummy_MHD_EVOLVE_primitives_solver = sizeof( primitives_solver ) \
    , dummy_MHD_EVOLVE_enable_shocktest_primitive_mode = sizeof( enable_shocktest_primitive_mode ) \
    , dummy_MHD_EVOLVE_enable_HARM_energyvariable = sizeof( enable_HARM_energyvariable ) \
    , dummy_MHD_EVOLVE_rho_b_atm = sizeof( rho_b_atm ) \
    , dummy_MHD_EVOLVE_rho_fact = sizeof( rho_fact ) \
    , dummy_MHD_EVOLVE_tau_atm = sizeof( tau_atm ) \
    , dummy_MHD_EVOLVE_tau_fact = sizeof( tau_fact ) \
    , dummy_MHD_EVOLVE_pfloor = sizeof( pfloor ) \
    , dummy_MHD_EVOLVE_P_fact = sizeof( P_fact ) \
    , dummy_MHD_EVOLVE_neos = sizeof( neos ) \
    , dummy_MHD_EVOLVE_gamma_th = sizeof( gamma_th ) \
    , dummy_MHD_EVOLVE_eps_flag = sizeof( eps_flag ) \
    , dummy_MHD_EVOLVE_Sym_Bz = sizeof( Sym_Bz ) \
    , dummy_MHD_EVOLVE_K_poly = sizeof( K_poly ) \
    , dummy_MHD_EVOLVE_M_B = sizeof( M_B ) \
    , dummy_MHD_EVOLVE_rad_evolve_enable = sizeof( rad_evolve_enable ) \
    , dummy_MHD_EVOLVE_rad_opacity_abs = sizeof( rad_opacity_abs ) \
    , dummy_MHD_EVOLVE_rad_opacity_sct = sizeof( rad_opacity_sct ) \
    , dummy_MHD_EVOLVE_enable_OS_collapse = sizeof( enable_OS_collapse ) \
    , dummy_MHD_EVOLVE_rad_closure_scheme = sizeof( rad_closure_scheme ) \
    , dummy_MHD_EVOLVE_Erad_cut = sizeof( Erad_cut ) \
    , dummy_MHD_EVOLVE_Erad_atm_cut = sizeof( Erad_atm_cut ) \
    , dummy_MHD_EVOLVE_gamma_OS = sizeof( gamma_OS ) \
  };

#endif  /* _OS_TOY_PARAMETERS_H_ */
