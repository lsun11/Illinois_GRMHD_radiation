/*@@
   @file    mhd_evolve.c
   @author  Automatically generated by GridFuncStuff.pl
   @desc
            Creates the CCTK variables for thorn mhd_evolve
   @enddesc
 @@*/


#define THORN_IS_mhd_evolve 1

#include <stddef.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameter.h"
#include "cctk_WarnLevel.h"
#include "cctki_Groups.h"
#include "cctki_FortranWrappers.h"

int CactusBindingsVariables_mhd_evolve_Initialise(void);
static int CCTKi_BindingsFortranWrappermhd_evolve(void *_GH, void *fpointer);

static int CCTKi_BindingsFortranWrappermhd_evolve(void *_GH, void *fpointer)
{
  cGH *GH = _GH;
  const int _cctk_zero = 0;
#ifndef CCTK_DEBUG
  CCTK_COMPLEX _cctk_dummy_var[4];
#endif
  void (*function)(MHD_EVOLVE_C2F_PROTO);
  DECLARE_MHD_EVOLVE_C2F
  INITIALISE_MHD_EVOLVE_C2F
  (void) (_cctk_zero + 0);
#ifndef CCTK_DEBUG
  (void) (_cctk_dummy_var + 0);
#endif

  function = (void (*) (MHD_EVOLVE_C2F_PROTO)) fpointer;
  function (PASS_MHD_EVOLVE_C2F (GH));

  return (0);
}

int CactusBindingsVariables_mhd_evolve_Initialise(void)
{
  const char * warn_mixeddim_gfs = "";
  int warn_mixeddim = 0;
  const CCTK_INT *allow_mixeddim_gfs;


  allow_mixeddim_gfs = CCTK_ParameterGet ("allow_mixeddim_gfs", "Cactus", 0);

  if (CCTKi_CreateGroup ("MHD_CONSERVATIVES", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 3,
                         "NONE", "DEFAULT",
                         "", "",
                         "prolongation=\"Lag3\"",
                         NULL,
                         5,
                         "rho_star",
                         "tau",
                         "mhd_st_x",
                         "mhd_st_y",
                         "mhd_st_z") == 1)
  {
    warn_mixeddim_gfs = "mhd_conservatives";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("RAD_CONSERVATIVES", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 3,
                         "NONE", "DEFAULT",
                         "", "",
                         "prolongation=\"Lag3\"",
                         NULL,
                         12,
                         "tau_rad",
                         "S_rad_x",
                         "S_rad_y",
                         "S_rad_z",
                         "tau_rad_nue",
                         "S_rad_x_nue",
                         "S_rad_y_nue",
                         "S_rad_z_nue",
                         "tau_rad_nux",
                         "S_rad_x_nux",
                         "S_rad_y_nux",
                         "S_rad_z_nux") == 1)
  {
    warn_mixeddim_gfs = "rad_conservatives";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("MICPHYS_CONSERVATIVES", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 3,
                         "NONE", "DEFAULT",
                         "", "",
                         "prolongation=\"Lag3\"",
                         NULL,
                         1,
                         "rhoYe") == 1)
  {
    warn_mixeddim_gfs = "micphys_conservatives";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("EM_CONSERVATIVEX", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 3,
                         "NONE", "DEFAULT",
                         "", "",
                         "prolongation=\"Lag3\"",
                         NULL,
                         1,
                         "Bxtilde") == 1)
  {
    warn_mixeddim_gfs = "em_conservativex";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("EM_CONSERVATIVEY", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 3,
                         "NONE", "DEFAULT",
                         "", "",
                         "prolongation=\"Lag3\"",
                         NULL,
                         1,
                         "Bytilde") == 1)
  {
    warn_mixeddim_gfs = "em_conservativey";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("EM_CONSERVATIVEZ", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 3,
                         "NONE", "DEFAULT",
                         "", "",
                         "prolongation=\"Lag3\"",
                         NULL,
                         1,
                         "Bztilde") == 1)
  {
    warn_mixeddim_gfs = "em_conservativez";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("EM_BLAGRANGEMULTIPLIER", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 3,
                         "NONE", "DEFAULT",
                         "", "",
                         "prolongation=\"Lag3\"",
                         NULL,
                         1,
                         "Blagrangemultiplier") == 1)
  {
    warn_mixeddim_gfs = "em_Blagrangemultiplier";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("EM_AX", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 3,
                         "NONE", "DEFAULT",
                         "", "",
                         "Prolongation=\"STAGGER011\"",
                         NULL,
                         1,
                         "Ax") == 1)
  {
    warn_mixeddim_gfs = "em_Ax";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("EM_AY", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 3,
                         "NONE", "DEFAULT",
                         "", "",
                         "Prolongation=\"STAGGER101\"",
                         NULL,
                         1,
                         "Ay") == 1)
  {
    warn_mixeddim_gfs = "em_Ay";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("EM_AZ", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 3,
                         "NONE", "DEFAULT",
                         "", "",
                         "Prolongation=\"STAGGER110\"",
                         NULL,
                         1,
                         "Az") == 1)
  {
    warn_mixeddim_gfs = "em_Az";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("EM_PHI", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 3,
                         "NONE", "DEFAULT",
                         "", "",
                         "Prolongation=\"STAGGER111\"",
                         NULL,
                         1,
                         "psi6phi") == 1)
  {
    warn_mixeddim_gfs = "em_Phi";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("EM_PHI_RHS", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         1,
                         "psi6phi_rhs") == 1)
  {
    warn_mixeddim_gfs = "em_Phi_rhs";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("STAGGER_BS", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "InterpNumTimelevels=1",
                         NULL,
                         3,
                         "Bx_stagger",
                         "By_stagger",
                         "Bz_stagger") == 1)
  {
    warn_mixeddim_gfs = "Stagger_Bs";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("MHD_RHS", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         5,
                         "rho_star_rhs",
                         "tau_rhs",
                         "mhd_st_x_rhs",
                         "mhd_st_y_rhs",
                         "mhd_st_z_rhs") == 1)
  {
    warn_mixeddim_gfs = "mhd_rhs";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("RAD_CONSERVATIVES_RHS", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         12,
                         "tau_rad_rhs",
                         "S_rad_x_rhs",
                         "S_rad_y_rhs",
                         "S_rad_z_rhs",
                         "tau_rad_nue_rhs",
                         "S_rad_x_nue_rhs",
                         "S_rad_y_nue_rhs",
                         "S_rad_z_nue_rhs",
                         "tau_rad_nux_rhs",
                         "S_rad_x_nux_rhs",
                         "S_rad_y_nux_rhs",
                         "S_rad_z_nux_rhs") == 1)
  {
    warn_mixeddim_gfs = "rad_conservatives_rhs";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("MICPHYS_CONSERVATIVES_RHS", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         1,
                         "rhoYe_rhs") == 1)
  {
    warn_mixeddim_gfs = "micphys_conservatives_rhs";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("EM_RHSX", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         1,
                         "Bxtilde_or_Ax_rhs") == 1)
  {
    warn_mixeddim_gfs = "em_rhsx";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("EM_RHSY", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         1,
                         "Bytilde_or_Ay_rhs") == 1)
  {
    warn_mixeddim_gfs = "em_rhsy";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("EM_RHSZ", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         1,
                         "Bztilde_or_Az_rhs") == 1)
  {
    warn_mixeddim_gfs = "em_rhsz";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("EM_BLAGRANGEMULTIPLIER_RHS", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         1,
                         "Blagrangemultiplier_rhs") == 1)
  {
    warn_mixeddim_gfs = "em_Blagrangemultiplier_rhs";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("FIELD_LINE_VARIABLES", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 3,
                         "NONE", "DEFAULT",
                         "", "",
                         "",
                         NULL,
                         4,
                         "mhd_psi_line",
                         "mhd_u_psi",
                         "mhd_chi_line",
                         "mhd_u_chi") == 1)
  {
    warn_mixeddim_gfs = "field_line_variables";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("FIELD_LINE_VARIABLES_RHS", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         4,
                         "mhd_psi_line_rhs",
                         "mhd_u_psi_rhs",
                         "mhd_chi_line_rhs",
                         "mhd_u_chi_rhs") == 1)
  {
    warn_mixeddim_gfs = "field_line_variables_rhs";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("MHD_PRIMITIVES", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "InterpNumTimelevels=1 prolongation=\"none\"",
                         NULL,
                         19,
                         "h_p",
                         "h",
                         "P",
                         "rho_b",
                         "w",
                         "st_x",
                         "st_y",
                         "st_z",
                         "Bx",
                         "By",
                         "Bz",
                         "Ex",
                         "Ey",
                         "Ez",
                         "sbt",
                         "sbx",
                         "sby",
                         "sbz",
                         "smallb2") == 1)
  {
    warn_mixeddim_gfs = "mhd_primitives";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("RAD_PRIMITIVES", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "InterpNumTimelevels=1 prolongation=\"none\"",
                         NULL,
                         12,
                         "E_rad",
                         "F_radx",
                         "F_rady",
                         "F_radz",
                         "E_rad_nue",
                         "F_radx_nue",
                         "F_rady_nue",
                         "F_radz_nue",
                         "E_rad_nux",
                         "F_radx_nux",
                         "F_rady_nux",
                         "F_radz_nux") == 1)
  {
    warn_mixeddim_gfs = "rad_primitives";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("RAD_PRESSURE", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "InterpNumTimelevels=1 prolongation=\"none\"",
                         NULL,
                         37,
                         "F_rad0",
                         "P_radxx",
                         "P_radyy",
                         "P_radzz",
                         "P_radxy",
                         "P_radxz",
                         "P_radyz",
                         "chi_rad",
                         "zeta_rad",
                         "F_rad_scalar",
                         "FaFar",
                         "FaFal",
                         "F_rad0_nue",
                         "P_radxx_nue",
                         "P_radyy_nue",
                         "P_radzz_nue",
                         "P_radxy_nue",
                         "P_radxz_nue",
                         "P_radyz_nue",
                         "chi_rad_nue",
                         "zeta_rad_nue",
                         "F_rad_scalar_nue",
                         "FaFar_nue",
                         "FaFal_nue",
                         "F_rad0_nux",
                         "P_radxx_nux",
                         "P_radyy_nux",
                         "P_radzz_nux",
                         "P_radxy_nux",
                         "P_radxz_nux",
                         "P_radyz_nux",
                         "chi_rad_nux",
                         "zeta_rad_nux",
                         "F_rad_scalar_nux",
                         "FaFar_nux",
                         "FaFal_nux",
                         "eta_nue") == 1)
  {
    warn_mixeddim_gfs = "rad_pressure";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("MHD_VS", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "InterpNumTimelevels=1 prolongation=\"none\"",
                         NULL,
                         3,
                         "vx",
                         "vy",
                         "vz") == 1)
  {
    warn_mixeddim_gfs = "mhd_vs";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("DISK_ATMOSPHERE", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "",
                         NULL,
                         3,
                         "rho_b_atm_gf",
                         "pfloor_gf",
                         "Fontfix_tracker_gf") == 1)
  {
    warn_mixeddim_gfs = "disk_atmosphere";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("MHD_SYNC_NABLAS", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         8,
                         "drho_b_m",
                         "dP_m",
                         "dvx_m",
                         "dvy_m",
                         "dvz_m",
                         "dBx_m",
                         "dBy_m",
                         "dBz_m") == 1)
  {
    warn_mixeddim_gfs = "mhd_sync_nablas";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("MICPHYS_SYNC_NABLAS", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         2,
                         "drhoYe_m",
                         "dT_fluid_m") == 1)
  {
    warn_mixeddim_gfs = "micphys_sync_nablas";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("MHD_SYNC_NABLAS_DIAG", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         6,
                         "drho_b_m_x",
                         "dvx_m_x",
                         "drho_b_m_xp1",
                         "dvx_m_xp1",
                         "drhoYe_m_x",
                         "drhoYe_m_xp1") == 1)
  {
    warn_mixeddim_gfs = "mhd_sync_nablas_diag";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("MHD_SYNC_RHO_BR_RHO_BL", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         2,
                         "rho_br",
                         "rho_bl") == 1)
  {
    warn_mixeddim_gfs = "mhd_sync_rho_br_rho_bl";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("MHD_SYNC_METRIC_FACEVALS", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         17,
                         "lapm1_f",
                         "shiftx_f",
                         "shifty_f",
                         "shiftz_f",
                         "gxx_f",
                         "gxy_f",
                         "gxz_f",
                         "gyy_f",
                         "gyz_f",
                         "gzz_f",
                         "phi_f",
                         "gupxx_f",
                         "gupyy_f",
                         "gupzz_f",
                         "gupxy_f",
                         "gupxz_f",
                         "gupyz_f") == 1)
  {
    warn_mixeddim_gfs = "mhd_sync_metric_facevals";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("MHD_SYNC_LR_HYDRO_QUANTITIES", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         8,
                         "Pr",
                         "Pl",
                         "vxr",
                         "vxl",
                         "vyr",
                         "vyl",
                         "vzr",
                         "vzl") == 1)
  {
    warn_mixeddim_gfs = "mhd_sync_lr_hydro_quantities";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("MICPHYS_SYNC_LR_HYDRO_QUANTITIES", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         4,
                         "Y_er",
                         "Y_el",
                         "T_fluidle",
                         "T_fluidr") == 1)
  {
    warn_mixeddim_gfs = "micphys_sync_lr_hydro_quantities";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("MHD_SYNC_LR_B_QUANTITIES", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         6,
                         "Bxr",
                         "Bxl",
                         "Byr",
                         "Byl",
                         "Bzr",
                         "Bzl") == 1)
  {
    warn_mixeddim_gfs = "mhd_sync_lr_B_quantities";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("RAD_SYNC_LR", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         24,
                         "E_radr",
                         "E_radl",
                         "F_radxr",
                         "F_radxle",
                         "F_radyr",
                         "F_radyle",
                         "F_radzr",
                         "F_radzle",
                         "E_rad_nuer",
                         "E_rad_nuel",
                         "F_radx_nuer",
                         "F_radx_nuele",
                         "F_rady_nuer",
                         "F_rady_nuele",
                         "F_radz_nuer",
                         "F_radz_nuele",
                         "E_rad_nuxr",
                         "E_rad_nuxl",
                         "F_radx_nuxr",
                         "F_radx_nuxle",
                         "F_rady_nuxr",
                         "F_rady_nuxle",
                         "F_radz_nuxr",
                         "F_radz_nuxle") == 1)
  {
    warn_mixeddim_gfs = "rad_sync_lr";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("MHD_NOSYNC", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "InterpNumTimelevels=1 prolongation=\"none\"",
                         NULL,
                         19,
                         "v02r",
                         "v02l",
                         "cmin",
                         "cmax",
                         "cmin_rad",
                         "cmax_rad",
                         "cmin_rad_nue",
                         "cmax_rad_nue",
                         "cmin_rad_nux",
                         "cmax_rad_nux",
                         "u0",
                         "rhob_floor",
                         "P_floor",
                         "v02_radr",
                         "v02_radl",
                         "v02_rad_nuer",
                         "v02_rad_nuel",
                         "v02_rad_nuxr",
                         "v02_rad_nuxl") == 1)
  {
    warn_mixeddim_gfs = "mhd_nosync";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("EOS_PARAMS1", "MHD_EVOLVE", "MHD_EVOLVE",
                         "ARRAY", "REAL", "PROTECTED",
                         1, 1,
                         "NONE", "CONSTANT",
                         "10", "",
                         "",
                         NULL,
                         3,
                         "rho_tab",
                         "P_tab",
                         "eps_tab") == 1)
  {
    warn_mixeddim_gfs = "eos_params1";
    warn_mixeddim = 1;
  }
  if (CCTKi_CreateGroup ("EOS_PARAMS2", "MHD_EVOLVE", "MHD_EVOLVE",
                         "ARRAY", "REAL", "PROTECTED",
                         1, 1,
                         "NONE", "CONSTANT",
                         "11", "",
                         "",
                         NULL,
                         2,
                         "k_tab",
                         "gamma_tab") == 1)
  {
    warn_mixeddim_gfs = "eos_params2";
    warn_mixeddim = 1;
  }
  if (CCTKi_CreateGroup ("MHDSCALAR", "MHD_EVOLVE", "MHD_EVOLVE",
                         "SCALAR", "REAL", "PROTECTED",
                         0, 1,
                         "NONE", "CONSTANT",
                         "", "",
                         "",
                         NULL,
                         1,
                         "n_poly") == 1)
  {
    warn_mixeddim_gfs = "mhdscalar";
    warn_mixeddim = 0;
  }
  if (CCTKi_CreateGroup ("RHOVECS", "MHD_EVOLVE", "MHD_EVOLVE",
                         "ARRAY", "REAL", "PROTECTED",
                         1, 1,
                         "NONE", "CONSTANT",
                         "1000", "",
                         "",
                         NULL,
                         3,
                         "rhovec",
                         "Pvec",
                         "vvec") == 1)
  {
    warn_mixeddim_gfs = "rhovecs";
    warn_mixeddim = 1;
  }
  if (CCTKi_CreateGroup ("MHD_TEMPS", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\" InterpNumTimelevels=1 prolongation=\"none\"",
                         NULL,
                         29,
                         "temp1",
                         "temp2",
                         "temp3",
                         "temp4",
                         "temp5",
                         "temp6",
                         "temp7",
                         "temp8",
                         "temp_g00",
                         "temp9",
                         "temp10",
                         "temp11",
                         "MONOPOLE",
                         "P_thermal",
                         "temp12",
                         "temp13",
                         "temp14",
                         "temp15",
                         "temp16",
                         "temp17",
                         "temp18",
                         "temp19",
                         "temp20",
                         "temp21",
                         "temp22",
                         "temp23",
                         "temp24",
                         "temp25",
                         "temp26") == 1)
  {
    warn_mixeddim_gfs = "mhd_temps";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("MICROPHYS_PRIMITIVES", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PROTECTED",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "InterpNumTimelevels=1 prolongation=\"none\"",
                         NULL,
                         20,
                         "Y_e",
                         "ka_gf",
                         "ks_gf",
                         "emission_gf",
                         "ka_gf_nue",
                         "ks_gf_nue",
                         "emission_gf_nue",
                         "ka_gf_nux",
                         "ks_gf_nux",
                         "emission_gf_nux",
                         "mu_nu",
                         "eps_thermal",
                         "eps_cld",
                         "P_cld",
                         "optd_x",
                         "optd_y",
                         "optd_z",
                         "optd",
                         "eps_tot",
                         "T_fluid") == 1)
  {
    warn_mixeddim_gfs = "microphys_primitives";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("RADSCALAR", "MHD_EVOLVE", "MHD_EVOLVE",
                         "SCALAR", "REAL", "PROTECTED",
                         0, 1,
                         "NONE", "CONSTANT",
                         "", "",
                         "",
                         NULL,
                         1,
                         "rad_const") == 1)
  {
    warn_mixeddim_gfs = "radscalar";
    warn_mixeddim = 0;
  }
  if (CCTKi_CreateGroup ("OS_STELLAR_SURFACE", "MHD_EVOLVE", "MHD_EVOLVE",
                         "SCALAR", "REAL", "PROTECTED",
                         0, 1,
                         "NONE", "CONSTANT",
                         "", "",
                         "",
                         NULL,
                         1,
                         "OS_surf_rad") == 1)
  {
    warn_mixeddim_gfs = "OS_stellar_surface";
    warn_mixeddim = 0;
  }
  if (CCTKi_CreateGroup ("EM_FIJS", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PRIVATE",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         6,
                         "fxy",
                         "fxz",
                         "fyx",
                         "fyz",
                         "fzx",
                         "fzy") == 1)
  {
    warn_mixeddim_gfs = "em_fijs";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("EM_FTIJS", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PRIVATE",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         6,
                         "ftxy",
                         "ftxz",
                         "ftyx",
                         "ftyz",
                         "ftzx",
                         "ftzy") == 1)
  {
    warn_mixeddim_gfs = "em_ftijs";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("RAD_FIJS", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PRIVATE",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         12,
                         "tau_rad_flux",
                         "S_radx_flux",
                         "S_rady_flux",
                         "S_radz_flux",
                         "tau_rad_nue_flux",
                         "S_radx_nue_flux",
                         "S_rady_nue_flux",
                         "S_radz_nue_flux",
                         "tau_rad_nux_flux",
                         "S_radx_nux_flux",
                         "S_rady_nux_flux",
                         "S_radz_nux_flux") == 1)
  {
    warn_mixeddim_gfs = "rad_fijs";
    warn_mixeddim = 3;
  }
  if (CCTKi_CreateGroup ("RAD_FIJS_DIAG", "MHD_EVOLVE", "MHD_EVOLVE",
                         "GF", "REAL", "PRIVATE",
                         3, 1,
                         "NONE", "DEFAULT",
                         "", "",
                         "Checkpoint=\"no\"",
                         NULL,
                         6,
                         "tau_rad_flux_x",
                         "S_radx_flux_x",
                         "tau_rad_advect_diag",
                         "S_radx_flux_xp1",
                         "tau_rad_flux_diag",
                         "tau_rad_scalar_diag") == 1)
  {
    warn_mixeddim_gfs = "rad_fijs_diag";
    warn_mixeddim = 3;
  }

  if (*warn_mixeddim_gfs)
  {
    if (allow_mixeddim_gfs && *allow_mixeddim_gfs)
    {
      CCTK_VWarn (2, __LINE__, __FILE__, "Cactus",
                  "CCTKi_CreateGroup: Working dimension already set, "
                  "creating GF group '%s' with different dimension %d",
                  warn_mixeddim_gfs, warn_mixeddim);
    }
    else
    {
      CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                  "CCTKi_CreateGroup: Working dimension already set, "
                  "cannot create GF group '%s' with dimension %d",
                  warn_mixeddim_gfs, warn_mixeddim);
    }
 }

  CCTKi_RegisterFortranWrapper("mhd_evolve", CCTKi_BindingsFortranWrappermhd_evolve);

  return 0;
}
