/*@@
   @file    TwoPuncturesAEI_Parameters.c
   @author  Automatically generated by CreateParameterBindings.pl
   @desc
            Creates/extends parameters for this thorn
   @enddesc
 @@*/


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "cctk_Config.h"
#include "cctk_Constants.h"
#include "ParameterBindings.h"
#include "CParameterStructNames.h"
#include "ParameterCRestrictedDIAGNOSTICS_VACUUM.h"
#include "ParameterCRestrictedGRID.h"
#include "ParameterCRestrictedGW_EXTRACTION.h"
#include "ParameterCRestrictedIO.h"
#include "ParameterCRestrictedSTATICCONFORMAL.h"

/* structure containing all private parameters of thorn TwoPuncturesAEI */
struct
{
  CCTK_REAL Newton_tol;
  CCTK_REAL TP_Extend_Radius;
  CCTK_REAL TP_Tiny;
  CCTK_REAL TP_epsilon;
  CCTK_REAL adm_tol;
  CCTK_REAL center_offset[3];
  CCTK_REAL initial_lapse_psi_exponent;
  CCTK_REAL moncrief_radius_GW[100];
  CCTK_REAL par_P_minus[3];
  CCTK_REAL par_P_plus[3];
  CCTK_REAL par_S_minus[3];
  CCTK_REAL par_S_plus[3];
  CCTK_REAL par_b;
  CCTK_REAL par_m_minus;
  CCTK_REAL par_m_plus;
  CCTK_REAL target_M_minus;
  CCTK_REAL target_M_plus;
  const char * grid_setup_method;
  const char * initial_lapse;
  CCTK_INT Newton_maxit;
  CCTK_INT do_initial_debug_output;
  CCTK_INT do_residuum_debug_output;
  CCTK_INT give_bare_mass;
  CCTK_INT keep_u_around;
  CCTK_INT moncrief_gw_num_radii;
  CCTK_INT multiply_old_lapse;
  CCTK_INT npoints_A;
  CCTK_INT npoints_B;
  CCTK_INT npoints_phi;
  CCTK_INT rescale_sources;
  CCTK_INT schedule_in_ADMBase_InitialData;
  CCTK_INT solve_momentum_constraint;
  CCTK_INT swap_xz;
  CCTK_INT use_external_initial_guess;
  CCTK_INT use_sources;
  CCTK_INT verbose;
} PRIVATE_TWOPUNCTURESAEI_STRUCT;


/* structure containing all restricted parameters of thorn TwoPuncturesAEI */
struct
{
  int dummy_parameter;
} RESTRICTED_TWOPUNCTURESAEI_STRUCT;


int CCTKi_BindingsCreateTwoPuncturesAEIParameters(void);
int CCTKi_BindingsCreateTwoPuncturesAEIParameters(void)
{
  CCTKi_ParameterCreate("Newton_maxit",
                        "TwoPuncturesAEI",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Maximum number of Newton iterations",
                        "5",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.Newton_maxit),
                        0,
                        NULL,
                        1,
                        "0:*", "");

  CCTKi_ParameterCreate("Newton_tol",
                        "TwoPuncturesAEI",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Tolerance for Newton solver",
                        "1.0e-10",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.Newton_tol),
                        0,
                        NULL,
                        1,
                        "(0:*)", "");

  CCTKi_ParameterCreate("TP_Extend_Radius",
                        "TwoPuncturesAEI",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Radius of an extended spacetime instead of the puncture",
                        "0.0",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.TP_Extend_Radius),
                        0,
                        NULL,
                        1,
                        "0:*", "anything positive, should be smaller than the horizon");

  CCTKi_ParameterCreate("TP_Tiny",
                        "TwoPuncturesAEI",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Tiny number to avoid nans near or at the pucture locations",
                        "0.0",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.TP_Tiny),
                        0,
                        NULL,
                        1,
                        "0:*", "anything positive, usually very small");

  CCTKi_ParameterCreate("TP_epsilon",
                        "TwoPuncturesAEI",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "A small number to smooth out singularities at the puncture locations",
                        "0.0",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.TP_epsilon),
                        0,
                        NULL,
                        1,
                        "0:*", "");

  CCTKi_ParameterCreate("adm_tol",
                        "TwoPuncturesAEI",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Tolerance of ADM masses when give_bare_mass=no",
                        "1.0e-10",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.adm_tol),
                        0,
                        NULL,
                        1,
                        "(0:*)", "");

  CCTKi_ParameterCreate("center_offset",
                        "TwoPuncturesAEI",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "offset b=0 to position (x,y,z)",
                        "0.0",
                        (PRIVATE_TWOPUNCTURESAEI_STRUCT.center_offset),
                        3,
                        NULL,
                        1,
                        "(*:*)", "");

  CCTKi_ParameterCreate("do_initial_debug_output",
                        "TwoPuncturesAEI",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Output debug information about initial guess",
                        "no",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.do_initial_debug_output),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("do_residuum_debug_output",
                        "TwoPuncturesAEI",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Output debug information about the residuum",
                        "no",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.do_residuum_debug_output),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("give_bare_mass",
                        "TwoPuncturesAEI",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "User provides bare masses rather than target ADM masses",
                        "yes",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.give_bare_mass),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("grid_setup_method",
                        "TwoPuncturesAEI",
                        "KEYWORD",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "How to fill the 3D grid from the spectral grid",
                        "evaluation",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.grid_setup_method),
                        0,
                        NULL,
                        2,
                        "Taylor expansion", "use a Taylor expansion about the nearest collocation point (fast, but might be inaccurate)",
                        "evaluation", "evaluate using all spectral coefficients (slow)");

  CCTKi_ParameterCreate("initial_lapse",
                        "TwoPuncturesAEI",
                        "KEYWORD",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "initial lapse, duh.",
                        "twopunctures-averaged",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.initial_lapse),
                        0,
                        NULL,
                        4,
                        "twopunctures-antisymmetric", "antisymmetric lapse for two puncture black holes, -1 <= alpha <= +1",
                        "twopunctures-averaged", "averaged lapse for two puncture black holes, 0 <= alpha <= +1",
                        "psi^n", "Based on the initial conformal factor",
                        "brownsville", "See Phys. Rev. D 74, 041501 (2006)");

  CCTKi_ParameterCreate("initial_lapse_psi_exponent",
                        "TwoPuncturesAEI",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Exponent n for psi^-n initial lapse profile",
                        "-2.0",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.initial_lapse_psi_exponent),
                        0,
                        NULL,
                        1,
                        "(*:*)", "Should be negative");

  CCTKi_ParameterCreate("keep_u_around",
                        "TwoPuncturesAEI",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Keep the variable u around after solving",
                        "no",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.keep_u_around),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("moncrief_gw_num_radii",
                        "TwoPuncturesAEI",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "How many radii will to measure GW's",
                        "1",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.moncrief_gw_num_radii),
                        0,
                        NULL,
                        1,
                        "1:10", "Positive int <= 10");

  CCTKi_ParameterCreate("moncrief_radius_GW",
                        "TwoPuncturesAEI",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_RECOVER,
                        "Radii at which to measure GW's",
                        "1.0",
                        (PRIVATE_TWOPUNCTURESAEI_STRUCT.moncrief_radius_GW),
                        100,
                        NULL,
                        1,
                        "0:*", "zero or any positive number");

  CCTKi_ParameterCreate("multiply_old_lapse",
                        "TwoPuncturesAEI",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Multiply the old lapse with the new one",
                        "no",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.multiply_old_lapse),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("npoints_A",
                        "TwoPuncturesAEI",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Number of coefficients in the compactified radial direction",
                        "30",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.npoints_A),
                        0,
                        NULL,
                        1,
                        "4:*", "");

  CCTKi_ParameterCreate("npoints_B",
                        "TwoPuncturesAEI",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Number of coefficients in the angular direction",
                        "30",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.npoints_B),
                        0,
                        NULL,
                        1,
                        "4:*", "");

  CCTKi_ParameterCreate("npoints_phi",
                        "TwoPuncturesAEI",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Number of coefficients in the phi direction",
                        "16",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.npoints_phi),
                        0,
                        NULL,
                        1,
                        "4:*:2", "");

  CCTKi_ParameterCreate("par_P_minus",
                        "TwoPuncturesAEI",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "momentum of the m- puncture",
                        "0.0",
                        (PRIVATE_TWOPUNCTURESAEI_STRUCT.par_P_minus),
                        3,
                        NULL,
                        1,
                        "(*:*)", "");

  CCTKi_ParameterCreate("par_P_plus",
                        "TwoPuncturesAEI",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "momentum of the m+ puncture",
                        "0.0",
                        (PRIVATE_TWOPUNCTURESAEI_STRUCT.par_P_plus),
                        3,
                        NULL,
                        1,
                        "(*:*)", "");

  CCTKi_ParameterCreate("par_S_minus",
                        "TwoPuncturesAEI",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "spin of the m- puncture",
                        "0.0",
                        (PRIVATE_TWOPUNCTURESAEI_STRUCT.par_S_minus),
                        3,
                        NULL,
                        1,
                        "(*:*)", "");

  CCTKi_ParameterCreate("par_S_plus",
                        "TwoPuncturesAEI",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "spin of the m+ puncture",
                        "0.0",
                        (PRIVATE_TWOPUNCTURESAEI_STRUCT.par_S_plus),
                        3,
                        NULL,
                        1,
                        "(*:*)", "");

  CCTKi_ParameterCreate("par_b",
                        "TwoPuncturesAEI",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "x coordinate of the m+ puncture",
                        "1.0",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.par_b),
                        0,
                        NULL,
                        1,
                        "(0.0:*)", "");

  CCTKi_ParameterCreate("par_m_minus",
                        "TwoPuncturesAEI",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "mass of the m- puncture",
                        "1.0",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.par_m_minus),
                        0,
                        NULL,
                        1,
                        "0.0:*)", "");

  CCTKi_ParameterCreate("par_m_plus",
                        "TwoPuncturesAEI",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "mass of the m+ puncture",
                        "1.0",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.par_m_plus),
                        0,
                        NULL,
                        1,
                        "0.0:*)", "");

  CCTKi_ParameterCreate("rescale_sources",
                        "TwoPuncturesAEI",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "If sources are used - rescale them after solving?",
                        "yes",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.rescale_sources),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("schedule_in_ADMBase_InitialData",
                        "TwoPuncturesAEI",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Schedule in (instead of after) ADMBase_InitialData",
                        "yes",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.schedule_in_ADMBase_InitialData),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("solve_momentum_constraint",
                        "TwoPuncturesAEI",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Solve for momentum constraint?",
                        "no",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.solve_momentum_constraint),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("swap_xz",
                        "TwoPuncturesAEI",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Swap x and z coordinates when interpolating, so that the black holes are separated in the z direction",
                        "no",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.swap_xz),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("target_M_minus",
                        "TwoPuncturesAEI",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "target ADM mass for m-",
                        "0.5",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.target_M_minus),
                        0,
                        NULL,
                        1,
                        "0.0:*)", "");

  CCTKi_ParameterCreate("target_M_plus",
                        "TwoPuncturesAEI",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "target ADM mass for m+",
                        "0.5",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.target_M_plus),
                        0,
                        NULL,
                        1,
                        "0.0:*)", "");

  CCTKi_ParameterCreate("use_external_initial_guess",
                        "TwoPuncturesAEI",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Set initial guess by external function?",
                        "no",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.use_external_initial_guess),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("use_sources",
                        "TwoPuncturesAEI",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Use sources?",
                        "no",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.use_sources),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("verbose",
                        "TwoPuncturesAEI",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Print screen output while solving",
                        "no",
                        &(PRIVATE_TWOPUNCTURESAEI_STRUCT.verbose),
                        0,
                        NULL,
                        0);

  return 0;
}

int CCTKi_BindingsTwoPuncturesAEIParameterExtensions(void);
int CCTKi_BindingsTwoPuncturesAEIParameterExtensions(void)
{
  return 0;
}

