/*@@
   @file    TwoPunctures_Parameters.c
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

/* structure containing all private parameters of thorn TwoPunctures */
struct
{
  CCTK_REAL bh_mass_minus;
  CCTK_REAL bh_mass_plus;
  CCTK_REAL bh_px_minus;
  CCTK_REAL bh_px_plus;
  CCTK_REAL bh_py_minus;
  CCTK_REAL bh_py_plus;
  CCTK_REAL bh_spin_minus;
  CCTK_REAL bh_spin_plus;
  CCTK_REAL excis_radius;
  CCTK_REAL half_binary_separation;
  CCTK_REAL moncrief_radius_GW[100];
  CCTK_REAL x_offset;
  CCTK_INT fill_excision_enable;
  CCTK_INT genID_cmdline_output_enable;
  CCTK_INT moncrief_gw_num_radii;
} PRIVATE_TWOPUNCTURES_STRUCT;


/* structure containing all restricted parameters of thorn TwoPunctures */
struct
{
  int dummy_parameter;
} RESTRICTED_TWOPUNCTURES_STRUCT;


int CCTKi_BindingsCreateTwoPuncturesParameters(void);
int CCTKi_BindingsCreateTwoPuncturesParameters(void)
{
  CCTKi_ParameterCreate("bh_mass_minus",
                        "TwoPunctures",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Mass of the black hole on -x axis",
                        "0.0",
                        &(PRIVATE_TWOPUNCTURES_STRUCT.bh_mass_minus),
                        0,
                        NULL,
                        1,
                        "0:*", "Positive");

  CCTKi_ParameterCreate("bh_mass_plus",
                        "TwoPunctures",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Mass of the black hole on +x axis",
                        "0.5",
                        &(PRIVATE_TWOPUNCTURES_STRUCT.bh_mass_plus),
                        0,
                        NULL,
                        1,
                        "0:*", "Positive");

  CCTKi_ParameterCreate("bh_px_minus",
                        "TwoPunctures",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "x-component of BH momentum on -x axis",
                        "0.0",
                        &(PRIVATE_TWOPUNCTURES_STRUCT.bh_px_minus),
                        0,
                        NULL,
                        1,
                        "*:*", "Real");

  CCTKi_ParameterCreate("bh_px_plus",
                        "TwoPunctures",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "x-component of BH momentum on +x axis",
                        "0.0",
                        &(PRIVATE_TWOPUNCTURES_STRUCT.bh_px_plus),
                        0,
                        NULL,
                        1,
                        "*:*", "Real");

  CCTKi_ParameterCreate("bh_py_minus",
                        "TwoPunctures",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "y-component of BH momentum on -x axis",
                        "0.0",
                        &(PRIVATE_TWOPUNCTURES_STRUCT.bh_py_minus),
                        0,
                        NULL,
                        1,
                        "*:*", "Real");

  CCTKi_ParameterCreate("bh_py_plus",
                        "TwoPunctures",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "y-component of BH momentum on +x axis",
                        "0.0",
                        &(PRIVATE_TWOPUNCTURES_STRUCT.bh_py_plus),
                        0,
                        NULL,
                        1,
                        "*:*", "Real");

  CCTKi_ParameterCreate("bh_spin_minus",
                        "TwoPunctures",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Spin of the BH on -x axis",
                        "0.0",
                        &(PRIVATE_TWOPUNCTURES_STRUCT.bh_spin_minus),
                        0,
                        NULL,
                        1,
                        "*:*", "Real");

  CCTKi_ParameterCreate("bh_spin_plus",
                        "TwoPunctures",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Spin of the BH on +x axis",
                        "0.0",
                        &(PRIVATE_TWOPUNCTURES_STRUCT.bh_spin_plus),
                        0,
                        NULL,
                        1,
                        "*:*", "Real");

  CCTKi_ParameterCreate("excis_radius",
                        "TwoPunctures",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Excision Radius, used only if fill_excision_enable=1",
                        "0.5",
                        &(PRIVATE_TWOPUNCTURES_STRUCT.excis_radius),
                        0,
                        NULL,
                        1,
                        "0:*", "Positive");

  CCTKi_ParameterCreate("fill_excision_enable",
                        "TwoPunctures",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Should I fill an excision zone? (no=0; yes=1)",
                        "0",
                        &(PRIVATE_TWOPUNCTURES_STRUCT.fill_excision_enable),
                        0,
                        NULL,
                        1,
                        "0:*", "Positive");

  CCTKi_ParameterCreate("genID_cmdline_output_enable",
                        "TwoPunctures",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Output initial data commandline?",
                        "0",
                        &(PRIVATE_TWOPUNCTURES_STRUCT.genID_cmdline_output_enable),
                        0,
                        NULL,
                        1,
                        "0:1", "Zero (no) or One (yes)");

  CCTKi_ParameterCreate("half_binary_separation",
                        "TwoPunctures",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Binary separation / 2",
                        "0.0",
                        &(PRIVATE_TWOPUNCTURES_STRUCT.half_binary_separation),
                        0,
                        NULL,
                        1,
                        "0:*", "Positive");

  CCTKi_ParameterCreate("moncrief_gw_num_radii",
                        "TwoPunctures",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "How many radii will to measure GW's",
                        "1",
                        &(PRIVATE_TWOPUNCTURES_STRUCT.moncrief_gw_num_radii),
                        0,
                        NULL,
                        1,
                        "1:10", "Positive int <= 10");

  CCTKi_ParameterCreate("moncrief_radius_GW",
                        "TwoPunctures",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_RECOVER,
                        "Radii at which to measure GW's",
                        "1.0",
                        (PRIVATE_TWOPUNCTURES_STRUCT.moncrief_radius_GW),
                        100,
                        NULL,
                        1,
                        "0:*", "zero or any positive number");

  CCTKi_ParameterCreate("x_offset",
                        "TwoPunctures",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Parameter necessary for setting up initial data",
                        "0.0",
                        &(PRIVATE_TWOPUNCTURES_STRUCT.x_offset),
                        0,
                        NULL,
                        1,
                        "*:*", "Real");

  return 0;
}

int CCTKi_BindingsTwoPuncturesParameterExtensions(void);
int CCTKi_BindingsTwoPuncturesParameterExtensions(void)
{
  return 0;
}
