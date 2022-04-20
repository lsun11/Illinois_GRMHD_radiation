/*@@
   @file    gw_extraction_Parameters.c
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
#include "ParameterCRestrictedIO.h"

/* structure containing all private parameters of thorn gw_extraction */
struct
{
  CCTK_REAL compute_Psi4_max_radius;
  CCTK_REAL compute_Psi4_min_radius;
  CCTK_REAL radius_GW_Psi4[101];
  const char * psif_vec;
  CCTK_INT enable_interp_onepoint;
  CCTK_INT nmodes_Psi4;
  CCTK_INT nmodes_ZM;
  CCTK_INT num_extraction_radii;
  CCTK_INT radius_power;
  CCTK_INT scale_with_radius;
  CCTK_INT use_Rij_from_compute_ricci;
} PRIVATE_GW_EXTRACTION_STRUCT;


/* structure containing all restricted parameters of thorn gw_extraction */
struct
{
  CCTK_REAL dR_GW;
  CCTK_REAL ddR_GW;
  CCTK_REAL phi_GW;
  CCTK_REAL radius_GW;
  CCTK_REAL radius_GW_phys;
  CCTK_REAL theta_GW;
} RESTRICTED_GW_EXTRACTION_STRUCT;


int CCTKi_BindingsCreategw_extractionParameters(void);
int CCTKi_BindingsCreategw_extractionParameters(void)
{
  CCTKi_ParameterCreate("dR_GW",
                        "gw_extraction",
                        "REAL",
                        "RESTRICTED",
                        CCTK_STEERABLE_NEVER,
                        "dR/dr",
                        "-1.0",
                        &(RESTRICTED_GW_EXTRACTION_STRUCT.dR_GW),
                        0,
                        NULL,
                        1,
                        "*:*", "Anything, negative means set = to surf_radius converted to phys coords!");

  CCTKi_ParameterCreate("ddR_GW",
                        "gw_extraction",
                        "REAL",
                        "RESTRICTED",
                        CCTK_STEERABLE_NEVER,
                        "d^2R/dr^2",
                        "-1.0",
                        &(RESTRICTED_GW_EXTRACTION_STRUCT.ddR_GW),
                        0,
                        NULL,
                        1,
                        "*:*", "Anything, negative means set = to surf_radius converted to phys coords!");

  CCTKi_ParameterCreate("phi_GW",
                        "gw_extraction",
                        "REAL",
                        "RESTRICTED",
                        CCTK_STEERABLE_ALWAYS,
                        "The angle phi in radians in which GW's are measured",
                        "0.52",
                        &(RESTRICTED_GW_EXTRACTION_STRUCT.phi_GW),
                        0,
                        NULL,
                        1,
                        "0:*", "Positive");

  CCTKi_ParameterCreate("radius_GW",
                        "gw_extraction",
                        "REAL",
                        "RESTRICTED",
                        CCTK_STEERABLE_ALWAYS,
                        "The distance from the origin where GW's are measured",
                        "-1.0",
                        &(RESTRICTED_GW_EXTRACTION_STRUCT.radius_GW),
                        0,
                        NULL,
                        1,
                        "*:*", "Anything, negative means set = to surf_radius!");

  CCTKi_ParameterCreate("radius_GW_phys",
                        "gw_extraction",
                        "REAL",
                        "RESTRICTED",
                        CCTK_STEERABLE_ALWAYS,
                        "The distance from the origin where GW's are measured, PHYSICAL COORDINATES!",
                        "-1.0",
                        &(RESTRICTED_GW_EXTRACTION_STRUCT.radius_GW_phys),
                        0,
                        NULL,
                        1,
                        "*:*", "Anything, negative means set = to surf_radius converted to phys coords!");

  CCTKi_ParameterCreate("theta_GW",
                        "gw_extraction",
                        "REAL",
                        "RESTRICTED",
                        CCTK_STEERABLE_ALWAYS,
                        "The angle theta in radians in which GW's are measured",
                        "0.79",
                        &(RESTRICTED_GW_EXTRACTION_STRUCT.theta_GW),
                        0,
                        NULL,
                        1,
                        "0:*", "Positive");

  CCTKi_ParameterCreate("compute_Psi4_max_radius",
                        "gw_extraction",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Maximum radius at which to compute Psi4",
                        "100000.0",
                        &(PRIVATE_GW_EXTRACTION_STRUCT.compute_Psi4_max_radius),
                        0,
                        NULL,
                        1,
                        "0:*", "Zero or positive");

  CCTKi_ParameterCreate("compute_Psi4_min_radius",
                        "gw_extraction",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Minimum radius at which to compute Psi4",
                        "0.0",
                        &(PRIVATE_GW_EXTRACTION_STRUCT.compute_Psi4_min_radius),
                        0,
                        NULL,
                        1,
                        "0:*", "Zero or positive");

  CCTKi_ParameterCreate("enable_interp_onepoint",
                        "gw_extraction",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Interpolate Psi4re and Psi4im to one point on the grid?  Useful for linearized wave; uses radius_GW,theta_GW,phi_GW",
                        "0",
                        &(PRIVATE_GW_EXTRACTION_STRUCT.enable_interp_onepoint),
                        0,
                        NULL,
                        1,
                        "0:*", "Positive");

  CCTKi_ParameterCreate("nmodes_Psi4",
                        "gw_extraction",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "number of modes for GW extraction routine",
                        "21",
                        &(PRIVATE_GW_EXTRACTION_STRUCT.nmodes_Psi4),
                        0,
                        NULL,
                        1,
                        "0:*", "Number of spin-weight -2 l,m modes to output. Set this to 21 for all l<=3 modes.");

  CCTKi_ParameterCreate("nmodes_ZM",
                        "gw_extraction",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "number of modes for GW extraction routine",
                        "7",
                        &(PRIVATE_GW_EXTRACTION_STRUCT.nmodes_ZM),
                        0,
                        NULL,
                        1,
                        "0:*", "Number of spin-weight -2 l,m modes to output. Set this to 7 for all l<=3 modes.");

  CCTKi_ParameterCreate("num_extraction_radii",
                        "gw_extraction",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Number of extraction radii",
                        "1",
                        &(PRIVATE_GW_EXTRACTION_STRUCT.num_extraction_radii),
                        0,
                        NULL,
                        1,
                        "0:*", "Positive");

  CCTKi_ParameterCreate("psif_vec",
                        "gw_extraction",
                        "KEYWORD",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Where does the basis vector point?",
                        "radial",
                        &(PRIVATE_GW_EXTRACTION_STRUCT.psif_vec),
                        0,
                        NULL,
                        4,
                        "radial", "",
                        "cartesian", "",
                        "metric_diag", "",
                        "shock", "");

  CCTKi_ParameterCreate("radius_GW_Psi4",
                        "gw_extraction",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "extraction radii",
                        "0",
                        (PRIVATE_GW_EXTRACTION_STRUCT.radius_GW_Psi4),
                        101,
                        NULL,
                        1,
                        "0:*", "Zero or positive");

  CCTKi_ParameterCreate("radius_power",
                        "gw_extraction",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Power of radial coordinate to multiply by",
                        "1",
                        &(PRIVATE_GW_EXTRACTION_STRUCT.radius_power),
                        0,
                        NULL,
                        1,
                        "1:*", "");

  CCTKi_ParameterCreate("scale_with_radius",
                        "gw_extraction",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Multiply all scalars by power of radial coordinate",
                        "no",
                        &(PRIVATE_GW_EXTRACTION_STRUCT.scale_with_radius),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("use_Rij_from_compute_ricci",
                        "gw_extraction",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Compute Ricci using the built-in PsiKadelia routine instead of our own Ricci computation?  WARNING: Do not yet trust results when this is set to 1.",
                        "0",
                        &(PRIVATE_GW_EXTRACTION_STRUCT.use_Rij_from_compute_ricci),
                        0,
                        NULL,
                        1,
                        "0:*", "Positive");

  return 0;
}

int CCTKi_Bindingsgw_extractionParameterExtensions(void);
int CCTKi_Bindingsgw_extractionParameterExtensions(void)
{
  return 0;
}

