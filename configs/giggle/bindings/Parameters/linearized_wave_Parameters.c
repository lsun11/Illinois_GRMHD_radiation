/*@@
   @file    linearized_wave_Parameters.c
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

/* structure containing all private parameters of thorn linearized_wave */
struct
{
  CCTK_REAL Psi4imsumlw;
  CCTK_REAL Psi4resumlw;
  CCTK_REAL amplitude;
  CCTK_REAL corrector_iteration;
  CCTK_REAL time_shift;
  CCTK_REAL width;
  CCTK_INT mode;
} PRIVATE_LINEARIZED_WAVE_STRUCT;


/* structure containing all restricted parameters of thorn linearized_wave */
struct
{
  const char * bound;
  const char * initial_data;
} RESTRICTED_LINEARIZED_WAVE_STRUCT;


int CCTKi_BindingsCreatelinearized_waveParameters(void);
int CCTKi_BindingsCreatelinearized_waveParameters(void)
{
  CCTKi_ParameterCreate("bound",
                        "linearized_wave",
                        "KEYWORD",
                        "RESTRICTED",
                        CCTK_STEERABLE_NEVER,
                        "Type of boundary condition to use",
                        "none",
                        &(RESTRICTED_LINEARIZED_WAVE_STRUCT.bound),
                        0,
                        NULL,
                        6,
                        "none", "Apply no boundary condition",
                        "flat", "Flat (von Neumann, n grad phi = 0) boundary condition",
                        "static", "Static (Dirichlet, dphi/dt=0) boundary condition",
                        "radiation", "Radiation boundary condition",
                        "robin", "Robin (phi(r) = C/r) boundary condition",
                        "zero", "Zero (Dirichlet, phi=0) boundary condition");

  CCTKi_ParameterCreate("initial_data",
                        "linearized_wave",
                        "KEYWORD",
                        "RESTRICTED",
                        CCTK_STEERABLE_NEVER,
                        "Type of initial data",
                        "gaussian",
                        &(RESTRICTED_LINEARIZED_WAVE_STRUCT.initial_data),
                        0,
                        NULL,
                        4,
                        "plane", "Plane wave",
                        "gaussian", "Gaussian wave",
                        "box", "Box wave",
                        "none", "No initial data, zero phi");

  CCTKi_ParameterCreate("Psi4imsumlw",
                        "linearized_wave",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "l=2,m=2 value of Psi_4im at Radius_GW",
                        "0.D0",
                        &(PRIVATE_LINEARIZED_WAVE_STRUCT.Psi4imsumlw),
                        0,
                        NULL,
                        1,
                        "0:", "");

  CCTKi_ParameterCreate("Psi4resumlw",
                        "linearized_wave",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "l=2,m=2 value of Psi_4re at Radius_GW",
                        "0.D0",
                        &(PRIVATE_LINEARIZED_WAVE_STRUCT.Psi4resumlw),
                        0,
                        NULL,
                        1,
                        "0:", "");

  CCTKi_ParameterCreate("amplitude",
                        "linearized_wave",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "The amplitude of the waves",
                        "0.001",
                        &(PRIVATE_LINEARIZED_WAVE_STRUCT.amplitude),
                        0,
                        NULL,
                        1,
                        "*:*", "No restriction");

  CCTKi_ParameterCreate("corrector_iteration",
                        "linearized_wave",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "num of iterations in corrector loop",
                        "1.0",
                        &(PRIVATE_LINEARIZED_WAVE_STRUCT.corrector_iteration),
                        0,
                        NULL,
                        1,
                        "*:*", "No restriction");

  CCTKi_ParameterCreate("mode",
                        "linearized_wave",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "choose a mode",
                        "20",
                        &(PRIVATE_LINEARIZED_WAVE_STRUCT.mode),
                        0,
                        NULL,
                        1,
                        "20:22", "20: even 20 mode, 21: odd 21 mode, 22: even 22 mode" );

  CCTKi_ParameterCreate("time_shift",
                        "linearized_wave",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Shifting time origin for the initial data",
                        "0.0",
                        &(PRIVATE_LINEARIZED_WAVE_STRUCT.time_shift),
                        0,
                        NULL,
                        1,
                        "*:*", "No restriction");

  CCTKi_ParameterCreate("width",
                        "linearized_wave",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "The width of the wave",
                        "1.0",
                        &(PRIVATE_LINEARIZED_WAVE_STRUCT.width),
                        0,
                        NULL,
                        1,
                        "0:*", "Positive");

  return 0;
}

int CCTKi_Bindingslinearized_waveParameterExtensions(void);
int CCTKi_Bindingslinearized_waveParameterExtensions(void)
{
  return 0;
}

