/*
   @file      CheckParameters.c
   @date      
   @author    Gabrielle Allen
   @desc 
              Check parameters for the wave equation initial data
   @enddesc 
   @version $Header: /cactusdevcvs/CactusExamples/IDWaveMoL/src/CheckParameters.c,v 1.1.1.1 2004/06/26 14:49:20 hawke Exp $
 @@*/

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

static const char *rcsid = "$Header: /cactusdevcvs/CactusExamples/IDWaveMoL/src/CheckParameters.c,v 1.1.1.1 2004/06/26 14:49:20 hawke Exp $";

CCTK_FILEVERSION(CactusExamples_IDWaveMoL_CheckParameters_c)

void IDWaveMoL_CheckParameters(CCTK_ARGUMENTS);

 /*@@
   @routine    IDScalarWaveC_CheckParameters
   @date       
   @author     Gabrielle Allen
   @desc 
               Check parameters for the wave equation initial data
   @enddesc 
   @calls      
   @calledby   
   @history 
   @hdate Mon Oct 11 11:49:21 1999 @hauthor Tom Goodale
   @hdesc Converted to C++ 
   @hdate Thu Feb 17 09:20:41 2000 @hauthor Tom Goodale
   @hdesc Converted to C
   @endhistory 

@@*/

void IDWaveMoL_CheckParameters(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

}
