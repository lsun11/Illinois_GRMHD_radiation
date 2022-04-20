 /*@@
   @file      InitSymBound.c
   @date      Fri Nov  9 13:29:53 2001
   @author    Ian Hawke
   @desc 
   Sets the symmetries for Wave Toy - again copies the standard.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"

#include "Symmetry.h"

static const char *rcsid = "$Header: /cactusdevcvs/CactusExamples/WaveMoL/src/InitSymBound.c,v 1.1.1.1 2004/06/26 14:49:02 hawke Exp $";

CCTK_FILEVERSION(CactusExamples_WaveMoL_InitSymBound_c)

void WaveMoL_InitSymBound(CCTK_ARGUMENTS);
 /*@@
   @routine    WaveMoL_InitSymBound
   @date       Fri Nov  9 13:31:28 2001
   @author     Ian Hawke
   @desc 
   Another straight copy of a wave toy routine, however extended 
   for all the new GFs.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void WaveMoL_InitSymBound(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
      
  int sym[3];

  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;

  SetCartSymVN(cctkGH, sym,"wavemol::phi");
  SetCartSymVN(cctkGH, sym,"wavemol::phit");

  sym[0] = -1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym,"wavemol::phix");

  sym[0] = 1;
  sym[1] = -1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym,"wavemol::phiy");

  sym[0] = 1;
  sym[1] = 1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym,"wavemol::phiz");

  return;
}
