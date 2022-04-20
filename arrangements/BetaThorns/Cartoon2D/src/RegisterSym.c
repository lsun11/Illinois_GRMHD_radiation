/*@@
   @file      RegisterSym.c
   @date      February 2004
   @author    Erik Schnetter
   @desc 
              Register the symmetry boundary faces.
   @enddesc
 @@*/

#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "Cartoon2D.h"

static const char * const rcsid = "$Header: /cactusdevcvs/BetaThorns/Cartoon2D/src/RegisterSym.c,v 1.2 2004/05/28 15:59:40 schnetter Exp $";

CCTK_FILEVERSION(BetaThorns_Cartoon2D_SetSym_c);



void Cartoon2D_RegisterSymmetries (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  
  int f;
  CCTK_INT handle;
  CCTK_INT faces[6];
  CCTK_INT width[6];
  CCTK_INT ierr;
  
  faces[0] = 1;
  faces[1] = 0;
  faces[2] = 1;
  faces[3] = 1;
  faces[4] = 0;
  faces[5] = 0;
  
  for (f=0; f<6; ++f) {
    width[f] = cctk_nghostzones[f/2];
  }
  
  handle = SymmetryRegister ("cartoon");
  if (handle < 0) {
    CCTK_WARN (0, "Could not register Cartoon boundary condition");
  }
  
  ierr = SymmetryRegisterGrid (cctkGH, handle, faces, width);
  if (ierr < 0) {
    CCTK_WARN (0, "Could not register the Cartoon boundaries -- probably some other thorn has already registered the same boundary faces for a different symmetry");
  }
  
  ierr = SymmetryRegisterGridInterpolator
    (cctkGH, handle, Cartoon2D_SymmetryInterpolate);
  if (ierr < 0) {
    CCTK_WARN (0, "Could not register the symmetry interpolator");
  }
  
}
