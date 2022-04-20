// Update boundaries using Cactus outer boundary methods
  
#include "cctk.h" 
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

#include "stdio.h"

extern "C" void scalarwaveMoL_update_boundary(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CartSymVN(cctkGH,"scalarwaveMoL::phi");
  CartSymVN(cctkGH,"scalarwaveMoL::phidot");
  /*
  int index = CCTK_GFINDEX3D(cctkGH,10,10,10);
  int varindex = CCTK_VarIndex("scalarwaveMoL::phi");

  double dataold = cctkGH->data[varindex][1][index];
  double datanew = cctkGH->data[varindex][0][index];
  */
  //  printf("HELLO old = %.15e, new = %.15e\n",phi_p[CCTK_GFINDEX3D(cctkGH,10,10,10)],phi[CCTK_GFINDEX3D(cctkGH,10,10,10)]);

  int ierr=-1;

  /* Uses all default arguments, so invalid table handle -1 can be passed */
  ierr = Boundary_SelectGroupForBC
    (cctkGH, CCTK_ALL_FACES, 1, -1, "scalarwaveMoL::scalarMoLevolve", bound);

  if (ierr < 0) 
  {
    CCTK_WARN(0,"Boundary conditions not applied - giving up!");
  }

  return;
}
