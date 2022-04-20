/*
  Set the symmetries for the BSSN variables
*/

#include <stdio.h>
#include <stdlib.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

#include "GenericFD.h"

static char *rcsid="$Header: /peter/piper/picked/a/peck/of/pickled/whatever $";

CCTK_FILEVERSION(BSSN_InitSymBound)

  extern "C" void BSSN_InitSymBound(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

#ifdef FD_C6
    /*
    if(cctk_nghostzones[1]<4 && cowling_enable==0) {
      printf("ERROR: LOOKS LIKE YOU'RE TRYING TO DO A 6TH ORDER EVOLUTION WITH %d GHOSTZONES!\n",cctk_nghostzones[1]);
      printf("Please set the number of ghostzones to >= 4, or else move to lower order!\n");
      exit(1);
    }
    */
#endif

#ifdef FD_C4
    if(cctk_nghostzones[1]<3 && cowling_enable==0) {
      printf("ERROR: LOOKS LIKE YOU'RE TRYING TO DO A 4TH ORDER EVOLUTION WITH %d GHOSTZONES!\n",cctk_nghostzones[1]);
      printf("Please set the number of ghostzones to >= 3, or else move to 2nd order!\n");
      exit(1);
    }
#endif

#ifdef FD_C2
    if(cctk_nghostzones[1]>1) {
      // MAKE THE WARNING VISIBLE!
      for(int i=0;i<100;i++) {
	CCTK_WARN (1, "================================================================================================");
	printf("WARNING: YOU'RE DOING A 2ND ORDER EVOLUTION WITH %d GHOSTZONES, WHICH IS MORE THAN NEEDED!\n",cctk_nghostzones[1]);
	CCTK_WARN (1, "================================================================================================");
      }
    }
#endif

  int sym[3];
  printf("SYMMETRY = %d\n",Symmetry);
  // First octant symmetry case:
  if(Symmetry==2) {
    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::phi");
    SetCartSymVN(cctkGH, sym,"BSSN::chi");
    SetCartSymVN(cctkGH, sym,"BSSN::trK");

    SetCartSymVN(cctkGH, sym,"BSSN::gxx");
    SetCartSymVN(cctkGH, sym,"BSSN::gyy");
    SetCartSymVN(cctkGH, sym,"BSSN::gzz");
    sym[0] = -1; sym[1] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::gxy");
    sym[0] = -1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::gxz");
    sym[0]= 1; sym[1] = -1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::gyz");

    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::psi");
    SetCartSymVN(cctkGH, sym,"BSSN::Kxx");
    SetCartSymVN(cctkGH, sym,"BSSN::Kyy");
    SetCartSymVN(cctkGH, sym,"BSSN::Kzz");
    sym[0] = -1; sym[1] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Kxy");
    sym[0] = -1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Kxz");
    sym[0]= 1; sym[1] = -1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Kyz");


    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Axx");
    SetCartSymVN(cctkGH, sym,"BSSN::Ayy");
    SetCartSymVN(cctkGH, sym,"BSSN::Azz");
    sym[0] = -1; sym[1] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Axy");
    sym[0] = -1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Axz");
    sym[0]= 1; sym[1] = -1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Ayz");

    sym[0] = -1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Gammax");
    sym[0] = 1; sym[1] = -1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Gammay");
    sym[0] = 1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Gammaz");

    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::gupxx");
    SetCartSymVN(cctkGH, sym,"BSSN::gupyy");
    SetCartSymVN(cctkGH, sym,"BSSN::gupzz");
    sym[0] = -1; sym[1] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::gupxy");
    sym[0] = -1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::gupxz");
    sym[0]= 1; sym[1] = -1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::gupyz");

    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::rho");
    SetCartSymVN(cctkGH, sym,"BSSN::S");

    sym[0] = -1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Sx");
    sym[0] = 1; sym[1] = -1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Sy");
    sym[0] = 1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Sz");

    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Sxx");
    SetCartSymVN(cctkGH, sym,"BSSN::Syy");
    SetCartSymVN(cctkGH, sym,"BSSN::Szz");
    sym[0] = -1; sym[1] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Sxy");
    sym[0] = -1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Sxz");
    sym[0]= 1; sym[1] = -1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Syz");

  }

  // Next, axisymmetry case:
  if (Symmetry==4) {   
    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::phi");
    SetCartSymVN(cctkGH, sym,"BSSN::chi");
    SetCartSymVN(cctkGH, sym,"BSSN::trK");

    SetCartSymVN(cctkGH, sym,"BSSN::gxx");
    SetCartSymVN(cctkGH, sym,"BSSN::gyy");
    SetCartSymVN(cctkGH, sym,"BSSN::gzz");
    sym[0] = 1; sym[1] = -1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::gxy");
    sym[0] = -1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::gxz");
    sym[0]= -1; sym[1] = -1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::gyz");

    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::psi");

    SetCartSymVN(cctkGH, sym,"BSSN::Kxx");
    SetCartSymVN(cctkGH, sym,"BSSN::Kyy");
    SetCartSymVN(cctkGH, sym,"BSSN::Kzz");
    sym[0] = 1; sym[1] = -1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Kxy");
    sym[0] = -1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Kxz");
    sym[0]= -1; sym[1] = -1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Kyz");


    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Axx");
    SetCartSymVN(cctkGH, sym,"BSSN::Ayy");
    SetCartSymVN(cctkGH, sym,"BSSN::Azz");
    sym[0] = 1; sym[1] = -1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Axy");
    sym[0] = -1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Axz");
    sym[0]= -1; sym[1] = -1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Ayz");

    sym[0] = -1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Gammax");
    sym[0] = -1; sym[1] =-1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Gammay");
    sym[0] = 1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Gammaz");

    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::gupxx");
    SetCartSymVN(cctkGH, sym,"BSSN::gupyy");
    SetCartSymVN(cctkGH, sym,"BSSN::gupzz");
    sym[0] = 1; sym[1] = -1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::gupxy");
    sym[0] = -1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::gupxz");
    sym[0]= -1; sym[1] = -1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::gupyz");

    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::rho");
    SetCartSymVN(cctkGH, sym,"BSSN::S");

    sym[0] = -1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Sx");
    sym[0] = -1; sym[1] =-1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Sy");
    sym[0] = 1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Sz");

    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Sxx");
    SetCartSymVN(cctkGH, sym,"BSSN::Syy");
    SetCartSymVN(cctkGH, sym,"BSSN::Szz");
    sym[0] = 1; sym[1] = -1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Sxy");
    sym[0] = -1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Sxz");
    sym[0]= -1; sym[1] = -1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Syz");

  }
  // Next equatorial symmetry case:
  if(Symmetry==1) {
    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::phi");
    SetCartSymVN(cctkGH, sym,"BSSN::chi");
    SetCartSymVN(cctkGH, sym,"BSSN::trK");

    SetCartSymVN(cctkGH, sym,"BSSN::gxx");
    SetCartSymVN(cctkGH, sym,"BSSN::gyy");
    SetCartSymVN(cctkGH, sym,"BSSN::gzz");
    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::gxy");
    sym[0] = 1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::gxz");
    sym[0] = 1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::gyz");

    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::psi");
    SetCartSymVN(cctkGH, sym,"BSSN::Kxx");
    SetCartSymVN(cctkGH, sym,"BSSN::Kyy");
    SetCartSymVN(cctkGH, sym,"BSSN::Kzz");
    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Kxy");
    sym[0] = 1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Kxz");
    sym[0] = 1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Kyz");


    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Axx");
    SetCartSymVN(cctkGH, sym,"BSSN::Ayy");
    SetCartSymVN(cctkGH, sym,"BSSN::Azz");
    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Axy");
    sym[0] = 1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Axz");
    sym[0] = 1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Ayz");

    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Gammax");
    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Gammay");
    sym[0] = 1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Gammaz");

    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::gupxx");
    SetCartSymVN(cctkGH, sym,"BSSN::gupyy");
    SetCartSymVN(cctkGH, sym,"BSSN::gupzz");
    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::gupxy");
    sym[0] = 1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::gupxz");
    sym[0] = 1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::gupyz");

    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::rho");
    SetCartSymVN(cctkGH, sym,"BSSN::S");

    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Sx");
    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Sy");
    sym[0] = 1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Sz");

    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Sxx");
    SetCartSymVN(cctkGH, sym,"BSSN::Syy");
    SetCartSymVN(cctkGH, sym,"BSSN::Szz");
    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::Sxy");
    sym[0] = 1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Sxz");
    sym[0] = 1; sym[1] = 1; sym[2] = -1;
    SetCartSymVN(cctkGH, sym,"BSSN::Syz");

  }

  // Finally no symmetry case:
  if(Symmetry==0) {
    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymVN(cctkGH, sym,"BSSN::phi");
    SetCartSymVN(cctkGH, sym,"BSSN::chi");
    SetCartSymVN(cctkGH, sym,"BSSN::trK");

    SetCartSymVN(cctkGH, sym,"BSSN::gxx");
    SetCartSymVN(cctkGH, sym,"BSSN::gyy");
    SetCartSymVN(cctkGH, sym,"BSSN::gzz");
    SetCartSymVN(cctkGH, sym,"BSSN::gxy");
    SetCartSymVN(cctkGH, sym,"BSSN::gxz");
    SetCartSymVN(cctkGH, sym,"BSSN::gyz");

    SetCartSymVN(cctkGH, sym,"BSSN::psi");
    SetCartSymVN(cctkGH, sym,"BSSN::Kxx");
    SetCartSymVN(cctkGH, sym,"BSSN::Kyy");
    SetCartSymVN(cctkGH, sym,"BSSN::Kzz");
    SetCartSymVN(cctkGH, sym,"BSSN::Kxy");
    SetCartSymVN(cctkGH, sym,"BSSN::Kxz");
    SetCartSymVN(cctkGH, sym,"BSSN::Kyz");


    SetCartSymVN(cctkGH, sym,"BSSN::Axx");
    SetCartSymVN(cctkGH, sym,"BSSN::Ayy");
    SetCartSymVN(cctkGH, sym,"BSSN::Azz");
    SetCartSymVN(cctkGH, sym,"BSSN::Axy");
    SetCartSymVN(cctkGH, sym,"BSSN::Axz");
    SetCartSymVN(cctkGH, sym,"BSSN::Ayz");

    SetCartSymVN(cctkGH, sym,"BSSN::Gammax");
    SetCartSymVN(cctkGH, sym,"BSSN::Gammay");
    SetCartSymVN(cctkGH, sym,"BSSN::Gammaz");

    SetCartSymVN(cctkGH, sym,"BSSN::gupxx");
    SetCartSymVN(cctkGH, sym,"BSSN::gupyy");
    SetCartSymVN(cctkGH, sym,"BSSN::gupzz");
    SetCartSymVN(cctkGH, sym,"BSSN::gupxy");
    SetCartSymVN(cctkGH, sym,"BSSN::gupxz");
    SetCartSymVN(cctkGH, sym,"BSSN::gupyz");

    SetCartSymVN(cctkGH, sym,"BSSN::rho");
    SetCartSymVN(cctkGH, sym,"BSSN::S");

    SetCartSymVN(cctkGH, sym,"BSSN::Sx");
    SetCartSymVN(cctkGH, sym,"BSSN::Sy");
    SetCartSymVN(cctkGH, sym,"BSSN::Sz");

    SetCartSymVN(cctkGH, sym,"BSSN::Sxx");
    SetCartSymVN(cctkGH, sym,"BSSN::Syy");
    SetCartSymVN(cctkGH, sym,"BSSN::Szz");
    SetCartSymVN(cctkGH, sym,"BSSN::Sxy");
    SetCartSymVN(cctkGH, sym,"BSSN::Sxz");
    SetCartSymVN(cctkGH, sym,"BSSN::Syz");

  }

  return;
}
