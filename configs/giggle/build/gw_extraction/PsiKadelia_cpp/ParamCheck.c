 /*@@
   @file      ParamCheck.c
   @date      April 26 2002
   @author    Gabrielle Allen
   @desc 
      Check parameters for PsiKadelia
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header: /cactusdevcvs/CactusEinstein/PsiKadelia/src/ParamCheck.c,v 1.1 2002/04/28 22:00:25 allen Exp $";

CCTK_FILEVERSION(CactusEinstein_PsiKadelia_ParamCheck_c)

void PsiKadelia_ParamCheck(CCTK_ARGUMENTS);

 /*@@
   @routine    PsiKadelia_ParamCheck
   @date       April 26 2002
   @author     Gabrielle Allen
   @desc 
      Check parameters for PsiKadelia
   @enddesc 
   @calls     
   @history  
 
   @endhistory 

@@*/

void PsiKadelia_ParamCheck(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

    /*
  if(! CCTK_EQUALS(metric_type, "physical") &&
     ! CCTK_EQUALS(metric_type, "static conformal"))
  {
    CCTK_PARAMWARN("Unknown ADMBase::metric_type - known types are \"physical\" and \"static conformal\"");
  }
    */
}
