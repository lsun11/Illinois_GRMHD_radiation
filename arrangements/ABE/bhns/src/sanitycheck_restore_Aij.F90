!-----------------------------------------------------------------------------
! Enforce det(g_ij)==1 and tr(A_ij)==0 to improve stability, if desired
!-----------------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bhns_sanitycheck_restore_Aij(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer, dimension(3)                    :: ext

  ext = cctk_lsh

  if(trA_detg_enforce.ne.2) then
     !Perform sanity check.  Already done if trA_detg_enforce==2!
     call sanitycheck_restore_Aij(cctkGH,ext,gxx,gxy,gxz,gyy,gyz,gzz, &
          Axx,Axy,Axz,Ayy,Ayz,Azz)
  end if


end subroutine bhns_sanitycheck_restore_Aij
