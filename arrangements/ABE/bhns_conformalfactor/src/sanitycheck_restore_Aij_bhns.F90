!-----------------------------------------------------------------------------
! Enforce det(g_ij)==1 and tr(A_ij)==0 to improve stability, if desired
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
subroutine magnetar_sanitycheck_restore_Aij_bhns(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer, dimension(3)                    :: ext
  real*8                                   :: dT,dX,dY,dZ,xmax
  real*8                                   :: Kmin,Kmax,detmin,detmax
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  ext = cctk_lsh
  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  if(trA_detg_enforce.ne.2) then
     !Perform sanity check.  Already done if trA_detg_enforce==2!
     call sanitycheck_restore_Aij(cctkGH,cctk_lsh,gxx,gxy,gxz,gyy,gyz,gzz, &
          Axx,Axy,Axz,Ayy,Ayz,Azz)
  end if
  !following line needed when we restart from checkpoint with excision enabled...
!When excision disabled, causes discrepancy.  This might indicate a bug in the DAGH version..... (and by extension this version of the code)
!  if(Symmetry .eq. AXISYM .and. excision_enable==1) then
!     call CCTK_VarIndex(index,'BSSN::Axx')
!     call BndCartoon2DVI(dummy, cctkGH, 2, index)
!  end if

end subroutine magnetar_sanitycheck_restore_Aij_bhns
