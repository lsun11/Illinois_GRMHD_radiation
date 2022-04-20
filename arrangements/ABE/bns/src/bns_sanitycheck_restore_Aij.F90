!-----------------------------------------------------------------------------
!
!$Id: predict_mhd.F90  $
!
!-----------------------------------------------------------------------------
!
! Predict & Correct RHS routines for mhd variables
!
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
!
! Predictor
!
!-----------------------------------------------------------------------------
subroutine bns_sanitycheck_restore_Aij(CCTK_ARGUMENTS)
  
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
     call sanitycheck_restore_Aij(ext,gxx,gxy,gxz,gyy,gyz,gzz, &
       Axx,Axy,Axz,Ayy,Ayz,Azz,detmin,detmax,Kmin,Kmax)
  end if

end subroutine bns_sanitycheck_restore_Aij
