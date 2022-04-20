!-----------------------------------------------------------------------------
!
!$Id: Updaters.F90,v 1.2 2006/03/20 10:32:24 ytliu Exp $
!
!-----------------------------------------------------------------------------
!
! Fortran evolution schemes for interior (Derivs subroutines)
!
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!

!-----------------------------------------------------------------------------
!
! Compute time derivative of lapse and shift (might change it later)
!
!-----------------------------------------------------------------------------
subroutine BSSN_Gauge_Derivs(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  integer, dimension(3)                    :: ext
  real*8					:: dT
!
  dT = CCTK_DELTA_TIME
  ext(1) = cctk_lsh(1)
  ext(2) = cctk_lsh(2)
  ext(3) = cctk_lsh(3)

  ! Zach says: I've commented out the following function call for efficiency, since 
  !            there are currently no shift conditions that support it.
  if(1==0) then
     call gauge_derivs(ext,dT, lapm1_p,lapm1, lapset_timederiv, &
          shiftx_p,shiftx,shiftxt_timederiv, &
          shifty_p,shifty,shiftyt_timederiv, &
          shiftz_p,shiftz,shiftzt_timederiv)
  end if
end subroutine BSSN_Gauge_Derivs

subroutine Derivs(ex,X,Y,Z,dX,dY,dZ,f,fx,fy,fz,Symmetry)
  implicit none
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))	      :: f,X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: fx,fy,fz  
  real*8				      :: dX,dY,dZ, SYM
  integer				      :: Symmetry
  integer				      :: NO_SYMM,EQUATORIAL
  integer                                     :: OCTANT,PI_SYMM,AXISYM
  parameter ( NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter (SYM = 1.d0)
!
  if (Symmetry==OCTANT) then
     call gderivs_oct(ex,f,fx,fy,fz,dX,dY,dZ,SYM,SYM,SYM)
  elseif (Symmetry==EQUATORIAL) then
     call gderivs_eq(ex,f,fx,fy,fz,dX,dY,dZ,SYM,SYM,SYM)
  elseif (Symmetry==AXISYM) then
     call gderivs_axi(ex,f,fx,fy,fz,dX,dY,dZ,SYM,SYM,SYM)
  elseif (Symmetry==NO_SYMM) then 
     call gderivs_gen(ex,f,fx,fy,fz,dX,dY,dZ)
  else 
    write(*,*) 'Symmetry type not supported in Derivs!'
    stop
  end if
end subroutine Derivs
