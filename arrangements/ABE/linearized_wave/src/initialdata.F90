!-----------------------------------------------------------------------------
!
! $Id: lin_wave.f90,v 1.1.1.1 2006/02/23 17:48:41 zetienne Exp $
!
!-----------------------------------------------------------------------------
!
! Find linear wave solution (See S. A. Teukolsky, PRD 26, 745 (1982) )
! Here: F = Amp * exp ( - Width * x^2 ),
!       with x = t +- r  (solution timesymmetric at t = 0)
!
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine lin_wave_initialdata(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

!
! Other variables:
!  
  integer, dimension(3)                       :: ext
  real*8                                      :: T,dT,dX,dY,dZ
  integer                            :: i,j,l

  ext = cctk_lsh

  T = cctk_time + time_shift
  dT=CCTK_DELTA_TIME
  dX=CCTK_DELTA_SPACE(1)
  dY=CCTK_DELTA_SPACE(2)
  dZ=CCTK_DELTA_SPACE(3)

  !Look inside analytic_solution.F90 for lin_wave_analytic() source code.
  call lin_wave_analytic(ext, X, Y, Z, T, amplitude, width, &
       gxx,gxy,gxz,gyy,gyz,gzz,Kxx,Kxy,Kxz,Kyy,Kyz,Kzz, &
       PhysicalRadius,RadiusDerivative, mode)

  !Look inside convert_metric_dumpanalytic.F90 for linwave_Convert() source code.
  call linwave_Convert(CCTK_PASS_FTOF)

  Psi = exp(phi)
  trK = (Kxx*gupxx + 2.d0*Kxy*gupxy + 2.d0*Kxz*gupxz + Kyy*gupyy + &
	 2.d0*Kyz*gupyz + Kzz*gupzz)/Psi**4
  Axx = Kxx/Psi**4 - gxx*trK/3.d0
  Axy = Kxy/Psi**4 - gxy*trK/3.d0
  Axz = Kxz/Psi**4 - gxz*trK/3.d0
  Ayy = Kyy/Psi**4 - gyy*trK/3.d0
  Ayz = Kyz/Psi**4 - gyz*trK/3.d0
  Azz = Kzz/Psi**4 - gzz*trK/3.d0

!Next, set BSSN matter sources to zero
  rho = 0.d0	
  S = 0.d0
  Sx = 0.d0
  Sy = 0.d0
  Sz = 0.d0
  Sxx = 0.d0
  Sxy = 0.d0
  Sxz = 0.d0
  Syy = 0.d0
  Syz = 0.d0
  Szz = 0.d0

  shiftx = 0.D0
  shifty = 0.D0
  shiftz = 0.D0

  lapm1 = 0.D0

  write(*,*) "HELLO SETTING UP LINWAVE INITIAL DATA ON GRID WITH DX=",dX
  write(*,*) "extents=",ext,X(ext(1),1,1)

end subroutine lin_wave_initialdata

