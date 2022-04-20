/*@@
  @file      AHFinder.F
  @date      April 1998
  @author    Miguel Alcubierre
  @desc
             Find an apparent horizon using a spherical harmonic
             expansion and either a minimization or Gundlach
             fast flow algorithm.

             Cactus4 by Lars Nerger / Gerd Lanfermann
  @enddesc
  @version   $Header: /cactus/CactusEinstein/AHFinder/src/AHFinder.F,v 1.85 2004/02/07 10:48:21 schnetter Exp $
@@*/

#include "cctk.h"

#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine AHFinder_bbh(CCTK_ARGUMENTS)
  
  use AHFinder_dat
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  REAL*8  :: xbh1tmp,ybh1tmp,zbh1tmp
  INTEGER :: ahflmaxtmp

  if(ahf_lmax .lt. 4) then
     write(*,*) "BAD PIE!  YOU MUST SET ahf_lmax to at least 4 WITH PARAMTER ahf_bbh TURNED ON!"
     write(*,*) "This is needed since we compute 4 l-modes for each individual horizon,"
     write(*,*) "and out_c0, out_cc variables are of length ~ahf_lmax (avoids memory error)"
     write(*,*) ""
     stop
  end if

  xbh1tmp = xbh1
  ybh1tmp = ybh1
  zbh1tmp = zbh1

  !First search for common horizon!

  xbh1 = 0.D0
  ybh1 = 0.D0
  zbh1 = 0.D0

  ahflmaxtmp = ahf_lmax
  call AHFinder(CCTK_PASS_FTOF)

  !Next search for horizon centered on BH2
  xbh1 = xbh2
  ybh1 = ybh2
  zbh1 = zbh2

!  ahf_lmax = 4
  call AHFinder(CCTK_PASS_FTOF)

  xbh1 = xbh1tmp
  ybh1 = ybh1tmp
  zbh1 = zbh1tmp

  !Finally, search for horizon centered on BH2
  call AHFinder(CCTK_PASS_FTOF)

  ahf_lmax = ahflmaxtmp

end subroutine AHFinder_bbh
