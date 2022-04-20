/*@@
  @file      AHFinder_find3.F
  @date      March 1999
  @author    Lars Nerger
  @desc 
             Muliplicate gridfunctions for final grid function
  @enddesc 
@@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

!     To output single gridfunctions when searching for 3 horizons
!     the gridfunctions are here multiplied.  The result of this is
!     that the functions ahfgrid3 and ahf_exp3 will have zeroes at
!     the location of all 3 horizons.

      subroutine AHFinder_find3(CCTK_ARGUMENTS,mtype)

      use AHFinder_dat

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

      integer mtype

      CCTK_REAL one

!     Description of variables
!
!     mtype     Type of surface: see AHFinder.F


!     **********************
!     ***   DEFINE one   ***
!     **********************

      one  = 1.0D0

!     *************************
!     ***   START ROUTINE   ***
!     *************************

!     For first outer horizon initialize {ahfgrid3,ahf_exp3}.

      if ( mfind .eq. 1) then
         if ( mtype .eq. 1 ) then

            ahfgrid3 = ahfgrid
            ahf_exp3 = ahf_exp

         else

!           If we did not find an outer horizon, make {ahfgrid3,ahf_exp3}
!           equal to one.

            ahfgrid3 = one
            ahf_exp3 = one

         end if

!     For second and third horizons, if we found an outer
!     horizon multiply the old values of {ahfgrid3,ahf_exp3}
!     with the new {ahfgrid3,ahf_exp3}.

      else if (mtype.eq.1) then

         ahfgrid3 = ahfgrid*ahfgrid3
         ahf_exp3 = ahf_exp*ahf_exp3

      end if


!     ***************
!     ***   END   ***
!     ***************

      end subroutine AHFinder_find3

