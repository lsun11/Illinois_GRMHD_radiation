 /*@@
   @file    SetSym.F
   @date    March 2000
   @author  Miguel Alcubierre
   @desc
            Sets the symmetries for the AHF grid functions.
   @enddesc
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

 /*@@
   @routine    AHFinder_SetSym
   @date       March 2000
   @author     Miguel Alcubierre
   @desc 
               Sets the symmetries for the AHF grid functions.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

      subroutine AHFinder_SetSym(CCTK_ARGUMENTS)

      implicit none
      
      DECLARE_CCTK_ARGUMENTS 
      DECLARE_CCTK_PARAMETERS

      integer sym(3)
      CCTK_INT CCTK_Equals
      integer, PARAMETER :: one = 1
      INTEGER :: ierr

!     +,+,+

      sym(1)=+1
      sym(2)=+1
      sym(3)=+1

      call SetCartSymVN(ierr,cctkGH,sym,'ahfinder::ahfgrid')
      call SetCartSymVN(ierr,cctkGH,sym,'ahfinder::ahf_exp')
      call SetCartSymVN(ierr,cctkGH,sym,'ahfinder::ahfgradn')
      call SetCartSymVN(ierr,cctkGH,sym,'ahfinder::ahfgauss')
      call SetCartSymVN(ierr,cctkGH,sym,'ahfinder::ahmask')

!     -,+,+

      sym(1)=-1
      sym(2)=+1
      sym(3)=+1

      call SetCartSymVN(ierr,cctkGH,sym,'ahfinder::ahfgradx')

!     +,-,+

      sym(1)=+1
      sym(2)=-1
      sym(3)=+1

      call SetCartSymVN(ierr,cctkGH,sym,'ahfinder::ahfgrady')

!     +,+,-

      sym(1)=+1
      sym(2)=+1
      sym(3)=-1

      call SetCartSymVN(ierr,cctkGH,sym,'ahfinder::ahfgradz')

!     End.

      end subroutine AHFinder_SetSym


 /*@@
   @file    SetSym.F
   @date    Jan 2003
   @author  Denis Pollney
   @desc
            Sets the reflection symmetry flags.
   @enddesc
 @@*/

      subroutine AHFinder_SetReflections(CCTK_ARGUMENTS)
      use AHFinder_dat

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      CCTK_REAL, parameter :: zero = 0.0D0
      CCTK_REAL, parameter :: half = 0.5D0
      CCTK_REAL, parameter :: one  = 1.0D0

!     Find out reflection symmetries depending on
!     the values of the parameters.

      if ((CCTK_Equals(ahf_octant,"yes").eq.1).or. &
           (CCTK_Equals(ahf_octant,"high").eq.1)) then

         refx = .true.
         refy = .true.
         refz = .true.

      else

         if (ahf_refx.ne.0) then
            refx = .true.
         else
            refx = .false.
         end if

         if (ahf_refy.ne.0) then
            refy = .true.
         else
            refy = .false.
         end if

         if (ahf_refz.ne.0) then
            refz = .true.
         else
            refz = .false.
         end if

      end if

      if (ahf_phi.eq.0) then
         refx = .true.
         refy = .true.
      end if

!     Force grid symmetries.  If the grid has some
!     symmetry (octant, quadrant, bitant), and the
!     horizon is centered on one of the symmetry
!     planes, then the grid symmetries are forced
!     on it regardless of the values of the symmetry
!     parameters.

      if (CCTK_Equals(domain,"octant").eq.1) then

         if (xc.eq.zero) refx = .true.
         if (yc.eq.zero) refy = .true.
         if (zc.eq.zero) refz = .true.

      else if (CCTK_Equals(domain,"quadrant").eq.1) then

         if (xc.eq.zero) refx = .true.
         if (yc.eq.zero) refy = .true.

      else if (CCTK_Equals(domain,"bitant").eq.1) then

         if (CCTK_Equals(bitant_plane,"xy").eq.1) then
            if (zc.eq.zero) refz = .true.
         else if (CCTK_Equals(bitant_plane,"xz").eq.1) then
            if (yc.eq.zero) refy = .true.
         else
            if (xc.eq.zero) refx = .true.
         end if

      end if

!     Check if we are using cartoon.  If we are, then the
!     corresponding symmetries are forced on the horizon.

      if (ahf_cartoon.ne.0) then
         cartoon = .true.
      else
         cartoon = .false.
      end if

      if (cartoon) then
         nonaxi = .false.
         refx = .true.
         refy = .true.
         !added by Zach:
         write(*,*) "CARTOON METHOD USED!"
      end if

!     If desired, write values of symmetry flags to screen.

      if (veryver) then
         write(*,*)
         write(*,*) 'Symmetries used for horizon:'
         write(*,*) 'refx = ',refx
         write(*,*) 'refy = ',refy
         write(*,*) 'refz = ',refz
         write(*,*)
      end if

      end subroutine AHFinder_SetReflections
