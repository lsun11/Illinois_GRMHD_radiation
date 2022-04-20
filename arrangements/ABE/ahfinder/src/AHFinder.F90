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

      subroutine AHFinder(CCTK_ARGUMENTS)

      use AHFinder_dat

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      CCTK_INT, parameter :: izero = 0
      integer,  parameter :: ik = kind(izero)

      logical status,horizon,gauss,ahfconformal
      logical giveup,report

      integer handle,x_index,y_index,z_index,ierror
      integer ahfso,ahfafter
      integer i,j,l,m
      integer NN
      integer mtype
      integer getoutpfx
      integer istat
      integer nhorizon

      CCTK_REAL atime
      CCTK_REAL rmx,r0
      CCTK_REAL intexp_h,intexp2_h,intarea_h
      CCTK_REAL aux,expin,expout,small1,small2
      CCTK_REAL c0_temp(0:18)
      CCTK_REAL maxval(3),minval(3)
      CCTK_REAL rhor

      CCTK_REAL, parameter :: zero = 0.0D0
      CCTK_REAL, parameter :: half = 0.5D0
      CCTK_REAL, parameter :: one  = 1.0D0

      CCTK_REAL, dimension (1:3) :: ahfmass
      CCTK_INT :: save_conformal_state

      character(200) logf

      save ahfafter,ahfso,atime

!     Description of variables:
!
!     status    Surface found?
!
!     horizon   Horizon found?
!
!     gauss     Output Gaussian curvature?
!
!     ahfso     How often we want to find a horizon?
!     ahfafter  After how many iterations start looking for horizons?
!
!     i,l,m     Counters.
!
!     NN        Total number of terms in expansion.
!
!     rmx       Radius of largest sphere that fits in the grid.
!
!     mtype     Type of surface:
!                  0.  No surface found.
!                  1.  Integral of H^2 is positive out and negative in.
!                  2.  Integral of H^2 is negative out and positive in.
!                  3.  Integral of H^2 is negative out and in.
!                  4.  Integral of H^2 is positive out and in.
!
!     nhorizon  Number of horizons to look for.


!     ***************************
!     ***   INITIALIZE MASK   ***
!     ***************************

!     Initialize mask to one, but only on first iteration.
!     Notice that I do this even if we do not want to find
!     a horizon from the very beginning.  This is to make
!     sure that the mask is 1 and not 0 in case some other
!     thorn wants to use it.

      if (cctk_iteration.eq.0 .or. CCTK_EQUALS(ahf_mask, "off")) then
         ahmask = one
      end if


!     ***********************************************************
!     ***   CHECK IF WE WANT TO FIND A HORIZON AT THIS TIME   ***
!     ***********************************************************

      atime = ahf_findaftertime

      if (atime.le.zero) then
         ahfafter = ahf_findafter
         if (cctk_iteration.lt.ahfafter) return
      else
         if (cctk_time.lt.atime) return
         if (ahf_ncall.eq.0) ahfafter = cctk_iteration
      end if

      ahfso = ahf_findevery
      if (mod(int(cctk_iteration-ahfafter),ahfso).ne.0) return


!     *******************************************
!     ***   ADD TO COUNT OF CALLS TO FINDER   ***
!     *******************************************

!     Add to count of calls to finder.  Notice that this counter was
!     initialized to zero, so it must be increased on first call to 1.

      ahf_ncall = ahf_ncall + 1


!     *********************************************************
!     ***   FIND OUT IF WE WANT TO FIND A TRAPPED SURFACE   ***
!     *********************************************************

      if (ahf_trapped_surface.ne.0) then
         find_trapped_surface = .true.
      else
         find_trapped_surface = .false.
      end if


!     ***************************
!     ***   GRID PARAMETERS   ***
!     ***************************

      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

      dx = cctk_delta_space(1)
      dy = cctk_delta_space(2)
      dz = cctk_delta_space(3)


!     *********************************************
!     ***   SET SYMMETRIES FOR GRID FUNCTIONS   ***
!     *********************************************

      call AHFinder_SetSym(CCTK_ARGUMENTS)


!     **********************************************
!     ***   FIND OUT ON WHICH PROCESSOR WE ARE   ***
!     **********************************************

      myproc = CCTK_MyProc(cctkGH)


!     *******************************************
!     ***   FIND  {xmn,xmx,ymn,ymx,zmn,zmx}   ***
!     *******************************************

      call CCTK_CoordRange(ierr,cctkGH,xmn,xmx,-1,"x","cart3d")
      call CCTK_CoordRange(ierr,cctkGH,ymn,ymx,-1,"y","cart3d")
      call CCTK_CoordRange(ierr,cctkGH,zmn,zmx,-1,"z","cart3d")


!     **************************************************
!     ***   TRANSFORM CONFORMAL TO PHYSICAL METRIC   ***
!     **************************************************

      ahfconformal = .false.

!     The horizon finding routines know nothing about conformal
!     metrics, so we need to make sure that we are working with
!     the physical metric.  I change not only the metric, but
!     also the "conformal state" flag to make sure that any
!     macros I use have always the physical metric.  However,
!     I need to remember the original value of the flag to change
!     it back at the end if necesary.

      if ( CCTK_Equals(metric_type,"static conformal") == 1 ) then
         ahfconformal    = .true.
         save_conformal_state = conformal_state
         conformal_state = 0
         
         call ConfToPhysInPlace (int(nx,ik),int(ny,ik),int(nz,ik),psi,&
              gxx,gxy,gxz,gyy,gyz,gzz)
      end if


!     ************************
!     ***   FIND  {lmax}   ***
!     ************************

      lmax = ahf_lmax

      if (lmax.ge.20) then
         write(*,*)
         write(*,*) 'ahf_lmax must be smaller than 20.'
         write(*,*) 'STOPPED IN AHFinder.F'
         write(*,*)
         call CCTK_WARN(0,"ahf_lmax must be smaller than 20.")
      end if


!     ***********************************************
!     ***   CHECK IF WE WANT TO FIND 3 HORIZONS   ***
!     ***********************************************

      if (ahf_find3.eq.0) then
         find3 = .false.
         nhorizon = 1
      else
         find3 = .true.
         nhorizon = 3
      end if


!     ***************************************************
!     ***   ALLOCATE STORAGE FOR COEFFICIENT ARRAYS   ***
!     ***************************************************

      if (ahf_ncall.eq.1) then

         allocate(c0(0:lmax),c0_old(0:lmax))
         allocate(cc(lmax,lmax),cc_old(lmax,lmax))
         allocate(cs(lmax,lmax),cs_old(lmax,lmax))

         c0 = zero
         cc = zero
         cs = zero

         c0_old = zero
         cc_old = zero
         cs_old = zero

         allocate(c0_0(0:lmax),c0_1(0:lmax),c0_2(0:lmax))
         allocate(cc_0(lmax,lmax),cc_1(lmax,lmax),cc_2(lmax,lmax))
         allocate(cs_0(lmax,lmax),cs_1(lmax,lmax),cs_2(lmax,lmax))

         c0_0 = zero
         cc_0 = zero
         cs_0 = zero

         c0_1 = zero
         cc_1 = zero
         cs_1 = zero

         c0_2 = zero
         cc_2 = zero
         cs_2 = zero

      end if

!     ***************************************************
!     ***   START LOOP TO LOOK FOR SEVERAL HORIZONS   ***
!     ***************************************************

!     IMPORTANT:  This is the point where we start the loop
!     to look for more than one horizon. Because of this,
!     anything that is likely to be different for the different
!     horizons should NOT be initialized until after this point!

      ahfmass = 0.0D0

      do mfind=1,nhorizon

!     Say what we are doing.

      if (nhorizon.gt.1) then
         if (mfind.eq.1) then
            write(*,*)
            write(*,*)
            write(*,*) 'AHFinder: MULTIPLE HORIZON MODE!'
         end if
         write(*,*)
         write(*,*) 'Searching for horizon ',mfind-1
      else
         write(*,*)
         write(*,*)
         write(*,*) 'AHFinder: SINGLE HORIZON MODE!'
         write(*,*)
         write(*,*) 'Searching for horizon'
      end if


!     ********************************
!     ***   INITIALIZE  {status}   ***
!     ********************************

      status = .false.


!     *********************************
!     ***   INITIALIZE FIRSTCALLS   ***
!     *********************************

      firstleg = .true.
      firstint = .true.


!     ***************************************
!     ***   FIND OUT THE TYPE OF OUTPUT   ***
!     ***************************************

      call CCTK_FortranString(nfile,out_dir,filestr)
      if (ahf_logfile.ne.0) then
         logfile = .true.
         logf = filestr(1:nfile)//"/ahf_logfile"
         if (find3) then
            logf = trim(logf) // "_" // char(ichar('0') + mfind - 1)
         end if
      else
         logfile = .false.
      end if

      if (ahf_verbose.ne.0) then
         verbose = .true.
      else
         verbose = .false.
      end if

      if (ahf_veryverbose.ne.0) then
         veryver = .true.
      else
         veryver = .false.
      end if

      if (ahf_guessverbose.ne.0) then
         guessverbose = .true.
      else
         guessverbose = .false.
      end if

      if (logfile.and.(myproc.eq.0)) then
         if (logfile) then
            if (myproc.eq.0) then
               open(11,file=logf,form='formatted',status='replace')
               write(11,*)
               write(11,*) 'LOG FILE FOR APPARENT HORIZON FINDER'
               write(11,*)
               write(11,*) 'Note:  During an evolution, this file will only'
               write(11,*) 'contain information about the last time that the'
               write(11,*) 'finder was called.'
               write(11,*)
               close(11)
            end if
         end if
      end if


!     **********************************
!     ***   INITIALIZE  {xc,yc,zc}   ***
!     **********************************

      if (find3) then
         ahf_xc_0 = 0.D0
         ahf_yc_0 = 0.D0
         ahf_zc_0 = 0.D0

         ahf_xc_1 = Xbh1
         ahf_yc_1 = Ybh1
         ahf_zc_1 = Zbh1

         ahf_xc_2 = Xbh2
         ahf_yc_2 = Ybh2
         ahf_zc_2 = Zbh2
         if (mfind.eq.1) then
            xc = ahf_xc_0
            yc = ahf_yc_0
            zc = ahf_zc_0
         else if (mfind.eq.2) then
            xc = ahf_xc_1
            yc = ahf_yc_1
            zc = ahf_zc_1
         else
            xc = ahf_xc_2
            yc = ahf_yc_2
            zc = ahf_zc_2
         end if
         write(*,*)'Centered on: '
         write(*,*)'xc = ',xc ! 'ahf_xc:',ahf_xc,ahf_centered
         write(*,*)'yc = ',yc ! 'ahf_xc:',ahf_yc,ahf_centered
         write(*,*)'zc = ',zc ! 'ahf_xc:',ahf_zc,ahf_centered
         write(*,*)'ahf_lmax = ',ahf_lmax

      else
         write(*,*)'Centered on: '
         write(*,*)'xbh1 = ',xbh1 ! 'ahf_xc:',ahf_xc,ahf_centered
         write(*,*)'ybh1 = ',ybh1 ! 'ahf_xc:',ahf_yc,ahf_centered
         write(*,*)'zbh1 = ',zbh1 ! 'ahf_xc:',ahf_zc,ahf_centered
         write(*,*)'ahf_lmax = ',ahf_lmax
         if(ahf_centered.eq.1) then
            ahf_xc = Xbh1
            ahf_yc = Ybh1
            ahf_zc = Zbh1
         endif
         xc = ahf_xc
         yc = ahf_yc
         zc = ahf_zc
      end if

!     Check if the center of the expansion is
!     out of the grid.  If it is, give up.

      giveup = .false.

      if ((xc.gt.xmx).or.(xc.lt.xmn)) then
         write(*,*)
         write(*,*) 'xc is out of the grid, giving up.'
         giveup = .true.
      end if

      if ((yc.gt.ymx).or.(yc.lt.ymn)) then
         write(*,*)
         write(*,*) 'yc is out of the grid, giving up.'
         giveup = .true.
      end if

      if ((zc.gt.zmx).or.(zc.lt.zmn)) then
         write(*,*)
         write(*,*) 'zc is out of the grid, giving up.'
         giveup = .true.
      end if

!     The end of this if statement is very far down,
!     it basically jumps all of the code.  To find
!     it look for the word "giveup" on a comment.

      if (.not.giveup) then


!     **************************
!     ***  FIND SYMMETRIES   ***
!     **************************

!     Find out if we want to move the center.

      if (ahf_offset.ne.0) then
         offset = .true.
      else
         offset = .false.
      end if

      if (ahf_wander.ne.0) then
         wander = .true.
      else
         wander = .false.
      end if

      if ((xc.ne.zero).or.(yc.ne.zero).or.(zc.ne.zero).or.     &
           (wander)) offset = .true.

!     Find out if we are in non-axisymmetric mode.

      if (ahf_phi.ne.0) then
         nonaxi = .true.
      else
         nonaxi = .false.
      end if


      call AHFinder_SetReflections(CCTK_ARGUMENTS)

!     Find the values of {stepx,stepy,stepz} depending
!     on the symmetries.

      stepx = 0
      stepy = 0
      stepz = 0

      if (refx) stepx=1
      if (refy) stepy=1
      if (refz) stepz=1

      if (CCTK_Equals(ahf_octant,"high").eq.1) stepx=3


!     *******************************
!     ***   FIND  {ntheta,nphi}   ***
!     *******************************

      ntheta = ahf_ntheta
      nphi   = ahf_nphi

!     For cartoon, I use only 2 points in the phi direction.
!     I could use one, but using 2 minimizes the changes in
!     the integration routine.

      if (cartoon) then
         nphi = 1
      end if


!     ********************************************************
!     ***   FIND TOTAL NUMBER OF TERMS IN EXPANSION {NN}   ***
!     ********************************************************

      NN = 0

!     Count {xc,yc,zc}.

      if (wander) then
         if (.not.refx) NN = NN + 1
         if (.not.refy) NN = NN + 1
         if (.not.refz) NN = NN + 1
      end if

!     Count c0(l).

      do l=0,lmax,1+stepz
         NN = NN + 1
      end do

      if (nonaxi) then

!        Count cc(l,m).

         do l=1,lmax
            do m=1+stepx,l,1+stepx
               if (stepz*mod(l-m,2).eq.0) NN = NN +1
            end do
         end do

!        Count cs(l,m).

         if (.not.refy) then
            do l=1,lmax
               do m=1,l,1+stepx
                  if (stepz*mod(l-m,2).eq.0) NN = NN + 1
               end do
            end do
         end if

      end if


!     *****************************************
!     ***   FIND RADIUS OF INITIAL SPHERE   ***
!     *****************************************

!     Find largest sphere that fits into the grid, taking
!     into account the position of the center.  Then multiply
!     this number with 0.9 to leave a safety margin at the edge
!     of the grid.

!     Cartoon.

      if (cartoon) then

         xc = zero
         yc = zero

         rmx = min(xmx,zmx-zc,zc-zmn)
         !below 3 lines added by Zach Etienne, for axi+equatorial symmetry
         if (CCTK_Equals(domain,"bitant").eq.1) then
            rmx = min(xmx-xc,zmx)
         end if
         write(*,*) "biggest radius = ",rmx

!     Octant.

      else if (CCTK_Equals(domain,"octant").eq.1) then

         if ((xc.eq.zero).and.(yc.eq.zero).and.(zc.eq.zero)) then
            rmx = min(xmx,ymx,zmx)
         else if ((xc.eq.zero).and.(yc.eq.zero)) then
            rmx = min(xmx,ymx,zmx-zc,zc)
         else if ((xc.eq.zero).and.(zc.eq.zero)) then
            rmx = min(xmx,ymx-yc,zmx,yc)
         else if ((yc.eq.zero).and.(zc.eq.zero)) then
            rmx = min(xmx-xc,ymx,zmx,xc)
         else if (xc.eq.zero) then
            rmx = min(xmx,ymx-yc,zmx-zc,yc,zc)
         else if (yc.eq.zero) then
            rmx = min(xmx-xc,ymx,zmx-zc,xc,zc)
         else if (zc.eq.zero) then
            rmx = min(xmx-xc,ymx-yc,zmx,xc,yc)
         else
            rmx = min(xmx-xc,ymx-yc,zmx-zc,xc,yc,zc)
         end if

!     Quadrant.

      else if (CCTK_Equals(domain,"quadrant").eq.1) then
         if ((xc.eq.zero).and.(yc.eq.zero)) then
            rmx = min(xmx,ymx,zmx-zc,zc-zmn)
         else if (xc.eq.zero) then
            rmx = min(xmx,ymx-yc,zmx-zc,yc,zc-zmn)
         else if (yc.eq.zero) then
            rmx = min(xmx-xc,ymx,zmx-zc,xc,zc-zmn)
         else
            rmx = min(xmx-xc,ymx-yc,zmx-zc,xc,yc,zc-zmn)
         end if

!     Bitant.

      else if (CCTK_Equals(domain,"bitant").eq.1) then

         if (CCTK_Equals(bitant_plane,"xy").eq.1) then

            if (zc.eq.zero) then
               rmx = min(xmx-xc,ymx-yc,xc-xmn,yc-ymn,zmx)
            else
               rmx = min(xmx-xc,ymx-yc,xc-xmn,yc-ymn,zmx-zc,zc)
            end if

         else if (CCTK_Equals(bitant_plane,"xz").eq.1) then

            if (yc.eq.zero) then
               rmx = min(xmx-xc,zmx-zc,xc-xmn,zc-zmn,ymx)
            else
               rmx = min(xmx-xc,zmx-zc,xc-xmn,zc-zmn,ymx-yc,yc)
            end if

         else

            if (xc.eq.zero) then
               rmx = min(ymx-yc,zmx-zc,yc-ymn,zc-zmn,xmx)
            else
               rmx = min(ymx-yc,zmx-zc,yc-ymn,zc-zmn,xmx-xc,xc)
            end if

         end if

!     Full.

      else

         rmx = min(xmx-xc,ymx-yc,zmx-zc,xc-xmn,yc-ymn,zc-zmn)

      end if

      rmx = 0.9*rmx

!     Check if a special radius for the initial sphere
!     has been forced by the parameter file.

      if (find3) then
         if (mfind.eq.1) then
            r0 = ahf_r0_0
         else if (mfind.eq.2) then
            r0 = ahf_r0_1
         else
            r0 = ahf_r0_2
         end if
      else
         r0 = ahf_r0
      end if

      if (r0.ne.zero) then
         if (r0.le.rmx) then
            rmx = r0
         else
            write (*,*)
            write (*,*) 'ahf_r0 is too large, rescaling ...'
         end if
      end if


!     *****************************
!     ***   GET INITIAL GUESS   ***
!     *****************************

!     Initializing flags

      sloppy = (ahf_sloppyguess.eq.1)

      if (cartoon) sloppy = .true.

      flow = (ahf_flow.eq.1)

      inner = (ahf_inner.eq.1)

      if (ahf_guessold.ne.0) then
         guessold = .true.
      else
         guessold = .false.
      end if

      if (ahf_guess_absmin.ne.0) then
         guess_absmin = .true.
      else
         guess_absmin = .false.
      end if

      if (ahf_manual_guess.ne.0) then
         manual_guess = .true.
      else
         manual_guess = .false.
      end if

      if (ahf_minarea.ne.0) then
         minarea = .true.
      else
         minarea = .false.
      end if

!     Get initial guess.

      if (find3) then
         if (mfind.eq.1) then
            status_old = status_old_0
         else if (mfind.eq.2) then
            status_old = status_old_1
         else
            status_old = status_old_2
         end if
      end if

      if (guessold.and.status_old) then

!        Old horizon as initial guess.

         if ((myproc.eq.0).and.verbose) then
            write(*,*)
            write(*,*)  'Using old horizon as initial guess'
         end if

!        If we are only looking for 1 horizon, then we are OK
!        since the arrays (c0,cc,cs) where saved from the last
!        call.  If we are looking for 3 horizons, we saved all
!        3 sets of coefficients, but now we need to copy them
!        back to the right place.

         if (find3) then
            if (mfind.eq.1) then
               c0 = c0_0
               cc = cc_0
               cs = cs_0
            else if (mfind.eq.2) then
               c0 = c0_1
               cc = cc_1
               cs = cs_1
            else
               c0 = c0_2
               cc = cc_2
               cs = cs_2
            end if
         end if

      else if (manual_guess) then

!        Manual guess.

         if ((myproc.eq.0).and.verbose) then
            write(*,*)
            write(*,*)  'Using manual initial guess'
         end if

         c0_temp(0) = ahf_l0_guess
         c0_temp(2) = ahf_l2_guess
         c0_temp(4) = ahf_l4_guess
         c0_temp(6) = ahf_l6_guess
         c0_temp(8) = ahf_l8_guess
         c0_temp(10) = ahf_l10_guess
         c0_temp(12) = ahf_l12_guess
         c0_temp(14) = ahf_l14_guess
         c0_temp(16) = ahf_l16_guess
         c0_temp(18) = ahf_l18_guess

         if (lmax.le.18) then
            do i=0,lmax
               c0(i) = c0_temp(i)
            end do
         else
            do i=0,18
               c0(i) = c0_temp(i)
            end do
         end if

      else

!        Standard initial guess.

         if (.not.flow) then

!           For minimization, call initial guess routine.

            call AHFinder_initguess(CCTK_ARGUMENTS,rmx)

         else

!           For flow, take as initial guess a large sphere.

            c0 = zero
            cc = zero
            cs = zero

            c0_old = zero
            cc_old = zero
            cs_old = zero

            if ((myproc.eq.0).and.verbose) then
               write(*,*)
               write(*,"(A35,ES11.3)")  &
                    ' Initial guess sphere with radius: ',rmx
               write(*,"(A13,3ES11.3)") &
                    ' centered on:',real(xc),real(yc),real(zc)
            end if

            c0(0) = rmx

         end if

      end if


!     ***************************
!     ***   LOOK FOR HORIZON  ***
!     ***************************

      if (.not.flow) then

!        Minimization algorithm.

         call AHFinder_min(CCTK_ARGUMENTS,NN,status,logf)

      else

!        Flow algorithm.

         if ((ahf_ncall.eq.1).and.(mfind.eq.1)) then
            allocate(hflow0(0:lmax),cflow0(0:lmax),nflow0(0:lmax))
            allocate(hflowc(lmax,lmax),cflowc(lmax,lmax),nflowc(lmax,lmax))
            allocate(hflows(lmax,lmax),cflows(lmax,lmax),nflows(lmax,lmax))
         end if

         call AHFinder_flow(CCTK_ARGUMENTS,rmx,status,logf)
         flow = .false.

      end if


!     ***************************************
!     ***   WRITE MESSAGES AND LOG FILE   ***
!     ***************************************

      if (find3) then
         if (mfind.eq.1) then
            status_old_0 = status
         else if (mfind.eq.2) then
            status_old_1 = status
         else
            status_old_2 = status
         end if
      else
         status_old = status
      end if

      if (status) then

!        Check if the expansion on the surface is small enough
!        for a horizon.  How small is small enough is of course
!        pretty arbitrary, so to reduce the subjectivity I use
!        three separate criteria:
!
!        1) The integral of the square of the expansion at the surface
!           must be smaller than a given threshold (small1).
!
!        2) The mean value of the expansion on the surface must be
!           smaller in absolute value than its standard deviation.
!           This means that the mean value of the expansion is
!           consistent with zero.
!
!        3) The integral of the square of the expansion at the surface
!           must be smaller by a factor of "small2" than the same
!           integral slightly outside and slightly inside.

         horizon = .true.

         small1 = 1.0D-2
         small2 = 0.5D0

!        Surface found.

         call AHFinder_fun(CCTK_ARGUMENTS)
         call AHFinder_exp(CCTK_ARGUMENTS)
         call AHFinder_int(CCTK_ARGUMENTS)

         if (ahf_gaussout.ne.0) then
            call AHFinder_gau(CCTK_ARGUMENTS)
         else
            circ_eq = 0.0D0
            meri_p1 = 0.0D0
            meri_p2 = 0.0D0
         end if

         intexp_h  = intexp
         intexp2_h = intexp2
         intarea_h = intarea

         if (intexp2_h.gt.small1) horizon = .false.

!        Find standard deviation:  sqrt(<exp^2> - <exp>^2)
!        The term inside the parenthesis is supposed to be positive,
!        but may fail to be because of numerical errors in pathological
!        situations, so be careful.

         aux = intexp2_h/intarea_h - (intexp_h/intarea_h)**2

         if (aux.ge.zero) then
            aux = sqrt(aux)
            if (abs(intexp_h/intarea_h)-aux.gt.zero) horizon = .false.
         else
            horizon = .false.
         end if

!        Surface slightly inside

         aux = c0(0)

         c0(0) = aux - half*min(dx,dy,dz)
!         write(*,*) "hello!",c0(0),aux,half*min(dx,dy,dz),dx,dy,dz
         call AHFinder_int(CCTK_ARGUMENTS)

         expin = intexp

         if (intexp2_h.gt.small2*intexp2) horizon = .false.

         if (veryver) then
            write(*,*)
            write(*,*) 'Surface slightly inside the minimum surface: '
            write(*,"(A21,ES14.6)") ' Integral of H      =',intexp
            write(*,"(A21,ES14.6)") ' Integral of H^2    =',intexp2
            write(*,"(A21,ES14.6)") ' Mean value of H    =',intexp/intarea
            write(*,"(A21,ES14.6)") ' Mean value of H^2  =',intexp2/intarea
            write(*,*) 'Number of interpolated points:     ', &
                 inside_min_count
            write(*,*) 'Number of those that are negative: ', &
                 inside_min_neg_count
         endif

!        Surface slightly outside.

         c0(0) = aux + half*min(dx,dy,dz)
         call AHFinder_int(CCTK_ARGUMENTS)

         expout = intexp

         if (intexp2_h.gt.small2*intexp2) horizon = .false.

         if (veryver) then
            write(*,*)
            write(*,*) 'Surface slightly outside the minimum surface: '
            write(*,"(A21,ES14.6)") ' Integral of H      =',intexp
            write(*,"(A21,ES14.6)") ' Integral of H^2    =',intexp2
            write(*,"(A21,ES14.6)") ' Mean value of H    =',intexp/intarea
            write(*,"(A21,ES14.6)") ' Mean value of H^2  =',intexp2/intarea
            write(*,*) 'Number of interpolated points:     ', &
                 inside_min_count
            write(*,*) 'Number of those that are negative: ', &
                 inside_min_neg_count
         endif

         c0(0) = aux

!        What kind of horizon?

         if ((expin.le.zero).and.(expout.ge.zero)) then
            mtype = 1
         else if ((expin.ge.zero).and.(expout.le.zero)) then
            mtype = 2
         else if ((expin.le.zero).and.(expout.le.zero)) then
            mtype = 3
         else
            mtype = 4
         end if

!        Messages to screen and logfile

         if (ahf_ReportAlways.eq.1) then
            report = .true.
         else
            report = mtype.eq.1
         end if

!        Report on the details of the surface.

         if (report.and.(myproc.eq.0)) then

!           Open logfile.

            if (logfile) then
               open(11,file=logf,form='formatted', &
                    status='old',position='append')
            end if

!           Find inverse of area.  Be carefull not to divide by zero!
!           Also, a value of 1.0D10 indicates an error in the integration
!           routine, and all integrals end up with the same value.
!           So do not divide by it or I will be confused later.

            if ((intarea_h.ne.zero).and.(intarea.lt.1.0D10)) then
               aux = one/intarea_h
            else
               aux = one
            end if

!           Surface information.  Notice that here I define the
!           surface "mass" as:
!
!           M = sqrt( A / (16 pi) ) = 0.141047396 sqrt(A)

            ahfmass(mfind) = 0.141047396D0*sqrt(intarea_h)

            write(*,*)
            write(*,*) 'Surface found, details below.'
            write(*,*)
            write(*,"(A21,ES14.6)") ' Surface area       =',intarea_h
            write(*,"(A21,ES14.6)") ' Surface mass       =', &
                 0.141047396D0*sqrt(intarea_h)
            write(*,"(A21,ES14.6)") ' Mean value of H    =',aux*intexp_h
            write(*,"(A21,ES14.6)") ' Mean value of H^2  =',aux*intexp2_h
            write(*,"(A21,ES14.6)") ' Standard deviation =', &
                 sqrt(abs(aux*intexp2_h - (aux*intexp_h)**2))
            write(*,*)
            if (ahf_gaussout.ne.0) then
               write(*,"(A11,ES14.6)") ' circ_eq  =',circ_eq
               write(*,"(A11,ES14.6)") ' meri_p1  =',meri_p1
               write(*,"(A11,ES14.6)") ' meri_p2  =',meri_p2
               write(*,*)
            end if

            if (logfile) then
               write(11,*)
               write(11,*) 'Surface found, details below'
               write(11,*)
               write(11,"(A21,ES14.6)") ' Surface area       =',intarea_h
               write(11,"(A21,ES14.6)") ' Surface mass       =', &
                    0.141047396D0*sqrt(intarea_h)
               write(11,"(A21,ES14.6)") ' Mean value of H    =',aux*intexp_h
               write(11,"(A21,ES14.6)") ' Mean value of H^2  =',aux*intexp2_h
               write(11,"(A21,ES14.6)") ' Standard deviation =', &
                    sqrt(abs(aux*intexp2_h - (aux*intexp_h)**2))
               write(11,*)
               if (ahf_gaussout.ne.0) then
                  write(11,"(A11,ES14.6)") ' circ_eq  =',circ_eq
                  write(11,"(A11,ES14.6)") ' meri_p1  =',meri_p1
                  write(11,"(A11,ES14.6)") ' meri_p2  =',meri_p2
                  write(11,*)
               end if
            end if

            if (mtype.eq.1) then
               if (horizon) then
                 write(*,*) 'The surface seems to be an outer horizon.'
                 if (logfile) then
                    write(11,*) 'The surface seems to be an outer horizon.'
                 end if
               else
                 write(*,*) 'The surface is probably an outer horizon, but'
                 write(*,*) 'the mean of the square of the expansion is'
                 write(*,*) 'not small enough. Try a larger lmax, or a'
                 write(*,*) 'smaller grid size.'
                 if (logfile) then
                    write(11,*) 'The surface is probably an outer horizon, but'
                    write(11,*) 'the mean of the square of the expansion is'
                    write(11,*) 'not small enough. Try a larger lmax, or a'
                    write(11,*) 'smaller grid size.'
                 end if
               end if
            else if (mtype.eq.2) then
               if (horizon) then
                 write(*,*) 'The surface seems to be an inner horizon.'
                 if (logfile) then
                    write(11,*) 'The surface seems to be an inner horizon.'
                 end if
               else
                 write(*,*) 'The surface is probably an inner horizon, but'
                 write(*,*) 'the mean of the square of the expansion is'
                 write(*,*) 'not small enough. Try a larger lmax, or a'
                 write(*,*) 'smaller grid size.'
                 if (logfile) then
                    write(11,*) 'The surface is probably an inner horizon, but'
                    write(11,*) 'the mean of the square of the expansion is'
                    write(11,*) 'not small enough. Try a larger lmax, or a'
                    write(11,*) 'smaller grid size.'
                 end if
               end if
            else if (mtype.eq.3)then
               if (horizon) then
                 write(*,*) 'The surface seems to be a marginal horizon'
                 write(*,*) '(the integral of the expansion is negative'
                 write(*,*) 'on both sides.)'
                 if (logfile) then
                    write(11,*) 'The surface seems to be a marginal horizon'
                    write(11,*) '(the integral of the expansion is negative'
                    write(11,*) 'on both sides.)'
                 end if
               else
                 write(*,*) 'The surface seems to be a trapped surface.'
                 if (logfile) then
                    write(11,*) 'The surface seems to be a trapped surface.'
                 end if
               end if
            else if (mtype.eq.4)then
               if (horizon) then
                 write(*,*) 'The surface seems to be a marginal horizon'
                 write(*,*) '(the integral of the expansion is positive'
                 write(*,*) 'on both sides.)'
                 if (logfile) then
                    write(11,*) 'The surface seems to be a marginal horizon'
                    write(11,*) '(the integral of the expansion is positive'
                    write(11,*) 'on both sides.)'
                 end if
               else
                 write(*,*) 'The surface does not seem to be a horizon.'
                 if (logfile) then
                    write(11,*) 'The surface does not seem to be a horizon.'
                 end if
               end if
            end if

!           Shape.

            if (verbose) then

               write(*,*)
               write(*,*) 'Shape coefficients:'

               if (offset) then
                  write(*,*)
                  write(*,"(A6,ES14.6)") ' xc  =',xc
                  write(*,"(A6,ES14.6)") ' yc  =',yc
                  write(*,"(A6,ES14.6)") ' zc  =',zc
               end if

               write(*,*)

               write(*,"(A4,I2,A3,ES14.6)") ' c0(',0,') =',c0(0)

               do l=1+stepz,lmax,1+stepz
                  write(*,"(A4,I2,A3,ES14.6)") ' c0(',l,') =',c0(l)
               end do

            end if

            if (logfile) then

               write(11,*)
               write(11,*) 'Shape coefficients:'

               if (offset) then
                  write(11,*)
                  write(11,"(A6,ES14.6)") ' xc  =',xc
                  write(11,"(A6,ES14.6)") ' yc  =',yc
                  write(11,"(A6,ES14.6)") ' zc  =',zc
               end if

               write(11,*)

               write(11,"(A4,I2,A3,ES14.6)") ' c0(',0,') =',c0(0)

               do l=1+stepz,lmax,1+stepz
                  write(11,"(A4,I2,A3,ES14.6)") ' c0(',l,') =',c0(l)
               end do

            end if

            if (nonaxi) then

               if (.not.refy) then

                  if (verbose) then

                     write(*,*)
                     do l=1,lmax
                        do m=1,l
                           if (stepz*mod(l-m,2).eq.0) then
                              write(*,"(A4,I2,A1,I2,A4,ES14.6)") ' cc(', &
                                   l,',',m,')  =',cc(l,m)
                           end if
                        end do
                     end do
                     write(*,*)
                     do l=1,lmax
                        do m=1,l
                           if (stepz*mod(l-m,2).eq.0) then
                              write(*,"(A4,I2,A1,I2,A4,ES14.6)") ' cs(', &
                                   l,',',m,')  =',cs(l,m)
                           end if
                        end do
                     end do

                  end if

                  if (logfile) then
                     write(11,*)

                     do l=1,lmax
                        do m=1,l
                           if (stepz*mod(l-m,2).eq.0) then
                              write(11,"(A4,I2,A1,I2,A4,ES14.6)") ' cc(', &
                                   l,',',m,')  =',cc(l,m)
                           end if
                        end do
                     end do
                     write(11,*)
                     do l=1,lmax
                        do m=1,l
                           if (stepz*mod(l-m,2).eq.0) then
                              write(11,"(A4,I2,A1,I2,A4,ES14.6)") ' cs(', &
                                   l,',',m,')  =',cs(l,m)
                           end if
                        end do
                     end do

                  end if

               else

                  if (verbose) then

                     write(*,*)
                     do l=1,lmax
                        do m=1,l
                           if (stepz*mod(l-m,2).eq.0) then
                              write(*,"(A4,I2,A1,I2,A4,ES14.6)") ' cc(', &
                                   l,',',m,')  =',cc(l,m)
                           end if
                        end do
                     end do

                  end if

                  if (logfile) then

                     write(11,*)
                     do l=1,lmax
                        do m=1,l
                           if (stepz*mod(l-m,2).eq.0) then
                              write(11,"(A4,I2,A1,I2,A4,ES14.6)") ' cc(', &
                                   l,',',m,')  =',cc(l,m)
                           end if
                        end do
                     end do

                  end if

               end if
            end if

!           Close logfile

            write(*,*)

            if (logfile) then
               write(11,*)
               write(11,*) 'END OF LOG FILE'
               write(11,*)
               close(11)
            end if

!        If we found a surface that is not a horizon and we
!        do not want to report all its details, just say that
!        we did not find a horizon.

         else if (myproc.eq.0) then

!           Open logfile.

            if (logfile) then
               open(11,file=logf,form='formatted', &
                    status='old',position='append')
            end if

!           Write messages.

            write(*,*)
            write(*,*) 'No horizon found.'

            if (logfile) then
               write(11,*)
               write(11,*) 'No horizon found.'
               write(11,*)
            end if

!           Close logfile

            write(*,*)

            if (logfile) then
               write(11,*)
               write(11,*) 'END OF LOG FILE'
               write(11,*)
               close(11)
            end if

         end if

!     No surface found.

      else

         mtype = 0
         report = .false.

         if (myproc.eq.0) then

!           Open logfile.

            if (logfile) then
               open(11,file=logf,form='formatted', &
                    status='old',position='append')
            end if

!           Write messages.

            write(*,*)
            write(*,*) 'No horizon found.'

            if (logfile) then
               write(11,*)
               write(11,*) 'No horizon found.'
               write(11,*)
            end if

!           Close logfile

            write(*,*)

            if (logfile) then
               write(11,*)
               write(11,*) 'END OF LOG FILE'
               write(11,*)
               close(11)
            end if

         end if
      end if


!     ******************************************************
!     ***   PREPARING  {ahfgrid3,ahf_exp3}  FOR OUTPUT   ***
!     ******************************************************

      if (find3) then
         call AHFinder_find3(CCTK_ARGUMENTS,mtype)
      end if


!     ****************************
!     ***   WRITE DATA FILES   ***
!     ****************************

      call AHFinder_Output(CCTK_ARGUMENTS,report,status,horizon,mtype,intarea_h)


!     ****************
!     ***   MASK   ***
!     ****************

!     Initialize the sides of the excised region to the center
!     of the expansion.

      if (ahf_ncall.eq.1) then

         if (find3) then

            if (mfind.eq.1) then

               dhole1_xmin = xc
               dhole1_ymin = yc
               dhole1_zmin = zc

               dhole1_xmax = xc
               dhole1_ymax = yc
               dhole1_zmax = zc

            else if (mfind.eq.2) then

               dhole2_xmin = xc
               dhole2_ymin = yc
               dhole2_zmin = zc

               dhole2_xmax = xc
               dhole2_ymax = yc
               dhole2_zmax = zc

            else

               dhole3_xmin = xc
               dhole3_ymin = yc
               dhole3_zmin = zc

               dhole3_xmax = xc
               dhole3_ymax = yc
               dhole3_zmax = zc

            end if

         else

            dhole1_xmin = xc
            dhole1_ymin = yc
            dhole1_zmin = zc

            dhole1_xmax = xc
            dhole1_ymax = yc
            dhole1_zmax = zc

         end if

      end if

!     The mask should only be different from 1 inside a horizon.
!     I therefore need to check first if I do have a horizon, and
!     also if the mask is desired.

      if ((horizon.and.(CCTK_EQUALS(ahf_mask,'strong'))).or. &
           (status.and.(mtype.eq.1).and.(CCTK_EQUALS(ahf_mask,'weak')))) then

         rhor = 0.0D0

         if (find3) then
            if ((mfind.eq.1).and.(ahf_mask_0.eq.1)) then
               call AHFinder_mask(CCTK_ARGUMENTS,rhor)
            else if ((mfind.eq.2).and.(ahf_mask_1.eq.1)) then
               call AHFinder_mask(CCTK_ARGUMENTS,rhor)
            else if ((mfind.eq.3).and.(ahf_mask_2.eq.1)) then
               call AHFinder_mask(CCTK_ARGUMENTS,rhor)
            end if
         else
            call AHFinder_mask(CCTK_ARGUMENTS,rhor)
         end if

!        Here I alter the parameters for SimpleExcision.  Notice that at
!        the moment, I only allow the excised cube to grow, and never
!        to shrink.  Also, notice that the sides of the cube have been
!        set to 0.57*rhor.  This is because if we want the corners of
!        the cube to touch a sphere of radius rhor, we need to make the
!        sides equal to 1/sqrt(3), which is roughly 0.57.  If you worry
!        that this is too close to the horizon, remember that rhor in fact
!        has already been shrunk by a safety margin.

         aux = 0.57D0*rhor

         if (find3) then

            if (mfind.eq.1) then

!              First hole.

               if (dhole1_xmin.gt.(xc-aux)) dhole1_xmin = xc-aux
               if (dhole1_ymin.gt.(yc-aux)) dhole1_ymin = yc-aux
               if (dhole1_zmin.gt.(zc-aux)) dhole1_zmin = zc-aux

               if (dhole1_xmax.lt.(xc+aux)) dhole1_xmax = xc+aux
               if (dhole1_ymax.lt.(yc+aux)) dhole1_ymax = yc+aux
               if (dhole1_zmax.lt.(zc+aux)) dhole1_zmax = zc+aux

            else if (mfind.eq.2) then

!              Second hole.

               if (dhole2_xmin.gt.(xc-aux)) dhole2_xmin = xc-aux
               if (dhole2_ymin.gt.(yc-aux)) dhole2_ymin = yc-aux
               if (dhole2_zmin.gt.(zc-aux)) dhole2_zmin = zc-aux

               if (dhole2_xmax.lt.(xc+aux)) dhole2_xmax = xc+aux
               if (dhole2_ymax.lt.(yc+aux)) dhole2_ymax = yc+aux
               if (dhole2_zmax.lt.(zc+aux)) dhole2_zmax = zc+aux

            else

!              Third hole.

               if (dhole3_xmin.gt.(xc-aux)) dhole3_xmin = xc-aux
               if (dhole3_ymin.gt.(yc-aux)) dhole3_ymin = yc-aux
               if (dhole3_zmin.gt.(zc-aux)) dhole3_zmin = zc-aux

               if (dhole3_xmax.lt.(xc+aux)) dhole3_xmax = xc+aux
               if (dhole3_ymax.lt.(yc+aux)) dhole3_ymax = yc+aux
               if (dhole3_zmax.lt.(zc+aux)) dhole3_zmax = zc+aux

            end if

         else

            if (dhole1_xmin.gt.(xc-aux)) dhole1_xmin = xc-aux
            if (dhole1_ymin.gt.(yc-aux)) dhole1_ymin = yc-aux
            if (dhole1_zmin.gt.(zc-aux)) dhole1_zmin = zc-aux

            if (dhole1_xmax.lt.(xc+aux)) dhole1_xmax = xc+aux
            if (dhole1_ymax.lt.(yc+aux)) dhole1_ymax = yc+aux
            if (dhole1_zmax.lt.(zc+aux)) dhole1_zmax = zc+aux

         end if

      end if


!     ********************************************************
!     ***   END OF LOOP FOR LOOKING FOR SEVERAL HORIZONS   ***
!     ********************************************************

!     End of "giveup" if statement.

      end if

!     If we are looking for 3 horizons, save values of coefficients.

      if (find3) then
         if (mfind.eq.1) then
            c0_0 = c0
            cc_0 = cc
            cs_0 = cs
         else if (mfind.eq.2) then
            c0_1 = c0
            cc_1 = cc
            cs_1 = cs
         else
            c0_2 = c0
            cc_2 = cc
            cs_2 = cs
         end if
      end if

      write (*,*)

      end do


!     *******************************************************
!     ***   CALCULATE PROPER DISTANCES BETWEEN HORIZONS   ***
!     *******************************************************

      if (find3) then
         call AHFinder_dis(CCTK_ARGUMENTS)
      end if


!     **************************************************
!     ***   TRANSFORM PHYSICAL TO CONFORMAL METRIC   ***
!     **************************************************

!     If necesary, change back the metric to conformal
!     and change back the "conformal_state" flag.

      if (ahfconformal) then
         conformal_state = save_conformal_state
         call PhysToConfInPlace (int(nx,ik),int(ny,ik),int(nz,ik),psi,&
             gxx,gxy,gxz,gyy,gyz,gzz)
      end if
!     ***************
!     ***   END   ***
!     ***************

      end subroutine AHFinder
