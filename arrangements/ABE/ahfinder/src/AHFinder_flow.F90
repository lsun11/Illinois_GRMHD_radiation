/*@@
  @file      AHFinder_flow.F
  @date      November 1998
  @author    Miguel Alcubierre
  @desc 
             Master routine to control the iterations for
             the flow algorithm.

             The flow algorithm used here is that of Carsten Gunlach:
             Phys. Rev. D 57, 863 (1998).

  @enddesc 
@@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

      subroutine AHFinder_flow(CCTK_ARGUMENTS,rmx,status,logf)

      use AHFinder_dat

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS 

      logical status,new,update

      integer i,j,l,m
      integer itmax

      CCTK_REAL rmx,ahftol
      CCTK_REAL A,B
      CCTK_REAL ahfsum,ahfdiff,ahfdiff_old
      CCTK_REAL reldiff
      CCTK_REAL maxchange,minchange
      CCTK_REAL zero,one
      CCTK_REAL dd,aux

      character(len=200) :: logf

!     Description of variables:
!
!     new         Flag to indicate that stepsize has changed.
!     update      Flag to indicate if coefficients will be updated.
!
!     i,j,l,m     Counters.
!
!     itmax       Maximum number of iterations.
!
!     ahftol      Tolerance used to decide if we have converged.
!
!     A,B         Flow parameters.
!
!     hw          Weigth of `H' flow (see Carsten's paper).
!     cw          Weigth of `C' flow (see Carsten's paper).
!     nw          Weigth of `N' flow (see Carsten's paper).
!
!     ahfsum      Square root of the sum of the squares of the
!                 expansion coefficients.
!
!     ahfdiff     Square root of the sum of the squares of the
!                 change in the coefficients.
!
!     rediff      Relative difference between having done one large
!                 step and two small ones.
!
!     maxchange   A value of "reldiff" larger than this forces a
!                 reduction of the stepsize.
!
!     minchange   A value of "reldiff" smaller than this forces an
!                 increase of the stepsize.
!
!     dd          Maximum of (dx,dy,dz).


!     ***************************
!     ***   DEFINE  NUMBERS   ***
!     ***************************

      zero = 0.0D0
      one  = 1.0D0

      dd = max(dx,dy,dz)


!     ***************************
!     ***   FIND PARAMETERS   ***
!     ***************************

      A = ahf_flowa
      B = ahf_flowb

      hw = ahf_flowh
      cw = ahf_flowc
      nw = ahf_flown

      maxchange = ahf_maxchange
      minchange = ahf_minchange


!     *******************************
!     ***   MAIN ITERATION LOOP   ***
!     *******************************

!     Find itmax and ahftol.

      itmax  = ahf_flowiter
      ahftol = ahf_flowtol

!     Start the iterations.

      do i=1,itmax

!        Save the old values of the coefficients.

         c0_old = c0
         cc_old = cc
         cs_old = cs

!        In order to adapt the stepsize, I do two steps
!        with the current stepsize, and compare them with
!        one step with double the stepsize.

         do j=1,2

            ahfsum  = zero
            ahfdiff = zero


!           ************************************
!           ***   FIND SPECTRAL COMPONENTS   ***
!           ************************************

            call AHFinder_fun(CCTK_ARGUMENTS)
            call AHFinder_exp(CCTK_ARGUMENTS)
            call AHFinder_int(CCTK_ARGUMENTS)


!           *****************************
!           ***   CHECK ERROR FLAGS   ***
!           *****************************

            update = .true.

!           If interror is true there was an error in the integral.
!           There are two possibilities here:

!           1) This is the first of the two small steps (j=1).
!              We then went out of bounds in the previous iteration,
!              so there is nothing we can do to fix it.

            if (interror.and.(j.eq.1)) then

!              We can be out of bounds for three reasons:

!              a) The radius is negative somewhere.  This can mean
!                 one of two things:
!
!                 * There is no horizon, and the flow has reached
!                   the origin.
!
!                 * The origin for our expansion is outside the horizon.
!                   This can happen if the black-hole is not centered,
!                   or if we have more than one horizon.

               if (interror1.ne.0) then

                  write(*,*)
                  write(*,*) 'Negative radius.'
                  return

!              b) We are out of the computational domain, pitty.

               else if (interror2.ne.0) then

                  write(*,*)
                  write(*,*) 'Out of bounds, giving up.'
                  return

!              c) We are inside the mask.

               else if (interror3.ne.0) then

                  write(*,*)
                  write(*,*) 'Inside mask, giving up.'
                  return

!              We should never get here!

               else

                  write(*,*)
                  write(*,*) 'PANIC! Something is very wrong.'
                  return

               end if

!           2) This is the second of the two small steps (j=2).
!              We can still try to reduce the stepsize and start
!              the iteration again.  We then return the coefficients
!              to their original values and, and do not update them.
!
!              Since I will not update the coefficients, I need to
!              assign values to "afhsum" and "ahfdiff". In order
!              to force a reduction of the step size, I assign
!              a very large number to "ahfdiff", which will make
!              "reldiff" close to 1.0, far too large.  I also
!              assign a value of 1.0 to "ahfsum", so it is much
!              smaller than "ahfdiff".  This is important in order
!              to prevent the algorithm from thinking we are done
!              (see below for termination criteria).

            else if (interror) then

               interror = .false.

               c0 = c0_old
               cc = cc_old
               cs = cs_old

               ahfsum   = one
               ahfdiff  = 1.0D10

               update = .false.

            end if


!           ****************************************
!           ***   UPDATE SPECTRAL COEFFICIENTS   ***
!           ****************************************

            if (update) then

!              Update c0.

               do l=0,lmax,1+stepz
                  aux = c0(l)
                  c0(l) = aux - A/(one + B*dble(l)*(dble(l) + one)) &
                       *(hw*hflow0(l) + cw*cflow0(l) + nw*nflow0(l))
                  ahfsum  = ahfsum + c0_old(l)**2
                  ahfdiff = ahfdiff + (c0(l) - c0_old(l))**2
               end do

!              Update {cc,cs}.

               if (nonaxi) then

!                 Update cc(l,m).

                  do l=1,lmax
                     do m=1+stepx,l,1+stepx
                        if (stepz*mod(l-m,2).eq.0) then
                           aux = cc(l,m)
                           cc(l,m) = aux - A/(one &
                                + B*dble(l)*(dble(l) + one)) &
                                *(hw*hflowc(l,m) + cw*cflowc(l,m) &
                                + nw*nflowc(l,m))
                           ahfsum  = ahfsum + cc_old(l,m)**2
                           ahfdiff = ahfdiff + (cc(l,m)-cc_old(l,m))**2
                        end if
                     end do
                  end do

!                 Update cs(l,m).

                  if (.not.refy) then
                     do l=1,lmax
                        do m=1,l,1+stepx
                           if (stepz*mod(l-m,2).eq.0) then
                              cs(l,m) = cs(l,m) - A/(one &
                                   + B*dble(l)*(dble(l) + one)) &
                                   *(hw*hflows(l,m) + cw*cflows(l,m) &
                                   + nw*nflows(l,m))
                              ahfsum  = ahfsum + cs_old(l,m)**2
                              ahfdiff = ahfdiff + (cs(l,m)-cs_old(l,m))**2
                           end if
                        end do
                     end do
                  end if

               end if

               ahfsum  = sqrt(ahfsum)
               ahfdiff = sqrt(ahfdiff)

!              If this is the first of the two small steps, save
!              twice the difference.  Why twice?  Because this would
!              have been the resulting difference if we had done only
!              one step that was twice as big.

               if (j.eq.1) then
                  ahfdiff_old = 2.0D0*ahfdiff
               end if

            end if

         end do


!        ************************************
!        ***   FIND RELATIVE DIFFERENCE   ***
!        ************************************

!        Here we find the relative difference between having done two
!        small steps and one large step.

         reldiff = abs(ahfdiff - ahfdiff_old) &
              / (abs(ahfdiff) + abs(ahfdiff_old))


!        Open logfile.

         logf = filestr(1:nfile)//"/ahf_logfile"

         if (logfile.and.(myproc.eq.0)) then
            open(11,file=logf,form='formatted',status='old', &
                 position='append')
         end if

!        ***************************
!        ***   ADJUST STEPSIZE   ***
!        ***************************

         new = .false.

!        If the relative difference is too large, we reduce the
!        stepsize for the next iteration.

         if ((reldiff.gt.maxchange).or. &
              (ahfdiff.gt.maxchange*ahfsum)) then

            new = .true.
            A = 0.5D0*A

            if (veryver) then
               write(*,*)
               write(*,"(A15,ES14.6)") ' New stepsize =',A
            end if

            if (logfile) then
               write(11,*)
               write(11,"(A15,ES14.6)") ' New stepsize =',A
            end if

!        If the diference is too small, we can safely increase
!        the stepsize for the next iteration.

         else if (reldiff.lt.minchange) then

            new = .true.
            A = 2.0D0*A

            if (veryver) then
               write(*,*)
               write(*,"(A15,ES14.6)") ' New stepsize =',A
            end if

            if (logfile) then
               write(11,*)
               write(11,"(A15,ES14.6)") ' New stepsize =',A
            end if

         end if
         

!        ******************************************
!        ***   LOGFILE AND MESSAGES TO SCREEN   ***
!        ******************************************

!        Write messages.

         if (myproc.eq.0) then

            if (veryver) then

               write(*,*)
               write(*,*)
               write(*,"(A16,I3)") ' FLOW ITERATION ',i

               if (intarea.ne.zero) then
                  aux = one/intarea
               else
                  aux = one
               end if

               write(*,*)
               write(*,"(A21,ES14.6)") ' Surface area       =',intarea
               write(*,"(A21,ES14.6)") ' Mean value of H    =',aux*intexp
               write(*,"(A21,ES14.6)") ' Mean value of H^2  =',aux*intexp2
               write(*,"(A21,ES14.6)") ' ahfdiff/ahfsum     =',ahfdiff/ahfsum

               if (offset) then
                  write(*,*)
                  write(*,"(A6,ES14.6)") ' xc  =',xc
                  write(*,"(A6,ES14.6)") ' yc  =',yc
                  write(*,"(A6,ES14.6)") ' zc  =',zc
               end if

               write(*,*)
               write(*,"(A20)") ' Shape coefficients:'
               write(*,*)

               write(*,"(A4,I2,A3,ES14.6)") ' c0(',0,') =',c0(0)

               do l=1+stepz,lmax,1+stepz
                  write(*,"(A4,I2,A3,ES14.6)") ' c0(',l,') =',c0(l)
               end do

            end if

            if (logfile) then

               write(11,*)
               write(11,"(A16,I3)") ' FLOW ITERATION ',i

               if (intarea.ne.zero) then
                  aux = one/intarea
               else
                  aux = one
               end if

               write(11,*)
               write(11,"(A21,ES14.6)") ' Surface area       =',intarea
               write(11,"(A21,ES14.6)") ' Mean value of H    =',aux*intexp
               write(11,"(A21,ES14.6)") ' Mean value of H^2  =',aux*intexp2
               write(11,"(A21,ES14.6)") ' ahfdiff/ahfsum     =',ahfdiff/ahfsum

               if (offset) then
                  write(11,*)
                  write(11,"(A6,ES14.6)") ' xc  =',xc
                  write(11,"(A6,ES14.6)") ' yc  =',yc
                  write(11,"(A6,ES14.6)") ' zc  =',zc
               end if

               write(11,*)
               write(11,"(A20)") ' Shape coefficients:'
               write(11,*)

               write(11,"(A4,I2,A3,ES14.6)") ' c0(',0,') =',c0(0)

               do l=1+stepz,lmax,1+stepz
                  write(11,"(A4,I2,A3,ES14.6)") ' c0(',l,') =',c0(l)
               end do

            end if

            if (nonaxi) then

               if (.not.refy) then

                  if (veryver) then

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

!   10             format(I4,I4,A2,ES14.6,A1,ES14.6)

               else

                  if (veryver) then

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

!   20             format(I4,I4,A2,ES14.6)

               end if
            end if
         end if

!        Close logfile.

         if (logfile.and.(myproc.eq.0)) then
            close(11)
         end if


!        ***********************************
!        ***   CHECK IF WE HAVE FAILED   ***
!        ***********************************

!        The finder needs a way to give up if there is no
!        horizon.  At the moment, the way in which I do this
!        is that I just give up if the monopole term is of
!        the order of 2 grid points in size.  Maybe this
!        could be improved on?

         if (c0(0).lt.2.0D0*dd) then
            write(*,*)
            write(*,*) 'Surface radius is too small, giving up.'
            return
         end if


!        ********************************************
!        ***   CHECK IF WE HAVE FOUND A HORIZON   ***
!        ********************************************

!        For this I use a straigthforward test to check if the
!        change in the surface is below a given tolerance.  This
!        usually works well, but if the tolerance is too small then
!        the convergence becomes very slow.  Also, since I use an
!        adaptive stepsize, I run the risk of having a small change
!        simply because the stepsize itself is very small, and not
!        because we found a horizon.  This can happen is the tolerance
!        is not small enough.  So we have a balancing act here.

         if (ahfdiff.lt.ahftol*ahfsum) then
            status = .true.
            return
         end if

      end do


!     ******************************
!     ***   TOO MANY IERATIONS   ***
!     ******************************

!     If we ever get here, we did too many iterations and
!     failed to converge.

      status = .true.

      if (veryver) then
         write(*,*)
         write(*,*) 'Too many iterations'
      end if


!     ***************
!     ***   END   ***
!     ***************

      end subroutine AHFinder_flow
