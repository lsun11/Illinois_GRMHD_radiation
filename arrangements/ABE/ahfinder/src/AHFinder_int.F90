      /*@@
      @file      AHFinder_int.F
      @date      April 1998
      @author    Miguel Alcubierre
      @desc
      Find surface integrals. The integrals are done
      in parallel, but the number of processors is
      assumed to be either the square of an integer
      or a power of two (if this is not the case,
      I just use less processors).
      @enddesc
      @version   $Header: /cactus/CactusEinstein/AHFinder/src/AHFinder_int.F,v 1.50 2003/11/05 20:31:57 diener Exp $
      @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

      subroutine AHFinder_int(CCTK_ARGUMENTS)

      use AHFinder_dat

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      integer i,j,l,m
      integer npt,npp
      integer l_ntheta,l_nphi,theta0,phiangle0
      integer npoints
      integer auxi,vindex
      integer ierror

      integer interp_handle,coord_system_handle
      integer reduction_handle,param_table_handle
      character(30) options_string
      character(128) operator
      CCTK_INT red_tmp, nchars

      CCTK_POINTER, dimension(3) :: interp_coords
      CCTK_INT, dimension(9) :: in_array_indices
      CCTK_POINTER, dimension(9) :: out_arrays
      CCTK_INT, dimension(9) :: out_array_type_codes


      CCTK_REAL LEGEN
      CCTK_REAL xp,yp,zp,rp
      CCTK_REAL xw,yw,zw
      CCTK_REAL cost,sint,cosp,sinp
      CCTK_REAL trr,ttt,tpp,trt,trp,ttp,ft,fp
      CCTK_REAL dtheta,dphiangle,dtp,idtheta,idphiangle,phianglestart
      CCTK_REAL det,idet
      CCTK_REAL zero,quarter,half,one,two,three,four,pi
      CCTK_REAL sina,cosa,intw,grad,sigma,h0
      CCTK_REAL theta,phiangle
      CCTK_REAL aux,aux1,aux2

      CCTK_REAL tempv(0:lmax),tempm(lmax,lmax)

      CCTK_REAL, dimension(3,3) :: gupij

      CCTK_REAL, allocatable, dimension(:) :: costheta,sintheta
      CCTK_REAL, allocatable, dimension(:) :: cosphi,sinphi

      CCTK_REAL, allocatable, dimension(:,:) :: rr,xa,ya,za,da,exp,gradn
      CCTK_REAL, allocatable, dimension(:,:) :: txx,tyy,tzz,txy,txz,tyz
      CCTK_REAL, allocatable, dimension(:,:) :: intmask

!      CCTK_REAL dX,dY,dZ
      CCTK_REAL, allocatable, dimension(:,:) :: pointcoords
      CCTK_REAL, allocatable, dimension(:) :: txxout,tyyout,tzzout,txyout,txzout,tyzout
      CCTK_REAL, allocatable, dimension(:) :: intmaskout
      CCTK_REAL, allocatable, dimension(:) :: expout

!     Variables to be saved for next call.

      save dtheta,dphiangle,dtp,idtheta,idphiangle,phianglestart
      save npt,npp,l_ntheta,l_nphi,theta0,phiangle0

!     Description of variables:
!     
!     i,j,l,m    Counters.
!     
!     theta      Latitude.
!     phiangle        Longitude.
!     
!     npt        number of processors in theta direction.
!     npp        number of processors in phiangle direction.
!     
!     l_ntheta   Local number of grid points in theta direction.
!     l_nphi     Local number of grid points in phiangle direction.
!     
!     theta0     Local origin in theta direction.
!     phiangle0       Local origin in phiangle direction.
!     
!     npoints    Number of points to interpolate .
!     
!     xp         x coordinate from surface centre.
!     yp         y coordinate from surface centre.
!     zp         z coordinate from surface centre.
!     
!     rp         Radius.
!     rr         Radius array.
!     
!     xa         Array with x coordinates from grid centre.
!     ya         Array with y coordinates from grid centre.
!     za         Array with z coordinates from grid centre.
!     
!     exp        Interpolated expansion.
!     
!     txx        Array with interpolated gxx metric component.
!     tyy        Array with interpolated gyy metric component.
!     tzz        Array with interpolated gzz metric component.
!     txy        Array with interpolated gxy metric component.
!     txz        Array with interpolated gxz metric component.
!     tyz        Array with interpolated gyz metric component.
!     
!     gradn      Array with interpolated norm of gradient of horizon function.
!     
!     intmask    Array with interpolated mask.
!     
!     cost       cos(theta)
!     sint       sin(theta)
!     cosp       cos(phiangle)
!     sinp       sin(phiangle)
!     
!     trr        Metric component {r,r}.
!     ttt        Metric component {theta,theta}.
!     tpp        Metric component {phiangle,phiangle}.
!     trt        Metric component {r,theta}.
!     trp        Metric component {r,phiangle}.
!     ttp        Metric component {theta,phiangle}.
!     
!     ft         dr/dtheta
!     fp         dr/dphiangle
!     
!     dtheta     Grid spacing in theta.
!     dphiangle       Grid spacing in phiangle.
!     dtp        dtheta*dphiangle.
!     
!     idtheta    1/(2 dtheta)
!     idphiangle      1/(2 dphiangle)
!     
!     phianglestart   Origin for phiangle (normally 0, but for some symmetries -pi/2).
!     
!     da         Area element.


!     *******************
!     ***   NUMBERS   ***
!     *******************

      zero    = 0.0D0
      quarter = 0.25D0
      half    = 0.5D0
      one     = 1.0D0
      two     = 2.0D0
      three   = 3.0D0
      four    = 4.0D0

      pi = acos(-one)


!     *****************************************************************
!     ***   TOTAL NUMBER OF PROCESSORS AND LOCAL PROCESSOR NUMBER   ***
!     *****************************************************************

      myproc = CCTK_MyProc(cctkGH)
      nprocs = CCTK_nProcs(cctkGH)


!     ****************************************
!     ***   GET REDUCTION HANDLE FOR SUM   ***
!     ****************************************

      call CCTK_ReductionArrayHandle(reduction_handle,"sum")

      if (reduction_handle.lt.0) then
         call CCTK_WARN(1,"Cannot get handle for sum reduction ! Forgot to activate an implementation providing reduction operators ??")
      end if


!     ******************************************************
!     ***   THINGS TO DO ON FIRST CALL TO THIS ROUTINE   ***
!     ******************************************************

      if (firstint) then

         firstint = .false.


!     ******************************************
!     ***   FIND  phianglestart,dtheta,dphiangle,dtp   ***
!     ******************************************

         phianglestart = zero

         dtheta = pi/dble(ntheta)

         if (cartoon) then

            dphiangle    = zero
            dtp     = dtheta
            idphiangle   = one
            idtheta = half/dtheta

            if (refz.or.CCTK_Equals(domain,"bitant").eq.1) then
               dtp = dtheta*half
               idtheta = half*half/dtheta
!     Following line added by Zach to fix axisymmetry + bitant (i.e., equatorial) symmetry:
               dtheta = half*dtheta
            end if

         else

            dphiangle = two*pi/dble(nphi)

            if (refz) dtheta = half*dtheta

            if (refx.and.refy) then
               dphiangle = quarter*dphiangle
            else if (refx) then
               dphiangle = half*dphiangle
               phianglestart = - half*pi
            else if (refy) then
               dphiangle = half*dphiangle
            end if

            dtp = dtheta*dphiangle

            idtheta = half/dtheta
            idphiangle   = half/dphiangle

         end if


!     *********************************************
!     ***   DIVIDE 2D DOMAIN AMONG PROCESSORS   ***
!     *********************************************

!     For cartoon this is trivial.

         if (cartoon) then

            npt = nprocs
            npp = 1

!     Otherwise it is more complicated.  At the moment I
!     only consider numbers of processors that are either
!     the square of an integer, or a power of two.
!     If neither of this is the case, then I use fewer
!     processors to make sure that I always deal with
!     either a square or a power of two.  This can
!     certainly be improved, but it might not be that
!     urgent since in most cases the number of processors
!     will be a power of two.

         else

            aux = sqrt(dble(nprocs))

!     One processor.

            if (nprocs.eq.1) then

               npt = 1
               npp = 1

!     Square of an integer.

            else if (aux.eq.int(aux)) then

               npt = int(aux)
               npp = int(aux)

!     The number of processors is not the square of an integer.
!     We will now check if it is a power of two.

            else

               l = 1
               m = 1

               i = 1

               do while (l*m.lt.nprocs)

                  npt = l
                  npp = m

                  if (mod(i,2).ne.0) then
                     m = 2*m
                  else
                     l = 2*l
                  end if

                  i = i + 1

               end do

!     "nprocs" is a power of two.

               if (l*m.eq.nprocs) then

                  npt = l
                  npp = m

!     "nprocs" is not a power of two.  This is where I
!     would need to do something clever.  At the moment,
!     I will just use fewer processors.  The number of
!     processors I use is the power of two or the square
!     of an integer that is closest to "nprocs".

               else

                  if (npt*npp.lt.int(aux)**2) then
                     npt = int(aux)
                     npp = int(aux)
                  end if

               end if

            end if

         end if

!     Now be careful not to use too many procs and end up
!     with very few points per proc.

         if (ntheta/npt.lt.3) npt = ntheta/3

         if (.not.cartoon) then
            if (nphi/npp.lt.3) npp = nphi/3
         end if

!     Figure out the number of grid points per processor,
!     and the "origin" for a given processor.

         if (myproc.ge.npt*npp) then

            l_ntheta = 0
            l_nphi   = 0

            theta0 = 0
            phiangle0   = 0

         else

!     First take the number of grid points per processor
!     in a given direction to be equal to the total number
!     of grid points divided by the number of processors
!     in that direction.

            l_ntheta = ntheta/npt
            l_nphi   = nphi/npp

!     And take the "origin" as the corresponding displacement
!     from the real origin.  Notice that mod(myproc,npt)
!     tells my on which "theta" bin I am, while myproc/npt
!     tells me on which "phiangle" bin I am.

            i = mod(myproc,npt)
            j = myproc/npt

            theta0 = i*l_ntheta
            phiangle0   = j*l_nphi

!     Now correct all this if the total number of grid points
!     in a given direction was not exactly divisible by the
!     number of processors on that direction.

            if (l_ntheta*npt.ne.ntheta) then

!     Find residue on theta direction.

               l = ntheta - l_ntheta*npt

!     Distribute residue in first few processors, and
!     displace theta0 accordingly.

               if (i.lt.l) then
                  l_ntheta = l_ntheta + 1
                  theta0 = theta0 + i
               else
                  theta0 = theta0 + l
               end if

            end if

            if (l_nphi*npp.ne.nphi) then

!     Find residue on phiangle direction.

               l = nphi - l_nphi*npp

!     Distribute residue in first few processors, and
!     displace phiangle0 accordingly.

               if (j.lt.l) then
                  l_nphi = l_nphi + 1
                  phiangle0 = phiangle0 + j
               else
                  phiangle0 = phiangle0 + l
               end if

            end if

         end if

      end if


!     **************************************
!     ***   ALLOCATE MEMORY FOR ARRAYS   ***
!     **************************************

      allocate(costheta(1:l_ntheta+1))
      allocate(sintheta(1:l_ntheta+1))

      allocate(cosphi(1:l_nphi+1),sinphi(1:l_nphi+1))

      allocate(rr(1:l_ntheta+1,1:l_nphi+1))

      allocate(xa(1:l_ntheta+1,1:l_nphi+1))
      allocate(ya(1:l_ntheta+1,1:l_nphi+1))
      allocate(za(1:l_ntheta+1,1:l_nphi+1))

      allocate(da(1:l_ntheta+1,1:l_nphi+1))
      allocate(exp(1:l_ntheta+1,1:l_nphi+1))
      allocate(gradn(1:l_ntheta+1,1:l_nphi+1))

      allocate(txx(1:l_ntheta+1,1:l_nphi+1))
      allocate(tyy(1:l_ntheta+1,1:l_nphi+1))
      allocate(tzz(1:l_ntheta+1,1:l_nphi+1))
      allocate(txy(1:l_ntheta+1,1:l_nphi+1))
      allocate(txz(1:l_ntheta+1,1:l_nphi+1))
      allocate(tyz(1:l_ntheta+1,1:l_nphi+1))

      allocate(intmask(1:l_ntheta+1,1:l_nphi+1))

      allocate(pointcoords((l_ntheta+1)*(l_nphi+1),3))
      allocate(txxout((l_ntheta+1)*(l_nphi+1)))
      allocate(txyout((l_ntheta+1)*(l_nphi+1)))
      allocate(txzout((l_ntheta+1)*(l_nphi+1)))
      allocate(tyyout((l_ntheta+1)*(l_nphi+1)))
      allocate(tyzout((l_ntheta+1)*(l_nphi+1)))
      allocate(tzzout((l_ntheta+1)*(l_nphi+1)))
      allocate(expout((l_ntheta+1)*(l_nphi+1)))
      allocate(intmaskout((l_ntheta+1)*(l_nphi+1)))


!     ********************************
!     ***   INITIALIZE VARIABLES   ***
!     ********************************

!     Initialize error flags.

      interror1 = 0
      interror2 = 0
      interror3 = 0

!     Initialize surface integrals.

      intexp  = zero
      intexp2 = zero
      intarea = zero
      intexpdel2 = zero

!     For flow algorithm, initialize spectral components.

      if (flow) then

         hflow0 = zero
         hflowc = zero
         hflows = zero

         cflow0 = zero
         cflowc = zero
         cflows = zero

         nflow0 = zero
         nflowc = zero
         nflows = zero

      end if


!     *************************************
!     ***   FIND INTERPOLATING POINTS   ***
!     *************************************

      if (myproc.ge.npt*npp) then

         xa = half*(xmx+xmn)
         ya = half*(ymx+ymn)
         za = half*(zmx+zmn)

      else

         if (refz.or.CCTK_Equals(domain,"bitant").eq.1) then
                                !Zach again... fix l_ntheta so that the interpolation stencil stays in-bounds!
!     if(abs(cos(dtheta*dble(l_ntheta+theta0))-0.D0) .lt. 0.0001) then
!     l_ntheta = l_ntheta - 1
!     write(*,*) "Fixed l_ntheta."
!     end if
         end if



         do j=1,l_nphi+1

!     Find phiangle.

            phiangle = dphiangle*dble(j-1+phiangle0) + phianglestart

!     Find sines and cosines of phiangle.

            cosp = cos(phiangle)
            sinp = sin(phiangle)

            cosphi(j) = cosp
            sinphi(j) = sinp

            do i=1,l_ntheta+1

!     Find theta.

               theta = dtheta*dble(i-1+theta0)

!     The following fixes theta so that when theta=0, we don't get interpolation errors with parameter Interpolation_Order = 2.
               if(interpolation_order .gt. 1) then
                  if(abs(theta) .lt. dtheta*0.00001) theta = theta+dtheta*0.00001
               end if

!     end if

!     Find sines and cosines of theta.

               cost = cos(theta)
               sint = sin(theta)

               costheta(i) = cost
               sintheta(i) = sint

!     Find radius rp.

               rp = c0(0)

               do l=1+stepz,lmax,1+stepz
                  rp = rp + c0(l)*LEGEN(l,0,cost)
               end do

!     Notice how the sum over m is first.  This will allow
!     me to use the recursion relations to avoid having to
!     start from scratch every time.  Also, I sum over all
!     l even if I do not want some terms.  This is because
!     in order to use the recursion relations I need all
!     polynomials.

               if (nonaxi) then
                  do m=1,lmax
                     aux = dble(m)*phiangle
                     sina = sin(aux)
                     cosa = cos(aux)
                     do l=m,lmax
                        aux = LEGEN(l,m,cost)
                        rp = rp + aux*cc(l,m)*cosa
                        if (.not.refy) then
                           rp = rp + aux*cs(l,m)*sina
                           
                        end if
                     end do
                  end do
               end if

!     Check for negative radius.

               if (rp.le.zero) then
                  interror1 = 1
                  write(6,*)'interror1!',rp
               end if

!     Find cartesian coordinates.

               aux = rp*sint

               xp = aux*cosp + xc
               yp = aux*sinp + yc
               zp = rp*cost + zc

!     Save arrays.

               rr(i,j) = rp

               xa(i,j) = xp
               ya(i,j) = yp
               za(i,j) = zp

!     Check if we are within bounds.

               if ((xp.gt.xmx).or.(xp.lt.xmn).or. &
                    (yp.gt.ymx).or.(yp.lt.ymn).or. &
                    (zp.gt.zmx).or.(zp.lt.zmn)) then
                  interror2 = 1
                  write(*,*) "interror",theta,theta0,rp,xp,yp,zp,xmx,ymx,zmx,&
                       rp,cost,zc,rp*cost

!     $                 xmn,ymn,zmn,theta,xc,yc,zc,aux,
!     $                 sint,cost,sinp,cosp,interror1
               end if
!     write(*,*) "nointerror",theta,theta0,xp,yp,zp,xmx,ymx,zmx, 
!     $                 xmn,ymn,zmn,theta,xc,yc,zc,aux,rp,
!     $                 sint,cost,sinp,cosp,interror1

            end do
         end do

      end if


!     *******************************************
!     ***   REDUCE ERRORS ACROSS PROCESSORS   ***
!     *******************************************

      call CCTK_ReduceLocalScalar(ierror,cctkGH,-1,reduction_handle, &
           interror1,red_tmp,CCTK_VARIABLE_INT)
      if (ierror.ne.0) then
         call CCTK_WARN(1,"Reduction of norm failed !")
      end if
      interror1 = red_tmp

      call CCTK_ReduceLocalScalar(ierror,cctkGH,-1,reduction_handle, &
           interror2,red_tmp,CCTK_VARIABLE_INT)
      if (ierror.ne.0) then
         call CCTK_WARN(1,"Reduction of norm failed !")
      end if
      interror2 = red_tmp

!     If there was an error on any processor then assign
!     large values to the integrals and return from the
!     subroutine (but remember to deallocate the arrays
!     first!).

      interror = .false.

      if ((interror1.ne.0).or.(interror2.ne.0)) then
         intexp  = 1.0D10
         intexp2 = 1.0D10
         intarea = 1.0D10
         intexpdel2 = 1.0D10
         interror = .true.
         goto 10
      end if


!     ***********************
!     ***   INTERPOLATE   ***
!     ***********************

!     Number of points to interpolate.

      npoints = (l_ntheta+1)*(l_nphi+1)

!     parameter, local interpolator, and coordinate system handle.
      param_table_handle = -1
      interp_handle = -1
      coord_system_handle = -1

      options_string = "order = " // char(ichar('0') + interpolation_order)
      call Util_TableCreateFromString (param_table_handle, options_string)
      if (param_table_handle .lt. 0) then
         call CCTK_WARN(0,"Cannot create parameter table for interpolator")
      endif

      call CCTK_FortranString (nchars, interpolation_operator, operator)
      call CCTK_InterpHandle (interp_handle, operator)
      if (interp_handle .lt. 0) then
         call CCTK_WARN(0,"Cannot get handle for interpolation ! Forgot to activate an implementation providing interpolation operators ??")
      endif

      call CCTK_CoordSystemHandle (coord_system_handle, "cart3d")
      if (coord_system_handle .lt. 0) then
         call CCTK_WARN(0,"Cannot get handle for cart3d coordinate system ! Forgot to activate an implementation providing coordinates ??")
      endif

!     fill in the input/output arrays for the interpolator
      interp_coords(1) = CCTK_PointerTo(xa)
      interp_coords(2) = CCTK_PointerTo(ya)
      interp_coords(3) = CCTK_PointerTo(za)

      call CCTK_VarIndex (vindex, "bssn::gxx")
      in_array_indices(1) = vindex
      call CCTK_VarIndex (vindex, "bssn::gyy")
      in_array_indices(2) = vindex
      call CCTK_VarIndex (vindex, "bssn::gzz")
      in_array_indices(3) = vindex
      call CCTK_VarIndex (vindex, "bssn::gxy")
      in_array_indices(4) = vindex
      call CCTK_VarIndex (vindex, "bssn::gxz")
      in_array_indices(5) = vindex
      call CCTK_VarIndex (vindex, "bssn::gyz")
      in_array_indices(6) = vindex
      call CCTK_VarIndex (vindex, "ahfinder::ahf_exp")
      in_array_indices(7) = vindex
      call CCTK_VarIndex (vindex, "ahfinder::ahmask")
      in_array_indices(8) = vindex
      call CCTK_VarIndex (vindex, "ahfinder::ahfgradn")
      in_array_indices(9) = vindex

      out_arrays(1) = CCTK_PointerTo(txx)
      out_arrays(2) = CCTK_PointerTo(tyy)
      out_arrays(3) = CCTK_PointerTo(tzz)
      out_arrays(4) = CCTK_PointerTo(txy)
      out_arrays(5) = CCTK_PointerTo(txz)
      out_arrays(6) = CCTK_PointerTo(tyz)
      out_arrays(7) = CCTK_PointerTo(exp)
      out_arrays(8) = CCTK_PointerTo(intmask)
      out_arrays(9) = CCTK_PointerTo(gradn)

      out_array_type_codes = CCTK_VARIABLE_REAL
      
!     write(*,*) "hewwo ahfinder_int1: ",exp(1  ,1  ),ahf_exp
!     write(*,*) "hewwo ahfinder_int2: ",exp(1  ,1  ),ahmask
!     write(*,*) "hewwo ahfinder_int3: ",exp(1  ,1  ),ahfgradn


!     For minimization we need the interpolated metric,
!     the expansion and the mask.

      if (.not.flow) then
!     When following if statement is true, use interpolate_pointset global interpolator instead of CCTK_InterpGridArrays "local" interpolator
         if(1==0) then
            dX = CCTK_DELTA_SPACE(1)
            dY = CCTK_DELTA_SPACE(2)
            dZ = CCTK_DELTA_SPACE(3)

            do j=1,l_nphi+1
               do i=1,l_ntheta+1         
                  pointcoords((j-1)*(l_ntheta+1) + i,1) = xa(i,j)
                  pointcoords((j-1)*(l_ntheta+1) + i,2) = ya(i,j)
                  pointcoords((j-1)*(l_ntheta+1) + i,3) = za(i,j)
               end do
            end do

            call interpolate_pointset(cctkGH,cctk_nghostzones,npoints,Interpolation_Order,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gxx,txxout)
            call interpolate_pointset(cctkGH,cctk_nghostzones,npoints,Interpolation_Order,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gxy,txyout)
            call interpolate_pointset(cctkGH,cctk_nghostzones,npoints,Interpolation_Order,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gxz,txzout)
            call interpolate_pointset(cctkGH,cctk_nghostzones,npoints,Interpolation_Order,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gyy,tyyout)
            call interpolate_pointset(cctkGH,cctk_nghostzones,npoints,Interpolation_Order,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gyz,tyzout)
            call interpolate_pointset(cctkGH,cctk_nghostzones,npoints,Interpolation_Order,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gzz,tzzout)
            call interpolate_pointset(cctkGH,cctk_nghostzones,npoints,Interpolation_Order,cctk_lsh,cctk_gsh,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,ahf_exp,expout)
            intmask = 1.D0

            do j=1,l_nphi+1
               do i=1,l_ntheta+1
                  txx(i,j) = txxout((j-1)*(l_ntheta+1) + i)
                  txy(i,j) = txyout((j-1)*(l_ntheta+1) + i)
                  txz(i,j) = txzout((j-1)*(l_ntheta+1) + i)
                  tyy(i,j) = tyyout((j-1)*(l_ntheta+1) + i)
                  tyz(i,j) = tyzout((j-1)*(l_ntheta+1) + i)
                  tzz(i,j) = tzzout((j-1)*(l_ntheta+1) + i)
                  exp(i,j) = expout((j-1)*(l_ntheta+1) + i)
                                !intmask(i,j) = intmaskout((j-1)*(l_ntheta+1) + i)

!     write(*,*) "hello0",xa(i,j),ya(i,j),za(i,j)
!     write(*,*) "hello",i,j,txx(i,j),txy(i,j),txz(i,j),tyy(i,j),tyz(i,j),tzz(i,j),exp(i,j),intmask(i,j)
!     write(*,*) "hello",i,j,npoints,pointcoords((j-1)*(l_ntheta+1) + i,2)-ya(i,j)
!     write(*,*) "hello",i,j,npoints,pointcoords((j-1)*(l_ntheta+1) + i,3)-za(i,j)
               end do
            end do
         else
            call CCTK_InterpGridArrays (ierror, cctkGH, 3, interp_handle, &
                 param_table_handle, coord_system_handle, &
                 npoints, CCTK_VARIABLE_REAL, interp_coords, &
                 8, in_array_indices, &
                 8, out_array_type_codes, out_arrays)

         end if


!     do j=1,l_nphi+1
!     do i=1,l_ntheta+1
!     write(*,*) "hello0",xa(i,j),ya(i,j),za(i,j)
!     write(*,*) "hello",i,j,txx(i,j),txy(i,j),txz(i,j),tyy(i,j),tyz(i,j),tzz(i,j),exp(i,j),intmask(i,j)
!     end do
!     end do
         if (ierror < 0) then
            call CCTK_WARN (1, "AHFinder_int.F1: interpolator call returned an error code");
         endif


!     For N flow, we need the interpolated metric, the expansion,
!     the norm of the gradient of the horizon function, and the mask.

      else if (nw.ne.zero) then

         call CCTK_InterpGridArrays (ierror, cctkGH, 3, interp_handle, &
              param_table_handle, coord_system_handle, &
              npoints, CCTK_VARIABLE_REAL, interp_coords, &
              9, in_array_indices, &
              9, out_array_type_codes, out_arrays)
         if (ierror < 0) then
            call CCTK_WARN (1, "AHFinder_int.F2: interpolator call returned an error code");
         endif

!     For H or C flows, we need the interpolated expansion, the
!     norm of the gradient of the horizon function, and the mask.

      else

         call CCTK_InterpGridArrays (ierror, cctkGH, 3, interp_handle, &
              param_table_handle, coord_system_handle, &
              npoints, CCTK_VARIABLE_REAL, interp_coords, &
              3, in_array_indices(7), &
              3, out_array_type_codes(7), out_arrays(7))
         if (ierror < 0) then
            call CCTK_WARN (1, "AHFinder_int.F3: interpolator call returned an error code");
         endif

      end if

!     release parameter table
      call Util_TableDestroy (ierror, param_table_handle)

!     ***************************************
!     ***   CHECK IF WE ARE INSIDE MASK   ***
!     ***************************************

!     Check if we are either inside mask, or to close for comfort.
!     The way I do this is very simple, I just check the interpolated
!     values of the mask.  If any of these values is not 1.0, it means
!     that a point inside the mask is contaminating the interpolation,
!     and this means in turn that we are too close to the mask.

!     Notice that in practice I do not ask for the interpolated values
!     to be exactly 1.0, since round off error can modify the value
!     slightly without really meaning that we are close to the mask.

      if (.not.CCTK_EQUALS(ahf_mask,'off')) then

         if (myproc.lt.npt*npp) then
            do j=1,l_nphi+1
               do i=1,l_ntheta+1
                  if (intmask(i,j).lt.0.99D0) then
                     interror3 = 1
                  end if
               end do
            end do
         end if

         call CCTK_ReduceLocalScalar(ierror,cctkGH,-1,reduction_handle, &
              interror3,red_tmp,CCTK_VARIABLE_INT)
         interror3 = red_tmp

         if (interror3.ne.0) then
            intexp  = 1.0D10
            intexp2 = 1.0D10
            intarea = 1.0D10
            intexpdel2 = 1.0D10
            interror = .true.
            goto 10
         end if

      end if


!     *************************
!     ***   AREA ELEMENTS   ***
!     *************************

!     Find 2D array for area elements.  Notice that
!     this is not needed for the flow algorithm!

      if (.not.flow) then

         if (myproc.ge.npt*npp) then

            da = zero

         else

            do j=1,l_nphi+1

!     Find sines and cosines of phiangle.

               cosp = cosphi(j)
               sinp = sinphi(j)

               do i=1,l_ntheta+1

!     Find {rp}.

                  rp = rr(i,j)

!     Find sines and cosines of theta.

                  cost = costheta(i)
                  sint = sintheta(i)

!     Find spherical metric components.

                  trr = sint**2*(txx(i,j)*cosp**2 &
                       + tyy(i,j)*sinp**2) + tzz(i,j)*cost**2 &
                       + two*sint*(txy(i,j)*sint*cosp*sinp &
                       + cost*(txz(i,j)*cosp + tyz(i,j)*sinp))

                  ttt = rp**2*(cost**2*(txx(i,j)*cosp**2 &
                       + tyy(i,j)*sinp**2) + tzz(i,j)*sint**2 &
                       + two*cost*(txy(i,j)*cost*cosp*sinp &
                       - sint*(txz(i,j)*cosp + tyz(i,j)*sinp)))

                  tpp = (rp*sint)**2*(txx(i,j)*sinp**2 + tyy(i,j)*cosp**2 &
                       - two*txy(i,j)*cosp*sinp)

                  trt = rp*(sint*cost*(txx(i,j)*cosp**2 + tyy(i,j)*sinp**2 &
                       - tzz(i,j) + two*txy(i,j)*sinp*cosp) &
                       + (cost**2-sint**2)*(txz(i,j)*cosp + tyz(i,j)*sinp))

                  trp = rp*sint*(sint*(cosp*sinp*(tyy(i,j) - txx(i,j)) &
                       + txy(i,j)*(cosp**2 - sinp**2)) &
                       + cost*(cosp*tyz(i,j) - sinp*txz(i,j)))

                  ttp = rp**2*sint*(cost*(sinp*cosp*(tyy(i,j) - txx(i,j)) &
                       + (cosp**2 - sinp**2)*txy(i,j)) &
                       + sint*(sinp*txz(i,j) - cosp*tyz(i,j)))

!     Find derivatives  {ft,fp}.  For interior points I
!     use centered differences, and for the boundary points
!     I use second order one sided differences. Notice that
!     since this is done even in inter-processor boundaries,
!     the final result will depend on the number of processors,
!     but the differences will converge away at high order.

                  if ((i.ne.1).and.(i.ne.l_ntheta+1)) then
                     ft = idtheta*(rr(i+1,j) - rr(i-1,j))
                  else if (i.eq.1) then
                     ft = - idtheta*(three*rr(1,j) &
                          - four*rr(2,j) + rr(3,j))
                  else
                     ft = + idtheta*(three*rr(l_ntheta+1,j) &
                          - four*rr(l_ntheta,j) + rr(l_ntheta-1,j))
                  end if

                  if (nonaxi) then
                     if ((j.ne.1).and.(j.ne.l_nphi+1)) then
                        fp = idphiangle*(rr(i,j+1) - rr(i,j-1))
                     else if (j.eq.1) then
                        fp = - idphiangle*(three*rr(i,1) &
                             - four*rr(i,2) + rr(i,3))
                     else
                        fp = + idphiangle*(three*rr(i,l_nphi+1) &
                             - four*rr(i,l_nphi) + rr(i,l_nphi-1))
                     end if
                  else
                     fp = zero
                  end if

!     Find area element.  For this I use the fact that if
!     we have  r = f(theta,phiangle), then the area element is:
!     
!     /                  2
!     da  =  | [ gtt + grr (f,t)  +  2 grt f,t ]
!     \
!     2
!     [ gpp + grr (f,p)  +  2 grp f,p ]  -  [ gtp
!     
!     2 \ 1/2
!     + grr f,t f,p  +  grt f,p + grp f,t ]  |     dtheta dphiangle
!     /
!     
!     where  f,t = df/dtheta  and  f,p = df/dphiangle.

                  aux = (ttt + trr*ft**2 + two*trt*ft) &
                       * (tpp + trr*fp**2 + two*trp*fp) &
                       - (ttp + ft*(trr*fp + trp) + trt*fp)**2

                  if (aux.gt.zero) then
                     aux = sqrt(aux)
                  else
                     aux = zero
                  end if

!     Multiply with dtheta*dphiangle.

                  da(i,j) = aux*dtp

               end do
            end do

         end if

      end if


!     *******************************************
!     ***   AREA AND INTEGRALS OF EXPANSION   ***
!     *******************************************

!     To find the integrals I use a second order method:
!     I approximate the value of the integrand at each
!     cell centre by an average of the values at the
!     four cell corners, and then add up all cells.

!     Notice that for the flow algorithm I do not need to
!     calculate the area nor the integrals of the expansion
!     except at the very end, but by then the flag 'flow' has
!     been switched off.

      inside_min_count     = 0
      inside_min_neg_count = 0

      if (.not.flow) then

         if (myproc.lt.npt*npp) then

            do j=1,l_nphi
               do i=1,l_ntheta

                  intarea = intarea + quarter &
                       *(da(i,j  ) + da(i+1,j  ) &
                       + da(i,j+1) + da(i+1,j+1))

                  intexp  = intexp  + quarter &
                       *(da(i  ,j  )*exp(i  ,j  ) &
                       + da(i+1,j  )*exp(i+1,j  ) &
                       + da(i  ,j+1)*exp(i  ,j+1) &
                       + da(i+1,j+1)*exp(i+1,j+1))
                  intexp2 = intexp2 + quarter &
                       *(da(i  ,j  )*exp(i  ,j  )**2 &
                       + da(i+1,j  )*exp(i+1,j  )**2 &
                       + da(i  ,j+1)*exp(i  ,j+1)**2 &
                       + da(i+1,j+1)*exp(i+1,j+1)**2)

                  intexpdel2 = intexpdel2 + quarter &
                       *(da(i  ,j  )*(exp(i  ,j  )+trapped_surface_delta)**2 &
                       + da(i+1,j  )*(exp(i+1,j  )+trapped_surface_delta)**2 &
                       + da(i  ,j+1)*(exp(i  ,j+1)+trapped_surface_delta)**2 &
                       + da(i+1,j+1)*(exp(i+1,j+1)+trapped_surface_delta)**2)

                  inside_min_count = inside_min_count + 1

!     write(*,*) "ahfinder_int: ",i,j,exp(i  ,j  )
!     write(*,*) "ahfinder_int: ",i,j,exp(i+1,j  )
!     write(*,*) "ahfinder_int: ",i,j,exp(i  ,j+1)
!     write(*,*) "ahfinder_int: ",i,j,exp(i+1,j+1)
                  if ((exp(i  ,j  ).lt.0.0D0).and. &
                       (exp(i+1,j  ).lt.0.0D0).and. &
                       (exp(i  ,j+1).lt.0.0D0).and. &
                       (exp(i+1,j+1).lt.0.0D0)) then
                     inside_min_neg_count = inside_min_neg_count + 1
                  end if

               end do
            end do

         end if

!     Reduce integrals.

         call CCTK_ReduceLocalScalar(ierror,cctkGH,-1,reduction_handle, &
              intarea,aux,CCTK_VARIABLE_REAL)
         intarea=aux

         call CCTK_ReduceLocalScalar(ierror,cctkGH,-1,reduction_handle, &
              intexp,aux,CCTK_VARIABLE_REAL)
         intexp=aux

         call CCTK_ReduceLocalScalar(ierror,cctkGH,-1,reduction_handle, &
              intexp2,aux,CCTK_VARIABLE_REAL)
         intexp2=aux

         call CCTK_ReduceLocalScalar(ierror,cctkGH,-1,reduction_handle, &
              intexpdel2,aux,CCTK_VARIABLE_REAL)
         intexpdel2=aux

         call CCTK_ReduceLocalScalar(ierror,cctkGH,-1,reduction_handle, &
              inside_min_neg_count,red_tmp,CCTK_VARIABLE_INT)
         inside_min_neg_count=red_tmp

         call CCTK_ReduceLocalScalar(ierror,cctkGH,-1,reduction_handle, &
              inside_min_count,red_tmp,CCTK_VARIABLE_INT)
         inside_min_count=red_tmp

      end if


!     *******************************
!     ***   SPECTRAL COMPONENTS   ***
!     *******************************

!     For flow algorithm, find spectral components.
!     For this I need to recalculate the Legendre polynomials,
!     but I see no way around it.  Also, to avoid having to
!     define lots of 2D arrays I do the integrals inside
!     the loop.  This means that I have to take care to see
!     if I am on an edge or corner point. Notice also that
!     these integrals do not involve the area element.

      if (flow) then

         if (myproc.lt.npt*npp) then

!     Find h0.

            if (find_trapped_surface) then
               h0 = - trapped_surface_delta
            else
               h0 = zero
            end if

!     Find integrals.

            do j=1,l_nphi+1

!     Find phiangle.

               phiangle = dphiangle*dble(j+phiangle0) + phianglestart

               do i=1,l_ntheta+1

!     Find sines and cosines of theta.

                  cost = costheta(i)
                  sint = sintheta(i)

!     Find norm of gradient of horizon function.

                  grad = gradn(i,j)

!     Find weight factor.

                  intw = dtp*(exp(i,j)-h0)*sint**2

!     Modify weight factor at corners and edges.

                  if (((j.eq.1).or.(j.eq.l_nphi+1)).and. &
                       ((i.eq.1).or.(i.eq.l_ntheta+1))) then
                     intw = quarter*intw
                  else if (((j.eq.1).or.(j.eq.l_nphi+1)).or. &
                       ((i.eq.1).or.(i.eq.l_ntheta+1))) then
                     intw = half*intw
                  end if

!     Find "sigma" for Nakamura flow (see Carstens paper).

                  if (nw.eq.0.0) then

                     sigma = zero

                  else

                     xw = xa(i,j) - xc
                     yw = ya(i,j) - yc
                     zw = za(i,j) - zc

!     Calculate metric on surface.

                     det = txx(i,j)*tyy(i,j)*tzz(i,j) &
                          + two*txy(i,j)*txz(i,j)*tyz(i,j) &
                          - txx(i,j)*tyz(i,j)**2 - tyy(i,j)*txz(i,j)**2 &
                          - tzz(i,j)*txy(i,j)**2

                     idet = one/det

                     gupij(1,1) = idet*(tyy(i,j)*tzz(i,j)-tyz(i,j)**2)
                     gupij(2,2) = idet*(txx(i,j)*tzz(i,j)-txz(i,j)**2)
                     gupij(3,3) = idet*(txx(i,j)*tyy(i,j)-txy(i,j)**2)

                     gupij(1,2) = idet*(txz(i,j)*tyz(i,j)-tzz(i,j)*txy(i,j))
                     gupij(2,1) = gupij(1,2)
                     gupij(1,3) = idet*(txy(i,j)*tyz(i,j)-tyy(i,j)*txz(i,j))
                     gupij(3,1) = gupij(1,3)
                     gupij(2,3) = idet*(txy(i,j)*txz(i,j)-txx(i,j)*tyz(i,j))
                     gupij(3,2) = gupij(2,3)

                     call AHFinder_calcsigma(CCTK_ARGUMENTS, &
                          xw,yw,zw,gupij,sigma)

                  end if

!     Monopole components.

                  aux = intw*grad

                  hflow0(0) = hflow0(0) + intw
                  cflow0(0) = cflow0(0) + aux
                  nflow0(0) = nflow0(0) + aux*sigma

!     Axisymmetric components.

                  do l=1+stepz,lmax,1+stepz

                     aux  = intw*LEGEN(l,0,cost)
                     aux1 = aux*grad

                     hflow0(l) = hflow0(l) + aux
                     cflow0(l) = cflow0(l) + aux1
                     nflow0(l) = nflow0(l) + aux1*sigma

                  end do

!     Non-axisymmetric components.

                  if (nonaxi) then

                     do m=1,lmax

                        aux = dble(m)*phiangle
                        sina = sin(aux)
                        cosa = cos(aux)

                        do l=m,lmax

                           aux  = intw*LEGEN(l,m,cost)

                           aux1 = aux*cosa
                           aux2 = aux1*grad
                           hflowc(l,m) = hflowc(l,m) + aux1
                           cflowc(l,m) = cflowc(l,m) + aux2
                           nflowc(l,m) = nflowc(l,m) + aux2*sigma

                           if (.not.refy) then
                              aux1 = aux*sina
                              aux2 = aux1*grad
                              hflows(l,m) = hflows(l,m) + aux1
                              cflows(l,m) = cflows(l,m) + aux2
                              nflows(l,m) = nflows(l,m) + aux2*sigma
                           end if

                        end do
                     end do

                  end if

               end do
            end do

         end if

!     Reduce integrals.

!     Axisymmetric components.

         auxi = lmax+1

         if (hw.ne.zero) then
            call CCTK_ReduceLocalArray1D(ierror,cctkGH,-1,reduction_handle, &
                 hflow0,tempv,auxi,CCTK_VARIABLE_REAL)
            hflow0 = tempv
         end if

         if (cw.ne.zero) then
            call CCTK_ReduceLocalArray1D(ierror,cctkGH,-1,reduction_handle, &
                 cflow0,tempv,auxi,CCTK_VARIABLE_REAL)
            cflow0 = tempv
         end if

         if (nw.ne.zero) then
            call CCTK_ReduceLocalArray1D(ierror,cctkGH,-1,reduction_handle, &
                 nflow0,tempv,auxi,CCTK_VARIABLE_REAL)
            nflow0 = tempv
         end if

!     Non-axisymmetric components.

         if (nonaxi) then

            auxi = lmax*lmax

            if (hw.ne.zero) then
               call CCTK_ReduceLocalArray1D(ierror,cctkGH,-1, &
                    reduction_handle,hflowc,tempm,auxi,CCTK_VARIABLE_REAL)
               hflowc = tempm
            end if

            if (cw.ne.zero) then
               call CCTK_ReduceLocalArray1D(ierror,cctkGH,-1, &
                    reduction_handle,cflowc,tempm,auxi,CCTK_VARIABLE_REAL)
               cflowc = tempm
            end if

            if (nw.ne.zero) then
               call CCTK_ReduceLocalArray1D(ierror,cctkGH,-1, &
                    reduction_handle,nflowc,tempm,auxi,CCTK_VARIABLE_REAL)
               nflowc = tempm
            end if

            if (.not.refy) then

               if (hw.ne.zero) then
                  call CCTK_ReduceLocalArray1D(ierror,cctkGH,-1, &
                       reduction_handle,hflows,tempm,auxi,CCTK_VARIABLE_REAL)
                  hflows = tempm
               end if

               if (cw.ne.zero) then
                  call CCTK_ReduceLocalArray1D(ierror,cctkGH,-1, &
                       reduction_handle,cflows,tempm,auxi,CCTK_VARIABLE_REAL)
                  cflows = tempm
               end if

               if (nw.ne.zero) then
                  call CCTK_ReduceLocalArray1D(ierror,cctkGH,-1, &
                       reduction_handle,nflows,tempm,auxi,CCTK_VARIABLE_REAL)
                  nflows = tempm
               end if

            end if

         end if

      end if


!     *****************************************************
!     ***   RESCALE INTEGRALS ACCORDING TO SYMMETRIES   ***
!     *****************************************************

!     For cartoon multiply the integrals with 2*pi.

      if (cartoon) then

         intw = 6.283185307D0

         if (refz.or.CCTK_Equals(domain,"bitant").eq.1) then
            intw = intw*2.D0
         end if

         if (.not.flow) then
            intexp  = intw*intexp
            intexp2 = intw*intexp2
            intarea = intw*intarea
            intexpdel2 = intw*intexpdel2
         end if

         if (flow) then
            hflow0 = intw*hflow0
            cflow0 = intw*cflow0
            nflow0 = intw*nflow0
         end if

!     Other symmetries.

      else if (refx.or.refy.or.refz) then

         intw = one

         if (refx) intw = two*intw
         if (refy) intw = two*intw
         if (refz) intw = two*intw

         if (.not.flow) then
            intexp  = intw*intexp
            intexp2 = intw*intexp2
            intarea = intw*intarea
            intexpdel2 = intw*intexpdel2
         end if

!     Spectral components for flow.

         if (flow) then

!     Axisymmetric.

            hflow0 = intw*hflow0
            cflow0 = intw*cflow0
            nflow0 = intw*nflow0

!     Non-axisymmetric.

            if (nonaxi) then

               hflowc = intw*hflowc
               cflowc = intw*cflowc
               nflowc = intw*nflowc

               if (.not.refy) then
                  hflows = intw*hflows
                  cflows = intw*cflows
                  nflows = intw*nflows
               end if

            end if

         end if

      end if


!     *****************************
!     ***   DEALLOCATE MEMORY   ***
!     *****************************

 10   continue

      deallocate(costheta,sintheta)
      deallocate(cosphi,sinphi)
      deallocate(rr,xa,ya,za)
      deallocate(da,exp,gradn)
      deallocate(txx,tyy,tzz,txy,txz,tyz)
      deallocate(intmask)
      deallocate(pointcoords)
      deallocate(txxout)
      deallocate(txyout)
      deallocate(txzout)
      deallocate(tyyout)
      deallocate(tyzout)
      deallocate(tzzout)
      deallocate(expout)
      deallocate(intmaskout)


!     ***************
!     ***   END   ***
!     ***************

      end subroutine AHFinder_int
