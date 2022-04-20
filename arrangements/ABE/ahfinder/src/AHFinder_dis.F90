/*@@
  @file      AHFinder_dis.F
  @date      September 2001
  @author    Miguel Alcubierre
  @desc
             Routine to calculate proper distances between horizons.
  @enddesc
@@*/

#include "cctk.h"

#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

      subroutine AHFinder_dis(CCTK_ARGUMENTS)

      use AHFinder_dat

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      integer i,l,m
      integer np,npoints,vindex,ierror
      integer param_table_handle,interp_handle,coord_system_handle
      character(30) options_string
      character(128) operator
      CCTK_INT nchars

      CCTK_POINTER, dimension(3) :: interp_coords
      CCTK_INT, dimension(6) :: in_array_indices
      CCTK_POINTER, dimension(6) :: out_arrays
      CCTK_INT, dimension(6) :: out_array_type_codes

      CCTK_REAL LEGEN
      CCTK_REAL zero,half,one,two,pi
      CCTK_REAL x12,y12,z12,r12,rho12
      CCTK_REAL theta,phiangle,cost,sint,cosp,sinp
      CCTK_REAL x1,y1,z1,r1
      CCTK_REAL x2,y2,z2,r2
      CCTK_REAL ddx,ddy,ddz,d12
      CCTK_REAL aux,cosa,sina,factor

      CCTK_REAL, allocatable, dimension(:) :: xa,ya,za
      CCTK_REAL, allocatable, dimension(:) :: txx,tyy,tzz,txy,txz,tyz

      character(200) disf
      logical disf_exists

!     **************************
!     ***   DEFINE NUMBERS   ***
!     **************************

      zero  = 0.0D0
      half  = 0.5D0
      one   = 1.0D0
      two   = 2.0D0

      pi = acos(-one)

      factor = one


!     ****************************************
!     ***   FIND LINE CONNECTING CENTERS   ***
!     ****************************************

!     For the moment, horizons 1 and 2 only.

      write(*,*)
      write(*,*) 'Integrating distances between horizons 1 and 2 ...'

      x1 = ahf_xc_1
      y1 = ahf_yc_1
      z1 = ahf_zc_1

      x2 = ahf_xc_2
      y2 = ahf_yc_2
      z2 = ahf_zc_2

!     Find line connecting centers.

!     Neither horizon exists.

      if ((.not.status_old_1).and.(.not.status_old_2)) then

         write(*,*)
         write(*,*) 'Neither horizon 1 nor horizon 2 exist,'
         write(*,*) 'giving up.'
         write(*,*)

         d12 = zero
         goto 100

!     Just horizon 1 exists.

      else if (status_old_1.and.(.not.status_old_2)) then

         write(*,*)
         write(*,*) 'Only horizon 1 exists.'
         write(*,*)

!        Check symmetries.

         if (CCTK_Equals(domain,"octant").eq.1) then

            write(*,*) 'Octant symmetry: calculating twice'
            write(*,*) 'twice the distance from horizon 1 to origin.'

            factor = two

            x2 = zero
            y2 = zero
            z2 = zero

         else if (CCTK_Equals(domain,"quadrant").eq.1) then

            if (CCTK_Equals(quadrant_direction,"x").eq.1) then

               write(*,*) 'Quadrant symmetry, x direction: calculating'
               write(*,*) 'twice the distance from horizon 1 to x axis.'

               factor = two

               x2 = x1
               y2 = zero
               z2 = zero

            else if (CCTK_Equals(quadrant_direction,"y").eq.1) then

               write(*,*) 'Quadrant symmetry, y direction: calculating'
               write(*,*) 'twice the distance from horizon 1 to y axis.'

               factor = two

               x2 = zero
               y2 = y1
               z2 = zero

            else if (CCTK_Equals(quadrant_direction,"z").eq.1) then

               write(*,*) 'Quadrant symmetry, z direction: calculating'
               write(*,*) 'twice the distance from horizon 1 to z axis.'

               factor = two

               x2 = zero
               y2 = zero
               z2 = z1

            end if

         else if (CCTK_Equals(domain,"bitant").eq.1) then

            if (CCTK_Equals(bitant_plane,"xy").eq.1) then

               write(*,*) 'Bitant symmetry, xy plane: calculating twice'
               write(*,*) 'the distance from horizon 1 to xy plane.'

               factor = two

               x2 = x1
               y2 = y1
               z2 = zero

            else if (CCTK_Equals(bitant_plane,"xz").eq.1) then

               write(*,*) 'Bitant symmetry, xz plane: calculating twice'
               write(*,*) 'the distance from horizon 1 to xz plane.'

               factor = two

               x2 = x1
               y2 = zero
               z2 = z1

            else if (CCTK_Equals(bitant_plane,"yz").eq.1) then

               write(*,*) 'Bitant symmetry, yz plane: calculating twice'
               write(*,*) 'the distance from horizon 1 to yz plane.'

               factor = two

               x2 = zero
               y2 = y1
               z2 = z1

            end if

         else

            write(*,*) 'No grid symmetries: no distance to calculate,'
            write(*,*) 'giving up.'

            d12 = zero
            goto 100

         end if

!     Just horizon 2 exists.

      else if (status_old_2.and.(.not.status_old_1)) then

         write(*,*)
         write(*,*) 'Only horizon 2 exists.'
         write(*,*)

!        Check symmetries.

         if (CCTK_Equals(domain,"octant").eq.1) then

            write(*,*) 'Octant symmetry: calculating'
            write(*,*) 'twice the distance from horizon 2 to origin.'

            factor = two

            x1 = zero
            y1 = zero
            z1 = zero

         else if (CCTK_Equals(domain,"quadrant").eq.1) then

            if (CCTK_Equals(quadrant_direction,"x").eq.1) then

               write(*,*) 'Quadrant symmetry, x direction: calculating'
               write(*,*) 'twice the distance from horizon 2 to x axis.'

               factor = two

               x1 = x2
               y1 = zero
               z1 = zero

            else if (CCTK_Equals(quadrant_direction,"y").eq.1) then

               write(*,*) 'Quadrant symmetry, y direction: calculating'
               write(*,*) 'twice the distance from horizon 2 to y axis.'

               factor = two

               x1 = zero
               y1 = y2
               z1 = zero

            else if (CCTK_Equals(quadrant_direction,"z").eq.1) then

               write(*,*) 'Quadrant symmetry, z direction: calculating'
               write(*,*) 'twice the distance from horizon 2 to z axis.'

               factor = two

               x1 = zero
               y1 = zero
               z1 = z2

            end if

         else if (CCTK_Equals(domain,"bitant").eq.1) then

            if (CCTK_Equals(bitant_plane,"xy").eq.1) then

               write(*,*) 'Bitant symmetry, xy plane: calculating twice'
               write(*,*) 'the distance from horizon 2 to xy plane.'

               factor = two

               x1 = x2
               y1 = y2
               z1 = zero

            else if (CCTK_Equals(bitant_plane,"xz").eq.1) then

               write(*,*) 'Bitant symmetry, xz plane: calculating twice'
               write(*,*) 'the distance from horizon 2 to xz plane.'

               factor = two

               x1 = x2
               y1 = zero
               z1 = z2

            else if (CCTK_Equals(bitant_plane,"yz").eq.1) then

               write(*,*) 'Bitant symmetry, yz plane: calculating twice'
               write(*,*) 'the distance from horizon 2 to yz plane.'

               factor = two

               x1 = zero
               y1 = y2
               z1 = z2

            end if

         else

            write(*,*) 'No grid symmetries: no distance to calculate.'
            write(*,*) 'giving up.'

            d12 = zero
            goto 100

         end if

      end if

!     Find connecting line.

      x12 = x2 - x1
      y12 = y2 - y1
      z12 = z2 - z1

      rho12 = sqrt(x12**2 + y12**2)
      r12   = sqrt(x12**2 + y12**2 + z12**2)


!     ***********************************************************
!     ***   FIND LOCATION OF HORIZONS ALONG CONNECTING LINE   ***
!     ***********************************************************

!     Find angles for horizon 1.

      if (r12.ne.zero) then
         cost = z12/r12
      else
         write(*,*)
         write(*,*) 'Coordinate distance is zero,'
         write(*,*) 'giving up.'
         write(*,*)
         d12 = zero
         goto 100
      end if

      sint = sqrt(one - cost**2)

      if (y12.ne.zero) then
         phiangle = atan2(y12,x12)
      else
         phiangle = zero
      end if

      sinp = sin(phiangle)
      cosp = cos(phiangle)

!     Find point in horizon 1.

      if (status_old_1) then

         r1 = c0(0)

         do l=1+stepz,lmax,1+stepz
            r1 = r1 + c0(l)*LEGEN(l,0,cost)
         end do

         if (nonaxi) then
            do m=1,lmax
               aux = dble(m)*phiangle
               sina = sin(aux)
               cosa = cos(aux)
               do l=m,lmax
                  aux = LEGEN(l,m,cost)
                  r1 = r1 + aux*cc(l,m)*cosa
                  if (.not.refy) then
                     r1 = r1 + aux*cs(l,m)*sina
                  end if
               end do
           end do
         end if

         x1 = x1 + r1*sint*cosp
         y1 = y1 + r1*sint*sinp
         z1 = z1 + r1*cost

      end if

!     Find angles for horizon 2.

      cost = -cost
      sint = +sint

      cosp = -cosp
      sinp = -sinp

      phiangle = phiangle + pi

!     Find point in horizon 2.

      if (status_old_2) then

         r2 = c0(0)

         do l=1+stepz,lmax,1+stepz
            r2 = r2 + c0(l)*LEGEN(l,0,cost)
         end do

         if (nonaxi) then
            do m=1,lmax
               aux = dble(m)*phiangle
               sina = sin(aux)
               cosa = cos(aux)
               do l=m,lmax
                  aux = LEGEN(l,m,cost)
                  r2 = r2 + aux*cc(l,m)*cosa
                  if (.not.refy) then
                     r2 = r2 + aux*cs(l,m)*sina
                  end if
               end do
            end do
         end if

         x2 = x2 + r2*sint*cosp
         y2 = y2 + r2*sint*sinp
         z2 = z2 + r2*cost

      end if


!     **************************************
!     ***   ALLOCATE MEMORY FOR ARRAYS   ***
!     **************************************

      aux = sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)/dx

      np = int(aux)

      if (myproc.eq.0) then

         allocate(xa(0:np))
         allocate(ya(0:np))
         allocate(za(0:np))

         allocate(txx(0:np))
         allocate(tyy(0:np))
         allocate(tzz(0:np))
         allocate(txy(0:np))
         allocate(txz(0:np))
         allocate(tyz(0:np))

      else

         allocate(xa(1))
         allocate(ya(1))
         allocate(za(1))

         allocate(txx(1))
         allocate(tyy(1))
         allocate(tzz(1))
         allocate(txy(1))
         allocate(txz(1))
         allocate(tyz(1))

      end if


!     *************************************
!     ***   FIND INTERPOLATING POINTS   ***
!     *************************************

      if (myproc.eq.0) then

         npoints = np + 1

         do i=0,np

            aux = dble(i)/dble(np)

            xa(i) = x1 + aux*(x2-x1)
            ya(i) = y1 + aux*(y2-y1)
            za(i) = z1 + aux*(z2-z1)

         end do

      else

         npoints = 1

         xa = half*(x1 + x2)
         ya = half*(y1 + y2)
         za = half*(z1 + z2)

      end if


!     ******************************
!     ***   INTERPOLATE METRIC   ***
!     ******************************

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
      if (interp_handle.lt.0) then
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

      out_arrays(1) = CCTK_PointerTo(txx)
      out_arrays(2) = CCTK_PointerTo(tyy)
      out_arrays(3) = CCTK_PointerTo(tzz)
      out_arrays(4) = CCTK_PointerTo(txy)
      out_arrays(5) = CCTK_PointerTo(txz)
      out_arrays(6) = CCTK_PointerTo(tyz)

      out_array_type_codes = CCTK_VARIABLE_REAL


!     Interpolation.

      call CCTK_InterpGridArrays (ierror, cctkGH, 3, interp_handle, &
           param_table_handle, coord_system_handle, &
           npoints, CCTK_VARIABLE_REAL, interp_coords, &
           6, in_array_indices, &
           6, out_array_type_codes, out_arrays)
      if (ierror < 0) then
        call CCTK_WARN (1, "AHFinder_dis.F: interpolator call returned an error code");
      endif

!     release parameter table
      call Util_TableDestroy (ierror, param_table_handle)


!     ******************************
!     ***   INTEGRATE DISTANCE   ***
!     ******************************

      ddx = (x2-x1)/dble(np)
      ddy = (y2-y1)/dble(np)
      ddz = (z2-z1)/dble(np)

      if (myproc.eq.0) then

!        Integrate using trapezoidal rule.

         d12 = zero

         do i=1,np

            d12 = d12 + half &
                 *(sqrt(txx(i-1)*ddx**2  + tyy(i-1)*ddy**2  + tzz(i-1)*ddz**2 &
                 +  two*(txy(i-1)*ddx*ddy + txz(i-1)*ddx*ddz + tyz(i-1)*ddy*ddz)) &
                 + sqrt(txx(i  )*ddx**2  + tyy(i  )*ddy**2  + tzz(i  )*ddz**2 &
                 +  two*(txy(i  )*ddx*ddy + txz(i  )*ddx*ddz + tyz(i  )*ddy*ddz)))

         end do

      end if

      d12 = factor*d12


!     *****************************
!     ***   DEALLOCATE MEMORY   ***
!     *****************************

      deallocate(xa,ya,za)
      deallocate(txx,tyy,tzz,txy,txz,tyz)


!     ******************
!     ***   OUTPUT   ***
!     ******************

  100 continue

      disf = filestr(1:nfile) // "/ahf_d12" // ".tl"
      if (myproc.eq.0) then
        inquire(file=disf, exist=disf_exists)
        if (disf_exists) then
          open(11,file=disf,form='formatted',status='old', position='append')
        else
          open(11,file=disf,form='formatted',status='replace')
          write(11,"(A18)") '" Horizon distance'
        end if
        write(11,"(2ES14.6)") cctk_time,d12
        close(11)
      end if


!     ***************
!     ***   END   ***
!     ***************

      end subroutine AHFinder_dis
