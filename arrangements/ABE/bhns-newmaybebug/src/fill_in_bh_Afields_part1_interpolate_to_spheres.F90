#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
! Fill interior of AH with smoothly extrapolated A^i data.
!-----------------------------------------------------------------------------
subroutine fill_in_bh_Afields_part1_interpolate_to_spheres(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  ! Variables needed for interpolation
  CCTK_POINTER, dimension(3)               :: interp_coords
  character(60)                            :: options_string
  integer                                  :: interpolation_order,nchars
  integer                                  :: ierr,N_dims,interp_handle,param_table_handle,coord_system_handle,N_interp_points
  integer                                  :: N_input_arrays,N_output_arrays
  real*8, dimension(INPUTARRAY_THETASIZE*INPUTARRAY_PHISIZE)         :: xinterp,yinterp,zinterp
  real*8, dimension(INPUTARRAY_THETASIZE*INPUTARRAY_PHISIZE)         :: Axinterp,Ayinterp,Azinterp
  real*8, dimension(INPUTARRAY_THETASIZE*INPUTARRAY_PHISIZE)         :: rhobinterp,Pinterp
  integer,dimension(5)                    :: input_array_type_codes,input_array_varindices,output_array_type_codes
  CCTK_POINTER,dimension(5)               :: output_array_pointers
  ! Output arrays, dummy indices, parameters
  real*8                                   :: xcenter,ycenter,zcenter
  real*8                                   :: horizdirn_x,horizdirn_y,horizdirn_z
  real*8                                   :: theta,phiangle,costheta,sintheta,PI
  integer                                  :: i,j,n,whichradius,order

  if(MOD(cctk_iteration,refill_horizons_magfields_every)==0 .and. cctk_iteration.gt.0) then
     xcenter = bh_posn_x(1)
     ycenter = bh_posn_y(1)
     zcenter = bh_posn_z(1)

!     Get horizon radius in direction where it is likely to be minimized:
     horizdirn_x = 0.D0
     horizdirn_y = 0.D0
     horizdirn_z = 100000.D0
     call get_ah_radius_in_dirn(cctkGH,horizdirn_x,horizdirn_y,horizdirn_z,horiz_radius);

     write(*,*) "HIWWW.",xcenter,ycenter,zcenter,horiz_radius

     PI = 3.14159265358979323846D0

     N_dims = 3
     interpolation_order = 2
     N_interp_points = INPUTARRAY_THETASIZE*INPUTARRAY_PHISIZE

     interp_handle = -1 
     call CCTK_InterpHandle (interp_handle, "Lagrange polynomial interpolation")
     !  call CCTK_InterpHandle (interp_handle, "uniform cartesian")
     if (interp_handle .lt. 0) then
        call CCTK_WARN(0,"Cannot get handle for interpolation ! Forgot to activate an implementation providing interpolation operators ??")
     endif

     param_table_handle = -1
     options_string = "order = " // char(ichar('0') + interpolation_order)
     call Util_TableCreateFromString (param_table_handle, options_string)

     if (param_table_handle .lt. 0) then
        call CCTK_WARN(0,"Cannot create parameter table for interpolator")
     endif

     coord_system_handle = -1
     call CCTK_CoordSystemHandle (coord_system_handle, "cart3d")
     if (coord_system_handle .lt. 0) then
        call CCTK_WARN(0,"Cannot get handle for cart3d coordinate system ! Forgot to activate an implementation providing coordinates ??")
     endif

     input_array_type_codes = CCTK_VARIABLE_REAL
     output_array_type_codes = CCTK_VARIABLE_REAL

     ! Specify interpolation input arrays, output arrays:
     N_input_arrays = 5
     N_output_arrays = 5

     call CCTK_VarIndex (input_array_varindices(1), "mhd_evolve::Ax")
     call CCTK_VarIndex (input_array_varindices(2), "mhd_evolve::Ay")
     call CCTK_VarIndex (input_array_varindices(3), "mhd_evolve::Az")
     call CCTK_VarIndex (input_array_varindices(4), "mhd_evolve::rho_b")
     call CCTK_VarIndex (input_array_varindices(5), "mhd_evolve::P")

     do whichradius=NUM_ZERO_PTS+1,RADIAL_INTERP_ORDER
        ! Set up interpolation coordinate arrays:
        n = 1
        do i=1,INPUTARRAY_THETASIZE
           if(Symmetry.eq.1) then
              theta = (i-1.D0)/(INPUTARRAY_THETASIZE-1.D0)*PI*0.5D0
           else if(Symmetry.eq.0) then
              theta = (i-1.D0)/(INPUTARRAY_THETASIZE-1.D0)*PI
           else
              write(*,*) "Sorry, Symmetry=",Symmetry," not supported in fill_in_bh_Afields_part1_interpolate_to_spheres.F90"
              stop
           end if
           costheta = cos(theta)
           sintheta = sin(theta)
           do j=1,INPUTARRAY_PHISIZE
              phiangle = (j-1.D0)/(INPUTARRAY_PHISIZE-1.D0)*2.D0*PI-PI
              xinterp(n) = horiz_radius*bhns_bh_filling_radius(whichradius)*sintheta*cos(phiangle)+xcenter
              yinterp(n) = horiz_radius*bhns_bh_filling_radius(whichradius)*sintheta*sin(phiangle)+ycenter
              zinterp(n) = horiz_radius*bhns_bh_filling_radius(whichradius)*costheta+zcenter
              n = n + 1
           end do
        end do

        interp_coords(1) = CCTK_PointerTo(xinterp)
        interp_coords(2) = CCTK_PointerTo(yinterp)
        interp_coords(3) = CCTK_PointerTo(zinterp)

        output_array_pointers(1) = CCTK_PointerTo(Axinterp)
        output_array_pointers(2) = CCTK_PointerTo(Ayinterp)
        output_array_pointers(3) = CCTK_PointerTo(Azinterp)
        output_array_pointers(4) = CCTK_PointerTo(rhobinterp)
        output_array_pointers(5) = CCTK_PointerTo(Pinterp)

        ! Perform interpolation:
        call CCTK_InterpGridArrays(ierr,cctkGH,N_dims,interp_handle, &
             param_table_handle,coord_system_handle, &
             N_interp_points,input_array_type_codes,interp_coords, &
             N_input_arrays,input_array_varindices, &
             N_output_arrays, output_array_type_codes, output_array_pointers)

        n = 1
        do i=1,INPUTARRAY_THETASIZE
           do j=1,INPUTARRAY_PHISIZE
              Axintr(whichradius-NUM_ZERO_PTS,i,j) = Axinterp(n)
              Ayintr(whichradius-NUM_ZERO_PTS,i,j) = Ayinterp(n)
              Azintr(whichradius-NUM_ZERO_PTS,i,j) = Azinterp(n)
              rhobintr(whichradius-NUM_ZERO_PTS,i,j) = rhobinterp(n)
              Pintr(whichradius-NUM_ZERO_PTS,i,j) = rhobinterp(n)
              n = n + 1
           end do
        end do
     end do
  end if
end subroutine fill_in_bh_Afields_part1_interpolate_to_spheres
