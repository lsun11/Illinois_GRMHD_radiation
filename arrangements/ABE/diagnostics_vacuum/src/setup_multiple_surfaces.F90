!-----------------------------------------------------------------------------
!
!$Id: setup_multiple_surfaces.F90  $
!
!-----------------------------------------------------------------------------
!
! setup surface integrals 
!
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
subroutine setup_multiple_surfaces(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,global_ext
  real*8                                   :: dX,dY,dZ
  real*8, dimension(1,3)                   :: pointcoords
  real*8                                   :: PI,costheta,sintheta,phiangle,zmax
  integer                                  :: handle,ierr,vindex
  integer                                  :: i,n,j,int_order
  integer                                  :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  real*8 :: maxsurf,dsurf,rr,rbr2,drbdr2,ddrbddr2

  if(cctk_iteration==0) then

     write(*,*) "HI... STARTING setup_multiple_surfaces"

     dX = CCTK_DELTA_SPACE(1)
     dY = CCTK_DELTA_SPACE(2)
     dZ = CCTK_DELTA_SPACE(3)

     ext = cctk_lsh
     global_ext = cctk_gsh

     PI = acos(-1.D0)

     int_order=2
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex(vindex,"grid::Z")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,zmax,1,vindex)

     maxsurf=0.5d0*dZ+zmax
     if(scaledsurf==1) then
        rsurf2=rsurf2*maxsurf
     endif

     write(6,*)'nsurf:',nsurf,' maxsurf:',maxsurf,' scaledsurf:',scaledsurf,surf_radius,rsurf2
     if(nsurf.gt.20) then
        write(6,*)'nsurf cannot be bigger than 20, since we cannot dynamically'
        write(6,*)'allocate the vectors in surf_params.  If you need more space'
        write(6,*)'edit moving_puncture/interface.ccl'
        stop
     else if(nsurf.eq.0) then
        return
     else if (nsurf.ge.1) then
        surfvec(1)=surf_radius
        rbrvec(1)=rbr
        drbrvec(1)=drbdr
        ddrbrvec(1)=ddrbddr
        if(nsurf.eq.1) then
           write(*,*) "HI... finished with setup_multiple_surfaces"
           return
        end if
        write(6,*)'surf 1:r,rb,drb,ddrb:',surfvec(1),rbrvec(1),drbrvec(1),ddrbrvec(1)

     endif

     if(arithsurf==1) then
        dsurf=(rsurf2-surf_radius)/(nsurf-1.0d0)
     else
        dsurf=exp((log(rsurf2)-log(surf_radius))/(nsurf-1.0d0))
     endif

     do i=2,nsurf

        if(arithsurf==1) then
           rr=surf_radius+(i-1)*dsurf
        else
           rr=surf_radius*dsurf**(i-1)
        endif

        pointcoords(1,1)=rr
        pointcoords(1,2)=0.
        pointcoords(1,3)=0.

        if(fisheye_enable.eq.1) then
           call CCTK_VarIndex(vindex,"fisheye::PhysicalRadius")
           call interp_driver_carp(cctkGH,1,pointcoords,vindex,rbr2)
           call CCTK_VarIndex(vindex,"fisheye::RadiusDerivative")
           call interp_driver_carp(cctkGH,1,pointcoords,vindex,drbdr2)
           call CCTK_VarIndex(vindex,"fisheye::RadiusDerivative2")
           call interp_driver_carp(cctkGH,1,pointcoords,vindex,ddrbddr2)
           !        call interpolate_pointset(cctkGH,cctk_nghostzones,1,int_order,ext,global_ext,dX,dY,dZ, &
           !             Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,PhysicalRadius,rbr2) 
           !        call interpolate_pointset(cctkGH,cctk_nghostzones,1,int_order,ext,global_ext,dX,dY,dZ, &
           !             Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,RadiusDerivative,drbdr2) 
           !        call interpolate_pointset(cctkGH,cctk_nghostzones,1,int_order,ext,global_ext,dX,dY,dZ, &
           !             Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,RadiusDerivative2,ddrbddr2) 
           rbr2=rbr2/rr
           !     drbdr=1.0/drbdr
        else
           rbr2=1.0
           drbdr2=1.0
           ddrbddr2=0.0
        endif

        surfvec(i)=rr
        rbrvec(i)=rbr2
        drbrvec(i)=drbdr2
        ddrbrvec(i)=ddrbddr2
        write(6,*)'surf ',i,':r,rb,drb,ddrb:',surfvec(i),rbrvec(i),drbrvec(i),ddrbrvec(i)

     enddo
  end if
end subroutine setup_multiple_surfaces
