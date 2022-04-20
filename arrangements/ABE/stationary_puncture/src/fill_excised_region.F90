#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine fill_excised_region(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer                                  :: i,j,k
  integer, dimension(3)                    :: ext,global_ext
  real*8                                   :: dT,dX,dY,dZ
  real*8, dimension(1,3)                   :: pointcoords
  real*8                                   :: phib,lapm1b,rb,PI,actual_excis_radius,excis_radius_fish
  real*8                                   :: shiftxb,shiftyb,shiftzb,shiftrb
  real*8                                   :: gxxb,gxyb,gxzb,gyyb,gyzb,gzzb,Gammaxb,Gammayb,Gammazb
  real*8                                   :: Axxb,Axyb,Axzb,Ayyb,Ayzb,Azzb
  real*8                                   :: gupxxb,gupxyb,gupxzb,gupyyb,gupyzb,gupzzb


  integer :: index,dummy,foundit_flag
  integer :: handle,ierr
  real*8  :: ONE,ZERO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(ONE = 1.D0, ZERO = 0.D0)

  ext = cctk_lsh

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  if(fill_excision_enable==2) then
     do k=1,ext(3)
        do j=1,ext(2)
           do i=1,ext(1)
              if(PhysicalRadius(i,j,k) .lt. excis_radius) then
                 psi(i,j,k) = 2.875D0 - 5.D0 * PhysicalRadius(i,j,k)**2  + 6.D0 * PhysicalRadius(i,j,k)**4 !2.875-5(r/M)**2+6(r/M)**4
                 phi(i,j,k) = log(psi(i,j,k))
                 
                 lapm1(i,j,k) = psi(i,j,k)**(-2) - 1.D0
              end if
           end do
        end do
     end do

  else if(fill_excision_enable==1 .and. asymptotic_bh==0) then
     do k=1,ext(3)
        do j=1,ext(2)
           do i=1,ext(1)
              if(PhysicalRadius(i,j,k) .lt. excis_radius) then
                 phi(i,j,k) = log(1.D0 + 1.D0/(2.D0*excis_radius))
                 psi(i,j,k) = exp(phi(i,j,k))
                 
                 lapm1(i,j,k) = psi(i,j,k)**(-2) - 1.D0
                 
              end if
           end do
        end do
     end do
     

  else if(fill_excision_enable==1 .and. asymptotic_bh==1) then
     !Fill in missing data from excised initial data.

     !First setup Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,and Zglobal, although they should already be set...
     call setup_global_coord_arrays(CCTK_PASS_FTOF)

     pointcoords(1,1) = excis_radius
     global_ext = cctk_gsh

     !set the y and z-coordinates
     if(Symmetry.eq.EQUATORIAL) then
        pointcoords(1,2) = dY*0.5D0
        pointcoords(1,3) = dZ*0.5D0
     else if(Symmetry.eq.AXISYM) then
        pointcoords(1,2) = 0.D0
        pointcoords(1,3) = dZ*0.5D0
     else if(Symmetry.eq.OCTANT) then
        pointcoords(1,2) = dY*0.5D0
        pointcoords(1,3) = dZ*0.5D0
     else if(Symmetry.eq.NO_SYMM) then
        pointcoords(1,2) = dY*0.5D0
        pointcoords(1,3) = dZ*0.5D0
     end if

     !Next find excis_radius_fish (= excis_radius in grid-native fisheye coordinates).  This works even if fisheye is turned off.
     foundit_flag = 0
     excis_radius_fish = 1.D307
     do i=1,global_ext(1)
        if(Xglobal(i).gt.0.D0) then
           pointcoords(1,1) = Xglobal(i)
           call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,PhysicalRadius,rb)
           if(rb.gt.excis_radius .and. foundit_flag.eq.0) then
              !           write(*,*) "FOUND THE EXCISION RADIUS (in grid-native coordinates)!"
              !           write(*,*) "excis_radius_fish =",Xglobal(i),rb,excis_radius,i
              !write(*,*) "PhysicalRadius",PhysicalRadius(:,2,2)
              excis_radius_fish = Xglobal(i)
              foundit_flag = 1
           end if
        end if
     end do
     call CCTK_ReductionHandle(handle,"minimum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,excis_radius_fish,excis_radius_fish,CCTK_VARIABLE_REAL)

     !Make sure the specified excision radius is in the grid!
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,foundit_flag,foundit_flag,CCTK_VARIABLE_INT)
     if(foundit_flag.eq.0) then
        write(*,*) "Could not find the excision radius specified (excis_radius) on the grid! :("
        stop
     end if

     pointcoords(1,1) = excis_radius_fish

     !Now fill excision region with zeros, foo (whom I pity)
     do i=1,ext(1)
        do j=1,ext(2)
           do k=1,ext(3)
              if(r(i,j,k) .lt. excis_radius_fish) then
                 psi(i,j,k) = 0.D0
                 phi(i,j,k) = 0.D0
                 lapm1(i,j,k) = 0.D0
                 shiftx(i,j,k) = 0.D0
                 shifty(i,j,k) = 0.D0
                 shiftz(i,j,k) = 0.D0
                 gxx(i,j,k) = 0.D0
                 gxy(i,j,k) = 0.D0
                 gxz(i,j,k) = 0.D0
                 gyy(i,j,k) = 0.D0
                 gyz(i,j,k) = 0.D0
                 gzz(i,j,k) = 0.D0
                 gupxx(i,j,k) = 0.D0
                 gupxy(i,j,k) = 0.D0
                 gupxz(i,j,k) = 0.D0
                 gupyy(i,j,k) = 0.D0
                 gupyz(i,j,k) = 0.D0
                 gupzz(i,j,k) = 0.D0
                 Axx(i,j,k) = 0.D0
                 Axy(i,j,k) = 0.D0
                 Axz(i,j,k) = 0.D0
                 Ayy(i,j,k) = 0.D0
                 Ayz(i,j,k) = 0.D0
                 Azz(i,j,k) = 0.D0
                 Gammax(i,j,k) = 0.D0
                 Gammay(i,j,k) = 0.D0
                 Gammaz(i,j,k) = 0.D0
              end if
           end do
        end do
     end do

     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,phi,phib)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,lapm1,lapm1b)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,shiftx,shiftxb)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,shifty,shiftyb)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,shiftz,shiftzb)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gxx,gxxb)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gxy,gxyb)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gxz,gxzb)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gyy,gyyb)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gyz,gyzb)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gzz,gzzb)

     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gupxx,gupxxb)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gupxy,gupxyb)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gupxz,gupxzb)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gupyy,gupyyb)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gupyz,gupyzb)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,gupzz,gupzzb)

     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,Axx,Axxb)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,Axy,Axyb)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,Axz,Axzb)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,Ayy,Ayyb)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,Ayz,Ayzb)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,Azz,Azzb)

     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,Gammax,Gammaxb)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,Gammay,Gammayb)
     call statpunc_return_funcvals_at_points(cctkGH,cctk_nghostzones,1,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,Gammaz,Gammazb)

     !IMPORTANT: in all cases, actual_excis_radius >= excis_radius, but only by a small amount ( <= dX)
     actual_excis_radius = sqrt(pointcoords(1,1)**2 + pointcoords(1,2)**2 + pointcoords(1,3)**2)

     do i=1,ext(1)
        do j=1,ext(2)
           do k=1,ext(3)
              if(r(i,j,k) .lt. actual_excis_radius) then
                 phi(i,j,k) = phib
                 psi(i,j,k) = exp(phib)

                 lapm1(i,j,k) = lapm1b

                 shiftrb = sqrt(shiftxb*shiftxb + shiftyb*shiftyb + shiftzb*shiftzb)

                 shiftx(i,j,k) = (X(i,j,k)/actual_excis_radius) * shiftrb
                 shifty(i,j,k) = (Y(i,j,k)/actual_excis_radius) * shiftrb
                 shiftz(i,j,k) = (Z(i,j,k)/actual_excis_radius) * shiftrb

                 gxx(i,j,k) = gxxb
                 gxy(i,j,k) = gxyb
                 gxz(i,j,k) = gxzb
                 gyy(i,j,k) = gyyb
                 gyz(i,j,k) = gyzb
                 gzz(i,j,k) = gzzb

                 gupxx(i,j,k) = gupxxb
                 gupxy(i,j,k) = gupxyb
                 gupxz(i,j,k) = gupxzb
                 gupyy(i,j,k) = gupyyb
                 gupyz(i,j,k) = gupyzb
                 gupzz(i,j,k) = gupzzb

                 Axx(i,j,k) = Axxb
                 Axy(i,j,k) = Axyb
                 Axz(i,j,k) = Axzb
                 Ayy(i,j,k) = Ayyb
                 Ayz(i,j,k) = Ayzb
                 Azz(i,j,k) = Azzb

                 Gammax(i,j,k) = Gammaxb
                 Gammay(i,j,k) = Gammayb
                 Gammaz(i,j,k) = Gammazb
              end if
           end do
        end do
     end do
  end if
end subroutine fill_excised_region

!-------------------------------------------------------------------------------
!
! Simplified version of interpolate_pointset: simply returns set of function 
!   values on stencil.
!
!------------------------------------------------------------------------------
subroutine statpunc_return_funcvals_at_points(cctkGH,nghostzones,numpoints,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,func,outputvals)
  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_POINTER :: cctkgh
  integer, dimension(3)                    :: ext,global_ext,nghostzones
  integer                                  :: numpoints
  real*8, dimension(3)                     :: offset
  real*8, dimension(numpoints,3)           :: pointcoords
  real*8, dimension(3)                     :: localpointcoord
  real*8,dimension(ext(1),ext(2),ext(3))   :: func
  real*8,dimension(ext(1))                 :: Xlocal1d
  real*8,dimension(ext(2))                 :: Ylocal1d
  real*8,dimension(ext(3))                 :: Zlocal1d
  real*8,dimension(global_ext(1))          :: Xglobal
  real*8,dimension(global_ext(2))          :: Yglobal
  real*8,dimension(global_ext(3))          :: Zglobal
  real*8                                   :: dX,dY,dZ,teeny,oneodX,oneodY,oneodZ
  real*8                                   :: interpvaluej,interpvaluek
  real*8,dimension(numpoints)              :: outputvals,funcvalues
  integer                                  :: n,ii,jj,kk
  integer, dimension(1)                    :: ivalue,jvalue,kvalue
  integer, dimension(1)                    :: ivalueglobal,jvalueglobal,kvalueglobal
  integer                                  :: imin,jmin,kmin,imax,jmax,kmax
  integer                                  :: handle,index,ierr
  integer                                  :: nsymghostzonesx,nsymghostzonesy,nsymghostzonesz
  integer                                  :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  imin = 0
  jmin = 0
  kmin = 0
  imax = ext(1)
  jmax = ext(2)
  kmax = ext(3)

  ivalue = 1
  jvalue = 1
  kvalue = 1
  ivalueglobal = 1
  jvalueglobal = 1
  kvalueglobal = 1

  teeny = dX*1.D-5
  oneodx = 1.D0/dX
  oneody = 1.D0/dY
  oneodz = 1.D0/dZ

  nsymghostzonesx = 0
  nsymghostzonesy = 0
  nsymghostzonesz = 0

  if(Symmetry.eq.EQUATORIAL) then
!     nsymghostzonesz=1
     nsymghostzonesz=nghostzones(3)
  else if(Symmetry.eq.AXISYM) then
     nsymghostzonesx=nghostzones(1)-1
     nsymghostzonesz=nghostzones(3)
  else if(Symmetry.eq.OCTANT) then
     nsymghostzonesx=nghostzones(1)
     nsymghostzonesy=nghostzones(2)
     nsymghostzonesz=nghostzones(3)
  end if

  !By default, assume upper bound is processor boundary and not outer boundary.
  !We subtract # of ghostzones*2 from ext(i) to avoid overlap in our stencil points
  imax = ext(1)-nghostzones(1)*2
  jmax = ext(2)-nghostzones(2)*2
  kmax = ext(3)-nghostzones(3)*2

  !check if upper proc bound is global boundary and fix ijkmax
  if(abs(Xlocal1d(ext(1))-Xglobal(global_ext(1))).lt.teeny) imax = ext(1)
  if(abs(Ylocal1d(ext(2))-Yglobal(global_ext(2))).lt.teeny) jmax = ext(2)
  if(abs(Zlocal1d(ext(3))-Zglobal(global_ext(3))).lt.teeny) kmax = ext(3)

  do n=1,numpoints
     !We define ijkvalueglobal(1) = corner of stencil closest to (-inf,-inf,-inf)
     ivalueglobal(1) = (pointcoords(n,1)-Xglobal(1))*oneodx+1
     jvalueglobal(1) = (pointcoords(n,2)-Yglobal(1))*oneody+1
     kvalueglobal(1) = (pointcoords(n,3)-Zglobal(1))*oneodz+1

     !If point touches outside of global domain, move it back into global domain
     if(ivalueglobal(1).gt.global_ext(1)) then
        ivalueglobal(1) = global_ext(1)
        write(*,*) "OUTSIDE X BOUND"
     end if
     if(jvalueglobal(1).gt.global_ext(2)) then
        jvalueglobal(1) = global_ext(2)
        write(*,*) "OUTSIDE Y BOUND"
     end if
     if(kvalueglobal(1).gt.global_ext(3)) then
        kvalueglobal(1) = global_ext(3)
        write(*,*) "OUTSIDE Z BOUND"
     end if

     !Now that we know global point, hunt for that point on local processor!
     localpointcoord(1) = Xglobal(ivalueglobal(1))+teeny
     localpointcoord(2) = Yglobal(jvalueglobal(1))+teeny
     localpointcoord(3) = Zglobal(kvalueglobal(1))+teeny

     !Next search for coordinate on local domain!
     ivalue(1) = (localpointcoord(1)-Xlocal1d(1))*oneodx+1
     jvalue(1) = (localpointcoord(2)-Ylocal1d(1))*oneody+1
     kvalue(1) = (localpointcoord(3)-Zlocal1d(1))*oneodz+1

     !If the glove doesn't fit, you must acquit!
     !I.e., when search returns <0 or >ijkmax, that means the point is off the local processor grid!
     if(ivalue(1).lt.1 .or. ivalue(1).gt.imax) ivalue(1) = -1
     if(jvalue(1).lt.1 .or. jvalue(1).gt.jmax) jvalue(1) = -1
     if(kvalue(1).lt.1 .or. kvalue(1).gt.kmax) kvalue(1) = -1

     !     write(*,*) "hix.",ivalue(1),n,pointcoords(n,1),Xglobal(ivalueglobal(1))
     !     write(*,*) "hiy.",jvalue(1),n,pointcoords(n,2),Yglobal(jvalueglobal(1))
     !     write(*,*) "hiz.",kvalue(1),n,pointcoords(n,3),Zglobal(kvalueglobal(1))

     !following shouldn't be a problem since 1.D307 is very close to infinity for double precision
     funcvalues(n) = 1.D307
     if(ivalue(1).gt.0 .and. jvalue(1).gt.0 .and. kvalue(1).gt.0) then
        ii = ivalue(1)
        jj = jvalue(1)
        kk = kvalue(1)
        funcvalues(n) = func(ii,jj,kk)
     end if
  end do
  outputvals = 0.D0
  call CCTK_ReductionHandle(handle,"minimum")
  call CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,&
       funcvalues,outputvals,numpoints,CCTK_VARIABLE_REAL)
end subroutine statpunc_return_funcvals_at_points
