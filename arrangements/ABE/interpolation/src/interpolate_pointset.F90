!
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-------------------------------------------------------------------------------
!
! Interpolate a set of points to 1st, 2nd, or 3rd order, 
!       over multiple procs if necessary.  
! Performs Lagrangian polynomial interpolation.
!
! By Zach Etienne, based on Cactus' LocalInterp interpolator
!   Special thanks to Josh Faber for useful suggestions.
!
! Unlike the interpolator built-in to PUGHInterp/LocalInterp, 
!  this is a global interpolator (i.e., will interpolate to _any_
!  point in the grid boundaries, regardless of which processor it is on)
!
! TODO, in decreasing priority 0->4:
! 0) Further debugging
! 1) Make warning statements appear on proc 0.
! 2) Better center 2nd order interpolation?  Currently goes x0 _y___ x1 _____ x2
!      where y is the interpolation point and x0<x1<x2.  y can never appear 
!      between x1 and x2 currently, unless x2 is the outer boundary point and 
!      x1 < y < x2
! 3) (0,0,0) in axisymmetry: force stencil centering, applying symmetry BCs?
! 4) Make it faster still.
!
!------------------------------------------------------------------------------
subroutine interpolate_pointset(cctkGH,nghostzones,numpoints,order,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal,pointcoords,func,outputvals)
  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_POINTER :: cctkgh
  integer, dimension(3)                    :: ext,nghostzones,global_ext
  integer                                  :: numpoints,order
  real*8, dimension(3,order+1)             :: coeff
  real*8, dimension(3)                     :: offset
  real*8, dimension(numpoints,3)           :: pointcoords
  real*8, dimension(3)                     :: stencpointcoord
  real*8,dimension(ext(1),ext(2),ext(3))   :: func
  real*8,dimension(ext(1))                 :: Xlocal1d
  real*8,dimension(ext(2))                 :: Ylocal1d
  real*8,dimension(ext(3))                 :: Zlocal1d
  real*8,dimension(global_ext(1))          :: Xglobal
  real*8,dimension(global_ext(2))          :: Yglobal
  real*8,dimension(global_ext(3))          :: Zglobal
  real*8                                   :: dX,dY,dZ,teeny,oneodX,oneodY,oneodZ
  real*8                                   :: interpvaluej,interpvaluek
  real*8,dimension(numpoints)              :: interpvalues,outputvals
  integer                                  :: n,i,j,k,ll,mm,nn
  integer, dimension(order+1)              :: ivalue,jvalue,kvalue
  integer, dimension(order+1)              :: ivalueglobal,jvalueglobal,kvalueglobal
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

  !Initialize interpvalues to zero
  interpvalues = 0.D0

  !The following ensures that stencil does not invade symmetry ghost zone region of grid.
  !Thus near a symmetry axis, extrapolation is performed instead of interpolation.
  ! This is for DAGH compatibility purposes only.
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

  if(Symmetry.ne.AXISYM) then
     !By default, assume upper bound is processor boundary and not outer boundary.
     !We subtract # of ghostzones*2 from ext(i) to avoid overlap in our stencil points
!     imax = ext(1)-2
!     jmax = ext(2)-2
!     kmax = ext(3)-2

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

        !Below if statement ensures that 3rd order stencil is centered.
        if(order==3) then
           ivalueglobal(1) = ivalueglobal(1) - 1
           jvalueglobal(1) = jvalueglobal(1) - 1
           kvalueglobal(1) = kvalueglobal(1) - 1
        end if

        !If stencil touches outside of global domain, move it back into global domain
        if(ivalueglobal(1)+(order).gt.global_ext(1)) then
           write(*,*) "WARNING OUTSIDE GLOBAL GRID in X direction!",ivalueglobal(1),ivalueglobal(1)+(order),global_ext(1),oneodx,pointcoords(n,1),Xglobal(1)
           ivalueglobal(1) = global_ext(1)-(order)
        end if
        if(jvalueglobal(1)+(order).gt.global_ext(2)) then
           write(*,*) "WARNING OUTSIDE GLOBAL GRID in Y direction!",ivalueglobal(2),ivalueglobal(2)+(order),global_ext(2),oneody,pointcoords(n,2),Yglobal(1)
           jvalueglobal(1) = global_ext(2)-(order)
        end if
        if(kvalueglobal(1)+(order).gt.global_ext(3)) then
           write(*,*) "WARNING OUTSIDE GLOBAL GRID in Z direction!",ivalueglobal(3),ivalueglobal(3)+(order),global_ext(3),oneodz,pointcoords(n,3),Zglobal(1)
           kvalueglobal(1) = global_ext(3)-(order)
        end if

        !setup global stencil array indices
        if(ivalueglobal(1).lt.1+nsymghostzonesx) ivalueglobal(1)=1+nsymghostzonesx
        if(jvalueglobal(1).lt.1+nsymghostzonesy) jvalueglobal(1)=1+nsymghostzonesy
        if(kvalueglobal(1).lt.1+nsymghostzonesz) kvalueglobal(1)=1+nsymghostzonesz
        do i=2,order+1
           ivalueglobal(i) = ivalueglobal(i-1)+1
           jvalueglobal(i) = jvalueglobal(i-1)+1
           kvalueglobal(i) = kvalueglobal(i-1)+1
        end do

        !Now that we know global stencil, hunt for all the stencil points on local processor!
        do i=1,order+1
           stencpointcoord(1) = Xglobal(ivalueglobal(i))+teeny
           stencpointcoord(2) = Yglobal(jvalueglobal(i))+teeny
           stencpointcoord(3) = Zglobal(kvalueglobal(i))+teeny
           
           !Next search for coordinate on local domain!
           ivalue(i) = (stencpointcoord(1)-Xlocal1d(1))*oneodx+1
           jvalue(i) = (stencpointcoord(2)-Ylocal1d(1))*oneody+1
           kvalue(i) = (stencpointcoord(3)-Zlocal1d(1))*oneodz+1

           !If the glove doesn't fit, you must acquit!
           !I.e., when search returns <0 or >ijkmax, that means the point is off the local processor grid!
           if(ivalue(i).lt.1 .or. ivalue(i).gt.imax) ivalue(i) = -1
           if(jvalue(i).lt.1 .or. jvalue(i).gt.jmax) jvalue(i) = -1
           if(kvalue(i).lt.1 .or. kvalue(i).gt.kmax) kvalue(i) = -1
        end do

!        if(n.eq.3121) then
!           write(*,*) "hi0",ivalueglobal,jvalueglobal,kvalueglobal
!           write(*,*) "hi1",ivalue,jvalue,kvalue
!        end if

        ! offset from grid point given by ijkvalueglobal(1), in fractions of grid points 
        !In language of coefficient equations below, offset(1) = (x-x0)/dX
        offset(1) = (pointcoords(n,1)-Xglobal(ivalueglobal(1)))*oneodx
        offset(2) = (pointcoords(n,2)-Yglobal(jvalueglobal(1)))*oneody
        offset(3) = (pointcoords(n,3)-Zglobal(kvalueglobal(1)))*oneodz

!        if(abs(offset(1)).gt.order) write(*,*) "bad pointx:",pointcoords(n,1),Xglobal(ivalueglobal(1))
!        if(abs(offset(2)).gt.order) write(*,*) "bad pointy:",pointcoords(n,2),Yglobal(jvalueglobal(1))
!        if(abs(offset(3)).gt.order) write(*,*) "bad pointz:",pointcoords(n,3),Zglobal(kvalueglobal(1))
!        if(n.eq.3121) then
!           write(*,*) "hi2",offset
!        end if

        ! *** compute the interpolation coefficients according to the order ***
        !
        ! (Thanks to Erik Schnetter/Numerical Recipes for formulating these so nicely.) 
        !
        ! These formulas are "just" the coefficients of the classical
        ! Lagrange interpolation polynomials along each dimension.
        ! For example, in 1 dimension the unique quadratic passing
        ! through the 3 points {(x0,y0), (x1,y1), (x2,y2)} is: 
        ! ( x-x1)( x-x2)        ( x-x0)( x-x2)        ( x-x0)( x-x1) 
        ! -------------- y0  +  -------------- y1  +  -------------- y2 
        ! (x0-x1)(x0-x2)        (x1-x0)(x1-x2)        (x2-x0)(x2-x1) 
        ! (It's easy to see this: each of the terms is yi if x=xi, or
        ! zero if x=any other xj.)  To get the formulas below, just negate 
        ! each (x-x) factor, and substitute the values xi=i.
        if(order==1) then
           do i=1,3
              !            (x-x1)/(x0-x1) = (offset(i)-1)/(-1)
              coeff(i,1) = 1 - offset(i)
              !            (x-x0)/(x1-x0) = offset(i)
              coeff(i,2) =     offset(i)
           end do
        else if(order==2) then
           do i=1,3
              ! It goes: x0 _y__ x1 ____ x2, where y is the point were interpolating to.
              !recall that x-x0 = offset
              !y0 coefficient: ((x-x1)    *  (x-x2))     / ((x2-x0)*(x1-x0))
              coeff(i,1) = (offset(i)-1) * (offset(i)-2) / (  2  *   1 )
              !y1 coefficient: ((x-x0)   *  (x-x2))      / ((x1-x0)*(x1-x2))
              coeff(i,2) =   (offset(i)) * (offset(i)-2) / (  1  * (-1))
              !y2 coefficient: ((x-x0)  *   ( x-x1))    / ((x2-x1)(x2-x0))
              coeff(i,3) =  (offset(i)) * (offset(i)-1) / ((1) * (2))
           end do
        else if(order==3) then
           do i=1,3
              !Below is basically the same scheme as above, but multiplied by -1 in num & denom.  
              !All of this could be optimized.
              coeff(i,1) = (1-offset(i)) * (2-offset(i)) * (3-offset(i)) / (  3  *   2  *   1 )
              coeff(i,2) = ( -offset(i)) * (2-offset(i)) * (3-offset(i)) / (  2  *   1  * (-1))
              coeff(i,3) = ( -offset(i)) * (1-offset(i)) * (3-offset(i)) / (  1  * (-1) * (-2))
              coeff(i,4) = ( -offset(i)) * (1-offset(i)) * (2-offset(i)) / ((-1) * (-2) * (-3))
           end do
        end if

        !if(ll.gt. ... statements ensure that points on interp stencil are not double counted!
        do k=1,order+1
           nn=kvalue(k)
           if(nn.gt.0) then
              interpvaluek = 0.D0
              do j=1,order+1
                 mm=jvalue(j)
                 if(mm.gt.0) then
                    interpvaluej = 0.D0
                    do i=1,order+1
                       ll=ivalue(i)
                       if(ll.gt.0) then
                          interpvaluej = interpvaluej + func(ll,mm,nn) * coeff(1,i)
!                          write(*,*) "hi.",ll,mm,nn,func(ll,mm,nn),interpvaluej
                       end if
                    end do
                    interpvaluek = interpvaluek + interpvaluej * coeff(2,j)
                 end if
              end do
              interpvalues(n) = interpvalues(n) + interpvaluek * coeff(3,k)
           end if
        end do
     end do
  else if(Symmetry.eq.AXISYM) then
     !By default, assume upper bound is processor boundary and not outer boundary.
     !We subtract # of ghostzones+1 from ext(i) to avoid overlap in our stencil points
     imax = ext(1)-4
     kmax = ext(3)-2
     !check if upper proc bound is global boundary and fix ijkmax
     if(abs(Xlocal1d(ext(1))-Xglobal(global_ext(1))).lt.teeny) imax = ext(1)
     if(abs(Zlocal1d(ext(3))-Zglobal(global_ext(3))).lt.teeny) kmax = ext(3)
     
     do n=1,numpoints
        !We define ijkvalueglobal(1) = corner of stencil closest to (-inf,-inf,-inf)
        ivalueglobal(1) = (pointcoords(n,1)-Xglobal(1))*oneodx+1
        kvalueglobal(1) = (pointcoords(n,3)-Zglobal(1))*oneodz+1

        !Below if statement ensures that 3rd order stencil is centered.
        if(order==3) then
           ivalueglobal(1) = ivalueglobal(1) - 1
           kvalueglobal(1) = kvalueglobal(1) - 1
        end if

        !If stencil touches outside of global domain, move it back into global domain
        if(ivalueglobal(1)+(order).gt.global_ext(1)) then
           ivalueglobal(1) = global_ext(1)-(order)
           write(*,*) "WARNING OUTSIDE GLOBAL GRID in X direction!",global_ext(1),ivalueglobal(1)+(order)
        end if
        if(kvalueglobal(1)+(order).gt.global_ext(3)) then
           kvalueglobal(1) = global_ext(3)-(order)
           write(*,*) "WARNING OUTSIDE GLOBAL GRID in Z direction!",global_ext(3),kvalueglobal(1)+(order)
        end if

        !setup global stencil array indices
        if(ivalueglobal(1).lt.1+nsymghostzonesx) ivalueglobal(1)=1+nsymghostzonesx
        if(kvalueglobal(1).lt.1+nsymghostzonesz) kvalueglobal(1)=1+nsymghostzonesz
        do i=2,order+1
           ivalueglobal(i) = ivalueglobal(i-1)+1
           kvalueglobal(i) = kvalueglobal(i-1)+1
        end do

        !Now that we know global stencil, hunt for all the stencil points on local processor!
        do i=1,order+1
           stencpointcoord(1) = Xglobal(ivalueglobal(i))+teeny
           stencpointcoord(3) = Zglobal(kvalueglobal(i))+teeny

           !Next search for coordinate on local domain!
           ivalue(i) = (stencpointcoord(1)-Xlocal1d(1))*oneodx+1
           kvalue(i) = (stencpointcoord(3)-Zlocal1d(1))*oneodz+1

           !If the glove doesn't fit, you must acquit!
           !I.e., when search returns <0 or >ijkmax, that means the point is off the local processor grid!
           if(ivalue(i).lt.1 .or. ivalue(i).gt.imax) ivalue(i) = -1
           if(kvalue(i).lt.1 .or. kvalue(i).gt.kmax) kvalue(i) = -1
        end do

        !write(*,*) "hi0",ivalueglobal,kvalueglobal
        !write(*,*) "hi1",ivalue,kvalue

        ! offset from grid point given by ijkvalueglobal(1), in fractions of grid points 
        !In language of coefficient equations below, offset(1) = (x-x0)/dX
        offset(1) = (pointcoords(n,1)-Xglobal(ivalueglobal(1)))*oneodx
        offset(2) = 0.D0
        offset(3) = (pointcoords(n,3)-Zglobal(kvalueglobal(1)))*oneodz

!        if(abs(offset(1)).gt.order) write(*,*) "bad pointx:",pointcoords(n,1),Xglobal(ivalueglobal(1))
!        if(abs(offset(2)).gt.order) write(*,*) "bad pointy:",pointcoords(n,2),Yglobal(jvalueglobal(1))
!        if(abs(offset(3)).gt.order) write(*,*) "bad pointz:",pointcoords(n,3),Zglobal(kvalueglobal(1))

!        write(*,*) "hi2",offset(1),offset(3)

        ! *** compute the interpolation coefficients according to the order ***
        !
        ! (Thanks to Erik Schnetter/Numerical Recipes for formulating these so nicely.) 
        !
        ! These formulas are "just" the coefficients of the classical
        ! Lagrange interpolation polynomials along each dimension.
        ! For example, in 1 dimension the unique quadratic passing
        ! through the 3 points {(x0,y0), (x1,y1), (x2,y2)} is: 
        ! ( x-x1)( x-x2)        ( x-x0)( x-x2)        ( x-x0)( x-x1) 
        ! -------------- y0  +  -------------- y1  +  -------------- y2 
        ! (x0-x1)(x0-x2)        (x1-x0)(x1-x2)        (x2-x0)(x2-x1) 
        ! (It's easy to see this: each of the terms is yi if x=xi, or
        ! zero if x=any other xj.)  To get the formulas below, just negate 
        ! each (x-x) factor, and substitute the values xi=i.

        if(order==1) then
!           write(*,*) "HELLO!",numpoints,pointcoords(n,1),pointcoords(n,2),pointcoords(n,3)
           do i=1,3
              !            (x-x1)/(x0-x1) = (offset(i)-1)/(-1)
              coeff(i,1) = 1 - offset(i)
              !            (x-x0)/(x1-x0) = offset(i)
              coeff(i,2) =     offset(i)
           end do
        else if(order==2) then
           do i=1,3
              ! It goes: x0 _y__ x1 ____ x2, where y is the point were interpolating to.
              !recall that x-x0 = offset
              !y0 coefficient: ((x-x1)    *  (x-x2))     / ((x2-x0)*(x1-x0))
              coeff(i,1) = (offset(i)-1) * (offset(i)-2) / (  2  *   1 )
              !y1 coefficient: ((x-x0)   *  (x-x2))      / ((x1-x0)*(x1-x2))
              coeff(i,2) =   (offset(i)) * (offset(i)-2) / (  1  * (-1))
              !y2 coefficient: ((x-x0)  *   ( x-x1))    / ((x2-x1)(x2-x0))
              coeff(i,3) =  (offset(i)) * (offset(i)-1) / ((1) * (2))
           end do
        else if(order==3) then
           do i=1,3
              !Below is basically the same scheme as above, but multiplied by -1 in num & denom.  
              !All of this could be optimized.
              coeff(i,1) = (1-offset(i)) * (2-offset(i)) * (3-offset(i)) / (  3  *   2  *   1 )
              coeff(i,2) = ( -offset(i)) * (2-offset(i)) * (3-offset(i)) / (  2  *   1  * (-1))
              coeff(i,3) = ( -offset(i)) * (1-offset(i)) * (3-offset(i)) / (  1  * (-1) * (-2))
              coeff(i,4) = ( -offset(i)) * (1-offset(i)) * (2-offset(i)) / ((-1) * (-2) * (-3))
           end do
        end if
        
        !if(ll.gt. ... statements ensure that points on interp stencil are not double counted!
        mm = 2
        do k=1,order+1
           nn=kvalue(k)
           if(nn.gt.0) then
              interpvaluek = 0.D0
              do j=1,order+1
                 interpvaluej = 0.D0
                 do i=1,order+1
                    ll=ivalue(i)
                    if(ll.gt.0) then
                       interpvaluej = interpvaluej + func(ll,mm,nn) * coeff(1,i)
                    end if
                 end do
                 interpvaluek = interpvaluek + interpvaluej * coeff(2,j)
              end do
              interpvalues(n) = interpvalues(n) + interpvaluek * coeff(3,k)
           end if
        end do
     end do
  end if


  outputvals = 0.D0
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,&
       interpvalues,outputvals,numpoints,CCTK_VARIABLE_REAL)
end subroutine interpolate_pointset
