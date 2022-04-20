!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!-----------------------------------------------------------------------------
subroutine find_equatorial_radius(ext,cctkgh,X,Y,Z,dX,dY,dZ,rho_tiny,rho_star,phi,gxx,gxy,gyy,Symmetry,finalvalue,myproc)
  
  implicit none
  
  
  integer, dimension(3)                    :: ext
  real*8                                   :: dT,dX,dY,dZ,finalvalue0,finalvalue
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in):: X,Y,Z
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in):: rho_star,phi,gxx,gxy,gyy

  real*8 :: xradiusvalue,finalxradiusvalue,gphiphi,sintheta,x_r,r_iso,phi_r,gxx_r,gxy_r,gyy_r,Psi4_r
  real*8 :: rho_tiny
  CCTK_POINTER :: cctkgh
  integer :: header_flag,handle,dummy,ivalue,oivalue
  integer :: index,i,j,k,sumoverprocs
  integer :: ierr,myproc,Symmetry
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  xradiusvalue = 1.D111
  finalvalue0 = 0.D0
  finalvalue = 0.D0
  ivalue = 0

  if(Symmetry==OCTANT) then
     j = 2
     k = 2
  else if(Symmetry==EQUATORIAL) then
     j = 0
     do i=1,ext(2)
        if(abs(Y(2,i,2)-dY*0.5).lt.dY*1.D-8) then
           j = i
        end if
     end do
     k = 2
  else if(Symmetry == AXISYM) then
     j = 1
     k = 2
  end if
  
  if(j .gt. 0 .and. abs(Z(2,2,k)-dZ*0.5).lt.dZ*1.D-8) then
     do i=ext(1),1,-1
!           write(*,*) "hi0.",i,j,k,rho_star(i,j,k),rho_tiny,X(i,j,k)
        if(rho_star(i,j,k) .lt. rho_tiny .and. X(i,j,k) .gt. 0.D0) then
           xradiusvalue = X(i,j,k) ! minus dX needed for DAGH compatibility...
!           xradiusvalue = X(i,j,k) - dX ! minus dX needed for DAGH compatibility...
        end if
     end do
  end if
  
!  xradiusvalue = abs(xradiusvalue)

  call CCTK_ReductionHandle(handle,"minimum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,xradiusvalue,finalxradiusvalue,CCTK_VARIABLE_REAL)  

  if(j .gt. 0 .and. abs(Z(2,2,k)-dZ*0.5).lt.dZ*1.D-8) then   
     do i=1,ext(1)
        if(X(i,j,k) .gt. 0.D0 .and. abs(X(i,j,k)-finalxradiusvalue).lt.abs(1.D-6*finalxradiusvalue)) then
           phi_r = phi(i,j,k)
           Psi4_r = exp(4.D0*phi_r)
           gxx_r  = Psi4_r * gxx(i,j,k)
           gxy_r  = Psi4_r * gxy(i,j,k)
           gyy_r  = Psi4_r * gyy(i,j,k)
           r_iso = sqrt(X(i,j,k)*X(i,j,k) + Y(i,j,k)*Y(i,j,k) + Z(i,j,k)*Z(i,j,k))
           sintheta = X(i,j,k)/r_iso
           gphiphi = Y(i,j,k)*Y(i,j,k)*gxx_r + X(i,j,k)*X(i,j,k)*gyy_r - 2.D0*X(i,j,k)*Y(i,j,k)*gxy_r
           finalvalue0 = sqrt(gphiphi)/sintheta
           
           ivalue = ivalue + 1
        end if
     end do
  end if
  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,finalvalue0,finalvalue,CCTK_VARIABLE_REAL)  
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,ivalue,oivalue,CCTK_VARIABLE_INT)
  
  if(myproc .eq. 0) then
     if(oivalue .eq. 0) then
        oivalue = oivalue + 1
        write(*,*) "owwweee!"
     end if
     finalvalue = finalvalue/oivalue
  end if

end subroutine find_equatorial_radius
