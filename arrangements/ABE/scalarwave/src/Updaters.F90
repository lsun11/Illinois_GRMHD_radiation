!-----------------------------------------------------------------------------
!
!$Id: Updaters.F90,v 1.2 2006/03/10 22:32:24 zetienne Exp $
!
!-----------------------------------------------------------------------------
!
! Fortran evolution schemes for interior (Predictor & Corrector subroutines)
!
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
!
! Predictor
!
!-----------------------------------------------------------------------------
subroutine ABE_Predictor(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
 
  real*8                                       :: dX, dY, dZ, dT
  integer                    :: i, j, k,temp
  integer                    :: imin, jmin, kmin, imax, jmax, kmax 
  CCTK_REAL                     :: phi_rhs, phidot_rhs
  CCTK_REAL                     :: ONE, TWO, HALF
  integer, dimension(3)         :: ext
  parameter(ONE = 1.D0, TWO = 2.D0, HALF = 0.5D0)

! Next, properly set the min. array values, depending on whether we are near an 
! axis.  This properly accounts for the axis ghost zones.
!  if(x(lbound(phi_p,1),1,1) < 0) then
!     imin = lbound(phi_p,1)+1
!  else 
!     imin = lbound(phi_p,1)
!  endif
!  if(y(1,lbound(phi_p,2),1) < 0) then
!     jmin = lbound(phi_p,2)+1
!  else 
!     jmin = lbound(phi_p,2)
!  endif
!  if(z(1,1,lbound(phi_p,3)) < 0) then
!     kmin = lbound(phi_p,3)+1
!  else 
!     kmin = lbound(phi_p,3)
!  endif

  imin = lbound(phi_p,1)
  jmin = lbound(phi_p,2)
  kmin = lbound(phi_p,3)

  imax=cctk_lsh(1)-1
  jmax=cctk_lsh(2)-1
  kmax=cctk_lsh(3)-1

  ext(1) = cctk_lsh(1)
  ext(2) = cctk_lsh(2)
  ext(3) = cctk_lsh(3)
      
  dT = CCTK_DELTA_TIME
  dX=X(2,1,1)-X(1,1,1)
  dY=Y(1,2,1)-Y(1,1,1)
  dZ=Z(1,1,2)-Z(1,1,1)

  write (*,*) "okayX... ",X(1,1,1),X(2,1,1),X(3,1,1)
  write (*,*) "okayY... ",Y(1,1,1),Y(1,2,1),Y(1,3,1)
  write (*,*) "okayZ... ",Z(1,1,1),Z(1,1,2),Z(1,1,3)

!  write (*,*) "okay1... ", dT,dX,dY,dZ,CCTK_DELTA_TIME,CCTK_DELTA_SPACE(1),CCTK_DELTA_SPACE(2),CCTK_DELTA_SPACE(3)
!  write (*,*) "okay2... ", dT,dX,dY,dZ,cctk_delta_time,cctk_delta_space(1),cctk_delta_space(2),cctk_delta_space(3)

! Compute Laplacian of phi:
!
  call Laplace(phi_p,DD_phi_p,ext,X,Y,Z,dX,dY,dZ,1,1,1);

  do k = kmin,kmax
     do j = jmin,jmax
        do i = imin,imax
           phi(i,j,k) = phi_p(i,j,k) + dT*phidot_p(i,j,k)
           phidot(i,j,k) = phidot_p(i,j,k) + dT*DD_phi_p(i,j,k)
        enddo
     enddo
  enddo

! If domain is set to octant, set octant properly.  Otherwise use axisymmetry.
  if (CCTK_Equals(domain,"octant").eq.1) then
      Symmetry = 0
  else
      Symmetry = 1
  endif
 
  if (Symmetry.eq.0) then
     call scalarwave_update_boundary(ext,X,Y,Z,dT,phi_p,phi,phidot_p,phidot,Symmetry)
  endif
  if (Symmetry.eq.1) then
      call BndCartoon2DVN(temp, cctkGH, 0,'ABE_scalar_wave::phi') 
      call BndCartoon2DVN(temp, cctkGH, 0,'ABE_scalar_wave::phidot') 
      call scalarwave_update_boundary(ext,X,Y,Z,dT,phi_p,phi,phidot_p,phidot,Symmetry)
  endif
  call CartSymVN(temp,cctkGH,'ABE_scalar_wave::phi')
  call CartSymVN(temp,cctkGH,'ABE_scalar_wave::phidot')


end subroutine ABE_Predictor

!===============================================================================
! Next, the Corrector step...
!===============================================================================

subroutine ABE_Corrector(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
 
  real*8                                       :: dX, dY, dZ, dT
  integer                    :: i, j, k, temp
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  CCTK_REAL                     :: ONE, TWO, HALF
  CCTK_REAL			:: rr,t, analy
  integer, dimension(3)         :: ext
  parameter(ONE = 1.D0, TWO = 2.D0, HALF = 0.5D0)

! Next, properly set the min. array values, depending on whether we are near an 
! axis.  This properly accounts for the axis ghost zones.
!  if(x(lbound(phi_p,1),1,1) < 0) then
!     imin = lbound(phi_p,1)+1
!  else 
!     imin = lbound(phi_p,1)
!  endif
!  if(y(1,lbound(phi_p,2),1) < 0) then
!     jmin = lbound(phi_p,2)+1
!  else 
!     jmin = lbound(phi_p,2)
!  endif
!  if(z(1,1,lbound(phi_p,3)) < 0) then
!     kmin = lbound(phi_p,3)+1
!  else 
!     kmin = lbound(phi_p,3)
!  endif

  imin = lbound(phi_p,1)
  jmin = lbound(phi_p,2)
  kmin = lbound(phi_p,3)

  imax=cctk_lsh(1)-1 
  jmax=cctk_lsh(2)-1
  kmax=cctk_lsh(3)-1

  ext(1) = cctk_lsh(1)
  ext(2) = cctk_lsh(2)
  ext(3) = cctk_lsh(3)

  dT = CCTK_DELTA_TIME
  dX=X(2,1,1)-X(1,1,1)
  dY=Y(1,2,1)-Y(1,1,1)
  dZ=Z(1,1,2)-Z(1,1,1)

! If domain is set to octant, set octant properly.  Otherwise use axisymmetry.
  if (CCTK_Equals(domain,"octant").eq.1) then
      Symmetry = 0
  else
      Symmetry = 1
  endif

! Compute Laplacian of phi_p (before predictor phi)
  call Laplace(phi_p,DD_phi_p,ext,X,Y,Z,dX,dY,dZ,1,1,1);

  if(corrector_iteration.eq.1) then
     call Laplace(phi,DD_phi,ext,X,Y,Z,dX,dY,dZ,1,1,1);
     do k = kmin,kmax
        do j = jmin,jmax
           do i = imin,imax
              phi_t(i,j,k) = phi_p(i,j,k) + dT*(HALF*(phidot_p(i,j,k) + phidot(i,j,k)))
              phidot_t(i,j,k) = phidot_p(i,j,k) + dT*(HALF*(DD_phi_p(i,j,k) + DD_phi(i,j,k)))
           enddo
        enddo
     enddo
     if (Symmetry.eq.0) then
        call scalarwave_update_boundary(ext,X,Y,Z,dT,phi_p,phi_t,phidot_p,phidot_t,Symmetry)
     endif
     if (Symmetry.eq.1) then
        call BndCartoon2DVN(temp, cctkGH, 0,'ABE_scalar_wave::phi_t') 
        call BndCartoon2DVN(temp, cctkGH, 0,'ABE_scalar_wave::phidot_t') 
        call scalarwave_update_boundary(ext,X,Y,Z,dT,phi_p,phi_t,phidot_p,phidot_t,Symmetry)
     endif
     call CartSymVN(temp,cctkGH,'ABE_scalar_wave::phi_t')
     call CartSymVN(temp,cctkGH,'ABE_scalar_wave::phidot_t')
  else
     call Laplace(phi_t,DD_phi,ext,X,Y,Z,dX,dY,dZ,1,1,1);
     do k = kmin,kmax
        do j = jmin,jmax
           do i = imin,imax
              phi(i,j,k) = phi_p(i,j,k) + dT*(HALF*(phidot_p(i,j,k) + phidot_t(i,j,k)))
              phidot(i,j,k) = phidot_p(i,j,k) + dT*(HALF*(DD_phi_p(i,j,k) + DD_phi(i,j,k)))
           enddo
        enddo
     enddo
     if (Symmetry.eq.0) then
        call scalarwave_update_boundary(ext,X,Y,Z,dT,phi_p,phi,phidot_p,phidot,Symmetry)
     endif
     if (Symmetry.eq.1) then
        call BndCartoon2DVN(temp, cctkGH, 0,'ABE_scalar_wave::phi') 
        call BndCartoon2DVN(temp, cctkGH, 0,'ABE_scalar_wave::phidot') 
        call scalarwave_update_boundary(ext,X,Y,Z,dT,phi_p,phi,phidot_p,phidot,Symmetry)
     endif
     call CartSymVN(temp,cctkGH,'ABE_scalar_wave::phi')
     call CartSymVN(temp,cctkGH,'ABE_scalar_wave::phidot')
  endif

! Increment corrector iteration.
  corrector_iteration = corrector_iteration + 1

end subroutine ABE_Corrector

! Reset the corrector iteration at the end of each timestep:
subroutine ABE_EndStep(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  corrector_iteration = 1

end subroutine ABE_EndStep
