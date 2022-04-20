!------------------------------------------------------------
! Master update boundary condition (bc) driver for the lapse
!------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine lapse_postbc(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer                                       :: dummy, i,j,k
  integer, dimension(3)                         :: ext
  real*8                                        :: dX,dY,dZ
!
  dX=CCTK_DELTA_SPACE(1)
  dY=CCTK_DELTA_SPACE(2)
  dZ=CCTK_DELTA_SPACE(3)


  ext = cctk_lsh
  if(cctk_iteration.gt.0) then
     call fill_lapse_symmetry_gz(ext,X,Y,Z,Symmetry,lapm1)
     call fill_lapsederivs_symmetry_gz(ext,X,Y,Z,Symmetry,lapsex,lapsey,lapsez)
     call Derivs(ext,X,Y,Z,dX,dY,dZ,lapm1,lapsex,lapsey,lapsez,Symmetry)
     call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,lapm1,lapsex,lapsey,lapsez)
     !FIXME: POSSIBLE BUG HERE.  Check output with following line commented/uncommented:
     !call fill_lapsederivs_symmetry_gz(ext,X,Y,Z,Symmetry,lapsex,lapsey,lapsez)
  end if

if(1==0) then
  do k=1, cctk_lsh(3)
     do j=1, cctk_lsh(2)
          do i=1, cctk_lsh(1)
if (i==78.and.j==70.and.k==2) then
    write(*,*) " In lapse post_bcs checkpoint <<2>>  tau(i,j,k) is ", tau(i,j,k)
    write(*,*) " In lapse post_bcs checkpoint <<2>>  tau_rhs(i,j,k) is ", tau_rhs(i,j,k)
    write(*,*) " In lapse post_bcs checkpoint <<2>>  tau_rad(i,j,k) is ", tau_rad(i,j,k)
    write(*,*) " In lapse post_bcs checkpoint <<2>>  tau_rad_rhs(i,j,k) is ", tau_rad_rhs(i,j,k)
end if
        end do
    end do
  end do
end if


end subroutine lapse_postbc
