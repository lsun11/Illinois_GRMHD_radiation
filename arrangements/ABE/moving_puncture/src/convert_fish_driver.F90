!-----------------------------------------------------------------------------
!
! Convert various quantities from physical coords to fisheye
!
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine movpunc_convert_fish(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext
  real*8 			           :: dT,dX,dY,dZ
  real*8 :: xmin,ymin,zmin,rho_fail_max_step,m_fail_step,reduction_value
  real*8 :: rho_max,xx,yy,zz,rr,xmax2,ymax2,zmax2
  integer :: i,j,l
  integer :: handle,dummy,proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax,glob_imax,glob_jmax,glob_kmax
  integer :: index
  integer :: ierr
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  ext = cctk_lsh

  if(fisheye_enable.eq.1) then
     call trans_phys_fish_tensor_flat(ext,X,Y,Z, &
          PhysicalRadius,RadiusDerivative, &
          gxx,gxy,gxz,gyy,gyz,gzz,Symmetry)
     call trans_phys_fish_tensor_inv(ext,X,Y,Z, &
          PhysicalRadius,RadiusDerivative, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,Symmetry)
     call trans_phys_fish_tensor(ext,X,Y,Z, &
          PhysicalRadius,RadiusDerivative, &
          Axx,Axy,Axz,Ayy,Ayz,Azz,Symmetry)
     call trans_phys_fish_phi(ext,X,Y,Z, &
          PhysicalRadius,RadiusDerivative,phi,Symmetry)
     call trans_phys_fish_gamt_flat(ext,X,Y,Z, &
          PhysicalRadius,RadiusDerivative,RadiusDerivative2, &
          Gammax,Gammay,Gammaz,Symmetry)
  endif

end subroutine movpunc_convert_fish
