#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine em_extraction_postregrid(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: i,j,k,ierr
  integer,dimension(3) :: ones
  real*8  :: dX,dY,dZ
  real*8  :: psi_4thpower, gxxL, gyyL, gzzL, gxyL, gxzL, gyzL, xL, yL, zL

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  write(*,*) "howdy...in em_extraction_postregrid",dX,cctk_lsh


  call em_extraction_lowlevel(cctkGH,cctk_lsh,cctk_nghostzones, &
       dx,dy,dz,r,x,y,z,gxx,gxy,gxz,gyy,gyz,gzz,lapm1,shiftx,shifty,shiftz, &
       phi,psi,Bx,By,Bz,Ex,Ey,Ez,phi0re,phi0im,phi2re,phi2im,set_up_spherical_EM_wave_in_flat_spacetime,compute_phi2_min_radius,compute_phi2_max_radius,scale_with_radius_em,radius_power_em)
  
  ones = 1

  ! Vasilis says: I 'm not so sure about the following, but it's probably OK

  call CartSymGN(ierr,cctkGH,"em_extraction::phi2_em")


end subroutine em_extraction_postregrid
