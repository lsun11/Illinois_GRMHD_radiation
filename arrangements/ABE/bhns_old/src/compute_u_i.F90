#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!------------------------------------------
!  Particle tracer for BHNS thorn.
!------------------------------------------

! compute lower components of u_i used to compute the 
! circulation on the z=0 plane
subroutine bhns_diagnostics_tracer_compute_u_i(CCTK_ARGUMENTS)
 
   implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  

! compute variables only when need them
  
  if(mod(cctk_iteration,out_every)==0) then

     ! store  u_x
     temp16  = u0*(gxx*(vx + shiftx) + &
           gxy*(vy + shifty) + &
           gxz*(vz + shiftz))
     ! store  u_y
     temp17  = u0*(gxy*(vx + shiftx) + &
           gyy*(vy + shifty) + &
           gyz*(vz + shiftz))

     temp18  = 1.0d0/u0
  end if
  
end subroutine bhns_diagnostics_tracer_compute_u_i
