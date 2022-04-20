#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
! Diagnostic output for bhns thorn
!-----------------------------------------------------------------------------
subroutine bhns_diagnostics_local(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer                                  :: i,j,k

  if(MOD(cctk_iteration,out_every)==0) then

     call bhns_compute_b2_cpp(cctkGH,cctk_lsh, phi, lapm1, &
          shiftx,shifty,shiftz,vx,vy,vz,Bx,By,Bz, & 
          gxx, gxy, gxz, gyy, gyz, gzz, temp8)

     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)

              ! Set b2 to zero inside horizon:
              if(abs(emask(i,j,k)-1.D0) .ge. 1.D-8) then
                 temp8(i,j,k) = 0.D0
              end if
           end do
        end do
     end do
     !$omp end parallel do
  end if
end subroutine bhns_diagnostics_local
