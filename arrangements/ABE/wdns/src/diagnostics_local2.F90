#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
! Diagnostic output for wdns thorn
!-----------------------------------------------------------------------------
subroutine wdns_diagnostics_local2(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer                                  :: i,j,k

  if(MOD(cctk_iteration,out_every)==0) then

     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              temp9(i,j,k)  = -1.D200
              temp10(i,j,k) = -1.D200
              temp11(i,j,k) = -1.D200
           end do
        end do
     end do
     !$omp end parallel do

     write(*,*) "wdns_max_b2 = ",wdns_max_b2

     if(wdns_max_b2.gt.0.D0) then
        !$omp parallel do
        do k=1,cctk_lsh(3)
           do j=1,cctk_lsh(2)
              do i=1,cctk_lsh(1)
                 if(abs((smallb2(i,j,k)-wdns_max_b2)/wdns_max_b2).lt.1e-8) then
                    ! make sure we're out of a ghostzone here:
                    temp9(7,7,7)  = X(i,j,k)
                    temp10(7,7,7) = Y(i,j,k)
                    temp11(7,7,7) = Z(i,j,k)
                 end if
              end do
           end do
        end do
        !$omp end parallel do
     else
        write(*,*) "It looks like the magnetic fields are zero."
     end if
  end if
end subroutine wdns_diagnostics_local2
