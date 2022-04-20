!----------------------------
! Set up rest mass integrand
!----------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine set_rho_b_old(CCTK_ARGUMENTS)
  implicit none

  !~~~~~> Input parameters:
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !~~~~~> Other variables:
  integer                            :: i,j,k


  if(MOD(cctk_iteration,out_every)==0) then
     !-----------------------------------------------------------------------------
     !
     !-----------------------------------------------------------------------------
     write(*,*) "*******************************"
     write(*,*) "setting rho_b_old"
     write(*,*) "*******************************"
     !Reset quickly, using OpenMP
     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
             rho_b_old(i,j,k)=rho_b(i,j,k)
           end do
        end do
     end do
     !$omp end parallel do
  endif
END subroutine set_rho_b_old
