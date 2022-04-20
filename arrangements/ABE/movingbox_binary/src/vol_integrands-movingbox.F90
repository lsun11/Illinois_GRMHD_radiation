!----------------------------
! Set up rest mass integrand
!----------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine compute_movingbox_vol_integrands_bin(CCTK_ARGUMENTS)
  implicit none

  !~~~~~> Input parameters:
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !~~~~~> Other variables:
  real*8                             :: xL,yL,zL,rho_starL, xo, yo, zo
  real*8                             :: dX, dY, dZ
  integer                            :: i,j,k,l
  integer                            :: AXISYM
  parameter(AXISYM = 4)

  !This could be useful if we want to make multiple boxes track the NS matter independently.
  !  WARNING: Try to avoid placing refinement boxes in such a way that their movement is erratic: 
  !           You'll end up with no improvement in your results, due to time prolongation.
!  WhichIntegral = 1351

  if(.not.(num_matter_refinement_box_centers.eq.1 .or. num_matter_refinement_box_centers.eq.0)) then
     write(*,*) "ERROR.  You must set num_matter_refinement_box_centers.eq.1 or num_matter_refinement_box_centers.eq.0 when the movingbox thorn is enabled!"
     stop
  end if

  if(MOD(cctk_iteration,Compute_VolIntegrands_Every)==0 .and. cctk_iteration .gt. 1) then

     write(*,*) "HOWDY INSIDE compute_movingbox_vol_integrands_bin!"

     !-----------------------------------------------------------------------------
     ! Compute integrands
     !-----------------------------------------------------------------------------
     !Reset VolIntegrand's quickly, using OpenMP
     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              box_x(i,j,k) = 0.D0
              box_y(i,j,k) = 0.D0
              box_z(i,j,k) = 0.D0
              box_denom(i,j,k) = 0.D0

              temp1(i,j,k) = 0.D0
              temp2(i,j,k) = 0.D0
              temp3(i,j,k) = 0.D0
              temp4(i,j,k) = 0.D0
           end do
        end do
     end do
     !$omp end parallel do


     xo = position_x(1)
     yo = position_y(1)
     zo = position_z(1)


     !do k = 1,cctk_lsh(3)
     !   do j = 1,cctk_lsh(2)
     !      do i = 1,cctk_lsh(1)
     !$omp parallel do
     do k=cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
        do j=cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
           do i=cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
              xL = X(i,j,k)
              yL = Y(i,j,k)
              zL = Z(i,j,k)
              
              if ( sqrt( (xL - xo)*(xL - xo)+(yL - yo)*(yL - yo)+(zL - zo)*(zL - zo)) <= radius_star1) then
                 rho_starL = rho_star(i,j,k)
                 
                 box_x(i,j,k) = rho_starL*xL
                 box_y(i,j,k) = rho_starL*yL
!                 box_z(i,j,k) = rho_starL*zL
                 box_denom(i,j,k) = rho_starL                               
              endif



              if(rho_star(i,j,k) .ne. rho_star(i,j,k)) then
                 write(*,*) "NAN FOUND IN RHO_STAR AT i,j,k=",i,j,k,rho_star(i,j,k)
              end if
              
              if(x(i,j,k) .ne. x(i,j,k)) then
                 write(*,*) "NAN FOUND IN X AT i,j,k=",i,j,k,x(i,j,k)
              end if

           end do ! i-loop
        end do ! j-loop
     end do ! k-loop
     !$omp end parallel do


     xo = position_x(2)
     yo = position_y(2)
     zo = position_z(2)
     !$omp parallel do
     do k=cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
        do j=cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
           do i=cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
              xL = X(i,j,k)
              yL = Y(i,j,k)
              zL = Z(i,j,k)
              
              if ( sqrt( (xL - xo)*(xL - xo)+(yL - yo)*(yL - yo)+(zL - zo)*(zL - zo) ) <= radius_star2) then
                 rho_starL = rho_star(i,j,k)
                 
                 temp1(i,j,k) = rho_starL*xL
                 temp2(i,j,k) = rho_starL*yL
!                 temp3(i,j,k) = rho_starL*zL
                 temp4(i,j,k) = rho_starL                 
              endif

           end do ! i-loop
        end do ! j-loop
     end do ! k-loop
     !$omp end parallel do

  end if
end subroutine compute_movingbox_vol_integrands_bin
