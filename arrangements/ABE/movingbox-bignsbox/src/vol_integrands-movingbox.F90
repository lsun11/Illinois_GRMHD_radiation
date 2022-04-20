!----------------------------
! Set up rest mass integrand
!----------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine compute_movingbox_vol_integrands(CCTK_ARGUMENTS)
  implicit none

  !~~~~~> Input parameters:
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !~~~~~> Other variables:
  real*8                             :: xL,yL,zL,rho_starL
  real*8                             :: dX, dY, dZ
  integer                            :: i,j,k,l
  integer                            :: AXISYM
  parameter(AXISYM = 4)

  WhichIntegral = 1351+current_center_to_track

  if(num_matter_refinement_box_centers.eq.0) then
     write(*,*) "ERROR.  You must disable the movingbox thorn if num_matter_refinement_box_centers==0!"
     stop
  end if

  if(MOD(cctk_iteration,Compute_VolIntegrands_Every)==0) then

     !First reset the current_center_to_track, if needed.
     if(current_center_to_track==num_matter_refinement_box_centers+1) current_center_to_track=1

     !-----------------------------------------------------------------------------
     ! Compute integrand
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
           end do
        end do
     end do
     !$omp end parallel do

     write(*,*) "current center = ",current_center_to_track,position_x(current_center_to_track),position_y(current_center_to_track),position_z(current_center_to_track)
     

!!!$omp parallel do
     do k = 1,cctk_lsh(3)
        do j = 1,cctk_lsh(2)
           do i = 1,cctk_lsh(1)
              xL = X(i,j,k)
              yL = Y(i,j,k)
              zL = Z(i,j,k)

              if( current_center_to_track.ne.1 .or. &
                   (current_center_to_track==1.and. &
                    sqrt((position_x(1)-position_x(2))**2 + (position_y(1)-position_y(2))**2 + (position_z(1)-position_z(2))**2) .le. Finest_Refinement_Box_Radius*2.D0)) then
                 
                 !Note that the first matter box is number 2.
                 if(  xL.ge.position_x(current_center_to_track)-3.D0*Finest_Refinement_Box_Radius .and. &
                      xL.le.position_x(current_center_to_track)+3.D0*Finest_Refinement_Box_Radius .and. &
                      yL.ge.position_y(current_center_to_track)-3.D0*Finest_Refinement_Box_Radius .and. &
                      yL.le.position_y(current_center_to_track)+3.D0*Finest_Refinement_Box_Radius .and. &
                      zL.ge.position_z(current_center_to_track)-3.D0*Finest_Refinement_Box_Radius .and. &
                      zL.le.position_z(current_center_to_track)+3.D0*Finest_Refinement_Box_Radius) then
                    
                    rho_starL = rho_star(i,j,k)
                    
                    box_x(i,j,k) = rho_starL*xL
                    box_y(i,j,k) = rho_starL*yL
                    box_z(i,j,k) = rho_starL*zL
                    box_denom(i,j,k) = rho_starL

!                    write(*,*) "hi",rho_starL
                 end if
              else
                    rho_starL = rho_star(i,j,k)
                    
                    box_x(i,j,k) = rho_starL*xL
                    box_y(i,j,k) = rho_starL*yL
                    box_z(i,j,k) = rho_starL*zL
                    box_denom(i,j,k) = rho_starL
              end if
           end do ! i-loop
        end do ! j-loop
     end do ! k-loop
!!!$omp end parallel do


!!!$omp parallel do
     do k = 1,cctk_lsh(3)
        do j = 1,cctk_lsh(2)
           do i = 1,cctk_lsh(1)
              xL = X(i,j,k)
              yL = Y(i,j,k)
              zL = Z(i,j,k)

              !bh1_distance=sqrt((x(i,j,k)-bh_posn_x(1))**2+(y(i,j,k)-bh_posn_y(1))**2+(z(i,j,k)-bh_posn_z(1))**2) 

              do l = 1,current_center_to_track
                 ! 
                 if(  xL.ge.position_x(l)-Finest_Refinement_Box_Radius .and. &
                      xL.le.position_x(l)+Finest_Refinement_Box_Radius .and. &
                      yL.ge.position_y(l)-Finest_Refinement_Box_Radius .and. &
                      yL.le.position_y(l)+Finest_Refinement_Box_Radius .and. &
                      zL.ge.position_z(l)-Finest_Refinement_Box_Radius .and. &
                      zL.le.position_z(l)+Finest_Refinement_Box_Radius) then
                    box_x(i,j,k) = 0.D0
                    box_y(i,j,k) = 0.D0
                    box_z(i,j,k) = 0.D0
                    box_denom(i,j,k) = 0.D0
                 end if
              end do ! l-loop
           end do ! i-loop
        end do ! j-loop
     end do ! k-loop
!!!$omp end parallel do

  end if
end subroutine Compute_movingbox_vol_integrands
