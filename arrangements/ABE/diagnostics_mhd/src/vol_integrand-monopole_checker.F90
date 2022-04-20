!---------------------------------------------------------
! Monopole checker: Computes Div(B)
!---------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "../../GenericFD.h"

subroutine monopole_integrand(CCTK_ARGUMENTS)
  implicit none

  !~~~~~> Input parameters:
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  !~~~~~> Other variables:
  real*8                             :: dX, dY, dZ, f1odX, f1odY, f1odZ, Psi6
  real*8                             :: dV,pi,FOUR,ONE
  integer                            :: i,j,k
  integer                            :: AXISYM
  parameter(AXISYM = 4)
  parameter(ONE = 1.D0, FOUR = 4.D0)

  ! 1st of 2 needed #includes for GenericFD.h:
#include "../../GenFD_decl_varF90.h"

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  !     Initialise finite differencing variables
  ! 2nd of 2 needed #includes for GenericFD.h:
#include "../../GenFD_set_varF90.h"

  WhichIntegral = 108

  if(Symmetry.ne.1 .and. Symmetry .ne. 0) then
     do i=1,10
        write(*,*) "WARNING: MONOPOLE CHECKER DOESN'T WORK WITH Symmetry.ne.1 and Symmetry.ne.0!"
     end do
  end if

  if(MOD(cctk_iteration,Compute_VolIntegrands_Every)==0) then
     if (constrained_transport_scheme==3) then
        call compute_absdivb_ct3(cctk_lsh,dX,dY,dZ,phi,Bx_stagger,By_stagger,Bz_stagger,& 
             vxr,vxl,vyr,VolIntegrand)

        ! Exclude region inside BH, as an extra diagnostic:
	if (num_BHs .gt. 0) then 
           !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)

                    if(abs(emask(i,j,k)-1.D0) .lt. 1.D-8) then
                       VolIntegrand2(i,j,k) = VolIntegrand(i,j,k)
                    else
                       VolIntegrand2(i,j,k) = 0.D0
                    end if
                 end do
              end do
           end do
           !$omp end parallel do
         end if
     else
        
        f1odX = 1.D0/dX
        f1odY = 1.D0/dY
        f1odZ = 1.D0/dZ
        
        !-----------------------------------------------------------------------------
        ! Compute integrand
        !-----------------------------------------------------------------------------
        !Reset VolIntegrand quickly, using OpenMP
        !$omp parallel do
        do k=1,cctk_lsh(3)
           do j=1,cctk_lsh(2)
              do i=1,cctk_lsh(1)
                 temp1(i,j,k) = 0.D0
                 temp2(i,j,k) = 0.D0
                 temp3(i,j,k) = 0.D0
                 temp4(i,j,k) = 0.D0
                 temp5(i,j,k) = 0.D0
                 temp6(i,j,k) = 0.D0
              end do
           end do
        end do
        !$omp end parallel do
        
        if(hyperbolic_divergence_cleaning_enable==0 .or. hyperbolic_divergence_cleaning_centered_differencing==0) then

           !In the following lines of code, we ignore the i=1, j=1, and k=1 planes.
           !  This should have no impact on the final result; the only zones
           !  that are actually integrated are ijk>=nghostzonesijk

           !First do derivatives in x,y,z-directions
!!!$omp parallel do
           do k = 2,cctk_lsh(3)
              do j = 2,cctk_lsh(2)
                 do i = 2,cctk_lsh(1)
                    Psi6 = exp(6.d0*phi(i,j,k))
                    temp1(i,j,k) = f1odX*(Psi6*Bx(i,j,k) - exp(6.d0*phi(i-1,j,k))*Bx(i-1,j,k))
                    temp2(i,j,k) = f1odY*(Psi6*By(i,j,k) - exp(6.d0*phi(i,j-1,k))*By(i,j-1,k))
                    temp3(i,j,k) = f1odZ*(Psi6*Bz(i,j,k) - exp(6.d0*phi(i,j,k-1))*Bz(i,j,k-1))
                 end do ! i-loop
              end do ! j-loop
           end do ! k-loop
!!!$omp end parallel do

           !Next, compute averages in y for Bxx and x for Byy and Bzz
           !$omp parallel do
           do k = 2,cctk_lsh(3)
              do j = 2,cctk_lsh(2)
                 do i = 2,cctk_lsh(1)
                    temp4(i,j,k) = 0.5D0*(temp1(i,j,k) + temp1(i,j-1,k))
                    temp5(i,j,k) = 0.5D0*(temp2(i,j,k) + temp2(i-1,j,k))
                    temp6(i,j,k) = 0.5D0*(temp3(i,j,k) + temp3(i-1,j,k))
                 end do ! i-loop
              end do ! j-loop
           end do ! k-loop
           !$omp end parallel do

           !Finally, compute averages in z for Bxx and Byy, and y for Bzz
           !$omp parallel do
           do k = 2,cctk_lsh(3)
              do j = 2,cctk_lsh(2)
                 do i = 2,cctk_lsh(1)
                    temp1(i,j,k) = 0.5D0*(temp4(i,j,k) + temp4(i,j,k-1))
                    temp2(i,j,k) = 0.5D0*(temp5(i,j,k) + temp5(i,j,k-1))
                    temp3(i,j,k) = 0.5D0*(temp6(i,j,k) + temp6(i,j-1,k))
                 end do ! i-loop
              end do ! j-loop
           end do ! k-loop
           !$omp end parallel do

           !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
                    VolIntegrand(i,j,k) = abs(temp1(i,j,k) + temp2(i,j,k) + temp3(i,j,k))
                 end do
              end do
           end do
           !$omp end parallel do

           !here we use vxr,vxl,vyl,vyr, and vzl as temporary storage, and store the measured monopole violation inside VolIntegrand
           call compute_absdivb(cctk_lsh,dX,dY,dZ,X,Z,Bx,By,Bz,phi,VolIntegrand,Sym_Bz,Symmetry,vxr,vxl,vyl,vyr,vzl)
        else
           !$omp parallel do
           do k=2,cctk_lsh(3)-1
              do j=2,cctk_lsh(2)-1
                 do i=2,cctk_lsh(1)-1
                    VolIntegrand(i,j,k) = abs(D1_c2(Bxtilde,i,j,k) + D2_c2(Bytilde,i,j,k) + D3_c2(Bztilde,i,j,k))
                 end do
              end do
           end do
           !$omp end parallel do
        end if
        if(Symmetry == 1 .or. Symmetry == 0) then
           !$omp parallel do
           do k=1,cctk_lsh(3)
              do j=1,cctk_lsh(2)
                 do i=1,cctk_lsh(1)
                    if( &
                         i.le.1 .or. &
                         j.le.1 .or. &
                         k.le.1 .or. &
                         i.gt.cctk_lsh(1) .or. &
                         j.gt.cctk_lsh(2) .or. &
                         k.gt.cctk_lsh(3)) then

!!$                      i.le.cctk_nghostzones(1)+1 .or. &
!!$                      j.le.cctk_nghostzones(2)+1 .or. &
!!$                      k.le.cctk_nghostzones(3)+1 .or. &
!!$                      i.gt.cctk_lsh(1)-cctk_nghostzones(1)-1 .or. &
!!$                      j.gt.cctk_lsh(2)-cctk_nghostzones(2)-1 .or. &
!!$                      k.gt.cctk_lsh(3)-cctk_nghostzones(3)-1) then

                       VolIntegrand(i,j,k) = 0.D0
                    end if
                 end do
              end do
           end do
           !$omp end parallel do
        end if

        !$omp parallel do
        do k=1,cctk_lsh(3)
           do j=1,cctk_lsh(2)
              do i=1,cctk_lsh(1)

                 if(abs(emask(i,j,k)-1.D0) .lt. 1.D-8) then
                    VolIntegrand2(i,j,k) = VolIntegrand(i,j,k)
                 else
                    VolIntegrand2(i,j,k) = 0.D0
                 end if
              end do
           end do
        end do
        !$omp end parallel do


        if(1==0) then
           if(abs(log(real(cctk_levfac(1)))/log(2.D0)) .lt. 1.D-10) then
              !write(*,*) "USING MONOPOLE DATA ON LEVEL ",log(real(cctk_levfac(1)))/log(2.D0)+1.D0
           else
              VolIntegrand = 0.D0
              write(*,*) "WARNING: MONOPOLE DIAGNOSTIC NOT TRUSTWORTHY FOR AMR/FMR SIMULATIONS:"
              write(*,*) "TODO: Exclude all refinement boundary points, plus one more since "
              write(*,*) "the diagnostic stencil reaches into the refinement boundary"
           end if
        end if

        if (Symmetry == AXISYM) then
           VolIntegrand = X*VolIntegrand
        end if
     end if
  end if
end subroutine monopole_integrand
