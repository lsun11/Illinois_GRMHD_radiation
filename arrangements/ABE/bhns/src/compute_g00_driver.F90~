#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


subroutine compute_g00_driver(CCTK_ARGUMENTS)
   implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  integer               :: i,j,k

  real*8                :: psiL,psiL4
  real*8                :: lapseL,betaxL,betayL,betazL
  real*8                :: gxxL,gxyL,gxzL,gyyL,gyzL,gzzL


     do k=1,ext(3)
        do j=1,ext(2)
           do i=1,ext(1)
              psiL  =  psi(i,j,k)
              psiL4 =  psiL*psiL*psiL*psiL
              
              lapseL  =  lapm1(i,j,k) + 1.0d0
              betaxL  =  shiftx(i,j,k)
              betayL  =  shifty(i,j,k)
              betazL  =  shiftz(i,j,k)

               ! Full 3-metric
               gxxL  =  psiL4*gxx(i,j,k)
               gxyL  =  psiL4*gxy(i,j,k)
               gxzL  =  psiL4*gxz(i,j,k)
               gyyL  =  psiL4*gyy(i,j,k)
               gyzL  =  psiL4*gyz(i,j,k)
               gzzL  =  psiL4*gzz(i,j,k)


               ! g00
               temp_g00(i,j,k) = psiL4*(gxxL*betaxL*betaxL  +         &
                                 gyyL*betayL*betayL  +         &
                                 gzzL*betazL*betazL  +         &
                          2.0d0*(gxyL*betaxL*betayL  +         &
                                 gxzL*betaxL*betazL  +         &
                                 gyzL*betayL*betazL)) - (lapseL*lapseL)


           end do
        end do
     end do

  
  end subroutine compute_g00_driver


