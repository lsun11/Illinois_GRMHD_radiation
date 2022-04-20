#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bhns_initia_velpert(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext

  real*8 			           :: dT,dX,dY,dZ
  real*8 			           :: n1,n2,n3,mf

  integer :: vindex,handle,ierr
  real*8  :: rho_max,tau_max,valuetosetvz
  real*8  :: ONE,ZERO,W_L,U2, gxxL, gxyL, gxzL, gyyL, gyzL, gzzL,psiL,psiL4
  real*8  :: cn2, cn3, an1, an2, an3
  real*8  :: gupxxL, gupxyL, gupxzL, gupyyL, gupyzL, gupzzL, u0L
  real*8  :: P_deplete, RADEQUAT, delta_bar, Omega_Frame
  real*8  :: rr, re2, ww, phia, pai=3.141592653589793d0
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(ONE = 1.D0, ZERO = 0.D0)
  
  real*8 :: Omega_r

  INTEGER :: i,j,k
  character :: filename*30,c2*2,c1
  CCTK_REAL :: reduction_value

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ext = cctk_lsh


! *********************
! ****   Cook ID   ****
! *********************

! add perturbation in v's

   write(6,*) "****************************** ADDING A PERTURBATION IN VELOCITY ***********************"
   write(6,*) "nperturb=", nperturb
   write(6,*) "lambda_perturb=", lambda_perturb
   write(6,*) "a2oa1_perturb=", a2oa1_perturb

   if ((use_new_bhns_initial_data==5).and.((nperturb==5).or.(nperturb==3))) then
      
      print *, "in bhns_initia_velpert"
      
      do k=1,ext(3)
         do j=1,ext(2)
            do i=1,ext(1)
               vx(i,j,k) = vx(i,j,k) + y(i,j,k)*lambda_perturb/a2oa1_perturb
               vy(i,j,k) = vy(i,j,k) - x(i,j,k)*lambda_perturb*a2oa1_perturb
            end do
         end do
      end do
   end if

end subroutine bhns_initia_velpert
