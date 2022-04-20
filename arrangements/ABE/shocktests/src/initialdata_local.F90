!-------------------------------------------------------------------
! Here we set up all required variables for shock tests, including:
!  emfields, BSSN variables, primitives, etc.
!-------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine shocktests_initialdata_local(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  real*8                :: dT,dX,dY,dZ,P_max,rho_max,tau_max,rhos_max,rho_fail_max_step,M_fail_step
  real*8                :: xmin,ymin,zmin,xmax,ymax,zmax
  integer               :: proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax
  integer               :: ierr,index,handle,dummy,glob_imax,glob_jmax,glob_kmax
  CCTK_REAL             :: reduction_value

  integer               :: i,j,k
  real*8                :: norm,w_l,w_r,h_l,h_r,e_s_l,e_s_r,ut_l,ut_r

  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  integer :: ONE,ZERO
  parameter(ONE = 1.D0, ZERO = 0.D0)

  ext = cctk_lsh

  if(cctk_lsh(1).le.2*cctk_nghostzones(1)) then
     write(*,*) "ERROR: MUST SET NX>NGHOSTZONESX*2!"
     stop
  else if(cctk_lsh(2).le.2*cctk_nghostzones(2)) then
     write(*,*) "ERROR: MUST SET NY>NGHOSTZONESY*2!"
     stop
  else if(cctk_lsh(3).le.2*cctk_nghostzones(3)) then
     write(*,*) "ERROR: MUST SET NZ>NGHOSTZONESZ*2!"
     stop
  end if

  write(*,*) "Okay... Looks like you set up the grid correctly EXTENTS:",cctk_lsh
  write(*,*) "Okay... Looks like you set up the grid correctly GZS:",cctk_nghostzones

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  w_l = 1.0/sqrt(1.0 - v_l*v_l)
  w_r = 1.0/sqrt(1.0 - v_r*v_r)
  h_l = 1.0 + ( npolyshock + 1.0 ) * P_l / rho_l
  h_r = 1.0 + ( npolyshock + 1.0 ) * P_r / rho_r
  !   e_s_l = rho_l*w_l*(1.0 - h_l*w_l) + P_l
  !   e_s_r = rho_r*w_r*(1.0 - h_r*w_r) + P_r
  !  rhosd(i,j,k)*(enthd(i,j,k)*u_T-1.0)-Pd(i,j,k)
  e_s_l = rho_l*w_l*(h_l*w_l - 1.0) - P_l
  e_s_r = rho_r*w_r*(h_r*w_r - 1.0) - P_r
  ut_l = w_l * h_l * v_l
  ut_r = w_r * h_r * v_r

  ! just to make sure, normalize n:
  !
  norm = sqrt(nx*nx + ny*ny + nz*nz)
  nx = nx/norm
  ny = ny/norm
  nz = nz/norm

  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
           if (X(i,1,1)*nx + Y(1,j,1)*ny + Z(1,1,k)*nz .le. 0.D0) then
              !
              ! "left" side of partition
              !
              P(i,j,k)    = P_l
              rho_b(i,j,k) = rho_l
              rho_star(i,j,k) = rho_l * w_l
              tau(i,j,k)   = e_s_l
              h(i,j,k) = h_l
              w(i,j,k)    = rho_star(i,j,k) *w_l
              st_x(i,j,k)  = nx * rho_star(i,j,k) *ut_l 
              st_y(i,j,k)  = ny * rho_star(i,j,k) *ut_l 
              st_z(i,j,k)  = nz * rho_star(i,j,k) *ut_l 
              vx(i,j,k)   = nx * v_l
              vy(i,j,k)   = ny * v_l
              vz(i,j,k)   = nz * v_l
           else
              !
              ! "right" side of partition
              !
              P(i,j,k)    = P_r
              rho_b(i,j,k) = rho_r
              rho_star(i,j,k) = rho_r * w_r
              tau(i,j,k)   = e_s_r
              h(i,j,k) = h_r
              w(i,j,k)    = rho_star(i,j,k) *w_r
              st_x(i,j,k)  = nx * rho_star(i,j,k) *ut_r 
              st_y(i,j,k)  = ny * rho_star(i,j,k) *ut_r 
              st_z(i,j,k)  = nz * rho_star(i,j,k) *ut_r
              vx(i,j,k)   = nx * v_r
              vy(i,j,k)   = ny * v_r
              vz(i,j,k)   = nz * v_r
           end if
           ! 	  vxd(i,j,k) = 0.5
           ! 	  vyd(i,j,k) = 0.0
           ! 	  vzd(i,j,k) = 0.0
           ! 	  rhosd(i,j,k) = exp(-pow(Xd(i)/0.4,4)) + 1.0
           ! 	  ut_l = 1.0/sqrt(1.0 - pow(vxd(i,j,k),2))
           ! 	  wd(i,j,k)  = rhosd(i,j,k)*ut_l
           ! 	  rhobd(i,j,k) = rhosd(i,j,k) / ut_l
           ! 	  Pd(i,j,k)  = pow(rhobd(i,j,k),Gam)
           ! 	  enthd(i,j,k) = 1.0 + Gam*pow(rhobd(i,j,k),Gam-1.0)
           ! 	  stxd(i,j,k) = rhosd(i,j,k)*enthd(i,j,k)*ut_l*vxd(i,j,k)
           ! 	  styd(i,j,k) = 0.0
           ! 	  stzd(i,j,k) = 0.0
           ! 	  esd(i,j,k)  = rhosd(i,j,k)*(enthd(i,j,k)*ut_l - 1.0) - Pd(i,j,k)
        end do
     end do
  end do

  !tau_p = tau

  neos = 1
  write(*,*) "Assuming polytropic EOS..."
  write(*,*) "Polytropic index = ",npolyshock
  write(*,*) "Polytropic constant K = ",K_poly
  rho_tab(1) = 0.03
  P_tab(1) = K_poly * rho_tab(1)**(1.0+1.0/npolyshock)
  eps_tab(1) = npolyshock*P_tab(1)/rho_tab(1)
  gamma_tab(1) = 1.0 + 1.0/npolyshock
  k_tab(1) = K_poly
  gamma_tab(2) = 1.0 + 1.0/npolyshock 
  k_tab(2) = K_poly

  gamma_th = 1.0 + 1.0/npolyshock 

  u0=1.D0/sqrt(1.D0-(vx*vx+vy*vy+vz*vz))

  mhd_st_x = st_x
  mhd_st_y = st_y
  mhd_st_z = st_z

  !======================================
  ! Set metric data to Minkowski
  !======================================
  phi = ZERO
  psi = exp(phi)
  phix = ZERO
  phiy = ZERO
  phiz = ZERO
  lapm1 = ZERO
  gxx = ONE
  gyy = ONE
  gzz = ONE
  gxy = ZERO
  gxz = ZERO
  gyz = ZERO
  gupxx = ONE
  gupyy = ONE
  gupzz = ONE
  gupxy = ZERO
  gupxz = ZERO
  gupyz = ZERO
  Gammax = ZERO
  Gammay = ZERO
  Gammaz = ZERO 
  !======================================
  ! Set K_ij to Zero!
  !======================================
  trK = ZERO
  Axx = ZERO
  Axy = ZERO
  Axz = ZERO
  Ayy = ZERO
  Ayz = ZERO
  Azz = ZERO
  !======================================
  ! Set everything else to Zero!
  !======================================
  shiftx = ZERO
  shifty = ZERO
  shiftz = ZERO

  !======================================
  ! Set initial EM fields to zero:
  !======================================
  Bx = ZERO
  By = ZERO
  Bz = ZERO
  Bxtilde = ZERO
  Bytilde = ZERO
  Bztilde = ZERO

  Ex = ZERO
  Ey = ZERO
  Ez = ZERO

  sbt = ZERO
  sbx = ZERO
  sby = ZERO
  sbz = ZERO

  !FIXME: read note:
  ! NOTE: The following lines are not strictly correct, but this should NOT matter, since Cowling is turned on!
  Sxx = ZERO
  Sxy = ZERO
  Sxz = ZERO
  Syy = ZERO
  Syz = ZERO
  Szz = ZERO

  Sx = ZERO
  Sy = ZERO
  Sz = ZERO

  rho = ZERO

  write(*,*) "LOCAL INITIALDATA FINISHED!",dx,mhd_st_x(4,4,4)
  write(*,*) "LOCAL INITIALDATA FINISHED! y=",y(1,:,1)

end subroutine shocktests_initialdata_local
