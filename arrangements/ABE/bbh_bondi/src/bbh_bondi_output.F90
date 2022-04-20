#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!------------------------------------
! Do bbh_bondi diagnostics
!------------------------------------
subroutine bbh_bondi_output(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  INTEGER            :: i,j,k
  INTEGER            :: ntot,n,vindex
  REAL*8             :: output_logradmin_rot_interp_1,output_logradmin_rot_interp_2,output_logradmin_rot_interp_3
  REAL*8             :: output_logradmax_rot_interp_1,output_logradmax_rot_interp_2,output_logradmax_rot_interp_3
  REAL*8             :: PI,thetamin,thetamax,phimin,phimax
  REAL*8             :: output_lograd1L,output_lograd2L,output_lograd3L
  REAL*8             :: phiL,thetaL,x1L,y1L,z1L,x2L,y2L,z2L,x3L,y3L,z3L
  REAL*8             :: BH1_posx,BH1_posy,BH2_posx,BH2_posy

  REAL*8, allocatable, dimension(:,:)       :: pointcoords1,pointcoords2,pointcoords3
  real*8, allocatable, dimension(:)         :: rho_b_int1,rho_b_int2,rho_b_int3
  real*8, allocatable, dimension(:)         :: P_int1,P_int2,P_int3
  real*8, allocatable, dimension(:)         :: vx_int1,vx_int2,vx_int3
  real*8, allocatable, dimension(:)         :: vy_int1,vy_int2,vy_int3
  real*8, allocatable, dimension(:)         :: vz_int1,vz_int2,vz_int3
  real*8, allocatable, dimension(:)         :: u0_int1,u0_int2,u0_int3
  real*8, allocatable, dimension(:)         :: lapse_rhs_int1,lapse_rhs_int2,lapse_rhs_int3
  real*8, allocatable, dimension(:)         :: shiftx_rhs_int1,shiftx_rhs_int2,shiftx_rhs_int3
  real*8, allocatable, dimension(:)         :: shifty_rhs_int1,shifty_rhs_int2,shifty_rhs_int3
  real*8, allocatable, dimension(:)         :: shiftz_rhs_int1,shiftz_rhs_int2,shiftz_rhs_int3
  real*8, allocatable, dimension(:)         :: phi_rhs_int1,phi_rhs_int2,phi_rhs_int3
  real*8, allocatable, dimension(:)         :: gxx_rhs_int1,gxx_rhs_int2,gxx_rhs_int3
  real*8, allocatable, dimension(:)         :: gxy_rhs_int1,gxy_rhs_int2,gxy_rhs_int3
  real*8, allocatable, dimension(:)         :: gxz_rhs_int1,gxz_rhs_int2,gxz_rhs_int3
  real*8, allocatable, dimension(:)         :: gyy_rhs_int1,gyy_rhs_int2,gyy_rhs_int3
  real*8, allocatable, dimension(:)         :: gyz_rhs_int1,gyz_rhs_int2,gyz_rhs_int3
  real*8, allocatable, dimension(:)         :: gzz_rhs_int1,gzz_rhs_int2,gzz_rhs_int3
  real*8, allocatable, dimension(:)         :: lapse_int1,lapse_int2,lapse_int3
  real*8, allocatable, dimension(:)         :: lapsex_int1,lapsex_int2,lapsex_int3
  real*8, allocatable, dimension(:)         :: lapsey_int1,lapsey_int2,lapsey_int3
  real*8, allocatable, dimension(:)         :: lapsez_int1,lapsez_int2,lapsez_int3
  real*8, allocatable, dimension(:)         :: shiftx_int1,shiftx_int2,shiftx_int3
  real*8, allocatable, dimension(:)         :: shifty_int1,shifty_int2,shifty_int3
  real*8, allocatable, dimension(:)         :: shiftz_int1,shiftz_int2,shiftz_int3
  real*8, allocatable, dimension(:)         :: shiftxx_int1,shiftxx_int2,shiftxx_int3
  real*8, allocatable, dimension(:)         :: shiftyx_int1,shiftyx_int2,shiftyx_int3
  real*8, allocatable, dimension(:)         :: shiftzx_int1,shiftzx_int2,shiftzx_int3
  real*8, allocatable, dimension(:)         :: shiftxy_int1,shiftxy_int2,shiftxy_int3
  real*8, allocatable, dimension(:)         :: shiftyy_int1,shiftyy_int2,shiftyy_int3
  real*8, allocatable, dimension(:)         :: shiftzy_int1,shiftzy_int2,shiftzy_int3
  real*8, allocatable, dimension(:)         :: shiftxz_int1,shiftxz_int2,shiftxz_int3
  real*8, allocatable, dimension(:)         :: shiftyz_int1,shiftyz_int2,shiftyz_int3
  real*8, allocatable, dimension(:)         :: shiftzz_int1,shiftzz_int2,shiftzz_int3
  real*8, allocatable, dimension(:)         :: phi_int1,phi_int2,phi_int3
  real*8, allocatable, dimension(:)         :: phix_int1,phiy_int1,phiz_int1
  real*8, allocatable, dimension(:)         :: phix_int2,phiy_int2,phiz_int2
  real*8, allocatable, dimension(:)         :: phix_int3,phiy_int3,phiz_int3
  real*8, allocatable, dimension(:)         :: gxx_int1,gxx_int2,gxx_int3
  real*8, allocatable, dimension(:)         :: gxy_int1,gxy_int2,gxy_int3
  real*8, allocatable, dimension(:)         :: gxz_int1,gxz_int2,gxz_int3
  real*8, allocatable, dimension(:)         :: gyy_int1,gyy_int2,gyy_int3
  real*8, allocatable, dimension(:)         :: gyz_int1,gyz_int2,gyz_int3
  real*8, allocatable, dimension(:)         :: gzz_int1,gzz_int2,gzz_int3
  real*8, allocatable, dimension(:)         :: gupxx_int1,gupxx_int2,gupxx_int3
  real*8, allocatable, dimension(:)         :: gupxy_int1,gupxy_int2,gupxy_int3
  real*8, allocatable, dimension(:)         :: gupxz_int1,gupxz_int2,gupxz_int3
  real*8, allocatable, dimension(:)         :: gupyy_int1,gupyy_int2,gupyy_int3
  real*8, allocatable, dimension(:)         :: gupyz_int1,gupyz_int2,gupyz_int3
  real*8, allocatable, dimension(:)         :: gupzz_int1,gupzz_int2,gupzz_int3
  real*8, allocatable, dimension(:)         :: gxxx_int1,gxxy_int1,gxxz_int1
  real*8, allocatable, dimension(:)         :: gxyx_int1,gxyy_int1,gxyz_int1
  real*8, allocatable, dimension(:)         :: gxzx_int1,gxzy_int1,gxzz_int1
  real*8, allocatable, dimension(:)         :: gyyx_int1,gyyy_int1,gyyz_int1
  real*8, allocatable, dimension(:)         :: gyzx_int1,gyzy_int1,gyzz_int1
  real*8, allocatable, dimension(:)         :: gzzx_int1,gzzy_int1,gzzz_int1
  real*8, allocatable, dimension(:)         :: gxxx_int2,gxxy_int2,gxxz_int2
  real*8, allocatable, dimension(:)         :: gxyx_int2,gxyy_int2,gxyz_int2
  real*8, allocatable, dimension(:)         :: gxzx_int2,gxzy_int2,gxzz_int2
  real*8, allocatable, dimension(:)         :: gyyx_int2,gyyy_int2,gyyz_int2
  real*8, allocatable, dimension(:)         :: gyzx_int2,gyzy_int2,gyzz_int2
  real*8, allocatable, dimension(:)         :: gzzx_int2,gzzy_int2,gzzz_int2
  real*8, allocatable, dimension(:)         :: gxxx_int3,gxxy_int3,gxxz_int3
  real*8, allocatable, dimension(:)         :: gxyx_int3,gxyy_int3,gxyz_int3
  real*8, allocatable, dimension(:)         :: gxzx_int3,gxzy_int3,gxzz_int3
  real*8, allocatable, dimension(:)         :: gyyx_int3,gyyy_int3,gyyz_int3
  real*8, allocatable, dimension(:)         :: gyzx_int3,gyzy_int3,gyzz_int3
  real*8, allocatable, dimension(:)         :: gzzx_int3,gzzy_int3,gzzz_int3

  if(MOD(cctk_iteration,out_every)==0) then 

     BH1_posx = 5.d0
     BH1_posy = 0.d0
     BH2_posx = -5.d0
     BH2_posy = 0.d0

     ntot = output_Nlograd*output_Nphi*output_Ntheta


     allocate(pointcoords1(ntot,3))
     allocate(pointcoords2(ntot,3))
     allocate(pointcoords3(ntot,3))

     allocate(rho_b_int1(ntot))
     allocate(rho_b_int2(ntot))
     allocate(rho_b_int3(ntot))

     allocate(P_int1(ntot))
     allocate(P_int2(ntot))
     allocate(P_int3(ntot))

     allocate(vx_int1(ntot))
     allocate(vx_int2(ntot))
     allocate(vx_int3(ntot))

     allocate(vy_int1(ntot))
     allocate(vy_int2(ntot))
     allocate(vy_int3(ntot))

     allocate(vz_int1(ntot))
     allocate(vz_int2(ntot))
     allocate(vz_int3(ntot))

     allocate(u0_int1(ntot))
     allocate(u0_int2(ntot))
     allocate(u0_int3(ntot))

     allocate(lapse_rhs_int1(ntot))
     allocate(lapse_rhs_int2(ntot))
     allocate(lapse_rhs_int3(ntot))

     allocate(shiftx_rhs_int1(ntot))
     allocate(shiftx_rhs_int2(ntot))
     allocate(shiftx_rhs_int3(ntot))
     allocate(shifty_rhs_int1(ntot))
     allocate(shifty_rhs_int2(ntot))
     allocate(shifty_rhs_int3(ntot))
     allocate(shiftz_rhs_int1(ntot))
     allocate(shiftz_rhs_int2(ntot))
     allocate(shiftz_rhs_int3(ntot))

     allocate(phi_rhs_int1(ntot))
     allocate(phi_rhs_int2(ntot))
     allocate(phi_rhs_int3(ntot))

     allocate(gxx_rhs_int1(ntot))
     allocate(gxx_rhs_int2(ntot))
     allocate(gxx_rhs_int3(ntot))
     allocate(gxy_rhs_int1(ntot))
     allocate(gxy_rhs_int2(ntot))
     allocate(gxy_rhs_int3(ntot))
     allocate(gxz_rhs_int1(ntot))
     allocate(gxz_rhs_int2(ntot))
     allocate(gxz_rhs_int3(ntot))
     allocate(gyy_rhs_int1(ntot))
     allocate(gyy_rhs_int2(ntot))
     allocate(gyy_rhs_int3(ntot))
     allocate(gyz_rhs_int1(ntot))
     allocate(gyz_rhs_int2(ntot))
     allocate(gyz_rhs_int3(ntot))
     allocate(gzz_rhs_int1(ntot))
     allocate(gzz_rhs_int2(ntot))
     allocate(gzz_rhs_int3(ntot))

     allocate(gxxx_int1(ntot))
     allocate(gxxy_int1(ntot))
     allocate(gxxz_int1(ntot))
     allocate(gxyx_int1(ntot))
     allocate(gxyy_int1(ntot))
     allocate(gxyz_int1(ntot))
     allocate(gxzx_int1(ntot))
     allocate(gxzy_int1(ntot))
     allocate(gxzz_int1(ntot))
     allocate(gyyx_int1(ntot))
     allocate(gyyy_int1(ntot))
     allocate(gyyz_int1(ntot))
     allocate(gyzx_int1(ntot))
     allocate(gyzy_int1(ntot))
     allocate(gyzz_int1(ntot))
     allocate(gzzx_int1(ntot))
     allocate(gzzy_int1(ntot))
     allocate(gzzz_int1(ntot))

     allocate(gxxx_int2(ntot))
     allocate(gxxy_int2(ntot))
     allocate(gxxz_int2(ntot))
     allocate(gxyx_int2(ntot))
     allocate(gxyy_int2(ntot))
     allocate(gxyz_int2(ntot))
     allocate(gxzx_int2(ntot))
     allocate(gxzy_int2(ntot))
     allocate(gxzz_int2(ntot))
     allocate(gyyx_int2(ntot))
     allocate(gyyy_int2(ntot))
     allocate(gyyz_int2(ntot))
     allocate(gyzx_int2(ntot))
     allocate(gyzy_int2(ntot))
     allocate(gyzz_int2(ntot))
     allocate(gzzx_int2(ntot))
     allocate(gzzy_int2(ntot))
     allocate(gzzz_int2(ntot))

     allocate(gxxx_int3(ntot))
     allocate(gxxy_int3(ntot))
     allocate(gxxz_int3(ntot))
     allocate(gxyx_int3(ntot))
     allocate(gxyy_int3(ntot))
     allocate(gxyz_int3(ntot))
     allocate(gxzx_int3(ntot))
     allocate(gxzy_int3(ntot))
     allocate(gxzz_int3(ntot))
     allocate(gyyx_int3(ntot))
     allocate(gyyy_int3(ntot))
     allocate(gyyz_int3(ntot))
     allocate(gyzx_int3(ntot))
     allocate(gyzy_int3(ntot))
     allocate(gyzz_int3(ntot))
     allocate(gzzx_int3(ntot))
     allocate(gzzy_int3(ntot))
     allocate(gzzz_int3(ntot))

     allocate(lapse_int1(ntot))
     allocate(lapse_int2(ntot))
     allocate(lapse_int3(ntot))

     allocate(lapsex_int1(ntot))
     allocate(lapsex_int2(ntot))
     allocate(lapsex_int3(ntot))
     allocate(lapsey_int1(ntot))
     allocate(lapsey_int2(ntot))
     allocate(lapsey_int3(ntot))
     allocate(lapsez_int1(ntot))
     allocate(lapsez_int2(ntot))
     allocate(lapsez_int3(ntot)) 

     allocate(shiftx_int1(ntot))
     allocate(shiftx_int2(ntot))
     allocate(shiftx_int3(ntot))
     allocate(shifty_int1(ntot))
     allocate(shifty_int2(ntot))
     allocate(shifty_int3(ntot))
     allocate(shiftz_int1(ntot))
     allocate(shiftz_int2(ntot))
     allocate(shiftz_int3(ntot))

     allocate(shiftxx_int1(ntot))
     allocate(shiftxx_int2(ntot))
     allocate(shiftxx_int3(ntot))
     allocate(shiftyx_int1(ntot))
     allocate(shiftyx_int2(ntot))
     allocate(shiftyx_int3(ntot))
     allocate(shiftzx_int1(ntot))
     allocate(shiftzx_int2(ntot))
     allocate(shiftzx_int3(ntot))
     allocate(shiftxy_int1(ntot))
     allocate(shiftxy_int2(ntot))
     allocate(shiftxy_int3(ntot))
     allocate(shiftyy_int1(ntot))
     allocate(shiftyy_int2(ntot))
     allocate(shiftyy_int3(ntot))
     allocate(shiftzy_int1(ntot))
     allocate(shiftzy_int2(ntot))
     allocate(shiftzy_int3(ntot))
     allocate(shiftxz_int1(ntot))
     allocate(shiftxz_int2(ntot))
     allocate(shiftxz_int3(ntot))
     allocate(shiftyz_int1(ntot))
     allocate(shiftyz_int2(ntot))
     allocate(shiftyz_int3(ntot))
     allocate(shiftzz_int1(ntot))
     allocate(shiftzz_int2(ntot))
     allocate(shiftzz_int3(ntot))

     allocate(phi_int1(ntot))
     allocate(phi_int2(ntot))
     allocate(phi_int3(ntot))

     allocate(phix_int1(ntot))
     allocate(phix_int2(ntot))
     allocate(phix_int3(ntot))
     allocate(phiy_int1(ntot))
     allocate(phiy_int2(ntot))
     allocate(phiy_int3(ntot))
     allocate(phiz_int1(ntot))
     allocate(phiz_int2(ntot))
     allocate(phiz_int3(ntot))

     allocate(gxx_int1(ntot))
     allocate(gxx_int2(ntot))
     allocate(gxx_int3(ntot))
     allocate(gxy_int1(ntot))
     allocate(gxy_int2(ntot))
     allocate(gxy_int3(ntot))
     allocate(gxz_int1(ntot))
     allocate(gxz_int2(ntot))
     allocate(gxz_int3(ntot))
     allocate(gyy_int1(ntot))
     allocate(gyy_int2(ntot))
     allocate(gyy_int3(ntot))
     allocate(gyz_int1(ntot))
     allocate(gyz_int2(ntot))
     allocate(gyz_int3(ntot))
     allocate(gzz_int1(ntot))
     allocate(gzz_int2(ntot))
     allocate(gzz_int3(ntot))

     allocate(gupxx_int1(ntot))
     allocate(gupxx_int2(ntot))
     allocate(gupxx_int3(ntot))
     allocate(gupxy_int1(ntot))
     allocate(gupxy_int2(ntot))
     allocate(gupxy_int3(ntot))
     allocate(gupxz_int1(ntot))
     allocate(gupxz_int2(ntot))
     allocate(gupxz_int3(ntot))
     allocate(gupyy_int1(ntot))
     allocate(gupyy_int2(ntot))
     allocate(gupyy_int3(ntot))
     allocate(gupyz_int1(ntot))
     allocate(gupyz_int2(ntot))
     allocate(gupyz_int3(ntot))
     allocate(gupzz_int1(ntot))
     allocate(gupzz_int2(ntot))
     allocate(gupzz_int3(ntot))

     PI=acos(-1.0)
     thetamin=-PI/100000.
     thetamax=PI/2.0 + PI/100000.
     phimin=-PI/100000.
     phimax=2.0*PI+PI/100000.

     output_logradmin_rot_interp_1 = log(output_radmin_rot_interp_1)
     output_logradmin_rot_interp_2 = log(output_radmin_rot_interp_2)
     output_logradmin_rot_interp_3 = log(output_radmin_rot_interp_3)
     output_logradmax_rot_interp_1 = log(output_radmax_rot_interp_1)
     output_logradmax_rot_interp_2 = log(output_radmax_rot_interp_2)
     output_logradmax_rot_interp_3 = log(output_radmax_rot_interp_3)

     n=1
     do k=1,output_Nlograd
        do j=1,output_Nphi
           do i=1,output_Ntheta
              output_lograd1L=output_logradmin_rot_interp_1+(output_logradmax_rot_interp_1-output_logradmin_rot_interp_1)/(output_Nlograd-1)*(k-1)
              output_lograd2L=output_logradmin_rot_interp_2+(output_logradmax_rot_interp_2-output_logradmin_rot_interp_2)/(output_Nlograd-1)*(k-1)
              output_lograd3L=output_logradmin_rot_interp_3+(output_logradmax_rot_interp_3-output_logradmin_rot_interp_3)/(output_Nlograd-1)*(k-1)
              phiL=phimin+(phimax-phimin)/(output_Nphi-1)*(j-1)
              thetaL = thetamin+(thetamax-thetamin)/(output_Ntheta-1)*(i-1)
              pointcoords1(n,1) = exp(output_lograd1L)*sin(thetaL)*cos(phiL)+BH1_posx
              pointcoords1(n,2) = exp(output_lograd1L)*sin(thetaL)*sin(phiL)+BH1_posy
              pointcoords1(n,3) = exp(output_lograd1L)*cos(thetaL)
              pointcoords2(n,1) = exp(output_lograd2L)*sin(thetaL)*cos(phiL)+BH2_posx
              pointcoords2(n,2) = exp(output_lograd2L)*sin(thetaL)*sin(phiL)+BH2_posy
              pointcoords2(n,3) = exp(output_lograd2L)*cos(thetaL)
              pointcoords3(n,1) = exp(output_lograd3L)*sin(thetaL)*cos(phiL)
              pointcoords3(n,2) = exp(output_lograd3L)*sin(thetaL)*sin(phiL)
              pointcoords3(n,3) = exp(output_lograd3L)*cos(thetaL)    
              !if ((i.eq.24).and.(j.eq.24).and.(k.eq.11)) then
              !write(*,*) "X1: ",pointcoords1(n,1) 
              !write(*,*) "Y1: ",pointcoords1(n,2)
              !write(*,*) "Z1: ",pointcoords1(n,3)
                 !endif
              n = n + 1
           enddo
        enddo
     enddo



     ! Interpolate the grid functions
     call CCTK_VarIndex(vindex,"mhd_evolve::rho_b")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,rho_b_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,rho_b_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,rho_b_int3)

     call CCTK_VarIndex(vindex,"mhd_evolve::P")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,P_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,P_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,P_int3)

     call CCTK_VarIndex(vindex,"mhd_evolve::vx")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,vx_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,vx_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,vx_int3)
     call CCTK_VarIndex(vindex,"mhd_evolve::vy")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,vy_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,vy_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,vy_int3)
     call CCTK_VarIndex(vindex,"mhd_evolve::vz")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,vz_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,vz_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,vz_int3)
     call CCTK_VarIndex(vindex,"mhd_evolve::u0")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,u0_int1) 
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,u0_int2) 
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,u0_int3)

     call CCTK_VarIndex(vindex,"lapse::lapm1_rhs")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,lapse_rhs_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,lapse_rhs_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,lapse_rhs_int3)

     call CCTK_VarIndex(vindex,"shift::shiftx_rhs")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,shiftx_rhs_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,shiftx_rhs_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,shiftx_rhs_int3)
     call CCTK_VarIndex(vindex,"shift::shifty_rhs")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,shifty_rhs_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,shifty_rhs_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,shifty_rhs_int3)
     call CCTK_VarIndex(vindex,"shift::shiftz_rhs")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,shiftz_rhs_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,shiftz_rhs_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,shiftz_rhs_int3)

     call CCTK_VarIndex(vindex,"bssn::phi_rhs")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,phi_rhs_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,phi_rhs_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,phi_rhs_int3)

     call CCTK_VarIndex(vindex,"bssn::gxx_rhs")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gxx_rhs_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gxx_rhs_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gxx_rhs_int3)
     call CCTK_VarIndex(vindex,"bssn::gxy_rhs")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gxy_rhs_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gxy_rhs_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gxy_rhs_int3)
     call CCTK_VarIndex(vindex,"bssn::gxz_rhs")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gxz_rhs_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gxz_rhs_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gxz_rhs_int3)
     call CCTK_VarIndex(vindex,"bssn::gyy_rhs")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gyy_rhs_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gyy_rhs_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gyy_rhs_int3)
     call CCTK_VarIndex(vindex,"bssn::gyz_rhs")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gyz_rhs_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gyz_rhs_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gyz_rhs_int3)
     call CCTK_VarIndex(vindex,"bssn::gzz_rhs")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gzz_rhs_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gzz_rhs_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gzz_rhs_int3)

     call CCTK_VarIndex(vindex,"lapse::lapm1")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,lapse_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,lapse_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,lapse_int3)

     call CCTK_VarIndex(vindex,"lapse::lapsex")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,lapsex_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,lapsex_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,lapsex_int3)
     call CCTK_VarIndex(vindex,"lapse::lapsey")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,lapsey_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,lapsey_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,lapsey_int3)
     call CCTK_VarIndex(vindex,"lapse::lapsez")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,lapsez_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,lapsez_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,lapsez_int3)

     call CCTK_VarIndex(vindex,"shift::shiftx")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,shiftx_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,shiftx_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,shiftx_int3)
     call CCTK_VarIndex(vindex,"shift::shifty")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,shifty_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,shifty_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,shifty_int3)
     call CCTK_VarIndex(vindex,"shift::shiftz")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,shiftz_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,shiftz_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,shiftz_int3)

     call CCTK_VarIndex(vindex,"shift::shiftxx")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,shiftxx_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,shiftxx_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,shiftxx_int3)
     call CCTK_VarIndex(vindex,"shift::shiftyx")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,shiftyx_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,shiftyx_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,shiftyx_int3)
     call CCTK_VarIndex(vindex,"shift::shiftzx")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,shiftzx_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,shiftzx_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,shiftzx_int3)
     call CCTK_VarIndex(vindex,"shift::shiftxy")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,shiftxy_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,shiftxy_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,shiftxy_int3)
     call CCTK_VarIndex(vindex,"shift::shiftyy")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,shiftyy_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,shiftyy_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,shiftyy_int3)
     call CCTK_VarIndex(vindex,"shift::shiftzy")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,shiftzy_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,shiftzy_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,shiftzy_int3)
     call CCTK_VarIndex(vindex,"shift::shiftxz")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,shiftxz_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,shiftxz_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,shiftxz_int3)
     call CCTK_VarIndex(vindex,"shift::shiftyz")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,shiftyz_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,shiftyz_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,shiftyz_int3)
     call CCTK_VarIndex(vindex,"shift::shiftzz")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,shiftzz_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,shiftzz_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,shiftzz_int3)

     call CCTK_VarIndex(vindex,"bssn::phi")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,phi_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,phi_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,phi_int3)
     call CCTK_VarIndex(vindex,"bssn::phix")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,phix_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,phix_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,phix_int3)
     call CCTK_VarIndex(vindex,"bssn::phiy")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,phiy_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,phiy_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,phiy_int3)
     call CCTK_VarIndex(vindex,"bssn::phiz")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,phiz_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,phiz_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,phiz_int3)

     call CCTK_VarIndex(vindex,"bssn::gxx")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gxx_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gxx_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gxx_int3)
     call CCTK_VarIndex(vindex,"bssn::gxy")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gxy_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gxy_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gxy_int3)
     call CCTK_VarIndex(vindex,"bssn::gxz")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gxz_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gxz_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gxz_int3)
     call CCTK_VarIndex(vindex,"bssn::gyy")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gyy_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gyy_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gyy_int3)
     call CCTK_VarIndex(vindex,"bssn::gyz")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gyz_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gyz_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gyz_int3)
     call CCTK_VarIndex(vindex,"bssn::gzz")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gzz_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gzz_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gzz_int3)

     call CCTK_VarIndex(vindex,"bssn::gupxx")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gupxx_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gupxx_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gupxx_int3)
     call CCTK_VarIndex(vindex,"bssn::gupxy")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gupxy_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gupxy_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gupxy_int3)
     call CCTK_VarIndex(vindex,"bssn::gupxz")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gupxz_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gupxz_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gupxz_int3)
     call CCTK_VarIndex(vindex,"bssn::gupyy")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gupyy_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gupyy_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gupyy_int3)
     call CCTK_VarIndex(vindex,"bssn::gupyz")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gupyz_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gupyz_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gupyz_int3)
     call CCTK_VarIndex(vindex,"bssn::gupzz")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gupzz_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gupzz_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gupzz_int3)

     call CCTK_VarIndex(vindex,"bssn::gxxx")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gxxx_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gxxx_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gxxx_int3)
     call CCTK_VarIndex(vindex,"bssn::gxxy")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gxxy_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gxxy_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gxxy_int3)

     call CCTK_VarIndex(vindex,"bssn::gxxz")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gxxz_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gxxz_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gxxz_int3)
     call CCTK_VarIndex(vindex,"bssn::gxyx")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gxyx_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gxyx_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gxyx_int3)
     call CCTK_VarIndex(vindex,"bssn::gxyy")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gxyy_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gxyy_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gxyy_int3)
     call CCTK_VarIndex(vindex,"bssn::gxyz")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gxyz_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gxyz_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gxyz_int3)
     call CCTK_VarIndex(vindex,"bssn::gxzx")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gxzx_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gxzx_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gxzx_int3)
     call CCTK_VarIndex(vindex,"bssn::gxzy")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gxzy_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gxzy_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gxzy_int3)
     call CCTK_VarIndex(vindex,"bssn::gxzz")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gxzz_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gxzz_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gxzz_int3)
     call CCTK_VarIndex(vindex,"bssn::gyyx")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gyyx_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gyyx_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gyyx_int3)
     call CCTK_VarIndex(vindex,"bssn::gyyy")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gyyy_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gyyy_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gyyy_int3)
     call CCTK_VarIndex(vindex,"bssn::gyyz")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gyyz_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gyyz_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gyyz_int3)
     call CCTK_VarIndex(vindex,"bssn::gyzx")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gyzx_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gyzx_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gyzx_int3)
     call CCTK_VarIndex(vindex,"bssn::gyzy")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gyzy_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gyzy_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gyzy_int3)
     call CCTK_VarIndex(vindex,"bssn::gyzz")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gyzz_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gyzz_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gyzz_int3)
     call CCTK_VarIndex(vindex,"bssn::gzzx")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gzzx_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gzzx_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gzzx_int3)
     call CCTK_VarIndex(vindex,"bssn::gzzy")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gzzy_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gzzy_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gzzy_int3)
     call CCTK_VarIndex(vindex,"bssn::gzzz")
     call interp_driver_carp(cctkGH,ntot,pointcoords1,vindex,gzzz_int1)
     call interp_driver_carp(cctkGH,ntot,pointcoords2,vindex,gzzz_int2)
     call interp_driver_carp(cctkGH,ntot,pointcoords3,vindex,gzzz_int3)


     call bbh_bondi_IO(cctkGH,cctk_time,cctk_iteration,1,&
          output_radmin_rot_interp_1,output_radmax_rot_interp_1,&
          output_Nlograd,output_Nphi,output_Ntheta,&
          rho_b_int1,P_int1,vx_int1,vy_int1,vz_int1,u0_int1,&
          lapse_int1,shiftx_int1,shifty_int1,shiftz_int1,&
          lapsex_int1,shiftxx_int1,shiftyx_int1,shiftzx_int1,&
          lapsey_int1,shiftxy_int1,shiftyy_int1,shiftzy_int1,&
          lapsez_int1,shiftxz_int1,shiftyz_int1,shiftzz_int1,&
          lapse_rhs_int1,shiftx_rhs_int1,shifty_rhs_int1,shiftz_rhs_int1,&
          phi_int1,gxx_int1,gxy_int1,gxz_int1,&
          gyy_int1,gyz_int1, gzz_int1,&
          gupxx_int1,gupxy_int1,gupxz_int1,&
          gupyy_int1,gupyz_int1, gupzz_int1,&
          phix_int1,gxxx_int1,gxyx_int1,gxzx_int1,&
          gyyx_int1,gyzx_int1, gzzx_int1,&
          phiy_int1,gxxy_int1,gxyy_int1,gxzy_int1,&
          gyyy_int1,gyzy_int1, gzzy_int1,&
          phiz_int1,gxxz_int1,gxyz_int1,gxzz_int1,&
          gyyz_int1,gyzz_int1, gzzz_int1,&
          phi_rhs_int1,gxx_rhs_int1,gxy_rhs_int1,gxz_rhs_int1,&
          gyy_rhs_int1,gyz_rhs_int1,gzz_rhs_int1,&
          BH1_posx,BH1_posy)

     call bbh_bondi_IO(cctkGH,cctk_time,cctk_iteration,2,&
          output_radmin_rot_interp_2,output_radmax_rot_interp_2,&
          output_Nlograd,output_Nphi,output_Ntheta,&
          rho_b_int2,P_int2,vx_int2,vy_int2,vz_int2,u0_int2,&
          lapse_int2,shiftx_int2,shifty_int2,shiftz_int2,&
          lapsex_int2,shiftxx_int2,shiftyx_int2,shiftzx_int2,&
          lapsey_int2,shiftxy_int2,shiftyy_int2,shiftzy_int2,&
          lapsez_int2,shiftxz_int2,shiftyz_int2,shiftzz_int2,&
          lapse_rhs_int2,shiftx_rhs_int2,shifty_rhs_int2,shiftz_rhs_int2,&
          phi_int2,gxx_int2,gxy_int2,gxz_int2,&
          gyy_int2,gyz_int2, gzz_int2,&
          gupxx_int2,gupxy_int2,gupxz_int2,&
          gupyy_int2,gupyz_int2, gupzz_int2,&
          phix_int2,gxxx_int2,gxyx_int2,gxzx_int2,&
          gyyx_int2,gyzx_int2, gzzx_int2,&
          phiy_int2,gxxy_int2,gxyy_int2,gxzy_int2,&
          gyyy_int2,gyzy_int2, gzzy_int2,&
          phiz_int2,gxxz_int2,gxyz_int2,gxzz_int2,&
          gyyz_int2,gyzz_int2, gzzz_int2,&
          phi_rhs_int2,gxx_rhs_int2,gxy_rhs_int2,gxz_rhs_int2,&
          gyy_rhs_int2,gyz_rhs_int2,gzz_rhs_int2,&
          BH2_posx,BH2_posy)

     call bbh_bondi_IO(cctkGH,cctk_time,cctk_iteration,3,&
          output_radmin_rot_interp_3,output_radmax_rot_interp_3,&
          output_Nlograd,output_Nphi,output_Ntheta,&
          rho_b_int3,P_int3,vx_int3,vy_int3,vz_int3,u0_int3,&
          lapse_int3,shiftx_int3,shifty_int3,shiftz_int3,&
          lapsex_int3,shiftxx_int3,shiftyx_int3,shiftzx_int3,&
          lapsey_int3,shiftxy_int3,shiftyy_int3,shiftzy_int3,&
          lapsez_int3,shiftxz_int3,shiftyz_int3,shiftzz_int3,&
          lapse_rhs_int3,shiftx_rhs_int3,shifty_rhs_int3,shiftz_rhs_int3,&
          phi_int3,gxx_int3,gxy_int3,gxz_int3,&
          gyy_int3,gyz_int3, gzz_int3,&
          gupxx_int3,gupxy_int3,gupxz_int3,&
          gupyy_int3,gupyz_int3, gupzz_int3,&
          phix_int3,gxxx_int3,gxyx_int3,gxzx_int3,&
          gyyx_int3,gyzx_int3, gzzx_int3,&
          phiy_int3,gxxy_int3,gxyy_int3,gxzy_int3,&
          gyyy_int3,gyzy_int3, gzzy_int3,&
          phiz_int3,gxxz_int3,gxyz_int3,gxzz_int3,&
          gyyz_int3,gyzz_int3, gzzz_int3,&
          phi_rhs_int3,gxx_rhs_int3,gxy_rhs_int3,gxz_rhs_int3,&
          gyy_rhs_int3,gyz_rhs_int3,gzz_rhs_int3,&
          0.d0,0.d0)
    


     !  write(*,*) "just called output lowlevel"
     deallocate(pointcoords1,pointcoords2,pointcoords3)
     deallocate(rho_b_int1,rho_b_int2,rho_b_int3)
     deallocate(P_int1,P_int2,P_int3)
     deallocate(vx_int1,vx_int2,vx_int3)
     deallocate(vy_int1,vy_int2,vy_int3)
     deallocate(vz_int1,vz_int2,vz_int3)
     deallocate(u0_int1,u0_int2,u0_int3)
     deallocate(lapse_rhs_int1,lapse_rhs_int2,lapse_rhs_int3)
     deallocate(shiftx_rhs_int1,shiftx_rhs_int2,shiftx_rhs_int3)
     deallocate(shifty_rhs_int1,shifty_rhs_int2,shifty_rhs_int3)
     deallocate(shiftz_rhs_int1,shiftz_rhs_int2,shiftz_rhs_int3)
     deallocate(phi_rhs_int1,phi_rhs_int2,phi_rhs_int3)
     deallocate(gxx_rhs_int1,gxx_rhs_int2,gxx_rhs_int3)
     deallocate(gxy_rhs_int1,gxy_rhs_int2,gxy_rhs_int3)
     deallocate(gxz_rhs_int1,gxz_rhs_int2,gxz_rhs_int3)
     deallocate(gyy_rhs_int1,gyy_rhs_int2,gyy_rhs_int3)
     deallocate(gyz_rhs_int1,gyz_rhs_int2,gyz_rhs_int3)
     deallocate(gzz_rhs_int1,gzz_rhs_int2,gzz_rhs_int3)
     deallocate(lapse_int1,lapse_int2,lapse_int3)
     deallocate(lapsex_int1,lapsex_int2,lapsex_int3)
     deallocate(lapsey_int1,lapsey_int2,lapsey_int3)
     deallocate(lapsez_int1,lapsez_int2,lapsez_int3)
     deallocate(shiftx_int1,shiftx_int2,shiftx_int3)
     deallocate(shifty_int1,shifty_int2,shifty_int3)
     deallocate(shiftz_int1,shiftz_int2,shiftz_int3)
     deallocate(shiftxx_int1,shiftxx_int2,shiftxx_int3)
     deallocate(shiftyx_int1,shiftyx_int2,shiftyx_int3)
     deallocate(shiftzx_int1,shiftzx_int2,shiftzx_int3)
     deallocate(shiftxy_int1,shiftxy_int2,shiftxy_int3)
     deallocate(shiftyy_int1,shiftyy_int2,shiftyy_int3)
     deallocate(shiftzy_int1,shiftzy_int2,shiftzy_int3)
     deallocate(shiftxz_int1,shiftxz_int2,shiftxz_int3)
     deallocate(shiftyz_int1,shiftyz_int2,shiftyz_int3)
     deallocate(shiftzz_int1,shiftzz_int2,shiftzz_int3)
     deallocate(phi_int1,phi_int2,phi_int3)
     deallocate(phix_int1,phix_int2,phix_int3)
     deallocate(phiy_int1,phiy_int2,phiy_int3)
     deallocate(phiz_int1,phiz_int2,phiz_int3)
     deallocate(gxx_int1,gxx_int2,gxx_int3)
     deallocate(gxy_int1,gxy_int2,gxy_int3)
     deallocate(gxz_int1,gxz_int2,gxz_int3)
     deallocate(gyy_int1,gyy_int2,gyy_int3)
     deallocate(gyz_int1,gyz_int2,gyz_int3)
     deallocate(gzz_int1,gzz_int2,gzz_int3)
     deallocate(gupxx_int1,gupxx_int2,gupxx_int3)
     deallocate(gupxy_int1,gupxy_int2,gupxy_int3)
     deallocate(gupxz_int1,gupxz_int2,gupxz_int3)
     deallocate(gupyy_int1,gupyy_int2,gupyy_int3)
     deallocate(gupyz_int1,gupyz_int2,gupyz_int3)
     deallocate(gupzz_int1,gupzz_int2,gupzz_int3)
     deallocate(gxxx_int1,gxxy_int1,gxxz_int1)
     deallocate(gxyx_int1,gxyy_int1,gxyz_int1)
     deallocate(gxzx_int1,gxzy_int1,gxzz_int1)
     deallocate(gyyx_int1,gyyy_int1,gyyz_int1)
     deallocate(gyzx_int1,gyzy_int1,gyzz_int1)
     deallocate(gzzx_int1,gzzy_int1,gzzz_int1)
     deallocate(gxxx_int2,gxxy_int2,gxxz_int2)
     deallocate(gxyx_int2,gxyy_int2,gxyz_int2)
     deallocate(gxzx_int2,gxzy_int2,gxzz_int2)
     deallocate(gyyx_int2,gyyy_int2,gyyz_int2)
     deallocate(gyzx_int2,gyzy_int2,gyzz_int2)
     deallocate(gzzx_int2,gzzy_int2,gzzz_int2)
     deallocate(gxxx_int3,gxxy_int3,gxxz_int3)
     deallocate(gxyx_int3,gxyy_int3,gxyz_int3)
     deallocate(gxzx_int3,gxzy_int3,gxzz_int3)
     deallocate(gyyx_int3,gyyy_int3,gyyz_int3)
     deallocate(gyzx_int3,gyzy_int3,gyzz_int3)
     deallocate(gzzx_int3,gzzy_int3,gzzz_int3)
     !stop
  endif
end subroutine bbh_bondi_output
