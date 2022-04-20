
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bhns_regridder_reader_driver_part2_primitives_etc(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  ! Set primitives & metric source terms.

  if(primitives_solver==11) then
     !Here we use the HARM primitives solver, with a new prescription that minimizes changes to
     !  conservatives without applying the "Font" fix.
     !We intend to make the below function truly generic, but currently it only supports Gamma-law EOS.
     !We hope the below function will eventually be used in place of other primitives solver options,
     !  since it can be easily extended to use your desired technique for solving primitives.
     call primitives_generic(cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z, &
          phi,gxx,gxy,gxz,gyy,gyz,gzz, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
          lapm1,shiftx,shifty,shiftz, &
          Bx,By,Bz, &
          mhd_st_x,mhd_st_y,mhd_st_z,tau,rho_star, &
          vx,vy,vz,P,rho_b,h,u0, &
          rho_b_atm,tau_atm, &
          neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
          temp1,temp2,temp3,temp4,temp5, &
          primitives_debug,Psi6threshold)

     !The goal of the below function is to just update the metric source terms.  We should
     !  really get rid of the old h,w,st_i,Ei gridfunctions...
     !Currently this only works for Gamma-law EOS's but could be extended for arbitary EOSs
     !  with some work.
     call metric_source_terms_and_misc_vars(cctkGH,cctk_lsh, &
          rho,Sx,Sy,Sz, &
          Sxx,Sxy,Sxz,Syy,Syz,Szz, &
          h,w,st_x,st_y,st_z, &
          Ex,Ey,Ez, &
          sbt,sbx,sby,sbz, &
          lapm1,shiftx,shifty,shiftz, &
          phi,gxx,gxy,gxz,gyy,gyz,gzz, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
          P,rho_b,u0,vx,vy,vz, &
          Bx,By,Bz, &
          rho_star,mhd_st_x,mhd_st_y,mhd_st_z)
  else
     write(*,*) "ERROR. Only primitives_solver==1 is supported in the regridder!"
     stop
  end if

  
  !Following lines are (empirically speaking) ESSENTIAL for excision runs:
  call Derivs(cctk_lsh,X,Y,Z,dX,dY,dZ,phi,phix,phiy,phiz,Symmetry) 
  call Derivs(cctk_lsh,X,Y,Z,dX,dY,dZ,lapm1,lapsex,lapsey,lapsez,Symmetry) 
  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,phi,phix,phiy,phiz)
  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,lapm1,lapsex,lapsey,lapsez)

  ! Initialize t_last, M0dot_last and int_M0dot: variables necessary to compute the
  !  time integrated M0 flux across BH
  t_last = 0.d0
  M0dot_last = 0.d0
  int_M0dot = 0.d0

end subroutine bhns_regridder_reader_driver_part2_primitives_etc
