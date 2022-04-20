#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bondi_initial_data(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext
  real*8 			           :: dT,dX,dY,dZ
  real*8 :: xmin,ymin,zmin,rho_fail_max_step,m_fail_step,reduction_value
  real*8 :: rho_max,xmax2,ymax2,zmax2
  integer :: i,j,l
  integer :: handle,dummy,proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax,glob_imax,glob_jmax,glob_kmax
  integer :: index
  integer :: ierr
  real*8 :: ZERO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  parameter(ZERO = 0.D0)

  ext(1) = cctk_lsh(1)
  ext(2) = cctk_lsh(2)
  ext(3) = cctk_lsh(3)  
  dT = CCTK_DELTA_TIME
  dX = X(2,1,1) - X(1,1,1)
  dY = Y(1,2,1) - Y(1,1,1)
  dZ = Z(1,1,2) - Z(1,1,1)

  mhd_r_crit = r_crit
  mhd_mdot = mdot
  mhd_mbh = mbh
  mhd_bigP = bigP

  if(bondi_enable.eq.1) then
     call moving_puncture_bondi(ext,X,Y,Z,PhysicalRadius, &
          Mbh,rho_star,tau,st_x,st_y,st_z,w,h,u0,rho_b,P,vx,vy,vz, &
          gamma_th,Mdot,r_crit, &
          K_poly,phi,lapm1, &
          xbh1,zbh1,bigP,Symmetry)

     n_poly=1.0d0/(gamma_th-1.0d0)
     do i=1,2
        k_tab(i)=k_poly
        gamma_tab(i)=gamma_th
     enddo
     rho_tab(1)=1.0d0
     P_tab(1)=K_poly
     eps_tab(1)=n_poly*P_tab(1)/rho_tab(1)
     
     call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
     call CCTK_VarIndex(index,"mhd_evolve::rho_b")
     call CCTK_ReductionHandle(handle,"maximum")
     if (handle .gt. 0) then
        call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,reduction_value,1,index)
        if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH).eq.0) then
           print *,"Maximum value of rho_b is ",reduction_value
        end if
     else
        call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
     end if
     rho_max = reduction_value
     
     rho_b_atm = rho_fact*rho_max

     write(6,*)'Maxrho:',rho_max,' x:',rho_fact,' rhoatm:',rho_b_atm

     if (Matter_BC==7) call Setup_Rhovec(CCTK_PASS_FTOF)

     rho = h*w*exp(-6.0*phi)-P
     Sx = st_x*exp(-6.0*phi)
     Sy = st_y*exp(-6.0*phi)
     Sz = st_z*exp(-6.0*phi)
     Sxx = st_x*st_x/w/h*exp(-6.0*phi) + P*gxx*exp(4.0*phi)
     Sxy = st_x*st_y/w/h*exp(-6.0*phi) + P*gxy*exp(4.0*phi)
     Sxz = st_x*st_z/w/h*exp(-6.0*phi) + P*gxz*exp(4.0*phi)
     Syy = st_y*st_y/w/h*exp(-6.0*phi) + P*gyy*exp(4.0*phi)
     Syz = st_y*st_z/w/h*exp(-6.0*phi) + P*gyz*exp(4.0*phi)
     Szz = st_z*st_z/w/h*exp(-6.0*phi) + P*gzz*exp(4.0*phi)

!     call primitive_vars_hybrid2(ext,cctk_nghostzones,X,Y,Z, &
!          rho_star,tau,st_x,st_y,st_z, &
!          mhd_st_x,mhd_st_y,mhd_st_z,neos, &
!          rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, &
!          w,w_p,rho_b,rho,P,h,Sx,Sy,Sz, &
!          Sxx,Sxy,Sxz,Syy,Syz,Szz, &
!          phi,lapm1,shiftx,shifty,shiftz, &
!          gxx,gxy,gxz,gyy,gyz,gzz, &
!          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
!          h_p,u0,rho_b_max,rho_b_atm, &
!          rho_fail_max_step,M_fail_step,rhos_max, &
!          Bx,By,Bz,Ex,Ey,Ez, &
!          vx,vy,vz, &
!          sbt,sbx,sby,sbz, &
!          proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
!          glob_imax,glob_jmax,glob_kmax,Symmetry,pfloor,excision_enable, &
!          excision_zone_gf, tau_stildefix_enable,0)
  else
  ! Set matter variables to zero!
     rho = ZERO
     S = ZERO
     Sx = ZERO
     Sy = ZERO
     Sz = ZERO
     Sxx = ZERO
     Sxy = ZERO
     Sxz = ZERO
     Syy = ZERO
     Syz = ZERO
     Szz = ZERO
  endif

  mhd_st_x=st_x
  mhd_st_y=st_y
  mhd_st_z=st_z
  sbt=ZERO
  sbx=ZERO
  sby=ZERO
  sbz=ZERO
  Bx=ZERO
  By=ZERO
  Bz=ZERO
  Ex=ZERO
  Ey=ZERO
  Ez=ZERO
  
  glob_imax = ubound(rho_star,1)
  glob_jmax = ubound(rho_star,2)
  glob_kmax = ubound(rho_star,3)
  
  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"grid::X")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xmax2,1,index)
  call CCTK_VarIndex(index,"grid::Y")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,ymax2,1,index)
  call CCTK_VarIndex(index,"grid::Z")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,zmax2,1,index)
  call CCTK_ReductionHandle(handle,"minimum")
  call CCTK_VarIndex(index,"grid::X")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xmin,1,index)
  call CCTK_VarIndex(index,"grid::Y")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,ymin,1,index)
  call CCTK_VarIndex(index,"grid::Z")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,zmin,1,index)
  
  if(X(cctk_lsh(1),1,1) .eq. xmax2) then
     proc_imax = glob_imax
  else
     proc_imax = -1
  end if
  if(Y(1,cctk_lsh(2),1) .eq. ymax2) then
     proc_jmax = glob_jmax
  else
     proc_jmax = -1
  end if
  if(Z(1,1,cctk_lsh(3)) .eq. zmax2) then
     proc_kmax = glob_kmax
  else
     proc_kmax = -1
  end if
  ! Possible bug if xmin,ymin,zmin not set (properly):
  !  if(abs(X(1,1,1)) .eq. dX*0.5D0) then
  ! Possible bug if xmin,ymin,zmin not set (properly):
  if(X(1,1,1) .eq. xmin) then
     proc_imin = 0
  else
     proc_imin = -100
  end if
  if(Y(1,1,1) .eq. ymin) then
     proc_jmin = 0
  else
     proc_jmin = -100
  end if
  if(Z(1,1,1) .eq. zmin) then
     proc_kmin = 0
  else
     proc_kmin = -100
  end if
  
  dtrhob=rho_b
  
  write(6,*)'rhos3:',rho_star(2,2,2),rho_star(3,2,2)
  write(6,*)'rhob3:',rho_b(2,2,2),rho_b(3,2,2)
  write(6,*)'vx3:',vx(2,2,2),vx(3,2,2)
  write(6,*)'rho3:',rho(2,2,2),rho(3,2,2)
  write(6,*)'sx3:',sx(2,2,2),sx(3,2,2)
  write(6,*)'sy3:',sy(2,2,2),sy(3,2,2)
  write(6,*)'sz3:',sz(2,2,2),sz(3,2,2)
  write(6,*)'sxx3:',sxx(2,2,2),sxx(3,2,2)
  write(6,*)'sxy3:',sxy(2,2,2),sxy(3,2,2)
  write(6,*)'sxz3:',sxz(2,2,2),sxz(3,2,2)
  write(6,*)'syy3:',syy(2,2,2),syy(3,2,2)
  write(6,*)'syz3:',syz(2,2,2),syz(3,2,2)
  write(6,*)'szz3:',szz(2,2,2),szz(3,2,2)
  write(6,*)'h3:',h(2,2,2),h(3,2,2)
  write(6,*)'u03:',u0(2,2,2),u0(3,2,2)
  
  
end subroutine bondi_initial_data
