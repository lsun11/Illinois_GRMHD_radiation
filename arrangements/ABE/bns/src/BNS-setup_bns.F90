#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine BNS_finalize_initialdata(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext
  real*8 			           :: dT,dX,dY,dZ,P_max,q_max,rho_max,rhos_max,rho_fail_max_step,M_fail_step,q_atmos
  real*8 			           :: detmin_l,detmax_l,delta_bar,xmin,ymin,zmin,xmax,ymax,zmax
  integer :: n1,n2,n3,mf
  integer :: handle,dummy,proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax,glob_imax,glob_jmax,glob_kmax
  integer :: i,j,k,index
  integer :: ierr,ONE,ZERO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
!  real*8, dimension(Nxh*2) :: Xh
!  real*8, dimension(Nyh*2) :: Yh
!  real*8, dimension(Nzh*2) :: Zh
!  double * Xh; Xh = new double[Nxh*2];
!  double * Yh; Yh = new double[Nyh*2];
!  double * Zh; Zh = new double[Nzh*2];

  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  CCTK_REAL reduction_value

  parameter(ONE = 1.D0, ZERO = 0.D0)


  ext(1) = cctk_lsh(1)
  ext(2) = cctk_lsh(2)
  ext(3) = cctk_lsh(3)  
  dT = CCTK_DELTA_TIME
  dX = X(2,1,1) - X(1,1,1)
  dY = Y(1,2,1) - Y(1,1,1)
  dZ = Z(1,1,2) - Z(1,1,1)

  glob_imax = ubound(rho_star,1)
  glob_jmax = ubound(rho_star,2)
  glob_kmax = ubound(rho_star,3)

  xpc1_fish = vxl(2,2,2)
  xpc2_fish = vyl(2,2,2)

! *** TEST ***
  write(*,*) 'xpc1_fish, xpc2_fish = ',xpc1_fish,xpc2_fish
! ************

  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"grid::X")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xmax,1,index)
  call CCTK_VarIndex(index,"grid::Y")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,ymax,1,index)
  call CCTK_VarIndex(index,"grid::Z")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,zmax,1,index)
  call CCTK_ReductionHandle(handle,"minimum")
  call CCTK_VarIndex(index,"grid::X")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xmin,1,index)
  call CCTK_VarIndex(index,"grid::Y")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,ymin,1,index)
  call CCTK_VarIndex(index,"grid::Z")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,zmin,1,index)

  if(X(cctk_lsh(1),1,1) .eq. xmax) then
     proc_imax = glob_imax
  else 
     proc_imax = -1
  end if
  if(Y(1,cctk_lsh(2),1) .eq. ymax) then
     proc_jmax = glob_jmax
  else 
     proc_jmax = -1
  end if
  if(Z(1,1,cctk_lsh(3)) .eq. zmax) then
     proc_kmax = glob_kmax
  else 
     proc_kmax = -1
  end if
  ! Possible bug if xmin,ymin,zmin not set (properly): 
  if(X(1,1,1) .eq. CCTK_ORIGIN_SPACE(1)) then
     proc_imin = 0
  else 
     proc_imin = -100
  end if
  if(Y(1,1,1) .eq. CCTK_ORIGIN_SPACE(2)) then
     proc_jmin = 0
  else 
     proc_jmin = -100
  end if
  if(Z(1,1,1) .eq. CCTK_ORIGIN_SPACE(3)) then
     proc_kmin = 0
  else 
     proc_kmin = -100
  end if

!First fix the symmetries on all the initial data that has been set.
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vars')
  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CartSymGN(dummy,cctkGH,'lapse::lapse_evolve')
  call CartSymGN(dummy,cctkGH,'shift::shift_evolve')

  !BAH: above lines have bugs in 'em...  Following works on one processor only.
  !  proc_imax = glob_imax
  !  proc_jmax = glob_jmax
  !  proc_kmax = glob_kmax
  !  proc_imin = 0.D0
  !  proc_jmin = 0.D0
  !  proc_kmin = 0.D0

  !Setting trivial functions:
  gupxx = ONE
  gupyy = ONE
  gupzz = ONE
  gupxy = ZERO
  gupxz = ZERO
  gupyz = ZERO
  gxx = ONE
  gxy = ZERO
  gxz = ZERO
  gyy = ONE
  gyz = ZERO
  gzz = ONE
  trK = ZERO
  Gammax = ZERO
  Gammay = ZERO
  Gammaz = ZERO
  if(fisheye_enable==1) then
     call trans_phys_fish_tensor_flat(ext,X,Y,Z, &
          PhysicalRadius,RadiusDerivative, &
          gxx,gxy,gxz,gyy,gyz,gzz,Symmetry)
     call trans_phys_fish_tensor_inv(ext,X,Y,Z, &
          PhysicalRadius,RadiusDerivative, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,Symmetry)
     call trans_phys_fish_gamt_flat(ext,X,Y,Z, &
          PhysicalRadius,RadiusDerivative,RadiusDerivative2, &
          Gammax,Gammay,Gammaz,Symmetry)
  endif
  MRsx = ZERO
  MRsy = ZERO
  MRsz = ZERO


  !write(*,*) "PROCS: ",cctk_gsh(1),cctk_lsh(1),proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax,glob_imax,glob_jmax,glob_kmax

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CartSymGN(dummy,cctkGH,'lapse::lapse_evolve')
  call CartSymGN(dummy,cctkGH,'shift::shift_evolve')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_rho_br_rho_bl')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')

!  write(*,*) "hewwo!",Nxh,Nyh,Nzh
!
!  do i = 1, Nxh
!     Xh(i+Nxh) = xminh + (i+0.5) * dxh;
!     Xh(Nxh-i-1) = xminh - (i+0.5) * dxh;
!  end do
!  for(i=0; i<Nyh; i++) {
!  do i = 1, Nyh
!     Yh(i+Nyh) = yminh + (i+0.5) * dyh;
!     Yh(Nyh-i-1) = yminh - (i+0.5) * dyh;
!  end do
!  do i = 1, Nzh
     !  for(i=0; i<Nzh; i++) {
!     Zh(i+Nzh) = zminh + (i+0.5) * dzh;
!     Zh(Nzh-i-1) = zminh - (i+0.5) * dzh;
!  end do

  !=========================================================
  ! Copy reference arrays to initial conditions   
  !========================================================= 
!  call translate_irr(ext,Nxh,Nyh,Nzh,nx_adm,ny_adm,n_tot_adm, &
!       Xh,Yh,Zh,X,Y,Z, &
!       shiftx_i,shifty_i,shiftz_i, &
!       shiftx,shifty,shiftz, &
!       q_i,q,lapse_i,lapm1, &
!       phi_i,phi, &
!       ux_i,uy_i,uz_i,ux,uy,uz)
!       lbx,lby,lbz)
!  write(*,*) "meowblauh0: ",q(2,2,2)
  !set symmetry on q:
  !rho_br(:,:,1) = rho_br(:,:,2)
  !set symmetry on u's (doesn't have to be perfect):
  !Pr(:,:,1) = Pr(:,:,2)
  !vxr(:,:,1) = vxr(:,:,2)
  !vyr(:,:,1) = vyr(:,:,2)
  !vzr(:,:,1) = -vzr(:,:,2)

  call CCTK_VarIndex(index,"mhd_evolve::rho_br")
  call CCTK_ReductionHandle(handle,"maximum")
  if (handle .gt. 0) then
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,reduction_value,1,index)
     if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH).eq.0) then
        print *,"Maximum value of q is ",reduction_value
     end if
  else
     call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
  end if
  q_max = reduction_value
  
  q_atmos = rho_fact*q_max
!may be incorrect:
  rho_b_atm = rho_fact*q_max

 write(*,*) "rho_b_atm = ",rho_b_atm

  !  TEMPORARY: Polytropic Index hardwired here!
  n_poly = 1.D0
  call compute_ibwh(ext,Omega_Frame,q_atmos,n_poly, &
       X,Y,Z, &
       rho_br,lapm1,phi, &
       shiftx,shifty,shiftz, &
       vxr,vyr,vzr,Pr, &
       rho,S,Sx,Sy,Sz, &
       Sxx,Sxy,Sxz,Syy,Syz,Szz, &
       rho_star,tau,st_x,st_y,st_z, &
       P,w,vx,vy,vz,rho_b, &
       u0,h)
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vars')
  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')

! Add B-field to the system
  Bx = ZERO
  By = ZERO
  Bz = ZERO
  sbt = ZERO
  sbx = ZERO
  sby = ZERO
  sbz = ZERO
  Ex = ZERO
  Ey = ZERO
  Ez = ZERO
  mhd_st_x = st_x
  mhd_st_y = st_y
  mhd_st_z = st_z
  call BNS_setup_emfield(CCTK_PASS_FTOF)

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

  call CCTK_VarIndex(index,"mhd_evolve::rho_star")
  call CCTK_ReductionHandle(handle,"maximum")
  if (handle .gt. 0) then
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,reduction_value,1,index)
     if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH).eq.0) then
        print *,"Maximum value of rho_star is ",reduction_value
     end if
  else
     call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
  end if
  rhos_max = reduction_value 

  call CCTK_VarIndex(index,"mhd_evolve::P")
  call CCTK_ReductionHandle(handle,"maximum")
  if (handle .gt. 0) then
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,reduction_value,1,index)
  else
     call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
  end if
  pfloor = P_fact*P_max

  call primitive_vars_hybrid2(ext,cctk_nghostzones,X,Y,Z, &
       rho_star,tau,st_x,st_y,st_z, &
       mhd_st_x,mhd_st_y,mhd_st_z,neos, &
       rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, &
       w,w,rho_b,rho,P,h,Sx,Sy,Sz, &
       Sxx,Sxy,Sxz,Syy,Syz,Szz, &
       phi,lapm1,shiftx,shifty,shiftz, &
       gxx,gxy,gxz,gyy,gyz,gzz, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
       h,u0,rho_max,rho_b_atm, &
       rho_fail_max_step,M_fail_step,rhos_max, &
       Bx,By,Bz,Ex,Ey,Ez, &
       vx,vy,vz, &
       sbt,sbx,sby,sbz, temp4, &
       proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
       glob_imax,glob_jmax,glob_kmax,Symmetry,pfloor,excision_enable, &
       excision_zone_gf,tau_stildefix_enable,0)

  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vars')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vs')
  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CartSymGN(dummy,cctkGH,'lapse::lapse_evolve')
  call CartSymGN(dummy,cctkGH,'shift::shift_evolve')

  if (kt_id==0) then 
     call kset_bns(ext,X,Y,Z, Axx,Axy,Axz,Ayy,Ayz,Azz,  & 
       shiftx,shifty,shiftz,phi,lapm1, Symmetry)

     call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')

     n1 = 0
     n2 = 0
     n3 = 0
     mf = 1
     call update_falloff(ext,X,Y,Z,Axx,n1,n2,n3,mf,Symmetry)
     call update_falloff(ext,X,Y,Z,Axy,n1,n2,n3,mf,Symmetry)
     call update_falloff(ext,X,Y,Z,Axz,n1,n2,n3,mf,Symmetry)
     call update_falloff(ext,X,Y,Z,Ayy,n1,n2,n3,mf,Symmetry)
     call update_falloff(ext,X,Y,Z,Ayz,n1,n2,n3,mf,Symmetry)
     call update_falloff(ext,X,Y,Z,Azz,n1,n2,n3,mf,Symmetry)
  end if

  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')
  call CCTK_SyncGroup(dummy,cctkGH,'lapse::lapse_evolve')
  call CCTK_SyncGroup(dummy,cctkGH,'shift::shift_evolve')
  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CartSymGN(dummy,cctkGH,'lapse::lapse_evolve')
  call CartSymGN(dummy,cctkGH,'shift::shift_evolve')

 !Compute K_ij's for apparent horizon finder:
 Psi = exp(phi)
 Kxx = Psi*Psi*Psi*Psi * (Axx + (1.D0/3.D0) * gxx * trK)
 Kxy = Psi*Psi*Psi*Psi * (Axy + (1.D0/3.D0) * gxy * trK)
 Kxz = Psi*Psi*Psi*Psi * (Axz + (1.D0/3.D0) * gxz * trK)
 Kyy = Psi*Psi*Psi*Psi * (Ayy + (1.D0/3.D0) * gyy * trK)
 Kyz = Psi*Psi*Psi*Psi * (Ayz + (1.D0/3.D0) * gyz * trK)
 Kzz = Psi*Psi*Psi*Psi * (Azz + (1.D0/3.D0) * gzz * trK)

 lapset = 0.D0
 shiftxt = 0.d0
 shiftyt = 0.d0
 shiftzt = 0.d0

end subroutine BNS_finalize_initialdata
