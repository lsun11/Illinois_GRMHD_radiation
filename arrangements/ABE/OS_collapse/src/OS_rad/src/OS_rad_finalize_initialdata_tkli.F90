#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine OS_rad_finalize_initialdata(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  integer, dimension(3)                    :: ext
  real*8                                   :: dT,dX,dY,dZ,P_max,rho_max,rhos_max,rho_fail_max_step,M_fail_step
  real*8                                   :: detmin_l,detmax_l,delta_bar,xmin,ymin,zmin,xmax,ymax,zmax               
  integer :: n1,n2,n3,mf
  integer :: handle,dummy,proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax,glob_imax,glob_jmax,glob_kmax
  integer :: index
  integer :: ierr,ONE,ZERO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  real*8 :: r_edge
  real*8 :: rho_b_0

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

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CartSymGN(dummy,cctkGH,'lapse::lapse_evolve')
  call CartSymGN(dummy,cctkGH,'shift::shift_evolve')

  !not sure why this is being done
  !trK = 0.D0

  !write(*,*) "SETUP Axx y parts: ",Axx(17,1,2),Axx(2,2,16),Axx(17,3,2)
  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')
  if(Symmetry .eq. AXISYM) then
!     call axibc_tensor(ext,X,Y,Z,Axx,Axy,Axz,Ayy,Ayz,Azz)
!     call axibc_scalar(ext,X,Y,Z,K)
!     call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
     call CCTK_VarIndex(index,'BSSN::Axx')
     call BndCartoon2DVI(dummy, cctkGH, 2, index)
  end if
  !write(*,*) "AFTER ROT SETUP Axx y parts: ",Axx(17,1,2),Axx(2,2,16),Axx(17,3,\2)

  !detmin_l = 0.D0
  !detmax_l = 0.D0

 !====================================================
  ! Convert to tilde metric and invert tilde metric
  !====================================================
 !call convert(ext,phi,gxx,gxy,gxz,gyy,gyz,gzz, &
  !     gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
   !    detmin_l,detmax_l)
! Invert metric...                                                                              
!    
              
!  gupxx =   ( gyy * gzz - gyz * gyz )   
!  gupxy = - ( gxy * gzz - gyz * gxz )   
!  gupxz =   ( gxy * gyz - gyy * gxz )  
!  gupyy =   ( gxx * gzz - gxz * gxz )  
!  gupyz = - ( gxx * gyz - gxy * gxz )  
!  gupzz =   ( gxx * gyy - gxy * gxy )  
   
  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')
  if(Symmetry .eq. AXISYM) then
!     call CCTK_VarIndex(index,'BSSN::gxx')
!     call BndCartoon2DVI(dummy, cctkGH, 2, index)
!     call CCTK_VarIndex(index,'BSSN::gupxx')
!     call BndCartoon2DVI(dummy, cctkGH, 2, index)
!     call CCTK_VarIndex(index,'BSSN::phi')
!     call BndCartoon2DVI(dummy, cctkGH, 0, index)
  end if

 !====================================================
  ! Get Gamma^i
  !====================================================
!  call setgamma(ext,X,Y,Z, &
!       phi,gxx,gxy,gxz,gyy,gyz,gzz, &
!       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
!       Gammax,Gammay,Gammaz, &
!       Symmetry)
  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')
  if(Symmetry .eq. AXISYM) then
     call CCTK_VarIndex(index,'BSSN::Gammax')
     call BndCartoon2DVI(dummy, cctkGH, 1, index)
  end if
 call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')
 
  n1 = 0
  n2 = 0
  n3 = 0
  mf = 1
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammax,n1,n2,n3,mf,Symmetry)
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammay,n1,n2,n3,mf,Symmetry)
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammaz,n1,n2,n3,mf,Symmetry)
  !Next line NEEDED!
  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')

!  Axx = Axx*exp(-4.d0*phi)
!  Axy = Axy*exp(-4.d0*phi)
!  Axz = Axz*exp(-4.d0*phi)
!  Ayy = Ayy*exp(-4.d0*phi)
!  Ayz = Ayz*exp(-4.d0*phi)
!  Azz = Azz*exp(-4.d0*phi)

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
  write(*,*) "rho_b_atm: ",rho_b_atm

  n_poly = Pr(1,1,1)
  k_poly = Pr(1,1,2)
  kappaa = Pr(1,1,3)
  kappas = Pr(1,1,4)
  aRmB4  = Pr(1,1,5)
  r_edge = Pr(1,1,6)
  rho_b_0 = Pr(1,1,7)
 gamma_th = 1.d0 + 1.d0/n_poly


  write(*,*) "n_poly is set to: ",n_poly
  write(*,*) "r_edge: ",r_edge

  !
  ! Compute the initial matter profile
  !

  call compute_OS_rad_hybrid(ext,X,Y,Z, &
       neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th, &
       lapm1,phi,shiftx,shifty,shiftz, &
       gxx,gxy,gxz,gyy,gyz,gzz, &
       rho,S,Sx,Sy,Sz, &
       Sxx,Sxy,Sxz,Syy,Syz,Szz, &
       rho_star,tau,st_x,st_y,st_z, mhd_st_x, mhd_st_y, mhd_st_z, &
       P,w,vx,vy,vz,rho_b, &
       u0,h,rho_b_atm,PhysicalRadius,eps_flag, K_poly, n_poly, &
       E_rad,F_rad0,F_radx,F_rady,F_radz,F,tau_rad,S_rad_x,S_rad_y,S_rad_z,Po4PiB,PoRho,rho_b_0,r_edge,width,rad_evolve_enable)
  write(*,*) "test it"

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vars')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vs')
  call CartSymGN(dummy,cctkGH,'rad_evolve::rad_vars')
  call CartSymGN(dummy,cctkGH,'rad_evolve::conserved_rad_vars')

  call CCTK_VarIndex(index,"mhd_evolve::P")
  call CCTK_ReductionHandle(handle,"maximum")
  if (handle .gt. 0) then
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,reduction_value,1,index)
     if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH).eq.0) then
        print *,"Maximum value of P is ",reduction_value
     end if
  else
     call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
  end if
  P_max = reduction_value
  pfloor = P_fact*P_max


  !
  ! Setup initial EM fields (from a vector potential A_phi)
  !
  Bx = 0.d0
  By = 0.d0
  Bz = 0.d0
 

  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vars')
  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_vars')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vs')
  
  !is this correct?
  call CartSymGN(dummy,cctkGH,'rad_evolve::rad_vars')
  call CCTK_SyncGroup(dummy,cctkGH,'rad_evolve::rad_vars')
  call CartSymGN(dummy,cctkGH,'rad_evolve::conserved_rad_vars')
  call CCTK_SyncGroup(dummy,cctkGH,'rad_evolve::conserved_rad_vars')

  
 if(Symmetry .eq. AXISYM) then
     call CCTK_VarIndex(index,'mhd_evolve::st_x')
     call BndCartoon2DVI(dummy, cctkGH, 1, index)
     call CCTK_VarIndex(index,'mhd_evolve::mhd_st_x')
     call BndCartoon2DVI(dummy, cctkGH, 1, index)
     call CCTK_VarIndex(index,'mhd_evolve::Bx')
     call BndCartoon2DVI(dummy, cctkGH, 1, index)
     call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::tau')

     !Is this correct?
     call CCTK_VarIndex(index,'rad_evolve::S_rad_x')
     call BndCartoon2DVI(dummy, cctkGH, 1, index)
     call BndCartoon2DVN(dummy, cctkGH, 0, 'rad_evolve::tau_rad')
  end if

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
 
 ! Set excision_zone_gf to avoid valgrind memory errors
  excision_zone_gf = 0


 call primitive_vars_hybrid2(ext,X,Y,Z, &
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
       sbt,sbx,sby,sbz, &
       proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
       glob_imax,glob_jmax,glob_kmax,Symmetry,pfloor,excision_enable,excision_zone_gf,tau_stildefix_enable)

  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vars')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vs')
  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CartSymGN(dummy,cctkGH,'lapse::lapse_evolve')
  call CartSymGN(dummy,cctkGH,'shift::shift_evolve')
  
  !is this correct?
  call CartSymGN(dummy,cctkGH,'rad_evolve::rad_vars')
  call CartSymGN(dummy,cctkGH,'rad_evolve::conserved_rad_vars')

  !Compute K_ij's for apparent horizon finder:
  Psi = exp(phi)
  Kxx = Psi*Psi*Psi*Psi * (Axx + (1.D0/3.D0) * gxx * trK)
  Kxy = Psi*Psi*Psi*Psi * (Axy + (1.D0/3.D0) * gxy * trK)
  Kxz = Psi*Psi*Psi*Psi * (Axz + (1.D0/3.D0) * gxz * trK)
  Kyy = Psi*Psi*Psi*Psi * (Ayy + (1.D0/3.D0) * gyy * trK)
  Kyz = Psi*Psi*Psi*Psi * (Ayz + (1.D0/3.D0) * gyz * trK)
  Kzz = Psi*Psi*Psi*Psi * (Azz + (1.D0/3.D0) * gzz * trK)

  lapset = 0.D0

    !particle tracer stuff
!  do index=1,narr
!     coord(index,1)=0.d0
!     coord(index,2)=index*(.95*r_edge/narr)
!     coord(index,3)=0.d0
!     coord(index,4)=0.d0
!  end do
  


  !initialize convergence testing stuff
  rho_b_conv=0.d0
end subroutine OS_rad_finalize_initialdata
