#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine disk_initialdata(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext
  real*8 			           :: dT,dX,dY,dZ,rhos_max
  real*8 			           :: xmin,ymin,zmin,xmax,ymax,zmax
  real*8                                   :: hc_mask_radius, fac, max_b2oP, maxb2, maxP, tau_max
  real*8                                   :: detmin_l, detmax_l
  real*8                                   :: rho_fail_max_step,M_fail_step
  logical :: fish_to_phys
  integer :: n1, n2, n3, mf
  integer :: handle,dummy,proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax
  integer :: proc_kmax,glob_imax,glob_jmax,glob_kmax
  integer :: index, Nfont, Nfont_l
  integer :: ierr,ONE,ZERO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  CCTK_REAL reduction_value

  parameter(ONE = 1.D0, ZERO = 0.D0)

  ext = cctk_lsh

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  glob_imax = ext(1)
  glob_jmax = ext(2)
  glob_kmax = ext(3)

  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"grid::X")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xmax,1,index)
  if(X(cctk_lsh(1),1,1) .eq. xmax) then
     proc_imax = glob_imax
  else 
     proc_imax = -1
  end if
  call CCTK_VarIndex(index,"grid::Y")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,ymax,1,index)
  if(Y(1,cctk_lsh(2),1) .eq. ymax) then
     proc_jmax = glob_jmax
  else 
     proc_jmax = -1
  end if
  call CCTK_VarIndex(index,"grid::Z")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,zmax,1,index)
  if(Z(1,1,cctk_lsh(3)) .eq. zmax) then
     proc_kmax = glob_kmax
  else 
     proc_kmax = -1
  end if

  call CCTK_ReductionHandle(handle,"minimum")
  call CCTK_VarIndex(index,"grid::X")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xmin,1,index)
  if(X(1,1,1) .eq. xmin) then
     proc_imin = 0
  else 
     proc_imin = -100
  end if
  call CCTK_VarIndex(index,"grid::Y")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,ymin,1,index)
  if(Y(1,1,1) .eq. ymin) then
     proc_jmin = 0
  else 
     proc_jmin = -100
  end if
  call CCTK_VarIndex(index,"grid::Z")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,zmin,1,index)
  if(Z(1,1,1) .eq. zmin) then
     proc_kmin = 0
  else 
     proc_kmin = -100
  end if

  ! Setup atmosphere:
  call disk_Setup_alt_Atmosphere(CCTK_PASS_FTOF)

  if (CCTK_MyProc(CCTKGH) == 0) write(*,*)  " Setting up Kerr-Schild Initial Data from file " 

  !|~~~~~> Set up *untilded* metric, extrinsic curvature, lapse,shifts
  if (radial_coordinate_type == 0) then
     call ks_initial_metric_fisheye(ext,X,Y,Z,PhysicalRadius,RadiusDerivative, &
          gxx,gxy,gxz,gyy,gyz,gzz,trK, &
          Axx,Axy,Axz,Ayy,Ayz,Azz, &
          lapm1,shiftx,shifty,shiftz,sam,Symmetry)

     fish_to_phys = .false.
     call trans_fish_phys_ttij2_fullgridfunction(ext,X,Y,Z,PhysicalRadius,RadiusDerivative, &
          gxx,gxy,gxz,gyy,gyz,gzz,fish_to_phys)

     call trans_fish_phys_ttij2_fullgridfunction(ext,X,Y,Z,PhysicalRadius,RadiusDerivative, &
          Axx,Axy,Axz,Ayy,Ayz,Azz,fish_to_phys)

     call trans_fish_phys_vector_fullgridfunction(ext,X,Y,Z,PhysicalRadius,RadiusDerivative, &
          shiftx,shifty,shiftz,fish_to_phys)

  else
     call ks_initial_metric_bl(ext,X,Y,Z,PhysicalRadius,RadiusDerivative, &
          gxx,gxy,gxz,gyy,gyz,gzz,trK, &
          Axx,Axy,Axz,Ayy,Ayz,Azz, &
          lapm1,shiftx,shifty,shiftz,sam,Symmetry)
  end if

  ! Sync before you call setgamma (which takes derivatives)
  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')
  call CartSymGN(dummy,cctkGH,'shift::shift_vars')
  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')

  if(Symmetry .eq. AXISYM) then 
     call CCTK_VarIndex(index,'BSSN::Axx')
     call BndCartoon2DVI(dummy, cctkGH, 2, index)
  end if

  !|~~~~~> convert to tilded metric and initialize \phi and \Gamma^a
  !Convert to tilde metric and invert tilde metric:
  !  call convert_cpp(cctkGH,cctk_lsh,phi, &
  !       gxx,gxy,gxz,gyy,gyz,gzz, &
  !       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)
  call convert(ext,phi,gxx,gxy,gxz,gyy,gyz,gzz, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,gxx_p, &
       detmin_l,detmax_l)
  gxx_p=gxx

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')

  !Get Gamma^i:
  !  call setgamma_cpp(cctkGH,cctk_lsh, &
  !       dx,dy,dz, &
  !       gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
  !       Gammax, Gammay, Gammaz)

  ! Here we use gxx_p as temporary storage:
  call setgamma(ext,X,Y,Z,gxx_p, &
       phi,gxx,gxy,gxz,gyy,gyz,gzz, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
       Gammax,Gammay,Gammaz, &
       Symmetry)
  gxx_p = gxx

  n1 = 0 
  n2 = 0 
  n3 = 0 
  mf = 1    
  call update_falloff(ext,X,Y,Z,PhysicalRadius,Gammax,n1,n2,n3,mf,Symmetry);
  call update_falloff(ext,X,Y,Z,PhysicalRadius,Gammay,n1,n2,n3,mf,Symmetry);
  call update_falloff(ext,X,Y,Z,PhysicalRadius,Gammaz,n1,n2,n3,mf,Symmetry);

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')
  call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')
  call CCTK_SyncGroup(dummy,cctkGH,'lapse::lapse_vars')

  call Derivs(ext,X,Y,Z,dX,dY,dZ,phi,phix,phiy,phiz,Symmetry)
  call Derivs(ext,X,Y,Z,dX,dY,dZ,lapm1,lapsex,lapsey,lapsez,Symmetry)
  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')
  call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')
  call CCTK_SyncGroup(dummy,cctkGH,'lapse::lapse_vars')

  !|~~~~~> covert traceless K_{ij} to tilded traceless K_{ij}
  Axx = Axx*exp(-4.D0*phi)
  Axy = Axy*exp(-4.D0*phi)
  Axz = Axz*exp(-4.D0*phi)
  Ayy = Ayy*exp(-4.D0*phi)
  Ayz = Ayz*exp(-4.D0*phi)
  Azz = Azz*exp(-4.D0*phi)

  !|~~~~~> Set everything else to Zero!
  rho = 0.d0
  S = 0.d0
  Sx = 0.d0
  Sy = 0.d0
  Sz = 0.d0
  Sxx = 0.d0
  Sxy = 0.d0
  Sxz = 0.d0
  Syy = 0.d0
  Syz = 0.d0
  Szz = 0.d0

  Pl = 0.d0

  !!
  !! Compute the initial matter profile
  !!
  ! This is ugly.  Is there another way?
  n_poly = 3.0

  ! Initialization, just in case.
  u0 = 1.d0
  h = 1.d0
  rho_star = 0.d0

  if (radial_coordinate_type == 0) then
     call setup_fishbone_moncrief_disk_spin(ext,X,Y,Z, &
          rho_star,tau,st_x,st_y,st_z,w,h,rho_b,P, &
          vx,vy,vz,u0, n_poly, ell, r_in, rho_b_atm_gf, K_poly, &
          sam,PhysicalRadius,RadiusDerivative, &
          pfloor_gf,lapm1,phi,shiftx,shifty,shiftz, &
          gxx,gxy,gxz,gyy,gyz,gzz,enable_HARM_energyvariable, proc_imin, proc_kmin)
  else
     call setup_fishbone_moncrief_disk_spin_bl(ext,X,Y,Z, &
          rho_star,tau,st_x,st_y,st_z,w,h,rho_b,P, &
          vx,vy,vz,u0, n_poly, ell, r_in, rho_b_atm_gf, K_poly, &
          sam,PhysicalRadius,RadiusDerivative, &
          pfloor_gf,lapm1,phi,shiftx,shifty,shiftz, &
          gxx,gxy,gxz,gyy,gyz,gzz,enable_HARM_energyvariable, proc_imin, proc_kmin)
  end if

  write(*,*) "PPPPPP:",P(2,2,2),tau(2,2,2)

  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_conservatives')
  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_primitives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_conservatives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_primitives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vs')

  ! Note that the regular old axibc stuff will not work on multiple processors.
  if (Symmetry == AXISYM) then
     call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::rho_b')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::P')
  end if

  !  Setup initial EM fields (from a vector potential A_phi)
  fac = 1.0
  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"mhd_evolve::rho_b")
  if (handle .gt. 0) then
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,reduction_value,1,index)
     if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH).eq.0) then
        print *,"Maximum value of rho_b is ",reduction_value
     end if
  else
     call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
  end if
  rho_b_max = reduction_value 

  ! b^2/P ---> Pr
  ! d2's ----> Aphi's an B^i phys's
  ! Warning:  the MRI wavelength thing will not work for spherical coords.
  call initial_emfields_disk(ext,r_in,rho_b_max,p_c,betam1,X,Y,Z, &
       phi,lapm1,shiftx,shifty,shiftz, &
       gxx,gxy,gxz,gyy,gyz,gzz, &
       rho_b,P,Bx,By,Bz, Pr, &
       st_x,st_y,st_z, &
       mhd_st_x,mhd_st_y,mhd_st_z, &
       Pl,u0,vx,vy,vz, &
       Symmetry,Sym_Bz,fac,enable_HARM_energyvariable, &
       PhysicalRadius,RadiusDerivative, &
       d2rho_b_bck,d2P_bck,d2vx_bck,d2vy_bck, &
       d2rho_b_fwd,d2P_fwd,d2vx_fwd)

  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_conservatives')
  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_primitives')
  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')

  ! Note that the regular old axibc stuff will not work on multiple processors.
  if (Symmetry == AXISYM) then
     call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::Pr')
  end if

  ! To set the strength of the field by choosing the minimum
  ! value of the plasma beta parameter:
  call CCTK_VarIndex(index,"mhd_evolve::Pr")
  call CCTK_ReductionHandle(handle,"maximum")
  if (handle .gt. 0) then
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,reduction_value,1,index)
     if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH).eq.0) then
        print *,"Maximum value of Pr is ",reduction_value
     end if
  else
     call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
  end if
  !  max_b2oP = reduction_value 

  maxb2 = reduction_value

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

  maxP = reduction_value
  max_b2oP = maxb2/maxP

  fac = sqrt(betam1/(max_b2oP/2.D0))

  ! To set the field strength according to the average 
  ! value of the beta parameter:
!!$// e^(6 phi) ---> rho_br
!!$      call f_volume_element_disk(ext,X,Y,Z,phi,rho_br,rho_b,
!!$				rho_bl,rho_tiny)
!!$      PV = weighted_integral(P,rho_br)
!!$      PemV = weighted_integral(Pr,rho_br)*0.5
!!$
!!$      if (betam1 .eq. 0.0) then
!!$	  fac = 0.0;
!!$      else 
!!$	  fac = sqrt(betam1*(PV/PemV))
!!$      end if 

  !  Renormalize the B-fields
  call initial_emfields_disk(ext,r_in,rho_b_max,p_c,betam1,X,Y,Z, &
       phi,lapm1,shiftx,shifty,shiftz, &
       gxx,gxy,gxz,gyy,gyz,gzz, &
       rho_b,P,Bx,By,Bz,Pr, &
       st_x,st_y,st_z, &
       mhd_st_x,mhd_st_y,mhd_st_z, &
       tau,u0,vx,vy,vz, &
       Symmetry,Sym_Bz,fac,enable_HARM_energyvariable, &
       PhysicalRadius,RadiusDerivative)
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_conservatives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_primitives')
  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_conservatives')
  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_primitives')

  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_sync_lr_hydro_quantities')

  ! Note that the regular old axibc stuff will not work on multiple processors.
  if (Symmetry == AXISYM) then
     call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::Pr')
  end if

  ! To set the strength of the field by choosing the minimum
  ! value of the plasma beta parameter:

  call CCTK_VarIndex(index,"mhd_evolve::Pr")
  call CCTK_ReductionHandle(handle,"maximum")
  if (handle .gt. 0) then
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,reduction_value,1,index)
     if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH).eq.0) then
        print *,"Maximum value of Pr is ",reduction_value
     end if
  else
     call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
  end if

  write(*,*) 'before axibc, u0(2,2,2) = ', u0(2,2,2)

  if (Symmetry==AXISYM) then
     call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::h')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::w')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::u0')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::rho_star')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::P')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::sbt')
     call CCTK_VarIndex(index,'mhd_evolve::sbx')
     call BndCartoon2DVI(dummy, cctkGH, 1, index)
     call CCTK_VarIndex(index,'mhd_evolve::vx')
     call BndCartoon2DVI(dummy, cctkGH, 1, index)
  end if

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

  call CCTK_VarIndex(index,"mhd_evolve::tau")
  call CCTK_ReductionHandle(handle,"maximum")
  if (handle .gt. 0) then
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,reduction_value,1,index)
     if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH).eq.0) then
        print *,"Maximum value of tau is ",reduction_value
     end if
  else
     call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
  end if
  tau_max = reduction_value 
  tau_atm = tau_fact*tau_max

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

  rho_b_atm = rho_fact * rho_b_max

  ! Let's initialize the Fontfix_tracker_gf array
  Fontfix_tracker_gf = 0.d0

  h_p = h

  if(enable_primitives_disk==0) then
     call primitive_vars_hybrid2(ext,cctk_nghostzones,X,Y,Z, &
          rho_star,tau,st_x,st_y,st_z, &
          mhd_st_x,mhd_st_y,mhd_st_z,neos, &
          rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, &
          w,w,rho_b,rho,P,h,Sx,Sy,Sz, &
          Sxx,Sxy,Sxz,Syy,Syz,Szz, &
          phi,lapm1,shiftx,shifty,shiftz, &
          gxx,gxy,gxz,gyy,gyz,gzz, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
          h,u0,rho_b_max,rho_b_atm, &
          rho_fail_max_step,M_fail_step,rhos_max, &
          Bx,By,Bz,Ex,Ey,Ez, &
          vx,vy,vz, &
          sbt,sbx,sby,sbz,temp4, &
          proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
          glob_imax,glob_jmax,glob_kmax,Symmetry,pfloor,excision_enable, &
          excision_zone_gf,tau_stildefix_enable,tau_atm,0)
  else if(enable_primitives_disk==1) then
     if(use_HARM_primitives==0) then
        call primitive_vars_alt_disk(ext,X,Y,Z, &
             rho_star,tau,st_x,st_y,st_z, &
             mhd_st_x,mhd_st_y,mhd_st_z,neos, &
             rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, &
             w,rho_b,rho,P,h,Sx,Sy,Sz, &
             Sxx,Sxy,Sxz,Syy,Syz,Szz, &
             phi,lapm1,shiftx,shifty,shiftz, &
             gxx,gxy,gxz,gyy,gyz,gzz, &
             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
             h,u0,rho_b_max,rho_b_atm_gf, &
             rho_fail_max_step,M_fail_step, &
             Bx,By,Bz,Ex,Ey,Ez, &
             vx,vy,vz, &
             sbt,sbx,sby,sbz, &
             proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax, &
             glob_imax,glob_jmax,glob_kmax,Symmetry, &
             Fontfix_tracker_gf,pfloor_gf, Nfont_l, &
             enable_HARM_energyvariable, excision_zone_gf, force_font_fix_fail, excision_enable)
     end if
!!$     else if(use_HARM_primitives==1) then
!!$        ! call the C function--a wrapper for the Noble et al code.
!!$
!!$        ! At one time, the second line of the function call was:
!!$        ! rho_b, rho_b_p, P, P_p, vx, vy, vz, &
!!$        call primitive_vars_harm(NFont_L, CCTKGH, X, Y, Z, ext, &
!!$             rho_star, tau, mhd_st_x, mhd_st_y, mhd_st_z, &
!!$             rho_b, rho_b, P, P, vx, vy, vz, &
!!$             u0, Bx, By, Bz, lapm1, phi, &
!!$             shiftx, shifty, shiftz, Fontfix_tracker_gf, rho_b_atm_gf, &
!!$             gxx, gxy, gxz, gyy, gyz, gzz, &
!!$             gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
!!$             Symmetry,  excision_radius, &
!!$             glob_imax, glob_jmax, glob_kmax, & 
!!$             proc_imax, proc_jmax, proc_kmax, & 
!!$             proc_imin, proc_jmin, proc_kmin)
!!$
!!$        call harm_impose_floors(ex,rho_b_atm_gf,pfloor_gf,rho_b,P, &
!!$             vx,vy,vz,Bx,By,Bz, &
!!$             rho_star,tau,mhd_st_x,mhd_st_y,mhd_st_z, &
!!$             lapm1,phi,shiftx,shifty,shiftz, &
!!$             gxx,gxy,gxz,gyy,gyz,gzz, &
!!$             neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th)
!!$
!!$        ! make sure you have the failures array stored in Fontfix_tracker_gf.
!!$        call repair_failures_mhd_disk(ex,X,Z,gamma_th,Fontfix_tracker_gf,rho_b,P, &
!!$             vx,vy,vz,u0,w,h,rho_star,tau, &
!!$             st_x,st_y,st_z,mhd_st_x,mhd_st_y,mhd_st_z, &
!!$             rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz, &
!!$             lapm1,shiftx,shifty,shiftz,phi, &
!!$             gxx,gxy,gxz,gyy,gyz,gzz, &
!!$             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,Bx,By,Bz, &
!!$             sbt,sbx,sby,sbz,rho_b_atm_gf, &
!!$             neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
!!$             proc_imin,proc_jmin,proc_kmin, &
!!$             proc_imax,proc_jmax,proc_kmax, & 
!!$             glob_imax,glob_jmax,glob_kmax, &
!!$             Symmetry,pfloor_gf,enable_HARM_energyvariable,excision_radius)
!!$
!!$        !!&     call repair_failures_alt(ex,X,Z,gamma,Fontfix_tracker_gf,rho_b,P, &
!!$        !!&                          vx,vy,vz,u0,w,h,rho_star,tau, &
!!$        !!&                          st_x,st_y,st_z,mhd_st_x,mhd_st_y,mhd_st_z, &
!!$        !!&                          rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz, &
!!$        !!&                          lapm1,shiftx,shifty,shiftz,phi, &
!!$        !!&                          gxx,gxy,gxz,gyy,gyz,gzz, &
!!$        !!&                          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,Bx,By,Bz, &
!!$        !!&                          sbt,sbx,sby,sbz,rho_b_atm_gf, &
!!$        !!&                          neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
!!$        !!&                          proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax,  &
!!$        !!&                          glob_imax,glob_jmax,glob_kmax,Symmetry,pfloor_gf,enable_HARM_energyvariable,excision_radius)
!!$
!!$        call calculate_aux_harm(ex,gamma_th,rho_b,P, &
!!$             vx,vy,vz,u0,w,h,rho_star,tau, &
!!$             st_x,st_y,st_z,mhd_st_x,mhd_st_y,mhd_st_z, &
!!$             rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz, &
!!$             lapm1,shiftx,shifty,shiftz,phi, &
!!$             gxx,gxy,gxz,gyy,gyz,gzz, &
!!$             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,Bx,By,Bz, &
!!$             sbt,sbx,sby,sbz, &
!!$             neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,proc_kmin)
!!$     end if

     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_VarIndex(index,"disk::Nfont_l")
     call CCTK_ReduceLocScalar (ierr,cctkGH,-1,handle,Nfont_l,Nfont,CCTK_VARIABLE_INT)

     if (CCTK_MyProc(cctkGH) == 0) write(*,*) 'Fixed ', Nfont, " zones in initial data."
  end if

  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_conservatives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_primitives')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vs')
  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_conservatives')
  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_primitives')
  call CCTK_SyncGroup(dummy,cctkGH,'mhd_evolve::mhd_vs')

  call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::w')
  call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::h')
  call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::u0')
  call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::rho_star')
  call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::P')
  call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::sbt')
  call CCTK_VarIndex(index,'mhd_evolve::sbx')
  call BndCartoon2DVI(dummy, cctkGH, 1, index)
  call CCTK_VarIndex(index,'mhd_evolve::vx')
  call BndCartoon2DVI(dummy, cctkGH, 1, index)

  ! You should probably convert these diagnostics
!!$  adm::AH_Manager.Find_Horizon(level,0.0,PhysTime);
!!$  double J_dyn = 0.0;
!!$//   adm::AH_Manager.dynamical_horizon_J(J_dyn,phi,
!!$// 				      gupxx,gupxy,gupxz,
!!$// 				      gupyy,gupyz,gupzz,
!!$// 				      Axx,Axy,Axz,
!!$// 				      Ayy,Ayz,Azz,trK,
!!$// 				      level,IDENT);
!!$  adm::AH_Manager.dynamical_horizon_J_ana(J_dyn,
!!$					  n_tot_horizon, N_theta, N_phi,
!!$					  xn, yn, zn, nx_g, ny_g, nz_g, nn,
!!$					  phi,
!!$					  gupxx,gupxy,gupxz,
!!$					  gupyy,gupyz,gupzz,
!!$					  Axx,Axy,Axz,
!!$					  Ayy,Ayz,Azz,trK,
!!$					  level,IDENT);
!!$  if (I_Print) cout << "J_dynamical_horizon = " << J_dyn << endl;

  write(*,*) "END OF disk_initialdata()"

end subroutine disk_initialdata
