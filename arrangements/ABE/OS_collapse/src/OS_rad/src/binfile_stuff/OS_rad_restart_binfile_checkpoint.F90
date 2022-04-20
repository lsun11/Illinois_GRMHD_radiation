
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine OS_rad_restart_binfile_checkpoint(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,local_origin
  real*8, dimension(3)                     :: xyz_glob_min,xyz_glob_max
  real*8 			           :: dT,dX,dY,dZ,P_max,rho_max,rhos_max,rho_fail_max_step,M_fail_step
  real*8 			           :: detmin_l,detmax_l,Kmin_l,Kmax_l,xmin,ymin,zmin,xmax,ymax,zmax
  integer :: i,n1,n2,n3,mf
  integer :: symm_iadj_glob,symm_jadj_glob,symm_kadj_glob
  integer :: handle,dummy,proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax,glob_imax,glob_jmax,glob_kmax
  integer :: index
  integer :: ierr,ONE,ZERO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  CCTK_REAL reduction_value

  character                                :: filename*25
  

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

  xyz_glob_min(1) = xmin
  xyz_glob_min(2) = ymin
  xyz_glob_min(3) = zmin

  xyz_glob_max(1) = xmax
  xyz_glob_max(2) = ymax
  xyz_glob_max(3) = zmax

  call setup_global_coord_arrays(CCTK_PASS_FTOF)

  do i=1,cctk_gsh(1)
     if(abs(X(1,1,1)-Xglobal(i)) .lt. dX*0.0001) then 
        local_origin(1) = i-1 ! -1 offset because outputbinfile() is a C file
     end if
  end do
  do i=1,cctk_gsh(2)
     if(abs(Y(1,1,1)-Yglobal(i)) .lt. dY*0.0001) then 
        local_origin(2) = i-1 ! -1 offset because outputbinfile() is a C file
     end if
  end do
  do i=1,cctk_gsh(3)
     if(abs(Z(1,1,1)-Zglobal(i)) .lt. dZ*0.0001) then 
        local_origin(3) = i-1 ! -1 offset because outputbinfile() is a C file
     end if
  end do

  symm_iadj_glob = 0
  symm_jadj_glob = 0
  symm_kadj_glob = 0
  if(Symmetry == AXISYM) then
     !adjustment here for symmetry considerations is necessary
     symm_iadj_glob = 1
     symm_kadj_glob = 1
  end if
  if(Symmetry == EQUATORIAL) then
     symm_kadj_glob = 1
  end if

  !Initialize quantities to zero so that ghost zone values are filled.  read_binfile does NOT fill ghostzones.
  call OS_initialize_to_zero(CCTK_PASS_FTOF)

  cctk_iteration = binfile_checkpoint_iteration

  filename = 'phi'
  call read_binfile(filename,binfile_checkpoint_iteration,phi,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'trK'
  call read_binfile(filename,binfile_checkpoint_iteration,trK,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'gxx'
  call read_binfile(filename,binfile_checkpoint_iteration,gxx,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'gxy'
  call read_binfile(filename,binfile_checkpoint_iteration,gxy,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'gxz'
  call read_binfile(filename,binfile_checkpoint_iteration,gxz,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'gyy'
  call read_binfile(filename,binfile_checkpoint_iteration,gyy,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'gyz'
  call read_binfile(filename,binfile_checkpoint_iteration,gyz,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'gzz'
  call read_binfile(filename,binfile_checkpoint_iteration,gzz,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);

  filename = 'Axx'
  call read_binfile(filename,binfile_checkpoint_iteration,Axx,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'Axy'
  call read_binfile(filename,binfile_checkpoint_iteration,Axy,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'Axz'
  call read_binfile(filename,binfile_checkpoint_iteration,Axz,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'Ayy'
  call read_binfile(filename,binfile_checkpoint_iteration,Ayy,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'Ayz'
  call read_binfile(filename,binfile_checkpoint_iteration,Ayz,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'Azz'
  call read_binfile(filename,binfile_checkpoint_iteration,Azz,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);

  filename = 'Gammax'
  call read_binfile(filename,binfile_checkpoint_iteration,Gammax,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'Gammay'
  call read_binfile(filename,binfile_checkpoint_iteration,Gammay,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'Gammaz'
  call read_binfile(filename,binfile_checkpoint_iteration,Gammaz,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);

  filename = 'lapm1'
  call read_binfile(filename,binfile_checkpoint_iteration,lapm1,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);

  filename = 'shiftx'
  call read_binfile(filename,binfile_checkpoint_iteration,shiftx,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'shifty'
  call read_binfile(filename,binfile_checkpoint_iteration,shifty,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'shiftz'
  call read_binfile(filename,binfile_checkpoint_iteration,shiftz,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);

  filename = 'lapset'
  call read_binfile(filename,binfile_checkpoint_iteration,lapset,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);

  filename = 'shiftxt'
  call read_binfile(filename,binfile_checkpoint_iteration,shiftxt,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'shiftyt'
  call read_binfile(filename,binfile_checkpoint_iteration,shiftyt,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'shiftzt'
  call read_binfile(filename,binfile_checkpoint_iteration,shiftzt,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);

  filename = 'rho_star'
  call read_binfile(filename,binfile_checkpoint_iteration,rho_star,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'tau'
  call read_binfile(filename,binfile_checkpoint_iteration,tau,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);

  filename = 'mhd_st_x'
  call read_binfile(filename,binfile_checkpoint_iteration,mhd_st_x,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'mhd_st_y'
  call read_binfile(filename,binfile_checkpoint_iteration,mhd_st_y,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'mhd_st_z'
  call read_binfile(filename,binfile_checkpoint_iteration,mhd_st_z,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);

  filename = 'Bx'
  call read_binfile(filename,binfile_checkpoint_iteration,Bx,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'By'
  call read_binfile(filename,binfile_checkpoint_iteration,By,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'Bz'
  call read_binfile(filename,binfile_checkpoint_iteration,Bz,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);

  filename = 'rho_b'
  call read_binfile(filename,binfile_checkpoint_iteration,rho_b,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'P'
  call read_binfile(filename,binfile_checkpoint_iteration,P,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'h'
  call read_binfile(filename,binfile_checkpoint_iteration,h,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  

  !radiation stuff
  filename = 'E_rad'
  call read_binfile(filename,binfile_checkpoint_iteration,E_rad,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'F_radx'
  call read_binfile(filename,binfile_checkpoint_iteration,F_radx,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
   filename = 'F_rady'
  call read_binfile(filename,binfile_checkpoint_iteration,F_rady,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
   filename = 'F_radz'
  call read_binfile(filename,binfile_checkpoint_iteration,F_radz,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
   filename = 'tau_rad'
  call read_binfile(filename,binfile_checkpoint_iteration,tau_rad,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
  filename = 'S_rad_x'
  call read_binfile(filename,binfile_checkpoint_iteration,S_rad_x,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
   filename = 'S_rad_y'
  call read_binfile(filename,binfile_checkpoint_iteration,S_rad_y,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
   filename = 'S_rad_z'
  call read_binfile(filename,binfile_checkpoint_iteration,S_rad_z,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);
 
  
  !Set up all the EOS parameters:
  filename = 'misc'
  call read_binfile(filename,binfile_checkpoint_iteration,vxr,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);

  K_poly = 0.D0
  neos = 0.D0
  rho_tab = 0.D0
  P_tab = 0.D0
  eps_tab = 0.D0
  k_tab = 0.D0
  gamma_tab = 0.D0
  
  !set up miscellaneous radiation stuff
 filename = 'misc_rad'
  call read_binfile(filename,binfile_checkpoint_iteration,vxl,cctkGH,cctk_time,symm_iadj_glob,symm_jadj_glob,symm_kadj_glob,local_origin,cctk_lsh);

  aRmB4 = 0.D0
  kappaa = 0.D0
  kappas = 0.D0

  if((Symmetry==AXISYM .and. abs(X(1,1,1) - Xglobal(1)) .lt. dX*0.0001 .and. abs(Z(1,1,1) - Zglobal(1)) .lt. dZ*0.0001) &
       .or. (Symmetry==EQUATORIAL .and. abs(X(1,1,1) - Xglobal(1)) .lt. dX*0.0001 .and. abs(Y(1,1,1) - Yglobal(1)) .lt. dY*0.0001 &
       .and. abs(Z(1,1,1) - Zglobal(1)) .lt. dZ*0.0001) ) then
     K_poly = vxr(2,1,2)
     neos = vxr(2,1,3)
     
     do i=1,neos
        rho_tab(i) = vxr(2,2,i+1)
        P_tab(i) = vxr(3,2,i+1)
        eps_tab(i) = vxr(4,2,i+1)
        k_tab(i) = vxr(5,2,i+1)
        gamma_tab(i) = vxr(6,2,i+1)
     end do
     
     k_tab(neos+1) = vxr(5,2,neos+2)
     gamma_tab(neos+1) = vxr(6,2,neos+2)
 
     !miscellaneous radiation stuff
     aRmB4 = vxl(2,1,2)
     kappaa = vxl(2,1,3)
     kappas = vxl(2,1,4)
     pfloor = vxl(2,1,5)
     rho_b_atm = vxl(2,1,6)

  end if

  call CCTK_ReductionHandle(handle,"sum")
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,K_poly,K_poly,CCTK_VARIABLE_REAL)
  call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,neos,neos,CCTK_VARIABLE_INT)

  call CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,rho_tab,rho_tab,neos,CCTK_VARIABLE_REAL);
  call CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,P_tab,P_tab,neos,CCTK_VARIABLE_REAL);
  call CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,eps_tab,eps_tab,neos,CCTK_VARIABLE_REAL);
  call CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,k_tab,k_tab,neos+1,CCTK_VARIABLE_REAL);
  call CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,gamma_tab,gamma_tab,neos+1,CCTK_VARIABLE_REAL);
  call CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,aRmB4,aRmB4,1,CCTK_VARIABLE_REAL);
  call CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,kappaa,kappaa,1,CCTK_VARIABLE_REAL);
  call CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,kappas,kappas,1,CCTK_VARIABLE_REAL);
 call CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,pfloor,pfloor,1,CCTK_VARIABLE_REAL);
 call CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,rho_b_atm,rho_b_atm,1,CCTK_VARIABLE_REAL);



  write(*,*) "Finished reading binfiles!"

  write(*,*) "eos params K_poly:",K_poly
  write(*,*) "eos params neos:",neos
  write(*,*) "eos params rho_tab:",rho_tab
  write(*,*) "eos params eps_tab:",eps_tab
  write(*,*) "eos params P_tab:",P_tab
  write(*,*) "eos params k_tab:",k_tab
  write(*,*) "eos params gamma_tab:",gamma_tab

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

  

  write(*,*) "HELLO INSIDE OS_rad-restart_binfile_checkpoint.F90!"
  
  if(Symmetry==AXISYM) then
     call BndCartoon2DVN(dummy, cctkGH, 0, 'bssn::phi')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'bssn::trK')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'lapse::lapm1')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::tau')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'lapse::lapset')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::rho_star')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::rho_b')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::P')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::h')
     call CCTK_VarIndex(index,'mhd_evolve::mhd_st_x')
     call BndCartoon2DVI(dummy, cctkGH, 1, index) 
     call CCTK_VarIndex(index,'shift::shiftx')
     call BndCartoon2DVI(dummy, cctkGH, 1, index) 
     call CCTK_VarIndex(index,'mhd_evolve::Bx')
     call BndCartoon2DVI(dummy, cctkGH, 1, index) 
     call CCTK_VarIndex(index,'BSSN::Gammax')
     call BndCartoon2DVI(dummy, cctkGH, 1, index) 
     call CCTK_VarIndex(index,'shift::shiftxt')
     call BndCartoon2DVI(dummy, cctkGH, 1, index) 
     call CCTK_VarIndex(index,'BSSN::Axx')
     call BndCartoon2DVI(dummy, cctkGH, 2, index) 
     call CCTK_VarIndex(index,'BSSN::gxx')
     call BndCartoon2DVI(dummy, cctkGH, 2, index) 
     !radiation stuff
     call BndCartoon2DVN(dummy, cctkGH, 0, 'rad_evolve::E_rad')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'rad_evolve::F_rad0')
     call BndCartoon2DVN(dummy, cctkGH, 0, 'rad_evolve::tau_rad')
     call CCTK_VarIndex(index,'rad_evolve::F_radx')
     call BndCartoon2DVI(dummy, cctkGH, 1, index) 
     call CCTK_VarIndex(index,'rad_evolve::S_rad_x')
     call BndCartoon2DVI(dummy, cctkGH, 1, index) 
  else
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vars')
     call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
     call CartSymGN(dummy,cctkGH,'lapse::lapse_evolve')
     call CartSymGN(dummy,cctkGH,'shift::shift_evolve')
     !radition stuff
     call CartSymGN(dummy,cctkGH,'rad_evolve::rad_vars')
     call CartSymGN(dummy,cctkGH,'rad_evolve::conserved_rad_vars')
  end if
  
  !First, set excision_zone_gf:
  if(excision_enable==1) then
     call find_excision_zone(ext,X,Y,Z,excision_zone_gf,excision_radius,Symmetry)
  else
     excision_zone_gf = 0
  end if
  

  ! Next, compute primitives:
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
 
!Brian commented this out because we are now saving pfloor in a binfile
! pfloor = P_max*P_fact 
  
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

!  rho_b_atm = rho_fact*rho_max

  !Brian commented this out because we are now saving rho_b_atm in a binfile
 !rho_b_atm = 6.82930999999999986D-09

  w = 1.D0

  call invert_g(ext,gxx,gxy,gxz,gyy,gyz,gzz, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
       detmin_l,detmax_l)
  write(*,*) "inside primitives0", mhd_st_x(2,2,30),mhd_st_y(2,2,30),mhd_st_z(2,2,30),h(2,2,30),gupxy(2,2,30)

!Following gridfunctions must be set for primitive solver.  Otherwise there are memory errors.
  vx = 0.D0 
  vy = 0.D0 
  vz = 0.D0 
  u0 = 1.D0

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
       glob_imax,glob_jmax,glob_kmax,Symmetry,pfloor,excision_enable,excision_zone_gf)

  if(excision_enable==1) then
     !Next set up the excision masks
     call set_masks(ext,X,Y,Z,mskf,hcmskf,excision_radius,Symmetry)
     
     
     !Finally remove points interior to the mask
     call remove_interior2(ext,X,Y,Z,lapm1,excision_zone_gf,Symmetry)
     call remove_interior2(ext,X,Y,Z,phi,excision_zone_gf,Symmetry)
     call remove_interior2(ext,X,Y,Z,trK,excision_zone_gf,Symmetry)
     call remove_interior2(ext,X,Y,Z,rho_star,excision_zone_gf,Symmetry)
     call remove_interior2(ext,X,Y,Z,tau,excision_zone_gf,Symmetry)
     call remove_interior2(ext,X,Y,Z,mhd_st_x,excision_zone_gf,Symmetry)
     call remove_interior2(ext,X,Y,Z,mhd_st_y,excision_zone_gf,Symmetry)
     call remove_interior2(ext,X,Y,Z,mhd_st_z,excision_zone_gf,Symmetry)
     
  end if

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
       glob_imax,glob_jmax,glob_kmax,Symmetry,pfloor,excision_enable,excision_zone_gf)
  
  !Does this need to be here?
  if (rad_evolve_enable==1) then
     call primitive_radiation_cpp(cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,&
          X,Y,Z,gupxx, gupxy, gupxz, &
          gupyy, gupyz, gupzz,gxx,gxy,gxz,gyy,gyz,gzz,shiftx,shifty,shiftz, &
          S_rad_x, S_rad_y, S_rad_z, tau_rad, E_rad, F_radx, F_rady, F_radz, &
          F_rad0,vx,vy,vz,u0,lapm1,phi,rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, &
          Szz)
  endif

  !what about boundary condition stuff?

  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vars')
  call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vs')

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CartSymGN(dummy,cctkGH,'lapse::lapse_evolve')
  call CartSymGN(dummy,cctkGH,'shift::shift_evolve')
  !Following lines are (empirically speaking) ESSENTIAL for excision runs:
  call Derivs(ext,X,Y,Z,dX,dY,dZ,phi,phix,phiy,phiz,Symmetry)
  call Derivs(ext,X,Y,Z,dX,dY,dZ,lapm1,lapsex,lapsey,lapsez,Symmetry)

  if(Symmetry==AXISYM) then
     call CCTK_VarIndex(index,'lapse::lapsex')
     call BndCartoon2DVI(dummy, cctkGH, 1, index) 
     call CCTK_VarIndex(index,'bssn::phix')
     call BndCartoon2DVI(dummy, cctkGH, 1, index) 
  end if

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars') 

!following is needed.
  call sanitycheck_restore_Aij(ext,gxx,gxy,gxz,gyy,gyz,gzz, &
       Axx,Axy,Axz,Ayy,Ayz,Azz,detmin_l,detmax_l,Kmin_l,Kmax_l) 

  if(Symmetry .eq. AXISYM) then    

  
  call CCTK_VarIndex(index,'BSSN::Axx')
  call BndCartoon2DVI(dummy, cctkGH, 2, index) 

!     call ezaxibc_scalar(ext,x,y,z,phi,excision_zone_gf,Symmetry) 
!     call ezaxibc_scalar(ext,x,y,z,trK,excision_zone_gf,Symmetry) 
!     call ezaxibc_vector(ext,x,y,z,Gammax,Gammay,Gammaz,excision_zone_gf,Symmetry)
!     call ezaxibc_tensor(ext,x,y,z,gxx,gxy,gxz,gyy,gyz,gzz,excision_zone_gf,Symmetry) 
!     call ezaxibc_tensor(ext,x,y,z,Axx,Axy,Axz,Ayy,Ayz,Azz,excision_zone_gf,Symmetry) 
  end if


  !call ezaxibc_tensor(ext,x,y,z,Axx,Axy,Axz,Ayy,Ayz,Azz,0.D0,Symmetry)


  !Compute K_ij's for apparent horizon finder:
  Psi = exp(phi)
  Kxx = Psi*Psi*Psi*Psi * ( Axx + (1.D0/3.D0) * gxx * trK )
  Kxy = Psi*Psi*Psi*Psi * ( Axy + (1.D0/3.D0) * gxy * trK )
  Kxz = Psi*Psi*Psi*Psi * ( Axz + (1.D0/3.D0) * gxz * trK )
  Kyy = Psi*Psi*Psi*Psi * ( Ayy + (1.D0/3.D0) * gyy * trK )
  Kyz = Psi*Psi*Psi*Psi * ( Ayz + (1.D0/3.D0) * gyz * trK )
  Kzz = Psi*Psi*Psi*Psi * ( Azz + (1.D0/3.D0) * gzz * trK )
  
  
  !call magnetar_copy_to_prev_timelevel(CCTK_PASS_FTOF)

  !lapset = 0.D0
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

end subroutine OS_rad_restart_binfile_checkpoint



subroutine OS_initialize_to_zero(CCTK_ARGUMENTS) 
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  rho_b = 0.D0
  phi = 0.D0
  trK = 0.D0
  lapm1 = 0.D0
  tau = 0.D0
  lapset = 0.D0
  rho_star = 0.D0
  P = 0.D0
  h = 0.D0

  mhd_st_x = 0.D0
  mhd_st_y = 0.D0
  mhd_st_z = 0.D0
  
  Gammax = 0.D0
  Gammay = 0.D0
  Gammaz = 0.D0
  
  shiftx = 0.D0
  shifty = 0.D0
  shiftz = 0.D0
  
  Bx = 0.D0
  By = 0.D0
  Bz = 0.D0
  
  shiftxt = 0.D0
  shiftyt = 0.D0
  shiftzt = 0.D0
  
  gxx = 0.D0
  gxy = 0.D0
  gxz = 0.D0
  gyy = 0.D0
  gyz = 0.D0
  gzz = 0.D0
  
  Axx = 0.D0
  Axy = 0.D0
  Axz = 0.D0
  Ayy = 0.D0
  Ayz = 0.D0
  Azz = 0.D0

  !initialize radiation stuff
  E_rad = 0.D0
  F_rad0 = 0.D0
  F_radx = 0.D0
  F_rady = 0.D0
  F_radz = 0.D0
  tau_rad = 0.D0
  S_rad_x = 0.D0
  S_rad_y = 0.D0
  S_rad_z = 0.D0

end subroutine OS_initialize_to_zero
