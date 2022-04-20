#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine binfile_check_do_checkpoint(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,local_origin
  real*8, dimension(3)                     :: xyz_glob_min,xyz_glob_max
  real*8 			           :: dT,dX,dY,dZ,P_max,rho_max,rhos_max,rho_fail_max_step,M_fail_step
  real*8 			           :: detmin_l,detmax_l,Kmin_l,Kmax_l,xmin,ymin,zmin,xmax,ymax,zmax
  real*8 			           :: lapsemin
  integer :: i,n1,n2,n3,mf
  integer :: iadj,jadj,kadj
  integer :: handle,dummy,proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax,glob_imax,glob_jmax,glob_kmax
  integer :: index
  integer :: ierr,ONE,ZERO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  CCTK_REAL reduction_value
  character                                :: filename*25
  parameter(ONE = 1.D0, ZERO = 0.D0)


  call CCTK_ReductionHandle(handle,"minimum")
  call CCTK_VarIndex(index,"lapse::lapm1")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,lapsemin,1,index)
  lapsemin = lapsemin + 1

  write(*,*) "Minimum Lapse: ",lapsemin


!I changed this here.
  if(cctk_time .lt. 0.0) then

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

     iadj = 0
     jadj = 0
     kadj = 0
     if(Symmetry == AXISYM) then
        if(X(1,1,1) .lt. 0.D0) then
           iadj = iadj + 1
        end if
        if(Z(1,1,1) .lt. 0.D0) then
           kadj = kadj + 1
        end if
     end if
     if(Symmetry == EQUATORIAL) then
        if(Z(1,1,1) .lt. 0.D0) then
           kadj = kadj + 1
        end if
     end if

     filename = 'phi'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,phi,filename,Symmetry)
     filename = 'trK'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,trK,filename,Symmetry)
     filename = 'gxx'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,gxx,filename,Symmetry)
     filename = 'gxy'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,gxy,filename,Symmetry)
     filename = 'gxz'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,gxz,filename,Symmetry)
     filename = 'gyy'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,gyy,filename,Symmetry)
     filename = 'gyz'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,gyz,filename,Symmetry)
     filename = 'gzz'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,gzz,filename,Symmetry)

     filename = 'Axx'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,Axx,filename,Symmetry)
     filename = 'Axy'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,Axy,filename,Symmetry)
     filename = 'Axz'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,Axz,filename,Symmetry)
     filename = 'Ayy'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,Ayy,filename,Symmetry)
     filename = 'Ayz'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,Ayz,filename,Symmetry)
     filename = 'Azz'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,Azz,filename,Symmetry)

     filename = 'Gammax'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,Gammax,filename,Symmetry)
     filename = 'Gammay'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,Gammay,filename,Symmetry)
     filename = 'Gammaz'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,Gammaz,filename,Symmetry)

     filename = 'lapm1'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,lapm1,filename,Symmetry)

     filename = 'shiftx'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,shiftx,filename,Symmetry)
     filename = 'shifty'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,shifty,filename,Symmetry)
     filename = 'shiftz'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,shiftz,filename,Symmetry)

     filename = 'lapset'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,lapset,filename,Symmetry)

     filename = 'shiftxt'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,shiftxt,filename,Symmetry)
     filename = 'shiftyt'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,shiftyt,filename,Symmetry)
     filename = 'shiftzt'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,shiftzt,filename,Symmetry)

     filename = 'rho_star'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,rho_star,filename,Symmetry)
     filename = 'tau'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,tau,filename,Symmetry)

     filename = 'mhd_st_x'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,mhd_st_x,filename,Symmetry)
     filename = 'mhd_st_y'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,mhd_st_y,filename,Symmetry)
     filename = 'mhd_st_z'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,mhd_st_z,filename,Symmetry)

     filename = 'Bx'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,Bx,filename,Symmetry)
     filename = 'By'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,By,filename,Symmetry)
     filename = 'Bz'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,Bz,filename,Symmetry)

     filename = 'rho_b'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,rho_b,filename,Symmetry)
     filename = 'P'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,P,filename,Symmetry)
     filename = 'h'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,h,filename,Symmetry)


!add radiation stuff
     filename = 'E_rad'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,E_rad,filename,Symmetry)
     filename = 'F_rad0'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,F_rad0,filename,Symmetry)
     filename = 'F_radx'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,F_radx,filename,Symmetry)
     filename = 'F_rady'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,F_rady,filename,Symmetry)
     filename = 'F_radz'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,F_radz,filename,Symmetry)
      filename = 'tau_rad'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,tau_rad,filename,Symmetry)
      filename = 'S_rad_x'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,S_rad_x,filename,Symmetry)
     filename = 'S_rad_y'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,S_rad_y,filename,Symmetry)
     filename = 'S_rad_z'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,S_rad_z,filename,Symmetry)

     !Write miscellaneous data required for checkpointing:
     vxr(2,1,2) = K_poly
     vxr(2,1,3) = neos
     do i=1,neos
        vxr(2,2,i+1) = rho_tab(i)
        vxr(3,2,i+1) = P_tab(i)
        vxr(4,2,i+1) = eps_tab(i)
        vxr(5,2,i+1) = k_tab(i)
        vxr(6,2,i+1) = gamma_tab(i)
     end do

     vxr(5,2,neos+2) = k_tab(neos+1)
     vxr(6,2,neos+2) = gamma_tab(neos+1)

     write(*,*) rho_tab
     write(*,*) P_tab
     write(*,*) eps_tab
     write(*,*) k_tab
     write(*,*) gamma_tab
     write(*,*) neos
     write(*,*) K_poly

     filename = 'misc'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,vxr,filename,Symmetry)


     !Write miscellaneous radiation data required for checkpointing:
     vxl(2,1,2) = aRmB4
     vxl(2,1,3) = kappaa
     vxl(2,1,4) = kappas
     
     !add some more miscellaneous data
     vxl(2,1,5) = pfloor
     vxl(2,1,6) = rho_b_atm
 
     filename = 'misc_rad'
     call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,vxl,filename,Symmetry)


     write(*,*) "wrote binfiles!"


     !  call outputbinfile(ierr,cctkGH,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,rho_b,filename,Symmetry)
     !  write(*,*) "wrote binfile!"
     stop

  end if
end subroutine binfile_check_do_checkpoint
