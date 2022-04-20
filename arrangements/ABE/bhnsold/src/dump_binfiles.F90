
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bhns_dump_binfiles(CCTK_ARGUMENTS)
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
  integer :: ierr
  real*8  :: ONE,ZERO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  CCTK_REAL reduction_value
  character           :: filename*25
  parameter(ONE = 1.D0, ZERO = 0.D0)

  if (mod(cctk_iteration,dump_every)==0) then
     ext = cctk_lsh
     dT = CCTK_DELTA_TIME
     dX = CCTK_DELTA_SPACE(1)
     dY = CCTK_DELTA_SPACE(2)
     dZ = CCTK_DELTA_SPACE(3)

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

!     call setup_global_coord_arrays(CCTK_PASS_FTOF)

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

     if(Symmetry==EQUATORIAL) then
        xyz_glob_min(3) = Zglobal(cctk_nghostzones(3)+1)
     end if

    ! Dump bin files
     open(19,file='DUMP_Input',status='old',ERR=110)
     if(1==0) then
110     write(*,*) "ERROR OPENING DUMP_Input!"; stop
     end if
     do 
        read(19,*,IOSTAT=ierr) filename 
        if (ierr .ne. 0) exit

    ! This is ugly. Need to fix it later.
       if (filename=='phi') then 
          call bhns_outputbinfile(ierr,cctkGH,cctk_nghostzones,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,phi,filename,Symmetry)
       elseif (filename=='rho_b') then
          call bhns_outputbinfile(ierr,cctkGH,cctk_nghostzones,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,rho_b,filename,Symmetry)
       elseif (filename=='rho_star') then
          call bhns_outputbinfile(ierr,cctkGH,cctk_nghostzones,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,rho_star,filename,Symmetry)
       elseif (filename=='P') then
          call bhns_outputbinfile(ierr,cctkGH,cctk_nghostzones,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,P,filename,Symmetry)
       elseif (filename=='gxx') then
          call bhns_outputbinfile(ierr,cctkGH,cctk_nghostzones,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,gxx,filename,Symmetry)
       elseif (filename=='vx') then
          call bhns_outputbinfile(ierr,cctkGH,cctk_nghostzones,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,vx,filename,Symmetry)
       elseif (filename=='vy') then
          call bhns_outputbinfile(ierr,cctkGH,cctk_nghostzones,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,vy,filename,Symmetry)
       elseif (filename=='vz') then
          call bhns_outputbinfile(ierr,cctkGH,cctk_nghostzones,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,vz,filename,Symmetry)
       elseif (filename=='lapm1') then
          call bhns_outputbinfile(ierr,cctkGH,cctk_nghostzones,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,lapm1,filename,Symmetry)
       elseif (filename=='Bx') then
          call bhns_outputbinfile(ierr,cctkGH,cctk_nghostzones,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,Bx,filename,Symmetry)
       elseif (filename=='By') then
          call bhns_outputbinfile(ierr,cctkGH,cctk_nghostzones,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,By,filename,Symmetry)
       elseif (filename=='Bz') then
          call bhns_outputbinfile(ierr,cctkGH,cctk_nghostzones,cctk_gsh,cctk_lsh,local_origin,cctk_time,xyz_glob_min,xyz_glob_max,cctk_iteration,Bz,filename,Symmetry)
       else 
	  write(*,*) "Grid function ",filename," can't be dumped..." 
	  write(*,*) "Either the grid function doesn't exist or it is not "
	  write(*,*) "recognized in dump_binfiles."
       end if
       
     end do
     close(19)

   end if

 end subroutine bhns_dump_binfiles
