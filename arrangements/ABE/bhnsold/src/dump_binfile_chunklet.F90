!=====================================================================================
! Binfile chunklet output.  Use the binfile unchunker on the chunklets to get
!   ordinary, DAGH-era binfiles.
!=====================================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bhns_dump_binfile_chunklet(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: local_origin,xyzmin_idx
  real*8 			           :: dX,dY,dZ
  real*8, dimension(3)                     :: xyzmin_coord,xyzmin_glob,deltaxyz
  integer :: handle,index,ierr
  integer :: i, Symm, iter, NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  CCTK_REAL reduction_value
  character           :: fname*25
  integer :: whichproc,num_procs,myproc,dummy

  call CCTK_ReductionHandle(handle,"sum")
  dummy=1

  if (mod(cctk_iteration,chunklet_dump_every)==0) then
     !=====================================================================================
     ! First fill in the metadata so that the chunklets can be unchunked.
     Symm=Symmetry
     iter=cctk_iteration
     deltaxyz = CCTK_DELTA_SPACE

     xyzmin_idx = 2000
     do i=1,cctk_gsh(1)
        if(abs(X(1,1,1)-Xglobal(i)) .lt. deltaxyz(1)*0.0001) then 
           xyzmin_idx(1) = i-1 ! -1 offset because binfile_chunklet() is a C++ function
        end if
     end do
     do i=1,cctk_gsh(2)
        if(abs(Y(1,1,1)-Yglobal(i)) .lt. deltaxyz(2)*0.0001) then 
           xyzmin_idx(2) = i-1 ! -1 offset because binfile_chunklet() is a C++ function
        end if
     end do
     do i=1,cctk_gsh(3)
        if(abs(Z(1,1,1)-Zglobal(i)) .lt. deltaxyz(3)*0.0001) then 
           xyzmin_idx(3) = i-1 ! -1 offset because binfile_chunklet() is a C++ function
        end if
     end do
     if(xyzmin_idx(1).eq.2000) then
        write(*,*) "INSIDE dump_binfile_chunklet: GRID ERROR..."; stop
     end if

     xyzmin_coord(1) = X(1,1,1)
     xyzmin_coord(2) = Y(1,1,1)
     xyzmin_coord(3) = Z(1,1,1)

     call CCTK_ReductionHandle(handle,"minimum")
     call CCTK_VarIndex(index,"grid::X")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xyzmin_glob(1),1,index)
     call CCTK_VarIndex(index,"grid::Y")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xyzmin_glob(2),1,index)
     call CCTK_VarIndex(index,"grid::Z")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xyzmin_glob(3),1,index)

     !=====================================================================================
     ! Now dump binfile chunklets, but only write N*chunklet_procs_at_a_time chunklets
     !    at a time, where N is the number of gridfunctions output at a time.
     ! This is necessary so that for large runs, we don't overflow the filesystem with a 
     !    bajillion file output requests at once.

     num_procs = CCTK_nProcs(cctkGH)
     myproc = CCTK_MyProc(cctkGH)
     do whichproc=0,num_procs,chunklet_procs_at_a_time

        ! The following function call forces all processors to wait while the file-writing processors are writing to disk.
        call CCTK_ReduceLocalScalar(ierr,cctkGH,-1,handle,dummy,dummy,CCTK_VARIABLE_INT)

        write(*,*) "HI",whichproc,myproc,whichproc+chunklet_procs_at_a_time

        if(myproc.ge.whichproc .and. myproc.le.whichproc+chunklet_procs_at_a_time) then

           ! Dump bin files
           open(19,file='DUMP_Input',status='old',ERR=110)
           if(1==0) then
110           write(*,*) "ERROR OPENING DUMP_Input!"; stop
           end if
           do 
              read(19,*,IOSTAT=ierr) fname 
              if (ierr .ne. 0) exit

              ! This is ugly. Need to fix it later.
              if (fname=='phi') then 
                 call bhns_outbinfile_chunklet(cctkGH,fname,cctk_lsh,xyzmin_idx,xyzmin_coord,cctk_time,cctk_gsh,xyzmin_glob,deltaxyz,cctk_nghostzones,Symm,phi,iter)
              elseif (fname=='rho_b') then
                 call bhns_outbinfile_chunklet(cctkGH,fname,cctk_lsh,xyzmin_idx,xyzmin_coord,cctk_time,cctk_gsh,xyzmin_glob,deltaxyz,cctk_nghostzones,Symm,rho_b,iter)
              elseif (fname=='rho_star') then
                 call bhns_outbinfile_chunklet(cctkGH,fname,cctk_lsh,xyzmin_idx,xyzmin_coord,cctk_time,cctk_gsh,xyzmin_glob,deltaxyz,cctk_nghostzones,Symm,rho_star,iter)
              elseif (fname=='P') then
                 call bhns_outbinfile_chunklet(cctkGH,fname,cctk_lsh,xyzmin_idx,xyzmin_coord,cctk_time,cctk_gsh,xyzmin_glob,deltaxyz,cctk_nghostzones,Symm,P,iter)
              elseif (fname=='gxx') then
                 call bhns_outbinfile_chunklet(cctkGH,fname,cctk_lsh,xyzmin_idx,xyzmin_coord,cctk_time,cctk_gsh,xyzmin_glob,deltaxyz,cctk_nghostzones,Symm,gxx,iter)
              elseif (fname=='vx') then
                 call bhns_outbinfile_chunklet(cctkGH,fname,cctk_lsh,xyzmin_idx,xyzmin_coord,cctk_time,cctk_gsh,xyzmin_glob,deltaxyz,cctk_nghostzones,Symm,vx,iter)
              elseif (fname=='vy') then
                 call bhns_outbinfile_chunklet(cctkGH,fname,cctk_lsh,xyzmin_idx,xyzmin_coord,cctk_time,cctk_gsh,xyzmin_glob,deltaxyz,cctk_nghostzones,Symm,vy,iter)
              elseif (fname=='vz') then
                 call bhns_outbinfile_chunklet(cctkGH,fname,cctk_lsh,xyzmin_idx,xyzmin_coord,cctk_time,cctk_gsh,xyzmin_glob,deltaxyz,cctk_nghostzones,Symm,vz,iter)
              elseif (fname=='lapm1') then
                 call bhns_outbinfile_chunklet(cctkGH,fname,cctk_lsh,xyzmin_idx,xyzmin_coord,cctk_time,cctk_gsh,xyzmin_glob,deltaxyz,cctk_nghostzones,Symm,lapm1,iter)
              elseif (fname=='Bx') then
                 call bhns_outbinfile_chunklet(cctkGH,fname,cctk_lsh,xyzmin_idx,xyzmin_coord,cctk_time,cctk_gsh,xyzmin_glob,deltaxyz,cctk_nghostzones,Symm,Bx,iter)
              elseif (fname=='By') then
                 call bhns_outbinfile_chunklet(cctkGH,fname,cctk_lsh,xyzmin_idx,xyzmin_coord,cctk_time,cctk_gsh,xyzmin_glob,deltaxyz,cctk_nghostzones,Symm,By,iter)
              elseif (fname=='Bz') then
                 call bhns_outbinfile_chunklet(cctkGH,fname,cctk_lsh,xyzmin_idx,xyzmin_coord,cctk_time,cctk_gsh,xyzmin_glob,deltaxyz,cctk_nghostzones,Symm,Bz,iter)
              else 
                 write(*,*) "Grid function ",fname," can't be dumped..." 
                 write(*,*) "Either the grid function doesn't exist or it is not "
                 write(*,*) "recognized in dump_binfile_chunklet."
              end if

           end do
           close(19)
        end if
     end do
  end if

end subroutine bhns_dump_binfile_chunklet
