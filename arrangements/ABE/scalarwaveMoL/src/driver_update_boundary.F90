!-----------------------------------------------------------------------------
! Update outer boundaries
!-----------------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine scalarwaveMoL_update_boundary(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_FUNCTIONS;

  integer, dimension(3)                    :: ext,fake_ext,local_have_bdry_min,local_have_bdry_max, gridfunc_syms
  real*8                                   :: dT,dx,dy,dz,Xglobmin,Yglobmin,Zglobmin,Xglobmax,Yglobmax,Zglobmax,levelnumber
  integer                                  :: i,AXISYM,ierr,handle,vindex
  integer 				   :: stagger_x,stagger_y,stagger_z

  AXISYM = 4

  ext = cctk_lsh

  dT = CCTK_DELTA_TIME
  dx = CCTK_DELTA_SPACE(1)
  dy = CCTK_DELTA_SPACE(2)
  dz = CCTK_DELTA_SPACE(3)

  levelnumber = cctk_levfac(1)
  levelnumber = log(levelnumber)/log(2.D0)+1.D0

  if(CCTK_nProcs(CCTKGH)+1.gt.16384) then
     write(*,*) "Error: Total number of processors exceeds hard-coded limit of 16384 inside update_boundary!"
     write(*,*) "       Please edit scalarwaveMoL/interface.ccl so that scawave_have_global_boundary arrays have > 16384 elements in dimension 1."
     stop
  end if

  !First we fill in symmetry ghostzones
  call fill_scalarwave_symmetry_gz(ext,X,Y,Z,scalarwave_Symmetry,phi,phidot)
  ! *** TEST ***
  stagger_x = 1
  stagger_y = 1
  stagger_z = 0
  gridfunc_syms(1) = 1
  gridfunc_syms(2) = 1
  gridfunc_syms(3) = 1
  !!call fill_scalarwave_symmetry_gz(ext,X,Y,Z,scalarwave_Symmetry,phi_stagger,phidot_stagger)
  call scalarwave_do_symmetry_gz_one_gridfunction_stagger(ext,X,Y,Z,Symmetry,phi_stagger, & 
        					gridfunc_syms,stagger_x,stagger_y,stagger_z)
  call scalarwave_do_symmetry_gz_one_gridfunction_stagger(ext,X,Y,Z,Symmetry,phidot_stagger, &
                                                gridfunc_syms,stagger_x,stagger_y,stagger_z)
  ! ************

  if(cctk_iteration.eq.0 .and. levelnumber==1) then
     !Next we find the global coordinate maxima/minima on the grid
     call CCTK_ReductionHandle(handle,"minimum")
     call CCTK_VarIndex(vindex,"grid::X")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,Xglobmin,1,vindex)
     call CCTK_VarIndex(vindex,"grid::Y")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,Yglobmin,1,vindex)
     call CCTK_VarIndex(vindex,"grid::Z")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,Zglobmin,1,vindex)

     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex(vindex,"grid::X")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,Xglobmax,1,vindex)
     call CCTK_VarIndex(vindex,"grid::Y")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,Yglobmax,1,vindex)
     call CCTK_VarIndex(vindex,"grid::Z")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,Zglobmax,1,vindex)

     !Then we determine whether any of the faces of our local grid corresponds to a global max/min:
     ! This way, we only update actual outer boundaries (and not refinement or symmetry boundaries)x
     call scalarwaveMol_global_edges(dX,dY,dZ, &
          Xglobmin,Yglobmin,Zglobmin, &
          X(1,1,1),Y(1,1,1),Z(1,1,1), &
          Xglobmax,Yglobmax,Zglobmax, &
          X(cctk_lsh(1),1,1),Y(1,cctk_lsh(2),1),Z(1,1,cctk_lsh(3)), &
          local_have_bdry_min,local_have_bdry_max,scalarwave_Symmetry)

     scawave_have_bdry_min(CCTK_MyProc(CCTKGH)+1,:) = local_have_bdry_min
     scawave_have_bdry_max(CCTK_MyProc(CCTKGH)+1,:) = local_have_bdry_max

     write(*,*) "Just set GLOBAL EDGES for the first time: levnum, blah=",levelnumber,local_have_bdry_min,local_have_bdry_max
  end if

  if(cctk_iteration.gt.0 .and. levelnumber==1) then
     local_have_bdry_min = scawave_have_bdry_min(CCTK_MyProc(CCTKGH)+1,:)
     local_have_bdry_max = scawave_have_bdry_max(CCTK_MyProc(CCTKGH)+1,:)

     write(*,*) "Update boundary: levnum, GLOBAL EDGES:",levelnumber,local_have_bdry_min,local_have_bdry_max
     
     !Then we update boundaries
     do i=1,cctk_nghostzones(2)
        fake_ext = cctk_lsh - cctk_nghostzones + i
        if(scalarwave_Symmetry==AXISYM) fake_ext(1) = fake_ext(1) + 1
	! *** TEST ***
        call scalarwaveMoL_update_boundary_lowlevel(ext, fake_ext, X, Y, Z, &
             dX,dY,dZ,dT, &
             phi_p,phi,phidot_p,phidot, phi_stagger_p,phi_stagger, & 
             phidot_stagger_p, phidot_stagger, scalarwave_Symmetry, &
             local_have_bdry_min,local_have_bdry_max)
        ! *************
     end do
  end if

end subroutine scalarwaveMoL_update_boundary
