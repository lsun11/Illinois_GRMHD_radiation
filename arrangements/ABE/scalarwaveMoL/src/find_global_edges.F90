!-----------------------------------------------------------------------------
! Update outer boundaries
!-----------------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine scalarwaveMoL_find_global_edges(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_FUNCTIONS;
  
  integer, dimension(3)                    :: ext,fake_ext
  real*8                                   :: dT,dx,dy,dz,Xglobmin,Yglobmin,Zglobmin,Xglobmax,Yglobmax,Zglobmax
  integer                                  :: i,AXISYM,ierr,handle,vindex

  AXISYM = 4

  ext = cctk_lsh

  dT = CCTK_DELTA_TIME
  dx = CCTK_DELTA_SPACE(1)
  dy = CCTK_DELTA_SPACE(2)
  dz = CCTK_DELTA_SPACE(3)

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
  call scalarwaveMol_global_edges(dX,dY,dZ, &
       Xglobmin,Yglobmin,Zglobmin, &
       X(1,1,1),Y(1,1,1),Z(1,1,1), &
       Xglobmax,Yglobmax,Zglobmax, &
       X(cctk_lsh(1),1,1),Y(1,cctk_lsh(2),1),Z(1,1,cctk_lsh(3)), &
       scawave_have_bdry_min,scawave_have_bdry_max,scalarwave_Symmetry)

  write(*,*) "GLOBAL EDGES:",Xglobmin,Yglobmin,Zglobmin,Xglobmax,Yglobmax,Zglobmax
  write(*,*) "GLOBAL EDGES:",scawave_have_bdry_min,scawave_have_bdry_max

end subroutine scalarwaveMoL_find_global_edges
