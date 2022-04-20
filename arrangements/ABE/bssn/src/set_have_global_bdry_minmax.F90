#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bssn_set_have_global_bdry_minmax(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext,thisproc_have_global_bdry_min,thisproc_have_global_bdry_max
  real*8                :: dX,dY,dZ,levelnumber
  real*8                :: Xglobmin,Yglobmin,Zglobmin,Xglobmax,Yglobmax,Zglobmax
  integer               :: dummy,vindex,handle,ierr

  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  !Here, we determine whether any of the faces of our local grid corresponds to a global max/min:
  ! This way, we only update actual outer boundaries (and not refinement or symmetry boundaries)x
  levelnumber = cctk_levfac(1)
  levelnumber = log(levelnumber)/log(2.D0)+1.D0

  have_global_bdry_min(int(levelnumber),:) = 0
  have_global_bdry_max(int(levelnumber),:) = 0

  if(levelnumber==1.D0) then
     have_global_bdry_min = 0
     have_global_bdry_max = 0

     !First we find the global maximum/minimum of our grid:
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

     call global_edges(dX,dY,dZ, &
          Xglobmin,Yglobmin,Zglobmin, &
          X(1,1,1),Y(1,1,1),Z(1,1,1), &
          Xglobmax,Yglobmax,Zglobmax, &
          X(cctk_lsh(1),1,1),Y(1,cctk_lsh(2),1),Z(1,1,cctk_lsh(3)), &
          thisproc_have_global_bdry_min,thisproc_have_global_bdry_max,Symmetry)

     write(*,*) "blah",Xglobmax,X(cctk_lsh(1),1,1)

     have_global_bdry_min(levelnumber,:) = thisproc_have_global_bdry_min
     have_global_bdry_max(levelnumber,:) = thisproc_have_global_bdry_max

     write(*,*) "SETTING thisproc_have_global_bdry_minmaxs:",thisproc_have_global_bdry_min,thisproc_have_global_bdry_max
  end if

  write(*,*) "BSSN MIN?:",levelnumber,have_global_bdry_min(1,:)

end subroutine bssn_set_have_global_bdry_minmax
