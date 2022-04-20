!! check this routine
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine OS_initialdata_local(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext

  real*8 			           :: dT,dX,dY,dZ
  real*8 			           :: n1,n2,n3,mf

  integer :: vindex,handle,ierr
  real*8  :: rho_max,tau_max
  real*8  :: ONE,ZERO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(ONE = 1.D0, ZERO = 0.D0)

  INTEGER :: i,j,k
  character :: filename*30,c2*2,c1
  CCTK_REAL :: reduction_value

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ext = cctk_lsh

  !initialize all this stuff
  P=ZERO
  rho_star=ZERO
  h=ZERO
  w=ZERO
  tau=ZERO
  st_x=ZERO
  st_y=ZERO
  st_z=ZERO

  sbt=ZERO
  sbx=ZERO
  sby=ZERO
  sbz=ZERO
  Bx=ZERO
  By=ZERO
  Bz=ZERO
  Ex=ZERO
  Ey=ZERO
  Ez=ZERO
  mhd_st_x=ZERO
  mhd_st_y=ZERO
  mhd_st_z=ZERO

  ! Set excision_zone_gf to avoid valgrind memory errors
  excision_zone_gf = 0
 
 
  !what is this?
  shiftxt = 0.D0
  shiftyt = 0.D0
  shiftzt = 0.D0
  lapset = 0.D0


  !Look, Gamma^i was set to zero,
  ! so next lines not strictly needed for brandt-brugmann ID, 
  ! but doesn't hurt in case this code is ever recycled!
  n1 = 0
  n2 = 0
  n3 = 0
  mf = 1
  !use PsiTau as temporary storage here:
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammax,PsiTau,n1,n2,n3,mf,Symmetry)
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammay,PsiTau,n1,n2,n3,mf,Symmetry)
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammaz,PsiTau,n1,n2,n3,mf,Symmetry)

  !Next line NEEDED!
!  call CCTK_SyncGroup(ierr,cctkGH,'BSSN::BSSN_vars')


  if(Symmetry==AXISYM) then
     write(*,*) "I PITY DA FOO TRIES TO USE AXISYMMETRY FOR ORBITING BLACK HOLES!"
     stop
  end if

end subroutine OS_initialdata_local
