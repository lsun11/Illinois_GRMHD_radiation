#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine stationary_puncture_initial_data_part1(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,global_ext
  real*8 			           :: dT,dX,dY,dZ

  real*8                                   :: phib,lapm1b,rb,PI,actual_excis_radius,excis_radius_fish
  real*8                                   :: shiftxb,shiftyb,shiftzb,shiftrb
  real*8                                   :: gxxb,gxyb,gxzb,gyyb,gyzb,gzzb,Gammaxb,Gammayb,Gammazb
  real*8                                   :: Axxb,Axyb,Axzb,Ayyb,Ayzb,Azzb
  real*8                                   :: gupxxb,gupxyb,gupxzb,gupyyb,gupyzb,gupzzb
  real*8, dimension(1,3)                   :: pointcoords
  ! integer                                  :: N_theta_sm,N_phi_sm,ntot_sm,numtheta_sm,numphi_sm
  !  real*8                                   :: dphi_sm,dcostheta_sm

  integer :: index,dummy,foundit_flag
  integer :: handle,ierr
  real*8  :: ONE,ZERO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(ONE = 1.D0, ZERO = 0.D0)

  INTEGER :: i,j,k !,iter,iterloop,ITERMAX=150,lo_or_hi
  !  REAL*8 :: xlo,xhi,r_iso,outputdiff,outputdiffold,mass,r_s,betar,dbdr,drdbr
  !  REAL*8 :: EPS=1.D-14

  ext = cctk_lsh

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  PI = 3.14159265358979323846D0

  if(asymptotic_bh==0) then
     call stationary_puncture_id_nonasymptoticisotropic(CCTK_PASS_FTOF)
!  else if(asymptotic_bh==1) then
!     call stationary_puncture_id_asymptoticBaumgarte(CCTK_PASS_FTOF)
  else if(asymptotic_bh==2) then
     call stationary_puncture_id_nonKeq0_soln_set_psi_gauge(cctkGH,cctk_nghostzones,cctk_lsh, &
          x,y,z, PhysicalRadius,r, &
          lapm1, shiftx, shifty, shiftz, Sz, psi)
     call stationary_puncture_id_nonKeq0_soln_finalize(CCTK_PASS_FTOF)
  end if
  if(asymptotic_bh==2 .and. fill_excision_enable.ne.0) then
     write(*,*) "Error: asymptotic_bh==2 automatically fills the excised BH with junk data!"
     stop
  end if
  if(fill_excision_enable.ne.0) then
! TODO: FIXME: 
!     call fill_excised_region(CCTK_PASS_FTOF)
  end if

  ! Set excision_zone_gf to avoid valgrind memory errors
  excision_zone_gf = 0

end subroutine stationary_puncture_initial_data_part1
