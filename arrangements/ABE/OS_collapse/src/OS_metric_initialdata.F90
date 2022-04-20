#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine OS_metric_initialdata(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,global_ext
  real*8 			           :: dT,dX,dY,dZ
  real*8 			           :: n1,n2,n3,mf

  real*8                                   :: PI,dmass,detmin,detmax

  integer :: index,dummy,foundit_flag
  integer :: handle,ierr
  real*8  :: ONE,ZERO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(ONE = 1.D0, ZERO = 0.D0)

  INTEGER :: i,j,k

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  PI = 3.14159265358979323846D0

  ext = cctk_lsh

  gupxx = 1.D0
  gupxy = 0.D0
  gupxz = 0.D0
  gupyy = 1.D0
  gupyz = 0.D0
  gupzz = 1.D0

  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
	   Kxx(i,j,k) = 0.d0
           Kxx(i,j,k) = 0.d0

           Kxy(i,j,k) = 0.d0
           Kxy(i,j,k) = 0.d0

           Kxz(i,j,k) = 0.d0
           Kxz(i,j,k) = 0.d0

           Kyy(i,j,k) = 0.d0
           Kyy(i,j,k) = 0.d0

           Kyz(i,j,k) = 0.d0
           Kyz(i,j,k) = 0.d0

           Kzz(i,j,k) = 0.d0
           Kzz(i,j,k) = 0.d0

        end do
     end do
  end do

  !  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_AH')
  !  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_AH')

  ! next compute phi:
  phi = log(psi)

  gxx = 1.D0
  gxy = 0.D0
  gxz = 0.D0
  gyy = 1.D0
  gyz = 0.D0
  gzz = 1.D0

  shiftx = 0.D0
  shifty = 0.D0
  shiftz = 0.D0

  ! Set excision_zone_gf to avoid valgrind memory errors
  excision_zone_gf = 0

  !Compute Atilde_ij's
  Axx = 0.d0
  Axy = 0.d0
  Axz = 0.d0
  Ayy = 0.d0
  Ayz = 0.d0
  Azz = 0.d0

  shiftxt = 0.D0
  shiftyt = 0.D0
  shiftzt = 0.D0
  lapset = 0.D0

  trK = 0.D0

  Gammax = 0.D0
  Gammay = 0.D0
  Gammaz = 0.D0

  !Look, Gamma^i was just set to zero,
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
!  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')

  !Set derivatives to zero
  gxxx = 0.D0
  gxxy = 0.D0 
  gxxz = 0.D0
  gxyx = 0.D0 
  gxyy = 0.D0 
  gxyz = 0.D0
  gxzx = 0.D0 
  gxzy = 0.D0 
  gxzz = 0.D0
  gyyx = 0.D0 
  gyyy = 0.D0 
  gyyz = 0.D0
  gyzx = 0.D0 
  gyzy = 0.D0 
  gyzz = 0.D0
  gzzx = 0.D0 
  gzzy = 0.D0 
  gzzz = 0.D0 

  Gammaxxx = 0.D0 
  Gammaxxy = 0.D0 
  Gammaxxz = 0.D0 
  Gammaxyy = 0.D0 
  Gammaxyz = 0.D0 
  Gammaxzz = 0.D0
  Gammayxx = 0.D0 
  Gammayxy = 0.D0 
  Gammayxz = 0.D0 
  Gammayyy = 0.D0 
  Gammayyz = 0.D0 
  Gammayzz = 0.D0
  Gammazxx = 0.D0 
  Gammazxy = 0.D0 
  Gammazxz = 0.D0 
  Gammazyy = 0.D0 
  Gammazyz = 0.D0 
  Gammazzz = 0.D0


  ! Set BSSN matter sources to zero!
  rho = ZERO
  S = ZERO
  Sx = ZERO
  Sy = ZERO
  Sz = ZERO
  Sxx = ZERO
  Sxy = ZERO
  Sxz = ZERO
  Syy = ZERO
  Syz = ZERO
  Szz = ZERO

  ! call find_apparent_horizon
  !  adm::AH_Manager.Find_Horizon(level,0.0,0.0);

  if(fisheye_enable .eq. 1) then
    write(*,*) 'Fisheye not supported in TwoPunctures! Bye...' 
    stop
  end if

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')
  call CartSymGN(dummy,cctkGH,'shift::shift_vars')

  if(Symmetry==AXISYM) then
     write(*,*) "I PITY DA FOO TRIES TO USE AXISYMMETRY FOR ORBITING BLACK HOLES!"
     stop
  end if

  !Following lines are (empirically speaking) ESSENTIAL for excision runs:
  call Derivs(ext,X,Y,Z,dX,dY,dZ,phi,phix,phiy,phiz,Symmetry) 
  call Derivs(ext,X,Y,Z,dX,dY,dZ,lapm1,lapsex,lapsey,lapsez,Symmetry) 

  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,phi,phix,phiy,phiz)
  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,lapm1,lapsex,lapsey,lapsez)

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  !  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars') 

end subroutine OS_metric_initialdata
