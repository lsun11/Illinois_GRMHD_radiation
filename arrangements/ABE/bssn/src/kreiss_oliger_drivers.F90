!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!
!-----------------------------------------------------------------------------
!
! Apply Kreiss-Oliger dissipation to Aij & trK.  Placing this file into the
!   excision thorn will cause problems since it uses Aij_rhs, private bssn
!   variables
!
!-----------------------------------------------------------------------------
subroutine kreiss_oliger_timestepping(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8 				:: dT
  integer, dimension(3)			:: ext
  integer            			:: dummy
  real*8, parameter                     :: SYM = 1.D0, ANT= - 1.D0
  integer              :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM,AXISYM
  parameter (NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  ext = cctk_lsh
  dT = CCTK_DELTA_TIME
  
  if(Symmetry==AXISYM) then
     !here we use Kij and psi as dummy variables.  
     ! Don't panic. Kij and psi are overwritten at the end of each predictor/corrector.
     call flat_lap(ext,X,Z,Axx,Kxx,SYM,SYM,Symmetry)
     call flat_lap(ext,X,Z,Axy,Kxy,SYM,SYM,Symmetry)
     call flat_lap(ext,X,Z,Axz,Kxz,ANT,ANT,Symmetry)
     call flat_lap(ext,X,Z,Ayy,Kyy,SYM,SYM,Symmetry)
     call flat_lap(ext,X,Z,Ayz,Kyz,ANT,ANT,Symmetry)
     call flat_lap(ext,X,Z,Azz,Kzz,SYM,SYM,Symmetry)
     call flat_lap(ext,X,Z,trK,psi,SYM,SYM,Symmetry)

     ! apply copy excision boundary conditions to the laplacians
     call scalar_excision_bc(ext,X,Y,Z,psi,Symmetry,excision_zone_gf)
     call vector_excision_bc(ext,X,Y,Z,Kxx,Kxy,Kxz, &
          Symmetry,excision_zone_gf)
     call vector_excision_bc(ext,X,Y,Z,Kyy,Kyz,Kzz, &
          Symmetry,excision_zone_gf)

     !FIXME (next two lines, for Carpet compatibility):
     call CartSymGN(dummy,cctkGH,'BSSN::BSSN_AH')
     call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_AH')

     call flat_lap(ext,X,Z,Kxx,gxxy,SYM,SYM,Symmetry)
     call flat_lap(ext,X,Z,Kxy,gxyy,SYM,SYM,Symmetry)
     call flat_lap(ext,X,Z,Kxz,gxzy,ANT,ANT,Symmetry)
     call flat_lap(ext,X,Z,Kyy,gyyy,SYM,SYM,Symmetry)
     call flat_lap(ext,X,Z,Kyz,gyzy,ANT,ANT,Symmetry)
     call flat_lap(ext,X,Z,Kzz,gzzy,SYM,SYM,Symmetry)
     call flat_lap(ext,X,Z,psi,gxyx,SYM,SYM,Symmetry)
      
     ! apply copy excision boundary conditions to the laplacians
     call scalar_excision_bc(ext,X,Y,Z,gxyx,Symmetry,excision_zone_gf)
     call vector_excision_bc(ext,X,Y,Z,gxxy,gxyy,gxzy, &
          Symmetry,excision_zone_gf)
     call vector_excision_bc(ext,X,Y,Z,gyyy,gyzy,gzzy, &
          Symmetry,excision_zone_gf)

     call kreiss_oliger3d(ext,X,Y,Z,Axx_rhs,gxxy,dT,C_ko)
     call kreiss_oliger3d(ext,X,Y,Z,Axy_rhs,gxyy,dT,C_ko)
     call kreiss_oliger3d(ext,X,Y,Z,Axz_rhs,gxzy,dT,C_ko)
     call kreiss_oliger3d(ext,X,Y,Z,Ayy_rhs,gyyy,dT,C_ko)
     call kreiss_oliger3d(ext,X,Y,Z,Ayz_rhs,gyzy,dT,C_ko)
     call kreiss_oliger3d(ext,X,Y,Z,Azz_rhs,gzzy,dT,C_ko)
     call kreiss_oliger3d(ext,X,Y,Z,trK_rhs,gxyx,dT,C_ko)
  else if(Symmetry==EQUATORIAL) then
     !here we use Kij and psi as dummy variables.  
     ! Don't panic. Kij and psi are overwritten at the end of each predictor/corrector.

     !When Symmetry==EQUATORIAL, first 2 symmetry arguments are IGNORED!
     call flat_lap3d(ext,X,Y,Z,Axx,Kxx,SYM,SYM,SYM,Symmetry)
     call flat_lap3d(ext,X,Y,Z,Axy,Kxy,SYM,SYM,SYM,Symmetry)
     call flat_lap3d(ext,X,Y,Z,Axz,Kxz,SYM,SYM,ANT,Symmetry)
     call flat_lap3d(ext,X,Y,Z,Ayy,Kyy,SYM,SYM,SYM,Symmetry)
     call flat_lap3d(ext,X,Y,Z,Ayz,Kyz,SYM,SYM,ANT,Symmetry)
     call flat_lap3d(ext,X,Y,Z,Azz,Kzz,SYM,SYM,SYM,Symmetry)
     call flat_lap3d(ext,X,Y,Z,trK,psi,SYM,SYM,SYM,Symmetry)

     ! apply copy excision boundary conditions to the laplacians
     call scalar_excision_bc(ext,X,Y,Z,psi,Symmetry,excision_zone_gf)
     call vector_excision_bc(ext,X,Y,Z,Kxx,Kxy,Kxz, &
          Symmetry,excision_zone_gf)
     call vector_excision_bc(ext,X,Y,Z,Kyy,Kyz,Kzz, &
          Symmetry,excision_zone_gf)

     call CartSymGN(dummy,cctkGH,'BSSN::BSSN_AH')
     call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_AH')

     !When Symmetry==EQUATORIAL, first 2 symmetry arguments are IGNORED!
     call flat_lap3d(ext,X,Y,Z,Kxx,gxxy,SYM,SYM,SYM,Symmetry)
     call flat_lap3d(ext,X,Y,Z,Kxy,gxyy,SYM,SYM,SYM,Symmetry)
     call flat_lap3d(ext,X,Y,Z,Kxz,gxzy,SYM,SYM,ANT,Symmetry)
     call flat_lap3d(ext,X,Y,Z,Kyy,gyyy,SYM,SYM,SYM,Symmetry)
     call flat_lap3d(ext,X,Y,Z,Kyz,gyzy,SYM,SYM,ANT,Symmetry)
     call flat_lap3d(ext,X,Y,Z,Kzz,gzzy,SYM,SYM,SYM,Symmetry)
     call flat_lap3d(ext,X,Y,Z,psi,gxyx,SYM,SYM,SYM,Symmetry)

     ! apply copy excision boundary conditions to the laplacians
     call scalar_excision_bc(ext,X,Y,Z,gxyx,Symmetry,excision_zone_gf)
     call vector_excision_bc(ext,X,Y,Z,gxxy,gxyy,gxzy, &
          Symmetry,excision_zone_gf)
     call vector_excision_bc(ext,X,Y,Z,gyyy,gyzy,gzzy, &
          Symmetry,excision_zone_gf)

     call kreiss_oliger3d(ext,X,Y,Z,Axx_rhs,gxxy,dT,C_ko)
     call kreiss_oliger3d(ext,X,Y,Z,Axy_rhs,gxyy,dT,C_ko)
     call kreiss_oliger3d(ext,X,Y,Z,Axz_rhs,gxzy,dT,C_ko)
     call kreiss_oliger3d(ext,X,Y,Z,Ayy_rhs,gyyy,dT,C_ko)
     call kreiss_oliger3d(ext,X,Y,Z,Ayz_rhs,gyzy,dT,C_ko)
     call kreiss_oliger3d(ext,X,Y,Z,Azz_rhs,gzzy,dT,C_ko)
     call kreiss_oliger3d(ext,X,Y,Z,trK_rhs,gxyx,dT,C_ko)

  end if

end subroutine kreiss_oliger_timestepping
