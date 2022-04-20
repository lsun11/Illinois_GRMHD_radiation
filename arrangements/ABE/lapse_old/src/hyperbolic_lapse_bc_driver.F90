!---------------------------------------------
! Update boundary driver for hyperbolic lapse
!---------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine hyperbolic_lapse_bc_driver(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS   

  interface
     subroutine hyper_lapse_bc(ext, fake_ext, dT, x, y, z, PhysicalRadius, RadiusDerivative,&
          lapse_new, dtlapse_new, &
          lapse, dtlapse, phi, &
          Xcenter,Ycenter,Zcenter,Symmetry, &
          gxx, gxy, gxz, gyy, gyz, gzz, &
          thisproc_have_bdry_min,thisproc_have_bdry_max)
       integer, dimension(3)                   :: ext,fake_ext,thisproc_have_bdry_min,thisproc_have_bdry_max
       real*8                                  :: dT,Xcenter,Ycenter,Zcenter
       real*8, dimension(ext(1),ext(2),ext(3)) :: PhysicalRadius, RadiusDerivative
       real*8, dimension(ext(1),ext(2),ext(3)) :: X,Y,Z
       real*8, dimension(ext(1),ext(2),ext(3)) :: lapse_new,dtlapse_new,lapse, dtlapse, phi
       real*8, dimension(ext(1),ext(2),ext(3)) :: gxx, gxy, gxz, gyy, gyz, gzz
       integer                                 :: Symmetry
     end subroutine hyper_lapse_bc
  end interface
  integer                                   :: n1,n2,n3,m,i,dummy
  integer, dimension(3)                     :: ext,fake_ext,thisproc_have_global_bdry_min,thisproc_have_global_bdry_max
  real*8                                    :: dT,dX,levelnumber
  real*8                                    :: HALF,ld_eps,ld_c
  integer                    :: PI_SYMM, AXISYM
  parameter(PI_SYMM = 3, AXISYM = 4)
  parameter ( HALF = 0.5D0 )
  !
  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dT = CCTK_DELTA_TIME

  levelnumber = cctk_levfac(1)
  levelnumber = log(levelnumber)/log(2.D0)+1.D0

  if(cctk_iteration.gt.0 .and. levelnumber==1.D0) then
     thisproc_have_global_bdry_min = have_global_bdry_min(int(levelnumber),:)
     thisproc_have_global_bdry_max = have_global_bdry_max(int(levelnumber),:)
     do i=1,cctk_nghostzones(2)
        fake_ext = cctk_lsh - cctk_nghostzones + i
        if(Symmetry==AXISYM) fake_ext(1) = fake_ext(1) + 1
        call hyper_lapse_bc(ext,fake_ext,dT,X,Y,Z,PhysicalRadius, RadiusDerivative,lapm1,lapset, &
             lapm1_p,lapset_p,phi_p, &
             Xcenter,Ycenter,Zcenter,Symmetry, &
             gxx_p, gxy_p, gxz_p, gyy_p, gyz_p, gzz_p, &
             thisproc_have_global_bdry_min,thisproc_have_global_bdry_max)
     end do

     if(Symmetry == AXISYM) then
        !first rotate, then compute values at the boundary. 
        !(needed since BndCartoon2D can only(?) be called inside a (CCTK_ARGUMENTS) function)
        call BndCartoon2DVN(dummy, cctkGH, 0, 'lapse::lapm1')
     end if

  end if

  !  call CCTK_SyncGroup(dummy,cctkGH,'lapse::lapse_vars')
  !  call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')

end subroutine hyperbolic_lapse_bc_driver
