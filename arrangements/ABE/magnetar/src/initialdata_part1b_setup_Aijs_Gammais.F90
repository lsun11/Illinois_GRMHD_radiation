!--------------------------------------------------------
! Okay, we've read in the initial data from files.
! Now we set up all other required variables, including:
!  emfields, BSSN variables, primitives, etc.
!--------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine magnetar_setup_Aijs_Gammais(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  real*8                :: detmin_l,detmax_l
  integer               :: i,j,k,n1,n2,n3,mf
  integer               :: ierr,index,handle,dummy,glob_imax,glob_jmax,glob_kmax
  CCTK_REAL             :: reduction_value

  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  integer :: ONE,ZERO
  parameter(ONE = 1.D0, ZERO = 0.D0)

  ext = cctk_lsh

  write(*,*) "Setting up Aij's, Gammai's"

  ! Note:  kset_c uses *physical* metric as input, and outputs A_{ij} using (d/dt)\gamma_{ij} = 0
  !we use all the _bck's and _fwd's as temporary storage, as well as Pl.
  call kset_c_v2(ext,X,Y,Z, &
       Axx,Axy,Axz,Ayy,Ayz,Azz,trK, &
       gxx,gxy,gxz,gyy,gyz,gzz, &
       shiftx,shifty,shiftz,lapm1, &
       Symmetry, &
       temp1,temp2,temp3,temp4,temp5, &
       temp6,temp7,temp8,temp9,temp10, &
       gxxx,gxxy,gxxz,gxyx,gxyy,gxyz, &
       gxzx,gxzy,gxzz,gyyx,gyyy,gyyz, &
       gyzx,gyzy,gyzz,gzzx,gzzy,gzzz, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
       Pl)

  !  call kset_c_cpp(cctkGH,cctk_lsh,cctk_nghostzones,Symmetry, &
  !       dx,dy,dz, &
  !       Axx,Axy,Axz,Ayy,Ayz,Azz,trK,  &
  !       gxx,gxy,gxz,gyy,gyz,gzz,  &
  !       shiftx,shifty,shiftz,lapm1)
  !  call kset_c(ext,X,Y,Z, &
  !         Axx,Axy,Axz,Ayy,Ayz,Azz,trK, &
  !         gxx,gxy,gxz,gyy,gyz,gzz, &
  !         shiftx,shifty,shiftz,lapm1, &
  !         Symmetry)

  trK = 0.D0

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  if(Symmetry .eq. AXISYM) then 
     call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')
     call CCTK_VarIndex(index,'BSSN::Axx')
     call BndCartoon2DVI(dummy, cctkGH, 2, index)
  end if

  detmin_l = 0.D0
  detmax_l = 0.D0

  !====================================================
  ! Convert to tilde metric and invert tilde metric 
  !====================================================
  !  call convert_cpp(cctkGH,cctk_lsh,phi, &
  !       gxx,gxy,gxz,gyy,gyz,gzz, &
  !       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)


  call convert(ext,phi,gxx,gxy,gxz,gyy,gyz,gzz, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,gxx_p, &
       detmin_l,detmax_l)
  gxx_p = gxx

  !====================================================
  ! Get Gamma^i
  !====================================================
  !  call setgamma_cpp(cctkGH,cctk_lsh, &
  !       dx,dy,dz, &
  !       gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
  !       Gammax, Gammay, Gammaz)
  ! Use gxx_p as temporary storage:
  call setgamma(ext,X,Y,Z,gxx_p, &
       phi,gxx,gxy,gxz,gyy,gyz,gzz, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
       Gammax,Gammay,Gammaz, &
       Symmetry)
  gxx_p = gxx

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  if(Symmetry .eq. AXISYM) then 
     call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')
     call CCTK_VarIndex(index,'BSSN::Gammax')
     call BndCartoon2DVI(dummy, cctkGH, 1, index)
  end if

  n1 = 0
  n2 = 0
  n3 = 0
  mf = 1
  !use gxx_p as temporary storage here:
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammax,gxx_p,n1,n2,n3,mf,Symmetry)
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammay,gxx_p,n1,n2,n3,mf,Symmetry)
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammaz,gxx_p,n1,n2,n3,mf,Symmetry)
  gxx_p=gxx

  !Scale Aij by appropriate conformal factor to make it BSSN-compatible.
  Axx = Axx*exp(-4.D0*phi)
  Axy = Axy*exp(-4.D0*phi)
  Axz = Axz*exp(-4.D0*phi)
  Ayy = Ayy*exp(-4.D0*phi)
  Ayz = Ayz*exp(-4.D0*phi)
  Azz = Azz*exp(-4.D0*phi)

  !Set psi
  psi = exp(phi)

end subroutine magnetar_setup_Aijs_Gammais
