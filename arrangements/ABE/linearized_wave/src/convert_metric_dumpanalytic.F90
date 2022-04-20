!-----------------------------------------------------------------------------
! Linearized wave driver routines
!-----------------------------------------------------------------------------

!======================================================
! Convert metric to tilded metric, initialize \Gamma^a
!======================================================
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "GenericFD.h"

subroutine linwave_Convert(CCTK_ARGUMENTS)
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  integer                    :: dummy,local_spatial_order
  CCTK_REAL                     :: n1,n2,n3,mf,ONE, TWO, HALF, detmin_l, detmax_l
  CCTK_REAL                     :: dX,dY,dZ
  integer, parameter     	:: AXISYM = 4
  integer, dimension(3)         :: ext
! 1st of 2 needed #includes for GenericFD.h:
#include "../../GenFD_decl_varF90.h"
  parameter(ONE = 1.D0, TWO = 2.D0, HALF = 0.5D0)

  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  !     Initialise finite differencing variables
  ! 2nd of 2 needed #includes for GenericFD.h:
#include "../../GenFD_set_varF90.h"

  !===================================================================
  ! convert to tilded metric,initialize \phi, and invert tilde metric 
  !===================================================================
  detmin_l=0.D0
  detmax_l=0.D0

  ! You'll find the source code for convert() in arrangements/ABE/bssn/src/adm_to_bssn.F90
  !  Note that we use PsiTau as temporary storage here.
  call convert(ext,phi,gxx,gxy,gxz,gyy,gyz,gzz, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,PsiTau, &
       detmin_l,detmax_l)

  !
  ! initialize \Gamma^a
  !
  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
           local_spatial_order=4
           if(i==1.or.j==1.or.k==1.or.i==cctk_lsh(1).or.j==cctk_lsh(2).or.k==cctk_lsh(3)) then
              local_spatial_order=0
           else if(i==2.or.j==2.or.k==2.or.i==cctk_lsh(1)-1.or.j==cctk_lsh(2)-1.or.k==cctk_lsh(3)-1) then
              local_spatial_order=2
           end if

           if(local_spatial_order==4) then
              Gammax(i,j,k) = -(D1_c4(gupxx,i,j,k) + D2_c4(gupxy,i,j,k) + D3_c4(gupxz,i,j,k))
              Gammay(i,j,k) = -(D1_c4(gupxy,i,j,k) + D2_c4(gupyy,i,j,k) + D3_c4(gupyz,i,j,k))
              Gammaz(i,j,k) = -(D1_c4(gupxz,i,j,k) + D2_c4(gupyz,i,j,k) + D3_c4(gupzz,i,j,k))
           else if(local_spatial_order==2) then
              Gammax(i,j,k) = -(D1_c2(gupxx,i,j,k) + D2_c2(gupxy,i,j,k) + D3_c2(gupxz,i,j,k))
              Gammay(i,j,k) = -(D1_c2(gupxy,i,j,k) + D2_c2(gupyy,i,j,k) + D3_c2(gupyz,i,j,k))
              Gammaz(i,j,k) = -(D1_c2(gupxz,i,j,k) + D2_c2(gupyz,i,j,k) + D3_c2(gupzz,i,j,k))
           else
              Gammax(i,j,k) = 0.D0
              Gammay(i,j,k) = 0.D0
              Gammaz(i,j,k) = 0.D0
           end if
        end do
     end do
  end do

!  call setgamma_v2(ext,cctk_nghostzones,dX,dY,dZ, &
!       phi,gxx,gxy,gxz,gyy,gyz,gzz, &
!       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
!       Gammax,Gammay,Gammaz,PsiTau, &
!       Symmetry)

!  call setgamma(ext,X,Y,Z, &
!       phi,gxx,gxy,gxz,gyy,gyz,gzz, &
!       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
!       Gammax,Gammay,Gammaz, &
!       Symmetry)

!  call CCTK_SyncGroup(dummy,cctkGH,"BSSN::BSSN_vars")
  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')

  !Okay, setgamma sets Gamma^i at all points EXCEPT outer boundaries.  Next lines correct this.
  ! FIXME: make this work with arbitrary number of ghostzones.
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

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')   
  call fill_bssn_symmetry_gz_bssn_vars(ext,X,Y,Z,Symmetry,phi,chi,trK,Gammax,Gammay,Gammaz,gxx,gxy,gxz,gyy,gyz,gzz,Axx,Axy,Axz,Ayy,Ayz,Azz)
  if (Symmetry == AxiSym) then
     call CCTK_VarIndex(index,'BSSN::Gammax')
     call BndCartoon2DVI(dummy, cctkGH, 1, index)
  end if

  phi_p = phi

  gxx_p = gxx
  gxy_p = gxy
  gxz_p = gxz
  gyy_p = gyy
  gyz_p = gyz
  gzz_p = gzz

  Gammax_p = Gammax
  Gammay_p = Gammay
  Gammaz_p = Gammaz

  ! Only need to sync Gamma^i?
  !    SyncList(Sync_List,t);
end subroutine linwave_Convert


!======================================================
!
! Dump Analytic solution
!
!======================================================

subroutine linwave_DumpAnalytic(CCTK_ARGUMENTS)
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  integer                    :: dummy
  CCTK_REAL                     :: ONE, TWO, HALF, detmin_l, detmax_l,PhysTime
  integer, dimension(3)         :: ext
  parameter(ONE = 1.D0, TWO = 2.D0, HALF = 0.5D0)
  do dummy=1,3
     ext(dummy) = cctk_lsh(dummy)
  end do

  !==================================================
  ! Find analytic result and store into phi_tilde 
  !     on coarse level!
  !==================================================
  PhysTime=cctk_time;
  call lin_wave_analytic(ext,X,Y,Z,PhysTime,amplitude,width, &
       gxx_anal,gxy_anal,gxz_anal,gyy_anal,gyz_anal,gzz_anal, & 
       PhysicalRadius, RadiusDerivative)

  !====================================================
  ! Invert tilde metric 
  !====================================================
  detmin_l=0.D0
  detmax_l=0.D0
  call convert(ext,phi_anal,gxx_anal,gxy_anal,gxz_anal,gyy_anal,gyz_anal,gzz_anal, &
       gupxx_anal,gupxy_anal,gupxz_anal,gupyy_anal,gupyz_anal,gupzz_anal,gxx_rhs, &
       detmin_l,detmax_l);
  gxx_rhs = 0.D0

  call CCTK_SyncGroup(dummy,cctkGH,"linearized_wave::linearized_wave_analytic");

  call setgamma(ext,X,Y,Z, &
       phi_anal,gxx_anal,gxy_anal,gxz_anal,gyy_anal,gyz_anal,gzz_anal, &
       gupxx_anal,gupxy_anal,gupxz_anal,gupyy_anal,gupyz_anal,gupzz_anal, &
       Gammax_anal,Gammay_anal,Gammaz_anal, &
       Symmetry);

end subroutine linwave_DumpAnalytic
