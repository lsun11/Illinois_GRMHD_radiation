!----------------------------------------------------------
! Driver routine for computing RHS and upwind differencing
!----------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine BSSN_timestepping(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8 					:: dT, dX,dY,dZ,detmin,detmax
  integer, dimension(3)			:: ext,fake_ext
  integer					:: dummy,i,j,k,index
  integer, parameter				:: AXISYM = 4
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  !

  dT = CCTK_DELTA_TIME
  dX=CCTK_DELTA_SPACE(1)
  dY=CCTK_DELTA_SPACE(2)
  dZ=CCTK_DELTA_SPACE(3)

  ext = cctk_lsh
  

     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
           if(isnan(tau_rhs(i,j,k))) then
                 write(*,*) "bssn driver timestepping:  Checkpoint << 1 >>, tau_rhs NAN FOUND!"
           end if
           end do
        end do
     end do



  !zero out source terms if zero_out_matter_source_terms==1
  if (zero_out_matter_source_terms .eq. 1) then
     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              rho(i,j,k)=0.d0
              Sx(i,j,k)=0.d0
              Sy(i,j,k)=0.d0
              Sz(i,j,k)=0.d0
              Sxx(i,j,k)=0.d0 
              Sxy(i,j,k)=0.d0
              Sxz(i,j,k)=0.d0
              Syy(i,j,k)=0.d0
              Syz(i,j,k)=0.d0
              Szz(i,j,k)=0.d0
              S(i,j,k)=0.d0
           end do
        end do
     end do
     !$omp end parallel do
  endif


  call BSSN_compute_gupij(cctkGH,cctk_lsh, &
       gxx,gxy,gxz,gyy,gyz,gzz, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)
  call BSSN_rhs(cctkGH, dT, dx, dy, dz, &
       cctk_nghostzones,cctk_lsh, &
       PhysicalRadius, RadiusDerivative,r, &
       gxx,gxy,gxz,gyy,gyz,gzz, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
       Axx,Axy,Axz,Ayy,Ayz,Azz, &
       Rxx,Rxy,Rxz,Ryy,Ryz,Rzz, &
       trRtilde,phi,chi,PsiRes,trK,lapm1, &
       lapsex,lapsey,lapsez, &
       shiftx,shifty,shiftz, &
       Gammax,Gammay,Gammaz, &
       S,rho, &
       Sx,Sy,Sz, &
       Sxx,Sxy,Sxz,Syy,Syz,Szz, &
       gxx_rhs,gxy_rhs,gxz_rhs,gyy_rhs,gyz_rhs,gzz_rhs, &
       Axx_rhs,Axy_rhs,Axz_rhs,Ayy_rhs,Ayz_rhs,Azz_rhs, &
       phi_rhs,chi_rhs,trK_rhs, &
       Gammax_rhs,Gammay_rhs,Gammaz_rhs, &
       phix, phiy, phiz, &
       gxxx, gxxy, gxxz, &
       gxyx, gxyy, gxyz, &
       gxzx, gxzy, gxzz, &
       gyyx, gyyy, gyyz, &
       gyzx, gyzy, gyzz, &
       gzzx, gzzy, gzzz, &
       Gammaxxx, Gammaxxy, Gammaxxz, Gammaxyy, Gammaxyz, Gammaxzz, &
       Gammayxx, Gammayxy, Gammayxz, Gammayyy, Gammayyz, Gammayzz, &
       Gammazxx, Gammazxy, Gammazxz, Gammazyy, Gammazyz, Gammazzz, &
       x,y,z, &
       bh_posn_x,bh_posn_y,bh_posn_z)

  if(enable_lower_order_at_boundaries==1) then
     write(*,*) "Although enable_lower_order_at_boundaries==1 is supported, the code has been disabled, due to long compile times."
     write(*,*) "You can re-enable this code by uncommenting all BSSN_ricci_and_constraints_{2,4} function calls and re-adding compute_rhsnew_2ndorderbdry.C and compute_rhsnew_4thorderbdry.C to make.code.defn inside the bssn/src directory."
     stop


!!$     call BSSN_rhs_4(cctkGH, dT, dx, dy, dz, &
!!$          cctk_nghostzones,cctk_lsh, &
!!$          PhysicalRadius, RadiusDerivative,r, &
!!$          gxx,gxy,gxz,gyy,gyz,gzz, &
!!$          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
!!$          Axx,Axy,Axz,Ayy,Ayz,Azz, &
!!$          Rxx,Rxy,Rxz,Ryy,Ryz,Rzz, &
!!$          trRtilde,phi,chi,PsiRes,trK,lapm1, &
!!$          lapsex,lapsey,lapsez, &
!!$          shiftx,shifty,shiftz, &
!!$          Gammax,Gammay,Gammaz, &
!!$          S,rho, &
!!$          Sx,Sy,Sz, &
!!$          Sxx,Sxy,Sxz,Syy,Syz,Szz, &
!!$          gxx_rhs,gxy_rhs,gxz_rhs,gyy_rhs,gyz_rhs,gzz_rhs, &
!!$          Axx_rhs,Axy_rhs,Axz_rhs,Ayy_rhs,Ayz_rhs,Azz_rhs, &
!!$          phi_rhs,chi_rhs,trK_rhs, &
!!$          Gammax_rhs,Gammay_rhs,Gammaz_rhs, &
!!$          phix, phiy, phiz, &
!!$          gxxx, gxxy, gxxz, &
!!$          gxyx, gxyy, gxyz, &
!!$          gxzx, gxzy, gxzz, &
!!$          gyyx, gyyy, gyyz, &
!!$          gyzx, gyzy, gyzz, &
!!$          gzzx, gzzy, gzzz, &
!!$          Gammaxxx, Gammaxxy, Gammaxxz, Gammaxyy, Gammaxyz, Gammaxzz, &
!!$          Gammayxx, Gammayxy, Gammayxz, Gammayyy, Gammayyz, Gammayzz, &
!!$          Gammazxx, Gammazxy, Gammazxz, Gammazyy, Gammazyz, Gammazzz, &
!!$          x,y,z, &
!!$          bh_posn_x,bh_posn_y,bh_posn_z)
!!$
!!$     call BSSN_rhs_2(cctkGH, dT, dx, dy, dz, &
!!$          cctk_nghostzones,cctk_lsh, &
!!$          PhysicalRadius, RadiusDerivative,r, &
!!$          gxx,gxy,gxz,gyy,gyz,gzz, &
!!$          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
!!$          Axx,Axy,Axz,Ayy,Ayz,Azz, &
!!$          Rxx,Rxy,Rxz,Ryy,Ryz,Rzz, &
!!$          trRtilde,phi,chi,PsiRes,trK,lapm1, &
!!$          lapsex,lapsey,lapsez, &
!!$          shiftx,shifty,shiftz, &
!!$          Gammax,Gammay,Gammaz, &
!!$          S,rho, &
!!$          Sx,Sy,Sz, &
!!$          Sxx,Sxy,Sxz,Syy,Syz,Szz, &
!!$          gxx_rhs,gxy_rhs,gxz_rhs,gyy_rhs,gyz_rhs,gzz_rhs, &
!!$          Axx_rhs,Axy_rhs,Axz_rhs,Ayy_rhs,Ayz_rhs,Azz_rhs, &
!!$          phi_rhs,chi_rhs,trK_rhs, &
!!$          Gammax_rhs,Gammay_rhs,Gammaz_rhs, &
!!$          phix, phiy, phiz, &
!!$          gxxx, gxxy, gxxz, &
!!$          gxyx, gxyy, gxyz, &
!!$          gxzx, gxzy, gxzz, &
!!$          gyyx, gyyy, gyyz, &
!!$          gyzx, gyzy, gyzz, &
!!$          gzzx, gzzy, gzzz, &
!!$          Gammaxxx, Gammaxxy, Gammaxxz, Gammaxyy, Gammaxyz, Gammaxzz, &
!!$          Gammayxx, Gammayxy, Gammayxz, Gammayyy, Gammayyz, Gammayzz, &
!!$          Gammazxx, Gammazxy, Gammazxz, Gammazyy, Gammazyz, Gammazzz, &
!!$          x,y,z, &
!!$          bh_posn_x,bh_posn_y,bh_posn_z)

  end if

  !This function is designed to stabilize puncture BH evolutions when the puncture sits on a gridpoint
  if(enable_second_order_inside_ah==1) then
     call BSSN_rhs_AHinterior(cctkGH, dT, dx, dy, dz, &
          cctk_nghostzones,cctk_lsh, &
          PhysicalRadius, RadiusDerivative,r, &
          gxx,gxy,gxz,gyy,gyz,gzz, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
          Axx,Axy,Axz,Ayy,Ayz,Azz, &
          Rxx,Rxy,Rxz,Ryy,Ryz,Rzz, &
          trRtilde,phi,chi,PsiRes,trK,lapm1, &
          lapsex,lapsey,lapsez, &
          shiftx,shifty,shiftz, &
          Gammax,Gammay,Gammaz, &
          S,rho, &
          Sx,Sy,Sz, &
          Sxx,Sxy,Sxz,Syy,Syz,Szz, &
          gxx_rhs,gxy_rhs,gxz_rhs,gyy_rhs,gyz_rhs,gzz_rhs, &
          Axx_rhs,Axy_rhs,Axz_rhs,Ayy_rhs,Ayz_rhs,Azz_rhs, &
          phi_rhs,chi_rhs,trK_rhs, &
          Gammax_rhs,Gammay_rhs,Gammaz_rhs, &
          phix, phiy, phiz, &
          gxxx, gxxy, gxxz, &
          gxyx, gxyy, gxyz, &
          gxzx, gxzy, gxzz, &
          gyyx, gyyy, gyyz, &
          gyzx, gyzy, gyzz, &
          gzzx, gzzy, gzzz, &
          Gammaxxx, Gammaxxy, Gammaxxz, Gammaxyy, Gammaxyz, Gammaxzz, &
          Gammayxx, Gammayxy, Gammayxz, Gammayyy, Gammayyz, Gammayzz, &
          Gammazxx, Gammazxy, Gammazxz, Gammazyy, Gammazyz, Gammazzz)
  end if

  if (bssn_enable_shift_upwind == 1) then
     if(chi_evolution_enable==0) then
        call shift_upwind2(cctkGH,cctk_lsh,cctk_nghostzones, &
             dT,dX,dY,dZ, &
             gxx,gxy,gxz,gyy,gyz,gzz, &
             gxx_rhs,gxy_rhs,gxz_rhs, &
             gyy_rhs,gyz_rhs,gzz_rhs, &
             Axx,Axy,Axz,Ayy,Ayz,Azz, &
             Axx_rhs,Axy_rhs,Axz_rhs, &
             Ayy_rhs,Ayz_rhs,Azz_rhs, &
             phi, phi_rhs, &
             trK,trK_rhs, &
             Gammax,Gammay,Gammaz, &
             Gammax_rhs,Gammay_rhs,Gammaz_rhs, &
             shiftx,shifty,shiftz,Symmetry)
     else 
        call shift_upwind2(cctkGH,cctk_lsh,cctk_nghostzones, &
             dT,dX,dY,dZ, &
             gxx,gxy,gxz,gyy,gyz,gzz, &
             gxx_rhs,gxy_rhs,gxz_rhs, &
             gyy_rhs,gyz_rhs,gzz_rhs, &
             Axx,Axy,Axz,Ayy,Ayz,Azz, &
             Axx_rhs,Axy_rhs,Axz_rhs, &
             Ayy_rhs,Ayz_rhs,Azz_rhs, &
             chi, chi_rhs, &
             trK,trK_rhs, &
             Gammax,Gammay,Gammaz, &
             Gammax_rhs,Gammay_rhs,Gammaz_rhs, &
             shiftx,shifty,shiftz,Symmetry)
     end if
  end if

  if(excision_enable == 1) then
     call metric_excision_bc(ext,X,Y,Z, &
          chi_rhs,   trK_rhs, &
          gxx_rhs,   gxy_rhs,   gxz_rhs, &
          gyy_rhs,   gyz_rhs,   gzz_rhs, &
          Axx_rhs,   Axy_rhs,   Axz_rhs, &
          Ayy_rhs,   Ayz_rhs,   Azz_rhs, &
          Gammax_rhs,Gammay_rhs,Gammaz_rhs, &
          Symmetry,excision_zone_gf)
     call kreiss_oliger_timestepping(CCTK_PASS_FTOF)
  end if

  !Need the following line so that gupij_f's are correctly computed in mhd_evolve!
  call fill_bssn_symmetry_gz_gupij(ext,X,Y,Z,Symmetry,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)



     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
           if(isnan(tau_rhs(i,j,k))) then
                 write(*,*) " bssn driver timestepping:  Checkpoint << 2 >>, tau_rhs NAN FOUND!"
           end if

          end do
	end do
     end do


end subroutine BSSN_timestepping
