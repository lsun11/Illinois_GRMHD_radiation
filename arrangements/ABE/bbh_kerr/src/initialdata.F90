#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bbh_kerr_initial_data(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,global_ext
  real*8 			           :: dT,dX,dY,dZ
  real*8, parameter  :: ZERO = 0.D0
  integer :: n1,n2,n3,mf, dummy
  real*8, parameter :: f1o3 = 1.d0/3.d0

  ext = cctk_lsh

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ! Initialize lapse, shift, metric and extrinsic curvature
  ! Note: this must be done before calling any one of 
  !       the initial metric subroutine for one BH
  lapm1 = 0.d0
  shiftx = 0.d0
  shifty = 0.d0
  shiftz = 0.d0
  gxx = 1.d0
  gxy = 0.d0
  gxz = 0.d0
  gyy = 1.d0
  gyz = 0.d0
  gzz = 1.d0
  Axx = 0.d0
  Axy = 0.d0
  Axz = 0.d0
  Ayy = 0.d0
  Ayz = 0.d0
  Azz = 0.d0
  psi = 1.d0

  if (coordinate_type .ne. 3) then 
     write(*,*) 'Coordinate_type ',coordinate_type,' no longer supported!'
     stop
  end if

  !!if (coordinate_type==1) then 

  !!   call ks_initial_metric_plus_junk_one_BH(ext, x, y, z,       &
  !!                  gxx, gxy, gxz, gyy, gyz, gzz,             &
  !!                  Axx, Axy, Axz, Ayy, Ayz, Azz, lapm1, shiftx,shifty,shiftz, &
  !!                  M0_BH1, sam1, r_junk1, lambda1, x0_BH1,y0_BH1,z0_BH1, &
  !!      	    x0_BH2,y0_BH2,z0_BH2,r2_shield)

  !!   call ks_initial_metric_plus_junk_one_BH(ext, x, y, z,       &
  !!                  gxx, gxy, gxz, gyy, gyz, gzz,             &
  !!                  Axx, Axy, Axz, Ayy, Ayz, Azz, lapm1, shiftx,shifty,shiftz, &
  !!                  M0_BH2, sam2, r_junk2, lambda2, x0_BH2,y0_BH2,z0_BH2, &
  !!                  x0_BH1,y0_BH1,z0_BH1,r1_shield)

  !!else if (coordinate_type==2) then 

  !!   call ks_initial_metric_bl_junk_one_BH(ext, x, y, z,             &
  !!     gxx, gxy, gxz, gyy, gyz, gzz,              &
  !!     Axx, Axy, Axz, Ayy, Ayz, Azz, lapm1, shiftx,shifty,shiftz, &
  !!     x0_BH1,y0_BH1,z0_BH1,x0_BH2,y0_BH2,z0_BH2,M0_BH1, sam1, r_junk1, lambda1, &
  !!                  r2_shield)

  !!   call ks_initial_metric_bl_junk_one_BH(ext, x, y, z,             &
  !!     gxx, gxy, gxz, gyy, gyz, gzz,              &
  !!     Axx, Axy, Axz, Ayy, Ayz, Azz, lapm1, shiftx,shifty,shiftz, &
  !!     x0_BH2,y0_BH2,z0_BH2,x0_BH1,y0_BH1,z0_BH1,M0_BH2, sam2, r_junk2, lambda2, &
  !!                  r1_shield)

  !!else

  if (ux_BH1==0.d0 .and. ux_BH2==0.d0) then
     call kerr_initial_metric_puncture_one_BH(ext, x, y, z,   &
       x0_BH1,y0_BH1,z0_BH1,x0_BH2,y0_BH2,z0_BH2,psi,gxx, gxy, gxz, gyy, gyz, gzz,     &
       Axx, Axy, Axz, Ayy, Ayz, Azz, lapm1, shiftx,shifty,shiftz, &
       sam1,M0_BH1,r2_shield)

     call kerr_initial_metric_puncture_one_BH(ext, x, y, z,   &
       x0_BH2,y0_BH2,z0_BH2,x0_BH1,y0_BH1,z0_BH1,psi,gxx, gxy, gxz, gyy, gyz, gzz,     &
       Axx, Axy, Axz, Ayy, Ayz, Azz, lapm1, shiftx,shifty,shiftz, &
       sam2,M0_BH2,r1_shield)

     psi = psi**4
     gxx = gxx*psi
     gxy = gxy*psi
     gxz = gxz*psi
     gyy = gyy*psi
     gyz = gyz*psi
     gzz = gzz*psi
     Axx = Axx*psi
     Axy = Axy*psi
     Axz = Axz*psi
     Ayy = Ayy*psi
     Ayz = Ayz*psi
     Azz = Azz*psi
  else
     call compute_boost_gij_Kij(M0_BH1,sam1,x0_BH1,y0_BH1,z0_BH1,ux_BH1,ext,x,y,z, &
                        lapm1,shiftx,shifty,shiftz,gxx,gxy,gxz, &
                        gyy,gyz,gzz,Axx,Axy,Axz,Ayy,Ayz,Azz)
! *** TEST ***
  !!psi = (gxx * gyy * gzz + &
  !!          gxy * gyz * gxz + &
  !!          gxz * gxy * gyz &
  !!          - gxz * gyy * gxz &
  !!          - gxy * gxy * gzz &
  !!          - gxx * gyz * gyz)

  !!gupxx =   ( gyy * gzz - gyz * gyz ) / psi
  !!gupxy = - ( gxy * gzz - gyz * gxz ) / psi
  !!gupxz =   ( gxy * gyz - gyy * gxz ) / psi
  !!gupyy =   ( gxx * gzz - gxz * gxz ) / psi
  !!gupyz = - ( gxx * gyz - gxy * gxz ) / psi 
  !!gupzz =   ( gxx * gyy - gxy * gxy ) / psi
  !!psi = psi**(1.d0/3.d0)
  !!Kxx = (gupxx*gupxx*Axx + 2.d0*gupxx*gupxy*Axy + 2.d0*gupxx*gupxz*Axz + &
  !!      gupxy*gupxy*Ayy + 2.d0*gupxy*gupxz*Ayz + gupxz*gupxz*Azz)*psi
  !!Kxy = (gupxx*gupxy*Axx + (gupxx*gupyy+gupxy*gupxy)*Axy + (gupxx*gupyz+gupxz*gupxy)*Axz + & 
  !!      gupxy*gupyy*Ayy + (gupxy*gupyz+gupxz*gupyy)*Ayz + gupxz*gupyz*Azz)*psi
  !!Kxz = (gupxx*gupxz*Axx + (gupxx*gupyz+gupxy*gupxz)*Axy + (gupxx*gupzz+gupxz*gupxz)*Axz + & 
  !!      gupxy*gupyz*Ayy + (gupxy*gupzz+gupxz*gupyz)*Ayz + gupxz*gupzz*Azz)*psi
  !!Kyy = (gupxy*gupxy*Axx + 2.d0*gupxy*gupyy*Axy + 2.d0*gupxy*gupyz*Axz + &
  !!      gupyy*gupyy*Ayy + 2.d0*gupyy*gupyz*Ayz + gupyz*gupyz*Azz)*psi
  !!Kyz = (gupxy*gupxz*Axx + (gupxy*gupyz+gupyy*gupxz)*Axy + (gupxy*gupzz+gupyz*gupxz)*Axz + &
  !!      gupyy*gupyz*Ayy + (gupyy*gupzz+gupyz*gupyz)*Ayz + gupyz*gupzz*Azz)*psi
  !!Kzz = (gupxz*gupxz*Axx + 2.d0*gupxz*gupyz*Axy + 2.d0*gupxz*gupzz*Axz + &
  !!      gupyz*gupyz*Ayy + 2.d0*gupyz*gupzz*Ayz + gupzz*gupzz*Azz)*psi

  !!Sxx = gxx/psi - 1.d0
  !!Sxy = gxy/psi 
  !!Sxz = gxz/psi
  !!Syy = gyy/psi - 1.d0
  !!Syz = gyz/psi
  !!Szz = gzz/psi - 1.d0
  !!psi = psi**(0.25d0) - 1.d0
  !!!!Sxx = Axx*psi**10
  !!!!Sxy = Axy*psi**10
  !!!!Sxz = Sxz*psi**10
  !!!!Syy = Syy*psi**10
  !!!!Syz = Syz*psi**10
  !!!!Szz = Szz*psi**10
  !!!!psi = psi -1.d0

  !!gxx = 1.d0
  !!gyy = 1.d0
  !!gzz = 1.d0
  !!gxy = 0.d0
  !!gxz = 0.d0
  !!gyz = 0.d0
  !!Axx = 0.d0
  !!Axy = 0.d0
  !!Axz = 0.d0
  !!Ayy = 0.d0
  !!Ayz = 0.d0
  !!Azz = 0.d0
! ************

     call compute_boost_gij_Kij(M0_BH2,sam2,x0_BH2,y0_BH2,z0_BH2,ux_BH2,ext,x,y,z, &
                        lapm1,shiftx,shifty,shiftz,gxx,gxy,gxz, &
                        gyy,gyz,gzz,Axx,Axy,Axz,Ayy,Ayz,Azz)
! *** TEST ***
  !!phi = (gxx * gyy * gzz + &
  !!          gxy * gyz * gxz + &
  !!          gxz * gxy * gyz &
  !!          - gxz * gyy * gxz &
  !!          - gxy * gxy * gzz &
  !!          - gxx * gyz * gyz)

  !!gupxx =   ( gyy * gzz - gyz * gyz ) / phi
  !!gupxy = - ( gxy * gzz - gyz * gxz ) / phi
  !!gupxz =   ( gxy * gyz - gyy * gxz ) / phi
  !!gupyy =   ( gxx * gzz - gxz * gxz ) / phi
  !!gupyz = - ( gxx * gyz - gxy * gxz ) / phi
  !!gupzz =   ( gxx * gyy - gxy * gxy ) / phi
  !!phi = phi**(1.d0/3.d0)
  !!!!! +K^i_j
  !!!!gxxx = gupxx*Axx + gupxy*Axy + gupxz*Axz + gxxx
  !!!!gxxy = gupxx*Axy + gupxy*Ayy + gupxz*Ayz + gxxy
  !!!!gxxz = gupxx*Axz + gupxy*Ayz + gupxz*Azz + gxxz
  !!!!gxyx = gupxy*Axx + gupyy*Axy + gupyz*Axz + gxyx
  !!!!gxyy = gupxy*Axy + gupyy*Ayy + gupyz*Ayz + gxyy
  !!!!gxyz = gupxy*Axz + gupyy*Ayz + gupyz*Azz + gxyz
  !!!!gxzx = gupxz*Axx + gupyz*Axy + gupzz*Axz + gxzx
  !!!!gxzy = gupxz*Axy + gupyz*Ayy + gupzz*Ayz + gxzy
  !!!!gxzz = gupxz*Axz + gupyz*Ayz + gupzz*Azz + gxzz
  !!Kxx = (gupxx*gupxx*Axx + 2.d0*gupxx*gupxy*Axy + 2.d0*gupxx*gupxz*Axz + &
  !!      gupxy*gupxy*Ayy + 2.d0*gupxy*gupxz*Ayz + gupxz*gupxz*Azz)*phi + Kxx
  !!Kxy = (gupxx*gupxy*Axx + (gupxx*gupyy+gupxy*gupxy)*Axy + (gupxx*gupyz+gupxz*gupxy)*Axz + &
  !!      gupxy*gupyy*Ayy + (gupxy*gupyz+gupxz*gupyy)*Ayz + gupxz*gupyz*Azz)*phi + Kxy
  !!Kxz = (gupxx*gupxz*Axx + (gupxx*gupyz+gupxy*gupxz)*Axy + (gupxx*gupzz+gupxz*gupxz)*Axz + &
  !!      gupxy*gupyz*Ayy + (gupxy*gupzz+gupxz*gupyz)*Ayz + gupxz*gupzz*Azz)*phi + Kxz
  !!Kyy = (gupxy*gupxy*Axx + 2.d0*gupxy*gupyy*Axy + 2.d0*gupxy*gupyz*Axz + &
  !!      gupyy*gupyy*Ayy + 2.d0*gupyy*gupyz*Ayz + gupyz*gupyz*Azz)*phi + Kyy
  !!Kyz = (gupxy*gupxz*Axx + (gupxy*gupyz+gupyy*gupxz)*Axy + (gupxy*gupzz+gupyz*gupxz)*Axz + &
  !!      gupyy*gupyz*Ayy + (gupyy*gupzz+gupyz*gupyz)*Ayz + gupyz*gupzz*Azz)*phi + Kyz
  !!Kzz = (gupxz*gupxz*Axx + 2.d0*gupxz*gupyz*Axy + 2.d0*gupxz*gupzz*Axz + &
  !!      gupyz*gupyz*Ayy + 2.d0*gupyz*gupzz*Ayz + gupzz*gupzz*Azz)*phi + Kzz

  !!gxx = gxx/phi + Sxx
  !!gxy = gxy/phi + Sxy
  !!gxz = gxz/phi + Sxz
  !!gyy = gyy/phi + Syy
  !!gyz = gyz/phi + Syz
  !!gzz = gzz/phi + Szz

  !!phi = phi**(0.25d0)

  !!psi = (psi + phi )**4
  !!gxx = gxx*psi
  !!gxy = gxy*psi
  !!gxz = gxz*psi
  !!gyy = gyy*psi
  !!gyz = gyz*psi
  !!gzz = gzz*psi

  !!psi = (gxx * gyy * gzz + &
  !!          gxy * gyz * gxz + &
  !!          gxz * gxy * gyz &
  !!          - gxz * gyy * gxz &
  !!          - gxy * gxy * gzz &
  !!          - gxx * gyz * gyz)**(-1.d0/3.d0)

  !!Axx = (gxx*gxx*Kxx + 2.d0*gxx*gxy*Kxy + 2.d0*gxx*gxz*Kxz + &
  !!      gxy*gxy*Kyy + 2.d0*gxy*gxz*Kyz + gxz*gxz*Kzz)*psi
  !!Axy = (gxx*gxy*Kxx + (gxx*gyy+gxy*gxy)*Kxy + (gxx*gyz+gxz*gxy)*Kxz + &
  !!      gxy*gyy*Kyy + (gxy*gyz+gxz*gyy)*Kyz + gxz*gyz*Kzz)*psi
  !!Axz = (gxx*gxz*Kxx + (gxx*gyz+gxy*gxz)*Kxy + (gxx*gzz+gxz*gxz)*Kxz + &
  !!      gxy*gyz*Kyy + (gxy*gzz+gxz*gyz)*Kyz + gxz*gzz*Kzz)*psi
  !!Ayy = (gxy*gxy*Kxx + 2.d0*gxy*gyy*Kxy + 2.d0*gxy*gyz*Kxz + &
  !!      gyy*gyy*Kyy + 2.d0*gyy*gyz*Kyz + gyz*gyz*Kzz)*psi
  !!Ayz = (gxy*gxz*Kxx + (gxy*gyz+gyy*gxz)*Kxy + (gxy*gzz+gyz*gxz)*Kxz + &
  !!      gyy*gyz*Kyy + (gyy*gzz+gyz*gyz)*Kyz + gyz*gzz*Kzz)*psi
  !!Azz = (gxz*gxz*Kxx + 2.d0*gxz*gyz*Kxy + 2.d0*gxz*gzz*Kxz + &
  !!      gyz*gyz*Kyy + 2.d0*gyz*gzz*Kyz + gzz*gzz*Kzz)*psi

  !!Kxx = 0.d0
  !!Kxy = 0.d0
  !!Kxz = 0.d0
  !!Kyy = 0.d0
  !!Kyz = 0.d0
  !!Kzz = 0.d0

  !!!!Axx = gxx*gxxx + gxy*gxyx + gxz*gxzx
  !!!!Axy = gxx*gxxy + gxy*gxyy + gxz*gxzy
  !!!!Axz = gxx*gxxz + gxy*gxyz + gxz*gxzz
  !!!!Ayy = gxy*gxxy + gyy*gxyy + gyz*gxzy
  !!!!Ayz = gxy*gxxz + gyy*gxyz + gyz*gxzz
  !!!!Azz = gxz*gxxz + gyz*gxyz + gzz*gxzz
! ************
  end if


  !!end if

  ! Compute trK

  ! Compute psi^4 and store it in psi
  psi = (gxx * gyy * gzz + &
            gxy * gyz * gxz + &
            gxz * gxy * gyz &
            - gxz * gyy * gxz &
            - gxy * gxy * gzz &
            - gxx * gyz * gyz)**(1.d0/3.d0)
  ! Convert g_ij to tilde{g}_ij
  gxx = gxx/psi
  gxy = gxy/psi
  gxz = gxz/psi
  gyy = gyy/psi
  gyz = gyz/psi
  gzz = gzz/psi 
  ! Compute tilde{g}^ij
  gupxx =   ( gyy * gzz - gyz * gyz )
  gupxy = - ( gxy * gzz - gyz * gxz )
  gupxz =   ( gxy * gyz - gyy * gxz )
  gupyy =   ( gxx * gzz - gxz * gxz )
  gupyz = - ( gxx * gyz - gxy * gxz )
  gupzz =   ( gxx * gyy - gxy * gxy )
  ! Compute trK
  trK = (gupxx*Axx + 2.d0*gupxy*Axy + 2.d0*gupxz*Axz + & 
	  	gupyy*Ayy + 2.d0*gupyz*Ayz + gupzz*Azz)/psi
  ! Convert K_ij to tilde{A}_ij
  Axx = Axx/psi - f1o3*trK*gxx
  Axy = Axy/psi - f1o3*trK*gxy
  Axz = Axz/psi - f1o3*trK*gxz
  Ayy = Ayy/psi - f1o3*trK*gyy
  Ayz = Ayz/psi - f1o3*trK*gyz
  Azz = Azz/psi - f1o3*trK*gzz
  ! Compute psi and phi
  psi = psi**(0.25d0)
  phi = log(psi)

 ! Use rho as temporary storage here:
  call setgamma_v2(ext, cctk_nghostzones, dX, dY, dZ, &
     phi, gxx, gxy, gxz, gyy, gyz, gzz, &
     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Gammax, Gammay, Gammaz, rho, &
     Symmetry)

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  
  n1 = 0
  n2 = 0
  n3 = 0
  mf = 1
  !use rho as temporary storage here:
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammax,rho,n1,n2,n3,mf,Symmetry)
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammay,rho,n1,n2,n3,mf,Symmetry)
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammaz,rho,n1,n2,n3,mf,Symmetry)

  ! Setup puncture initial lapse and shift
  lapm1 = psi**(-2) - 1.d0
  !!shiftx = 0.d0
  !!shifty = 0.d0
  !!shiftz = 0.d0

  shiftxt = 0.D0
  shiftyt = 0.D0
  shiftzt = 0.D0
  lapset = 0.D0


  ! Set matter variables to zero!
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

  !Following lines are (empirically speaking) ESSENTIAL for excision runs:
  call Derivs(ext,X,Y,Z,dX,dY,dZ,phi,phix,phiy,phiz,Symmetry)
  call Derivs(ext,X,Y,Z,dX,dY,dZ,lapm1,lapsex,lapsey,lapsez,Symmetry)

  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,phi,phix,phiy,phiz)
  call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,lapm1,lapsex,lapsey,lapsez)

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')

end subroutine bbh_kerr_initial_data
