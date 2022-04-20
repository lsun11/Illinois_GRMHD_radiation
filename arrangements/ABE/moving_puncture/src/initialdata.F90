#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine movpunc_initial_data(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer, dimension(3)                    :: ext
  real*8 			           :: dT,dX,dY,dZ
  real*8 :: xmin,ymin,zmin,rho_fail_max_step,m_fail_step,reduction_value
  real*8 :: rho_max,xx,yy,zz,rr,xmax2,ymax2,zmax2
  integer :: i,j,l
  integer :: handle,dummy,proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax,glob_imax,glob_jmax,glob_kmax
  integer :: index
  integer :: ierr
  real*8 :: ONE,ZERO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  parameter(ONE = 1.D0, ZERO = 0.D0)

  ext(1) = cctk_lsh(1)
  ext(2) = cctk_lsh(2)
  ext(3) = cctk_lsh(3)  
  dT = CCTK_DELTA_TIME
  dX = X(2,1,1) - X(1,1,1)
  dY = Y(1,2,1) - Y(1,1,1)
  dZ = Z(1,1,2) - Z(1,1,1)

  write(6,*)'puncture_enable:',puncture_enable

  if(puncture_enable >= 1) then 
     
     write(6,*)'calling laguna!',bhpos
     
     call laguna_initial_metric(ext,X,Y,Z,PhysicalRadius, &
          Axx,Axy,Axz,Ayy,Ayz,Azz, &
          lapm1,shiftx,shifty,shiftz, &
          phi,bigP,bhpos,xbh1,zbh1,Mbh,Symmetry,puncture_enable)

     !     ahf_xc=xbh1
     ! ahf_yx=ybh1
     !     ahf_zc=zbh1
     
     write(6,*)'init data: bhx,z',xbh1,zbh1 !,ahf_xc,ahf_zc

     !======================================
     ! Set everything else to flat data
     !======================================
     gxx = ONE
     gxy = ZERO
     gxz = ZERO
     gyy = ONE
     gyz = ZERO
     gzz = ONE
     trK = ZERO
     gupxx = ONE
     gupxy = ZERO
     gupxz = ZERO
     gupyy = ONE
     gupyz = ZERO
     gupzz = ONE

     Gammax = ZERO
     Gammay = ZERO
     Gammaz = ZERO
     
     !======================================
     ! Set everything else to Zero!
     !======================================
     Gammaxxx = ZERO
     Gammaxxy = ZERO
     Gammaxxz = ZERO
     Gammaxyy = ZERO
     Gammaxyz = ZERO
     Gammaxzz = ZERO
     Gammayxx = ZERO
     Gammayxy = ZERO
     Gammayxz = ZERO
     Gammayyy = ZERO
     Gammayyz = ZERO
     Gammayzz = ZERO
     Gammazxx = ZERO
     Gammazxy = ZERO
     Gammazxz = ZERO
     Gammazyy = ZERO
     Gammazyz = ZERO
     Gammazzz = ZERO
     rho = ZERO
  else ! Setting up K-S metric!
     
     call ks_initial_metric(ext,X,Y,Z, &
          gxx,gxy,gxz,gyy,gyz,gzz,trK, &
          Axx,Axy,Axz,Ayy,Ayz,Azz, &
          lapm1,shiftx,shifty,shiftz,ZERO,Symmetry)
     
     !|~~~~~> convert to tilded metric and initialize \phi and \Gamma^a
     call convert_a(ext,X,Y,Z,phi, &
          gxx,  gxy,  gxz,  gyy,  gyz,  gzz, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
          gxxx,gxxy,gxxz,gxyx,gxyy,gxyz, &
          gxzx,gxzy,gxzz,gyyx,gyyy,gyyz, &
          gyzx,gyzy,gyzz,gzzx,gzzy,gzzz, &
          Gammax,Gammay,Gammaz, &
          phix,phiy,phiz,lapsex,lapsey,lapsez, &
          Symmetry)

     call anagam(ext, &
          gxx,gxy,gxz,gyy,gyz,gzz, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
          Gammaxxx, Gammaxxy, Gammaxxz, &
          Gammaxyy, Gammaxyz, Gammaxzz, &
          Gammayxx, Gammayxy, Gammayxz, &
          Gammayyy, Gammayyz, Gammayzz, &
          Gammazxx, Gammazxy, Gammazxz, &
          Gammazyy, Gammazyz, Gammazzz, &
          gxxx,gxyx,gxzx,gyyx,gyzx,gzzx, &
          gxxy,gxyy,gxzy,gyyy,gyzy,gzzy, &
          gxxz,gxyz,gxzz,gyyz,gyzz,gzzz)
     
     !|~~~~~> covert traceless K_{ij} to tilded traceless K_{ij}
     !use trK_p as a temp variable here for psi4:
     trK_p = exp(-4.D0*phi)
     Axx = trK_p*Axx
     Axy = trK_p*Axy
     Axz = trK_p*Axz
     Ayy = trK_p*Ayy
     Ayz = trK_p*Ayz
     Azz = trK_p*Azz
     !|~~~~~> Set everything else to Zero!
     lapset = ZERO
     shiftxt = ZERO
     shiftyt = ZERO
     shiftzt = ZERO
     
     rho = ZERO     

  end if

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
    
  glob_imax = ubound(rho,1)
  glob_jmax = ubound(rho,2)
  glob_kmax = ubound(rho,3)
  
  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"grid::X")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xmax2,1,index)
  call CCTK_VarIndex(index,"grid::Y")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,ymax2,1,index)
  call CCTK_VarIndex(index,"grid::Z")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,zmax2,1,index)
  call CCTK_ReductionHandle(handle,"minimum")
  call CCTK_VarIndex(index,"grid::X")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xmin,1,index)
  call CCTK_VarIndex(index,"grid::Y")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,ymin,1,index)
  call CCTK_VarIndex(index,"grid::Z")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,zmin,1,index)
  
  if(X(cctk_lsh(1),1,1) .eq. xmax2) then
     proc_imax = glob_imax
  else
     proc_imax = -1
  end if
  if(Y(1,cctk_lsh(2),1) .eq. ymax2) then
     proc_jmax = glob_jmax
  else
     proc_jmax = -1
  end if
  if(Z(1,1,cctk_lsh(3)) .eq. zmax2) then
     proc_kmax = glob_kmax
  else
     proc_kmax = -1
  end if
  ! Possible bug if xmin,ymin,zmin not set (properly):
  !  if(abs(X(1,1,1)) .eq. dX*0.5D0) then
  ! Possible bug if xmin,ymin,zmin not set (properly):
  if(X(1,1,1) .eq. xmin) then
     proc_imin = 0
  else
     proc_imin = -100
  end if
  if(Y(1,1,1) .eq. ymin) then
     proc_jmin = 0
  else
     proc_jmin = -100
  end if
  if(Z(1,1,1) .eq. zmin) then
     proc_kmin = 0
  else
     proc_kmin = -100
  end if

end subroutine movpunc_initial_data
