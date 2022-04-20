#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bbhlorene_initialdata(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,global_ext
  !  real*8,dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))   :: alpha,beta,u,nablau,ux,uy,uz,uxx,uxy,uxz,uyy,uyz,uzz,integrand
  real*8 :: rbh1,rbh2,xphys,yphys,zphys
  real*8 :: nx1,ny1,nz1,nx2,ny2,nz2
  real*8 :: gupxxL,gupxyL,gupxzL,gupyyL,gupyzL,gupzzL
  real*8 			           :: dT,dX,dY,dZ
  real*8 			           :: n1,n2,n3,mf

  real*8                                   :: PI,dmass,detmin,detmax
  real*8                                   :: px1,py1,pz1,px2,py2,pz2
  real*8                                   :: sx1,sy1,sz1,sx2,sy2,sz2
  real*8                                   :: esx1,esy1,esz1,esx2,esy2,esz2
  real*8                                   :: m1,m2,y1,y2,bigM
  real*8, dimension(1,3)                   :: pointcoords

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

  PI = 3.141592653589793238462D0

  ext = cctk_lsh

  write(*,*) "inside initialdata!"

  ! OVERTURE:

  write(*,*) "hello!",each_black_hole_mass,each_black_hole_momentum,each_black_hole_spin,total_binary_separation/2.D0

  bigM = 1.D0

  m1 = each_black_hole_mass*bigM
  m2 = each_black_hole_mass*bigM

  px1 = each_black_hole_momentum*bigM
  px2 = -1.0*each_black_hole_momentum*bigM
  py1 = 0.D0
  py2 = 0.D0
  pz1 = 0.D0
  pz2 = 0.D0

  sx1=0.
  sx2=0.
  sy1=0.
  sy2=0.
  sz1=each_black_hole_spin*bigM**2
  sz2=each_black_hole_spin*bigM**2

  y1 = -1.D0*total_binary_separation/2.D0*bigM
  y2 = total_binary_separation/2.D0*bigM

  !  r = sqrt(x*x + y*y + z*z)
  write(*,*) "Check it!:",m1,px1,sz1,y1,y2

  ! THE MEATY PART:

  gupxx = 1.D0
  gupxy = 0.D0
  gupxz = 0.D0
  gupyy = 1.D0
  gupyz = 0.D0
  gupzz = 1.D0

  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)

           if(abs(r(i,j,k)).gt.dX*1.D-10) then
              xphys = x(i,j,k)*PhysicalRadius(i,j,k)/r(i,j,k)
              yphys = y(i,j,k)*PhysicalRadius(i,j,k)/r(i,j,k)
              zphys = z(i,j,k)*PhysicalRadius(i,j,k)/r(i,j,k)
           else
              xphys = 0.D0
              yphys = 0.D0
              zphys = 0.D0
              write(*,*) "hitasodfapsdofijaspdfoij:",r(i,j,k),x(i,j,k),y(i,j,k),z(i,j,k),abs(r(i,j,k)),dX*1.D-10
              !stop
           end if

           rbh1 = sqrt(xphys**2 + (yphys-y1)**2 + zphys**2)
           rbh2 = sqrt(xphys**2 + (yphys-y2)**2 + zphys**2)

           if(rbh1==0.D0) rbh1 = 1.D-10*dX
           if(rbh2==0.D0) rbh2 = 1.D-10*dX

           ny1 = (yphys-y1)/rbh1
           ny2 = (yphys-y2)/rbh2

           nx1 = xphys/rbh1
           nx2 = xphys/rbh2

           nz1 = zphys/rbh1
           nz2 = zphys/rbh2

           esx1=sy1*nz1-sz1*ny1
           esx2=sy2*nz2-sz2*ny2
           esy1=sz1*nx1-sx1*nz1
           esy2=sz2*nx2-sx2*nz2
           esz1=sx1*ny1-sy1*nx1
           esz2=sx2*ny2-sy2*nx2

           gupxxL = gupxx(i,j,k)
           gupxyL = gupxy(i,j,k)
           gupxzL = gupxz(i,j,k)
           gupyyL = gupyy(i,j,k)
           gupyzL = gupyz(i,j,k)
           gupzzL = gupzz(i,j,k)

           ! In the below lines, K_ij actually == Abar_ij
           Kxx(i,j,k) =       3.D0/(2.D0*rbh1*rbh1) * ((2.D0 * px1*nx1) - (gupxxL - nx1*nx1)*(px1*nx1 + py1*ny1 + pz1*nz1))
           Kxx(i,j,k) = Kxx(i,j,k) + 3.D0/(2.D0*rbh2*rbh2) * ((2.D0 * px2*nx2) - (gupxxL - nx2*nx2)*(px2*nx2 + py2*ny2 + pz2*nz2))

           Kxy(i,j,k) =       3.D0/(2.D0*rbh1*rbh1) * ((px1*ny1 + py1*nx1) - (gupxyL- nx1*ny1)*(px1*nx1 + py1*ny1 + pz1*nz1))
           Kxy(i,j,k) = Kxy(i,j,k) + 3.D0/(2.D0*rbh2*rbh2) * ((px2*ny2 + py2*nx2) - (gupxyL- nx2*ny2)*(px2*nx2 + py2*ny2 + pz2*nz2))

           Kxz(i,j,k) =       3.D0/(2.D0*rbh1*rbh1) * ((px1*nz1 + pz1*nx1) - (gupxzL- nx1*nz1)*(px1*nx1 + py1*ny1 + pz1*nz1))
           Kxz(i,j,k) = Kxz(i,j,k) + 3.D0/(2.D0*rbh2*rbh2) * ((px2*nz2 + pz2*nx2) - (gupxzL- nx2*nz2)*(px2*nx2 + py2*ny2 + pz2*nz2))

           Kyy(i,j,k) =       3.D0/(2.D0*rbh1*rbh1) * ((2.D0 * py1*ny1) - (gupyyL - ny1*ny1)*(px1*nx1 + py1*ny1 + pz1*nz1))
           Kyy(i,j,k) = Kyy(i,j,k) + 3.D0/(2.D0*rbh2*rbh2) * ((2.D0 * py2*ny2) - (gupyyL - ny2*ny2)*(px2*nx2 + py2*ny2 + pz2*nz2))

           Kyz(i,j,k) =       3.D0/(2.D0*rbh1*rbh1) * ((py1*nz1 + pz1*ny1) - (gupyzL- ny1*nz1)*(px1*nx1 + py1*ny1 + pz1*nz1))
           Kyz(i,j,k) = Kyz(i,j,k) + 3.D0/(2.D0*rbh2*rbh2) * ((py2*nz2 + pz2*ny2) - (gupyzL- ny2*nz2)*(px2*nx2 + py2*ny2 + pz2*nz2))

           Kzz(i,j,k) =       3.D0/(2.D0*rbh1*rbh1) * ((2.D0 * pz1*nz1) - (gupzzL - nz1*nz1)*(px1*nx1 + py1*ny1 + pz1*nz1))
           Kzz(i,j,k) = Kzz(i,j,k) + 3.D0/(2.D0*rbh2*rbh2) * ((2.D0 * pz2*nz2) - (gupzzL - nz2*nz2)*(px2*nx2 + py2*ny2 + pz2*nz2))

           if(each_black_hole_spin.ne.0) then
              Kxx(i,j,k)=Kxx(i,j,k) + 6.0d0/rbh1**3*esx1*nx1 + 6.0/rbh2**3*esx2*nx2
              Kxy(i,j,k)=Kxy(i,j,k) + 3.0d0/rbh1**3*(esx1*ny1+esy1*nx1) + 3.0d0/rbh2**3*(esx2*ny2+esy2*nx2)
              Kxz(i,j,k)=Kxz(i,j,k) + 3.0d0/rbh1**3*(esx1*nz1+esz1*nx1) + 3.0d0/rbh2**3*(esx2*nz2+esz2*nx2)
              Kyy(i,j,k)=Kyy(i,j,k) + 6.0d0/rbh1**3*esy1*ny1 + 6.0/rbh2**3*esy2*ny2
              Kyz(i,j,k)=Kyz(i,j,k) + 3.0d0/rbh1**3*(esy1*nz1+esz1*ny1) + 3.0d0/rbh2**3*(esy2*nz2+esz2*ny2)
              Kzz(i,j,k)=Kzz(i,j,k) + 6.0d0/rbh1**3*esz1*nz1 + 6.0/rbh2**3*esz2*nz2
           endif
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

  !  write(*,*) "start!"
  !  call invert_g(ext,gxx,gxy,gxz,gyy,gyz,gzz,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,detmin,detmax)
  !  write(*,*) "finished!"

  ! Convert Abar_ij to K_ij
  Kxx = Kxx / (Psi*Psi)
  Kxy = Kxy / (Psi*Psi)
  Kxz = Kxz / (Psi*Psi)
  Kyy = Kyy / (Psi*Psi)
  Kyz = Kyz / (Psi*Psi)
  Kzz = Kzz / (Psi*Psi)


  shiftx = 0.D0
  shifty = 0.D0
  shiftz = 0.D0

  lapm1 = psi**(-2) - 1.D0
  where(lapm1+1.D0.lt.1.D-10)
     lapm1 = (1.D-8)-1.D0
  end where

  ! Set excision_zone_gf to avoid valgrind memory errors
  excision_zone_gf = 0

  !Compute Atilde_ij's
  Axx = Kxx / (Psi*Psi*Psi*Psi)
  Axy = Kxy / (Psi*Psi*Psi*Psi)
  Axz = Kxz / (Psi*Psi*Psi*Psi)
  Ayy = Kyy / (Psi*Psi*Psi*Psi)
  Ayz = Kyz / (Psi*Psi*Psi*Psi)
  Azz = Kzz / (Psi*Psi*Psi*Psi)

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

  write(*,*) "setup id:",Axx(31,31,3),phi(31,31,3)

  if(fisheye_enable .eq. 1) then
     call trans_phys_fish_tensor_flat(cctkGH,cctk_lsh,Symmetry, &
          x, y, z,  &
          PhysicalRadius,RadiusDerivative, &
          gxx,gxy,gxz,gyy,gyz,gzz)
     call trans_phys_fish_tensor_inv(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)
     call trans_phys_fish_tensor(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative, &
          Axx,Axy,Axz,Ayy,Ayz,Azz)
     call trans_phys_fish_phi(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative,phi)
     call trans_phys_fish_gamt_flat(cctkGH,cctk_lsh,Symmetry, & 
          x, y, z,  &
          PhysicalRadius,RadiusDerivative,RadiusDerivative2, &
          Gammax,Gammay,Gammaz)
  end if

  write(*,*) "setup id:",Axx(31,31,3),phi(31,31,3)

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

end subroutine bbhlorene_initialdata
