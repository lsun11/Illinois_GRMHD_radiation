#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bhns_initialdata_local(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext

  real*8 			           :: dT,dX,dY,dZ
  real*8 			           :: n1,n2,n3,mf

  integer :: vindex,handle,ierr
  real*8  :: rho_max,tau_max
  real*8  :: ONE,ZERO
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(ONE = 1.D0, ZERO = 0.D0)

  INTEGER :: i,j,k
  character :: filename*30,c2*2,c1
  CCTK_REAL :: reduction_value

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ext = cctk_lsh

  ! set psi and phi correctly, after the symmetry operator
  psi=phi
  phi=log(phi)

  if(use_new_bhns_initial_data==0) then
     Axx=Axx/(lapm1+ONE)
     Axy=Axy/(lapm1+ONE)
     Axz=Axz/(lapm1+ONE)
     Ayy=Ayy/(lapm1+ONE)
     Ayz=Ayz/(lapm1+ONE)
     Azz=-(Axx+Ayy)
  else 
     Axx=Axx*psi**(-4)
     Axy=Axy*psi**(-4)
     Axz=Axz*psi**(-4)
     Ayy=Ayy*psi**(-4)
     Ayz=Ayz*psi**(-4)
     !     Azz=-(Axx+Ayy)
     Azz=Azz*psi**(-4)
  end if

  !The following chunk of code MUST come AFTER Aij is set (when use_new_bhns_initial_data==0, at least)!
  if(reset_shift_lapse==1) then
     !     lapm1=psi**(-2)-1.0d0
     !     lapm1=psi**(-1)-1.0d0
     !     lapm1 = 0.D0

     shiftx=0.
     shifty=0.
     shiftz=0.
  endif

  write(6,*)'phi,axx-1:a',phi(1,1,1)
  write(6,*)'phi,axx-1:b',Axx(1,1,1)
  write(6,*)'phi,axx-1:c',Ayy(1,1,1)
  write(6,*)'phi,axx-1:d',rho_b(1,1,1)

  !call CartSymGN(ierr,cctkGH,'mhd_evolve::mhd_primitives')
  !$omp parallel do
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)

           P(i,j,k)=ZERO
           rho_star(i,j,k)=ZERO
           h(i,j,k)=ZERO
           w(i,j,k)=ZERO
           tau(i,j,k)=ZERO
           st_x(i,j,k)=ZERO
           st_y(i,j,k)=ZERO
           st_z(i,j,k)=ZERO

           sbt(i,j,k)=ZERO
           sbx(i,j,k)=ZERO
           sby(i,j,k)=ZERO
           sbz(i,j,k)=ZERO
           Bx(i,j,k)=ZERO
           By(i,j,k)=ZERO
           Bz(i,j,k)=ZERO
           Ex(i,j,k)=ZERO
           Ey(i,j,k)=ZERO
           Ez(i,j,k)=ZERO
           mhd_st_x(i,j,k)=ZERO
           mhd_st_y(i,j,k)=ZERO
           mhd_st_z(i,j,k)=ZERO

           gupxx(i,j,k)=ONE
           gupyy(i,j,k)=ONE
           gupzz(i,j,k)=ONE
           gupxy(i,j,k)=ZERO
           gupxz(i,j,k)=ZERO
           gupyz(i,j,k)=ZERO

           gxx(i,j,k)=ONE
           gyy(i,j,k)=ONE
           gzz(i,j,k)=ONE
           gxy(i,j,k)=ZERO
           gxz(i,j,k)=ZERO
           gyz(i,j,k)=ZERO

           ! Set excision_zone_gf to avoid valgrind memory errors
           excision_zone_gf = 0

           shiftxt(i,j,k)= 0.D0
           shiftyt(i,j,k)= 0.D0
           shiftzt(i,j,k)= 0.D0
           lapset(i,j,k)= 0.D0

           trK(i,j,k)= 0.D0

           Gammax(i,j,k)= 0.D0
           Gammay(i,j,k)= 0.D0
           Gammaz(i,j,k)= 0.D0

           !======================================
           ! Set everything else to Zero! 
           !======================================
           !Set derivatives to zero
           gxxx(i,j,k) = ZERO
           gxxy(i,j,k) = ZERO 
           gxxz(i,j,k) = ZERO
           gxyx(i,j,k) = ZERO 
           gxyy(i,j,k) = ZERO 
           gxyz(i,j,k) = ZERO
           gxzx(i,j,k) = ZERO 
           gxzy(i,j,k) = ZERO 
           gxzz(i,j,k) = ZERO
           gyyx(i,j,k) = ZERO 
           gyyy(i,j,k) = ZERO 
           gyyz(i,j,k) = ZERO
           gyzx(i,j,k) = ZERO 
           gyzy(i,j,k) = ZERO 
           gyzz(i,j,k) = ZERO
           gzzx(i,j,k) = ZERO 
           gzzy(i,j,k) = ZERO 
           gzzz(i,j,k) = ZERO 
           Gammaxxx(i,j,k) = ZERO
           Gammaxxy(i,j,k) = ZERO
           Gammaxxz(i,j,k) = ZERO
           Gammaxyy(i,j,k) = ZERO
           Gammaxyz(i,j,k) = ZERO
           Gammaxzz(i,j,k) = ZERO
           Gammayxx(i,j,k) = ZERO
           Gammayxy(i,j,k) = ZERO
           Gammayxz(i,j,k) = ZERO
           Gammayyy(i,j,k) = ZERO
           Gammayyz(i,j,k) = ZERO
           Gammayzz(i,j,k) = ZERO
           Gammazxx(i,j,k) = ZERO
           Gammazxy(i,j,k) = ZERO
           Gammazxz(i,j,k) = ZERO
           Gammazyy(i,j,k) = ZERO
           Gammazyz(i,j,k) = ZERO
           Gammazzz(i,j,k) = ZERO

           sbt(i,j,k)=ZERO
           sbx(i,j,k)=ZERO
           sby(i,j,k)=ZERO
           sbz(i,j,k)=ZERO
           if(em_evolve_enable==0) then
              Bx(i,j,k)=ZERO
              By(i,j,k)=ZERO
              Bz(i,j,k)=ZERO
              Ex(i,j,k)=ZERO
              Ey(i,j,k)=ZERO
              Ez(i,j,k)=ZERO
           end if
        end do
     end do
  end do
  !$omp end parallel do

  call fill_bssn_symmetry_gz_bssn_vars(ext,X,Y,Z,Symmetry,phi,chi,trK,Gammax,Gammay,Gammaz,gxx,gxy,gxz,gyy,gyz,gzz,Axx,Axy,Axz,Ayy,Ayz,Azz)

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
  !  call CCTK_SyncGroup(ierr,cctkGH,'BSSN::BSSN_vars')



  if(Symmetry==AXISYM) then
     write(*,*) "I PITY DA FOO TRIES TO USE AXISYMMETRY FOR ORBITING BLACK HOLES!"
     stop
  end if

end subroutine bhns_initialdata_local
