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

  integer :: vindex,handle,ierr, dummy
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


!  b2_gridfunction = 0.D0
!  b2oP_gridfunction = 0.D0

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
  else if (use_new_bhns_initial_data==1) then
     Axx=Axx*psi**(-4)
     Axy=Axy*psi**(-4)
     Axz=Axz*psi**(-4)
     Ayy=Ayy*psi**(-4)
     Ayz=Ayz*psi**(-4)
     !     Azz=-(Axx+Ayy)
     Azz=Azz*psi**(-4)
  else if (use_new_bhns_initial_data==2) then
     ! kij is the extrinsic curvature Kij. We now convert it to the BSSN traceless, conformal extrinsic curvature
     Axx=kxx*psi**(-4)-gxx*trK/3.d0
     Axy=kxy*psi**(-4)-gxy*trK/3.d0
     Axz=kxz*psi**(-4)-gxz*trK/3.d0
     Ayy=kyy*psi**(-4)-gyy*trK/3.d0
     Ayz=kyz*psi**(-4)-gyz*trK/3.d0
     Azz=kzz*psi**(-4)-gzz*trK/3.d0
     ! Also set the inverse tilde metric
     gupxx =   ( gyy * gzz - gyz * gyz )
     gupxy = - ( gxy * gzz - gyz * gxz )
     gupxz =   ( gxy * gyz - gyy * gxz )
     gupyy =   ( gxx * gzz - gxz * gxz )
     gupyz = - ( gxx * gyz - gxy * gxz )
     gupzz =   ( gxx * gyy - gxy * gxy )

     ! Sync because we are about to take derivatives to calculate Gamma^i
     call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')

     ! Use P as temporary storage here:
     call setgamma_v2(ext, cctk_nghostzones, dX, dY, dZ, &
          phi, gxx, gxy, gxz, gyy, gyz, gzz, &
          gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Gammax, Gammay, Gammaz, P, &
          Symmetry)
     
     ! Sync again
     call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  else if ((use_new_bhns_initial_data==3).or.(use_new_bhns_initial_data==4)) then

     Axx= Axx*psi**(-4)    ! The lorene/Cocal data give K_ij which assuming K=0 is K_ij = A_ij, the BSSN \tilde A_ij = \psi^{-4} A_{ij}
     Axy= Axy*psi**(-4)
     Axz= Axz*psi**(-4)
     Ayy= Ayy*psi**(-4)
     Ayz= Ayz*psi**(-4)

     Azz=-(Axx+Ayy)     
  else

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


           Ax(i,j,k)=ZERO
           Ay(i,j,k)=ZERO
           Az(i,j,k)=ZERO
           psi6phi(i,j,k)=ZERO

           Bx(i,j,k)=ZERO
           By(i,j,k)=ZERO
           Bz(i,j,k)=ZERO
           Ex(i,j,k)=ZERO
           Ey(i,j,k)=ZERO
           Ez(i,j,k)=ZERO
           mhd_st_x(i,j,k)=ZERO
           mhd_st_y(i,j,k)=ZERO
           mhd_st_z(i,j,k)=ZERO

           ! Set excision_zone_gf to avoid valgrind memory errors
           excision_zone_gf(i,j,k) = 0

           if (use_new_bhns_initial_data.ne.2) then
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
              trK(i,j,k)= 0.D0
              Gammax(i,j,k)= 0.D0
              Gammay(i,j,k)= 0.D0
              Gammaz(i,j,k)= 0.D0              
           end if
              
           shiftxt(i,j,k)= 0.D0
           shiftyt(i,j,k)= 0.D0
           shiftzt(i,j,k)= 0.D0
           lapset(i,j,k)= 0.D0
           
           
           
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

  ! Check for NaNs in the initial data
  if (genID_cmdline_output_enable.ne.1) then
     call check_for_nans(cctk_lsh,lapm1,shiftx,shifty,shiftz,gxx,gxy,gxz,gyy,gyz,gzz,Axx,Axy,Axz,Ayy,Ayz,Azz,rho_b,vx,vy,vz)
  end if


end subroutine bhns_initialdata_local


subroutine check_for_nans(ext,lapse,shiftx,shifty,shiftz,gxx,gxy,gxz,gyy,gyz,gzz,Axx,Axy,Axz,Ayy,Ayz,Azz,rhob,ux,uy,uz)
  implicit none

  integer, dimension(3)                    :: ext
  INTEGER :: i,j,k
  real*8, dimension(ext(1),ext(2),ext(3))  :: gxx,gxy,gxz,gyy,gyz,gzz,Axx,Axy,Axz,Ayy,Ayz,Azz,lapse,shiftx,shifty,shiftz,rhob,ux,uy,uz

  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
           if(isnan(lapse(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3) 
              write(*,*)  "lapse= ",lapse(i,j,k),gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
              stop
           end if

           if(isnan(shiftx(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3) 
              write(*,*)  "shiftx= ",shiftx(i,j,k),gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
              stop
           end if

           if(isnan(shifty(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3) 
              write(*,*)  "shifty= ",shifty(i,j,k),gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
              stop
           end if
           
           if(isnan(shiftz(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3) 
              write(*,*)  "shiftz= ",shiftz(i,j,k),gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
              stop
           end if

           if(isnan(gxx(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3) 
              write(*,*)  "gxx= ",gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
              stop
           end if

           if(isnan(gxy(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3) 
              write(*,*)  "gxy= ",gxy(i,j,k),gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
              stop
           end if

           if(isnan(gxz(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3) 
              write(*,*)  "gxz= ",gxz(i,j,k),gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
              stop
           end if

           if(isnan(gyy(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3) 
              write(*,*)  "gyy= ",gyy(i,j,k),gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
              stop
           end if

           if(isnan(gyz(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3) 
              write(*,*)  "gyz= ",gyz(i,j,k),gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
              stop
           end if

           if(isnan(gzz(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3) 
              write(*,*)  "gzz= ",gzz(i,j,k),gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
              stop
           end if

           if(isnan(Axx(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3) 
              write(*,*)  "Axx= ",Axx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
              stop
           end if

           if(isnan(Axy(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3) 
              write(*,*)  "Axy= ",Axy(i,j,k),gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
              stop
           end if

           if(isnan(Axz(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3) 
              write(*,*)  "Axz= ",Axz(i,j,k),gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
              stop
           end if

           if(isnan(Ayy(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3) 
              write(*,*)  "Ayy= ",Ayy(i,j,k),gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
              stop
           end if

           if(isnan(Ayz(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3) 
              write(*,*)  "Ayz= ",Ayz(i,j,k),gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
              stop
           end if

           if(isnan(Azz(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3) 
              write(*,*)  "Azz= ",Azz(i,j,k),gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
              stop
           end if

           if(isnan(rhob(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3) 
              write(*,*)  "rhob= ",rhob(i,j,k),gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
              stop
           end if

           if(isnan(ux(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3) 
              write(*,*)  "ux= ",ux(i,j,k),gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
              stop
           end if

           if(isnan(uy(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3) 
              write(*,*)  "uy= ",uy(i,j,k),gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
              stop
           end if
           
           if(isnan(uz(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3) 
              write(*,*)  "uz= ",uz(i,j,k),gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
              stop
           end if


        end do
     end do
  end do





end subroutine check_for_nans
