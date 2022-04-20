#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

! Equation numbers referenced from Baumgarte & Naculich: gr-qc/0701037

subroutine stationary_puncture_id_asymptoticBaumgarte(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext
  real*8 			           :: dT,dX,dY,dZ
  real*8 			           :: scaled_x,scaled_y,scaled_z
  real*8 			           :: n1,n2,n3,mf
  real*8 			           :: detmin_l,detmax_l

  integer :: index,dummy,ierr
  real*8  :: ONE,ZERO
  parameter(ONE = 1.D0, ZERO = 0.D0)
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  INTEGER :: i,j,k,iter,iterloop,ITERMAX=150,lo_or_hi
  REAL*8 :: xlo,xhi,r_iso,outputdiff,outputdiffold,mass,r_s,betar,dbdr,drdbr
  !  REAL*8 :: EPS=2.D-15
  !REAL*8 :: EPS=1.D-14
  REAL*8 :: EPS=1.D-12

  INTERFACE 
     FUNCTION funcv_bisect_step(xlo,xhi,fval,outputdiff,lo_or_hi)
       IMPLICIT NONE
       REAL*8, INTENT(INOUT) :: xlo,xhi,fval,outputdiff
       INTEGER :: funcv_bisect_step,lo_or_hi
     END FUNCTION funcv_bisect_step
  END INTERFACE

  ext(1) = cctk_lsh(1)
  ext(2) = cctk_lsh(2)
  ext(3) = cctk_lsh(3)  
  dT = CCTK_DELTA_TIME
  dX = X(2,1,1) - X(1,1,1)
  dY = Y(1,2,1) - Y(1,1,1)
  dZ = Z(1,1,2) - Z(1,1,1)


  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
           r_iso=PhysicalRadius(i,j,k)
           scaled_x = X(i,j,k)*r_iso/r(i,j,k)
           scaled_y = Y(i,j,k)*r_iso/r(i,j,k)
           scaled_z = Z(i,j,k)*r_iso/r(i,j,k)
           mass = 1.D0
           xlo = 3.D0*mass*0.5D0+dX*0.01
           xhi = 1000.D0

           outputdiff = 1.D0
           outputdiffold = 2.D0
           lo_or_hi = -1
           iter = 0
           iterloop = 1
           do while(outputdiff .gt. EPS)
              if(funcv_bisect_step(xlo,xhi,r_iso,outputdiff,lo_or_hi) == 1) then
                 xlo = xlo + (xhi - xlo)*0.5
              else
                 xhi = xhi - (xhi - xlo)*0.5
              end if

              if((outputdiffold - outputdiff) .lt. EPS .and. iter.gt.ITERMAX) then
                 if(iterloop.gt.10) write(*,*) "Having problems... CORRECTING!",r_iso,xlo,xhi,outputdiff
                 xlo = xlo - dX*0.8*iterloop
                 xhi = xhi + dX*iterloop
                 lo_or_hi = -1
                 iter = 0
                 iterloop = iterloop + 1
              end if

              outputdiffold = outputdiff
              iter = iter + 1
           end do

           if(MOD(k,2)==0 .and. j==1 .and. i==1) then
              index = 100.D0*k/ext(3)
              write(*,*) "Doing bisection to generate initial data:",index, "percent finished..."
           end if
           r_s = xlo
!           r_schw(i,j,k) = xlo
!           r_s = r_schw(i,j,k)

           !           write(*,*) sqrt(r_s/r_iso) - sqrt(4.D0*r_s / (2.D0*r_s + mass + sqrt(4.D0*r_s*r_s + 4.D0*mass*r_s + 3.D0*mass*mass))) * &
           !                ((8.D0*r_s + 6.D0*mass + 3.D0*sqrt(8.D0*r_s*r_s + 8.D0*mass*r_s + 6.D0*mass*mass)) / &
           !                ((4.D0+3.D0*sqrt(2.D0)) * (2.D0*r_s - 3.D0*mass)))**(1.D0/(sqrt(8.D0)))
           ! eq. 5
           psi(i,j,k) = sqrt(r_s/r_iso)
           phi(i,j,k) = log(psi(i,j,k))

           ! eq. 3b,d - 1.D0
           lapm1(i,j,k) = sqrt(1.D0 - 2.D0*mass/r_s + 27.D0*mass**4/(16.D0*r_s**4)) - 1.D0
           ! eq. 15
           betar = sqrt(27.D0)*mass*mass*r_iso/(4.D0*r_s*r_s*r_s)

           !fill in coord transformations here:
           shiftx(i,j,k) = betar*scaled_x/r_iso
           shifty(i,j,k) = betar*scaled_y/r_iso
           shiftz(i,j,k) = betar*scaled_z/r_iso

           !write(*,*) "hello.",r_iso,betar,sqrt(shiftx(i,j,k)*shiftx(i,j,k) + shifty(i,j,k)*shifty(i,j,k) + shiftz(i,j,k)*shiftz(i,j,k))

           ! eq. 9 (drdbr = dR/dr)
           drdbr=r_iso*r_s/(r_s-1.5d0)/sqrt(r_s**2+r_s+0.75d0)

           ! derivative of betar w.r.t. r
           dbdr=sqrt(27.0d0)/4.0d0*(r_s**(-3)-3.0d0*r_iso/r_s**4/drdbr)

           ! A_ij = 1/(2 alpha)* (d^i beta^j + d^j beta^i 
           !                      -2./3. delta_{ij} d_k beta^k)
           !
           ! raising and lowering is done with flat metric, 
           ! so is trivial (no factors of psi)
           ! the results is:
           !
           ! \tilde{A}_{ij}=(xhat_i xhat_j -1./3.*delta_ij)*
           !                   (d beta^r/dr-beta^r/r)/alpha
           ! 
           Axx(i,j,k)=(1.0d0*(scaled_x/r_iso)**2-1.0d0/3.0d0)*(dbdr-betar/r_iso) / (lapm1(i,j,k) + 1.D0)
           Ayy(i,j,k)=(1.0d0*(scaled_y/r_iso)**2-1.0d0/3.0d0)*(dbdr-betar/r_iso) / (lapm1(i,j,k) + 1.D0)
           Azz(i,j,k)=(1.0d0*(scaled_z/r_iso)**2-1.0d0/3.0d0)*(dbdr-betar/r_iso) / (lapm1(i,j,k) + 1.D0)
           Axy(i,j,k)=(scaled_x/r_iso)*(scaled_y/r_iso)*(dbdr-betar/r_iso)/(lapm1(i,j,k) + 1.D0)
           Axz(i,j,k)=(scaled_x/r_iso)*(scaled_z/r_iso)*(dbdr-betar/r_iso)/(lapm1(i,j,k) + 1.D0)
           Ayz(i,j,k)=(scaled_y/r_iso)*(scaled_z/r_iso)*(dbdr-betar/r_iso)/(lapm1(i,j,k) + 1.D0)
        end do
     end do
  end do

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CartSymGN(dummy,cctkGH,'lapse::lapse_vars')
  call CartSymGN(dummy,cctkGH,'shift::shift_vars')
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


  if(1==0) then
     ! Note:  kset_c uses *physical* metric as input, so convert gij to physical metric
     Sx = psi**(4.D0)
     gxx = gxx*Sx
     gxy = gxy*Sx
     gxz = gxz*Sx
     gyy = gyy*Sx
     gyz = gyz*Sx
     gzz = gzz*Sx

     if(Symmetry==AXISYM) then
        call BndCartoon2DVN(dummy, cctkGH, 0, 'bssn::phi')
        call BndCartoon2DVN(dummy, cctkGH, 0, 'bssn::trK')
        call BndCartoon2DVN(dummy, cctkGH, 0, 'lapse::lapm1')
        call CCTK_VarIndex(index,'shift::shiftx')
        call BndCartoon2DVI(dummy, cctkGH, 1, index) 
        call CCTK_VarIndex(index,'BSSN::gxx')
        call BndCartoon2DVI(dummy, cctkGH, 2, index) 
     end if

     call kset_c(ext,X,Y,Z, &
          Axx,Axy,Axz,Ayy,Ayz,Azz,trK, &
          gxx,gxy,gxz,gyy,gyz,gzz, &
          shiftx,shifty,shiftz,lapm1, &
          Symmetry)

     !  write(*,*) ""
     trK = 0.D0


     !use Sx for temporary array allocation
     Sx = psi**(-4.D0)
     Axx = Axx*Sx
     Axy = Axy*Sx
     Axz = Axz*Sx
     Ayy = Ayy*Sx
     Ayz = Ayz*Sx
     Azz = Azz*Sx
  end if

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')
  if(Symmetry .eq. AXISYM) then 
     call CCTK_VarIndex(index,'BSSN::Axx')
     call BndCartoon2DVI(dummy, cctkGH, 2, index)
  end if

  !======================================================
  ! Convert BACK to tilde metric and invert tilde metric 
  !======================================================
  gxx = ONE
  gxy = ZERO
  gxz = ZERO
  gyy = ONE
  gyz = ZERO
  gzz = ONE

  !====================================================
  ! Get Gamma^i
  !====================================================
  !Here we use gxx_p as temporary storage:
  call setgamma(ext,X,Y,Z,gxx_p, &
       phi,gxx,gxy,gxz,gyy,gyz,gzz, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
       Gammax,Gammay,Gammaz, &
       Symmetry)
  gxx_p = gxx

  call CartSymGN(dummy,cctkGH,'BSSN::BSSN_vars')
  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')
  if(Symmetry .eq. AXISYM) then 
     call CCTK_VarIndex(index,'BSSN::Gammax')
     call BndCartoon2DVI(dummy, cctkGH, 1, index)
  end if
  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')

  n1 = 0
  n2 = 0
  n3 = 0
  mf = 1
  !use gxx_p as temporary storage here:
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammax,gxx_p,n1,n2,n3,mf,Symmetry)
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammay,gxx_p,n1,n2,n3,mf,Symmetry)
  call update_falloff_fisheye(ext,X,Y,Z,PhysicalRadius,Gammaz,gxx_p,n1,n2,n3,mf,Symmetry)
  gxx_p=gxx

  !Next line NEEDED!
  call CCTK_SyncGroup(dummy,cctkGH,'BSSN::BSSN_vars')
 
end subroutine stationary_puncture_id_asymptoticBaumgarte
