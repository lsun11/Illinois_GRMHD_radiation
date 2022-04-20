#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"


subroutine ang_freq_ofBfield_driver(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: i,j,k,ierr
  integer,dimension(3) :: ones
  real*8  :: dX,dY,dZ
  real*8  :: bhx,bhy,bhz

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  !-----------------------------------------!
  !track the position of the BH
  bhx = bh_posn_x(1)
  bhy = bh_posn_y(1)
  bhz = bh_posn_z(1)

  print *,"bhx =",bhx
  print *,"bhy =",bhy
  print *,"bhz =",bhz


  if(MOD(cctk_iteration,out_every).eq.0) then
     
     write(*,*) "howdy...in ang_freq_ofBfield driver",dX,cctk_lsh

     Bfreq1 = 1.d100
     Bfreq2 = 1.d100

     
     call ang_freq_of_mag_field_lines(cctkGH,cctk_lsh,cctk_nghostzones, &
          dX,dY,dZ,r,x,y,z, &
          gxx,gxy,gxz,gyy,gyz,gzz,lapm1,shiftx,shifty,shiftz, &
          phi,psi, &
          Bx,By,Bz,Ex,Ey,Ez, Bfreq1, Bfreq2,bhx,bhy,bhz)
    
     ones = 1
     
     ! Vasilis says: I 'm not so sure about the following, but it's probably OK
     
     call CartSymGN(ierr,cctkGH,"ang_freq_ofBfield::OmegF")
     
     
     ! Finally, we set the previous level stuff
!!$     do k=cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
!!$        do j=cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
!!$           do i=cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)

!!$     do k=cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
!!$        do j=cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
!!$           do i=cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)
     !$omp parallel do
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
!!$     do k=cctk_nghostzones(3)+1,cctk_lsh(3)-cctk_nghostzones(3)
!!$        do j=cctk_nghostzones(2)+1,cctk_lsh(2)-cctk_nghostzones(2)
!!$           do i=cctk_nghostzones(1)+1,cctk_lsh(1)-cctk_nghostzones(1)

              Bfreq1_p(i,j,k) = Bfreq1(i,j,k)
              Bfreq1_p_p(i,j,k) = Bfreq1(i,j,k)

              Bfreq2_p(i,j,k) = Bfreq2(i,j,k)
              Bfreq2_p_p(i,j,k) = Bfreq2(i,j,k)

           end do
        end do
     end do
     !$omp end parallel do         

  end if

end subroutine ang_freq_ofBfield_driver
