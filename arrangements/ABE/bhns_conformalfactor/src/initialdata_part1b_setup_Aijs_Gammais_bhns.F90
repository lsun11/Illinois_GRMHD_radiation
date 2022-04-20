!--------------------------------------------------------
! Okay, we've read in the initial data from files.
! Now we set up all other required variables, including:
!  emfields, BSSN variables, primitives, etc.
!--------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine magnetar_setup_Aijs_Gammais_bhns(CCTK_ARGUMENTS)
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

  !$omp parallel do
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)

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
           
           sbt(i,j,k)=ZERO
           sbx(i,j,k)=ZERO
           sby(i,j,k)=ZERO
           sbz(i,j,k)=ZERO
        end do
     end do
  end do
  !$omp end parallel do
  
  call fill_bssn_symmetry_gz_bssn_vars(ext,X,Y,Z,Symmetry,phi,chi,trK,Gammax,Gammay,Gammaz,gxx,gxy,gxz,gyy,gyz,gzz,Axx,Axy,Axz,Ayy,Ayz,Azz)  


  if(Symmetry==AXISYM) then
     write(*,*) "I PITY DA FOO TRIES TO USE AXISYMMETRY FOR ORBITING BLACK HOLES!"
     stop
  end if

  ! Check for NaNs in the initial data
  if (genID_cmdline_output_enable.ne.1) then
     call check_for_nans(cctk_lsh,lapm1,shiftx,shifty,shiftz,gxx,gxy,gxz,gyy,gyz,gzz,Axx,Axy,Axz,Ayy,Ayz,Azz,rho_b,vx,vy,vz)
  end if

end subroutine magnetar_setup_Aijs_Gammais_bhns

 





subroutine check_for_nans_bhns(ext,lapse,shiftx,shifty,shiftz,gxx,gxy,gxz,gyy,gyz,gzz,Axx,Axy,Axz,Ayy,Ayz,Azz,rhob,ux,uy,uz)
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





end subroutine check_for_nans_bhns
