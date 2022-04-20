!--------------------------------------------
! Driver routine for updating outer boundary
!--------------------------------------------
#include "cctk.h"
#include "GenericFD.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine driver_bssn_update_boundary(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8                  :: dT, dX,dY,dZ,chi_floor,levelnumber
  integer, dimension(3)   :: ext,fake_ext,thisproc_have_global_bdry_min,thisproc_have_global_bdry_max
  integer                 :: i,ierr,handle,vindex,j,k
  integer, parameter      :: AXISYM = 4
  !
  dT = CCTK_DELTA_TIME
  dX=CCTK_DELTA_SPACE(1)
  dY=CCTK_DELTA_SPACE(2)
  dZ=CCTK_DELTA_SPACE(3)

  ext = cctk_lsh




  levelnumber = cctk_levfac(1)
  levelnumber = log(levelnumber)/log(2.D0)+1.D0

  if(cctk_iteration.gt.0) then
     if(chi_evolution_enable==1) then
        if(chi_floor_enable==1) then
           !Make sure chi remains positive!
           !This is in the same vein as the floor imposed by the Rochester (formerly Brownsville) 
           !    group.  See: Phys. Rev. Lett. 96, 111101 (2006).
           ! We use this floor in tandem with the floor imposed in compute_rhsnew.
           chi_floor = exp(chi_exponent*phi_cap);

           write(*,*) "chi floor: ",chi_floor,log(chi_floor)/chi_exponent
           !$omp parallel do
           do k=1,ext(3)
              do j=1,ext(2)
                 do i=1,ext(1)
                    if(chi(i,j,k).lt.0.D0) chi(i,j,k) = chi_floor*1.D-2
                 end do
              end do
           end do
           !$omp end parallel do
        end if
        !$omp parallel do
        do k=1,ext(3)
           do j=1,ext(2)
              do i=1,ext(1)
                 phi(i,j,k) = log(chi(i,j,k))/chi_exponent
                 psi(i,j,k) = exp(phi(i,j,k))
              end do
           end do
        end do
        !$omp end parallel do
     end if

     if(phi_cap_enable==1) then
        !$omp parallel do
        do k=1,ext(3)
           do j=1,ext(2)
              do i=1,ext(1)
                 if(phi(i,j,k).gt.phi_cap) then
                    phi(i,j,k) = phi_cap
                    psi(i,j,k) = exp(phi(i,j,k))
                 end if
              end do
           end do
        end do
        !$omp end parallel do
     end if


     if(cctk_iteration.gt.0 .and. levelnumber==1.D0) then
        !thisproc_have_global_bdry_min = have_global_bdry_min(int(levelnumber),:)
        !thisproc_have_global_bdry_max = have_global_bdry_max(int(levelnumber),:)

        ! This is a quick and dirty fix: set thisproc_have_global_bdry_min/max so that
        !    ALL boundaries are updated, EXCEPT in symmetry ghostzones!
        !    We must apply this fix, since have_global_bdry_min/max does not save all the
        !    necessary information when being saved to checkpoint.
        !    In other words, have_global_bdry_min/max is BROKEN.
        call global_edges(dX,dY,dZ, &
             X(1,1,1),Y(1,1,1),Z(1,1,1), &
             X(1,1,1),Y(1,1,1),Z(1,1,1), &
             X(cctk_lsh(1),1,1),Y(1,cctk_lsh(2),1),Z(1,1,cctk_lsh(3)), &
             X(cctk_lsh(1),1,1),Y(1,cctk_lsh(2),1),Z(1,1,cctk_lsh(3)), &
             thisproc_have_global_bdry_min,thisproc_have_global_bdry_max,Symmetry)

        write(*,*) "EDGES:",thisproc_have_global_bdry_min,thisproc_have_global_bdry_max

        do i=1,cctk_nghostzones(2)
           fake_ext = cctk_lsh - cctk_nghostzones + i
           if(Symmetry==AXISYM) fake_ext(1) = fake_ext(1) + 1
           call update_boundary(ext,fake_ext,X,Y,Z,r,PhysicalRadius, RadiusDerivative,fisheye_enable,&
                dT,dX,dY,dZ,lapm1_p,phi_p,phi,trK_p,trK,  &
                gxx_p,gxy_p,gxz_p,gyy_p,gyz_p,gzz_p, &
                gxx,gxy,gxz,gyy,gyz,gzz, &
                Axx_p,Axy_p,Axz_p,Ayy_p,Ayz_p,Azz_p, &
                Axx,Axy,Axz,Ayy,Ayz,Azz, &
                Gammax_p,Gammay_p,Gammaz_p, &
                Gammax,Gammay,Gammaz,Symmetry, &
                Spatial_Gauge,Xcenter, Ycenter, Zcenter,bc_type, &
                gambound,adm_mass,thisproc_have_global_bdry_min,thisproc_have_global_bdry_max,use_trans_fish_phys_new)
        end do

        !need Sync's for axisymmetric, non-Carpet runs
        if(Symmetry == AXISYM) then
           if(excision_enable.ne.1) then
              call BndCartoon2DVN(ierr, cctkGH, 0, 'BSSN::phi')
              call BndCartoon2DVN(ierr, cctkGH, 0, 'BSSN::trK')
              call CCTK_VarIndex(vindex,'BSSN::gxx')
              call BndCartoon2DVI(ierr, cctkGH, 2, vindex)
              call CCTK_VarIndex(vindex,'BSSN::Axx')
              call BndCartoon2DVI(ierr, cctkGH, 2, vindex)
              call CCTK_VarIndex(vindex,'BSSN::Gammax')
              call BndCartoon2DVI(ierr, cctkGH, 1, vindex)
           else
              call CCTK_SyncGroup(ierr,cctkGH,'BSSN::BSSN_vars')
              call ezaxibc_scalar(ext,x,y,z,phi,excision_zone_gf,Symmetry)
              call ezaxibc_scalar(ext,x,y,z,trK,excision_zone_gf,Symmetry)
              call ezaxibc_vector(ext,x,y,z,Gammax,Gammay,Gammaz,excision_zone_gf,Symmetry)
              call ezaxibc_tensor(ext,x,y,z,gxx,gxy,gxz,gyy,gyz,gzz,excision_zone_gf,Symmetry)
              call ezaxibc_tensor(ext,x,y,z,Axx,Axy,Axz,Ayy,Ayz,Azz,excision_zone_gf,Symmetry)
           end if
           call CCTK_SyncGroup(ierr,cctkGH,'BSSN::BSSN_vars')
        end if

     end if

     !$omp parallel do
     do k=1,ext(3)
        do j=1,ext(2)
           do i=1,ext(1)
              chi(i,j,k) = exp(chi_exponent*phi(i,j,k))

           end do
        end do
     end do
     !$omp end parallel do
     !     chi = exp(chi_exponent*phi)
  end if



end subroutine driver_bssn_update_boundary
