!------------------------------------------------------------------------------
! Driver routine for updating outer boundary, computing Ricci, and constraints
!------------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine driver_post_update_boundary(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8                  :: dX,dY,dZ,psi4
  real*8                  :: radius_old, rbr_old, drbdr_old, ddrbddr_old
  integer, dimension(3)   :: ext
  integer                 :: i,j,k
  !
  dX=CCTK_DELTA_SPACE(1)
  dY=CCTK_DELTA_SPACE(2)
  dZ=CCTK_DELTA_SPACE(3)
  
  ext = cctk_lsh

! write(*,*) "Start bssn_postbc!!!!!!"

if (1==0) then
  do k=1, ext(3)
     do j=1, ext(2)
     do i=1, ext(1)
       if (ABS(phi(i,j,k)).ge.10.D0) then
       write(*,*) "Inside bssn driver_post_bcs, phi(i,j,k) is unphysical, NOT try to fix it", i,j,k, phi(i,j,k), phi(i-1,j,k), phi(i,j-1,k),phi(i,j,k-1)
       end if
       end do
     end do
  end do



  do k=1, ext(3)
     do j=1, ext(2)
        do i=1, ext(1)

if (i==27.and.j==24.and.k==16) then
   write(*,*) "1.Inside bssn driver_post_bcs, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==25.and.j==14.and.k==19) then
   write(*,*) "2.Inside bssn driver_post_bcs, phi(i,j,k) is ", phi(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
   write(*,*) "Inside bssn driver_post_bcs, tau(i,j,k) is ", tau(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
   write(*,*) "Inside bssn driver_post_bcs, tau_rhs(i,j,k) is ", tau_rhs(i,j,k)
end if


if (i==50.and.j==2.and.k==2) then
   write(*,*) "Inside bssn driver_post_bcs, tau_rad(i,j,k) is ", tau_rad(i,j,k)
end if

if (i==50.and.j==2.and.k==2) then
   write(*,*) "Inside bssn driver_post_bcs, tau_rad_rhs(i,j,k) is ", tau_rad_rhs(i,j,k)
end if


	 end do
    end do
end do

end if 
  
  if(cctk_iteration.gt.0) then
     !Note: already computed these at the end of driver_ricci...
     !     call Derivs(ext,X,Y,Z,dX,dY,dZ,phi,phix,phiy,phiz,Symmetry)
     !     call Derivs_interior(cctkGH,dx,dy,dz,cctk_nghostzones,cctk_lsh,Symmetry,phi,phix,phiy,phiz)
     
     !Finally, compute K_ij's for apparent horizon finder and PsiKadelia:
     !$omp parallel do private(psi4)
     do k=1,ext(3)
        do j=1,ext(2)
           do i=1,ext(1)
              Psi(i,j,k) = exp(phi(i,j,k))
              psi4 = Psi(i,j,k)*Psi(i,j,k)*Psi(i,j,k)*Psi(i,j,k)
              Kxx(i,j,k) = psi4 * (Axx(i,j,k) + (1.D0/3.D0) * gxx(i,j,k) * trK(i,j,k))
              Kxy(i,j,k) = psi4 * (Axy(i,j,k) + (1.D0/3.D0) * gxy(i,j,k) * trK(i,j,k))
              Kxz(i,j,k) = psi4 * (Axz(i,j,k) + (1.D0/3.D0) * gxz(i,j,k) * trK(i,j,k))
              Kyy(i,j,k) = psi4 * (Ayy(i,j,k) + (1.D0/3.D0) * gyy(i,j,k) * trK(i,j,k))
              Kyz(i,j,k) = psi4 * (Ayz(i,j,k) + (1.D0/3.D0) * gyz(i,j,k) * trK(i,j,k))
              Kzz(i,j,k) = psi4 * (Azz(i,j,k) + (1.D0/3.D0) * gzz(i,j,k) * trK(i,j,k))

              if(isnan(Psi(i,j,k))) then
                 write(*,*) "Psi:driver_post_update_boundary: FOUND NAN AT ", i,cctk_lsh(1),j,cctk_lsh(2),k,cctk_lsh(3) 
                 write(*,*) "dx =", dX, " \txyz = ", CCTK_TIME, x(i,j,k), y(i,j,k), z(i,j,k), "gij_old= ",gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
                 write(*,*) "Psi=", Psi(i,j,k)
		 stop
              end if
              if(isnan(Kxx(i,j,k))) then
                 write(*,*) "Psi:driver_post_update_boundary: FOUND NAN AT ", i,cctk_lsh(1),j,cctk_lsh(2),k,cctk_lsh(3) 
                 write(*,*) "dx =", dX, " \txyz = ", CCTK_TIME, x(i,j,k), y(i,j,k), z(i,j,k), "gij_old= ",gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
                 write(*,*) "Kxx, psi4, Axx, gxx, trK", Kxx(i,j,k), psi4, Axx(i,j,k), gxx(i,j,k), trK(i,j,k)
		 stop
              end if
              if(isnan(Kxy(i,j,k))) then
                 write(*,*) "Psi:driver_post_update_boundary: FOUND NAN AT ", i,cctk_lsh(1),j,cctk_lsh(2),k,cctk_lsh(3) 
                 write(*,*) "dx =", dX, " \txyz = ", CCTK_TIME, x(i,j,k), y(i,j,k), z(i,j,k), "gij_old= ",gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
                 write(*,*) "Kxy=", Kxy(i,j,k)
		 stop
              end if
              if(isnan(Kxz(i,j,k))) then
                 write(*,*) "Psi:driver_post_update_boundary: FOUND NAN AT ", i,cctk_lsh(1),j,cctk_lsh(2),k,cctk_lsh(3) 
                 write(*,*) "dx =", dX, " \txyz = ", CCTK_TIME, x(i,j,k), y(i,j,k), z(i,j,k), "gij_old= ",gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
                 write(*,*) "Kxz=", Kxz(i,j,k)
		 stop
              end if
              if(isnan(Kyy(i,j,k))) then
                 write(*,*) "Psi:driver_post_update_boundary: FOUND NAN AT ", i,cctk_lsh(1),j,cctk_lsh(2),k,cctk_lsh(3) 
                 write(*,*) "dx =", dX, " \txyz = ", CCTK_TIME, x(i,j,k), y(i,j,k), z(i,j,k), "gij_old= ",gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
                 write(*,*) "Kyy=", Kyy(i,j,k)
		 stop
              end if
              if(isnan(Kyz(i,j,k))) then
                 write(*,*) "Psi:driver_post_update_boundary: FOUND NAN AT ", i,cctk_lsh(1),j,cctk_lsh(2),k,cctk_lsh(3) 
                 write(*,*) "dx =", dX, " \txyz = ", CCTK_TIME, x(i,j,k), y(i,j,k), z(i,j,k), "gij_old= ",gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
                 write(*,*) "Kyz=", Kyz(i,j,k)
		 stop
              end if
              if(isnan(Kzz(i,j,k))) then
                 write(*,*) "Psi:driver_post_update_boundary: FOUND NAN AT ", i,cctk_lsh(1),j,cctk_lsh(2),k,cctk_lsh(3) 
                 write(*,*) "dx =", dX, " \txyz = ", CCTK_TIME, x(i,j,k), y(i,j,k), z(i,j,k), "gij_old= ",gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
                 write(*,*) "Kzz=", Kzz(i,j,k)
		 stop
              end if
           end do
        end do
     end do
     !$omp end parallel do
  end if
  
  ! The case bc_type==7 requires that we compute the ADM mass:
  ! I.e., bc_type==7 imposes phi = M_ADM/2r + O(1/r^2) at the outer boundary
  if (bc_type==7 .and. iter_count==number_of_mol_ministeps .and. cctk_iteration.gt.0) then
     radius_old = surf_radius
     rbr_old = rbr
     drbdr_old = drbdr
     ddrbddr_old = ddrbddr

     surf_radius = radius_adm
     rbr = rbr_adm
     drbdr = drbdr_adm
     ddrbddr = ddrbddr_adm

     call M_surf_integral(cctkGH,adm_mass)

     surf_radius = radius_old
     rbr = rbr_old
     drbdr = drbdr_old
     ddrbddr = ddrbddr_old
     !call surf_Mass3_fisheye(cctkGH,ext,cctk_gsh,cctk_nghostzones,X,Y,Z,dX,dY,dZ,ntot, &
     !     radius_adm,N_theta,N_phi,Xlocal1d,Ylocal1d,Zlocal1d, &
     !     Xglobal,Yglobal,Zglobal,phi,gupxx,gupxy,gupxz,gupyy, &
     !     gupyz,gupzz,phix,phiy,phiz,Gammax,Gammay,Gammaz,rbr_adm,drbdr_adm, &
     !     ddrbddr_adm,adm_mass,sym_factor,dcostheta,dphi)
     
     !write(*,*) "Computed ADM Mass for this system:",adm_mass,rbr_adm,drbdr_adm,ddrbddr_adm,N_theta,N_phi 
  end if

end subroutine driver_post_update_boundary
