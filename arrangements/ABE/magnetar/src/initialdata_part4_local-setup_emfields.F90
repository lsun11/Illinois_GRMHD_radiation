!--------------------------------------------------------
! Okay, we've read in the initial data from files.
! Now we set up all other required variables, including:
!  emfields, BSSN variables, primitives, etc.
!--------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine magnetar_setup_emfields(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  real*8                :: detmin_l,detmax_l,delta_bar
  real*8                :: Aphi_scaling_factor
  integer               :: proc_imin,proc_jmin,proc_kmin,proc_imax,proc_jmax,proc_kmax
  integer               :: ierr,index,handle,dummy,glob_imax,glob_jmax,glob_kmax
  integer               :: Nfont, Nfont_l
  CCTK_REAL             :: reduction_value

  integer               :: i,j,k

  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  integer :: ONE,ZERO
  parameter(ONE = 1.D0, ZERO = 0.D0)

  ext = cctk_lsh

  write(*,*) "part4: setup emfields..."

  ! Setup initial EM fields (from a vector potential A_phi)
  ! Note: Using rho_bl,rho_br,Pl,Pr,gxx_p,gxy_p,gxz_p,gyy_p as temporary storage!
  if(em_field_type == 0 .or. em_field_type == 2) then
     Aphi_scaling_factor = 1.D0
     call setup_poloidal_emfields(ext,RADEQUAT,magnetar_P_max,p_c,betam1,X,Y,Z,PhysicalRadius, &
          phi,lapm1,shiftx,shifty,shiftz, &
          gxx,gxy,gxz,gyy,gyz,gzz, &
          rho_b,P,Bx,By,Bz, Ax,Ay,Az,  &
          Bx_stagger, By_stagger, Bz_stagger, &
          st_x,st_y,st_z,mhd_st_x,mhd_st_y,mhd_st_z, &
          tau,u0,vx,vy,vz,Symmetry,Sym_Bz, &
          rho_bl,rho_br,Pl,Pr,Aphi_scaling_factor, & 
          constrained_transport_scheme, em_field_type, &
          gxx_p,gxy_p,gxz_p,gyy_p)

     !$omp parallel do 
     do k=1,ext(3)
	do j=1,ext(2)
	   do i=1,ext(1)
              gxx_p(i,j,k)=gxx(i,j,k) 
              gxy_p(i,j,k)=gxy(i,j,k)
              gxz_p(i,j,k)=gxz(i,j,k)
              gyy_p(i,j,k)=gyy(i,j,k)
	   end do
	end do
     end do
     !$omp end parallel do
  else 
     if(fisheye_enable==1) then
        write(*,*) "TOROIDAL FIELDS NOT FISHEYE COMPATIBLE.  PLEASE TRY AGAIN!"
        stop
     end if
     call initial_toroidal_emfields(ext,magnetar_P_max,p_c,betam1,X,Y,Z, &
          phi,lapm1, shiftx,shifty,shiftz, &
          gxx,gxy,gxz,gyy,gyz,gzz, &
          rho_b,P,Bx,By,Bz, &
          st_x,st_y,st_z,mhd_st_x,mhd_st_y,mhd_st_z, &
          tau,u0,vx,vy,vz,Symmetry,Sym_Bz)
  end if

  if(Symmetry .eq. AXISYM) then 
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_conservatives')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_primitives')
     call CartSymGN(dummy,cctkGH,'mhd_evolve::mhd_vs')

     call CCTK_VarIndex(index,'mhd_evolve::st_x')
     call BndCartoon2DVI(dummy, cctkGH, 1, index)
     call CCTK_VarIndex(index,'mhd_evolve::mhd_st_x')
     call BndCartoon2DVI(dummy, cctkGH, 1, index)
     call CCTK_VarIndex(index,'mhd_evolve::Bx')
     call BndCartoon2DVI(dummy, cctkGH, 1, index)
     call BndCartoon2DVN(dummy, cctkGH, 0, 'mhd_evolve::tau')
  end if

end subroutine magnetar_setup_emfields
