!|----------------------------------------------------+
!| Function to compute fluxes of EM fields 
!!   across a spherical surface of radius surf_radius.
!!   and centered at the origin
!|----------------------------------------------------+
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine em_radiative_luminosity(cctkGH,L_em_out,L_em_in,L_em_out_nosqrtmg,L_em_in_nosqrtmp, L_em_out_Poynting,L_em_out_Poynting2,L_em_out_T0r,ifile,cctktime)
  implicit none
  
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  ! Variables in the function call
  CCTK_POINTER                             :: cctkGH
  real*8                                   :: F_em_out, F_em_in, L_em_out, L_em_in, L_em_out_nosqrtmg,L_em_in_nosqrtmp, F_em_out_nosqrtmg,F_em_in_nosqrtmp, L_em_out_Poynting, F_em_out_Poynting
  real*8                                   :: F_em_out_T_0r, F_em_out_Poynting2, L_em_out_Poynting2, L_em_out_T0r
  ! Variables needed for interpolation
  real*8, dimension(N_theta*N_phi,3)       :: pointcoords
  integer                                  :: vindex, ifile
  ! Output arrays, dummy indices, parameters
!!  real*8                                   :: phi2reint(:), phi2imint(:), phi0reint(:), phi0imint(:), pointcoords(:)
  real*8, dimension(N_theta*N_phi)         :: Psi6int, lapm1int, phi2reint, phi2imint, phi0reint, phi0imint, stxint, styint, stzint
  real*8, dimension(N_theta*N_phi)         :: SPxint, SPyint, SPzint, T_0xint, T_0yint, T_0zint
  integer                                  :: i,j,n
  integer                                  :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  integer                                  :: num_cols2, header_flag2, numrows2
  real*8, dimension(N_theta*N_phi,50)      :: export_data2
  character, dimension(50)                 :: data_headers2*20
  character                                :: filename2*50
  real*8                                   :: cctktime



  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  real*8 PI, costheta, sintheta, phiangle, nx, ny, nz, F_em_x, F_em_y, F_em_z, stx, sty, stz

  PI = 3.14159265358979323846D0

  !| Psi6int <~~~ exp(6.D0*phi)
  !| lapm1int   <~~~ lapm1,
  !| phi2reint <~~~ phi2re
  !| phi2imint <~~~ phi2im
  !| phi0reint <~~~ phi0re
  !| phi0imint <~~~ phi0im
 
  n = 1
  do i=1,N_theta
     costheta = 1.0 - (i - 0.5)*dcostheta
     sintheta = sqrt(1.0 - costheta*costheta)
     do j=1,N_phi
        phiangle = (j - 0.5)*dphi
        if(N_phi==1) phiangle = 0.D0
        pointcoords(n,1) = surf_radius*sintheta*cos(phiangle)
        pointcoords(n,2) = surf_radius*sintheta*sin(phiangle)
        pointcoords(n,3) = surf_radius*costheta
        n = n + 1
     end do
  end do


  call CCTK_VarIndex(vindex,"mhd_evolve::mhd_st_x")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,stxint)

  call CCTK_VarIndex(vindex,"mhd_evolve::mhd_st_y")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,styint)

  call CCTK_VarIndex(vindex,"mhd_evolve::mhd_st_z")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,stzint)


  call CCTK_VarIndex(vindex,"em_extraction::SPx")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,SPxint)

  call CCTK_VarIndex(vindex,"em_extraction::SPy")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,SPyint)

  call CCTK_VarIndex(vindex,"em_extraction::SPz")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,SPzint)

  call CCTK_VarIndex(vindex,"em_extraction::T_0x")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,T_0xint)

  call CCTK_VarIndex(vindex,"em_extraction::T_0y")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,T_0yint)

  call CCTK_VarIndex(vindex,"em_extraction::T_0z")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,T_0zint)



  call CCTK_VarIndex(vindex,"bssn::phi")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,Psi6int)

  call CCTK_VarIndex(vindex,"lapse::lapm1")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,lapm1int)

  do i=1,N_theta*N_phi
     if (set_up_spherical_EM_wave_in_flat_spacetime==1) then
        Psi6int(i) = 1.d0
        lapm1int(i) = 1.D0
     else
        Psi6int(i) = exp(6.D0*Psi6int(i))
        lapm1int(i) = lapm1int(i) + 1.D0
     end if
  end do

  call CCTK_VarIndex(vindex,"em_extraction::NPphi2re")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,phi2reint)

  call CCTK_VarIndex(vindex,"em_extraction::NPphi2im")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,phi2imint)

  call CCTK_VarIndex(vindex,"em_extraction::NPphi0re")
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,phi0reint)

  call CCTK_VarIndex(vindex,"em_extraction::NPphi0im")                            
  call interp_driver_carp(cctkGH,N_theta*N_phi,pointcoords,vindex,phi0imint)


  ! Compute the luminosities-------------------------------------

  !Zero out our flux sum first.
  F_em_out    = 0.D0 
  F_em_in    = 0.D0 
  F_em_out_nosqrtmg = 0.d0
  F_em_in_nosqrtmp = 0.d0

  F_em_out_poynting = 0.d0

  ! Now do the surface integration
  do i=1,N_theta*N_phi
    F_em_out = F_em_out + (phi2reint(i)*phi2reint(i) + phi2imint(i)*phi2imint(i))*lapm1int(i)*Psi6int(i) 

    F_em_in = F_em_in + (phi0reint(i)*phi0reint(i) + phi0imint(i)*phi0imint(i))*lapm1int(i)*Psi6int(i) 

    F_em_out_nosqrtmg = F_em_out_nosqrtmg + (phi2reint(i)*phi2reint(i) + phi2imint(i)*phi2imint(i))
    F_em_in_nosqrtmp = F_em_in_nosqrtmp + (phi0reint(i)*phi0reint(i) + phi0imint(i)*phi0imint(i))

    ! All the above uses the Newman Penrose phi2 scalar to compute fluxes. Below we use the Poynting vector
    ! This Poynting vector diagnostic is valid only the asymptotically flat region (far from the sources)
    ! Note that the vector is not unit here, because below in Luminosity calculation we multiply by r instead of r^2.
    nx = pointcoords(i,1)
    ny = pointcoords(i,2)
    nz = pointcoords(i,3)

    F_em_out_Poynting = F_em_out_Poynting + ( nx*stxint(i) + ny*styint(i) + nz*stzint(i) ) 
    F_em_out_Poynting2 = F_em_out_Poynting2 + ( nx*SPxint(i) + ny*SPyint(i) + nz*SPzint(i) ) 
    F_em_out_T_0r = F_em_out_T_0r + ( nx*T_0xint(i) + ny*T_0yint(i) + nz*T_0zint(i) ) 
    
  end do

  ! L_em is the Poynting luminosity through the surface at surf_radius !
  L_em_out    = F_em_out *  surf_radius * surf_radius * dphi * dcostheta*sym_factor/(4.d0*PI)
  L_em_in    = F_em_in *  surf_radius * surf_radius * dphi * dcostheta*sym_factor/(4.d0*PI)
  L_em_out_nosqrtmg = F_em_out_nosqrtmg *  surf_radius * surf_radius * dphi * dcostheta*sym_factor/(4.d0*PI)
  L_em_in_nosqrtmp = F_em_in_nosqrtmp *  surf_radius * surf_radius * dphi * dcostheta*sym_factor/(4.d0*PI)

  ! And now the Poynting luminosity
  L_em_out_Poynting = F_em_out_Poynting * surf_radius * dphi * dcostheta*sym_factor
  L_em_out_Poynting2 = F_em_out_Poynting2 * surf_radius * dphi * dcostheta*sym_factor/(4.d0*PI)
  L_em_out_T0r = -F_em_out_T_0r * surf_radius * dphi * dcostheta*sym_factor


  ! ------- END OF COMPUTING LUMINOSITIES---------------------------------


  ! Output 2D data     ---------------------------------

  if(CCTKTIME == 0.D0) then
     header_flag2 = 1
  else
     header_flag2 = 0
  end if

  write(filename2,31)(ifile-1)
31 FORMAT("phi2_2D.data.",I1)

  num_cols2 = 8
  numrows2 = N_theta*N_phi

  data_headers2(1) = '# Time' 
  data_headers2(2) = "X"
  data_headers2(3) = "Y"
  data_headers2(4) = "Z"
  data_headers2(5) = 'phi2_Re' 
  data_headers2(6) = "phi2_Im"
  data_headers2(7) = "phi0_Re"
  data_headers2(8) = "phi0_Im"
  
  do i=1,N_theta*N_phi 
     export_data2(i,1) = cctktime
     export_data2(i,2) = pointcoords(i,1)
     export_data2(i,3) = pointcoords(i,2)
     export_data2(i,4) = pointcoords(i,3)
     export_data2(i,5) = phi2reint(i)
     export_data2(i,6) = phi2imint(i)
     export_data2(i,7) = phi0reint(i)
     export_data2(i,8) = phi0imint(i)
  end do

  if(CCTK_MyProc(CCTKGH)==0) call output_multiline_data_to_file(filename2,num_cols2,data_headers2,numrows2,export_data2,header_flag2) 


end subroutine em_radiative_luminosity
