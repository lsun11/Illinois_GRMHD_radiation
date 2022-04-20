#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
! PsiKadelia diagnostic output
!-----------------------------------------------------------------------------
subroutine psikad_diagnostics(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8, dimension(100)                   :: export_data
  character, dimension(100)                :: data_headers*20
  character                                :: filename*50

  !Interpolation stuff:
  real*8, allocatable, dimension(:,:)      :: pointcoords
  real*8, allocatable, dimension(:)        :: Psi4reint,Psi4imint
  real*8, allocatable, dimension(:)        :: phiint,lapm1int
  real*8, allocatable, dimension(:)        :: gupxxint,gupxyint,gupxzint
  real*8, allocatable, dimension(:)        :: gupyyint,gupyzint,gupzzint
  real*8, allocatable, dimension(:)        :: betaxint,betayint,betazint
  real*8, allocatable, dimension(:)        :: betarint
  integer                                  :: N_theta,N_phi,ntot,n,j
  real*8                                   :: sintheta,costheta,phiangle,dcostheta,PI,dphi

  real*8                                   :: Psi4relm,Psi4imlm
  real*8,dimension(1,3)                    :: pointcoords1pt
  real*8,dimension(1)                      :: outputblah
  real*8                                   :: psi4r,psi4i
  real*8                                   :: radius_GWorig
  real*8                                   :: r_areal,avg_gutt,avg_gutr,avg_gurr
  integer :: num_cols,i,l,m,radiusindex,psikad_disable
  integer, dimension(nmodes_Psi4)               :: lmode,mmode
  integer :: header_flag,vindex
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  PI = 3.14159265358979323846D0

  !  if(MOD(cctk_iteration,out_every).eq.0) then

  if(MOD(cctk_iteration,out_every).eq.0) then

     !header_flag == 1 -> write file headers
     if(CCTK_TIME == 0.D0) then
        header_flag = 1
     else
        header_flag = 0
     end if

     psikad_disable=0
     !-----------------------------------------------------------------------------------!
     ! Setup initial parameters

     if (Symmetry==AXISYM) then
        m=1
        do i=1,nmodes_Psi4
           lmode(i) = m+i
           mmode(i) = 0
        end do
     else if (Symmetry==EQUATORIAL .or. Symmetry==NO_SYMM) then
        i = 1
        l = 2
        do
           m = l
           do
              lmode(i) = l
              mmode(i) = m
              m = m - 1
              i = i + 1
              if (m .lt. -l .or. i .gt. nmodes_Psi4) exit
           end do
           if (i .gt. nmodes_Psi4) exit
           l = l+1
        end do
     else
        psikad_disable=1
        !        write(*,*) 'Symmetry not supported in Psi4 diagnostics'
        !        stop
     end if

     if(psikad_disable==0) then
        !--------------------------------------------------------------------------------------------
        ! Output spin-weight -2 decomposition of Psi4 at different radii

        do radiusindex=1,num_extraction_radii
           if (radius_GW_Psi4(radiusindex) .gt. 0.d0) then 

              if(num_extraction_radii.eq.1) then
                 filename = 'Psi4_rad.mon'
              else
                 write(filename,31)radiusindex
31               FORMAT('Psi4_rad.mon.',I1)
              endif
              num_cols = 2*nmodes_Psi4 + 5

              data_headers(1) = "# Time"
              export_data(1) = CCTK_TIME

              !| Psi4reint <~~~ psi4re
              !| Psi4imint <~~~ psi4im
              N_theta = numtheta
              N_phi = numphi

              dphi=2.D0 * PI / N_phi
              dcostheta = 2.D0 / N_theta
              if (Symmetry==AXISYM) dphi=2.D0*PI

              if (Symmetry==OCTANT) then 
                 N_theta = N_theta/2
                 N_phi = N_phi/4
              else if (Symmetry==EQUATORIAL) then
                 N_theta = N_theta/2
              else if (Symmetry==PI_SYMM) then 
                 N_theta = N_theta/2
                 N_phi = N_phi/2
              else if (Symmetry==AXISYM) then
                 !  FIXME: This only works if equatorial symmetry is also assumed.
                 N_theta = N_theta/2
                 N_phi = 1
              end if

              ntot = N_theta*N_phi

              allocate(pointcoords(ntot,3))
              allocate(lapm1int(ntot))
    	      allocate(betaxint(ntot))
              allocate(betayint(ntot))
              allocate(betazint(ntot))
              allocate(betarint(ntot))

              n = 1
              do i=1,N_theta
                 costheta = 1.D0 - (i - 0.5D0)*dcostheta
                 sintheta = sqrt(1.D0 - costheta*costheta)
                 do j=1,N_phi
                    phiangle = (j - 0.5D0)*dphi
                    if(N_phi==1) phiangle = 0.D0
                    pointcoords(n,1) = radius_GW_Psi4(radiusindex)*sintheta*cos(phiangle)
                    pointcoords(n,2) = radius_GW_Psi4(radiusindex)*sintheta*sin(phiangle)
                    pointcoords(n,3) = radius_GW_Psi4(radiusindex)*costheta
                    n = n + 1
                 end do
              end do

              call CCTK_VarIndex(vindex,"lapse::lapm1")
              call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,lapm1int)
              call CCTK_VarIndex(vindex,"shift::shiftx")
              call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,betaxint)
              call CCTK_VarIndex(vindex,"shift::shifty")
              call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,betayint)
              call CCTK_VarIndex(vindex,"shift::shiftz")
              call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,betazint)
                          
              n = 1
              do i=1,N_theta
                 do j=1,N_phi
	  	    betarint(n) = (pointcoords(n,1)*betaxint(n) + &
				  pointcoords(n,2)*betayint(n) + &
				  pointcoords(n,3)*betazint(n) )/radius_GW_Psi4(radiusindex)
                    n = n + 1
                 end do
              end do

	      deallocate(betaxint,betayint,betazint)

              allocate(phiint(ntot))
              allocate(gupxxint(ntot))
              allocate(gupxyint(ntot))
              allocate(gupxzint(ntot))
              allocate(gupyyint(ntot))
              allocate(gupyzint(ntot))
              allocate(gupzzint(ntot))

	      call CCTK_VarIndex(vindex,"bssn::phi")
              call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,phiint)
              call CCTK_VarIndex(vindex,"bssn::gupxx")
              call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gupxxint)
              call CCTK_VarIndex(vindex,"bssn::gupxy")
              call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gupxyint)
              call CCTK_VarIndex(vindex,"bssn::gupxz")
              call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gupxzint)
              call CCTK_VarIndex(vindex,"bssn::gupyy")
              call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gupyyint)
              call CCTK_VarIndex(vindex,"bssn::gupyz")
              call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gupyzint)
              call CCTK_VarIndex(vindex,"bssn::gupzz")
              call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,gupzzint)

	      ! Compute areal radius and averaged 1-lapse
	      call areal_radius_averaged_gutt_gutr_guprr(cctkGH,Symmetry,N_theta,N_phi,ntot, &
			dphi,dcostheta,radius_GW_Psi4(radiusindex), & 
			phiint,gupxxint,gupxyint,gupxzint,gupyyint,gupyzint, & 
			gupzzint, lapm1int, betarint, r_areal,avg_gutt, & 
			avg_gutr,avg_gurr)

              deallocate(phiint,gupxxint,gupxyint,gupxzint,gupyyint)
              deallocate(gupyzint,gupzzint,lapm1int,betarint)

              allocate(Psi4reint(ntot))
              allocate(Psi4imint(ntot))
              call CCTK_VarIndex(vindex,"gw_extraction::psi0re")
              call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,Psi4reint)
              call CCTK_VarIndex(vindex,"gw_extraction::psi0im")
              call interp_driver_carp(cctkGH,ntot,pointcoords,vindex,Psi4imint)

              do i=1,nmodes_Psi4
                 l = lmode(i)
                 m = mmode(i)
                 write(*,*) "outputting l mode =",l,"m mode = ",m,radius_GW_Psi4(radiusindex)
                 call lm_modes(cctkGH,Symmetry,l,m,N_theta,N_phi,ntot, &
                      dphi,dcostheta,radius_GW_Psi4(radiusindex),Psi4reint,Psi4imint,Psi4relm,Psi4imlm)
                 data_headers(2*i) = "Re(psi4_lm)"
                 data_headers(2*i+1) = "Im(psi4_lm)"
                 export_data(2*i) = Psi4relm
                 export_data(2*i+1) = Psi4imlm
              end do
	      i = nmodes_Psi4
              data_headers(2*i+2) = "r_areal"
              data_headers(2*i+3) = "Averaged g^tt"
	      data_headers(2*i+4) = "Averaged g^tr"
              data_headers(2*i+5) = "Averaged g^rr"
              export_data(2*i+2) = r_areal
              export_data(2*i+3) = avg_gutt
              export_data(2*i+4) = avg_gutr
              export_data(2*i+5) = avg_gurr

              if(CCTK_MyProc(CCTKGH).eq.0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

              deallocate(pointcoords)
              deallocate(Psi4reint)
              deallocate(Psi4imint)
          end if
        end do

     	if(num_extraction_radii.gt.0) write(*,*) 'In Psikadelia, N_theta, N_phi = ',N_theta,N_phi

        !--------------------------------------------------------------------------------------------
        ! Output Psi4 interpolated to a given point.  Useful for linearized wave analyses.

        if (enable_interp_onepoint==1) then
           filename = 'Psi4_1pt_rad.mon'

           num_cols = 3
           data_headers(1) = "# Time"
           export_data(1) = CCTK_TIME

           pointcoords1pt(1,1) = radius_GW*sin(theta_GW)*cos(phi_GW)
           pointcoords1pt(1,2) = radius_GW*sin(theta_GW)*sin(phi_GW)
           pointcoords1pt(1,3) = radius_GW*cos(theta_GW)
           call CCTK_VarIndex(vindex,"gw_extraction::psi0re")
           call interp_driver_carp(cctkGH,1,pointcoords1pt,vindex,outputblah)
           data_headers(2) = "psi4_re"
           export_data(2) = outputblah(1)
           call CCTK_VarIndex(vindex,"gw_extraction::psi0im")
           call interp_driver_carp(cctkGH,1,pointcoords1pt,vindex,outputblah)
           data_headers(3) = "psi4_im"
           export_data(3) = outputblah(1)

           if(CCTK_MyProc(CCTKGH).eq.0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
        end if
        !--------------------------------------------------------------------------------------------
     end if
  end if
end subroutine psikad_diagnostics
