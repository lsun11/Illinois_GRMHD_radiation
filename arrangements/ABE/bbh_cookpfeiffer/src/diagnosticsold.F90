#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!------------------------------------
! Do stationary puncture diagnostics
!------------------------------------
subroutine bbhcp_diagnostics(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,ext_noghost
  character                                :: filename*25
  real*8                                   :: dT,dX,dY,dZ,Zmin
  real*8, dimension(20)                    :: export_data
  character, dimension(20)                 :: data_headers*20
  integer                                  :: vindex
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  !Needed for GW extraction:
  real*8                                   :: Mass_GW,hplus,hcross,moncrief_gw_extr_radius
  real*8                                   :: psi4r,psi4i
  integer                                  :: jj,nmodes_moncrief
  real*8,dimension(2)                      :: gw_amplitudes
  !End GW extraction variables  
  real*8                                   :: old_surf_radius,multfactor
  CCTK_POINTER :: reduction_value_pointer
  integer :: header_flag,handle
  integer :: index,num_cols
  integer :: ierr,myproc_rank,ii

  ext = cctk_lsh

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  if(trA_detg_enforce.ne.2) then
     !First perform sanity check.  Already done if trA_detg_enforce==2!
     call sanitycheck_restore_Aij(cctkGH,cctk_lsh,gxx,gxy,gxz,gyy,gyz,gzz, &
          Axx,Axy,Axz,Ayy,Ayz,Azz)
  end if

  if(MOD(cctk_iteration,out_every)==0) then
     !header_flag == 1 -> write file headers
     if(CCTK_TIME == 0.D0) then
        header_flag = 1
     else
        header_flag = 0
     end if

     !--------------------------------------------------------------------------------!

     filename = 'statpunc.mon'
     num_cols = 12
     export_data = 0.D0
     data_headers = 'unknowndude'

     data_headers(1) = 'Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'Phicent'
     call CCTK_VarIndex (vindex, "bssn::phi")
     call interp_center_carp(cctkGH,vindex,export_data(2))

     data_headers(3) = 'Kcent'
     call CCTK_VarIndex (vindex, "bssn::trK")
     call interp_center_carp(cctkGH,vindex,export_data(3))

     data_headers(4) = 'lapsecent'
     call CCTK_VarIndex (vindex, "lapse::lapm1")
     call interp_center_carp(cctkGH,vindex,export_data(4))
     export_data(4) = export_data(4) + 1.D0

     data_headers(5) = 'psicent'
     Psi = exp(phi)
     call CCTK_VarIndex (vindex, "bssn::psi")
     call interp_center_carp(cctkGH,vindex,export_data(5))

     data_headers(6) = 'lapsemin'
     call CCTK_ReductionHandle(handle,"minimum")
     call CCTK_VarIndex(index,"lapse::lapm1")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(6),1,index)
     export_data(6) = export_data(6) + 1.D0

     data_headers(7) = 'phimax'
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex(index,"bssn::phi")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data(7),1,index)

        data_headers(9) = 'Jsurfpunc1.1M'
        old_surf_radius = surf_radius
        surf_radius=1.1D0*BigM
        call J_surf_integral_offcenter(cctkGH,0.D0,0.D0,0.D0,export_data(9))
        surf_radius = old_surf_radius

        data_headers(11) = "M_ADM1.1Surf"
        old_surf_radius = surf_radius
        surf_radius=1.1D0*BigM
        call M_surf_integral(cctkGH,export_data(11))
        surf_radius = old_surf_radius

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if
     write(*,*) "hiafter5:",vindex

     !--------------------------------------------------------------------------------!

     filename = 'statpunc.kon'
     num_cols = 4
     export_data = 0.D0

     data_headers(1) = 'Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'Xbh1'
     export_data(2) = bh_posn_x(1)

     data_headers(3) = 'Ybh1'
     export_data(3) = bh_posn_y(1)

     data_headers(4) = 'Zbh1'
     export_data(4) = bh_posn_z(1)

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if

     !--------------------------------------------------------------------------------!

     do ii=1,nsurf

        surf_radius=surfvec(ii)
        rbr=rbrvec(ii)
        drbdr=drbrvec(ii)
        ddrbddr=ddrbrvec(ii)

        if(nsurf.eq.1) then
           filename = 'statpunc.jon'
        else
           write(filename,31)ii
31         FORMAT('statpunc.jon.',I1)
        endif

        num_cols = 9

        data_headers(1) = "Time"
        export_data(1) = CCTK_TIME

        data_headers(2) = "fisheye radius"
        export_data(2) = surf_radius

        data_headers(3) = "phys radius"
        export_data(3) = surf_radius*rbr

        data_headers(4) = "Mass_sur"
        call M_surf_integral(cctkGH,export_data(4))

        data_headers(5) = "Ang_mom_sur"
        call J_surf_integral_offcenter(cctkGH,0.D0,0.D0,0.D0,export_data(5))

        data_headers(6) = "Komar Mass"
        call Komar_surf_integralv2(cctkGH,export_data(6))

        data_headers(7) = "P_x2"
        data_headers(8) = "P_y2"
        data_headers(9) = "P_z2"
        call P_surf_integralv2(cctkGH,export_data(7),export_data(8),export_data(9))

        if(CCTK_MyProc(CCTKGH) .eq. 0) then
           call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
        end if
     end do

     !--------------------------------------------------------------------------------!
     call CCTK_ReductionHandle(handle,"minimum")
     call CCTK_VarIndex(index,"grid::Z")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,Zmin,1,index)

     do jj=1,moncrief_gw_num_radii

        moncrief_gw_extr_radius = moncrief_radius_GW(jj)        

        ! Gravitational wave extraction!
        Mass_GW = 0.0D0

        ! Initialize several variables at t=0
        if (cctk_iteration==0) then
           momentsr_old = 0.d0
           momentsi_old = 0.d0
           odd_momentsr_old = 0.d0
           odd_momentsi_old = 0.d0
           int_momentsr = 0.d0
           int_momentsi = 0.d0
           momentsr = 0.D0
           momentsi = 0.D0
        end if

        !gw_amplitudes = 0.D0
        hplus = 0.d0
        hcross = 0.d0

        nmodes_moncrief = ubound(momentsr,1)

        call gw_wave_flux(cctkGH,moncrief_gw_extr_radius,Mass_GW,cctk_lsh, &
             dT,phi,gxx,gxy,gxz,gyy,gyz,gzz,momentsr,momentsi, &
             momentsr_old,momentsi_old,odd_momentsr,odd_momentsi, &
             odd_momentsr_old,odd_momentsi_old,int_momentsr,int_momentsi, &
             hplus,hcross,theta_GW,phi_GW, &
             PhysicalRadius,Radiusderivative,Zmin,dZ)

        if(moncrief_gw_num_radii.eq.1) then
           filename = 'statpunc-evenmodes.wave'
        else
           write(filename,41)jj
41         FORMAT('statpunc-evenmodes.wave.',I1)
        endif

        num_cols = 2*nmodes_moncrief + 1

        data_headers(1) = "# Time"
        export_data(1) = CCTK_TIME

        do ii=1,nmodes_moncrief
           data_headers(2*ii) = " " 
           data_headers(2*ii+1) = " " 
           export_data(2*ii) = momentsr(ii)
           export_data(2*ii+1) = momentsi(ii)
        end do

        if(CCTK_MyProc(CCTKGH) .eq. 0) then
           call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
        end if
        !--------------------------------------------------------------------------------!

        if(moncrief_gw_num_radii.eq.1) then
           filename = 'statpunc-oddmodes.wave'
        else
           write(filename,51)jj
51         FORMAT('statpunc-oddmodes.wave.',I1)
        endif
        num_cols = 2*nmodes_moncrief + 1

        data_headers(1) = "# Time"
        export_data(1) = CCTK_TIME

        do ii=1,nmodes_moncrief
           data_headers(2*ii) = " "
           data_headers(2*ii+1) = " "
           export_data(2*ii) = odd_momentsr(ii)
           export_data(2*ii+1) = odd_momentsi(ii)
        end do

        if(CCTK_MyProc(CCTKGH) .eq. 0) then
           call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
        end if

     end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end if

!  PsiTau = (shiftx*lapsex + shifty*lapsey + shiftz*lapsez)/(2.D0 * (1.D0 + lapm1)) - trK

  write(*,*) "FINISHED WITH bbhcp_diagnostics!"

end subroutine bbhcp_diagnostics
