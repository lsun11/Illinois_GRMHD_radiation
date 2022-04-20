#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
! Diagnostic output for bondi thorn
!-----------------------------------------------------------------------------
subroutine OS_diagnostics(CCTK_ARGUMENTS)
  
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  integer                                  :: ii,jj,vindex,num_cols,header_flag
  real*8                                   :: rho_analytic
  real*8                                   :: dT,lapse_center,rho_center,rho_center_analytic,P_center,trK_center,Q_analytic
  character                                :: filename*50
  real*8, dimension(50)                    :: export_data
  character, dimension(50)                 :: data_headers*20
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !particle tracer stuff for analytic solutions

 ! if(CCTK_MyProc(CCTKGH) .eq. 0) then
 !    do ii=1,3
 !       do jj=1,narr
 !          pos(jj,ii) = coord(jj,ii+1)
 !       end do
 !    end do  
     
 ! endif
    
  
 ! filename = 'analytic_profiles.xg'
 ! if (CCTK_TIME .eq. 0) then
 !    open(UNIT=15,FILE=filename,STATUS="REPLACE")
 !    write(15,*) ""
 !    write(15,*) '"analytic profiles'
 ! else
 !    open(UNIT=15,FILE=filename,STATUS="OLD",POSITION="APPEND")
 ! endif
  

  dT = CCTK_DELTA_TIME

 ! if (CCTK_TIME .gt. 0) then
 !    call CCTK_VarIndex (vindex, "lapse::lapm1")
 !    call interp_center_carp(cctkGH,vindex,lapse_center)
 !    call CCTK_VarIndex (vindex, "mhd_evolve::rho_b")
 !    call interp_center_carp(cctkGH,vindex,rho_center)
     
 !    lapse_center=lapse_center+1.d0
 !    tau_center_OS = tau_center_OS+dT*lapse_center
 !    call OS_analytic(tau_center_OS,rho_center_analytic,R_edge,1.d0)
 !    write(*,*) "dT: ",dT
 !    write(*,*) "lapse_center: ",lapse_center
 !    write(*,*) "rho_center: ", rho_center
 !    write(*,*) "tau_center_OS: ",tau_center_OS
 !    write(*,*) "rho_center_analytic: ",rho_center_analytic
 ! endif
 ! write(*,*) "t: ",CCTK_TIME, " tau_center: ",tau_center_OS
  
 ! if (mod(CCTK_iteration,out_every) .eq. 0) then
 !    write(15,*) ""
 !    write(15,*) ""
 !    write(15,*) '"Time =',CCTK_TIME
     
     ! write(15,'(t2,a6,t17,a1,t31,a5,t45,a1,t59,a1,t73,a1,t87,a2,t101,a2)') '"coord','z','rho_b','P','E','F','Ec','Fc'
  !   do ii=1,narr
        !           write(15,'(7e18.10)') coord(ii,2),z_anal(ii),rho_b_anal(ii)*rho_b_initial,&
        !                rho_b_anal(ii)*E_anal(ii)**0.25d0*PoRho*rho_b_initial,&
        !                E_anal(ii)*rho_b_initial*PoRho/Po4PiB,&
        !                F_anal(ii)*rho_b_initial*PoRho/Po4PiB,&
        !                Ec_anal(ii)*rho_b_initial*PoRho/Po4PiB,Fc_anal(ii)*rho_b_initial*PoRho/Po4PiB
!        call OS_analytic(coord(ii,1),rho_analytic,R_edge,1.d0) 
!        write(15,'(2e18.10)') coord(ii,2),rho_analytic
  !   end do
 ! endif
 ! close(15)
  

  write(*,*) "before interp_center"
  !--------------------------------------------------------------------------------!
  call CCTK_VarIndex (vindex, "lapse::lapm1")
  call interp_center_carp(cctkGH,vindex,lapse_center)
  call CCTK_VarIndex (vindex, "mhd_evolve::rho_b")
  call interp_center_carp(cctkGH,vindex,rho_center)
  call CCTK_VarIndex (vindex, "mhd_evolve::P") 
  call interp_center_carp(cctkGH,vindex,P_center)
  call CCTK_VarIndex (vindex, "bssn::trK")
  call interp_center_carp(cctkGH,vindex,trK_center)
  write(*,*) "after interp_center"

  lapse_center=lapse_center+1.d0
  if (CCTK_TIME .gt. 0) then 
     tau_center_OS = tau_center_OS+dT*lapse_center
  endif
  if (MOD(cctk_iteration,out_every)==0) then  

     !header_flag == 1 -> write file headers
     if(CCTK_TIME == 0.D0) then 
        header_flag = 1 
     else
        header_flag = 0 
     end if
 
     call OS_analytic(tau_center_OS,rho_center_analytic,Q_analytic,R_edge,1.d0)
     
     filename = 'OS_collapse.mon'
     num_cols = 9
     export_data = 0.d0
     data_headers(1) = '# Time' 
     export_data(1) = CCTK_TIME 
     data_headers(2) = 'tau_center'
     export_data(2) = tau_center_os
     data_headers(3) = 'lapse_center'
     export_data(3) = lapse_center
     data_headers(4) = 'rho_center'
     export_data(4) = rho_center
     data_headers(5) = 'rho_center_analytic'
     export_data(5) = rho_center_analytic
     data_headers(6) = 'Q_analytic'
     export_data(6) = Q_analytic
     data_headers(7) = 'baumgarte_analytic_lapse'
     export_data(7) = 1.d0+6.d0*log(Q_analytic)
     data_headers(8) = 'P_center'
     export_data(8) = P_center
     data_headers(9) = 'trK_center'
     export_data(9) = trK_center
     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  endif
  !--------------------------------------------------------------------------------!
     do ii=1,nsurf

        surf_radius=surfvec(ii)
        rbr=rbrvec(ii)
        drbdr=drbrvec(ii)
        ddrbddr=ddrbrvec(ii)

        if(nsurf.eq.1) then
           filename = 'OS_collapse.jon'
        else
           write(filename,31)ii
31         FORMAT("OS_collapse.jon.",I1)
        endif

        num_cols = 9

        data_headers(1) = '# Time'
        export_data(1) = CCTK_TIME

        data_headers(2) = "fisheye radius"
        export_data(2) = surf_radius

        data_headers(3) = "phys radius"
        export_data(3) = surf_radius*rbr

        data_headers(4) = "Mass_sur"
        call M_surf_integral(cctkGH,export_data(4))

        data_headers(5) = "Ang_mom_surf"
        call J_surf_integral_offcenter(cctkGH,0.D0,0.D0,0.D0,export_data(5))

        data_headers(6) = "Komar Mass"
        call Komar_surf_integralv2(cctkGH,export_data(6))

        data_headers(7) = "P_x2"
        data_headers(8) = "P_y2"
        data_headers(9) = "P_z2"
        call P_surf_integralv2(cctkGH,export_data(7),export_data(8),export_data(9))

        if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end do

     !--------------------------------------------------------------------------------!

end subroutine OS_diagnostics
