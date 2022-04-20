#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!------------------------------------
! Do stationary puncture diagnostics
!------------------------------------
subroutine kerr_diagnostics(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,ext_noghost
  character                                :: filename*50
  real*8                                   :: dT,dX,dY,dZ,Zmin
  real*8, dimension(50)                    :: export_data
  character, dimension(50)                 :: data_headers*20
  integer                                  :: vindex
  integer				   :: npts
  real*8, allocatable, dimension(:,:)      :: pointcoords
  real*8, allocatable, dimension(:)	   :: phiint
  real*8 				   :: xstart,xend,d_spacing
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  !Needed for GW extraction:
  real*8                                   :: Mass_GW,hplus,hcross,moncrief_gw_extr_radius
  real*8                                   :: psi4r,psi4i
  integer                                  :: jj,nmodes_moncrief
  real*8,dimension(2)                      :: gw_amplitudes
  !End GW extraction variables  

  real*8                                   :: old_surf_radius,SymmFactor
  CCTK_POINTER :: reduction_value_pointer
  integer :: header_flag,handle
  integer :: index,num_cols
  integer :: ierr,myproc_rank,ii
  ! Variables used for BH local diagnotsics:
  integer :: found_horizon, horizon_number
  integer :: N_theta_bh_diagnostics, N_phi_bh_diagnostics
  real*8  :: Mirr, M_BH, J_BH, xi_err

  ext = cctk_lsh

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  !--------------------------------------------------------------------------------!
  ! ******** TEST ****************
  ! Update movingbox center
  !---------------------------------------------------------------------------------
  !!if (cctk_time .gt. 140.d0) stop
  !!if (cctk_iteration .gt. 0) then 
  !!   position_y(1) = 0.d0
  !!   position_z(1) = 0.d0
  !!   ii = 2
  !!   do while (cctk_time .gt. t_tab(ii))
  !!      ii = ii+1
  !!   end do
  !!   position_x(1) = bh_x_tab(ii-1) + (cctk_time-t_tab(ii-1))/(t_tab(ii)-t_tab(ii-1)) * (bh_x_tab(ii)-bh_x_tab(ii-1))
  !!end if
! !***************************

  if(trA_detg_enforce.ne.2) then
     !First perform sanity check.  Already done if trA_detg_enforce==2!
     call sanitycheck_restore_Aij(cctkGH,cctk_lsh,gxx,gxy,gxz,gyy,gyz,gzz, &
          Axx,Axy,Axz,Ayy,Ayz,Azz)
  end if

  if(MOD(cctk_iteration,out_every)==0) then
     call set_SymmFactor(Symmetry,SymmFactor)

     !header_flag == 1 -> write file headers
     if(CCTK_TIME == 0.D0) then
        header_flag = 1
     else
        header_flag = 0
     end if

     !--------------------------------------------------------------------------------!

     filename = 'ks.mon'
     num_cols = 13
     export_data = 0.D0
     data_headers = 'unknowndude'

     data_headers(1) = '# Time'
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

     data_headers(8) = 'Jsurf_inner'
     old_surf_radius = surf_radius
     surf_radius=BH_Vol_Excise_Radius
     call J_surf_integral_offcenter(cctkGH,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1),export_data(8))
     surf_radius = old_surf_radius

     ! J_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(9) = "J_exVol"
     export_data(9) = J_ADM_VolInt*SymmFactor

     data_headers(10) = "J_Surf+Vol"
     export_data(10) = J_ADM_VolInt*SymmFactor + export_data(8)

     data_headers(11) = "M_ADM_inner_surf"
     old_surf_radius = surf_radius
     surf_radius=BH_Vol_Excise_Radius
     call M_surf_integral_offcenter_nofish(cctkGH,bh_posn_x(1),bh_posn_y(1),bh_posn_z(1),export_data(11))
     surf_radius = old_surf_radius

     ! M_ADM volume integral is computed in scheduler of diagnostics_vacuum
     data_headers(12) = "M_ADM_exVol"
     export_data(12) = M_ADM_VolInt*SymmFactor

     data_headers(13) = "M_Surf_+_Vol"
     export_data(13) = M_ADM_VolInt*SymmFactor + export_data(11)

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     write(*,*) "hiafter5:",vindex

     !--------------------------------------------------------------------------------!
     ! BH local diagnostics
     !
     filename = 'ks_BH_diagnostics.mon'
     num_cols = 5
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     N_theta_bh_diagnostics = 50
     N_phi_bh_diagnostics = 50

     Mirr = 0.d0
     M_BH = 0.d0
     J_BH = 0.d0
     xi_err = 0.d0

     horizon_number = 1
     if (Symmetry==EQUATORIAL) then
        call BH_local_diagnostics(cctkGH,Mirr,M_BH,J_BH, xi_err, &
                       horizon_number,N_theta_bh_diagnostics, &
                       N_phi_bh_diagnostics,Symmetry,found_horizon)
     else 
        Mirr = 0.d0
	M_BH = 0.d0
	J_BH = 0.d0
	xi_err = 1.9718d300
     end if

     data_headers(2) = 'Mirr_bh1'
     export_data(2) = Mirr

     data_headers(3) = 'M_bh1'
     export_data(3) = M_BH

     data_headers(4) = 'J_bh1'
     export_data(4) = J_BH

     data_headers(5) = 'xi_err1'
     export_data(5) = xi_err

     if(CCTK_MyProc(CCTKGH)==0 .and. found_horizon==1)  &
       call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!

     filename = 'ks.kon'
     num_cols = 4
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME
     data_headers(2) = 'Xbh1'
     data_headers(3) = 'Ybh1'
     data_headers(4) = 'Zbh1'
     export_data(2) = bh_posn_x(1)
     export_data(3) = bh_posn_y(1)
     export_data(4) = bh_posn_z(1)

     if(CCTK_MyProc(CCTKGH)==0) & 
       call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     !************ HAMILTONIAN CONSTRAINT ************!

     filename = 'ks-ham.con'
     num_cols = 5
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = 'Haml ResN'
     export_data(2) = Ham_const_VolIntN*SymmFactor

     data_headers(3) = 'Haml ResD'
     export_data(3) = Ham_const_VolIntD*SymmFactor

     data_headers(4) = "Ham Res2P1.1N"
     export_data(4) = Ham_const_excised_VolIntN*SymmFactor

     data_headers(5) = "Ham Res2P1.1D"
     export_data(5) = Ham_const_excised_VolIntD*SymmFactor

     if (enable_M_constraint==1) then
        data_headers(6) = "M_constraint"
        export_data(6) = M_constraint*SymmFactor
        data_headers(7) = "M_constraint_inner"
        export_data(7) = M_constraint_inner*SymmFactor
        num_cols = 7
     end if

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     !************ MOMENTUM CONSTRAINT ************!

     filename = 'ks-mom.con'
     num_cols = 9
     export_data = 0.D0

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = "Momx ResN"
     export_data(2) = momx_const_VolIntN*SymmFactor
     data_headers(3) = "Momy ResN"
     export_data(3) = momy_const_VolIntN*SymmFactor
     data_headers(4) = "Momz ResN"
     export_data(4) = momz_const_VolIntN*SymmFactor
     data_headers(5) = "Mom ResD"
     export_data(5) = mom_const_VolIntD*SymmFactor

     data_headers(6) = "Momx Res1.1N"
     export_data(6) = momx_const_excised_VolIntN*SymmFactor
     data_headers(7) = "Momy Res1.1N"
     export_data(7) = momy_const_excised_VolIntN*SymmFactor
     data_headers(8) = "Momz Res1.1N"
     export_data(8) = momz_const_excised_VolIntN*SymmFactor
     data_headers(9) = "Mom Res1.1D"
     export_data(9) = mom_const_excised_VolIntD*SymmFactor

     if (enable_Jz_constraint==1) then
        data_headers(num_cols+1) = "Jz_constraint"
        export_data(num_cols+1) = Jz_constraint*SymmFactor
        data_headers(num_cols+2) = "Jz_constraint_inner"
        export_data(num_cols+2) = Jz_constraint_inner*SymmFactor
        num_cols = num_cols+2
     end if

     if (enable_P_constraint==1) then
        data_headers(num_cols+1) = "P_constraint"
        export_data(num_cols+1) = P_constraint*SymmFactor
        data_headers(num_cols+2) = "P_constraint_inner"
        export_data(num_cols+2) = P_constraint_inner*SymmFactor
        num_cols = num_cols+2
     end if

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!
     !*********** GAMMA CONSTRAINT ************!
     filename = 'ks-gam.con'
     num_cols = 4

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     data_headers(2) = "Gamx Res"
     export_data(2) = Gamx_const_VolInt*SymmFactor
     data_headers(3) = "Gamy Res"
     export_data(3) = Gamy_const_VolInt*SymmFactor
     data_headers(4) = "Gamz Res"
     export_data(4) = Gamz_const_VolInt*SymmFactor

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

     !--------------------------------------------------------------------------------!

     do ii=1,nsurf

        surf_radius=surfvec(ii)
        rbr=rbrvec(ii)
        drbdr=drbrvec(ii)
        ddrbddr=ddrbrvec(ii)

        if(nsurf.eq.1) then
           filename = 'ks.jon'
        else
           write(filename,31)ii
31         FORMAT('ks.jon.',I1)
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

        data_headers(5) = "Ang_mom_sur"
        call J_surf_integral_offcenter(cctkGH,0.D0,0.D0,0.D0,export_data(5))

        data_headers(6) = "Komar Mass"
        call Komar_surf_integralv2(cctkGH,export_data(6))

        data_headers(7) = "P_x2"
        data_headers(8) = "P_y2"
        data_headers(9) = "P_z2"
        call P_surf_integralv2(cctkGH,export_data(7),export_data(8),export_data(9))

        if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end do

     ! Output a grid function along x-axis
     !!npts = 100
     !!xstart = 0.1414d0
     !!xend = 1.d0
     !!d_spacing = (xend-xstart)/(npts-1)
     !!allocate(pointcoords(npts,3))
     !!allocate(phiint(npts))

     !!do jj=1,npts
     !!   pointcoords(jj,1) = xstart + (jj-1)*d_spacing
     !!   pointcoords(jj,2) = xstart + (jj-1)*d_spacing
     !!   pointcoords(jj,3) = xstart + (jj-1)*d_spacing
     !!end do
     !!
     !!call CCTK_VarIndex(vindex,"bssn::phi")
     !!call interp_driver_carp(cctkGH,npts,pointcoords,vindex,phiint)

     !!open(18,file='phi_x.dat',position='append')
     !!write(18,*) '# t = ',cctk_time
     !!do jj=1,npts
     !!   write(18,*) pointcoords(jj,1)*sqrt(2.d0),phiint(jj)
     !!end do
     !!write(18,*) ' '
     !!close(18)

     !!deallocate(pointcoords,phiint)
  end if

end subroutine kerr_diagnostics
