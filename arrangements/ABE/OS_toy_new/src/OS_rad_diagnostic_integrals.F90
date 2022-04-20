!-----------------------------------------------------------------------------
!
!$Id: predict_mhd.F90  $
!
!-----------------------------------------------------------------------------
!
! Predict & Correct RHS routines for mhd variables
!
!-----------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
!
! Predictor
!
!-----------------------------------------------------------------------------
subroutine OS_rad_diag_integrals(CCTK_ARGUMENTS)

  implicit none
  interface
     subroutine OS_rad_analytic_new(Fcoeff,z,R_a,tau,R_OS,M,Eo4PiB,tau_a,tau_s,eta,Q,rho_b,Ec,Fc,E,F)
       use nr
       use nrtype
       implicit none

       real(DP), intent(in) :: R_a,tau,R_OS,M,Eo4PiB,tau_a,tau_s
       real(DP), intent(out) :: E,F,rho_b,z,Fcoeff
       real(DP) :: eta,tau_t,sinchi0,chi0,Q,Ec,Fc,x1,x2,xacc,term
       integer :: n
       interface
          FUNCTION func(eta,tau,R_OS,M)
            USE nrtype
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: eta,tau,R_OS,M
            REAL(DP) :: func
          END FUNCTION func
       end interface
     end subroutine OS_rad_analytic_new
  end interface

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext,ext_noghost,global_ext
  real*8, dimension(20)                     :: export_data
  character, dimension(20)                     :: data_headers*20
  character                                :: filename*50
  real*8                                   :: dT,dX,dY,dZ,xmax
  real*8                                   :: RestMass,ADMMass,Jtot,temp,resid,resid_norm
  real*8                                   :: resid_x,resid_y,resid_z,resid_mom_norm
  real*8                                   :: dIntegral
  real*8                                   :: gxres,gyres,gzres
  real*8                                   :: multfactor
  real*8                                   :: rho_tiny
  integer                                  :: interpolate_order

  !Needed for GW extraction:
  real*8                                   :: Mass_GW,E_GW,hplus,hcross,J_GW
  real*8                                   :: psi4r,psi4i
  real*8,dimension(2)                      :: gw_amplitudes
  !End GW extraction variables  

  integer                                  :: adjimin, adjjmin, adjkmin
  integer                                  :: adjimax, adjjmax, adjkmax
  CCTK_REAL :: reduction_value
  CCTK_INT :: red_tmp
  integer :: header_flag,handle,dummy
  integer :: index,num_cols
  integer :: ierr,myproc_rank
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  real*8                                     :: PI
  real*8                                     :: rho_b_initial,SymmFactor, E_initial
  real*8, allocatable,dimension(:)           :: gxx_pcle2, gxz_pcle2, gzz_pcle2, gyy_pcle2, gxy_pcle2, gyz_pcle2, phi_pcle2
  real*8, allocatable,dimension(:)           :: rho_b_pcle2, P_pcle2, E_rad_pcle2, F_rad_scalar_pcle2 
  integer :: ii,jj, vindex

  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  write(*,*) "Start OS_rad_diagnostic_integrals.F90!!!!"
  
  PI = acos(-1.d0)

  ext(1) = cctk_lsh(1)
  ext(2) = cctk_lsh(2)
  ext(3) = cctk_lsh(3)
  global_ext = cctk_gsh

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  if(mod(cctk_iteration,out_every)==0) then
     call set_SymmFactor(Symmetry,SymmFactor)
  end if

  write(*,*) "SymmFactor =", SymmFactor

  call CCTK_ReductionHandle(handle,"maximum")
  call CCTK_VarIndex(index,"grid::X")
  call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xmax,1,index)


  !header_flag == 1 -> write file headers
  if(CCTK_TIME == 0.D0) then
     header_flag = 1
  else
     header_flag = 0
  end if

  !--------------------------------------------------------------------------------!

 if (mod(CCTK_iteration,out_every) .eq. 0) then

  filename = 'OS_rad.mon'
  num_cols = 7
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

  data_headers(4) = 'gxxcent-1'
  call CCTK_VarIndex (vindex, "bssn::gxx")
  call interp_center_carp(cctkGH,vindex,export_data(4))
  export_data(4) = export_data(4) - 1

  data_headers(5) = 'rhocent'
  call CCTK_VarIndex (vindex, "mhd_evolve::rho_b")
  call interp_center_carp(cctkGH,vindex,export_data(5))
  write(*,*) "rho center = ",export_data(5)

  data_headers(6) = 'lapsecent'
  call CCTK_VarIndex (vindex, "lapse::lapm1")
  call interp_center_carp(cctkGH,vindex,export_data(6))
  export_data(6) = export_data(6) + 1.D0


  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
  end if

  !--------------------------------------------------------------------------------!


     rho_b_initial=M_OS/(4.d0*PI*R_OS**3/3.d0)

     write(*,*) "rho_b_initial = ", rho_b_initial

     do ii=1,3
        do jj=1,narr
           pos(jj,ii) = coord(jj,ii+1)
        end do
     end do

     allocate(gxx_pcle2(narr))
     allocate(gzz_pcle2(narr))
     allocate(gxz_pcle2(narr))
     allocate(gyy_pcle2(narr))
     allocate(gyz_pcle2(narr))
     allocate(gxy_pcle2(narr))
     allocate(phi_pcle2(narr))
     allocate(rho_b_pcle2(narr))
     allocate(P_pcle2(narr))
     allocate(E_rad_pcle2(narr))
     allocate(F_rad_scalar_pcle2(narr))


     call CCTK_VarIndex(vindex,"bssn::gxx")       
     call interp_driver_carp(cctkGH,narr,pos,vindex,gxx_pcle2)
     call CCTK_VarIndex(vindex,"bssn::gxz")
     call interp_driver_carp(cctkGH,narr,pos,vindex,gxz_pcle2)
     call CCTK_VarIndex(vindex,"bssn::gzz")
     call interp_driver_carp(cctkGH,narr,pos,vindex,gzz_pcle2)
     call CCTK_VarIndex(vindex,"bssn::gyy")
     call interp_driver_carp(cctkGH,narr,pos,vindex,gyy_pcle2)
     call CCTK_VarIndex(vindex,"bssn::gxy")
     call interp_driver_carp(cctkGH,narr,pos,vindex,gxy_pcle2)
     call CCTK_VarIndex(vindex,"bssn::gyz")
     call interp_driver_carp(cctkGH,narr,pos,vindex,gyz_pcle2)     
     call CCTK_VarIndex(vindex,"bssn::phi")
     call interp_driver_carp(cctkGH,narr,pos,vindex,phi_pcle2)
     call CCTK_VarIndex(vindex,"mhd_evolve::rho_b")
     call interp_driver_carp(cctkGH,narr,pos,vindex,rho_b_pcle2)
     call CCTK_VarIndex(vindex,"mhd_evolve::P")
     call interp_driver_carp(cctkGH,narr,pos,vindex,P_pcle2)
     call CCTK_VarIndex(vindex,"mhd_evolve::E_rad")
     call interp_driver_carp(cctkGH,narr,pos,vindex,E_rad_pcle2)
     call CCTK_VarIndex(vindex,"mhd_evolve::F_rad_scalar")
     call interp_driver_carp(cctkGH,narr,pos,vindex,F_rad_scalar_pcle2)


     areal_radius=exp(2.d0*phi_pcle2)*sqrt(coord(:,4)*coord(:,4)*gxx_pcle2-2.d0*coord(:,2)*coord(:,4)*gxz_pcle2+coord(:,2)*coord(:,2)*gzz_pcle2)
 
!     areal_radius=exp(2.d0*phi_pcle2)*sqrt(coord(:,3)*coord(:,3)*gxx_pcle2-2.d0*coord(:,2)*coord(:,3)*gxy_pcle2+coord(:,2)*coord(:,2)*gyy_pcle2) 

     E_initial = E_over_rho*rho_b_initial**(gamma_OS)
!     E_initial = E_over_rho*rho_b_initial
     write(*,*) "E_initial = ", E_initial
     if (rad_closure_scheme.eq.0) then
        do ii=1,narr
           call OS_rad_analytic_new(F_coeff(ii),z_anal(ii),areal_radius(ii),&
                coord(ii,1),R_OS,1.d0,E_initial,opt_depth_a,opt_depth_s,eta(ii),Q(ii),&
                rho_b_anal(ii),Ec_anal(ii),Fc_anal(ii),E_anal(ii),F_anal(ii))
           
           Ec_pcle(ii) = E_rad_pcle2(ii)*Q(ii)**4
           Fc_pcle(ii) = F_rad_scalar_pcle2(ii)*Q(ii)**4
        end do
     end if
     write(*,*) "End OS_rad_analytic_new "
  !--------------------------------------------------------------------------------!                                              
  !--------------------------------------------------------------------------------!
     filename = 'OS_rad_cent.bon'

     num_cols = 16

     data_headers(1) = "Coord time"
     export_data(1)=CCTK_TIME
     data_headers(2) = "tau"
     export_data(2)=coord(1,1)
     data_headers(3) = "tc/tcoll"
     export_data(3)=eta(1)/PI
     data_headers(4) = "R_a"
     export_data(4)=areal_radius(1)
     data_headers(5) = "rho_b"
     export_data(5)=rho_b_pcle2(1)
     data_headers(6) = "P"
     export_data(6)=P_pcle2(1)
     data_headers(7) = "E_rad"
     export_data(7)=E_rad_pcle2(1)
     data_headers(8) = "F_rad_scalar"
     export_data(8)=F_rad_scalar_pcle2(1)
     data_headers(9) = "Ec"
     export_data(9)=Ec_pcle(1)
     data_headers(10) = "Fc"
     export_data(10)=Fc_pcle(1)
     data_headers(11) = "rho_b_anal"
     export_data(11)=rho_b_anal(1)*rho_b_initial
     data_headers(12) = "P_anal"
     export_data(12)=rho_b_anal(1)*E_anal(1)**0.25d0*P_over_rho*rho_b_initial
     data_headers(13) = "E_rad_anal"
     export_data(13)=E_anal(1)*rho_b_initial*P_over_rho/Po4PiB
     data_headers(14) = "F_anal"
     export_data(14)=F_anal(1)*rho_b_initial*P_over_rho/Po4PiB
     data_headers(15) = "Ec_anal"
     export_data(15)=Ec_anal(1)*rho_b_initial*P_over_rho/Po4PiB
     data_headers(16) = "Fc_anal"
     export_data(16)=Fc_anal(1)*rho_b_initial*P_over_rho/Po4PiB

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if

     filename = 'OS_rad_mid.bon'

     num_cols = 16

     data_headers(1) = "Coord time"
     export_data(1)=CCTK_TIME
     data_headers(2) = "tau"
     export_data(2)=coord(narr/2,1)
     data_headers(3) = "tc/tcoll"
     export_data(3)=eta(narr/2)/PI
     data_headers(4) = "R_a"
     export_data(4)=areal_radius(narr/2)
     data_headers(5) = "rho_b"
     export_data(5)=rho_b_pcle2(narr/2)
     data_headers(6) = "P"
     export_data(6)=P_pcle2(narr/2)
     data_headers(7) = "E_rad"
     export_data(7)=E_rad_pcle2(narr/2)
     data_headers(8) = "F"
     export_data(8)=F_rad_scalar_pcle2(narr/2)
     data_headers(9) = "Ec"
     export_data(9)=Ec_pcle(narr/2)
     data_headers(10) = "Fc"
     export_data(10)=Fc_pcle(narr/2)
     data_headers(11) = "rho_b_anal"
     export_data(11)=rho_b_anal(narr/2)*rho_b_initial
     data_headers(12) = "P_anal"
     export_data(12)=rho_b_anal(narr/2)*E_anal(narr/2)**0.25d0*P_over_rho*rho_b_initial
     data_headers(13) = "E_rad_anal"
     export_data(13)=E_anal(narr/2)*rho_b_initial*P_over_rho/Po4PiB
     data_headers(14) = "F_anal"
     export_data(14)=F_anal(narr/2)*rho_b_initial*P_over_rho/Po4PiB
     data_headers(15) = "Ec_anal"
     export_data(15)=Ec_anal(narr/2)*rho_b_initial*P_over_rho/Po4PiB
     data_headers(16) = "Fc_anal"
     export_data(16)=Fc_anal(narr/2)*rho_b_initial*P_over_rho/Po4PiB

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if




     filename = 'OS_rad_mid.bon'
 
     num_cols = 16

     data_headers(1) = "Coord time"                                             
     export_data(1)=CCTK_TIME
     data_headers(2) = "tau"                        
     export_data(2)=coord(narr/2,1)
     data_headers(3) = "tc/tcoll"        
     export_data(3)=eta(narr/2)/PI
     data_headers(4) = "R_a"
     export_data(4)=areal_radius(narr/2)
     data_headers(5) = "rho_b"   
     export_data(5)=rho_b_pcle2(narr/2)
     data_headers(6) = "P"   
     export_data(6)=P_pcle2(narr/2)
     data_headers(7) = "E_rad"   
     export_data(7)=E_rad_pcle2(narr/2)
     data_headers(8) = "F"   
     export_data(8)=F_rad_scalar_pcle2(narr/2)
     data_headers(9) = "Ec"   
     export_data(9)=Ec_pcle(narr/2)
     data_headers(10) = "Fc"   
     export_data(10)=Fc_pcle(narr/2)
     data_headers(11) = "rho_b_anal"   
     export_data(11)=rho_b_anal(narr/2)*rho_b_initial
     data_headers(12) = "P_anal"   
     export_data(12)=rho_b_anal(narr/2)*E_anal(narr/2)**0.25d0*P_over_rho*rho_b_initial
     data_headers(13) = "E_rad_anal"   
     export_data(13)=E_anal(narr/2)*rho_b_initial*P_over_rho/Po4PiB
     data_headers(14) = "F_anal"   
     export_data(14)=F_anal(narr/2)*rho_b_initial*P_over_rho/Po4PiB
     data_headers(15) = "Ec_anal"   
     export_data(15)=Ec_anal(narr/2)*rho_b_initial*P_over_rho/Po4PiB
     data_headers(16) = "Fc_anal"   
     export_data(16)=Fc_anal(narr/2)*rho_b_initial*P_over_rho/Po4PiB

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if
     !--------------------------------------------------------------------------------!
     !--------------------------------------------------------------------------------!

     filename = 'OS_rad_edge.bon'

     num_cols = 16

     data_headers(1) = "Coord time"                                             
     export_data(1)=CCTK_TIME
     data_headers(2) = "tau"                        
     export_data(2)=coord(narr,1)
     data_headers(3) = "tc/tcoll"        
     export_data(3)=eta(narr)/PI
     data_headers(4) = "R_a"
     export_data(4)=areal_radius(narr)
     data_headers(5) = "rho_b"   
     export_data(5)=rho_b_pcle2(narr)
     data_headers(6) = "P"   
     export_data(6)=P_pcle2(narr)
     data_headers(7) = "E_rad"   
     export_data(7)=E_rad_pcle2(narr)
     data_headers(8) = "F"   
     export_data(8)=F_rad_scalar_pcle2(narr)
     data_headers(9) = "Ec"   
     export_data(9)=Ec_pcle(narr)
     data_headers(10) = "Fc"   
     export_data(10)=Fc_pcle(narr)
     data_headers(11) = "rho_b_anal"   
     export_data(11)=rho_b_anal(narr)*rho_b_initial
     data_headers(12) = "P_anal"   
     export_data(12)=rho_b_anal(narr)*E_anal(narr)**0.25d0*P_over_rho*rho_b_initial
     data_headers(13) = "E_rad_anal"   
     export_data(13)=E_anal(narr)*rho_b_initial*P_over_rho/Po4PiB
     data_headers(14) = "F_anal"   
     export_data(14)=F_anal(narr)*rho_b_initial*P_over_rho/Po4PiB
     data_headers(15) = "Ec_anal"   
     export_data(15)=Ec_anal(narr)*rho_b_initial*P_over_rho/Po4PiB
     data_headers(16) = "Fc_anal"   
     export_data(16)=Fc_anal(narr)*rho_b_initial*P_over_rho/Po4PiB

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if


     filename = 'OS_rad.don'

     num_cols = 5

     data_headers(1) = '# Time'
     export_data(1) = CCTK_TIME

     ! M0 volume integral is computed in scheduler of diagnostics_mhd                                                               
     data_headers(2) = "M0"
     export_data(2) = M0_VolInt*SymmFactor

     ! Eint (minternal_hybrid) volume integral is computed in scheduler of diagnostics_mhd                                           
     data_headers(3) = 'Eint'
     export_data(3) = Minternal_VolInt*SymmFactor

     ! Kinetic energy (T) volume integral is computed in scheduler of diagnostics_mhd 
     data_headers(4) = 'T'
     export_data(4) = T_VolInt*SymmFactor
 
    ! M_ADM volume integral is computed in scheduler of diagnostics_vacuum                                                   
     data_headers(5) = "M_ADM_Vol"
     export_data(5) = M_ADM_VolInt*SymmFactor

     if(CCTK_MyProc(CCTKGH) .eq. 0) then
        call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)
     end if

     !--------------------------------------------------------------------------------!                                                                                                                      
     !************ HAMILTONIAN CONSTRAINT ************!                                                                                                                                                      

     filename = "OS-toy-ham.con"
     num_cols = 11
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

     data_headers(6) = "HR1.1-NoD"
     export_data(6) = Ham_const_excised_VolIntN/Ham_const_excised_VolIntD

     data_headers(7) = "HRes_outAH_N"
     !export_data(7) = Ham_const_outside_AH_mask_VolIntN*SymmFactor                                                                                                                                          

     data_headers(8) = "HRes_outAH_D"
     !export_data(8) = Ham_const_outside_AH_mask_VolIntD*SymmFactor                                                                                                                                          

     data_headers(9) = "HRes_outAHNoD"
     !export_data(9) = Ham_const_outside_AH_mask_VolIntN/Ham_const_outside_AH_mask_VolIntD   

     data_headers(10) = "HRes2P1.1NNR"
     export_data(10) = Ham_const_excised_innerregion_VolIntN*SymmFactor

     data_headers(11) = "HRes2P1.1DNR"
     export_data(11) = Ham_const_excised_innerregion_VolIntD*SymmFactor


     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)



     !--------------------------------------------------------------------------------!                                                                                                                      
     !************ MOMENTUM CONSTRAINT ************!                                                                                                                                                         
     filename = 'OS-toy-mom.con'
     num_cols = 10
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

     data_headers(10) = "MRx 1.1NoD"
     export_data(10) = momx_const_excised_VolIntN/mom_const_excised_VolIntD

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data,header_flag)

end if


  !--------------------------------------------------------------------------------!
  !particle tracer stuff for analytic solutions

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     filename = 'analytic_profiles.xg'
     if (CCTK_TIME .eq. 0) then
        open(UNIT=15,FILE=filename,STATUS="REPLACE")
        write(15,*) ""
        write(15,*) '"analytic profiles'
     else
        open(UNIT=15,FILE=filename,STATUS="OLD",POSITION="APPEND")
     endif
     write(*,*) "iteration: ",CCTK_iteration
     if (mod(CCTK_iteration,out_every) .eq. 0) then
        write(15,*) ""
        write(15,*) ""
        write(15,*) '"Time =',CCTK_TIME

        write(15,'(t2,a6,t17,a3,t31,a3,t45,a5,t59,a1,t73,a1,t87,a1,t101,a2,t115,a2,t129,a3,t143,a1)') '"coord','tau','R_a','rho_b','P','E','F','Ec','Fc','eta','z'
        do ii=1,narr
           write(15,'(11e18.10)') coord(ii,2),coord(ii,1),areal_radius(ii),rho_b_anal(ii)*rho_b_initial,&
                rho_b_anal(ii)*rho_b_initial*(E_anal(ii)/rad_const)**(0.25d0),&
                E_anal(ii),&
                F_anal(ii),&
                Ec_anal(ii),Fc_anal(ii),eta(ii),z_anal(ii)
        end do
     endif
     close(15)
  endif

  !--------------------------------------------------------------------------------!

  !numerical solutions at the positions of particles 

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     filename = 'numerical_profiles.xg'
     if (CCTK_TIME .eq. 0) then
        open(UNIT=17,FILE=filename,STATUS="REPLACE")
        write(17,*) ""
        write(17,*) '"numerical profiles'
     else
        open(UNIT=17,FILE=filename,STATUS="OLD",POSITION="APPEND")
     endif
     if (mod(CCTK_iteration,out_every) .eq. 0) then
        write(17,*) ""
        write(17,*) ""
        write(17,*) '"Time =',CCTK_TIME

        write(17,'(t2,a6,t17,a3,t31,a3,t45,a5,t59,a1,t73,a1,t87,a1,t101,a2,t115,a2,t129,a3,t143,a1)') '"coord','tau','R_a','rho_b','P','E','F','Ec','Fc','eta','z'
        do ii=1,narr
           write(17,'(11e18.10)') coord(ii,2),coord(ii,1),areal_radius(ii),rho_b_pcle2(ii),P_pcle2(ii),E_rad_pcle2(ii),F_rad_scalar_pcle2(ii),Ec_pcle(ii),Fc_pcle(ii),eta(ii),z_anal(ii)
        end do
     endif
     close(17)
  endif

  !----------------------------------------------------------------------------------!

  !coord files

  if(CCTK_MyProc(CCTKGH) .eq. 0) then
     filename = 'coord_profiles.xg'
     if (CCTK_TIME .eq. 0) then
        open(UNIT=19,FILE=filename,STATUS="REPLACE")
        write(19,*) ""
        write(19,*) '"coordinate profiles'
     else
        open(UNIT=19,FILE=filename,STATUS="OLD",POSITION="APPEND")
     endif
     if (mod(CCTK_iteration,out_every) .eq. 0) then
        write(19,*) ""
        write(19,*) ""
        write(19,*) '"Time =',CCTK_TIME

        write(19,'(t2,a11,t17,a11,t31,a11,t45,a11,t59,a12,t73,a12,t87,a12,t101,a12,t115)') '"coord(ii,2)','coord(ii,1)','slope(ii,2)','slope(ii,1)','slopet(ii,2)','slopet(ii,1)','slopem(ii,2)','slopem(ii,1)'
        do ii=1,narr
           write(19,'(11e18.10)') coord(ii,2),coord(ii,1),coord(ii,1),slope(ii,2),slope(ii,1),slopet(ii,2),slopet(ii,1),slopem(ii,2),slopem(ii,1)
        end do
     endif
     close(19)
  endif


write(*,*) "END OF OS_rad_diag_integrals"

end subroutine OS_rad_diag_integrals
