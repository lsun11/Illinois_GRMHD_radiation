#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine BHNS_setup_emfield_part2_local_setBi_after_checkpoint(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext
  real*8                                   :: dX,dY,dZ, xc,yc,zc,rad
  real*8                                   :: x_NS_CoM_coord,y_NS_CoM_coord,fs4pi
  real*8                                   :: x_NS_CoM_coord2,y_NS_CoM_coord2
  real*8                                   :: vx_CM,vy_CM
  real*8,allocatable,dimension(:,:,:)      :: A_phi,A_phix,A_phiy,A_phiz
  real*8,dimension(1,3)                    :: points_CM
  integer                                  :: AXISYM,EQUATORIAL, i,j,k
  integer                                  :: OCTANT
  integer                                  :: vindex
  real*8, parameter                        :: SYM = 1.d0, ANTI = -1.d0
  CCTK_REAL reduction_value
  parameter(EQUATORIAL = 1, OCTANT = 2, AXISYM = 4)
  !

  write(*,*) "INSIDE BHNS_setup_emfield_part2_local_setBi_after_checkpoint!",bhnsinsertBNOW

  ! Note that if ITERATION_TO_INSERT_MAGNETIC_FIELDS is set to 0, we will instead just use the regular BHNS_setup_emfield_part2_local_setBi to set the magnetic fields.
  if(CCTK_ITERATION .eq. ITERATION_TO_INSERT_MAGNETIC_FIELDS .and. ITERATION_TO_INSERT_MAGNETIC_FIELDS.gt.0) then
     write(*,*) "INSERTING B FIELDS AT T>0!"

     fs4pi = sqrt(4.d0*acos(-1.d0))
     ext = cctk_lsh

     allocate(A_phi(ext(1),ext(2),ext(3)))
     allocate(A_phix(ext(1),ext(2),ext(3)))
     allocate(A_phiy(ext(1),ext(2),ext(3)))
     allocate(A_phiz(ext(1),ext(2),ext(3)))

     ! Compute magnetic vector potential A_{\phi} with Ab=1
     ! Note that A_{\phi} is calculated from the center of each NS

     !     x_NS_CoM_coord = CoMx_VolInt/CoM_VolInt_denominator
     !     y_NS_CoM_coord = CoMy_VolInt/CoM_VolInt_denominator

   !  x_NS_CoM_coord = CoMx_VolInt/CoM_VolInt_denominator !initial_ns_coord_x
   !  y_NS_CoM_coord = CoMy_VolInt/CoM_VolInt_denominator !initial_ns_coord_y


     if(use_new_bhns_initial_data.eq.3) then
        
        ! compute the position of the NSs
        x_NS_CoM_coord = Box1X_VolInt1/Box1denom_VolInt1
        y_NS_CoM_coord = Box1Y_VolInt1/Box1denom_VolInt1

        x_NS_CoM_coord2 = Box1X_VolInt2/Box1denom_VolInt2
        y_NS_CoM_coord2 = Box1Y_VolInt2/Box1denom_VolInt2
     else
        ! compute the position of the NS
        x_NS_CoM_coord = Box1X_VolInt/Box1denom_VolInt
        y_NS_CoM_coord = Box1Y_VolInt/Box1denom_VolInt
     end if

     
     if(center_Bfields_around_BH_instead==1) then
        x_NS_CoM_coord = CoMx_VolInt/CoM_VolInt_denominator !initial_ns_coord_x
        y_NS_CoM_coord = CoMy_VolInt/CoM_VolInt_denominator !initial_ns_coord_y
     end if
        

     write(*,*) "INSIDE EMFIELDS SETUP: x,y of ns:",x_NS_CoM_coord,y_NS_CoM_coord
     write(*,*) "INSIDE EMFIELDS SETUP: Pmax:",bhns_P_max,p_c

     if(unequalmass.eq.1) then
        write(*,*) "INSIDE EMFIELDS SETUP: Pmax2:",bhns_P_max2
     end if


     print *, " driver_setup_emfields_part2-local-setBi_after_checkpoint.F90 (x_NS,y_NS)=",x_NS_CoM_coord,y_NS_CoM_coord

     if (em_field_type==1 .and. Sym_Bz .gt. 0.d0) then 
        write(*,*) "Sorry, em_field_type==1 (toroidal fields) not yet supported with Sym_Bz==1.  You can set Sym_Bz=-1 !"
        stop
     end if

     if (enable_trace_field_line==1 .and. em_field_type==1) then 
        write(*,*) "Sorry, the initial data for the field line tracer variables currently only supports the poloidal initial data. You have to write your own initial data for the toroidal configuration."
     end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! FIRST, WE SET UP MAGNETIC FIELDS ON _P LEVEL:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! The magnetic fields are set based on pressure. We don't have P_p or P_p_p, so we must compute that using the primitives solver.
     ! The primitives solver depends on gupij, so first compute gupij, based on gij_p:
     temp1 = 1.D0/(gxx_p * gyy_p * gzz_p + gxy_p * gyz_p * gxz_p + gxz_p * gxy_p * gyz_p &
          - gxz_p * gyy_p * gxz_p - gxy_p * gxy_p * gzz_p - gxx_p * gyz_p * gyz_p)

      gupxx =   ( gyy_p * gzz_p - gyz_p * gyz_p )* temp1
      gupxy = - ( gxy_p * gzz_p - gyz_p * gxz_p )* temp1
      gupxz =   ( gxy_p * gyz_p - gyy_p * gxz_p )* temp1
      gupyy =   ( gxx_p * gzz_p - gxz_p * gxz_p )* temp1
      gupyz = - ( gxx_p * gyz_p - gxy_p * gxz_p )* temp1
      gupzz =   ( gxx_p * gyy_p - gxy_p * gxy_p )* temp1

     write(*,*) "HELLO BEFORE _P primitives!"

     ! Next compute primitives, based on _p conservatives:
     if(primitives_solver==11) then
        !SKIP INITIAL GUESS.  THIS IS NOT NEEDED ANYWAY.

        !Here we use the HARM primitives solver, with a new prescription that minimizes changes to
        !  conservatives without applying the "Font" fix.
        !We intend to make the below function truly generic, but currently it only supports Gamma=2.
        !  It can be trivially extended for arbitrary Gamma-law EOS's, but will require work for more
        !  generic EOS's.
        !We hope the below function will eventually be used in place of other primitives solver options,
        !  since it can be easily extended to use your desired technique for solving primitives.
        call primitives_generic(cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z, &
             phi_p,gxx_p,gxy_p,gxz_p,gyy_p,gyz_p,gzz_p, &
             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
             lapm1_p,shiftx_p,shifty_p,shiftz_p, &
             Bx,By,Bz, &
             mhd_st_x_p,mhd_st_y_p,mhd_st_z_p,tau_p,rho_star_p, &
             vx,vy,vz,P,rho_b,h,u0, &
             rho_b_atm,tau_atm, &
             neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
             temp1,temp2,temp3,temp4,temp5, &
             primitives_debug,Psi6threshold)

        !SKIP COMPUTATION OF METRIC SOURCE TERMS.  THAT WOULD BE UNNECESSARY!
     else
        write(*,*) "ERROR. Primitives solver != 11 NOT ALLOWED when setting B fields from checkpoint."
        stop
     end if

     !
     ! Compute the vector potential A_phi
     !
     if (em_field_type .ne. 1) then 
        if(Symmetry==1 .and. angle_to_tilt_magnetic_fields.ne.0) then
           write(*,*) "CANNOT TILT THE MAGNETIC FIELDS WITH Symmetry==1."
           stop
        end if
        call BHNS_compute_Aphi(ext,X,Y,Z,PhysicalRadius,P,A_phi,Ax_p,Ay_p,Az_p,  & 
             mhd_psi_line, mhd_u_psi, mhd_chi_line, mhd_u_chi,                   &
             betam1,p_c,x_NS_CoM_coord,y_NS_CoM_coord,bhns_P_max,Sym_Bz,         & 
             enable_trace_field_line,constrained_transport_scheme,em_field_type, &
             angle_to_tilt_magnetic_fields,Aphi_power,rloop,Iloop,               &
             use_new_bhns_initial_data,x_NS_CoM_coord2,y_NS_CoM_coord2,          &
             bhns_P_max2,unequalmass,betam1_2,rloop2,Iloop2)
     end if

     !
     ! Compute the vector potential Ax and Ay for near toroidal configuration
     !
     if (em_field_type==1) then 
        if(constrained_transport_scheme==3) then
           call BHNS_compute_A_toroidal_constrained_transport_eq_3(ext,X,Y,Z,P,Ax_p,Ay_p,Az_p,phi_p,betam1, &
                p_c,x_NS_CoM_coord,y_NS_CoM_coord,r0,bhns_P_max)
        else if(constrained_transport_scheme==1) then
           call BHNS_compute_A_toroidal(ext,X,Y,Z,P,Ax_p,Ay_p,Az_p,phi_p,betam1, &
                p_c,x_NS_CoM_coord,y_NS_CoM_coord,r0,bhns_P_max)
        else
           write(*,*) "Sorry, em_field_type==1 (toroidal fields) only supported for constrained_transport_scheme == 1 or 3."
           stop
        end if
     end if

     ! Compute B^i's based on values of A_i's.  
     !   Note that the conservatives up to this point have assumed B=0.  Here we
     !   update the conservatives to include the now nonzero magnetic terms.
     call BHNS_compute_Bi(ext,X,Y,Z,PhysicalRadius, &
          phi_p, &
          A_phi,A_phix,A_phiy,A_phiz, & 
          Ax_p,Ay_p,Az_p, & 
          Bx,By,Bz, &
          x_NS_CoM_coord,y_NS_CoM_coord,rho_b_atm, &
          constrained_transport_scheme,em_field_type)


     ! Reset the atmospheric density to be variable and not constant,
     ! falling off according to some prescription specified by parameter atm_type
     if (RESET_RHO_B_ATM.eq.1) then 

!        if (atm_type.eq.1) then
        call bhns_compute_b2_cpp(cctkGH,cctk_lsh, phi_p, lapm1_p, &
              shiftx_p,shifty_p,shiftz_p,vx,vy,vz,Bx,By,Bz, & 
              gxx_p, gxy_p, gxz_p, gyy_p, gyz_p, gzz_p, temp2)
!        end if

        call reset_rhob_atm_f90(cctk_lsh,X,Y,Z,temp2,rho_b,P,vx,vy,vz,shiftx_p,shifty_p,bhns_B_v,atm_type, &
             x_NS_CoM_coord,y_NS_CoM_coord,bhns_vx_CM,bhns_vy_CM,bhns_vz_CM,bhns_rhob_fac,rad_rhob_fac,    &
             bhns_rhob_max,rhobatm_falloff_power,rhob_o_b2,bhns_R_NS,rho_b_atm,bhns_fac_atm,beta_ext,k_tab,&
             gamma_tab,neos,x_NS_CoM_coord2,y_NS_CoM_coord2,use_new_bhns_initial_data)
        
     end if

     ! Set the st_i's:
     st_x = mhd_st_x_p
     st_y = mhd_st_y_p
     st_z = mhd_st_z_p

     ! Recompute conservatives to be consistent with current magnetic fields and atmospheric density
     call recompute_conservatives_fast_standalone_gf(cctkGH,cctk_lsh, &
          rho_b,P,vx,vy,vz, &
          phi_p,lapm1_p, &
          shiftx_p,shifty_p,shiftz_p, &
          gxx_p,gxy_p,gxz_p,gyy_p,gyz_p,gzz_p, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
          Bx,By,Bz, &
          rho_star_p,tau_p,mhd_st_x_p,mhd_st_y_p,mhd_st_z_p)

     
     ! Reset magnetic fields to zero, so that in the below primitives solver call, the conservatives & primitives are consistent with B=0.
     !   Remember that the magnetic fields are assumed zero until they are computed inside BHNS_compute_Bi() on each
     !   timelevel. 
     Bx = 0.D0
     By = 0.D0
     Bz = 0.D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! NEXT, _P_P LEVEL:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! The magnetic fields are set based on pressure. We don't have P_p or P_p_p, so we must compute that using the primitives solver.
     ! The primitives solver depends on gupij, so first compute gupij, based on gij_p:
     temp1 = 1.D0/(gxx_p_p * gyy_p_p * gzz_p_p + gxy_p_p * gyz_p_p * gxz_p_p + gxz_p_p * gxy_p_p * gyz_p_p &
          - gxz_p_p * gyy_p_p * gxz_p_p - gxy_p_p * gxy_p_p * gzz_p_p - gxx_p_p * gyz_p_p * gyz_p_p)

      gupxx =   ( gyy_p_p * gzz_p_p - gyz_p_p * gyz_p_p )* temp1
      gupxy = - ( gxy_p_p * gzz_p_p - gyz_p_p * gxz_p_p )* temp1
      gupxz =   ( gxy_p_p * gyz_p_p - gyy_p_p * gxz_p_p )* temp1
      gupyy =   ( gxx_p_p * gzz_p_p - gxz_p_p * gxz_p_p )* temp1
      gupyz = - ( gxx_p_p * gyz_p_p - gxy_p_p * gxz_p_p )* temp1
      gupzz =   ( gxx_p_p * gyy_p_p - gxy_p_p * gxy_p_p )* temp1

     write(*,*) "HELLO BEFORE _P_P primitives!"

     ! Next compute primitives, based on _p_p conservatives:
     if(primitives_solver==11) then
        !SKIP INITIAL GUESS.  THIS IS NOT NEEDED ANYWAY.

        !Here we use the HARM primitives solver, with a new prescription that minimizes changes to
        !  conservatives without applying the "Font" fix.
        !We intend to make the below function truly generic, but currently it only supports Gamma=2.
        !  It can be trivially extended for arbitrary Gamma-law EOS's, but will require work for more
        !  generic EOS's.
        !We hope the below function will eventually be used in place of other primitives solver options,
        !  since it can be easily extended to use your desired technique for solving primitives.
        call primitives_generic(cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z, &
             phi_p_p,gxx_p_p,gxy_p_p,gxz_p_p,gyy_p_p,gyz_p_p,gzz_p_p, &
             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
             lapm1_p_p,shiftx_p_p,shifty_p_p,shiftz_p_p, &
             Bx,By,Bz, &
             mhd_st_x_p_p,mhd_st_y_p_p,mhd_st_z_p_p,tau_p_p,rho_star_p_p, &
             vx,vy,vz,P,rho_b,h,u0, &
             rho_b_atm,tau_atm, &
             neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
             temp1,temp2,temp3,temp4,temp5, &
             primitives_debug,Psi6threshold)

        !SKIP COMPUTATION OF METRIC SOURCE TERMS.  THAT WOULD BE UNNECESSARY!
     else
        write(*,*) "ERROR. Primitives solver != 11 NOT ALLOWED when setting B fields from checkpoint."
        stop
     end if

     ! Set the st_i's:
     st_x = mhd_st_x_p_p
     st_y = mhd_st_y_p_p
     st_z = mhd_st_z_p_p

     !
     ! Compute the vector potential A_phi
     !
     if (em_field_type .ne. 1) then 
        call BHNS_compute_Aphi(ext,X,Y,Z,PhysicalRadius,P,A_phi,Ax_p_p,Ay_p_p,Az_p_p, & 
             mhd_psi_line, mhd_u_psi, mhd_chi_line, mhd_u_chi,                        &
             betam1,p_c,x_NS_CoM_coord,y_NS_CoM_coord,bhns_P_max,Sym_Bz,              & 
             enable_trace_field_line,constrained_transport_scheme,em_field_type,      &
             angle_to_tilt_magnetic_fields,Aphi_power,rloop,Iloop,                    &
             use_new_bhns_initial_data,x_NS_CoM_coord2,y_NS_CoM_coord2,               &
             bhns_P_max2,unequalmass,betam1_2,rloop2,Iloop2)
     end if

     !
     ! Compute the vector potential Ax and Ay for near toroidal configuration
     !
     if (em_field_type==1) then 
        if(constrained_transport_scheme==3) then
           call BHNS_compute_A_toroidal_constrained_transport_eq_3(ext,X,Y,Z,P,Ax_p_p,Ay_p_p,Az_p_p,phi_p_p,betam1, &
                p_c,x_NS_CoM_coord,y_NS_CoM_coord,r0,bhns_P_max)
        else if(constrained_transport_scheme==1) then
           call BHNS_compute_A_toroidal(ext,X,Y,Z,P,Ax_p_p,Ay_p_p,Az_p_p,phi_p_p,betam1, &
                p_c,x_NS_CoM_coord,y_NS_CoM_coord,r0,bhns_P_max)
        else
           write(*,*) "Sorry, em_field_type==1 (toroidal fields) only supported for constrained_transport_scheme == 1 or 3."
           stop
        end if
     end if

     ! Compute B^i's based on values of A_i's.  
     !   Note that the conservatives up to this point have assumed B=0.  Here we
     !   update the conservatives to include the now nonzero magnetic terms.
     call BHNS_compute_Bi(ext,X,Y,Z,PhysicalRadius, &
          phi_p_p, &
          A_phi,A_phix,A_phiy,A_phiz, & 
          Ax_p_p,Ay_p_p,Az_p_p, & 
          Bx,By,Bz, &
          x_NS_CoM_coord,y_NS_CoM_coord,rho_b_atm, &
          constrained_transport_scheme,em_field_type)

     ! Reset the atmospheric density to be variable and not constant, falling off according to some prescription
     if (RESET_RHO_B_ATM.eq.1) then 
!        if (atm_type.eq.1) then
        call bhns_compute_b2_cpp(cctkGH,cctk_lsh, phi_p_p, lapm1_p_p, &
              shiftx_p_p,shifty_p_p,shiftz_p_p,vx,vy,vz,Bx,By,Bz, & 
              gxx_p_p, gxy_p_p, gxz_p_p, gyy_p_p, gyz_p_p, gzz_p_p, temp2)
!        end if

        call reset_rhob_atm_f90(cctk_lsh,X,Y,Z,temp2,rho_b,P,vx,vy,vz,shiftx,shifty,bhns_B_v,atm_type,      &
              x_NS_CoM_coord,y_NS_CoM_coord,bhns_vx_CM,bhns_vy_CM,bhns_vz_CM,bhns_rhob_fac,rad_rhob_fac,    &
              bhns_rhob_max,rhobatm_falloff_power,rhob_o_b2,bhns_R_NS,rho_b_atm,bhns_fac_atm,beta_ext,k_tab,&
              gamma_tab,neos,x_NS_CoM_coord2,y_NS_CoM_coord2,use_new_bhns_initial_data)

     end if

     ! Recompute conservatives to be consistent with current magnetic fields
     call recompute_conservatives_fast_standalone_gf(cctkGH,cctk_lsh, &
          rho_b,P,vx,vy,vz, &
          phi_p_p,lapm1_p_p, &
          shiftx_p_p,shifty_p_p,shiftz_p_p, &
          gxx_p_p,gxy_p_p,gxz_p_p,gyy_p_p,gyz_p_p,gzz_p_p, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
          Bx,By,Bz, &
          rho_star_p_p,tau_p_p,mhd_st_x_p_p,mhd_st_y_p_p,mhd_st_z_p_p)

     ! Reset magnetic fields to zero, so that in the below primitives solver call, the conservatives & primitives are consistent with B=0.
     !   Remember that the magnetic fields are assumed zero until they are computed inside BHNS_compute_Bi() on each
     !   timelevel.
     Bx = 0.D0
     By = 0.D0
     Bz = 0.D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! FINALLY, LATEST LEVEL (neither _p nor _p_p):
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! The magnetic fields are set based on pressure. We don't have P_p or P_p_p, so we must compute that using the primitives solver.
     ! The primitives solver depends on gupij, so first compute gupij, based on gij:
     temp1 = 1.D0/(gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz &
          - gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz)

      gupxx =   ( gyy * gzz - gyz * gyz )* temp1
      gupxy = - ( gxy * gzz - gyz * gxz )* temp1
      gupxz =   ( gxy * gyz - gyy * gxz )* temp1
      gupyy =   ( gxx * gzz - gxz * gxz )* temp1
      gupyz = - ( gxx * gyz - gxy * gxz )* temp1
      gupzz =   ( gxx * gyy - gxy * gxy )* temp1

     write(*,*) "HELLO BEFORE  primitives!"

     ! Next compute primitives, based on  conservatives:
     if(primitives_solver==11) then
        !SKIP INITIAL GUESS.  THIS IS NOT NEEDED ANYWAY.

        !Here we use the HARM primitives solver, with a new prescription that minimizes changes to
        !  conservatives without applying the "Font" fix.
        !We intend to make the below function truly generic, but currently it only supports Gamma=2.
        !  It can be trivially extended for arbitrary Gamma-law EOS's, but will require work for more
        !  generic EOS's.
        !We hope the below function will eventually be used in place of other primitives solver options,
        !  since it can be easily extended to use your desired technique for solving primitives.
        call primitives_generic(cctkGH,cctk_lsh,cctk_nghostzones,X,Y,Z, &
             phi,gxx,gxy,gxz,gyy,gyz,gzz, &
             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
             lapm1,shiftx,shifty,shiftz, &
             Bx,By,Bz, &
             mhd_st_x,mhd_st_y,mhd_st_z,tau,rho_star, &
             vx,vy,vz,P,rho_b,h,u0, &
             rho_b_atm,tau_atm, &
             neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &
             temp1,temp2,temp3,temp4,temp5, &
             primitives_debug,Psi6threshold)

        !SKIP COMPUTATION OF METRIC SOURCE TERMS.  THAT WOULD BE UNNECESSARY!
     else
        write(*,*) "ERROR. Primitives solver != 11 NOT ALLOWED when setting B fields from checkpoint."
        stop
     end if

     ! Set the st_i's:
     st_x = mhd_st_x
     st_y = mhd_st_y
     st_z = mhd_st_z

     !
     ! Compute the vector potential A_phi
     !
     if (em_field_type .ne. 1) then 
        call BHNS_compute_Aphi(ext,X,Y,Z,PhysicalRadius,P,A_phi,Ax,Ay,Az,        & 
             mhd_psi_line, mhd_u_psi, mhd_chi_line, mhd_u_chi,                   &
             betam1,p_c,x_NS_CoM_coord,y_NS_CoM_coord,bhns_P_max,Sym_Bz,         & 
             enable_trace_field_line,constrained_transport_scheme,em_field_type, &
             angle_to_tilt_magnetic_fields,Aphi_power,rloop,Iloop,               &
             use_new_bhns_initial_data,x_NS_CoM_coord2,y_NS_CoM_coord2,          &
             bhns_P_max2,unequalmass,betam1_2,rloop2,Iloop2)
     end if

     !
     ! Compute the vector potential Ax and Ay for near toroidal configuration
     !
     if (em_field_type==1) then 
        if(constrained_transport_scheme==3) then
           call BHNS_compute_A_toroidal_constrained_transport_eq_3(ext,X,Y,Z,P,Ax,Ay,Az,phi,betam1, &
                p_c,x_NS_CoM_coord,y_NS_CoM_coord,r0,bhns_P_max)
        else if(constrained_transport_scheme==1) then
           call BHNS_compute_A_toroidal(ext,X,Y,Z,P,Ax,Ay,Az,phi,betam1, &
                p_c,x_NS_CoM_coord,y_NS_CoM_coord,r0,bhns_P_max)
        else
           write(*,*) "Sorry, em_field_type==1 (toroidal fields) only supported for constrained_transport_scheme == 1 or 3."
           stop
        end if
     end if

     ! Compute B^i's based on values of A_i's.  
     !   Note that the conservatives up to this point have assumed B=0.  Here we
     !   update the conservatives to include the now nonzero magnetic terms.
     call BHNS_compute_Bi(ext,X,Y,Z,PhysicalRadius, &
          phi, &
          A_phi,A_phix,A_phiy,A_phiz, & 
          Ax,Ay,Az, & 
          Bx,By,Bz, & 
          x_NS_CoM_coord,y_NS_CoM_coord,rho_b_atm, &
          constrained_transport_scheme,em_field_type)

     ! Reset the atmospheric density to be variable and not constant, falling off according to some prescription
     if (RESET_RHO_B_ATM.eq.1) then 
!        if (atm_type.eq.1) then
        call bhns_compute_b2_cpp(cctkGH,cctk_lsh, phi, lapm1, &
              shiftx,shifty,shiftz,vx,vy,vz,Bx,By,Bz, & 
              gxx, gxy, gxz, gyy, gyz, gzz, temp2)
!        end if

        call reset_rhob_atm_f90(cctk_lsh,X,Y,Z,temp2,rho_b,P,vx,vy,vz,shiftx,shifty,bhns_B_v,atm_type,      &
              x_NS_CoM_coord,y_NS_CoM_coord,bhns_vx_CM,bhns_vy_CM,bhns_vz_CM,bhns_rhob_fac,rad_rhob_fac,    &
              bhns_rhob_max,rhobatm_falloff_power,rhob_o_b2,bhns_R_NS,rho_b_atm,bhns_fac_atm,beta_ext,k_tab,&
              gamma_tab,neos,x_NS_CoM_coord2,y_NS_CoM_coord2,use_new_bhns_initial_data)

     end if

     ! Recompute conservatives to be consistent with current magnetic fields
     call recompute_conservatives_fast_standalone_gf(cctkGH,cctk_lsh, &
          rho_b,P,vx,vy,vz, &
          phi,lapm1, &
          shiftx,shifty,shiftz, &
          gxx,gxy,gxz,gyy,gyz,gzz, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
          Bx,By,Bz, &
          rho_star,tau,mhd_st_x,mhd_st_y,mhd_st_z)

     ! Set the gauge variable to zero.
     psi6phi = 0.D0
     psi6phi_p = 0.D0
     psi6phi_p_p = 0.D0

     ! Compute B from A to set Stagger_Bs (on current timelevel) !
     call bhns_compute_B_from_A(CCTK_PASS_FTOF)

     ! recompute u0 velocity
     call reset_u0(cctk_lsh,phi,lapm1,shiftx,shifty,shiftz,gxx,gxy,gxz,gyy,gyz,gzz,vx,vy,vz,u0,X,Y,Z)
                  
     ! Free up memory:
     deallocate(A_phi, A_phix, A_phiy, A_phiz)
     
  end if
end subroutine BHNS_setup_emfield_part2_local_setBi_after_checkpoint



subroutine reset_rhob_atm_f90(cctk_lsh,X,Y,Z,b2,rho_b,P,vx,vy,vz,shiftx,shifty,B_v,atm_type, &
     x_NS_CoM_coord,y_NS_CoM_coord,vx_CM,vy_CM,vz_CM,bhns_rhob_fac,rad_rhob_fac,             &
     bhns_rhob_max,rhobatm_falloff_power,rhob_o_b2,R_NS,rho_b_atm,fac_atm,beta_ext,k_tab,    &
     gamma_tab,neos,x_NS_CoM_coord2,y_NS_CoM_coord2,initial_data)
  
  implicit none
  integer, dimension(3)                    :: cctk_lsh
  integer                                  :: neos
  real*8, dimension(neos+1)                :: gamma_tab, k_tab 
  real*8                                   :: xc,yc,zc,rad,bhns_rhob_fac,rhob_o_b2
  real*8                                   :: rad_rhob_fac,bhns_rhob_max,rhobatm_falloff_power
  real*8                                   :: x_NS_CoM_coord,y_NS_CoM_coord,beta_ext,rho_match
  real*8                                   :: vx_CM,vy_CM,vz_CM,R_NS,rho_b_atm,fac_atm, kloc, Pcold
  real*8                                   :: x_NS_CoM_coord2,y_NS_CoM_coord2,xc_b,yc_b,zc_b,rad_b
  integer                                  :: i,j,k, atm_type,B_v,initial_data   
  real*8,dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))      :: rho_b,b2,X,Y,Z,P
  real*8,dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))      :: vx,vy,vz
  real*8,dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))      :: shiftx,shifty
  
  if(beta_ext>0.0d0) then
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              xc = X(i,j,k)-x_NS_CoM_coord
              yc = Y(i,j,k)-y_NS_CoM_coord
              zc = Z(i,j,k)
              rad = sqrt(xc*xc+yc*yc)
              
              if(initial_data.eq.3) then
                 xc_b  =  X(i,j,k)-x_NS_CoM_coord2
                 yc_b  =  Y(i,j,k)-y_NS_CoM_coord2
                 zc_b  =  Z(i,j,k)
                 rad_b =  sqrt(xc_b*xc_b+yc_b*yc_b)
              end if
                            
              !Compute local cold pressure (may need to call the Cold EOS routine in the future)
              Pcold = k_tab(1)*rho_b(i,j,k)**gamma_tab(1)
              !Compute local K
              kloc=k_tab(1)*P(i,j,k)/Pcold
              if (kloc > 100.d0) then
                 kloc=100.d0
              end if
              
              if (kloc < 1.d0) then
                 kloc=1.d0
              end if

              ! Reset the atmospheric density
              rho_match  = sqrt(0.5d0*b2(i,j,k)*beta_ext/kloc)
              
              if (rho_b(i,j,k) <= rho_match) then
                 ! Now reset the atmosphere
                 if (atm_type==0) then
                    rho_b(i,j,k) = rho_match*(rad_rhob_fac/rad)**rhobatm_falloff_power 
                    P(i,j,k)     = kloc*rho_b(i,j,k)*rho_b(i,j,k)
                 else if (atm_type==1) then                
                    rho_b(i,j,k) = sqrt(beta_ext*0.5d0*b2(i,j,k)/kloc)
                    P(i,j,k)     = kloc*rho_b(i,j,k)*rho_b(i,j,k)
                 end if
                 ! Reset the velocity in the atmopshere      

                 if (initial_data.ne.3) then
                    if(rad<2.0D0*R_NS) then
                       vx(i,j,k)    = vx_CM
                       vy(i,j,k)    = vy_CM
                       vz(i,j,k)    = 0.0d0
                    else if(rad.ge.2.0D0*R_NS) then
                       vx(i,j,k)    = vx_CM*(2.0d0*R_NS/rad)**4 - shiftx(i,j,k)*(1.0d0 - (2.0d0*R_NS/rad)**4)
                       vy(i,j,k)    = vy_CM*(2.0d0*R_NS/rad)**4 - shifty(i,j,k)*(1.0d0 - (2.0d0*R_NS/rad)**4)
                       vz(i,j,k)    = 0.0d0
                    end if
                 end if
              end if
           end do
        end do
     end do
  end if
  print *,"(vx_CM,vy_CM,vz_CM,R_NS)=", vx_CM, vy_CM,vz_CM,R_NS,rho_b_atm
  
end subroutine reset_rhob_atm_f90
 
 
subroutine reset_u0(cctk_lsh,phi,lapm1,shiftx,shifty,shiftz,gxx,gxy,gxz,gyy,gyz,gzz,vx,vy,vz,u0,X,Y,Z)
   
   implicit none
   integer, dimension(3)                    :: cctk_lsh
   integer                                  :: i,j,k
   real*8,dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))      :: phi,lapm1
   real*8,dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))      :: shiftx,shifty,shiftz
   real*8,dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))      :: vx,vy,vz,u0,X,Y,Z
   real*8,dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))      :: gxx,gxy,gxz,gyy,gyz,gzz
   real*8                                    :: gxx_physL,gxy_physL,gxz_physL
   real*8                                    :: gyy_physL,gyz_physL,gzz_physL
   real*8                                    :: alpha,alpha_inv,alp_u02,u02        
   real*8                                    :: Psi4,v2,max,A


   do k=1,cctk_lsh(3)
      do j=1,cctk_lsh(2)
         do i=1,cctk_lsh(1)

            alpha     =  lapm1(i,j,k) + 1.0d0
            alpha_inv =  1.0/alpha 
            
            Psi4     =   exp(4.0d0*phi(i,j,k))
            
            gxx_physL  =  gxx(i,j,k)*Psi4
            gxy_physL  =  gxy(i,j,k)*Psi4
            gxz_physL  =  gxz(i,j,k)*Psi4
            gyy_physL  =  gyy(i,j,k)*Psi4
            gyz_physL  =  gyz(i,j,k)*Psi4
            gzz_physL  =  gzz(i,j,k)*Psi4
  

            v2 = (      gxx_physL*(vx(i,j,k) + shiftx(i,j,k))*(vx(i,j,k) + shiftx(i,j,k)) &
               +  2.0d0*gxy_physL*(vx(i,j,k) + shiftx(i,j,k))*(vy(i,j,k) + shifty(i,j,k)) &
               +  2.0d0*gxz_physL*(vx(i,j,k) + shiftx(i,j,k))*(vz(i,j,k) + shiftz(i,j,k)) &
               +        gyy_physL*(vy(i,j,k) + shifty(i,j,k))*(vy(i,j,k) + shifty(i,j,k)) &
               +  2.0d0*gyz_physL*(vy(i,j,k) + shifty(i,j,k))*(vz(i,j,k) + shiftz(i,j,k)) &
               +        gzz_physL*(vz(i,j,k) + shiftz(i,j,k))*(vz(i,j,k) + shiftz(i,j,k)))*alpha_inv*alpha_inv


            if(v2.ge.1.0d0) then
               max = 100.0d0
               A   =  (1.0d0 -  1.0d0/(max*max))/v2
               v2  =  A*v2
            end if

             alp_u02    =  1.0d0/(1.0d0 - v2)
             u02        =  alp_u02*alpha_inv*alpha_inv
             u0(i,j,k)  =  sqrt(u02)

          end do
       end do
    end do
    
 end subroutine reset_u0

