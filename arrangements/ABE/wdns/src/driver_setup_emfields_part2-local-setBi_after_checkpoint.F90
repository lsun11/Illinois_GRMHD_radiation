#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine WDNS_setup_emfield_part2_local_setBi_after_checkpoint(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext
  real*8                                   :: dX,dY,dZ
  real*8                                   :: x_NS_CoM_coord,y_NS_CoM_coord,fs4pi
  real*8,allocatable,dimension(:,:,:)      :: A_phi,A_phix,A_phiy,A_phiz
  integer                                  :: AXISYM,EQUATORIAL
  integer                                  :: OCTANT
  real*8, parameter                        :: SYM = 1.d0, ANTI = -1.d0
  CCTK_REAL reduction_value
  parameter(EQUATORIAL = 1, OCTANT = 2, AXISYM = 4)
  !

  if(CCTK_ITERATION .eq. ITERATION_TO_RESET_MAGNETIC_FIELDS) then
     write(*,*) "RESETTING B FIELDS AT CHECKPOINT!"

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

     x_NS_CoM_coord = CoMx_VolInt/CoM_VolInt_denominator !initial_ns_coord_x
     y_NS_CoM_coord = CoMy_VolInt/CoM_VolInt_denominator !initial_ns_coord_y

     write(*,*) "INSIDE EMFIELDS SETUP: x,y of ns:",x_NS_CoM_coord,y_NS_CoM_coord
     write(*,*) "INSIDE EMFIELDS SETUP: Pmax:",wdns_P_max,p_c

     if(em_field_type==1 .and. constrained_transport_scheme .ne. 3) then
        write(*,*) "Sorry, em_field_type==1 (toroidal fields) not yet supported for constrained_transport_scheme != 3."
        stop
     end if

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

     write(*,*) "HELLO BEFORE _P gupij_wdns calculation!",cctk_lsh,gxx(9,9,9)
     write(*,*) "HELLO BEFORE _P gupij_wdns calculation: _P LEVEL!",cctk_lsh,gxx_p(9,9,9)
     write(*,*) "HELLO BEFORE _P gupij_wdns calculation: _P_P LEVEL!",cctk_lsh,gxx_p_p(9,9,9)

     ! First compute gupij, based on gij_p:
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

     ! Set the st_i's:
     st_x = mhd_st_x_p
     st_y = mhd_st_y_p
     st_z = mhd_st_z_p

     !
     ! Compute the vector potential A_phi
     !
     if (em_field_type .ne. 1) then 
        call WDNS_compute_Aphi(ext,X,Y,Z,PhysicalRadius,P,A_phi,Ax_p,Ay_p,Az_p, & 
             mhd_psi_line, mhd_u_psi, mhd_chi_line, mhd_u_chi, &
             betam1,p_c,x_NS_CoM_coord,y_NS_CoM_coord,wdns_P_max,Sym_Bz, & 
             enable_trace_field_line,constrained_transport_scheme,em_field_type)
     end if

     !
     ! Compute the vector potential Ax and Ay for near toroidal configuration
     !
     if (em_field_type==1) then 
        call WDNS_compute_A_toroidal(ext,X,Y,Z,P,Ax_p,Ay_p,Az_p,phi_p,betam1, &
             p_c,x_NS_CoM_coord,y_NS_CoM_coord,r0,wdns_P_max)
     end if

     write(*,*) "HELLO BEFORE _P WDNS_compute_Bi_and_update_conservatives"

     ! Compute B^i's based on values of A_i's.  
     !   Note that the conservatives up to this point have assumed B=0.  Here we
     !   update the conservatives to include the now nonzero magnetic terms.
     call WDNS_compute_Bi_and_update_conservatives(ext,X,Y,Z,PhysicalRadius, &
          P,rho_b,vx,vy,vz,u0, &
          phi_p,gxx_p,gxy_p,gxz_p,gyy_p,gyz_p,gzz_p,shiftx_p,shifty_p,shiftz_p,lapm1_p, &
          tau_p,mhd_st_x_p,mhd_st_y_p,mhd_st_z_p,st_x,st_y,st_z, &
          A_phi,A_phix,A_phiy,A_phiz, & 
          Ax_p,Ay_p,Az_p, & 
          Bx,By,Bz, & 
          temp1,temp2,temp3, & 
          x_NS_CoM_coord,y_NS_CoM_coord,rho_b_atm, &
          constrained_transport_scheme,em_field_type)
     
     ! Reset magnetic fields to zero, so that in the below primitives solver call, the conservatives & primitives are consistent with B=0.
     !   Remember that the magnetic fields are assumed zero until they are computed inside WDNS_compute_Bi_and_update_conservatives() on each
     !   timelevel. 
     Bx = 0.D0
     By = 0.D0
     Bz = 0.D0

     write(*,*) "HELLO AFTER _P WDNS_compute_Bi_and_update_conservatives"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! NEXT, _P_P LEVEL:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! Compute gupij, based on gij_p_p:
     temp1 = 1.D0/(gxx_p_p * gyy_p_p * gzz_p_p + gxy_p_p * gyz_p_p * gxz_p_p + gxz_p_p * gxy_p_p * gyz_p_p &
          - gxz_p_p * gyy_p_p * gxz_p_p - gxy_p_p * gxy_p_p * gzz_p_p - gxx_p_p * gyz_p_p * gyz_p_p)

      gupxx =   ( gyy_p_p * gzz_p_p - gyz_p_p * gyz_p_p )* temp1
      gupxy = - ( gxy_p_p * gzz_p_p - gyz_p_p * gxz_p_p )* temp1
      gupxz =   ( gxy_p_p * gyz_p_p - gyy_p_p * gxz_p_p )* temp1
      gupyy =   ( gxx_p_p * gzz_p_p - gxz_p_p * gxz_p_p )* temp1
      gupyz = - ( gxx_p_p * gyz_p_p - gxy_p_p * gxz_p_p )* temp1
      gupzz =   ( gxx_p_p * gyy_p_p - gxy_p_p * gxy_p_p )* temp1

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
        call WDNS_compute_Aphi(ext,X,Y,Z,PhysicalRadius,P,A_phi,Ax_p_p,Ay_p_p,Az_p_p, & 
             mhd_psi_line, mhd_u_psi, mhd_chi_line, mhd_u_chi, &
             betam1,p_c,x_NS_CoM_coord,y_NS_CoM_coord,wdns_P_max,Sym_Bz, & 
             enable_trace_field_line,constrained_transport_scheme,em_field_type)
     end if

     !
     ! Compute the vector potential Ax and Ay for near toroidal configuration
     !
     if (em_field_type==1) then 
        call WDNS_compute_A_toroidal(ext,X,Y,Z,P,Ax_p_p,Ay_p_p,Az_p_p,phi_p_p,betam1, &
             p_c,x_NS_CoM_coord,y_NS_CoM_coord,r0,wdns_P_max)
     end if

     ! Compute B^i's based on values of A_i's.  
     !   Note that the conservatives up to this point have assumed B=0.  Here we
     !   update the conservatives to include the now nonzero magnetic terms.
     call WDNS_compute_Bi_and_update_conservatives(ext,X,Y,Z,PhysicalRadius, &
          P,rho_b,vx,vy,vz,u0, &
          phi_p_p,gxx_p_p,gxy_p_p,gxz_p_p,gyy_p_p,gyz_p_p,gzz_p_p,shiftx_p_p,shifty_p_p,shiftz_p_p,lapm1_p_p, &
          tau_p_p,mhd_st_x_p_p,mhd_st_y_p_p,mhd_st_z_p_p,st_x,st_y,st_z, &
          A_phi,A_phix,A_phiy,A_phiz, & 
          Ax_p_p,Ay_p_p,Az_p_p, & 
          Bx,By,Bz, & 
          temp1,temp2,temp3, & 
          x_NS_CoM_coord,y_NS_CoM_coord,rho_b_atm, &
          constrained_transport_scheme,em_field_type)

     ! Reset magnetic fields to zero, so that in the below primitives solver call, the conservatives & primitives are consistent with B=0.
     !   Remember that the magnetic fields are assumed zero until they are computed inside WDNS_compute_Bi_and_update_conservatives() on each
     !   timelevel.
     Bx = 0.D0
     By = 0.D0
     Bz = 0.D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! FINALLY, LATEST LEVEL (neither _p nor _p_p):
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! Update gupij and primitives on the latest level (neither _p nor _p_p):
     temp1 = 1.D0/(gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz &
          - gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz)

      gupxx =   ( gyy * gzz - gyz * gyz )* temp1
      gupxy = - ( gxy * gzz - gyz * gxz )* temp1
      gupxz =   ( gxy * gyz - gyy * gxz )* temp1
      gupyy =   ( gxx * gzz - gxz * gxz )* temp1
      gupyz = - ( gxx * gyz - gxy * gxz )* temp1
      gupzz =   ( gxx * gyy - gxy * gxy )* temp1

     if(primitives_solver==11) then

        !SKIP INITIAL GUESS.  IT'S NOT NEEDED ANYWAY.

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
     end if

     ! Set the st_i's:
     st_x = mhd_st_x
     st_y = mhd_st_y
     st_z = mhd_st_z

     !
     ! Compute the vector potential A_phi
     !
     if (em_field_type .ne. 1) then 
        call WDNS_compute_Aphi(ext,X,Y,Z,PhysicalRadius,P,A_phi,Ax,Ay,Az, & 
             mhd_psi_line, mhd_u_psi, mhd_chi_line, mhd_u_chi, &
             betam1,p_c,x_NS_CoM_coord,y_NS_CoM_coord,wdns_P_max,Sym_Bz, & 
             enable_trace_field_line,constrained_transport_scheme,em_field_type)
     end if

     !
     ! Compute the vector potential Ax and Ay for near toroidal configuration
     !
     if (em_field_type==1) then 
        call WDNS_compute_A_toroidal(ext,X,Y,Z,P,Ax,Ay,Az,phi,betam1, &
             p_c,x_NS_CoM_coord,y_NS_CoM_coord,r0,wdns_P_max)
     end if

     ! Compute B^i's based on values of A_i's.  
     !   Note that the conservatives up to this point have assumed B=0.  Here we
     !   update the conservatives to include the now nonzero magnetic terms.
     call WDNS_compute_Bi_and_update_conservatives(ext,X,Y,Z,PhysicalRadius, &
          P,rho_b,vx,vy,vz,u0, &
          phi,gxx,gxy,gxz,gyy,gyz,gzz,shiftx,shifty,shiftz,lapm1, &
          tau,mhd_st_x,mhd_st_y,mhd_st_z,st_x,st_y,st_z, &
          A_phi,A_phix,A_phiy,A_phiz, & 
          Ax,Ay,Az, & 
          Bx,By,Bz, & 
          temp1,temp2,temp3, & 
          x_NS_CoM_coord,y_NS_CoM_coord,rho_b_atm, &
          constrained_transport_scheme,em_field_type)

     ! Set the gauge variable to zero.
     psi6phi = 0.D0
     psi6phi_p = 0.D0
     psi6phi_p_p = 0.D0

     ! Compute B from A
     call wdns_compute_B_from_A(CCTK_PASS_FTOF)

     ! Free up memory:
     deallocate(A_phi, A_phix, A_phiy, A_phiz)

  end if
end subroutine WDNS_setup_emfield_part2_local_setBi_after_checkpoint



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine WDNS_compute_Bi_and_update_conservatives(ext,X,Y,Z,PhysicalRadius, &
     P,rho_b,vx,vy,vz,u0, &
     phi,gxx,gxy,gxz,gyy,gyz,gzz,shiftx,shifty,shiftz,lapm1, &
     tau,mhd_st_x,mhd_st_y,mhd_st_z,st_x,st_y,st_z, &
     A_phi,A_phix,A_phiy,A_phiz, & 
     Ax,Ay,Az, & 
     Bx,By,Bz, & 
     temp1,temp2,temp3, & 
     x_NS_CoM_coord,y_NS_CoM_coord,rho_b_atm, &
     constrained_transport_scheme,em_field_type)

  ! INPUT/OUTPUT PARAMETERS
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: X,Y,Z,PhysicalRadius
  real*8, dimension(ext(1),ext(2),ext(3))  :: P,rho_b,vx,vy,vz,u0
  real*8, dimension(ext(1),ext(2),ext(3))  :: phi,gxx,gxy,gxz,gyy,gyz,gzz,shiftx,shifty,shiftz,lapm1
  real*8, dimension(ext(1),ext(2),ext(3))  :: tau,mhd_st_x,mhd_st_y,mhd_st_z,st_x,st_y,st_z
  real*8, dimension(ext(1),ext(2),ext(3))  :: A_phi,A_phix,A_phiy,A_phiz
  real*8, dimension(ext(1),ext(2),ext(3))  :: Ax,Ay,Az
  real*8, dimension(ext(1),ext(2),ext(3))  :: Bx,By,Bz
  real*8, dimension(ext(1),ext(2),ext(3))  :: temp1,temp2,temp3
  real*8                                   :: x_NS_CoM_coord,y_NS_CoM_coord,rho_b_atm
  integer 				   :: constrained_transport_scheme,em_field_type

  ! INTERNAL PARAMETERS
  integer                                  :: i,j,k
  integer                                  :: imin,imax,jmin,jmax,kmin,kmax
  real*8                                   :: dX,dY,dZ,psim6,psim6_s,pomega2,al,sqrtg,sqrtg4,B2s
  real*8                                   :: fs4pi,B_xs,B_ys,B_zs,sb0,sb2,xn,yn
  real*8                                   :: sb_x,sb_y,sb_z,psi4,u_x,u_y,u_z
  real*8                                   :: psin
  real*8                                   :: Yijk,Yijkp1,Yijp1k,Yijp1kp1,Yip1jk,Yip1jkp1,Yip1jp1k,Yip1jp1kp1
  real*8                                   :: Xijk,Xijkp1,Xijp1k,Xijp1kp1,Xip1jk,Xip1jkp1,Xip1jp1k,Xip1jp1kp1
  real*8				   :: Ax000,Ax001,Ax010,Ax011,Ax100,Ax101,Ax110,Ax111
  real*8                                   :: Ay000,Ay001,Ay010,Ay011,Ay100,Ay101,Ay110,Ay111

  integer                                  :: AXISYM,EQUATORIAL,OCTANT
  parameter(EQUATORIAL = 1, OCTANT = 2, AXISYM = 4)
  real*8, parameter                        :: SYM = 1.d0, ANTI = -1.d0

  dX = X(2,1,1) - X(1,1,1)
  dY = Y(1,2,1) - Y(1,1,1)
  dZ = Z(1,1,2) - Z(1,1,1)

  imin = 1
  jmin = 1
  kmin = 1
  imax = ext(1)
  jmax = ext(2)
  kmax = ext(3)

  fs4pi = sqrt(4.d0*acos(-1.d0))

  ! Compute the derivatives of A_phi, if we are not using constrained_transport_scheme==3
  !
  if (constrained_transport_scheme .ne. 3) then 
     write(*,*) 'SORRY, CONSTRAINED_TRANSPORT_SCHEME != 3 NOT SUPPORTED WHEN SETTING B FIELDS AFTER CHECKPOINT RESTART.'
     write(*,*) 'You will need to compute Bitildes at the end of the ...setBi_after_checkpoint.F90 routine!'
     stop
  end if

  ! Compute B^i on staggered grid and temporarily store them in A_phix,A_phiy,A_phiz
  !
  if (constrained_transport_scheme==3) then 
     do k=1,ext(3)
        do j=1,ext(2)
           do i=1,ext(1)
              im1 = max(i-1,1)
              jm1 = max(j-1,1)
              km1 = max(k-1,1)
              ip1 = min(i+1,ext(1))
              jp1 = min(j+1,ext(2))
              kp1 = min(k+1,ext(3))

              psim6_s = exp(-3.d0 * (phi(i,j,k) + phi(ip1,j,k)) )
              A_phix(i,j,k) = ( (Az(i,j,k)-Az(i,jm1,k))/dY   & 
                   - (Ay(i,j,k)-Ay(i,j,km1))/dZ ) * psim6_s

              psim6_s = exp(-3.d0 * (phi(i,j,k) + phi(i,jp1,k)) )
              A_phiy(i,j,k) = ( (Ax(i,j,k)-Ax(i,j,km1))/dZ & 
                   - (Az(i,j,k)-Az(im1,j,k))/dX ) * psim6_s

              psim6_s = exp(-3.d0 * (phi(i,j,k) + phi(i,j,kp1)) )
              A_phiz(i,j,k) = ( (Ay(i,j,k)-Ay(im1,j,k))/dX &
                   - (Ax(i,j,k)-Ax(i,jm1,k))/dY ) * psim6_s

           end do
        end do
     end do

     ! Now compute B^i on unstaggered grid by simple averge
     do k=1,ext(3)
        do j=1,ext(2)
           do i=1,ext(1)
              im1 = max(i-1,1)
              jm1 = max(j-1,1)
              km1 = max(k-1,1)
              Bx(i,j,k) = 0.5d0* (A_phix(i,j,k) + A_phix(im1,j,k))
              By(i,j,k) = 0.5d0* (A_phiy(i,j,k) + A_phiy(i,jm1,k))
              Bz(i,j,k) = 0.5d0* (A_phiz(i,j,k) + A_phiz(i,j,km1))
           end do
        end do
     end do
  end if

  !
  ! Now compute B^i according to (exercise for the readers)
  !  B^x = (-x/pomega^2) e^(-6 phi) * A_{phi,z}; 
  !  B^y = (-y/pomega^2) e^(-6 phi) * A_{phi,z};
  !  B^z = e^(-6 phi) * (x A_{phi,x} + y A_{phi,y})/pomega^2; 
  !  pomega^2 = x^2 + y^2
  !
  ! and then calculate mhd_st_i and tau
  !
  do k = kmin,kmax
     do j = jmin,jmax
        do i = imin,imax
           xn = X(i,j,k) - x_NS_CoM_coord
           yn = Y(i,j,k) - y_NS_CoM_coord

           psim6 = exp(-6.d0*phi(i,j,k))
           pomega2 = xn**2 + yn**2
	   if (constrained_transport_scheme==3) then 
	      ! do nothing since B^i has been computed 
   	   elseif (em_field_type==0 .or. i==imax .or. j==jmax .or. k==kmax) then 
              Bx(i,j,k) = -xn/pomega2 * psim6 * A_phiz(i,j,k)
              By(i,j,k) = -yn/pomega2 * psim6 * A_phiz(i,j,k)
              Bz(i,j,k) = psim6/pomega2 * (xn*A_phix(i,j,k) +  & 
                   yn*A_phiy(i,j,k))
	   else
	      ! Compute Bx, By, Bz that satisfy div(B)=0 to machine precision.
  	      ! Here is the recipe, which has been verified by Mathematica: 
	      !
	      ! \tilde{B}^x(i,j,k)=-(Ay(i,j,k+1)-Ay(i,j,k) + Ay(i+1,j,k+1)-Ay(i+1,j,k) + Ay(i,j+1,k+1)-Ay(i,j+1,k) + Ay(i+1,j+1,k+1)-Ay(i+1,j+1,k))/(4*dz) + (Az(i,j+1,k)-Az(i,j,k) + Az(i+1,j+1,k)-Az(i+1,j,k) + Az(i,j+1,k+1)-Az(i,j,k+1) + Az(i+1,j+1,k+1)-Az(i+1,j,k+1))/(4*dy)
	      ! \tilde{B}^y(i,j,k)=(Ax(i,j,k+1)-Ax(i,j,k) + Ax(i+1,j,k+1)-Ax(i+1,j,k) + Ax(i,j+1,k+1)-Ax(i,j+1,k) + Ax(i+1,j+1,k+1)-Ax(i+1,j+1,k))/(4*dz) - (Az(i+1,j,k)-Az(i,j,k) + Az(i+1,j+1,k)-Az(i,j+1,k) + Az(i+1,j,k+1)-Az(i,j,k+1) + Az(i+1,j+1,k+1)-Az(i,j+1,k+1))/(4*dx)
	      ! \tilde{B}^z(i,j,k):=(Ay(i+1,j,k)-Ay(i,j,k) + Ay(i+1,j+1,k)-Ay(i,j+1,k) + Ay(i+1,j,k+1)-Ay(i,j,k+1) + Ay(i+1,j+1,k+1)-Ay(i,j+1,k+1))/(4*dx) - (Ax(i,j+1,k)-Ax(i,j,k) + Ax(i+1,j+1,k)-Ax(i+1,j,k) + Ax(i,j+1,k+1)-Ax(i,j,k+1) + Ax(i+1,j+1,k+1)-Ax(i+1,j,k+1))/(4*dy)
              ! Here \tilde{B}^i = psi^6 B^i, and 
	      ! Ax(i,j,k), Ay(i,j,k), and Az(i,j,k) are the 
              ! 3 (covariant) components of the vector potential at the corner 
	      ! point (x_i-dx/2, y_j-dy/2, z_k-dz/2).
    	      ! In the present case, we set Ax = -A_phi * y/(x^2+y^2), 
	      !  Ay = A_phi * x/(x^2+y^2), and Az = 0.
	      ! Note that B^i can't be set this way at the boundary points
	      ! becuase of the sentcil structure. Presumably, this 
	      ! is not a problem when we have enough ghostzones.
	      !
	      ! In the following, Axabc denotes Ax(i+a,j+b,k+c) and 
	      ! similarly for Ayabc.
	      ! 

              ! Here we set the coordinates, relative to the origin, which in this case is the center of the NS.
              Yijk      = yn
              Yijkp1    = yn
              Yijp1k    = yn + dY
              Yijp1kp1  = yn + dY
              Yip1jk    = yn
              Yip1jkp1  = yn
              Yip1jp1k  = yn + dY
              Yip1jp1kp1= yn + dY

              Xijk      = xn
              Xijkp1    = xn
              Xijp1k    = xn
              Xijp1kp1  = xn
              Xip1jk    = xn + dX
              Xip1jkp1  = xn + dX
              Xip1jp1k  = xn + dX
              Xip1jp1kp1= xn + dX

              Ax000 = -A_phi(i,j,k)* (Yijk-0.5d0*dY)/ & 
                   ((Xijk-0.5d0*dX)**2 + (Yijk-0.5d0*dY)**2)
              Ax001 = -A_phi(i,j,k+1)*(Yijkp1-0.5d0*dY)/ & 
                   ((Xijkp1-0.5d0*dX)**2 + (Yijkp1-0.5d0*dY)**2)
              Ax010 = -A_phi(i,j+1,k)* (Yijp1k-0.5d0*dY)/ & 
                   ((Xijp1k-0.5d0*dX)**2 + (Yijp1k-0.5d0*dY)**2)
              Ax011 = -A_phi(i,j+1,k+1)*(Yijp1kp1-0.5d0*dY)/ & 
                   ((Xijp1kp1-0.5d0*dX)**2 + (Yijp1kp1-0.5d0*dY)**2)
              Ax100 = -A_phi(i+1,j,k)*(Yip1jk-0.5d0*dY)/ & 
                   ((Xip1jk-0.5d0*dX)**2 + (Yip1jk-0.5d0*dY)**2)
              Ax101 = -A_phi(i+1,j,k+1)*(Yip1jkp1-0.5d0*dY)/ & 
                   ((Xip1jkp1-0.5d0*dX)**2 + (Yip1jkp1-0.5d0*dY)**2)
              Ax110 = -A_phi(i+1,j+1,k)*(Yip1jp1k-0.5d0*dY)/ & 
                   ((Xip1jp1k-0.5d0*dX)**2 + (Yip1jp1k-0.5d0*dY)**2)
              Ax111 = -A_phi(i+1,j+1,k+1)*(Yip1jp1kp1-0.5d0*dY)/ & 
                   ((Xip1jp1kp1-0.5d0*dX)**2 + (Yip1jp1kp1-0.5d0*dY)**2)
              Ay000 = A_phi(i,j,k)* (Xijk-0.5d0*dX)/ &
                   ((Xijk-0.5d0*dX)**2 + (Yijk-0.5d0*dY)**2)
              Ay001 = A_phi(i,j,k+1)*(Xijkp1-0.5d0*dX)/ &
                   ((Xijkp1-0.5d0*dX)**2 + (Yijkp1-0.5d0*dY)**2)
              Ay010 = A_phi(i,j+1,k)* (Xijp1k-0.5d0*dX)/ &
                   ((Xijp1k-0.5d0*dX)**2 + (Yijp1k-0.5d0*dY)**2)
              Ay011 = A_phi(i,j+1,k+1)*(Xijp1kp1-0.5d0*dX)/ &
                   ((Xijp1kp1-0.5d0*dX)**2 + (Yijp1kp1-0.5d0*dY)**2)
              Ay100 = A_phi(i+1,j,k)*(Xip1jk-0.5d0*dX)/ &
                   ((Xip1jk-0.5d0*dX)**2 + (Yip1jk-0.5d0*dY)**2)
              Ay101 = A_phi(i+1,j,k+1)*(Xip1jkp1-0.5d0*dX)/ &
                   ((Xip1jkp1-0.5d0*dX)**2 + (Yip1jkp1-0.5d0*dY)**2)
              Ay110 = A_phi(i+1,j+1,k)*(Xip1jp1k-0.5d0*dX)/ &
                   ((Xip1jp1k-0.5d0*dX)**2 + (Yip1jp1k-0.5d0*dY)**2)
              Ay111 = A_phi(i+1,j+1,k+1)*(Xip1jp1kp1-0.5d0*dX)/ &
                   ((Xip1jp1kp1-0.5d0*dX)**2 + (Yip1jp1kp1-0.5d0*dY)**2)

	      Bx(i,j,k) = -( (Ay001-Ay000) + (Ay101-Ay100) + (Ay011-Ay010) + (Ay111-Ay110) ) * 0.25d0/dZ * psim6
	      By(i,j,k) = ( (Ax001-Ax000) + (Ax101-Ax100) + (Ax011-Ax010) + (Ax111-Ax110) ) * 0.25d0/dZ * psim6
	      Bz(i,j,k) = ( (Ay100-Ay000) + (Ay110-Ay010) + (Ay101-Ay001) + (Ay111-Ay011) ) * 0.25d0/dX * psim6 - ( (Ax010-Ax000) + (Ax110-Ax100) + (Ax011-Ax001) + (Ax111-Ax101) ) * 0.25d0/dY * psim6
	   end if
           ! Compute b^0 and b_i
           al = 1.d0 + lapm1(i,j,k)
           sqrtg = 1.d0/psim6
           sqrtg4 = al * sqrtg
           B2s = exp(4.d0*phi(i,j,k))*(gxx(i,j,k)*Bx(i,j,k)**2 + & 
		2.d0*gxy(i,j,k)*Bx(i,j,k)*By(i,j,k) + & 
		2.d0*gxz(i,j,k)*Bx(i,j,k)*Bz(i,j,k) + &
        	gyy(i,j,k)*By(i,j,k)**2 + 2.d0*gyz(i,j,k)*By(i,j,k)*Bz(i,j,k) + & 
		gzz(i,j,k)*Bz(i,j,k)**2)/(fs4pi*al)**2
           psin = exp(4.d0*phi(i,j,k))/al/fs4pi
           B_xs  = psin * (gxx(i,j,k) * Bx(i,j,k) + gxy(i,j,k) * By(i,j,k) + & 
                gxz(i,j,k) * Bz(i,j,k))
           B_ys  = psin * (gxy(i,j,k) * Bx(i,j,k) + gyy(i,j,k) * By(i,j,k) + & 
                gyz(i,j,k) * Bz(i,j,k))
           B_zs  = psin * (gxz(i,j,k) * Bx(i,j,k) + gyz(i,j,k) * By(i,j,k) + & 
                gzz(i,j,k) * Bz(i,j,k))

           psin = psi4*u0(i,j,k)
           u_x = ( gxx(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) +  &
                gxy(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) +  & 
                gxz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) )*psin
           u_y = ( gxy(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) +  & 
                gyy(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) +  & 
                gyz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) )*psin
           u_z = ( gxz(i,j,k)*(shiftx(i,j,k)+vx(i,j,k)) +  &
                gyz(i,j,k)*(shifty(i,j,k)+vy(i,j,k)) +  & 
                gzz(i,j,k)*(shiftz(i,j,k)+vz(i,j,k)) )*psin


           sb0 = (u_x*Bx(i,j,k) + u_y*By(i,j,k) + &
                u_z*Bz(i,j,k))/fs4pi/al
           sb2 = (B2s + sb0**2)/u0(i,j,k)**2
           !Branson's way of ensuring b^2/P is what we want: 
           !  Call this routine once to calibrate sb2/P, then average it.
           !   Scale sb2 accordingly, then call this function again to obtain the desired sb2/P
           ! Here's the necessary line of code.  You'll need to add the input parameters as well.
           !bsq(i,j,k) = sb2/P(i,j,k)

           sb_x = (B_xs + u_x*sb0)/u0(i,j,k)
           sb_y = (B_ys + u_y*sb0)/u0(i,j,k)
           sb_z = (B_zs + u_z*sb0)/u0(i,j,k)
           ! Now compute mhd_st_i and tau
           mhd_st_x(i,j,k) = st_x(i,j,k)+sqrtg4*(u0(i,j,k)*sb2*u_x-sb0*sb_x)
           mhd_st_y(i,j,k) = st_y(i,j,k)+sqrtg4*(u0(i,j,k)*sb2*u_y-sb0*sb_y)
           mhd_st_z(i,j,k) = st_z(i,j,k)+sqrtg4*(u0(i,j,k)*sb2*u_z-sb0*sb_z)
           tau(i,j,k) = tau(i,j,k) + sqrtg*( sb2*(al*u0(i,j,k))**2 &
                - sb2*0.5d0 - (al*sb0)**2 )

           ! CHECK NEXT FUNCTION IN SCHEDULE.CCL for reason why temp1 and temp2 are defined.
           ! Temporarily store e^(6 phi)*b^2/2 (b^2/2 = P_mag) to temp1, and e^(6 phi)*P(i,j,k) to temp2, e^(6 phi) to temp3 
           temp1(i,j,k) = sb2*0.5d0*sqrtg
           !if(sb2.gt.0.D0) write(*,*) sb2
           temp2(i,j,k) = P(i,j,k)*sqrtg
           if (rho_b(i,j,k) .gt. rho_b_atm*1.d5) then
              temp3(i,j,k) = 1.d0
           else
              temp3(i,j,k) = 0.d0
           end if
        end do
     end do
  end do

end subroutine WDNS_compute_Bi_and_update_conservatives
