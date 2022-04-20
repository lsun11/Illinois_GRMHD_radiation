#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


subroutine bhns_initialdata_read_binfile_driver(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext

  real*8 		:: dX,dY,dZ,superposition_factor,eta_falloffr,eta_falloffdr,inner_value,outer_value, psiL
  real*8                :: gamma1, gamma2, gamma3
  real*8                :: kappa1, kappa2, kappa3
  real*8                :: rhoo1, rhoo2
  integer, pointer      :: nphase(:)
  real(8), pointer      :: abc(:), abi(:), rhoi(:), tedi(:)
  real(8), pointer      :: loc_rho(:), loc_pre(:), loc_eps(:), loc_kappa(:), loc_gamma(:)
  real*8                :: an1, an2, an3, c2
  real*8                :: gxxL, gxyL, gxzL, gyyL, gyzL, gzzL, psiLm4, detgij
  real*8                :: gupxxL, gupxyL, gupxzL, gupyyL, gupyzL, gupzzL
  real*8                :: shiftxL,shiftyL,shiftzL,lapl
  real*8                :: aaa, bbb, ccc
  integer               :: i,ierr,j,k, ia
  integer               :: n1, n2
  write(*,*) "If Valgrind gives an error right before this line, ignore the error."

  Allocate(nphase(0:0));   nphase(0) = 2

  write(6,*) "Size of gamma_tab :", size(gamma_tab)
  write(6,*) "Size of k_tab     :", size(k_tab)
  write(6,*) "Size of rho_tab   :", size(rho_tab)
  write(6,*) "Size of P_tab     :", size(P_tab)
  write(6,*) "Size of eps_tab   :", size(eps_tab)

  n1 = size(gamma_tab)
  n2 = size(rho_tab)

  Allocate(abc(1:n1));   Allocate(abi(1:n1));  Allocate(rhoi(1:n2));   Allocate(tedi(1:n2))
  abc=0.0d0;             abi=0.0d0;            rhoi=0.0d0;             tedi=0.0d0

  Allocate(loc_rho(1:n2));       loc_rho = 0.0d0
  Allocate(loc_pre(1:n2));       loc_pre = 0.0d0
  Allocate(loc_eps(1:n2));       loc_eps = 0.0d0
  Allocate(loc_kappa(1:n1));   loc_kappa = 0.0d0
  Allocate(loc_gamma(1:n1));   loc_gamma = 0.0d0



  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  ext = cctk_lsh

  !Set initial BH position to zero.
  bh_posn_x = 0.D0
  bh_posn_y = 0.D0
  bh_posn_z = 0.D0
  
 
  write(*,*) "n_ijk:",cctk_lsh
  write(*,*) "ZZZ:",Z(1,1,:)
  !  write(*,*) "X:",X(:,1,1)

  if(use_new_bhns_initial_data==0) then
     call read_inputfile_bhns(cctkGH,genID_cmdline_output_enable,fisheye_enable,reset_shift_lapse, &
          X(1,1,1),Y(1,1,1),Z(1,1,1), &
          dX,dY,dZ,cctk_lsh, &
          lapm1,shiftx,shifty,shiftz, &
          phi,Axx,Axy,Axz,Ayy,Ayz, &
          rho_b,vx,vy,vz, &
          gamma_th,bh_posn_x(1),BigM)
  else if(use_new_bhns_initial_data==1) then
     call read_inputfile_bhns_bhspin(cctkGH,genID_cmdline_output_enable,fisheye_enable, &
          X(1,1,1),Y(1,1,1),Z(1,1,1), &
          dX,dY,dZ,cctk_lsh, &
          lapm1,shiftx,shifty,shiftz, &
          phi,Axx,Axy,Axz,Ayy,Ayz,Azz, &
          rho_b,vx,vy,vz, &
          gamma_th,bh_posn_x(1),BigM,initial_ns_coord_x,initial_ns_coord_y,0)

     ! Call the following if you want a rotated metric superposition:
     if(1==0) then
        call read_inputfile_bhns_bhspin(cctkGH,genID_cmdline_output_enable,fisheye_enable, &
             X(1,1,1),Y(1,1,1),Z(1,1,1), &
             dX,dY,dZ,cctk_lsh, &
             temp1,temp2,temp3,temp4, &
             temp5,Kxx,Kxy,Kxz,Kyy,Kyz,Kzz, &
             temp6,temp7,temp8,temp9, &
             gamma_th,bh_posn_x(1),BigM,initial_ns_coord_x,initial_ns_coord_y,1)

        temp11 = 1.D0

        !superposition_factor = 0.80D0
        !write(*,*) "SUPERPOSITION FACTOR=",superposition_factor

        ! Note that z -> y, and y -> -z in the rotated data
        lapm1 = temp11*lapm1 + (1.D0 - temp11)*temp1
        shiftx = temp11*shiftx + (1.D0 - temp11)*temp2
        shifty = temp11*shifty + (1.D0 - temp11)*(-1.D0)*temp4
        shiftz = temp11*shiftz + (1.D0 - temp11)*temp3

        phi = temp11*phi + (1.D0 - temp11)*temp5

        Axx = temp11*Axx + (1.D0 - temp11)*Kxx
        Axy = temp11*Axy + (1.D0 - temp11)*(-1.D0)*Kxz
        Axz = temp11*Axz + (1.D0 - temp11)*Kxy
        Ayy = temp11*Ayy + (1.D0 - temp11)*Kzz
        Ayz = temp11*Ayz + (1.D0 - temp11)*(-1.D0)*Kyz
        Azz = temp11*Azz + (1.D0 - temp11)*Kyy

        rho_b = temp11*rho_b + (1.D0 - temp11)*temp6
        vx = temp11*vx + (1.D0 - temp11)*temp7
        vy = temp11*vy + (1.D0 - temp11)*(-1.D0)*temp9
        vz = temp11*vz + (1.D0 - temp11)*temp8
     end if

! read Lorene ID for NSNS
  else if(use_new_bhns_initial_data==3) then
     call read_inputfile_wdns_lorene2(cctkGH,genID_cmdline_output_enable,fisheye_enable,reset_shift_lapse, &
           X(1,1,1),Y(1,1,1),Z(1,1,1), &
           dX,dY,dZ,cctk_lsh, &
           lapm1,shiftx,shifty,shiftz, &
           phi,Axx,Axy,Axz,Ayy,Ayz, &
           rho_b,vx,vy,vz, &
           gamma1, gamma2, gamma3, kappa1, kappa2, kappa3, rhoo1, rhoo2)

     print *,"WDNS:: done with Lorene. We need to compute more stuff ... "

! read COCAL ID for NSNS
  else if(use_new_bhns_initial_data==4) then
     call read_inputfile_cocal(cctkGH, num_CO, piecewise,  &
           genID_cmdline_output_enable,fisheye_enable,reset_shift_lapse, &
           X(1,1,1),Y(1,1,1),Z(1,1,1), &
           dX,dY,dZ,cctk_lsh, &
           lapm1,shiftx,shifty,shiftz, &
           phi,Axx,Axy,Axz,Ayy,Ayz, &
           rho_b,vx,vy,vz, &
           nphase, abc, abi, rhoi, tedi)

     write(6,*) "~~~~~~~~~~~~~~~~~~   FROM read_inputfile_driver.F90   ~~~~~~~~~~~~~~~~~~~~~~~~~"
     write(6,*) "Number of phases nphase=", nphase(0) 
     write(6,'(a4,2a18)') "i", "Gamma", "Kappa"
     do ia=1,n1
       write(6,'(i4,1p,2e18.10)') ia, abi(ia), abc(ia)
     end do
     write(6,'(a4,3a18)') "i", "rho0", "rho"
     do ia=1,n2
       write(6,'(i4,1p,3e18.10)') ia, rhoi(ia), tedi(ia)
     end do
     write(6,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

     kappa1 = abc(1)      !  used below: k_poly=kappa1    gamma_th=gamma1
     gamma1 = abi(1)

     print *,"kappa1= ",kappa1
     print *,"gamma1= ",gamma1

     print *,"Cocal ID. We need to compute more stuff ... "


! read COCAL ID for black hole-disk systems
  else if(use_new_bhns_initial_data==6) then

!  Notice that the 3-metric is the conformal metric
!  kij is the full extrinsic curvature
!  and the metric is not conformally flat and trk is not zero
  
     call read_inputfile_BH_disk(cctkGH, num_CO, piecewise,          &
           genID_cmdline_output_enable,fisheye_enable,reset_shift_lapse,  &
           X(1,1,1),Y(1,1,1),Z(1,1,1),                                    &
           dX,dY,dZ,cctk_lsh,                                             &
           lapm1,shiftx,shifty,shiftz,                                    &
           gxx,gxy,gxz,gyy,gyz,gzz,                                       &
           gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,                           &
           kxx,kxy,kxz,kyy,kyz,kzz,                                       &
           phi,trK,rho_b,vx,vy,vz,                                        &
           nphase, abc, abi, rhoi, tedi)

     write(6,*) "~~~~~~~~~~~~~~~~~~   FROM read_inputfile_driver.F90   ~~~~~~~~~~~~~~~~~~~~~~~~~"
     write(6,*) "Number of phases nphase=", nphase(0) 
     write(6,'(a4,2a18)') "i", "Gamma", "Kappa"
     do ia=1,n1
       write(6,'(i4,1p,2e18.10)') ia, abi(ia), abc(ia)
     end do
     write(6,'(a4,3a18)') "i", "rho0", "rho"
     do ia=1,n2
       write(6,'(i4,1p,3e18.10)') ia, rhoi(ia), tedi(ia)
     end do
     write(6,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

     kappa1 = abc(1)      !  used below: k_poly=kappa1    gamma_th=gamma1
     gamma1 = abi(1)

     print *,"kappa1= ",kappa1
     print *,"gamma1= ",gamma1

        
     print *,"**********************************************************************"
     print *,"BH+disk Cocal ID has been read in. We need to compute more stuff ... "
     print *,"**********************************************************************"

  else if(use_new_bhns_initial_data==2) then
     call read_inputfile_SpEC_bhns(cctkGH,genID_cmdline_output_enable,fisheye_enable, &
          X(1,1,1),Y(1,1,1),Z(1,1,1), &
          dX,dY,dZ,cctk_lsh, &
          lapm1,shiftx,shifty,shiftz, &
          gxx,gxy,gxz,gyy,gyz,gzz, phi, &
          kxx,kxy,kxz,kyy,kyz,kzz, &
          rho_b,vx,vy,vz, &
          gamma_th,k_poly,bh_posn_x(1),bh_posn_y(1),BigM,initial_ns_coord_x,initial_ns_coord_y,0)  

     ! Set the conformal factor now
     ! Note that here we store the oconformal factor to phi
     if (genID_cmdline_output_enable.ne.1) then
        do k=1,cctk_lsh(3)
           do j=1,cctk_lsh(2)
              do i=1,cctk_lsh(1)
                 gxxL = gxx(i,j,k)
                 gxyL = gxy(i,j,k)
                 gxzL = gxz(i,j,k)
                 gyyL = gyy(i,j,k)
                 gyzL = gyz(i,j,k)
                 gzzL = gzz(i,j,k)
                 detgij =  gxxL * gyyL * gzzL + gxyL * gyzL * gxzL + gxzL * gxyL * gyzL &
                      - gxzL * gyyL * gxzL - gxyL * gxyL * gzzL - gxxL * gyzL * gyzL
                 psiL = detgij**(1.d0/12.d0)
                 phi(i,j,k) = psiL
                 if (detgij <= 0.d0) then
                    write(*,*) "Stopping: Negative 3-metric determinant. Something's wrong!"
                    write(*,*) "X, Y, Z=", X(i,j,k), Y(i,j,k), Z(i,j,k)
                    write(*,*) "gxx, gyy, gzz=", gxxL, gyyL, gzzL
                    write(*,*) "gxy, gxz, gyz=", gxyL, gxzL, gyzL
                    stop
                 end if
                 psiLm4=1.d0/(psiL*psiL*psiL*psiL)
                 ! We do not know whether the data have trK=0 so we compute trK here, remember Aij is the full extrinsic curvature so far
                 trK(i,j,k) = gxxL*kxx(i,j,k) + gyyL*kyy(i,j,k) + gzzL*kzz(i,j,k) + & 
                      2.d0*(gxyL*kxy(i,j,k) + gxzL*kxz(i,j,k) + gyzL*kyz(i,j,k) ) 
                 
                 ! Also the 3 metric here is not conformally flat, so we set the BSSN metric, now
                 gxx(i,j,k) = gxxL*psiLm4
                 gxy(i,j,k) = gxyL*psiLm4
                 gxz(i,j,k) = gxzL*psiLm4
                 gyy(i,j,k) = gyyL*psiLm4
                 gyz(i,j,k) = gyzL*psiLm4
                 gzz(i,j,k) = gzzL*psiLm4
              end do
           end do
        end do
     end if
        
  end if

! define lapse -1 array

  if(genID_cmdline_output_enable.ne.1) then
     lapm1 = lapm1 - 1.D0
  else
     lapm1=0.D0

     shiftx=0.D0
     shifty=0.D0
     shiftz=0.D0

     phi = 1.D0

     gxx = 1.D0
     gxy = 0.D0
     gxz = 0.D0
     gyy = 1.D0
     gyz = 0.D0
     gzz = 1.D0

     Axx = 0.D0
     Axy = 0.D0
     Axz = 0.D0
     Ayy = 0.D0
     Ayz = 0.D0
     Azz = 0.D0

     rho_b = 1.D-7
     vx = 0.D0
     vy = 0.D0
     vz = 0.D0

     gamma_th = 1.3D0
       
  end if
     
  ! K_poly=1 because Keisuke's datafiles are in polytropic units already!!!
  if (use_new_bhns_initial_data<2) then
     k_poly=1.0d0
  end if
 
  if((use_new_bhns_initial_data.eq.3).or.(use_new_bhns_initial_data.eq.4).or.(use_new_bhns_initial_data.eq.6)) then
     k_poly=kappa1
     gamma_th=1.+1./gamma1
     ! note that I'm assuming that Cocal gives us gamma no n
     if((use_new_bhns_initial_data.eq.4).or.(use_new_bhns_initial_data.eq.6)) gamma_th=gamma1
  endif
  write(6,*)'gamma_th,BigM,bh_posn_x(1),k_poly = ',gamma_th,BigM,bh_posn_x(1), k_poly     
  
  n_poly=1.0d0/(gamma_th-1.0d0)
  write(6,*)'BH irreducible mass scaling factor:',BigM

  if((use_new_bhns_initial_data.eq.0).or.(use_new_bhns_initial_data.eq.1).or.(use_new_bhns_initial_data.eq.2)) then
     write(6,*)'BH position x=',bh_posn_x(1),' xbh/mirr=',bh_posn_x(1)/BigM
  end if

  if(piecewise.eq.0) then  
  
     !single polytrope
     neos=1

     rho_tab(1)=1.0

     P_tab(1)=K_poly

     eps_tab(1)=K_poly/(gamma_th-1.0d0)

     do i=1,2
        k_tab(i)=K_poly
        gamma_tab(i)=gamma_th
     enddo
 
  else
     ! piecewise polytrope
     gamma_th= 1.66d0

!     neos=2
     neos = nphase(0)    ! **** nphase = neos  ****

     write(6,*) "=============", neos
     write(6,*) "======= neos, k_tab, gamma_tab, rho_tab, P_tab, eps_tab ======"
     do ia=1, neos+1
       write(6,'(i4,1p,5e18.10)') neos, k_tab(ia), gamma_tab(ia), rho_tab(ia), P_tab(ia), eps_tab(ia)
     end do

     do ia=1, neos+1
       loc_gamma(ia) = abi(ia)   !  abi(0) at surface
       loc_kappa(ia) = abc(ia)
!       loc_rho(ia) = rhoi(ia)
     enddo

     do ia=1, neos
       ! Transition densities, pressures
       loc_rho(ia) = rhoi(ia)
       loc_pre(ia) = loc_kappa(ia)*(loc_rho(ia)**loc_gamma(ia))
     enddo
     
!     loc_pre(neos+1) =  loc_kappa(neos)*(loc_rho(neos+1)**loc_gamma(neos))

     loc_eps(1) = loc_pre(1)/loc_rho(1)/(loc_gamma(1)-1.0d0)
!     print  *,"loc_eps1,2 = before loop",loc_eps(1),loc_eps(2), neos

     do ia=2,neos
       ! Transition specific internal energies
       loc_eps(ia) = loc_eps(ia-1) + (loc_pre(ia)/loc_rho(ia) - loc_pre(ia-1)/loc_rho(ia-1))/(loc_gamma(ia)-1.0d0)
     enddo
    
    
!     print *,"after loop", loc_eps(1), loc_eps(2),neos
     write(6,*) "======= neos, loc_kappa, loc_gamma, loc_rho, loc_pre, loc_eps ======"
     do ia=1, neos+1
       write(6,'(i4,1p,5e18.10)') neos, loc_kappa(ia), loc_gamma(ia), loc_rho(ia), loc_pre(ia), loc_eps(ia)
     end do

     gamma_tab = loc_gamma
     k_tab     = loc_kappa
     P_tab     = loc_pre
     rho_tab   = loc_rho
     eps_tab   = loc_eps


     print *, "Diagnostic for piecewise polytropes."
     write(6,'(2a18)') "Gamma", "Kappa"
     do ia=1, n1
       write(6,'(1p,2e18.10)') gamma_tab(ia), k_tab(ia) 
     enddo
     write(6,'(3a18)') "rho0", "P", "epsilon"
     do ia=1, n2
       write(6,'(1p,3e18.10)') rho_tab(ia), P_tab(ia), eps_tab(ia)                               
     enddo
  end if

  print *,"Done with bhns_initialdata_read_binfile_driver"

  deallocate(nphase)
  deallocate(abc)  
  deallocate(abi)  
  deallocate(rhoi)  
  deallocate(tedi)  

  deallocate(loc_rho)
  deallocate(loc_pre)
  deallocate(loc_eps)
  deallocate(loc_kappa)
  deallocate(loc_gamma)
  
end subroutine bhns_initialdata_read_binfile_driver


