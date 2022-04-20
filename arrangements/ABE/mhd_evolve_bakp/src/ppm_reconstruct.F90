!-----------------------------------------------------------------!
! Calculate a_{j-1/2}, set a_R, a_L equal to this.
!-----------------------------------------------------------------!
!
subroutine ppm_find_face_vals(ext,X,Y,Z,rho,delta_rho,rhor,rhol,P,delta_P,Pr,Pl, &
     vx,delta_vx,vxr,vxl,vy,delta_vy,vyr,vyl,vz,delta_vz,vzr,vzl,m,Symmetry)
 implicit none
 integer, dimension(3)                          :: ext
 real*8,  dimension(ext(1),ext(2),ext(3))       :: X,Y,Z
 real*8, dimension(ext(1),ext(2),ext(3))        :: rho,delta_rho,rhor,rhol
 real*8, dimension(ext(1),ext(2),ext(3))        :: P,delta_P,Pr,Pl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vx,delta_vx,vxr,vxl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vy,delta_vy,vyr,vyl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vz,delta_vz,vzr,vzl
 integer                                        :: m,Symmetry
!
 integer                                        :: j,jmin,jmax
 real*8                     :: SYM, ANTI, HALF, SIXTH, EIGHTH, ZERO
 parameter(HALF = 1.0D0/2.D0, SIXTH = 1.D0/6.D0, EIGHTH = 1.D0/8.D0, SYM = 1.D0, ANTI = -1.D0, ZERO = 0.D0)
 integer                            :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
 parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

 jmin = lbound(rho,m)
 jmax = ubound(rho,m)
 if(Symmetry == OCTANT) then 
    if(m==1 .and. X(1,1,1).lt.ZERO) then
       jmin = lbound(rho,1)+1
    else if(m==2 .and. Y(1,1,1).lt.ZERO) then
       jmin = lbound(rho,2)+1
    else if(m==3 .and. Z(1,1,1).lt.ZERO) then
       jmin = lbound(rho,3)+1     
    end if
 else if(Symmetry == AXISYM) then 
    if(m==1 .and. X(1,1,1).lt.ZERO) then
       jmin = lbound(rho,1)+1
    else if(m==3 .and. Z(1,1,1).lt.ZERO) then
       jmin = lbound(rho,3)+1     
    end if    
 else if(Symmetry == EQUATORIAL) then 
    if(m==3 .and. Z(1,1,1).lt.ZERO) then
       jmin = lbound(rho,3)+1     
    end if    
 end if

 if(m==1) then
    rhol(jmin,:,:) = SYM*rho(jmin,:,:) + HALF*(rho(jmin,:,:)-SYM*rho(jmin,:,:)) &
            + EIGHTH*(ANTI*delta_rho(jmin,:,:) - delta_rho(jmin,:,:))
    Pl(jmin,:,:) = SYM*P(jmin,:,:) + HALF*(P(jmin,:,:)-SYM*P(jmin,:,:)) &
            + EIGHTH*(ANTI*delta_P(jmin,:,:) - delta_P(jmin,:,:))
    if(Symmetry==NO_SYMM.or.Symmetry==EQUATORIAL) then
       vxl(jmin,:,:) = SYM*vx(jmin,:,:) + HALF*(vx(jmin,:,:)-SYM*vx(jmin,:,:)) &
            + EIGHTH*(ANTI*delta_vx(jmin,:,:) - delta_vx(jmin,:,:))
    else
       vxl(jmin,:,:) = ANTI*vx(jmin,:,:) + HALF*(vx(jmin,:,:)-ANTI*vx(jmin,:,:)) &
            + EIGHTH*(SYM*delta_vx(jmin,:,:) - delta_vx(jmin,:,:))
    end if
    if(Symmetry==AXISYM) then
       vyl(jmin,:,:) = ANTI*vy(jmin,:,:) + HALF*(vy(jmin,:,:)-ANTI*vy(jmin,:,:)) &
            + EIGHTH*(SYM*delta_vy(jmin,:,:) - delta_vy(jmin,:,:))
    else
       vyl(jmin,:,:) = SYM*vy(jmin,:,:) + HALF*(vy(jmin,:,:)-SYM*vy(jmin,:,:)) &
            + EIGHTH*(ANTI*delta_vy(jmin,:,:) - delta_vy(jmin,:,:))
    end if
    vzl(jmin,:,:) = SYM*vz(jmin,:,:) + HALF*(vz(jmin,:,:)-SYM*vz(jmin,:,:)) &
            + EIGHTH*(ANTI*delta_vz(jmin,:,:) - delta_vz(jmin,:,:))
    do j=jmin+1,jmax
       rhol(j,:,:) = rho(j-1,:,:) + HALF*(rho(j,:,:)-rho(j-1,:,:)) &
            + EIGHTH*(delta_rho(j-1,:,:) - delta_rho(j,:,:))
       Pl(j,:,:) = P(j-1,:,:) + HALF*(P(j,:,:)-P(j-1,:,:)) &
            + EIGHTH*(delta_P(j-1,:,:) - delta_P(j,:,:))
       vxl(j,:,:) = vx(j-1,:,:) + HALF*(vx(j,:,:)-vx(j-1,:,:)) &
            + EIGHTH*(delta_vx(j-1,:,:) - delta_vx(j,:,:))
       vyl(j,:,:) = vy(j-1,:,:) + HALF*(vy(j,:,:)-vy(j-1,:,:)) &
            + EIGHTH*(delta_vy(j-1,:,:) - delta_vy(j,:,:))
       vzl(j,:,:) = vz(j-1,:,:) + HALF*(vz(j,:,:)-vz(j-1,:,:)) &
            + EIGHTH*(delta_vz(j-1,:,:) - delta_vz(j,:,:))
    end do
 
    do j=jmin, jmax-1
       rhor(j,:,:) = rhol(j+1,:,:)
       Pr(j,:,:)   = Pl(j+1,:,:)
       vxr(j,:,:)  = vxl(j+1,:,:)
       vyr(j,:,:)  = vyl(j+1,:,:)
       vzr(j,:,:)  = vzl(j+1,:,:)
    end do
    rhor(jmax,:,:) = rho(jmax,:,:)
    Pr(jmax,:,:)   = P(jmax,:,:)
    vxr(jmax,:,:)  = vx(jmax,:,:)
    vyr(jmax,:,:)  = vy(jmax,:,:)
    vzr(jmax,:,:)  = vz(jmax,:,:)

 else if(m==2) then
    rhol(:,jmin,:) = SYM*rho(:,jmin,:) + HALF*(rho(:,jmin,:)-SYM*rho(:,jmin,:)) &
            + EIGHTH*(ANTI*delta_rho(:,jmin,:) - delta_rho(:,jmin,:))
    Pl(:,jmin,:) = SYM*P(:,jmin,:) + HALF*(P(:,jmin,:)-SYM*P(:,jmin,:)) &
            + EIGHTH*(ANTI*delta_P(:,jmin,:) - delta_P(:,jmin,:))
    vxl(:,jmin,:) = SYM*vx(:,jmin,:) + HALF*(vx(:,jmin,:)-SYM*vx(:,jmin,:)) &
         + EIGHTH*(ANTI*delta_vx(:,jmin,:) - delta_vx(:,jmin,:))
    if(Symmetry==NO_SYMM.or.Symmetry==EQUATORIAL) then
       vyl(:,jmin,:) = SYM*vy(:,jmin,:) + HALF*(vy(:,jmin,:)-SYM*vy(:,jmin,:)) &
            + EIGHTH*(ANTI*delta_vy(:,jmin,:) - delta_vy(:,jmin,:))
    else
       vyl(:,jmin,:) = ANTI*vy(:,jmin,:) + HALF*(vy(:,jmin,:)-ANTI*vy(:,jmin,:)) &
            + EIGHTH*(SYM*delta_vy(:,jmin,:) - delta_vy(:,jmin,:))
    end if
    vzl(:,jmin,:) = SYM*vz(:,jmin,:) + HALF*(vz(:,jmin,:)-SYM*vz(:,jmin,:)) &
            + EIGHTH*(ANTI*delta_vz(:,jmin,:) - delta_vz(:,jmin,:))
    do j=jmin+1,jmax
       rhol(:,j,:) = rho(:,j-1,:) + HALF*(rho(:,j,:)-rho(:,j-1,:)) &
            + EIGHTH*(delta_rho(:,j-1,:) - delta_rho(:,j,:))
       Pl(:,j,:) = P(:,j-1,:) + HALF*(P(:,j,:)-P(:,j-1,:)) &
            + EIGHTH*(delta_P(:,j-1,:) - delta_P(:,j,:))
       vxl(:,j,:) = vx(:,j-1,:) + HALF*(vx(:,j,:)-vx(:,j-1,:)) &
            + EIGHTH*(delta_vx(:,j-1,:) - delta_vx(:,j,:))
       vyl(:,j,:) = vy(:,j-1,:) + HALF*(vy(:,j,:)-vy(:,j-1,:)) &
            + EIGHTH*(delta_vy(:,j-1,:) - delta_vy(:,j,:))
       vzl(:,j,:) = vz(:,j-1,:) + HALF*(vz(:,j,:)-vz(:,j-1,:)) &
            + EIGHTH*(delta_vz(:,j-1,:) - delta_vz(:,j,:))
    end do
    
    do j=jmin, jmax-1
       rhor(:,j,:) = rhol(:,j+1,:)
       Pr(:,j,:)   = Pl(:,j+1,:)
       vxr(:,j,:)  = vxl(:,j+1,:)
       vyr(:,j,:)  = vyl(:,j+1,:)
       vzr(:,j,:)  = vzl(:,j+1,:)
    end do
    rhor(:,jmax,:) = rho(:,jmax,:)
    Pr(:,jmax,:)   = P(:,jmax,:)
    vxr(:,jmax,:)  = vx(:,jmax,:)
    vyr(:,jmax,:)  = vy(:,jmax,:)
    vzr(:,jmax,:)  = vz(:,jmax,:)

 else
    rhol(:,:,jmin) = SYM*rho(:,:,jmin) + HALF*(rho(:,:,jmin)-SYM*rho(:,:,jmin)) &
            + EIGHTH*(ANTI*delta_rho(:,:,jmin) - delta_rho(:,:,jmin))
    Pl(:,:,jmin) = SYM*P(:,:,jmin) + HALF*(P(:,:,jmin)-SYM*P(:,:,jmin)) &
            + EIGHTH*(ANTI*delta_P(:,:,jmin) - delta_P(:,:,jmin))
    vxl(:,:,jmin) = SYM*vx(:,:,jmin) + HALF*(vx(:,:,jmin)-SYM*vx(:,:,jmin)) &
         + EIGHTH*(ANTI*delta_vx(:,:,jmin) - delta_vx(:,:,jmin))
    vyl(:,:,jmin) = SYM*vy(:,:,jmin) + HALF*(vy(:,:,jmin)-SYM*vy(:,:,jmin)) &
         + EIGHTH*(ANTI*delta_vy(:,:,jmin) - delta_vy(:,:,jmin))
    if(Symmetry==NO_SYMM .or. Z(1,1,jmin).lt.ZERO) then
       vzl(:,:,jmin) = SYM*vz(:,:,jmin) + HALF*(vz(:,:,jmin)-SYM*vz(:,:,jmin)) &
            + EIGHTH*(ANTI*delta_vz(:,:,jmin) - delta_vz(:,:,jmin))
    else
       vzl(:,:,jmin) = ANTI*vz(:,:,jmin) + HALF*(vz(:,:,jmin)-ANTI*vz(:,:,jmin)) &
            + EIGHTH*(SYM*delta_vz(:,:,jmin) - delta_vz(:,:,jmin))
    end if
    do j=jmin+1,jmax
       rhol(:,:,j) = rho(:,:,j-1) + HALF*(rho(:,:,j)-rho(:,:,j-1)) &
            + EIGHTH*(delta_rho(:,:,j-1) - delta_rho(:,:,j))
       Pl(:,:,j) = P(:,:,j-1) + HALF*(P(:,:,j)-P(:,:,j-1)) &
            + EIGHTH*(delta_P(:,:,j-1) - delta_P(:,:,j))
       vxl(:,:,j) = vx(:,:,j-1) + HALF*(vx(:,:,j)-vx(:,:,j-1)) &
            + EIGHTH*(delta_vx(:,:,j-1) - delta_vx(:,:,j))
       vyl(:,:,j) = vy(:,:,j-1) + HALF*(vy(:,:,j)-vy(:,:,j-1)) &
            + EIGHTH*(delta_vy(:,:,j-1) - delta_vy(:,:,j))
       vzl(:,:,j) = vz(:,:,j-1) + HALF*(vz(:,:,j)-vz(:,:,j-1)) &
            + EIGHTH*(delta_vz(:,:,j-1) - delta_vz(:,:,j))
    end do
    
    do j=jmin, jmax-1
       rhor(:,:,j) = rhol(:,:,j+1)
       Pr(:,:,j)   = Pl(:,:,j+1)
       vxr(:,:,j)  = vxl(:,:,j+1)
       vyr(:,:,j)  = vyl(:,:,j+1)
       vzr(:,:,j)  = vzl(:,:,j+1)
    end do
    rhor(:,:,jmax) = rho(:,:,jmax)
    Pr(:,:,jmax)   = P(:,:,jmax)
    vxr(:,:,jmax)  = vx(:,:,jmax)
    vyr(:,:,jmax)  = vy(:,:,jmax)
    vzr(:,:,jmax)  = vz(:,:,jmax)

 end if
 

end subroutine ppm_find_face_vals

!-----------------------------------------------------------------!
! steepen rho_b
!-----------------------------------------------------------------!
!
subroutine ppm_steepen_rho(ext,X,Y,Z,rho,d1rho,d2rho,rhor_ppm,rhol_ppm, &
     rhor_mc,rhol_mc,P,n,rho_max,m,Symmetry)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8, dimension(ext(1),ext(2),ext(3))        :: rho,d1rho,d2rho,P
 real*8, dimension(ext(1),ext(2),ext(3))        :: rhor_ppm,rhol_ppm
 real*8, dimension(ext(1),ext(2),ext(3))        :: rhor_mc,rhol_mc
 integer                                        :: m,Symmetry
 real*8                                         :: n,rho_max
!
 integer                                        :: i,j,k
 integer                    :: imin,jmin,kmin,imax,jmax,kmax
 real*8                     :: SYM, ANTI, HALF, SIXTH, ONE, ZERO
 real*8                     :: K0, eta1, eta2, epsilon, Gamma, TINY
 real*8                     :: contact_disc_check,second_dir_check,relative_change_check
 real*8                     :: eta_tilde,eta
 integer                            :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
 parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
 parameter(ONE = 1.D0, ZERO = 0.D0)
 parameter(HALF = ONE/2.D0, SIXTH = ONE/6.D0, SYM = ONE, ANTI = -ONE)
! standard Colella-Woodward parameters:
!    K0 = 0.1d0, eta1 = 20.d0, eta2 = 0.05d0, epsilon = 0.01d0
! our parameters:
 parameter(K0 = 0.1d0, eta1 = 20.d0, eta2 = 0.05d0, epsilon = 0.01d0)
! The following values are Font-like:
! parameter(K0 = 1.d0, eta1 = 5.d0, eta2 = 0.05d0, epsilon = 0.1d0)
 imax = ubound(rho,1)
 imin = lbound(rho,1)
 jmax = ubound(rho,2)
 jmin = lbound(rho,2)
 kmax = ubound(rho,3)
 kmin = lbound(rho,3)
 Gamma = ONE + ONE/n
 TINY = 1.d-15*rho_max

 if(Symmetry == OCTANT) then 
    if(X(1,1,1).lt.ZERO) then
       imin = imin + 1
    end if
    if(Y(1,1,1).lt.ZERO) then
       jmin = jmin + 1
    end if
    if(Z(1,1,1).lt.ZERO) then
       kmin = kmin + 1
    end if
 else if(Symmetry == AXISYM) then 
    if(X(1,1,1).lt.ZERO) then
       imin = imin + 1
    end if
    if(Z(1,1,1).lt.ZERO) then
       kmin = kmin + 1
    end if
 else if(Symmetry == EQUATORIAL) then 
    if(Z(1,1,1).lt.ZERO) then
       kmin = kmin + 1
    end if
 end if


 do k=kmin,kmax
    do j=jmin,jmax
       do i=imin,imax
          
          if(m==1) then
             if(i==imin .or. i==imax) then
                contact_disc_check = -ONE
                second_dir_check = -ONE
                relative_change_check = -ONE
             else
                if ((min(rho(i+1,j,k),rho(i-1,j,k)) .lt. TINY) .or. &
                     (min(P(i+1,j,k),P(i-1,j,k)) .lt. TINY)) then
                   contact_disc_check = -ONE
                else
                   contact_disc_check = Gamma*K0*abs(rho(i+1,j,k) &
                        -rho(i-1,j,k))/min(rho(i+1,j,k),rho(i-1,j,k)) &
                        - abs(P(i+1,j,k)-P(i-1,j,k))/min(P(i+1,j,k),P(i-1,j,k))
                end if
                second_dir_check = -d2rho(i+1,j,k)*d2rho(i-1,j,k)
                relative_change_check = abs(2.d0*d1rho(i,j,k)) &
                     - epsilon*min(rho(i+1,j,k),rho(i-1,j,k))
             end if
             
             if(contact_disc_check.ge.ZERO .and. second_dir_check.ge.ZERO &
                  .and. relative_change_check.ge.ZERO) then
                if (abs(d1rho(i,j,k)) .gt. ZERO) then
                   eta_tilde = -SIXTH*(d2rho(i+1,j,k)-d2rho(i-1,j,k))/d1rho(i,j,k)
                else
                   eta_tilde = ZERO
                end if
                eta = max(ZERO,min(eta1*(eta_tilde - eta2),ONE))
                rhol_ppm(i,j,k) = rhol_ppm(i,j,k)*(ONE-eta) &
                     + rhor_mc(i,j,k)*eta
                rhor_ppm(i,j,k) = rhor_ppm(i,j,k)*(ONE-eta) &
                     + rhol_mc(i+1,j,k)*eta
             end if

          else if(m==2) then
             if(j==jmin .or. j==jmax) then
                contact_disc_check = -ONE
                second_dir_check = -ONE
                relative_change_check = -ONE
             else
                if ((min(rho(i,j+1,k),rho(i,j-1,k)) .lt. TINY) .or. &
                     (min(P(i,j+1,k),P(i,j-1,k)) .lt. TINY)) then
                   contact_disc_check = -ONE
                else
                   contact_disc_check = Gamma*K0*abs(rho(i,j+1,k) &
                        -rho(i,j-1,k))/min(rho(i,j+1,k),rho(i,j-1,k)) &
                        - abs(P(i,j+1,k)-P(i,j-1,k))/min(P(i,j+1,k),P(i,j-1,k))
                end if
                second_dir_check = -d2rho(i,j+1,k)*d2rho(i,j-1,k)
                relative_change_check = abs(2.d0*d1rho(i,j,k)) &
                     - epsilon*min(rho(i,j+1,k),rho(i,j-1,k))
             end if
             
             if(contact_disc_check.ge.ZERO .and. second_dir_check.ge.ZERO &
                  .and. relative_change_check.ge.ZERO) then
                if (abs(d1rho(i,j,k)) .gt. ZERO) then
                   eta_tilde = -SIXTH*(d2rho(i,j+1,k)-d2rho(i,j-1,k))/d1rho(i,j,k)
                else
                   eta_tilde = ZERO
                end if
                eta = max(ZERO,min(eta1*(eta_tilde - eta2),ONE))
                rhol_ppm(i,j,k) = rhol_ppm(i,j,k)*(ONE-eta) &
                     + rhor_mc(i,j,k)*eta
                rhor_ppm(i,j,k) = rhor_ppm(i,j,k)*(ONE-eta) &
                     + rhol_mc(i,j+1,k)*eta
             end if

          else if(m==3) then
             if(k==kmin .or. k==kmax) then
                contact_disc_check = -ONE
                second_dir_check = -ONE
                relative_change_check = -ONE
             else
                if ((min(rho(i,j,k+1),rho(i,j,k-1)) .lt. TINY) .or. &
                     (min(P(i,j,k+1),P(i,j,k-1)) .lt. TINY)) then
                   contact_disc_check = -ONE
                else
                   contact_disc_check = Gamma*K0*abs(rho(i,j,k+1) &
                        -rho(i,j,k-1))/min(rho(i,j,k+1),rho(i,j,k-1)) &
                        - abs(P(i,j,k+1)-P(i,j,k-1))/min(P(i,j,k+1),P(i,j,k-1))
                end if
                second_dir_check = -d2rho(i,j,k+1)*d2rho(i,j,k-1)
                relative_change_check = abs(2.d0*d1rho(i,j,k)) &
                     - epsilon*min(rho(i,j,k+1),rho(i,j,k-1))
             end if
             
             if(contact_disc_check.ge.ZERO .and. second_dir_check.ge.ZERO &
                  .and. relative_change_check.ge.ZERO) then
                if (abs(d1rho(i,j,k)) .gt. ZERO) then
                   eta_tilde = -SIXTH*(d2rho(i,j,k+1)-d2rho(i,j,k-1))/d1rho(i,j,k)
                else
                   eta_tilde = ZERO
                end if
                eta = max(ZERO,min(eta1*(eta_tilde - eta2),ONE))
                rhol_ppm(i,j,k) = rhol_ppm(i,j,k)*(ONE-eta) &
                     + rhor_mc(i,j,k)*eta
                rhor_ppm(i,j,k) = rhor_ppm(i,j,k)*(ONE-eta) &
                     + rhol_mc(i,j,k+1)*eta
             end if

          end if

       end do
    end do
 end do

end subroutine ppm_steepen_rho
!-----------------------------------------------------------------!
! ensure monotonicity of the interpolating polynomial             !
!-----------------------------------------------------------------!
!
subroutine ppm_monotonize(ext,X,Y,Z,rhomax,rho_b,rho_br,rho_bl, &
     d02rho_b_m,d2rho_b_bck,d2rho_b_fwd, P,Pr,Pl,d02P_m,d2P_bck,d2P_fwd, &
     vx,vxr,vxl,d02vx_m,d2vx_bck,d2vx_fwd, vy,vyr,vyl,d02vy_m, &
     d2vy_bck,d2vy_fwd, vz,vzr,vzl,d02vz_m,d2vz_bck,d2vz_fwd, &
     m,Symmetry, Reconstruction)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: rho_b,rho_br,rho_bl
 real*8, dimension(ext(1),ext(2),ext(3))        :: P,Pr,Pl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vx,vxr,vxl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vy,vyr,vyl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vz,vzr,vzl
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8, dimension(ext(1),ext(2),ext(3)),intent(in) :: d02rho_b_m,d02P_m
 real*8, dimension(ext(1),ext(2),ext(3)),intent(in) :: d02vx_m,d02vy_m,d02vz_m
 real*8, dimension(ext(1),ext(2),ext(3)),intent(in) :: d2rho_b_bck,d2rho_b_fwd
 real*8, dimension(ext(1),ext(2),ext(3)),intent(in) :: d2P_bck,d2P_fwd
 real*8, dimension(ext(1),ext(2),ext(3)),intent(in) :: d2vx_bck,d2vx_fwd
 real*8, dimension(ext(1),ext(2),ext(3)),intent(in) :: d2vy_bck,d2vy_fwd
 real*8, dimension(ext(1),ext(2),ext(3)),intent(in) :: d2vz_bck,d2vz_fwd
 integer, intent(in)                            :: m,Symmetry,Reconstruction
 integer                                        :: i,j,k
 integer                                        :: imin,imax,jmin,jmax,kmin,kmax
 real*8                                         :: SYM, ANTI, rhomax, rhomax_eff
 parameter(SYM = 1.D0, ANTI = -1.D0)
 integer :: PPM_PLUS, PPM
 integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
 parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
 parameter(PPM_PLUS = 1, PPM = 2)
 integer, dimension(ext(1),ext(2),ext(3))  :: instar
 kmax = ubound(rho_b,3)
 kmin = lbound(rho_b,3)

 if(Symmetry == OCTANT) then 
    if(m==1 .and. X(1,1,1) .lt. 0.D0) then
       kmin = kmin + 1
    end if
    if(m==2 .and. Y(1,1,1) .lt. 0.D0) then
       kmin = kmin + 1
    end if
    if(m==3 .and. Z(1,1,1) .lt. 0.D0) then
       kmin = kmin + 1
    end if
 else if(Symmetry == AXISYM) then
    if(m==1 .and. X(1,1,1) .lt. 0.D0) then
       kmin = kmin + 1
    end if
    if(m==3 .and. Z(1,1,1) .lt. 0.D0) then
       kmin = kmin + 1
    end if    
 else if(Symmetry == EQUATORIAL) then
    if(m==3 .and. Z(1,1,1) .lt. 0.D0) then
       kmin = kmin + 1
    end if    
 end if

 if(Reconstruction == PPM_PLUS) then
    where(rho_b.gt.0.05d0*rhomax)
       instar = 1
    elsewhere
       instar = 0
    end where
    rhomax_eff = rhomax
 else
    instar = 0
    rhomax_eff = 2.d0 * rhomax
 end if

 if(m==1) then
    call monotonize_rho(ext,X,Y,Z,rhomax_eff,rho_b,rho_br,rho_bl, &
         d02rho_b_m,d2rho_b_bck,d2rho_b_fwd,m,SYM,instar,Symmetry)
    call monotonize_P(ext,X,Y,Z,rhomax_eff,rho_b,P,Pr,Pl,d02P_m,d2P_bck,d2P_fwd,m,SYM,instar,Symmetry)
    call monotonize(ext,X,Y,Z,vz,vzr,vzl,d02vz_m,d2vz_bck,d2vz_fwd,m,SYM,instar,Symmetry)
    if(Symmetry==OCTANT) then
       call monotonize(ext,X,Y,Z,vx,vxr,vxl,d02vx_m,d2vx_bck,d2vx_fwd,m,ANTI,instar,Symmetry)
       call monotonize(ext,X,Y,Z,vy,vyr,vyl,d02vy_m,d2vy_bck,d2vy_fwd,m,SYM,instar,Symmetry)
    else if(Symmetry==AXISYM) then
       call monotonize(ext,X,Y,Z,vx,vxr,vxl,d02vx_m,d2vx_bck,d2vx_fwd,m,ANTI,instar,Symmetry)
       call monotonize(ext,X,Y,Z,vy,vyr,vyl,d02vy_m,d2vy_bck,d2vy_fwd,m,ANTI,instar,Symmetry)
    else
       call monotonize(ext,X,Y,Z,vx,vxr,vxl,d02vx_m,d2vx_bck,d2vx_fwd,m,SYM,instar,Symmetry)
       call monotonize(ext,X,Y,Z,vy,vyr,vyl,d02vy_m,d2vy_bck,d2vy_fwd,m,SYM,instar,Symmetry)
    end if
 else if(m==2) then
    if(Symmetry.ne.AXISYM) then
       call monotonize_rho(ext,X,Y,Z,rhomax_eff,rho_b,rho_br,rho_bl, &
            d02rho_b_m,d2rho_b_bck,d2rho_b_fwd,m,SYM,instar,Symmetry)
       call monotonize_P(ext,X,Y,Z,rhomax_eff,rho_b,P,Pr,Pl,d02P_m,d2P_bck,d2P_fwd,m,SYM,instar,Symmetry)
       call monotonize(ext,X,Y,Z,vx,vxr,vxl,d02vx_m,d2vx_bck,d2vx_fwd,m,SYM,instar,Symmetry)
       call monotonize(ext,X,Y,Z,vz,vzr,vzl,d02vz_m,d2vz_bck,d2vz_fwd,m,SYM,instar,Symmetry)
       if(Symmetry==OCTANT) then
          call monotonize(ext,X,Y,Z,vy,vyr,vyl,d02vy_m,d2vy_bck,d2vy_fwd,m,ANTI,instar,Symmetry)
       else
          call monotonize(ext,X,Y,Z,vy,vyr,vyl,d02vy_m,d2vy_bck,d2vy_fwd,m,SYM,instar,Symmetry)
       end if
    end if
 else if(m==3) then
    call monotonize_rho(ext,X,Y,Z,rhomax_eff,rho_b,rho_br,rho_bl, &
         d02rho_b_m,d2rho_b_bck,d2rho_b_fwd,m,SYM,instar,Symmetry)
    call monotonize_P(ext,X,Y,Z,rhomax_eff,rho_b,P,Pr,Pl,d02P_m,d2P_bck,d2P_fwd,m,SYM,instar,Symmetry)
    call monotonize(ext,X,Y,Z,vx,vxr,vxl,d02vx_m,d2vx_bck,d2vx_fwd,m,SYM,instar,Symmetry)
    call monotonize(ext,X,Y,Z,vy,vyr,vyl,d02vy_m,d2vy_bck,d2vy_fwd,m,SYM,instar,Symmetry)
    if(Symmetry.eq.NO_SYMM .or. (Z(1,1,kmin) .lt. 0.d0)) then
       call monotonize(ext,X,Y,Z,vz,vzr,vzl,d02vz_m,d2vz_bck,d2vz_fwd,m,SYM,instar,Symmetry)
    else
       call monotonize(ext,X,Y,Z,vz,vzr,vzl,d02vz_m,d2vz_bck,d2vz_fwd,m,ANTI,instar,Symmetry)
    end if
 else
    write(*,*) "m not in proper range",m
 end if

! write(*,*) "inside ppm_monot",m,vxl(2,2,2),vxr(2,2,2),vyl(2,2,2),vyr(2,2,2),vzl(2,2,2),vzr(2,2,2),d02vx_m(2,2,2),d2vx_bck(2,2,2), &
!      d2vx_fwd(2,2,2),rho_bl(2,2,2),rho_br(2,2,2),rho_br(2,3,2),rhomax_eff,rho_b(2,2,2),rho_b(2,3,2),d02rho_b_m(2,2,2),d2rho_b_bck(2,2,2),d2rho_b_fwd(2,2,2),rhomax,vx(2,2,2),vy(2,2,2),vz(2,2,2),vx(2,3,2),vy(2,3,2),vz(2,3,2)
 

end subroutine ppm_monotonize
!-----------------------------------------------------------------!
! ensure monotonicity of the interpolating polynomial (Clark fix) !
!-----------------------------------------------------------------!
!
subroutine monotonize(ext,X,Y,Z,a,ar,al,d02a_m,d2a_bck,d2a_fwd,m,sym,instar,Symmetry)
 implicit none
 integer, dimension(3)                     :: ext
 real*8, dimension(ext(1),ext(2),ext(3))   :: a,ar,al
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8, dimension(ext(1),ext(2),ext(3)), intent(in) :: d02a_m,d2a_bck,d2a_fwd
 integer, intent(in)                       :: m
 real*8                                    :: sym
 integer                                   :: i,j,k
 integer                                   :: imin,imax,jmin,jmax,kmin,kmax
 real*8, parameter                         :: ZERO = 0.d0, HALF = 0.5d0, SIX = 6.d0
 real*8, parameter                         :: THREE = 3.d0, TWO = 2.d0
 real*8                                    :: q1, q2   !temporary storage variables
 integer, dimension(ext(1),ext(2),ext(3))  :: instar
 integer :: Symmetry, NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
 parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

 imax = ubound(a,1)
 imin = lbound(a,1)
 jmax = ubound(a,2)
 jmin = lbound(a,2)
 kmax = ubound(a,3)
 kmin = lbound(a,3) 

 if(Symmetry == OCTANT) then 
    if(X(1,1,1) .lt. ZERO) then
       imin = imin + 1
    end if
    if(Y(1,1,1) .lt. ZERO) then
       jmin = jmin + 1
    end if
    if(Z(1,1,1) .lt. ZERO) then
       kmin = kmin + 1
    end if
 else if(Symmetry == AXISYM) then
    if(X(1,1,1) .lt. ZERO) then
       imin = imin + 1
    end if
    if(Z(1,1,1) .lt. ZERO) then
       kmin = kmin + 1
    end if  
 else if(Symmetry == EQUATORIAL) then
    if(Z(1,1,1) .lt. ZERO) then
       kmin = kmin + 1
    end if  
 end if


 if (m==1) then

    do k=kmin,kmax
       do j=jmin,jmax
          do i=imin,imax-1

             if (i==imin) then
                if ( ( (sym*d02a_m(imin+1,j,k))*(sym*d02a_m(imin,j,k)) .gt. 0.d0 ) .and. &
                     ( (sym*d02a_m(imin,j,k))  * d02a_m(imin,j,k)      .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,k)        * d02a_m(i+1,j,k)    .gt. 0.d0 ) .and. &
                     ( d02a_m(i+1,j,k)      * d02a_m(i+2,j,k)    .gt. 0.d0 ) .and.instar(i,j,k)==1) then
                   ! do nothing
                else
                   if ((ar(i,j,k)-a(i,j,k))*(a(i,j,k)-al(i,j,k)) < ZERO) then
                      al(i,j,k) = a(i,j,k)
                      ar(i,j,k) = a(i,j,k)
                   else 
                      ! Otherwise, ensure monotonicity 
                      q1 = (ar(i,j,k)-al(i,j,k)) &
                           *(a(i,j,k) - HALF*(ar(i,j,k)+al(i,j,k)))
                      q2 = (ar(i,j,k)-al(i,j,k))**2/SIX
                      
                      if (q2 < q1) then
                         al(i,j,k) = THREE*a(i,j,k) - TWO*ar(i,j,k)
                      else if (-q2 > q1) then
                         ar(i,j,k) = THREE*a(i,j,k) - TWO*al(i,j,k)
                      end if
                   end if
                end if
                
             else
                if ( ( d2a_bck(i-1,j,k) * d02a_m(i-1,j,k)  .gt. 0.d0 ) .and. &
                     ( d02a_m(i-1,j,k) * d02a_m(i,j,k)    .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,k)   * d02a_m(i+1,j,k)  .gt. 0.d0 ) .and. &
                     ( d02a_m(i+1,j,k) * d2a_fwd(i+1,j,k)  .gt. 0.d0 ) .and.instar(i,j,k)==1) then
                   ! do nothing
                else
                   if ((ar(i,j,k)-a(i,j,k))*(a(i,j,k)-al(i,j,k)) < ZERO) then
                      al(i,j,k) = a(i,j,k)
                      ar(i,j,k) = a(i,j,k)
                   else 
                      ! Otherwise, ensure monotonicity 
                      q1 = (ar(i,j,k)-al(i,j,k)) &
                           *(a(i,j,k) - HALF*(ar(i,j,k)+al(i,j,k)))
                      q2 = (ar(i,j,k)-al(i,j,k))**2/SIX

                      if (q2 < q1) then
                         al(i,j,k) = THREE*a(i,j,k) - TWO*ar(i,j,k)
                      else if (-q2 > q1) then
                         ar(i,j,k) = THREE*a(i,j,k) - TWO*al(i,j,k)
                      end if
                   end if
                end if
             end if
             
             if (i==imin) then
                if (sym < ZERO) then
                   al(i,j,k) = ZERO
                end if
             end if
             
          end do
       end do
    end do

 else if (m==2) then

    do k=kmin,kmax
       do j=jmin,jmax-1
          do i=imin,imax

             if (j==jmin) then
                if ( ( (sym*d02a_m(i,jmin+1,k))*(sym*d02a_m(i,jmin,k)) .gt. 0.d0 ) .and. &
                     ( (sym*d02a_m(i,jmin,k))  * d02a_m(i,jmin,k)      .gt. 0.d0 ) .and. &
                     ( d02a_m(i,jmin,k)        * d02a_m(i,jmin+1,k)    .gt. 0.d0 .and.instar(i,j,k)==1) .and. &
                     ( d02a_m(i,jmin+1,k)      * d02a_m(i,jmin+2,k)    .gt. 0.d0 ) ) then
                   ! do nothing
                else
                   if ((ar(i,j,k)-a(i,j,k))*(a(i,j,k)-al(i,j,k)) < ZERO) then
                      al(i,j,k) = a(i,j,k)
                      ar(i,j,k) = a(i,j,k)
                   else 
                      ! Otherwise, ensure monotonicity 
                      q1 = (ar(i,j,k)-al(i,j,k))*(a(i,j,k) - HALF*(ar(i,j,k)+al(i,j,k)))
                      q2 = (ar(i,j,k)-al(i,j,k))**2/SIX
                      
                      if (q2 < q1) then
                         al(i,j,k) = THREE*a(i,j,k) - TWO*ar(i,j,k)
                      else if (-q2 > q1) then
                         ar(i,j,k) = THREE*a(i,j,k) - TWO*al(i,j,k)
                      end if
                   end if
                end if                
             else
                if ( ( d2a_bck(i,j-1,k) * d02a_m(i,j-1,k)  .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j-1,k) * d02a_m(i,j,k)    .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,k)   * d02a_m(i,j+1,k)  .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j+1,k) * d2a_fwd(i,j+1,k)  .gt. 0.d0 ) .and.instar(i,j,k)==1) then
                   ! do nothing
                else
                   if ((ar(i,j,k)-a(i,j,k))*(a(i,j,k)-al(i,j,k)) < ZERO) then
                      al(i,j,k) = a(i,j,k)
                      ar(i,j,k) = a(i,j,k)
                   else 
                      ! Otherwise, ensure monotonicity 
                      q1 = (ar(i,j,k)-al(i,j,k))*(a(i,j,k) - HALF*(ar(i,j,k)+al(i,j,k)))
                      q2 = (ar(i,j,k)-al(i,j,k))**2/SIX
                      
!                      if ((j /= 32) .and. (j /= 33)) then
                         if (q2 < q1) then
                            al(i,j,k) = THREE*a(i,j,k) - TWO*ar(i,j,k)
                         else if (-q2 > q1) then
                            ar(i,j,k) = THREE*a(i,j,k) - TWO*al(i,j,k)
                         end if
!                      end if
                   end if
                end if
             end if
             
             if (j==jmin) then
                if (sym < ZERO) then
                   al(i,j,k) = ZERO
                end if
             end if
             
          end do
       end do
    end do

 else

    do k=kmin,kmax-1
       do j=jmin,jmax
          do i=imin,imax

             if (k==kmin) then
                
                if ( ( (sym*d02a_m(i,j,kmin+1))*(sym*d02a_m(i,j,kmin)) .gt. 0.d0 ) .and. &
                     ( (sym*d02a_m(i,j,kmin))  * d02a_m(i,j,kmin)      .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,kmin)        * d02a_m(i,j,kmin+1)    .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,kmin+1)      * d02a_m(i,j,kmin+2)    .gt. 0.d0 ) .and.instar(i,j,k)==1) then
                   ! do nothing
                else
                   if ((ar(i,j,k)-a(i,j,k))*(a(i,j,k)-al(i,j,k)) < ZERO) then
                      al(i,j,k) = a(i,j,k)
                      ar(i,j,k) = a(i,j,k)
                   else
                      ! Otherwise, ensure monotonicity 
                      q1 = (ar(i,j,k)-al(i,j,k)) &
                           *(a(i,j,k) - HALF*(ar(i,j,k)+al(i,j,k)))
                      q2 = (ar(i,j,k)-al(i,j,k))**2/SIX
                      
                      if (q2 < q1) then
                         al(i,j,k) = THREE*a(i,j,k) - TWO*ar(i,j,k)
                      else if (-q2 > q1) then
                         ar(i,j,k) = THREE*a(i,j,k) - TWO*al(i,j,k)
                      end if
                   end if
                end if
                
             else
                if ( ( d2a_bck(i,j,k-1) * d02a_m(i,j,k-1)  .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,k-1) * d02a_m(i,j,k)    .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,k)   * d02a_m(i,j,k+1)  .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,k+1) * d2a_fwd(i,j,k+1)  .gt. 0.d0 ) .and.instar(i,j,k)==1) then
                   ! do nothing
                else
                   if ((ar(i,j,k)-a(i,j,k))*(a(i,j,k)-al(i,j,k)) < ZERO) then
                      al(i,j,k) = a(i,j,k)
                      ar(i,j,k) = a(i,j,k)
                   else
                      ! Otherwise, ensure monotonicity 
                      q1 = (ar(i,j,k)-al(i,j,k)) &
                           *(a(i,j,k) - HALF*(ar(i,j,k)+al(i,j,k)))
                      q2 = (ar(i,j,k)-al(i,j,k))**2/SIX
                      
                      if (q2 < q1) then
                         al(i,j,k) = THREE*a(i,j,k) - TWO*ar(i,j,k)
                      else if (-q2 > q1) then
                         ar(i,j,k) = THREE*a(i,j,k) - TWO*al(i,j,k)
                      end if
                   end if
                end if
             end if
             
             if (k==kmin) then
                if (sym < ZERO) then
                   al(i,j,k) = 0.d0
                end if
             end if
          end do
       end do
    end do

 end if

end subroutine monotonize

!-----------------------------------------------------------------!
! ensure monotonicity of the interpolating polynomial (Clark fix) !
! do not apply monotonization when rho > 0.85 rho_max
!-----------------------------------------------------------------!
!
subroutine monotonize_rho(ext,X,Y,Z,amax,a,ar,al,d02a_m,d2a_bck,d2a_fwd,m,sym,instar,Symmetry)
 implicit none
 integer, dimension(3)                     :: ext
 real*8, dimension(ext(1),ext(2),ext(3))   :: a,ar,al
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8, dimension(ext(1),ext(2),ext(3)), intent(in) :: d02a_m,d2a_bck,d2a_fwd
 integer, intent(in)                       :: m
 real*8                                    :: sym, amax
 integer                                   :: i,j,k
 integer                                   :: imin,imax,jmin,jmax,kmin,kmax
 real*8, parameter                         :: ZERO = 0.d0, HALF = 0.5d0, SIX = 6.d0
 real*8, parameter                         :: THREE = 3.d0, TWO = 2.d0
 real*8                                    :: q1, q2   !temporary storage variables
 integer, dimension(ext(1),ext(2),ext(3))  :: instar
 integer :: Symmetry, NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
 parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

 imax = ubound(a,1)
 imin = lbound(a,1)
 jmax = ubound(a,2)
 jmin = lbound(a,2)
 kmax = ubound(a,3)
 kmin = lbound(a,3) 

 if(Symmetry == OCTANT) then 
    if(X(1,1,1) .lt. ZERO) then
       imin = imin + 1
    end if
    if(Y(1,1,1) .lt. ZERO) then
       jmin = jmin + 1
    end if
    if(Z(1,1,1) .lt. ZERO) then
       kmin = kmin + 1
    end if
 else if(Symmetry == AXISYM) then
    if(X(1,1,1) .lt. ZERO) then
       imin = imin + 1
    end if
    if(Z(1,1,1) .lt. ZERO) then
       kmin = kmin + 1
    end if  
 else if(Symmetry == EQUATORIAL) then
    if(Z(1,1,1) .lt. ZERO) then
       kmin = kmin + 1
    end if  
 end if

 if (m==1) then

    do k=kmin,kmax
       do j=jmin,jmax
          do i=imin,imax-1

             if(abs(a(i,j,k)-amax)/amax .gt. 0.15d0) then

             if (i==imin) then
                if ( ( (sym*d02a_m(imin+1,j,k))*(sym*d02a_m(imin,j,k)) .gt. 0.d0 ) .and. &
                     ( (sym*d02a_m(imin,j,k))  * d02a_m(imin,j,k)      .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,k)        * d02a_m(i+1,j,k)    .gt. 0.d0 ) .and. &
                     ( d02a_m(i+1,j,k)      * d02a_m(i+2,j,k)    .gt. 0.d0 ) .and.instar(i,j,k)==1) then
                   ! do nothing
                else
                   if ((ar(i,j,k)-a(i,j,k))*(a(i,j,k)-al(i,j,k)) < ZERO) then
                      al(i,j,k) = a(i,j,k)
                      ar(i,j,k) = a(i,j,k)
                   else 
                      ! Otherwise, ensure monotonicity 
                      q1 = (ar(i,j,k)-al(i,j,k)) &
                           *(a(i,j,k) - HALF*(ar(i,j,k)+al(i,j,k)))
                      q2 = (ar(i,j,k)-al(i,j,k))**2/SIX
                      
                      if (q2 < q1) then
                         al(i,j,k) = THREE*a(i,j,k) - TWO*ar(i,j,k)
                      else if (-q2 > q1) then
                         ar(i,j,k) = THREE*a(i,j,k) - TWO*al(i,j,k)
                      end if
                   end if
                end if
                
             else
                if ( ( d2a_bck(i-1,j,k) * d02a_m(i-1,j,k)  .gt. 0.d0 ) .and. &
                     ( d02a_m(i-1,j,k) * d02a_m(i,j,k)    .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,k)   * d02a_m(i+1,j,k)  .gt. 0.d0 ) .and. &
                     ( d02a_m(i+1,j,k) * d2a_fwd(i+1,j,k)  .gt. 0.d0 ) .and.instar(i,j,k)==1) then
                   ! do nothing
                else
                   if ((ar(i,j,k)-a(i,j,k))*(a(i,j,k)-al(i,j,k)) < ZERO) then
                      al(i,j,k) = a(i,j,k)
                      ar(i,j,k) = a(i,j,k)
                   else 
                      ! Otherwise, ensure monotonicity 
                      q1 = (ar(i,j,k)-al(i,j,k)) &
                           *(a(i,j,k) - HALF*(ar(i,j,k)+al(i,j,k)))
                      q2 = (ar(i,j,k)-al(i,j,k))**2/SIX

                      if (q2 < q1) then
                         al(i,j,k) = THREE*a(i,j,k) - TWO*ar(i,j,k)
                      else if (-q2 > q1) then
                         ar(i,j,k) = THREE*a(i,j,k) - TWO*al(i,j,k)
                      end if
                   end if
                end if
             end if
             
             if (i==imin) then
                if (sym < ZERO) then
                   al(i,j,k) = ZERO
                end if
             end if
          end if
             
          end do
       end do
    end do

 else if (m==2) then

    do k=kmin,kmax
       do j=jmin,jmax-1
          do i=imin,imax

             if(abs(a(i,j,k)-amax)/amax .gt. 0.15d0) then

             if (j==jmin) then
                if ( ( (sym*d02a_m(i,jmin+1,k))*(sym*d02a_m(i,jmin,k)) .gt. 0.d0 ) .and. &
                     ( (sym*d02a_m(i,jmin,k))  * d02a_m(i,jmin,k)      .gt. 0.d0 ) .and. &
                     ( d02a_m(i,jmin,k)        * d02a_m(i,jmin+1,k)    .gt. 0.d0 ) .and. &
                     ( d02a_m(i,jmin+1,k)      * d02a_m(i,jmin+2,k)    .gt. 0.d0 ) .and.instar(i,j,k)==1) then
                   ! do nothing
                else
                   if ((ar(i,j,k)-a(i,j,k))*(a(i,j,k)-al(i,j,k)) < ZERO) then
                      al(i,j,k) = a(i,j,k)
                      ar(i,j,k) = a(i,j,k)
                   else 
                      ! Otherwise, ensure monotonicity 
                      q1 = (ar(i,j,k)-al(i,j,k))*(a(i,j,k) - HALF*(ar(i,j,k)+al(i,j,k)))
                      q2 = (ar(i,j,k)-al(i,j,k))**2/SIX
                      
                      if (q2 < q1) then
                         al(i,j,k) = THREE*a(i,j,k) - TWO*ar(i,j,k)
                      else if (-q2 > q1) then
                         ar(i,j,k) = THREE*a(i,j,k) - TWO*al(i,j,k)
                      end if
                   end if
                end if                
             else
                if ( ( d2a_bck(i,j-1,k) * d02a_m(i,j-1,k)  .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j-1,k) * d02a_m(i,j,k)    .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,k)   * d02a_m(i,j+1,k)  .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j+1,k) * d2a_fwd(i,j+1,k)  .gt. 0.d0 ) .and.instar(i,j,k)==1) then
                   ! do nothing
                else
                   if ((ar(i,j,k)-a(i,j,k))*(a(i,j,k)-al(i,j,k)) < ZERO) then
                      al(i,j,k) = a(i,j,k)
                      ar(i,j,k) = a(i,j,k)
                   else 
                      ! Otherwise, ensure monotonicity 
                      q1 = (ar(i,j,k)-al(i,j,k))*(a(i,j,k) - HALF*(ar(i,j,k)+al(i,j,k)))
                      q2 = (ar(i,j,k)-al(i,j,k))**2/SIX
                      
!                      if ((j /= 32) .and. (j /= 33)) then
                         if (q2 < q1) then
                            al(i,j,k) = THREE*a(i,j,k) - TWO*ar(i,j,k)
                         else if (-q2 > q1) then
                            ar(i,j,k) = THREE*a(i,j,k) - TWO*al(i,j,k)
                         end if
!                      end if
                   end if
                end if
             end if
             
             if (j==jmin) then
                if (sym < ZERO) then
                   al(i,j,k) = ZERO
                end if
             end if
             
          end if
          end do
       end do
    end do

 else

    do k=kmin,kmax-1
       do j=jmin,jmax
          do i=imin,imax

             if(abs(a(i,j,k)-amax)/amax .gt. 0.15d0) then

             if (k==kmin) then
                
                if ( ( (sym*d02a_m(i,j,kmin+1))*(sym*d02a_m(i,j,kmin)) .gt. 0.d0 ) .and. &
                     ( (sym*d02a_m(i,j,kmin))  * d02a_m(i,j,kmin)      .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,kmin)        * d02a_m(i,j,kmin+1)    .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,kmin+1)      * d02a_m(i,j,kmin+2)    .gt. 0.d0 ) .and.instar(i,j,k)==1) then
                   ! do nothing
                else
                   if ((ar(i,j,k)-a(i,j,k))*(a(i,j,k)-al(i,j,k)) < ZERO) then
                      al(i,j,k) = a(i,j,k)
                      ar(i,j,k) = a(i,j,k)
                   else
                      ! Otherwise, ensure monotonicity 
                      q1 = (ar(i,j,k)-al(i,j,k)) &
                           *(a(i,j,k) - HALF*(ar(i,j,k)+al(i,j,k)))
                      q2 = (ar(i,j,k)-al(i,j,k))**2/SIX
                      
                      if (q2 < q1) then
                         al(i,j,k) = THREE*a(i,j,k) - TWO*ar(i,j,k)
                      else if (-q2 > q1) then
                         ar(i,j,k) = THREE*a(i,j,k) - TWO*al(i,j,k)
                      end if
                   end if
                end if
                
             else
                if ( ( d2a_bck(i,j,k-1) * d02a_m(i,j,k-1)  .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,k-1) * d02a_m(i,j,k)    .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,k)   * d02a_m(i,j,k+1)  .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,k+1) * d2a_fwd(i,j,k+1)  .gt. 0.d0 ) .and.instar(i,j,k)==1) then
                   ! do nothing
                else
                   if ((ar(i,j,k)-a(i,j,k))*(a(i,j,k)-al(i,j,k)) < ZERO) then
                      al(i,j,k) = a(i,j,k)
                      ar(i,j,k) = a(i,j,k)
                   else
                      ! Otherwise, ensure monotonicity 
                      q1 = (ar(i,j,k)-al(i,j,k)) &
                           *(a(i,j,k) - HALF*(ar(i,j,k)+al(i,j,k)))
                      q2 = (ar(i,j,k)-al(i,j,k))**2/SIX
                      
                      if (q2 < q1) then
                         al(i,j,k) = THREE*a(i,j,k) - TWO*ar(i,j,k)
                      else if (-q2 > q1) then
                         ar(i,j,k) = THREE*a(i,j,k) - TWO*al(i,j,k)
                      end if
                   end if
                end if
             end if
             
             if (k==kmin) then
                if (sym < ZERO) then
                   al(i,j,k) = 0.d0
                end if
             end if

          end if
          end do
       end do
    end do

 end if

end subroutine monotonize_rho
 
!-----------------------------------------------------------------!
! ensure monotonicity of the interpolating polynomial (Clark fix) !
! do not apply monotonization when rho > 0.85 rho_max
!-----------------------------------------------------------------!
!
subroutine monotonize_P(ext,X,Y,Z,rhomax,rho,a,ar,al,d02a_m,d2a_bck,d2a_fwd,m,sym,instar,Symmetry)
 implicit none
 integer, dimension(3)                     :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8, dimension(ext(1),ext(2),ext(3))   :: a,ar,al,rho
 real*8, dimension(ext(1),ext(2),ext(3)), intent(in) :: d02a_m,d2a_bck,d2a_fwd
 integer, intent(in)                       :: m
 real*8                                    :: sym, rhomax
 integer                                   :: i,j,k
 integer                                   :: imin,imax,jmin,jmax,kmin,kmax
 real*8, parameter                         :: ZERO = 0.d0, HALF = 0.5d0, SIX = 6.d0
 real*8, parameter                         :: THREE = 3.d0, TWO = 2.d0
 real*8                                    :: q1, q2   !temporary storage variables
 integer, dimension(ext(1),ext(2),ext(3))  :: instar
 integer :: Symmetry, NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
 parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

 imax = ubound(a,1)
 imin = lbound(a,1)
 jmax = ubound(a,2)
 jmin = lbound(a,2)
 kmax = ubound(a,3)
 kmin = lbound(a,3) 

 if(Symmetry == OCTANT) then 
    if(X(1,1,1) .lt. ZERO) then
       imin = imin + 1
    end if
    if(Y(1,1,1) .lt. ZERO) then
       jmin = jmin + 1
    end if
    if(Z(1,1,1) .lt. ZERO) then
       kmin = kmin + 1
    end if
 else if(Symmetry == AXISYM) then
    if(X(1,1,1) .lt. ZERO) then
       imin = imin + 1
    end if
    if(Z(1,1,1) .lt. ZERO) then
       kmin = kmin + 1
    end if  
 else if(Symmetry == EQUATORIAL) then
    if(Z(1,1,1) .lt. ZERO) then
       kmin = kmin + 1
    end if  
 end if

 if (m==1) then

    do k=kmin,kmax
       do j=jmin,jmax
          do i=imin,imax-1

             if(abs(rho(i,j,k)-rhomax)/rhomax .gt. 0.15d0) then

             if (i==imin) then
                if ( ( (sym*d02a_m(imin+1,j,k))*(sym*d02a_m(imin,j,k)) .gt. 0.d0 ) .and. &
                     ( (sym*d02a_m(imin,j,k))  * d02a_m(imin,j,k)      .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,k)        * d02a_m(i+1,j,k)    .gt. 0.d0 ) .and. &
                     ( d02a_m(i+1,j,k)      * d02a_m(i+2,j,k)    .gt. 0.d0 ) .and.instar(i,j,k)==1) then
                   ! do nothing
                else
                   if ((ar(i,j,k)-a(i,j,k))*(a(i,j,k)-al(i,j,k)) < ZERO) then
                      al(i,j,k) = a(i,j,k)
                      ar(i,j,k) = a(i,j,k)
                      
                   else 
                      ! Otherwise, ensure monotonicity 
                      q1 = (ar(i,j,k)-al(i,j,k)) &
                           *(a(i,j,k) - HALF*(ar(i,j,k)+al(i,j,k)))
                      q2 = (ar(i,j,k)-al(i,j,k))**2/SIX
                      
                      if (q2 < q1) then
                         al(i,j,k) = THREE*a(i,j,k) - TWO*ar(i,j,k)
                      else if (-q2 > q1) then
                         ar(i,j,k) = THREE*a(i,j,k) - TWO*al(i,j,k)
                      end if
                   end if
                end if
                
             else
                if ( ( d2a_bck(i-1,j,k) * d02a_m(i-1,j,k)  .gt. 0.d0 ) .and. &
                     ( d02a_m(i-1,j,k) * d02a_m(i,j,k)    .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,k)   * d02a_m(i+1,j,k)  .gt. 0.d0 ) .and. &
                     ( d02a_m(i+1,j,k) * d2a_fwd(i+1,j,k)  .gt. 0.d0 ) .and.instar(i,j,k)==1) then
                   ! do nothing
                else
                   if ((ar(i,j,k)-a(i,j,k))*(a(i,j,k)-al(i,j,k)) < ZERO) then
                      al(i,j,k) = a(i,j,k)
                      ar(i,j,k) = a(i,j,k)
                   else 
                      ! Otherwise, ensure monotonicity 
                      q1 = (ar(i,j,k)-al(i,j,k)) &
                           *(a(i,j,k) - HALF*(ar(i,j,k)+al(i,j,k)))
                      q2 = (ar(i,j,k)-al(i,j,k))**2/SIX

                      if (q2 < q1) then
                         al(i,j,k) = THREE*a(i,j,k) - TWO*ar(i,j,k)
                      else if (-q2 > q1) then
                         ar(i,j,k) = THREE*a(i,j,k) - TWO*al(i,j,k)
                      end if
                   end if
                end if
             end if
             
             if (i==imin) then
                if (sym < ZERO) then
                   al(i,j,k) = ZERO
                end if
             end if
          end if
             
          end do
       end do
    end do

 else if (m==2) then

    do k=kmin,kmax
       do j=jmin,jmax-1
          do i=imin,imax

             if(abs(rho(i,j,k)-rhomax)/rhomax .gt. 0.15d0) then

             if (j==jmin) then
                if ( ( (sym*d02a_m(i,jmin+1,k))*(sym*d02a_m(i,jmin,k)) .gt. 0.d0 ) .and. &
                     ( (sym*d02a_m(i,jmin,k))  * d02a_m(i,jmin,k)      .gt. 0.d0 ) .and. &
                     ( d02a_m(i,jmin,k)        * d02a_m(i,jmin+1,k)    .gt. 0.d0 ) .and. &
                     ( d02a_m(i,jmin+1,k)      * d02a_m(i,jmin+2,k)    .gt. 0.d0 ) .and.instar(i,j,k)==1) then
                   ! do nothing
                else
                   if ((ar(i,j,k)-a(i,j,k))*(a(i,j,k)-al(i,j,k)) < ZERO) then
                      al(i,j,k) = a(i,j,k)
                      ar(i,j,k) = a(i,j,k)
                   else 
                      ! Otherwise, ensure monotonicity 
                      q1 = (ar(i,j,k)-al(i,j,k))*(a(i,j,k) - HALF*(ar(i,j,k)+al(i,j,k)))
                      q2 = (ar(i,j,k)-al(i,j,k))**2/SIX
                      
                      if (q2 < q1) then
                         al(i,j,k) = THREE*a(i,j,k) - TWO*ar(i,j,k)
                      else if (-q2 > q1) then
                         ar(i,j,k) = THREE*a(i,j,k) - TWO*al(i,j,k)
                      end if
                   end if
                end if                
             else
                if ( ( d2a_bck(i,j-1,k) * d02a_m(i,j-1,k)  .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j-1,k) * d02a_m(i,j,k)    .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,k)   * d02a_m(i,j+1,k)  .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j+1,k) * d2a_fwd(i,j+1,k)  .gt. 0.d0 ) .and.instar(i,j,k)==1) then
                   ! do nothing
                else
                   if ((ar(i,j,k)-a(i,j,k))*(a(i,j,k)-al(i,j,k)) < ZERO) then
                      al(i,j,k) = a(i,j,k)
                      ar(i,j,k) = a(i,j,k)
                   else 
                      ! Otherwise, ensure monotonicity 
                      q1 = (ar(i,j,k)-al(i,j,k))*(a(i,j,k) - HALF*(ar(i,j,k)+al(i,j,k)))
                      q2 = (ar(i,j,k)-al(i,j,k))**2/SIX
                      
!                      if ((j /= 32) .and. (j /= 33)) then
                         if (q2 < q1) then
                            al(i,j,k) = THREE*a(i,j,k) - TWO*ar(i,j,k)
                         else if (-q2 > q1) then
                            ar(i,j,k) = THREE*a(i,j,k) - TWO*al(i,j,k)
                         end if
!                      end if
                   end if
                end if
             end if
             
             if (j==jmin) then
                if (sym < ZERO) then
                   al(i,j,k) = ZERO
                end if
             end if
             
          end if
          end do
       end do
    end do

 else

    do k=kmin,kmax-1
       do j=jmin,jmax
          do i=imin,imax

             if(abs(rho(i,j,k)-rhomax)/rhomax .gt. 0.15d0) then

             if (k==kmin) then
                
                if ( ( (sym*d02a_m(i,j,kmin+1))*(sym*d02a_m(i,j,kmin)) .gt. 0.d0 ) .and. &
                     ( (sym*d02a_m(i,j,kmin))  * d02a_m(i,j,kmin)      .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,kmin)        * d02a_m(i,j,kmin+1)    .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,kmin+1)      * d02a_m(i,j,kmin+2)    .gt. 0.d0 ) .and.instar(i,j,k)==1) then
                   ! do nothing
                else
                   if ((ar(i,j,k)-a(i,j,k))*(a(i,j,k)-al(i,j,k)) < ZERO) then
                      al(i,j,k) = a(i,j,k)
                      ar(i,j,k) = a(i,j,k)
                   else
                      ! Otherwise, ensure monotonicity 
                      q1 = (ar(i,j,k)-al(i,j,k)) &
                           *(a(i,j,k) - HALF*(ar(i,j,k)+al(i,j,k)))
                      q2 = (ar(i,j,k)-al(i,j,k))**2/SIX
                      
                      if (q2 < q1) then
                         al(i,j,k) = THREE*a(i,j,k) - TWO*ar(i,j,k)
                      else if (-q2 > q1) then
                         ar(i,j,k) = THREE*a(i,j,k) - TWO*al(i,j,k)
                      end if
                   end if
                end if
                
             else
                if ( ( d2a_bck(i,j,k-1) * d02a_m(i,j,k-1)  .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,k-1) * d02a_m(i,j,k)    .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,k)   * d02a_m(i,j,k+1)  .gt. 0.d0 ) .and. &
                     ( d02a_m(i,j,k+1) * d2a_fwd(i,j,k+1)  .gt. 0.d0 ) .and.instar(i,j,k)==1) then
                   ! do nothing
                else
                   if ((ar(i,j,k)-a(i,j,k))*(a(i,j,k)-al(i,j,k)) < ZERO) then
                      al(i,j,k) = a(i,j,k)
                      ar(i,j,k) = a(i,j,k)
                   else
                      ! Otherwise, ensure monotonicity 
                      q1 = (ar(i,j,k)-al(i,j,k)) &
                           *(a(i,j,k) - HALF*(ar(i,j,k)+al(i,j,k)))
                      q2 = (ar(i,j,k)-al(i,j,k))**2/SIX
                      
                      if (q2 < q1) then
                         al(i,j,k) = THREE*a(i,j,k) - TWO*ar(i,j,k)
                      else if (-q2 > q1) then
                         ar(i,j,k) = THREE*a(i,j,k) - TWO*al(i,j,k)
                      end if
                   end if
                end if
             end if
             
             if (k==kmin) then
                if (sym < ZERO) then
                   al(i,j,k) = 0.d0
                end if
             end if

          end if
          end do
       end do
    end do

 end if

end subroutine monotonize_P

!-----------------------------------------------------------------!
! calculate \tilde{f}_j                                           !
!-----------------------------------------------------------------!
!
subroutine ppm_ftilde(ext,X,Y,Z,ftilde,P,d0P_m,vx,vy,vz,P_max,m,Symmetry)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8, dimension(ext(1),ext(2),ext(3))        :: ftilde,P,d0P_m
 real*8, dimension(ext(1),ext(2),ext(3))        :: vx,vy,vz
 integer, intent(in)                            :: m,Symmetry
 integer                                        :: i,j,k
 integer                                        :: imin,imax,jmin,jmax,kmin,kmax
 real*8                             :: ONE, TWO, SYM, ANTI, ZERO, P_max
 parameter(ONE = 1.D0, TWO = 2.D0, SYM = 1.D0, ANTI = -1.D0, ZERO = 0.d0)
 integer                            :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM
 integer                            :: AXISYM
 parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
 real*8, parameter                  :: omega1 = 0.75d0, omega2 = 10.d0
 real*8, parameter                  :: epsilon2 = 0.33d0
! The following parameter choices are Font-like.
! real*8, parameter                  :: omega1 = 0.52d0, omega2 = 10.d0
! real*8, parameter                  :: epsilon2 = 0.5d0
 real*8                             :: w, q1, TINY

 imax = ubound(ftilde,1)
 imin = lbound(ftilde,1)
 jmax = ubound(ftilde,2)
 jmin = lbound(ftilde,2)
 kmax = ubound(ftilde,3)
 kmin = lbound(ftilde,3) 
 TINY = 1.d-15*P_max

 if(Symmetry == OCTANT) then 
    if(X(1,1,1).lt.ZERO) then
       imin = imin + 1
    end if
    if(Y(1,1,1).lt.ZERO) then
       jmin = jmin + 1
    end if
    if(Z(1,1,1).lt.ZERO) then
       kmin = kmin + 1
    end if
 else if(Symmetry == AXISYM) then
    if(X(1,1,1).lt.ZERO) then
       imin = imin + 1
    end if
    if(Z(1,1,1).lt.ZERO) then
       kmin = kmin + 1
    end if
 else if(Symmetry == EQUATORIAL) then
    if(Z(1,1,1).lt.ZERO) then
       kmin = kmin + 1
    end if
 end if

 if (m==1) then

    do k=kmin,kmax
       do j=jmin,jmax
          do i=imin+1,imax-1

             if (min(P(i+1,j,k),P(i-1,j,k)) .gt. TINY) then
                if (abs(P(i+1,j,k) - P(i-1,j,k))/min(P(i+1,j,k),P(i-1,j,k)) > epsilon2 .and. &
                     ((P(i+1,j,k)-P(i-1,j,k))*(vx(i-1,j,k)-vx(i+1,j,k)) > ZERO) ) then
                   w = ONE
                else
                   w = ZERO
                end if
                
                if ((d0P_m(i+1,j,k)+d0P_m(i-1,j,k)).gt.TINY) then
                   q1 = omega2*(d0P_m(i,j,k)/(d0P_m(i+1,j,k)+d0P_m(i-1,j,k)) - omega1)
                   ftilde(i,j,k) = min(ONE,w*max(ZERO,q1))
                else
                   ftilde(i,j,k) = ZERO
                end if
             else
                ftilde(i,j,k) = ZERO
             end if
          end do
       end do
    end do
    
    ftilde(imin,:,:) = ZERO
    ftilde(imax,:,:) = ZERO

 else if (m==2) then

    do k=kmin,kmax
       do j=jmin+1,jmax-1
          do i=imin,imax

             if (min(P(i,j+1,k),P(i,j-1,k)) .gt. TINY) then
                if (abs(P(i,j+1,k) - P(i,j-1,k))/min(P(i,j+1,k),P(i,j-1,k)) > epsilon2 .and. &
                     ((P(i,j+1,k)-P(i,j-1,k))*(vy(i,j-1,k)-vy(i,j+1,k)) > ZERO) ) then
                   w = ONE
                else
                   w = ZERO
                end if
                
                if ((d0P_m(i,j+1,k)+d0P_m(i,j-1,k)) .gt. TINY) then
                   q1 = omega2*(d0P_m(i,j,k)/(d0P_m(i,j+1,k)+d0P_m(i,j-1,k)) - omega1)
                   ftilde(i,j,k) = min(ONE,w*max(ZERO,q1))
                else
                   ftilde(i,j,k) = ZERO
                end if
             else
                ftilde(i,j,k) = ZERO
             end if
          end do
       end do
    end do
    
    ftilde(:,jmin,:) = ZERO
    ftilde(:,jmax,:) = ZERO
   
 else

    do k=kmin+1,kmax-1
       do j=jmin,jmax
          do i=imin,imax

             if (min(P(i,j,k+1),P(i,j,k-1)) .gt. TINY) then
                if (abs(P(i,j,k+1) - P(i,j,k-1))/min(P(i,j,k+1),P(i,j,k-1)) > epsilon2 .and. &
                     ((P(i,j,k+1)-P(i,j,k-1))*(vz(i,j,k-1)-vz(i,j,k+1)) > ZERO) ) then
                   w = ONE
                else
                   w = ZERO
                end if

                if ((d0P_m(i,j,k+1)+d0P_m(i,j,k-1)) .gt. TINY) then
                   q1 = omega2*(d0P_m(i,j,k)/(d0P_m(i,j,k+1)+d0P_m(i,j,k-1)) - omega1)
                   ftilde(i,j,k) = min(ONE,w*max(ZERO,q1))
                else
                   ftilde(i,j,k) = ZERO
                end if

             else
                ftilde(i,j,k) = ZERO
             end if
          end do
       end do
    end do
    
    ftilde(:,:,kmin) = ZERO
    ftilde(:,:,kmax) = ZERO

 end if

end subroutine ppm_ftilde

!-----------------------------------------------------------------!
! apply flattening                                                !
!-----------------------------------------------------------------!
!
subroutine ppm_flatten(ext,X,Y,Z,ftilde,rho_b,rho_br,rho_bl,P,Pr,Pl, &
     vx,vxr,vxl,vy,vyr,vyl,vz,vzr,vzl,m,Symmetry)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8, dimension(ext(1),ext(2),ext(3))        :: rho_b, rho_br, rho_bl
 real*8, dimension(ext(1),ext(2),ext(3))        :: ftilde,P,Pr,Pl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vx,vxr,vxl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vy,vyr,vyl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vz,vzr,vzl
 integer, intent(in)                            :: m,Symmetry
 real*8                             :: ONE, TWO, SYM, ANTI, ZERO, f
 parameter(ONE = 1.D0, TWO = 2.D0, SYM = 1.D0, ANTI = -1.D0, ZERO = 0.d0)
 integer                            :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM
 integer                            :: AXISYM
 parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
 integer                            :: i,j,k, imin,imax,jmin,jmax,kmin,kmax
 integer                            :: s

 imax = ubound(ftilde,1)
 imin = lbound(ftilde,1)
 jmax = ubound(ftilde,2)
 jmin = lbound(ftilde,2)
 kmax = ubound(ftilde,3)
 kmin = lbound(ftilde,3) 

 if(Symmetry == OCTANT) then 
    if(X(1,1,1) .lt. ZERO) then
       imin = imin + 1
    end if
    if(Y(1,1,1) .lt. ZERO) then
       jmin = jmin + 1
    end if
    if(Z(1,1,1) .lt. ZERO) then
       kmin = kmin + 1
    end if
 else if(Symmetry == AXISYM) then
    if(X(1,1,1) .lt. ZERO) then
       imin = imin + 1
    end if
    if(Z(1,1,1) .lt. ZERO) then
       kmin = kmin + 1
    end if
 else if(Symmetry == EQUATORIAL) then
    if(Z(1,1,1) .lt. ZERO) then
       kmin = kmin + 1
    end if
 end if


 if (m==1) then

    do k=kmin,kmax
       do j=jmin,jmax
          do i=imin+1,imax-1
             ! Calculate s
             if ((P(i+1,j,k) - P(i-1,j,k)) <= ZERO) then
                s = 1
             else
                s = -1
             end if

             ! Calculate f
             f = max(ftilde(i,j,k),ftilde(i+s,j,k))

             ! Flatten
             rho_br(i,j,k)   = rho_b(i,j,k)*f + rho_br(i,j,k)*(ONE-f)
             rho_bl(i+1,j,k) = rho_b(i,j,k)*f + rho_bl(i+1,j,k)*(ONE-f)
             Pr(i,j,k)       = P(i,j,k)*f     + Pr(i,j,k)*(ONE-f)
             Pl(i+1,j,k)     = P(i,j,k)*f     + Pl(i+1,j,k)*(ONE-f)
             vxr(i,j,k)      = vx(i,j,k)*f    + vxr(i,j,k)*(ONE-f)
             vxl(i+1,j,k)    = vx(i,j,k)*f    + vxl(i+1,j,k)*(ONE-f)
             vyr(i,j,k)      = vy(i,j,k)*f    + vyr(i,j,k)*(ONE-f)
             vyl(i+1,j,k)    = vy(i,j,k)*f    + vyl(i+1,j,k)*(ONE-f)
             vzr(i,j,k)      = vz(i,j,k)*f    + vzr(i,j,k)*(ONE-f)
             vzl(i+1,j,k)    = vz(i,j,k)*f    + vzl(i+1,j,k)*(ONE-f)

          end do
       end do
    end do

 else if (m==2) then

    do k=kmin,kmax
       do j=jmin+1,jmax-1
          do i=imin,imax
             ! Calculate s
             if ((P(i,j+1,k) - P(i,j-1,k)) <= ZERO) then
                s = 1
             else
                s = -1
             end if

             ! Calculate f
             f = max(ftilde(i,j,k),ftilde(i,j+s,k))

             ! Flatten
             rho_br(i,j,k)   = rho_b(i,j,k)*f + rho_br(i,j,k)*(ONE-f)
             rho_bl(i,j+1,k) = rho_b(i,j,k)*f + rho_bl(i,j+1,k)*(ONE-f)
             Pr(i,j,k)       = P(i,j,k)*f     + Pr(i,j,k)*(ONE-f)
             Pl(i,j+1,k)     = P(i,j,k)*f     + Pl(i,j+1,k)*(ONE-f)
             vxr(i,j,k)      = vx(i,j,k)*f    + vxr(i,j,k)*(ONE-f)
             vxl(i,j+1,k)    = vx(i,j,k)*f    + vxl(i,j+1,k)*(ONE-f)
             vyr(i,j,k)      = vy(i,j,k)*f    + vyr(i,j,k)*(ONE-f)
             vyl(i,j+1,k)    = vy(i,j,k)*f    + vyl(i,j+1,k)*(ONE-f)
             vzr(i,j,k)      = vz(i,j,k)*f    + vzr(i,j,k)*(ONE-f)
             vzl(i,j+1,k)    = vz(i,j,k)*f    + vzl(i,j+1,k)*(ONE-f)

          end do
       end do
    end do

 else

    do k=kmin+1,kmax-1
       do j=jmin,jmax
          do i=imin,imax
             ! Calculate s
             if ((P(i,j,k+1) - P(i,j,k-1)) <= ZERO) then
                s = 1
             else
                s = -1
             end if

             ! Calculate f
             f = max(ftilde(i,j,k),ftilde(i,j,k+s))

             ! Flatten
             rho_br(i,j,k)   = rho_b(i,j,k)*f + rho_br(i,j,k)*(ONE-f)
             rho_bl(i,j,k+1) = rho_b(i,j,k)*f + rho_bl(i,j,k+1)*(ONE-f)
             Pr(i,j,k)       = P(i,j,k)*f     + Pr(i,j,k)*(ONE-f)
             Pl(i,j,k+1)     = P(i,j,k)*f     + Pl(i,j,k+1)*(ONE-f)
             vxr(i,j,k)      = vx(i,j,k)*f    + vxr(i,j,k)*(ONE-f)
             vxl(i,j,k+1)    = vx(i,j,k)*f    + vxl(i,j,k+1)*(ONE-f)
             vyr(i,j,k)      = vy(i,j,k)*f    + vyr(i,j,k)*(ONE-f)
             vyl(i,j,k+1)    = vy(i,j,k)*f    + vyl(i,j,k+1)*(ONE-f)
             vzr(i,j,k)      = vz(i,j,k)*f    + vzr(i,j,k)*(ONE-f)
             vzl(i,j,k+1)    = vz(i,j,k)*f    + vzl(i,j,k+1)*(ONE-f)

          end do
       end do
    end do

 end if

end subroutine ppm_flatten

!-----------------------------------------------------------------!
!  find f_b(i) = f_c(i-1),  f_f(i) = f_c(i+1)                     !
!-----------------------------------------------------------------!
subroutine hydro_forwardback(ext,X,Y,Z,d2rho_c,d2rho_b,d2rho_f,d2P_c,d2P_b,d2P_f, &
     d2vx_c,d2vx_b,d2vx_f,d2vy_c,d2vy_b,d2vy_f,d2vz_c,d2vz_b,d2vz_f, &
     m,Symmetry)
 implicit none
 interface
    subroutine forwardback(ext,X,Y,Z,m,f_c,f_b,f_f,sym,Symmetry)
      implicit none
      integer, dimension(3)                    :: ext
      real*8, dimension(ext(1),ext(2),ext(3))  :: X,Y,Z
      real*8, dimension(ext(1),ext(2),ext(3))  :: f_c,f_b,f_f
      real*8                                   :: sym
      integer                                  :: m,Symmetry
    end subroutine forwardback
 end interface
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8, dimension(ext(1),ext(2),ext(3))        :: d2rho_c,d2rho_b,d2rho_f
 real*8, dimension(ext(1),ext(2),ext(3))        :: d2P_c,d2P_b,d2P_f
 real*8, dimension(ext(1),ext(2),ext(3))        :: d2vx_c,d2vx_b,d2vx_f
 real*8, dimension(ext(1),ext(2),ext(3))        :: d2vy_c,d2vy_b,d2vy_f
 real*8, dimension(ext(1),ext(2),ext(3))        :: d2vz_c,d2vz_b,d2vz_f
 integer                                        :: m, Symmetry
 real*8                     :: ONE, TWO, SYM, ANTI
 parameter(ONE = 1.0D0, TWO = 2.D0, SYM = 1.D0, ANTI = -1.D0)
 integer                            :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
 parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)


 call forwardback(ext,X,Y,Z,m,d2rho_c,d2rho_b,d2rho_f,SYM,Symmetry)
 call forwardback(ext,X,Y,Z,m,d2P_c,  d2P_b,  d2P_f,  SYM,Symmetry)
 if(m==1) then
    call forwardback(ext,X,Y,Z,m,d2vz_c, d2vz_b, d2vz_f, SYM,Symmetry)
    if(Symmetry==2) then
       call forwardback(ext,X,Y,Z,m,d2vx_c,d2vx_b,d2vx_f,ANTI,Symmetry)
       call forwardback(ext,X,Y,Z,m,d2vy_c,d2vy_b,d2vy_f,SYM,Symmetry)
    else if(Symmetry==4) then
       call forwardback(ext,X,Y,Z,m,d2vx_c,d2vx_b,d2vx_f,ANTI,Symmetry)
       call forwardback(ext,X,Y,Z,m,d2vy_c,d2vy_b,d2vy_f,ANTI,Symmetry)
    else
       call forwardback(ext,X,Y,Z,m,d2vx_c,d2vx_b,d2vx_f,SYM,Symmetry)
       call forwardback(ext,X,Y,Z,m,d2vy_c,d2vy_b,d2vy_f,SYM,Symmetry)
    end if

 else if(m==2) then
    if(Symmetry.ne.4) then
       call forwardback(ext,X,Y,Z,m,d2vz_c, d2vz_b, d2vz_f, SYM,Symmetry)
       call forwardback(ext,X,Y,Z,m,d2vx_c, d2vx_b, d2vx_f, SYM,Symmetry)
       if(Symmetry==2) then
          call forwardback(ext,X,Y,Z,m,d2vy_c,d2vy_b,d2vy_f,ANTI,Symmetry)
       else
          call forwardback(ext,X,Y,Z,m,d2vy_c,d2vy_b,d2vy_f,SYM,Symmetry)
       end if
    end if
 else if(m==3) then
    call forwardback(ext,X,Y,Z,m,d2vx_c, d2vx_b, d2vx_f, SYM,Symmetry)
    call forwardback(ext,X,Y,Z,m,d2vy_c, d2vy_b, d2vy_f, SYM,Symmetry)
    if(Symmetry.ne.0) then
       call forwardback(ext,X,Y,Z,m,d2vz_c,d2vz_b,d2vz_f,ANTI,Symmetry)
    else
       call forwardback(ext,X,Y,Z,m,d2vz_c,d2vz_b,d2vz_f,SYM,Symmetry)
    end if
 else
    write(*,*) "m not in proper range"
 end if

! write(*,*) "inside hydro_fb",m,d2vx_c(2,2,2),d2vx_b(2,2,2),d2vx_f(2,2,2),d2vx_c(2,2,4),d2vz_c(2,2,2),d2vz_b(2,2,2),d2vz_f(2,2,2),d2vz_c(2,2,2),d2vz_b(2,2,2),d2vz_f(2,2,2)

end subroutine hydro_forwardback

subroutine forwardback(ext,X,Y,Z,m,f_c,f_b,f_f,sym)
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: X,Y,Z
  real*8, dimension(ext(1),ext(2),ext(3))  :: f_c,f_b,f_f
  real*8                                   :: sym
  integer                                  :: m
  integer                                  :: k,kmin,kmax
  integer                            :: Symmetry, NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  kmin = lbound(f_c,m)
  kmax = ubound(f_c,m)

  if(Symmetry == OCTANT) then 
     if(m==1 .and. X(1,1,1).lt.0.D0) then
        kmin = lbound(f_c,1)+1
     else if(m==2 .and. Y(1,1,1).lt.0.D0) then
        kmin = lbound(f_c,2)+1
     else if(m==3 .and. Z(1,1,1).lt.0.D0) then
        kmin = lbound(f_c,3)+1
     end if
  else if(Symmetry == AXISYM) then 
     if(m==1 .and. X(1,1,1).lt.0.D0) then
        kmin = lbound(f_c,1)+1
     else if(m==3 .and. Z(1,1,1).lt.0.D0) then
        kmin = lbound(f_c,3)+1
     end if
  end if

  if(m==1) then
     f_b(kmin,:,:) = sym*f_c(kmin,:,:)
     f_f(kmin,:,:) =     f_c(kmin+1,:,:)
     do k = kmin+1, kmax-1
        f_b(k,:,:) = f_c(k-1,:,:)
        f_f(k,:,:) = f_c(k+1,:,:)
     end do
     f_b(kmax,:,:) = f_c(kmax-1,:,:)
     f_f(kmax,:,:) = f_c(kmax,:,:)

  else if(m==2) then
     f_b(:,kmin,:) = sym*f_c(:,kmin,:)
     f_f(:,kmin,:) =     f_c(:,kmin+1,:)
     do k = kmin+1, kmax-1
        f_b(:,k,:) = f_c(:,k-1,:)
        f_f(:,k,:) = f_c(:,k+1,:)
     end do
     f_b(:,kmax,:) = f_c(:,kmax-1,:)
     f_f(:,kmax,:) = f_c(:,kmax,:)
     
  else if(m==3) then
     f_b(:,:,kmin) = sym*f_c(:,:,kmin)
     f_f(:,:,kmin) =     f_c(:,:,kmin+1)
     do k = kmin+1, kmax-1
        f_b(:,:,k) = f_c(:,:,k-1)
        f_f(:,:,k) = f_c(:,:,k+1)
     end do
     f_b(:,:,kmax) = f_c(:,:,kmax-1)
     f_f(:,:,kmax) = f_c(:,:,kmax)
     
  else
     write(*,*)  "m not in proper range"
  end if

end subroutine forwardback

subroutine ppm_shift(ext,X,Y,Z,rhor,rhol,Pr,Pl,vxr,vxl,vyr,vyl,vzr,vzl,m,Symmetry)
  implicit none
   integer, dimension(3)                        :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8, dimension(ext(1),ext(2),ext(3))        :: rhor,rhol,Pr,Pl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vxr,vxl,vyr,vyl,vzr,vzl
 real*8, dimension(ext(1),ext(2),ext(3))        :: rhot,Pt,vxt,vyt,vzt
 integer                                        :: m, Symmetry
 integer                                        :: k, kmin, kmax
  integer                            :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

 kmin = lbound(rhor,m)
 kmax = ubound(rhor,m)
 
 if(Symmetry == OCTANT) then 
    if(m==1 .and. X(1,1,1).lt.0.D0) then
       kmin = lbound(rhor,1)+1
    else if(m==2 .and. Y(1,1,1).lt.0.D0) then
       kmin = lbound(rhor,2)+1
    else if(m==3 .and. Z(1,1,1).lt.0.D0) then
       kmin = lbound(rhor,3)+1     
    end if
 else if(Symmetry == AXISYM) then
    if(m==1 .and. X(1,1,1).lt.0.D0) then
       kmin = lbound(rhor,1)+1
    else if(m==3 .and. Z(1,1,1).lt.0.D0) then
       kmin = lbound(rhor,3)+1     
    end if    
 end if

 rhot = rhor
 Pt   = Pr
 vxt  = vxr
 vyt  = vyr
 vzt  = vzr

 rhor = rhol
 Pr   = Pl
 vxr  = vxl
 vyr  = vyl
 vzr  = vzl

 if(m==1) then
    
    do k=kmin+1, kmax
       rhol(k,:,:) = rhot(k-1,:,:)
       Pl(k,:,:)   = Pt(k-1,:,:)
       vxl(k,:,:)  = vxt(k-1,:,:)
       vyl(k,:,:)  = vyt(k-1,:,:)
       vzl(k,:,:)  = vzt(k-1,:,:)
    end do

 else if (m==2) then

    do k=kmin+1, kmax
       rhol(:,k,:) = rhot(:,k-1,:)
       Pl(:,k,:)   = Pt(:,k-1,:)
       vxl(:,k,:)  = vxt(:,k-1,:)
       vyl(:,k,:)  = vyt(:,k-1,:)
       vzl(:,k,:)  = vzt(:,k-1,:)
    end do

 else if (m==3) then

    do k=kmin+1, kmax
       rhol(:,:,k) = rhot(:,:,k-1)
       Pl(:,:,k)   = Pt(:,:,k-1)
       vxl(:,:,k)  = vxt(:,:,k-1)
       vyl(:,:,k)  = vyt(:,:,k-1)
       vzl(:,:,k)  = vzt(:,:,k-1)
    end do

 else
    write(*,*) "m not in proper range"
 end if

end subroutine ppm_shift
