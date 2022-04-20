!-----------------------------------------------------------------------------
!
! routines for evolving the conducting fluid
! in the presence of a magnetic field
!
!-----------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------
!
! Compute the flux for advecting rho_star, tau (Font's energy variable), 
!  and S_i .
!
!-----------------------------------------------------------------------------
!
subroutine mhdflux_hybrid(ext,X,f_rho_i,f_tau_i,f_Sx_i,f_Sy_i,f_Sz_i,Pr,Pl,  &
     rho_br,rho_bl, vxr,vxl, vyr,vyl, vzr,vzl, Bxr,Bxl,Byr,Byl,Bzr,Bzl, &
     v02r,v02l,gupxx_f,gupyy_f,gupzz_f,cmax,cmin,i, &
     alpha_f, betax_f, betay_f, betaz_f, &
     gxx_f, gxy_f, gxz_f, gyy_f, gyz_f, gzz_f, phi_f, Symmetry, & 
     neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th,enable_HARM_energyvariable)
  implicit none

  interface
     subroutine find_cplus_cminus(ext,cplus,cminus,v02,u0,vi,lapse,shifti,phi,gupii)
       implicit none
       integer, dimension(3)                      :: ext
       real*8, dimension(ext(1),ext(2),ext(3))    :: cplus,cminus,v02,u0,vi
       real*8, dimension(ext(1),ext(2),ext(3))    :: lapse,shifti,phi,gupii
     end subroutine find_cplus_cminus
     subroutine flux_hll(ext, qr,ql, Fr,Fl, tau_ij, cmax,cmin)
       implicit none
       integer, dimension(3)                      :: ext
       real*8, dimension(ext(1),ext(2),ext(3))    :: qr,ql, Fr,Fl, tau_ij, cmax,cmin
     end subroutine flux_hll
     subroutine compute_pcold_epscold(rhob, P_cold, eps_cold, &
          neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
       implicit none
       integer :: neos,ergo_star
       real*8  :: rhob, P_cold, eps_cold, ergo_sigma
       real*8, dimension(neos) :: rho_tab,P_tab,eps_tab
       real*8, dimension(neos+1) :: k_tab, gamma_tab
     end subroutine compute_pcold_epscold
  end interface

  integer, dimension(3)					:: ext
  integer						:: i,imax,imin,m
  integer                                                :: Symmetry
  real*8, dimension(ext(1),ext(2),ext(3))		:: X
  real*8, dimension(ext(1),ext(2),ext(3))		:: f_rho_i, f_tau_i
  real*8, dimension(ext(1),ext(2),ext(3))                :: f_Sx_i,f_Sy_i,f_Sz_i
  real*8, dimension(ext(1),ext(2),ext(3))                :: rhor,rhol,er,el
  real*8, dimension(ext(1),ext(2),ext(3))                :: Pr,Pl,hr1,hl1
  real*8, dimension(ext(1),ext(2),ext(3))                :: rho_br,rho_bl
  real*8, dimension(ext(1),ext(2),ext(3))                :: vxr,vxl,vyr,vyl
  real*8, dimension(ext(1),ext(2),ext(3))                :: vzr,vzl
  real*8, dimension(ext(1),ext(2),ext(3))                :: v02l, v02r
  real*8, dimension(ext(1),ext(2),ext(3))                :: Bxr,Bxl,Byr,Byl
  real*8, dimension(ext(1),ext(2),ext(3))                :: Bzr,Bzl
  real*8, dimension(ext(1),ext(2),ext(3))                :: gupxx_f, gupyy_f, gupzz_f
  real*8					                :: gamma_th,ergo_sigma
  integer                                                :: neos,ergo_star
  real*8, dimension(neos)				:: rho_tab,P_tab,eps_tab
  real*8, dimension(neos+1)                              :: k_tab, gamma_tab
  !
  real*8, dimension(ext(1),ext(2),ext(3))                :: cplusl,cminusl
  real*8, dimension(ext(1),ext(2),ext(3))                :: cplusr,cminusr
  real*8, dimension(ext(1),ext(2),ext(3))                :: cmin,cmax,Fr,Fl
  real*8, dimension(ext(1),ext(2),ext(3))                :: au0r1,au0l1, X_f
  real*8, dimension(ext(1),ext(2),ext(3))                :: u0r,u0l,sb2l
  real*8, dimension(ext(1),ext(2),ext(3))                :: sb2r,sbxr,sbyr
  real*8, dimension(ext(1),ext(2),ext(3))                :: sbzr,sbtr,sbtl
  real*8, dimension(ext(1),ext(2),ext(3))                :: sbxl,sbyl,sbzl
  real*8							:: dX2
  !
  real*8, dimension(ext(1),ext(2),ext(3))                :: alpha_f, betax_f
  real*8, dimension(ext(1),ext(2),ext(3))                :: betay_f, betaz_f
  real*8, dimension(ext(1),ext(2),ext(3))                :: gxx_f, gxy_f, gxz_f
  real*8, dimension(ext(1),ext(2),ext(3))                :: gyy_f, gyz_f, gzz_f
  real*8, dimension(ext(1),ext(2),ext(3))                :: phi_f, al, psi4,psi6
  real*8, dimension(ext(1),ext(2),ext(3))                :: u_xl,u_yl,u_zl
  real*8, dimension(ext(1),ext(2),ext(3))                :: sb_xl,sb_yl,sb_zl
  real*8, dimension(ext(1),ext(2),ext(3))                :: u_xr,u_yr,u_zr
  real*8, dimension(ext(1),ext(2),ext(3))                :: sb_xr,sb_yr,sb_zr
  real*8, dimension(ext(1),ext(2),ext(3))                :: dPcold_drho,eps_th
  real*8                                                 :: f1os4pi,rhob
  real*8                                                 :: P_cold,eps_cold,eps
  integer						:: jmin,jmax,kmin,kmax,j,k
  integer						:: ii,jj,kk, nn, enable_HARM_energyvariable
  integer, parameter					:: AXISYM = 4
  real*8, parameter					:: fac = 0.99d0
  logical						:: exit_do

  f1os4pi = 1.d0/sqrt(4.d0*acos(-1.d0))

  imax = ubound(vxr,1)
  imin = lbound(vxr,1)
  jmax = ubound(vxr,2)
  jmin = lbound(vxr,2)
  kmax = ubound(vxr,3)
  kmin = lbound(vxr,3)

  if (Symmetry==AXISYM .and. i==1) then
     dX2 = 0.5d0*(X(imin+1,1,1)-X(imin,1,1))
     do m=imin,imax
        X_f(m,:,:) = X(m,1,1) - dX2
     end do
  end if

  al = 1.d0 + alpha_f
  psi4 = exp(4.d0 * phi_f)

  ! Compute al*u0-1
  er = psi4*(gxx_f*(vxr + betax_f)**2 + &
       2.d0*gxy_f*(vxr + betax_f)*(vyr + betay_f) +         &
       2.d0*gxz_f*(vxr + betax_f)*(vzr + betaz_f) +         &
       gyy_f*(vyr + betay_f)**2 +                           &
       2.d0*gyz_f*(vyr + betay_f)*(vzr + betaz_f) +         &
       gzz_f*(vzr + betaz_f)**2 )/al**2

  ! *** Check for superluminal velocity ***
  where (er .gt. 1.d0)
     vxr = (vxr + betax_f)*sqrt(fac/er)-betax_f
     vyr = (vyr + betay_f)*sqrt(fac/er)-betay_f
     vzr = (vzr + betaz_f)*sqrt(fac/er)-betaz_f
     er = fac
  end where
  ! ***************************************
  el = sqrt(1.d0-er)
  au0r1 = er/el/(1.d0+el)
  u0r = (au0r1+1.d0)/al

  er = psi4*(gxx_f*(vxl + betax_f)**2 + &
       2.d0*gxy_f*(vxl + betax_f)*(vyl + betay_f) +         &
       2.d0*gxz_f*(vxl + betax_f)*(vzl + betaz_f) +         &
       gyy_f*(vyl + betay_f)**2 +                           &
       2.d0*gyz_f*(vyl + betay_f)*(vzl + betaz_f) +         &
       gzz_f*(vzl + betaz_f)**2 )/al**2
!  write(*,*) "er",i,er(5,2,2),vxl(5,2,2),vyl(5,2,2),vzl(5,2,2)
  ! *** Check for superluminal velocity ***
  where (er .gt. 1.d0)
     vxl = (vxl + betax_f)*sqrt(fac/er)-betax_f
     vyl = (vyl + betay_f)*sqrt(fac/er)-betay_f
     vzl = (vzl + betaz_f)*sqrt(fac/er)-betaz_f
     er = fac
  end where
  ! ***************************************
  el = sqrt(1.d0-er)
  au0l1 = er/el/(1.d0+el)
  u0l = (au0l1+1.d0)/al


  !  do kk=kmin,kmax
  !     do jj=jmin,jmax
  !	do ii=imin,imax
  !	   rhob = rho_br(ii,jj,kk) 
  !
  !           nn = 1
  !           exit_do = .FALSE.
  !           do
  !             if (rhob .lt. rho_tab(nn)) then
  !                exit_do = .TRUE.
  !                P_cold = k_tab(nn)*rhob**gamma_tab(nn)
  !		if (rhob .ne. 0.d0) then
  !                   dPcold_drho(ii,jj,kk) = gamma_tab(nn)*P_cold/rhob
  !		else
  !		   dPcold_drho(ii,jj,kk) = 0.d0
  !		   eps_cold = 0.d0
  !		end if
  !                if (nn==1) then
  !		   if (rhob .ne. 0.d0) &
  !                      eps_cold = P_cold/rhob/(gamma_tab(nn)-1.d0)
  !                else
  !                   eps_cold = eps_tab(nn) +  &
  !                 (P_cold/rhob - P_tab(nn-1)/rho_tab(nn-1))/(gamma_tab(nn)-1.d0)
  !                end if
  !             end if
  !             if (nn==neos .or. exit_do) exit
  !             nn = nn + 1
  !           end do
  !           if (rhob .gt. rho_tab(neos)) then
  !              P_cold = k_tab(neos+1)*rhob**gamma_tab(neos+1)
  !              eps_cold = eps_tab(neos) + (P_cold/rhob - P_tab(neos) &
  !                        /rho_tab(neos))/(gamma_tab(neos+1)-1.d0)
  !              dPcold_drho(ii,jj,kk) = gamma_tab(neos+1)*P_cold/rhob
  !           end if
  !	   eps_th(ii,jj,kk) = (Pr(ii,jj,kk)-P_cold)/(gamma_th-1.d0)/rhob
  !
  !	   eps = eps_cold + eps_th(ii,jj,kk)
  !	   if (rhob .ne. 0.d0) then
  !	      hr1(ii,jj,kk) = Pr(ii,jj,kk)/rhob + eps
  !	   else
  !	      hr1(ii,jj,kk) = 0.d0
  !	   end if
  !
  !	end do
  !     end do
  !  end do 

  !  write(*,*) "mhd_flux0:",k_tab,rho_br(16,2,2),gamma_tab
  where (rho_br .le. rho_tab(1)) 
     where (rho_br .eq. 0.d0)
        dPcold_drho = 0.d0
        eps_th = 0.d0
        hr1 = 0.d0
     elsewhere
        v02l = k_tab(1)*rho_br**gamma_tab(1)           ! v02l -> P_cold
        cplusl = v02l/rho_br/(gamma_tab(1)-1.d0)       ! cplusl -> eps_cold
        dPcold_drho = gamma_tab(1)*v02l/rho_br
        eps_th = (Pr - v02l)/(gamma_th-1.d0)/rho_br
        hr1 = cplusl + eps_th + Pr/rho_br
     end where
  end where

  do nn=2,neos 
     where (rho_br .le. rho_tab(nn) .and. rho_br .gt. rho_tab(nn-1)) 
        v02l = k_tab(nn)*rho_br**gamma_tab(nn)         ! v02l -> P_cold
        cplusl = eps_tab(nn-1) + (v02l/rho_br - P_tab(nn-1)/rho_tab(nn-1))/ & 
             (gamma_tab(nn)-1.d0)		  ! cplusl -> eps_cold
        dPcold_drho = gamma_tab(nn)*v02l/rho_br
        eps_th = (Pr - v02l)/(gamma_th-1.d0)/rho_br
        hr1 = cplusl + eps_th + Pr/rho_br
     end where
  end do

if (ergo_star .eq. 0) then
 where (rho_br .gt. rho_tab(neos))
       v02l = k_tab(neos+1)*rho_br**gamma_tab(neos+1)   ! v02 -> P_cold
       cplusl = eps_tab(neos) + (v02l/rho_br - P_tab(neos)/rho_tab(neos))/(gamma_tab(neos+1)-1.d0)               ! va2 -> eps_cold
       dPcold_drho = gamma_tab(neos+1)*v02l/rho_br
       eps_th = (Pr-v02l)/(gamma_th-1.d0)/rho_br
       hr1 = cplusl + eps_th + Pr/rho_br
 end where
else
 where (rho_br .gt. rho_tab(neos))
       v02l = ((ergo_sigma* (1+eps_tab(neos)+P_tab(neos)/rho_tab(neos))/rho_tab(neos)**ergo_sigma) * rho_br**(ergo_sigma+1) + P_tab(neos)- ergo_sigma*(1+eps_tab(neos))*rho_tab(neos))/(ergo_sigma+1)
       cplusl = (((1+eps_tab(neos)+P_tab(neos)/rho_tab(neos))/rho_tab(neos)**ergo_sigma) * rho_br**(ergo_sigma+1) -P_tab(neos) + ergo_sigma*((1+eps_tab(neos))*rho_tab(neos)))/((ergo_sigma+1)*rho_br)-1
       dPcold_drho = ergo_sigma*(1+eps_tab(neos)+P_tab(neos)/rho_tab(neos))/rho_tab(neos)**ergo_sigma*rho_br
       eps_th = (Pr-v02l)/(gamma_th-1.d0)/rho_br
       hr1 = cplusl + eps_th + Pr/rho_br
 end where
end if





  sbtr = f1os4pi*psi4*u0r/al*( gxx_f*Bxr*(vxr+betax_f) + &
       gxy_f*( Bxr*(vyr+betay_f) + Byr*(vxr+betax_f) ) + &
       gxz_f*( Bxr*(vzr+betaz_f) + Bzr*(vxr+betax_f) ) + &
       gyy_f*Byr*(vyr+betay_f) + &
       gyz_f*( Byr*(vzr+betaz_f) + Bzr*(vyr+betay_f) ) + &
       gzz_f*Bzr*(vzr+betaz_f) )
  sbxr = f1os4pi*Bxr/(al*u0r) + sbtr*vxr
  sbyr = f1os4pi*Byr/(al*u0r) + sbtr*vyr
  sbzr = f1os4pi*Bzr/(al*u0r) + sbtr*vzr
  sb2r = -(al*sbtr)**2 + psi4*( gxx_f*(sbxr+betax_f*sbtr)**2 + &
       2.d0*gxy_f*(sbxr+betax_f*sbtr)*(sbyr+betay_f*sbtr) + &
       2.d0*gxz_f*(sbxr+betax_f*sbtr)*(sbzr+betaz_f*sbtr) + &
       gyy_f*(sbyr+betay_f*sbtr)**2 + &
       2.d0*gyz_f*(sbyr+betay_f*sbtr)*(sbzr+betaz_f*sbtr) + &
       gzz_f*(sbzr+betaz_f*sbtr)**2 )

  ! Compute v02r
  where (rho_br .gt. 0.d0)
     cplusl = (dPcold_drho + gamma_th*(gamma_th-1.d0)*eps_th)/(1.d0+hr1) 
     cminusl = sb2r/(sb2r + rho_br*(1.d0+hr1))
     v02r = cminusl + cplusl*(1.d0-cminusl)
  elsewhere
     v02r = sb2r/(sb2r+1.d-300)
  end where

  !  do kk=kmin,kmax
  !     do jj=jmin,jmax
  !        do ii=imin,imax
  !           rhob = rho_bl(ii,jj,kk)
  !
  !           nn = 1
  !           exit_do = .FALSE.
  !           do
  !             if (rhob .lt. rho_tab(nn)) then
  !                exit_do = .TRUE.
  !                P_cold = k_tab(nn)*rhob**gamma_tab(nn)
  !                if (rhob .ne. 0.d0) then
  !                   dPcold_drho(ii,jj,kk) = gamma_tab(nn)*P_cold/rhob
  !                else
  !                   dPcold_drho(ii,jj,kk) = 0.d0
  !                   eps_cold = 0.d0
  !                end if
  !                if (nn==1) then
  !                   if (rhob .ne. 0.d0) &
  !                      eps_cold = P_cold/rhob/(gamma_tab(nn)-1.d0)
  !                else
  !                   eps_cold = eps_tab(nn) +  &
  !                 (P_cold/rhob - P_tab(nn-1)/rho_tab(nn-1))/(gamma_tab(nn)-1.d0)
  !                end if
  !             end if
  !             if (nn==neos .or. exit_do) exit
  !             nn = nn + 1
  !           end do
  !           if (rhob .gt. rho_tab(neos)) then
  !              P_cold = k_tab(neos+1)*rhob**gamma_tab(neos+1)
  !              eps_cold = eps_tab(neos) + (P_cold/rhob - P_tab(neos) &
  !                        /rho_tab(neos))/(gamma_tab(neos+1)-1.d0)
  !              dPcold_drho(ii,jj,kk) = gamma_tab(neos+1)*P_cold/rhob
  !           end if
  !           eps_th(ii,jj,kk) = (Pl(ii,jj,kk)-P_cold)/(gamma_th-1.d0)/rhob
  !
  !           eps = eps_cold + eps_th(ii,jj,kk)
  !	   if (rhob .ne. 0.d0) then
  !              hl1(ii,jj,kk) = Pl(ii,jj,kk)/rhob + eps
  !           else
  !              hl1(ii,jj,kk) = 0.d0
  !           end if
  !        end do
  !     end do
  !  end do

  where (rho_bl .le. rho_tab(1)) 
     where (rho_bl .eq. 0.d0)
        dPcold_drho = 0.d0
        eps_th = 0.d0
        hl1 = 0.d0
     elsewhere
        v02l = k_tab(1)*rho_bl**gamma_tab(1)           ! v02l -> P_cold
        cplusl = v02l/rho_bl/(gamma_tab(1)-1.d0)       ! cplusl -> eps_cold
        dPcold_drho = gamma_tab(1)*v02l/rho_bl
        eps_th = (Pl - v02l)/(gamma_th-1.d0)/rho_bl
        hl1 = cplusl + eps_th + Pl/rho_bl
     end where
  end where

  do nn=2,neos 
     where (rho_bl .le. rho_tab(nn) .and. rho_bl .gt. rho_tab(nn-1)) 
        v02l = k_tab(nn)*rho_bl**gamma_tab(nn)         ! v02l -> P_cold
        cplusl = eps_tab(nn-1) + (v02l/rho_bl - P_tab(nn-1)/rho_tab(nn-1))/ & 
             (gamma_tab(nn)-1.d0)		  ! cplusl -> eps_cold
        dPcold_drho = gamma_tab(nn)*v02l/rho_bl
        eps_th = (Pl - v02l)/(gamma_th-1.d0)/rho_bl
        hl1 = cplusl + eps_th + Pl/rho_bl
     end where
  end do


if (ergo_star .eq. 0) then
 where (rho_bl .gt. rho_tab(neos))
       v02l = k_tab(neos+1)*rho_bl**gamma_tab(neos+1)   ! v02 -> P_cold
       cplusl = eps_tab(neos) + (v02l/rho_bl - P_tab(neos)/rho_tab(neos))/(gamma_tab(neos+1)-1.d0)               ! va2 -> eps_cold
       dPcold_drho = gamma_tab(neos+1)*v02l/rho_bl
       eps_th = (Pl-v02l)/(gamma_th-1.d0)/rho_bl
       hl1 = cplusl + eps_th + Pl/rho_bl
 end where
else
 where (rho_bl .gt. rho_tab(neos))
       v02l = ((ergo_sigma* (1+eps_tab(neos)+P_tab(neos)/rho_tab(neos))/rho_tab(neos)**ergo_sigma) * rho_bl**(ergo_sigma+1) + P_tab(neos) - ergo_sigma*(1+eps_tab(neos))*rho_tab(neos))/(ergo_sigma+1)
       cplusl = (((1+eps_tab(neos)+P_tab(neos)/rho_tab(neos))/rho_tab(neos)**ergo_sigma) * rho_bl**(ergo_sigma+1) -P_tab(neos) + ergo_sigma*((1+eps_tab(neos))*rho_tab(neos)))/((ergo_sigma+1)*rho_bl)-1
       dPcold_drho = ergo_sigma*(1+eps_tab(neos)+P_tab(neos)/rho_tab(neos))/rho_tab(neos)**ergo_sigma*rho_bl
       eps_th = (Pl-v02l)/(gamma_th-1.d0)/rho_bl
       hl1 = cplusl + eps_th + Pl/rho_bl
 end where
end if






  ! Here sba = b^a defined in Gammie's paper
  sbtl = f1os4pi*psi4*u0l/al*( gxx_f*Bxl*(vxl+betax_f) + &
       gxy_f*( Bxl*(vyl+betay_f) + Byl*(vxl+betax_f) ) + &
       gxz_f*( Bxl*(vzl+betaz_f) + Bzl*(vxl+betax_f) ) + &
       gyy_f*Byl*(vyl+betay_f) + &
       gyz_f*( Byl*(vzl+betaz_f) + Bzl*(vyl+betay_f) ) + &
       gzz_f*Bzl*(vzl+betaz_f) )
  sbxl = f1os4pi*Bxl/(al*u0l) + sbtl*vxl
  sbyl = f1os4pi*Byl/(al*u0l) + sbtl*vyl
  sbzl = f1os4pi*Bzl/(al*u0l) + sbtl*vzl
  ! Compute b^2
  sb2l = -(al*sbtl)**2 + psi4*( gxx_f*(sbxl+betax_f*sbtl)**2 + &
       2.d0*gxy_f*(sbxl+betax_f*sbtl)*(sbyl+betay_f*sbtl) + &
       2.d0*gxz_f*(sbxl+betax_f*sbtl)*(sbzl+betaz_f*sbtl) + &
       gyy_f*(sbyl+betay_f*sbtl)**2 + &
       2.d0*gyz_f*(sbyl+betay_f*sbtl)*(sbzl+betaz_f*sbtl) + &
       gzz_f*(sbzl+betaz_f*sbtl)**2 )
  !  write(*,*) "inside mhdflux1", sbtl(1,2,2),psi4(1,2,2),u0l(1,2,2),al(1,2,2),gyy_f(1,2,2),Byl(1,2,2),vyl(1,2,2),betay_f(1,2,2)

  ! Compute v02l
  where (rho_bl .gt. 0.d0)
     cplusl = (dPcold_drho + gamma_th*(gamma_th-1.d0)*eps_th)/(1.d0+hl1)
     cminusl = sb2l/(sb2l + rho_bl*(1.d0+hl1))
     v02l = cminusl + cplusl*(1.d0-cminusl)
  elsewhere
     v02l = sb2l/(sb2l+1.d-300)
  end where
  !  write(*,*) "inside mhdflux1a", v02l(1,2,2),cplusl(1,2,2),cminusl(1,2,2),sb2l(1,2,2),dPcold_drho(1,2,2),gamma_th,hl1(1,2,2),rho_bl(1,2,2)

  if (i==1) then
     call find_cplus_cminus(ext,cplusl,cminusl,v02l,u0l,  &
          vxl,alpha_f,betax_f,phi_f,gupxx_f)

     !write(*,*) "cpluslx",u0l(10,2,10),vxl(10,2,10),v02l(10,2,10),phi_f(10,2,10),gupxx_f(10,2,10),betax_f(10,2,10),alpha_f(10,2,10)

     call find_cplus_cminus(ext,cplusr,cminusr,v02r,u0r,  &
          vxr,alpha_f,betax_f,phi_f,gupxx_f)
  else if (i==2) then
     call find_cplus_cminus(ext,cplusl,cminusl,v02l,u0l,  &
          vyl,alpha_f,betay_f,phi_f,gupyy_f)
     call find_cplus_cminus(ext,cplusr,cminusr,v02r,u0r,  &
          vyr,alpha_f,betay_f,phi_f,gupyy_f)
  else
     call find_cplus_cminus(ext,cplusl,cminusl,v02l,u0l,  &
          vzl,alpha_f,betaz_f,phi_f,gupzz_f)
     call find_cplus_cminus(ext,cplusr,cminusr,v02r,u0r,  &
          vzr,alpha_f,betaz_f,phi_f,gupzz_f)
  end if

  !write(*,*) "cplusl",i,cplusl(30,30,2),v02l(30,30,2),u0l(30,30,2),vxl(30,30,2),alpha_f(30,30,2),betax_f(30,30,2),phi_f(30,30,2),gupxx_f(30,30,2)

  cmax = max(0.d0,cplusl,cplusr)
  cmin = -min(0.d0,cminusl,cminusr)

  !  write(*,*) "cmaxcmin",cmin(2,2,2),cminusl(2,2,2),cminusr(2,2,2)

  !Above lines of code will produce NaN's in the HLL step if cmax=cmin=0...
  !  This problem occurs in axisymmetry when we are in the x<0 ghostzones and have mhd_evolve::em_evolve_enable = "no".
  !   Apparently with the B field evolution turned off, cplus's and cminus's are not set properly.
  !  if(Symmetry==AXISYM) then
  !     where(cmin==0 .and. cmax==0)
  !        cmax = 1.D-100
  !     end where
  !  end if

  psi6 = exp(6.d0 * phi_f)

  ! density flux
  if (i==1) then
     Fr = al*psi6*rho_br*vxr*u0r
     Fl = al*psi6*rho_bl*vxl*u0l
  elseif (i==2) then
     Fr = al*psi6*rho_br*vyr*u0r
     Fl = al*psi6*rho_bl*vyl*u0l
  elseif (i==3) then
     Fr = al*psi6*rho_br*vzr*u0r
     Fl = al*psi6*rho_bl*vzl*u0l
  end if

  rhor = al*psi6*rho_br*u0r
  rhol = al*psi6*rho_bl*u0l

  ! write(*,*) "inside mhd_flux_hybrid, before drho_b_m comp: ",rhor(1,2,2),rhol(1,2,2),Fr(1,2,2),Fl(1,2,2),cmax(1,2,2),cmin(1,2,2)
  ! write(*,*) "hi.",(cmin(1,2,2)*Fr(1,2,2) + cmax(1,2,2)*Fl(1,2,2) - cmin(1,2,2)*cmax(1,2,2)*(rhor(1,2,2)-rhol(1,2,2)) ),(cmax(1,2,2) + cmin(1,2,2))
  call flux_hll(ext, rhor,rhol, Fr,Fl, f_rho_i, cmax,cmin)
  ! write(*,*) "frhoi",i,f_rho_i(10,2,10)

  ! write(*,*) "inside mhd_flux_hybrid, after drho_b_m comp: ",f_rho_i(1,2,2)

  ! energy flux
  u_xl = u0l*psi4*( gxx_f*(vxl+betax_f) + gxy_f*(vyl+betay_f) + &
       gxz_f*(vzl+betaz_f) )
  u_yl = u0l*psi4*( gxy_f*(vxl+betax_f) + gyy_f*(vyl+betay_f) + &
       gyz_f*(vzl+betaz_f) )
  u_zl = u0l*psi4*( gxz_f*(vxl+betax_f) + gyz_f*(vyl+betay_f) + &
       gzz_f*(vzl+betaz_f) )

!  write(*,*) "u_zl:",u0l(5,2,2)*psi4(5,2,2),gxz_f(5,2,2)*(vxl(5,2,2)+betax_f(5,2,2)),vxl(5,2,2),vyl(5,2,2),vzl(5,2,2)

  

  sb_xl = psi4*( gxx_f*(sbxl+sbtl*betax_f) + gxy_f*(sbyl+sbtl*betay_f) + &
       gxz_f*(sbzl+sbtl*betaz_f) )
  sb_yl = psi4*( gxy_f*(sbxl+sbtl*betax_f) + gyy_f*(sbyl+sbtl*betay_f) + &
       gyz_f*(sbzl+sbtl*betaz_f) )
  sb_zl = psi4*( gxz_f*(sbxl+sbtl*betax_f) + gyz_f*(sbyl+sbtl*betay_f) + &
       gzz_f*(sbzl+sbtl*betaz_f) )

  u_xr = u0r*psi4*( gxx_f*(vxr+betax_f) + gxy_f*(vyr+betay_f) + &
       gxz_f*(vzr+betaz_f) )
  u_yr = u0r*psi4*( gxy_f*(vxr+betax_f) + gyy_f*(vyr+betay_f) + &
       gyz_f*(vzr+betaz_f) )
  u_zr = u0r*psi4*( gxz_f*(vxr+betax_f) + gyz_f*(vyr+betay_f) + &
       gzz_f*(vzr+betaz_f) )
  sb_xr = psi4*( gxx_f*(sbxr+sbtr*betax_f) + gxy_f*(sbyr+sbtr*betay_f) + &
       gxz_f*(sbzr+sbtr*betaz_f) )
  sb_yr = psi4*( gxy_f*(sbxr+sbtr*betax_f) + gyy_f*(sbyr+sbtr*betay_f) + &
       gyz_f*(sbzr+sbtr*betaz_f) )
  sb_zr = psi4*( gxz_f*(sbxr+sbtr*betax_f) + gyz_f*(sbyr+sbtr*betay_f) + &
       gzz_f*(sbzr+sbtr*betaz_f) )

  if (enable_HARM_energyvariable==1) then
     ! Temporarily store 1+u_0r and 1+u_0l in er and el
     el = -alpha_f - al*au0l1 + u_xl*betax_f + u_yl*betay_f + u_zl*betaz_f
     er = -alpha_f - al*au0r1 + u_xr*betax_f + u_yr*betay_f + u_zr*betaz_f
     ! Temporarily store b_tr and b_tl to cplusr and cplusl
     cplusl = -al*al*sbtl + sb_xl*betax_f + sb_yl*betay_f + sb_zl*betaz_f
     cplusr = -al*al*sbtr + sb_xr*betax_f + sb_yr*betay_f + sb_zr*betaz_f
     ! Now compute energy flux
     if (i==1) then
        Fl = -rhol*vxl*(el + (el-1.d0)*hl1) - al*psi6*( & 
             sb2l*u0l*vxl*(el-1.d0) - sbxl*cplusl )
        Fr = -rhor*vxr*(er + (er-1.d0)*hr1) - al*psi6*( &
             sb2r*u0r*vxr*(er-1.d0) - sbxr*cplusr )
     elseif (i==2) then
        Fl = -rhol*vyl*(el + (el-1.d0)*hl1) - al*psi6*( &
             sb2l*u0l*vyl*(el-1.d0) - sbyl*cplusl )
        Fr = -rhor*vyr*(er + (er-1.d0)*hr1) - al*psi6*( &
             sb2r*u0r*vyr*(er-1.d0) - sbyr*cplusr )
     else
        Fl = -rhol*vzl*(el + (el-1.d0)*hl1) - al*psi6*( &
             sb2l*u0l*vzl*(el-1.d0) - sbzl*cplusl )
        Fr = -rhor*vzr*(er + (er-1.d0)*hr1) - al*psi6*( &
             sb2r*u0r*vzr*(er-1.d0) - sbzr*cplusr )
     end if
     ! Finally, compute taur and taul and store them in el and er
     el = -rhol*(el + (el-1.d0)*hl1)-al*psi6*(Pl + (u0l*(el-1.d0)+0.5d0)*sb2l & 
          -sbtl*cplusl )
     er = -rhor*(er + (er-1.d0)*hr1)-al*psi6*(Pr + (u0r*(er-1.d0)+0.5d0)*sb2r &
          -sbtr*cplusr )
  else
     el = (au0l1+hl1+au0l1*hl1)*rhol-psi6*Pl
     if (i==1) then
        Fl = el*vxl + Pl*psi6*(vxl+betax_f) + psi6*( (al*u0l)**2 *sb2l*vxl &
             + 0.5d0*betax_f*sb2l - al**2*sbtl*sbxl )
     elseif (i==2) then
        Fl = el*vyl + Pl*psi6*(vyl+betay_f) + psi6*( (al*u0l)**2 *sb2l*vyl &
             + 0.5d0*betay_f*sb2l - al**2*sbtl*sbyl )
     else
        Fl = el*vzl + Pl*psi6*(vzl+betaz_f) + psi6*( (al*u0l)**2 *sb2l*vzl &
             + 0.5d0*betaz_f*sb2l - al**2*sbtl*sbzl )
     end if
     el = el + psi6*(sb2l*(al*u0l)**2 - 0.5d0*sb2l - (al*sbtl)**2)

     er = (au0r1+hr1+au0r1*hr1)*rhor-psi6*Pr
     if (i==1) then
        Fr = er*vxr + Pr*psi6*(vxr+betax_f) + psi6*( (al*u0r)**2 *sb2r*vxr &
             + 0.5d0*betax_f*sb2r - al**2*sbtr*sbxr )
     elseif (i==2) then
        Fr = er*vyr + Pr*psi6*(vyr+betay_f) + psi6*( (al*u0r)**2 *sb2r*vyr &
             + 0.5d0*betay_f*sb2r - al**2*sbtr*sbyr )
     else
        Fr = er*vzr + Pr*psi6*(vzr+betaz_f) + psi6*( (al*u0r)**2 *sb2r*vzr &
             + 0.5d0*betaz_f*sb2r - al**2*sbtr*sbzr )
     end if
     er = er + psi6*(sb2r*(al*u0r)**2 - 0.5d0*sb2r - (al*sbtr)**2)
  end if


  call flux_hll(ext, er,el, Fr,Fl, f_tau_i, cmax,cmin)
  ! write(*,*) "inside mhd_flux_hybrid: ",er(1,2,2),el(1,2,2),Fr(1,2,2),Fl(1,2,2),cmax(1,2,2),cmin(1,2,2),f_tau_i(1,2,2)

  if (Symmetry==AXISYM .and. i==1) then
     f_rho_i = f_rho_i * X_f
     f_tau_i = f_tau_i * X_f
  end if

  ! Flux for S_j

  ! Flux for S_x
  if (i==1) then
     Fl = al*psi6*( (rho_bl*(1.d0+hl1)+sb2l)*u0l*vxl*u_xl  &
          + Pl+0.5d0*sb2l - sbxl*sb_xl )
     Fr = al*psi6*( (rho_br*(1.d0+hr1)+sb2r)*u0r*vxr*u_xr  &
          + Pr+0.5d0*sb2r - sbxr*sb_xr )
  elseif (i==2) then
     Fl = al*psi6*( (rho_bl*(1.d0+hl1)+sb2l)*u0l*vyl*u_xl  &
          - sbyl*sb_xl )
     Fr = al*psi6*( (rho_br*(1.d0+hr1)+sb2r)*u0r*vyr*u_xr  &
          - sbyr*sb_xr )
  else
     Fl = al*psi6*( (rho_bl*(1.d0+hl1)+sb2l)*u0l*vzl*u_xl  &
          - sbzl*sb_xl )
     Fr = al*psi6*( (rho_br*(1.d0+hr1)+sb2r)*u0r*vzr*u_xr  &
          - sbzr*sb_xr )
  end if
  el = al*psi6*( (rho_bl*(1.d0+hl1)+sb2l)*u0l*u_xl - sbtl*sb_xl )
  er = al*psi6*( (rho_br*(1.d0+hr1)+sb2r)*u0r*u_xr - sbtr*sb_xr )

  ! write(*,*) "before x flux hll:",psi4(17,2,15),gxx_f(17,2,15),sbxl(17,2,15),sbtl(17,2,15),betax_f(17,2,15),gxy_f(17,2,15),sbyl(17,2,15),sbtl(17,2,15),betay_f(17,2,15),gxz_f(17,2,15),sbzl(17,2,15),sbtl(17,2,15),betaz_f(17,2,15),vxl(17,2,15),u0l(17,2,15),Bxl(17,2,15)
  ! write(*,*) "before x flux hll:",al(17,2,15),psi6(17,2,15),rho_bl(17,2,15),hl1(17,2,15),sb2l(17,2,15),u0l(17,2,15),vzl(17,2,15),u_xl(17,2,15),sbzl(17,2,15),sb_xl(17,2,15)

  ! write(*,*) "before x flux hll:",er(17,2,15),el(17,2,15),Fr(17,2,15),Fl(17,2,15),cmax(17,2,15),cmin(17,2,15)
  call flux_hll(ext, er,el, Fr,Fl, f_Sx_i, cmax,cmin)
  ! write(*,*) "post x flux hll:",f_Sx_i(17,2,15)

  ! Flux for S_y
  if (i==1) then
     Fl = al*psi6*( (rho_bl*(1.d0+hl1)+sb2l)*u0l*vxl*u_yl  &
          - sbxl*sb_yl )
     Fr = al*psi6*( (rho_br*(1.d0+hr1)+sb2r)*u0r*vxr*u_yr  &
          - sbxr*sb_yr )
  elseif (i==2) then
     Fl = al*psi6*( (rho_bl*(1.d0+hl1)+sb2l)*u0l*vyl*u_yl  &
          + Pl+0.5d0*sb2l - sbyl*sb_yl )
     Fr = al*psi6*( (rho_br*(1.d0+hr1)+sb2r)*u0r*vyr*u_yr  &
          + Pr+0.5d0*sb2r - sbyr*sb_yr )
  else
     Fl = al*psi6*( (rho_bl*(1.d0+hl1)+sb2l)*u0l*vzl*u_yl  &
          - sbzl*sb_yl )
     Fr = al*psi6*( (rho_br*(1.d0+hr1)+sb2r)*u0r*vzr*u_yr  &
          - sbzr*sb_yr )
  end if
  el = al*psi6*( (rho_bl*(1.d0+hl1)+sb2l)*u0l*u_yl - sbtl*sb_yl )
  er = al*psi6*( (rho_br*(1.d0+hr1)+sb2r)*u0r*u_yr - sbtr*sb_yr )

  call flux_hll(ext, er,el, Fr,Fl, f_Sy_i, cmax,cmin)

  ! Flux for S_z
  if (i==1) then
     Fl = al*psi6*( (rho_bl*(1.d0+hl1)+sb2l)*u0l*vxl*u_zl  &
          - sbxl*sb_zl )
     Fr = al*psi6*( (rho_br*(1.d0+hr1)+sb2r)*u0r*vxr*u_zr  &
          - sbxr*sb_zr )
  elseif (i==2) then
     Fl = al*psi6*( (rho_bl*(1.d0+hl1)+sb2l)*u0l*vyl*u_zl  &
          - sbyl*sb_zl )
     Fr = al*psi6*( (rho_br*(1.d0+hr1)+sb2r)*u0r*vyr*u_zr  &
          - sbyr*sb_zr )
  else
     Fl = al*psi6*( (rho_bl*(1.d0+hl1)+sb2l)*u0l*vzl*u_zl  &
          + Pl+0.5d0*sb2l - sbzl*sb_zl )
     Fr = al*psi6*( (rho_br*(1.d0+hr1)+sb2r)*u0r*vzr*u_zr  &
          + Pr+0.5d0*sb2r - sbzr*sb_zr )
  end if
  el = al*psi6*( (rho_bl*(1.d0+hl1)+sb2l)*u0l*u_zl - sbtl*sb_zl )
!  write(*,*) "el:",rho_bl(30,30,2),hl1(30,30,2),sb2l(30,30,2),u0l(30,30,2),u_zl(30,30,2)

  er = al*psi6*( (rho_br*(1.d0+hr1)+sb2r)*u0r*u_zr - sbtr*sb_zr )

  call flux_hll(ext, er,el, Fr,Fl, f_Sz_i, cmax,cmin)
  !write(*,*) "FLUXES:",cmin(30,30,2),Fr(30,30,2),cmax(30,30,2),er(30,30,2),el(30,30,2)

  !bug here?????? (.and. i==2?) answer: No.  i is a direction here.  Duh.
  if (Symmetry==AXISYM .and. i==1) then
     f_Sx_i = f_Sx_i * X_f
     f_Sy_i = f_Sy_i * X_f**2
     f_Sz_i = f_Sz_i * X_f
  end if

  ! write(*,*) "FLUXES:",f_rho_i(10,2,10),f_tau_i(10,2,10),f_Sx_i(10,2,10),f_Sy_i(10,2,10),f_Sz_i(10,2,10)

end subroutine mhdflux_hybrid
