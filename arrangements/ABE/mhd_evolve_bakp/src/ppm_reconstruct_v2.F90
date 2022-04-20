!-----------------------------------------------------------------!
! Calculate center differences df = (f_{j+1} - f_{j-1} )/2
!                           d^2 f = f_{j+1} - 2 f_j + f_{j-1}
!-----------------------------------------------------------------!
!
subroutine find_centderivs_f(ext,nghostzones,f,df,d2f,m,sym, & 
                              Symmetry,cell_centering_enabled)
   implicit none
   integer, dimension(3)                          :: ext,nghostzones
   real*8, dimension(ext(1),ext(2),ext(3))        :: f,df,d2f
   integer                                        :: imax,jmax,kmax,i,j,k,m,jmin
   integer                                        :: Symmetry
   integer                                        :: sym_index,cell_centering_enabled
   integer, parameter :: AXISYM=4, OCTANT=2
   real*8, dimension(3) 			  :: sym
!
   imax = ubound(f,1)
   jmax = ubound(f,2)
   kmax = ubound(f,3)
   jmin = 1
   if (Symmetry==AXISYM) then 
      jmin = 2
      jmax = 2
   end if

   if (cell_centering_enabled==1) then
      sym_index = 2*nghostzones(m) + 1
      ! In axisymmetry, nghostzones=2, though [# symmetry ghostzones] = 1.
      if (Symmetry==4 .and. m==1) sym_index = 3
   else 
      sym_index = 2*nghostzones(m) + 2
   end if

   if (m==1) then 
      !$omp parallel do
      do k=1,kmax
	 do j=jmin,jmax
            i = 1
	    if (Symmetry==AXISYM .or. Symmetry==OCTANT) then 
	       df(i,j,k) = 0.5d0*(f(i+1,j,k) - sym(1)*f(sym_index,j,k))
	       d2f(i,j,k) = f(i+1,j,k) - 2.d0*f(i,j,k) + sym(1)*f(sym_index,j,k)
 	    else
               df(i,j,k) = -1.5d0*f(i,j,k) + 2.d0*f(i+1,j,k) - 0.5d0*f(i+2,j,k)
               d2f(i,j,k) = f(i,j,k) - 2.d0*f(i+1,j,k) + f(i+2,j,k)
	    end if
	    do i=2,imax-1
	       df(i,j,k) = 0.5d0*(f(i+1,j,k) - f(i-1,j,k))
	       d2f(i,j,k) = f(i+1,j,k) - 2.d0*f(i,j,k) + f(i-1,j,k)
	    end do
	    i = imax
	    df(i,j,k) = 1.5d0*f(i,j,k) - 2.d0*f(i-1,j,k) + 0.5d0*f(i-2,j,k)
	    d2f(i,j,k) = f(i,j,k) - 2.d0*f(i-1,j,k) + f(i-2,j,k)
	 end do
      end do
      !$omp end parallel do
   else if (m==2) then 
      !$omp parallel do
      do k=1,kmax
         do i=1,imax
            j = 1
	    if (Symmetry==OCTANT) then
               df(i,j,k) = 0.5d0*(f(i,j+1,k) - sym(2)*f(i,sym_index,k))
               d2f(i,j,k) = f(i,j+1,k) - 2.d0*f(i,j,k) + sym(2)*f(i,sym_index,k)
	    else
               df(i,j,k) = -1.5d0*f(i,j,k) + 2.d0*f(i,j+1,k) - 0.5d0*f(i,j+2,k)
               d2f(i,j,k) = f(i,j,k) - 2.d0*f(i,j+1,k) + f(i,j+2,k)
	    end if
            do j=2,jmax-1                
               df(i,j,k) = 0.5d0*(f(i,j+1,k) - f(i,j-1,k))
               d2f(i,j,k) = f(i,j+1,k) - 2.d0*f(i,j,k) + f(i,j-1,k)
            end do
            j = jmax
            df(i,j,k) = 1.5d0*f(i,j,k) - 2.d0*f(i,j-1,k) + 0.5d0*f(i,j-2,k)
            d2f(i,j,k) = f(i,j,k) - 2.d0*f(i,j-1,k) + f(i,j-2,k)
         end do
      end do
      !$omp end parallel do
   else
      !$omp parallel do
      do j=jmin,jmax
         do i=1,imax
            k = 1
	    if (Symmetry .ne. 0.d0) then 
               df(i,j,k) = 0.5d0*(f(i,j,k+1) - sym(3)*f(i,j,sym_index))
               d2f(i,j,k) = f(i,j,k+1) - 2.d0*f(i,j,k) + sym(3)*f(i,j,sym_index)
	    else 
               df(i,j,k) = -1.5d0*f(i,j,k) + 2.d0*f(i,j,k+1) - 0.5d0*f(i,j,k+2)
               d2f(i,j,k) = f(i,j,k) - 2.d0*f(i,j,k+1) + f(i,j,k+2)
	    end if
            do k=2,kmax-1
               df(i,j,k) = 0.5d0*(f(i,j,k+1) - f(i,j,k-1))
               d2f(i,j,k) = f(i,j,k+1) - 2.d0*f(i,j,k) + f(i,j,k-1)
            end do
            k = kmax
            df(i,j,k) = 1.5d0*f(i,j,k) - 2.d0*f(i,j,k-1) + 0.5d0*f(i,j,k-2)
            d2f(i,j,k) = f(i,j,k) - 2.d0*f(i,j,k-1) + f(i,j,k-2)
         end do
      end do
      !$omp end parallel do
   end if

end subroutine find_centderivs_f

!----------------------------------------------------------------------!
! Calculate a_{j+1/2}, set a_R = a_{j+1/2}, a_L  = a_{j-1/2} = a_R(j-1)
!----------------------------------------------------------------------!
!
subroutine ppm_find_face_vals_v2(ext,nghostzones,a,delta_a,ar,al,m,sym, & 
          	Symmetry,cell_centering_enabled)
 implicit none
 integer, dimension(3)                          :: ext,nghostzones
 integer, dimension(3)                          :: sym
 real*8, dimension(ext(1),ext(2),ext(3))        :: a,delta_a,ar,al
 integer                                        :: m,Symmetry
 integer                                        :: cell_centering_enabled
 integer					:: i,j,k,imax,jmax,kmax,sym_index,jmin
 integer, parameter 			        :: AXISYM = 4
 real*8, parameter :: sixth = 0.16666666666666666666666666667d0,eighth=0.125d0
!
 imax = ubound(a,1)
 jmax = ubound(a,2)
 kmax = ubound(a,3)
 jmin = 1
 if (Symmetry==AXISYM) then 
    jmin = 2
    jmax = 2
 end if

   if (cell_centering_enabled==1) then
      sym_index = 2*nghostzones(m) + 1
      ! In axisymmetry, nghostzones=2, though [# symmetry ghostzones] = 1.
      if (Symmetry==4 .and. m==1) sym_index = 3
   else
      sym_index = 2*nghostzones(m) + 2
   end if

   if (m==1) then
      !$omp parallel do
      do k=1,kmax
	 do j=jmin,jmax
	    do i=1,imax-1
	       ar(i,j,k) = a(i,j,k) + 0.5d0*(a(i+1,j,k)-a(i,j,k)) + &
                 sixth*(delta_a(i,j,k)-delta_a(i+1,j,k))
	    end do
	    ar(imax,j,k) = a(imax,j,k)
	 end do
      end do
      !$omp end parallel do
      !$omp parallel do
      do k=1,kmax
    	 do j=jmin,jmax
	    al(1,j,k) = sym(1)*a(sym_index,j,k) +  &
                  0.5d0*(a(1,j,k)-sym(1)*a(sym_index,j,k)) + &
                 sixth*(-sym(1)*delta_a(sym_index,j,k) - delta_a(1,j,k))
            do i=2,imax
               al(i,j,k) = ar(i-1,j,k)
            end do
	 end do
      end do
      !$omp end parallel do
   else if (m==2) then 
      !$omp parallel do
      do k=1,kmax
         do i=1,imax
            do j=1,jmax-1
               ar(i,j,k) = a(i,j,k) + 0.5d0*(a(i,j+1,k)-a(i,j,k)) + &
                 sixth*(delta_a(i,j,k)-delta_a(i,j+1,k))
            end do
            ar(i,jmax,k) = a(i,jmax,k)
         end do
      end do
      !$omp end parallel do
      !$omp parallel do
      do k=1,kmax
         do i=1,imax
            al(i,1,k) = sym(2)*a(i,sym_index,k) +  &
                   0.5d0*(a(i,1,k)-sym(2)*a(i,sym_index,k)) + &
                 sixth*(-sym(2)*delta_a(i,sym_index,k) - delta_a(i,1,k))
            do j=2,jmax
               al(i,j,k) = ar(i,j-1,k)
            end do
         end do
      end do
      !$omp end parallel do
   else 
      !$omp parallel do
      do j=jmin,jmax
         do i=1,imax
            do k=1,kmax-1
               ar(i,j,k) = a(i,j,k) + 0.5d0*(a(i,j,k+1)-a(i,j,k)) + &
                 sixth*(delta_a(i,j,k)-delta_a(i,j,k+1))
            end do
            ar(i,j,kmax) = a(i,j,kmax)
         end do
      end do
      !$omp end parallel do
      !$omp parallel do
      do j=jmin,jmax
         do i=1,imax
            al(i,j,1) = sym(3)*a(i,j,sym_index) +  & 
                 0.5d0*(a(i,j,1)-sym(3)*a(i,j,sym_index)) + &
                 sixth*(-sym(3)*delta_a(i,j,sym_index) - delta_a(i,j,1))
            do k=2,kmax
               al(i,j,k) = ar(i,j,k-1)
            end do
         end do
      end do
      !$omp end parallel do
   end if   

end subroutine ppm_find_face_vals_v2

!-----------------------------------------------------------------!
! steepen rho_b
! Note from YT: This subroutine is a major modification to 
!  the DAGH's version, and I believe there are bugs in 
!  the DAGH's version.
!-----------------------------------------------------------------!
!
subroutine ppm_steepen_rho_v2(ext,rho,d1rho,d2rho,rhor_ppm,rhol_ppm, &
     rhor_mc,rhol_mc,P,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th, & 
     neos, ergo_star, ergo_sigma, m,Symmetry)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: rho,d1rho,d2rho,P
 real*8, dimension(ext(1),ext(2),ext(3))        :: rhor_ppm,rhol_ppm
 real*8, dimension(ext(1),ext(2),ext(3))        :: rhor_mc,rhol_mc
 integer                                        :: m,neos,ergo_star,Symmetry
 real*8, dimension(neos) 			:: rho_tab,P_tab,eps_tab
 real*8, dimension(neos+1)			:: gamma_tab,k_tab
!
 integer                                        :: i,j,k,n
 integer                    :: imax,jmax,kmax,jmin,exit_do
 real*8, parameter          :: SIXTH = 0.1666666666666666666666667d0
 real*8                     :: K0, eta1, eta2, epsilon, Gamma, gamma_th,ergo_sigma
 real*8                     :: gamma_cold,P_cold,eps_cold
 real*8                     :: contact_disc_check,second_dir_check,relative_change_check
 real*8                     :: eta_tilde,eta
! standard Colella-Woodward parameters:
!    K0 = 0.1d0, eta1 = 20.d0, eta2 = 0.05d0, epsilon = 0.01d0
! our parameters:
 parameter(K0 = 0.1d0, eta1 = 20.d0, eta2 = 0.05d0, epsilon = 0.01d0)
! The following values are Font-like:
! parameter(K0 = 1.d0, eta1 = 5.d0, eta2 = 0.05d0, epsilon = 0.1d0)
!
 imax = ubound(rho,1)
 jmax = ubound(rho,2)
 kmax = ubound(rho,3)
 jmin = 1
 if (Symmetry ==4) then 
    jmin = 2
    jmax = 2
 end if

 if (m==1) then 
    !$omp parallel do
    do k=1,kmax
       do j=jmin,jmax
	  do i=2,imax-1
	     gamma_cold = gamma_tab(neos+1)
	     exit_do = 0
	     do n=1,neos
	        if (rho(i,j,k) .lt. rho_tab(n)) then 
                   gamma_cold = gamma_tab(n)
	           exit_do = 1
	        end if
	        if (exit_do == 1) exit
	     end do
	     call compute_pcold_epscold(rho(i,j,k), P_cold, eps_cold, &
                     neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
	     ! Compute effective Gamma = (partial P / partial rho0)_s /(P/rho0)
	     Gamma = gamma_th + (gamma_cold-gamma_th)*P_cold/P(i,j,k)
	     contact_disc_check = Gamma*K0*abs(rho(i+1,j,k)-rho(i-1,j,k))* & 
                         min(P(i+1,j,k),P(i-1,j,k)) & 
                        -abs(P(i+1,j,k)-P(i-1,j,k))*min(rho(i+1,j,k),rho(i-1,j,k))
	     second_dir_check = -d2rho(i+1,j,k)*d2rho(i-1,j,k)
             relative_change_check = abs(2.d0*d1rho(i,j,k)) &
                     - epsilon*min(rho(i+1,j,k),rho(i-1,j,k))

	     if(contact_disc_check.ge.0.d0 .and. second_dir_check.ge.0.d0 &
                  .and. relative_change_check.ge.0.d0) then
               if (abs(d1rho(i,j,k)) .gt. 0.d0) then 
	          eta_tilde = -SIXTH*(d2rho(i+1,j,k)-d2rho(i-1,j,k))/d1rho(i,j,k)
	       else 
	          eta_tilde = 0.d0
	       end if
	       eta = max(0.d0,min(eta1*(eta_tilde - eta2),1.d0))
               rhol_ppm(i,j,k) = rhol_ppm(i,j,k)*(1.d0-eta) &
                    + rhol_mc(i,j,k)*eta
               rhor_ppm(i,j,k) = rhor_ppm(i,j,k)*(1.d0-eta) &
                    + rhor_mc(i+1,j,k)*eta
	     end if
	  end do
       end do
    end do
    !$omp end parallel do  
 else if(m==2) then 
    !$omp parallel do
    do k=1,kmax
       do j=jmin+1,jmax-1
          do i=1,imax
             gamma_cold = gamma_tab(neos+1)
             exit_do = 0
             do n=1,neos
                if (rho(i,j,k) .lt. rho_tab(n)) then
                   gamma_cold = gamma_tab(n)
                   exit_do = 1
                end if
                if (exit_do == 1) exit
             end do
             call compute_pcold_epscold(rho(i,j,k), P_cold, eps_cold, &
                     neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
             ! Compute effective Gamma = (partial P / partial rho0)_s /(P/rho0)
             Gamma = gamma_th + (gamma_cold-gamma_th)*P_cold/P(i,j,k)
             contact_disc_check = Gamma*K0*abs(rho(i,j+1,k)-rho(i,j-1,k))* & 
			  min(P(i,j+1,k),P(i,j-1,k)) &
                        - abs(P(i,j+1,k)-P(i,j-1,k))*min(rho(i,j+1,k),rho(i,j-1,k))
             second_dir_check = -d2rho(i,j+1,k)*d2rho(i,j-1,k)
             relative_change_check = abs(2.d0*d1rho(i,j,k)) &
                     - epsilon*min(rho(i,j+1,k),rho(i,j-1,k))

             if(contact_disc_check.ge.0.d0 .and. second_dir_check.ge.0.d0 &
                  .and. relative_change_check.ge.0.d0) then
               if (abs(d1rho(i,j,k)) .gt. 0.d0) then
                  eta_tilde = -SIXTH*(d2rho(i,j+1,k)-d2rho(i,j-1,k))/d1rho(i,j,k)
               else
                  eta_tilde = 0.d0
               end if
               eta = max(0.d0,min(eta1*(eta_tilde - eta2),1.d0))
               rhol_ppm(i,j,k) = rhol_ppm(i,j,k)*(1.d0-eta) &
                    + rhol_mc(i,j,k)*eta
               rhor_ppm(i,j,k) = rhor_ppm(i,j,k)*(1.d0-eta) &
                    + rhor_mc(i,j+1,k)*eta
             end if
          end do
       end do
    end do
    !omp end parallel do
 else
    !$omp parallel do
    do k=2,kmax-1
       do j=jmin,jmax
          do i=1,imax
             gamma_cold = gamma_tab(neos+1)
             exit_do = 0
             do n=1,neos
                if (rho(i,j,k) .lt. rho_tab(n)) then
                   gamma_cold = gamma_tab(n)
                   exit_do = 1
                end if
                if (exit_do == 1) exit
             end do
             call compute_pcold_epscold(rho(i,j,k), P_cold, eps_cold, &
                     neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
             ! Compute effective Gamma = (partial P / partial rho0)_s / (P/rho0)
             Gamma = gamma_th + (gamma_cold-gamma_th)*P_cold/P(i,j,k)
             contact_disc_check = Gamma*K0*abs(rho(i,j,k+1)-rho(i,j,k-1))* & 
			  min(P(i,j,k+1),P(i,j,k-1)) & 
                        - abs(P(i,j,k+1)-P(i,j,k-1))*min(rho(i,j,k+1),rho(i,j,k-1))
             second_dir_check = -d2rho(i,j,k+1)*d2rho(i,j,k-1)
             relative_change_check = abs(2.d0*d1rho(i,j,k)) &
                     - epsilon*min(rho(i,j,k+1),rho(i,j,k-1))

             if(contact_disc_check.ge.0.d0 .and. second_dir_check.ge.0.d0 &
                  .and. relative_change_check.ge.0.d0) then
               if (abs(d1rho(i,j,k)) .gt. 0.d0) then
                  eta_tilde = -SIXTH*(d2rho(i,j,k+1)-d2rho(i,j,k-1))/d1rho(i,j,k)
               else
                  eta_tilde = 0.d0
               end if
               eta = max(0.d0,min(eta1*(eta_tilde - eta2),1.d0))
               rhol_ppm(i,j,k) = rhol_ppm(i,j,k)*(1.d0-eta) &
                    + rhol_mc(i,j,k)*eta
               rhor_ppm(i,j,k) = rhor_ppm(i,j,k)*(1.d0-eta) &
                    + rhor_mc(i,j,k+1)*eta
             end if
          end do
       end do
    end do
    !omp end parallel do
 end if

end subroutine ppm_steepen_rho_v2

!-----------------------------------------------------------------!
! ensure monotonicity of the interpolating polynomial (remove the PPM plus stuff)
!-----------------------------------------------------------------!
!
subroutine ppm_monotonize_v2(ext,a,ar,al,Symmetry)
 implicit none
 integer, dimension(3)                     :: ext
 real*8, dimension(ext(1),ext(2),ext(3))   :: a,ar,al
 real*8					   :: aijk,arijk,alijk,daijk,maijk
 integer                                   :: i,j,k
 integer                                   :: imax,jmin,jmax,kmax
 integer :: Symmetry
 real*8, parameter 			   :: sixth = 0.16666666666666666667d0
!
 imax = ubound(a,1)
 jmax = ubound(a,2)
 kmax = ubound(a,3)
 jmin = 1
 if (Symmetry == 4) then 
    jmin = 2
    jmax = 2
 end if

 !omp parallel do
 do k=1,kmax
    do j=jmin,jmax
       do i=1,imax
	  arijk=ar(i,j,k)
	  alijk=al(i,j,k)
	  aijk = a(i,j,k)
	  daijk = arijk - alijk
	  maijk = 0.5d0*(arijk+alijk)
	  if ( (arijk-aijk)*(aijk-alijk) .le. 0.d0) then 
	     ar(i,j,k) = aijk
	     al(i,j,k) = aijk
	  else if ( daijk*(aijk-maijk) .gt. sixth*daijk**2) then 
	     al(i,j,k) = 3.d0*aijk - 2.d0*arijk
	  else if ( daijk*(aijk-maijk) .lt. -sixth*daijk**2) then
	     ar(i,j,k) = 3.d0*aijk - 2.d0*alijk
	  end if
       end do
    end do
 end do
 !omp end parallel do

 end subroutine ppm_monotonize_v2

!-----------------------------------------------------------------!
! calculate \tilde{f}_j                                           !
!-----------------------------------------------------------------!
!
subroutine ppm_ftilde_v2(ext,nghostzones,ftilde,P,vx,vy,vz,symP,symvx,symvy,symvz, & 
                       m,Symmetry,cell_centering_enabled)
 implicit none
 integer, dimension(3)                          :: ext,nghostzones
 real*8, dimension(ext(1),ext(2),ext(3))        :: ftilde,P
 real*8, dimension(ext(1),ext(2),ext(3))        :: vx,vy,vz
 integer, dimension(3)				:: symP,symvx,symvy,symvz
 integer, intent(in)                            :: m,Symmetry,cell_centering_enabled
 integer                                        :: i,j,k, sym_index1,sym_index2
 integer                                        :: imax,jmin,jmax,kmax
 integer, parameter :: NO_SYMM=0, EQUATORIAL=1, OCTANT=2, PI_SYMM=3, AXISYM=4
 real*8, parameter                  :: omega1 = 0.75d0, omega2 = 10.d0
 real*8, parameter                  :: epsilon2 = 0.33d0
! The following parameter choices are Font-like.
! real*8, parameter                  :: omega1 = 0.52d0, omega2 = 10.d0
! real*8, parameter                  :: epsilon2 = 0.5d0
 real*8                             :: w, q1, dP1,dP2,dP1odP2,q2,symv

 imax = ubound(ftilde,1)
 jmax = ubound(ftilde,2)
 kmax = ubound(ftilde,3)
 jmin = 1
 if (Symmetry==AXISYM) then 
    jmin = 2
    jmax = 2
 end if

 if (cell_centering_enabled==1) then
    sym_index1 = 2*nghostzones(m) + 1
    ! In axisymmetry, nghostzones=2, though [# symmetry ghostzones] = 1.
    if (Symmetry==4 .and. m==1) sym_index1 = 3
 else
    sym_index1 = 2*nghostzones(m) + 2
 end if
 sym_index2 = sym_index1 + 1

 if (m==1) then 
    symv = symvx(1)
    ! First compute ftilde(i,j,k) for i=1 and i=2 
    if (Symmetry==NO_SYMM .or. Symmetry==EQUATORIAL .or. Symmetry==PI_SYMM) then 
       !$omp parallel do
       do k=1,kmax
	  do j=1,jmax
             ftilde(1,j,k) = 1.d0
             ftilde(2,j,k) = 1.d0
	  end do
       end do
       !$omp end parallel do
    else 
       !$omp parallel do private (dP1,dP2,dP1odP2,q1,q2,w)
       do k=1,kmax
          do j=jmin,jmax

             dP1 = P(2,j,k) - symP(m)*P(sym_index1,j,k)
	     dP2 = P(3,j,k) - symP(m)*P(sym_index2,j,k)
	     if (dP2 .ne. 0.d0) then 
	        dP1odP2 = dP1/dP2 
   	     else
	        dP1odP2 = 1.d0
             end if    
	     q1 = (dP1odP2-omega1)*omega2
	     q2 = abs(dP1) / min(P(2,j,k),symP(m)*P(sym_index1,j,k))
	     if (q2 .gt. epsilon2 .and. & 
	   	  (symv*vx(sym_index1,j,k)-vx(2,j,k)) .gt. 0.d0) then 
         ! "inside a shock":
	        w = 1.d0
	     else
         ! NOT "inside a shock":
	        w = 0.d0
	     end if
	     ftilde(1,j,k) = min(1.d0, w*max(0.d0,q1) )

	     dP1 = P(3,j,k) - P(1,j,k)
	     dP2 = P(4,j,k) - symP(m)*P(sym_index1,j,k)
             if (dP2 .ne. 0.d0) then
                dP1odP2 = dP1/dP2
             else
                dP1odP2 = 1.d0
             end if
             q1 = (dP1odP2-omega1)*omega2
	     q2 = abs(dP1)/min(P(3,j,k),P(1,j,k))
             if (q2 .gt. epsilon2 .and. &
                  (vx(1,j,k)-vx(3,j,k)) .gt. 0.d0) then
                w = 1.d0
             else
                w = 0.d0
             end if
             ftilde(2,j,k) = min(1.d0, w*max(0.d0,q1) )
	  end do
       end do
       !$omp end parallel do
     end if

     ! Next set ftilde(i,j,k) = 1 for i = imax and i = imax-1
     !$omp parallel do
     do k=1,kmax
        do j=jmin,jmax
           ftilde(imax-1,j,k) = 1.d0
	   ftilde(imax,j,k) =  1.d0
	end do
     end do
     !$omp end parallel do

     ! Now compute ftilde(i,j,k) for the rest of i
     !$omp parallel do private (dP1,dP2,dP1odP2,q1,q2,w)
     do k=1,kmax
        do j=jmin,jmax
           do i=3,imax-2
	      dP1 = P(i+1,j,k)-P(i-1,j,k)
	      dP2 = P(i+2,j,k)-P(i-2,j,k)
              if (dP2 .ne. 0.d0) then
                 dP1odP2 = dP1/dP2
              else
                 dP1odP2 = 1.d0
              end if
              q1 = (dP1odP2-omega1)*omega2
              q2 = abs(dP1)/min(P(i+1,j,k),P(i-1,j,k))
              if (q2 .gt. epsilon2 .and. &
                (vx(i-1,j,k)-vx(i+1,j,k)) .gt. 0.d0) then
                 ! "inside a shock":
                 w = 1.d0
              else
                 ! NOT "inside a shock":
                 w = 0.d0
              end if
              ftilde(i,j,k) = min(1.d0, w*max(0.d0,q1) )
	   end do
	end do
     end do
     !omp end parallel do
 else if (m==2) then 
    symv = symvy(2)
    ! First compute filde(i,j,k) for j=1 and j=2
    if (Symmetry==NO_SYMM .or. Symmetry==EQUATORIAL .or. Symmetry==PI_SYMM) then 
       !$omp parallel do
       do k=1,kmax
	  do i=1,imax
             ftilde(i,1,k) = 1.d0
             ftilde(i,2,k) = 1.d0
	  end do
       end do
       !$omp end parallel do
    else 
       !$omp parallel do private (dP1,dP2,dP1odP2,q1,q2,w)
       do k=1,kmax
          do i=1,imax

             dP1 = P(i,2,k) - symP(m)*P(i,sym_index1,k)
	     dP2 = P(i,3,k) - symP(m)*P(i,sym_index2,k)
	     if (dP2 .ne. 0.d0) then 
	        dP1odP2 = dP1/dP2 
   	     else
	        dP1odP2 = 1.d0
             end if    
	     q1 = (dP1odP2-omega1)*omega2
	     q2 = abs(dP1)/min(P(i,2,k),symP(m)*P(i,sym_index1,k))
	     if (q2 .gt. epsilon2 .and. & 
	   	  (symv*vy(i,sym_index1,k)-vy(i,2,k)) .gt. 0.d0) then 
	        w = 1.d0
	     else
	        w = 0.d0
	     end if
	     ftilde(i,1,k) = min(1.d0, w*max(0.d0,q1) )

	     dP1 = P(i,3,k) - P(i,1,k)
	     dP2 = P(i,4,k) - symP(m)*P(i,sym_index1,k)
             if (dP2 .ne. 0.d0) then
                dP1odP2 = dP1/dP2
             else
                dP1odP2 = 1.d0
             end if
             q1 = (dP1odP2-omega1)*omega2
	     q2 = abs(dP1)/min(P(i,3,k),P(i,1,k))
             if (q2 .gt. epsilon2 .and. &
                  (vy(i,1,k)-vy(i,3,k)) .gt. 0.d0) then
                 ! "inside a shock":
                w = 1.d0
             else
                 ! NOT "inside a shock":
                w = 0.d0
             end if
             ftilde(i,2,k) = min(1.d0, w*max(0.d0,q1) )
	  end do
       end do
       !$omp end parallel do
    end if

    ! Next set ftilde(i,j,k) = 1 for j = jmax and j = jmax-1
    !$omp parallel do
    do k=1,kmax
       do i=1,imax
          ftilde(i,jmax-1,k) = 1.d0
	  ftilde(i,jmax,k) =  1.d0
       end do
    end do
    !$omp end parallel do

    ! Now compute ftilde(i,j,k) for the rest of j
    !$omp parallel do private (dP1,dP2,dP1odP2,q1,q2,w)
    do k=1,kmax
       do j=3,jmax-2
          do i=1,imax
	     dP1 = P(i,j+1,k)-P(i,j-1,k)
	     dP2 = P(i,j+2,k)-P(i,j-2,k)
             if (dP2 .ne. 0.d0) then
                dP1odP2 = dP1/dP2
             else
                dP1odP2 = 1.d0
             end if
             q1 = (dP1odP2-omega1)*omega2
             q2 = abs(dP1)/min(P(i,j+1,k),P(i,j-1,k))
             if (q2 .gt. epsilon2 .and. &
                (vy(i,j-1,k)-vy(i,j+1,k)) .gt. 0.d0) then
                 ! "inside a shock":
                w = 1.d0
             else
                 ! NOT "inside a shock":
                w = 0.d0
             end if
             ftilde(i,j,k) = min(1.d0, w*max(0.d0,q1) )
           end do
	end do
     end do
     !$omp end parallel do
 else 
    symv = symvz(3)
    ! First compute filde(i,j,k) for k=1 and k=2
    if (Symmetry==NO_SYMM) then 
       !$omp parallel do
       do j=1,jmax
	  do i=1,imax
             ftilde(i,j,1) = 1.d0
             ftilde(i,j,2) = 1.d0
	  end do
       end do
       !$omp end parallel do
    else 
       !$omp parallel do private (dP1,dP2,dP1odP2,q1,q2,w)
       do j=jmin,jmax
          do i=1,imax

             dP1 = P(i,j,2) - symP(m)*P(i,j,sym_index1)
	     dP2 = P(i,j,3) - symP(m)*P(i,j,sym_index2)
	     if (dP2 .ne. 0.d0) then 
	        dP1odP2 = dP1/dP2 
   	     else
	        dP1odP2 = 1.d0
             end if    
	     q1 = (dP1odP2-omega1)*omega2
	     q2 = abs(dP1)/min(P(i,j,2),symP(m)*P(i,j,sym_index1))
	     if (q2 .gt. epsilon2 .and. & 
	   	  (symv*vz(i,j,sym_index1)-vz(i,j,2)) .gt. 0.d0) then 
                 ! "inside a shock":
	        w = 1.d0
	     else
                 ! NOT "inside a shock":
	        w = 0.d0
	     end if
	     ftilde(i,j,1) = min(1.d0, w*max(0.d0,q1) )

	     dP1 = P(i,j,3) - P(i,j,1)
	     dP2 = P(i,j,4) - symP(m)*P(i,j,sym_index1)
             if (dP2 .ne. 0.d0) then
                dP1odP2 = dP1/dP2
             else
                dP1odP2 = 1.d0
             end if
             q1 = (dP1odP2-omega1)*omega2
	     q2 = abs(dP1)/min(P(i,j,3),P(i,j,1))
             if (q2 .gt. epsilon2 .and. &
                  (vz(i,j,1)-vz(i,j,3)) .gt. 0.d0) then
                ! "inside a shock":
                w = 1.d0
             else
                ! NOT "inside a shock":
                w = 0.d0
             end if
             ftilde(i,j,2) = min(1.d0, w*max(0.d0,q1) )
	  end do
       end do
       !$omp end parallel do
    end if

    ! Next set ftilde(i,j,k) = 1 for k = kmax and k = kmax-1
    !$omp parallel do
    do j=jmin,jmax
       do i=1,imax
          ftilde(i,j,kmax-1) = 1.d0
	  ftilde(i,j,kmax) =  1.d0
       end do
    end do
    !$omp end parallel do

    ! Now compute ftilde(i,j,k) for the rest of k
    !$omp parallel do private (dP1,dP2,dP1odP2,q1,q2,w)
    do k=3,kmax-2
       do j=jmin,jmax
          do i=1,imax
	     dP1 = P(i,j,k+1)-P(i,j,k-1)
	     dP2 = P(i,j,k+2)-P(i,j,k-2)
             if (dP2 .ne. 0.d0) then
                dP1odP2 = dP1/dP2
             else
                dP1odP2 = 1.d0
             end if
             q1 = (dP1odP2-omega1)*omega2
             q2 = abs(dP1)/min(P(i,j,k+1),P(i,j,k-1))
             if (q2 .gt. epsilon2 .and. &
                (vz(i,j,k-1)-vz(i,j,k+1)) .gt. 0.d0) then
                ! "inside a shock":
                w = 1.d0
             else
                ! NOT "inside a shock":
                w = 0.d0
             end if
             ftilde(i,j,k) = min(1.d0, w*max(0.d0,q1) )
	  end do
        end do
    end do
    !$omp end parallel do
 end if

end subroutine ppm_ftilde_v2

!-----------------------------------------------------------------!
! apply flattening                                                !
!-----------------------------------------------------------------!
!
subroutine ppm_flatten_v2(ext,ftilde,rho_b,rho_br,rho_bl,P,Pr,Pl, &
     vx,vxr,vxl,vy,vyr,vyl,vz,vzr,vzl,Bx,Bxr,Bxl,By,Byr,Byl,Bz,Bzr,Bzl, & 
     em_evolve_enable,m,Symmetry)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: rho_b, rho_br, rho_bl
 real*8, dimension(ext(1),ext(2),ext(3))        :: ftilde,P,Pr,Pl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vx,vxr,vxl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vy,vyr,vyl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vz,vzr,vzl
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bx,Bxr,Bxl
 real*8, dimension(ext(1),ext(2),ext(3))        :: By,Byr,Byl
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bz,Bzr,Bzl
 integer, intent(in)                            :: m,Symmetry,em_evolve_enable
 integer                            :: i,j,k, imax,jmin,jmax,kmax
 integer                            :: s
 real*8				    :: f,omf
!
 imax = ubound(ftilde,1)
 jmax = ubound(ftilde,2)
 kmax = ubound(ftilde,3)
 jmin = 1
 if (Symmetry==4) then 
    jmin = 2
    jmax = 2
 end if

 if (m==1) then

    !omp parallel do
    do k=1,kmax
       do j=jmin,jmax
          do i=2,imax-1
             ! Calculate s
             if ((P(i+1,j,k) - P(i-1,j,k)) .le. 0.d0) then
                s = 1
             else
                s = -1
             end if

             ! Calculate f
             f = max(ftilde(i,j,k),ftilde(i+s,j,k))
  	     omf = 1.d0-f

             ! Flatten
             rho_br(i,j,k)   = rho_b(i,j,k)*f + rho_br(i,j,k)*omf
             rho_bl(i,j,k)   = rho_b(i,j,k)*f + rho_bl(i,j,k)*omf
             Pr(i,j,k)       = P(i,j,k)*f     + Pr(i,j,k)*omf
             Pl(i,j,k)       = P(i,j,k)*f     + Pl(i,j,k)*omf
             vxr(i,j,k)      = vx(i,j,k)*f    + vxr(i,j,k)*omf
             vxl(i,j,k)      = vx(i,j,k)*f    + vxl(i,j,k)*omf
             vyr(i,j,k)      = vy(i,j,k)*f    + vyr(i,j,k)*omf
             vyl(i,j,k)      = vy(i,j,k)*f    + vyl(i,j,k)*omf
             vzr(i,j,k)      = vz(i,j,k)*f    + vzr(i,j,k)*omf
             vzl(i,j,k)      = vz(i,j,k)*f    + vzl(i,j,k)*omf
	    
	     if (em_evolve_enable ==1) then 
                Bxr(i,j,k)      = Bx(i,j,k)*f    + Bxr(i,j,k)*omf
                Bxl(i,j,k)      = Bx(i,j,k)*f    + Bxl(i,j,k)*omf
                Byr(i,j,k)      = By(i,j,k)*f    + Byr(i,j,k)*omf
                Byl(i,j,k)      = By(i,j,k)*f    + Byl(i,j,k)*omf
                Bzr(i,j,k)      = Bz(i,j,k)*f    + Bzr(i,j,k)*omf
                Bzl(i,j,k)      = Bz(i,j,k)*f    + Bzl(i,j,k)*omf
	     end if

          end do
       end do
    end do
    !omp end parallel do

 else if (m==2) then

    !omp parallel do
    do k=1,kmax
       do j=jmin+1,jmax-1
          do i=1,imax
             ! Calculate s
             if ((P(i,j+1,k) - P(i,j-1,k)) .le. 0.d0) then
                s = 1
             else
                s = -1
             end if

             ! Calculate f
             f = max(ftilde(i,j,k),ftilde(i,j+s,k))
  	     omf = 1.d0-f

             ! Flatten
             rho_br(i,j,k)   = rho_b(i,j,k)*f + rho_br(i,j,k)*omf
             rho_bl(i,j,k)   = rho_b(i,j,k)*f + rho_bl(i,j,k)*omf
             Pr(i,j,k)       = P(i,j,k)*f     + Pr(i,j,k)*omf
             Pl(i,j,k)       = P(i,j,k)*f     + Pl(i,j,k)*omf
             vxr(i,j,k)      = vx(i,j,k)*f    + vxr(i,j,k)*omf
             vxl(i,j,k)      = vx(i,j,k)*f    + vxl(i,j,k)*omf
             vyr(i,j,k)      = vy(i,j,k)*f    + vyr(i,j,k)*omf
             vyl(i,j,k)      = vy(i,j,k)*f    + vyl(i,j,k)*omf
             vzr(i,j,k)      = vz(i,j,k)*f    + vzr(i,j,k)*omf
             vzl(i,j,k)      = vz(i,j,k)*f    + vzl(i,j,k)*omf

             if (em_evolve_enable ==1) then
                Bxr(i,j,k)      = Bx(i,j,k)*f    + Bxr(i,j,k)*omf
                Bxl(i,j,k)      = Bx(i,j,k)*f    + Bxl(i,j,k)*omf
                Byr(i,j,k)      = By(i,j,k)*f    + Byr(i,j,k)*omf
                Byl(i,j,k)      = By(i,j,k)*f    + Byl(i,j,k)*omf
                Bzr(i,j,k)      = Bz(i,j,k)*f    + Bzr(i,j,k)*omf
                Bzl(i,j,k)      = Bz(i,j,k)*f    + Bzl(i,j,k)*omf
             end if

          end do
       end do
    end do
    !omp end parallel do

 else

    !omp parallel do
    do k=2,kmax-1
       do j=jmin,jmax
          do i=1,imax
             ! Calculate s
             if ((P(i,j,k+1) - P(i,j,k-1)) .le. 0.d0) then
                s = 1
             else
                s = -1
             end if

             ! Calculate f
             f = max(ftilde(i,j,k),ftilde(i,j,k+s))
             omf = 1.d0-f

             ! Flatten
             rho_br(i,j,k)   = rho_b(i,j,k)*f + rho_br(i,j,k)*omf
             rho_bl(i,j,k)   = rho_b(i,j,k)*f + rho_bl(i,j,k)*omf
             Pr(i,j,k)       = P(i,j,k)*f     + Pr(i,j,k)*omf
             Pl(i,j,k)       = P(i,j,k)*f     + Pl(i,j,k)*omf
             vxr(i,j,k)      = vx(i,j,k)*f    + vxr(i,j,k)*omf
             vxl(i,j,k)      = vx(i,j,k)*f    + vxl(i,j,k)*omf
             vyr(i,j,k)      = vy(i,j,k)*f    + vyr(i,j,k)*omf
             vyl(i,j,k)      = vy(i,j,k)*f    + vyl(i,j,k)*omf
             vzr(i,j,k)      = vz(i,j,k)*f    + vzr(i,j,k)*omf
             vzl(i,j,k)      = vz(i,j,k)*f    + vzl(i,j,k)*omf

             if (em_evolve_enable ==1) then
                Bxr(i,j,k)      = Bx(i,j,k)*f    + Bxr(i,j,k)*omf
                Bxl(i,j,k)      = Bx(i,j,k)*f    + Bxl(i,j,k)*omf
                Byr(i,j,k)      = By(i,j,k)*f    + Byr(i,j,k)*omf
                Byl(i,j,k)      = By(i,j,k)*f    + Byl(i,j,k)*omf
                Bzr(i,j,k)      = Bz(i,j,k)*f    + Bzr(i,j,k)*omf
                Bzl(i,j,k)      = Bz(i,j,k)*f    + Bzl(i,j,k)*omf
             end if

          end do
       end do
    end do
    !omp end parallel do

 end if

end subroutine ppm_flatten_v2

subroutine ppm_shift_v2(ext,ar,al,temp,m,Symmetry)
  implicit none
   integer, dimension(3)                        :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: ar,al,temp
 integer                                        :: m, Symmetry
 integer                                        :: i,j,k,imax,jmax,kmax,jmin
!
 imax = ubound(ar,1)
 jmax = ubound(ar,2)
 kmax = ubound(ar,3)
 jmin = 1
 if (Symmetry ==4) then 
    jmin = 2
    jmax = 2
 end if

 !omp parallel do
 do k=1,kmax
    do j=jmin,jmax
       do i=1,imax
	  temp(i,j,k) = ar(i,j,k)
       end do
    end do
 end do
 !omp end parallel do

 !omp parallel do
 do k=1,kmax
    do j=jmin,jmax
       do i=1,imax
          ar(i,j,k) = al(i,j,k)
       end do
    end do
 end do
 !omp end parallel do

 if (m==1) then 
    !omp parallel do
    do k=1,kmax
       do j=jmin,jmax
	  do i=2,imax
	     al(i,j,k) = temp(i-1,j,k)
	  end do
       end do
    end do   
    !omp end parallel do
 else if (m==2) then 
    !omp parallel do
    do k=1,kmax
       do j=jmin+1,jmax
          do i=1,imax
             al(i,j,k) = temp(i,j-1,k)
          end do
       end do
    end do
    !omp end parallel do
 else
    !omp parallel do
    do k=2,kmax
       do j=jmin,jmax
          do i=1,imax
             al(i,j,k) = temp(i,j,k-1)
          end do
       end do
    end do
    !omp end parallel do
 end if
  
end subroutine ppm_shift_v2

! Use PPM scheme to reconstruct a smooth function 
!
subroutine compute_face_avg_ppm(m,ext,nghostzones,cell_centering_enabled,Symmetry,sym,f,ff)
    implicit none
    integer, dimension(3) :: ext,sym,nghostzones
    integer :: m,Symmetry,cell_centering_enabled
    integer :: sym_index1,sym_index2
    real*8, dimension(ext(1),ext(2),ext(3)) :: f,ff
    integer :: imax,jmax,kmax,jmin,i,j,k
    real*8, parameter :: am2=-1.d0/12.d0, am1=7.d0/12.d0, a0=7.d0/12.d0, a1=-1.d0/12.d0
    real*8, parameter :: bm1=0.25d0, b0=13.d0/12d0, b1=-5.d0/12.d0, b2=1.d0/12.d0
    real*8, parameter :: c0=25.d0/12.d0, c1=-23.d0/12.d0, c2=13.d0/12.d0, c3=-0.25d0
    real*8, parameter :: dm3=1.d0/12.d0, dm2=-5.d0/12.d0, dm1=13.d0/12.d0, d0=0.25d0
    integer, parameter :: NO_SYM=0,EQUATORIAL=1,OCTANT=2,PI_SYMM=3,AXISYM=4
!
    imax = ubound(f,1)
    jmax = ubound(f,2)
    kmax = ubound(f,3)
    jmin = 1
    if (Symmetry==AXISYM) then
       jmin = 2
       jmax = 2
    end if

    if (cell_centering_enabled==1) then
       sym_index1 = 2*nghostzones(m) + 1
       ! In axisymmetry, nghostzones=2, though [# symmetry ghostzones] = 1.
       if (Symmetry==4 .and. m==1) sym_index1 = 3
    else
       sym_index1 = 2*nghostzones(m) + 2
    end if
    sym_index2 = sym_index1 + 1

    if (m==1) then 

       ! i = 1
       if (Symmetry==AXISYM .or. Symmetry==OCTANT) then 
	  !$omp parallel do
	  do k=1,kmax
	     do j=jmin,jmax
	        ff(1,j,k) = am2*sym(m)*f(sym_index2,j,k)+am1*sym(m)*f(sym_index1,j,k) + a0*f(1,j,k) + a1*f(2,j,k)
	     end do
	  end do
	  !$omp end parallel do  
       else
	  !$omp parallel do
	  do k=1,kmax
             do j=jmin,jmax
                ff(1,j,k) = c0*f(1,j,k) + c1*f(2,j,k) + c2*f(3,j,k) + c3*f(4,j,k)
             end do
          end do
          !$omp end parallel do
       end if

       ! i = 2
       if (Symmetry==AXISYM .or. Symmetry==OCTANT) then
	  !$omp parallel do
          do k=1,kmax
             do j=jmin,jmax
	        ff(2,j,k) = am2*sym(m)*f(sym_index1,j,k) + am1*f(1,j,k) + a0*f(2,j,k) + a1*f(3,j,k)
	     end do
          end do
          !$omp end parallel do
       else 
          !$omp parallel do
          do k=1,kmax
             do j=jmin,jmax
	        ff(2,j,k) = bm1*f(1,j,k) + b0*f(2,j,k) + b1*f(3,j,k) + b2*f(4,j,k)
             end do
          end do
          !$omp end parallel do
       end if
      
       ! i = 3 -- imax-1
       !$omp parallel do
       do k=1,kmax
	  do j=jmin,jmax
	     do i=3,imax-1
	        ff(i,j,k) = am2*f(i-2,j,k) + am1*f(i-1,j,k) + a0*f(i,j,k) + a1*f(i+1,j,k)
	     end do
	  end do
       end do
       !$omp end parallel do

       !i = imax
       !$omp parallel do
       do k=1,kmax
          do j=jmin,jmax
	     ff(imax,j,k) = dm3*f(imax-3,j,k) + dm2*f(imax-2,j,k) + dm1*f(imax-1,j,k) + d0*f(imax,j,k)
          end do
       end do
       !$omp end parallel do

    else if (m==2) then 

       ! j = 1
       if (Symmetry==OCTANT) then 
	  !$omp parallel do
	  do k=1,kmax
	     do i=1,imax
	        ff(i,1,k) = am2*sym(m)*f(i,sym_index2,k)+am1*sym(m)*f(i,sym_index1,k) + a0*f(i,1,k) + a1*f(i,2,k)
	     end do
	  end do
	  !$omp end parallel do  
       else
	  !$omp parallel do
	  do k=1,kmax
             do i=1,imax
                ff(i,1,k) = c0*f(i,1,k) + c1*f(i,2,k) + c2*f(i,3,k) + c3*f(i,4,k)
             end do
          end do
          !$omp end parallel do
       end if

       ! j = 2
       if (Symmetry==OCTANT) then
	  !$omp parallel do
          do k=1,kmax
             do i=1,imax
	        ff(i,2,k) = am2*sym(m)*f(i,sym_index1,k) + am1*f(i,1,k) + a0*f(i,2,k) + a1*f(i,3,k)
	     end do
          end do
          !$omp end parallel do
       else 
          !$omp parallel do
          do k=1,kmax
             do i=1,imax
	        ff(i,2,k) = bm1*f(i,1,k) + b0*f(i,2,k) + b1*f(i,3,k) + b2*f(i,4,k)
             end do
          end do
          !$omp end parallel do
       end if
      
       ! j = 3 -- jmax-1
       !$omp parallel do
       do k=1,kmax
	  do j=3,jmax-1
	     do i=1,imax
	        ff(i,j,k) = am2*f(i,j-2,k) + am1*f(i,j-1,k) + a0*f(i,j,k) + a1*f(i,j+1,k)
	     end do
	  end do
       end do
       !$omp end parallel do

       !j = jmax
       !$omp parallel do
       do k=1,kmax
          do i=1,imax
	     ff(i,jmax,k) = dm3*f(i,jmax-3,k) + dm2*f(i,jmax-2,k) + dm1*f(i,jmax-1,k) + d0*f(i,jmax,k)
          end do
       end do
       !$omp end parallel do

    else

       ! k = 1
       if (Symmetry .ne. NO_SYM) then 
	  !$omp parallel do
	  do j=jmin,jmax
	     do i=1,imax
	        ff(i,j,1) = am2*sym(m)*f(i,j,sym_index2)+am1*sym(m)*f(i,j,sym_index1) + a0*f(i,j,1) + a1*f(i,j,2)
	     end do
	  end do
	  !$omp end parallel do  
       else
	  !$omp parallel do
	  do j=jmin,jmax
             do i=1,imax
                ff(i,j,1) = c0*f(i,j,1) + c1*f(i,j,2) + c2*f(i,j,3) + c3*f(i,j,4)
             end do
          end do
          !$omp end parallel do
       end if

       ! k = 2
       if (Symmetry .ne. NO_SYM) then
	  !$omp parallel do
          do j=jmin,jmax
             do i=1,imax
	        ff(i,j,2) = am2*sym(m)*f(i,j,sym_index1) + am1*f(i,j,1) + a0*f(i,j,2) + a1*f(i,j,3)
	     end do
          end do
          !$omp end parallel do
       else 
          !$omp parallel do
          do j=jmin,jmax
             do i=1,imax
	        ff(i,j,2) = bm1*f(i,j,1) + b0*f(i,j,2) + b1*f(i,j,3) + b2*f(i,j,4)
             end do
          end do
          !$omp end parallel do
       end if
      
       ! k = 3 -- kmax-1
       !$omp parallel do
       do k=3,kmax-1
	  do j=jmin,jmax
	     do i=1,imax
	        ff(i,j,k) = am2*f(i,j,k-2) + am1*f(i,j,k-1) + a0*f(i,j,k) + a1*f(i,j,k+1)
	     end do
	  end do
       end do
       !$omp end parallel do

       !k = kmax
       !$omp parallel do
       do j=jmin,jmax
          do i=1,imax
	     ff(i,j,kmax) = dm3*f(i,j,kmax-3) + dm2*f(i,j,kmax-2) + dm1*f(i,j,kmax-1) + d0*f(i,j,kmax)
          end do
       end do
       !$omp end parallel do

    end if
end subroutine compute_face_avg_ppm
