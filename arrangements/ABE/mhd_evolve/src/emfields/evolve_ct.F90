!---------------------------------------------------------------------!
! fij = v^i B^j - v^j B^i					      !
!---------------------------------------------------------------------!
!
subroutine flux_induction(ext,X,fij,vir,vjr,Bir,Bjr,vil,vjl,Bil,Bjl, & 
     cmax,cmin,pow_axi,Symmetry,Fr,Fl, X_f)
  implicit none
  integer, dimension(3)                          :: ext
  real*8, dimension(ext(1),ext(2),ext(3))	:: X
  real*8, dimension(ext(1),ext(2),ext(3))        :: vir,vjr,Bir,Bjr,fij
  real*8, dimension(ext(1),ext(2),ext(3))        :: vil,vjl,Bil,Bjl,v02
  real*8, dimension(ext(1),ext(2),ext(3))        :: cmin,cmax
  integer					:: pow_axi,Symmetry
  integer                                        :: imin,imax,jmin,jmax
  integer                                        :: kmin,kmax, i,j,k
  real*8						:: dX2
  integer, parameter 				:: AXISYM = 4
  ! Auxiliary arrays
  real*8, dimension(ext(1),ext(2),ext(3))        :: Fr,Fl, X_f
  !
  Fr = vir*Bjr - vjr*Bir
  Fl = vil*Bjl - vjl*Bil
  call flux_hll(ext, Bjr,Bjl, Fr,Fl, fij, cmax,cmin)

  if (Symmetry==AXISYM .and. pow_axi .ne. 0) then
     imin = lbound(fij,1)
     jmin = lbound(fij,2)
     kmin = lbound(fij,3)
     imax = ubound(fij,1)
     jmax = ubound(fij,2)
     kmax = ubound(fij,3)
     dX2 = 0.5d0*(X(imin+1,1,1)-X(imin,1,1))
     do k=kmin,kmax
        do j=jmin,jmax
           do i=imin,imax
              if (pow_axi == 1) then
                 X_f(i,j,k) = X(i,j,k) - dX2
              else if (pow_axi == 2) then
                 !pow_axi == 2 for fzx!
                 X_f(i,j,k) = X(i,j,k)
              end if
           end do
        end do
     end do
     fij = fij*X_f
  end if

end subroutine flux_induction

!-----------------------------------------------------------!
! Constraint Transport: f^{ij} -> \tilde{f}^{ij}	    !
!-----------------------------------------------------------!
!
subroutine flux_ct(ex,X,Y,Z,fxy,fxz,fyx,fyz,fzx,fzy, ftxy,ftxz,ftyx,ftyz, &
     ftzx,ftzy,Sym_Bz,Symmetry)
  implicit none
  integer, dimension(3)				:: ex
  real*8, dimension(ex(1),ex(2),ex(3))		:: fxy,fxz,fyx,fyz,fzx,fzy
  real*8, dimension(ex(1),ex(2),ex(3))           :: ftxy,ftxz,ftyx,ftyz,ftzx,ftzy
  real*8, dimension(ex(1),ex(2),ex(3))      	:: X,Y,Z
  integer                                        :: imin,imax,jmin,jmax,kmin,kmax
  integer                                        :: i,j,k,Symmetry
  real*8						:: Sym_Bz,Symx, Symz
  integer                                        :: NO_SYMM, EQUATORIAL
  integer                                        :: PI_SYMM, AXISYM, OCTANT
  PARAMETER (NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  !
  ! Input translation
  !
  imin = lbound(fxy,1)
  jmin = lbound(fxy,2)
  kmin = lbound(fxy,3)
  imax = ubound(fxy,1)
  jmax = ubound(fxy,2)
  kmax = ubound(fxy,3)

  if(Symmetry == OCTANT) then 
     if(X(1,1,1).lt.0.D0) then
        imin = imin + 1
     end if
     if(Y(1,1,1).lt.0.D0) then
        jmin = jmin + 1
     end if
     if(Z(1,1,1).lt.0.D0) then
        kmin = kmin + 1
     end if
  else if(Symmetry == AXISYM) then
     if(X(1,1,1).lt.0.D0) then
        imin = imin + 1
     end if
     if(Z(1,1,1).lt.0.D0) then
        kmin = kmin + 1
     end if
  else if(Symmetry == EQUATORIAL) then
     if(Z(1,1,1).lt.0.D0) then
        kmin = kmin + 1
     end if
  end if

  if (Symmetry==NO_SYMM .or. Z(1,1,kmin) .lt. 0.d0) then
     Symz = 1.d0
  else
     Symz = Sym_Bz
     !     Symz = -Sym_Bz
  end if

  !
  ! Calculate \tilde{f}^{xy}(i-1/2,j,k) and store in ftxy(i,j,k),
  !	    \tilde{f}^{yx}(i,j-1/2,k) and store in ftyx(i,j,k)
  !

  ! Testing:  Initialize everything to zero first:
  ftxy = 0.d0
  ftxz = 0.d0
  ftyx = 0.d0
  ftyz = 0.d0
  ftzx = 0.d0 
  ftzy = 0.d0

  ftxy = fxy
  ftxz = fxz
  ftyx = fyx
  ftyz = fyz
  ftzx = fzx
  ftzy = fzy

  !Enforce fzx, fzy = 0 on equatorial plane!
  if (Sym_Bz .lt. 0.d0) then
     fzx(:,:,kmin) = 0.d0
     if (Symmetry .ne. AXISYM) fzy(:,:,kmin) = 0.d0
  end if


  if (Symmetry .ne. AXISYM) then
     do k=kmin,kmax
	do j=jmin,jmax
	   do i=imin,imax
	      ftxy(i,j,k)=2.d0*fxy(i,j,k)-fyx(i,j,k)
	      if (j .ne. jmax) then
		 ftxy(i,j,k)=ftxy(i,j,k)+fxy(i,j+1,k)-fyx(i,j+1,k)
	      else
		 ftxy(i,j,k)=ftxy(i,j,k)+fxy(i,j,k)-fyx(i,j,k)
	      end if
	      if (j==jmin .and. Symmetry==OCTANT) then
	         ftxy(i,j,k)=ftxy(i,j,k)-fxy(i,j,k)
	      elseif (j==jmin) then
		 ftxy(i,j,k)=ftxy(i,j,k)+fxy(i,j,k)
	      elseif (j .ne. jmin) then 
		 ftxy(i,j,k)=ftxy(i,j,k)+fxy(i,j-1,k)
	      end if
	      if (i==imin .and. Symmetry==OCTANT) then
	         ftxy(i,j,k)=ftxy(i,j,k)+fyx(i,j,k) 
	    	 if (j .ne. jmax) then 
		    ftxy(i,j,k)=ftxy(i,j,k)+fyx(i,j+1,k)
	         else
		    ftxy(i,j,k)=ftxy(i,j,k)+fyx(i,j,k)
		 end if
	      elseif (i==imin) then
	         ftxy(i,j,k)=ftxy(i,j,k)-fyx(i,j,k)
		 if (j .ne. jmax) then
                    ftxy(i,j,k)=ftxy(i,j,k)-fyx(i,j+1,k)
                 else
                    ftxy(i,j,k)=ftxy(i,j,k)-fyx(i,j,k)
                 end if
	      elseif (i .ne. imin) then
		 ftxy(i,j,k)=ftxy(i,j,k)-fyx(i-1,j,k)
	         if (j .ne. jmax) then 
		    ftxy(i,j,k)=ftxy(i,j,k)-fyx(i-1,j+1,k)
	         else
		    ftxy(i,j,k)=ftxy(i,j,k)-fyx(i-1,j,k)
		 end if
	      end if

	      ftyx(i,j,k)=2.d0*fyx(i,j,k)-fxy(i,j,k)
	      if (i .ne. imax) then
		 ftyx(i,j,k)=ftyx(i,j,k)+fyx(i+1,j,k)-fxy(i+1,j,k)
	      else
		 ftyx(i,j,k)=ftyx(i,j,k)+fyx(i,j,k)-fxy(i,j,k)
	      end if
	      if (i==imin .and. Symmetry==OCTANT) then
		 ftyx(i,j,k)=ftyx(i,j,k)-fyx(i,j,k)
	      elseif (i==imin) then
		 ftyx(i,j,k)=ftyx(i,j,k)+fyx(i,j,k)
	      elseif (i .ne. imin) then
		 ftyx(i,j,k)=ftyx(i,j,k)+fyx(i-1,j,k)
              end if
	      if (j==jmin .and. Symmetry==OCTANT) then
                 ftyx(i,j,k)=ftyx(i,j,k)+fxy(i,j,k)
	      elseif (j==jmin) then
                 ftyx(i,j,k)=ftyx(i,j,k)-fxy(i,j,k)
	      elseif (j .ne. jmin) then
                 ftyx(i,j,k)=ftyx(i,j,k)-fxy(i,j-1,k)
	      end if
	      if (i .ne. imax) then
		 if (j==jmin .and. Symmetry==OCTANT) then
		    ftyx(i,j,k)=ftyx(i,j,k)+fxy(i+1,j,k)
		 elseif (j .ne. jmin) then
		    ftyx(i,j,k)=ftyx(i,j,k)-fxy(i+1,j-1,k)
		 end if
	      else
	         if (j==jmin .and. Symmetry==OCTANT) then
                    ftyx(i,j,k)=ftyx(i,j,k)+fxy(i,j,k)
                 elseif (j .ne. jmin) then
                    ftyx(i,j,k)=ftyx(i,j,k)-fxy(i,j-1,k)
                 end if
	      end if
	   end do
	end do
     end do
     ftxy=0.125d0*ftxy
     ftyx=0.125d0*ftyx
  else
     ftxy=fxy
  end if
  !
  ! Calculate \tilde{f}^{yz}(i,j-1/2,k) and store in ftyz(i,j,k),
  !           \tilde{f}^{zy}(i,j,k-1/2) and store in ftzy(i,j,k)
  !
  if (Symmetry .ne. AXISYM) then
     do k=kmin,kmax
	do j=jmin,jmax
	   do i=imin,imax
	      ftyz(i,j,k)=2.d0*fyz(i,j,k)-fzy(i,j,k)
	      if (k .ne. kmax) then 
		 ftyz(i,j,k)=ftyz(i,j,k)+fyz(i,j,k+1)-fzy(i,j,k+1)
		 if (j==jmin .and. Symmetry==OCTANT) then
		    ftyz(i,j,k)=ftyz(i,j,k)+fzy(i,j,k+1)
		 elseif (j==jmin) then
		    ftyz(i,j,k)=ftyz(i,j,k)-fzy(i,j,k+1)
	         elseif (j .ne. jmin) then
		    ftyz(i,j,k)=ftyz(i,j,k)-fzy(i,j-1,k+1)
		 end if
	      else    ! k==kmax
		 ftyz(i,j,k)=ftyz(i,j,k)+fyz(i,j,k)-fzy(i,j,k)
                 if (j==jmin .and. Symmetry==OCTANT) then
                    ftyz(i,j,k)=ftyz(i,j,k)+fzy(i,j,k)
		 elseif (j==jmin) then
                    ftyz(i,j,k)=ftyz(i,j,k)-fzy(i,j,k)
                 elseif (j .ne. jmin) then
                    ftyz(i,j,k)=ftyz(i,j,k)-fzy(i,j-1,k)
                 end if
	      end if
	      if (k==kmin) then
	         ftyz(i,j,k)=ftyz(i,j,k)+Symz*fyz(i,j,k)
	      else
	         ftyz(i,j,k)=ftyz(i,j,k)+fyz(i,j,k-1)
	      end if
	      if (j==jmin .and. Symmetry==OCTANT) then
		 ftyz(i,j,k)=ftyz(i,j,k)+fzy(i,j,k)
	      elseif (j==jmin) then
		 ftyz(i,j,k)=ftyz(i,j,k)-fzy(i,j,k)
	      elseif (j .ne. jmin) then
		 ftyz(i,j,k)=ftyz(i,j,k)-fzy(i,j-1,k)
	      end if

	      ftzy(i,j,k)=2.d0*fzy(i,j,k)-fyz(i,j,k)
	      if (j .ne. jmax) then 
		 ftzy(i,j,k)=ftzy(i,j,k)+fzy(i,j+1,k)-fyz(i,j+1,k)
	      else
                 ftzy(i,j,k)=ftzy(i,j,k)+fzy(i,j,k)-fyz(i,j,k)
	      end if
	      if (j==jmin .and. Symmetry==OCTANT) then
		 ftzy(i,j,k)=ftzy(i,j,k)-fzy(i,j,k)
	      elseif (j==jmin) then
		 ftzy(i,j,k)=ftzy(i,j,k)+fzy(i,j,k)
	      elseif (j .ne. jmin) then
		 ftzy(i,j,k)=ftzy(i,j,k)+fzy(i,j-1,k)
	      end if
	      if (k==kmin) then
		 ftzy(i,j,k)=ftzy(i,j,k)-fyz(i,j,k)*Symz
		 if (j .ne. jmax)  then
		    ftzy(i,j,k)=ftzy(i,j,k)-fyz(i,j+1,k)*Symz
	         else
                    ftzy(i,j,k)=ftzy(i,j,k)-fyz(i,j,k)*Symz
		 end if
	      else
		 ftzy(i,j,k)=ftzy(i,j,k)-fyz(i,j,k-1)
		 if (j .ne. jmax) then 
		    ftzy(i,j,k)=ftzy(i,j,k)-fyz(i,j+1,k-1)
		 else
                    ftzy(i,j,k)=ftzy(i,j,k)-fyz(i,j,k-1)
		 end if
	      end if
	   end do
	end do
     end do
     ftyz=0.125d0*ftyz
     ftzy=0.125d0*ftzy
  else
     ftzy=fzy
  end if
  !
  ! Calculate \tilde{f}^{xz}(i-1/2,j,k) and store in ftxz(i,j,k),
  !           \tilde{f}^{zx}(i,j,k-1/2) and store in ftzx(i,j,k)
  !
  if (Symmetry==OCTANT) then
     Symx=-1.d0
  else
     Symx=1.d0
  end if
  do k=kmin,kmax
     do j=jmin,jmax
	do i=imin,imax
	   ftxz(i,j,k)=2.d0*fxz(i,j,k)-fzx(i,j,k)
	   if (k .ne. kmax) then
	      ftxz(i,j,k)=ftxz(i,j,k)+fxz(i,j,k+1)-fzx(i,j,k+1)
	      if (i==imin) then
		 ftxz(i,j,k)=ftxz(i,j,k)-Symx*fzx(i,j,k+1)
	      else
	    	 ftxz(i,j,k)=ftxz(i,j,k)-fzx(i-1,j,k+1)
	      end if
	   else  ! k==kmax
	      ftxz(i,j,k)=ftxz(i,j,k)+fxz(i,j,k)-fzx(i,j,k)
              if (i==imin) then
                 ftxz(i,j,k)=ftxz(i,j,k)-Symx*fzx(i,j,k)
              else
                 ftxz(i,j,k)=ftxz(i,j,k)-fzx(i-1,j,k)
              end if
	   end if
	   if (k==kmin) then
	      ftxz(i,j,k)=ftxz(i,j,k)+Symz*fxz(i,j,k)
	   else  
	      ftxz(i,j,k)=ftxz(i,j,k)+fxz(i,j,k-1)
	   end if
	   if (i==imin) then
	      ftxz(i,j,k)=ftxz(i,j,k)-Symx*fzx(i,j,k)
	   else
	      ftxz(i,j,k)=ftxz(i,j,k)-fzx(i-1,j,k)
	   end if

	   ftzx(i,j,k)=2.d0*fzx(i,j,k)-fxz(i,j,k)
	   if (i .ne. imax) then 
	      ftzx(i,j,k)=ftzx(i,j,k)+fzx(i+1,j,k)-fxz(i+1,j,k)
              if (k==kmin) then
	         ftzx(i,j,k)=ftzx(i,j,k)-Symz*fxz(i+1,j,k)
	      else
	         ftzx(i,j,k)=ftzx(i,j,k)-fxz(i+1,j,k-1)
	      end if
	   else  ! i==imax
	      ftzx(i,j,k)=ftzx(i,j,k)+fzx(i,j,k)-fxz(i,j,k)
              if (k==kmin) then
                 ftzx(i,j,k)=ftzx(i,j,k)-Symz*fxz(i,j,k)
              else
                 ftzx(i,j,k)=ftzx(i,j,k)-fxz(i,j,k-1)
              end if
	   end if
	   if (i==imin) then
	      ftzx(i,j,k)=ftzx(i,j,k)+Symx*fzx(i,j,k)
	   else
	      ftzx(i,j,k)=ftzx(i,j,k)+fzx(i-1,j,k)
	   end if
	   if (k==kmin) then
	      ftzx(i,j,k)=ftzx(i,j,k)-Symz*fxz(i,j,k)
	   else
	      ftzx(i,j,k)=ftzx(i,j,k)-fxz(i,j,k-1)
	   end if
	end do
     end do
  end do
  ftxz=0.125d0*ftxz
  ftzx=0.125d0*ftzx
  !ftzx(2,2,2),ftzy(2,2,2)
  !write(*,*) "fluxct",ftxy(2,2,2),ftxz(2,2,2),ftyx(2,2,2),ftyz(2,2,2),ftzx(2,2,2),ftzy(2,2,2),fxy(2,2,2),fxz(2,2,2),fyx(2,2,2),fyz(2,2,2),fzx(2,2,2),fzy(2,2,2)


  ! *** TEST ***
  !  if (Symmetry==AXISYM) then 
  !     ftxz(1,:,:) = 0.d0
  !     ftzx(1,:,:) = 0.d0
  !  end if
  ! ************

end subroutine flux_ct

!------------------------------------------------------!
! Advect B^i					       !
!------------------------------------------------------!
!
subroutine advect_b_ct(ex,X,Y,Z, Bx_rhs, By_rhs, Bz_rhs,  &
     fxy,fxz,fyx,fyz,fzx,fzy, Symmetry)
  implicit none
  integer, dimension(3)				:: ex
  real*8, dimension(ex(1),ex(2),ex(3))		:: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))		:: fxy,fxz,fyx,fyz,fzx,fzy
  real*8, dimension(ex(1),ex(2),ex(3))		:: Bx_rhs, By_rhs, Bz_rhs
  integer                                        :: imin,imax,jmin,jmax,kmin,kmax
  integer                                        :: i,j,k,Symmetry
  integer                                        :: glob_imin,glob_kmin,Nx,Nz
  real*8						:: ddx,ddy,ddz
  integer                                        :: NO_SYMM, EQUATORIAL
  integer                                        :: PI_SYMM, AXISYM, OCTANT
  PARAMETER (NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  !
  ! Input translation
  !
  imin = lbound(fxy,1)
  jmin = lbound(fxy,2)
  kmin = lbound(fxy,3)
  imax = ubound(fxy,1)
  jmax = ubound(fxy,2)
  kmax = ubound(fxy,3)
  if(Symmetry == OCTANT) then 
     if(X(1,1,1).lt.0.D0) then
        imin = imin + 1
     end if
     if(Y(1,1,1).lt.0.D0) then
        jmin = jmin + 1
     end if
     if(Z(1,1,1).lt.0.D0) then
        kmin = kmin + 1
     end if
  else if(Symmetry == AXISYM) then
     if(X(1,1,1).lt.0.D0) then
        imin = imin + 1
     end if
     if(Z(1,1,1).lt.0.D0) then
        kmin = kmin + 1
     end if
  else if(Symmetry == EQUATORIAL) then
     if(Z(1,1,1).lt.0.D0) then
        kmin = kmin + 1
     end if
  end if
  !
  ddx=1.d0/(X(imin+1,1,1)-X(imin,1,1))
  ddy=1.d0/(Y(1,jmin+1,1)-Y(1,jmin,1))
  ddz=1.d0/(Z(1,1,kmin+1)-Z(1,1,kmin))
  !
  do k=kmin,kmax-1
     do j=jmin,jmax-1
	do i=imin,imax-1
	   if (Symmetry==AXISYM) then
	      Bx_rhs(i,j,k)=ddz*(fzx(i,j,k)-fzx(i,j,k+1))/X(i,1,1)
	      By_rhs(i,j,k)=ddx*(fxy(i,j,k)-fxy(i+1,j,k)) + &
                   ddz*(fzy(i,j,k)-fzy(i,j,k+1))
	      Bz_rhs(i,j,k)=ddx*(fxz(i,j,k)-fxz(i+1,j,k))/X(i,1,1)
	   else
	      Bx_rhs(i,j,k)=ddy*(fyx(i,j,k)-fyx(i,j+1,k)) + &
                   ddz*(fzx(i,j,k)-fzx(i,j,k+1))
	      By_rhs(i,j,k)=ddx*(fxy(i,j,k)-fxy(i+1,j,k)) + &
                   ddz*(fzy(i,j,k)-fzy(i,j,k+1))
	      Bz_rhs(i,j,k)=ddx*(fxz(i,j,k)-fxz(i+1,j,k)) + &
                   ddy*(fyz(i,j,k)-fyz(i,j+1,k))

	   end if
	end do
     end do
  end do
  !write(*,*) "advect_b_ct:",Bx_rhs(2,2,2),fyx(2,2,2),fyx(2,2+1,2),fzx(2,2,2),fzx(2,2,2+1)

end subroutine advect_b_ct
