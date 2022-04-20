!------------------------------------------------------------
! 
! Given a "physical" radius R, this subroutine calculates 
! the fisheye radius rbar. 
!
!------------------------------------------------------------
!
subroutine R_to_rfish(R, rbar, rbar_max, aFE, kappaFE, r0FE, sFE, n)
 implicit none
 integer :: n,i,j,k
 integer, parameter :: nd = 20
 real*8 :: a(0:nd),kappa(nd),r0(nd),s(nd)
 real*8 :: aFE(n+1),kappaFE(n+1),r0FE(n+1),sFE(n+1)
 real*8 :: R, drbar, rbar_max, rbar1, rbar2, rbart
 real*8, dimension(1) :: rbar
 integer, parameter :: ntab = 20
 real*8  :: rtab(0:ntab)
 real*8  :: func_FE,zbrent
 external func_FE
 real*8, parameter :: tol = 1.d-14
 logical :: exit_do
!
 if (R==0.d0) then
    rbar = 0.d0
    return
 end if

 if (n .gt. nd) then
    write(*,*) 'n > nd in subroutine R_to_rfish'
    stop
 end if
 a(0) = aFE(1)
 do i=1,n
    a(i) = aFE(i+1)
    kappa(i) = kappaFE(i+1)
    r0(i) = r0FE(i+1)
    s(i) = sFE(i+1)
 end do

 drbar = rbar_max/dble(ntab)
 do k=0,ntab
    rbart = dble(k)*drbar
    call fisheye_transform(rtab(k),rbart,a,kappa,r0,s,n,nd)
 end do

 rbar1 = -1.d0
 j=0
 exit_do = .FALSE.
 do
    if (R .ge. rtab(j) .and. R .lt. rtab(j+1)) then
       exit_do = .TRUE.
       rbar1 = dble(j)*drbar
       rbar2 = dble(j+1)*drbar
    end if
    if (j==ntab-1 .or. exit_do) exit
    j = j + 1
 end do

 if (rbar1 .lt. 0.d0) then
    write(*,*) 'Cannot bracket the root in subroutine R_to_rfish!'
    write(*,*) 'R = ',R
    write(*,*) 'Allowed min and max of R: ',rtab(0),rtab(ntab)
    stop
 end if
  
 rbar = zbrent(func_FE,rbar1,rbar2,tol,R,a,kappa,s,r0,n,nd)

end subroutine R_to_rfish

subroutine Set_xpc_fish(ext,xpc1_fish_gf,xpc2_fish_gf,xpc1_fish,xpc2_fish)
 implicit none
 integer, dimension(3) :: ext
 real*8, dimension(ext(1),ext(2),ext(3)) :: xpc1_fish_gf,xpc2_fish_gf
 real*8 :: xpc1_fish,xpc2_fish
!
 xpc1_fish_gf = xpc1_fish
 xpc2_fish_gf = xpc2_fish
end subroutine Set_xpc_fish

! Fisheye transformation
!
   subroutine fisheye_transform(r,rbar,a,kappa,r0,s,n,nd)
   implicit none
   integer :: n,j,nd
   real*8 :: a(0:nd),kappa(nd),r0(nd),s(nd)
   real*8 :: r,rbar
!
   r = a(n)*rbar
   do j=1,n
      if (rbar-r0(j) .gt. 0.d0) then
         r = r + kappa(j)*(2.d0*r0(j)/s(j) + log( &
                (1.d0+exp(-2.d0/s(j)*(rbar+r0(j))) )/ &
                (1.d0+exp(-2.d0/s(j)*(rbar-r0(j))) ) ) )
      else
         r = r + kappa(j)*(2.d0*rbar/s(j) + log( &
                (1.d0+exp(-2.d0/s(j)*(rbar+r0(j))) )/ &
                (1.d0+exp(-2.d0/s(j)*(r0(j)-rbar)) ) ) )
      end if
   end do
   end subroutine fisheye_transform

! Define the function for the root finder
!
   real*8 function func_FE(rbar,Rp,a,kappa,s,r0,n,nd)
   implicit none
   real*8 :: rbar,r,Rp
   integer :: n,nd
   real*8 :: a(0:nd),kappa(nd),r0(nd),s(nd)
!
   call fisheye_transform(r,rbar,a,kappa,r0,s,n,nd)
   func_FE = (r-Rp)/(Rp + 1.d-306)
   end

! *************************************************************
!
! Numerical root-finding subroutine
!
! *************************************************************
      FUNCTION zbrent(func,x1,x2,tol,Rp,a_fish,kappa,s_fish,r0,n,nd)
      implicit none
      INTEGER :: ITMAX
      REAL*8  :: zbrent,tol,x1,x2,func,EPS
      EXTERNAL func
      PARAMETER (ITMAX=100,EPS=1.d-15)
      INTEGER :: iter
      REAL*8 :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      integer :: n,nd
      real*8 :: Rp
      real*8 :: a_fish(0:nd),kappa(nd),r0(nd),s_fish(nd)
!
      a=x1
      b=x2
      fa=func(a,Rp,a_fish,kappa,s_fish,r0,n,nd)
      fb=func(b,Rp,a_fish,kappa,s_fish,r0,n,nd)
      if((fa.gt.0.d0.and.fb.gt.0.d0).or.(fa.lt.0.d0.and.fb.lt.0.d0)) then
        write(*,*) 'root must be bracketed for zbrent'
        stop
      end if
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0.d0.and.fc.gt.0.d0).or.(fb.lt.0.d0.and.fc.lt.0.d0))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.d0*EPS*abs(b)+0.5d0*tol
        xm=0.5d0*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.d0)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.d0*xm*s
            q=1.d0-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.d0*xm*q*(q-r)-(b-a)*(r-1.d0))
            q=(q-1.d0)*(r-1.d0)*(s-1.d0)
          endif
          if(p.gt.0.d0) q=-q
          p=abs(p)
          if(2.d0*p .lt. min(3.d0*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b,Rp,a_fish,kappa,s_fish,r0,n,nd)
11    continue
      write(*,*) 'zbrent exceeding maximum iterations'
      zbrent=b
      return
      END
