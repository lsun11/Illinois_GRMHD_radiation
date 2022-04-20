!-----------------------------------------------------------------------------
!
! $Id: laguna_initial_metric.f90,v 1.1.1.1 2006/02/23 17:48:41 zetienne Exp $
!
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!
! Find Laguna PRD 69, 104020 (2004) metric initial data for boosted BH
!
!-----------------------------------------------------------------------------
  subroutine laguna_initial_metric(ex, x, y, z,RP,            &
                    Axx, Axy, Axz, Ayy, Ayz, Azz, alpha, Bx,By,Bz, phi, &
                    bigP,bhpos,xbh,zbh,mbh,Symmetry,puncture_enable)
  implicit none

!
! Input parameters:
!
  integer, dimension(3),     intent(in) :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in) :: x,y,z,RP
  integer, intent(in)                   :: Symmetry,puncture_enable
  real*8, intent(in)                    :: bigP,bhpos,mbh
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: alpha,Bx,By,Bz,phi

!
! Other variables:
!
  integer :: i,j,k,m,n
  integer :: imin,jmin,kmin,imax,jmax,kmax
  real*8  :: el(0:3),boo(0:3),phh(0:3),eta(0:3,0:3),pel(0:3,0:3)
  real*8  :: r,phii,bigR,bigS,rp1,st,ct,sp,cp,P0,P2,rp0
  real*8  :: u0,u2,v0,v2,nx,ny,nz,pknk,xi2,zk2,rfact,xp,yp,zp,rpp
  real*8  :: betar,betat,betax,betay,betaz,u,v,lapsei,xbh,zbh
!  real*8,  dimension(ex(1),ex(2),ex(3)) :: qxx,qxy,qxz,qyy,qyz,qzz,q
  real*8, parameter :: zeo=0.d0,one=1.d0,two=2.d0,thr=3.d0
  real*8, parameter :: f1o3=1.d0/3.d0,inbd=1.d0,hlf=0.5d0,fou=4.d0

!~~~~~> Interface

  interface
    real*8 function radius(x,y,z,Symmetry)
      implicit none
      real*8,  intent(in) :: x,y,z
      integer, intent(in) :: Symmetry
    end function radius
  end interface

!
! Input translation
!
  imin = lbound(phi,1)
  jmin = lbound(phi,2)
  kmin = lbound(phi,3)
  imax = ubound(phi,1)
  jmax = ubound(phi,2)
  kmax = ubound(phi,3)


  write(*,*)'parameters: Sym',Symmetry,' bhpos:',bhpos,' bigP:',bigP 
!
! go to each gridpoint...
!
  do k = kmin, kmax
   do j = jmin, jmax
    do i = imin, imax
!
       r=sqrt(x(i,j,k)**2+y(i,j,k)**2+z(i,j,k)**2)
       rp0=RP(i,j,k)
       if(r.gt.0)then
          rfact=rp0/r
       else
          rfact=1.0
       endif
       xp=x(i,j,k)*rfact
       yp=y(i,j,k)*rfact
       zp=z(i,j,k)*rfact

       if(Symmetry.eq.4) then
          zbh=bhpos
          xbh=0.
          zk2=zp-bhpos
          xi2=xp
       else if (Symmetry.eq.1) then
          xbh=bhpos
          zbh=0.
          xi2=xp-bhpos
          zk2=zp
       else
          write(6,*)'We have not programmed that symmetry yet'
          write(6,*)'for the laguna initial data'
          stop
       endif
       rpp=sqrt(xi2**2+yp**2+zk2**2)
       rp1=rpp+1.0
       
       if(Symmetry.eq.4) then
          ! Axisymmetry: momentum in z direction          
          st = sqrt(1.0-zk2**2/rpp**2)
          ct = zk2/rpp
          cp = xp/sqrt(xp**2+yp**2)
          sp = yp/sqrt(xp**2+yp**2)
          P0=1.0
          P2=1.5*zk2**2/rpp**2-0.5d0
          
       else if (Symmetry.eq.1) then
          ! Equatorial: momentum in x-direction
          st = sqrt(1.0-xi2**2/rpp**2)
          ct = xi2/rpp
          cp = zp/sqrt(zp**2+yp**2)
          sp = yp/sqrt(zp**2+yp**2)
          P0=1.0
          P2=1.5*xi2**2/rpp**2-0.5d0
          
       endif
          
          ! DBP eq. 20
       u0=mbh**4+10.0d0*mbh**3*rpp+40.0d0*mbh**2*rpp**2+80.0d0*mbh*rpp**3+80.0d0*rpp**4

	! DBP eq. 21
	u2=mbh/5.0d0/rpp**3*(42.0*mbh**5*rpp+378.0*mbh**4*rpp**2+1316.0*mbh**3*rpp**3+2156.0*mbh**2*rpp**4 &
		+1536.0*mbh*rpp**5+240.0*rpp**6+21.0d0*mbh*(mbh+2.0*rpp)**5*log(mbh/(mbh+2.0*rpp)))

       u=bigP**2*mbh/8.0d0/(mbh+2.0*rpp)**5*(P0*u0+P2*u2)

       if(puncture_enable == 1) then
          phii=1.0d0+u+0.5d0*mbh/rpp
       else if (puncture_enable == 2 .or. puncture_enable == 3) then
          if(rpp.gt.0.5d0*mbh) then
             phii=1.0d0+u+0.5d0*mbh/rpp
          else
             phii=2.875d0-5.0d0*(rpp/mbh)**2+6.0d0*(rpp/mbh)**4+u
          endif
       else if (puncture_enable == 4 .or. puncture_enable == 5) then
          if(rpp.gt.0.5d0*mbh) then
             phii=1.0d0+u+0.5d0*mbh/rpp
          else
             phii=2.0d0+u
          endif
       else if (puncture_enable == 6 .or. puncture_enable == 7) then
          if(rpp.gt.0.25d0*mbh) then
             phii=1.0d0+u+0.5d0*mbh/rpp
          else
             phii=3.0d0+u
          endif
       endif
       phi(i,j,k)=log(phii)
       if(i.eq.2.and.j.eq.2.and.k.eq.2)write(6,*)'phi0',phi(i,j,k),phii,u,rpp,x(2,2,2),z(2,2,2)
!       alpha(i,j,k)=(1.0+v-cl/r)/phii
       alpha(i,j,k)=phii**(-2)
!
!-- the 3-metric
!
     nx=xi2/rpp
     ny=yp/rpp
     nz=zk2/rpp

     if(puncture_enable == 3 .or. puncture_enable == 5 .or. puncture_enable == 7)then
	if(rpp.lt.mbh/2.0)then
	    rpp=0.25d0+(rpp/mbh)**2
	endif
     endif			

     if(Symmetry.eq.4)then
        ! Axisymmetry: momentum in z-direction
        pknk=bigP*nz
        Azz(i,j,k) = 1.5*phii**(-6)*(2.0*bigP*nz - (1.0 - nz**2) *pknk)/rpp**2
        Ayy(i,j,k) = 1.5*phii**(-6)*(-1.0)*(1.0-ny**2)*pknk/rpp**2
        Axx(i,j,k) = 1.5*phii**(-6)*(-1.0)*(1.0-nx**2)*pknk/rpp**2
        Ayz(i,j,k) = 1.5*phii**(-6)*(    bigP*ny +        nz*ny  *pknk)/rpp**2
        Axz(i,j,k) = 1.5*phii**(-6)*(    bigP*nx +        nz*nx  *pknk)/rpp**2
        Axy(i,j,k) = 1.5*phii**(-6)*nx*ny*pknk/rpp**2

     else if (Symmetry.eq.1)then
        ! Equatorial: momentum in x-direction
        pknk=bigP*nx
        Axx(i,j,k) = 1.5*phii**(-6)*(2.0*bigP*nx - (1.0 - nx**2) *pknk)/rpp**2
        Ayy(i,j,k) = 1.5*phii**(-6)*(-1.0)*(1.0-ny**2)*pknk/rpp**2
        Azz(i,j,k) = 1.5*phii**(-6)*(-1.0)*(1.0-nz**2)*pknk/rpp**2
        Axy(i,j,k) = 1.5*phii**(-6)*(    bigP*ny +        nx*ny  *pknk)/rpp**2
        Axz(i,j,k) = 1.5*phii**(-6)*(    bigP*nz +        nx*nz  *pknk)/rpp**2
        Ayz(i,j,k) = 1.5*phii**(-6)*ny*nz*pknk/rpp**2
     endif

    end do
   end do
  end do

  alpha = alpha - one

 Bx=0
 By=0
 Bz=0

! write(*,*)'Initial data-phi:',phi(1,1,2),phi(10,10,11),' Axx:',Axx(1,1,2),Axx(10,10,11),x(10,10,11),y(10,10,11),z(10,10,11),' lap:',alpha(1,1,2),alpha(10,10,11),x(10,10,11),y(10,10,11),z(10,10,11)

 return
 end subroutine laguna_initial_metric
