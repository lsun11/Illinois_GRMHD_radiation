!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     calculates the special relativistic shock tube problem          c
!     (see Thompson 1986)                                             c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine shocktests_rieman_analytic(t,Rmin,Rmax,gamma,rho_l,rho_r,nsg,R_out,rho_out,P_out,v_out)
  implicit none
  real*8 :: gamma,rho_l,rho_r

  real*8,dimension(nsg) :: R_out,rho_out,P_out,v_out

  double precision pr,pl,nr,nl,gam,cl,cr,zbrent,Rstep
  double precision pm,nml,nmr,vm,vt,vs,help,t,func,um1,um2
  double precision um,yl,wm,f,c,y,vss,esmax,nsmax,esmin
  character*1 chr
  logical succes

  integer ngrid,i,nmax,sm_device,pflag,nsg
  parameter(ngrid=300)
  real U(ngrid),R(ngrid),n(ngrid),p(ngrid),e(ngrid)
  real Rs(ngrid),ps(ngrid),ns(ngrid),Us(ngrid),es(ngrid)

  common/root/pr,gam,cl,cr,pl,nl,nr,help
  external func
  double precision Rmin, Rmax, v

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !     get data                                                        c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  nr = rho_r
  nl = rho_l
  gam = gamma

  help=gam/(gam-1.)

  pl = nl**gam
  pr = nr**gam

  write(*,500) t,Rmin,Rmax,gam
  write(*,600) nl,pl,nr,pr
500 format('# t : ',e12.4,'  Rmin : ',e12.4,'  Rmax : ',e12.4,'  Gamma : ',e12.4)
600 format('# nl : ',e12.4,' pl : ',e12.4,' nr : ',e12.4,' pr : ',e12.4)

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !     get parameters                                                  c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  cl=sqrt(gam*pl/(nl+gam/(gam-1.)*pl))
  cr=sqrt(gam*pr/(nr+gam/(gam-1.)*pr))

  write(*,*) 'cl =',cl,'  cr =',cr

  um1=1.d-5
  um2=1.d-3

  call zbrac(func,um1,um2,succes)
  write(*,*) "HIosdifasdopfi",um1,um2,succes
  um = zbrent(func,um1,um2,1.d-10)
  wm=sqrt(1.+um**2)
  vm=um/wm
  yl=sqrt(help*pl/nl)
  pm=pl*(((f(yl))**2-(f(um))**sqrt(gam-1.))/(2.*f(yl)*(f(um))**(0.5*sqrt(gam-1.))*yl))**(2.*help)
  vs=(pm-pr+help*pm*um**2)/((nr+help*wm*pm)*um)
  nmr=nr*(pm-pr+help*pm*um**2)/(wm*(pm-pr)-nr*um**2)
  nml=nl*(pm/pl)**(1./gam)
  vt=vm-c(nml,pm)/(wm**2*(1.-c(nml,pm)*vm))

  write(*,*) 'um =',um,'  vm =',vm
  write(*,*) 'vs =',vs,'  vt =',vt

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !     assign data to gridpoints                                       c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Rs(1) = Rmin
  ns(1) = nl
  ps(1) = pl
  Us(1) = 0.

  Rs(2) = -cl*t
  ns(2) = nl
  ps(2) = pl
  Us(2) = 0.

  Rs(nsg) = Rmax
  ns(nsg) = nr
  ps(nsg) = pr
  Us(nsg) = 0.

  Rs(nsg-1) = vs*t
  ns(nsg-1) = nr
  ps(nsg-1) = pr
  Us(nsg-1) = 0.

  Rs(nsg-2) = vs*t
  ns(nsg-2) = nmr
  ps(nsg-2) = pm
  Us(nsg-2) = um

  Rs(nsg-3) = vm*t
  ns(nsg-3) = nmr
  ps(nsg-3) = pm
  Us(nsg-3) = um

  Rs(nsg-4) = vm*t
  ns(nsg-4) = nml
  ps(nsg-4) = pm
  Us(nsg-4) = um

  Rs(nsg-5) = vt*t
  ns(nsg-5) = nml
  ps(nsg-5) = pm
  Us(nsg-5) = um

  do i=3,nsg-6
     Us(i) = real(i-2)/real(nsg-7)*um
     y=((f(yl))**2-(f(dble(Us(i))))**sqrt(gam-1.))/(2.*f(yl)*(f(dble(Us(i))))**(0.5*sqrt(gam-1.)))
     ns(i) = nl*(y/yl)**(2./(gam-1))
     ps(i) = pl*(ns(i)/nl)**gam 
     vss = Us(i)/sqrt(1. + Us(i)**2)
     Rs(i) = t*(vss - c(dble(ns(i)),dble(ps(i)))/((1.+Us(i)**2)*(1.- c(dble(ns(i)),dble(ps(i)))*vss)))
!     Rs(i) = t*(vss - sqrt(gam*ps(i)/(ns(i)+help*ps(i)))/((1.+Us(i)**2)*(1.- sqrt(gam*ps(i)/(ns(i)+help*ps(i)))*vss)))
!     write(*,*) "Rs = ",i,Rs(i),dble(ns(i)),dble(ps(i)),sqrt(gam*ps(i)/(ns(i)+help*ps(i)))
  enddo

  esmax = 0.
  esmin = 1.d10
  nsmax = 0.
  do i=1,nsg
     es(i) = ps(i)/(ns(i)*(gam-1.))
     if (es(i).gt.esmax) esmax = es(i)
     if (es(i).lt.esmin) esmin = es(i)
     if (ns(i).gt.nsmax) nsmax = ns(i)
  enddo

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !     plot                                                            c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  do i = 1, nsg
     R_out(i) = Rs(i)
     if(i>1) then
        if(Rs(i)==Rs(i-1)) R_out(i) = R_out(i-1) + 1e-15
     end if
     rho_out(i) = ns(i)
     v_out(i) = Us(i)/sqrt(1.D0 + Us(i)*Us(i))
     P_out(i) = ps(i)
  enddo
end subroutine shocktests_rieman_analytic


double precision function func(um)
  implicit none
  double precision pr,gam,cl,cr,pl,nl,um,pm,vs,nmr,yl,f,help,nr
  double precision wm
  common/root/pr,gam,cl,cr,pl,nl,nr,help

  wm=sqrt(1.+um**2)
  yl=sqrt(help*pl/nl)
  pm=pl*(((f(yl))**2-(f(um))**sqrt(gam-1.))/(2.*f(yl)*(f(um))**(0.5*sqrt(gam-1.))*yl))**(2.*help)
  vs=(pm-pr+help*pm*um**2)/((nr+help*wm*pm)*um)
  nmr=nr*(pm-pr+help*pm*um**2)/(wm*(pm-pr)-nr*um**2)
  func=vs*(wm**2*(nmr+help*pm)-pm-nr-pr/(gam-1.)) - (nmr+help*pm)*wm*um
  return
end function func

double precision function f(x)
  implicit none
  double precision x

  f=x+sqrt(1.+x**2)

  return
end function f

double precision function c(n,p)
  implicit none
  double precision pr,gam,cl,cr,pl,nl,nr,help,n,p
  common/root/pr,gam,cl,cr,pl,nl,nr,help

  c=sqrt(gam*p/(n+help*p))

  return
end function c
