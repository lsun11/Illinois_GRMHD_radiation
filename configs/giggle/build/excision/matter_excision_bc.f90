function fr2(x,y,z)
  implicit none
  real*8 :: x,y,z
  real*8 :: fr2
  fr2 = x**2 + y**2 + z**2
end function fr2
function fr1(x,y,z)
  implicit none
  real*8 :: x,y,z,fr1
  fr1 = sqrt(x**2 + y**2 + z**2)
end function fr1
function frm1(x,y,z)
  implicit none
  real*8 :: x,y,z,frm1
  frm1 = 1.d0/sqrt(x**2 + y**2 + z**2)
end function frm1
!-------------------------------------------------------------------------+
!
! Use Alcubierres simple excision method for fluid
!
! excision region : sphere
!
!-------------------------------------------------------------------------+
subroutine hydro_ezbc(ex,x,y,z,rhos,tau,mhdstx,mhdsty,mhdstz,stx,sty,stz, &
     rhob,P,enth,vx,vy,vz,w,sbt,sbx,sby,sbz,Bx,By,Bz,alp,betx,bety,betz,phi, &
     gxx,gxy,gxz,gyy,gyz,gzz,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
     symmetry,n,dm,kpoly_init)
  implicit none
!~~~~~~> Input arguments
  integer, intent(in) :: ex(1:3),symmetry
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   ::x,y,z
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::rhos,tau,mhdstx,mhdsty,mhdstz
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::stx,sty,stz
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::rhob,P,enth,vx,vy,vz,w
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::sbt,sbx,sby,sbz
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::Bx,By,Bz
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   ::alp,betx,bety,betz,phi
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   ::gxx,gxy,gxz,gyy,gyz,gzz
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   ::gupxx,gupxy,gupxz
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   ::gupyy,gupyz,gupzz
  integer :: dm(1:ex(1),1:ex(2),1:ex(3))
  real*8    :: n,kpoly_init
!~~~~~~> local variable
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k,flax,flay,flaz,ni
  integer :: ksym,jsym,isym,diax,diay,diaz
  real*8  :: sf(1:3),rr,xs,zs,rxs,rzs,rxzs,c1,c2
  real*8  :: al,ephi,u0,au0,u_x,u_y,u_z,sb_x,sb_y,sb_z,sb2,v2,gam,f1os4pi
  real*8  :: r1,r2,fac1,fac2,fac3,r,sinthet,costhet,sinphi,cosphi,Pmin
  integer :: extrap_method
  integer, parameter :: inEZ=0, onbd=1, outEZ=2
  real*8,  parameter :: zeo=0.d0, one=1.d0, two=2.d0
  integer, parameter :: nosymmetry=0, equatorial=1, octant=2
  integer, parameter :: pi_symm=3, axisym=4
  integer, parameter :: FLAT=0, LINEAR=1
  real*8,dimension(ex(1),ex(2),ex(3))   :: vr,vt,vp,Br,Bt,Bp
!~~~~~~> interface
  interface
!     subroutine find_excision_zone(ex,x,y,z,dm,innerbd,Symmetry)
!       implicit none
!       integer, intent(in)  :: ex(1:3)
!       real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   ::x,y,z
!       integer, intent(out) :: dm(1:ex(1),1:ex(2),1:ex(3))
!       real*8               :: innerbd
!       integer, intent(in)  :: Symmetry
!     end subroutine find_excision_zone
  end interface
  real*8  :: fr2,fr1,frm1
  EXTERNAL fr2,fr1,frm1
  extrap_method = FLAT
!  extrap_method = LINEAR
  imin = lbound(rhos,1)
  jmin = lbound(rhos,2)
  kmin = lbound(rhos,3)
  imax = ubound(rhos,1)
  jmax = ubound(rhos,2)
  kmax = ubound(rhos,3)
 f1os4pi = 1.d0/sqrt(4.d0*acos(-1.d0))
!!$  c1 = 2.d0
  c1 = 1.d0
  c2 = c1 - 1.d0
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
!~~~~~~> Find the grid points on the boundary
!  call find_excision_zone(ex,x,y,z,dm,innerbd,Symmetry)
!~~~~~~> inner boundary condition
  if( symmetry == nosymmetry ) then
     isym = imin + 1
     jsym = jmin + 1
     ksym = kmin + 1
  else if( symmetry == equatorial ) then
     isym = imin + 1
     jsym = jmin + 1
     ksym = kmin
  else if( symmetry == octant ) then
     isym = imin
     jsym = jmin
     ksym = kmin
  else if( symmetry == pi_symm ) then
     isym = imin + 1
     jsym = jmin
     ksym = kmin
  endif
  do k = kmin, kmax, 1
     do j = jmin, jmax, 1
        do i = imin, imax, 1
           r = sqrt(x(i,1,1)**2 + y(1,j,1)**2 + z(1,1,k)**2)
           costhet = z(1,1,k)/r
           sinthet = sqrt(1.d0 - costhet**2)
           cosphi  = x(i,1,1)/(r*sinthet)
           sinphi  = y(1,j,1)/(r*sinthet)
           vr(i,j,k) = sinthet*cosphi*vx(i,j,k) + &
                sinthet*sinphi*vy(i,j,k) + costhet*vz(i,j,k)
           vt(i,j,k) = costhet*cosphi/r*vx(i,j,k) + &
                costhet*sinphi/r*vy(i,j,k) - sinthet/r*vz(i,j,k)
           vp(i,j,k) = -sinphi/sinthet/r*vx(i,j,k) + &
                cosphi/sinthet/r*vy(i,j,k)
           Br(i,j,k) = sinthet*cosphi*Bx(i,j,k) + &
                sinthet*sinphi*By(i,j,k) + costhet*Bz(i,j,k)
           Bt(i,j,k) = costhet*cosphi/r*Bx(i,j,k) + &
                costhet*sinphi/r*By(i,j,k) - sinthet/r*Bz(i,j,k)
           Bp(i,j,k) = -sinphi/sinthet/r*Bx(i,j,k) + &
                cosphi/sinthet/r*By(i,j,k)
        end do
     end do
  end do
  if(symmetry .ne. axisym) then
     do k = ksym, kmax-1, 1
        do j = jsym, jmax-1, 1
           do i = isym, imax-1, 1
              if( dm(i,j,k) == onbd ) then
                 ni = 2
                 flax = 1
                 flay = 1
                 flaz = 1
                 if(x(i,1,1) < zeo) flax = -1
                 if(y(1,j,1) < zeo) flay = -1
                 if(z(1,1,k) < zeo) flaz = -1
                 sf(1) = x(i,1,1) * flax
                 sf(2) = y(1,j,1) * flay
                 sf(3) = z(1,1,k) * flaz
                 if( dm(i+flax,j     ,k     ) == onbd ) sf(1) = zeo
                 if( dm(i     ,j+flay,k     ) == onbd ) sf(2) = zeo
                 if( dm(i     ,j     ,k+flaz) == onbd ) sf(3) = zeo
                 if((dm(i+flax,j,k) == outEZ) .or. (dm(i,j+flay,k) == outEZ)) &
                      sf(3) = zeo
                 if(( dm(i+flax,j     ,k     ) == onbd ) .and. &
                      ( dm(i     ,j+flay,k     ) == onbd ) .and. &
                      ( dm(i     ,j     ,k+flaz) == onbd )) then
                    ni = 1
                    if( dm(i     ,j+flay,k+flaz) == outEZ ) sf(2) = one
                    if( dm(i+flax,j     ,k+flaz) == outEZ ) sf(3) = one
                    if( dm(i+flax,j+flay,k     ) == outEZ ) then
                       sf(1) = one
                       sf(2) = zeo
                       sf(3) = zeo
                    endif
                    if(( dm(i+flax,j+flay,k) /= outEZ ) .and. &
                         ( dm(i,j+flay,k+flaz) /= outEZ ) .and. &
                         ( dm(i+flax,j,k+flaz) /= outEZ )) then
                         write(*,*) 'something wrong in find_excision_zone'
                    end if
                 endif
                 if( ni == 2 ) then
                    rr = sf(1) + sf(2) + sf(3)
                    sf(1) = sf(1) / rr
                    sf(2) = sf(2) / rr
                    sf(3) = sf(3) / rr
                    rhob(i,j,k) = sf(1) * rhob(i+flax,j,k) + &
                         sf(2) * rhob(i,j+flay,k) + &
                         sf(3) * rhob(i,j,k+flaz)
                    P(i,j,k) = sf(1) * P(i+flax,j,k) + &
                         sf(2) * P(i,j+flay,k) + &
                         sf(3) * P(i,j,k+flaz)
                    vx(i,j,k) = sf(1) * vx(i+flax,j,k) + &
                         sf(2) * vx(i,j+flay,k) + &
                         sf(3) * vx(i,j,k+flaz)
                    vy(i,j,k) = sf(1) * vy(i+flax,j,k) + &
                         sf(2) * vy(i,j+flay,k) + &
                         sf(3) * vy(i,j,k+flaz)
                    vz(i,j,k) = sf(1) * vz(i+flax,j,k) + &
                         sf(2) * vz(i,j+flay,k) + &
                         sf(3) * vz(i,j,k+flaz)
                    Bx(i,j,k) = sf(1) * Bx(i+flax,j,k) + &
                         sf(2) * Bx(i,j+flay,k) + &
                         sf(3) * Bx(i,j,k+flaz)
                    By(i,j,k) = sf(1) * By(i+flax,j,k) + &
                         sf(2) * By(i,j+flay,k) + &
                         sf(3) * By(i,j,k+flaz)
                    Bz(i,j,k) = sf(1) * Bz(i+flax,j,k) + &
                         sf(2) * Bz(i,j+flay,k) + &
                         sf(3) * Bz(i,j,k+flaz)
                 else if( ni == 1 ) then
                    rr = sf(1) + sf(2) +sf(3)
                    sf(1) = sf(1) / rr
                    sf(2) = sf(2) / rr
                    sf(3) = sf(3) / rr
                    rhob(i,j,k) = sf(1) * rhob(i+flax,j+flay,k) + &
                         sf(2) * rhob(i,j+flay,k+flaz) + &
                         sf(3) * rhob(i+flax,j,k+flaz)
                    P(i,j,k) = sf(1) * P(i+flax,j+flay,k) + &
                         sf(2) * P(i,j+flay,k+flaz) + &
                         sf(3) * P(i+flax,j,k+flaz)
                    vx(i,j,k) = sf(1) * vx(i+flax,j+flay,k) + &
                         sf(2) * vx(i,j+flay,k+flaz) + &
                         sf(3) * vx(i+flax,j,k+flaz)
                    vy(i,j,k) = sf(1) * vy(i+flax,j+flay,k) + &
                         sf(2) * vy(i,j+flay,k+flaz) + &
                         sf(3) * vy(i+flax,j,k+flaz)
                    vz(i,j,k) = sf(1) * vz(i+flax,j+flay,k) + &
                         sf(2) * vz(i,j+flay,k+flaz) + &
                         sf(3) * vz(i+flax,j,k+flaz)
                    Bx(i,j,k) = sf(1) * Bx(i+flax,j+flay,k) + &
                         sf(2) * Bx(i,j+flay,k+flaz) + &
                         sf(3) * Bx(i+flax,j,k+flaz)
                    By(i,j,k) = sf(1) * By(i+flax,j+flay,k) + &
                         sf(2) * By(i,j+flay,k+flaz) + &
                         sf(3) * By(i+flax,j,k+flaz)
                    Bz(i,j,k) = sf(1) * Bz(i+flax,j+flay,k) + &
                         sf(2) * Bz(i,j+flay,k+flaz) + &
                         sf(3) * Bz(i+flax,j,k+flaz)
                 endif
                    al = 1.d0 + alp(i,j,k)
                    ephi = exp(4.d0 * phi(i,j,k))
                    ! Compute al*u0-1
                    v2 = ephi*(gxx(i,j,k)*(vx(i,j,k) + betx(i,j,k))**2 + &
                         2.d0*gxy(i,j,k)*(vx(i,j,k) + betx(i,j,k)) &
                         *(vy(i,j,k) + bety(i,j,k)) +  &
                         2.d0*gxz(i,j,k)*(vx(i,j,k) + betx(i,j,k)) &
                         *(vz(i,j,k) + betz(i,j,k)) +  &
                         gyy(i,j,k)*(vy(i,j,k) + bety(i,j,k))**2 + &
                         2.d0*gyz(i,j,k)*(vy(i,j,k) + bety(i,j,k)) &
                         *(vz(i,j,k) + betz(i,j,k)) +  &
                         gzz(i,j,k)*(vz(i,j,k) + betz(i,j,k))**2 )/al**2
                    ! *** Check for superluminal velocity ***
                    if (v2 .gt. 1.d0) then
                       vx(i,j,k)=(vx(i,j,k)+betx(i,j,k))*sqrt(0.99d0/v2)-betx(i,j,k)
                       vy(i,j,k)=(vy(i,j,k)+bety(i,j,k))*sqrt(0.99d0/v2)-bety(i,j,k)
                       vz(i,j,k)=(vz(i,j,k)+betz(i,j,k))*sqrt(0.99d0/v2)-betz(i,j,k)
                       v2 = 0.99d0
                    end if
                    ! ***************************************
                    gam = sqrt(1.d0-v2)
                    au0 = 1.d0/gam
                    u0  = au0/al
                    sbt(i,j,k) = f1os4pi*ephi*u0/al*( &
                         gxx(i,j,k)*Bx(i,j,k)*(vx(i,j,k)+betx(i,j,k)) + &
                         gxy(i,j,k)*( Bx(i,j,k)*(vy(i,j,k)+bety(i,j,k)) &
                         + By(i,j,k)*(vx(i,j,k)+betx(i,j,k)) ) + &
                         gxz(i,j,k)*( Bx(i,j,k)*(vz(i,j,k)+betz(i,j,k)) &
                         + Bz(i,j,k)*(vx(i,j,k)+betx(i,j,k)) ) + &
                         gyy(i,j,k)*By(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gyz(i,j,k)*( By(i,j,k)*(vz(i,j,k)+betz(i,j,k)) &
                         + Bz(i,j,k)*(vy(i,j,k)+bety(i,j,k)) ) + &
                         gzz(i,j,k)*Bz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    sbx(i,j,k) = f1os4pi*Bx(i,j,k)/(al*u0) &
                         + sbt(i,j,k)*vx(i,j,k)
                    sby(i,j,k) = f1os4pi*By(i,j,k)/(al*u0) &
                         + sbt(i,j,k)*vy(i,j,k)
                    sbz(i,j,k) = f1os4pi*Bz(i,j,k)/(al*u0) &
                         + sbt(i,j,k)*vz(i,j,k)
                    sb2 = -(al*sbt(i,j,k))**2 &
                         + ephi*( gxx(i,j,k)*(&
                         sbx(i,j,k)+betx(i,j,k)*sbt(i,j,k))**2 + &
                         2.d0*gxy(i,j,k)*(sbx(i,j,k)+betx(i,j,k)*sbt(i,j,k)) &
                         *(sby(i,j,k)+bety(i,j,k)*sbt(i,j,k)) + &
                         2.d0*gxz(i,j,k)*(sbx(i,j,k)+betx(i,j,k)*sbt(i,j,k)) &
                         *(sbz(i,j,k)+betz(i,j,k)*sbt(i,j,k)) + &
                         gyy(i,j,k)*(sby(i,j,k)+bety(i,j,k)*sbt(i,j,k))**2 + &
                         2.d0*gyz(i,j,k)*(sby(i,j,k)+bety(i,j,k)*sbt(i,j,k)) &
                         *(sbz(i,j,k)+betz(i,j,k)*sbt(i,j,k)) + &
                         gzz(i,j,k)*(sbz(i,j,k)+betz(i,j,k)*sbt(i,j,k))**2 )
                    u_x = u0*ephi*( gxx(i,j,k)*(vx(i,j,k)+betx(i,j,k)) &
                         + gxy(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gxz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    u_y = u0*ephi*( gxy(i,j,k)*(vx(i,j,k)+betx(i,j,k)) &
                         + gyy(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gyz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    u_z = u0*ephi*( gxz(i,j,k)*(vx(i,j,k)+betx(i,j,k)) &
                         + gyz(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gzz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    sb_x = ephi*( gxx(i,j,k)*(sbx(i,j,k) &
                         +sbt(i,j,k)*betx(i,j,k)) &
                         + gxy(i,j,k)*(sby(i,j,k)+sbt(i,j,k)*bety(i,j,k)) + &
                         gxz(i,j,k)*(sbz(i,j,k)+sbt(i,j,k)*betz(i,j,k)) )
                    sb_y = ephi*( gxy(i,j,k)*(sbx(i,j,k) &
                         +sbt(i,j,k)*betx(i,j,k)) &
                         + gyy(i,j,k)*(sby(i,j,k)+sbt(i,j,k)*bety(i,j,k)) + &
                         gyz(i,j,k)*(sbz(i,j,k)+sbt(i,j,k)*betz(i,j,k)) )
                    sb_z = ephi*( gxz(i,j,k)*(sbx(i,j,k) &
                         +sbt(i,j,k)*betx(i,j,k)) &
                         + gyz(i,j,k)*(sby(i,j,k)+sbt(i,j,k)*bety(i,j,k)) + &
                         gzz(i,j,k)*(sbz(i,j,k)+sbt(i,j,k)*betz(i,j,k)) )
                    ephi = exp(6.d0 * phi(i,j,k))
                    rhos(i,j,k) = rhob(i,j,k)*ephi*au0
                    w(i,j,k) = au0*rhos(i,j,k)
                    enth(i,j,k) = 1.d0 + (1.d0+n)*P(i,j,k)/rhob(i,j,k)
                    tau(i,j,k) = w(i,j,k)*enth(i,j,k) - ephi*P(i,j,k) &
                         - rhos(i,j,k) &
                         + ephi*(sb2*au0**2 - 0.5d0*sb2 - (al*sbt(i,j,k))**2)
                    mhdstx(i,j,k) = al*ephi*((rhob(i,j,k)*enth(i,j,k)+sb2)*u0*u_x &
                         - sbt(i,j,k)*sb_x )
                    mhdsty(i,j,k) = al*ephi*((rhob(i,j,k)*enth(i,j,k)+sb2)*u0*u_y &
                         - sbt(i,j,k)*sb_y )
                    mhdstz(i,j,k) = al*ephi*((rhob(i,j,k)*enth(i,j,k)+sb2)*u0*u_z &
                         - sbt(i,j,k)*sb_z )
                    stx(i,j,k) = al*ephi*rhob(i,j,k)*enth(i,j,k)*u0*u_x
                    sty(i,j,k) = al*ephi*rhob(i,j,k)*enth(i,j,k)*u0*u_y
                    stz(i,j,k) = al*ephi*rhob(i,j,k)*enth(i,j,k)*u0*u_z
                    sbt(i,j,k) = sbt(i,j,k)/f1os4pi
                    sbx(i,j,k) = sbx(i,j,k)/f1os4pi
                    sby(i,j,k) = sby(i,j,k)/f1os4pi
                    sbz(i,j,k) = sbz(i,j,k)/f1os4pi
                 endif
              enddo
           enddo
        enddo
     else
!!$        j = 2
        do j = 1, 3, 1
        do k = kmin, kmax-1, 1
              do i = imin, imax-1, 1
              if( dm(i,j,k) == onbd ) then
                 ni = 2
                 flax = 1
                 flaz = 1
                 if(x(i,1,1) < zeo) flax = -1
                 if(z(1,1,k) < zeo) flaz = -1
                 sf(1) = x(i,1,1) * flax
                 sf(3) = z(1,1,k) * flaz
                 if (k==kmin) flaz = max(0,flaz)
                 if( dm(i+flax,j     ,k     ) == onbd ) sf(1) = zeo
                 if( dm(i     ,j     ,k+flaz) == onbd ) sf(3) = zeo
!!$                 if(dm(i+flax,j,k) == outEZ) sf(3) = zeo
                 if((dm(i+flax,j,k) == onbd).and.(dm(i,j,k+flaz) == onbd)) then
                    ni = 1
                    if( dm(i+flax,j,k+flaz) /= outEZ .and. flaz .ne. 0) then
                         write(*,*) 'something wrong in find_excision_zone'
                         stop
                    end if
                 endif
                 if( ni == 2 ) then
                    rr = sf(1) + sf(3)
                    sf(1) = sf(1) / rr
                    sf(3) = sf(3) / rr
                    if(extrap_method==FLAT) then
                       r1 =   (sf(1)*x(i+flax,1,1)+(1.d0-sf(1))*x(i,1,1))**2 &
                            + y(1,j,1)**2 &
                            + (sf(3)*z(1,1,k+flaz)+(1.d0-sf(3))*z(1,1,k))**2
                       r2 = x(i,1,1)**2 + y(1,j,1)**2 + z(1,1,k)**2
                       fac1 = r1/r2
                       fac2 = 1.d0 - &
                            2.d0*(sqrt(r1)-sqrt(r2))/(sqrt(r1)+sqrt(r2))
                       fac3 = 1.d0 + &
                            2.d0*(sqrt(r1)-sqrt(r2))/(sqrt(r1)+sqrt(r2))
                       rhob(i,j,k) = fac1 * ( &
                            sf(1)*rhob(i+flax,j,k) + sf(3)*rhob(i,j,k+flaz) )
                       P(i,j,k) = fac1 * ( &
                            sf(1)*P(i+flax,j,k) + sf(3)*P(i,j,k+flaz) )
                       vr(i,j,k) = fac3 * ( &
                            sf(1)*vr(i+flax,j,k) + sf(3)*vr(i,j,k+flaz) )
                       vt(i,j,k) = fac2 * ( &
                            sf(1)*vt(i+flax,j,k) + sf(3)*vt(i,j,k+flaz) )
                       vp(i,j,k) = fac2 * ( &
                            sf(1)*vp(i+flax,j,k) + sf(3)*vp(i,j,k+flaz) )
                       Br(i,j,k) = fac1 * ( &
                            sf(1)*Br(i+flax,j,k) + sf(3)*Br(i,j,k+flaz) )
                       Bt(i,j,k) = fac2 * ( &
                            sf(1)*Bt(i+flax,j,k) + sf(3)*Bt(i,j,k+flaz) )
                       Bp(i,j,k) = fac2 * ( &
                            sf(1)*Bp(i+flax,j,k) + sf(3)*Bp(i,j,k+flaz) )
                    else if(extrap_method==LINEAR) then
                       rhob(i,j,k) = 1.d0/fr2(x(i,1,1),y(1,j,1),z(1,1,k)) * ( &
                            sf(1)* (2.d0*fr2(x(i+flax,1,1),y(1,j,1),z(1,1,k))*rhob(i+flax,j,k) &
                                       - fr2(x(i+2*flax,1,1),y(1,j,1),z(1,1,k))*rhob(i+2*flax,j,k)) &
                          + sf(3)* (2.d0*fr2(x(i,1,1),y(1,j,1),z(1,1,k+flaz))*rhob(i,j,k+flaz) &
                                       - fr2(x(i,1,1),y(1,j,1),z(1,1,k+2*flaz))*rhob(i,j,k+2*flaz)) )
                       P(i,j,k) = 1.d0/fr2(x(i,1,1),y(1,j,1),z(1,1,k)) * ( &
                            sf(1)* (2.d0*fr2(x(i+flax,1,1),y(1,j,1),z(1,1,k))*P(i+flax,j,k) &
                                       - fr2(x(i+2*flax,1,1),y(1,j,1),z(1,1,k))*P(i+2*flax,j,k)) &
                          + sf(3)* (2.d0*fr2(x(i,1,1),y(1,j,1),z(1,1,k+flaz))*P(i,j,k+flaz) &
                                       - fr2(x(i,1,1),y(1,j,1),z(1,1,k+2*flaz))*P(i,j,k+2*flaz)) )
                       vr(i,j,k) = 1.d0/fr1(x(i,1,1),y(1,j,1),z(1,1,k)) * ( &
                            sf(1)* (2.d0*fr1(x(i+flax,1,1),y(1,j,1),z(1,1,k))*vr(i+flax,j,k) &
                                       - fr1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k))*vr(i+2*flax,j,k)) &
                          + sf(3)* (2.d0*fr1(x(i,1,1),y(1,j,1),z(1,1,k+flaz))*vr(i,j,k+flaz) &
                                       - fr1(x(i,1,1),y(1,j,1),z(1,1,k+2*flaz))*vr(i,j,k+2*flaz)) )
                       vt(i,j,k) = 1.d0/frm1(x(i,1,1),y(1,j,1),z(1,1,k)) * ( &
                            sf(1)* (2.d0*frm1(x(i+flax,1,1),y(1,j,1),z(1,1,k))*vt(i+flax,j,k) &
                                       - frm1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k))*vt(i+2*flax,j,k)) &
                          + sf(3)* (2.d0*frm1(x(i,1,1),y(1,j,1),z(1,1,k+flaz))*vt(i,j,k+flaz) &
                                       - frm1(x(i,1,1),y(1,j,1),z(1,1,k+2*flaz))*vt(i,j,k+2*flaz)) )
                       vp(i,j,k) = 1.d0/frm1(x(i,1,1),y(1,j,1),z(1,1,k)) * ( &
                            sf(1)* (2.d0*frm1(x(i+flax,1,1),y(1,j,1),z(1,1,k))*vp(i+flax,j,k) &
                                       - frm1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k))*vp(i+2*flax,j,k)) &
                          + sf(3)* (2.d0*frm1(x(i,1,1),y(1,j,1),z(1,1,k+flaz))*vp(i,j,k+flaz) &
                                       - frm1(x(i,1,1),y(1,j,1),z(1,1,k+2*flaz))*vp(i,j,k+2*flaz)) )
                       Br(i,j,k) = 1.d0/fr2(x(i,1,1),y(1,j,1),z(1,1,k)) * ( &
                            sf(1)* (2.d0*fr2(x(i+flax,1,1),y(1,j,1),z(1,1,k))*Br(i+flax,j,k) &
                                       - fr2(x(i+2*flax,1,1),y(1,j,1),z(1,1,k))*Br(i+2*flax,j,k)) &
                          + sf(3)* (2.d0*fr2(x(i,1,1),y(1,j,1),z(1,1,k+flaz))*Br(i,j,k+flaz) &
                                       - fr2(x(i,1,1),y(1,j,1),z(1,1,k+2*flaz))*Br(i,j,k+2*flaz)) )
                       Bt(i,j,k) = 1.d0/frm1(x(i,1,1),y(1,j,1),z(1,1,k)) * ( &
                            sf(1)* (2.d0*frm1(x(i+flax,1,1),y(1,j,1),z(1,1,k))*Bt(i+flax,j,k) &
                                       - frm1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k))*Bt(i+2*flax,j,k)) &
                          + sf(3)* (2.d0*frm1(x(i,1,1),y(1,j,1),z(1,1,k+flaz))*Bt(i,j,k+flaz) &
                                       - frm1(x(i,1,1),y(1,j,1),z(1,1,k+2*flaz))*Bt(i,j,k+2*flaz)) )
                       Bp(i,j,k) = 1.d0/frm1(x(i,1,1),y(1,j,1),z(1,1,k)) * ( &
                            sf(1)* (2.d0*frm1(x(i+flax,1,1),y(1,j,1),z(1,1,k))*Bp(i+flax,j,k) &
                                       - frm1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k))*Bp(i+2*flax,j,k)) &
                          + sf(3)* (2.d0*frm1(x(i,1,1),y(1,j,1),z(1,1,k+flaz))*Bp(i,j,k+flaz) &
                                       - frm1(x(i,1,1),y(1,j,1),z(1,1,k+2*flaz))*Bp(i,j,k+2*flaz)) )
                    end if
                 else if( ni == 1 ) then
                    if(extrap_method==FLAT) then
                       r1 = x(i+flax,1,1)**2 + y(1,j,1)**2 + z(1,1,k+flaz)**2
                       r2 = x(i,1,1)**2 + y(1,j,1)**2 + z(1,1,k)**2
                       fac1 = r1/r2
                       fac2 = 1.d0 - 2.d0*(sqrt(r1)-sqrt(r2))/(sqrt(r1)+sqrt(r2))
                       fac3 = 1.d0 + 2.d0*(sqrt(r1)-sqrt(r2))/(sqrt(r1)+sqrt(r2))
                       rhob(i,j,k) = fac1 * rhob(i+flax,j,k+flaz)
                       P(i,j,k)    = fac1 * P(i+flax,j,k+flaz)
                       vr(i,j,k)   = fac3 * vr(i+flax,j,k+flaz)
                       vt(i,j,k)   = fac2 * vt(i+flax,j,k+flaz)
                       vp(i,j,k)   = fac2 * vp(i+flax,j,k+flaz)
                       Br(i,j,k)   = fac1 * Br(i+flax,j,k+flaz)
                       Bt(i,j,k)   = fac2 * Bt(i+flax,j,k+flaz)
                       Bp(i,j,k)   = fac2 * Bp(i+flax,j,k+flaz)
                    else if(extrap_method==LINEAR) then
                       rhob(i,j,k) = 1.d0/fr2(x(i,1,1),y(1,j,1),z(1,1,k)) * &
                            (2.d0*fr2(x(i+flax,1,1),y(1,j,1),z(1,1,k+flaz))*rhob(i+flax,j,k+flaz) &
                                - fr2(x(i+2*flax,1,1),y(1,j,1),z(1,1,k+2*flaz))*rhob(i+2*flax,j,k+2*flaz))
                       P(i,j,k) = 1.d0/fr2(x(i,1,1),y(1,j,1),z(1,1,k)) * &
                            (2.d0*fr2(x(i+flax,1,1),y(1,j,1),z(1,1,k+flaz))*P(i+flax,j,k+flaz) &
                                - fr2(x(i+2*flax,1,1),y(1,j,1),z(1,1,k+2*flaz))*P(i+2*flax,j,k+2*flaz))
                       vr(i,j,k) = 1.d0/fr1(x(i,1,1),y(1,j,1),z(1,1,k)) * &
                            (2.d0*fr1(x(i+flax,1,1),y(1,j,1),z(1,1,k+flaz))*vr(i+flax,j,k+flaz) &
                                - fr1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k+2*flaz))*vr(i+2*flax,j,k+2*flaz))
                       vt(i,j,k) = 1.d0/frm1(x(i,1,1),y(1,j,1),z(1,1,k)) * &
                            (2.d0*frm1(x(i+flax,1,1),y(1,j,1),z(1,1,k+flaz))*vt(i+flax,j,k+flaz) &
                                - frm1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k+2*flaz))*vt(i+2*flax,j,k+2*flaz))
                       vp(i,j,k) = 1.d0/frm1(x(i,1,1),y(1,j,1),z(1,1,k)) * &
                            (2.d0*frm1(x(i+flax,1,1),y(1,j,1),z(1,1,k+flaz))*vp(i+flax,j,k+flaz) &
                                - frm1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k+2*flaz))*vp(i+2*flax,j,k+2*flaz))
                       Br(i,j,k) = 1.d0/fr2(x(i,1,1),y(1,j,1),z(1,1,k)) * &
                            (2.d0*fr2(x(i+flax,1,1),y(1,j,1),z(1,1,k+flaz))*Br(i+flax,j,k+flaz) &
                                - fr2(x(i+2*flax,1,1),y(1,j,1),z(1,1,k+2*flaz))*Br(i+2*flax,j,k+2*flaz))
                       Bt(i,j,k) = 1.d0/frm1(x(i,1,1),y(1,j,1),z(1,1,k)) * &
                            (2.d0*frm1(x(i+flax,1,1),y(1,j,1),z(1,1,k+flaz))*Bt(i+flax,j,k+flaz) &
                                - frm1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k+2*flaz))*Bt(i+2*flax,j,k+2*flaz))
                       Bp(i,j,k) = 1.d0/frm1(x(i,1,1),y(1,j,1),z(1,1,k)) * &
                            (2.d0*frm1(x(i+flax,1,1),y(1,j,1),z(1,1,k+flaz))*Bp(i+flax,j,k+flaz) &
                                - frm1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k+2*flaz))*Bp(i+2*flax,j,k+2*flaz))
                    end if
              endif
!!$              P(i,j,k) = 0.439377566268456810*rhob(i,j,k)**(1.d0 + 1.d0/n)
              P(i,j,k) = kpoly_init*rhob(i,j,k)**(1.d0 + 1.d0/n)
              r = sqrt(x(i,1,1)**2 + y(1,j,1)**2 + z(1,1,k)**2)
              costhet = z(1,1,k)/r
              sinthet = sqrt(1.d0 - costhet**2)
              cosphi  = x(i,1,1)/(r*sinthet)
              sinphi  = y(1,j,1)/(r*sinthet)
              vx(i,j,k) = sinthet*cosphi*vr(i,j,k) + &
                   r*costhet*cosphi*vt(i,j,k) - r*sinthet*sinphi*vp(i,j,k)
              vy(i,j,k) = sinthet*sinphi*vr(i,j,k) + &
                   r*costhet*sinphi*vt(i,j,k) + r*sinthet*cosphi*vp(i,j,k)
              vz(i,j,k) = costhet*vr(i,j,k) - r*sinthet*vt(i,j,k)
              Bx(i,j,k) = sinthet*cosphi*Br(i,j,k) + &
                   r*costhet*cosphi*Bt(i,j,k) - r*sinthet*sinphi*Bp(i,j,k)
              By(i,j,k) = sinthet*sinphi*Br(i,j,k) + &
                   r*costhet*sinphi*Bt(i,j,k) + r*sinthet*cosphi*Bp(i,j,k)
              Bz(i,j,k) = costhet*Br(i,j,k) - r*sinthet*Bt(i,j,k)
                 al = 1.d0 + alp(i,j,k)
                    ephi = exp(4.d0 * phi(i,j,k))
                    ! Compute al*u0-1
                    v2 = ephi*(gxx(i,j,k)*(vx(i,j,k) + betx(i,j,k))**2 + &
                         2.d0*gxy(i,j,k)*(vx(i,j,k) + betx(i,j,k)) &
                         *(vy(i,j,k) + bety(i,j,k)) +  &
                         2.d0*gxz(i,j,k)*(vx(i,j,k) + betx(i,j,k)) &
                         *(vz(i,j,k) + betz(i,j,k)) +  &
                         gyy(i,j,k)*(vy(i,j,k) + bety(i,j,k))**2 + &
                         2.d0*gyz(i,j,k)*(vy(i,j,k) + bety(i,j,k)) &
                         *(vz(i,j,k) + betz(i,j,k)) +  &
                         gzz(i,j,k)*(vz(i,j,k) + betz(i,j,k))**2 )/al**2
                    ! *** Check for superluminal velocity ***
                    if (v2 .gt. 1.d0) then
                       vx(i,j,k)=(vx(i,j,k)+betx(i,j,k))*sqrt(0.99d0/v2)-betx(i,j,k)
                       vy(i,j,k)=(vy(i,j,k)+bety(i,j,k))*sqrt(0.99d0/v2)-bety(i,j,k)
                       vz(i,j,k) =(vz(i,j,k)+betz(i,j,k))*sqrt(0.99d0/v2)-betz(i,j,k)
                       v2 = 0.99d0
                    end if
                    ! ***************************************
                    gam = sqrt(1.d0-v2)
                    au0 = 1.d0/gam
                    u0  = au0/al
                    sbt(i,j,k) = f1os4pi*ephi*u0/al*( &
                         gxx(i,j,k)*Bx(i,j,k)*(vx(i,j,k)+betx(i,j,k)) + &
                         gxy(i,j,k)*( Bx(i,j,k)*(vy(i,j,k)+bety(i,j,k)) &
                         + By(i,j,k)*(vx(i,j,k)+betx(i,j,k)) ) + &
                         gxz(i,j,k)*( Bx(i,j,k)*(vz(i,j,k)+betz(i,j,k)) &
                         + Bz(i,j,k)*(vx(i,j,k)+betx(i,j,k)) ) + &
                         gyy(i,j,k)*By(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gyz(i,j,k)*( By(i,j,k)*(vz(i,j,k)+betz(i,j,k)) &
                         + Bz(i,j,k)*(vy(i,j,k)+bety(i,j,k)) ) + &
                         gzz(i,j,k)*Bz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    sbx(i,j,k) = f1os4pi*Bx(i,j,k)/(al*u0) &
                         + sbt(i,j,k)*vx(i,j,k)
                    sby(i,j,k) = f1os4pi*By(i,j,k)/(al*u0) &
                         + sbt(i,j,k)*vy(i,j,k)
                    sbz(i,j,k) = f1os4pi*Bz(i,j,k)/(al*u0) &
                         + sbt(i,j,k)*vz(i,j,k)
                    sb2 = -(al*sbt(i,j,k))**2 &
                         + ephi*( gxx(i,j,k)*(&
                         sbx(i,j,k)+betx(i,j,k)*sbt(i,j,k))**2 + &
                         2.d0*gxy(i,j,k)*(sbx(i,j,k)+betx(i,j,k)*sbt(i,j,k)) &
                         *(sby(i,j,k)+bety(i,j,k)*sbt(i,j,k)) + &
                         2.d0*gxz(i,j,k)*(sbx(i,j,k)+betx(i,j,k)*sbt(i,j,k)) &
                         *(sbz(i,j,k)+betz(i,j,k)*sbt(i,j,k)) + &
                         gyy(i,j,k)*(sby(i,j,k)+bety(i,j,k)*sbt(i,j,k))**2 + &
                         2.d0*gyz(i,j,k)*(sby(i,j,k)+bety(i,j,k)*sbt(i,j,k)) &
                         *(sbz(i,j,k)+betz(i,j,k)*sbt(i,j,k)) + &
                         gzz(i,j,k)*(sbz(i,j,k)+betz(i,j,k)*sbt(i,j,k))**2 )
                    u_x = u0*ephi*( gxx(i,j,k)*(vx(i,j,k)+betx(i,j,k)) &
                         + gxy(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gxz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    u_y = u0*ephi*( gxy(i,j,k)*(vx(i,j,k)+betx(i,j,k)) &
                         + gyy(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gyz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    u_z = u0*ephi*( gxz(i,j,k)*(vx(i,j,k)+betx(i,j,k)) &
                         + gyz(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gzz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    sb_x = ephi*( gxx(i,j,k)*(sbx(i,j,k) &
                         +sbt(i,j,k)*betx(i,j,k)) &
                         + gxy(i,j,k)*(sby(i,j,k)+sbt(i,j,k)*bety(i,j,k)) + &
                         gxz(i,j,k)*(sbz(i,j,k)+sbt(i,j,k)*betz(i,j,k)) )
                    sb_y = ephi*( gxy(i,j,k)*(sbx(i,j,k) &
                         +sbt(i,j,k)*betx(i,j,k)) &
                         + gyy(i,j,k)*(sby(i,j,k)+sbt(i,j,k)*bety(i,j,k)) + &
                         gyz(i,j,k)*(sbz(i,j,k)+sbt(i,j,k)*betz(i,j,k)) )
                    sb_z = ephi*( gxz(i,j,k)*(sbx(i,j,k) &
                         +sbt(i,j,k)*betx(i,j,k)) &
                         + gyz(i,j,k)*(sby(i,j,k)+sbt(i,j,k)*bety(i,j,k)) + &
                         gzz(i,j,k)*(sbz(i,j,k)+sbt(i,j,k)*betz(i,j,k)) )
                    ephi = exp(6.d0 * phi(i,j,k))
                    rhos(i,j,k) = rhob(i,j,k)*ephi*au0
                    w(i,j,k) = au0*rhos(i,j,k)
                    enth(i,j,k) = 1.d0 + (1.d0+n)*P(i,j,k)/rhob(i,j,k)
                    tau(i,j,k) = w(i,j,k)*enth(i,j,k) - ephi*P(i,j,k) &
                         - rhos(i,j,k) &
                         + ephi*(sb2*au0**2 - 0.5d0*sb2 - (al*sbt(i,j,k))**2)
                    mhdstx(i,j,k) = al*ephi*((rhob(i,j,k)*enth(i,j,k)+sb2)*u0*u_x &
                         - sbt(i,j,k)*sb_x )
                    mhdsty(i,j,k) = al*ephi*((rhob(i,j,k)*enth(i,j,k)+sb2)*u0*u_y &
                         - sbt(i,j,k)*sb_y )
                    mhdstz(i,j,k) = al*ephi*((rhob(i,j,k)*enth(i,j,k)+sb2)*u0*u_z &
                         - sbt(i,j,k)*sb_z )
                    stx(i,j,k) = al*ephi*rhob(i,j,k)*enth(i,j,k)*u0*u_x
                    sty(i,j,k) = al*ephi*rhob(i,j,k)*enth(i,j,k)*u0*u_y
                    stz(i,j,k) = al*ephi*rhob(i,j,k)*enth(i,j,k)*u0*u_z
                    sbt(i,j,k) = sbt(i,j,k)/f1os4pi*al
                    sbx(i,j,k) = sbx(i,j,k)/f1os4pi*al
                    sby(i,j,k) = sby(i,j,k)/f1os4pi*al
                    sbz(i,j,k) = sbz(i,j,k)/f1os4pi*al
                 end if
           end do
        end do
        end do
  end if
  return
end subroutine hydro_ezbc
!-------------------------------------------------------------------------+
!
! Use Alcubierres simple excision method for fluid (for hybrid EOS)
!
! excision region : sphere
!
!-------------------------------------------------------------------------+
subroutine hydro_ezbc_hybrid(ex,x,y,z,rhos,tau,mhdstx,mhdsty,mhdstz,stx,sty,stz, &
     rhob,P,enth,vx,vy,vz,w,sbt,sbx,sby,sbz,Bx,By,Bz,alp,betx,bety,betz,phi, &
     gxx,gxy,gxz,gyy,gyz,gzz,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
     symmetry,dm,gamma_th,neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
  implicit none
!~~~~~~> Input arguments
  integer, intent(in) :: ex(1:3),symmetry
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   ::x,y,z
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::rhos,tau,mhdstx,mhdsty,mhdstz
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::stx,sty,stz
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::rhob,P,enth,vx,vy,vz,w
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::sbt,sbx,sby,sbz
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::Bx,By,Bz
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   ::alp,betx,bety,betz,phi
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   ::gxx,gxy,gxz,gyy,gyz,gzz
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   ::gupxx,gupxy,gupxz
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   ::gupyy,gupyz,gupzz
  real*8    :: gamma_th,ergo_sigma
  integer :: neos,ergo_star
  real*8, dimension(neos) :: rho_tab,P_tab,eps_tab
  real*8, dimension(neos+1) :: k_tab, gamma_tab
  integer :: dm(1:ex(1),1:ex(2),1:ex(3))
!~~~~~~> local variable
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k,flax,flay,flaz,ni
  integer :: ksym,jsym,isym,diax,diay,diaz
  real*8  :: sf(1:3),rr,xs,zs,rxs,rzs,rxzs,c1,c2
  real*8  :: al,ephi,u0,au0,u_x,u_y,u_z,sb_x,sb_y,sb_z,sb2,v2,gam,f1os4pi
  real*8  :: r1,r2,fac1,fac2,fac3,r,sinthet,costhet,sinphi,cosphi,Pmin
  integer :: extrap_method
  real*8  :: P_cold, eps_cold, eps
  integer, parameter :: inEZ=0, onbd=1, outEZ=2
  real*8,  parameter :: zeo=0.d0, one=1.d0, two=2.d0
  integer, parameter :: nosymmetry=0, equatorial=1, octant=2
  integer, parameter :: pi_symm=3, axisym=4
  integer, parameter :: FLAT=0, LINEAR=1
  real*8,dimension(ex(1),ex(2),ex(3))   :: vr,vt,vp,Br,Bt,Bp
!~~~~~~> interface
  interface
     !subroutine find_excision_zone(ex,x,y,z,dm,innerbd,Symmetry)
     !  implicit none
     !  integer, intent(in)  :: ex(1:3)
     !  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   ::x,y,z 
     !  integer, intent(out) :: dm(1:ex(1),1:ex(2),1:ex(3))
     !  real*8               :: innerbd
     !  integer, intent(in)  :: Symmetry
     !end subroutine find_excision_zone
     subroutine compute_pcold_epscold(rhob, P_cold, eps_cold, &
                neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
       implicit none
       integer :: neos, ergo_star
       real*8  :: rhob, P_cold, eps_cold, ergo_sigma
       real*8, dimension(neos) :: rho_tab,P_tab,eps_tab
       real*8, dimension(neos+1) :: k_tab, gamma_tab
     end subroutine compute_pcold_epscold
  end interface
  real*8  :: fr2,fr1,frm1
  EXTERNAL fr2,fr1,frm1
  extrap_method = FLAT
!  extrap_method = LINEAR
  imin = lbound(rhos,1)
  jmin = lbound(rhos,2)
  kmin = lbound(rhos,3)
  imax = ubound(rhos,1)
  jmax = ubound(rhos,2)
  kmax = ubound(rhos,3)
 f1os4pi = 1.d0/sqrt(4.d0*acos(-1.d0))
!!$  c1 = 2.d0
  c1 = 1.d0
  c2 = c1 - 1.d0
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
!~~~~~~> Find the grid points on the boundary
!  call find_excision_zone(ex,x,y,z,dm,innerbd,Symmetry)
!~~~~~~> inner boundary condition
  if( symmetry == nosymmetry ) then
     isym = imin + 1
     jsym = jmin + 1
     ksym = kmin + 1
  else if( symmetry == equatorial ) then
     isym = imin + 1
     jsym = jmin + 1
     ksym = kmin
  else if( symmetry == octant ) then
     isym = imin
     jsym = jmin
     ksym = kmin
  else if( symmetry == pi_symm ) then
     isym = imin + 1
     jsym = jmin
     ksym = kmin
  endif
  do k = kmin, kmax, 1
     do j = jmin, jmax, 1
        do i = imin, imax, 1
           r = sqrt(x(i,1,1)**2 + y(1,j,1)**2 + z(1,1,k)**2)
           costhet = z(1,1,k)/r
           sinthet = sqrt(1.d0 - costhet**2)
           cosphi  = x(i,1,1)/(r*sinthet)
           sinphi  = y(1,j,1)/(r*sinthet)
           vr(i,j,k) = sinthet*cosphi*vx(i,j,k) + &
                sinthet*sinphi*vy(i,j,k) + costhet*vz(i,j,k)
           vt(i,j,k) = costhet*cosphi/r*vx(i,j,k) + &
                costhet*sinphi/r*vy(i,j,k) - sinthet/r*vz(i,j,k)
           vp(i,j,k) = -sinphi/sinthet/r*vx(i,j,k) + &
                cosphi/sinthet/r*vy(i,j,k)
           Br(i,j,k) = sinthet*cosphi*Bx(i,j,k) + &
                sinthet*sinphi*By(i,j,k) + costhet*Bz(i,j,k)
           Bt(i,j,k) = costhet*cosphi/r*Bx(i,j,k) + &
                costhet*sinphi/r*By(i,j,k) - sinthet/r*Bz(i,j,k)
           Bp(i,j,k) = -sinphi/sinthet/r*Bx(i,j,k) + &
                cosphi/sinthet/r*By(i,j,k)
        end do
     end do
  end do
  if(symmetry .ne. axisym) then
     do k = ksym, kmax-1, 1
        do j = jsym, jmax-1, 1
           do i = isym, imax-1, 1
              if( dm(i,j,k) == onbd ) then
                 ni = 2
                 flax = 1
                 flay = 1
                 flaz = 1
                 if(x(i,1,1) < zeo) flax = -1
                 if(y(1,j,1) < zeo) flay = -1
                 if(z(1,1,k) < zeo) flaz = -1
                 sf(1) = x(i,1,1) * flax
                 sf(2) = y(1,j,1) * flay
                 sf(3) = z(1,1,k) * flaz
                 if( dm(i+flax,j     ,k     ) == onbd ) sf(1) = zeo
                 if( dm(i     ,j+flay,k     ) == onbd ) sf(2) = zeo
                 if( dm(i     ,j     ,k+flaz) == onbd ) sf(3) = zeo
                 if((dm(i+flax,j,k) == outEZ) .or. (dm(i,j+flay,k) == outEZ)) &
                      sf(3) = zeo
                 if(( dm(i+flax,j     ,k     ) == onbd ) .and. &
                      ( dm(i     ,j+flay,k     ) == onbd ) .and. &
                      ( dm(i     ,j     ,k+flaz) == onbd )) then
                    ni = 1
                    if( dm(i     ,j+flay,k+flaz) == outEZ ) sf(2) = one
                    if( dm(i+flax,j     ,k+flaz) == outEZ ) sf(3) = one
                    if( dm(i+flax,j+flay,k     ) == outEZ ) then
                       sf(1) = one
                       sf(2) = zeo
                       sf(3) = zeo
                    endif
                    if(( dm(i+flax,j+flay,k) /= outEZ ) .and. &
                         ( dm(i,j+flay,k+flaz) /= outEZ ) .and. &
                         ( dm(i+flax,j,k+flaz) /= outEZ )) then
                         write(*,*) 'something wrong in find_excision_zone'
                    end if
                 endif
                 if( ni == 2 ) then
                    rr = sf(1) + sf(2) + sf(3)
                    sf(1) = sf(1) / rr
                    sf(2) = sf(2) / rr
                    sf(3) = sf(3) / rr
                    rhob(i,j,k) = sf(1) * rhob(i+flax,j,k) + &
                         sf(2) * rhob(i,j+flay,k) + &
                         sf(3) * rhob(i,j,k+flaz)
                    P(i,j,k) = sf(1) * P(i+flax,j,k) + &
                         sf(2) * P(i,j+flay,k) + &
                         sf(3) * P(i,j,k+flaz)
                    vx(i,j,k) = sf(1) * vx(i+flax,j,k) + &
                         sf(2) * vx(i,j+flay,k) + &
                         sf(3) * vx(i,j,k+flaz)
                    vy(i,j,k) = sf(1) * vy(i+flax,j,k) + &
                         sf(2) * vy(i,j+flay,k) + &
                         sf(3) * vy(i,j,k+flaz)
                    vz(i,j,k) = sf(1) * vz(i+flax,j,k) + &
                         sf(2) * vz(i,j+flay,k) + &
                         sf(3) * vz(i,j,k+flaz)
                    Bx(i,j,k) = sf(1) * Bx(i+flax,j,k) + &
                         sf(2) * Bx(i,j+flay,k) + &
                         sf(3) * Bx(i,j,k+flaz)
                    By(i,j,k) = sf(1) * By(i+flax,j,k) + &
                         sf(2) * By(i,j+flay,k) + &
                         sf(3) * By(i,j,k+flaz)
                    Bz(i,j,k) = sf(1) * Bz(i+flax,j,k) + &
                         sf(2) * Bz(i,j+flay,k) + &
                         sf(3) * Bz(i,j,k+flaz)
                 else if( ni == 1 ) then
                    rr = sf(1) + sf(2) +sf(3)
                    sf(1) = sf(1) / rr
                    sf(2) = sf(2) / rr
                    sf(3) = sf(3) / rr
                    rhob(i,j,k) = sf(1) * rhob(i+flax,j+flay,k) + &
                         sf(2) * rhob(i,j+flay,k+flaz) + &
                         sf(3) * rhob(i+flax,j,k+flaz)
                    P(i,j,k) = sf(1) * P(i+flax,j+flay,k) + &
                         sf(2) * P(i,j+flay,k+flaz) + &
                         sf(3) * P(i+flax,j,k+flaz)
                    vx(i,j,k) = sf(1) * vx(i+flax,j+flay,k) + &
                         sf(2) * vx(i,j+flay,k+flaz) + &
                         sf(3) * vx(i+flax,j,k+flaz)
                    vy(i,j,k) = sf(1) * vy(i+flax,j+flay,k) + &
                         sf(2) * vy(i,j+flay,k+flaz) + &
                         sf(3) * vy(i+flax,j,k+flaz)
                    vz(i,j,k) = sf(1) * vz(i+flax,j+flay,k) + &
                         sf(2) * vz(i,j+flay,k+flaz) + &
                         sf(3) * vz(i+flax,j,k+flaz)
                    Bx(i,j,k) = sf(1) * Bx(i+flax,j+flay,k) + &
                         sf(2) * Bx(i,j+flay,k+flaz) + &
                         sf(3) * Bx(i+flax,j,k+flaz)
                    By(i,j,k) = sf(1) * By(i+flax,j+flay,k) + &
                         sf(2) * By(i,j+flay,k+flaz) + &
                         sf(3) * By(i+flax,j,k+flaz)
                    Bz(i,j,k) = sf(1) * Bz(i+flax,j+flay,k) + &
                         sf(2) * Bz(i,j+flay,k+flaz) + &
                         sf(3) * Bz(i+flax,j,k+flaz)
                 endif
                    al = 1.d0 + alp(i,j,k)
                    ephi = exp(4.d0 * phi(i,j,k))
                    ! Compute al*u0-1
                    v2 = ephi*(gxx(i,j,k)*(vx(i,j,k) + betx(i,j,k))**2 + &
                         2.d0*gxy(i,j,k)*(vx(i,j,k) + betx(i,j,k)) &
                         *(vy(i,j,k) + bety(i,j,k)) +  &
                         2.d0*gxz(i,j,k)*(vx(i,j,k) + betx(i,j,k)) &
                         *(vz(i,j,k) + betz(i,j,k)) +  &
                         gyy(i,j,k)*(vy(i,j,k) + bety(i,j,k))**2 + &
                         2.d0*gyz(i,j,k)*(vy(i,j,k) + bety(i,j,k)) &
                         *(vz(i,j,k) + betz(i,j,k)) +  &
                         gzz(i,j,k)*(vz(i,j,k) + betz(i,j,k))**2 )/al**2
                    ! *** Check for superluminal velocity ***
                    if (v2 .gt. 1.d0) then
                       vx(i,j,k)=(vx(i,j,k)+betx(i,j,k))*sqrt(0.99d0/v2)-betx(i,j,k)
                       vy(i,j,k)=(vy(i,j,k)+bety(i,j,k))*sqrt(0.99d0/v2)-bety(i,j,k)
                       vz(i,j,k)=(vz(i,j,k)+betz(i,j,k))*sqrt(0.99d0/v2)-betz(i,j,k)
                       v2 = 0.99d0
                    end if
                    ! ***************************************
                    gam = sqrt(1.d0-v2)
                    au0 = 1.d0/gam
                    u0  = au0/al
                    sbt(i,j,k) = f1os4pi*ephi*u0/al*( &
                         gxx(i,j,k)*Bx(i,j,k)*(vx(i,j,k)+betx(i,j,k)) + &
                         gxy(i,j,k)*( Bx(i,j,k)*(vy(i,j,k)+bety(i,j,k)) &
                         + By(i,j,k)*(vx(i,j,k)+betx(i,j,k)) ) + &
                         gxz(i,j,k)*( Bx(i,j,k)*(vz(i,j,k)+betz(i,j,k)) &
                         + Bz(i,j,k)*(vx(i,j,k)+betx(i,j,k)) ) + &
                         gyy(i,j,k)*By(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gyz(i,j,k)*( By(i,j,k)*(vz(i,j,k)+betz(i,j,k)) &
                         + Bz(i,j,k)*(vy(i,j,k)+bety(i,j,k)) ) + &
                         gzz(i,j,k)*Bz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    sbx(i,j,k) = f1os4pi*Bx(i,j,k)/(al*u0) &
                         + sbt(i,j,k)*vx(i,j,k)
                    sby(i,j,k) = f1os4pi*By(i,j,k)/(al*u0) &
                         + sbt(i,j,k)*vy(i,j,k)
                    sbz(i,j,k) = f1os4pi*Bz(i,j,k)/(al*u0) &
                         + sbt(i,j,k)*vz(i,j,k)
                    sb2 = -(al*sbt(i,j,k))**2 &
                         + ephi*( gxx(i,j,k)*(&
                         sbx(i,j,k)+betx(i,j,k)*sbt(i,j,k))**2 + &
                         2.d0*gxy(i,j,k)*(sbx(i,j,k)+betx(i,j,k)*sbt(i,j,k)) &
                         *(sby(i,j,k)+bety(i,j,k)*sbt(i,j,k)) + &
                         2.d0*gxz(i,j,k)*(sbx(i,j,k)+betx(i,j,k)*sbt(i,j,k)) &
                         *(sbz(i,j,k)+betz(i,j,k)*sbt(i,j,k)) + &
                         gyy(i,j,k)*(sby(i,j,k)+bety(i,j,k)*sbt(i,j,k))**2 + &
                         2.d0*gyz(i,j,k)*(sby(i,j,k)+bety(i,j,k)*sbt(i,j,k)) &
                         *(sbz(i,j,k)+betz(i,j,k)*sbt(i,j,k)) + &
                         gzz(i,j,k)*(sbz(i,j,k)+betz(i,j,k)*sbt(i,j,k))**2 )
                    u_x = u0*ephi*( gxx(i,j,k)*(vx(i,j,k)+betx(i,j,k)) &
                         + gxy(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gxz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    u_y = u0*ephi*( gxy(i,j,k)*(vx(i,j,k)+betx(i,j,k)) &
                         + gyy(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gyz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    u_z = u0*ephi*( gxz(i,j,k)*(vx(i,j,k)+betx(i,j,k)) &
                         + gyz(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gzz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    sb_x = ephi*( gxx(i,j,k)*(sbx(i,j,k) &
                         +sbt(i,j,k)*betx(i,j,k)) &
                         + gxy(i,j,k)*(sby(i,j,k)+sbt(i,j,k)*bety(i,j,k)) + &
                         gxz(i,j,k)*(sbz(i,j,k)+sbt(i,j,k)*betz(i,j,k)) )
                    sb_y = ephi*( gxy(i,j,k)*(sbx(i,j,k) &
                         +sbt(i,j,k)*betx(i,j,k)) &
                         + gyy(i,j,k)*(sby(i,j,k)+sbt(i,j,k)*bety(i,j,k)) + &
                         gyz(i,j,k)*(sbz(i,j,k)+sbt(i,j,k)*betz(i,j,k)) )
                    sb_z = ephi*( gxz(i,j,k)*(sbx(i,j,k) &
                         +sbt(i,j,k)*betx(i,j,k)) &
                         + gyz(i,j,k)*(sby(i,j,k)+sbt(i,j,k)*bety(i,j,k)) + &
                         gzz(i,j,k)*(sbz(i,j,k)+sbt(i,j,k)*betz(i,j,k)) )
                    call compute_pcold_epscold(rhob(i,j,k), P_cold, eps_cold, &
                        neos, ergo_star, ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
write(*,*) "hi0.",i,j,k,rhob(i,j,k), P_cold, eps_cold,neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab
                    eps = eps_cold + (P(i,j,k)-P_cold)/(gamma_th-1.d0)/rhob(i,j,k)
                    ephi = exp(6.d0 * phi(i,j,k))
                    rhos(i,j,k) = rhob(i,j,k)*ephi*au0
                    w(i,j,k) = au0*rhos(i,j,k)
                    enth(i,j,k) = 1.d0 + P(i,j,k)/rhob(i,j,k) + eps
                    tau(i,j,k) = w(i,j,k)*enth(i,j,k) - ephi*P(i,j,k) &
                         - rhos(i,j,k) &
                         + ephi*(sb2*au0**2 - 0.5d0*sb2 - (al*sbt(i,j,k))**2)
                    mhdstx(i,j,k) = al*ephi*((rhob(i,j,k)*enth(i,j,k)+sb2)*u0*u_x &
                         - sbt(i,j,k)*sb_x )
                    mhdsty(i,j,k) = al*ephi*((rhob(i,j,k)*enth(i,j,k)+sb2)*u0*u_y &
                         - sbt(i,j,k)*sb_y )
                    mhdstz(i,j,k) = al*ephi*((rhob(i,j,k)*enth(i,j,k)+sb2)*u0*u_z &
                         - sbt(i,j,k)*sb_z )
                    stx(i,j,k) = al*ephi*rhob(i,j,k)*enth(i,j,k)*u0*u_x
                    sty(i,j,k) = al*ephi*rhob(i,j,k)*enth(i,j,k)*u0*u_y
                    stz(i,j,k) = al*ephi*rhob(i,j,k)*enth(i,j,k)*u0*u_z
                    sbt(i,j,k) = sbt(i,j,k)/f1os4pi
                    sbx(i,j,k) = sbx(i,j,k)/f1os4pi
                    sby(i,j,k) = sby(i,j,k)/f1os4pi
                    sbz(i,j,k) = sbz(i,j,k)/f1os4pi
                 endif
              enddo
           enddo
        enddo
     else
!!$        j = 2
        do j = 1, 3, 1
        do k = kmin, kmax-1, 1
              do i = imin, imax-1, 1
              if( dm(i,j,k) == onbd ) then
                 ni = 2
                 flax = 1
                 flaz = 1
                 if(x(i,1,1) < zeo) flax = -1
                 if(z(1,1,k) < zeo) flaz = -1
                 sf(1) = x(i,1,1) * flax
                 sf(3) = z(1,1,k) * flaz
                 if (k==kmin) flaz = max(0,flaz)
                 if( dm(i+flax,j     ,k     ) == onbd ) sf(1) = zeo
                 if( dm(i     ,j     ,k+flaz) == onbd ) sf(3) = zeo
!!$                 if(dm(i+flax,j,k) == outEZ) sf(3) = zeo
                 if((dm(i+flax,j,k) == onbd).and.(dm(i,j,k+flaz) == onbd)) then
                    ni = 1
                    if( dm(i+flax,j,k+flaz) /= outEZ .and. flaz .ne. 0) then
                         write(*,*) 'something wrong in find_excision_zone'
                         stop
                    end if
                 endif
                 if( ni == 2 ) then
                    rr = sf(1) + sf(3)
                    sf(1) = sf(1) / rr
                    sf(3) = sf(3) / rr
                    if(extrap_method==FLAT) then
                       r1 =   (sf(1)*x(i+flax,1,1)+(1.d0-sf(1))*x(i,1,1))**2 &
                            + y(1,j,1)**2 &
                            + (sf(3)*z(1,1,k+flaz)+(1.d0-sf(3))*z(1,1,k))**2
                       r2 = x(i,1,1)**2 + y(1,j,1)**2 + z(1,1,k)**2
                       fac1 = r1/r2
                       fac2 = 1.d0 - &
                            2.d0*(sqrt(r1)-sqrt(r2))/(sqrt(r1)+sqrt(r2))
                       fac3 = 1.d0 + &
                            2.d0*(sqrt(r1)-sqrt(r2))/(sqrt(r1)+sqrt(r2))
                       rhob(i,j,k) = fac1 * ( &
                            sf(1)*rhob(i+flax,j,k) + sf(3)*rhob(i,j,k+flaz) )
                       P(i,j,k) = fac1 * ( &
                            sf(1)*P(i+flax,j,k) + sf(3)*P(i,j,k+flaz) )
                       vr(i,j,k) = fac3 * ( &
                            sf(1)*vr(i+flax,j,k) + sf(3)*vr(i,j,k+flaz) )
                       vt(i,j,k) = fac2 * ( &
                            sf(1)*vt(i+flax,j,k) + sf(3)*vt(i,j,k+flaz) )
                       vp(i,j,k) = fac2 * ( &
                            sf(1)*vp(i+flax,j,k) + sf(3)*vp(i,j,k+flaz) )
                       Br(i,j,k) = fac1 * ( &
                            sf(1)*Br(i+flax,j,k) + sf(3)*Br(i,j,k+flaz) )
                       Bt(i,j,k) = fac2 * ( &
                            sf(1)*Bt(i+flax,j,k) + sf(3)*Bt(i,j,k+flaz) )
                       Bp(i,j,k) = fac2 * ( &
                            sf(1)*Bp(i+flax,j,k) + sf(3)*Bp(i,j,k+flaz) )
                    else if(extrap_method==LINEAR) then
                       rhob(i,j,k) = 1.d0/fr2(x(i,1,1),y(1,j,1),z(1,1,k)) * ( &
                            sf(1)* (2.d0*fr2(x(i+flax,1,1),y(1,j,1),z(1,1,k))*rhob(i+flax,j,k) &
                                       - fr2(x(i+2*flax,1,1),y(1,j,1),z(1,1,k))*rhob(i+2*flax,j,k)) &
                          + sf(3)* (2.d0*fr2(x(i,1,1),y(1,j,1),z(1,1,k+flaz))*rhob(i,j,k+flaz) &
                                       - fr2(x(i,1,1),y(1,j,1),z(1,1,k+2*flaz))*rhob(i,j,k+2*flaz)) )
                       P(i,j,k) = 1.d0/fr2(x(i,1,1),y(1,j,1),z(1,1,k)) * ( &
                            sf(1)* (2.d0*fr2(x(i+flax,1,1),y(1,j,1),z(1,1,k))*P(i+flax,j,k) &
                                       - fr2(x(i+2*flax,1,1),y(1,j,1),z(1,1,k))*P(i+2*flax,j,k)) &
                          + sf(3)* (2.d0*fr2(x(i,1,1),y(1,j,1),z(1,1,k+flaz))*P(i,j,k+flaz) &
                                       - fr2(x(i,1,1),y(1,j,1),z(1,1,k+2*flaz))*P(i,j,k+2*flaz)) )
                       vr(i,j,k) = 1.d0/fr1(x(i,1,1),y(1,j,1),z(1,1,k)) * ( &
                            sf(1)* (2.d0*fr1(x(i+flax,1,1),y(1,j,1),z(1,1,k))*vr(i+flax,j,k) &
                                       - fr1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k))*vr(i+2*flax,j,k)) &
                          + sf(3)* (2.d0*fr1(x(i,1,1),y(1,j,1),z(1,1,k+flaz))*vr(i,j,k+flaz) &
                                       - fr1(x(i,1,1),y(1,j,1),z(1,1,k+2*flaz))*vr(i,j,k+2*flaz)) )
                       vt(i,j,k) = 1.d0/frm1(x(i,1,1),y(1,j,1),z(1,1,k)) * ( &
                            sf(1)* (2.d0*frm1(x(i+flax,1,1),y(1,j,1),z(1,1,k))*vt(i+flax,j,k) &
                                       - frm1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k))*vt(i+2*flax,j,k)) &
                          + sf(3)* (2.d0*frm1(x(i,1,1),y(1,j,1),z(1,1,k+flaz))*vt(i,j,k+flaz) &
                                       - frm1(x(i,1,1),y(1,j,1),z(1,1,k+2*flaz))*vt(i,j,k+2*flaz)) )
                       vp(i,j,k) = 1.d0/frm1(x(i,1,1),y(1,j,1),z(1,1,k)) * ( &
                            sf(1)* (2.d0*frm1(x(i+flax,1,1),y(1,j,1),z(1,1,k))*vp(i+flax,j,k) &
                                       - frm1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k))*vp(i+2*flax,j,k)) &
                          + sf(3)* (2.d0*frm1(x(i,1,1),y(1,j,1),z(1,1,k+flaz))*vp(i,j,k+flaz) &
                                       - frm1(x(i,1,1),y(1,j,1),z(1,1,k+2*flaz))*vp(i,j,k+2*flaz)) )
                       Br(i,j,k) = 1.d0/fr2(x(i,1,1),y(1,j,1),z(1,1,k)) * ( &
                            sf(1)* (2.d0*fr2(x(i+flax,1,1),y(1,j,1),z(1,1,k))*Br(i+flax,j,k) &
                                       - fr2(x(i+2*flax,1,1),y(1,j,1),z(1,1,k))*Br(i+2*flax,j,k)) &
                          + sf(3)* (2.d0*fr2(x(i,1,1),y(1,j,1),z(1,1,k+flaz))*Br(i,j,k+flaz) &
                                       - fr2(x(i,1,1),y(1,j,1),z(1,1,k+2*flaz))*Br(i,j,k+2*flaz)) )
                       Bt(i,j,k) = 1.d0/frm1(x(i,1,1),y(1,j,1),z(1,1,k)) * ( &
                            sf(1)* (2.d0*frm1(x(i+flax,1,1),y(1,j,1),z(1,1,k))*Bt(i+flax,j,k) &
                                       - frm1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k))*Bt(i+2*flax,j,k)) &
                          + sf(3)* (2.d0*frm1(x(i,1,1),y(1,j,1),z(1,1,k+flaz))*Bt(i,j,k+flaz) &
                                       - frm1(x(i,1,1),y(1,j,1),z(1,1,k+2*flaz))*Bt(i,j,k+2*flaz)) )
                       Bp(i,j,k) = 1.d0/frm1(x(i,1,1),y(1,j,1),z(1,1,k)) * ( &
                            sf(1)* (2.d0*frm1(x(i+flax,1,1),y(1,j,1),z(1,1,k))*Bp(i+flax,j,k) &
                                       - frm1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k))*Bp(i+2*flax,j,k)) &
                          + sf(3)* (2.d0*frm1(x(i,1,1),y(1,j,1),z(1,1,k+flaz))*Bp(i,j,k+flaz) &
                                       - frm1(x(i,1,1),y(1,j,1),z(1,1,k+2*flaz))*Bp(i,j,k+2*flaz)) )
                    end if
                 else if( ni == 1 ) then
                    if(extrap_method==FLAT) then
                       r1 = x(i+flax,1,1)**2 + y(1,j,1)**2 + z(1,1,k+flaz)**2
                       r2 = x(i,1,1)**2 + y(1,j,1)**2 + z(1,1,k)**2
                       fac1 = r1/r2
                       fac2 = 1.d0 - 2.d0*(sqrt(r1)-sqrt(r2))/(sqrt(r1)+sqrt(r2))
                       fac3 = 1.d0 + 2.d0*(sqrt(r1)-sqrt(r2))/(sqrt(r1)+sqrt(r2))
                       rhob(i,j,k) = fac1 * rhob(i+flax,j,k+flaz)
                       P(i,j,k)    = fac1 * P(i+flax,j,k+flaz)
                       vr(i,j,k)   = fac3 * vr(i+flax,j,k+flaz)
                       vt(i,j,k)   = fac2 * vt(i+flax,j,k+flaz)
                       vp(i,j,k)   = fac2 * vp(i+flax,j,k+flaz)
                       Br(i,j,k)   = fac1 * Br(i+flax,j,k+flaz)
                       Bt(i,j,k)   = fac2 * Bt(i+flax,j,k+flaz)
                       Bp(i,j,k)   = fac2 * Bp(i+flax,j,k+flaz)
                    else if(extrap_method==LINEAR) then
                       rhob(i,j,k) = 1.d0/fr2(x(i,1,1),y(1,j,1),z(1,1,k)) * &
                            (2.d0*fr2(x(i+flax,1,1),y(1,j,1),z(1,1,k+flaz))*rhob(i+flax,j,k+flaz) &
                                - fr2(x(i+2*flax,1,1),y(1,j,1),z(1,1,k+2*flaz))*rhob(i+2*flax,j,k+2*flaz))
                       P(i,j,k) = 1.d0/fr2(x(i,1,1),y(1,j,1),z(1,1,k)) * &
                            (2.d0*fr2(x(i+flax,1,1),y(1,j,1),z(1,1,k+flaz))*P(i+flax,j,k+flaz) &
                                - fr2(x(i+2*flax,1,1),y(1,j,1),z(1,1,k+2*flaz))*P(i+2*flax,j,k+2*flaz))
                       vr(i,j,k) = 1.d0/fr1(x(i,1,1),y(1,j,1),z(1,1,k)) * &
                            (2.d0*fr1(x(i+flax,1,1),y(1,j,1),z(1,1,k+flaz))*vr(i+flax,j,k+flaz) &
                                - fr1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k+2*flaz))*vr(i+2*flax,j,k+2*flaz))
                       vt(i,j,k) = 1.d0/frm1(x(i,1,1),y(1,j,1),z(1,1,k)) * &
                            (2.d0*frm1(x(i+flax,1,1),y(1,j,1),z(1,1,k+flaz))*vt(i+flax,j,k+flaz) &
                                - frm1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k+2*flaz))*vt(i+2*flax,j,k+2*flaz))
                       vp(i,j,k) = 1.d0/frm1(x(i,1,1),y(1,j,1),z(1,1,k)) * &
                            (2.d0*frm1(x(i+flax,1,1),y(1,j,1),z(1,1,k+flaz))*vp(i+flax,j,k+flaz) &
                                - frm1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k+2*flaz))*vp(i+2*flax,j,k+2*flaz))
                       Br(i,j,k) = 1.d0/fr2(x(i,1,1),y(1,j,1),z(1,1,k)) * &
                            (2.d0*fr2(x(i+flax,1,1),y(1,j,1),z(1,1,k+flaz))*Br(i+flax,j,k+flaz) &
                                - fr2(x(i+2*flax,1,1),y(1,j,1),z(1,1,k+2*flaz))*Br(i+2*flax,j,k+2*flaz))
                       Bt(i,j,k) = 1.d0/frm1(x(i,1,1),y(1,j,1),z(1,1,k)) * &
                            (2.d0*frm1(x(i+flax,1,1),y(1,j,1),z(1,1,k+flaz))*Bt(i+flax,j,k+flaz) &
                                - frm1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k+2*flaz))*Bt(i+2*flax,j,k+2*flaz))
                       Bp(i,j,k) = 1.d0/frm1(x(i,1,1),y(1,j,1),z(1,1,k)) * &
                            (2.d0*frm1(x(i+flax,1,1),y(1,j,1),z(1,1,k+flaz))*Bp(i+flax,j,k+flaz) &
                                - frm1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k+2*flaz))*Bp(i+2*flax,j,k+2*flaz))
                    end if
              endif
!write(*,*) hi1a.,i,j,k,rhob(i,j,k)
              call compute_pcold_epscold(rhob(i,j,k), P_cold, eps_cold, &
                        neos, ergo_star, ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab)
              P(i,j,k) = P_cold
              eps = eps_cold
              r = sqrt(x(i,1,1)**2 + y(1,j,1)**2 + z(1,1,k)**2)
              costhet = z(1,1,k)/r
              sinthet = sqrt(1.d0 - costhet**2)
              cosphi  = x(i,1,1)/(r*sinthet)
              sinphi  = y(1,j,1)/(r*sinthet)
              vx(i,j,k) = sinthet*cosphi*vr(i,j,k) + &
                   r*costhet*cosphi*vt(i,j,k) - r*sinthet*sinphi*vp(i,j,k)
              vy(i,j,k) = sinthet*sinphi*vr(i,j,k) + &
                   r*costhet*sinphi*vt(i,j,k) + r*sinthet*cosphi*vp(i,j,k)
              vz(i,j,k) = costhet*vr(i,j,k) - r*sinthet*vt(i,j,k)
              Bx(i,j,k) = sinthet*cosphi*Br(i,j,k) + &
                   r*costhet*cosphi*Bt(i,j,k) - r*sinthet*sinphi*Bp(i,j,k)
              By(i,j,k) = sinthet*sinphi*Br(i,j,k) + &
                   r*costhet*sinphi*Bt(i,j,k) + r*sinthet*cosphi*Bp(i,j,k)
              Bz(i,j,k) = costhet*Br(i,j,k) - r*sinthet*Bt(i,j,k)
                 al = 1.d0 + alp(i,j,k)
                    ephi = exp(4.d0 * phi(i,j,k))
                    ! Compute al*u0-1
                    v2 = ephi*(gxx(i,j,k)*(vx(i,j,k) + betx(i,j,k))**2 + &
                         2.d0*gxy(i,j,k)*(vx(i,j,k) + betx(i,j,k)) &
                         *(vy(i,j,k) + bety(i,j,k)) +  &
                         2.d0*gxz(i,j,k)*(vx(i,j,k) + betx(i,j,k)) &
                         *(vz(i,j,k) + betz(i,j,k)) +  &
                         gyy(i,j,k)*(vy(i,j,k) + bety(i,j,k))**2 + &
                         2.d0*gyz(i,j,k)*(vy(i,j,k) + bety(i,j,k)) &
                         *(vz(i,j,k) + betz(i,j,k)) +  &
                         gzz(i,j,k)*(vz(i,j,k) + betz(i,j,k))**2 )/al**2
                    ! *** Check for superluminal velocity ***
                    if (v2 .gt. 1.d0) then
                       vx(i,j,k)=(vx(i,j,k)+betx(i,j,k))*sqrt(0.99d0/v2)-betx(i,j,k)
                       vy(i,j,k)=(vy(i,j,k)+bety(i,j,k))*sqrt(0.99d0/v2)-bety(i,j,k)
                       vz(i,j,k) =(vz(i,j,k)+betz(i,j,k))*sqrt(0.99d0/v2)-betz(i,j,k)
                       v2 = 0.99d0
                    end if
                    ! ***************************************
                    gam = sqrt(1.d0-v2)
                    au0 = 1.d0/gam
                    u0  = au0/al
                    sbt(i,j,k) = f1os4pi*ephi*u0/al*( &
                         gxx(i,j,k)*Bx(i,j,k)*(vx(i,j,k)+betx(i,j,k)) + &
                         gxy(i,j,k)*( Bx(i,j,k)*(vy(i,j,k)+bety(i,j,k)) &
                         + By(i,j,k)*(vx(i,j,k)+betx(i,j,k)) ) + &
                         gxz(i,j,k)*( Bx(i,j,k)*(vz(i,j,k)+betz(i,j,k)) &
                         + Bz(i,j,k)*(vx(i,j,k)+betx(i,j,k)) ) + &
                         gyy(i,j,k)*By(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gyz(i,j,k)*( By(i,j,k)*(vz(i,j,k)+betz(i,j,k)) &
                         + Bz(i,j,k)*(vy(i,j,k)+bety(i,j,k)) ) + &
                         gzz(i,j,k)*Bz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    sbx(i,j,k) = f1os4pi*Bx(i,j,k)/(al*u0) &
                         + sbt(i,j,k)*vx(i,j,k)
                    sby(i,j,k) = f1os4pi*By(i,j,k)/(al*u0) &
                         + sbt(i,j,k)*vy(i,j,k)
                    sbz(i,j,k) = f1os4pi*Bz(i,j,k)/(al*u0) &
                         + sbt(i,j,k)*vz(i,j,k)
                    sb2 = -(al*sbt(i,j,k))**2 &
                         + ephi*( gxx(i,j,k)*(&
                         sbx(i,j,k)+betx(i,j,k)*sbt(i,j,k))**2 + &
                         2.d0*gxy(i,j,k)*(sbx(i,j,k)+betx(i,j,k)*sbt(i,j,k)) &
                         *(sby(i,j,k)+bety(i,j,k)*sbt(i,j,k)) + &
                         2.d0*gxz(i,j,k)*(sbx(i,j,k)+betx(i,j,k)*sbt(i,j,k)) &
                         *(sbz(i,j,k)+betz(i,j,k)*sbt(i,j,k)) + &
                         gyy(i,j,k)*(sby(i,j,k)+bety(i,j,k)*sbt(i,j,k))**2 + &
                         2.d0*gyz(i,j,k)*(sby(i,j,k)+bety(i,j,k)*sbt(i,j,k)) &
                         *(sbz(i,j,k)+betz(i,j,k)*sbt(i,j,k)) + &
                         gzz(i,j,k)*(sbz(i,j,k)+betz(i,j,k)*sbt(i,j,k))**2 )
                    u_x = u0*ephi*( gxx(i,j,k)*(vx(i,j,k)+betx(i,j,k)) &
                         + gxy(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gxz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    u_y = u0*ephi*( gxy(i,j,k)*(vx(i,j,k)+betx(i,j,k)) &
                         + gyy(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gyz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    u_z = u0*ephi*( gxz(i,j,k)*(vx(i,j,k)+betx(i,j,k)) &
                         + gyz(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gzz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    sb_x = ephi*( gxx(i,j,k)*(sbx(i,j,k) &
                         +sbt(i,j,k)*betx(i,j,k)) &
                         + gxy(i,j,k)*(sby(i,j,k)+sbt(i,j,k)*bety(i,j,k)) + &
                         gxz(i,j,k)*(sbz(i,j,k)+sbt(i,j,k)*betz(i,j,k)) )
                    sb_y = ephi*( gxy(i,j,k)*(sbx(i,j,k) &
                         +sbt(i,j,k)*betx(i,j,k)) &
                         + gyy(i,j,k)*(sby(i,j,k)+sbt(i,j,k)*bety(i,j,k)) + &
                         gyz(i,j,k)*(sbz(i,j,k)+sbt(i,j,k)*betz(i,j,k)) )
                    sb_z = ephi*( gxz(i,j,k)*(sbx(i,j,k) &
                         +sbt(i,j,k)*betx(i,j,k)) &
                         + gyz(i,j,k)*(sby(i,j,k)+sbt(i,j,k)*bety(i,j,k)) + &
                         gzz(i,j,k)*(sbz(i,j,k)+sbt(i,j,k)*betz(i,j,k)) )
                    ephi = exp(6.d0 * phi(i,j,k))
                    rhos(i,j,k) = rhob(i,j,k)*ephi*au0
                    w(i,j,k) = au0*rhos(i,j,k)
                    enth(i,j,k) = 1.d0 + P(i,j,k)/rhob(i,j,k) + eps
                    tau(i,j,k) = w(i,j,k)*enth(i,j,k) - ephi*P(i,j,k) &
                         - rhos(i,j,k) &
                         + ephi*(sb2*au0**2 - 0.5d0*sb2 - (al*sbt(i,j,k))**2)
                    mhdstx(i,j,k) = al*ephi*((rhob(i,j,k)*enth(i,j,k)+sb2)*u0*u_x &
                         - sbt(i,j,k)*sb_x )
                    mhdsty(i,j,k) = al*ephi*((rhob(i,j,k)*enth(i,j,k)+sb2)*u0*u_y &
                         - sbt(i,j,k)*sb_y )
                    mhdstz(i,j,k) = al*ephi*((rhob(i,j,k)*enth(i,j,k)+sb2)*u0*u_z &
                         - sbt(i,j,k)*sb_z )
                    stx(i,j,k) = al*ephi*rhob(i,j,k)*enth(i,j,k)*u0*u_x
                    sty(i,j,k) = al*ephi*rhob(i,j,k)*enth(i,j,k)*u0*u_y
                    stz(i,j,k) = al*ephi*rhob(i,j,k)*enth(i,j,k)*u0*u_z
                    sbt(i,j,k) = sbt(i,j,k)/f1os4pi*al
                    sbx(i,j,k) = sbx(i,j,k)/f1os4pi*al
                    sby(i,j,k) = sby(i,j,k)/f1os4pi*al
                    sbz(i,j,k) = sbz(i,j,k)/f1os4pi*al
                 end if
           end do
        end do
     end do
  end if
  return
end subroutine hydro_ezbc_hybrid
!-------------------------------------------------------------------------+
!
! Use Alcubierres simple excision method for fluid, pure hydro case
!
! excision region : sphere
!
!-------------------------------------------------------------------------+
subroutine hydro_ezbc2(ex,x,y,z,rhos,tau,stx,sty,stz, &
     rhob,P,enth,vx,vy,vz,w,alp,betx,bety,betz,phi, &
     gxx,gxy,gxz,gyy,gyz,gzz,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
     symmetry,n,dm,kpoly_init)
  implicit none
!~~~~~~> Input arguments
  integer, intent(in) :: ex(1:3),symmetry
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   ::x,y,z
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::rhos,tau
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::stx,sty,stz
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::rhob,P,enth,vx,vy,vz,w
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   ::alp,betx,bety,betz,phi
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   ::gxx,gxy,gxz,gyy,gyz,gzz
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   ::gupxx,gupxy,gupxz
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   ::gupyy,gupyz,gupzz
  real*8    :: n,kpoly_init
  integer :: dm(1:ex(1),1:ex(2),1:ex(3))
!~~~~~~> local variable
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k,flax,flay,flaz,ni
  integer :: ksym,jsym,isym,diax,diay,diaz
  real*8  :: sf(1:3),rr,xs,zs,rxs,rzs,rxzs,c1,c2
  real*8  :: al,ephi,u0,au0,u_x,u_y,u_z,v2,gam,f1os4pi
  real*8  :: r1,r2,fac1,fac2,fac3,r,sinthet,costhet,sinphi,cosphi,Pmin
  integer :: extrap_method
  integer, parameter :: inEZ=0, onbd=1, outEZ=2
  real*8,  parameter :: zeo=0.d0, one=1.d0, two=2.d0
  integer, parameter :: nosymmetry=0, equatorial=1, octant=2
  integer, parameter :: pi_symm=3, axisym=4
  integer, parameter :: FLAT=0, LINEAR=1
  real*8,dimension(ex(1),ex(2),ex(3))   :: vr,vt,vp
!~~~~~~> interface
  interface
 !    subroutine find_excision_zone(ex,x,y,z,dm,innerbd,Symmetry)
 !      implicit none
 !      integer, intent(in)  :: ex(1:3)
 !      real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   ::x,y,z
 !      integer, intent(out) :: dm(1:ex(1),1:ex(2),1:ex(3))
 !      real*8               :: innerbd
 !      integer, intent(in)  :: Symmetry
 !    end subroutine find_excision_zone
  end interface
  real*8  :: fr2,fr1,frm1
  EXTERNAL fr2,fr1,frm1
!  extrap_method = FLAT
  extrap_method = LINEAR
  imin = lbound(rhos,1)
  jmin = lbound(rhos,2)
  kmin = lbound(rhos,3)
  imax = ubound(rhos,1)
  jmax = ubound(rhos,2)
  kmax = ubound(rhos,3)
 f1os4pi = 1.d0/sqrt(4.d0*acos(-1.d0))
!!$  c1 = 2.d0
  c1 = 1.d0
  c2 = c1 - 1.d0
!~~~~~~> Find the grid points on the boundary
!  call find_excision_zone(ex,x,y,z,dm,innerbd,Symmetry)
!~~~~~~> inner boundary condition
  if( symmetry == nosymmetry ) then
     isym = imin + 1
     jsym = jmin + 1
     ksym = kmin + 1
  else if( symmetry == equatorial ) then
     isym = imin + 1
     jsym = jmin + 1
     ksym = kmin
  else if( symmetry == octant ) then
     isym = imin
     jsym = jmin
     ksym = kmin
  else if( symmetry == pi_symm ) then
     isym = imin + 1
     jsym = jmin
     ksym = kmin
  endif
  do k = kmin, kmax, 1
     do j = jmin, jmax, 1
        do i = imin, imax, 1
           r = sqrt(x(i,1,1)**2 + y(1,j,1)**2 + z(1,1,k)**2)
           costhet = z(1,1,k)/r
           sinthet = sqrt(1.d0 - costhet**2)
           cosphi  = x(i,1,1)/(r*sinthet)
           sinphi  = y(1,j,1)/(r*sinthet)
           vr(i,j,k) = sinthet*cosphi*vx(i,j,k) + &
                sinthet*sinphi*vy(i,j,k) + costhet*vz(i,j,k)
           vt(i,j,k) = costhet*cosphi/r*vx(i,j,k) + &
                costhet*sinphi/r*vy(i,j,k) - sinthet/r*vz(i,j,k)
           vp(i,j,k) = -sinphi/sinthet/r*vx(i,j,k) + &
                cosphi/sinthet/r*vy(i,j,k)
        end do
     end do
  end do
  if(symmetry .ne. axisym) then
     do k = ksym, kmax-1, 1
        do j = jsym, jmax-1, 1
           do i = isym, imax-1, 1
              if( dm(i,j,k) == onbd ) then
                 ni = 2
                 flax = 1
                 flay = 1
                 flaz = 1
                 if(x(i,1,1) < zeo) flax = -1
                 if(y(1,j,1) < zeo) flay = -1
                 if(z(1,1,k) < zeo) flaz = -1
                 sf(1) = x(i,1,1) * flax
                 sf(2) = y(1,j,1) * flay
                 sf(3) = z(1,1,k) * flaz
                 if( dm(i+flax,j     ,k     ) == onbd ) sf(1) = zeo
                 if( dm(i     ,j+flay,k     ) == onbd ) sf(2) = zeo
                 if( dm(i     ,j     ,k+flaz) == onbd ) sf(3) = zeo
                 if((dm(i+flax,j,k) == outEZ) .or. (dm(i,j+flay,k) == outEZ)) &
                      sf(3) = zeo
                 if(( dm(i+flax,j     ,k     ) == onbd ) .and. &
                      ( dm(i     ,j+flay,k     ) == onbd ) .and. &
                      ( dm(i     ,j     ,k+flaz) == onbd )) then
                    ni = 1
                    if( dm(i     ,j+flay,k+flaz) == outEZ ) sf(2) = one
                    if( dm(i+flax,j     ,k+flaz) == outEZ ) sf(3) = one
                    if( dm(i+flax,j+flay,k     ) == outEZ ) then
                       sf(1) = one
                       sf(2) = zeo
                       sf(3) = zeo
                    endif
                    if(( dm(i+flax,j+flay,k) /= outEZ ) .and. &
                         ( dm(i,j+flay,k+flaz) /= outEZ ) .and. &
                         ( dm(i+flax,j,k+flaz) /= outEZ )) then
                         write(*,*) 'something wrong in find_excision_zone'
                    end if
                 endif
                 if( ni == 2 ) then
                    rr = sf(1) + sf(2) + sf(3)
                    sf(1) = sf(1) / rr
                    sf(2) = sf(2) / rr
                    sf(3) = sf(3) / rr
                    rhob(i,j,k) = sf(1) * rhob(i+flax,j,k) + &
                         sf(2) * rhob(i,j+flay,k) + &
                         sf(3) * rhob(i,j,k+flaz)
                    P(i,j,k) = sf(1) * P(i+flax,j,k) + &
                         sf(2) * P(i,j+flay,k) + &
                         sf(3) * P(i,j,k+flaz)
                    vx(i,j,k) = sf(1) * vx(i+flax,j,k) + &
                         sf(2) * vx(i,j+flay,k) + &
                         sf(3) * vx(i,j,k+flaz)
                    vy(i,j,k) = sf(1) * vy(i+flax,j,k) + &
                         sf(2) * vy(i,j+flay,k) + &
                         sf(3) * vy(i,j,k+flaz)
                    vz(i,j,k) = sf(1) * vz(i+flax,j,k) + &
                         sf(2) * vz(i,j+flay,k) + &
                         sf(3) * vz(i,j,k+flaz)
                 else if( ni == 1 ) then
                    rr = sf(1) + sf(2) +sf(3)
                    sf(1) = sf(1) / rr
                    sf(2) = sf(2) / rr
                    sf(3) = sf(3) / rr
                    rhob(i,j,k) = sf(1) * rhob(i+flax,j+flay,k) + &
                         sf(2) * rhob(i,j+flay,k+flaz) + &
                         sf(3) * rhob(i+flax,j,k+flaz)
                    P(i,j,k) = sf(1) * P(i+flax,j+flay,k) + &
                         sf(2) * P(i,j+flay,k+flaz) + &
                         sf(3) * P(i+flax,j,k+flaz)
                    vx(i,j,k) = sf(1) * vx(i+flax,j+flay,k) + &
                         sf(2) * vx(i,j+flay,k+flaz) + &
                         sf(3) * vx(i+flax,j,k+flaz)
                    vy(i,j,k) = sf(1) * vy(i+flax,j+flay,k) + &
                         sf(2) * vy(i,j+flay,k+flaz) + &
                         sf(3) * vy(i+flax,j,k+flaz)
                    vz(i,j,k) = sf(1) * vz(i+flax,j+flay,k) + &
                         sf(2) * vz(i,j+flay,k+flaz) + &
                         sf(3) * vz(i+flax,j,k+flaz)
                 endif
                    al = 1.d0 + alp(i,j,k)
                    ephi = exp(4.d0 * phi(i,j,k))
                    ! Compute al*u0-1
                    v2 = ephi*(gxx(i,j,k)*(vx(i,j,k) + betx(i,j,k))**2 + &
                         2.d0*gxy(i,j,k)*(vx(i,j,k) + betx(i,j,k)) &
                         *(vy(i,j,k) + bety(i,j,k)) +  &
                         2.d0*gxz(i,j,k)*(vx(i,j,k) + betx(i,j,k)) &
                         *(vz(i,j,k) + betz(i,j,k)) +  &
                         gyy(i,j,k)*(vy(i,j,k) + bety(i,j,k))**2 + &
                         2.d0*gyz(i,j,k)*(vy(i,j,k) + bety(i,j,k)) &
                         *(vz(i,j,k) + betz(i,j,k)) +  &
                         gzz(i,j,k)*(vz(i,j,k) + betz(i,j,k))**2 )/al**2
                    ! *** Check for superluminal velocity ***
                    if (v2 .gt. 1.d0) then
                       vx(i,j,k)=(vx(i,j,k)+betx(i,j,k))*sqrt(0.99d0/v2)-betx(i,j,k)
                       vy(i,j,k)=(vy(i,j,k)+bety(i,j,k))*sqrt(0.99d0/v2)-bety(i,j,k)
                       vz(i,j,k)=(vz(i,j,k)+betz(i,j,k))*sqrt(0.99d0/v2)-betz(i,j,k)
                       v2 = 0.99d0
                    end if
                    ! ***************************************
                    gam = sqrt(1.d0-v2)
                    au0 = 1.d0/gam
                    u0  = au0/al
                    u_x = u0*ephi*( gxx(i,j,k)*(vx(i,j,k)+betx(i,j,k)) &
                         + gxy(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gxz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    u_y = u0*ephi*( gxy(i,j,k)*(vx(i,j,k)+betx(i,j,k)) &
                         + gyy(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gyz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    u_z = u0*ephi*( gxz(i,j,k)*(vx(i,j,k)+betx(i,j,k)) &
                         + gyz(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gzz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    ephi = exp(6.d0 * phi(i,j,k))
                    rhos(i,j,k) = rhob(i,j,k)*ephi*au0
                    w(i,j,k) = au0*rhos(i,j,k)
                    enth(i,j,k) = 1.d0 + (1.d0+n)*P(i,j,k)/rhob(i,j,k)
                    tau(i,j,k) = w(i,j,k)*enth(i,j,k) - ephi*P(i,j,k) &
                         - rhos(i,j,k)
                    stx(i,j,k) = al*ephi*rhob(i,j,k)*enth(i,j,k)*u0*u_x
                    sty(i,j,k) = al*ephi*rhob(i,j,k)*enth(i,j,k)*u0*u_y
                    stz(i,j,k) = al*ephi*rhob(i,j,k)*enth(i,j,k)*u0*u_z
                 endif
              enddo
           enddo
        enddo
     else
!!$        j = 2
        do j = 1, 3, 1
        do k = kmin, kmax-1, 1
              do i = imin, imax-1, 1
              if( dm(i,j,k) == onbd ) then
                 ni = 2
                 flax = 1
                 flaz = 1
                 if(x(i,1,1) < zeo) flax = -1
                 if(z(1,1,k) < zeo) flaz = -1
                 sf(1) = x(i,1,1) * flax
                 sf(3) = z(1,1,k) * flaz
                 if( dm(i+flax,j     ,k     ) == onbd ) sf(1) = zeo
                 if( dm(i     ,j     ,k+flaz) == onbd ) sf(3) = zeo
!!$                 if(dm(i+flax,j,k) == outEZ) sf(3) = zeo
                 if((dm(i+flax,j,k) == onbd).and.(dm(i,j,k+flaz) == onbd)) then
                    ni = 1
                    if( dm(i+flax,j,k+flaz) /= outEZ ) then
                         write(*,*) 'something wrong in find_excision_zone'
                      end if
                 endif
                 if( ni == 2 ) then
                    rr = sf(1) + sf(3)
                    sf(1) = sf(1) / rr
                    sf(3) = sf(3) / rr
                    if(extrap_method==FLAT) then
                       r1 =   (sf(1)*x(i+flax,1,1)+(1.d0-sf(1))*x(i,1,1))**2 &
                            + y(1,j,1)**2 &
                            + (sf(3)*z(1,1,k+flaz)+(1.d0-sf(3))*z(1,1,k))**2
                       r2 = x(i,1,1)**2 + y(1,j,1)**2 + z(1,1,k)**2
                       fac1 = r1/r2
                       fac2 = 1.d0 - &
                            2.d0*(sqrt(r1)-sqrt(r2))/(sqrt(r1)+sqrt(r2))
                       fac3 = 1.d0 + &
                            2.d0*(sqrt(r1)-sqrt(r2))/(sqrt(r1)+sqrt(r2))
                       rhob(i,j,k) = fac1 * ( &
                            sf(1)*rhob(i+flax,j,k) + sf(3)*rhob(i,j,k+flaz) )
                       P(i,j,k) = fac1 * ( &
                            sf(1)*P(i+flax,j,k) + sf(3)*P(i,j,k+flaz) )
                       vr(i,j,k) = fac3 * ( &
                            sf(1)*vr(i+flax,j,k) + sf(3)*vr(i,j,k+flaz) )
                       vt(i,j,k) = fac2 * ( &
                            sf(1)*vt(i+flax,j,k) + sf(3)*vt(i,j,k+flaz) )
                       vp(i,j,k) = fac2 * ( &
                            sf(1)*vp(i+flax,j,k) + sf(3)*vp(i,j,k+flaz) )
                    else if(extrap_method==LINEAR) then
                       rhob(i,j,k) = 1.d0/fr2(x(i,1,1),y(1,j,1),z(1,1,k)) * ( &
                            sf(1)* (2.d0*fr2(x(i+flax,1,1),y(1,j,1),z(1,1,k))*rhob(i+flax,j,k) &
                                       - fr2(x(i+2*flax,1,1),y(1,j,1),z(1,1,k))*rhob(i+2*flax,j,k)) &
                          + sf(3)* (2.d0*fr2(x(i,1,1),y(1,j,1),z(1,1,k+flaz))*rhob(i,j,k+flaz) &
                                       - fr2(x(i,1,1),y(1,j,1),z(1,1,k+2*flaz))*rhob(i,j,k+2*flaz)) )
                       P(i,j,k) = 1.d0/fr2(x(i,1,1),y(1,j,1),z(1,1,k)) * ( &
                            sf(1)* (2.d0*fr2(x(i+flax,1,1),y(1,j,1),z(1,1,k))*P(i+flax,j,k) &
                                       - fr2(x(i+2*flax,1,1),y(1,j,1),z(1,1,k))*P(i+2*flax,j,k)) &
                          + sf(3)* (2.d0*fr2(x(i,1,1),y(1,j,1),z(1,1,k+flaz))*P(i,j,k+flaz) &
                                       - fr2(x(i,1,1),y(1,j,1),z(1,1,k+2*flaz))*P(i,j,k+2*flaz)) )
                       vr(i,j,k) = 1.d0/fr1(x(i,1,1),y(1,j,1),z(1,1,k)) * ( &
                            sf(1)* (2.d0*fr1(x(i+flax,1,1),y(1,j,1),z(1,1,k))*vr(i+flax,j,k) &
                                       - fr1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k))*vr(i+2*flax,j,k)) &
                          + sf(3)* (2.d0*fr1(x(i,1,1),y(1,j,1),z(1,1,k+flaz))*vr(i,j,k+flaz) &
                                       - fr1(x(i,1,1),y(1,j,1),z(1,1,k+2*flaz))*vr(i,j,k+2*flaz)) )
                       vt(i,j,k) = 1.d0/frm1(x(i,1,1),y(1,j,1),z(1,1,k)) * ( &
                            sf(1)* (2.d0*frm1(x(i+flax,1,1),y(1,j,1),z(1,1,k))*vt(i+flax,j,k) &
                                       - frm1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k))*vt(i+2*flax,j,k)) &
                          + sf(3)* (2.d0*frm1(x(i,1,1),y(1,j,1),z(1,1,k+flaz))*vt(i,j,k+flaz) &
                                       - frm1(x(i,1,1),y(1,j,1),z(1,1,k+2*flaz))*vt(i,j,k+2*flaz)) )
                       vp(i,j,k) = 1.d0/frm1(x(i,1,1),y(1,j,1),z(1,1,k)) * ( &
                            sf(1)* (2.d0*frm1(x(i+flax,1,1),y(1,j,1),z(1,1,k))*vp(i+flax,j,k) &
                                       - frm1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k))*vp(i+2*flax,j,k)) &
                          + sf(3)* (2.d0*frm1(x(i,1,1),y(1,j,1),z(1,1,k+flaz))*vp(i,j,k+flaz) &
                                       - frm1(x(i,1,1),y(1,j,1),z(1,1,k+2*flaz))*vp(i,j,k+2*flaz)) )
                    end if
                 else if( ni == 1 ) then
                    if(extrap_method==FLAT) then
                       r1 = x(i+flax,1,1)**2 + y(1,j,1)**2 + z(1,1,k+flaz)**2
                       r2 = x(i,1,1)**2 + y(1,j,1)**2 + z(1,1,k)**2
                       fac1 = r1/r2
                       fac2 = 1.d0 - 2.d0*(sqrt(r1)-sqrt(r2))/(sqrt(r1)+sqrt(r2))
                       fac3 = 1.d0 + 2.d0*(sqrt(r1)-sqrt(r2))/(sqrt(r1)+sqrt(r2))
                       rhob(i,j,k) = fac1 * rhob(i+flax,j,k+flaz)
                       P(i,j,k)    = fac1 * P(i+flax,j,k+flaz)
                       vr(i,j,k)   = fac3 * vr(i+flax,j,k+flaz)
                       vt(i,j,k)   = fac2 * vt(i+flax,j,k+flaz)
                       vp(i,j,k)   = fac2 * vp(i+flax,j,k+flaz)
                    else if(extrap_method==LINEAR) then
                       rhob(i,j,k) = 1.d0/fr2(x(i,1,1),y(1,j,1),z(1,1,k)) * &
                            (2.d0*fr2(x(i+flax,1,1),y(1,j,1),z(1,1,k+flaz))*rhob(i+flax,j,k+flaz) &
                                - fr2(x(i+2*flax,1,1),y(1,j,1),z(1,1,k+2*flaz))*rhob(i+2*flax,j,k+2*flaz))
                       P(i,j,k) = 1.d0/fr2(x(i,1,1),y(1,j,1),z(1,1,k)) * &
                            (2.d0*fr2(x(i+flax,1,1),y(1,j,1),z(1,1,k+flaz))*P(i+flax,j,k+flaz) &
                                - fr2(x(i+2*flax,1,1),y(1,j,1),z(1,1,k+2*flaz))*P(i+2*flax,j,k+2*flaz))
                       vr(i,j,k) = 1.d0/fr1(x(i,1,1),y(1,j,1),z(1,1,k)) * &
                            (2.d0*fr1(x(i+flax,1,1),y(1,j,1),z(1,1,k+flaz))*vr(i+flax,j,k+flaz) &
                                - fr1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k+2*flaz))*vr(i+2*flax,j,k+2*flaz))
                       vt(i,j,k) = 1.d0/frm1(x(i,1,1),y(1,j,1),z(1,1,k)) * &
                            (2.d0*frm1(x(i+flax,1,1),y(1,j,1),z(1,1,k+flaz))*vt(i+flax,j,k+flaz) &
                                - frm1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k+2*flaz))*vt(i+2*flax,j,k+2*flaz))
                       vp(i,j,k) = 1.d0/frm1(x(i,1,1),y(1,j,1),z(1,1,k)) * &
                            (2.d0*frm1(x(i+flax,1,1),y(1,j,1),z(1,1,k+flaz))*vp(i+flax,j,k+flaz) &
                                - frm1(x(i+2*flax,1,1),y(1,j,1),z(1,1,k+2*flaz))*vp(i+2*flax,j,k+2*flaz))
                    end if
              endif
!!$              P(i,j,k) = kpoly_init*rhob(i,j,k)**(1.d0 + 1.d0/n)
              r = sqrt(x(i,1,1)**2 + y(1,j,1)**2 + z(1,1,k)**2)
              costhet = z(1,1,k)/r
              sinthet = sqrt(1.d0 - costhet**2)
              cosphi  = x(i,1,1)/(r*sinthet)
              sinphi  = y(1,j,1)/(r*sinthet)
              vx(i,j,k) = sinthet*cosphi*vr(i,j,k) + &
                   r*costhet*cosphi*vt(i,j,k) - r*sinthet*sinphi*vp(i,j,k)
              vy(i,j,k) = sinthet*sinphi*vr(i,j,k) + &
                   r*costhet*sinphi*vt(i,j,k) + r*sinthet*cosphi*vp(i,j,k)
              vz(i,j,k) = costhet*vr(i,j,k) - r*sinthet*vt(i,j,k)
                 al = 1.d0 + alp(i,j,k)
                    ! Compute al*u0-1
                    v2 = ephi*(gxx(i,j,k)*(vx(i,j,k) + betx(i,j,k))**2 + &
                         2.d0*gxy(i,j,k)*(vx(i,j,k) + betx(i,j,k)) &
                         *(vy(i,j,k) + bety(i,j,k)) +  &
                         2.d0*gxz(i,j,k)*(vx(i,j,k) + betx(i,j,k)) &
                         *(vz(i,j,k) + betz(i,j,k)) +  &
                         gyy(i,j,k)*(vy(i,j,k) + bety(i,j,k))**2 + &
                         2.d0*gyz(i,j,k)*(vy(i,j,k) + bety(i,j,k)) &
                         *(vz(i,j,k) + betz(i,j,k)) +  &
                         gzz(i,j,k)*(vz(i,j,k) + betz(i,j,k))**2 )/al**2
                    ! *** Check for superluminal velocity ***
                    if (v2 .gt. 1.d0) then
                       vx(i,j,k)=(vx(i,j,k)+betx(i,j,k))*sqrt(0.99d0/v2)-betx(i,j,k)
                       vy(i,j,k)=(vy(i,j,k)+bety(i,j,k))*sqrt(0.99d0/v2)-bety(i,j,k)
                       vz(i,j,k) =(vz(i,j,k)+betz(i,j,k))*sqrt(0.99d0/v2)-betz(i,j,k)
                       v2 = 0.99d0
                    end if
                    ! ***************************************
                    gam = sqrt(1.d0-v2)
                    au0 = 1.d0/gam
                    u0  = au0/al
                    u_x = u0*ephi*( gxx(i,j,k)*(vx(i,j,k)+betx(i,j,k)) &
                         + gxy(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gxz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    u_y = u0*ephi*( gxy(i,j,k)*(vx(i,j,k)+betx(i,j,k)) &
                         + gyy(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gyz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    u_z = u0*ephi*( gxz(i,j,k)*(vx(i,j,k)+betx(i,j,k)) &
                         + gyz(i,j,k)*(vy(i,j,k)+bety(i,j,k)) + &
                         gzz(i,j,k)*(vz(i,j,k)+betz(i,j,k)) )
                    ephi = exp(6.d0 * phi(i,j,k))
                    rhos(i,j,k) = rhob(i,j,k)*ephi*au0
                    w(i,j,k) = au0*rhos(i,j,k)
                    enth(i,j,k) = 1.d0 + (1.d0+n)*P(i,j,k)/rhob(i,j,k)
                    tau(i,j,k) = w(i,j,k)*enth(i,j,k) - ephi*P(i,j,k) &
                         - rhos(i,j,k)
                    stx(i,j,k) = al*ephi*rhob(i,j,k)*enth(i,j,k)*u0*u_x
                    sty(i,j,k) = al*ephi*rhob(i,j,k)*enth(i,j,k)*u0*u_y
                    stz(i,j,k) = al*ephi*rhob(i,j,k)*enth(i,j,k)*u0*u_z
                 end if
           end do
        end do
        end do
  end if
  return
end subroutine hydro_ezbc2
