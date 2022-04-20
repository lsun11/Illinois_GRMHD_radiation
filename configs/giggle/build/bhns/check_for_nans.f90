subroutine check_for_nans(ext,rho_b,P,vx,vy,vz, phi,lapm1,shiftx,shifty,shiftz, gxx,gxy,gxz,gyy,gyz,gzz,gupxx,gupxy,gupxz,gupyy,gup&
  &yz,gupzz, &
          Bx,By,Bz,rho_star,tau,st_x,st_y,st_z)
  implicit none
  integer, dimension(3)                    :: ext
  INTEGER :: i,j,k
  real*8, dimension(ext(1),ext(2),ext(3))  :: rho_b,P,vx,vy,vz, phi,lapm1,shiftx,shifty,shiftz, gxx,gxy,gxz,gyy,gyz,gzz,gupxx,gupxy&
  &,gupxz,gupyy,gupyz,gupzz,&
        Bx,By,Bz,rho_star,tau,st_x,st_y,st_z
  do k=1,ext(3)
     do j=1,ext(2)
        do i=1,ext(1)
           if(isnan(rho_b(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "rho_b= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(P(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "P= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(vx(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "vx= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(vy(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "vy= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(vz(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "vz= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(phi(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "phi= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(lapm1(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "lapm1= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(shiftx(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "shiftx= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(shifty(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "shifty= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(shiftz(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "shiftz= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(gxx(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "gxx= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(gxy(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "gxy= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(gxz(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "gxz= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(gyy(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "gyy= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(gyz(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "gyz= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(gzz(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "gzz= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(gupxx(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "gxx= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(gupxy(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "gxy= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(gupxz(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "gxz= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(gupyy(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "gyy= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(gupyz(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "gyz= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(gupzz(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "gzz= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(Bx(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "Bx= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(By(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "By= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(Bz(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "Bz= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(rho_star(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "rho_star= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(tau(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "tau= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(st_x(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "st_x= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(st_y(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "st_y= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
           if(isnan(st_z(i,j,k))) then
              write(*,*) "initialdata_local: FOUND NAN AT ", i,ext(1),j,ext(2),k,ext(3)
              write(*,*)  "st_z= ",rho_b(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),phi(i,j,k),Bx(i,j,k),by(i,j,k),bz(i,j,k)
              stop
           end if
        end do
     end do
  end do
end subroutine check_for_nans
