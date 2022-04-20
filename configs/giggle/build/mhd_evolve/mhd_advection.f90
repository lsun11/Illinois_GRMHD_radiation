!-----------------------------------------------------------!
!-----------------------------------------------------------!
subroutine find_ur_ul_emfields(ext,X,Y,Z,Bx,By,Bz, Bxr,Bxl, Byr,Byl, &
                Bzr,Bzl, dBx,dBy,dBz, m,Symmetry,Sym_Bz,Reconstruction)
 implicit none
  interface
     subroutine find_ur_ul(ext,X,Y,Z,U,Ur,Ul,nabla_U,m,Sym_arg,Reconstruction,Symmetry)
       implicit none
       integer, dimension(3)                      :: ext
       real*8, dimension(ext(1),ext(2),ext(3))    :: X,Y,Z
       integer                                    :: m,Reconstruction,Symmetry
       real*8                                     :: Sym_arg
       real*8, dimension(ext(1),ext(2),ext(3))    :: U, Ur, Ul, nabla_u
     end subroutine find_ur_ul
  end interface
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bx,By,Bz
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bxr,Bxl, Byr,Byl, Bzr,Bzl
 real*8, dimension(ext(1),ext(2),ext(3))        :: dBx,dBy,dBz
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8                                                :: Sym_Bz,SYM,ANTI
 integer                                        :: Symmetry, m, kmin
 integer                                        :: Reconstruction
 integer                                        :: NO_SYMM,EQUATORIAL,OCTANT
 integer                                        :: PI_SYMM,AXISYM
 parameter (NO_SYMM=0, EQUATORIAL=1, OCTANT=2, PI_SYMM=3, AXISYM=4)
 parameter (SYM = 1.d0, ANTI = -1.d0)
!
 if (m==1) then
    if (Symmetry==OCTANT) then
       call find_ur_ul(ext,X,Y,Z,Bx,Bxr,Bxl,dBx,m,ANTI,Reconstruction,Symmetry)
       call find_ur_ul(ext,X,Y,Z,By,Byr,Byl,dBy,m,SYM,Reconstruction,Symmetry)
    else if (Symmetry==AXISYM) then
       call find_ur_ul(ext,X,Y,Z,Bx,Bxr,Bxl,dBx,m,ANTI,Reconstruction,Symmetry)
       call find_ur_ul(ext,X,Y,Z,By,Byr,Byl,dBy,m,ANTI,Reconstruction,Symmetry)
    else
       call find_ur_ul(ext,X,Y,Z,Bx,Bxr,Bxl,dBx,m,SYM,Reconstruction,Symmetry)
       call find_ur_ul(ext,X,Y,Z,By,Byr,Byl,dBy,m,SYM,Reconstruction,Symmetry)
    end if
    call find_ur_ul(ext,X,Y,Z,Bz,Bzr,Bzl,dBz,m,SYM,Reconstruction,Symmetry)
 else if (m==2) then
    if (Symmetry==OCTANT) then
       call find_ur_ul(ext,X,Y,Z,By,Byr,Byl,dBy,m,ANTI,Reconstruction,Symmetry)
    else
       call find_ur_ul(ext,X,Y,Z,By,Byr,Byl,dBy,m,SYM,Reconstruction,Symmetry)
    end if
    call find_ur_ul(ext,X,Y,Z,Bx,Bxr,Bxl,dBx,m,SYM,Reconstruction,Symmetry)
    call find_ur_ul(ext,X,Y,Z,Bz,Bzr,Bzl,dBz,m,SYM,Reconstruction,Symmetry)
 else if (m==3) then
    if(Symmetry == AXISYM .or. Symmetry==OCTANT .or. Symmetry==EQUATORIAL) then
       kmin = lbound(Bx,3)+1
    end if
    if (Symmetry==NO_SYMM .or. Z(1,1,kmin) .lt. 0.d0) then
       call find_ur_ul(ext,X,Y,Z,Bz,Bzr,Bzl,dBz,m,SYM,Reconstruction,Symmetry)
       call find_ur_ul(ext,X,Y,Z,Bx,Bxr,Bxl,dBx,m,SYM,Reconstruction,Symmetry)
       call find_ur_ul(ext,X,Y,Z,By,Byr,Byl,dBy,m,SYM,Reconstruction,Symmetry)
    else
       call find_ur_ul(ext,X,Y,Z,Bz,Bzr,Bzl,dBz,m,Sym_Bz,Reconstruction,Symmetry)
       call find_ur_ul(ext,X,Y,Z,Bx,Bxr,Bxl,dBx,m,ANTI*Sym_Bz,Reconstruction,Symmetry)
       call find_ur_ul(ext,X,Y,Z,By,Byr,Byl,dBy,m,ANTI*Sym_Bz,Reconstruction,Symmetry)
    end if
 end if
end subroutine find_ur_ul_emfields
!--------------------------------------------------!
!
!--------------------------------------------------!
!
subroutine compute_nabla_emfields(ext,dX,dY,dZ,Bx,By,Bz,dBx,dBy,dBz, &
                        m,Symmetry,Sym_Bz,X,Y,Z,Reconstruction)
 implicit none
 interface
    subroutine compute_nabla(ext,dX,dY,dZ,U,nabla_U,m,Sym_arg,X,Y,Z,LIMITER,Symmetry)
      implicit none
      integer, dimension(3)                      :: ext
      integer                                    :: m,LIMITER,Symmetry
      real*8                                     :: Sym_arg
      real*8, dimension(ext(1),ext(2),ext(3))    :: U, nabla_U
      real*8, dimension(ext(1),ext(2),ext(3))         :: X,Y,Z
      real*8                                     :: dX, dY, dZ
    end subroutine compute_nabla
 end interface
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bx,By,Bz,dBx,dBy,dBz
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 integer                                        :: m, Symmetry, Reconstruction
 real*8                                                :: dX,dY,dZ, Sym_Bz, SYM,ANTI
 integer                                        :: NO_SYMM,EQUATORIAL,OCTANT
 integer                                        :: PI_SYMM,AXISYM
  integer                                            :: SPPM, LIMITER,dU,MC
 parameter (NO_SYMM=0, EQUATORIAL=1, OCTANT=2, PI_SYMM=3, AXISYM=4)
 parameter (SPPM=5, MC=1, dU=4)
 parameter (SYM = 1.d0, ANTI = -1.d0)
! 
 if (Reconstruction == SPPM) then
    LIMITER = dU
 else
    LIMITER = MC
 end if
 if (m==1) then
    if (Symmetry==OCTANT) then
       call compute_nabla(ext,dX,dY,dZ,Bx,dBx,m,ANTI,X,Y,Z,LIMITER,Symmetry)
       call compute_nabla(ext,dX,dY,dZ,By,dBy,m,SYM,X,Y,Z,LIMITER,Symmetry)
    else if (Symmetry==AXISYM) then
       call compute_nabla(ext,dX,dY,dZ,Bx,dBx,m,ANTI,X,Y,Z,LIMITER,Symmetry)
       call compute_nabla(ext,dX,dY,dZ,By,dBy,m,ANTI,X,Y,Z,LIMITER,Symmetry)
    else
       call compute_nabla(ext,dX,dY,dZ,Bx,dBx,m,SYM,X,Y,Z,LIMITER,Symmetry)
       call compute_nabla(ext,dX,dY,dZ,By,dBy,m,SYM,X,Y,Z,LIMITER,Symmetry)
    end if
    call compute_nabla(ext,dX,dY,dZ,Bz,dBz,m,SYM,X,Y,Z,LIMITER,Symmetry)
 else if (m==2) then
    if (Symmetry==OCTANT) then
       call compute_nabla(ext,dX,dY,dZ,By,dBy,m,ANTI,X,Y,Z,LIMITER,Symmetry)
    else
       call compute_nabla(ext,dX,dY,dZ,By,dBy,m,SYM,X,Y,Z,LIMITER,Symmetry)
    end if
    call compute_nabla(ext,dX,dY,dZ,Bx,dBx,m,SYM,X,Y,Z,LIMITER,Symmetry)
    call compute_nabla(ext,dX,dY,dZ,Bz,dBz,m,SYM,X,Y,Z,LIMITER,Symmetry)
 else if (m==3) then
    if (Symmetry==NO_SYMM) then
       call compute_nabla(ext,dX,dY,dZ,Bz,dBz,m,SYM,X,Y,Z,LIMITER,Symmetry)
       call compute_nabla(ext,dX,dY,dZ,Bx,dBx,m,SYM,X,Y,Z,LIMITER,Symmetry)
       call compute_nabla(ext,dX,dY,dZ,By,dBy,m,SYM,X,Y,Z,LIMITER,Symmetry)
    else
       call compute_nabla(ext,dX,dY,dZ,Bz,dBz,m,Sym_Bz,X,Y,Z,LIMITER,Symmetry)
       call compute_nabla(ext,dX,dY,dZ,Bx,dBx,m,ANTI*Sym_Bz,X,Y,Z,LIMITER,Symmetry)
       call compute_nabla(ext,dX,dY,dZ,By,dBy,m,ANTI*Sym_Bz,X,Y,Z,LIMITER,Symmetry)
    end if
 end if
! write(*,*) NABLA_MHD: ,Bz(11,2,2),Bz(10,2,2),Bz(15,2,2),dBz(11,2,2)
end subroutine compute_nabla_emfields
!
!-----------------------------------------------------------------------------
! Advect Z_i
!-----------------------------------------------------------------------------
!
subroutine advect_z(ext, X, Z_y, vy, P,alpha,phi, betax,betay,betaz, &
     gxx,gxy,gxz,gyy,gyz,gzz, sbt,sbx,sby,sbz, st_x_rhs, st_y_rhs, st_z_rhs, &
     fmx, fmy, fmz, m, dX, dY, dZ,Symmetry)
  implicit none
  integer, dimension(3)                    :: ext
  real*8, dimension(ext(1),ext(2),ext(3))  :: X
  real*8, dimension(ext(1),ext(2),ext(3))  :: st_x_rhs, st_y_rhs, st_z_rhs
  real*8, dimension(ext(1),ext(2),ext(3))  :: P,alpha,phi
  real*8, dimension(ext(1),ext(2),ext(3))  :: fmx, fmy, fmz
  real*8, dimension(ext(1),ext(2),ext(3))  :: Z_y, vy
  real*8, dimension(ext(1),ext(2),ext(3))  :: betax,betay,betaz
  real*8, dimension(ext(1),ext(2),ext(3))  :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ext(1),ext(2),ext(3))  :: sbt,sbx,sby,sbz
  real*8, dimension(ext(1),ext(2),ext(3))  :: sb_y,b2
  real*8  :: dX, dY, dZ
  integer                                   :: Symmetry
!
! Internal variables
!
  real*8  :: PI, f1o4p, f1o8p, ONE  = 1.D0, TWO = 2.D0, HALF = 0.5D0
  real*8  :: FOUR = 4.D0, SIX = 6.D0, SYM  = 1.D0, ANTI = -1.D0
  integer              :: i,j,k,mmax,mmin,m,imin,imax,jmin,jmax,kmin,kmax
  real*8               :: Gamma
  integer, parameter   :: AXISYM = 4
 mmax = ubound(fmx,m)
 mmin = lbound(fmx,m)
! write(*,*) m1, st_z_rhs(30,30,2),fmz(30,30,2),fmz(31,30,2),fmz(30,31,2),fmz(30,30,3)
 if (m==1) then
    if (Symmetry==AXISYM) then
       f1o4p = 0.25d0/(acos(-1.d0))
       sb_y = exp(4.d0*phi)*( gxy*(sbt*betax + sbx) + gyy*(sbt*betay + sby) + &
            gyz*(sbt*betaz + sbz) )
       b2 = - ((1.d0+alpha)*sbt)**2 + exp(4.d0*phi)*( gxx*(sbx + betax*sbt)**2 + &
            2.d0*gxy*(sbx + betax*sbt)*(sby + betay*sbt) + &
            2.d0*gxz*(sbx + betax*sbt)*(sbz + betaz*sbt) + &
            gyy*(sby + betay*sbt)**2 + &
            2.d0*gyz*(sby + betay*sbt)*(sbz + betaz*sbt) + &
            gzz*(sbz + betaz*sbt)**2 )
       do i=mmin,mmax-1
          st_x_rhs(i,:,:) = st_x_rhs(i,:,:) + &
               ( (fmx(i,:,:)-fmx(i+1,:,:))/dX + Z_y(i,:,:)*vy(i,:,:) + &
               (1.d0+alpha(i,:,:))*exp(6.d0*phi(i,:,:))*(P(i,:,:)+ &
               f1o4p/(1.d0+alpha(i,:,:))**2 *(HALF*b2(i,:,:) + sb_y(i,:,:) * &
               (sbt(i,:,:)*vy(i,:,:) - sby(i,:,:)) ) ))/X(i,1,1)
          st_y_rhs(i,:,:) = st_y_rhs(i,:,:) + &
               (fmy(i,:,:)-fmy(i+1,:,:))/dX/X(i,1,1)**2
          st_z_rhs(i,:,:) = st_z_rhs(i,:,:) + (fmz(i,:,:)-fmz(i+1,:,:))/dX/X(i,1,1)
       end do
!       write(*,*) hiblah,st_x_rhs(5,2,2),(fmx(5,2,2)-fmx(6,2,2))/dX,Z_y(5,2,2)*vy(5,2,2),P(5,2,2),f1o4p/(alpha(5,2,2)+1.D0)**2 *(0.
    else
       do i=mmin,mmax-1
          st_x_rhs(i,:,:) = st_x_rhs(i,:,:) + (fmx(i,:,:)-fmx(i+1,:,:))/dX
          st_y_rhs(i,:,:) = st_y_rhs(i,:,:) + (fmy(i,:,:)-fmy(i+1,:,:))/dX
          st_z_rhs(i,:,:) = st_z_rhs(i,:,:) + (fmz(i,:,:)-fmz(i+1,:,:))/dX
       end do
    end if
    st_x_rhs(mmax,:,:) = 0.d0
    st_y_rhs(mmax,:,:) = 0.d0
    st_z_rhs(mmax,:,:) = 0.d0
!    write(*,*) after advect_z: ,m,st_x_rhs(2,2,2),st_y_rhs(2,2,2),st_z_rhs(2,2,2),fmx(2,2,2),fmy(2,2,2),fmz(2,2,2),fmx(17,2,15),fmy
    return
 elseif (m==2) then
  if (Symmetry .ne. AXISYM) then
    do j=mmin,mmax-1
       st_x_rhs(:,j,:) = st_x_rhs(:,j,:) + (fmx(:,j,:)-fmx(:,j+1,:))/dY
       st_y_rhs(:,j,:) = st_y_rhs(:,j,:) + (fmy(:,j,:)-fmy(:,j+1,:))/dY
       st_z_rhs(:,j,:) = st_z_rhs(:,j,:) + (fmz(:,j,:)-fmz(:,j+1,:))/dY
    end do
  end if
    st_x_rhs(:,mmax,:) = 0.d0
    st_y_rhs(:,mmax,:) = 0.d0
    st_z_rhs(:,mmax,:) = 0.d0
!    write(*,*) after advect_z: ,m,st_x_rhs(2,2,2),st_y_rhs(2,2,2),st_z_rhs(2,2,2),fmx(2,2,2),fmy(2,2,2),fmz(2,2,2),fmx(16,3,15),fmy
    return
 elseif (m==3) then
    do k=mmin,mmax-1
       st_x_rhs(:,:,k) = st_x_rhs(:,:,k) + (fmx(:,:,k)-fmx(:,:,k+1))/dZ
       st_y_rhs(:,:,k) = st_y_rhs(:,:,k) + (fmy(:,:,k)-fmy(:,:,k+1))/dZ
       st_z_rhs(:,:,k) = st_z_rhs(:,:,k) + (fmz(:,:,k)-fmz(:,:,k+1))/dZ
    end do
    st_x_rhs(:,:,mmax) = 0.d0
    st_y_rhs(:,:,mmax) = 0.d0
    st_z_rhs(:,:,mmax) = 0.d0
!    write(*,*) after advect_z: ,m,st_x_rhs(2,2,2),st_y_rhs(2,2,2),st_z_rhs(2,2,2),fmx(2,2,2),fmy(2,2,2),fmz(2,2,2),fmx(16,2,16),fmy
    return
 end if
end subroutine advect_z
!-----------------------------------------------------------------!
! Calculate an estimate of maximum velocity in fluids rest frame !
!-----------------------------------------------------------------!
!
subroutine find_v02(ext,v02,rho_b,P,n,u0,vx,vy,vz, sqBx,sqBy,sqBz, &
        lapse,shiftx,shifty,shiftz, phi,gxx,gxy,gxz,gyy,gyz,gzz)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: vx,vy,vz,u0
 real*8, dimension(ext(1),ext(2),ext(3))        :: sqBx,sqBy,sqBz
 real*8, dimension(ext(1),ext(2),ext(3))        :: rho_b,P,v02,h
 real*8, dimension(ext(1),ext(2),ext(3))        :: lapse,shiftx,shifty,shiftz
 real*8, dimension(ext(1),ext(2),ext(3))        :: gxx,gxy,gxz,gyy,gyz,gzz,phi
 real*8, dimension(ext(1),ext(2),ext(3))        :: sbt,sbx,sby,sbz
 real*8                                         :: n,gam
 real*8, parameter                              :: fac = 0.99d0
! 
! Auxiliary arrays 
!
 real*8, dimension(ext(1),ext(2),ext(3))        :: cs2,va2,b2
!
 gam = 1.d0 + 1.d0/n
 where (rho_b .gt. 0.d0)
   h = 1.d0 + (1.d0+n)*P/rho_b
 elsewhere
   h = 1.d0
 end where
 u0 = exp(4.d0*phi)*(gxx*(vx + shiftx)**2 + &
       2.d0*gxy*(vx + shiftx)*(vy + shifty) +         &
       2.d0*gxz*(vx + shiftx)*(vz + shiftz) +         &
       gyy*(vy + shifty)**2 +                           &
       2.d0*gyz*(vy + shifty)*(vz + shiftz) +         &
       gzz*(vz + shiftz)**2 )/(1.d0+lapse)**2
! *** Check for superluminal velocity ***
  where (u0 .gt. 1.d0)
    vx = (vx + shiftx)*sqrt(fac/u0)-shiftx
    vy = (vy + shifty)*sqrt(fac/u0)-shifty
    vz = (vz + shiftz)*sqrt(fac/u0)-shiftz
    u0 = fac**2
  end where
  u0 = 1.d0/(1.d0+lapse)/sqrt(1.d0-u0)
 sbt = u0*exp(4.d0*phi)*( gxx*(vx+shiftx)*sqBx + &
        gxy*( (vx+shiftx)*sqBy + (vy+shifty)*sqBx ) + &
        gxz*( (vx+shiftx)*sqBz + (vz+shiftz)*sqBx ) + &
        gyy*(vy+shifty)*sqBy + &
        gyz*( (vy+shifty)*sqBz + (vz+shiftz)*sqBy ) + &
        gzz*(vz+shiftz)*sqBz )
 sbx = sqBx/u0 + sbt*vx
 sby = sqBy/u0 + sbt*vy
 sbz = sqBz/u0 + sbt*vz
 b2 = - ((1.d0+lapse)*sbt)**2 + exp(4.d0*phi)*(gxx*(sbx + shiftx*sbt)**2 + &
        2.d0*gxy*(sbx + shiftx*sbt)*(sby + shifty*sbt) + &
        2.d0*gxz*(sbx + shiftx*sbt)*(sbz + shiftz*sbt) + &
        gyy*(sby + shifty*sbt)**2 + &
        2.d0*gyz*(sby + shifty*sbt)*(sbz + shiftz*sbt) + &
        gzz*(sbz + shiftz*sbt)**2 )
 b2 = b2/(4.d0*acos(-1.d0))*exp(-12.d0*phi)/(1.d0+lapse)**4
 where (rho_b .gt. 0.d0)
   cs2 = gam*P/(rho_b*h)
   va2 = b2/(b2 + rho_b*h)
   v02 = va2+cs2*(1.d0-va2)
 elsewhere
   v02 = b2/(b2+1.d-300)
 end where
end subroutine find_v02
subroutine sppm_shift_emfields_ul(ext,Z,Bxr,Bxl,Byr,Byl,Bzr,Bzl,m,Symmetry, &
                Sym_Bz)
 implicit none
 integer, dimension(3)                                        :: ext
 real*8, dimension(ext(1),ext(2),ext(3))                                :: Z
 real*8, dimension(ext(1),ext(2),ext(3))                :: Bxr,Byr,Bzr
 real*8, dimension(ext(1),ext(2),ext(3))                :: Bxl,Byl,Bzl
 real*8                                                        :: Sym_Bz
 integer                                                :: m, Symmetry
 integer                                                :: mmin,mmax
 integer                                                :: i,j,k
 integer, parameter  :: NO_SYMM=0, EQUATORIAL=1, OCTANT=2, PI_SYMM=3, AXISYM=4
!
 mmin = lbound(Bxl,m)
 mmax = ubound(Bxl,m)
 if (m==1) then
    do i=mmax,mmin+1,-1
       Bxl(i,:,:) = Bxl(i-1,:,:)
       Byl(i,:,:) = Byl(i-1,:,:)
       Bzl(i,:,:) = Bzl(i-1,:,:)
    end do
    if (Symmetry==OCTANT) then
       Bxl(mmin,:,:) = -Bxr(mmin,:,:)
       Byl(mmin,:,:) = Byr(mmin,:,:)
       Bzl(mmin,:,:) = Bzr(mmin,:,:)
    else if (Symmetry==AXISYM) then
       Bxl(mmin,:,:) = -Bxr(mmin,:,:)
       Byl(mmin,:,:) = -Byr(mmin,:,:)
       Bzl(mmin,:,:) = Bzr(mmin,:,:)
    end if
    return
 end if
 if (m==2) then
    do j=mmax,mmin+1,-1
       Bxl(:,j,:) = Bxl(:,j-1,:)
       Byl(:,j,:) = Byl(:,j-1,:)
       Bzl(:,j,:) = Bzl(:,j-1,:)
    end do
    if (Symmetry==OCTANT) then
       Bxl(:,mmin,:) = Bxr(:,mmin,:)
       Byl(:,mmin,:) = -Byr(:,mmin,:)
       Bzl(:,mmin,:) = Bzr(:,mmin,:)
    end if
    return
 end if
 if (m==3) then
    do k=mmax,mmin+1,-1
       Bxl(:,:,k) = Bxl(:,:,k-1)
       Byl(:,:,k) = Byl(:,:,k-1)
       Bzl(:,:,k) = Bzl(:,:,k-1)
    end do
    if (Symmetry==OCTANT .or. (Symmetry==AXISYM .and. Z(1,1,mmin) .gt. 0.d0)) then
       Bxl(:,:,mmin) = -Sym_Bz*Bxr(:,:,mmin)
       Byl(:,:,mmin) = -Sym_Bz*Byr(:,:,mmin)
       Bzl(:,:,mmin) =  Sym_Bz*Bzr(:,:,mmin)
    end if
    return
 end if
end subroutine sppm_shift_emfields_ul
!-----------------------------------------------------------------!
! Calculate an estimate of maximum velocity in fluids rest frame !
! for hybrid EOS                                                  !
!-----------------------------------------------------------------!
!
subroutine find_v02_hybrid(ext,v02,rho_b,P,u0,vx,vy,vz, &
        neos,ergo_star,ergo_sigma,rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, &
        sqBx,sqBy,sqBz, &
        lapse,shiftx,shifty,shiftz, phi,gxx,gxy,gxz,gyy,gyz,gzz, &
        cs2,va2,b2,eps_th, dPcold_drho)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: vx,vy,vz,u0
 real*8, dimension(ext(1),ext(2),ext(3))        :: sqBx,sqBy,sqBz
 real*8, dimension(ext(1),ext(2),ext(3))        :: rho_b,P,v02,h
 real*8, dimension(ext(1),ext(2),ext(3))        :: lapse,shiftx,shifty,shiftz
 real*8, dimension(ext(1),ext(2),ext(3))        :: gxx,gxy,gxz,gyy,gyz,gzz,phi
 real*8, dimension(ext(1),ext(2),ext(3))        :: sbt,sbx,sby,sbz
 integer                                        :: neos,ergo_star
 real*8, dimension(neos)                        :: rho_tab, P_tab, eps_tab
 real*8, dimension(neos+1)                        :: k_tab, gamma_tab
 real*8                                         :: gamma_th,ergo_sigma
 real*8                                         :: P_cold,eps_cold,rhob,eps
 real*8                                         :: depscold_drho
 real*8, parameter                              :: fac = 0.99d0
 logical                                        :: exit_do
 integer                                        :: i,j,k, imin,imax,jmin,jmax
 integer                                        :: kmin,kmax, ii
! 
! Auxiliary arrays 
!
 real*8, dimension(ext(1),ext(2),ext(3))        :: cs2,va2,b2
 real*8, dimension(ext(1),ext(2),ext(3))        :: eps_th, dPcold_drho
!
!  imin = lbound(rho_b,1)
!  imax = ubound(rho_b,1)
!  jmin = lbound(rho_b,2)
!  jmax = ubound(rho_b,2)
!  kmin = lbound(rho_b,3)
!  kmax = ubound(rho_b,3)
!  
!  do k=kmin,kmax
!     do j=jmin,jmax
!        do i=imin,imax
!           rhob = rho_b(i,j,k)
!           ii = 1
!           exit_do = .FALSE.
!           do
!             if (rhob .lt. rho_tab(ii)) then
!                exit_do = .TRUE.
!                P_cold = k_tab(ii)*rhob**gamma_tab(ii)
!                dPcold_drho(i,j,k) = gamma_tab(ii)*P_cold/rhob
!                if (ii==1) then
!                      eps_cold = P_cold/rhob/(gamma_tab(ii)-1.d0)
!                else
!                      eps_cold = eps_tab(ii) +  &
!                 (P_cold/rhob - P_tab(ii-1)/rho_tab(ii-1))/(gamma_tab(ii)-1.d0)
!                end if
!                  end if
!                  if (ii==neos .or. exit_do) exit
!                  ii = ii + 1
!           end do   
!           if (rhob .gt. rho_tab(neos)) then
!              P_cold = k_tab(neos+1)*rhob**gamma_tab(neos+1)
!              eps_cold = eps_tab(neos) + (P_cold/rhob - P_tab(neos) & 
!                        /rho_tab(neos))/(gamma_tab(neos+1)-1.d0)
!                   dPcold_drho(i,j,k) = gamma_tab(neos+1)*P_cold/rhob
!           end if
!           eps_th(i,j,k) = (P(i,j,k)-P_cold)/(gamma_th-1.d0)/rhob
!           eps = eps_cold + eps_th(i,j,k)
!           h(i,j,k) = 1.d0 + eps + P(i,j,k)/rhob
!        end do ! i-loop
!     end do    ! j-loop
!  end do       ! k-loop
 where (rho_b .le. rho_tab(1))
       where (rho_b .le. 0.d0)
          v02 = 0.d0
          dPcold_drho = 0.d0
          eps_th = 0.d0
          h = 1.d0
       elsewhere
          v02 = k_tab(1)*rho_b**gamma_tab(1)
          va2 = v02/rho_b/(gamma_tab(1)-1.d0)
          dPcold_drho = gamma_tab(1)*v02/rho_b
          eps_th = (P-v02)/(gamma_th-1.d0)/rho_b
          h = 1.d0 + eps_th + va2 + P/rho_b
       end where
 end where
 do ii=2,neos
    where (rho_b .le. rho_tab(ii) .and. rho_b .gt. rho_tab(ii-1))
          v02 = k_tab(ii)*rho_b**gamma_tab(ii)
          va2 = eps_tab(ii-1) + (v02/rho_b - P_tab(ii-1)/rho_tab(ii-1))/ &
                (gamma_tab(ii)-1.d0)
          dPcold_drho = gamma_tab(ii)*v02/rho_b
          eps_th = (P-v02)/(gamma_th-1.d0)/rho_b
          h = 1.d0 + eps_th + va2 + P/rho_b
    end where
 end do
if (ergo_star .eq. 0) then
 where (rho_b .gt. rho_tab(neos))
       v02 = k_tab(neos+1)*rho_b**gamma_tab(neos+1)
       va2 = eps_tab(neos) + (v02/rho_b - P_tab(neos)/rho_tab(neos))/ &
               (gamma_tab(neos+1)-1.d0)
       dPcold_drho = gamma_tab(neos+1)*v02/rho_b
       eps_th = (P-v02)/(gamma_th-1.d0)/rho_b
       h = 1.d0 + eps_th + va2 + P/rho_b
 end where
else
 where (rho_b .gt. rho_tab(neos))
       v02 = ((ergo_sigma* (1+eps_tab(neos)+P_tab(neos)/rho_tab(neos))/rho_tab(neos)**ergo_sigma) * rho_b**(ergo_sigma+1) + P_tab(n&
  &eos) - ergo_sigma*(1+eps_tab(neos))*rho_tab(neos))/(ergo_sigma+1)
       va2 = (((1+eps_tab(neos)+P_tab(neos)/rho_tab(neos))/rho_tab(neos)**ergo_sigma) * rho_b**(ergo_sigma+1) - P_tab(neos) + ergo_&
  &sigma*((1+eps_tab(neos))*rho_tab(neos)))/((ergo_sigma+1)*rho_b)-1
       dPcold_drho = ergo_sigma*(1+eps_tab(neos)+P_tab(neos)/rho_tab(neos))/rho_tab(neos)**ergo_sigma*rho_b
       eps_th = (P-v02)/(gamma_th-1.d0)/rho_b
       h = 1.d0 + eps_th + va2 + P/rho_b
 end where
end if
 u0 = exp(4.d0*phi)*(gxx*(vx + shiftx)**2 + &
       2.d0*gxy*(vx + shiftx)*(vy + shifty) +         &
       2.d0*gxz*(vx + shiftx)*(vz + shiftz) +         &
       gyy*(vy + shifty)**2 +                           &
       2.d0*gyz*(vy + shifty)*(vz + shiftz) +         &
       gzz*(vz + shiftz)**2 )/(1.d0+lapse)**2
! *** Check for superluminal velocity ***
  where (u0 .gt. 1.d0)
    vx = (vx + shiftx)*sqrt(fac/u0)-shiftx
    vy = (vy + shifty)*sqrt(fac/u0)-shifty
    vz = (vz + shiftz)*sqrt(fac/u0)-shiftz
    u0 = fac**2
  end where
  u0 = 1.d0/(1.d0+lapse)/sqrt(1.d0-u0)
 sbt = u0*exp(4.d0*phi)*( gxx*(vx+shiftx)*sqBx + &
        gxy*( (vx+shiftx)*sqBy + (vy+shifty)*sqBx ) + &
        gxz*( (vx+shiftx)*sqBz + (vz+shiftz)*sqBx ) + &
        gyy*(vy+shifty)*sqBy + &
        gyz*( (vy+shifty)*sqBz + (vz+shiftz)*sqBy ) + &
        gzz*(vz+shiftz)*sqBz )
 sbx = sqBx/u0 + sbt*vx
 sby = sqBy/u0 + sbt*vy
 sbz = sqBz/u0 + sbt*vz
 b2 = - ((1.d0+lapse)*sbt)**2 + exp(4.d0*phi)*(gxx*(sbx + shiftx*sbt)**2 + &
        2.d0*gxy*(sbx + shiftx*sbt)*(sby + shifty*sbt) + &
        2.d0*gxz*(sbx + shiftx*sbt)*(sbz + shiftz*sbt) + &
        gyy*(sby + shifty*sbt)**2 + &
        2.d0*gyz*(sby + shifty*sbt)*(sbz + shiftz*sbt) + &
        gzz*(sbz + shiftz*sbt)**2 )
 b2 = b2/(4.d0*acos(-1.d0))*exp(-12.d0*phi)/(1.d0+lapse)**4
 where (rho_b .gt. 0.d0)
   cs2 = (dPcold_drho + gamma_th*(gamma_th-1.d0)*eps_th)/h
   va2 = b2/(b2 + rho_b*h)
   v02 = va2+cs2*(1.d0-va2)
 elsewhere
   v02 = b2/(b2+1.d-300)
 end where
end subroutine find_v02_hybrid
!-----------------------------------------------------------------------------
!
! Update outer boundaries
!
!-----------------------------------------------------------------------------
subroutine matter_bc_mhd(ext, X,Y,Z, &
     rho_star, tau, mhd_st_x, mhd_st_y, mhd_st_z, &
     tau_rad, S_rad_x, S_rad_y, S_rad_z, &
     st_x, st_y, st_z, &
     rho_b, P, vx, vy, vz, Bx, By, Bz, &
     w, h, u0, Ex, Ey, Ez, sbt, sbx, sby, sbz, &
     rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
     phi, lapse, shiftx, shifty, shiftz, &
     gxx, gxy, gxz, gyy, gyz, gzz, &
     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
     Symmetry,bc, &
     neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th, rho_b_atm, rad_const, &
     have_bdry_min,have_bdry_max, enable_OS_collapse, rad_evolve_enable, rad_closure_scheme, Erad_atm_cut)
  implicit none
  integer, dimension(3)                   :: ext,have_bdry_min,have_bdry_max
  real*8, dimension(ext(1),ext(2),ext(3)) :: X,Y,Z
  real*8, dimension(ext(1),ext(2),ext(3)) :: rho_star, tau, mhd_st_x, mhd_st_y, mhd_st_z
  real*8, dimension(ext(1),ext(2),ext(3)) :: tau_rad, S_rad_x, S_rad_y, S_rad_z
  real*8, dimension(ext(1),ext(2),ext(3)) :: st_x, st_y, st_z
  real*8, dimension(ext(1),ext(2),ext(3)) :: rho_b, P, vx, vy, vz, Bx, By, Bz
  real*8, dimension(ext(1),ext(2),ext(3)) :: w, h, u0, Ex, Ey, Ez, sbt, sbx, sby, sbz
  real*8, dimension(ext(1),ext(2),ext(3)) :: rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz
  real*8, dimension(ext(1),ext(2),ext(3)) :: phi, lapse, shiftx, shifty, shiftz
  real*8, dimension(ext(1),ext(2),ext(3)) :: gxx, gxy, gxz, gyy, gyz, gzz
  real*8, dimension(ext(1),ext(2),ext(3)) :: gupxx, gupxy, gupxz, gupyy, gupyz, gupzz
  integer                                 :: Symmetry
  integer                                 :: bc
  integer                                 :: neos,ergo_star,enable_OS_collapse, rad_evolve_enable, rad_closure_scheme,compute_micro&
  &physics
  real*8, dimension(neos)                 :: rho_tab,P_tab,eps_tab
  real*8, dimension(neos+1)               :: k_tab, gamma_tab
  real*8                                  :: rho_b_atm, rad_const, gamma_th, ergo_sigma, Erad_atm_cut
!
! Other variables:
! 
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin, imax, jmax, kmax
  real*8, parameter                  :: ZERO = 0.D0
  real*8, parameter                  :: TWO = 2.D0, THREE = 3.D0
  real*8                             :: nx,ny,nz
  integer                            :: NO_SYMM, EQUATORIAL, OCTANT
  integer                            :: PI_SYMM, AXISYM
  integer                            :: FREEZE, EXTRAP, PERIODIC, OUTF, COPY, QUAD
  integer                            :: PLANAR
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(OUTF = 1, FREEZE = 2, COPY = 3, EXTRAP = 4, QUAD = 5, PLANAR = 6)
! Storage variables:
  integer, parameter                 :: m = 17
  real*8, dimension(m)               :: AUX
  real*8                             :: rho_bl, Pl, vxl, vyl, vzl, Bxl, Byl, Bzl
  real*8                             :: rho_sl, taul, mhd_st_xl, mhd_st_yl, mhd_st_zl, rho_Yel, T_fluidl
  real*8                             :: st_xl, st_yl, st_zl, wl, hl, u0l, Exl, Eyl, Ezl
  real*8                             :: sbtl, sbxl, sbyl, sbzl, rhol, Sxl, Syl, Szl
  real*8                                 :: E_radl,F_radxl,F_radyl,F_radzl, Y_el
  real*8                             :: F_rad0l, P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl,eps_totl, eps_cldl, eps&
  &_thermall, P_cldl
  real*8                                 :: tau_radl, S_rad_xl, S_rad_yl, S_rad_zl
  real*8                             :: Sxxl, Sxyl, Sxzl, Syyl, Syzl, Szzl
!
! Input translation
!
  imin = lbound(tau,1)
  jmin = lbound(tau,2)
  kmin = lbound(tau,3)
!  if(Symmetry .eq. OCTANT) then
!     if(X(1,1,1) .lt. ZERO) then
!        imin = imin + 1
!     end if
!     if(Y(1,1,1) .lt. ZERO) then
!        jmin = jmin + 1
!     end if
!     if(Z(1,1,1) .lt. ZERO) then
!        kmin = kmin + 1
!     end if
!  end if
  imax = ubound(tau,1)
  jmax = ubound(tau,2)
  kmax = ubound(tau,3)
 write(*,*) "Inside MATTER_BC_MHD!!!!!"
  !first apply boundary conditions on the hydro primitive variables
  if(bc==FREEZE) then
     write(*,*) 'Warning:  Freeze boundary conditions are not supported!'
     stop
  else if(bc==OUTF) then
     if (have_bdry_max(1)==1) then
        rho_b(imax,:,:) = rho_b(imax-1,:,:)
        P(imax,:,:) = P(imax-1,:,:)
        vx(imax,:,:) = vx(imax-1,:,:)
        vy(imax,:,:) = vy(imax-1,:,:)
        vz(imax,:,:) = vz(imax-1,:,:)
        where (vx(imax,:,:) .lt. ZERO)
              vx(imax,:,:) = ZERO
        end where
     end if
     if (have_bdry_min(1)==1) then
        rho_b(imin,:,:) = rho_b(imin+1,:,:)
        P(imin,:,:) = P(imin+1,:,:)
        vx(imin,:,:) = vx(imin+1,:,:)
        vy(imin,:,:) = vy(imin+1,:,:)
        vz(imin,:,:) = vz(imin+1,:,:)
       where (vx(imin,:,:) .gt. ZERO)
             vx(imin,:,:) = ZERO
       end where
     end if
     if (have_bdry_max(2)==1) then
        rho_b(:,jmax,:) = rho_b(:,jmax-1,:)
        P(:,jmax,:) = P(:,jmax-1,:)
        vx(:,jmax,:) = vx(:,jmax-1,:)
        vy(:,jmax,:) = vy(:,jmax-1,:)
        vz(:,jmax,:) = vz(:,jmax-1,:)
        where (vy(:,jmax,:) .lt. ZERO)
              vy(:,jmax,:) = ZERO
        end where
     end if
     if (have_bdry_min(2)==1) then
        rho_b(:,jmin,:) = rho_b(:,jmin+1,:)
        P(:,jmin,:) = P(:,jmin+1,:)
        vx(:,jmin,:) = vx(:,jmin+1,:)
        vy(:,jmin,:) = vy(:,jmin+1,:)
        vz(:,jmin,:) = vz(:,jmin+1,:)
        where (vy(:,jmin,:) .gt. ZERO)
              vy(:,jmin,:) = ZERO
        end where
     end if
     if (have_bdry_max(3)==1) then
        rho_b(:,:,kmax) = rho_b(:,:,kmax-1)
        P(:,:,kmax) = P(:,:,kmax-1)
        vx(:,:,kmax) = vx(:,:,kmax-1)
        vy(:,:,kmax) = vy(:,:,kmax-1)
        vz(:,:,kmax) = vz(:,:,kmax-1)
        where (vz(:,:,kmax) .lt. ZERO)
              vz(:,:,kmax) = ZERO
        end where
     end if
     if (have_bdry_min(3)==1) then
        rho_b(:,:,kmin) = rho_b(:,:,kmin+1)
        P(:,:,kmin) = P(:,:,kmin+1)
        vx(:,:,kmin) = vx(:,:,kmin+1)
        vy(:,:,kmin) = vy(:,:,kmin+1)
        vz(:,:,kmin) = vz(:,:,kmin+1)
        where (vz(:,:,kmin) .gt. ZERO)
              vz(:,:,kmin) = ZERO
        end where
     end if
  else if(bc==EXTRAP) then
     if (have_bdry_max(1)==1) then
        rho_b(imax,:,:) = TWO*rho_b(imax-1,:,:) &
             - rho_b(imax-2,:,:)
        P(imax,:,:) = TWO*P(imax-1,:,:) &
             - P(imax-2,:,:)
        vx(imax,:,:) = TWO*vx(imax-1,:,:) &
             - vx(imax-2,:,:)
        vy(imax,:,:) = TWO*vy(imax-1,:,:) &
             - vy(imax-2,:,:)
        vz(imax,:,:) = TWO*vz(imax-1,:,:) &
             - vz(imax-2,:,:)
     end if
     if (have_bdry_max(2)==1) then
        rho_b(:,jmax,:) = TWO*rho_b(:,jmax-1,:) &
             - rho_b(:,jmax-2,:)
        P(:,jmax,:) = TWO*P(:,jmax-1,:) &
             - P(:,jmax-2,:)
        vx(:,jmax,:) = TWO*vx(:,jmax-1,:) &
             - vx(:,jmax-2,:)
        vy(:,jmax,:) = TWO*vy(:,jmax-1,:) &
             - vy(:,jmax-2,:)
        vz(:,jmax,:) = TWO*vz(:,jmax-1,:) &
             - vz(:,jmax-2,:)
     end if
     if (have_bdry_max(3)==1) then
        rho_b(:,:,kmax) = TWO*rho_b(:,:,kmax-1) &
             - rho_b(:,:,kmax-2)
        P(:,:,kmax) = TWO*P(:,:,kmax-1) &
             - P(:,:,kmax-2)
        vx(:,:,kmax) = TWO*vx(:,:,kmax-1) &
             - vx(:,:,kmax-2)
        vy(:,:,kmax) = TWO*vy(:,:,kmax-1) &
             - vy(:,:,kmax-2)
        vz(:,:,kmax) = TWO*vz(:,:,kmax-1) &
             - vz(:,:,kmax-2)
     end if
     if (have_bdry_min(1)==1) then
        rho_b(imin,:,:) = TWO*rho_b(imin+1,:,:) &
             - rho_b(imin+2,:,:)
        P(imin,:,:) = TWO*P(imin+1,:,:) &
             - P(imin+2,:,:)
        vx(imin,:,:) = TWO*vx(imin+1,:,:) &
             - vx(imin+2,:,:)
        vy(imin,:,:) = TWO*vy(imin+1,:,:) &
             - vy(imin+2,:,:)
        vz(imin,:,:) = TWO*vz(imin+1,:,:) &
             - vz(imin+2,:,:)
     end if
     if (have_bdry_min(2)==1) then
        rho_b(:,jmin,:) = TWO*rho_b(:,jmin+1,:) &
             - rho_b(:,jmin+2,:)
        P(:,jmin,:) = TWO*P(:,jmin+1,:) &
             - P(:,jmin+2,:)
        vx(:,jmin,:) = TWO*vx(:,jmin+1,:) &
             - vx(:,jmin+2,:)
        vy(:,jmin,:) = TWO*vy(:,jmin+1,:) &
             - vy(:,jmin+2,:)
        vz(:,jmin,:) = TWO*vz(:,jmin+1,:) &
             - vz(:,jmin+2,:)
     end if
     if (have_bdry_min(3)==1) then
        rho_b(:,:,kmin) = TWO*rho_b(:,:,kmin+1) &
             - rho_b(:,:,kmin+2)
        P(:,:,kmin) = TWO*P(:,:,kmin+1) &
             - P(:,:,kmin+2)
        vx(:,:,kmin) = TWO*vx(:,:,kmin+1) &
             - vx(:,:,kmin+2)
        vy(:,:,kmin) = TWO*vy(:,:,kmin+1) &
             - vy(:,:,kmin+2)
        vz(:,:,kmin) = TWO*vz(:,:,kmin+1) &
             - vz(:,:,kmin+2)
     end if
  elseif(bc==PLANAR) then
     do k=kmin+1,kmax-1
        do j=jmin,jmax
           do i=imin,imax
              if(i==imax) then
                 if(j==jmax) then
                    rho_b(i,j,k) = rho_b(i-1,j-1,k)
                    P(i,j,k) = P(i-1,j-1,k)
                    vx(i,j,k) = vx(i-1,j-1,k)
                    vy(i,j,k) = vy(i-1,j-1,k)
                    vz(i,j,k) = vz(i-1,j-1,k)
                 else if(j==jmax-1) then
                    rho_b(i,j,k) = rho_b(i-1,j,k)
                    P(i,j,k) = P(i-1,j,k)
                    vx(i,j,k) = vx(i-1,j,k)
                    vy(i,j,k) = vy(i-1,j,k)
                    vz(i,j,k) = vz(i-1,j,k)
                 else
                    rho_b(i,j,k) = rho_b(i-1,j+1,k)
                    P(i,j,k) = P(i-1,j+1,k)
                    vx(i,j,k) = vx(i-1,j+1,k)
                    vy(i,j,k) = vy(i-1,j+1,k)
                    vz(i,j,k) = vz(i-1,j+1,k)
                 end if
              else if(i==imin) then
                 if(j==jmin) then
                    rho_b(i,j,k) = rho_b(i+1,j+1,k)
                    P(i,j,k) = P(i+1,j+1,k)
                    vx(i,j,k) = vx(i+1,j+1,k)
                    vy(i,j,k) = vy(i+1,j+1,k)
                    vz(i,j,k) = vz(i+1,j+1,k)
                 else if(j==jmin+1) then
                    rho_b(i,j,k) = rho_b(i+1,j,k)
                    P(i,j,k) = P(i+1,j,k)
                    vx(i,j,k) = vx(i+1,j,k)
                    vy(i,j,k) = vy(i+1,j,k)
                    vz(i,j,k) = vz(i+1,j,k)
                 else
                    rho_b(i,j,k) = rho_b(i+1,j-1,k)
                    P(i,j,k) = P(i+1,j-1,k)
                    vx(i,j,k) = vx(i+1,j-1,k)
                    vy(i,j,k) = vy(i+1,j-1,k)
                    vz(i,j,k) = vz(i+1,j-1,k)
                 end if
              else if(j==jmin) then
                 if(i==imin+1) then
                    rho_b(i,j,k) = rho_b(i,j+1,k)
                    P(i,j,k) = P(i,j+1,k)
                    vx(i,j,k) = vx(i,j+1,k)
                    vy(i,j,k) = vy(i,j+1,k)
                    vz(i,j,k) = vz(i,j+1,k)
                 else
                    rho_b(i,j,k) = rho_b(i-1,j+1,k)
                    P(i,j,k) = P(i-1,j+1,k)
                    vx(i,j,k) = vx(i-1,j+1,k)
                    vy(i,j,k) = vy(i-1,j+1,k)
                    vz(i,j,k) = vz(i-1,j+1,k)
                 end if
              else if(j==jmax) then
                 if(i==imax-1) then
                    rho_b(i,j,k) = rho_b(i,j-1,k)
                    P(i,j,k) = P(i,j-1,k)
                    vx(i,j,k) = vx(i,j-1,k)
                    vy(i,j,k) = vy(i,j-1,k)
                    vz(i,j,k) = vz(i,j-1,k)
                 else
                    rho_b(i,j,k) = rho_b(i+1,j-1,k)
                    P(i,j,k) = P(i+1,j-1,k)
                    vx(i,j,k) = vx(i+1,j-1,k)
                    vy(i,j,k) = vy(i+1,j-1,k)
                    vz(i,j,k) = vz(i+1,j-1,k)
                 end if
              end if
           end do
        end do
     end do
  else if (bc==COPY) then
     if (have_bdry_max(1)==1) then
        rho_b(imax,:,:) = rho_b(imax-1,:,:)
        P(imax,:,:) = P(imax-1,:,:)
        vx(imax,:,:) = vx(imax-1,:,:)
        vy(imax,:,:) = vy(imax-1,:,:)
        vz(imax,:,:) = vz(imax-1,:,:)
     end if
     if (have_bdry_min(1)==1) then
        rho_b(imin,:,:) = rho_b(imin+1,:,:)
        P(imin,:,:) = P(imin+1,:,:)
        vx(imin,:,:) = vx(imin+1,:,:)
        vy(imin,:,:) = vy(imin+1,:,:)
        vz(imin,:,:) = vz(imin+1,:,:)
     end if
     if (have_bdry_max(2)==1) then
        rho_b(:,jmax,:) = rho_b(:,jmax-1,:)
        P(:,jmax,:) = P(:,jmax-1,:)
        vx(:,jmax,:) = vx(:,jmax-1,:)
        vy(:,jmax,:) = vy(:,jmax-1,:)
        vz(:,jmax,:) = vz(:,jmax-1,:)
     end if
     if (have_bdry_min(2)==1) then
        rho_b(:,jmin,:) = rho_b(:,jmin+1,:)
        P(:,jmin,:) = P(:,jmin+1,:)
        vx(:,jmin,:) = vx(:,jmin+1,:)
        vy(:,jmin,:) = vy(:,jmin+1,:)
        vz(:,jmin,:) = vz(:,jmin+1,:)
     end if
     if (have_bdry_max(3)==1) then
        rho_b(:,:,kmax) = rho_b(:,:,kmax-1)
        P(:,:,kmax) = P(:,:,kmax-1)
        vx(:,:,kmax) = vx(:,:,kmax-1)
        vy(:,:,kmax) = vy(:,:,kmax-1)
        vz(:,:,kmax) = vz(:,:,kmax-1)
     end if
     if (have_bdry_min(3)==1) then
        rho_b(:,:,kmin) = rho_b(:,:,kmin+1)
        P(:,:,kmin) = P(:,:,kmin+1)
        vx(:,:,kmin) = vx(:,:,kmin+1)
        vy(:,:,kmin) = vy(:,:,kmin+1)
        vz(:,:,kmin) = vz(:,:,kmin+1)
     end if
  else if(bc==QUAD) then
     if (have_bdry_max(1)==1) then
        rho_b(imax,:,:) = THREE*rho_b(imax-1,:,:) &
             - THREE*rho_b(imax-2,:,:) + rho_b(imax-3,:,:)
        P(imax,:,:) = THREE*P(imax-1,:,:) &
             - THREE*P(imax-2,:,:) + P(imax-3,:,:)
        vx(imax,:,:) = THREE*vx(imax-1,:,:) &
             - THREE*vx(imax-2,:,:) + vx(imax-3,:,:)
        vy(imax,:,:) = THREE*vy(imax-1,:,:) &
             - THREE*vy(imax-2,:,:) + vy(imax-3,:,:)
        vz(imax,:,:) = THREE*vz(imax-1,:,:) &
             - THREE*vz(imax-2,:,:) + vz(imax-3,:,:)
     end if
     if (have_bdry_min(1)==1) then
        rho_b(imin,:,:) = THREE*rho_b(imin+1,:,:) &
             - THREE*rho_b(imin+2,:,:) + rho_b(imin+3,:,:)
        P(imin,:,:) = THREE*P(imin+1,:,:) &
             - THREE*P(imin+2,:,:) + P(imin+3,:,:)
        vx(imin,:,:) = THREE*vx(imin+1,:,:) &
             - THREE*vx(imin+2,:,:) + vx(imin+3,:,:)
        vy(imin,:,:) = THREE*vy(imin+1,:,:) &
             - THREE*vy(imin+2,:,:) + vy(imin+3,:,:)
        vz(imin,:,:) = THREE*vz(imin+1,:,:) &
             - THREE*vz(imin+2,:,:) + vz(imin+3,:,:)
     end if
     if (have_bdry_max(2)==1) then
        rho_b(:,jmax,:) = THREE*rho_b(:,jmax-1,:) &
             - THREE*rho_b(:,jmax-2,:) + rho_b(:,jmax-3,:)
        P(:,jmax,:) = THREE*P(:,jmax-1,:) &
             - THREE*P(:,jmax-2,:) + P(:,jmax-3,:)
        vx(:,jmax,:) = THREE*vx(:,jmax-1,:) &
             - THREE*vx(:,jmax-2,:) + vx(:,jmax-3,:)
        vy(:,jmax,:) = THREE*vy(:,jmax-1,:) &
             - THREE*vy(:,jmax-2,:) + vy(:,jmax-3,:)
        vz(:,jmax,:) = THREE*vz(:,jmax-1,:) &
             - THREE*vz(:,jmax-2,:) + vz(:,jmax-3,:)
     end if
     if (have_bdry_min(2)==1) then
        rho_b(:,jmin,:) = THREE*rho_b(:,jmin+1,:) &
             - THREE*rho_b(:,jmin+2,:) + rho_b(:,jmin+3,:)
        P(:,jmin,:) = THREE*P(:,jmin+1,:) &
             - THREE*P(:,jmin+2,:) + P(:,jmin+3,:)
        vx(:,jmin,:) = THREE*vx(:,jmin+1,:) &
             - THREE*vx(:,jmin+2,:) + vx(:,jmin+3,:)
        vy(:,jmin,:) = THREE*vy(:,jmin+1,:) &
             - THREE*vy(:,jmin+2,:) + vy(:,jmin+3,:)
        vz(:,jmin,:) = THREE*vz(:,jmin+1,:) &
             - THREE*vz(:,jmin+2,:) + vz(:,jmin+3,:)
     end if
     if (have_bdry_max(3)==1) then
        rho_b(:,:,kmax) = THREE*rho_b(:,:,kmax-1) &
             - THREE*rho_b(:,:,kmax-2) + rho_b(:,:,kmax-3)
        P(:,:,kmax) = THREE*P(:,:,kmax-1) &
             - THREE*P(:,:,kmax-2) + P(:,:,kmax-3)
        vx(:,:,kmax) = THREE*vx(:,:,kmax-1) &
             - THREE*vx(:,:,kmax-2) + vx(:,:,kmax-3)
        vy(:,:,kmax) = THREE*vy(:,:,kmax-1) &
             - THREE*vy(:,:,kmax-2) + vy(:,:,kmax-3)
        vz(:,:,kmax) = THREE*vz(:,:,kmax-1) &
             - THREE*vz(:,:,kmax-2) + vz(:,:,kmax-3)
     end if
     if (have_bdry_min(3)==1) then
        rho_b(:,:,kmin) = THREE*rho_b(:,:,kmin+1) &
             - THREE*rho_b(:,:,kmin+2) + rho_b(:,:,kmin+3)
        P(:,:,kmin) = THREE*P(:,:,kmin+1) &
             - THREE*P(:,:,kmin+2) + P(:,:,kmin+3)
        vx(:,:,kmin) = THREE*vx(:,:,kmin+1) &
             - THREE*vx(:,:,kmin+2) + vx(:,:,kmin+3)
        vy(:,:,kmin) = THREE*vy(:,:,kmin+1) &
             - THREE*vy(:,:,kmin+2) + vy(:,:,kmin+3)
        vz(:,:,kmin) = THREE*vz(:,:,kmin+1) &
             - THREE*vz(:,:,kmin+2) + vz(:,:,kmin+3)
     endif
  end if
  !write(*,*) midway bc rho_star = ,rho_b(15,15,5),P(15,15,5),vx(15,15,5),vy(15,15,5),vz(15,15,5),have_bdry_max(3)==1,proc_kmin,Z(1,
  ! now recompute the conserved variables and the auxiliary variables.
  ! i=imax, k=kmax
  if (Symmetry == AXISYM) then
     if (have_bdry_max(1)==1) then
        i=imax
        j=2
        do k=kmin,kmax
           rho_bl = rho_b(i,j,k)
           Pl = P(i,j,k)
           vxl = vx(i,j,k)
           vyl = vy(i,j,k)
           vzl = vz(i,j,k)
           Bxl = Bx(i,j,k)
           Byl = By(i,j,k)
           Bzl = Bz(i,j,k)
           AUX(1)  = phi(i,j,k)
           AUX(2)  = lapse(i,j,k)
           AUX(3)  = shiftx(i,j,k)
           AUX(4)  = shifty(i,j,k)
           AUX(5)  = shiftz(i,j,k)
           AUX(6)  = gxx(i,j,k)
           AUX(7)  = gxy(i,j,k)
           AUX(8)  = gxz(i,j,k)
           AUX(9)  = gyy(i,j,k)
           AUX(10) = gyz(i,j,k)
           AUX(11) = gzz(i,j,k)
           AUX(12) = gupxx(i,j,k)
           AUX(13) = gupxy(i,j,k)
           AUX(14) = gupxz(i,j,k)
           AUX(15) = gupyy(i,j,k)
           AUX(16) = gupyz(i,j,k)
           AUX(17) = gupzz(i,j,k)
           call recompute_vars(rho_bl,Pl,vxl,vyl,vzl,Bxl,Byl,Bzl,AUX,m,rho_sl,taul, &
                mhd_st_xl,mhd_st_yl,mhd_st_zl,st_xl,st_yl,st_zl,wl,hl,u0l,Exl,Eyl,Ezl,&
                E_radl,F_radxl,F_radyl,F_radzl, Y_el, T_fluidl, tau_radl, S_rad_xl, S_rad_yl, S_rad_zl, rho_Yel,&
                F_rad0l, P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl, &
                eps_totl, eps_cldl, eps_thermall, P_cldl, &
                sbtl,sbxl,sbyl,sbzl,rhol,Sxl,Syl,Szl,Sxxl,Sxyl,Sxzl,Syyl,Syzl,Szzl, &
                neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, &
                gamma_tab, gamma_th,rho_b_atm, rad_const, compute_microphysics, &
                enable_OS_collapse, rad_evolve_enable, rad_closure_scheme, Erad_atm_cut)
           ! Assign primitive variables which may have been modified
           vx(i,j,k) = vxl
           vy(i,j,k) = vyl
           vz(i,j,k) = vzl
           ! Assign conserved variables
           rho_star(i,j,k) = rho_sl
           tau(i,j,k) = taul
           mhd_st_x(i,j,k) = mhd_st_xl
           mhd_st_y(i,j,k) = mhd_st_yl
           mhd_st_z(i,j,k) = mhd_st_zl
           ! Assign raditaion variables
           tau_rad(i,j,k) = tau_radl
           S_rad_x(i,j,k) = S_rad_xl
           S_rad_y(i,j,k) = S_rad_yl
           S_rad_z(i,j,k) = S_rad_zl
           ! Assign other variables
           st_x(i,j,k) = st_xl
           st_y(i,j,k) = st_yl
           st_z(i,j,k) = st_zl
           w(i,j,k) = wl
           h(i,j,k) = hl
           u0(i,j,k) = u0l
           Ex(i,j,k) = Exl
           Ey(i,j,k) = Eyl
           Ez(i,j,k) = Ezl
           sbt(i,j,k) = sbtl
           sbx(i,j,k) = sbxl
           sby(i,j,k) = sbyl
           sbz(i,j,k) = sbzl
           rho(i,j,k) = rhol
           Sx(i,j,k) = Sxl
           Sy(i,j,k) = Syl
           Sz(i,j,k) = Szl
           Sxx(i,j,k) = Sxxl
           Sxy(i,j,k) = Sxyl
           Sxz(i,j,k) = Sxzl
           Syy(i,j,k) = Syyl
           Syz(i,j,k) = Syzl
           Szz(i,j,k) = Szzl
        end do
     end if
     if (have_bdry_max(3)==1) then
        k=kmax
        j=2
        do i=imin,imax
           rho_bl = rho_b(i,j,k)
           Pl = P(i,j,k)
           vxl = vx(i,j,k)
           vyl = vy(i,j,k)
           vzl = vz(i,j,k)
           Bxl = Bx(i,j,k)
           Byl = By(i,j,k)
           Bzl = Bz(i,j,k)
           AUX(1)  = phi(i,j,k)
           AUX(2)  = lapse(i,j,k)
           AUX(3)  = shiftx(i,j,k)
           AUX(4)  = shifty(i,j,k)
           AUX(5)  = shiftz(i,j,k)
           AUX(6)  = gxx(i,j,k)
           AUX(7)  = gxy(i,j,k)
           AUX(8)  = gxz(i,j,k)
           AUX(9)  = gyy(i,j,k)
           AUX(10) = gyz(i,j,k)
           AUX(11) = gzz(i,j,k)
           AUX(12) = gupxx(i,j,k)
           AUX(13) = gupxy(i,j,k)
           AUX(14) = gupxz(i,j,k)
           AUX(15) = gupyy(i,j,k)
           AUX(16) = gupyz(i,j,k)
           AUX(17) = gupzz(i,j,k)
           call recompute_vars(rho_bl,Pl,vxl,vyl,vzl,Bxl,Byl,Bzl,AUX,m,rho_sl,taul, &
                mhd_st_xl,mhd_st_yl,mhd_st_zl,st_xl,st_yl,st_zl,wl,hl,u0l,Exl,Eyl,Ezl,&
                E_radl,F_radxl,F_radyl,F_radzl, Y_el, T_fluidl, tau_radl, S_rad_xl, S_rad_yl, S_rad_zl, rho_Yel,&
                F_rad0l, P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl, &
                eps_totl, eps_cldl, eps_thermall, P_cldl, &
                sbtl,sbxl,sbyl,sbzl,rhol,Sxl,Syl,Szl,Sxxl,Sxyl,Sxzl,Syyl,Syzl,Szzl, &
                neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, &
                gamma_tab, gamma_th,rho_b_atm, rad_const, &
                compute_microphysics, enable_OS_collapse, rad_evolve_enable, rad_closure_scheme, Erad_atm_cut)
           ! Assign primitive variables which may have been modified
           vx(i,j,k) = vxl
           vy(i,j,k) = vyl
           vz(i,j,k) = vzl
           ! Assign conserved variables
           rho_star(i,j,k) = rho_sl
           tau(i,j,k) = taul
           mhd_st_x(i,j,k) = mhd_st_xl
           mhd_st_y(i,j,k) = mhd_st_yl
           mhd_st_z(i,j,k) = mhd_st_zl
           ! Assign raditaion variables                                                                                             
           tau_rad(i,j,k) = tau_radl
           S_rad_x(i,j,k) = S_rad_xl
           S_rad_y(i,j,k) = S_rad_yl
           S_rad_z(i,j,k) = S_rad_zl
           ! Assign other variables
           st_x(i,j,k) = st_xl
           st_y(i,j,k) = st_yl
           st_z(i,j,k) = st_zl
           w(i,j,k) = wl
           h(i,j,k) = hl
           u0(i,j,k) = u0l
           Ex(i,j,k) = Exl
           Ey(i,j,k) = Eyl
           Ez(i,j,k) = Ezl
           sbt(i,j,k) = sbtl
           sbx(i,j,k) = sbxl
           sby(i,j,k) = sbyl
           sbz(i,j,k) = sbzl
           rho(i,j,k) = rhol
           Sx(i,j,k) = Sxl
           Sy(i,j,k) = Syl
           Sz(i,j,k) = Szl
           Sxx(i,j,k) = Sxxl
           Sxy(i,j,k) = Sxyl
           Sxz(i,j,k) = Sxzl
           Syy(i,j,k) = Syyl
           Syz(i,j,k) = Syzl
           Szz(i,j,k) = Szzl
        end do
     end if
     if (have_bdry_min(3)==1) then
        ! k = kmin
        k=kmin
        j=2
        do i=imin,imax
           rho_bl = rho_b(i,j,k)
           Pl = P(i,j,k)
           vxl = vx(i,j,k)
           vyl = vy(i,j,k)
           vzl = vz(i,j,k)
           Bxl = Bx(i,j,k)
           Byl = By(i,j,k)
           Bzl = Bz(i,j,k)
           AUX(1)  = phi(i,j,k)
           AUX(2)  = lapse(i,j,k)
           AUX(3)  = shiftx(i,j,k)
           AUX(4)  = shifty(i,j,k)
           AUX(5)  = shiftz(i,j,k)
           AUX(6)  = gxx(i,j,k)
           AUX(7)  = gxy(i,j,k)
           AUX(8)  = gxz(i,j,k)
           AUX(9)  = gyy(i,j,k)
           AUX(10) = gyz(i,j,k)
           AUX(11) = gzz(i,j,k)
           AUX(12) = gupxx(i,j,k)
           AUX(13) = gupxy(i,j,k)
           AUX(14) = gupxz(i,j,k)
           AUX(15) = gupyy(i,j,k)
           AUX(16) = gupyz(i,j,k)
           AUX(17) = gupzz(i,j,k)
           call recompute_vars(rho_bl,Pl,vxl,vyl,vzl,Bxl,Byl,Bzl,AUX,m,rho_sl,taul, &
                mhd_st_xl,mhd_st_yl,mhd_st_zl,st_xl,st_yl,st_zl,wl,hl,u0l,Exl,Eyl,Ezl,&
                E_radl,F_radxl,F_radyl,F_radzl, Y_el, T_fluidl, tau_radl, S_rad_xl, S_rad_yl, S_rad_zl, rho_Yel,&
                F_rad0l, P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl, &
                eps_totl, eps_cldl, eps_thermall, P_cldl, &
                sbtl,sbxl,sbyl,sbzl,rhol,Sxl,Syl,Szl,Sxxl,Sxyl,Sxzl,Syyl,Syzl,Szzl, &
                neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, &
                gamma_tab, gamma_th,rho_b_atm, rad_const, &
                compute_microphysics, enable_OS_collapse, rad_evolve_enable, rad_closure_scheme, Erad_atm_cut)
           ! Assign primitive variables which may have been modified
           vx(i,j,k) = vxl
           vy(i,j,k) = vyl
           vz(i,j,k) = vzl
           ! Assign conserved variables
           rho_star(i,j,k) = rho_sl
           tau(i,j,k) = taul
           mhd_st_x(i,j,k) = mhd_st_xl
           mhd_st_y(i,j,k) = mhd_st_yl
           mhd_st_z(i,j,k) = mhd_st_zl
           ! Assign raditaion variables                                                                                             
           tau_rad(i,j,k) = tau_radl
           S_rad_x(i,j,k) = S_rad_xl
           S_rad_y(i,j,k) = S_rad_yl
           S_rad_z(i,j,k) = S_rad_zl
           ! Assign other variables
           st_x(i,j,k) = st_xl
           st_y(i,j,k) = st_yl
           st_z(i,j,k) = st_zl
           w(i,j,k) = wl
           h(i,j,k) = hl
           u0(i,j,k) = u0l
           Ex(i,j,k) = Exl
           Ey(i,j,k) = Eyl
           Ez(i,j,k) = Ezl
           sbt(i,j,k) = sbtl
           sbx(i,j,k) = sbxl
           sby(i,j,k) = sbyl
           sbz(i,j,k) = sbzl
           rho(i,j,k) = rhol
           Sx(i,j,k) = Sxl
           Sy(i,j,k) = Syl
           Sz(i,j,k) = Szl
           Sxx(i,j,k) = Sxxl
           Sxy(i,j,k) = Sxyl
           Sxz(i,j,k) = Sxzl
           Syy(i,j,k) = Syyl
           Syz(i,j,k) = Syzl
           Szz(i,j,k) = Szzl
        end do
     end if
  else
     if (have_bdry_max(1)==1) then
        i=imax
        do k=kmin,kmax
           do j=jmin,jmax
              rho_bl = rho_b(i,j,k)
              Pl = P(i,j,k)
              vxl = vx(i,j,k)
              vyl = vy(i,j,k)
              vzl = vz(i,j,k)
              Bxl = Bx(i,j,k)
              Byl = By(i,j,k)
              Bzl = Bz(i,j,k)
              AUX(1)  = phi(i,j,k)
              AUX(2)  = lapse(i,j,k)
              AUX(3)  = shiftx(i,j,k)
              AUX(4)  = shifty(i,j,k)
              AUX(5)  = shiftz(i,j,k)
              AUX(6)  = gxx(i,j,k)
              AUX(7)  = gxy(i,j,k)
              AUX(8)  = gxz(i,j,k)
              AUX(9)  = gyy(i,j,k)
              AUX(10) = gyz(i,j,k)
              AUX(11) = gzz(i,j,k)
              AUX(12) = gupxx(i,j,k)
              AUX(13) = gupxy(i,j,k)
              AUX(14) = gupxz(i,j,k)
              AUX(15) = gupyy(i,j,k)
              AUX(16) = gupyz(i,j,k)
              AUX(17) = gupzz(i,j,k)
              call recompute_vars(rho_bl,Pl,vxl,vyl,vzl,Bxl,Byl,Bzl,AUX,m,rho_sl,taul, &
                   mhd_st_xl,mhd_st_yl,mhd_st_zl,st_xl,st_yl,st_zl,wl,hl,u0l,Exl,Eyl,Ezl,&
                   E_radl,F_radxl,F_radyl,F_radzl, Y_el, T_fluidl, tau_radl, S_rad_xl, S_rad_yl, S_rad_zl, rho_Yel,&
                   F_rad0l, P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl, &
                   eps_totl, eps_cldl, eps_thermall, P_cldl, &
                   sbtl,sbxl,sbyl,sbzl,rhol,Sxl,Syl,Szl,Sxxl,Sxyl,Sxzl,Syyl,Syzl,Szzl, &
                   neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, &
                   gamma_tab, gamma_th,rho_b_atm, rad_const, &
                   compute_microphysics, enable_OS_collapse, rad_evolve_enable, rad_closure_scheme, Erad_atm_cut)
              ! Assign primitive variables which may have been modified
              vx(i,j,k) = vxl
              vy(i,j,k) = vyl
              vz(i,j,k) = vzl
              ! Assign conserved variables
              rho_star(i,j,k) = rho_sl
              tau(i,j,k) = taul
              mhd_st_x(i,j,k) = mhd_st_xl
              mhd_st_y(i,j,k) = mhd_st_yl
              mhd_st_z(i,j,k) = mhd_st_zl
              ! Assign raditaion variables                                                                                          
              tau_rad(i,j,k) = tau_radl
              S_rad_x(i,j,k) = S_rad_xl
              S_rad_y(i,j,k) = S_rad_yl
              S_rad_z(i,j,k) = S_rad_zl
              ! Assign other variables
              st_x(i,j,k) = st_xl
              st_y(i,j,k) = st_yl
              st_z(i,j,k) = st_zl
              w(i,j,k) = wl
              h(i,j,k) = hl
              u0(i,j,k) = u0l
              Ex(i,j,k) = Exl
              Ey(i,j,k) = Eyl
              Ez(i,j,k) = Ezl
              sbt(i,j,k) = sbtl
              sbx(i,j,k) = sbxl
              sby(i,j,k) = sbyl
              sbz(i,j,k) = sbzl
              rho(i,j,k) = rhol
              Sx(i,j,k) = Sxl
              Sy(i,j,k) = Syl
              Sz(i,j,k) = Szl
              Sxx(i,j,k) = Sxxl
              Sxy(i,j,k) = Sxyl
              Sxz(i,j,k) = Sxzl
              Syy(i,j,k) = Syyl
              Syz(i,j,k) = Syzl
              Szz(i,j,k) = Szzl
           end do
        end do
     end if
     if (have_bdry_min(1)==1) then
        i=imin
        do k=kmin,kmax
           do j=jmin,jmax
              rho_bl = rho_b(i,j,k)
              Pl = P(i,j,k)
              vxl = vx(i,j,k)
              vyl = vy(i,j,k)
              vzl = vz(i,j,k)
              Bxl = Bx(i,j,k)
              Byl = By(i,j,k)
              Bzl = Bz(i,j,k)
              AUX(1)  = phi(i,j,k)
              AUX(2)  = lapse(i,j,k)
              AUX(3)  = shiftx(i,j,k)
              AUX(4)  = shifty(i,j,k)
              AUX(5)  = shiftz(i,j,k)
              AUX(6)  = gxx(i,j,k)
              AUX(7)  = gxy(i,j,k)
              AUX(8)  = gxz(i,j,k)
              AUX(9)  = gyy(i,j,k)
              AUX(10) = gyz(i,j,k)
              AUX(11) = gzz(i,j,k)
              AUX(12) = gupxx(i,j,k)
              AUX(13) = gupxy(i,j,k)
              AUX(14) = gupxz(i,j,k)
              AUX(15) = gupyy(i,j,k)
              AUX(16) = gupyz(i,j,k)
              AUX(17) = gupzz(i,j,k)
              call recompute_vars(rho_bl,Pl,vxl,vyl,vzl,Bxl,Byl,Bzl,AUX,m,rho_sl,taul, &
                   mhd_st_xl,mhd_st_yl,mhd_st_zl,st_xl,st_yl,st_zl,wl,hl,u0l,Exl,Eyl,Ezl,&
                   E_radl,F_radxl,F_radyl,F_radzl, Y_el, T_fluidl, tau_radl,S_rad_xl, S_rad_yl, S_rad_zl, rho_Yel,&
                   F_rad0l, P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl, &
                   eps_totl, eps_cldl, eps_thermall, P_cldl, &
                   sbtl,sbxl,sbyl,sbzl,rhol,Sxl,Syl,Szl,Sxxl,Sxyl,Sxzl,Syyl,Syzl,Szzl, &
                   neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, &
                   gamma_tab, gamma_th,rho_b_atm, rad_const, &
                   compute_microphysics, enable_OS_collapse, rad_evolve_enable, rad_closure_scheme, Erad_atm_cut)
              ! Assign primitive variables which may have been modified
              vx(i,j,k) = vxl
              vy(i,j,k) = vyl
              vz(i,j,k) = vzl
              ! Assign conserved variables
              rho_star(i,j,k) = rho_sl
              tau(i,j,k) = taul
              mhd_st_x(i,j,k) = mhd_st_xl
              mhd_st_y(i,j,k) = mhd_st_yl
              mhd_st_z(i,j,k) = mhd_st_zl
              ! Assign raditaion variables                                                                                          
              tau_rad(i,j,k) = tau_radl
              S_rad_x(i,j,k) = S_rad_xl
              S_rad_y(i,j,k) = S_rad_yl
              S_rad_z(i,j,k) = S_rad_zl
              ! Assign other variables
              st_x(i,j,k) = st_xl
              st_y(i,j,k) = st_yl
              st_z(i,j,k) = st_zl
              w(i,j,k) = wl
              h(i,j,k) = hl
              u0(i,j,k) = u0l
              Ex(i,j,k) = Exl
              Ey(i,j,k) = Eyl
              Ez(i,j,k) = Ezl
              sbt(i,j,k) = sbtl
              sbx(i,j,k) = sbxl
              sby(i,j,k) = sbyl
              sbz(i,j,k) = sbzl
              rho(i,j,k) = rhol
              Sx(i,j,k) = Sxl
              Sy(i,j,k) = Syl
              Sz(i,j,k) = Szl
              Sxx(i,j,k) = Sxxl
              Sxy(i,j,k) = Sxyl
              Sxz(i,j,k) = Sxzl
              Syy(i,j,k) = Syyl
              Syz(i,j,k) = Syzl
              Szz(i,j,k) = Szzl
           end do
        end do
     end if
     if (have_bdry_max(2)==1) then
        j=jmax
        do k=kmin,kmax
           do i=imin,imax
              rho_bl = rho_b(i,j,k)
              Pl = P(i,j,k)
              vxl = vx(i,j,k)
              vyl = vy(i,j,k)
              vzl = vz(i,j,k)
              Bxl = Bx(i,j,k)
              Byl = By(i,j,k)
              Bzl = Bz(i,j,k)
              AUX(1)  = phi(i,j,k)
              AUX(2)  = lapse(i,j,k)
              AUX(3)  = shiftx(i,j,k)
              AUX(4)  = shifty(i,j,k)
              AUX(5)  = shiftz(i,j,k)
              AUX(6)  = gxx(i,j,k)
              AUX(7)  = gxy(i,j,k)
              AUX(8)  = gxz(i,j,k)
              AUX(9)  = gyy(i,j,k)
              AUX(10) = gyz(i,j,k)
              AUX(11) = gzz(i,j,k)
              AUX(12) = gupxx(i,j,k)
              AUX(13) = gupxy(i,j,k)
              AUX(14) = gupxz(i,j,k)
              AUX(15) = gupyy(i,j,k)
              AUX(16) = gupyz(i,j,k)
              AUX(17) = gupzz(i,j,k)
              call recompute_vars(rho_bl,Pl,vxl,vyl,vzl,Bxl,Byl,Bzl,AUX,m,rho_sl,taul, &
                   mhd_st_xl,mhd_st_yl,mhd_st_zl,st_xl,st_yl,st_zl,wl,hl,u0l,Exl,Eyl,Ezl,&
                   E_radl,F_radxl,F_radyl,F_radzl, Y_el, T_fluidl,tau_radl, S_rad_xl, S_rad_yl, S_rad_zl, rho_Yel,&
                   F_rad0l, P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl, &
                   eps_totl, eps_cldl, eps_thermall, P_cldl, &
                   sbtl,sbxl,sbyl,sbzl,rhol,Sxl,Syl,Szl,Sxxl,Sxyl,Sxzl,Syyl,Syzl,Szzl, &
                   neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab,&
                   gamma_tab, gamma_th,rho_b_atm, rad_const, &
                   compute_microphysics, enable_OS_collapse, rad_evolve_enable, rad_closure_scheme, Erad_atm_cut)
              ! Assign primitive variables which may have been modified
              vx(i,j,k) = vxl
              vy(i,j,k) = vyl
              vz(i,j,k) = vzl
              ! Assign conserved variables
              rho_star(i,j,k) = rho_sl
              tau(i,j,k) = taul
              mhd_st_x(i,j,k) = mhd_st_xl
              mhd_st_y(i,j,k) = mhd_st_yl
              mhd_st_z(i,j,k) = mhd_st_zl
              ! Assign raditaion variables                                                                                          
              tau_rad(i,j,k) = tau_radl
              S_rad_x(i,j,k) = S_rad_xl
              S_rad_y(i,j,k) = S_rad_yl
              S_rad_z(i,j,k) = S_rad_zl
              ! Assign other variables
              st_x(i,j,k) = st_xl
              st_y(i,j,k) = st_yl
              st_z(i,j,k) = st_zl
              w(i,j,k) = wl
              h(i,j,k) = hl
              u0(i,j,k) = u0l
              Ex(i,j,k) = Exl
              Ey(i,j,k) = Eyl
              Ez(i,j,k) = Ezl
              sbt(i,j,k) = sbtl
              sbx(i,j,k) = sbxl
              sby(i,j,k) = sbyl
              sbz(i,j,k) = sbzl
              rho(i,j,k) = rhol
              Sx(i,j,k) = Sxl
              Sy(i,j,k) = Syl
              Sz(i,j,k) = Szl
              Sxx(i,j,k) = Sxxl
              Sxy(i,j,k) = Sxyl
              Sxz(i,j,k) = Sxzl
              Syy(i,j,k) = Syyl
              Syz(i,j,k) = Syzl
              Szz(i,j,k) = Szzl
           end do
        end do
     end if
     if (have_bdry_min(2)==1) then
        j=jmin
        do k=kmin,kmax
           do i=imin,imax
              rho_bl = rho_b(i,j,k)
              Pl = P(i,j,k)
              vxl = vx(i,j,k)
              vyl = vy(i,j,k)
              vzl = vz(i,j,k)
              Bxl = Bx(i,j,k)
              Byl = By(i,j,k)
              Bzl = Bz(i,j,k)
              AUX(1)  = phi(i,j,k)
              AUX(2)  = lapse(i,j,k)
              AUX(3)  = shiftx(i,j,k)
              AUX(4)  = shifty(i,j,k)
              AUX(5)  = shiftz(i,j,k)
              AUX(6)  = gxx(i,j,k)
              AUX(7)  = gxy(i,j,k)
              AUX(8)  = gxz(i,j,k)
              AUX(9)  = gyy(i,j,k)
              AUX(10) = gyz(i,j,k)
              AUX(11) = gzz(i,j,k)
              AUX(12) = gupxx(i,j,k)
              AUX(13) = gupxy(i,j,k)
              AUX(14) = gupxz(i,j,k)
              AUX(15) = gupyy(i,j,k)
              AUX(16) = gupyz(i,j,k)
              AUX(17) = gupzz(i,j,k)
              call recompute_vars(rho_bl,Pl,vxl,vyl,vzl,Bxl,Byl,Bzl,AUX,m,rho_sl,taul, &
                   mhd_st_xl,mhd_st_yl,mhd_st_zl,st_xl,st_yl,st_zl,wl,hl,u0l,Exl,Eyl,Ezl,&
                   E_radl,F_radxl,F_radyl,F_radzl, Y_el, T_fluidl,tau_radl, S_rad_xl, S_rad_yl, S_rad_zl, rho_Yel,&
                   F_rad0l, P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl, &
                   eps_totl, eps_cldl, eps_thermall, P_cldl, &
                   sbtl,sbxl,sbyl,sbzl,rhol,Sxl,Syl,Szl,Sxxl,Sxyl,Sxzl,Syyl,Syzl,Szzl, &
                   neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab,k_tab, &
                   gamma_tab, gamma_th,rho_b_atm, rad_const, &
                   compute_microphysics, enable_OS_collapse, rad_evolve_enable, rad_closure_scheme, Erad_atm_cut)
              ! Assign primitive variables which may have been modified
              vx(i,j,k) = vxl
              vy(i,j,k) = vyl
              vz(i,j,k) = vzl
              ! Assign conserved variables
              rho_star(i,j,k) = rho_sl
              tau(i,j,k) = taul
              mhd_st_x(i,j,k) = mhd_st_xl
              mhd_st_y(i,j,k) = mhd_st_yl
              mhd_st_z(i,j,k) = mhd_st_zl
              ! Assign raditaion variables                                                                                          
              tau_rad(i,j,k) = tau_radl
              S_rad_x(i,j,k) = S_rad_xl
              S_rad_y(i,j,k) = S_rad_yl
              S_rad_z(i,j,k) = S_rad_zl
              ! Assign other variables
              st_x(i,j,k) = st_xl
              st_y(i,j,k) = st_yl
              st_z(i,j,k) = st_zl
              w(i,j,k) = wl
              h(i,j,k) = hl
              u0(i,j,k) = u0l
              Ex(i,j,k) = Exl
              Ey(i,j,k) = Eyl
              Ez(i,j,k) = Ezl
              sbt(i,j,k) = sbtl
              sbx(i,j,k) = sbxl
              sby(i,j,k) = sbyl
              sbz(i,j,k) = sbzl
              rho(i,j,k) = rhol
              Sx(i,j,k) = Sxl
              Sy(i,j,k) = Syl
              Sz(i,j,k) = Szl
              Sxx(i,j,k) = Sxxl
              Sxy(i,j,k) = Sxyl
              Sxz(i,j,k) = Sxzl
              Syy(i,j,k) = Syyl
              Syz(i,j,k) = Syzl
              Szz(i,j,k) = Szzl
           end do
        end do
     end if
     if (have_bdry_min(3)==1) then
        k=kmin
        do j=jmin,jmax
           do i=imin,imax
              rho_bl = rho_b(i,j,k)
              Pl = P(i,j,k)
              vxl = vx(i,j,k)
              vyl = vy(i,j,k)
              vzl = vz(i,j,k)
              Bxl = Bx(i,j,k)
              Byl = By(i,j,k)
              Bzl = Bz(i,j,k)
              AUX(1)  = phi(i,j,k)
              AUX(2)  = lapse(i,j,k)
              AUX(3)  = shiftx(i,j,k)
              AUX(4)  = shifty(i,j,k)
              AUX(5)  = shiftz(i,j,k)
              AUX(6)  = gxx(i,j,k)
              AUX(7)  = gxy(i,j,k)
              AUX(8)  = gxz(i,j,k)
              AUX(9)  = gyy(i,j,k)
              AUX(10) = gyz(i,j,k)
              AUX(11) = gzz(i,j,k)
              AUX(12) = gupxx(i,j,k)
              AUX(13) = gupxy(i,j,k)
              AUX(14) = gupxz(i,j,k)
              AUX(15) = gupyy(i,j,k)
              AUX(16) = gupyz(i,j,k)
              AUX(17) = gupzz(i,j,k)
              call recompute_vars(rho_bl,Pl,vxl,vyl,vzl,Bxl,Byl,Bzl,AUX,m,rho_sl,taul, &
                   mhd_st_xl,mhd_st_yl,mhd_st_zl,st_xl,st_yl,st_zl,wl,hl,u0l,Exl,Eyl,Ezl,&
                   E_radl,F_radxl,F_radyl,F_radzl, Y_el, T_fluidl, tau_radl, S_rad_xl, S_rad_yl, S_rad_zl, rho_Yel,&
                   F_rad0l, P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl, &
                   eps_totl, eps_cldl, eps_thermall, P_cldl, &
                   sbtl,sbxl,sbyl,sbzl,rhol,Sxl,Syl,Szl,Sxxl,Sxyl,Sxzl,Syyl,Syzl,Szzl, &
                   neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab,k_tab, &
                   gamma_tab, gamma_th,rho_b_atm, rad_const, &
                   compute_microphysics, enable_OS_collapse, rad_evolve_enable, rad_closure_scheme, Erad_atm_cut)
              ! Assign primitive variables which may have been modified
              vx(i,j,k) = vxl
              vy(i,j,k) = vyl
              vz(i,j,k) = vzl
              ! Assign conserved variables
              rho_star(i,j,k) = rho_sl
              tau(i,j,k) = taul
              mhd_st_x(i,j,k) = mhd_st_xl
              mhd_st_y(i,j,k) = mhd_st_yl
              mhd_st_z(i,j,k) = mhd_st_zl
              ! Assign raditaion variables                                                                                          
              tau_rad(i,j,k) = tau_radl
              S_rad_x(i,j,k) = S_rad_xl
              S_rad_y(i,j,k) = S_rad_yl
              S_rad_z(i,j,k) = S_rad_zl
              ! Assign other variables
              st_x(i,j,k) = st_xl
              st_y(i,j,k) = st_yl
              st_z(i,j,k) = st_zl
              w(i,j,k) = wl
              h(i,j,k) = hl
              u0(i,j,k) = u0l
              Ex(i,j,k) = Exl
              Ey(i,j,k) = Eyl
              Ez(i,j,k) = Ezl
              sbt(i,j,k) = sbtl
              sbx(i,j,k) = sbxl
              sby(i,j,k) = sbyl
              sbz(i,j,k) = sbzl
              rho(i,j,k) = rhol
              Sx(i,j,k) = Sxl
              Sy(i,j,k) = Syl
              Sz(i,j,k) = Szl
              Sxx(i,j,k) = Sxxl
              Sxy(i,j,k) = Sxyl
              Sxz(i,j,k) = Sxzl
              Syy(i,j,k) = Syyl
              Syz(i,j,k) = Syzl
              Szz(i,j,k) = Szzl
           end do
        end do
     end if
  end if
end subroutine matter_bc_mhd
subroutine recompute_vars(rho_bl,Pl,vxl,vyl,vzl,Bxl,Byl,Bzl,AUX,m,rho_sl,taul, &
                mhd_st_xl,mhd_st_yl,mhd_st_zl,st_xl,st_yl,st_zl,wl,hl,u0l,Exl,Eyl,Ezl,&
                E_radl,F_radxl,F_radyl,F_radzl, Y_el, T_fluidl, tau_radl, S_rad_xl, S_rad_yl, S_rad_zl, rho_Yel,&
                F_rad0l, P_radxxl, P_radyyl, P_radzzl, P_radxyl,P_radxzl,P_radyzl, &
                E_rad_nuel,F_radx_nuel,F_rady_nuel,F_radz_nuel, tau_rad_nuel, S_rad_x_nuel, S_rad_y_nuel, S_rad_z_nuel, &
                F_rad0_nuel, P_radxx_nuel, P_radyy_nuel, P_radzz_nuel, P_radxy_nuel, P_radxz_nuel, P_radyz_nuel, &
                E_rad_nuxl,F_radx_nuxl,F_rady_nuxl,F_radz_nuxl, tau_rad_nuxl, S_rad_x_nuxl, S_rad_y_nuxl, S_rad_z_nuxl, &
                F_rad0_nuxl, P_radxx_nuxl, P_radyy_nuxl, P_radzz_nuxl, P_radxy_nuxl, P_radxz_nuxl, P_radyz_nuxl, &
                eps_totl, eps_cldl, eps_thermall, P_cldl, &
                sbtl,sbxl,sbyl,sbzl,rhol,Sxl,Syl,Szl,Sxxl,Sxyl,Sxzl,Syyl,Syzl,Szzl, &
                neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, &
                gamma_tab, gamma_th,rho_b_atm, rad_const, &
                compute_microphysics, enable_OS_collapse, rad_evolve_enable, rad_closure_scheme, microphysics_scheme, Erad_atm_cut)
  implicit none
  integer                   :: m
  integer                   :: neos,ergo_star,enable_OS_collapse, rad_evolve_enable, rad_closure_scheme, microphysics_scheme, compu&
  &te_microphysics
  real*8, dimension(neos)   :: rho_tab,P_tab,eps_tab
  real*8, dimension(neos+1) :: k_tab, gamma_tab
  real*8                    :: gamma_th,ergo_sigma,rad_const
  real*8, dimension(m)      :: AUX
  real*8                    :: rho_bl, Pl, vxl, vyl, vzl, Bxl, Byl, Bzl
  real*8                    :: rho_sl, taul, mhd_st_xl, mhd_st_yl, mhd_st_zl
  real*8                    :: st_xl, st_yl, st_zl
  real*8                    :: E_rad_denom, E_rad_num, E_radl,F_radxl,F_radyl,F_radzl, Y_el, T_fluidl
  real*8                    :: E_rad_nuel,F_radx_nuel,F_rady_nuel,F_radz_nuel, E_rad_nuxl,F_radx_nuxl,F_rady_nuxl,F_radz_nuxl
  real*8                    :: S_rad_xl, S_rad_yl, S_rad_zl, tau_radl, rho_Yel
  real*8                    :: S_rad_x_nuel, S_rad_y_nuel, S_rad_z_nuel, tau_rad_nuel
  real*8                    :: S_rad_x_nuxl, S_rad_y_nuxl, S_rad_z_nuxl, tau_rad_nuxl
  real*8                    :: shift_x, shift_y, shift_z, v_xl, v_yl, v_zl
  real*8                    :: wl, hl, u0l, Exl, Eyl, Ezl, sbtl, sbxl, sbyl, sbzl
  real*8                    :: rhol, Sxl, Syl, Szl, Sxxl, Sxyl, Sxzl, Syyl, Syzl, Szzl
  real*8                    :: rho_b_atm, Erad_atm_cut
  !other variables
  real*8                   :: phi, lapse, shiftx, shifty, shiftz
  real*8                   :: gxx, gxy, gxz, gyy, gyz, gzz
  real*8                   :: gupxx, gupxy, gupxz, gupyy, gupyz, gupzz
  real*8                   :: Psi6, Psi4, eps, el, er, temp3, au0r1, alpn1
  real*8                   :: P_cold, eps_cold, Pmax, Pmin
  real*8                   :: P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl, P_radl, F_rad0l
  real*8                   :: P_radxx_nuel, P_radyy_nuel, P_radzz_nuel, P_radxy_nuel, P_radxz_nuel, P_radyz_nuel, P_rad_nuel, F_rad&
  &0_nuel
  real*8                   :: P_radxx_nuxl, P_radyy_nuxl, P_radzz_nuxl, P_radxy_nuxl, P_radxz_nuxl, P_radyz_nuxl, P_rad_nuxl, F_rad&
  &0_nuxl
  real*8                   :: P_rad0xl, P_rad0yl, P_rad0zl, P_rad00l
  real*8                   :: eps_totl, eps_cldl, eps_thermall, P_cldl
  real*8                   :: sb_x, sb_y, sb_z
  real*8                   :: E_x, E_y, E_z, u_x, u_y, u_z, uxl, uyl, uzl
  real*8                   :: temp_rad, temp_rad1, F_rad_xl, F_rad_yl, F_rad_zl, F_rad_0l, Fksq, Fasq, zeta, zeta_temp, zeta_cut,chi
  real*8                   :: beta2, udotbeta, g_00l, u_0l
  real*8, parameter        :: fac2 = 0.99d0
  real*8                   :: f1o4p, f1o8p
  integer                  :: i,j,k
  real*8                   :: f1o4pa,B_x,B_y,B_z,B2,sb2, fac
  real*8, parameter        :: one = 1.d0
  real*8, parameter        :: two = 2.d0
  f1o4p = 1.d0/(4.d0*acos(-1.d0))
  f1o8p = 1.d0/(8.d0*acos(-1.d0))
  ! Just in case, initialize output variables to zero
  rho_sl    = 0.d0
  taul      = 0.d0
  rho_Yel   = 0.d0
  mhd_st_xl = 0.d0
  mhd_st_yl = 0.d0
  mhd_st_zl = 0.d0
  st_xl     = 0.d0
  st_yl     = 0.d0
  st_zl     = 0.d0
  wl        = 0.d0
  hl        = 0.d0
  u0l       = 0.d0
  Exl       = 0.d0
  Eyl       = 0.d0
  Ezl       = 0.d0
  sbtl      = 0.d0
  sbxl      = 0.d0
  sbyl      = 0.d0
  sbzl      = 0.d0
  rhol      = 0.d0
  tau_radl  = 0.d0
  S_rad_xl  = 0.d0
  S_rad_yl  = 0.d0
  S_rad_zl  = 0.d0
  Sxl       = 0.d0
  Syl       = 0.d0
  Szl       = 0.d0
  Sxxl      = 0.d0
  Sxyl      = 0.d0
  Sxzl      = 0.d0
  Syyl      = 0.d0
  Syzl      = 0.d0
  Szzl      = 0.d0
  ! assign metric variables
  phi    = AUX(1)
  lapse  = AUX(2)
  shiftx = AUX(3)
  shifty = AUX(4)
  shiftz = AUX(5)
  gxx    = AUX(6)
  gxy    = AUX(7)
  gxz    = AUX(8)
  gyy    = AUX(9)
  gyz    = AUX(10)
  gzz    = AUX(11)
  gupxx  = AUX(12)
  gupxy  = AUX(13)
  gupxz  = AUX(14)
  gupyy  = AUX(15)
  gupyz  = AUX(16)
  gupzz  = AUX(17)
  alpn1 = 1.d0 + lapse
  Psi6 = exp(6.d0*phi)
  Psi4 = exp(4.d0*phi)
  f1o4pa = 1.d0/(sqrt(4.d0*acos(-1.d0))*alpn1)
  ! First check that everything is within bounds
  if (rho_bl .lt. rho_b_atm) then
     rho_bl = rho_b_atm
     call compute_pcold_epscold(rho_bl, P_cold, eps_cold, &
          neos,ergo_star, ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,enable_OS_collapse)
     Pl = P_cold
     eps = eps_cold
     hl = 1.d0 + Pl/rho_bl + eps
  else
     call compute_pcold_epscold(rho_bl, P_cold, eps_cold, &
          neos,ergo_star, ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,enable_OS_collapse)
     if (compute_microphysics.eq.1) then
        call compute_eps_T_microphys(eps, T_fluidl, eps_cold, Pl, P_cold, rho_bl)
     else
         !! Just to calculate T_fluid
!         call compute_eps_T_microphys(eps, T_fluidl, eps_cold, Pl, P_cold, rho_bl)
         eps = eps_cold + (Pl-P_cold)/(gamma_th-1.d0)/rho_bl
     end if
     hl = 1.d0 + Pl/rho_bl + eps
  end if
  Pmax = 10.d0*P_cold
  if (Pl .gt. Pmax) then
     Pl = Pmax
     if (compute_microphysics.eq.1) then
        call compute_eps_T_microphys(eps, T_fluidl, eps_cold, Pmax, P_cold, rho_bl)
     else
        !! Just to calculate T_fluid
!        call compute_eps_T_microphys(eps, T_fluidl, eps_cold, Pmax, P_cold, rho_bl)
        eps = eps_cold + (Pmax-P_cold)/(gamma_th-1.d0)/rho_bl
     end if
     hl = 1.d0 + eps + Pmax/rho_bl
  end if
  Pmin = 0.5d0*P_cold
  if (Pl .lt. Pmin) then
     Pl = Pmin
     if (compute_microphysics.eq.1) then
        call compute_eps_T_microphys(eps, T_fluidl, eps_cold, Pmin, P_cold, rho_bl)
     else
        !! Just to calculate T_fluid
!        call compute_eps_T_microphys(eps, T_fluidl, eps_cold, Pmin, P_cold, rho_bl)
        eps = eps_cold + (Pmin-P_cold)/(gamma_th-1.d0)/rho_bl
     end if
     hl = 1.d0 + eps + Pmin/rho_bl
  end if
  eps_totl = eps
  eps_cldl = eps_cold
  eps_thermall = eps - eps_cold
  P_cldl = P_cold
  ! Compute al*u0-1
  er = Psi4*(gxx*(vxl + shiftx)**2 +             &
       2.d0*gxy*(vxl + shiftx)*(vyl + shifty) +  &
       2.d0*gxz*(vxl + shiftx)*(vzl + shiftz) +  &
       gyy*(vyl + shifty)**2 +                   &
       2.d0*gyz*(vyl + shifty)*(vzl + shiftz) +  &
       gzz*(vzl + shiftz)**2 )/alpn1**2
  ! *** Check for superluminal velocity ***
  if (er .gt. 1.d0) then
     vxl = (vxl + shiftx)*sqrt(fac2/er)-shiftx
     vyl = (vyl + shifty)*sqrt(fac2/er)-shifty
     vzl = (vzl + shiftz)*sqrt(fac2/er)-shiftz
     er = fac2
  end if
  ! ***************************************
  el = sqrt(1.d0-er)
  au0r1 = er/el/(1.d0+el)
  u0l  = (au0r1+1.d0)/alpn1
  rho_sl = alpn1 * Psi6 * rho_bl * u0l
  rho_Yel = rho_sl * Y_el
  wl = alpn1 * u0l * rho_sl
  u_x = u0l*Psi4*(gxx*(vxl + shiftx) + &
       gxy*(vyl + shifty) + &
       gxz*(vzl + shiftz))
  u_y = u0l*Psi4*(gxy*(vxl + shiftx) + &
       gyy*(vyl + shifty) + &
       gyz*(vzl + shiftz))
  u_z = u0l*Psi4*(gxz*(vxl + shiftx) + &
       gyz*(vyl + shifty) + &
       gzz*(vzl + shiftz))
  st_xl = rho_sl*hl*u_x
  st_yl = rho_sl*hl*u_y
  st_zl = rho_sl*hl*u_z
  B_x  = Psi4 * (gxx*Bxl + gxy*Byl + gxz*Bzl)
  B_y  = Psi4 * (gxy*Bxl + gyy*Byl + gyz*Bzl)
  B_z  = Psi4 * (gxz*Bxl + gyz*Byl + gzz*Bzl)
  B2 = Psi4*(gxx*Bxl*Bxl + two*gxy*Bxl*Byl + &
       TWO*gxz*Bxl*Bzl + gyy*Byl*Byl + &
       TWO*gyz*Byl*Bzl + gzz*Bzl*Bzl)
  sbtl = (u_x*Bxl + u_y*Byl + u_z*Bzl)*f1o4pa
  sb2 = ( B2*f1o4pa**2 + sbtl**2)/u0l**2
  sb_x = (B_x*f1o4pa + u_x*sbtl)/u0l
  sb_y = (B_y*f1o4pa + u_y*sbtl)/u0l
  sb_z = (B_z*f1o4pa + u_z*sbtl)/u0l
  mhd_st_xl = st_xl + &
       alpn1*Psi6*u0l*sb2*u_x - alpn1*Psi6*sbtl*sb_x
  mhd_st_yl = st_yl + &
       alpn1*Psi6*u0l*sb2*u_y - alpn1*Psi6*sbtl*sb_y
  mhd_st_zl = st_zl + &
       alpn1*Psi6*u0l*sb2*u_z - alpn1*Psi6*sbtl*sb_z
  taul = (au0r1 + (Pl/rho_bl + eps)*alpn1*u0l)*rho_sl + &
       Psi6*sb2*(alpn1*u0l)**2 &
       - Psi6*(Pl+sb2*0.5d0)-Psi6*(alpn1*sbtl)**2
  ! Now calculate the sources....
  fac  = one / ( Psi6 * wl*hl )
  rhol = hl * wl/Psi6 - Pl
  Sxl  = st_xl/Psi6
  Syl  = st_yl/Psi6
  Szl  = st_zl/Psi6
  Sxxl = fac * st_xl*st_xl + Psi4 * gxx * Pl
  Sxyl = fac * st_xl*st_yl + Psi4 * gxy * Pl
  Sxzl = fac * st_xl*st_zl + Psi4 * gxz * Pl
  Syyl = fac * st_yl*st_yl + Psi4 * gyy * Pl
  Syzl = fac * st_yl*st_zl + Psi4 * gyz * Pl
  Szzl = fac * st_zl*st_zl + Psi4 * gzz * Pl
  !
  ! MHD metric sources
  E_x = Psi6/alpn1 * ( Byl*(vzl+shiftz) - Bzl*(vyl+shifty) )
  E_y = Psi6/alpn1 * ( Bzl*(vxl+shiftx) - Bxl*(vzl+shiftz) )
  E_z = Psi6/alpn1 * ( Bxl*(vyl+shifty) - Byl*(vxl+shiftx) )
  Exl = (gupxx*E_x + gupxy*E_y + gupxz*E_z)/Psi4
  Eyl = (gupxy*E_x + gupyy*E_y + gupyz*E_z)/Psi4
  Ezl = (gupxz*E_x + gupyz*E_y + gupzz*E_z)/Psi4
  temp3 = f1o8p*(Exl*E_x + Eyl*E_y + Ezl*E_z + Bxl*B_x + Byl*B_y + Bzl*B_z)
  rhol = rhol + temp3
  Sxxl = Sxxl + temp3*Psi4*gxx - f1o4p*(E_x*E_x + B_x*B_x)
  Sxyl = Sxyl + temp3*Psi4*gxy - f1o4p*(E_x*E_y + B_x*B_y)
  Sxzl = Sxzl + temp3*Psi4*gxz - f1o4p*(E_x*E_z + B_x*B_z)
  Syyl = Syyl + temp3*Psi4*gyy - f1o4p*(E_y*E_y + B_y*B_y)
  Syzl = Syzl + temp3*Psi4*gyz - f1o4p*(E_y*E_z + B_y*B_z)
  Szzl = Szzl + temp3*Psi4*gzz - f1o4p*(E_z*E_z + B_z*B_z)
  Sxl  = Sxl  + f1o4p*Psi6*(Eyl*Bzl - Ezl*Byl)
  Syl  = Syl  + f1o4p*Psi6*(Ezl*Bxl - Exl*Bzl)
  Szl  = Szl  + f1o4p*Psi6*(Exl*Byl - Eyl*Bxl)
  sbtl = u_x*Bxl + u_y*Byl + u_z*Bzl
  sbxl = Bxl/u0l + vxl*sbtl
  sbyl = Byl/u0l + vyl*sbtl
  sbzl = Bzl/u0l + vzl*sbtl
  if (enable_OS_collapse.eq.1) then
     T_fluidl = Pl/rho_bl
  end if
  ! Add radiation terms in the source.
  if (rad_evolve_enable.eq.1) then
     shift_x = Psi4 *(gxx*shiftx + gxy*shifty + gxz*shiftz)
     shift_y = Psi4 *(gxy*shiftx + gyy*shifty + gyz*shiftz)
     shift_z = Psi4 *(gxz*shiftx + gyz*shifty + gzz*shiftz)
     v_xl = Psi4 *(gxx*vxl + gxy*vyl + gxz*vzl)
     v_yl = Psi4 *(gxy*vxl + gyy*vyl + gyz*vzl)
     v_zl = Psi4 *(gxz*vxl + gyz*vyl + gzz*vzl)
     beta2 = shiftx*shift_x + shifty*shift_y + shiftz*shift_z
     udotbeta = u0l*(vxl*shift_x + vyl*shift_y + vzl*shift_z)
     g_00l =beta2-alpn1*alpn1
     u_0l = g_00l*u0l + udotbeta
     uxl = u0l*vxl
     uyl = u0l*vyl
     uzl = u0l*vzl
     if (E_radl .lt. Erad_atm_cut) then
        E_radl = Erad_atm_cut
        F_radxl = 0.0
        F_radyl = 0.0
        F_radzl = 0.0
     end if
     F_rad0l = - (F_radxl*u_x + F_radyl*u_y + F_radzl*u_z)/u_0l
     F_rad_xl = Psi4 * (gxx * F_radxl + gxy * F_radyl + gxz * F_radzl) + shift_x* F_rad0l
     F_rad_yl = Psi4 * (gxy * F_radxl + gyy * F_radyl + gyz * F_radzl) + shift_y* F_rad0l
     F_rad_zl = Psi4 * (gxz * F_radxl + gyz * F_radyl + gzz * F_radzl) + shift_z* F_rad0l
     F_rad_0l = - (F_rad_xl*uxl + F_rad_yl*uyl + F_rad_zl*uzl)/u0l
     if (enable_OS_collapse.eq.1) then
        T_fluidl = Pl/rho_bl
     end if
     if (rad_closure_scheme.eq.0) then
        P_radl = E_radl/3.d0
        temp_rad = alpn1*u0l
        temp_rad1 = temp_rad*temp_rad*(E_radl+P_radl) - P_radl + 2.d0*alpn1*u0l*F_rad0l
        tau_radl = alpn1*alpn1*Psi6*(E_radl*u0l*u0l+2.0*F_rad0l*u0l+P_radl*u0l*u0l)-Psi6*P_radl
        S_rad_xl = alpn1*Psi6*((E_radl+P_radl)*u0l*u_x + F_rad0l*u_x + F_rad_xl * u0l)
        S_rad_yl = alpn1*Psi6*((E_radl+P_radl)*u0l*u_y + F_rad0l*u_y + F_rad_yl * u0l)
        S_rad_zl = alpn1*Psi6*((E_radl+P_radl)*u0l*u_z + F_rad0l*u_z + F_rad_zl * u0l)
        rhol = rhol + temp_rad1
        Sxl = Sxl + temp_rad*( ( (E_radl + P_radl) * u0l + F_rad0l) * (shift_x + v_xl) + F_rad_xl)
        Syl = Syl + temp_rad*( ( (E_radl + P_radl) * u0l + F_rad0l) * (shift_y + v_yl) + F_rad_yl)
        Szl = Szl + temp_rad*( ( (E_radl + P_radl) * u0l + F_rad0l) * (shift_z + v_zl) + F_rad_zl)
        Sxxl = Sxxl + (E_radl+P_radl)*(u0l*(shift_x + v_xl))**2 + 2.0*F_rad_xl*u0l*(shift_x + v_xl) + Psi4 * P_radl * gxx
        Syyl = Syyl + (E_radl+P_radl)*(u0l*(shift_y + v_yl))**2 + 2.0*F_rad_yl*u0l*(shift_y + v_yl) + Psi4 * P_radl * gyy
        Szzl = Szzl + (E_radl+P_radl)*(u0l*(shift_z + v_zl))**2 + 2.0*F_rad_zl*u0l*(shift_z + v_zl) + Psi4 * P_radl * gzz
        Sxyl = Sxyl + (E_radl+P_radl)*u0l**2*(shift_x + v_xl)*(shift_y + v_yl) + u0l*(F_rad_xl*(shift_y + v_yl)+F_rad_yl*(shift_x +&
  & v_xl)) + Psi4 * P_radl * gxy
        Sxzl = Sxzl + (E_radl+P_radl)*u0l**2*(shift_x + v_xl)*(shift_z + v_zl) + u0l*(F_rad_xl*(shift_z + v_zl)+F_rad_zl*(shift_x +&
  & v_xl)) + Psi4 * P_radl * gxz
        Syzl = Syzl + (E_radl+P_radl)*u0l**2*(shift_y + v_yl)*(shift_z + v_zl) + u0l*(F_rad_yl*(shift_z + v_zl)+F_rad_zl*(shift_y +&
  & v_yl)) + Psi4 * P_radl * gyz
     else
        Fasq = F_rad_0l*F_rad0l + F_rad_xl*F_radxl +  F_rad_yl*F_radyl +  F_rad_zl*F_radzl
        zeta_temp = sqrt(abs(F_rad_0l*F_rad0l + F_rad_xl*F_radxl +  F_rad_yl*F_radyl +  F_rad_zl*F_radzl)/E_radl**2)
        !zeta_cut = Erad_atm_cut*1.5
        zeta_cut = 1.0e-40
        if (E_radl.le.zeta_cut) then
           zeta = 1.0
        else
           zeta = zeta_temp
        end if
        if (zeta .gt. 1.0) then
           zeta = 1.0;
        end if
        chi = 1/3.0 + zeta**2*(6.0-2.0*zeta+6.0*zeta**2)/15.0
        if (E_radl .lt. Erad_atm_cut) then
           P_radxxl = 0.0
           P_radyyl = 0.0
           P_radzzl = 0.0
           P_radxyl = 0.0
           P_radxzl = 0.0
           P_radyzl = 0.0
        else
           if (Fasq .le. 0) then
              P_radxxl = E_radl*(gupxx/Psi4 - shiftx*shiftx/alpn1**2 + uxl**2)/2.0*(1.0-chi)
              P_radyyl = E_radl*(gupyy/Psi4 - shifty*shifty/alpn1**2 + uyl**2)/2.0*(1.0-chi)
              P_radzzl = E_radl*(gupzz/Psi4 - shiftz*shiftz/alpn1**2 + uzl**2)/2.0*(1.0-chi)
              P_radxyl = E_radl*(gupxy/Psi4 - shiftx*shifty/alpn1**2 + uxl*uyl)/2.0*(1.0-chi)
              P_radxzl = E_radl*(gupxz/Psi4 - shiftx*shiftz/alpn1**2 + uxl*uzl)/2.0*(1.0-chi)
              P_radyzl = E_radl*(gupyz/Psi4 - shifty*shiftz/alpn1**2 + uyl*uzl)/2.0*(1.0-chi)
           else
              P_radxxl = E_radl*((F_radxl**2/Fasq)*(3.0*chi -1.0)/2.0 + (gupxx/Psi4 - shiftx*shiftx/alpn1**2 + uxl**2)/2.0*(1.0-chi&
  &))
              P_radyyl = E_radl*((F_radyl**2/Fasq)*(3.0*chi -1.0)/2.0 + (gupyy/Psi4 - shifty*shifty/alpn1**2 + uyl**2)/2.0*(1.0-chi&
  &))
              P_radzzl = E_radl*((F_radzl**2/Fasq)*(3.0*chi -1.0)/2.0 + (gupzz/Psi4 - shiftz*shiftz/alpn1**2 + uzl**2)/2.0*(1.0-chi&
  &))
              P_radxyl = E_radl*((F_radxl*F_radyl/Fasq)*(3.0*chi -1.0)/2.0 + (gupxy/Psi4 - shiftx*shifty/alpn1**2 + uxl*uyl)/2.0*(1&
  &.0-chi))
              P_radxzl = E_radl*((F_radxl*F_radzl/Fasq)*(3.0*chi -1.0)/2.0 + (gupxz/Psi4 - shiftx*shiftz/alpn1**2 + uxl*uzl)/2.0*(1&
  &.0-chi))
              P_radyzl = E_radl*((F_radyl*F_radzl/Fasq)*(3.0*chi -1.0)/2.0 + (gupyz/Psi4 - shifty*shiftz/alpn1**2 ++ uyl*uzl)/2.0*(&
  &1.0-chi))
           end if
        end if
        P_rad0xl = - (P_radxxl * u_x + P_radxyl * u_y + P_radxzl * u_z)/u_0l
        P_rad0yl = - (P_radxyl * u_x + P_radyyl * u_y + P_radyzl * u_z)/u_0l
        P_rad0zl = - (P_radxzl * u_x + P_radyzl * u_y + P_radzzl * u_z)/u_0l
        P_rad00l = - (P_rad0xl * u_x + P_rad0yl * u_y + P_rad0zl * u_z)/u_0l
        tau_radl = alpn1*alpn1*Psi6*(E_radl*u0l*u0l+2.0*F_rad0l*u0l+P_rad00l)
        S_rad_xl = alpn1*Psi6*(E_radl*u0l*u_x + F_rad0l*u_x + F_rad_xl * u0l + P_rad00l*shift_x + Psi4*(P_rad0xl*gxx + P_rad0yl*gxy&
  & + P_rad0zl*gxz))
        S_rad_yl = alpn1*Psi6*(E_radl*u0l*u_y + F_rad0l*u_y + F_rad_yl * u0l + P_rad00l*shift_y + Psi4*(P_rad0xl*gxy + P_rad0yl*gyy&
  & + P_rad0zl*gyz))
        S_rad_zl = alpn1*Psi6*(E_radl*u0l*u_z + F_rad0l*u_z + F_rad_zl * u0l + P_rad00l*shift_z + Psi4*(P_rad0xl*gxz + P_rad0yl*gyz&
  & + P_rad0zl*gzz))
        rhol = rhol + alpn1**2.0*(E_radl*u0l**2 + 2.0 * F_rad0l * u0l + P_rad00l)
        Sxl = Sxl + alpn1*(u0l*E_radl*u_x + F_rad0l*u_x + u0l*F_rad_xl +&
             P_rad00l*shift_x + Psi4*(P_rad0xl*gxx + P_rad0yl*gxy + P_rad0zl*gxz))
        Syl = Syl + alpn1*(u0l*E_radl*u_y + F_rad0l*u_y + u0l*F_rad_yl +&
             P_rad00l*shift_y + Psi4*(P_rad0xl*gxy + P_rad0yl*gyy + P_rad0zl*gyz))
        Szl = Szl + alpn1*(u0l*E_radl*u_z + F_rad0l*u_z + u0l*F_rad_zl +&
             P_rad00l*shift_z + Psi4*(P_rad0xl*gxz + P_rad0yl*gyz + P_rad0zl*gzz))
        Sxxl = Sxxl + E_radl*u_x*u_x + 2.0*F_rad_xl*u_x +&
             shift_x**2.0*P_rad00l + shift_x*2.0*(gxx*P_rad0xl+gxy*P_rad0yl+gxz*P_rad0zl) +&
             Psi4**2.0*( gxx**2.0*P_radxxl + gxy**2.0*P_radyyl + gxz**2.0*P_radzzl +&
             2.0*(gxx*gxy*P_radxyl+gxx*gxz*P_radxzl+gxy*gxz*P_radyzl) )
        Syyl = Syyl + E_radl*u_y*u_y + 2.0*F_rad_yl*u_y +&
             shift_y**2.0*P_rad00l + shift_y*2.0*(gxy*P_rad0xl+gyy*P_rad0yl+gyz*P_rad0zl) +&
             Psi4**2.0*( gxy**2.0*P_radxxl + gyy**2.0*P_radyyl + gyz**2.0*P_radzzl +&
             2.0*(gxy*gyy*P_radxyl+gxy*gyz*P_radxzl+gyy*gyz*P_radyzl) )
        Szzl = Szzl + E_radl*u_z*u_z + 2.0*F_rad_zl*u_z +&
             shift_z**2.0*P_rad00l + shift_z*2.0*(gxz*P_rad0xl+gyz*P_rad0yl+gzz*P_rad0zl) +&
             Psi4**2.0*( gxz**2.0*P_radxxl + gyz**2.0*P_radyyl + gzz**2.0*P_radzzl +&
             2.0*(gxz*gyz*P_radxyl+gxz*gzz*P_radxzl+gyz*gzz*P_radyzl) )
        Sxyl = Sxyl + E_radl*u_x*u_y + F_rad_xl*u_y + F_rad_yl*u_x +&
             shift_x*shift_y*P_rad00l +&
             Psi4*(shift_x*(gxy*P_rad0xl + gyy*P_rad0yl + gyz*P_rad0zl) + (shift_y*(gxx*P_rad0xl + gxy*P_rad0yl + gxz*P_rad0zl)) )+&
             Psi4**2.0*(gxx*gxy*P_radxxl + gxy*gyy*P_radyyl + gxz*gyz*P_radzzl +&
             (gxx*gyy + gxy*gxy)*P_radxyl + (gxx*gyz + gxz*gxy)*P_radxzl + (gxy*gyz + gxz*gyy)*P_radyzl)
        Sxzl = Sxzl + E_radl*u_x*u_z + F_rad_xl*u_z + F_rad_zl*u_x +&
             shift_x*shift_z*P_rad00l +&
             Psi4*(shift_x*(gxz*P_rad0xl + gyz*P_rad0yl + gzz*P_rad0zl) + (shift_z*(gxx*P_rad0xl + gxy*P_rad0yl + gxz*P_rad0zl)) )+&
             Psi4**2.0*(gxx*gxz*P_radxxl + gxy*gyz*P_radyyl + gxz*gzz*P_radzzl +&
             (gxx*gyz + gxy*gxz)*P_radxyl + (gxx*gzz + gxz*gxz)*P_radxzl + (gxy*gzz + gxz*gyz)*P_radyzl)
        Syzl = Syzl + E_radl*u_y*u_z + F_rad_yl*u_z + F_rad_zl*u_y +&
             shift_y*shift_z*P_rad00l +&
             Psi4*(shift_y*(gxz*P_rad0xl + gyz*P_rad0yl + gzz*P_rad0zl) + (shift_z*(gxy*P_rad0xl + gyy*P_rad0yl + gyz*P_rad0zl)) )+&
             Psi4**2.0*(gxy*gxz*P_radxxl + gyy*gyz*P_radyyl + gyz*gzz*P_radzzl +&
             (gxy*gyz + gxz*gyy)*P_radxyl + (gxy*gzz + gyz*gxz)*P_radxzl + (gyy*gzz + gyz*gyz)*P_radyzl)
     end if
  end if
if(1==0) then
     call rad_sources(Psi4, Psi6, alpn1, shiftx, shifty, shiftz, vxl, vyl, vzl, &
          u0l, u_x, u_y, u_z, &
          gxx, gyy, gzz, gxy, gxz, gyz, &
          gupxx, gupyy, gupzz, gupxy, gupxz, gupyz,&
          E_radl, F_radxl, F_radyl, F_radzl, F_rad0l, &
          tau_radl, S_rad_xl, S_rad_yl, S_rad_zl, &
          P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl, &
          rhol, Sxl, Syl, Szl, Sxxl, Syyl, Szzl, Sxyl, Sxzl, Syzl, &
          rad_closure_scheme, Erad_atm_cut)
  end if
     if (microphysics_scheme==1) then
        call rad_sources(Psi4, Psi6, alpn1, shiftx, shifty, shiftz, vxl, vyl, vzl, &
             u0l, u_x, u_y, u_z, &
             gxx, gyy, gzz, gxy, gxz, gyz, &
             gupxx, gupyy, gupzz, gupxy, gupxz, gupyz,&
             E_rad_nuel, F_radx_nuel, F_rady_nuel, F_radz_nuel, F_rad0_nuel, &
             tau_rad_nuel, S_rad_x_nuel, S_rad_y_nuel, S_rad_z_nuel, &
             P_radxx_nuel, P_radyy_nuel, P_radzz_nuel, P_radxy_nuel, P_radxz_nuel, P_radyz_nuel, &
             rhol, Sxl, Syl, Szl, Sxxl, Syyl, Szzl, Sxyl, Sxzl, Syzl, &
             rad_closure_scheme, Erad_atm_cut)
        call rad_sources(Psi4, Psi6, alpn1, shiftx, shifty, shiftz, vxl, vyl, vzl, &
             u0l, u_x, u_y, u_z, &
             gxx, gyy, gzz, gxy, gxz, gyz, &
             gupxx, gupyy, gupzz, gupxy, gupxz, gupyz,&
             E_rad_nuxl, F_radx_nuxl, F_rady_nuxl, F_radz_nuxl, F_rad0_nuxl, &
             tau_rad_nuxl, S_rad_x_nuxl, S_rad_y_nuxl, S_rad_z_nuxl, &
             P_radxx_nuxl, P_radyy_nuxl, P_radzz_nuxl, P_radxy_nuxl, P_radxz_nuxl, P_radyz_nuxl, &
             rhol, Sxl, Syl, Szl, Sxxl, Syyl, Szzl, Sxyl, Sxzl, Syzl, &
             rad_closure_scheme, Erad_atm_cut)
     end if
     if(tau_radl.lt.0.0) then
!        write(*,*) Inside recompute_vars, tau_rad and S_rad are, tau_radl, S_rad_xl
!        write(*,*) E_radl, u0l, F_rad0l, P_rad00l, F_rad_xl, P_rad0xl, Fasq are, E_radl, u0l, F_rad0l, P_rad00l, F_rad_xl, P_rad0xl
        tau_radl = alpn1*alpn1*Psi6*(Erad_atm_cut*u0l*u0l)
     end if
     if(tau_rad_nuel.lt.0.0) then
        tau_rad_nuel = alpn1*alpn1*Psi6*(Erad_atm_cut*u0l*u0l)
     end if
     if(tau_rad_nuxl.lt.0.0) then
        tau_rad_nuxl = alpn1*alpn1*Psi6*(Erad_atm_cut*u0l*u0l)
     end if
end subroutine recompute_vars
subroutine rad_sources(Psi4, Psi6, alpn1, &
     shiftx, shifty, shiftz, vxl, vyl, vzl, &
     u0l, u_x, u_y, u_z, &
     gxx, gyy, gzz, gxy, gxz, gyz, &
     gupxx, gupyy, gupzz, gupxy, gupxz, gupyz,&
     E_radl, F_radxl, F_radyl, F_radzl, F_rad0l, &
     tau_radl, S_rad_xl, S_rad_yl, S_rad_zl, &
     P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl, &
     rhol, Sxl, Syl, Szl, Sxxl, Syyl, Szzl, Sxyl, Sxzl, Syzl, &
     rad_closure_scheme, Erad_atm_cut)
  implicit none
  !! Input:
  real*8                    :: Psi4, Psi6, alpn1, gxx, gyy, gzz, gxy, gxz, gyz
  real*8                    :: gupxx, gupyy, gupzz, gupxy, gupxz, gupyz
  real*8                    :: shiftx, shifty, shiftz, vxl, vyl, vzl, u0l, u_x, u_y, u_z
  real*8                    :: E_radl, F_radxl, F_radyl, F_radzl
  !! Output:
  real*8                    :: tau_radl, S_rad_xl, S_rad_yl, S_rad_zl, F_rad0l
  real*8                    :: P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl
  real*8                    :: rhol, Sxl, Syl, Szl, Sxxl, Syyl, Szzl, Sxyl, Sxzl, Syzl
  !! Others:
  real*8                    :: shift_x, shift_y, shift_z, v_xl, v_yl, v_zl
  real*8                    :: beta2, udotbeta, g_00l, u_0l
  real*8                    :: uxl, uyl, uzl
  real*8                    :: F_rad_xl, F_rad_yl, F_rad_zl, F_rad_0l, P_radl
  real*8                    :: temp_rad, temp_rad1, Fasq, zeta_temp, zeta_cut, zeta, chi
  real*8                    :: P_rad0xl, P_rad0yl, P_rad0zl, P_rad00l
  real*8                    :: Erad_atm_cut
  integer*8                 :: rad_closure_scheme
  shift_x = Psi4 *(gxx*shiftx + gxy*shifty + gxz*shiftz)
  shift_y = Psi4 *(gxy*shiftx + gyy*shifty + gyz*shiftz)
  shift_z = Psi4 *(gxz*shiftx + gyz*shifty + gzz*shiftz)
  v_xl = Psi4 *(gxx*vxl + gxy*vyl + gxz*vzl)
  v_yl = Psi4 *(gxy*vxl + gyy*vyl + gyz*vzl)
  v_zl = Psi4 *(gxz*vxl + gyz*vyl + gzz*vzl)
  beta2 = shiftx*shift_x + shifty*shift_y + shiftz*shift_z
  udotbeta = u0l*(vxl*shift_x + vyl*shift_y + vzl*shift_z)
  g_00l =beta2-alpn1*alpn1
  u_0l = g_00l*u0l + udotbeta
  uxl = u0l*vxl
  uyl = u0l*vyl
  uzl = u0l*vzl
  if (E_radl .lt. Erad_atm_cut) then
     E_radl = Erad_atm_cut
     F_radxl = 0.0
     F_radyl = 0.0
     F_radzl = 0.0
  end if
  F_rad0l = - (F_radxl*u_x + F_radyl*u_y + F_radzl*u_z)/u_0l
  F_rad_xl = Psi4 * (gxx * F_radxl + gxy * F_radyl + gxz * F_radzl) + shift_x* F_rad0l
  F_rad_yl = Psi4 * (gxy * F_radxl + gyy * F_radyl + gyz * F_radzl) + shift_y* F_rad0l
  F_rad_zl = Psi4 * (gxz * F_radxl + gyz * F_radyl + gzz * F_radzl) + shift_z* F_rad0l
  F_rad_0l = - (F_rad_xl*uxl + F_rad_yl*uyl + F_rad_zl*uzl)/u0l
  if (rad_closure_scheme.eq.0) then
     P_radl = E_radl/3.d0
     temp_rad = alpn1*u0l
     temp_rad1 = temp_rad*temp_rad*(E_radl+P_radl) - P_radl + 2.d0*alpn1*u0l*F_rad0l
     tau_radl = alpn1*alpn1*Psi6*(E_radl*u0l*u0l+2.0*F_rad0l*u0l+P_radl*u0l*u0l)-Psi6*P_radl
     S_rad_xl = alpn1*Psi6*((E_radl+P_radl)*u0l*u_x + F_rad0l*u_x + F_rad_xl * u0l)
     S_rad_yl = alpn1*Psi6*((E_radl+P_radl)*u0l*u_y + F_rad0l*u_y + F_rad_yl * u0l)
     S_rad_zl = alpn1*Psi6*((E_radl+P_radl)*u0l*u_z + F_rad0l*u_z + F_rad_zl * u0l)
     rhol = rhol + temp_rad1
     Sxl = Sxl + temp_rad*( ( (E_radl + P_radl) * u0l + F_rad0l) * (shift_x + v_xl) + F_rad_xl)
     Syl = Syl + temp_rad*( ( (E_radl + P_radl) * u0l + F_rad0l) * (shift_y + v_yl) + F_rad_yl)
     Szl = Szl + temp_rad*( ( (E_radl + P_radl) * u0l + F_rad0l) * (shift_z + v_zl) + F_rad_zl)
     Sxxl = Sxxl + (E_radl+P_radl)*(u0l*(shift_x + v_xl))**2 + 2.0*F_rad_xl*u0l*(shift_x + v_xl) + Psi4 * P_radl * gxx
     Syyl = Syyl + (E_radl+P_radl)*(u0l*(shift_y + v_yl))**2 + 2.0*F_rad_yl*u0l*(shift_y + v_yl) + Psi4 * P_radl * gyy
     Szzl = Szzl + (E_radl+P_radl)*(u0l*(shift_z + v_zl))**2 + 2.0*F_rad_zl*u0l*(shift_z + v_zl) + Psi4 * P_radl * gzz
     Sxyl = Sxyl + (E_radl+P_radl)*u0l**2*(shift_x + v_xl)*(shift_y + v_yl) + u0l*(F_rad_xl*(shift_y + v_yl)+F_rad_yl*(shift_x + v_&
  &xl)) + Psi4 * P_radl * gxy
     Sxzl = Sxzl + (E_radl+P_radl)*u0l**2*(shift_x + v_xl)*(shift_z + v_zl) + u0l*(F_rad_xl*(shift_z + v_zl)+F_rad_zl*(shift_x + v_&
  &xl)) + Psi4 * P_radl * gxz
     Syzl = Syzl + (E_radl+P_radl)*u0l**2*(shift_y + v_yl)*(shift_z + v_zl) + u0l*(F_rad_yl*(shift_z + v_zl)+F_rad_zl*(shift_y + v_&
  &yl)) + Psi4 * P_radl * gyz
  else
     Fasq = F_rad_0l*F_rad0l + F_rad_xl*F_radxl +  F_rad_yl*F_radyl +  F_rad_zl*F_radzl
     zeta_temp = sqrt(abs(F_rad_0l*F_rad0l + F_rad_xl*F_radxl +  F_rad_yl*F_radyl +  F_rad_zl*F_radzl)/E_radl**2)
     zeta_cut = Erad_atm_cut*1.5
     if (E_radl.le.zeta_cut) then
        zeta = 1.0
     else
        zeta = zeta_temp
     end if
     if (zeta .gt. 1.0) then
        zeta = 1.0
     end if
     chi = 1/3.0 + zeta**2*(6.0-2.0*zeta+6.0*zeta**2)/15.0
     if (E_radl .lt. Erad_atm_cut) then
        P_radxxl = 0.0
        P_radyyl = 0.0
        P_radzzl = 0.0
        P_radxyl = 0.0
        P_radxzl = 0.0
        P_radyzl = 0.0
     else
        if (Fasq .le. 0) then
           P_radxxl = E_radl*(gupxx/Psi4 - shiftx*shiftx/alpn1**2 + uxl**2)/2.0*(1.0-chi)
           P_radyyl = E_radl*(gupyy/Psi4 - shifty*shifty/alpn1**2 + uyl**2)/2.0*(1.0-chi)
           P_radzzl = E_radl*(gupzz/Psi4 - shiftz*shiftz/alpn1**2 + uzl**2)/2.0*(1.0-chi)
           P_radxyl = E_radl*(gupxy/Psi4 - shiftx*shifty/alpn1**2 + uxl*uyl)/2.0*(1.0-chi)
           P_radxzl = E_radl*(gupxz/Psi4 - shiftx*shiftz/alpn1**2 + uxl*uzl)/2.0*(1.0-chi)
           P_radyzl = E_radl*(gupyz/Psi4 - shifty*shiftz/alpn1**2 + uyl*uzl)/2.0*(1.0-chi)
        else
           P_radxxl = E_radl*((F_radxl**2/Fasq)*(3.0*chi -1.0)/2.0 + (gupxx/Psi4 - shiftx*shiftx/alpn1**2 + uxl**2)/2.0*(1.0-chi))
           P_radyyl = E_radl*((F_radyl**2/Fasq)*(3.0*chi -1.0)/2.0 + (gupyy/Psi4 - shifty*shifty/alpn1**2 + uyl**2)/2.0*(1.0-chi))
           P_radzzl = E_radl*((F_radzl**2/Fasq)*(3.0*chi -1.0)/2.0 + (gupzz/Psi4 - shiftz*shiftz/alpn1**2 + uzl**2)/2.0*(1.0-chi))
           P_radxyl = E_radl*((F_radxl*F_radyl/Fasq)*(3.0*chi -1.0)/2.0 + (gupxy/Psi4 - shiftx*shifty/alpn1**2 + uxl*uyl)/2.0*(1.0-&
  &chi))
           P_radxzl = E_radl*((F_radxl*F_radzl/Fasq)*(3.0*chi -1.0)/2.0 + (gupxz/Psi4 - shiftx*shiftz/alpn1**2 + uxl*uzl)/2.0*(1.0-&
  &chi))
           P_radyzl = E_radl*((F_radyl*F_radzl/Fasq)*(3.0*chi -1.0)/2.0 + (gupyz/Psi4 - shifty*shiftz/alpn1**2 ++ uyl*uzl)/2.0*(1.0&
  &-chi))
        end if
     end if
     P_rad0xl = - (P_radxxl * u_x + P_radxyl * u_y + P_radxzl * u_z)/u_0l
     P_rad0yl = - (P_radxyl * u_x + P_radyyl * u_y + P_radyzl * u_z)/u_0l
     P_rad0zl = - (P_radxzl * u_x + P_radyzl * u_y + P_radzzl * u_z)/u_0l
     P_rad00l = - (P_rad0xl * u_x + P_rad0yl * u_y + P_rad0zl * u_z)/u_0l
     tau_radl = alpn1*alpn1*Psi6*(E_radl*u0l*u0l+2.0*F_rad0l*u0l+P_rad00l)
     S_rad_xl = alpn1*Psi6*(E_radl*u0l*u_x + F_rad0l*u_x + F_rad_xl * u0l + P_rad00l*shift_x + Psi4*(P_rad0xl*gxx + P_rad0yl*gxy + &
  &P_rad0zl*gxz))
     S_rad_yl = alpn1*Psi6*(E_radl*u0l*u_y + F_rad0l*u_y + F_rad_yl * u0l + P_rad00l*shift_y + Psi4*(P_rad0xl*gxy + P_rad0yl*gyy + &
  &P_rad0zl*gyz))
     S_rad_zl = alpn1*Psi6*(E_radl*u0l*u_z + F_rad0l*u_z + F_rad_zl * u0l + P_rad00l*shift_z + Psi4*(P_rad0xl*gxz + P_rad0yl*gyz + &
  &P_rad0zl*gzz))
     rhol = rhol + alpn1**2.0*(E_radl*u0l**2 + 2.0 * F_rad0l * u0l + P_rad00l)
     Sxl = Sxl + alpn1*(u0l*E_radl*u_x + F_rad0l*u_x + u0l*F_rad_xl +&
          P_rad00l*shift_x + Psi4*(P_rad0xl*gxx + P_rad0yl*gxy + P_rad0zl*gxz))
     Syl = Syl + alpn1*(u0l*E_radl*u_y + F_rad0l*u_y + u0l*F_rad_yl +&
          P_rad00l*shift_y + Psi4*(P_rad0xl*gxy + P_rad0yl*gyy + P_rad0zl*gyz))
     Szl = Szl + alpn1*(u0l*E_radl*u_z + F_rad0l*u_z + u0l*F_rad_zl +&
          P_rad00l*shift_z + Psi4*(P_rad0xl*gxz + P_rad0yl*gyz + P_rad0zl*gzz))
     Sxxl = Sxxl + E_radl*u_x*u_x + 2.0*F_rad_xl*u_x +&
          shift_x**2.0*P_rad00l + shift_x*2.0*(gxx*P_rad0xl+gxy*P_rad0yl+gxz*P_rad0zl) +&
          Psi4**2.0*( gxx**2.0*P_radxxl + gxy**2.0*P_radyyl + gxz**2.0*P_radzzl +&
          2.0*(gxx*gxy*P_radxyl+gxx*gxz*P_radxzl+gxy*gxz*P_radyzl) )
     Syyl = Syyl + E_radl*u_y*u_y + 2.0*F_rad_yl*u_y +&
          shift_y**2.0*P_rad00l + shift_y*2.0*(gxy*P_rad0xl+gyy*P_rad0yl+gyz*P_rad0zl) +&
          Psi4**2.0*( gxy**2.0*P_radxxl + gyy**2.0*P_radyyl + gyz**2.0*P_radzzl +&
          2.0*(gxy*gyy*P_radxyl+gxy*gyz*P_radxzl+gyy*gyz*P_radyzl) )
     Szzl = Szzl + E_radl*u_z*u_z + 2.0*F_rad_zl*u_z +&
          shift_z**2.0*P_rad00l + shift_z*2.0*(gxz*P_rad0xl+gyz*P_rad0yl+gzz*P_rad0zl) +&
          Psi4**2.0*( gxz**2.0*P_radxxl + gyz**2.0*P_radyyl + gzz**2.0*P_radzzl +&
          2.0*(gxz*gyz*P_radxyl+gxz*gzz*P_radxzl+gyz*gzz*P_radyzl) )
     Sxyl = Sxyl + E_radl*u_x*u_y + F_rad_xl*u_y + F_rad_yl*u_x +&
          shift_x*shift_y*P_rad00l +&
          Psi4*(shift_x*(gxy*P_rad0xl + gyy*P_rad0yl + gyz*P_rad0zl) + (shift_y*(gxx*P_rad0xl + gxy*P_rad0yl + gxz*P_rad0zl)) )+&
          Psi4**2.0*(gxx*gxy*P_radxxl + gxy*gyy*P_radyyl + gxz*gyz*P_radzzl +&
          (gxx*gyy + gxy*gxy)*P_radxyl + (gxx*gyz + gxz*gxy)*P_radxzl + (gxy*gyz + gxz*gyy)*P_radyzl)
     Sxzl = Sxzl + E_radl*u_x*u_z + F_rad_xl*u_z + F_rad_zl*u_x +&
          shift_x*shift_z*P_rad00l +&
          Psi4*(shift_x*(gxz*P_rad0xl + gyz*P_rad0yl + gzz*P_rad0zl) + (shift_z*(gxx*P_rad0xl + gxy*P_rad0yl + gxz*P_rad0zl)) )+&
          Psi4**2.0*(gxx*gxz*P_radxxl + gxy*gyz*P_radyyl + gxz*gzz*P_radzzl +&
          (gxx*gyz + gxy*gxz)*P_radxyl + (gxx*gzz + gxz*gxz)*P_radxzl + (gxy*gzz + gxz*gyz)*P_radyzl)
     Syzl = Syzl + E_radl*u_y*u_z + F_rad_yl*u_z + F_rad_zl*u_y +&
          shift_y*shift_z*P_rad00l +&
          Psi4*(shift_y*(gxz*P_rad0xl + gyz*P_rad0yl + gzz*P_rad0zl) + (shift_z*(gxy*P_rad0xl + gyy*P_rad0yl + gyz*P_rad0zl)) )+&
          Psi4**2.0*(gxy*gxz*P_radxxl + gyy*gyz*P_radyyl + gyz*gzz*P_radzzl +&
          (gxy*gyz + gxz*gyy)*P_radxyl + (gxy*gzz + gyz*gxz)*P_radxzl + (gyy*gzz + gyz*gyz)*P_radyzl)
  end if
end subroutine rad_sources