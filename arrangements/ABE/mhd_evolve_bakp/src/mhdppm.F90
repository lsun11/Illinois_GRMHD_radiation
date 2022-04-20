!-----------------------------------------------------------------!
! Calculate a_{j-1/2}, set a_R, a_L equal to this.
!-----------------------------------------------------------------!
!
subroutine mhdppm_find_face_vals(ext,X,Y,Z,rho,delta_rho,rhor,rhol,P,delta_P,Pr,Pl, &
     vx,delta_vx,vxr,vxl,vy,delta_vy,vyr,vyl,vz,delta_vz,vzr,vzl, &
     Bx,delta_Bx,Bxr,Bxl,By,delta_By,Byr,Byl,Bz,delta_Bz,Bzr,Bzl, &
     m,Symmetry,Sym_Bz)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z,rho,delta_rho,rhor,rhol
 real*8, dimension(ext(1),ext(2),ext(3))        :: P,delta_P,Pr,Pl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vx,delta_vx,vxr,vxl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vy,delta_vy,vyr,vyl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vz,delta_vz,vzr,vzl
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bx,delta_Bx,Bxr,Bxl
 real*8, dimension(ext(1),ext(2),ext(3))        :: By,delta_By,Byr,Byl
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bz,delta_Bz,Bzr,Bzl
 integer                                        :: m,Symmetry
 real*8  :: Sym_Bz, mSym_Bz
!
 integer                                        :: j,jmin,jmax
 real*8                     :: SYM, ANTI, ZERO, HALF, SIXTH, EIGHTH
 parameter(ZERO = 0.D0, HALF = 1.0D0/2.D0, SIXTH = 1.D0/6.D0, EIGHTH = 1.D0/8.D0, SYM = 1.D0, ANTI = -1.D0)
 integer                            :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
 parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

 mSym_Bz = -Sym_Bz
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
       Bxl(jmin,:,:) = SYM*Bx(jmin,:,:) + HALF*(Bx(jmin,:,:)-SYM*Bx(jmin,:,:)) &
            + EIGHTH*(ANTI*delta_Bx(jmin,:,:) - delta_Bx(jmin,:,:))
    else
       vxl(jmin,:,:) = ANTI*vx(jmin,:,:) + HALF*(vx(jmin,:,:)-ANTI*vx(jmin,:,:)) &
            + EIGHTH*(SYM*delta_vx(jmin,:,:) - delta_vx(jmin,:,:))
       Bxl(jmin,:,:) = ANTI*Bx(jmin,:,:) + HALF*(Bx(jmin,:,:)-ANTI*Bx(jmin,:,:)) &
            + EIGHTH*(SYM*delta_Bx(jmin,:,:) - delta_Bx(jmin,:,:))
    end if
    if(Symmetry==AXISYM) then
       vyl(jmin,:,:) = ANTI*vy(jmin,:,:) + HALF*(vy(jmin,:,:)-ANTI*vy(jmin,:,:)) &
            + EIGHTH*(SYM*delta_vy(jmin,:,:) - delta_vy(jmin,:,:))
       Byl(jmin,:,:) = ANTI*By(jmin,:,:) + HALF*(By(jmin,:,:)-ANTI*By(jmin,:,:)) &
            + EIGHTH*(SYM*delta_By(jmin,:,:) - delta_By(jmin,:,:))
    else
       vyl(jmin,:,:) = SYM*vy(jmin,:,:) + HALF*(vy(jmin,:,:)-SYM*vy(jmin,:,:)) &
            + EIGHTH*(ANTI*delta_vy(jmin,:,:) - delta_vy(jmin,:,:))
       Byl(jmin,:,:) = SYM*By(jmin,:,:) + HALF*(By(jmin,:,:)-SYM*By(jmin,:,:)) &
            + EIGHTH*(ANTI*delta_By(jmin,:,:) - delta_By(jmin,:,:))
    end if
    vzl(jmin,:,:) = SYM*vz(jmin,:,:) + HALF*(vz(jmin,:,:)-SYM*vz(jmin,:,:)) &
            + EIGHTH*(ANTI*delta_vz(jmin,:,:) - delta_vz(jmin,:,:))
    Bzl(jmin,:,:) = SYM*Bz(jmin,:,:) + HALF*(Bz(jmin,:,:)-SYM*Bz(jmin,:,:)) &
            + EIGHTH*(ANTI*delta_Bz(jmin,:,:) - delta_Bz(jmin,:,:))
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
       Bxl(j,:,:) = Bx(j-1,:,:) + HALF*(Bx(j,:,:)-Bx(j-1,:,:)) &
            + EIGHTH*(delta_Bx(j-1,:,:) - delta_Bx(j,:,:))
       Byl(j,:,:) = By(j-1,:,:) + HALF*(By(j,:,:)-By(j-1,:,:)) &
            + EIGHTH*(delta_By(j-1,:,:) - delta_By(j,:,:))
       Bzl(j,:,:) = Bz(j-1,:,:) + HALF*(Bz(j,:,:)-Bz(j-1,:,:)) &
            + EIGHTH*(delta_Bz(j-1,:,:) - delta_Bz(j,:,:))
    end do
!    write(*,*) "inside mhdppm_ffv0",vxl(15,2,2),Bxl(15,2,2),Bx(15,2,2),delta_Bx(15,2,2),Bx(1,2,2),delta_Bx(1,2,2)
 
!    write(*,*) "inside mhdppm_ffv0", vxl(jmin+1,2,2),delta_vx(jmin+1,2,2),delta_vx(jmin+15,2,2),delta_vx(jmin,2,2),vx(jmin,2,2),jmin
    do j=jmin, jmax-1
       rhor(j,:,:) = rhol(j+1,:,:)
       Pr(j,:,:)   = Pl(j+1,:,:)
       vxr(j,:,:)  = vxl(j+1,:,:)
       vyr(j,:,:)  = vyl(j+1,:,:)
       vzr(j,:,:)  = vzl(j+1,:,:)
       Bxr(j,:,:)  = Bxl(j+1,:,:)
       Byr(j,:,:)  = Byl(j+1,:,:)
       Bzr(j,:,:)  = Bzl(j+1,:,:)
    end do
    rhor(jmax,:,:) = rho(jmax,:,:)
    Pr(jmax,:,:)   = P(jmax,:,:)
    vxr(jmax,:,:)  = vx(jmax,:,:)
    vyr(jmax,:,:)  = vy(jmax,:,:)
    vzr(jmax,:,:)  = vz(jmax,:,:)
    Bxr(jmax,:,:)  = Bx(jmax,:,:)
    Byr(jmax,:,:)  = By(jmax,:,:)
    Bzr(jmax,:,:)  = Bz(jmax,:,:)

 else if(m==2) then
    rhol(:,jmin,:) = SYM*rho(:,jmin,:) + HALF*(rho(:,jmin,:)-SYM*rho(:,jmin,:)) &
            + EIGHTH*(ANTI*delta_rho(:,jmin,:) - delta_rho(:,jmin,:))
    Pl(:,jmin,:) = SYM*P(:,jmin,:) + HALF*(P(:,jmin,:)-SYM*P(:,jmin,:)) &
            + EIGHTH*(ANTI*delta_P(:,jmin,:) - delta_P(:,jmin,:))
    vxl(:,jmin,:) = SYM*vx(:,jmin,:) + HALF*(vx(:,jmin,:)-SYM*vx(:,jmin,:)) &
         + EIGHTH*(ANTI*delta_vx(:,jmin,:) - delta_vx(:,jmin,:))
    Bxl(:,jmin,:) = SYM*Bx(:,jmin,:) + HALF*(Bx(:,jmin,:)-SYM*Bx(:,jmin,:)) &
         + EIGHTH*(ANTI*delta_Bx(:,jmin,:) - delta_Bx(:,jmin,:))
    if(Symmetry==NO_SYMM.or.Symmetry==EQUATORIAL) then
       vyl(:,jmin,:) = SYM*vy(:,jmin,:) + HALF*(vy(:,jmin,:)-SYM*vy(:,jmin,:)) &
            + EIGHTH*(ANTI*delta_vy(:,jmin,:) - delta_vy(:,jmin,:))
       Byl(:,jmin,:) = SYM*By(:,jmin,:) + HALF*(By(:,jmin,:)-SYM*By(:,jmin,:)) &
            + EIGHTH*(ANTI*delta_By(:,jmin,:) - delta_By(:,jmin,:))
    else
       vyl(:,jmin,:) = ANTI*vy(:,jmin,:) + HALF*(vy(:,jmin,:)-ANTI*vy(:,jmin,:)) &
            + EIGHTH*(SYM*delta_vy(:,jmin,:) - delta_vy(:,jmin,:))
       Byl(:,jmin,:) = ANTI*By(:,jmin,:) + HALF*(By(:,jmin,:)-ANTI*By(:,jmin,:)) &
            + EIGHTH*(SYM*delta_By(:,jmin,:) - delta_By(:,jmin,:))
    end if
    vzl(:,jmin,:) = SYM*vz(:,jmin,:) + HALF*(vz(:,jmin,:)-SYM*vz(:,jmin,:)) &
            + EIGHTH*(ANTI*delta_vz(:,jmin,:) - delta_vz(:,jmin,:))
    Bzl(:,jmin,:) = SYM*Bz(:,jmin,:) + HALF*(Bz(:,jmin,:)-SYM*Bz(:,jmin,:)) &
            + EIGHTH*(ANTI*delta_Bz(:,jmin,:) - delta_Bz(:,jmin,:))
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
       Bxl(:,j,:) = Bx(:,j-1,:) + HALF*(Bx(:,j,:)-Bx(:,j-1,:)) &
            + EIGHTH*(delta_Bx(:,j-1,:) - delta_Bx(:,j,:))
       Byl(:,j,:) = By(:,j-1,:) + HALF*(By(:,j,:)-By(:,j-1,:)) &
            + EIGHTH*(delta_By(:,j-1,:) - delta_By(:,j,:))
       Bzl(:,j,:) = Bz(:,j-1,:) + HALF*(Bz(:,j,:)-Bz(:,j-1,:)) &
            + EIGHTH*(delta_Bz(:,j-1,:) - delta_Bz(:,j,:))
    end do
    
    do j=jmin, jmax-1
       rhor(:,j,:) = rhol(:,j+1,:)
       Pr(:,j,:)   = Pl(:,j+1,:)
       vxr(:,j,:)  = vxl(:,j+1,:)
       vyr(:,j,:)  = vyl(:,j+1,:)
       vzr(:,j,:)  = vzl(:,j+1,:)
       Bxr(:,j,:)  = Bxl(:,j+1,:)
       Byr(:,j,:)  = Byl(:,j+1,:)
       Bzr(:,j,:)  = Bzl(:,j+1,:)
    end do
    rhor(:,jmax,:) = rho(:,jmax,:)
    Pr(:,jmax,:)   = P(:,jmax,:)
    vxr(:,jmax,:)  = vx(:,jmax,:)
    vyr(:,jmax,:)  = vy(:,jmax,:)
    vzr(:,jmax,:)  = vz(:,jmax,:)
    Bxr(:,jmax,:)  = Bx(:,jmax,:)
    Byr(:,jmax,:)  = By(:,jmax,:)
    Bzr(:,jmax,:)  = Bz(:,jmax,:)

 else
    rhol(:,:,jmin) = SYM*rho(:,:,jmin) + HALF*(rho(:,:,jmin)-SYM*rho(:,:,jmin)) &
            + EIGHTH*(ANTI*delta_rho(:,:,jmin) - delta_rho(:,:,jmin))
    Pl(:,:,jmin) = SYM*P(:,:,jmin) + HALF*(P(:,:,jmin)-SYM*P(:,:,jmin)) &
            + EIGHTH*(ANTI*delta_P(:,:,jmin) - delta_P(:,:,jmin))
    vxl(:,:,jmin) = SYM*vx(:,:,jmin) + HALF*(vx(:,:,jmin)-SYM*vx(:,:,jmin)) &
         + EIGHTH*(ANTI*delta_vx(:,:,jmin) - delta_vx(:,:,jmin))
    vyl(:,:,jmin) = SYM*vy(:,:,jmin) + HALF*(vy(:,:,jmin)-SYM*vy(:,:,jmin)) &
         + EIGHTH*(ANTI*delta_vy(:,:,jmin) - delta_vy(:,:,jmin))
    if(Symmetry==NO_SYMM) then
       vzl(:,:,jmin) = SYM*vz(:,:,jmin) + HALF*(vz(:,:,jmin)-SYM*vz(:,:,jmin)) &
            + EIGHTH*(ANTI*delta_vz(:,:,jmin) - delta_vz(:,:,jmin))
       Bzl(:,:,jmin) = SYM*Bz(:,:,jmin) + HALF*(Bz(:,:,jmin)-SYM*Bz(:,:,jmin)) &
            + EIGHTH*(ANTI*delta_Bz(:,:,jmin) - delta_Bz(:,:,jmin))
       Bxl(:,:,jmin) = SYM*Bx(:,:,jmin) + HALF*(Bx(:,:,jmin)-SYM*Bx(:,:,jmin)) &
            + EIGHTH*(ANTI*delta_Bx(:,:,jmin) - delta_Bx(:,:,jmin))
       Byl(:,:,jmin) = SYM*By(:,:,jmin) + HALF*(By(:,:,jmin)-SYM*By(:,:,jmin)) &
            + EIGHTH*(ANTI*delta_By(:,:,jmin) - delta_By(:,:,jmin))
    else
       vzl(:,:,jmin) = ANTI*vz(:,:,jmin) + HALF*(vz(:,:,jmin)-ANTI*vz(:,:,jmin)) &
            + EIGHTH*(SYM*delta_vz(:,:,jmin) - delta_vz(:,:,jmin))
       Bzl(:,:,jmin) = Sym_Bz*SYM*Bz(:,:,jmin) + HALF*(Bz(:,:,jmin)-Sym_Bz*SYM*Bz(:,:,jmin)) &
            + EIGHTH*(Sym_Bz*ANTI*delta_Bz(:,:,jmin) - delta_Bz(:,:,jmin))
       Bxl(:,:,jmin) = mSym_Bz*SYM*Bx(:,:,jmin) + HALF*(Bx(:,:,jmin)-mSym_Bz*SYM*Bx(:,:,jmin)) &
            + EIGHTH*(mSym_Bz*ANTI*delta_Bx(:,:,jmin) - delta_Bx(:,:,jmin))
       Byl(:,:,jmin) = mSym_Bz*SYM*By(:,:,jmin) + HALF*(By(:,:,jmin)-mSym_Bz*SYM*By(:,:,jmin)) &
            + EIGHTH*(mSym_Bz*ANTI*delta_By(:,:,jmin) - delta_By(:,:,jmin))
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
       Bxl(:,:,j) = Bx(:,:,j-1) + HALF*(Bx(:,:,j)-Bx(:,:,j-1)) &
            + EIGHTH*(delta_Bx(:,:,j-1) - delta_Bx(:,:,j))
       Byl(:,:,j) = By(:,:,j-1) + HALF*(By(:,:,j)-By(:,:,j-1)) &
            + EIGHTH*(delta_By(:,:,j-1) - delta_By(:,:,j))
       Bzl(:,:,j) = Bz(:,:,j-1) + HALF*(Bz(:,:,j)-Bz(:,:,j-1)) &
            + EIGHTH*(delta_Bz(:,:,j-1) - delta_Bz(:,:,j))
    end do
    
    do j=jmin, jmax-1
       rhor(:,:,j) = rhol(:,:,j+1)
       Pr(:,:,j)   = Pl(:,:,j+1)
       vxr(:,:,j)  = vxl(:,:,j+1)
       vyr(:,:,j)  = vyl(:,:,j+1)
       vzr(:,:,j)  = vzl(:,:,j+1)
       Bxr(:,:,j)  = Bxl(:,:,j+1)
       Byr(:,:,j)  = Byl(:,:,j+1)
       Bzr(:,:,j)  = Bzl(:,:,j+1)
    end do
    rhor(:,:,jmax) = rho(:,:,jmax)
    Pr(:,:,jmax)   = P(:,:,jmax)
    vxr(:,:,jmax)  = vx(:,:,jmax)
    vyr(:,:,jmax)  = vy(:,:,jmax)
    vzr(:,:,jmax)  = vz(:,:,jmax)
    Bxr(:,:,jmax)  = Bx(:,:,jmax)
    Byr(:,:,jmax)  = By(:,:,jmax)
    Bzr(:,:,jmax)  = Bz(:,:,jmax)

 end if
! write(*,*) "inside mhdppm_ffv1",vx(15,2,2),vy(15,2,2),vz(15,2,2),vxl(15,2,2),vxr(15,2,2),vyl(15,2,2),vyr(15,2,2),vzl(15,2,2),vzr(15,2,2),Bxl(15,2,2),Byl(15,2,2),Bzl(15,2,2),delta_Bx(15,2,2),delta_By(15,2,2),delta_Bz(15,2,2),rhol(15,2,2),rhol(2,3,2),rhor(15,2,2),rhor(2,3,2)
! write(*,*) "inside mhdppm_ffv1",vx(2,3,2),vy(2,3,2),vz(2,3,2),vxl(2,3,2),vxr(2,3,2),vyl(2,3,2),vyr(2,3,2),vzl(2,3,2),vzr(2,3,2),Bxl(2,3,2),Byl(2,3,2),Bzl(2,3,2),delta_Bx(2,3,2),delta_By(2,3,2),delta_Bz(2,3,2),rhol(2,3,2),rhol(2,3,2),rhor(2,3,2),rhor(2,3,2)

! write(*,*) "inside mhdppm_ffv1",rhor(15,2,2),Pr(15,2,2),vxr(15,2,2),Bxr(15,2,2),vyr(15,2,2),Byr(15,2,2),vzr(15,2,2),Bzr(15,2,2),delta_Bx(15,2,2),Bx(15,2,2)
! write(*,*) "inside mhdppm_ffv2",rhol(15,2,2),Pl(15,2,2),vxl(15,2,2),Bxl(15,2,2),vyl(15,2,2),Byl(15,2,2),vzl(15,2,2),Bzl(15,2,2),delta_Bx(15,2,2),Bx(15,2,2)
 

end subroutine mhdppm_find_face_vals
!-----------------------------------------------------------------!
! ensure monotonicity of the interpolating polynomial             !
!-----------------------------------------------------------------!
!
subroutine mhdppm_monotonize(ext, X,Y,Z, &
     Bx,Bxr,Bxl, By,Byr,Byl, Bz,Bzr,Bzl, &
     m, Symmetry, Sym_Bz)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bx,Bxr,Bxl
 real*8, dimension(ext(1),ext(2),ext(3))        :: By,Byr,Byl
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bz,Bzr,Bzl
 integer, intent(in)                            :: m,Symmetry
 integer                                        :: i,j,k
 integer                                        :: imin,imax,jmin,jmax,kmin,kmax
 real*8                                         :: SYM, ANTI, Sym_Bz, mSym_Bz
 parameter(SYM = 1.D0, ANTI = -1.D0)
 integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
 parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
 mSym_Bz = -Sym_Bz
 
 kmin = lbound(Bx,3)
 kmax = ubound(Bx,3)

 if(Symmetry == OCTANT) then 
    if(X(1,1,1).lt.0.D0) then
       kmin = lbound(Bx,1)+1
    else if(Y(1,1,1).lt.0.D0) then
       kmin = lbound(Bx,2)+1
    else if(Z(1,1,1).lt.0.D0) then
       kmin = lbound(Bx,3)+1     
    end if
 else if(Symmetry == AXISYM) then 
    if(X(1,1,1).lt.0.D0) then
       kmin = lbound(Bx,1)+1
    else if(Z(1,1,1).lt.0.D0) then
       kmin = lbound(Bx,3)+1     
    end if
 else if(Symmetry == EQUATORIAL) then 
    if(Z(1,1,1).lt.0.D0) then
       kmin = lbound(Bx,3)+1     
    end if
 end if

! write(*,*) "inside mhdppm_monot",Bxl(15,2,2),Bxr(15,2,2),Byl(15,2,2),Byr(15,2,2),Bzl(15,2,2),Bzr(15,2,2), &
!      Bx(15,2,2),By(15,2,2),Bz(15,2,2),Bx(2,3,2),By(2,3,2),Bz(2,3,2)
 if(m==1) then
    call monotonize_nock(ext,X,Y,Z,Bz,Bzr,Bzl,m,SYM,Symmetry)
    if(Symmetry==OCTANT) then
       call monotonize_nock(ext,X,Y,Z,Bx,Bxr,Bxl,m,ANTI,Symmetry)
       call monotonize_nock(ext,X,Y,Z,By,Byr,Byl,m,SYM,Symmetry)
    else if(Symmetry==AXISYM) then
       call monotonize_nock(ext,X,Y,Z,Bx,Bxr,Bxl,m,ANTI,Symmetry)
       call monotonize_nock(ext,X,Y,Z,By,Byr,Byl,m,ANTI,Symmetry)
    else
       call monotonize_nock(ext,X,Y,Z,Bx,Bxr,Bxl,m,SYM,Symmetry)
       call monotonize_nock(ext,X,Y,Z,By,Byr,Byl,m,SYM,Symmetry)
    end if
 else if(m==2) then
    if(Symmetry.ne.AXISYM) then
       call monotonize_nock(ext,X,Y,Z,Bx,Bxr,Bxl,m,SYM,Symmetry)
       call monotonize_nock(ext,X,Y,Z,Bz,Bzr,Bzl,m,SYM,Symmetry)
       if(Symmetry==OCTANT) then
          call monotonize_nock(ext,X,Y,Z,By,Byr,Byl,m,ANTI,Symmetry)
       else
          call monotonize_nock(ext,X,Y,Z,By,Byr,Byl,m,SYM,Symmetry)
       end if
    end if
 else if(m==3) then
    if((Symmetry.eq.NO_SYMM) .or. (Z(1,1,kmin+1) .lt. 0.d0)) then
       call monotonize_nock(ext,X,Y,Z,Bx,Bxr,Bxl,m,SYM,Symmetry)
       call monotonize_nock(ext,X,Y,Z,By,Byr,Byl,m,SYM,Symmetry)
       call monotonize_nock(ext,X,Y,Z,Bz,Bzr,Bzl,m,SYM,Symmetry)
    else
       call monotonize_nock(ext,X,Y,Z,Bx,Bxr,Bxl,m,mSym_Bz,Symmetry)
       call monotonize_nock(ext,X,Y,Z,By,Byr,Byl,m,mSym_Bz,Symmetry)
       call monotonize_nock(ext,X,Y,Z,Bz,Bzr,Bzl,m,Sym_Bz,Symmetry)
    end if
 else
    write(*,*) "mhdppm_monotonize: m not in proper range",m
 end if

! write(*,*) "inside mhdppm_monot",Bxl(15,2,2),Bxr(15,2,2),Byl(15,2,2),Byr(15,2,2),Bzl(15,2,2),Bzr(15,2,2), &
!      Bx(15,2,2),By(15,2,2),Bz(15,2,2),Bx(2,3,2),By(2,3,2),Bz(2,3,2)
end subroutine mhdppm_monotonize
!--------------------------------------------------------------------!
! ensure monotonicity of the interpolating polynomial (no Clark fix) !
!--------------------------------------------------------------------!
!
subroutine monotonize_nock(ext,X,Y,Z,a,ar,al,m,sym,Symmetry)
  implicit none
 integer, dimension(3)                     :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8, dimension(ext(1),ext(2),ext(3))   :: a,ar,al
 integer, intent(in)                       :: m,Symmetry
 real*8                                    :: sym
 integer                                   :: i,j,k
 integer                                   :: imin,imax,jmin,jmax,kmin,kmax
 real*8, parameter                         :: ZERO = 0.d0, HALF = 0.5d0, SIX = 6.d0
 real*8, parameter                         :: THREE = 3.d0, TWO = 2.d0
 real*8                                    :: q1, q2
 integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
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


 do k=kmin,kmax
    do j=jmin,jmax
       do i=imin,imax

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

          if(m==1 .and. i==imin .and. sym.lt.ZERO) al(i,j,k) = ZERO
          if(m==2 .and. j==jmin .and. sym.lt.ZERO) al(i,j,k) = ZERO
          if(m==3 .and. k==kmin .and. sym.lt.ZERO) al(i,j,k) = ZERO
             
       end do
    end do
 end do

end subroutine monotonize_nock

!-----------------------------------------------------------------!
!  redefine aL, aR                                                !
!-----------------------------------------------------------------!
subroutine mhdppm_shift(ext,X,Y,Z,rhor,rhol,Pr,Pl,vxr,vxl,vyr,vyl,vzr,vzl, &
     Bxr,Bxl,Byr,Byl,Bzr,Bzl,m,Symmetry)
  implicit none
   integer, dimension(3)                        :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8, dimension(ext(1),ext(2),ext(3))        :: rhor,rhol,Pr,Pl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vxr,vxl,vyr,vyl,vzr,vzl
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bxr,Bxl,Byr,Byl,Bzr,Bzl
 real*8, dimension(ext(1),ext(2),ext(3))        :: rhot,Pt,vxt,vyt,vzt
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bxt,Byt,Bzt
 integer                                        :: m, Symmetry
 integer                                        :: k, kmin, kmax
 integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
 parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

 kmin = lbound(rhor,m)
 kmax = ubound(rhor,m)

! write(*,*) "inside mhdppm_shift0",vzl(2,3,2),vzr(2,3,2),Bzl(2,3,2),Bzr(2,3,2),rhol(2,3,2),rhor(2,3,2)
! write(*,*) "inside mhdppm_shift1",vzl(15,2,2),vzr(15,2,2),Bzl(15,2,2),Bzr(15,2,2),rhol(15,2,2),rhor(15,2,2)

 rhot = rhor
 Pt   = Pr
 vxt  = vxr
 vyt  = vyr
 vzt  = vzr
 Bxt  = Bxr
 Byt  = Byr
 Bzt  = Bzr

 rhor = rhol
 Pr   = Pl
 vxr  = vxl
 vyr  = vyl
 vzr  = vzl
 Bxr  = Bxl
 Byr  = Byl
 Bzr  = Bzl

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
 else if(Symmetry == EQUATORIAL) then 
    if(m==3 .and. Z(1,1,1).lt.0.D0) then
       kmin = lbound(rhor,3)+1     
    end if
 end if



 if(m==1) then
    
    do k=kmin+1, kmax
       rhol(k,:,:) = rhot(k-1,:,:)
       Pl(k,:,:)   = Pt(k-1,:,:)
       vxl(k,:,:)  = vxt(k-1,:,:)
       vyl(k,:,:)  = vyt(k-1,:,:)
       vzl(k,:,:)  = vzt(k-1,:,:)
       Bxl(k,:,:)  = Bxt(k-1,:,:)
       Byl(k,:,:)  = Byt(k-1,:,:)
       Bzl(k,:,:)  = Bzt(k-1,:,:)
    end do

 else if (m==2) then

    do k=kmin+1, kmax
       rhol(:,k,:) = rhot(:,k-1,:)
       Pl(:,k,:)   = Pt(:,k-1,:)
       vxl(:,k,:)  = vxt(:,k-1,:)
       vyl(:,k,:)  = vyt(:,k-1,:)
       vzl(:,k,:)  = vzt(:,k-1,:)
       Bxl(:,k,:)  = Bxt(:,k-1,:)
       Byl(:,k,:)  = Byt(:,k-1,:)
       Bzl(:,k,:)  = Bzt(:,k-1,:)
    end do

 else if (m==3) then

    do k=kmin+1, kmax
       rhol(:,:,k) = rhot(:,:,k-1)
       Pl(:,:,k)   = Pt(:,:,k-1)
       vxl(:,:,k)  = vxt(:,:,k-1)
       vyl(:,:,k)  = vyt(:,:,k-1)
       vzl(:,:,k)  = vzt(:,:,k-1)
       Bxl(:,:,k)  = Bxt(:,:,k-1)
       Byl(:,:,k)  = Byt(:,:,k-1)
       Bzl(:,:,k)  = Bzt(:,:,k-1)
    end do

 else
    write(*,*) "mhdppm_shift: m not in proper range",m
 end if
! write(*,*) "inside mhdppm_shift0",vzl(2,3,2),vzr(2,3,2),Bzl(2,3,2),Bzr(2,3,2),rhol(2,3,2),rhor(2,3,2)
! write(*,*) "inside mhdppm_shift1",vzl(15,2,2),vzr(15,2,2),Bzl(15,2,2),Bzr(15,2,2),rhol(15,2,2),rhor(15,2,2)
! write(*,*) "inside mhdppm_shift2",rhor(15,2,2),Pr(15,2,2),vxr(15,2,2),Bxr(15,2,2),vyr(15,2,2),Byr(15,2,2),vzr(15,2,2),Bzr(15,2,2)

end subroutine mhdppm_shift
!-----------------------------------------------------------------!
! apply flattening                                                !
!-----------------------------------------------------------------!
!
subroutine mhdppm_flatten(ext,X,Y,Z,ftilde,rho_b,rho_br,rho_bl,P,Pr,Pl, &
     vx,vxr,vxl,vy,vyr,vyl,vz,vzr,vzl, &
     Bx,Bxr,Bxl,By,Byr,Byl,Bz,Bzr,Bzl, &
     m,Symmetry)
 implicit none
 integer, dimension(3)                          :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 real*8, dimension(ext(1),ext(2),ext(3))        :: rho_b, rho_br, rho_bl
 real*8, dimension(ext(1),ext(2),ext(3))        :: ftilde,P,Pr,Pl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vx,vxr,vxl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vy,vyr,vyl
 real*8, dimension(ext(1),ext(2),ext(3))        :: vz,vzr,vzl
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bx,Bxr,Bxl
 real*8, dimension(ext(1),ext(2),ext(3))        :: By,Byr,Byl
 real*8, dimension(ext(1),ext(2),ext(3))        :: Bz,Bzr,Bzl
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
    if(Y(1,1,1) .lt. ZERO) then
       jmin = jmin + 1
    end if
    if(Z(1,1,1) .lt. ZERO) then
       kmin = kmin + 1
    end if
 else if(Symmetry == EQUATORIAL) then 
    if(Z(1,1,1) .lt. ZERO) then
       kmin = kmin + 1
    end if
 end if


! write(*,*) "inside mhdppm_flatten0",rho_br(15,2,2),Pr(15,2,2),vxr(15,2,2),Bxr(15,2,2),vyr(15,2,2),Byr(15,2,2),vzr(15,2,2),Bzr(15,2,2),ftilde(15,2,2)

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
             Bxr(i,j,k)      = Bx(i,j,k)*f    + Bxr(i,j,k)*(ONE-f)
             Bxl(i+1,j,k)    = Bx(i,j,k)*f    + Bxl(i+1,j,k)*(ONE-f)
             Byr(i,j,k)      = By(i,j,k)*f    + Byr(i,j,k)*(ONE-f)
             Byl(i+1,j,k)    = By(i,j,k)*f    + Byl(i+1,j,k)*(ONE-f)
             Bzr(i,j,k)      = Bz(i,j,k)*f    + Bzr(i,j,k)*(ONE-f)
             Bzl(i+1,j,k)    = Bz(i,j,k)*f    + Bzl(i+1,j,k)*(ONE-f)

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
             Bxr(i,j,k)      = Bx(i,j,k)*f    + Bxr(i,j,k)*(ONE-f)
             Bxl(i,j+1,k)    = Bx(i,j,k)*f    + Bxl(i,j+1,k)*(ONE-f)
             Byr(i,j,k)      = By(i,j,k)*f    + Byr(i,j,k)*(ONE-f)
             Byl(i,j+1,k)    = By(i,j,k)*f    + Byl(i,j+1,k)*(ONE-f)
             Bzr(i,j,k)      = Bz(i,j,k)*f    + Bzr(i,j,k)*(ONE-f)
             Bzl(i,j+1,k)    = Bz(i,j,k)*f    + Bzl(i,j+1,k)*(ONE-f)

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
             Bxr(i,j,k)      = Bx(i,j,k)*f    + Bxr(i,j,k)*(ONE-f)
             Bxl(i,j,k+1)    = Bx(i,j,k)*f    + Bxl(i,j,k+1)*(ONE-f)
             Byr(i,j,k)      = By(i,j,k)*f    + Byr(i,j,k)*(ONE-f)
             Byl(i,j,k+1)    = By(i,j,k)*f    + Byl(i,j,k+1)*(ONE-f)
             Bzr(i,j,k)      = Bz(i,j,k)*f    + Bzr(i,j,k)*(ONE-f)
             Bzl(i,j,k+1)    = Bz(i,j,k)*f    + Bzl(i,j,k+1)*(ONE-f)

          end do
       end do
    end do

 end if

! write(*,*) "inside mhdppm_flatten1",rho_br(15,2,2),Pr(15,2,2),vxr(15,2,2),Bxr(15,2,2),vyr(15,2,2),Byr(15,2,2),vzr(15,2,2),Bzr(15,2,2),ftilde(15,2,2)
! write(*,*) "inside mhdppm_flatten2",rho_bl(15,2,2),Pl(15,2,2),vxl(15,2,2),Bxl(15,2,2),vyl(15,2,2),Byl(15,2,2),vzl(15,2,2),Bzl(15,2,2),ftilde(15,2,2)
end subroutine mhdppm_flatten
