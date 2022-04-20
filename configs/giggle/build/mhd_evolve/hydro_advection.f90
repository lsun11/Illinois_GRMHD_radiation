!
!-----------------------------------------------------------------------------
!
! Advect rho_star and e
!
!-----------------------------------------------------------------------------
!
subroutine advect_rho_e(ext,rho_star_rhs,e_rhs,frm,fem,m,X,Y,Z,Symmetry)
 implicit none
 integer, dimension(3)                                :: ext
 real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
 integer                                        :: m,mmax,mmin,i,j,k
 integer                                        :: Symmetry
 real*8                                                :: dX,dY,dZ
 real*8, dimension(ext(1),ext(2),ext(3))        :: rho_star_rhs,e_rhs
 real*8, dimension(ext(1),ext(2),ext(3))        :: frm,fem
 integer, parameter                                :: AXISYM = 4
!
 mmax = ubound(frm,m)
 mmin = lbound(frm,m)
! write(*,*) inside advect_rho_e, m,mmax,rho_star_rhs(15,2,2),e_rhs(15,2,2),frm(15,2,2),fem(15,2,2),fem(15,2,2),frm(15,2,2),frm(16,2
 if (m==1) then
    dX = X(mmin+1,1,1)-X(mmin,1,1)
  if (Symmetry==AXISYM) then
    do i=mmin,mmax-1
       rho_star_rhs(i,:,:) = rho_star_rhs(i,:,:) +  &
                (frm(i,:,:)-frm(i+1,:,:))/dX/X(i,1,1)
       e_rhs(i,:,:) = e_rhs(i,:,:) + (fem(i,:,:)-fem(i+1,:,:))/dX/X(i,1,1)
    end do
  else
    do i=mmin,mmax-1
       rho_star_rhs(i,:,:) = rho_star_rhs(i,:,:) + (frm(i,:,:)-frm(i+1,:,:))/dX
       e_rhs(i,:,:) = e_rhs(i,:,:) + (fem(i,:,:)-fem(i+1,:,:))/dX
    end do
  end if
    e_rhs(mmax,:,:) = 0.d0
    rho_star_rhs(mmax,:,:) = 0.d0
    !frm = drho_bm
    !fem = dP_m
 !   write(*,*) inside advect_rho_e, m,mmax,rho_star_rhs(15,2,2),e_rhs(15,2,2),frm(15,2,2),fem(15,2,2),fem(3,2,2),frm(16,2,2),frm(17
 elseif (m==2) then
    if (Symmetry .ne. AXISYM) then
       dY = Y(1,mmin+1,1)-Y(1,mmin,1)
       do j=mmin,mmax-1
          rho_star_rhs(:,j,:) = rho_star_rhs(:,j,:) + (frm(:,j,:)-frm(:,j+1,:))/dY
          e_rhs(:,j,:) = e_rhs(:,j,:) + (fem(:,j,:)-fem(:,j+1,:))/dY
       end do
    end if
    e_rhs(:,mmax,:) = 0.d0
    rho_star_rhs(:,mmax,:) = 0.d0
 !   write(*,*) inside advect_rho_e, m,mmax,rho_star_rhs(15,2,2),e_rhs(15,2,2),frm(15,2,2),fem(15,2,2),fem(2,3,2),frm(16,2,2),frm(17
 elseif (m==3) then
    dZ = Z(1,1,mmin+1)-Z(1,1,mmin)
    do k=mmin,mmax-1
       rho_star_rhs(:,:,k) = rho_star_rhs(:,:,k) + (frm(:,:,k)-frm(:,:,k+1))/dZ
       e_rhs(:,:,k) = e_rhs(:,:,k) + (fem(:,:,k)-fem(:,:,k+1))/dZ
    end do
    e_rhs(:,:,mmax) = 0.d0
    rho_star_rhs(:,:,mmax) = 0.d0
 !   write(*,*) inside advect_rho_e, m,mmax,rho_star_rhs(15,2,2),e_rhs(15,2,2),frm(15,2,2),fem(15,2,2),fem(15,2,2),frm(15,2,2),frm(1
 end if
 return
end subroutine advect_rho_e
!-----------------------------------------------------------------!
! Calculate face averages
!-----------------------------------------------------------------!
!
subroutine face_avg(ext,X,Y,Z,func,func_f,m,Sym_arg,Symmetry)
  implicit none
  integer, dimension(3)                               :: ext
  real*8, dimension(ext(1),ext(2),ext(3)), intent(in) :: func
  real*8, dimension(ext(1),ext(2),ext(3))             :: func_f
  real*8, dimension(ext(1),ext(2),ext(3))              :: X,Y,Z
  integer                                             :: m, i,j,k
  integer                                             :: mmax, mmin
  real*8, parameter                                   :: TWO = 2.d0, ZERO = 0.d0
  real*8                                              :: Sym_arg
  real*8                                              :: SYM, ANTI,Zmin
  integer                            :: Symmetry, NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(SYM = 1.d0, ANTI = -1.d0)
!
  mmax = ubound(func,m)
  mmin = lbound(func,m)
  if(Symmetry == OCTANT) then
     if(m==1 .and. X(1,1,1) .lt. ZERO) then
        mmin = lbound(func_f,1)+1
     else if(m==2 .and. Y(1,1,1) .lt. ZERO) then
        mmin = lbound(func_f,2)+1
     else if(m==3 .and. Z(1,1,1) .lt. ZERO) then
        mmin = lbound(func_f,3)+1
     end if
  else if(Symmetry == AXISYM) then
      if(m==1 .and. X(1,1,1) .lt. ZERO) then
        mmin = lbound(func_f,1)+1
     else if(m==3 .and. Z(1,1,1) .lt. ZERO) then
        mmin = lbound(func_f,3)+1
     end if
  else if(Symmetry == EQUATORIAL) then
     if(m==3 .and. Z(1,1,1) .lt. ZERO) then
        mmin = lbound(func_f,3)+1
     end if
  end if
  if (m==1) then
    i = mmin
    if (Sym_arg > 0) then
       func_f(i,:,:) = func(i,:,:)
    else
       func_f(i,:,:) = 0.d0
    end if
    do i=mmin + 1, mmax
      func_f(i,:,:) = (func(i,:,:) + func(i-1,:,:))/TWO
    end do
  end if
  if (m==2) then
    j = mmin
    if (Sym_arg > 0) then
       func_f(:,j,:) = func(:,j,:)
    else
       func_f(:,j,:) = 0.d0
    end if
    do j=mmin + 1, mmax
      func_f(:,j,:) = (func(:,j,:) + func(:,j-1,:))/TWO
    end do
  end if
  if (m==3) then
    k = mmin
    if (Sym_arg > 0) then
       func_f(:,:,k) = func(:,:,k)
    else
       func_f(:,:,k) = 0.d0
    end if
    do k=mmin + 1, mmax
      func_f(:,:,k) = (func(:,:,k) + func(:,:,k-1))/TWO
    end do
  end if
end subroutine face_avg
!-----------------------------------------------------------------!
! Calculate face averages to third order accuracy
!-----------------------------------------------------------------!
!
subroutine face_avg3(ext,func,func_f,m,Sym_arg)
  implicit none
  integer, dimension(3)                               :: ext
  real*8, dimension(ext(1),ext(2),ext(3)), intent(in) :: func
  real*8, dimension(ext(1),ext(2),ext(3))             :: func_f
  integer                                             :: m, i,j,k
  integer                                             :: mmax, mmin
  real*8                                              :: Sym_arg
  real*8                                              :: SYM, ANTI
  real*8, parameter        :: f3o8=0.375d0, f3o4=0.75d0, f1o8=0.125d0
  parameter(SYM = 1.d0, ANTI = -1.d0)
!
  mmax = ubound(func,m)
  mmin = lbound(func,m)
  if (m==1) then
    func_f(mmin,:,:) = (f3o4 + Sym_arg*f3o8)*func(mmin,:,:) &
                        - f1o8*func(mmin+1,:,:)
    do i=mmin + 1, mmax - 1
      func_f(i,:,:) = f3o4*func(i,:,:) + f3o8*func(i-1,:,:) &
                        - f1o8*func(i+1,:,:)
    end do
    func_f(mmax,:,:) = f3o8*func(mmax,:,:) + f3o4*func(mmax-1,:,:) &
                        - f1o8*func(mmax-2,:,:)
  end if
  if (m==2) then
    func_f(:,mmin,:) = (f3o4 + Sym_arg*f3o8)*func(:,mmin,:) &
                        - f1o8*func(:,mmin+1,:)
    do j=mmin + 1, mmax - 1
      func_f(:,j,:) = f3o4*func(:,j,:) + f3o8*func(:,j-1,:) &
                        - f1o8*func(:,j+1,:)
    end do
    func_f(:,mmax,:) = f3o8*func(:,mmax,:) + f3o4*func(:,mmax-1,:) &
                        - f1o8*func(:,mmax-2,:)
  end if
  if (m==3) then
    func_f(:,:,mmin) = (f3o4 + Sym_arg*f3o8)*func(:,:,mmin) &
                        - f1o8*func(:,:,mmin+1)
    do k=mmin + 1, mmax - 1
      func_f(:,:,k) = f3o4*func(:,:,k) + f3o8*func(:,:,k-1) &
                        - f1o8*func(:,:,k+1)
    end do
    func_f(:,:,mmax) = f3o8*func(:,:,mmax) + f3o4*func(:,:,mmax-1) &
                        - f1o8*func(:,:,mmax-2)
  end if
end subroutine face_avg3
!-----------------------------------------------------------------!
! Calculate nablas for hydro quantities
!-----------------------------------------------------------------!
!
subroutine compute_nabla_hydro(ext,dX,dY,dZ,vx,dvx_m,vy,dvy_m,vz,dvz_m, &
                       rho_b,drho_b_m,P,dP_m,m,Symmetry,X,Y,Z,Reconstruction)
  implicit none
 interface
    subroutine compute_nabla(ext,dX,dY,dZ,U,nabla_U,m,Sym_arg,X,Y,Z,LIMITER,Symmetry)
      implicit none
      integer, dimension(3)                      :: ext
      integer                                    :: m, LIMITER,Symmetry
      real*8                                     :: Sym_arg
      real*8, dimension(ext(1),ext(2),ext(3))    :: U, nabla_U
      real*8, dimension(ext(1),ext(2),ext(3))    :: X,Y,Z
      real*8                                     :: dX, dY, dZ
    end subroutine compute_nabla
 end interface
  integer, dimension(3)                          :: ext
  real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
  real*8, dimension(ext(1),ext(2),ext(3))        :: vx, vy, vz, rho_star
  real*8, dimension(ext(1),ext(2),ext(3))        :: dvx_m, dvy_m, dvz_m
  real*8, dimension(ext(1),ext(2),ext(3))        :: rho_b, drho_b_m, P, dP_m
  integer                                        :: m, Symmetry, Reconstruction
  real*8, parameter                              :: TWO = 2.d0
  integer                            :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM
  integer                            :: AXISYM
  integer                            :: SPPM, LIMITER,dU,MC
  real*8                             :: SYM, ANTI, dX, dY, dZ, rho_tiny
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(SPPM=5, MC=1, dU=4)
  parameter(SYM = 1.d0, ANTI = -1.d0)
!
  if (Reconstruction==SPPM) then
     LIMITER = dU
  else
     LIMITER = MC
  end if
  if (Symmetry==OCTANT) then
     if (m==1) then
        call compute_nabla(ext,dX,dY,dZ,vx,dvx_m,m,ANTI,X,Y,Z,LIMITER,Symmetry)
     else
        call compute_nabla(ext,dX,dY,dZ,vx,dvx_m,m,SYM,X,Y,Z,LIMITER,Symmetry)
     end if
     if (m==2) then
        call compute_nabla(ext,dX,dY,dZ,vy,dvy_m,m,ANTI,X,Y,Z,LIMITER,Symmetry)
     else
        call compute_nabla(ext,dX,dY,dZ,vy,dvy_m,m,SYM,X,Y,Z,LIMITER,Symmetry)
     end if
     if (m==3) then
        call compute_nabla(ext,dX,dY,dZ,vz,dvz_m,m,ANTI,X,Y,Z,LIMITER,Symmetry)
     else
        call compute_nabla(ext,dX,dY,dZ,vz,dvz_m,m,SYM,X,Y,Z,LIMITER,Symmetry)
     end if
  elseif (Symmetry==AXISYM) then
     if (m==1) then
        call compute_nabla(ext,dX,dY,dZ,vx,dvx_m,m,ANTI,X,Y,Z,LIMITER,Symmetry)
        call compute_nabla(ext,dX,dY,dZ,vy,dvy_m,m,ANTI,X,Y,Z,LIMITER,Symmetry)
     else
        call compute_nabla(ext,dX,dY,dZ,vx,dvx_m,m,SYM,X,Y,Z,LIMITER,Symmetry)
        call compute_nabla(ext,dX,dY,dZ,vy,dvy_m,m,SYM,X,Y,Z,LIMITER,Symmetry)
     end if
     if (m==3) then
        call compute_nabla(ext,dX,dY,dZ,vz,dvz_m,m,ANTI,X,Y,Z,LIMITER,Symmetry)
     else
        call compute_nabla(ext,dX,dY,dZ,vz,dvz_m,m,SYM,X,Y,Z,LIMITER,Symmetry)
     end if
  else if (Symmetry==EQUATORIAL .or. Symmetry==PI_SYMM) then
     call compute_nabla(ext,dX,dY,dZ,vx,dvx_m,m,SYM,X,Y,Z,LIMITER,Symmetry)
     call compute_nabla(ext,dX,dY,dZ,vy,dvy_m,m,SYM,X,Y,Z,LIMITER,Symmetry)
     if (m==3) then
        call compute_nabla(ext,dX,dY,dZ,vz,dvz_m,m,ANTI,X,Y,Z,LIMITER,Symmetry)
     else
        call compute_nabla(ext,dX,dY,dZ,vz,dvz_m,m,SYM,X,Y,Z,LIMITER,Symmetry)
     end if
  else
    call compute_nabla(ext,dX,dY,dZ,vx,dvx_m,m,SYM,X,Y,Z,LIMITER,Symmetry)
    call compute_nabla(ext,dX,dY,dZ,vy,dvy_m,m,SYM,X,Y,Z,LIMITER,Symmetry)
    call compute_nabla(ext,dX,dY,dZ,vz,dvz_m,m,SYM,X,Y,Z,LIMITER,Symmetry)
  end if
  call compute_nabla(ext,dX,dY,dZ,rho_b,drho_b_m,m,SYM,X,Y,Z,LIMITER,Symmetry)
  call compute_nabla(ext,dX,dY,dZ,P,dP_m,m,SYM,X,Y,Z,LIMITER,Symmetry)
!  write(*,*) NABLA_HYDRO: ,m,rho_b(15,2,2),drho_b_m(15,2,2),P(15,2,2),dP_m(15,2,2),dvx_m(15,2,2),vx(15,2,2),Symmetry,vx(15,2,2),vy(
end subroutine compute_nabla_hydro
!-----------------------------------------------------------------!
! Calculate nablas for hydro quantities
!-----------------------------------------------------------------!
!
subroutine find_ur_ul_hydro(ext, X,Y,Z, vx, vxr, vxl, dvx_m, &
                            vy, vyr, vyl, dvy_m, &
                            vz, vzr, vzl, dvz_m, &
                            rho_b,rho_br, rho_bl, drho_b_m, &
                            P, Pr, Pl, dP_m, &
                            m, Symmetry,Reconstruction)
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
  real*8, dimension(ext(1),ext(2),ext(3))        :: X,Y,Z
  real*8, dimension(ext(1),ext(2),ext(3))        :: vx, vy, vz,v
  real*8, dimension(ext(1),ext(2),ext(3))        :: vxl, vyl, vzl, vxr, vyr, vzr
  real*8, dimension(ext(1),ext(2),ext(3))        :: dvx_m, dvy_m, dvz_m
  real*8, dimension(ext(1),ext(2),ext(3))        :: rho_b, drho_b_m, rho_br, rho_bl
  real*8, dimension(ext(1),ext(2),ext(3))        :: P, Pr, Pl, dP_m
  integer                                        :: m, Symmetry, kmin
  integer                                        :: Reconstruction,SPPM
  real*8, parameter                              :: TWO = 2.d0
  integer                            :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM
  integer                            :: AXISYM
  real*8                             :: SYM, ANTI, fac, Zmin
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(SPPM = 5)
  parameter(SYM = 1.d0, ANTI = -1.d0)
!
  if (Symmetry==OCTANT) then
     if (m==1) then
        call find_ur_ul(ext,X,Y,Z,vx,vxr,vxl,dvx_m,m,ANTI,Reconstruction,Symmetry)
     else
        call find_ur_ul(ext,X,Y,Z,vx,vxr,vxl,dvx_m,m,SYM,Reconstruction,Symmetry)
     end if
     if (m==2) then
        call find_ur_ul(ext,X,Y,Z,vy,vyr,vyl,dvy_m,m,ANTI,Reconstruction,Symmetry)
     else
        call find_ur_ul(ext,X,Y,Z,vy,vyr,vyl,dvy_m,m,SYM,Reconstruction,Symmetry)
     end if
     if (m==3) then
        call find_ur_ul(ext,X,Y,Z,vz,vzr,vzl,dvz_m,m,ANTI,Reconstruction,Symmetry)
     else
        call find_ur_ul(ext,X,Y,Z,vz,vzr,vzl,dvz_m,m,SYM,Reconstruction,Symmetry)
     end if
  elseif (Symmetry==AXISYM) then
     if (m==1) then
        call find_ur_ul(ext,X,Y,Z,vx,vxr,vxl,dvx_m,m,ANTI,Reconstruction,Symmetry)
        call find_ur_ul(ext,X,Y,Z,vy,vyr,vyl,dvy_m,m,ANTI,Reconstruction,Symmetry)
     else
        call find_ur_ul(ext,X,Y,Z,vx,vxr,vxl,dvx_m,m,SYM,Reconstruction,Symmetry)
        call find_ur_ul(ext,X,Y,Z,vy,vyr,vyl,dvy_m,m,SYM,Reconstruction,Symmetry)
     end if
     kmin = lbound(vx,3)+1
     Zmin = Z(1,1,kmin)
     if (m==3 .and. Zmin .gt. 0.d0) then
        call find_ur_ul(ext,X,Y,Z,vz,vzr,vzl,dvz_m,m,ANTI,Reconstruction,Symmetry)
     else
        call find_ur_ul(ext,X,Y,Z,vz,vzr,vzl,dvz_m,m,SYM,Reconstruction,Symmetry)
     end if
  else if (Symmetry==EQUATORIAL .or. Symmetry==PI_SYMM) then
     if(Symmetry==EQUATORIAL) then
        kmin = lbound(vx,3)+1
        Zmin = Z(1,1,kmin)
     else
        kmin = lbound(vx,3)
        Zmin = Z(1,1,kmin)
     end if
     call find_ur_ul(ext,X,Y,Z,vx,vxr,vxl,dvx_m,m,SYM,Reconstruction,Symmetry)
     call find_ur_ul(ext,X,Y,Z,vy,vyr,vyl,dvy_m,m,SYM,Reconstruction,Symmetry)
     if (m==3 .and. Zmin .gt. 0.d0) then
        call find_ur_ul(ext,X,Y,Z,vz,vzr,vzl,dvz_m,m,ANTI,Reconstruction,Symmetry)
     else
        call find_ur_ul(ext,X,Y,Z,vz,vzr,vzl,dvz_m,m,SYM,Reconstruction,Symmetry)
     end if
  else
     call find_ur_ul(ext,X,Y,Z,vx,vxr,vxl,dvx_m,m,SYM,Reconstruction,Symmetry)
     call find_ur_ul(ext,X,Y,Z,vy,vyr,vyl,dvy_m,m,SYM,Reconstruction,Symmetry)
     call find_ur_ul(ext,X,Y,Z,vz,vzr,vzl,dvz_m,m,SYM,Reconstruction,Symmetry)
  end if
  call find_ur_ul(ext,X,Y,Z,rho_b,rho_br,rho_bl,drho_b_m,m,SYM,Reconstruction,Symmetry)
  !write(*,*) hewwo1.,ext
  call find_ur_ul(ext,X,Y,Z,P,Pr,Pl,dP_m,m,SYM,Reconstruction,Symmetry)
! write(*,*) ur_ul hydro: SYM,Reconstruction,Symmetry,ext,m,SYM,Reconstruction,Symmetry
  ! Check for negative pressure and density resulting from interpolation
  where (Pl .lt. 0.d0 .or. rho_bl .lt. 0.d0)
     Pl = 0.d0
     rho_bl = 0.d0
  end where
  where (Pr .lt. 0.d0 .or. rho_br .lt. 0.d0)
     Pr = 0.d0
     rho_br = 0.d0
  end where
!  write(*,*) find_ur_ul_hydro,rho_bl(15,2,2),Pl(15,2,2),vxr(16,2,2),vyr(16,2,2),vzr(16,2,2),vxr(15,2,2),vyr(15,2,2),vzr(15,2,2)
end subroutine find_ur_ul_hydro
!-----------------------------------------------------------------!
! Calculate nablas for hydro quantities
!-----------------------------------------------------------------!
!
subroutine face_avg_hydro(ext,X,Y,Z,alpha,alpha_f,betax,betax_f, &
                          betay, betay_f, betaz, betaz_f, &
                          gxx, gxx_f, gxy, gxy_f, gxz, gxz_f, &
                          gyy, gyy_f, gyz, gyz_f, gzz, gzz_f, &
                          phi, phi_f, gupxx, gupxx_f, gupyy, &
                          gupyy_f, gupzz, gupzz_f, m, Symmetry)
  implicit none
 interface
    subroutine face_avg(ext,X,Y,Z,func,func_f,m,Sym_arg,Symmetry)
      implicit none
      integer, dimension(3)                      :: ext
      integer                                    :: Symmetry
      real*8, dimension(ext(1),ext(2),ext(3))    :: X,Y,Z
      integer                                    :: m
      real*8                                     :: Sym_arg
      real*8, dimension(ext(1),ext(2),ext(3))    :: func, func_f
    end subroutine face_avg
 end interface
  integer, dimension(3)                          :: ext
  real*8, dimension(ext(1),ext(2),ext(3))              :: X,Y,Z
  real*8, dimension(ext(1),ext(2),ext(3))        :: alpha, alpha_f, betax, betax_f
  real*8, dimension(ext(1),ext(2),ext(3))        :: betay, betay_f, betaz, betaz_f
  real*8, dimension(ext(1),ext(2),ext(3))        :: gxx, gxx_f, gxy, gxy_f, gxz, gxz_f
  real*8, dimension(ext(1),ext(2),ext(3))        :: gyy, gyy_f, gyz, gyz_f, gzz, gzz_f
  real*8, dimension(ext(1),ext(2),ext(3))        :: gupxx, gupxx_f, gupyy, gupyy_f
  real*8, dimension(ext(1),ext(2),ext(3))        :: gupzz, gupzz_f
  real*8, dimension(ext(1),ext(2),ext(3))        :: phi, phi_f
  integer                                        :: m, Symmetry
  real*8, parameter                              :: TWO = 2.d0
  integer                            :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM
  integer                            :: AXISYM
  real*8                             :: SYM, ANTI, Zmin
  integer                             :: kmin
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(SYM = 1.d0, ANTI = -1.d0)
! 
  kmin = lbound(phi,3)+1
  Zmin = Z(1,1,kmin)
!
  if (Symmetry==OCTANT) then
     ! Calculate face averages of \beta^i
     if (m==1) then
        call face_avg(ext,X,Y,Z,betax,betax_f,m,ANTI,Symmetry)
     else
        call face_avg(ext,X,Y,Z,betax,betax_f,m,SYM,Symmetry)
     end if
     if (m==2) then
        call face_avg(ext,X,Y,Z,betay,betay_f,m,ANTI,Symmetry)
     else
        call face_avg(ext,X,Y,Z,betay,betay_f,m,SYM,Symmetry)
     end if
     if (m==3) then
        call face_avg(ext,X,Y,Z,betaz,betaz_f,m,ANTI,Symmetry)
     else
        call face_avg(ext,X,Y,Z,betaz,betaz_f,m,SYM,Symmetry)
     end if
     if (m==1 .or. m==2) then
        call face_avg(ext,X,Y,Z,gxy,gxy_f,m,ANTI,Symmetry)
     else
        call face_avg(ext,X,Y,Z,gxy,gxy_f,m,SYM,Symmetry)
     end if
     if (m==1 .or. m==3) then
        call face_avg(ext,X,Y,Z,gxz,gxz_f,m,ANTI,Symmetry)
     else
        call face_avg(ext,X,Y,Z,gxz,gxz_f,m,SYM,Symmetry)
     end if
     if (m==2 .or. m==3) then
        call face_avg(ext,X,Y,Z,gyz,gyz_f,m,ANTI,Symmetry)
     else
        call face_avg(ext,X,Y,Z,gyz,gyz_f,m,SYM,Symmetry)
     end if
  else if (Symmetry==AXISYM) then
     if (m==1) then
        call face_avg(ext,X,Y,Z,betax,betax_f,m,ANTI,Symmetry)
        call face_avg(ext,X,Y,Z,betay,betay_f,m,ANTI,Symmetry)
     else
        call face_avg(ext,X,Y,Z,betax,betax_f,m,SYM,Symmetry)
        call face_avg(ext,X,Y,Z,betay,betay_f,m,SYM,Symmetry)
     end if
     if (m==3 .and. Zmin .gt. 0.d0) then
        call face_avg(ext,X,Y,Z,betaz,betaz_f,m,ANTI,Symmetry)
     else
        call face_avg(ext,X,Y,Z,betaz,betaz_f,m,SYM,Symmetry)
     end if
     call face_avg(ext,X,Y,Z,gxy,gxy_f,m,SYM,Symmetry)
     call face_avg(ext,X,Y,Z,gxz,gxz_f,m,ANTI,Symmetry)
     call face_avg(ext,X,Y,Z,gyz,gyz_f,m,ANTI,Symmetry)
  else
     call face_avg(ext,X,Y,Z,betax,betax_f,m,SYM,Symmetry)
     call face_avg(ext,X,Y,Z,betay,betay_f,m,SYM,Symmetry)
     call face_avg(ext,X,Y,Z,betaz,betaz_f,m,SYM,Symmetry)
     call face_avg(ext,X,Y,Z,gxy,gxy_f,m,SYM,Symmetry)
     call face_avg(ext,X,Y,Z,gxz,gxz_f,m,SYM,Symmetry)
     call face_avg(ext,X,Y,Z,gyz,gyz_f,m,SYM,Symmetry)
  end if
  if(Symmetry==EQUATORIAL .and. m==3) then
     call face_avg(ext,X,Y,Z,betaz,betaz_f,m,ANTI,Symmetry)
     call face_avg(ext,X,Y,Z,gxz,gxz_f,m,ANTI,Symmetry)
     call face_avg(ext,X,Y,Z,gyz,gyz_f,m,ANTI,Symmetry)
  end if
  call face_avg(ext,X,Y,Z,gxx,gxx_f,m,SYM,Symmetry)
  call face_avg(ext,X,Y,Z,gyy,gyy_f,m,SYM,Symmetry)
  call face_avg(ext,X,Y,Z,gzz,gzz_f,m,SYM,Symmetry)
  call face_avg(ext,X,Y,Z,phi,phi_f,m,SYM,Symmetry)
  call face_avg(ext,X,Y,Z,alpha,alpha_f,m,SYM,Symmetry)
  call face_avg(ext,X,Y,Z,gupxx,gupxx_f,m,SYM,Symmetry)
  call face_avg(ext,X,Y,Z,gupyy,gupyy_f,m,SYM,Symmetry)
  call face_avg(ext,X,Y,Z,gupzz,gupzz_f,m,SYM,Symmetry)
!  write(*,*) face_avg,gxx(15,2,2),gxy(15,2,2),gxz(15,2,2),gyy(15,2,2),gyz(15,2,2),gzz(15,2,2),gxx_f(15,2,2),gxy_f(15,2,2),gxz_f(15,
end subroutine face_avg_hydro
subroutine sppm_shift_hydro_ul(ext,Z,rho_br,rho_bl,Pr,Pl, &
                vxr,vxl,vyr,vyl,vzr,vzl,m,Symmetry)
 implicit none
 integer, dimension(3)                                        :: ext
 real*8, dimension(ext(1),ext(2),ext(3))                :: rho_bl,Pl
 real*8, dimension(ext(1),ext(2),ext(3))                :: vxl,vyl,vzl
 real*8, dimension(ext(1),ext(2),ext(3))                :: rho_br,Pr
 real*8, dimension(ext(1),ext(2),ext(3))                :: vxr,vyr,vzr
 real*8, dimension(ext(1),ext(2),ext(3))                :: Z
 integer                                                :: m,Symmetry
 integer                                                :: mmin,mmax
 integer                                                :: i,j,k
 integer, parameter  :: NO_SYMM=0, EQUATORIAL=1, OCTANT=2, PI_SYMM=3, AXISYM=4
!
 mmin = lbound(Pl,m)
 mmax = ubound(Pl,m)
 if (m==1) then
    do i=mmax,mmin+1,-1
       rho_bl(i,:,:) = rho_bl(i-1,:,:)
       Pl(i,:,:) = Pl(i-1,:,:)
       vxl(i,:,:) = vxl(i-1,:,:)
       vyl(i,:,:) = vyl(i-1,:,:)
       vzl(i,:,:) = vzl(i-1,:,:)
    end do
    if (Symmetry==OCTANT) then
       rho_bl(mmin,:,:) = rho_br(mmin,:,:)
       Pl(mmin,:,:) = Pr(mmin,:,:)
       vxl(mmin,:,:) = -vxr(mmin,:,:)
       vyl(mmin,:,:) = vyr(mmin,:,:)
       vzl(mmin,:,:) = vzr(mmin,:,:)
    else if (Symmetry==AXISYM) then
       rho_bl(mmin,:,:) = rho_br(mmin,:,:)
       Pl(mmin,:,:) = Pr(mmin,:,:)
       vxl(mmin,:,:) = -vxr(mmin,:,:)
       vyl(mmin,:,:) = -vyr(mmin,:,:)
       vzl(mmin,:,:) = vzr(mmin,:,:)
    end if
    return
 end if
 if (m==2) then
    do j=mmax,mmin+1,-1
       rho_bl(:,j,:) = rho_bl(:,j-1,:)
       Pl(:,j,:) = Pl(:,j-1,:)
       vxl(:,j,:) = vxl(:,j-1,:)
       vyl(:,j,:) = vyl(:,j-1,:)
       vzl(:,j,:) = vzl(:,j-1,:)
    end do
    if (Symmetry==OCTANT) then
       rho_bl(:,mmin,:) = rho_br(:,mmin,:)
       Pl(:,mmin,:) = Pr(:,mmin,:)
       vxl(:,mmin,:) = vxr(:,mmin,:)
       vyl(:,mmin,:) = -vyr(:,mmin,:)
       vzl(:,mmin,:) = vzr(:,mmin,:)
    end if
    return
 end if
 if (m==3) then
    do k=mmax,mmin+1,-1
       rho_bl(:,:,k) = rho_bl(:,:,k-1)
       Pl(:,:,k) = Pl(:,:,k-1)
       vxl(:,:,k) = vxl(:,:,k-1)
       vyl(:,:,k) = vyl(:,:,k-1)
       vzl(:,:,k) = vzl(:,:,k-1)
    end do
    if (Symmetry==OCTANT .or. (Symmetry==AXISYM .and. Z(1,1,mmin) .gt. 0.d0)) then
       rho_bl(:,:,mmin) = rho_br(:,:,mmin)
       Pl(:,:,mmin) = Pr(:,:,mmin)
       vxl(:,:,mmin) = vxr(:,:,mmin)
       vyl(:,:,mmin) = vyr(:,:,mmin)
       vzl(:,:,mmin) = -vzr(:,:,mmin)
    end if
    return
 end if
end subroutine sppm_shift_hydro_ul
