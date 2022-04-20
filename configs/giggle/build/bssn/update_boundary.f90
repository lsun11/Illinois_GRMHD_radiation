!--------------------------------------------------------------------
! Update outer boundaries
! Note that we do not update Gammai if Spatial_Gauge=6 (GF shift).
!--------------------------------------------------------------------
subroutine update_boundary(ex,fake_ex,X,Y,Z,R_fish, PhysicalR, RadiusDerivative,Fisheye_Flag, &
     dT,dX,dY,dZ, alpha, phi_old, phi_new, &
     K_old, K_new, gxx_old,gxy_old,gxz_old,gyy_old,gyz_old,gzz_old,           &
     gxx_new,gxy_new,gxz_new,gyy_new,gyz_new,gzz_new, Axx_old, &
     Axy_old,Axz_old,Ayy_old,Ayz_old,Azz_old, Axx_new,Axy_new, &
     Axz_new,Ayy_new,Ayz_new,Azz_new, Gammax_old, Gammay_old,  &
     Gammaz_old, Gammax_new, Gammay_new, Gammaz_new, Symmetry, &
     Spatial_Gauge,Xcenter, Ycenter, Zcenter,bc_type,gambound, adm_mass, &
     have_bdry_min,have_bdry_max,use_trans_fish_phys_new)
  implicit none
  !
  ! Input parameters:
  !
  integer, dimension(3)                    :: ex,fake_ex,have_bdry_min,have_bdry_max
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z,R_fish
  real*8, dimension(ex(1),ex(2),ex(3))     :: PhysicalR,RadiusDerivative
  real*8                                   :: dT,dX,dY,dZ
  real*8                                   :: Xcenter,Ycenter,Zcenter
  real*8, dimension(ex(1),ex(2),ex(3))     :: alpha,phi_old,phi_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: K_old,K_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: gxx_old,gxy_old,gxz_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: gyy_old,gyz_old,gzz_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: gxx_new,gxy_new,gxz_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: gyy_new,gyz_new,gzz_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: Axx_old,Axy_old,Axz_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: Ayy_old,Ayz_old,Azz_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: Axx_new,Axy_new,Axz_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: Ayy_new,Ayz_new,Azz_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: Gammax_old, Gammay_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: Gammaz_old, Gammax_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: Gammay_new, Gammaz_new
  integer                                  :: Symmetry,Spatial_Gauge
  integer                                  :: Fisheye_Flag,gambound,use_trans_fish_phys_new
  !
  ! Other variables:
  !
  integer                    :: i, j, k
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  integer                    :: il,jl,kl,iu,ju,ku
  integer                    :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM&
       &, AXISYM
  integer                    :: WAVELIKE,FALLOFF,QE,bc_type,EXTRAP, COPY, QUAD
  !BCs for puncture (extrapr: r*f at bdry = 2r*f at bdry-1 - r*f at bdry-2 [good for 1/r falloff!]) :
  integer :: EXTRAPR, MoR, bc_type2
  real*8                     :: Courant,Dlx,Dly,Dlz,r,r1,Rp,Rp1,r2
  real*8                     :: facr,facRp,fac,fac1,fac2,rp2,adm_mass
  real*8                     :: f_int,f_int_phi,f_int_K,f_int_Gamx&
       &,f_int_Gamy,f_int_Gamz,f_int_dR
  real*8                     :: f_int_gxx,f_int_gxy,f_int_gxz&
       &,f_int_gyy,f_int_gyz,f_int_gzz
  real*8                     :: f_int_Axx,f_int_Axy,f_int_Axz&
       &,f_int_Ayy,f_int_Ayz,f_int_Azz, psi4_int,psi4
  real*8                     :: HALF, ONE, TWO, THREE, ZERO, sym_sa, sym_as
  real*8, parameter          :: f1o3=1.d0/3.d0,f1o6=1.d0/6.d0
  parameter(HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0, THREE = 3.D0, ZERO = 0.D0)
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3,&
       & AXISYM = 4)
  parameter(WAVELIKE = 0, FALLOFF = 1, QE = 2, EXTRAP = 3, COPY = 4, QUAD = 5)
  parameter(EXTRAPR = 6, MoR = 7)
  real*8, parameter          :: SYM = 1.D0, ANT = -1.D0
  ! The following were added for the purpose of axisymmetrization:
  integer                    :: jlb, jub
  logical                      :: fish_to_phys
  if(Symmetry==AXISYM) then
     sym_sa = ANT
     sym_as = SYM
  else
     sym_sa = SYM
     sym_as = ANT
  end if
  imin = ex(1)-fake_ex(1)+1
  jmin = ex(2)-fake_ex(2)+1
  kmin = ex(3)-fake_ex(3)+1
  imax = fake_ex(1)
  jmax = fake_ex(2)
  kmax = fake_ex(3)
if(1==0) then
do k=1, ex(3)
     do j=1, ex(2)
             do i=1, ex(1)
                if(isnan(Axx_old(i,j,k))) then
                    write(*,*) "bssn update bdy. F90:  Checkpoint << 1 >>, Axx_old = NAN FOUND!"
                       end if
                if(isnan(Axx_new(i,j,k))) then
                    write(*,*) "bssn update bdy. F90:  Checkpoint << 1 >>, Axx = NAN FOUND!"
                       end if
       end do
     end do
  end do
end if
  if (bc_type==MoR) then
     bc_type2 = WAVELIKE
  else
     bc_type2 = bc_type
  end if
  !Following is USEFUL when nghostzone > 1 !!!!
  write(*,*) "Calling upd_bdry with imin,imax =",imin,imax
  if(bc_type2==WAVELIKE) then
     if (bc_type==MoR) then
        ! Impose phi = M/2r + O(1/r^2) (in physical coordinates)
        phi_new = phi_new - log(PhysicalR/R_fish)/3.d0 - log(RadiusDerivative)/6.d0 &
             - 0.5d0*adm_mass/PhysicalR
        call update_falloff_fisheye2(ex, X, Y, Z, PhysicalR, phi_new, &
             0, 0, 0, 2 , Symmetry,imin,imax,jmin,jmax,kmin,kmax)
        phi_new = phi_new + log(PhysicalR/R_fish)/3.d0 + log(RadiusDerivative)/6.d0 &
             + 0.5d0*adm_mass/PhysicalR
     end if
     ! Transform components of gij and Aij from fisheye to physical coordinates
     if (Fisheye_Flag==1) then
        fish_to_phys = .TRUE.
        call trans_fish_phys(ex,fake_ex, X, Y, Z, PhysicalR, RadiusDerivative, &
             gxx_old,gxy_old,gxz_old,gyy_old,gyz_old,gzz_old, &
             fish_to_phys,2,have_bdry_min,have_bdry_max,use_trans_fish_phys_new)
        call trans_fish_phys(ex,fake_ex, X, Y, Z, PhysicalR, RadiusDerivative, &
             Axx_old,Axy_old,Axz_old,Ayy_old,Ayz_old,Azz_old, &
             fish_to_phys,2,have_bdry_min,have_bdry_max,use_trans_fish_phys_new)
     end if
     ! Change:  Set bounds on y loops
     if (Symmetry.ne.AXISYM) then
        jlb = jmin
        jub = jmax - 1
     else
        jlb = jmin + 1
        jub = jlb
     end if
     !
     ! x = max
     if(have_bdry_max(1)==1) then
        i = imax
!!!$omp parallel do
        do k = kmin, kmax - 1
           do j = jlb, jub
              call weights2(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),RadiusDerivative(i,j,k), &
                   & alpha(i,j,k),phi_old(i,j,k),&
                   & dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il&
                   &,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, Xcenter,&
                   & Ycenter, Zcenter,gxx_old(i,j,k),gxy_old(i,j,k),gxz_old(i,j,k),gyy_old(i,j,k),gyz_old(i,j,k),gzz_old(i,j,k))
              !
              ! Interpolate 
              !
              if (bc_type==WAVELIKE) then
                 call interp(ex,phi_old,f_int_phi,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,SYM,SYM,SYM)
              end if
              call interp(ex,K_old  ,f_int_K  ,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,gxx_old,f_int_gxx,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,gxy_old,f_int_gxy,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,sym_as,ANT,SYM)
              call interp(ex,gxz_old,f_int_gxz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,ANT,SYM,ANT)
              call interp(ex,gyy_old,f_int_gyy,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,gyz_old,f_int_gyz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,sym_sa,ANT,ANT)
              call interp(ex,gzz_old,f_int_gzz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,Axx_old,f_int_Axx,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,Axy_old,f_int_Axy,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,sym_as,ANT,SYM)
              call interp(ex,Axz_old,f_int_Axz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,ANT,SYM,ANT)
              call interp(ex,Ayy_old,f_int_Ayy,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,Ayz_old,f_int_Ayz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,sym_sa,ANT,ANT)
              call interp(ex,Azz_old,f_int_Azz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,RadiusDerivative,f_int_dR,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              if(gambound .ge. 1) then
                 call interp(ex,Gammax_old,f_int_Gamx,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,ANT,SYM,SYM)
                 call interp(ex,Gammay_old,f_int_Gamy,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,SYM,ANT,SYM)
                 call interp(ex,Gammaz_old,f_int_Gamz,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,SYM,SYM,ANT)
              endif
              r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
              r1 = facr*r
              Rp = PhysicalR(i,j,k)
              Rp1 = facRp*Rp
              if (bc_type==WAVELIKE) then
                 phi_new(i,j,k) = facRp*f_int_phi +  &
                      f1o3*(log(Rp/r)-facRp*log(Rp1/r1)) + &
                      f1o6*(log(RadiusDerivative(i,j,k))-facRp*log(f_int_dR))
              end if
              K_new(i,j,k)   = f_int_K * facRp
              gxx_new(i,j,k) = ONE + facRp*(f_int_gxx-ONE)
              gxy_new(i,j,k) = facRp*f_int_gxy
              gxz_new(i,j,k) = facRp*f_int_gxz
              gyy_new(i,j,k) = ONE + facRp*(f_int_gyy-ONE)
              gyz_new(i,j,k) = facRp*f_int_gyz
              gzz_new(i,j,k) = ONE + facRp*(f_int_gzz-ONE)
              Axx_new(i,j,k) = f_int_Axx * facRp
              Axy_new(i,j,k) = f_int_Axy * facRp
              Axz_new(i,j,k) = f_int_Axz * facRp
              Ayy_new(i,j,k) = f_int_Ayy * facRp
              Ayz_new(i,j,k) = f_int_Ayz * facRp
              Azz_new(i,j,k) = f_int_Azz * facRp
              if(gambound.eq.0) then
                 Gammax_new(i,j,k) = Gammax_old(i,j,k)
                 Gammay_new(i,j,k) = Gammay_old(i,j,k)
                 Gammaz_new(i,j,k) = Gammaz_old(i,j,k)
              else
                 Gammax_new(i,j,k) = f_int_Gamx * facRp**(gambound)
                 Gammay_new(i,j,k) = f_int_Gamy * facRp**(gambound)
                 Gammaz_new(i,j,k) = f_int_Gamz * facRp**(gambound)
              endif
           enddo
        enddo
!!!$omp end parallel do
     end if
     !
     ! y = ymax.  Do not do this for axisymmetric cases.
     !
     if(have_bdry_max(2)==1) then
        j = jmax
!!!$omp parallel do
        do k = kmin, kmax - 1
           do i = imin, imax
              call weights2(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),RadiusDerivative(i,j,k),alpha(i,j,k),phi_old(i,j,k),&
                   & dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il&
                   &,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, Xcenter,&
                   & Ycenter, Zcenter,gxx_old(i,j,k),gxy_old(i,j,k),gxz_old(i,j,k),gyy_old(i,j,k),gyz_old(i,j,k),gzz_old(i,j,k))
              ! 
              ! Interpolate 
              !
              if (bc_type==WAVELIKE) then
                 call interp(ex,phi_old,f_int_phi,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,SYM,SYM,SYM)
              end if
              call interp(ex,K_old  ,f_int_K  ,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,gxx_old,f_int_gxx,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,gxy_old,f_int_gxy,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,sym_as,ANT,SYM)
              call interp(ex,gxz_old,f_int_gxz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,ANT,SYM,ANT)
              call interp(ex,gyy_old,f_int_gyy,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,gyz_old,f_int_gyz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,sym_sa,ANT,ANT)
              call interp(ex,gzz_old,f_int_gzz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,Axx_old,f_int_Axx,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,Axy_old,f_int_Axy,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,sym_as,ANT,SYM)
              call interp(ex,Axz_old,f_int_Axz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,ANT,SYM,ANT)
              call interp(ex,Ayy_old,f_int_Ayy,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,Ayz_old,f_int_Ayz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,sym_sa,ANT,ANT)
              call interp(ex,Azz_old,f_int_Azz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,RadiusDerivative,f_int_dR,il,iu,jl,ju,kl,ku,Dlx, &
                   Dly,Dlz,SYM,SYM,SYM)
              if(gambound .ne. 0) then
                 call interp(ex,Gammax_old,f_int_Gamx,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,ANT,SYM,SYM)
                 call interp(ex,Gammay_old,f_int_Gamy,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,SYM,ANT,SYM)
                 call interp(ex,Gammaz_old,f_int_Gamz,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,SYM,SYM,ANT)
              endif
              r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
              r1 = facr*r
              Rp = PhysicalR(i,j,k)
              Rp1 = facRp*Rp
              if (bc_type==WAVELIKE) then
                 phi_new(i,j,k) = facRp*f_int_phi +  &
                      f1o3*(log(Rp/r)-facRp*log(Rp1/r1)) + &
                      f1o6*(log(RadiusDerivative(i,j,k))-facRp*log(f_int_dR))
              end if
              K_new(i,j,k)   = f_int_K * facRp
              gxx_new(i,j,k) = ONE + facRp*(f_int_gxx-ONE)
              gxy_new(i,j,k) = facRp*f_int_gxy
              gxz_new(i,j,k) = facRp*f_int_gxz
              gyy_new(i,j,k) = ONE + facRp*(f_int_gyy-ONE)
              gyz_new(i,j,k) = facRp*f_int_gyz
              gzz_new(i,j,k) = ONE + facRp*(f_int_gzz-ONE)
              Axx_new(i,j,k) = f_int_Axx * facRp
              Axy_new(i,j,k) = f_int_Axy * facRp
              Axz_new(i,j,k) = f_int_Axz * facRp
              Ayy_new(i,j,k) = f_int_Ayy * facRp
              Ayz_new(i,j,k) = f_int_Ayz * facRp
              Azz_new(i,j,k) = f_int_Azz * facRp
              if(gambound.eq.0) then
                 Gammax_new(i,j,k) = Gammax_old(i,j,k)
                 Gammay_new(i,j,k) = Gammay_old(i,j,k)
                 Gammaz_new(i,j,k) = Gammaz_old(i,j,k)
              else
                 Gammax_new(i,j,k) = f_int_Gamx * facRp**(gambound)
                 Gammay_new(i,j,k) = f_int_Gamy * facRp**(gambound)
                 Gammaz_new(i,j,k) = f_int_Gamz * facRp**(gambound)
              endif
           enddo
        enddo
!!!$omp end parallel do
     endif
     !
     ! z = zmax
     !
     if (Symmetry.ne.AXISYM) then
        jub = jub + 1
     end if
     if(have_bdry_max(3)==1) then
        k = kmax
!!!$omp parallel do
        do j = jlb, jub
           do i = imin, imax
              call weights2(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),RadiusDerivative(i,j,k),alpha(i,j,k),phi_old(i,j,k),&
                   & dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il&
                   &,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry, Xcenter,&
                   & Ycenter, Zcenter,gxx_old(i,j,k),gxy_old(i,j,k),gxz_old(i,j,k),gyy_old(i,j,k),gyz_old(i,j,k),gzz_old(i,j,k))
              !
              ! Interpolate 
              !
              if (bc_type==WAVELIKE) then
                 call interp(ex,phi_old,f_int_phi,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,SYM,SYM,SYM)
              end if
              call interp(ex,K_old  ,f_int_K  ,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,gxx_old,f_int_gxx,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,gxy_old,f_int_gxy,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,sym_as,ANT,SYM)
              call interp(ex,gxz_old,f_int_gxz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,ANT,SYM,ANT)
              call interp(ex,gyy_old,f_int_gyy,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,gyz_old,f_int_gyz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,sym_sa,ANT,ANT)
              call interp(ex,gzz_old,f_int_gzz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,Axx_old,f_int_Axx,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,Axy_old,f_int_Axy,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,sym_as,ANT,SYM)
              call interp(ex,Axz_old,f_int_Axz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,ANT,SYM,ANT)
              call interp(ex,Ayy_old,f_int_Ayy,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,Ayz_old,f_int_Ayz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,sym_sa,ANT,ANT)
              call interp(ex,Azz_old,f_int_Azz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,RadiusDerivative,f_int_dR,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              if(gambound .ne. 0) then
                 call interp(ex,Gammax_old,f_int_Gamx,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,ANT,SYM,SYM)
                 call interp(ex,Gammay_old,f_int_Gamy,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,SYM,ANT,SYM)
                 call interp(ex,Gammaz_old,f_int_Gamz,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,SYM,SYM,ANT)
              endif
              r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
              r1 = facr*r
              Rp = PhysicalR(i,j,k)
              Rp1 = facRp*Rp
              psi4_int = exp(4.d0*f_int_phi)
              if (bc_type==WAVELIKE) then
                 phi_new(i,j,k) = facRp*f_int_phi +  &
                      f1o3*(log(Rp/r)-facRp*log(Rp1/r1)) + &
                      f1o6*(log(RadiusDerivative(i,j,k))-facRp*log(f_int_dR))
              end if
              psi4 = exp(4.d0*phi_new(i,j,k))
              fac = psi4_int/psi4*facRp
              K_new(i,j,k)   = f_int_K * facRp
              gxx_new(i,j,k) = ONE + facRp*(f_int_gxx-ONE)
              gxy_new(i,j,k) = facRp*f_int_gxy
              gxz_new(i,j,k) = facRp*f_int_gxz
              gyy_new(i,j,k) = ONE + facRp*(f_int_gyy-ONE)
              gyz_new(i,j,k) = facRp*f_int_gyz
              gzz_new(i,j,k) = ONE + facRp*(f_int_gzz-ONE)
              Axx_new(i,j,k) = f_int_Axx * facRp
              Axy_new(i,j,k) = f_int_Axy * facRp
              Axz_new(i,j,k) = f_int_Axz * facRp
              Ayy_new(i,j,k) = f_int_Ayy * facRp
              Ayz_new(i,j,k) = f_int_Ayz * facRp
              Azz_new(i,j,k) = f_int_Azz * facRp
              if(gambound.eq.0) then
                 Gammax_new(i,j,k) = Gammax_old(i,j,k)
                 Gammay_new(i,j,k) = Gammay_old(i,j,k)
                 Gammaz_new(i,j,k) = Gammaz_old(i,j,k)
              else
                 Gammax_new(i,j,k) = f_int_Gamx * facRp**(gambound)
                 Gammay_new(i,j,k) = f_int_Gamy * facRp**(gambound)
                 Gammaz_new(i,j,k) = f_int_Gamz * facRp**(gambound)
              endif
           enddo
        enddo
!!!$omp end parallel do
     end if
     if(have_bdry_min(1)==1) then
        !
        ! x = xmin
        !
        i = imin
!!!$omp parallel do
        do k = kmin, kmax-1
           do j = jmin, jmax-1
              call weights2(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),RadiusDerivative(i,j,k),alpha(i,j,k),phi_old(i,j,k)&
                   &, dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin&
                   &,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry,&
                   & Xcenter, Ycenter, Zcenter,gxx_old(i,j,k),gxy_old(i,j,k),gxz_old(i,j,k),gyy_old(i,j,k),gyz_old(i,j,k),gzz_old(i&
  &,j,k))
              !
              ! Interpolate 
              !
              if (bc_type==WAVELIKE) then
                 call interp(ex,phi_old,f_int_phi,il,iu,jl,ju,kl,ku&
                      &,Dlx,Dly,Dlz,SYM,SYM,SYM)
              end if
              call interp(ex,K_old  ,f_int_K  ,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,gxx_old,f_int_gxx,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,gxy_old,f_int_gxy,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,ANT,ANT,SYM)
              call interp(ex,gxz_old,f_int_gxz,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,ANT,SYM,ANT)
              call interp(ex,gyy_old,f_int_gyy,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,gyz_old,f_int_gyz,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,SYM,ANT,ANT)
              call interp(ex,gzz_old,f_int_gzz,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,Axx_old,f_int_Axx,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,Axy_old,f_int_Axy,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,ANT,ANT,SYM)
              call interp(ex,Axz_old,f_int_Axz,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,ANT,SYM,ANT)
              call interp(ex,Ayy_old,f_int_Ayy,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,Ayz_old,f_int_Ayz,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,SYM,ANT,ANT)
              call interp(ex,Azz_old,f_int_Azz,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,RadiusDerivative,f_int_dR,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              if(gambound .ne. 0) then
                 call interp(ex,Gammax_old,f_int_Gamx,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,ANT,SYM,SYM)
                 call interp(ex,Gammay_old,f_int_Gamy,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,SYM,ANT,SYM)
                 call interp(ex,Gammaz_old,f_int_Gamz,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,SYM,SYM,ANT)
              endif
              r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
              r1 = facr*r
              Rp = PhysicalR(i,j,k)
              Rp1 = facRp*Rp
              if (bc_type==WAVELIKE) then
                 phi_new(i,j,k) = facRp*f_int_phi +  &
                      f1o3*(log(Rp/r)-facRp*log(Rp1/r1)) + &
                      f1o6*(log(RadiusDerivative(i,j,k))-facRp*log(f_int_dR))
              end if
              K_new(i,j,k)   = f_int_K * facRp
              gxx_new(i,j,k) = ONE + facRp*(f_int_gxx-ONE)
              gxy_new(i,j,k) = facRp*f_int_gxy
              gxz_new(i,j,k) = facRp*f_int_gxz
              gyy_new(i,j,k) = ONE + facRp*(f_int_gyy-ONE)
              gyz_new(i,j,k) = facRp*f_int_gyz
              gzz_new(i,j,k) = ONE + facRp*(f_int_gzz-ONE)
              Axx_new(i,j,k) = f_int_Axx * facRp
              Axy_new(i,j,k) = f_int_Axy * facRp
              Axz_new(i,j,k) = f_int_Axz * facRp
              Ayy_new(i,j,k) = f_int_Ayy * facRp
              Ayz_new(i,j,k) = f_int_Ayz * facRp
              Azz_new(i,j,k) = f_int_Azz * facRp
              if(gambound.eq.0) then
                 Gammax_new(i,j,k) = Gammax_old(i,j,k)
                 Gammay_new(i,j,k) = Gammay_old(i,j,k)
                 Gammaz_new(i,j,k) = Gammaz_old(i,j,k)
              else
                 Gammax_new(i,j,k) = f_int_Gamx * facRp**(gambound)
                 Gammay_new(i,j,k) = f_int_Gamy * facRp**(gambound)
                 Gammaz_new(i,j,k) = f_int_Gamz * facRp**(gambound)
              endif
           enddo
        enddo
!!!$omp end parallel do
     endif
     if(have_bdry_min(2)==1) then
        !
        ! y = ymin
        !
        j = jmin
!!!$omp parallel do
        do k = kmin, kmax-1
           do i = imin +1, imax-1
              call weights2(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),RadiusDerivative(i,j,k),alpha(i,j,k),phi_old(i,j,k)&
                   &, dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin&
                   &,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry,&
                   & Xcenter, Ycenter, Zcenter,gxx_old(i,j,k),gxy_old(i,j,k),gxz_old(i,j,k),gyy_old(i,j,k),gyz_old(i,j,k),gzz_old(i&
  &,j,k))
              !
              ! Interpolate 
              !
              if (bc_type==WAVELIKE) then
                 call interp(ex,phi_old,f_int_phi,il,iu,jl,ju,kl,ku&
                      &,Dlx,Dly,Dlz,SYM,SYM,SYM)
              end if
              call interp(ex,K_old  ,f_int_K  ,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,gxx_old,f_int_gxx,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,gxy_old,f_int_gxy,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,ANT,ANT,SYM)
              call interp(ex,gxz_old,f_int_gxz,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,ANT,SYM,ANT)
              call interp(ex,gyy_old,f_int_gyy,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,gyz_old,f_int_gyz,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,SYM,ANT,ANT)
              call interp(ex,gzz_old,f_int_gzz,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,Axx_old,f_int_Axx,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,Axy_old,f_int_Axy,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,ANT,ANT,SYM)
              call interp(ex,Axz_old,f_int_Axz,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,ANT,SYM,ANT)
              call interp(ex,Ayy_old,f_int_Ayy,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,Ayz_old,f_int_Ayz,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,SYM,ANT,ANT)
              call interp(ex,Azz_old,f_int_Azz,il,iu,jl,ju,kl,ku&
                   &,Dlx,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,RadiusDerivative,f_int_dR,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              if(gambound .ne. 0) then
                 call interp(ex,Gammax_old,f_int_Gamx,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,ANT,SYM,SYM)
                 call interp(ex,Gammay_old,f_int_Gamy,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,SYM,ANT,SYM)
                 call interp(ex,Gammaz_old,f_int_Gamz,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,SYM,SYM,ANT)
              endif
              r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
              r1 = facr*r
              Rp = PhysicalR(i,j,k)
              Rp1 = facRp*Rp
              if (bc_type==WAVELIKE) then
                 phi_new(i,j,k) = facRp*f_int_phi +  &
                      f1o3*(log(Rp/r)-facRp*log(Rp1/r1)) + &
                      f1o6*(log(RadiusDerivative(i,j,k))-facRp*log(f_int_dR))
              end if
              K_new(i,j,k)   = f_int_K * facRp
              gxx_new(i,j,k) = ONE + facRp*(f_int_gxx-ONE)
              gxy_new(i,j,k) = facRp*f_int_gxy
              gxz_new(i,j,k) = facRp*f_int_gxz
              gyy_new(i,j,k) = ONE + facRp*(f_int_gyy-ONE)
              gyz_new(i,j,k) = facRp*f_int_gyz
              gzz_new(i,j,k) = ONE + facRp*(f_int_gzz-ONE)
              Axx_new(i,j,k) = f_int_Axx * facRp
              Axy_new(i,j,k) = f_int_Axy * facRp
              Axz_new(i,j,k) = f_int_Axz * facRp
              Ayy_new(i,j,k) = f_int_Ayy * facRp
              Ayz_new(i,j,k) = f_int_Ayz * facRp
              Azz_new(i,j,k) = f_int_Azz * facRp
              if(gambound.eq.0) then
                 Gammax_new(i,j,k) = Gammax_old(i,j,k)
                 Gammay_new(i,j,k) = Gammay_old(i,j,k)
                 Gammaz_new(i,j,k) = Gammaz_old(i,j,k)
              else
                 Gammax_new(i,j,k) = f_int_Gamx * facRp**(gambound)
                 Gammay_new(i,j,k) = f_int_Gamy * facRp**(gambound)
                 Gammaz_new(i,j,k) = f_int_Gamz * facRp**(gambound)
              endif
           enddo
        enddo
!!!$omp end parallel do
     endif
     if(have_bdry_min(3)==1) then
        !
        ! z = zmin
        !
        k = kmin
        do j = jmin+1, jmax-1
           do i = imin+1, imax-1
              call weights2(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),RadiusDerivative(i,j,k),alpha(i,j,k),phi_old(i,j,k)&
                   &, dX,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin&
                   &,kmax,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry,&
                   & Xcenter, Ycenter, Zcenter,gxx_old(i,j,k),gxy_old(i,j,k),gxz_old(i,j,k),gyy_old(i,j,k),gyz_old(i,j,k),gzz_old(i&
  &,j,k))
              !
              ! Interpolate 
              !
              if (bc_type==WAVELIKE) then
                 call interp(ex,phi_old,f_int_phi,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,SYM,SYM,SYM)
              end if
              call interp(ex,K_old  ,f_int_K  ,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,gxx_old,f_int_gxx,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,gxy_old,f_int_gxy,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,ANT,ANT,SYM)
              call interp(ex,gxz_old,f_int_gxz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,ANT,SYM,ANT)
              call interp(ex,gyy_old,f_int_gyy,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,gyz_old,f_int_gyz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,ANT,ANT)
              call interp(ex,gzz_old,f_int_gzz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,Axx_old,f_int_Axx,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,Axy_old,f_int_Axy,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,ANT,ANT,SYM)
              call interp(ex,Axz_old,f_int_Axz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,ANT,SYM,ANT)
              call interp(ex,Ayy_old,f_int_Ayy,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,Ayz_old,f_int_Ayz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,ANT,ANT)
              call interp(ex,Azz_old,f_int_Azz,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              call interp(ex,RadiusDerivative,f_int_dR,il,iu,jl,ju,kl,ku,Dlx&
                   &,Dly,Dlz,SYM,SYM,SYM)
              if(gambound .ne. 0) then
                 call interp(ex,Gammax_old,f_int_Gamx,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,ANT,SYM,SYM)
                 call interp(ex,Gammay_old,f_int_Gamy,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,SYM,ANT,SYM)
                 call interp(ex,Gammaz_old,f_int_Gamz,il,iu,jl,ju,kl,ku,Dlx&
                      &,Dly,Dlz,SYM,SYM,ANT)
              endif
              r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
              r1 = facr*r
              Rp = PhysicalR(i,j,k)
              Rp1 = facRp*Rp
              if (bc_type==WAVELIKE) then
                 phi_new(i,j,k) = facRp*f_int_phi +  &
                      f1o3*(log(Rp/r)-facRp*log(Rp1/r1)) + &
                      f1o6*(log(RadiusDerivative(i,j,k))-facRp*log(f_int_dR))
              end if
              K_new(i,j,k)   = f_int_K * facRp
              gxx_new(i,j,k) = ONE + facRp*(f_int_gxx-ONE)
              gxy_new(i,j,k) = facRp*f_int_gxy
              gxz_new(i,j,k) = facRp*f_int_gxz
              gyy_new(i,j,k) = ONE + facRp*(f_int_gyy-ONE)
              gyz_new(i,j,k) = facRp*f_int_gyz
              gzz_new(i,j,k) = ONE + facRp*(f_int_gzz-ONE)
              Axx_new(i,j,k) = f_int_Axx * facRp
              Axy_new(i,j,k) = f_int_Axy * facRp
              Axz_new(i,j,k) = f_int_Axz * facRp
              Ayy_new(i,j,k) = f_int_Ayy * facRp
              Ayz_new(i,j,k) = f_int_Ayz * facRp
              Azz_new(i,j,k) = f_int_Azz * facRp
              if(gambound.eq.0) then
                 Gammax_new(i,j,k) = Gammax_old(i,j,k)
                 Gammay_new(i,j,k) = Gammay_old(i,j,k)
                 Gammaz_new(i,j,k) = Gammaz_old(i,j,k)
              else
                 Gammax_new(i,j,k) = f_int_Gamx * facRp**(gambound)
                 Gammay_new(i,j,k) = f_int_Gamy * facRp**(gambound)
                 Gammaz_new(i,j,k) = f_int_Gamz * facRp**(gambound)
              endif
           enddo
        enddo
     endif
     ! Transform the components of gij and Aij back to fisheye coordinates 
     if (Fisheye_Flag==1) then
        fish_to_phys = .FALSE.
        call trans_fish_phys(ex, fake_ex, X, Y, Z, PhysicalR, RadiusDerivative, &
             gxx_old,gxy_old,gxz_old,gyy_old,gyz_old,gzz_old, &
             fish_to_phys,2,have_bdry_min,have_bdry_max,use_trans_fish_phys_new)
        call trans_fish_phys(ex, fake_ex, X, Y, Z, PhysicalR, RadiusDerivative, &
             Axx_old,Axy_old,Axz_old,Ayy_old,Ayz_old,Azz_old, &
             fish_to_phys,2,have_bdry_min,have_bdry_max,use_trans_fish_phys_new)
        call trans_fish_phys(ex, fake_ex, X, Y, Z, PhysicalR, RadiusDerivative, &
             gxx_new,gxy_new,gxz_new,gyy_new,gyz_new,gzz_new, &
             fish_to_phys,0,have_bdry_min,have_bdry_max,use_trans_fish_phys_new)
        call trans_fish_phys(ex, fake_ex, X, Y, Z, PhysicalR, RadiusDerivative, &
             Axx_new,Axy_new,Axz_new,Ayy_new,Ayz_new,Azz_new, &
             fish_to_phys,0,have_bdry_min,have_bdry_max,use_trans_fish_phys_new)
     end if
  else if(bc_type==FALLOFF) then
     !
     ! x = max
     !
     write(*,*) "BSSN BC_TYPE: FALLOFF"
     if(have_bdry_max(1)==1) then
        i = imax
        do k = kmin, kmax - 1
           do j = jmin, jmax - 1
              r  = sqrt(X(i,j,k)**2   + Y(i,j,k)**2 + Z(i,j,k)**2)
              r1 = sqrt(x(i-1,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
              fac = r1/r
              phi_new(i,j,k) = phi_new(i-1,j,k) * fac
              K_new(i,j,k)   = K_new(i-1,j,k)   * fac
              gxx_new(i,j,k) = (gxx_new(i-1,j,k) - ONE) * fac + ONE
              gxy_new(i,j,k) = gxy_new(i-1,j,k) * fac
              gxz_new(i,j,k) = gxz_new(i-1,j,k) * fac
              gyy_new(i,j,k) = (gyy_new(i-1,j,k) - ONE) * fac + ONE
              gyz_new(i,j,k) = gyz_new(i-1,j,k) * fac
              gzz_new(i,j,k) = (gzz_new(i-1,j,k) - ONE) * fac + ONE
              Axx_new(i,j,k) = Axx_new(i-1,j,k) * fac
              Axy_new(i,j,k) = Axy_new(i-1,j,k) * fac
              Axz_new(i,j,k) = Axz_new(i-1,j,k) * fac
              Ayy_new(i,j,k) = Ayy_new(i-1,j,k) * fac
              Ayz_new(i,j,k) = Ayz_new(i-1,j,k) * fac
              Azz_new(i,j,k) = Azz_new(i-1,j,k) * fac
              if(gambound.eq.0) then
                 Gammax_new(i,j,k) = Gammax_old(i,j,k)
                 Gammay_new(i,j,k) = Gammay_old(i,j,k)
                 Gammaz_new(i,j,k) = Gammaz_old(i,j,k)
              else
                 Gammax_new(i,j,k) = Gammax_old(i,j,k)*fac**gambound
                 Gammay_new(i,j,k) = Gammay_old(i,j,k)*fac**gambound
                 Gammaz_new(i,j,k) = Gammaz_old(i,j,k)*fac**gambound
              end if
           enddo
        enddo
     end if
     !
     ! y = ymax
     !
     if(have_bdry_max(2)==1) then
        j = jmax
        do k = kmin, kmax - 1
           do i = imin, imax
              r  = sqrt(X(i,j,k)**2 + Y(i,j,k)**2   + Z(i,j,k)**2)
              r1 = sqrt(X(i,j,k)**2 + y(i,j-1,k)**2 + Z(i,j,k)**2)
              fac = r1/r
              phi_new(i,j,k) = phi_new(i,j-1,k) * fac
              K_new(i,j,k)   = K_new(i,j-1,k)   * fac
              gxx_new(i,j,k) = (gxx_new(i,j-1,k) - ONE) * fac + ONE
              gxy_new(i,j,k) = gxy_new(i,j-1,k) * fac
              gxz_new(i,j,k) = gxz_new(i,j-1,k) * fac
              gyy_new(i,j,k) = (gyy_new(i,j-1,k) - ONE) * fac + ONE
              gyz_new(i,j,k) = gyz_new(i,j-1,k) * fac
              gzz_new(i,j,k) = (gzz_new(i,j-1,k) - ONE) * fac + ONE
              Axx_new(i,j,k) = Axx_new(i,j-1,k) * fac
              Axy_new(i,j,k) = Axy_new(i,j-1,k) * fac
              Axz_new(i,j,k) = Axz_new(i,j-1,k) * fac
              Ayy_new(i,j,k) = Ayy_new(i,j-1,k) * fac
              Ayz_new(i,j,k) = Ayz_new(i,j-1,k) * fac
              Azz_new(i,j,k) = Azz_new(i,j-1,k) * fac
              if(gambound.eq.0) then
                 Gammax_new(i,j,k) = Gammax_old(i,j,k)
                 Gammay_new(i,j,k) = Gammay_old(i,j,k)
                 Gammaz_new(i,j,k) = Gammaz_old(i,j,k)
              else
                 Gammax_new(i,j,k) = Gammax_old(i,j,k)*fac**gambound
                 Gammay_new(i,j,k) = Gammay_old(i,j,k)*fac**gambound
                 Gammaz_new(i,j,k) = Gammaz_old(i,j,k)*fac**gambound
              end if
           enddo
        enddo
     end if
     !
     ! z = zmax
     !
     if(have_bdry_max(3)==1) then
        k = kmax
        do j = jmin, jmax
           do i = imin, imax
              r  = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
              r1 = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + z(i,j,k-1)**2)
              fac = r1/r
              phi_new(i,j,k) = phi_new(i,j,k-1) * fac
              K_new(i,j,k)   = K_new(i,j,k-1)   * fac
              gxx_new(i,j,k) = (gxx_new(i,j,k-1) - ONE) * fac + ONE
              gxy_new(i,j,k) = gxy_new(i,j,k-1) * fac
              gxz_new(i,j,k) = gxz_new(i,j,k-1) * fac
              gyy_new(i,j,k) = (gyy_new(i,j,k-1) - ONE) * fac + ONE
              gyz_new(i,j,k) = gyz_new(i,j,k-1) * fac
              gzz_new(i,j,k) = (gzz_new(i,j,k-1) - ONE) * fac + ONE
              Axx_new(i,j,k) = Axx_new(i,j,k-1) * fac
              Axy_new(i,j,k) = Axy_new(i,j,k-1) * fac
              Axz_new(i,j,k) = Axz_new(i,j,k-1) * fac
              Ayy_new(i,j,k) = Ayy_new(i,j,k-1) * fac
              Ayz_new(i,j,k) = Ayz_new(i,j,k-1) * fac
              Azz_new(i,j,k) = Azz_new(i,j,k-1) * fac
              if(gambound.eq.0) then
                 Gammax_new(i,j,k) = Gammax_old(i,j,k)
                 Gammay_new(i,j,k) = Gammay_old(i,j,k)
                 Gammaz_new(i,j,k) = Gammaz_old(i,j,k)
              else
                 Gammax_new(i,j,k) = Gammax_old(i,j,k)*fac**gambound
                 Gammay_new(i,j,k) = Gammay_old(i,j,k)*fac**gambound
                 Gammaz_new(i,j,k) = Gammaz_old(i,j,k)*fac**gambound
              end if
           enddo
        enddo
     end if
     if(have_bdry_min(1)==1) then
        !
        ! x = xmin
        !
        i = imin
        do k = kmin, kmax-1
           do j = jmin, jmax-1
              r  = sqrt(X(i,j,k)**2   + Y(i,j,k)**2 + Z(i,j,k)**2)
              r1 = sqrt(x(i+1,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
              fac = r1/r
              phi_new(i,j,k) = phi_new(i+1,j,k) * fac
              K_new(i,j,k)   = K_new(i+1,j,k)   * fac
              gxx_new(i,j,k) = (gxx_new(i+1,j,k) - ONE) * fac + ONE
              gxy_new(i,j,k) = gxy_new(i+1,j,k) * fac
              gxz_new(i,j,k) = gxz_new(i+1,j,k) * fac
              gyy_new(i,j,k) = (gyy_new(i+1,j,k) - ONE) * fac + ONE
              gyz_new(i,j,k) = gyz_new(i+1,j,k) * fac
              gzz_new(i,j,k) = (gzz_new(i+1,j,k) - ONE) * fac + ONE
              Axx_new(i,j,k) = Axx_new(i+1,j,k) * fac
              Axy_new(i,j,k) = Axy_new(i+1,j,k) * fac
              Axz_new(i,j,k) = Axz_new(i+1,j,k) * fac
              Ayy_new(i,j,k) = Ayy_new(i+1,j,k) * fac
              Ayz_new(i,j,k) = Ayz_new(i+1,j,k) * fac
              Azz_new(i,j,k) = Azz_new(i+1,j,k) * fac
              if(gambound.eq.0) then
                 Gammax_new(i,j,k) = Gammax_old(i,j,k)
                 Gammay_new(i,j,k) = Gammay_old(i,j,k)
                 Gammaz_new(i,j,k) = Gammaz_old(i,j,k)
              else
                 Gammax_new(i,j,k) = Gammax_old(i,j,k)*fac**gambound
                 Gammay_new(i,j,k) = Gammay_old(i,j,k)*fac**gambound
                 Gammaz_new(i,j,k) = Gammaz_old(i,j,k)*fac**gambound
              end if
           enddo
        enddo
     endif
     if(have_bdry_min(2)==1) then
        !
        ! y = ymin
        !
        j = jmin
        do k = kmin, kmax-1
           do i = imin +1, imax-1
              r  = sqrt(X(i,j,k)**2 + Y(i,j,k)**2   + Z(i,j,k)**2)
              r1 = sqrt(X(i,j,k)**2 + y(i,j+1,k)**2 + Z(i,j,k)**2)
              fac = r1/r
              phi_new(i,j,k) = phi_new(i,j+1,k) * fac
              K_new(i,j,k)   = K_new(i,j+1,k)   * fac
              gxx_new(i,j,k) = (gxx_new(i,j+1,k) - ONE) * fac + ONE
              gxy_new(i,j,k) = gxy_new(i,j+1,k) * fac
              gxz_new(i,j,k) = gxz_new(i,j+1,k) * fac
              gyy_new(i,j,k) = (gyy_new(i,j+1,k) - ONE) * fac + ONE
              gyz_new(i,j,k) = gyz_new(i,j+1,k) * fac
              gzz_new(i,j,k) = (gzz_new(i,j+1,k) - ONE) * fac + ONE
              Axx_new(i,j,k) = Axx_new(i,j+1,k) * fac
              Axy_new(i,j,k) = Axy_new(i,j+1,k) * fac
              Axz_new(i,j,k) = Axz_new(i,j+1,k) * fac
              Ayy_new(i,j,k) = Ayy_new(i,j+1,k) * fac
              Ayz_new(i,j,k) = Ayz_new(i,j+1,k) * fac
              Azz_new(i,j,k) = Azz_new(i,j+1,k) * fac
              if(gambound.eq.0) then
                 Gammax_new(i,j,k) = Gammax_old(i,j,k)
                 Gammay_new(i,j,k) = Gammay_old(i,j,k)
                 Gammaz_new(i,j,k) = Gammaz_old(i,j,k)
              else
                 Gammax_new(i,j,k) = Gammax_old(i,j,k)*fac**gambound
                 Gammay_new(i,j,k) = Gammay_old(i,j,k)*fac**gambound
                 Gammaz_new(i,j,k) = Gammaz_old(i,j,k)*fac**gambound
              end if
           enddo
        enddo
     endif
     if(have_bdry_min(3)==1) then
        !
        ! z = zmin
        !
        k = kmin
        do j = jmin+1, jmax-1
           do i = imin+1, imax-1
              r  = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
              r1 = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + z(i,j,k+1)**2)
              fac = r1/r
              phi_new(i,j,k) = phi_new(i,j,k+1) * fac
              K_new(i,j,k)   = K_new(i,j,k+1)   * fac
              gxx_new(i,j,k) = (gxx_new(i,j,k+1) - ONE) * fac + ONE
              gxy_new(i,j,k) = gxy_new(i,j,k+1) * fac
              gxz_new(i,j,k) = gxz_new(i,j,k+1) * fac
              gyy_new(i,j,k) = (gyy_new(i,j,k+1) - ONE) * fac + ONE
              gyz_new(i,j,k) = gyz_new(i,j,k+1) * fac
              gzz_new(i,j,k) = (gzz_new(i,j,k+1) - ONE) * fac + ONE
              Axx_new(i,j,k) = Axx_new(i,j,k+1) * fac
              Axy_new(i,j,k) = Axy_new(i,j,k+1) * fac
              Axz_new(i,j,k) = Axz_new(i,j,k+1) * fac
              Ayy_new(i,j,k) = Ayy_new(i,j,k+1) * fac
              Ayz_new(i,j,k) = Ayz_new(i,j,k+1) * fac
              Azz_new(i,j,k) = Azz_new(i,j,k+1) * fac
              if(gambound.eq.0) then
                 Gammax_new(i,j,k) = Gammax_old(i,j,k)
                 Gammay_new(i,j,k) = Gammay_old(i,j,k)
                 Gammaz_new(i,j,k) = Gammaz_old(i,j,k)
              else
                 Gammax_new(i,j,k) = Gammax_old(i,j,k)*fac**gambound
                 Gammay_new(i,j,k) = Gammay_old(i,j,k)*fac**gambound
                 Gammaz_new(i,j,k) = Gammaz_old(i,j,k)*fac**gambound
              end if
           enddo
        enddo
     endif
  else if(bc_type==QE) then
        write(*,*)     "BSSN BC_TYPE: QE"
     !
     ! x = max
     !
     if(have_bdry_max(1)==1) then
        i = imax
        do k = kmin, kmax - 1
           do j = jmin, jmax - 1
              phi_new(i,j,k) = phi_old(i,j,k)
              K_new(i,j,k)   = K_old(i,j,k)
              gxx_new(i,j,k) = gxx_old(i,j,k)
              gxy_new(i,j,k) = gxy_old(i,j,k)
              gxz_new(i,j,k) = gxz_old(i,j,k)
              gyy_new(i,j,k) = gyy_old(i,j,k)
              gyz_new(i,j,k) = gyz_old(i,j,k)
              gzz_new(i,j,k) = gzz_old(i,j,k)
              Axx_new(i,j,k) = Axx_old(i,j,k)
              Axy_new(i,j,k) = Axy_old(i,j,k)
              Axz_new(i,j,k) = Axz_old(i,j,k)
              Ayy_new(i,j,k) = Ayy_old(i,j,k)
              Ayz_new(i,j,k) = Ayz_old(i,j,k)
              Azz_new(i,j,k) = Azz_old(i,j,k)
              Gammax_new(i,j,k) = Gammax_old(i,j,k)
              Gammay_new(i,j,k) = Gammay_old(i,j,k)
              Gammaz_new(i,j,k) = Gammaz_old(i,j,k)
           enddo
        enddo
     end if
     !
     ! y = ymax
     !
     if(have_bdry_max(2)==1) then
        j = jmax
        do k = kmin, kmax - 1
           do i = imin, imax
              phi_new(i,j,k) = phi_old(i,j,k)
              K_new(i,j,k)   = K_old(i,j,k)
              gxx_new(i,j,k) = gxx_old(i,j,k)
              gxy_new(i,j,k) = gxy_old(i,j,k)
              gxz_new(i,j,k) = gxz_old(i,j,k)
              gyy_new(i,j,k) = gyy_old(i,j,k)
              gyz_new(i,j,k) = gyz_old(i,j,k)
              gzz_new(i,j,k) = gzz_old(i,j,k)
              Axx_new(i,j,k) = Axx_old(i,j,k)
              Axy_new(i,j,k) = Axy_old(i,j,k)
              Axz_new(i,j,k) = Axz_old(i,j,k)
              Ayy_new(i,j,k) = Ayy_old(i,j,k)
              Ayz_new(i,j,k) = Ayz_old(i,j,k)
              Azz_new(i,j,k) = Azz_old(i,j,k)
              Gammax_new(i,j,k) = Gammax_old(i,j,k)
              Gammay_new(i,j,k) = Gammay_old(i,j,k)
              Gammaz_new(i,j,k) = Gammaz_old(i,j,k)
           enddo
        enddo
     end if
     !
     ! z = zmax
     !
     if(have_bdry_max(3)==1) then
        k = kmax
        do j = jmin, jmax
           do i = imin, imax
              phi_new(i,j,k) = phi_old(i,j,k)
              K_new(i,j,k)   = K_old(i,j,k)
              gxx_new(i,j,k) = gxx_old(i,j,k)
              gxy_new(i,j,k) = gxy_old(i,j,k)
              gxz_new(i,j,k) = gxz_old(i,j,k)
              gyy_new(i,j,k) = gyy_old(i,j,k)
              gyz_new(i,j,k) = gyz_old(i,j,k)
              gzz_new(i,j,k) = gzz_old(i,j,k)
              Axx_new(i,j,k) = Axx_old(i,j,k)
              Axy_new(i,j,k) = Axy_old(i,j,k)
              Axz_new(i,j,k) = Axz_old(i,j,k)
              Ayy_new(i,j,k) = Ayy_old(i,j,k)
              Ayz_new(i,j,k) = Ayz_old(i,j,k)
              Azz_new(i,j,k) = Azz_old(i,j,k)
              Gammax_new(i,j,k) = Gammax_old(i,j,k)
              Gammay_new(i,j,k) = Gammay_old(i,j,k)
              Gammaz_new(i,j,k) = Gammaz_old(i,j,k)
           enddo
        enddo
     end if
     if(have_bdry_min(1)==1) then
        !
        ! x = xmin
        !
        i = imin
        do k = kmin, kmax-1
           do j = jmin, jmax-1
              phi_new(i,j,k) = phi_old(i,j,k)
              K_new(i,j,k)   = K_old(i,j,k)
              gxx_new(i,j,k) = gxx_old(i,j,k)
              gxy_new(i,j,k) = gxy_old(i,j,k)
              gxz_new(i,j,k) = gxz_old(i,j,k)
              gyy_new(i,j,k) = gyy_old(i,j,k)
              gyz_new(i,j,k) = gyz_old(i,j,k)
              gzz_new(i,j,k) = gzz_old(i,j,k)
              Axx_new(i,j,k) = Axx_old(i,j,k)
              Axy_new(i,j,k) = Axy_old(i,j,k)
              Axz_new(i,j,k) = Axz_old(i,j,k)
              Ayy_new(i,j,k) = Ayy_old(i,j,k)
              Ayz_new(i,j,k) = Ayz_old(i,j,k)
              Azz_new(i,j,k) = Azz_old(i,j,k)
              Gammax_new(i,j,k) = Gammax_old(i,j,k)
              Gammay_new(i,j,k) = Gammay_old(i,j,k)
              Gammaz_new(i,j,k) = Gammaz_old(i,j,k)
           enddo
        enddo
     endif
     if(have_bdry_min(2)==1) then
        !
        ! y = ymin
        !
        j = jmin
        do k = kmin, kmax-1
           do i = imin +1, imax-1
              phi_new(i,j,k) = phi_old(i,j,k)
              K_new(i,j,k)   = K_old(i,j,k)
              gxx_new(i,j,k) = gxx_old(i,j,k)
              gxy_new(i,j,k) = gxy_old(i,j,k)
              gxz_new(i,j,k) = gxz_old(i,j,k)
              gyy_new(i,j,k) = gyy_old(i,j,k)
              gyz_new(i,j,k) = gyz_old(i,j,k)
              gzz_new(i,j,k) = gzz_old(i,j,k)
              Axx_new(i,j,k) = Axx_old(i,j,k)
              Axy_new(i,j,k) = Axy_old(i,j,k)
              Axz_new(i,j,k) = Axz_old(i,j,k)
              Ayy_new(i,j,k) = Ayy_old(i,j,k)
              Ayz_new(i,j,k) = Ayz_old(i,j,k)
              Azz_new(i,j,k) = Azz_old(i,j,k)
              Gammax_new(i,j,k) = Gammax_old(i,j,k)
              Gammay_new(i,j,k) = Gammay_old(i,j,k)
              Gammaz_new(i,j,k) = Gammaz_old(i,j,k)
           enddo
        enddo
     endif
     if(have_bdry_min(3)==1) then
        !
        ! z = zmin
        !
        k = kmin
        do j = jmin+1, jmax-1
           do i = imin+1, imax-1
              phi_new(i,j,k) = phi_old(i,j,k)
              K_new(i,j,k)   = K_old(i,j,k)
              gxx_new(i,j,k) = gxx_old(i,j,k)
              gxy_new(i,j,k) = gxy_old(i,j,k)
              gxz_new(i,j,k) = gxz_old(i,j,k)
              gyy_new(i,j,k) = gyy_old(i,j,k)
              gyz_new(i,j,k) = gyz_old(i,j,k)
              gzz_new(i,j,k) = gzz_old(i,j,k)
              Axx_new(i,j,k) = Axx_old(i,j,k)
              Axy_new(i,j,k) = Axy_old(i,j,k)
              Axz_new(i,j,k) = Axz_old(i,j,k)
              Ayy_new(i,j,k) = Ayy_old(i,j,k)
              Ayz_new(i,j,k) = Ayz_old(i,j,k)
              Azz_new(i,j,k) = Azz_old(i,j,k)
              Gammax_new(i,j,k) = Gammax_old(i,j,k)
              Gammay_new(i,j,k) = Gammay_old(i,j,k)
              Gammaz_new(i,j,k) = Gammaz_old(i,j,k)
           enddo
        enddo
     endif
  else if(bc_type==EXTRAP) then
        write(*,*)           "BSSN BC_TYPE: EXTRAP"
     !
     ! x = max
     !
     if(have_bdry_max(1)==1) then
        i = imax
        do k = kmin, kmax
           do j = jmin, jmax
              phi_new(i,j,k) = TWO*phi_new(i-1,j,k) - phi_new(i-2,j,k)
              K_new(i,j,k)   = TWO*K_new(i-1,j,k) - K_new(i-2,j,k)
              gxx_new(i,j,k) = TWO*gxx_new(i-1,j,k) - gxx_new(i-2,j,k)
              gxy_new(i,j,k) = TWO*gxy_new(i-1,j,k) - gxy_new(i-2,j,k)
              gxz_new(i,j,k) = TWO*gxz_new(i-1,j,k) - gxz_new(i-2,j,k)
              gyy_new(i,j,k) = TWO*gyy_new(i-1,j,k) - gyy_new(i-2,j,k)
              gyz_new(i,j,k) = TWO*gyz_new(i-1,j,k) - gyz_new(i-2,j,k)
              gzz_new(i,j,k) = TWO*gzz_new(i-1,j,k) - gzz_new(i-2,j,k)
              Axx_new(i,j,k) = TWO*Axx_new(i-1,j,k) - Axx_new(i-2,j,k)
              Axy_new(i,j,k) = TWO*Axy_new(i-1,j,k) - Axy_new(i-2,j,k)
              Axz_new(i,j,k) = TWO*Axz_new(i-1,j,k) - Axz_new(i-2,j,k)
              Ayy_new(i,j,k) = TWO*Ayy_new(i-1,j,k) - Ayy_new(i-2,j,k)
              Ayz_new(i,j,k) = TWO*Ayz_new(i-1,j,k) - Ayz_new(i-2,j,k)
              Azz_new(i,j,k) = TWO*Azz_new(i-1,j,k) - Azz_new(i-2,j,k)
              Gammax_new(i,j,k) = TWO*Gammax_new(i-1,j,k) - Gammax_new(i-2,j,k)
              Gammay_new(i,j,k) = TWO*Gammay_new(i-1,j,k) - Gammay_new(i-2,j,k)
              Gammaz_new(i,j,k) = TWO*Gammaz_new(i-1,j,k) - Gammaz_new(i-2,j,k)
           enddo
        enddo
     end if
     !
     ! y = ymax
     !
     if(have_bdry_max(2)==1) then
        j = jmax
        do k = kmin, kmax
           do i = imin, imax
              phi_new(i,j,k) = TWO*phi_new(i,j-1,k) - phi_new(i,j-2,k)
              K_new(i,j,k)   = TWO*K_new(i,j-1,k) - K_new(i,j-2,k)
              gxx_new(i,j,k) = TWO*gxx_new(i,j-1,k) - gxx_new(i,j-2,k)
              gxy_new(i,j,k) = TWO*gxy_new(i,j-1,k) - gxy_new(i,j-2,k)
              gxz_new(i,j,k) = TWO*gxz_new(i,j-1,k) - gxz_new(i,j-2,k)
              gyy_new(i,j,k) = TWO*gyy_new(i,j-1,k) - gyy_new(i,j-2,k)
              gyz_new(i,j,k) = TWO*gyz_new(i,j-1,k) - gyz_new(i,j-2,k)
              gzz_new(i,j,k) = TWO*gzz_new(i,j-1,k) - gzz_new(i,j-2,k)
              Axx_new(i,j,k) = TWO*Axx_new(i,j-1,k) - Axx_new(i,j-2,k)
              Axy_new(i,j,k) = TWO*Axy_new(i,j-1,k) - Axy_new(i,j-2,k)
              Axz_new(i,j,k) = TWO*Axz_new(i,j-1,k) - Axz_new(i,j-2,k)
              Ayy_new(i,j,k) = TWO*Ayy_new(i,j-1,k) - Ayy_new(i,j-2,k)
              Ayz_new(i,j,k) = TWO*Ayz_new(i,j-1,k) - Ayz_new(i,j-2,k)
              Azz_new(i,j,k) = TWO*Azz_new(i,j-1,k) - Azz_new(i,j-2,k)
              Gammax_new(i,j,k) = TWO*Gammax_new(i,j-1,k) - Gammax_new(i,j-2,k)
              Gammay_new(i,j,k) = TWO*Gammay_new(i,j-1,k) - Gammay_new(i,j-2,k)
              Gammaz_new(i,j,k) = TWO*Gammaz_new(i,j-1,k) - Gammaz_new(i,j-2,k)
           enddo
        enddo
     end if
     !
     ! z = zmax
     !
     if(have_bdry_max(3)==1) then
        k = kmax
        do j = jmin, jmax
           do i = imin, imax
              phi_new(i,j,k) = TWO*phi_new(i,j,k-1) - phi_new(i,j,k-2)
              K_new(i,j,k)   = TWO*K_new(i,j,k-1) - K_new(i,j,k-2)
              gxx_new(i,j,k) = TWO*gxx_new(i,j,k-1) - gxx_new(i,j,k-2)
              gxy_new(i,j,k) = TWO*gxy_new(i,j,k-1) - gxy_new(i,j,k-2)
              gxz_new(i,j,k) = TWO*gxz_new(i,j,k-1) - gxz_new(i,j,k-2)
              gyy_new(i,j,k) = TWO*gyy_new(i,j,k-1) - gyy_new(i,j,k-2)
              gyz_new(i,j,k) = TWO*gyz_new(i,j,k-1) - gyz_new(i,j,k-2)
              gzz_new(i,j,k) = TWO*gzz_new(i,j,k-1) - gzz_new(i,j,k-2)
              Axx_new(i,j,k) = TWO*Axx_new(i,j,k-1) - Axx_new(i,j,k-2)
              Axy_new(i,j,k) = TWO*Axy_new(i,j,k-1) - Axy_new(i,j,k-2)
              Axz_new(i,j,k) = TWO*Axz_new(i,j,k-1) - Axz_new(i,j,k-2)
              Ayy_new(i,j,k) = TWO*Ayy_new(i,j,k-1) - Ayy_new(i,j,k-2)
              Ayz_new(i,j,k) = TWO*Ayz_new(i,j,k-1) - Ayz_new(i,j,k-2)
              Azz_new(i,j,k) = TWO*Azz_new(i,j,k-1) - Azz_new(i,j,k-2)
              Gammax_new(i,j,k) = TWO*Gammax_new(i,j,k-1) - Gammax_new(i,j,k-2)
              Gammay_new(i,j,k) = TWO*Gammay_new(i,j,k-1) - Gammay_new(i,j,k-2)
              Gammaz_new(i,j,k) = TWO*Gammaz_new(i,j,k-1) - Gammaz_new(i,j,k-2)
           enddo
        enddo
     end if
     if(have_bdry_min(1)==1) then
        !
        ! x = xmin
        !
        i = imin
        do k = kmin, kmax
           do j = jmin, jmax
              phi_new(i,j,k) = TWO*phi_new(imin+1,j,k) - phi_new(imin+2,j,k)
              K_new(i,j,k) = TWO*K_new(imin+1,j,k) - K_new(imin+2,j,k)
              gxx_new(i,j,k) = TWO*gxx_new(imin+1,j,k) - gxx_new(imin+2,j,k)
              gxy_new(i,j,k) = TWO*gxy_new(imin+1,j,k) - gxy_new(imin+2,j,k)
              gxz_new(i,j,k) = TWO*gxz_new(imin+1,j,k) - gxz_new(imin+2,j,k)
              gyy_new(i,j,k) = TWO*gyy_new(imin+1,j,k) - gyy_new(imin+2,j,k)
              gyz_new(i,j,k) = TWO*gyz_new(imin+1,j,k) - gyz_new(imin+2,j,k)
              gzz_new(i,j,k) = TWO*gzz_new(imin+1,j,k) - gzz_new(imin+2,j,k)
              Axx_new(i,j,k) = TWO*Axx_new(imin+1,j,k) - Axx_new(imin+2,j,k)
              Axy_new(i,j,k) = TWO*Axy_new(imin+1,j,k) - Axy_new(imin+2,j,k)
              Axz_new(i,j,k) = TWO*Axz_new(imin+1,j,k) - Axz_new(imin+2,j,k)
              Ayy_new(i,j,k) = TWO*Ayy_new(imin+1,j,k) - Ayy_new(imin+2,j,k)
              Ayz_new(i,j,k) = TWO*Ayz_new(imin+1,j,k) - Ayz_new(imin+2,j,k)
              Azz_new(i,j,k) = TWO*Azz_new(imin+1,j,k) - Azz_new(imin+2,j,k)
              Gammax_new(i,j,k) = TWO*Gammax_new(imin+1,j,k) - Gammax_new(imin+2,j,k)
              Gammay_new(i,j,k) = TWO*Gammay_new(imin+1,j,k) - Gammay_new(imin+2,j,k)
              Gammaz_new(i,j,k) = TWO*Gammaz_new(imin+1,j,k) - Gammaz_new(imin+2,j,k)
           enddo
        enddo
     endif
     if(have_bdry_min(2)==1) then
        !
        ! y = ymin
        !
        j = jmin
        do k = kmin, kmax
           do i = imin, imax
              phi_new(i,j,k) = TWO*phi_new(i,jmin+1,k) - phi_new(i,jmin+2,k)
              K_new(i,j,k) = TWO*K_new(i,jmin+1,k) - K_new(i,jmin+2,k)
              gxx_new(i,j,k) = TWO*gxx_new(i,jmin+1,k) - gxx_new(i,jmin+2,k)
              gxy_new(i,j,k) = TWO*gxy_new(i,jmin+1,k) - gxy_new(i,jmin+2,k)
              gxz_new(i,j,k) = TWO*gxz_new(i,jmin+1,k) - gxz_new(i,jmin+2,k)
              gyy_new(i,j,k) = TWO*gyy_new(i,jmin+1,k) - gyy_new(i,jmin+2,k)
              gyz_new(i,j,k) = TWO*gyz_new(i,jmin+1,k) - gyz_new(i,jmin+2,k)
              gzz_new(i,j,k) = TWO*gzz_new(i,jmin+1,k) - gzz_new(i,jmin+2,k)
              Axx_new(i,j,k) = TWO*Axx_new(i,jmin+1,k) - Axx_new(i,jmin+2,k)
              Axy_new(i,j,k) = TWO*Axy_new(i,jmin+1,k) - Axy_new(i,jmin+2,k)
              Axz_new(i,j,k) = TWO*Axz_new(i,jmin+1,k) - Axz_new(i,jmin+2,k)
              Ayy_new(i,j,k) = TWO*Ayy_new(i,jmin+1,k) - Ayy_new(i,jmin+2,k)
              Ayz_new(i,j,k) = TWO*Ayz_new(i,jmin+1,k) - Ayz_new(i,jmin+2,k)
              Azz_new(i,j,k) = TWO*Azz_new(i,jmin+1,k) - Azz_new(i,jmin+2,k)
              Gammax_new(i,j,k) = TWO*Gammax_new(i,jmin+1,k) - Gammax_new(i,jmin+2,k)
              Gammay_new(i,j,k) = TWO*Gammay_new(i,jmin+1,k) - Gammay_new(i,jmin+2,k)
              Gammaz_new(i,j,k) = TWO*Gammaz_new(i,jmin+1,k) - Gammaz_new(i,jmin+2,k)
           enddo
        enddo
     endif
     if(have_bdry_min(3)==1) then
        !
        ! z = zmin
        !
        k = kmin
        do j = jmin, jmax
           do i = imin, imax
              phi_new(i,j,k) = TWO*phi_new(i,j,kmin+1) - phi_new(i,j,kmin+2)
              K_new(i,j,k) = TWO*K_new(i,j,kmin+1) - K_new(i,j,kmin+2)
              gxx_new(i,j,k) = TWO*gxx_new(i,j,kmin+1) - gxx_new(i,j,kmin+2)
              gxy_new(i,j,k) = TWO*gxy_new(i,j,kmin+1) - gxy_new(i,j,kmin+2)
              gxz_new(i,j,k) = TWO*gxz_new(i,j,kmin+1) - gxz_new(i,j,kmin+2)
              gyy_new(i,j,k) = TWO*gyy_new(i,j,kmin+1) - gyy_new(i,j,kmin+2)
              gyz_new(i,j,k) = TWO*gyz_new(i,j,kmin+1) - gyz_new(i,j,kmin+2)
              gzz_new(i,j,k) = TWO*gzz_new(i,j,kmin+1) - gzz_new(i,j,kmin+2)
              Axx_new(i,j,k) = TWO*Axx_new(i,j,kmin+1) - Axx_new(i,j,kmin+2)
              Axy_new(i,j,k) = TWO*Axy_new(i,j,kmin+1) - Axy_new(i,j,kmin+2)
              Axz_new(i,j,k) = TWO*Axz_new(i,j,kmin+1) - Axz_new(i,j,kmin+2)
              Ayy_new(i,j,k) = TWO*Ayy_new(i,j,kmin+1) - Ayy_new(i,j,kmin+2)
              Ayz_new(i,j,k) = TWO*Ayz_new(i,j,kmin+1) - Ayz_new(i,j,kmin+2)
              Azz_new(i,j,k) = TWO*Azz_new(i,j,kmin+1) - Azz_new(i,j,kmin+2)
              Gammax_new(i,j,k) = TWO*Gammax_new(i,j,kmin+1) - Gammax_new(i,j,kmin+2)
              Gammay_new(i,j,k) = TWO*Gammay_new(i,j,kmin+1) - Gammay_new(i,j,kmin+2)
              Gammaz_new(i,j,k) = TWO*Gammaz_new(i,j,kmin+1) - Gammaz_new(i,j,kmin+2)
           enddo
        enddo
     endif
  else if(bc_type==COPY) then
        write(*,*)         "BSSN BC_TYPE: COPY"
     !
     ! x = max
     !
     if(have_bdry_max(1)==1) then
        i = imax
        do k = kmin, kmax
           do j = jmin, jmax
              phi_new(i,j,k) = phi_new(i-1,j,k)
              K_new(i,j,k)   = K_new(i-1,j,k)
              gxx_new(i,j,k) = gxx_new(i-1,j,k)
              gxy_new(i,j,k) = gxy_new(i-1,j,k)
              gxz_new(i,j,k) = gxz_new(i-1,j,k)
              gyy_new(i,j,k) = gyy_new(i-1,j,k)
              gyz_new(i,j,k) = gyz_new(i-1,j,k)
              gzz_new(i,j,k) = gzz_new(i-1,j,k)
              Axx_new(i,j,k) = Axx_new(i-1,j,k)
              Axy_new(i,j,k) = Axy_new(i-1,j,k)
              Axz_new(i,j,k) = Axz_new(i-1,j,k)
              Ayy_new(i,j,k) = Ayy_new(i-1,j,k)
              Ayz_new(i,j,k) = Ayz_new(i-1,j,k)
              Azz_new(i,j,k) = Azz_new(i-1,j,k)
              Gammax_new(i,j,k) = Gammax_new(i-1,j,k)
              Gammay_new(i,j,k) = Gammay_new(i-1,j,k)
              Gammaz_new(i,j,k) = Gammaz_new(i-1,j,k)
           enddo
        enddo
     end if
     !
     ! y = ymax
     !
     if(have_bdry_max(2)==1) then
        j = jmax
        do k = kmin, kmax
           do i = imin, imax
              phi_new(i,j,k) = phi_new(i,j-1,k)
              K_new(i,j,k)   = K_new(i,j-1,k)
              gxx_new(i,j,k) = gxx_new(i,j-1,k)
              gxy_new(i,j,k) = gxy_new(i,j-1,k)
              gxz_new(i,j,k) = gxz_new(i,j-1,k)
              gyy_new(i,j,k) = gyy_new(i,j-1,k)
              gyz_new(i,j,k) = gyz_new(i,j-1,k)
              gzz_new(i,j,k) = gzz_new(i,j-1,k)
              Axx_new(i,j,k) = Axx_new(i,j-1,k)
              Axy_new(i,j,k) = Axy_new(i,j-1,k)
              Axz_new(i,j,k) = Axz_new(i,j-1,k)
              Ayy_new(i,j,k) = Ayy_new(i,j-1,k)
              Ayz_new(i,j,k) = Ayz_new(i,j-1,k)
              Azz_new(i,j,k) = Azz_new(i,j-1,k)
              Gammax_new(i,j,k) = Gammax_new(i,j-1,k)
              Gammay_new(i,j,k) = Gammay_new(i,j-1,k)
              Gammaz_new(i,j,k) = Gammaz_new(i,j-1,k)
           enddo
        enddo
     end if
     !
     ! z = zmax
     !
     if(have_bdry_max(3)==1) then
        k = kmax
        do j = jmin, jmax
           do i = imin, imax
              phi_new(i,j,k) = phi_new(i,j,k-1)
              K_new(i,j,k)   = K_new(i,j,k-1)
              gxx_new(i,j,k) = gxx_new(i,j,k-1)
              gxy_new(i,j,k) = gxy_new(i,j,k-1)
              gxz_new(i,j,k) = gxz_new(i,j,k-1)
              gyy_new(i,j,k) = gyy_new(i,j,k-1)
              gyz_new(i,j,k) = gyz_new(i,j,k-1)
              gzz_new(i,j,k) = gzz_new(i,j,k-1)
              Axx_new(i,j,k) = Axx_new(i,j,k-1)
              Axy_new(i,j,k) = Axy_new(i,j,k-1)
              Axz_new(i,j,k) = Axz_new(i,j,k-1)
              Ayy_new(i,j,k) = Ayy_new(i,j,k-1)
              Ayz_new(i,j,k) = Ayz_new(i,j,k-1)
              Azz_new(i,j,k) = Azz_new(i,j,k-1)
              Gammax_new(i,j,k) = Gammax_new(i,j,k-1)
              Gammay_new(i,j,k) = Gammay_new(i,j,k-1)
              Gammaz_new(i,j,k) = Gammaz_new(i,j,k-1)
           enddo
        enddo
     end if
     if(have_bdry_min(1)==1) then
        !
        ! x = xmin
        !
        i = imin
        do k = kmin, kmax
           do j = jmin, jmax
              phi_new(i,j,k) = phi_new(imin+1,j,k)
              K_new(i,j,k) = K_new(imin+1,j,k)
              gxx_new(i,j,k) = gxx_new(imin+1,j,k)
              gxy_new(i,j,k) = gxy_new(imin+1,j,k)
              gxz_new(i,j,k) = gxz_new(imin+1,j,k)
              gyy_new(i,j,k) = gyy_new(imin+1,j,k)
              gyz_new(i,j,k) = gyz_new(imin+1,j,k)
              gzz_new(i,j,k) = gzz_new(imin+1,j,k)
              Axx_new(i,j,k) = Axx_new(imin+1,j,k)
              Axy_new(i,j,k) = Axy_new(imin+1,j,k)
              Axz_new(i,j,k) = Axz_new(imin+1,j,k)
              Ayy_new(i,j,k) = Ayy_new(imin+1,j,k)
              Ayz_new(i,j,k) = Ayz_new(imin+1,j,k)
              Azz_new(i,j,k) = Azz_new(imin+1,j,k)
              Gammax_new(i,j,k) = Gammax_new(imin+1,j,k)
              Gammay_new(i,j,k) = Gammay_new(imin+1,j,k)
              Gammaz_new(i,j,k) = Gammaz_new(imin+1,j,k)
           enddo
        enddo
     endif
     if(have_bdry_min(2)==1) then
        !
        ! y = ymin
        !
        j = jmin
        do k = kmin, kmax
           do i = imin, imax
              phi_new(i,j,k) = phi_new(i,jmin+1,k)
              K_new(i,j,k) = K_new(i,jmin+1,k)
              gxx_new(i,j,k) = gxx_new(i,jmin+1,k)
              gxy_new(i,j,k) = gxy_new(i,jmin+1,k)
              gxz_new(i,j,k) = gxz_new(i,jmin+1,k)
              gyy_new(i,j,k) = gyy_new(i,jmin+1,k)
              gyz_new(i,j,k) = gyz_new(i,jmin+1,k)
              gzz_new(i,j,k) = gzz_new(i,jmin+1,k)
              Axx_new(i,j,k) = Axx_new(i,jmin+1,k)
              Axy_new(i,j,k) = Axy_new(i,jmin+1,k)
              Axz_new(i,j,k) = Axz_new(i,jmin+1,k)
              Ayy_new(i,j,k) = Ayy_new(i,jmin+1,k)
              Ayz_new(i,j,k) = Ayz_new(i,jmin+1,k)
              Azz_new(i,j,k) = Azz_new(i,jmin+1,k)
              Gammax_new(i,j,k) = Gammax_new(i,jmin+1,k)
              Gammay_new(i,j,k) = Gammay_new(i,jmin+1,k)
              Gammaz_new(i,j,k) = Gammaz_new(i,jmin+1,k)
           enddo
        enddo
     endif
     if(have_bdry_min(3)==1) then
        !
        ! z = zmin
        !
        k = kmin
        do j = jmin, jmax
           do i = imin, imax
              phi_new(i,j,k) = phi_new(i,j,kmin+1)
              K_new(i,j,k) = K_new(i,j,kmin+1)
              gxx_new(i,j,k) = gxx_new(i,j,kmin+1)
              gxy_new(i,j,k) = gxy_new(i,j,kmin+1)
              gxz_new(i,j,k) = gxz_new(i,j,kmin+1)
              gyy_new(i,j,k) = gyy_new(i,j,kmin+1)
              gyz_new(i,j,k) = gyz_new(i,j,kmin+1)
              gzz_new(i,j,k) = gzz_new(i,j,kmin+1)
              Axx_new(i,j,k) = Axx_new(i,j,kmin+1)
              Axy_new(i,j,k) = Axy_new(i,j,kmin+1)
              Axz_new(i,j,k) = Axz_new(i,j,kmin+1)
              Ayy_new(i,j,k) = Ayy_new(i,j,kmin+1)
              Ayz_new(i,j,k) = Ayz_new(i,j,kmin+1)
              Azz_new(i,j,k) = Azz_new(i,j,kmin+1)
              Gammax_new(i,j,k) = Gammax_new(i,j,kmin+1)
              Gammay_new(i,j,k) = Gammay_new(i,j,kmin+1)
              Gammaz_new(i,j,k) = Gammaz_new(i,j,kmin+1)
           enddo
        enddo
     endif
  else if(bc_type==QUAD) then
        write(*,*)         "BSSN BC_TYPE: QUAD"
     !
     ! x = max
     !
     if(have_bdry_max(1)==1) then
        phi_new(imax,:,:) = THREE*phi_new(imax-1,:,:) - THREE*phi_new(imax-2,:,:) + phi_new(imax-3,:,:)
        K_new(imax,:,:)   = THREE*K_new(imax-1,:,:)   - THREE*K_new(imax-2,:,:)   + K_new(imax-3,:,:)
        gxx_new(imax,:,:) = THREE*gxx_new(imax-1,:,:) - THREE*gxx_new(imax-2,:,:) + gxx_new(imax-3,:,:)
        gxy_new(imax,:,:) = THREE*gxy_new(imax-1,:,:) - THREE*gxy_new(imax-2,:,:) + gxy_new(imax-3,:,:)
        gxz_new(imax,:,:) = THREE*gxz_new(imax-1,:,:) - THREE*gxz_new(imax-2,:,:) + gxz_new(imax-3,:,:)
        gyy_new(imax,:,:) = THREE*gyy_new(imax-1,:,:) - THREE*gyy_new(imax-2,:,:) + gyy_new(imax-3,:,:)
        gyz_new(imax,:,:) = THREE*gyz_new(imax-1,:,:) - THREE*gyz_new(imax-2,:,:) + gyz_new(imax-3,:,:)
        gzz_new(imax,:,:) = THREE*gzz_new(imax-1,:,:) - THREE*gzz_new(imax-2,:,:) + gzz_new(imax-3,:,:)
        Axx_new(imax,:,:) = THREE*Axx_new(imax-1,:,:) - THREE*Axx_new(imax-2,:,:) + Axx_new(imax-3,:,:)
        Axy_new(imax,:,:) = THREE*Axy_new(imax-1,:,:) - THREE*Axy_new(imax-2,:,:) + Axy_new(imax-3,:,:)
        Axz_new(imax,:,:) = THREE*Axz_new(imax-1,:,:) - THREE*Axz_new(imax-2,:,:) + Axz_new(imax-3,:,:)
        Ayy_new(imax,:,:) = THREE*Ayy_new(imax-1,:,:) - THREE*Ayy_new(imax-2,:,:) + Ayy_new(imax-3,:,:)
        Ayz_new(imax,:,:) = THREE*Ayz_new(imax-1,:,:) - THREE*Ayz_new(imax-2,:,:) + Ayz_new(imax-3,:,:)
        Azz_new(imax,:,:) = THREE*Azz_new(imax-1,:,:) - THREE*Azz_new(imax-2,:,:) + Azz_new(imax-3,:,:)
        Gammax_new(imax,:,:) = THREE*Gammax_new(imax-1,:,:) - THREE*Gammax_new(imax-2,:,:) + Gammax_new(imax-3,:,:)
        Gammay_new(imax,:,:) = THREE*Gammay_new(imax-1,:,:) - THREE*Gammay_new(imax-2,:,:) + Gammay_new(imax-3,:,:)
        Gammaz_new(imax,:,:) = THREE*Gammaz_new(imax-1,:,:) - THREE*Gammaz_new(imax-2,:,:) + Gammaz_new(imax-3,:,:)
     end if
     !
     ! y = ymax
     !
     if(have_bdry_max(2)==1) then
        phi_new(:,jmax,:) = THREE*phi_new(:,jmax-1,:) - THREE*phi_new(:,jmax-2,:) + phi_new(:,jmax-3,:)
        K_new(:,jmax,:)   = THREE*K_new(:,jmax-1,:)   - THREE*K_new(:,jmax-2,:)   + K_new(:,jmax-3,:)
        gxx_new(:,jmax,:) = THREE*gxx_new(:,jmax-1,:) - THREE*gxx_new(:,jmax-2,:) + gxx_new(:,jmax-3,:)
        gxy_new(:,jmax,:) = THREE*gxy_new(:,jmax-1,:) - THREE*gxy_new(:,jmax-2,:) + gxy_new(:,jmax-3,:)
        gxz_new(:,jmax,:) = THREE*gxz_new(:,jmax-1,:) - THREE*gxz_new(:,jmax-2,:) + gxz_new(:,jmax-3,:)
        gyy_new(:,jmax,:) = THREE*gyy_new(:,jmax-1,:) - THREE*gyy_new(:,jmax-2,:) + gyy_new(:,jmax-3,:)
        gyz_new(:,jmax,:) = THREE*gyz_new(:,jmax-1,:) - THREE*gyz_new(:,jmax-2,:) + gyz_new(:,jmax-3,:)
        gzz_new(:,jmax,:) = THREE*gzz_new(:,jmax-1,:) - THREE*gzz_new(:,jmax-2,:) + gzz_new(:,jmax-3,:)
        Axx_new(:,jmax,:) = THREE*Axx_new(:,jmax-1,:) - THREE*Axx_new(:,jmax-2,:) + Axx_new(:,jmax-3,:)
        Axy_new(:,jmax,:) = THREE*Axy_new(:,jmax-1,:) - THREE*Axy_new(:,jmax-2,:) + Axy_new(:,jmax-3,:)
        Axz_new(:,jmax,:) = THREE*Axz_new(:,jmax-1,:) - THREE*Axz_new(:,jmax-2,:) + Axz_new(:,jmax-3,:)
        Ayy_new(:,jmax,:) = THREE*Ayy_new(:,jmax-1,:) - THREE*Ayy_new(:,jmax-2,:) + Ayy_new(:,jmax-3,:)
        Ayz_new(:,jmax,:) = THREE*Ayz_new(:,jmax-1,:) - THREE*Ayz_new(:,jmax-2,:) + Ayz_new(:,jmax-3,:)
        Azz_new(:,jmax,:) = THREE*Azz_new(:,jmax-1,:) - THREE*Azz_new(:,jmax-2,:) + Azz_new(:,jmax-3,:)
        Gammax_new(:,jmax,:) = THREE*Gammax_new(:,jmax-1,:) - THREE*Gammax_new(:,jmax-2,:) + Gammax_new(:,jmax-3,:)
        Gammay_new(:,jmax,:) = THREE*Gammay_new(:,jmax-1,:) - THREE*Gammay_new(:,jmax-2,:) + Gammay_new(:,jmax-3,:)
        Gammaz_new(:,jmax,:) = THREE*Gammaz_new(:,jmax-1,:) - THREE*Gammaz_new(:,jmax-2,:) + Gammaz_new(:,jmax-3,:)
     end if
     !
     ! z = zmax
     !
     if(have_bdry_max(3)==1) then
        phi_new(:,:,kmax) = THREE*phi_new(:,:,kmax-1) - THREE*phi_new(:,:,kmax-2) + phi_new(:,:,kmax-3)
        K_new(:,:,kmax)   = THREE*K_new(:,:,kmax-1)   - THREE*K_new(:,:,kmax-2)   + K_new(:,:,kmax-3)
        gxx_new(:,:,kmax) = THREE*gxx_new(:,:,kmax-1) - THREE*gxx_new(:,:,kmax-2) + gxx_new(:,:,kmax-3)
        gxy_new(:,:,kmax) = THREE*gxy_new(:,:,kmax-1) - THREE*gxy_new(:,:,kmax-2) + gxy_new(:,:,kmax-3)
        gxz_new(:,:,kmax) = THREE*gxz_new(:,:,kmax-1) - THREE*gxz_new(:,:,kmax-2) + gxz_new(:,:,kmax-3)
        gyy_new(:,:,kmax) = THREE*gyy_new(:,:,kmax-1) - THREE*gyy_new(:,:,kmax-2) + gyy_new(:,:,kmax-3)
        gyz_new(:,:,kmax) = THREE*gyz_new(:,:,kmax-1) - THREE*gyz_new(:,:,kmax-2) + gyz_new(:,:,kmax-3)
        gzz_new(:,:,kmax) = THREE*gzz_new(:,:,kmax-1) - THREE*gzz_new(:,:,kmax-2) + gzz_new(:,:,kmax-3)
        Axx_new(:,:,kmax) = THREE*Axx_new(:,:,kmax-1) - THREE*Axx_new(:,:,kmax-2) + Axx_new(:,:,kmax-3)
        Axy_new(:,:,kmax) = THREE*Axy_new(:,:,kmax-1) - THREE*Axy_new(:,:,kmax-2) + Axy_new(:,:,kmax-3)
        Axz_new(:,:,kmax) = THREE*Axz_new(:,:,kmax-1) - THREE*Axz_new(:,:,kmax-2) + Axz_new(:,:,kmax-3)
        Ayy_new(:,:,kmax) = THREE*Ayy_new(:,:,kmax-1) - THREE*Ayy_new(:,:,kmax-2) + Ayy_new(:,:,kmax-3)
        Ayz_new(:,:,kmax) = THREE*Ayz_new(:,:,kmax-1) - THREE*Ayz_new(:,:,kmax-2) + Ayz_new(:,:,kmax-3)
        Azz_new(:,:,kmax) = THREE*Azz_new(:,:,kmax-1) - THREE*Azz_new(:,:,kmax-2) + Azz_new(:,:,kmax-3)
        Gammax_new(:,:,kmax) = THREE*Gammax_new(:,:,kmax-1) - THREE*Gammax_new(:,:,kmax-2) + Gammax_new(:,:,kmax-3)
        Gammay_new(:,:,kmax) = THREE*Gammay_new(:,:,kmax-1) - THREE*Gammay_new(:,:,kmax-2) + Gammay_new(:,:,kmax-3)
        Gammaz_new(:,:,kmax) = THREE*Gammaz_new(:,:,kmax-1) - THREE*Gammaz_new(:,:,kmax-2) + Gammaz_new(:,:,kmax-3)
     end if
     if(have_bdry_min(1)==1) then
        !
        ! x = xmin
        !
        phi_new(imin,:,:) = THREE*phi_new(imin+1,:,:) - THREE*phi_new(imin+2,:,:) + phi_new(imin+3,:,:)
        K_new(imin,:,:)   = THREE*K_new(imin+1,:,:)   - THREE*K_new(imin+2,:,:)   + K_new(imin+3,:,:)
        gxx_new(imin,:,:) = THREE*gxx_new(imin+1,:,:) - THREE*gxx_new(imin+2,:,:) + gxx_new(imin+3,:,:)
        gxy_new(imin,:,:) = THREE*gxy_new(imin+1,:,:) - THREE*gxy_new(imin+2,:,:) + gxy_new(imin+3,:,:)
        gxz_new(imin,:,:) = THREE*gxz_new(imin+1,:,:) - THREE*gxz_new(imin+2,:,:) + gxz_new(imin+3,:,:)
        gyy_new(imin,:,:) = THREE*gyy_new(imin+1,:,:) - THREE*gyy_new(imin+2,:,:) + gyy_new(imin+3,:,:)
        gyz_new(imin,:,:) = THREE*gyz_new(imin+1,:,:) - THREE*gyz_new(imin+2,:,:) + gyz_new(imin+3,:,:)
        gzz_new(imin,:,:) = THREE*gzz_new(imin+1,:,:) - THREE*gzz_new(imin+2,:,:) + gzz_new(imin+3,:,:)
        Axx_new(imin,:,:) = THREE*Axx_new(imin+1,:,:) - THREE*Axx_new(imin+2,:,:) + Axx_new(imin+3,:,:)
        Axy_new(imin,:,:) = THREE*Axy_new(imin+1,:,:) - THREE*Axy_new(imin+2,:,:) + Axy_new(imin+3,:,:)
        Axz_new(imin,:,:) = THREE*Axz_new(imin+1,:,:) - THREE*Axz_new(imin+2,:,:) + Axz_new(imin+3,:,:)
        Ayy_new(imin,:,:) = THREE*Ayy_new(imin+1,:,:) - THREE*Ayy_new(imin+2,:,:) + Ayy_new(imin+3,:,:)
        Ayz_new(imin,:,:) = THREE*Ayz_new(imin+1,:,:) - THREE*Ayz_new(imin+2,:,:) + Ayz_new(imin+3,:,:)
        Azz_new(imin,:,:) = THREE*Azz_new(imin+1,:,:) - THREE*Azz_new(imin+2,:,:) + Azz_new(imin+3,:,:)
        Gammax_new(imin,:,:) = THREE*Gammax_new(imin+1,:,:) - THREE*Gammax_new(imin+2,:,:) + Gammax_new(imin+3,:,:)
        Gammay_new(imin,:,:) = THREE*Gammay_new(imin+1,:,:) - THREE*Gammay_new(imin+2,:,:) + Gammay_new(imin+3,:,:)
        Gammaz_new(imin,:,:) = THREE*Gammaz_new(imin+1,:,:) - THREE*Gammaz_new(imin+2,:,:) + Gammaz_new(imin+3,:,:)
     endif
     if(have_bdry_min(2)==1) then
        !
        ! y = ymin
        !
        phi_new(:,jmin,:) = THREE*phi_new(:,jmin+1,:) - THREE*phi_new(:,jmin+2,:) + phi_new(:,jmin+3,:)
        K_new(:,jmin,:)   = THREE*K_new(:,jmin+1,:)   - THREE*K_new(:,jmin+2,:)   + K_new(:,jmin+3,:)
        gxx_new(:,jmin,:) = THREE*gxx_new(:,jmin+1,:) - THREE*gxx_new(:,jmin+2,:) + gxx_new(:,jmin+3,:)
        gxy_new(:,jmin,:) = THREE*gxy_new(:,jmin+1,:) - THREE*gxy_new(:,jmin+2,:) + gxy_new(:,jmin+3,:)
        gxz_new(:,jmin,:) = THREE*gxz_new(:,jmin+1,:) - THREE*gxz_new(:,jmin+2,:) + gxz_new(:,jmin+3,:)
        gyy_new(:,jmin,:) = THREE*gyy_new(:,jmin+1,:) - THREE*gyy_new(:,jmin+2,:) + gyy_new(:,jmin+3,:)
        gyz_new(:,jmin,:) = THREE*gyz_new(:,jmin+1,:) - THREE*gyz_new(:,jmin+2,:) + gyz_new(:,jmin+3,:)
        gzz_new(:,jmin,:) = THREE*gzz_new(:,jmin+1,:) - THREE*gzz_new(:,jmin+2,:) + gzz_new(:,jmin+3,:)
        Axx_new(:,jmin,:) = THREE*Axx_new(:,jmin+1,:) - THREE*Axx_new(:,jmin+2,:) + Axx_new(:,jmin+3,:)
        Axy_new(:,jmin,:) = THREE*Axy_new(:,jmin+1,:) - THREE*Axy_new(:,jmin+2,:) + Axy_new(:,jmin+3,:)
        Axz_new(:,jmin,:) = THREE*Axz_new(:,jmin+1,:) - THREE*Axz_new(:,jmin+2,:) + Axz_new(:,jmin+3,:)
        Ayy_new(:,jmin,:) = THREE*Ayy_new(:,jmin+1,:) - THREE*Ayy_new(:,jmin+2,:) + Ayy_new(:,jmin+3,:)
        Ayz_new(:,jmin,:) = THREE*Ayz_new(:,jmin+1,:) - THREE*Ayz_new(:,jmin+2,:) + Ayz_new(:,jmin+3,:)
        Azz_new(:,jmin,:) = THREE*Azz_new(:,jmin+1,:) - THREE*Azz_new(:,jmin+2,:) + Azz_new(:,jmin+3,:)
        Gammax_new(:,jmin,:) = THREE*Gammax_new(:,jmin+1,:) - THREE*Gammax_new(:,jmin+2,:) + Gammax_new(:,jmin+3,:)
        Gammay_new(:,jmin,:) = THREE*Gammay_new(:,jmin+1,:) - THREE*Gammay_new(:,jmin+2,:) + Gammay_new(:,jmin+3,:)
        Gammaz_new(:,jmin,:) = THREE*Gammaz_new(:,jmin+1,:) - THREE*Gammaz_new(:,jmin+2,:) + Gammaz_new(:,jmin+3,:)
     end if
     if(have_bdry_min(3)==1) then
        !
        ! z = zmin
        !
        phi_new(:,:,kmin) = THREE*phi_new(:,:,kmin+1) - THREE*phi_new(:,:,kmin+2) + phi_new(:,:,kmin+3)
        K_new(:,:,kmin)   = THREE*K_new(:,:,kmin+1)   - THREE*K_new(:,:,kmin+2)   + K_new(:,:,kmin+3)
        gxx_new(:,:,kmin) = THREE*gxx_new(:,:,kmin+1) - THREE*gxx_new(:,:,kmin+2) + gxx_new(:,:,kmin+3)
        gxy_new(:,:,kmin) = THREE*gxy_new(:,:,kmin+1) - THREE*gxy_new(:,:,kmin+2) + gxy_new(:,:,kmin+3)
        gxz_new(:,:,kmin) = THREE*gxz_new(:,:,kmin+1) - THREE*gxz_new(:,:,kmin+2) + gxz_new(:,:,kmin+3)
        gyy_new(:,:,kmin) = THREE*gyy_new(:,:,kmin+1) - THREE*gyy_new(:,:,kmin+2) + gyy_new(:,:,kmin+3)
        gyz_new(:,:,kmin) = THREE*gyz_new(:,:,kmin+1) - THREE*gyz_new(:,:,kmin+2) + gyz_new(:,:,kmin+3)
        gzz_new(:,:,kmin) = THREE*gzz_new(:,:,kmin+1) - THREE*gzz_new(:,:,kmin+2) + gzz_new(:,:,kmin+3)
        Axx_new(:,:,kmin) = THREE*Axx_new(:,:,kmin+1) - THREE*Axx_new(:,:,kmin+2) + Axx_new(:,:,kmin+3)
        Axy_new(:,:,kmin) = THREE*Axy_new(:,:,kmin+1) - THREE*Axy_new(:,:,kmin+2) + Axy_new(:,:,kmin+3)
        Axz_new(:,:,kmin) = THREE*Axz_new(:,:,kmin+1) - THREE*Axz_new(:,:,kmin+2) + Axz_new(:,:,kmin+3)
        Ayy_new(:,:,kmin) = THREE*Ayy_new(:,:,kmin+1) - THREE*Ayy_new(:,:,kmin+2) + Ayy_new(:,:,kmin+3)
        Ayz_new(:,:,kmin) = THREE*Ayz_new(:,:,kmin+1) - THREE*Ayz_new(:,:,kmin+2) + Ayz_new(:,:,kmin+3)
        Azz_new(:,:,kmin) = THREE*Azz_new(:,:,kmin+1) - THREE*Azz_new(:,:,kmin+2) + Azz_new(:,:,kmin+3)
        Gammax_new(:,:,kmin) = THREE*Gammax_new(:,:,kmin+1) - THREE*Gammax_new(:,:,kmin+2) + Gammax_new(:,:,kmin+3)
        Gammay_new(:,:,kmin) = THREE*Gammay_new(:,:,kmin+1) - THREE*Gammay_new(:,:,kmin+2) + Gammay_new(:,:,kmin+3)
        Gammaz_new(:,:,kmin) = THREE*Gammaz_new(:,:,kmin+1) - THREE*Gammaz_new(:,:,kmin+2) + Gammaz_new(:,:,kmin+3)
     endif
  else if(bc_type==EXTRAPR) then
   write(*,*)                    "BSSN BC_TYPE: EXTRAPR"
     ! Transform components of gij and Aij from fisheye to physical coordinates
     if (Fisheye_Flag==1) then
        fish_to_phys = .TRUE.
        call trans_fish_phys(ex, fake_ex, X, Y, Z, PhysicalR, RadiusDerivative, &
             gxx_new,gxy_new,gxz_new,gyy_new,gyz_new,gzz_new, &
             fish_to_phys,2,have_bdry_min,have_bdry_max,use_trans_fish_phys_new)
        call trans_fish_phys(ex, fake_ex, X, Y, Z, PhysicalR, RadiusDerivative, &
             Axx_new,Axy_new,Axz_new,Ayy_new,Ayz_new,Azz_new, &
             fish_to_phys,2,have_bdry_min,have_bdry_max,use_trans_fish_phys_new)
     end if
     !
     ! x = max
     !
     if(have_bdry_max(1)==1) then
        i = imax
        do k = kmin, kmax
           do j = jmin, jmax
              r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
              r1 = sqrt(X(i-1,j,k)**2 + Y(i-1,j,k)**2 + Z(i-1,j,k)**2)
              r2 = sqrt(X(i-2,j,k)**2 + Y(i-2,j,k)**2 + Z(i-2,j,k)**2)
              rp  = PhysicalR(i,j,k)
              rp1 = PhysicalR(i-1,j,k)
              rp2 = PhysicalR(i-2,j,k)
              fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
              fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i-1,j,k))
              fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i-2,j,k))
              rp1=rp1/rp
              rp2=rp2/rp
              phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i-1,j,k)-fac1)-rp2*(phi_new(i-2,j,k)-fac2))
              K_new(i,j,k)   = rp1*TWO*K_new(i-1,j,k) - rp2*K_new(i-2,j,k)
              gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i-1,j,k)-1.0d0) - rp2*(gxx_new(i-2,j,k)-1.0d0)
              gxy_new(i,j,k) = rp1*TWO*gxy_new(i-1,j,k) - rp2*gxy_new(i-2,j,k)
              gxz_new(i,j,k) = rp1*TWO*gxz_new(i-1,j,k) - rp2*gxz_new(i-2,j,k)
              gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i-1,j,k)-1.0d0) - rp2*(gyy_new(i-2,j,k)-1.0d0)
              gyz_new(i,j,k) = rp1*TWO*gyz_new(i-1,j,k) - rp2*gyz_new(i-2,j,k)
              gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i-1,j,k)-1.0d0) - rp2*(gzz_new(i-2,j,k)-1.0d0)
              Axx_new(i,j,k) = rp1*TWO*Axx_new(i-1,j,k) - rp2*Axx_new(i-2,j,k)
              Axy_new(i,j,k) = rp1*TWO*Axy_new(i-1,j,k) - rp2*Axy_new(i-2,j,k)
              Axz_new(i,j,k) = rp1*TWO*Axz_new(i-1,j,k) - rp2*Axz_new(i-2,j,k)
              Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i-1,j,k) - rp2*Ayy_new(i-2,j,k)
              Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i-1,j,k) - rp2*Ayz_new(i-2,j,k)
              Azz_new(i,j,k) = rp1*TWO*Azz_new(i-1,j,k) - rp2*Azz_new(i-2,j,k)
              Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i-1,j,k) - rp2*Gammax_new(i-2,j,k)
              Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i-1,j,k) - rp2*Gammay_new(i-2,j,k)
              Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i-1,j,k) - rp2*Gammaz_new(i-2,j,k)
           enddo
        enddo
     end if
     !
     ! y = ymax
     !
     if(have_bdry_max(2)==1) then
        j = jmax
        do k = kmin, kmax
           do i = imin, imax
              r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
              r1 = sqrt(X(i,j-1,k)**2 + Y(i,j-1,k)**2 + Z(i,j-1,k)**2)
              r2 = sqrt(X(i,j-2,k)**2 + Y(i,j-2,k)**2 + Z(i,j-2,k)**2)
              rp  = PhysicalR(i,j,k)
              rp1 = PhysicalR(i,j-1,k)
              rp2 = PhysicalR(i,j-2,k)
              fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
              fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i,j-1,k))
              fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i,j-2,k))
              rp1=rp1/rp
              rp2=rp2/rp
              phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i,j-1,k)-fac1)-rp2*(phi_new(i,j-2,k)-fac2))
              K_new(i,j,k)   = rp1*TWO*K_new(i,j-1,k) - rp2*K_new(i,j-2,k)
              gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i,j-1,k)-1.0d0) - rp2*(gxx_new(i,j-2,k)-1.0d0)
              gxy_new(i,j,k) = rp1*TWO*gxy_new(i,j-1,k) - rp2*gxy_new(i,j-2,k)
              gxz_new(i,j,k) = rp1*TWO*gxz_new(i,j-1,k) - rp2*gxz_new(i,j-2,k)
              gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i,j-1,k)-1.0d0) - rp2*(gyy_new(i,j-2,k)-1.0d0)
              gyz_new(i,j,k) = rp1*TWO*gyz_new(i,j-1,k) - rp2*gyz_new(i,j-2,k)
              gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i,j-1,k)-1.0d0) - rp2*(gzz_new(i,j-2,k)-1.0d0)
              Axx_new(i,j,k) = rp1*TWO*Axx_new(i,j-1,k) - rp2*Axx_new(i,j-2,k)
              Axy_new(i,j,k) = rp1*TWO*Axy_new(i,j-1,k) - rp2*Axy_new(i,j-2,k)
              Axz_new(i,j,k) = rp1*TWO*Axz_new(i,j-1,k) - rp2*Axz_new(i,j-2,k)
              Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i,j-1,k) - rp2*Ayy_new(i,j-2,k)
              Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i,j-1,k) - rp2*Ayz_new(i,j-2,k)
              Azz_new(i,j,k) = rp1*TWO*Azz_new(i,j-1,k) - rp2*Azz_new(i,j-2,k)
              Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i,j-1,k) - rp2*Gammax_new(i,j-2,k)
              Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i,j-1,k) - rp2*Gammay_new(i,j-2,k)
              Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i,j-1,k) - rp2*Gammaz_new(i,j-2,k)
           enddo
        enddo
     endif
     !
     ! z = zmax
     !
     if(have_bdry_max(3)==1) then
        k = kmax
        do j = jmin, jmax
           do i = imin, imax
              r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
              r1 = sqrt(X(i,j,k-1)**2 + Y(i,j,k-1)**2 + Z(i,j,k-1)**2)
              r2 = sqrt(X(i,j,k-2)**2 + Y(i,j,k-2)**2 + Z(i,j,k-2)**2)
              rp  = PhysicalR(i,j,k)
              rp1 = PhysicalR(i,j,k-1)
              rp2 = PhysicalR(i,j,k-2)
              fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
              fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i,j,k-1))
              fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i,j,k-2))
              rp1=rp1/rp
              rp2=rp2/rp
              phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i,j,k-1)-fac1)-rp2*(phi_new(i,j,k-2)-fac2))
              K_new(i,j,k)   = rp1*TWO*K_new(i,j,k-1) - rp2*K_new(i,j,k-2)
              gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i,j,k-1)-1.0d0) - rp2*(gxx_new(i,j,k-2)-1.0d0)
              gxy_new(i,j,k) = rp1*TWO*gxy_new(i,j,k-1) - rp2*gxy_new(i,j,k-2)
              gxz_new(i,j,k) = rp1*TWO*gxz_new(i,j,k-1) - rp2*gxz_new(i,j,k-2)
              gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i,j,k-1)-1.0d0) - rp2*(gyy_new(i,j,k-2)-1.0d0)
              gyz_new(i,j,k) = rp1*TWO*gyz_new(i,j,k-1) - rp2*gyz_new(i,j,k-2)
              gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i,j,k-1)-1.0d0) - rp2*(gzz_new(i,j,k-2)-1.0d0)
              Axx_new(i,j,k) = rp1*TWO*Axx_new(i,j,k-1) - rp2*Axx_new(i,j,k-2)
              Axy_new(i,j,k) = rp1*TWO*Axy_new(i,j,k-1) - rp2*Axy_new(i,j,k-2)
              Axz_new(i,j,k) = rp1*TWO*Axz_new(i,j,k-1) - rp2*Axz_new(i,j,k-2)
              Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i,j,k-1) - rp2*Ayy_new(i,j,k-2)
              Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i,j,k-1) - rp2*Ayz_new(i,j,k-2)
              Azz_new(i,j,k) = rp1*TWO*Azz_new(i,j,k-1) - rp2*Azz_new(i,j,k-2)
              Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i,j,k-1) - rp2*Gammax_new(i,j,k-2)
              Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i,j,k-1) - rp2*Gammay_new(i,j,k-2)
              Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i,j,k-1) - rp2*Gammaz_new(i,j,k-2)
           enddo
        enddo
     end if
     if(have_bdry_min(1)==1) then
        !
        ! x = xmin
        !
        i = imin
        do k = kmin, kmax
           do j = jmin, jmax
              r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
              r1 = sqrt(X(i+1,j,k)**2 + Y(i+1,j,k)**2 + Z(i+1,j,k)**2)
              r2 = sqrt(X(i+2,j,k)**2 + Y(i+2,j,k)**2 + Z(i+2,j,k)**2)
              rp  = PhysicalR(i,j,k)
              rp1 = PhysicalR(i+1,j,k)
              rp2 = PhysicalR(i+2,j,k)
              fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
              fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i+1,j,k))
              fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i+2,j,k))
              rp1=rp1/rp
              rp2=rp2/rp
              phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i+1,j,k)-fac1)-rp2*(phi_new(i+2,j,k)-fac2))
              K_new(i,j,k) = rp1*TWO*K_new(imin+1,j,k) - rp2*K_new(imin+2,j,k)
              gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(imin+1,j,k)-1.0d0) - rp2*(gxx_new(imin+2,j,k)-1.0d0)
              gxy_new(i,j,k) = rp1*TWO*gxy_new(imin+1,j,k) - rp2*gxy_new(imin+2,j,k)
              gxz_new(i,j,k) = rp1*TWO*gxz_new(imin+1,j,k) - rp2*gxz_new(imin+2,j,k)
              gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(imin+1,j,k)-1.0d0) - rp2*(gyy_new(imin+2,j,k)-1.0d0)
              gyz_new(i,j,k) = rp1*TWO*gyz_new(imin+1,j,k) - rp2*gyz_new(imin+2,j,k)
              gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(imin+1,j,k)-1.0d0) - rp2*(gzz_new(imin+2,j,k)-1.0d0)
              Axx_new(i,j,k) = rp1*TWO*Axx_new(imin+1,j,k) - rp2*Axx_new(imin+2,j,k)
              Axy_new(i,j,k) = rp1*TWO*Axy_new(imin+1,j,k) - rp2*Axy_new(imin+2,j,k)
              Axz_new(i,j,k) = rp1*TWO*Axz_new(imin+1,j,k) - rp2*Axz_new(imin+2,j,k)
              Ayy_new(i,j,k) = rp1*TWO*Ayy_new(imin+1,j,k) - rp2*Ayy_new(imin+2,j,k)
              Ayz_new(i,j,k) = rp1*TWO*Ayz_new(imin+1,j,k) - rp2*Ayz_new(imin+2,j,k)
              Azz_new(i,j,k) = rp1*TWO*Azz_new(imin+1,j,k) - rp2*Azz_new(imin+2,j,k)
              Gammax_new(i,j,k) = rp1*TWO*Gammax_new(imin+1,j,k) - rp2*Gammax_new(imin+2,j,k)
              Gammay_new(i,j,k) = rp1*TWO*Gammay_new(imin+1,j,k) - rp2*Gammay_new(imin+2,j,k)
              Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(imin+1,j,k) - rp2*Gammaz_new(imin+2,j,k)
           enddo
        enddo
     endif
     if(have_bdry_min(2)==1) then
        !
        ! y = ymin
        !
        j = jmin
        do k = kmin, kmax
           do i = imin, imax
              r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
              r1 = sqrt(X(i,j+1,k)**2 + Y(i,j+1,k)**2 + Z(i,j+1,k)**2)
              r2 = sqrt(X(i,j+2,k)**2 + Y(i,j+2,k)**2 + Z(i,j+2,k)**2)
              rp  = PhysicalR(i,j,k)
              rp1 = PhysicalR(i,j+1,k)
              rp2 = PhysicalR(i,j+2,k)
              fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
              fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i,j+1,k))
              fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i,j+2,k))
              rp1=rp1/rp
              rp2=rp2/rp
              phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i,j+1,k)-fac1)-rp2*(phi_new(i,j+2,k)-fac2))
              K_new(i,j,k) = rp1*TWO*K_new(i,jmin+1,k) - rp2*K_new(i,jmin+2,k)
              gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i,jmin+1,k)-1.0d0) - rp2*(gxx_new(i,jmin+2,k)-1.0d0)
              gxy_new(i,j,k) = rp1*TWO*gxy_new(i,jmin+1,k) - rp2*gxy_new(i,jmin+2,k)
              gxz_new(i,j,k) = rp1*TWO*gxz_new(i,jmin+1,k) - rp2*gxz_new(i,jmin+2,k)
              gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i,jmin+1,k)-1.0d0) - rp2*(gyy_new(i,jmin+2,k)-1.0d0)
              gyz_new(i,j,k) = rp1*TWO*gyz_new(i,jmin+1,k) - rp2*gyz_new(i,jmin+2,k)
              gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i,jmin+1,k)-1.0d0) - rp2*(gzz_new(i,jmin+2,k)-1.0d0)
              Axx_new(i,j,k) = rp1*TWO*Axx_new(i,jmin+1,k) - rp2*Axx_new(i,jmin+2,k)
              Axy_new(i,j,k) = rp1*TWO*Axy_new(i,jmin+1,k) - rp2*Axy_new(i,jmin+2,k)
              Axz_new(i,j,k) = rp1*TWO*Axz_new(i,jmin+1,k) - rp2*Axz_new(i,jmin+2,k)
              Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i,jmin+1,k) - rp2*Ayy_new(i,jmin+2,k)
              Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i,jmin+1,k) - rp2*Ayz_new(i,jmin+2,k)
              Azz_new(i,j,k) = rp1*TWO*Azz_new(i,jmin+1,k) - rp2*Azz_new(i,jmin+2,k)
              Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i,jmin+1,k) - rp2*Gammax_new(i,jmin+2,k)
              Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i,jmin+1,k) - rp2*Gammay_new(i,jmin+2,k)
              Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i,jmin+1,k) - rp2*Gammaz_new(i,jmin+2,k)
           enddo
        enddo
     endif
     if(have_bdry_min(3)==1) then
        !
        ! z = zmin
        !
        k = kmin
        do j = jmin, jmax
           do i = imin, imax
              r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
              r1 = sqrt(X(i,j,k+1)**2 + Y(i,j,k+1)**2 + Z(i,j,k+1)**2)
              r2 = sqrt(X(i,j,k+2)**2 + Y(i,j,k+2)**2 + Z(i,j,k+2)**2)
              rp  = PhysicalR(i,j,k)
              rp1 = PhysicalR(i,j,k+1)
              rp2 = PhysicalR(i,j,k+2)
              fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
              fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i,j,k+1))
              fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i,j,k+2))
              rp1=rp1/rp
              rp2=rp2/rp
              phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i,j,k+1)-fac1)-rp2*(phi_new(i,j,k+2)-fac2))
              K_new(i,j,k) = rp1*TWO*K_new(i,j,kmin+1) - rp2*K_new(i,j,kmin+2)
              gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i,j,kmin+1)-1.0d0) - rp2*(gxx_new(i,j,kmin+2)-1.0d0)
              gxy_new(i,j,k) = rp1*TWO*gxy_new(i,j,kmin+1) - rp2*gxy_new(i,j,kmin+2)
              gxz_new(i,j,k) = rp1*TWO*gxz_new(i,j,kmin+1) - rp2*gxz_new(i,j,kmin+2)
              gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i,j,kmin+1)-1.0d0) - rp2*(gyy_new(i,j,kmin+2)-1.0d0)
              gyz_new(i,j,k) = rp1*TWO*gyz_new(i,j,kmin+1) - rp2*gyz_new(i,j,kmin+2)
              gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i,j,kmin+1)-1.0d0) - rp2*(gzz_new(i,j,kmin+2)-1.0d0)
              Axx_new(i,j,k) = rp1*TWO*Axx_new(i,j,kmin+1) - rp2*Axx_new(i,j,kmin+2)
              Axy_new(i,j,k) = rp1*TWO*Axy_new(i,j,kmin+1) - rp2*Axy_new(i,j,kmin+2)
              Axz_new(i,j,k) = rp1*TWO*Axz_new(i,j,kmin+1) - rp2*Axz_new(i,j,kmin+2)
              Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i,j,kmin+1) - rp2*Ayy_new(i,j,kmin+2)
              Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i,j,kmin+1) - rp2*Ayz_new(i,j,kmin+2)
              Azz_new(i,j,k) = rp1*TWO*Azz_new(i,j,kmin+1) - rp2*Azz_new(i,j,kmin+2)
              Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i,j,kmin+1) - rp2*Gammax_new(i,j,kmin+2)
              Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i,j,kmin+1) - rp2*Gammay_new(i,j,kmin+2)
              Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i,j,kmin+1) - rp2*Gammaz_new(i,j,kmin+2)
           enddo
        enddo
     end if
     !
     ! Now for the edges, to avoid instabilities
     !
     ! x=xmax,z=zmax (used by all symmetries)
     !
     if(have_bdry_max(1)==1.and.have_bdry_max(3)==1) then
        i = imax
        k = kmax
        do j = jmin, jmax
           r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
           r1 = sqrt(X(i-1,j,k-1)**2 + Y(i-1,j,k-1)**2 + Z(i-1,j,k-1)**2)
           r2 = sqrt(X(i-2,j,k-2)**2 + Y(i-2,j,k-2)**2 + Z(i-2,j,k-2)**2)
           rp  = PhysicalR(i,j,k)
           rp1 = PhysicalR(i-1,j,k-1)
           rp2 = PhysicalR(i-2,j,k-2)
           fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
           fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i-1,j,k-1))
           fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i-2,j,k-2))
           rp1=rp1/rp
           rp2=rp2/rp
           phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i-1,j,k-1)-fac1)-rp2*(phi_new(i-2,j,k-2)-fac2))
           K_new(i,j,k)   = rp1*TWO*K_new(i-1,j,k-1) - rp2*K_new(i-2,j,k-2)
           gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i-1,j,k-1)-1.0d0) - rp2*(gxx_new(i-2,j,k-2)-1.0d0)
           gxy_new(i,j,k) = rp1*TWO*gxy_new(i-1,j,k-1) - rp2*gxy_new(i-2,j,k-2)
           gxz_new(i,j,k) = rp1*TWO*gxz_new(i-1,j,k-1) - rp2*gxz_new(i-2,j,k-2)
           gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i-1,j,k-1)-1.0d0) - rp2*(gyy_new(i-2,j,k-2)-1.0d0)
           gyz_new(i,j,k) = rp1*TWO*gyz_new(i-1,j,k-1) - rp2*gyz_new(i-2,j,k-2)
           gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i-1,j,k-1)-1.0d0) - rp2*(gzz_new(i-2,j,k-2)-1.0d0)
           Axx_new(i,j,k) = rp1*TWO*Axx_new(i-1,j,k-1) - rp2*Axx_new(i-2,j,k-2)
           Axy_new(i,j,k) = rp1*TWO*Axy_new(i-1,j,k-1) - rp2*Axy_new(i-2,j,k-2)
           Axz_new(i,j,k) = rp1*TWO*Axz_new(i-1,j,k-1) - rp2*Axz_new(i-2,j,k-2)
           Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i-1,j,k-1) - rp2*Ayy_new(i-2,j,k-2)
           Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i-1,j,k-1) - rp2*Ayz_new(i-2,j,k-2)
           Azz_new(i,j,k) = rp1*TWO*Azz_new(i-1,j,k-1) - rp2*Azz_new(i-2,j,k-2)
           Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i-1,j,k-1) - rp2*Gammax_new(i-2,j,k-2)
           Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i-1,j,k-1) - rp2*Gammay_new(i-2,j,k-2)
           Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i-1,j,k-1) - rp2*Gammaz_new(i-2,j,k-2)
        enddo
     endif
     if(have_bdry_max(1)==1.and.have_bdry_max(2)==1) then
        !
        ! x=xmax,y=ymax (used by all symmetries except AXISYM)
        !
        i = imax
        j = jmax
        do k = kmin, kmax
           r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
           r1 = sqrt(X(i-1,j-1,k)**2 + Y(i-1,j-1,k)**2 + Z(i-1,j-1,k)**2)
           r2 = sqrt(X(i-2,j-2,k)**2 + Y(i-2,j-2,k)**2 + Z(i-2,j-2,k)**2)
           rp  = PhysicalR(i,j,k)
           rp1 = PhysicalR(i-1,j-1,k)
           rp2 = PhysicalR(i-2,j-2,k)
           fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
           fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i-1,j-1,k))
           fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i-2,j-2,k))
           rp1=rp1/rp
           rp2=rp2/rp
           phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i-1,j-1,k)-fac1)-rp2*(phi_new(i-2,j-2,k)-fac2))
           K_new(i,j,k)   = rp1*TWO*K_new(i-1,j-1,k) - rp2*K_new(i-2,j-2,k)
           gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i-1,j-1,k)-1.0d0) - rp2*(gxx_new(i-2,j-2,k)-1.0d0)
           gxy_new(i,j,k) = rp1*TWO*gxy_new(i-1,j-1,k) - rp2*gxy_new(i-2,j-2,k)
           gxz_new(i,j,k) = rp1*TWO*gxz_new(i-1,j-1,k) - rp2*gxz_new(i-2,j-2,k)
           gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i-1,j-1,k)-1.0d0) - rp2*(gyy_new(i-2,j-2,k)-1.0d0)
           gyz_new(i,j,k) = rp1*TWO*gyz_new(i-1,j-1,k) - rp2*gyz_new(i-2,j-2,k)
           gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i-1,j-1,k)-1.0d0) - rp2*(gzz_new(i-2,j-2,k)-1.0d0)
           Axx_new(i,j,k) = rp1*TWO*Axx_new(i-1,j-1,k) - rp2*Axx_new(i-2,j-2,k)
           Axy_new(i,j,k) = rp1*TWO*Axy_new(i-1,j-1,k) - rp2*Axy_new(i-2,j-2,k)
           Axz_new(i,j,k) = rp1*TWO*Axz_new(i-1,j-1,k) - rp2*Axz_new(i-2,j-2,k)
           Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i-1,j-1,k) - rp2*Ayy_new(i-2,j-2,k)
           Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i-1,j-1,k) - rp2*Ayz_new(i-2,j-2,k)
           Azz_new(i,j,k) = rp1*TWO*Azz_new(i-1,j-1,k) - rp2*Azz_new(i-2,j-2,k)
           Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i-1,j-1,k) - rp2*Gammax_new(i-2,j-2,k)
           Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i-1,j-1,k) - rp2*Gammay_new(i-2,j-2,k)
           Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i-1,j-1,k) - rp2*Gammaz_new(i-2,j-2,k)
        enddo
        !        endif
        !
        ! y=ymax,z=zmax (used by all symmetries except AXISYM)
        !
        if(have_bdry_max(2)==1.and.have_bdry_max(3)==1) then
           j = jmax
           k = kmax
           do i = imin, imax
              r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
              r1 = sqrt(X(i,j-1,k-1)**2 + Y(i,j-1,k-1)**2 + Z(i,j-1,k-1)**2)
              r2 = sqrt(X(i,j-2,k-2)**2 + Y(i,j-2,k-2)**2 + Z(i,j-2,k-2)**2)
              rp  = PhysicalR(i,j,k)
              rp1 = PhysicalR(i,j-1,k-1)
              rp2 = PhysicalR(i,j-2,k-2)
              fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
              fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i,j-1,k-1))
              fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i,j-2,k-2))
              rp1=rp1/rp
              rp2=rp2/rp
              phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i,j-1,k-1)-fac1)-rp2*(phi_new(i,j-2,k-2)-fac2))
              K_new(i,j,k)   = rp1*TWO*K_new(i,j-1,k-1) - rp2*K_new(i,j-2,k-2)
              gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i,j-1,k-1)-1.0d0) - rp2*(gxx_new(i,j-2,k-2)-1.0d0)
              gxy_new(i,j,k) = rp1*TWO*gxy_new(i,j-1,k-1) - rp2*gxy_new(i,j-2,k-2)
              gxz_new(i,j,k) = rp1*TWO*gxz_new(i,j-1,k-1) - rp2*gxz_new(i,j-2,k-2)
              gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i,j-1,k-1)-1.0d0) - rp2*(gyy_new(i,j-2,k-2)-1.0d0)
              gyz_new(i,j,k) = rp1*TWO*gyz_new(i,j-1,k-1) - rp2*gyz_new(i,j-2,k-2)
              gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i,j-1,k-1)-1.0d0) - rp2*(gzz_new(i,j-2,k-2)-1.0d0)
              Axx_new(i,j,k) = rp1*TWO*Axx_new(i,j-1,k-1) - rp2*Axx_new(i,j-2,k-2)
              Axy_new(i,j,k) = rp1*TWO*Axy_new(i,j-1,k-1) - rp2*Axy_new(i,j-2,k-2)
              Axz_new(i,j,k) = rp1*TWO*Axz_new(i,j-1,k-1) - rp2*Axz_new(i,j-2,k-2)
              Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i,j-1,k-1) - rp2*Ayy_new(i,j-2,k-2)
              Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i,j-1,k-1) - rp2*Ayz_new(i,j-2,k-2)
              Azz_new(i,j,k) = rp1*TWO*Azz_new(i,j-1,k-1) - rp2*Azz_new(i,j-2,k-2)
              Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i,j-1,k-1) - rp2*Gammax_new(i,j-2,k-2)
              Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i,j-1,k-1) - rp2*Gammay_new(i,j-2,k-2)
              Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i,j-1,k-1) - rp2*Gammaz_new(i,j-2,k-2)
           enddo
        endif
     endif
     if(have_bdry_min(1)==1.and.have_bdry_max(2)==1) then
        !
        ! x=xmin,y=ymax
        !
        i = imin
        j = jmax
        do k = kmin, kmax
           r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
           r1 = sqrt(X(i+1,j-1,k)**2 + Y(i+1,j-1,k)**2 + Z(i+1,j-1,k)**2)
           r2 = sqrt(X(i+2,j-2,k)**2 + Y(i+2,j-2,k)**2 + Z(i+2,j-2,k)**2)
           rp  = PhysicalR(i,j,k)
           rp1 = PhysicalR(i+1,j-1,k)
           rp2 = PhysicalR(i+2,j-2,k)
           fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
           fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i+1,j-1,k))
           fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i+2,j-2,k))
           rp1=rp1/rp
           rp2=rp2/rp
           phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i+1,j-1,k)-fac1)-rp2*(phi_new(i+2,j-2,k)-fac2))
           K_new(i,j,k)   = rp1*TWO*K_new(i+1,j-1,k) - rp2*K_new(i+2,j-2,k)
           gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i+1,j-1,k)-1.0d0) - rp2*(gxx_new(i+2,j-2,k)-1.0d0)
           gxy_new(i,j,k) = rp1*TWO*gxy_new(i+1,j-1,k) - rp2*gxy_new(i+2,j-2,k)
           gxz_new(i,j,k) = rp1*TWO*gxz_new(i+1,j-1,k) - rp2*gxz_new(i+2,j-2,k)
           gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i+1,j-1,k)-1.0d0) - rp2*(gyy_new(i+2,j-2,k)-1.0d0)
           gyz_new(i,j,k) = rp1*TWO*gyz_new(i+1,j-1,k) - rp2*gyz_new(i+2,j-2,k)
           gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i+1,j-1,k)-1.0d0) - rp2*(gzz_new(i+2,j-2,k)-1.0d0)
           Axx_new(i,j,k) = rp1*TWO*Axx_new(i+1,j-1,k) - rp2*Axx_new(i+2,j-2,k)
           Axy_new(i,j,k) = rp1*TWO*Axy_new(i+1,j-1,k) - rp2*Axy_new(i+2,j-2,k)
           Axz_new(i,j,k) = rp1*TWO*Axz_new(i+1,j-1,k) - rp2*Axz_new(i+2,j-2,k)
           Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i+1,j-1,k) - rp2*Ayy_new(i+2,j-2,k)
           Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i+1,j-1,k) - rp2*Ayz_new(i+2,j-2,k)
           Azz_new(i,j,k) = rp1*TWO*Azz_new(i+1,j-1,k) - rp2*Azz_new(i+2,j-2,k)
           Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i+1,j-1,k) - rp2*Gammax_new(i+2,j-2,k)
           Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i+1,j-1,k) - rp2*Gammay_new(i+2,j-2,k)
           Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i+1,j-1,k) - rp2*Gammaz_new(i+2,j-2,k)
        enddo
     endif
     !
     ! x=xmin,z=zmax
     !
     if(have_bdry_min(1)==1.and.have_bdry_max(3)==1) then
        i = imin
        k = kmax
        do j = jmin, jmax
           r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
           r1 = sqrt(X(i+1,j,k-1)**2 + Y(i+1,j,k-1)**2 + Z(i+1,j,k-1)**2)
           r2 = sqrt(X(i+2,j,k-2)**2 + Y(i+2,j,k-2)**2 + Z(i+2,j,k-2)**2)
           rp  = PhysicalR(i,j,k)
           rp1 = PhysicalR(i+1,j,k-1)
           rp2 = PhysicalR(i+2,j,k-2)
           fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
           fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i+1,j,k-1))
           fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i+2,j,k-2))
           rp1=rp1/rp
           rp2=rp2/rp
           phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i+1,j,k-1)-fac1)-rp2*(phi_new(i+2,j,k-2)-fac2))
           K_new(i,j,k)   = rp1*TWO*K_new(i+1,j,k-1) - rp2*K_new(i+2,j,k-2)
           gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i+1,j,k-1)-1.0d0) - rp2*(gxx_new(i+2,j,k-2)-1.0d0)
           gxy_new(i,j,k) = rp1*TWO*gxy_new(i+1,j,k-1) - rp2*gxy_new(i+2,j,k-2)
           gxz_new(i,j,k) = rp1*TWO*gxz_new(i+1,j,k-1) - rp2*gxz_new(i+2,j,k-2)
           gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i+1,j,k-1)-1.0d0) - rp2*(gyy_new(i+2,j,k-2)-1.0d0)
           gyz_new(i,j,k) = rp1*TWO*gyz_new(i+1,j,k-1) - rp2*gyz_new(i+2,j,k-2)
           gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i+1,j,k-1)-1.0d0) - rp2*(gzz_new(i+2,j,k-2)-1.0d0)
           Axx_new(i,j,k) = rp1*TWO*Axx_new(i+1,j,k-1) - rp2*Axx_new(i+2,j,k-2)
           Axy_new(i,j,k) = rp1*TWO*Axy_new(i+1,j,k-1) - rp2*Axy_new(i+2,j,k-2)
           Axz_new(i,j,k) = rp1*TWO*Axz_new(i+1,j,k-1) - rp2*Axz_new(i+2,j,k-2)
           Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i+1,j,k-1) - rp2*Ayy_new(i+2,j,k-2)
           Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i+1,j,k-1) - rp2*Ayz_new(i+2,j,k-2)
           Azz_new(i,j,k) = rp1*TWO*Azz_new(i+1,j,k-1) - rp2*Azz_new(i+2,j,k-2)
           Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i+1,j,k-1) - rp2*Gammax_new(i+2,j,k-2)
           Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i+1,j,k-1) - rp2*Gammay_new(i+2,j,k-2)
           Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i+1,j,k-1) - rp2*Gammaz_new(i+2,j,k-2)
        enddo
     endif
     !
     ! x=xmax,y=ymin (used by EQ and NO)
     !
     if(have_bdry_max(1)==1.and.have_bdry_min(2)==1) then
        i = imax
        j = jmin
        do k = kmin, kmax
           r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
           r1 = sqrt(X(i-1,j+1,k)**2 + Y(i-1,j+1,k)**2 + Z(i-1,j+1,k)**2)
           r2 = sqrt(X(i-2,j+2,k)**2 + Y(i-2,j+2,k)**2 + Z(i-2,j+2,k)**2)
           rp  = PhysicalR(i,j,k)
           rp1 = PhysicalR(i-1,j+1,k)
           rp2 = PhysicalR(i-2,j+2,k)
           fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
           fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i-1,j+1,k))
           fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i-2,j+2,k))
           rp1=rp1/rp
           rp2=rp2/rp
           phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i-1,j+1,k)-fac1)-rp2*(phi_new(i-2,j+2,k)-fac2))
           K_new(i,j,k)   = rp1*TWO*K_new(i-1,j+1,k) - rp2*K_new(i-2,j+2,k)
           gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i-1,j+1,k)-1.0d0) - rp2*(gxx_new(i-2,j+2,k)-1.0d0)
           gxy_new(i,j,k) = rp1*TWO*gxy_new(i-1,j+1,k) - rp2*gxy_new(i-2,j+2,k)
           gxz_new(i,j,k) = rp1*TWO*gxz_new(i-1,j+1,k) - rp2*gxz_new(i-2,j+2,k)
           gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i-1,j+1,k)-1.0d0) - rp2*(gyy_new(i-2,j+2,k)-1.0d0)
           gyz_new(i,j,k) = rp1*TWO*gyz_new(i-1,j+1,k) - rp2*gyz_new(i-2,j+2,k)
           gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i-1,j+1,k)-1.0d0) - rp2*(gzz_new(i-2,j+2,k)-1.0d0)
           Axx_new(i,j,k) = rp1*TWO*Axx_new(i-1,j+1,k) - rp2*Axx_new(i-2,j+2,k)
           Axy_new(i,j,k) = rp1*TWO*Axy_new(i-1,j+1,k) - rp2*Axy_new(i-2,j+2,k)
           Axz_new(i,j,k) = rp1*TWO*Axz_new(i-1,j+1,k) - rp2*Axz_new(i-2,j+2,k)
           Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i-1,j+1,k) - rp2*Ayy_new(i-2,j+2,k)
           Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i-1,j+1,k) - rp2*Ayz_new(i-2,j+2,k)
           Azz_new(i,j,k) = rp1*TWO*Azz_new(i-1,j+1,k) - rp2*Azz_new(i-2,j+2,k)
           Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i-1,j+1,k) - rp2*Gammax_new(i-2,j+2,k)
           Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i-1,j+1,k) - rp2*Gammay_new(i-2,j+2,k)
           Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i-1,j+1,k) - rp2*Gammaz_new(i-2,j+2,k)
        enddo
     endif
     !
     ! x=xmin,y=ymin (used by EQ and NO)
     !
     if(have_bdry_min(1)==1.and.have_bdry_min(2)==1) then
        i = imin
        j = jmin
        do k = kmin, kmax
           r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
           r1 = sqrt(X(i+1,j+1,k)**2 + Y(i+1,j+1,k)**2 + Z(i+1,j+1,k)**2)
           r2 = sqrt(X(i+2,j+2,k)**2 + Y(i+2,j+2,k)**2 + Z(i+2,j+2,k)**2)
           rp  = PhysicalR(i,j,k)
           rp1 = PhysicalR(i+1,j+1,k)
           rp2 = PhysicalR(i+2,j+2,k)
           fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
           fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i+1,j+1,k))
           fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i+2,j+2,k))
           rp1=rp1/rp
           rp2=rp2/rp
           phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i+1,j+1,k)-fac1)-rp2*(phi_new(i+2,j+2,k)-fac2))
           K_new(i,j,k)   = rp1*TWO*K_new(i+1,j+1,k) - rp2*K_new(i+2,j+2,k)
           gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i+1,j+1,k)-1.0d0) - rp2*(gxx_new(i+2,j+2,k)-1.0d0)
           gxy_new(i,j,k) = rp1*TWO*gxy_new(i+1,j+1,k) - rp2*gxy_new(i+2,j+2,k)
           gxz_new(i,j,k) = rp1*TWO*gxz_new(i+1,j+1,k) - rp2*gxz_new(i+2,j+2,k)
           gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i+1,j+1,k)-1.0d0) - rp2*(gyy_new(i+2,j+2,k)-1.0d0)
           gyz_new(i,j,k) = rp1*TWO*gyz_new(i+1,j+1,k) - rp2*gyz_new(i+2,j+2,k)
           gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i+1,j+1,k)-1.0d0) - rp2*(gzz_new(i+2,j+2,k)-1.0d0)
           Axx_new(i,j,k) = rp1*TWO*Axx_new(i+1,j+1,k) - rp2*Axx_new(i+2,j+2,k)
           Axy_new(i,j,k) = rp1*TWO*Axy_new(i+1,j+1,k) - rp2*Axy_new(i+2,j+2,k)
           Axz_new(i,j,k) = rp1*TWO*Axz_new(i+1,j+1,k) - rp2*Axz_new(i+2,j+2,k)
           Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i+1,j+1,k) - rp2*Ayy_new(i+2,j+2,k)
           Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i+1,j+1,k) - rp2*Ayz_new(i+2,j+2,k)
           Azz_new(i,j,k) = rp1*TWO*Azz_new(i+1,j+1,k) - rp2*Azz_new(i+2,j+2,k)
           Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i+1,j+1,k) - rp2*Gammax_new(i+2,j+2,k)
           Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i+1,j+1,k) - rp2*Gammay_new(i+2,j+2,k)
           Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i+1,j+1,k) - rp2*Gammaz_new(i+2,j+2,k)
        enddo
     endif
     !
     ! y=ymin,z=zmax (used by EQ and NO)
     !
     if(have_bdry_min(2)==1.and.have_bdry_max(3)==1) then
        j = jmin
        k = kmax
        do i = imin, imax
           r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
           r1 = sqrt(X(i,j+1,k-1)**2 + Y(i,j+1,k-1)**2 + Z(i,j+1,k-1)**2)
           r2 = sqrt(X(i,j+2,k-2)**2 + Y(i,j+2,k-2)**2 + Z(i,j+2,k-2)**2)
           rp  = PhysicalR(i,j,k)
           rp1 = PhysicalR(i,j+1,k-1)
           rp2 = PhysicalR(i,j+2,k-2)
           fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
           fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i,j+1,k-1))
           fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i,j+2,k-2))
           rp1=rp1/rp
           rp2=rp2/rp
           phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i,j+1,k-1)-fac1)-rp2*(phi_new(i,j+2,k-2)-fac2))
           K_new(i,j,k)   = rp1*TWO*K_new(i,j+1,k-1) - rp2*K_new(i,j+2,k-2)
           gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i,j+1,k-1)-1.0d0) - rp2*(gxx_new(i,j+2,k-2)-1.0d0)
           gxy_new(i,j,k) = rp1*TWO*gxy_new(i,j+1,k-1) - rp2*gxy_new(i,j+2,k-2)
           gxz_new(i,j,k) = rp1*TWO*gxz_new(i,j+1,k-1) - rp2*gxz_new(i,j+2,k-2)
           gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i,j+1,k-1)-1.0d0) - rp2*(gyy_new(i,j+2,k-2)-1.0d0)
           gyz_new(i,j,k) = rp1*TWO*gyz_new(i,j+1,k-1) - rp2*gyz_new(i,j+2,k-2)
           gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i,j+1,k-1)-1.0d0) - rp2*(gzz_new(i,j+2,k-2)-1.0d0)
           Axx_new(i,j,k) = rp1*TWO*Axx_new(i,j+1,k-1) - rp2*Axx_new(i,j+2,k-2)
           Axy_new(i,j,k) = rp1*TWO*Axy_new(i,j+1,k-1) - rp2*Axy_new(i,j+2,k-2)
           Axz_new(i,j,k) = rp1*TWO*Axz_new(i,j+1,k-1) - rp2*Axz_new(i,j+2,k-2)
           Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i,j+1,k-1) - rp2*Ayy_new(i,j+2,k-2)
           Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i,j+1,k-1) - rp2*Ayz_new(i,j+2,k-2)
           Azz_new(i,j,k) = rp1*TWO*Azz_new(i,j+1,k-1) - rp2*Azz_new(i,j+2,k-2)
           Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i,j+1,k-1) - rp2*Gammax_new(i,j+2,k-2)
           Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i,j+1,k-1) - rp2*Gammay_new(i,j+2,k-2)
           Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i,j+1,k-1) - rp2*Gammaz_new(i,j+2,k-2)
        enddo
     endif
     !
     ! do all the z=zmin symmetries
     !
     k=kmin
     !
     ! x=xmax,z=zmin
     !
     if(have_bdry_max(1)==1.and.have_bdry_min(3)==1) then
        i=imax
        do j = jmin, jmax
           r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
           r1 = sqrt(X(i-1,j,k+1)**2 + Y(i-1,j,k+1)**2 + Z(i-1,j,k+1)**2)
           r2 = sqrt(X(i-2,j,k+2)**2 + Y(i-2,j,k+2)**2 + Z(i-2,j,k+2)**2)
           rp  = PhysicalR(i,j,k)
           rp1 = PhysicalR(i-1,j,k+1)
           rp2 = PhysicalR(i-2,j,k+2)
           fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
           fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i-1,j,k+1))
           fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i-2,j,k+2))
           rp1=rp1/rp
           rp2=rp2/rp
           phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i-1,j,k+1)-fac1)-rp2*(phi_new(i-2,j,k+2)-fac2))
           K_new(i,j,k)   = rp1*TWO*K_new(i-1,j,k+1) - rp2*K_new(i-2,j,k+2)
           gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i-1,j,k+1)-1.0d0) - rp2*(gxx_new(i-2,j,k+2)-1.0d0)
           gxy_new(i,j,k) = rp1*TWO*gxy_new(i-1,j,k+1) - rp2*gxy_new(i-2,j,k+2)
           gxz_new(i,j,k) = rp1*TWO*gxz_new(i-1,j,k+1) - rp2*gxz_new(i-2,j,k+2)
           gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i-1,j,k+1)-1.0d0) - rp2*(gyy_new(i-2,j,k+2)-1.0d0)
           gyz_new(i,j,k) = rp1*TWO*gyz_new(i-1,j,k+1) - rp2*gyz_new(i-2,j,k+2)
           gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i-1,j,k+1)-1.0d0) - rp2*(gzz_new(i-2,j,k+2)-1.0d0)
           Axx_new(i,j,k) = rp1*TWO*Axx_new(i-1,j,k+1) - rp2*Axx_new(i-2,j,k+2)
           Axy_new(i,j,k) = rp1*TWO*Axy_new(i-1,j,k+1) - rp2*Axy_new(i-2,j,k+2)
           Axz_new(i,j,k) = rp1*TWO*Axz_new(i-1,j,k+1) - rp2*Axz_new(i-2,j,k+2)
           Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i-1,j,k+1) - rp2*Ayy_new(i-2,j,k+2)
           Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i-1,j,k+1) - rp2*Ayz_new(i-2,j,k+2)
           Azz_new(i,j,k) = rp1*TWO*Azz_new(i-1,j,k+1) - rp2*Azz_new(i-2,j,k+2)
           Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i-1,j,k+1) - rp2*Gammax_new(i-2,j,k+2)
           Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i-1,j,k+1) - rp2*Gammay_new(i-2,j,k+2)
           Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i-1,j,k+1) - rp2*Gammaz_new(i-2,j,k+2)
        enddo
     endif
     !
     ! x=xmin,z=zmin
     !
     if(have_bdry_min(1)==1.and.have_bdry_min(3)==1) then
        i=imin
        do j = jmin, jmax
           r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
           r1 = sqrt(X(i+1,j,k+1)**2 + Y(i+1,j,k+1)**2 + Z(i+1,j,k+1)**2)
           r2 = sqrt(X(i+2,j,k+2)**2 + Y(i+2,j,k+2)**2 + Z(i+2,j,k+2)**2)
           rp  = PhysicalR(i,j,k)
           rp1 = PhysicalR(i+1,j,k+1)
           rp2 = PhysicalR(i+2,j,k+2)
           fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
           fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i+1,j,k+1))
           fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i+2,j,k+2))
           rp1=rp1/rp
           rp2=rp2/rp
           phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i+1,j,k+1)-fac1)-rp2*(phi_new(i+2,j,k+2)-fac2))
           K_new(i,j,k)   = rp1*TWO*K_new(i+1,j,k+1) - rp2*K_new(i+2,j,k+2)
           gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i+1,j,k+1)-1.0d0) - rp2*(gxx_new(i+2,j,k+2)-1.0d0)
           gxy_new(i,j,k) = rp1*TWO*gxy_new(i+1,j,k+1) - rp2*gxy_new(i+2,j,k+2)
           gxz_new(i,j,k) = rp1*TWO*gxz_new(i+1,j,k+1) - rp2*gxz_new(i+2,j,k+2)
           gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i+1,j,k+1)-1.0d0) - rp2*(gyy_new(i+2,j,k+2)-1.0d0)
           gyz_new(i,j,k) = rp1*TWO*gyz_new(i+1,j,k+1) - rp2*gyz_new(i+2,j,k+2)
           gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i+1,j,k+1)-1.0d0) - rp2*(gzz_new(i+2,j,k+2)-1.0d0)
           Axx_new(i,j,k) = rp1*TWO*Axx_new(i+1,j,k+1) - rp2*Axx_new(i+2,j,k+2)
           Axy_new(i,j,k) = rp1*TWO*Axy_new(i+1,j,k+1) - rp2*Axy_new(i+2,j,k+2)
           Axz_new(i,j,k) = rp1*TWO*Axz_new(i+1,j,k+1) - rp2*Axz_new(i+2,j,k+2)
           Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i+1,j,k+1) - rp2*Ayy_new(i+2,j,k+2)
           Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i+1,j,k+1) - rp2*Ayz_new(i+2,j,k+2)
           Azz_new(i,j,k) = rp1*TWO*Azz_new(i+1,j,k+1) - rp2*Azz_new(i+2,j,k+2)
           Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i+1,j,k+1) - rp2*Gammax_new(i+2,j,k+2)
           Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i+1,j,k+1) - rp2*Gammay_new(i+2,j,k+2)
           Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i+1,j,k+1) - rp2*Gammaz_new(i+2,j,k+2)
        enddo
     endif
     !
     ! y=ymax,z=zmin
     !
     if(have_bdry_max(2)==1.and.have_bdry_min(3)==1) then
        j = jmax
        do i = imin, imax
           r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
           r1 = sqrt(X(i,j-1,k+1)**2 + Y(i,j-1,k+1)**2 + Z(i,j-1,k+1)**2)
           r2 = sqrt(X(i,j-2,k+2)**2 + Y(i,j-2,k+2)**2 + Z(i,j-2,k+2)**2)
           rp  = PhysicalR(i,j,k)
           rp1 = PhysicalR(i,j-1,k+1)
           rp2 = PhysicalR(i,j-2,k+2)
           fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
           fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i,j-1,k+1))
           fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i,j-2,k+2))
           rp1=rp1/rp
           rp2=rp2/rp
           phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i,j-1,k+1)-fac1)-rp2*(phi_new(i,j-2,k+2)-fac2))
           K_new(i,j,k)   = rp1*TWO*K_new(i,j-1,k+1) - rp2*K_new(i,j-2,k+2)
           gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i,j-1,k+1)-1.0d0) - rp2*(gxx_new(i,j-2,k+2)-1.0d0)
           gxy_new(i,j,k) = rp1*TWO*gxy_new(i,j-1,k+1) - rp2*gxy_new(i,j-2,k+2)
           gxz_new(i,j,k) = rp1*TWO*gxz_new(i,j-1,k+1) - rp2*gxz_new(i,j-2,k+2)
           gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i,j-1,k+1)-1.0d0) - rp2*(gyy_new(i,j-2,k+2)-1.0d0)
           gyz_new(i,j,k) = rp1*TWO*gyz_new(i,j-1,k+1) - rp2*gyz_new(i,j-2,k+2)
           gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i,j-1,k+1)-1.0d0) - rp2*(gzz_new(i,j-2,k+2)-1.0d0)
           Axx_new(i,j,k) = rp1*TWO*Axx_new(i,j-1,k+1) - rp2*Axx_new(i,j-2,k+2)
           Axy_new(i,j,k) = rp1*TWO*Axy_new(i,j-1,k+1) - rp2*Axy_new(i,j-2,k+2)
           Axz_new(i,j,k) = rp1*TWO*Axz_new(i,j-1,k+1) - rp2*Axz_new(i,j-2,k+2)
           Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i,j-1,k+1) - rp2*Ayy_new(i,j-2,k+2)
           Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i,j-1,k+1) - rp2*Ayz_new(i,j-2,k+2)
           Azz_new(i,j,k) = rp1*TWO*Azz_new(i,j-1,k+1) - rp2*Azz_new(i,j-2,k+2)
           Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i,j-1,k+1) - rp2*Gammax_new(i,j-2,k+2)
           Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i,j-1,k+1) - rp2*Gammay_new(i,j-2,k+2)
           Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i,j-1,k+1) - rp2*Gammaz_new(i,j-2,k+2)
        enddo
     endif
     !
     ! y=ymin,z=zmin
     !
     if(have_bdry_min(2)==1.and.have_bdry_min(3)==1) then
        j = jmin
        do i = imin, imax
           r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
           r1 = sqrt(X(i,j+1,k+1)**2 + Y(i,j+1,k+1)**2 + Z(i,j+1,k+1)**2)
           r2 = sqrt(X(i,j+2,k+2)**2 + Y(i,j+2,k+2)**2 + Z(i,j+2,k+2)**2)
           rp  = PhysicalR(i,j,k)
           rp1 = PhysicalR(i,j+1,k+1)
           rp2 = PhysicalR(i,j+2,k+2)
           fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
           fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i,j+1,k+1))
           fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i,j+2,k+2))
           rp1=rp1/rp
           rp2=rp2/rp
           phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i,j+1,k+1)-fac1)-rp2*(phi_new(i,j+2,k+2)-fac2))
           K_new(i,j,k)   = rp1*TWO*K_new(i,j+1,k+1) - rp2*K_new(i,j+2,k+2)
           gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i,j+1,k+1)-1.0d0) - rp2*(gxx_new(i,j+2,k+2)-1.0d0)
           gxy_new(i,j,k) = rp1*TWO*gxy_new(i,j+1,k+1) - rp2*gxy_new(i,j+2,k+2)
           gxz_new(i,j,k) = rp1*TWO*gxz_new(i,j+1,k+1) - rp2*gxz_new(i,j+2,k+2)
           gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i,j+1,k+1)-1.0d0) - rp2*(gyy_new(i,j+2,k+2)-1.0d0)
           gyz_new(i,j,k) = rp1*TWO*gyz_new(i,j+1,k+1) - rp2*gyz_new(i,j+2,k+2)
           gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i,j+1,k+1)-1.0d0) - rp2*(gzz_new(i,j+2,k+2)-1.0d0)
           Axx_new(i,j,k) = rp1*TWO*Axx_new(i,j+1,k+1) - rp2*Axx_new(i,j+2,k+2)
           Axy_new(i,j,k) = rp1*TWO*Axy_new(i,j+1,k+1) - rp2*Axy_new(i,j+2,k+2)
           Axz_new(i,j,k) = rp1*TWO*Axz_new(i,j+1,k+1) - rp2*Axz_new(i,j+2,k+2)
           Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i,j+1,k+1) - rp2*Ayy_new(i,j+2,k+2)
           Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i,j+1,k+1) - rp2*Ayz_new(i,j+2,k+2)
           Azz_new(i,j,k) = rp1*TWO*Azz_new(i,j+1,k+1) - rp2*Azz_new(i,j+2,k+2)
           Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i,j+1,k+1) - rp2*Gammax_new(i,j+2,k+2)
           Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i,j+1,k+1) - rp2*Gammay_new(i,j+2,k+2)
           Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i,j+1,k+1) - rp2*Gammaz_new(i,j+2,k+2)
        enddo
     endif
     !
     !     finally the corners
     !
     !     x=xmax,y=ymax,z=zmax
     !
     if(have_bdry_max(1)==1.and.have_bdry_max(2)==1.and.have_bdry_max(3)==1) then
        i = imax
        j = jmax
        k = kmax
        r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
        r1 = sqrt(X(i-1,j-1,k-1)**2 + Y(i-1,j-1,k-1)**2 + Z(i-1,j-1,k-1)**2)
        r2 = sqrt(X(i-2,j-2,k-2)**2 + Y(i-2,j-2,k-2)**2 + Z(i-2,j-2,k-2)**2)
        rp  = PhysicalR(i,j,k)
        rp1 = PhysicalR(i-1,j-1,k-1)
        rp2 = PhysicalR(i-2,j-2,k-2)
        fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
        fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i-1,j-1,k-1))
        fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i-2,j-2,k-2))
        rp1=rp1/rp
        rp2=rp2/rp
        phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i-1,j-1,k-1)-fac1)-rp2*(phi_new(i-2,j-2,k-2)-fac2))
        K_new(i,j,k)   = rp1*TWO*K_new(i-1,j-1,k-1) - rp2*K_new(i-2,j-2,k-2)
        gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i-1,j-1,k-1)-1.0d0) - rp2*(gxx_new(i-2,j-2,k-2)-1.0d0)
        gxy_new(i,j,k) = rp1*TWO*gxy_new(i-1,j-1,k-1) - rp2*gxy_new(i-2,j-2,k-2)
        gxz_new(i,j,k) = rp1*TWO*gxz_new(i-1,j-1,k-1) - rp2*gxz_new(i-2,j-2,k-2)
        gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i-1,j-1,k-1)-1.0d0) - rp2*(gyy_new(i-2,j-2,k-2)-1.0d0)
        gyz_new(i,j,k) = rp1*TWO*gyz_new(i-1,j-1,k-1) - rp2*gyz_new(i-2,j-2,k-2)
        gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i-1,j-1,k-1)-1.0d0) - rp2*(gzz_new(i-2,j-2,k-2)-1.0d0)
        Axx_new(i,j,k) = rp1*TWO*Axx_new(i-1,j-1,k-1) - rp2*Axx_new(i-2,j-2,k-2)
        Axy_new(i,j,k) = rp1*TWO*Axy_new(i-1,j-1,k-1) - rp2*Axy_new(i-2,j-2,k-2)
        Axz_new(i,j,k) = rp1*TWO*Axz_new(i-1,j-1,k-1) - rp2*Axz_new(i-2,j-2,k-2)
        Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i-1,j-1,k-1) - rp2*Ayy_new(i-2,j-2,k-2)
        Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i-1,j-1,k-1) - rp2*Ayz_new(i-2,j-2,k-2)
        Azz_new(i,j,k) = rp1*TWO*Azz_new(i-1,j-1,k-1) - rp2*Azz_new(i-2,j-2,k-2)
        Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i-1,j-1,k-1) - rp2*Gammax_new(i-2,j-2,k-2)
        Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i-1,j-1,k-1) - rp2*Gammay_new(i-2,j-2,k-2)
        Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i-1,j-1,k-1) - rp2*Gammaz_new(i-2,j-2,k-2)
     endif
     !  endif
     !     
     !    x=xmin,y=ymax,z=zmax
     !
     if(have_bdry_min(1)==1.and.have_bdry_max(2)==1.and.have_bdry_max(3)==1) then
        i = imin
        j = jmax
        k = kmax
        r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
        r1 = sqrt(X(i+1,j-1,k-1)**2 + Y(i+1,j-1,k-1)**2 + Z(i+1,j-1,k-1)**2)
        r2 = sqrt(X(i+2,j-2,k-2)**2 + Y(i+2,j-2,k-2)**2 + Z(i+2,j-2,k-2)**2)
        rp  = PhysicalR(i,j,k)
        rp1 = PhysicalR(i+1,j-1,k-1)
        rp2 = PhysicalR(i+2,j-2,k-2)
        fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
        fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i+1,j-1,k-1))
        fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i+2,j-2,k-2))
        rp1=rp1/rp
        rp2=rp2/rp
        phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i+1,j-1,k-1)-fac1)-rp2*(phi_new(i+2,j-2,k-2)-fac2))
        K_new(i,j,k)   = rp1*TWO*K_new(i+1,j-1,k-1) - rp2*K_new(i+2,j-2,k-2)
        gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i+1,j-1,k-1)-1.0d0) - rp2*(gxx_new(i+2,j-2,k-2)-1.0d0)
        gxy_new(i,j,k) = rp1*TWO*gxy_new(i+1,j-1,k-1) - rp2*gxy_new(i+2,j-2,k-2)
        gxz_new(i,j,k) = rp1*TWO*gxz_new(i+1,j-1,k-1) - rp2*gxz_new(i+2,j-2,k-2)
        gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i+1,j-1,k-1)-1.0d0) - rp2*(gyy_new(i+2,j-2,k-2)-1.0d0)
        gyz_new(i,j,k) = rp1*TWO*gyz_new(i+1,j-1,k-1) - rp2*gyz_new(i+2,j-2,k-2)
        gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i+1,j-1,k-1)-1.0d0) - rp2*(gzz_new(i+2,j-2,k-2)-1.0d0)
        Axx_new(i,j,k) = rp1*TWO*Axx_new(i+1,j-1,k-1) - rp2*Axx_new(i+2,j-2,k-2)
        Axy_new(i,j,k) = rp1*TWO*Axy_new(i+1,j-1,k-1) - rp2*Axy_new(i+2,j-2,k-2)
        Axz_new(i,j,k) = rp1*TWO*Axz_new(i+1,j-1,k-1) - rp2*Axz_new(i+2,j-2,k-2)
        Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i+1,j-1,k-1) - rp2*Ayy_new(i+2,j-2,k-2)
        Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i+1,j-1,k-1) - rp2*Ayz_new(i+2,j-2,k-2)
        Azz_new(i,j,k) = rp1*TWO*Azz_new(i+1,j-1,k-1) - rp2*Azz_new(i+2,j-2,k-2)
        Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i+1,j-1,k-1) - rp2*Gammax_new(i+2,j-2,k-2)
        Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i+1,j-1,k-1) - rp2*Gammay_new(i+2,j-2,k-2)
        Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i+1,j-1,k-1) - rp2*Gammaz_new(i+2,j-2,k-2)
     endif
     !
     !    x=xmax,y=ymin,z=zmax
     !
     if(have_bdry_max(1)==1.and.have_bdry_min(2)==1.and.have_bdry_max(3)==1) then
        i = imax
        j = jmin
        k = kmax
        r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
        r1 = sqrt(X(i-1,j+1,k-1)**2 + Y(i-1,j+1,k-1)**2 + Z(i-1,j+1,k-1)**2)
        r2 = sqrt(X(i-2,j+2,k-2)**2 + Y(i-2,j+2,k-2)**2 + Z(i-2,j+2,k-2)**2)
        rp  = PhysicalR(i,j,k)
        rp1 = PhysicalR(i-1,j+1,k-1)
        rp2 = PhysicalR(i-2,j+2,k-2)
        fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
        fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i-1,j+1,k-1))
        fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i-2,j+2,k-2))
        rp1=rp1/rp
        rp2=rp2/rp
        phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i-1,j+1,k-1)-fac1)-rp2*(phi_new(i-2,j+2,k-2)-fac2))
        K_new(i,j,k)   = rp1*TWO*K_new(i-1,j+1,k-1) - rp2*K_new(i-2,j+2,k-2)
        gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i-1,j+1,k-1)-1.0d0) - rp2*(gxx_new(i-2,j+2,k-2)-1.0d0)
        gxy_new(i,j,k) = rp1*TWO*gxy_new(i-1,j+1,k-1) - rp2*gxy_new(i-2,j+2,k-2)
        gxz_new(i,j,k) = rp1*TWO*gxz_new(i-1,j+1,k-1) - rp2*gxz_new(i-2,j+2,k-2)
        gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i-1,j+1,k-1)-1.0d0) - rp2*(gyy_new(i-2,j+2,k-2)-1.0d0)
        gyz_new(i,j,k) = rp1*TWO*gyz_new(i-1,j+1,k-1) - rp2*gyz_new(i-2,j+2,k-2)
        gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i-1,j+1,k-1)-1.0d0) - rp2*(gzz_new(i-2,j+2,k-2)-1.0d0)
        Axx_new(i,j,k) = rp1*TWO*Axx_new(i-1,j+1,k-1) - rp2*Axx_new(i-2,j+2,k-2)
        Axy_new(i,j,k) = rp1*TWO*Axy_new(i-1,j+1,k-1) - rp2*Axy_new(i-2,j+2,k-2)
        Axz_new(i,j,k) = rp1*TWO*Axz_new(i-1,j+1,k-1) - rp2*Axz_new(i-2,j+2,k-2)
        Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i-1,j+1,k-1) - rp2*Ayy_new(i-2,j+2,k-2)
        Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i-1,j+1,k-1) - rp2*Ayz_new(i-2,j+2,k-2)
        Azz_new(i,j,k) = rp1*TWO*Azz_new(i-1,j+1,k-1) - rp2*Azz_new(i-2,j+2,k-2)
        Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i-1,j+1,k-1) - rp2*Gammax_new(i-2,j+2,k-2)
        Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i-1,j+1,k-1) - rp2*Gammay_new(i-2,j+2,k-2)
        Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i-1,j+1,k-1) - rp2*Gammaz_new(i-2,j+2,k-2)
     endif
     !
     !    x=xmin,y=ymin,z=zmax
     !
     if(have_bdry_min(1)==1.and.have_bdry_min(2)==1.and.have_bdry_max(3)==1) then
        i = imin
        r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
        r1 = sqrt(X(i+1,j+1,k-1)**2 + Y(i+1,j+1,k-1)**2 + Z(i+1,j+1,k-1)**2)
        r2 = sqrt(X(i+2,j+2,k-2)**2 + Y(i+2,j+2,k-2)**2 + Z(i+2,j+2,k-2)**2)
        rp  = PhysicalR(i,j,k)
        rp1 = PhysicalR(i+1,j+1,k-1)
        rp2 = PhysicalR(i+2,j+2,k-2)
        fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
        fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i+1,j+1,k-1))
        fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i+2,j+2,k-2))
        rp1=rp1/rp
        rp2=rp2/rp
        phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i+1,j+1,k-1)-fac1)-rp2*(phi_new(i+2,j+2,k-2)-fac2))
        K_new(i,j,k)   = rp1*TWO*K_new(i+1,j+1,k-1) - rp2*K_new(i+2,j+2,k-2)
        gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i+1,j+1,k-1)-1.0d0) - rp2*(gxx_new(i+2,j+2,k-2)-1.0d0)
        gxy_new(i,j,k) = rp1*TWO*gxy_new(i+1,j+1,k-1) - rp2*gxy_new(i+2,j+2,k-2)
        gxz_new(i,j,k) = rp1*TWO*gxz_new(i+1,j+1,k-1) - rp2*gxz_new(i+2,j+2,k-2)
        gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i+1,j+1,k-1)-1.0d0) - rp2*(gyy_new(i+2,j+2,k-2)-1.0d0)
        gyz_new(i,j,k) = rp1*TWO*gyz_new(i+1,j+1,k-1) - rp2*gyz_new(i+2,j+2,k-2)
        gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i+1,j+1,k-1)-1.0d0) - rp2*(gzz_new(i+2,j+2,k-2)-1.0d0)
        Axx_new(i,j,k) = rp1*TWO*Axx_new(i+1,j+1,k-1) - rp2*Axx_new(i+2,j+2,k-2)
        Axy_new(i,j,k) = rp1*TWO*Axy_new(i+1,j+1,k-1) - rp2*Axy_new(i+2,j+2,k-2)
        Axz_new(i,j,k) = rp1*TWO*Axz_new(i+1,j+1,k-1) - rp2*Axz_new(i+2,j+2,k-2)
        Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i+1,j+1,k-1) - rp2*Ayy_new(i+2,j+2,k-2)
        Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i+1,j+1,k-1) - rp2*Ayz_new(i+2,j+2,k-2)
        Azz_new(i,j,k) = rp1*TWO*Azz_new(i+1,j+1,k-1) - rp2*Azz_new(i+2,j+2,k-2)
        Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i+1,j+1,k-1) - rp2*Gammax_new(i+2,j+2,k-2)
        Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i+1,j+1,k-1) - rp2*Gammay_new(i+2,j+2,k-2)
        Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i+1,j+1,k-1) - rp2*Gammaz_new(i+2,j+2,k-2)
     endif
     !
     !    z = zmin symmetries
     !
     k=kmin
     !
     !    x=xmax,y=ymax,z=zmin
     !
     if(have_bdry_max(1)==1.and.have_bdry_max(2)==1.and.have_bdry_min(3)==1) then
        i=imax
        j=jmax
        r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
        r1 = sqrt(X(i-1,j-1,k+1)**2 + Y(i-1,j-1,k+1)**2 + Z(i-1,j-1,k+1)**2)
        r2 = sqrt(X(i-2,j-2,k+2)**2 + Y(i-2,j-2,k+2)**2 + Z(i-2,j-2,k+2)**2)
        rp  = PhysicalR(i,j,k)
        rp1 = PhysicalR(i-1,j-1,k+1)
        rp2 = PhysicalR(i-2,j-2,k+2)
        fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
        fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i-1,j-1,k+1))
        fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i-2,j-2,k+2))
        rp1=rp1/rp
        rp2=rp2/rp
        phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i-1,j-1,k+1)-fac1)-rp2*(phi_new(i-2,j-2,k+2)-fac2))
        K_new(i,j,k)   = rp1*TWO*K_new(i-1,j-1,k+1) - rp2*K_new(i-2,j-2,k+2)
        gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i-1,j-1,k+1)-1.0d0) - rp2*(gxx_new(i-2,j-2,k+2)-1.0d0)
        gxy_new(i,j,k) = rp1*TWO*gxy_new(i-1,j-1,k+1) - rp2*gxy_new(i-2,j-2,k+2)
        gxz_new(i,j,k) = rp1*TWO*gxz_new(i-1,j-1,k+1) - rp2*gxz_new(i-2,j-2,k+2)
        gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i-1,j-1,k+1)-1.0d0) - rp2*(gyy_new(i-2,j-2,k+2)-1.0d0)
        gyz_new(i,j,k) = rp1*TWO*gyz_new(i-1,j-1,k+1) - rp2*gyz_new(i-2,j-2,k+2)
        gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i-1,j-1,k+1)-1.0d0) - rp2*(gzz_new(i-2,j-2,k+2)-1.0d0)
        Axx_new(i,j,k) = rp1*TWO*Axx_new(i-1,j-1,k+1) - rp2*Axx_new(i-2,j-2,k+2)
        Axy_new(i,j,k) = rp1*TWO*Axy_new(i-1,j-1,k+1) - rp2*Axy_new(i-2,j-2,k+2)
        Axz_new(i,j,k) = rp1*TWO*Axz_new(i-1,j-1,k+1) - rp2*Axz_new(i-2,j-2,k+2)
        Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i-1,j-1,k+1) - rp2*Ayy_new(i-2,j-2,k+2)
        Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i-1,j-1,k+1) - rp2*Ayz_new(i-2,j-2,k+2)
        Azz_new(i,j,k) = rp1*TWO*Azz_new(i-1,j-1,k+1) - rp2*Azz_new(i-2,j-2,k+2)
        Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i-1,j-1,k+1) - rp2*Gammax_new(i-2,j-2,k+2)
        Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i-1,j-1,k+1) - rp2*Gammay_new(i-2,j-2,k+2)
        Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i-1,j-1,k+1) - rp2*Gammaz_new(i-2,j-2,k+2)
     endif
     !
     !    x=xmin,y=ymax,z=zmin
     !
     if(have_bdry_min(1)==1.and.have_bdry_max(2)==1.and.have_bdry_min(3)==1) then
        i=imin
        r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
        r1 = sqrt(X(i+1,j-1,k+1)**2 + Y(i+1,j-1,k+1)**2 + Z(i+1,j-1,k+1)**2)
        r2 = sqrt(X(i+2,j-2,k+2)**2 + Y(i+2,j-2,k+2)**2 + Z(i+2,j-2,k+2)**2)
        rp  = PhysicalR(i,j,k)
        rp1 = PhysicalR(i+1,j-1,k+1)
        rp2 = PhysicalR(i+2,j-2,k+2)
        fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
        fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i+1,j-1,k+1))
        fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i+2,j-2,k+2))
        rp1=rp1/rp
        rp2=rp2/rp
        phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i+1,j-1,k+1)-fac1)-rp2*(phi_new(i+2,j-2,k+2)-fac2))
        K_new(i,j,k)   = rp1*TWO*K_new(i+1,j-1,k+1) - rp2*K_new(i+2,j-2,k+2)
        gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i+1,j-1,k+1)-1.0d0) - rp2*(gxx_new(i+2,j-2,k+2)-1.0d0)
        gxy_new(i,j,k) = rp1*TWO*gxy_new(i+1,j-1,k+1) - rp2*gxy_new(i+2,j-2,k+2)
        gxz_new(i,j,k) = rp1*TWO*gxz_new(i+1,j-1,k+1) - rp2*gxz_new(i+2,j-2,k+2)
        gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i+1,j-1,k+1)-1.0d0) - rp2*(gyy_new(i+2,j-2,k+2)-1.0d0)
        gyz_new(i,j,k) = rp1*TWO*gyz_new(i+1,j-1,k+1) - rp2*gyz_new(i+2,j-2,k+2)
        gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i+1,j-1,k+1)-1.0d0) - rp2*(gzz_new(i+2,j-2,k+2)-1.0d0)
        Axx_new(i,j,k) = rp1*TWO*Axx_new(i+1,j-1,k+1) - rp2*Axx_new(i+2,j-2,k+2)
        Axy_new(i,j,k) = rp1*TWO*Axy_new(i+1,j-1,k+1) - rp2*Axy_new(i+2,j-2,k+2)
        Axz_new(i,j,k) = rp1*TWO*Axz_new(i+1,j-1,k+1) - rp2*Axz_new(i+2,j-2,k+2)
        Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i+1,j-1,k+1) - rp2*Ayy_new(i+2,j-2,k+2)
        Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i+1,j-1,k+1) - rp2*Ayz_new(i+2,j-2,k+2)
        Azz_new(i,j,k) = rp1*TWO*Azz_new(i+1,j-1,k+1) - rp2*Azz_new(i+2,j-2,k+2)
        Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i+1,j-1,k+1) - rp2*Gammax_new(i+2,j-2,k+2)
        Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i+1,j-1,k+1) - rp2*Gammay_new(i+2,j-2,k+2)
        Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i+1,j-1,k+1) - rp2*Gammaz_new(i+2,j-2,k+2)
     endif
     !
     !    x=xmax,y=ymin,z=zmin
     !
     if(have_bdry_max(1)==1.and.have_bdry_min(2)==1.and.have_bdry_min(3)==1) then
        i=imax
        j=jmin
        r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
        r1 = sqrt(X(i-1,j+1,k+1)**2 + Y(i-1,j+1,k+1)**2 + Z(i-1,j+1,k+1)**2)
        r2 = sqrt(X(i-2,j+2,k+2)**2 + Y(i-2,j+2,k+2)**2 + Z(i-2,j+2,k+2)**2)
        rp  = PhysicalR(i,j,k)
        rp1 = PhysicalR(i-1,j+1,k+1)
        rp2 = PhysicalR(i-2,j+2,k+2)
        fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
        fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i-1,j+1,k+1))
        fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i-2,j+2,k+2))
        rp1=rp1/rp
        rp2=rp2/rp
        phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i-1,j+1,k+1)-fac1)-rp2*(phi_new(i-2,j+2,k+2)-fac2))
        K_new(i,j,k)   = rp1*TWO*K_new(i-1,j+1,k+1) - rp2*K_new(i-2,j+2,k+2)
        gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i-1,j+1,k+1)-1.0d0) - rp2*(gxx_new(i-2,j+2,k+2)-1.0d0)
        gxy_new(i,j,k) = rp1*TWO*gxy_new(i-1,j+1,k+1) - rp2*gxy_new(i-2,j+2,k+2)
        gxz_new(i,j,k) = rp1*TWO*gxz_new(i-1,j+1,k+1) - rp2*gxz_new(i-2,j+2,k+2)
        gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i-1,j+1,k+1)-1.0d0) - rp2*(gyy_new(i-2,j+2,k+2)-1.0d0)
        gyz_new(i,j,k) = rp1*TWO*gyz_new(i-1,j+1,k+1) - rp2*gyz_new(i-2,j+2,k+2)
        gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i-1,j+1,k+1)-1.0d0) - rp2*(gzz_new(i-2,j+2,k+2)-1.0d0)
        Axx_new(i,j,k) = rp1*TWO*Axx_new(i-1,j+1,k+1) - rp2*Axx_new(i-2,j+2,k+2)
        Axy_new(i,j,k) = rp1*TWO*Axy_new(i-1,j+1,k+1) - rp2*Axy_new(i-2,j+2,k+2)
        Axz_new(i,j,k) = rp1*TWO*Axz_new(i-1,j+1,k+1) - rp2*Axz_new(i-2,j+2,k+2)
        Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i-1,j+1,k+1) - rp2*Ayy_new(i-2,j+2,k+2)
        Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i-1,j+1,k+1) - rp2*Ayz_new(i-2,j+2,k+2)
        Azz_new(i,j,k) = rp1*TWO*Azz_new(i-1,j+1,k+1) - rp2*Azz_new(i-2,j+2,k+2)
        Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i-1,j+1,k+1) - rp2*Gammax_new(i-2,j+2,k+2)
        Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i-1,j+1,k+1) - rp2*Gammay_new(i-2,j+2,k+2)
        Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i-1,j+1,k+1) - rp2*Gammaz_new(i-2,j+2,k+2)
     endif
     !
     !    x=xmin,y=ymin,z=zmin
     !
     if(have_bdry_min(1)==1.and.have_bdry_min(2)==1.and.have_bdry_min(3)==1) then
        i=imin
        r = sqrt(X(i,j,k)**2 + Y(i,j,k)**2 + Z(i,j,k)**2)
        r1 = sqrt(X(i+1,j+1,k+1)**2 + Y(i+1,j+1,k+1)**2 + Z(i+1,j+1,k+1)**2)
        r2 = sqrt(X(i+2,j+2,k+2)**2 + Y(i+2,j+2,k+2)**2 + Z(i+2,j+2,k+2)**2)
        rp  = PhysicalR(i,j,k)
        rp1 = PhysicalR(i+1,j+1,k+1)
        rp2 = PhysicalR(i+2,j+2,k+2)
        fac = f1o3*log(rp/r)+f1o6*log(RadiusDerivative(i,j,k))
        fac1 = f1o3*log(rp1/r1)+f1o6*log(RadiusDerivative(i+1,j+1,k+1))
        fac2 = f1o3*log(rp2/r2)+f1o6*log(RadiusDerivative(i+2,j+2,k+2))
        rp1=rp1/rp
        rp2=rp2/rp
        phi_new(i,j,k) = fac+(rp1*TWO*(phi_new(i+1,j+1,k+1)-fac1)-rp2*(phi_new(i+2,j+2,k+2)-fac2))
        K_new(i,j,k)   = rp1*TWO*K_new(i+1,j+1,k+1) - rp2*K_new(i+2,j+2,k+2)
        gxx_new(i,j,k) = 1.0d0+ rp1*TWO*(gxx_new(i+1,j+1,k+1)-1.0d0) - rp2*(gxx_new(i+2,j+2,k+2)-1.0d0)
        gxy_new(i,j,k) = rp1*TWO*gxy_new(i+1,j+1,k+1) - rp2*gxy_new(i+2,j+2,k+2)
        gxz_new(i,j,k) = rp1*TWO*gxz_new(i+1,j+1,k+1) - rp2*gxz_new(i+2,j+2,k+2)
        gyy_new(i,j,k) = 1.0d0+ rp1*TWO*(gyy_new(i+1,j+1,k+1)-1.0d0) - rp2*(gyy_new(i+2,j+2,k+2)-1.0d0)
        gyz_new(i,j,k) = rp1*TWO*gyz_new(i+1,j+1,k+1) - rp2*gyz_new(i+2,j+2,k+2)
        gzz_new(i,j,k) = 1.0d0+ rp1*TWO*(gzz_new(i+1,j+1,k+1)-1.0d0) - rp2*(gzz_new(i+2,j+2,k+2)-1.0d0)
        Axx_new(i,j,k) = rp1*TWO*Axx_new(i+1,j+1,k+1) - rp2*Axx_new(i+2,j+2,k+2)
        Axy_new(i,j,k) = rp1*TWO*Axy_new(i+1,j+1,k+1) - rp2*Axy_new(i+2,j+2,k+2)
        Axz_new(i,j,k) = rp1*TWO*Axz_new(i+1,j+1,k+1) - rp2*Axz_new(i+2,j+2,k+2)
        Ayy_new(i,j,k) = rp1*TWO*Ayy_new(i+1,j+1,k+1) - rp2*Ayy_new(i+2,j+2,k+2)
        Ayz_new(i,j,k) = rp1*TWO*Ayz_new(i+1,j+1,k+1) - rp2*Ayz_new(i+2,j+2,k+2)
        Azz_new(i,j,k) = rp1*TWO*Azz_new(i+1,j+1,k+1) - rp2*Azz_new(i+2,j+2,k+2)
        Gammax_new(i,j,k) = rp1*TWO*Gammax_new(i+1,j+1,k+1) - rp2*Gammax_new(i+2,j+2,k+2)
        Gammay_new(i,j,k) = rp1*TWO*Gammay_new(i+1,j+1,k+1) - rp2*Gammay_new(i+2,j+2,k+2)
        Gammaz_new(i,j,k) = rp1*TWO*Gammaz_new(i+1,j+1,k+1) - rp2*Gammaz_new(i+2,j+2,k+2)
     endif
     if (Fisheye_Flag==1) then
        fish_to_phys = .FALSE.
        ! FIXME: Shouldnt the ,2, below be ,0,????!!!  -Zach
        call trans_fish_phys(ex, fake_ex, X, Y, Z, PhysicalR, RadiusDerivative, &
             gxx_new,gxy_new,gxz_new,gyy_new,gyz_new,gzz_new, &
             fish_to_phys,2,have_bdry_min,have_bdry_max,use_trans_fish_phys_new)
        call trans_fish_phys(ex, fake_ex, X, Y, Z, PhysicalR, RadiusDerivative, &
             Axx_new,Axy_new,Axz_new,Ayy_new,Ayz_new,Azz_new, &
             fish_to_phys,2,have_bdry_min,have_bdry_max,use_trans_fish_phys_new)
     end if
  endif
if(1==0) then
do k=1, ex(3)
     do j=1, ex(2)
        do i=1, ex(1)
                if(isnan(Axx_old(i,j,k))) then
                    write(*,*) "bssn update bdy. F90:  Checkpoint << 2 >>, Axx_old = NAN FOUND!"
                end if
                if(isnan(Axx_new(i,j,k))) then
                    write(*,*) "bssn update bdy. F90:  Checkpoint << 2 >>, Axx = NAN FOUND! i,j,k =", i,j,k
                end if
       end do
     end do
  end do
end if
end subroutine update_boundary
!
!-------------------------------------------------------------------
! The following is a subroutine to rotate things off of the plane
! when evolving in axisymmetry (x<xmax)
!-------------------------------------------------------------------
!
subroutine rotate(ex,X,Y,Z,alpha,phi_new,K_new,gxx_new,&
     &gxy_new,gxz_new,gyy_new,&
     &gyz_new,gzz_new,Axx_new,&
     &Axy_new,Axz_new,Ayy_new,&
     &Ayz_new,Azz_new, Gammax_new,&
     &Gammay_new, Gammaz_new)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: alpha,phi_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: K_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: gxx_new,gxy_new,gxz_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: gyy_new,gyz_new,gzz_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: Axx_new,Axy_new,Axz_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: Ayy_new,Ayz_new,Azz_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: Gammax_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: Gammay_new, Gammaz_new
  real*8                     :: f_int,f_int_phi,f_int_K,f_int_Gamx&
       &,f_int_Gamy,f_int_Gamz
  real*8                     :: f_int_gxx,f_int_gxy,f_int_gxz&
       &,f_int_gyy,f_int_gyz,f_int_gzz
  real*8                     :: f_int_Axx,f_int_Axy,f_int_Axz&
       &,f_int_Ayy,f_int_Ayz,f_int_Azz
  real*8                     :: HALF, ONE, TWO, ZERO
  parameter(HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0)
  integer                                  :: i,j,k
  real*8                                   :: cosine, sine, rho, err
  integer                                  :: imin, imax, jmin, jmax
  integer                                  :: kmin, kmax
  real*8                                   :: dX,dY,dZ
  integer                                  :: im, l
  real*8, dimension(4)                     :: xa,ya1,ya2,ya3,ya4,ya5,ya6
  imin = lbound(gxx_new,1)
  jmin = lbound(gxx_new,2)
  kmin = lbound(gxx_new,3)
  imax = ubound(gxx_new,1)
  jmax = ubound(gxx_new,2)
  kmax = ubound(gxx_new,3)
  dX = X(imin+1,1,1) - X(imin,1,1)
  dY = Y(1,jmin+1,1) - Y(1,jmin,1)
  dZ = Z(1,1,kmin+1) - Z(1,1,kmin)
  j=jmin + 1
  do k=kmin, kmax
     do i=imin, imax-1
        rho = sqrt(X(i,j,k)*X(i,j,k) + dY*dY)
        cosine = X(i,j,k)/rho
        sine = dY/rho
        if (i-imin < 2) then
           im=imin
        elseif(imax-i < 3) then
           im=imax-3
        else
           im =i-1
        end if
        ! Handle the scalars first
        if (i==imin) then
           xa(1) = X(im,j,k)-dX
           ya1(1) = phi_new(im,j,k)
           ya2(1) = alpha(im,j,k)
           ya3(1) = K_new(im,j,k)
           do l=2,4
              xa(l) = X(im+l-2,j,k)
              ya1(l) = phi_new(im+l-2,j,k)
              ya2(l) = alpha(im+l-2,j,k)
              ya3(l) = K_new(im+l-2,j,k)
           enddo
        else
           do l=1,4
              xa(l) = X(im+l-1,j,k)
              ya1(l) = phi_new(im+l-1,j,k)
              ya2(l) = alpha(im+l-1,j,k)
              ya3(l) = K_new(im+l-1,j,k)
           enddo
        endif
        call polint(xa,ya1,4,rho,f_int_phi,err)
        phi_new(i,j+1,k) = f_int_phi
        phi_new(i,j-1,k) = f_int_phi
        call polint(xa,ya2,4,rho,f_int,err)
        alpha(i,j+1,k) = f_int
        alpha(i,j-1,k) = f_int
        call polint(xa,ya3,4,rho,f_int_K,err)
        K_new(i,j+1,k) = f_int_K
        K_new(i,j-1,k) = f_int_K
        !Get interpolated values of metric components
        if (i==imin) then
           ya1(1) = gxx_new(im,j,k)
           ya2(1) = -gxy_new(im,j,k)
           ya3(1) = -gxz_new(im,j,k)
           ya4(1) = gyy_new(im,j,k)
           ya5(1) = gyz_new(im,j,k)
           ya6(1) = gzz_new(im,j,k)
           do l=2,4
              ya1(l) = gxx_new(im+l-2,j,k)
              ya2(l) = gxy_new(im+l-2,j,k)
              ya3(l) = gxz_new(im+l-2,j,k)
              ya4(l) = gyy_new(im+l-2,j,k)
              ya5(l) = gyz_new(im+l-2,j,k)
              ya6(l) = gzz_new(im+l-2,j,k)
           enddo
        else
           do l=1,4
              ya1(l) = gxx_new(im+l-1,j,k)
              ya2(l) = gxy_new(im+l-1,j,k)
              ya3(l) = gxz_new(im+l-1,j,k)
              ya4(l) = gyy_new(im+l-1,j,k)
              ya5(l) = gyz_new(im+l-1,j,k)
              ya6(l) = gzz_new(im+l-1,j,k)
           enddo
        endif
        call polint(xa,ya1,4,rho,f_int_gxx,err)
        call polint(xa,ya2,4,rho,f_int_gxy,err)
        call polint(xa,ya3,4,rho,f_int_gxz,err)
        call polint(xa,ya4,4,rho,f_int_gyy,err)
        call polint(xa,ya5,4,rho,f_int_gyz,err)
        call polint(xa,ya6,4,rho,f_int_gzz,err)
        !Now rotate the metric off of the plane
        gxx_new(i,j+1,k) = cosine*cosine*f_int_gxx - TWO*sine*cosine*f_int_gxy &
             & + sine*sine*f_int_gyy
        gxx_new(i,j-1,k) = cosine*cosine*f_int_gxx + TWO*sine*cosine*f_int_gxy &
             & + sine*sine*f_int_gyy
        gxy_new(i,j+1,k) = cosine*sine*f_int_gxx + (cosine*cosine-sine*sine)*f_int_gxy &
             & - sine*cosine*f_int_gyy
        gxy_new(i,j-1,k) = -cosine*sine*f_int_gxx + (cosine*cosine-sine*sine)*f_int_gxy &
             & + sine*cosine*f_int_gyy
        gxz_new(i,j+1,k) = cosine*f_int_gxz - sine*f_int_gyz
        gxz_new(i,j-1,k) = cosine*f_int_gxz + sine*f_int_gyz
        gyy_new(i,j+1,k) = sine*sine*f_int_gxx + TWO*sine*cosine*f_int_gxy &
             & + cosine*cosine*f_int_gyy
        gyy_new(i,j-1,k) = sine*sine*f_int_gxx - TWO*sine*cosine*f_int_gxy &
             & + cosine*cosine*f_int_gyy
        gyz_new(i,j+1,k) = sine*f_int_gxz + cosine*f_int_gyz
        gyz_new(i,j-1,k) = -sine*f_int_gxz + cosine*f_int_gyz
        gzz_new(i,j+1,k) = f_int_gzz
        gzz_new(i,j-1,k) = f_int_gzz
        !Now get interpolated values of the Aij
        if (i==imin) then
           ya1(1) = Axx_new(im,j,k)
           ya2(1) = -Axy_new(im,j,k)
           ya3(1) = -Axz_new(im,j,k)
           ya4(1) = Ayy_new(im,j,k)
           ya5(1) = Ayz_new(im,j,k)
           ya6(1) = Azz_new(im,j,k)
           do l=2,4
              ya1(l) = Axx_new(im+l-2,j,k)
              ya2(l) = Axy_new(im+l-2,j,k)
              ya3(l) = Axz_new(im+l-2,j,k)
              ya4(l) = Ayy_new(im+l-2,j,k)
              ya5(l) = Ayz_new(im+l-2,j,k)
              ya6(l) = Azz_new(im+l-2,j,k)
           enddo
        else
           do l=1,4
              ya1(l) = Axx_new(im+l-1,j,k)
              ya2(l) = Axy_new(im+l-1,j,k)
              ya3(l) = Axz_new(im+l-1,j,k)
              ya4(l) = Ayy_new(im+l-1,j,k)
              ya5(l) = Ayz_new(im+l-1,j,k)
              ya6(l) = Azz_new(im+l-1,j,k)
           enddo
        endif
        call polint(xa,ya1,4,rho,f_int_Axx,err)
        call polint(xa,ya2,4,rho,f_int_Axy,err)
        call polint(xa,ya3,4,rho,f_int_Axz,err)
        call polint(xa,ya4,4,rho,f_int_Ayy,err)
        call polint(xa,ya5,4,rho,f_int_Ayz,err)
        call polint(xa,ya6,4,rho,f_int_Azz,err)
        !Now rotate the Aij off of the plane
        Axx_new(i,j+1,k) = cosine*cosine*f_int_Axx - TWO*sine*cosine*f_int_Axy &
             & + sine*sine*f_int_Ayy
        Axx_new(i,j-1,k) = cosine*cosine*f_int_Axx + TWO*sine*cosine*f_int_Axy &
             & + sine*sine*f_int_Ayy
        Axy_new(i,j+1,k) = cosine*sine*f_int_Axx + (cosine*cosine-sine*sine)*f_int_Axy &
             & - sine*cosine*f_int_Ayy
        Axy_new(i,j-1,k) = -cosine*sine*f_int_Axx + (cosine*cosine-sine*sine)*f_int_Axy &
             & + sine*cosine*f_int_Ayy
        Axz_new(i,j+1,k) = cosine*f_int_Axz - sine*f_int_Ayz
        Axz_new(i,j-1,k) = cosine*f_int_Axz + sine*f_int_Ayz
        Ayy_new(i,j+1,k) = sine*sine*f_int_Axx + TWO*sine*cosine*f_int_Axy &
             & + cosine*cosine*f_int_Ayy
        Ayy_new(i,j-1,k) = sine*sine*f_int_Axx - TWO*sine*cosine*f_int_Axy &
             & + cosine*cosine*f_int_Ayy
        Ayz_new(i,j+1,k) = sine*f_int_Axz + cosine*f_int_Ayz
        Ayz_new(i,j-1,k) = -sine*f_int_Axz + cosine*f_int_Ayz
        Azz_new(i,j+1,k) = f_int_Azz
        Azz_new(i,j-1,k) = f_int_Azz
        !Now get interpolated values of the Gamma^i
        if (i==imin) then
           ya1(1) = -Gammax_new(im,j,k)
           ya2(1) = Gammay_new(im,j,k)
           ya3(1) = Gammaz_new(im,j,k)
           do l=2,4
              ya1(l) = Gammax_new(im+l-2,j,k)
              ya2(l) = Gammay_new(im+l-2,j,k)
              ya3(l) = Gammaz_new(im+l-2,j,k)
           enddo
        else
           do l=1,4
              ya1(l) = Gammax_new(im+l-1,j,k)
              ya2(l) = Gammay_new(im+l-1,j,k)
              ya3(l) = Gammaz_new(im+l-1,j,k)
           enddo
        endif
        call polint(xa,ya1,4,rho,f_int_Gamx,err)
        call polint(xa,ya2,4,rho,f_int_Gamy,err)
        call polint(xa,ya3,4,rho,f_int_Gamz,err)
        !Finally, rotate the Gamma^i off of the plane
        Gammax_new(i,j+1,k) = cosine*f_int_Gamx - sine*f_int_Gamy
        Gammax_new(i,j-1,k) = cosine*f_int_Gamx + sine*f_int_Gamy
        Gammay_new(i,j+1,k) = sine*f_int_Gamx + cosine*f_int_Gamy
        Gammay_new(i,j-1,k) = -sine*f_int_Gamx + cosine*f_int_Gamy
        Gammaz_new(i,j+1,k) = f_int_Gamz
        Gammaz_new(i,j-1,k) = f_int_Gamz
     enddo
  end do
end subroutine rotate
!
!--------------------------------------------------------------------
!---------
!
! Interpolation for point i,j,k
!
!--------------------------------------------------------------------
!---------
!
subroutine weights(X,Y,Z,alpha,phi,dX,dY,dZ,dT,i,j,k, imin,imax,jmin&
     &,jmax,kmin,kmax, il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry,&
     & Xcenter, Ycenter, Zcenter)
  implicit none
  real*8, intent(in)   :: X,Y,Z,dX,dY,dZ,dT,alpha,phi
  real*8, intent(in)   :: Xcenter, Ycenter, Zcenter
  integer, intent(in)  :: i,j,k
  integer, intent(in)  :: imin,imax,jmin,jmax,kmin,kmax
  integer, intent(out) :: il,iu,jl,ju,kl,ku
  real*8, intent(out)  :: Dlx,Dly,Dlz,fac
  integer, intent(in)  :: Symmetry
  !
  real*8               :: xo,yo,zo,xl,xu,yl,yu,zl,zu,r,r1,r2
  real*8               :: ONE,TWO
  integer              :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM,AXISYM
  parameter ( ONE = 1.D0, TWO = 2.D0 )
  parameter (NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  !
  ! find intersection of characteristic with previous timeslice
  ! 
  r  = sqrt(X**2 + Y**2 + Z**2)
  !  r1 = sqrt((X-Xcenter)**2 + (Y-Ycenter)**2 + (Z-Zcenter)**2)
  !  r2 = sqrt((X+Xcenter)**2 + (Y+Ycenter)**2 + (Z+Zcenter)**2)
  fac = ONE - (alpha + ONE)*exp(-TWO*phi)*dT/r
  xo = X*fac
  yo = Y*fac
  zo = Z*fac
  !
  ! find indices for interpolation
  !        
  if (xo < X) then
     il = i - 1
     xl = X - dX
  else
     il = i
     xl = X
  end if
  if (yo < Y) then
     jl = j - 1
     yl = Y - dY
  else
     jl = j
     yl = Y
  end if
  if (zo < Z) then
     kl = k - 1
     zl = Z - dZ
  else
     kl = k
     zl = Z
  end if
  !
  ! make sure that interpolation molecule does not extent beyond grid:
  !
  if (il > imax - 1) then
     il = imax - 1
     xl = X - dX
  end if
  if (jl > jmax - 1) then
     jl = jmax - 1
     yl = Y - dY
  end if
  if (kl > kmax - 1) then
     kl = kmax - 1
     zl = Z - dZ
  end if
  if (Symmetry == EQUATORIAL) then
     if (il < imin) then
        il = imin
        xl = X
     end if
     if (jl < jmin) then
        jl = jmin
        yl = Y
     end if
     iu = il + 1
     ju = jl + 1
     ku = kl + 1
     xu = xl + dX
     yu = yl + dY
     zu = zl + dZ
     if (kl < kmin) then
        kl = kmin
        !       zl = Z
     end if
  endif
  if (Symmetry == PI_SYMM) then
     if (il < imin) then
        il = imin
        xl = X
     end if
     iu = il + 1
     ju = jl + 1
     ku = kl + 1
     xu = xl + dX
     yu = yl + dY
     zu = zl + dZ
     if (kl < kmin) then
        kl = kmin
        !        zl = Z
     end if
  endif
  ! a change: use the following case for both octant
  ! and axisymmetry
  if (Symmetry == OCTANT .or. Symmetry == AXISYM) then
     iu = il + 1
     ju = jl + 1
     ku = kl + 1
     xu = xl + dX
     yu = yl + dY
     zu = zl + dZ
     if (il < imin) then
        il = imin
        !        xl = X
     end if
     if (jl < jmin) then
        jl = jmin
        !        yl = Y 
     end if
     if (kl < kmin) then
        kl = kmin
        !        zl = Z
     end if
  endif
  if (Symmetry == NO_SYMM) then
     if (il < imin) then
        il = imin
        xl = X
     end if
     if (jl < jmin) then
        jl = jmin
        yl = Y
     end if
     if (kl < kmin) then
        kl = kmin
        zl = Z
     end if
     iu = il + 1
     ju = jl + 1
     ku = kl + 1
     xu = xl + dX
     yu = yl + dY
     zu = zl + dZ
  endif
  !
  !
  !!
  ! get upper values
  !
  ! find coefficients for interpolation
  !
  Dlx = (xo - xl)/(xu - xl)
  Dly = (yo - yl)/(yu - yl)
  Dlz = (zo - zl)/(zu - zl)
  return
end subroutine weights
!
!--------------------------------------------------------------------
!---------
!
! Same as above, but use \sqrt{g_{rr}} instead of e^{-2\phi}
!
!--------------------------------------------------------------------
!---------
!
subroutine weights2(X,Y,Z,Rp,dR,alpha,phi,dX,dY,dZ,dT,i,j,k, imin,imax,jmin&
     &,jmax,kmin,kmax, il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,facr,facRp,Symmetry,&
     & Xcenter, Ycenter, Zcenter, gxx,gxy,gxz,gyy,gyz,gzz)
  implicit none
  real*8, intent(in)   :: X,Y,Z,dX,dY,dZ,Rp,dT,alpha,phi,dR
  real*8, intent(in)   :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8, intent(in)   :: Xcenter, Ycenter, Zcenter
  integer, intent(in)  :: i,j,k
  integer, intent(in)  :: imin,imax,jmin,jmax,kmin,kmax
  integer, intent(out) :: il,iu,jl,ju,kl,ku
  real*8, intent(out)  :: Dlx,Dly,Dlz,facr, facRp
  integer, intent(in)  :: Symmetry
  !
  real*8               :: xo,yo,zo,xl,xu,yl,yu,zl,zu,r,r1,r2,sqrt_tilgrr,tmp
  real*8               :: ONE,TWO
  integer              :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM,AXISYM
  parameter ( ONE = 1.D0, TWO = 2.D0 )
  parameter (NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  !
  ! find intersection of characteristic with previous timeslice
  ! 
  r  = sqrt(X**2 + Y**2 + Z**2)
  sqrt_tilgrr = (gxx*X*X + gyy*Y*Y + gzz*Z*Z + TWO*(gxy*X*Y + gxz*X*Z + gyz*Y*Z))
  sqrt_tilgrr = sqrt(sqrt_tilgrr)/r
  !  r1 = sqrt((X-Xcenter)**2 + (Y-Ycenter)**2 + (Z-Zcenter)**2)
  !  r2 = sqrt((X+Xcenter)**2 + (Y+Ycenter)**2 + (Z+Zcenter)**2)
!!$  fac = ONE - (alpha + ONE)*exp(-TWO*phi)*dT/r
  facr = ONE - (alpha + ONE)*exp(-TWO*phi)*dT/(r*sqrt_tilgrr)
  xo = X*facr
  yo = Y*facr
  zo = Z*facr
  facRp = ONE - (alpha + ONE)*exp(-TWO*phi)*dT*dR/Rp/sqrt_tilgrr
  !
  ! find indices for interpolation
  !        
  if (xo < X) then
     il = i - 1
     xl = X - dX
  else
     il = i
     xl = X
  end if
  if (yo < Y) then
     jl = j - 1
     yl = Y - dY
  else
     jl = j
     yl = Y
  end if
  if (zo < Z) then
     kl = k - 1
     zl = Z - dZ
  else
     kl = k
     zl = Z
  end if
  !
  ! make sure that interpolation molecule does not extent beyond grid:
  !
  if (il > imax - 1) then
     il = imax - 1
     xl = X - dX
  end if
  if (jl > jmax - 1) then
     jl = jmax - 1
     yl = Y - dY
  end if
  if (kl > kmax - 1) then
     kl = kmax - 1
     zl = Z - dZ
  end if
  if (Symmetry == EQUATORIAL) then
     if (il < imin) then
        il = imin
        xl = X
     end if
     if (jl < jmin) then
        jl = jmin
        yl = Y
     end if
     iu = il + 1
     ju = jl + 1
     ku = kl + 1
     xu = xl + dX
     yu = yl + dY
     zu = zl + dZ
     if (kl < kmin) then
        kl = kmin
        !       zl = Z
     end if
  endif
  if (Symmetry == PI_SYMM) then
     if (il < imin) then
        il = imin
        xl = X
     end if
     iu = il + 1
     ju = jl + 1
     ku = kl + 1
     xu = xl + dX
     yu = yl + dY
     zu = zl + dZ
     if (kl < kmin) then
        kl = kmin
        !        zl = Z
     end if
  endif
  ! a change: use the following case for both octant
  ! and axisymmetry
  if (Symmetry == OCTANT .or. Symmetry == AXISYM) then
     iu = il + 1
     ju = jl + 1
     ku = kl + 1
     xu = xl + dX
     yu = yl + dY
     zu = zl + dZ
     if (il < imin) then
        il = imin
        !        xl = X
     end if
     if (jl < jmin) then
        jl = jmin
        !        yl = Y 
     end if
     if (kl < kmin) then
        kl = kmin
        !        zl = Z
     end if
  endif
  if (Symmetry == NO_SYMM) then
     if (il < imin) then
        il = imin
        xl = X
     end if
     if (jl < jmin) then
        jl = jmin
        yl = Y
     end if
     if (kl < kmin) then
        kl = kmin
        zl = Z
     end if
     iu = il + 1
     ju = jl + 1
     ku = kl + 1
     xu = xl + dX
     yu = yl + dY
     zu = zl + dZ
  endif
  !
  !
  !!
  ! get upper values
  !
  ! find coefficients for interpolation
  !
!  write(*,*) HI,xo,xl,xu,(xo - xl)/(xu - xl)
!!$  if(xo<0.D0) then
!!$     tmp = xl
!!$     xl = xu
!!$     xu = tmp
!!$  end if
!!$  if(yo<0.D0) then
!!$     tmp = yl
!!$     yl = yu
!!$     yu = tmp
!!$  end if
!!$  if(zo<0.D0) then
!!$     tmp = zl
!!$     zl = zu
!!$     zu = tmp
!!$  end if
  Dlx = (xo - xl)/(xu - xl)
  Dly = (yo - yl)/(yu - yl)
  Dlz = (zo - zl)/(zu - zl)
  return
end subroutine weights2
!
!--------------------------------------------------------------------
!---------
!
! Interpolate function f using weights Dlx, Dly and Dlz
!
!--------------------------------------------------------------------
!---------
!
subroutine interp(ex,f,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,xsym,ysym&
     &,zsym)
  implicit none
  !
  ! Input:
  !
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: f
  real*8                                   :: f_int
  integer                                  :: il,iu,jl,ju,kl,ku
  real*8                                   :: Dlx,Dly,Dlz
  real*8                                   :: xsym,ysym,zsym
  !
  ! Other variables:
  !
  real*8      :: xs, ys, zs
  real*8      :: l
  parameter ( l = 1.D0 )
  !
  ! set symmetries if necessary
  !
  xs = l
  ys = l
  zs = l
  if (il == iu) xs = xsym
  if (jl == ju) ys = ysym
  if (kl == ku) zs = zsym
  !
  f_int = xs*(l-Dlx)*(ys*(l-Dly)*(zs*(l-Dlz)*f(il,jl,kl) + Dlz*f(il&
       &,jl,ku)) +          Dly *(zs*(l-Dlz)*f(il,ju,kl) + Dlz*f(il&
       &,ju,ku)) ) +    Dlx*( ys*(l-Dly)*(zs*(l-Dlz)*f(iu,jl,kl) +&
       & Dlz*f(iu,jl,ku)) +          Dly *(zs*(l-Dlz)*f(iu,ju,kl) +&
       & Dlz*f(iu,ju,ku)) )
  return
end subroutine interp
!---------------------------------------------------------------------
! Boundary condition on the extrinsic curvature.
! Only to be used in the binary case.
!---------------------------------------------------------------------
!
subroutine update_k(ex, X, Y, Z, PhysicalR, RDeriv,dT, &
     lapse_old, phi, kxx_old,kxy_old&
     &,kxz_old,kyy_old,kyz_old,kzz_old, kxx_new,kxy_new,kxz_new&
     &,kyy_new,kyz_new,kzz_new, Axx_old,Axy_old,Axz_old,Ayy_old&
     &,Ayz_old,Azz_old, Axx_new,Axy_new,Axz_new,Ayy_new,Ayz_new&
     &,Azz_new, Symmetry)
  implicit none
  interface
     subroutine gderivs_eq(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8                                   :: dX,dY,dZ,SYM1,SYM2&
            &,SYM3
     end subroutine gderivs_eq
  end interface
  !
  ! Input parameters:
  !
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z
  real*8                                   :: dT
  real*8, dimension(ex(1),ex(2),ex(3))     :: PhysicalR,RDeriv
  real*8, dimension(ex(1),ex(2),ex(3))     :: lapse_old,phi
  real*8, dimension(ex(1),ex(2),ex(3))     :: kxx_old,kxy_old,kxz_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: kyy_old,kyz_old,kzz_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: kxx_new,kxy_new,kxz_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: kyy_new,kyz_new,kzz_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: Axx_old,Axy_old,Axz_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: Ayy_old,Ayz_old,Azz_old
  real*8, dimension(ex(1),ex(2),ex(3))     :: Axx_new,Axy_new,Axz_new
  real*8, dimension(ex(1),ex(2),ex(3))     :: Ayy_new,Ayz_new,Azz_new
  integer                                  :: Symmetry
  !
  ! Other variables:
  !
  integer                    :: i, j, k
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  integer                    :: il,jl,kl,iu,ju,ku
  integer                    :: NO_SYMM, EQUATORIAL, OCTANT
  real*8                     :: Courant,dX,dY,dZ,Dlx,Dly,Dlz,f_int,fac
  real*8                      :: HALF, ONE, TWO, ZERO, F2o3
  real*8                     :: SYM, ANT, Xcenter, Ycenter, Zcenter
  real*8, dimension(ex(1),ex(2),ex(3))  :: f_rad
  parameter(HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, F2o3 =&
       & 2.D0/3.D0 )
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2)
  parameter ( SYM = 1.D0, ANT = - 1.D0 )
  !
  ! where the action is...
  !
  imin = lbound(Axx_old,1)
  jmin = lbound(Axx_old,2)
  kmin = lbound(Axx_old,3)
  imax = ubound(Axx_old,1)
  jmax = ubound(Axx_old,2)
  kmax = ubound(Axx_old,3)
  dX = X(imin+1,1,1) - X(imin,1,1)
  dY = Y(1,jmin+1,1) - Y(1,jmin,1)
  dZ = Z(1,1,kmin+1) - Z(1,1,kmin)
  Xcenter = ZERO
  Ycenter = ZERO
  Zcenter = ZERO
  !
  !--------------------------------------------------------------------
  !---------
  ! Axx
  !--------------------------------------------------------------------
  !---------
  f_rad = Axx_old - kxx_old
  i = imax
  do j = jmin, jmax - 1
     do k = kmin, kmax - 1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,SYM&
             &,SYM,SYM)
        Axx_new(i,j,k) = f_int * fac + kxx_new(i,j,k)
     end do
  end do
  j = jmax
  do i = imin, imax
     do k = kmin, kmax - 1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,SYM&
             &,SYM,SYM)
        Axx_new(i,j,k) = f_int * fac + kxx_new(i,j,k)
     end do
  end do
  k = kmax
  do i = imin, imax
     do j = jmin, jmax
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,SYM&
             &,SYM,SYM)
        Axx_new(i,j,k) = f_int * fac + kxx_new(i,j,k)
     end do
  end do
  i = imin
  do j = jmin, jmax-1
     do k = kmin, kmax-1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,SYM&
             &,SYM,SYM)
        Axx_new(i,j,k) = f_int * fac + kxx_new(i,j,k)
     end do
  end do
  j = jmin
  do i = imin, imax-1
     do k = kmin, kmax-1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,SYM&
             &,SYM,SYM)
        Axx_new(i,j,k) = f_int * fac + kxx_new(i,j,k)
     end do
  end do
  !
  !--------------------------------------------------------------------
  !---------
  ! Axy
  !--------------------------------------------------------------------
  !---------
  f_rad = Axy_old - kxy_old
  i = imax
  do j = jmin, jmax - 1
     do k = kmin, kmax - 1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,ANT&
             &,ANT,SYM)
        Axy_new(i,j,k) = f_int * fac + kxy_new(i,j,k)
     end do
  end do
  j = jmax
  do i = imin, imax
     do k = kmin, kmax - 1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,ANT&
             &,ANT,SYM)
        Axy_new(i,j,k) = f_int * fac + kxy_new(i,j,k)
     end do
  end do
  k = kmax
  do i = imin, imax
     do j = jmin, jmax
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,ANT&
             &,ANT,SYM)
        Axy_new(i,j,k) = f_int * fac + kxy_new(i,j,k)
     end do
  end do
  i = imin
  do j = jmin, jmax-1
     do k = kmin, kmax-1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,ANT&
             &,ANT,SYM)
        Axy_new(i,j,k) = f_int * fac + kxy_new(i,j,k)
     end do
  end do
  j = jmin
  do i = imin, imax-1
     do k = kmin, kmax-1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,ANT&
             &,ANT,SYM)
        Axy_new(i,j,k) = f_int * fac + kxy_new(i,j,k)
     end do
  end do
  !
  !--------------------------------------------------------------------
  !---------
  ! Axz
  !--------------------------------------------------------------------
  !---------
  f_rad = Axz_old - kxz_old
  i = imax
  do j = jmin, jmax - 1
     do k = kmin, kmax - 1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,ANT&
             &,SYM,ANT)
        Axz_new(i,j,k) = f_int * fac + kxz_new(i,j,k)
     end do
  end do
  j = jmax
  do i = imin, imax
     do k = kmin, kmax - 1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,ANT&
             &,SYM,ANT)
        Axz_new(i,j,k) = f_int * fac + kxz_new(i,j,k)
     end do
  end do
  k = kmax
  do i = imin, imax
     do j = jmin, jmax
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,ANT&
             &,SYM,ANT)
        Axz_new(i,j,k) = f_int * fac + kxz_new(i,j,k)
     end do
  end do
  i = imin
  do j = jmin, jmax-1
     do k = kmin, kmax-1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,ANT&
             &,SYM,ANT)
        Axz_new(i,j,k) = f_int * fac + kxz_new(i,j,k)
     end do
  end do
  j = jmin
  do i = imin, imax-1
     do k = kmin, kmax-1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,ANT&
             &,SYM,ANT)
        Axz_new(i,j,k) = f_int * fac + kxz_new(i,j,k)
     end do
  end do
  !
  !--------------------------------------------------------------------
  !---------
  ! Ayy
  !--------------------------------------------------------------------
  !---------
  f_rad = Ayy_old - kyy_old
  i = imax
  do j = jmin, jmax - 1
     do k = kmin, kmax - 1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,SYM&
             &,SYM,SYM)
        Ayy_new(i,j,k) = f_int * fac + kyy_new(i,j,k)
     end do
  end do
  j = jmax
  do i = imin, imax
     do k = kmin, kmax - 1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,SYM&
             &,SYM,SYM)
        Ayy_new(i,j,k) = f_int * fac + kyy_new(i,j,k)
     end do
  end do
  k = kmax
  do i = imin, imax
     do j = jmin, jmax
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,SYM&
             &,SYM,SYM)
        Ayy_new(i,j,k) = f_int * fac + kyy_new(i,j,k)
     end do
  end do
  i = imin
  do j = jmin, jmax-1
     do k = kmin, kmax-1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,SYM&
             &,SYM,SYM)
        Ayy_new(i,j,k) = f_int * fac + kyy_new(i,j,k)
     end do
  end do
  j = jmin
  do i = imin, imax-1
     do k = kmin, kmax-1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,SYM&
             &,SYM,SYM)
        Ayy_new(i,j,k) = f_int * fac + kyy_new(i,j,k)
     end do
  end do
  !
  !--------------------------------------------------------------------
  !---------
  ! Ayz
  !--------------------------------------------------------------------
  !---------
  f_rad = Ayz_old - kyz_old
  i = imax
  do j = jmin, jmax - 1
     do k = kmin, kmax - 1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,SYM&
             &,ANT,ANT)
        Ayz_new(i,j,k) = f_int * fac + kyz_new(i,j,k)
     end do
  end do
  j = jmax
  do i = imin, imax
     do k = kmin, kmax - 1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,SYM&
             &,ANT,ANT)
        Ayz_new(i,j,k) = f_int * fac + kyz_new(i,j,k)
     end do
  end do
  k = kmax
  do i = imin, imax
     do j = jmin, jmax
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,SYM&
             &,ANT,ANT)
        Ayz_new(i,j,k) = f_int * fac + kyz_new(i,j,k)
     end do
  end do
  i = imin
  do j = jmin, jmax-1
     do k = kmin, kmax-1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,SYM&
             &,ANT,ANT)
        Ayz_new(i,j,k) = f_int * fac + kyz_new(i,j,k)
     end do
  end do
  j = jmin
  do i = imin, imax-1
     do k = kmin, kmax-1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,SYM&
             &,ANT,ANT)
        Ayz_new(i,j,k) = f_int * fac + kyz_new(i,j,k)
     end do
  end do
  !
  !--------------------------------------------------------------------
  !---------
  ! Azz
  !--------------------------------------------------------------------
  !---------
  f_rad = Azz_old - kzz_old
  i = imax
  do j = jmin, jmax - 1
     do k = kmin, kmax - 1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,SYM&
             &,SYM,SYM)
        Azz_new(i,j,k) = f_int * fac + kzz_new(i,j,k)
     end do
  end do
  j = jmax
  do i = imin, imax
     do k = kmin, kmax - 1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,SYM&
             &,SYM,SYM)
        Azz_new(i,j,k) = f_int * fac + kzz_new(i,j,k)
     end do
  end do
  k = kmax
  do i = imin, imax
     do j = jmin, jmax
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,SYM&
             &,SYM,SYM)
        Azz_new(i,j,k) = f_int * fac + kzz_new(i,j,k)
     end do
  end do
  i = imin
  do j = jmin, jmax-1
     do k = kmin, kmax-1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,SYM&
             &,SYM,SYM)
        Azz_new(i,j,k) = f_int * fac + kzz_new(i,j,k)
     end do
  end do
  j = jmin
  do i = imin, imax-1
     do k = kmin, kmax-1
        call weights(X(i,j,k),Y(i,j,k),Z(i,j,k),PhysicalR(i,j,k),lapse_old(i,j,k),phi(i,j,k), dX&
             &,dY,dZ,dT,i,j,k,imin,imax,jmin,jmax, kmin,kmax,il,iu,jl&
             &,ju,kl,ku,Dlx,Dly,Dlz,fac,Symmetry, Xcenter,Ycenter&
             &,Zcenter)
        call interp(ex,f_rad,f_int,il,iu,jl,ju,kl,ku,Dlx,Dly,Dlz,SYM&
             &,SYM,SYM)
        Azz_new(i,j,k) = f_int * fac + kzz_new(i,j,k)
     end do
  end do
end subroutine update_k
!
!-----------------------------------------------------------------------------
!
! c1/r + c2/r^3 boundary condition for \phi
!
!-----------------------------------------------------------------------------
!
subroutine phibc(ex,x,y,z,f,Symmetry)
  implicit none
  !
  ! Input parameters:
  !
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: f
  real*8                                   :: monopole,quadrupole
  integer                                  :: Symmetry
  !
  ! Other variables:
  !
  integer                    :: i, j, k
  integer                    :: imin, jmin, kmin, imax, jmax, kmax
  integer                    :: il,jl,kl,iu,ju,ku,bc
  integer                    :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  real*8                     :: Courant,dX,dY,dZ
  real*8                     :: xu,yu,zu,ru,rl,xl,yl,zl
  real*8                     :: f1,f2,r,P,r1,r2,P1,P2,d,c1,c2,xa,ya,za
  real*8                     :: HALF, ONE, TWO, THREE
  real*8, dimension(ex(1),ex(2),ex(3))     :: fold
  parameter(HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0, THREE = 3.D0)
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  !
  ! where the action is...
  !
  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)
  imax = ubound(f,1)
  jmax = ubound(f,2)
  kmax = ubound(f,3)
  dX = X(imin+1,1,1) - X(imin,1,1)
  dY = Y(1,jmin+1,1) - Y(1,jmin,1)
  dZ = Z(1,1,kmin+1) - Z(1,1,kmin)
  !
  ! x = max
  !
  i = imax
  do j = jmin, jmax
     do k = kmin,  kmax
        f1 = f(i-1,j,k)
        f2 = f(i-2,j,k)
        r = sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k) + z(i,j,k)*z(i,j,k))
        r1 = sqrt(x(i-1,j,k)*x(i-1,j,k) + y(i,j,k)*y(i,j,k) + z(i,j,k)*z(i,j,k))
        r2 = sqrt(x(i-2,j,k)*x(i-2,j,k) + y(i,j,k)*y(i,j,k) + z(i,j,k)*z(i,j,k))
        d = r2*r2 - r1*r1
        c1 = (f2*r2*r2*r2 - f1*r1*r1*r1)/d
        c2 = r1*r1*r2*r2*(f1*r1 - f2*r2)/d
        f(i,j,k) = c1/r + c2/(r*r*r)
     enddo
  enddo
  !
  ! y = max
  !
  if(Symmetry .ne. AXISYM) then
     j = jmax
     do i = imin, imax
        do k = kmin,  kmax
           f1 = f(i,j-1,k)
           f2 = f(i,j-2,k)
           r = sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k) + z(i,j,k)*z(i,j,k))
           r1 = sqrt(x(i,j,k)*x(i,j,k) + y(i,j-1,k)*y(i,j-1,k) + z(i,j,k)*z(i,j,k))
           r2 = sqrt(x(i,j,k)*x(i,j,k) + y(i,j-2,k)*y(i,j-2,k) + z(i,j,k)*z(i,j,k))
           d = r2*r2 - r1*r1
           c1 = (f2*r2*r2*r2 - f1*r1*r1*r1)/d
           c2 = r1*r1*r2*r2*(f1*r1 - f2*r2)/d
           f(i,j,k) = c1/r + c2/(r*r*r)
        enddo
     enddo
  end if
  !
  ! z = max
  !
  k = kmax
  do i = imin, imax - 1
     do j = jmin, jmax - 1
        f1 = f(i,j,k-1)
        f2 = f(i,j,k-2)
        r = sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k) + z(i,j,k)*z(i,j,k))
        r1 = sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k) + z(i,j,k-1)*z(i,j,k-1))
        r2 = sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k) + z(i,j,k-2)*z(i,j,k-2))
        d = r2*r2 - r1*r1
        c1 = (f2*r2*r2*r2 - f1*r1*r1*r1)/d
        c2 = r1*r1*r2*r2*(f1*r1 - f2*r2)/d
        f(i,j,k) = c1/r + c2/(r*r*r)
     enddo
  enddo
  if (Symmetry == EQUATORIAL.or.Symmetry == NO_SYMM.or.Symmetry == PI_SYMM) then
     !
     ! x = min
     !
     i = imin
     do j = jmin, jmax
        do k = kmin, kmax
           f1 = f(i+1,j,k)
           f2 = f(i+2,j,k)
           r = sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k) + z(i,j,k)*z(i,j,k))
           r1 = sqrt(x(i+1,j,k)*x(i+1,j,k) + y(i,j,k)*y(i,j,k) + z(i,j,k)*z(i,j,k))
           r2 = sqrt(x(i+2,j,k)*x(i+2,j,k) + y(i,j,k)*y(i,j,k) + z(i,j,k)*z(i,j,k))
           d = r2*r2 - r1*r1
           f(i,j,k) = f1 * (r1/r)
           c1 = (f2*r2*r2*r2 - f1*r1*r1*r1)/d
           c2 = r1*r1*r2*r2*(f1*r1 - f2*r2)/d
           f(i,j,k) = c1/r + c2/(r*r*r)
        enddo
     enddo
  endif
  if (Symmetry == EQUATORIAL.or.Symmetry == NO_SYMM) then
     !
     ! y = min
     !
     j = jmin
     do i = imin, imax
        do k = kmin, kmax
           f1 = f(i,j+1,k)
           f2 = f(i,j+2,k)
           r = sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k) + z(i,j,k)*z(i,j,k))
           r1 = sqrt(x(i,j,k)*x(i,j,k) + y(i,j+1,k)*y(i,j+1,k) + z(i,j,k)*z(i,j,k))
           r2 = sqrt(x(i,j,k)*x(i,j,k) + y(i,j+2,k)*y(i,j+2,k) + z(i,j,k)*z(i,j,k))
           d = r2*r2 - r1*r1
           c1 = (f2*r2*r2*r2 - f1*r1*r1*r1)/d
           c2 = r1*r1*r2*r2*(f1*r1 - f2*r2)/d
           f(i,j,k) = c1/r + c2/(r*r*r)
        enddo
     enddo
  endif
  if (Symmetry == NO_SYMM) then
     !
     ! z = min
     !
     k = kmin
     do i = imin + 1, imax
        do j = jmin + 1, jmax
           f1 = f(i,j,k+1)
           f2 = f(i,j,k+2)
           r = sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k) + z(i,j,k)*z(i,j,k))
           r1 = sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k) + z(i,j,k+1)*z(i,j,k+1))
           r2 = sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k) + z(i,j,k+2)*z(i,j,k+2))
           d = r2*r2 - r1*r1
           c1 = (f2*r2*r2*r2 - f1*r1*r1*r1)/d
           c2 = r1*r1*r2*r2*(f1*r1 - f2*r2)/d
           f(i,j,k) = c1/r + c2/(r*r*r)
        enddo
     enddo
  end if
end subroutine phibc
!
! Given a symmetric tensor tij = exp(4 phi) ttij, this subroutine transform 
! the components of ttij between fisheye coordinates and physical coordinates
! near the boundary of the grid
!
subroutine trans_fish_phys(ex,fake_ex, X, Y, Z, PhysicalR, RadiusDerivative, &
     ttxx,ttxy,ttxz,ttyy,ttyz,ttzz,fish_to_phys,lb, &
     have_bdry_min,have_bdry_max,use_trans_fish_phys_new)
  implicit none
  integer, dimension(3)                    :: ex,fake_ex,have_bdry_min,have_bdry_max
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: PhysicalR,RadiusDerivative
  real*8, dimension(ex(1),ex(2),ex(3))     :: ttxx,ttxy,ttxz,ttyy,ttyz,ttzz
  integer                                   :: lb,use_trans_fish_phys_new
  integer :: global_edge_x_max,global_edge_y_max,global_edge_z_max
  integer :: global_edge_x_min,global_edge_y_min,global_edge_z_min
  !
  integer                                   :: i,j,k,imax,jmax,kmax
  integer                                  :: imin,jmin,kmin
  logical                                    :: fish_to_phys
  logical, dimension(ex(1),ex(2),ex(3))    :: tree
  !
  imin = ex(1)-fake_ex(1)+1
  jmin = ex(2)-fake_ex(2)+1
  kmin = ex(3)-fake_ex(3)+1
  imax = fake_ex(1)
  jmax = fake_ex(2)
  kmax = fake_ex(3)
  ! We dont want to do the transformation more than once at each grid point.
  ! Thats why its necessary to mark the grid point where the
  ! transformation has been carried out. Tree is the logical variable
  ! to make the mark.
  !
  tree = .FALSE.
  ! xmax
  if(have_bdry_max(1)==1) then
     do k=kmin,kmax
        do j=jmin,jmax
           do i=imax-lb,imax
              ! Recognize the tree?
              if (tree(i,j,k) .eqv. .FALSE.) then
                 if(use_trans_fish_phys_new==0) then
                    call trans_fish_phys_ttij(ex,X,Y,Z,PhysicalR,RadiusDerivative, &
                         ttxx,ttxy,ttxz,ttyy,ttyz,ttzz,i,j,k,fish_to_phys)
                 else
                    call trans_fish_phys_ttij_new(ex,X,Y,Z,PhysicalR,RadiusDerivative, &
                         ttxx,ttxy,ttxz,ttyy,ttyz,ttzz,i,j,k,fish_to_phys)
                 end if
                 tree(i,j,k) = .TRUE.
              end if
           end do
        end do
     end do
  endif
  ! ymax
  if(have_bdry_max(2)==1) then
     do k=kmin,kmax
        do j=jmax-lb,jmax
           do i=imin,imax
              ! Recognize the tree?
              if (tree(i,j,k) .eqv. .FALSE.) then
                 if(use_trans_fish_phys_new==0) then
                    call trans_fish_phys_ttij(ex,X,Y,Z,PhysicalR,RadiusDerivative, &
                         ttxx,ttxy,ttxz,ttyy,ttyz,ttzz,i,j,k,fish_to_phys)
                 else
                    call trans_fish_phys_ttij_new(ex,X,Y,Z,PhysicalR,RadiusDerivative, &
                         ttxx,ttxy,ttxz,ttyy,ttyz,ttzz,i,j,k,fish_to_phys)
                 end if
                 tree(i,j,k) = .TRUE.
              end if
           end do
               end do
     end do
  end if
  ! zmax
  if(have_bdry_max(3)==1) then
     do k=kmax-lb,kmax
        do j=jmin,jmax
           do i=imin,imax
              ! Recognize the tree?
              if (tree(i,j,k) .eqv. .FALSE.) then
                 if(use_trans_fish_phys_new==0) then
                    call trans_fish_phys_ttij(ex,X,Y,Z,PhysicalR,RadiusDerivative, &
                         ttxx,ttxy,ttxz,ttyy,ttyz,ttzz,i,j,k,fish_to_phys)
                 else
                    call trans_fish_phys_ttij_new(ex,X,Y,Z,PhysicalR,RadiusDerivative, &
                         ttxx,ttxy,ttxz,ttyy,ttyz,ttzz,i,j,k,fish_to_phys)
                 end if
                 tree(i,j,k) = .TRUE.
              end if
           end do
        end do
     end do
  endif
  ! xmin
  if(have_bdry_min(1)==1) then
     do k=kmin,kmax
        do j=jmin,jmax
           do i=imin,imin+lb
              ! Recognize the tree?
              if (tree(i,j,k) .eqv. .FALSE.) then
                 if(use_trans_fish_phys_new==0) then
                    call trans_fish_phys_ttij(ex,X,Y,Z,PhysicalR,RadiusDerivative, &
                         ttxx,ttxy,ttxz,ttyy,ttyz,ttzz,i,j,k,fish_to_phys)
                 else
                    call trans_fish_phys_ttij_new(ex,X,Y,Z,PhysicalR,RadiusDerivative, &
                         ttxx,ttxy,ttxz,ttyy,ttyz,ttzz,i,j,k,fish_to_phys)
                 end if
                 tree(i,j,k) = .TRUE.
              end if
           end do
        end do
     end do
  end if
  ! ymin 
  if(have_bdry_min(2)==1) then
     do k=kmin,kmax
        do j=jmin,jmin+lb
           do i=imin,imax
              ! Recognize the tree?
              if (tree(i,j,k) .eqv. .FALSE.) then
                 if(use_trans_fish_phys_new==0) then
                    call trans_fish_phys_ttij(ex,X,Y,Z,PhysicalR,RadiusDerivative, &
                         ttxx,ttxy,ttxz,ttyy,ttyz,ttzz,i,j,k,fish_to_phys)
                 else
                    call trans_fish_phys_ttij_new(ex,X,Y,Z,PhysicalR,RadiusDerivative, &
                         ttxx,ttxy,ttxz,ttyy,ttyz,ttzz,i,j,k,fish_to_phys)
                 end if
                 tree(i,j,k) = .TRUE.
              end if
           end do
        end do
     end do
  end if
  ! zmin
  if(have_bdry_min(3)==1) then
     do k=kmin,kmin+lb
        do j=jmin,jmax
           do i=imin,imax
              ! Recognize the tree?
              if (tree(i,j,k) .eqv. .FALSE.) then
                 if(use_trans_fish_phys_new==0) then
                    call trans_fish_phys_ttij(ex,X,Y,Z,PhysicalR,RadiusDerivative, &
                         ttxx,ttxy,ttxz,ttyy,ttyz,ttzz,i,j,k,fish_to_phys)
                 else
                    call trans_fish_phys_ttij_new(ex,X,Y,Z,PhysicalR,RadiusDerivative, &
                         ttxx,ttxy,ttxz,ttyy,ttyz,ttzz,i,j,k,fish_to_phys)
                 end if
                 tree(i,j,k) = .TRUE.
              end if
           end do
        end do
     end do
  end if
end subroutine trans_fish_phys
subroutine trans_fish_phys_ttij(ex, X, Y, Z, PhysicalR, RadiusDerivative, &
     ttxx,ttxy,ttxz,ttyy,ttyz,ttzz,i,j,k,fish_to_phys)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: PhysicalR,RadiusDerivative
  real*8, dimension(ex(1),ex(2),ex(3))     :: ttxx,ttxy,ttxz,ttyy,ttyz,ttzz
  real*8                                   :: Rp,dr,r,xi,yi,zi,fac
  real*8                                   :: Jxx,Jxy,Jxz,Jyy,Jyz,Jzz
  real*8                                   :: ttxxi,ttxyi,ttxzi,ttyyi,ttyzi,ttzzi
  integer                                   :: i,j,k
  logical                                   :: fish_to_phys
  !
  xi = X(i,1,1)
  yi = Y(1,j,1)
  zi = Z(1,1,k)
  ttxxi = ttxx(i,j,k)
  ttxyi = ttxy(i,j,k)
  ttxzi = ttxz(i,j,k)
  ttyyi = ttyy(i,j,k)
  ttyzi = ttyz(i,j,k)
  ttzzi = ttzz(i,j,k)
  r = sqrt(xi**2 + yi**2 + zi**2)
  Rp = PhysicalR(i,j,k)
  dR = RadiusDerivative(i,j,k)
  ! transformation matrix
  if (fish_to_phys) then
     fac = 1.d0/dR - r/Rp
     Jxx = r/Rp + xi/r * xi/r * fac
     Jxy = xi/r * yi/r * fac
     Jxz = xi/r * zi/r * fac
     Jyy = r/Rp + yi/r * yi/r * fac
     Jyz = yi/r * zi/r * fac
     Jzz = r/Rp + zi/r * zi/r * fac
     fac = (Rp/r)**(4.d0/3.d0) * dR**(2.d0/3.d0)
  else
     fac = dR - Rp/r
     Jxx = Rp/r + xi/r * xi/r * fac
     Jxy = xi/r * yi/r * fac
     Jxz = xi/r * zi/r * fac
     Jyy = Rp/r + yi/r * yi/r * fac
     Jyz = yi/r * zi/r * fac
     Jzz = Rp/r + zi/r * zi/r * fac
     fac = (Rp/r)**(-4.d0/3.d0) * dR**(-2.d0/3.d0)
  end if
  ttxx(i,j,k) = (Jxx*Jxx*ttxxi + 2.d0*Jxx*Jxy*ttxyi + 2.d0*Jxx*Jxz*ttxzi + &
       Jxy*Jxy*ttyyi + 2.d0*Jxy*Jxz*ttyzi + Jxz*Jxz*ttzzi)*fac
  ttxy(i,j,k) = (Jxx*Jxy*ttxxi + (Jxx*Jyy+Jxy*Jxy)*ttxyi + &
       (Jxx*Jyz+Jxz*Jxy)*ttxzi + Jxy*Jyy*ttyyi + &
       (Jxy*Jyz+Jxz*Jyy)*ttyzi + Jxz*Jyz*ttzzi)*fac
  ttxz(i,j,k) = (Jxx*Jxz*ttxxi + (Jxx*Jyz+Jxy*Jxz)*ttxyi + &
       (Jxx*Jzz+Jxz*Jxz)*ttxzi + Jxy*Jyz*ttyyi + &
       (Jxy*Jzz+Jxz*Jyz)*ttyzi + Jxz*Jzz*ttzzi)*fac
  ttyy(i,j,k) = (Jxy*Jxy*ttxxi + 2.d0*Jxy*Jyy*ttxyi + &
       2.d0*Jxy*Jyz*ttxzi + Jyy*Jyy*ttyyi + &
       2.d0*Jyy*Jyz*ttyzi + Jyz*Jyz*ttzzi)*fac
  ttyz(i,j,k) = (Jxy*Jxz*ttxxi + (Jxy*Jyz+Jyy*Jxz)*ttxyi + &
       (Jxy*Jzz+Jyz*Jxz)*ttxzi + Jyy*Jyz*ttyyi + &
       (Jyy*Jzz+Jyz*Jyz)*ttyzi + Jyz*Jzz*ttzzi)*fac
  ttzz(i,j,k) = (Jxz*Jxz*ttxxi + 2.d0*Jxz*Jyz*ttxyi + &
       2.d0*Jxz*Jzz*ttxzi + Jyz*Jyz*ttyyi + &
       2.d0*Jyz*Jzz*ttyzi + Jzz*Jzz*ttzzi)*fac
end subroutine trans_fish_phys_ttij
subroutine trans_fish_phys_ttij_new(ex, X, Y, Z, PhysicalR, RadiusDerivative, &
     ttxx,ttxy,ttxz,ttyy,ttyz,ttzz,i,j,k,fish_to_phys)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: PhysicalR,RadiusDerivative
  real*8, dimension(ex(1),ex(2),ex(3))     :: ttxx,ttxy,ttxz,ttyy,ttyz,ttzz
  real*8                                   :: Rp,dr,r,xi,yi,zi,fac
  real*8                                   :: Jxx,Jxy,Jxz,Jyy,Jyz,Jzz
  real*8                                   :: ttxxi,ttxyi,ttxzi,ttyyi,ttyzi,ttzzi
  real*8                                   :: rpr,F2o3,F4o3,rinv,rpinv,drinv,nx,ny,nz
  integer                                  :: i,j,k
  logical                                  :: fish_to_phys
  !
  xi = X(i,1,1)
  yi = Y(1,j,1)
  zi = Z(1,1,k)
  ttxxi = ttxx(i,j,k)
  ttxyi = ttxy(i,j,k)
  ttxzi = ttxz(i,j,k)
  ttyyi = ttyy(i,j,k)
  ttyzi = ttyz(i,j,k)
  ttzzi = ttzz(i,j,k)
  r = sqrt(xi**2 + yi**2 + zi**2)
  rinv=1.0/r
  nx=xi*rinv
  ny=yi*rinv
  nz=zi*rinv
  Rp = PhysicalR(i,j,k)
  rpinv=1.0/Rp
  dR = RadiusDerivative(i,j,k)
  ! transformation matrix
  if (fish_to_phys) then
     drinv=1.0d0/dR
     rpr=r*rpinv
     fac = drinv-rpr
     Jxx = rpr + nx * nx * fac
     Jxy = nx * ny * fac
     Jxz = nx * nz * fac
     Jyy = rpr + ny * ny * fac
     Jyz = ny * nz * fac
     Jzz = rpr + nz * nz * fac
     fac = (rpr * rpr * dRinv)**(-2.d0/3.d0)
     !     write(*,*) hi transforming to phys.,dR,Rp,r,ttxxi,xi,fac,Jxx,Jxy,Jxz,Jyy,Jyz,Jzz,rpr,rinv,drinv,rpinv
     !      write(*,*) hi transforming to phys.,1./rpr,dR,rpr**(-4.d0/3.d0), dR**(2.d0/3.d0),fac
  else
     !write(*,*) hi2.
     rpr=Rp*rinv
     fac = dR - rpr
     Jxx = rpr + nx * nx * fac
     Jxy = nx * ny * fac
     Jxz = nx * nz * fac
     Jyy = rpr + ny * ny * fac
     Jyz = ny * nz * fac
     Jzz = rpr + nz * nz * fac
     fac = (rpr * rpr * dR)**(-2.d0/3.d0)
  end if
  ttxx(i,j,k) = (Jxx*Jxx*ttxxi + 2.d0*Jxx*Jxy*ttxyi + 2.d0*Jxx*Jxz*ttxzi + &
       Jxy*Jxy*ttyyi + 2.d0*Jxy*Jxz*ttyzi + Jxz*Jxz*ttzzi)*fac
  ttxy(i,j,k) = (Jxx*Jxy*ttxxi + (Jxx*Jyy+Jxy*Jxy)*ttxyi + &
       (Jxx*Jyz+Jxz*Jxy)*ttxzi + Jxy*Jyy*ttyyi + &
       (Jxy*Jyz+Jxz*Jyy)*ttyzi + Jxz*Jyz*ttzzi)*fac
  ttxz(i,j,k) = (Jxx*Jxz*ttxxi + (Jxx*Jyz+Jxy*Jxz)*ttxyi + &
       (Jxx*Jzz+Jxz*Jxz)*ttxzi + Jxy*Jyz*ttyyi + &
       (Jxy*Jzz+Jxz*Jyz)*ttyzi + Jxz*Jzz*ttzzi)*fac
  ttyy(i,j,k) = (Jxy*Jxy*ttxxi + 2.d0*Jxy*Jyy*ttxyi + &
       2.d0*Jxy*Jyz*ttxzi + Jyy*Jyy*ttyyi + &
       2.d0*Jyy*Jyz*ttyzi + Jyz*Jyz*ttzzi)*fac
  ttyz(i,j,k) = (Jxy*Jxz*ttxxi + (Jxy*Jyz+Jyy*Jxz)*ttxyi + &
       (Jxy*Jzz+Jyz*Jxz)*ttxzi + Jyy*Jyz*ttyyi + &
       (Jyy*Jzz+Jyz*Jyz)*ttyzi + Jyz*Jzz*ttzzi)*fac
  ttzz(i,j,k) = (Jxz*Jxz*ttxxi + 2.d0*Jxz*Jyz*ttxyi + &
       2.d0*Jxz*Jzz*ttxzi + Jyz*Jyz*ttyyi + &
       2.d0*Jyz*Jzz*ttyzi + Jzz*Jzz*ttzzi)*fac
end subroutine trans_fish_phys_ttij_new
!
!-----------------------------------------------------------------------------
!
!  Numerical Recipe interpolation subroutine 
!
!-----------------------------------------------------------------------------
!
SUBROUTINE polint(xa,ya,n,x,y,dy)
  IMPLICIT NONE
  integer                    :: n,NMAX,i,m,ns
  real*8                         :: dy,x,y
  real*8, dimension(n)         :: xa,ya
  PARAMETER (NMAX=10)
  real*8, dimension(NMAX)    :: c,d
  real*8                     :: den,dif,dift,ho,hp,w
  !
  ns=1
  dif=abs(x-xa(1))
  do 11 i=1,n
     dift=abs(x-xa(i))
     if (dift.lt.dif) then
        ns=i
        dif=dift
     endif
     c(i)=ya(i)
     d(i)=ya(i)
11   continue
     y=ya(ns)
     ns=ns-1
     do 13 m=1,n-1
        do 12 i=1,n-m
           ho=xa(i)-x
           hp=xa(i+m)-x
           w=c(i+1)-d(i)
           den=ho-hp
           if(den.eq.0.) then
              write(*,*) 'failure in polint'
              stop
           end if
           den=w/den
           d(i)=hp*den
           c(i)=ho*den
12         continue
           if (2*ns.lt.n-m)then
              dy=c(ns+1)
           else
              dy=d(ns)
              ns=ns-1
           endif
           y=y+dy
13         continue
           return
         END SUBROUTINE polint
