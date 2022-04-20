! Note:  in this version, it is not necessary to multiply by powers
! of the determinant.  (I.e. this applies to tensor-densities of weight
! zero.)
! Note2: this function is really not optimized.  Do NOT use outside of initial data routines.
subroutine trans_fish_phys_ttij2_fullgridfunction(ex, X, Y, Z, PhysicalR, RadiusDerivative, &
     ttxx,ttxy,ttxz,ttyy,ttyz,ttzz,fish_to_phys)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: X, Y, Z
  real*8, dimension(ex(1),ex(2),ex(3))     :: PhysicalR,RadiusDerivative
  real*8, dimension(ex(1),ex(2),ex(3))     :: ttxx,ttxy,ttxz,ttyy,ttyz,ttzz
  real*8                                   :: Rp,dr,r,xi,yi,zi,fac
  real*8                                   :: Jxx,Jxy,Jxz,Jyy,Jyz,Jzz
  real*8                                   :: ttxxi,ttxyi,ttxzi,ttyyi,ttyzi,ttzzi
  integer				   :: i,j,k
  logical				   :: fish_to_phys
  !

  do k=1,ex(3)
     do j=1,ex(2)
        do i=1,ex(1)

           xi = x(i,j,k)
           yi = y(i,j,k)
           zi = z(i,j,k)
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
              Jxz = xi/r * zi/r * fac;
              Jyy = r/Rp + yi/r * yi/r * fac
              Jyz = yi/r * zi/r * fac
              Jzz = r/Rp + zi/r * zi/r * fac
           else 
              fac = dR - Rp/r
              Jxx = Rp/r + xi/r * xi/r * fac
              Jxy = xi/r * yi/r * fac
              Jxz = xi/r * zi/r * fac;
              Jyy = Rp/r + yi/r * yi/r * fac
              Jyz = yi/r * zi/r * fac
              Jzz = Rp/r + zi/r * zi/r * fac
           end if


           ttxx(i,j,k) = (Jxx*Jxx*ttxxi + 2.d0*Jxx*Jxy*ttxyi + 2.d0*Jxx*Jxz*ttxzi + &
                Jxy*Jxy*ttyyi + 2.d0*Jxy*Jxz*ttyzi + Jxz*Jxz*ttzzi)
           ttxy(i,j,k) = (Jxx*Jxy*ttxxi + (Jxx*Jyy+Jxy*Jxy)*ttxyi + &
                (Jxx*Jyz+Jxz*Jxy)*ttxzi + Jxy*Jyy*ttyyi + &
                (Jxy*Jyz+Jxz*Jyy)*ttyzi + Jxz*Jyz*ttzzi)
           ttxz(i,j,k) = (Jxx*Jxz*ttxxi + (Jxx*Jyz+Jxy*Jxz)*ttxyi + &
                (Jxx*Jzz+Jxz*Jxz)*ttxzi + Jxy*Jyz*ttyyi + &
                (Jxy*Jzz+Jxz*Jyz)*ttyzi + Jxz*Jzz*ttzzi)
           ttyy(i,j,k) = (Jxy*Jxy*ttxxi + 2.d0*Jxy*Jyy*ttxyi + &
                2.d0*Jxy*Jyz*ttxzi + Jyy*Jyy*ttyyi + &
                2.d0*Jyy*Jyz*ttyzi + Jyz*Jyz*ttzzi)
           ttyz(i,j,k) = (Jxy*Jxz*ttxxi + (Jxy*Jyz+Jyy*Jxz)*ttxyi + &
                (Jxy*Jzz+Jyz*Jxz)*ttxzi + Jyy*Jyz*ttyyi + &
                (Jyy*Jzz+Jyz*Jyz)*ttyzi + Jyz*Jzz*ttzzi)
           ttzz(i,j,k) = (Jxz*Jxz*ttxxi + 2.d0*Jxz*Jyz*ttxyi + &
                2.d0*Jxz*Jzz*ttxzi + Jyz*Jyz*ttyyi + &
                2.d0*Jyz*Jzz*ttyzi + Jzz*Jzz*ttzzi)

        end do
     end do
  end do

end subroutine trans_fish_phys_ttij2_fullgridfunction

! Warning: this function is really not optimized.  Do NOT use outside of initial data routines.
subroutine trans_fish_phys_vector_fullgridfunction(ex, X, Y, Z, r_fish, PhysicalR, RadiusDerivative, &
     vx,vy,vz,fish_to_phys)
  implicit none
  integer, dimension(3)                    :: ex
  real*8, dimension(ex(1),ex(2),ex(3))     :: X,Y,Z,r_fish,PhysicalR,RadiusDerivative
  real*8, dimension(ex(1),ex(2),ex(3))     :: vx,vy,vz
  real*8                                   :: Rp,dr,r,xi,yi,zi,fac
  real*8                                   :: Jxx,Jxy,Jxz,Jyy,Jyz,Jzz
  real*8                                   :: vxi,vyi,vzi
  integer                                  :: i,j,k
  logical                                  :: fish_to_phys
!
  
  do k=1,ex(3)
     do j=1,ex(2)
        do i=1,ex(1)
           xi = X(i,1,1)
           yi = Y(1,j,1)
           zi = Z(1,1,k)
           vxi = vx(i,j,k)
           vyi = vy(i,j,k)
           vzi = vz(i,j,k)
           r = r_fish(i,j,k)
           Rp = PhysicalR(i,j,k)
           dR = RadiusDerivative(i,j,k)
           
           ! transformation matrix
           if (fish_to_phys) then
              fac = dR - Rp/r
              Jxx = Rp/r + xi/r * xi/r * fac
              Jxy = xi/r * yi/r * fac
              Jxz = xi/r * zi/r * fac;
              Jyy = Rp/r + yi/r * yi/r * fac
              Jyz = yi/r * zi/r * fac
              Jzz = Rp/r + zi/r * zi/r * fac
           else 
              fac = 1.d0/dR - r/Rp
              Jxx = r/Rp + xi/r * xi/r * fac
              Jxy = xi/r * yi/r * fac
              Jxz = xi/r * zi/r * fac;
              Jyy = r/Rp + yi/r * yi/r * fac
              Jyz = yi/r * zi/r * fac
              Jzz = r/Rp + zi/r * zi/r * fac
           end if
           
           
           vx(i,j,k) = Jxx*vxi + Jxy*vyi + Jxz*vzi
           vy(i,j,k) = Jxy*vxi + Jyy*vyi + Jyz*vzi
           vz(i,j,k) = Jxz*vxi + Jyz*vyi + Jzz*vzi
        end do
     end do
  end do
  
end subroutine trans_fish_phys_vector_fullgridfunction

