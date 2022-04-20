! Warning: this function is really not optimized.  Do NOT use outside of initial data routines.
subroutine statpunc_trans_fish_phys_vector_fullgridfunction(ex, X, Y, Z, r_fish, PhysicalR, RadiusDerivative, &
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
  
end subroutine statpunc_trans_fish_phys_vector_fullgridfunction
