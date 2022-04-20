#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine convert_fish_matter(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3)                    :: ext

  ext = cctk_lsh

  call trans_phys_fish_matter(ext,X,Y,Z, &
       PhysicalRadius,RadiusDerivative, &
       rho_star,tau,st_x,st_y,st_z,w,vx,vy,vz,Symmetry)

  mhd_st_x=st_x
  mhd_st_y=st_y
  mhd_st_z=st_z

  !  if(bondi_enable.eq.1) then
  !     rho = h*w*exp(-6.0*phi)-P
  !     Sx = st_x*exp(-6.0*phi)
  !     Sy = st_y*exp(-6.0*phi)
  !     Sz = st_z*exp(-6.0*phi)
  !     Sxx = st_x*st_x/w/h*exp(-6.0*phi) + P*gxx*exp(4.0*phi)
  !     Sxy = st_x*st_y/w/h*exp(-6.0*phi) + P*gxy*exp(4.0*phi)
  !     Sxz = st_x*st_z/w/h*exp(-6.0*phi) + P*gxz*exp(4.0*phi)
  !     Syy = st_y*st_y/w/h*exp(-6.0*phi) + P*gyy*exp(4.0*phi)
  !     Syz = st_y*st_z/w/h*exp(-6.0*phi) + P*gyz*exp(4.0*phi)
  !     Szz = st_z*st_z/w/h*exp(-6.0*phi) + P*gzz*exp(4.0*phi)
  !  endif

end subroutine convert_fish_matter

!-----------------------------------------------------------------------------
!
! Convert various quantities from physical coords to fisheye
!
!-----------------------------------------------------------------------------
subroutine trans_phys_fish_matter(ex, xG, yG, zG, &
     RPG,dRG,rho_s,tau,st_x,st_y,st_z,w,vx,vy,vz,Symmetry)

  implicit none  
  !
  ! Input parameters:
  !
  integer, dimension(3),     intent(in) :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in) :: xG,yG,zG,RPG,dRG
  real*8                                            :: x,y,z,RP,dR
  integer, intent(in)                   :: Symmetry
  integer                               :: i,j,k
  real*8,  dimension(ex(1),ex(2),ex(3)) :: rho_s,tau,st_x,st_y,st_z,w,vx,vy,vz
  real*8                                :: vxi,vyi,vzi,r
  real*8                                :: fac,jxx,jxy,jxz,jyy,jyz,jzz


  do k=1,ex(3)
     do j=1,ex(2)
        do i=1,ex(1)
           x = xG(i,j,k)
           y = yG(i,j,k)
           z = zG(i,j,k)

           RP = RPG(i,j,k)
           dR = dRG(i,j,k)

           r = sqrt(x*x + y*y + z*z)
           fac = 1.0/dR - r/Rp
           Jxx = r/Rp + x/r * x/r * fac
           Jxy = x/r * y/r * fac
           Jxz = x/r * z/r * fac
           Jyy = r/Rp + y/r * y/r * fac
           Jyz = y/r * z/r * fac
           Jzz = r/Rp + z/r * z/r * fac

           vxi=vx(i,j,k)
           vyi=vy(i,j,k)
           vzi=vz(i,j,k)

           vx(i,j,k)=jxx*vxi+jxy*vyi+jxz*vzi
           vy(i,j,k)=jxy*vxi+jyy*vyi+jyz*vzi
           vz(i,j,k)=jxz*vxi+jyz*vyi+jzz*vzi

           !  write(6,*)'Boundary0:',vxi(60,1,1),vyi(60,1,1),vzi(60,1,1),jxx(60,1,1),jxy(60,1,1),jxz(60,1,1),vx(60,1,1),vy(60,1,1),vz(60,1,1)

           fac = DR*(RP/r)**2 

           rho_s(i,j,k)=rho_s(i,j,k)*fac
           w(i,j,k)=w(i,j,k)*fac
           tau(i,j,k)=tau(i,j,k)*fac

           fac = dR - Rp/r
           Jxx = Rp/r + x/r * x/r * fac
           Jxy = x/r * y/r * fac
           Jxz = x/r * z/r * fac
           Jyy = Rp/r + y/r * y/r * fac
           Jyz = y/r * z/r * fac
           Jzz = Rp/r + z/r * z/r * fac
           fac = DR*(RP/r)**2 

           vxi=st_x(i,j,k)
           vyi=st_y(i,j,k)
           vzi=st_z(i,j,k)

           st_x(i,j,k)=fac*(jxx*vxi+jxy*vyi+jxz*vzi)
           st_y(i,j,k)=fac*(jxy*vxi+jyy*vyi+jyz*vzi)
           st_z(i,j,k)=fac*(jxz*vxi+jyz*vyi+jzz*vzi)

        end do
     end do
  end do

  return
end subroutine trans_phys_fish_matter
