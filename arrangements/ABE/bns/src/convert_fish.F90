!-----------------------------------------------------------------------------
!
! Convert various quantities from physical coords to fisheye
!
!-----------------------------------------------------------------------------
subroutine trans_phys_fish_tensor_flat(ex, x, y, z, &
     Rp,dR,gxx,gxy,gxz,gyy,gyz,gzz,Symmetry)
  
  implicit none
  
  !
  ! Input parameters:
  !
  integer, dimension(3),     intent(in) :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in) :: x,y,z,RP,dR
  integer, intent(in)                   :: Symmetry
  real*8,  dimension(ex(1),ex(2),ex(3)) :: gxx,gxy,gxz,gyy,gyz,gzz

  real*8,  dimension(ex(1),ex(2),ex(3)) :: jxx,jxy,jxz,jyy,jyz,jzz,r,fac

  r = sqrt(x**2 + y**2 + z**2)
  fac = dR - Rp/r
  Jxx = Rp/r + x/r * x/r * fac
  Jxy = x/r * y/r * fac
  Jxz = x/r * z/r * fac
  Jyy = Rp/r + y/r * y/r * fac
  Jyz = y/r * z/r * fac
  Jzz = Rp/r + z/r * z/r * fac
  fac = (Rp/r)**(-4.d0/3.d0) * dR**(-2.d0/3.d0)

  gxx = (Jxx*Jxx + Jxy*Jxy + Jxz*Jxz)*fac
  gxy = (Jxx*Jxy + Jxy*Jyy + Jxz*Jyz)*fac
  gxz = (Jxx*Jxz + Jxy*Jyz + Jxz*Jzz)*fac
  gyy = (Jxy*Jxy + Jyy*Jyy + Jyz*Jyz)*fac
  gyz = (Jxy*Jxz + Jyy*Jyz + Jyz*Jzz)*fac
  gzz = (Jxz*Jxz + Jyz*Jyz + Jzz*Jzz)*fac

  return
end subroutine trans_phys_fish_tensor_flat

subroutine trans_phys_fish_tensor(ex, x, y, z, &
     RP,dR,gxx,gxy,gxz,gyy,gyz,gzz,Symmetry)
  
  implicit none
  
  !
  ! Input parameters:
  !
  integer, dimension(3),     intent(in) :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in) :: x,y,z,RP,dR
  integer, intent(in)                   :: Symmetry
  real*8,  dimension(ex(1),ex(2),ex(3)) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8,  dimension(ex(1),ex(2),ex(3)) :: gxxi,gxyi,gxzi,gyyi,gyzi,gzzi
  real*8,  dimension(ex(1),ex(2),ex(3)) :: r,fac,jxx,jxy,jxz,jyy,jyz,jzz

  r = sqrt(x**2 + y**2 + z**2)
  fac = dR - Rp/r
  Jxx = Rp/r + x/r * x/r * fac
  Jxy = x/r * y/r * fac
  Jxz = x/r * z/r * fac
  Jyy = Rp/r + y/r * y/r * fac
  Jyz = y/r * z/r * fac
  Jzz = Rp/r + z/r * z/r * fac
  fac = (Rp/r)**(-4.d0/3.d0) * dR**(-2.d0/3.d0)

  gxxi=gxx
  gyyi=gyy
  gzzi=gzz
  gxyi=gxy
  gxzi=gxz
  gyzi=gyz

  gxx = (Jxx*Jxx*gxxi + 2.d0*Jxx*Jxy*gxyi + 2.d0*Jxx*Jxz*gxzi + &
                Jxy*Jxy*gyyi + 2.d0*Jxy*Jxz*gyzi + Jxz*Jxz*gzzi)*fac
  gxy = (Jxx*Jxy*gxxi + (Jxx*Jyy+Jxy*Jxy)*gxyi + &
                (Jxx*Jyz+Jxz*Jxy)*gxzi + Jxy*Jyy*gyyi + &
                (Jxy*Jyz+Jxz*Jyy)*gyzi + Jxz*Jyz*gzzi)*fac
  gxz = (Jxx*Jxz*gxxi + (Jxx*Jyz+Jxy*Jxz)*gxyi + &
                (Jxx*Jzz+Jxz*Jxz)*gxzi + Jxy*Jyz*gyyi + &
                (Jxy*Jzz+Jxz*Jyz)*gyzi + Jxz*Jzz*gzzi)*fac
  gyy = (Jxy*Jxy*gxxi + 2.d0*Jxy*Jyy*gxyi + &
                2.d0*Jxy*Jyz*gxzi + Jyy*Jyy*gyyi + &
                2.d0*Jyy*Jyz*gyzi + Jyz*Jyz*gzzi)*fac
  gyz = (Jxy*Jxz*gxxi + (Jxy*Jyz+Jyy*Jxz)*gxyi + &
                (Jxy*Jzz+Jyz*Jxz)*gxzi + Jyy*Jyz*gyyi + &
                (Jyy*Jzz+Jyz*Jyz)*gyzi + Jyz*Jzz*gzzi)*fac
  gzz = (Jxz*Jxz*gxxi + 2.d0*Jxz*Jyz*gxyi + &
                2.d0*Jxz*Jzz*gxzi + Jyz*Jyz*gyyi + &
                2.d0*Jyz*Jzz*gyzi + Jzz*Jzz*gzzi)*fac
  
  return
end subroutine trans_phys_fish_tensor

subroutine trans_phys_fish_tensor_inv(ex, x, y, z, &
     RP,dR,gxx,gxy,gxz,gyy,gyz,gzz,Symmetry)
  
  implicit none
  
  !
  ! Input parameters:
  !
  integer, dimension(3),     intent(in) :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in) :: x,y,z,RP,dR
  integer, intent(in)                   :: Symmetry
  real*8,  dimension(ex(1),ex(2),ex(3)) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8,  dimension(ex(1),ex(2),ex(3)) :: gxxi,gxyi,gxzi,gyyi,gyzi,gzzi
  real*8,  dimension(ex(1),ex(2),ex(3)) :: r,fac,jxx,jxy,jxz,jyy,jyz,jzz

  r = sqrt(x**2 + y**2 + z**2)
  fac = 1.0/dR - r/Rp
  Jxx = r/Rp + x/r * x/r * fac
  Jxy = x/r * y/r * fac
  Jxz = x/r * z/r * fac
  Jyy = r/Rp + y/r * y/r * fac
  Jyz = y/r * z/r * fac
  Jzz = r/Rp + z/r * z/r * fac
  fac = (r/Rp)**(-4.d0/3.d0) * dR**(2.d0/3.d0)

  gxxi=gxx
  gyyi=gyy
  gzzi=gzz
  gxyi=gxy
  gxzi=gxz
  gyzi=gyz

  gxx = (Jxx*Jxx*gxxi + 2.d0*Jxx*Jxy*gxyi + 2.d0*Jxx*Jxz*gxzi + &
                Jxy*Jxy*gyyi + 2.d0*Jxy*Jxz*gyzi + Jxz*Jxz*gzzi)*fac
  gxy = (Jxx*Jxy*gxxi + (Jxx*Jyy+Jxy*Jxy)*gxyi + &
                (Jxx*Jyz+Jxz*Jxy)*gxzi + Jxy*Jyy*gyyi + &
                (Jxy*Jyz+Jxz*Jyy)*gyzi + Jxz*Jyz*gzzi)*fac
  gxz = (Jxx*Jxz*gxxi + (Jxx*Jyz+Jxy*Jxz)*gxyi + &
                (Jxx*Jzz+Jxz*Jxz)*gxzi + Jxy*Jyz*gyyi + &
                (Jxy*Jzz+Jxz*Jyz)*gyzi + Jxz*Jzz*gzzi)*fac
  gyy = (Jxy*Jxy*gxxi + 2.d0*Jxy*Jyy*gxyi + &
                2.d0*Jxy*Jyz*gxzi + Jyy*Jyy*gyyi + &
                2.d0*Jyy*Jyz*gyzi + Jyz*Jyz*gzzi)*fac
  gyz = (Jxy*Jxz*gxxi + (Jxy*Jyz+Jyy*Jxz)*gxyi + &
                (Jxy*Jzz+Jyz*Jxz)*gxzi + Jyy*Jyz*gyyi + &
                (Jyy*Jzz+Jyz*Jyz)*gyzi + Jyz*Jzz*gzzi)*fac
  gzz = (Jxz*Jxz*gxxi + 2.d0*Jxz*Jyz*gxyi + &
                2.d0*Jxz*Jzz*gxzi + Jyz*Jyz*gyyi + &
                2.d0*Jyz*Jzz*gyzi + Jzz*Jzz*gzzi)*fac
  
  return
end subroutine trans_phys_fish_tensor_inv

subroutine trans_phys_fish_phi(ex, x, y, z, &
     RP,RD,phi,Symmetry)
  
  implicit none
  
  !
  ! Input parameters:
  !
  integer, dimension(3),     intent(in) :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in) :: x,y,z,RP,RD
  integer, intent(in)                   :: Symmetry
  real*8,  dimension(ex(1),ex(2),ex(3)) :: phi
  real*8,  dimension(ex(1),ex(2),ex(3)) :: r
  
  r=sqrt(x**2+y**2+z**2)

  phi=phi+1.d0/3.d0*log(Rp/r)+1.d0/6.d0*log(Rd)

  return
end subroutine trans_phys_fish_phi

subroutine trans_phys_fish_gamt_flat(ex, x, y, z, &
     RP,RD,RD2,gx,gy,gz,Symmetry)
  
  implicit none
  
  !
  ! Input parameters:
  !
  integer, dimension(3),     intent(in) :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in) :: x,y,z,RP,RD,RD2
  integer, intent(in)                   :: Symmetry
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: gx,gy,gz
 
  real*8,  dimension(ex(1),ex(2),ex(3)) :: r,nx,ny,nz,juxx,juyy,juzz,juxy,juxz,juyz
  real*8,  dimension(ex(1),ex(2),ex(3)) :: p1,dxp,dyp,dzp,g2,g3
  real*8,  dimension(ex(1),ex(2),ex(3)) :: jxxx,jxxy,jxxz,jyyx,jyyy,jyyz,jzzx,jzzy,jzzz
  real*8,  dimension(ex(1),ex(2),ex(3)) :: jxyx,jxyy,jxyz,jxzx,jxzy,jxzz,jyzx,jyzy,jyzz
!
  r=sqrt(x**2+y**2+z**2)
  nx=x/r
  ny=y/r
  nz=z/r
 
  !
  !   Now the Jacobian terms (we need to transform the downstairs version)
  !     
  juxx=r/RP+nx**2*(1.0d0/RD-r/RP)
  juyy=r/RP+ny**2*(1.0d0/RD-r/RP)
  juzz=r/RP+nz**2*(1.0d0/RD-r/RP)
  juxy=nx*ny*(1.0d0/RD-r/RP)
  juxz=nx*nz*(1.0d0/RD-r/RP)
  juyz=ny*nz*(1.0d0/RD-r/RP)
       
  ! (phii/phii_p)^4 = (RP/r)^4/3 * (dr_p/dr_fish)^2/3
  p1=(RP/r)**(4.d0/3.d0) * RD**(2.d0/3.d0)
  ! dxp is /partial_{x-bar} (RP/r)
  dxp=nx/r*(RD-RP/r)
  dyp=ny/r*(RD-RP/r)
  dzp=nz/r*(RD-RP/r)
  dxp=2.0d0/3.0d0*p1*(2.0*dxp*r/RP+RD2/RD*nx)
  dyp=2.0d0/3.0d0*p1*(2.0*dyp*r/RP+RD2/RD*ny)
  dzp=2.0d0/3.0d0*p1*(2.0*dzp*r/RP+RD2/RD*nz)

  !
  ! and now the Gamma's
  !                   =del_ij       +del_ik  +del_jk  -terms
  g2=RP/r-RD
  g3=1.0/RD-r/RP
  
  jxxx=nx*r/(RP**2)*g2+nx/r*g3+nx/r*g3-2*nx*nx*nx/r*g3 &
       -nx*nx*nx*r/(RP**2)*g2-nx*nx*nx*RD2/RD**2
  jxxy=ny*r/(RP**2)*g2                  -2*nx*nx*ny/r*g3 &
       -nx*nx*ny*r/(RP**2)*g2-nx*nx*ny*RD2/RD**2
  jxxz=nz*r/(RP**2)*g2                  -2*nx*nx*nz/r*g3 &
       -nx*nx*nz*r/(RP**2)*g2-nx*nx*nz*RD2/RD**2
  
  jyyx=nx*r/(RP**2)*g2                  -2*ny*ny*nx/r*g3 &
       -ny*ny*nx*r/(RP**2)*g2-ny*ny*nx*RD2/RD**2
  jyyy=ny*r/(RP**2)*g2+ny/r*g3+ny/r*g3-2*ny*ny*ny/r*g3 &
       -ny*ny*ny*r/(RP**2)*g2-ny*ny*ny*RD2/RD**2
  jyyz=nz*r/(RP**2)*g2                  -2*ny*ny*nz/r*g3 &
       -ny*ny*nz*r/(RP**2)*g2-ny*ny*nz*RD2/RD**2
  
  jzzx=nx*r/(RP**2)*g2                  -2*nz*nz*nx/r*g3 &
       -nz*nz*nx*r/(RP**2)*g2-nz*nz*nx*RD2/RD**2
  jzzy=ny*r/(RP**2)*g2                  -2*nz*nz*ny/r*g3 &
       -nz*nz*ny*r/(RP**2)*g2-nz*nz*ny*RD2/RD**2
  jzzz=nz*r/(RP**2)*g2+nz/r*g3+nz/r*g3-2*nz*nz*nz/r*g3 &
       -nz*nz*nz*r/(RP**2)*g2-nz*nz*nz*RD2/RD**2
  
  jxyx=                       ny/r*g3-2*nx*ny*nx/r*g3 &
       -nx*ny*nx*r/(RP**2)*g2-nx*ny*nx*RD2/RD**2 
  jxyy=                       nx/r*g3-2*nx*ny*ny/r*g3 &
       -nx*ny*ny*r/(RP**2)*g2-nx*ny*ny*RD2/RD**2 
  jxyz=                               -2*nx*ny*nz/r*g3 &
       -nx*ny*nz*r/(RP**2)*g2-nx*ny*nz*RD2/RD**2 
  
  jxzx=                       nz/r*g3-2*nx*nz*nx/r*g3 &
       -nx*nz*nx*r/(RP**2)*g2-nx*nz*nx*RD2/RD**2 
  jxzy=                               -2*nx*nz*ny/r*g3 &
       -nx*nz*ny*r/(RP**2)*g2-nx*nz*ny*RD2/RD**2 
  jxzz=                       nx/r*g3-2*nx*nz*nz/r*g3 &
       -nx*nz*nz*r/(RP**2)*g2-nx*nz*nz*RD2/RD**2 
  
  jyzx=                               -2*ny*nz*nx/r*g3 &
       -ny*nz*nx*r/(RP**2)*g2-ny*nz*nx*RD2/RD**2 
  jyzy=                       nz/r*g3-2*ny*nz*ny/r*g3 &
       -ny*nz*ny*r/(RP**2)*g2-ny*nz*ny*RD2/RD**2 
  jyzz=                       ny/r*g3-2*ny*nz*nz/r*g3 &
       -ny*nz*nz*r/(RP**2)*g2-ny*nz*nz*RD2/RD**2 
  
  gx=-1.d0*(dxp*(juxx**2+juxy**2+juxz**2)+ &
       2.0*p1*(juxx*jxxx+juxy*jxyx+juxz*jxzx)+ &
       dyp*(juxx*juxy+juxy*juyy+juxz*juyz)+ &
       p1*(juxx*jxyy+juxy*jyyy+juxz*jyzy+ &
       jxxy*juxy+jxyy*juyy+jxzy*juyz)+ &
       dzp*(juxx*juxz+juxy*juyz+juxz*juzz)+ &
       p1*(juxx*jxzz+juxy*jyzz+juxz*jzzz+ &
       jxxz*juxz+jxyz*juyz+jxzz*juzz))
  gy=-1.d0*(dxp*(juxx*juxy+juxy*juyy+juxz*juyz)+ &
       p1*(juxx*jxyx+juxy*jyyx+juxz*jyzx+ &
       jxxx*juxy+jxyx*juyy+jxzx*juyz)+ &
       dyp*(juxy**2+juyy**2+juyz**2)+ &
       2.0*p1*(juxy*jxyy+juyy*jyyy+juyz*jyzy)+ &
       dzp*(juxy*juxz+juyy*juyz+juyz*juzz)+ &
       p1*(juxy*jxzz+juyy*jyzz+juyz*jzzz+ &
       jxyz*juxz+jyyz*juyz+jyzz*juzz))
  gz=-1.d0*(dxp*(juxx*juxz+juxy*juyz+juxz*juzz)+ &
       p1*(juxx*jxzx+juxy*jyzx+juxz*jzzx+ &
       jxxx*juxz+jxyx*juyz+jxzx*juzz)+ &
       dyp*(juxy*juxz+juyy*juyz+juyz*juzz)+ &
       p1*(juxy*jxzy+juyy*jyzy+juyz*jzzy+ &
       jxyy*juxz+jyyy*juyz+jyzy*juzz)+ &
       dzp*(juxz**2+juyz**2+juzz**2)+ &
       2.0*p1*(juxz*jxzz+juyz*jyzz+juzz*jzzz)) 

  return

end subroutine trans_phys_fish_gamt_flat
