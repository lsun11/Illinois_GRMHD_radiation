!-----------------------------------------------------------------------------
!
! $Id
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
! Compute Aij using (d/dt)\gamma_{ij} = 0
!
!-----------------------------------------------------------------------------
subroutine kset_c_v2(ex, X, Y, Z,  &
     Axx, Axy, Axz, Ayy, Ayz, Azz, trK, &
     gxx, gxy, gxz, gyy, gyz, gzz, &
     betax, betay, betaz, lapse,  &
     Symmetry,betaxx, betaxy, betaxz, &
     betayx, betayy, betayz, &
     betazx, betazy, betazz, &
     gxxx,gxxy,gxxz,gxyx,gxyy,gxyz, &
     gxzx,gxzy,gxzz,gyyx,gyyy,gyyz, &
     gyzx,gyzy,gyzz,gzzx,gzzy,gzzz, &
     gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
     temp)
  implicit none
  interface
     subroutine gderivs_oct(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
     end subroutine gderivs_oct
     subroutine gderivs_eq(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
     end subroutine gderivs_eq
     subroutine gderivs_axi(ex,f,fx,fy,fz,dX,dY,dZ,SYM1,SYM2,SYM3)
       implicit none
       integer, dimension(3)                    :: ex
       real*8, dimension(ex(1),ex(2),ex(3))     :: f,fx,fy,fz
       real*8                                   :: dX,dY,dZ,SYM1,SYM2,SYM3
     end subroutine gderivs_axi
  end interface
!
! Input parameters:
!
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3))        :: Axx, Axy, Axz, Ayy, Ayz, Azz
  real*8, dimension(ex(1),ex(2),ex(3))        :: trK
  real*8, dimension(ex(1),ex(2),ex(3))        :: gxx, gxy, gxz, gyy, gyz, gzz  
  real*8, dimension(ex(1),ex(2),ex(3))        :: betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3))        :: lapse
  integer                              	      :: Symmetry
!
! Other input variables
!
  real*8, dimension(ex(1),ex(2),ex(3))        :: betaxx, betaxy, betaxz
  real*8, dimension(ex(1),ex(2),ex(3))        :: betayx, betayy, betayz
  real*8, dimension(ex(1),ex(2),ex(3))        :: betazx, betazy, betazz
  real*8, dimension(ex(1),ex(2),ex(3))        :: gxxx,gxxy,gxxz,gxyx,gxyy,gxyz
  real*8, dimension(ex(1),ex(2),ex(3))        :: gxzx,gxzy,gxzz,gyyx,gyyy,gyyz
  real*8, dimension(ex(1),ex(2),ex(3))        :: gyzx,gyzy,gyzz,gzzx,gzzy,gzzz
  real*8, dimension(ex(1),ex(2),ex(3))        :: gupxx,gupxy,gupxz,gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3))        :: temp  


!
! Other variables:
!
  integer                                     :: imin, jmin, kmin
  real*8                                      :: dX, dY, dZ
  real*8                                      :: psi, kxxout
  real*8                                      :: F1o3, F1o6, ONE, TWO, FOUR, ZERO
  real*8                                      :: F2o3, SIX, EIGHT, HALF, PI  
  real*8                                      :: SYM, ANTI, Sym_z
  integer    :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter ( ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0, F1o6 = 1.D0/6.D0 )
  parameter ( ZERO = 0.D0, F1o3 = 1.D0/3.D0, F2o3 = 2.D0/3.D0, SIX = 6.D0 )
  parameter ( SYM = 1.D0, ANTI = - 1.D0, EIGHT = 8.D0, HALF = 0.5D0 )
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  PI = acos(-ONE)
  imin = lbound(gxx,1)
  jmin = lbound(gxx,2)
  kmin = lbound(gxx,3)
  dX = X(imin + 1,1,1) - X(imin,1,1)
  dY = Y(1,jmin + 1,1) - Y(1,jmin,1)
  dZ = Z(1,1,kmin + 1) - Z(1,1,kmin)
  if (Z(1,1,kmin+1) .lt. 0.d0) then
     Sym_z = SYM
  else
     Sym_z = ANTI
  end if
!-----------------------------------------------------------------------------
! Compute first derivatives of the shift and physical metric
!-----------------------------------------------------------------------------
  if(Symmetry==OCTANT) then
     call gderivs_oct(ex,betax,betaxx,betaxy,betaxz,dX,dY,dZ,SYM,ANTI,SYM)
     call gderivs_oct(ex,betay,betayx,betayy,betayz,dX,dY,dZ,ANTI,SYM,SYM)
     call gderivs_oct(ex,betaz,betazx,betazy,betazz,dX,dY,dZ,SYM,SYM,ANTI)
     call gderivs_oct(ex,gxx,gxxx,gxxy,gxxz,dX,dY,dZ,SYM ,SYM ,SYM )
     call gderivs_oct(ex,gxy,gxyx,gxyy,gxyz,dX,dY,dZ,ANTI,ANTI,SYM )
     call gderivs_oct(ex,gxz,gxzx,gxzy,gxzz,dX,dY,dZ,ANTI,SYM ,ANTI)
     call gderivs_oct(ex,gyy,gyyx,gyyy,gyyz,dX,dY,dZ,SYM ,SYM ,SYM )  
     call gderivs_oct(ex,gyz,gyzx,gyzy,gyzz,dX,dY,dZ,SYM ,ANTI,ANTI)  
     call gderivs_oct(ex,gzz,gzzx,gzzy,gzzz,dX,dY,dZ,SYM ,SYM ,SYM )
  else if (Symmetry == EQUATORIAL) then
     call gderivs_eq(ex,betax,betaxx,betaxy,betaxz,dX,dY,dZ,SYM,ANTI,SYM)
     call gderivs_eq(ex,betay,betayx,betayy,betayz,dX,dY,dZ,ANTI,SYM,SYM)
     call gderivs_eq(ex,betaz,betazx,betazy,betazz,dX,dY,dZ,SYM,SYM,ANTI)
     call gderivs_eq(ex,gxx,gxxx,gxxy,gxxz,dX,dY,dZ,SYM ,SYM ,SYM )
     call gderivs_eq(ex,gxy,gxyx,gxyy,gxyz,dX,dY,dZ,ANTI,ANTI,SYM )
     call gderivs_eq(ex,gxz,gxzx,gxzy,gxzz,dX,dY,dZ,ANTI,SYM ,ANTI)
     call gderivs_eq(ex,gyy,gyyx,gyyy,gyyz,dX,dY,dZ,SYM ,SYM ,SYM )  
     call gderivs_eq(ex,gyz,gyzx,gyzy,gyzz,dX,dY,dZ,SYM ,ANTI,ANTI)  
     call gderivs_eq(ex,gzz,gzzx,gzzy,gzzz,dX,dY,dZ,SYM ,SYM ,SYM )
  else if (Symmetry == AXISYM) then
     call gderivs_axi(ex,betax,betaxx,betaxy,betaxz,dX,dY,dZ,ANTI,SYM,SYM)
     call gderivs_axi(ex,betay,betayx,betayy,betayz,dX,dY,dZ,ANTI,ANTI,SYM)
     call gderivs_axi(ex,betaz,betazx,betazy,betazz,dX,dY,dZ,SYM,SYM,Sym_z)
!     call gderivs_axi(ex,betax,betaxx,betaxy,betaxz,dX,dY,dZ,ANTI,ANTI,SYM)
!     call gderivs_axi(ex,betay,betayx,betayy,betayz,dX,dY,dZ,ANTI,SYM,SYM)
!     call gderivs_axi(ex,betaz,betazx,betazy,betazz,dX,dY,dZ,SYM,SYM,Sym_z)
     call gderivs_axi(ex,gxx,gxxx,gxxy,gxxz,dX,dY,dZ,SYM ,SYM ,SYM )
     call gderivs_axi(ex,gxy,gxyx,gxyy,gxyz,dX,dY,dZ,SYM,ANTI,SYM )
!     write(*,*) "KSET_C:",gxyy(2,2,2)
     call gderivs_axi(ex,gxz,gxzx,gxzy,gxzz,dX,dY,dZ,ANTI,SYM ,Sym_z)
     call gderivs_axi(ex,gyy,gyyx,gyyy,gyyz,dX,dY,dZ,SYM ,SYM ,SYM )  
     call gderivs_axi(ex,gyz,gyzx,gyzy,gyzz,dX,dY,dZ,ANTI ,ANTI,Sym_z)
     call gderivs_axi(ex,gzz,gzzx,gzzy,gzzz,dX,dY,dZ,SYM ,SYM ,SYM )
  else if (Symmetry == NO_SYMM) then 
     call gderivs_gen(ex,betax,betaxx,betaxy,betaxz,dX,dY,dZ)
     call gderivs_gen(ex,betay,betayx,betayy,betayz,dX,dY,dZ)
     call gderivs_gen(ex,betaz,betazx,betazy,betazz,dX,dY,dZ)
     call gderivs_gen(ex,gxx,gxxx,gxxy,gxxz,dX,dY,dZ)
     call gderivs_gen(ex,gxy,gxyx,gxyy,gxyz,dX,dY,dZ)
     call gderivs_gen(ex,gxz,gxzx,gxzy,gxzz,dX,dY,dZ)
     call gderivs_gen(ex,gyy,gyyx,gyyy,gyyz,dX,dY,dZ)
     call gderivs_gen(ex,gyz,gyzx,gyzy,gyzz,dX,dY,dZ)
     call gderivs_gen(ex,gzz,gzzx,gzzy,gzzz,dX,dY,dZ)
  end if
!----------------------------------------------------------
! Output
!----------------------------------------------------------
  Axx = ( betax*gxxx + betay*gxxy + betaz*gxxz &
       + gxx*betaxx + gxy*betayx + gxz*betazx &
       + gxx*betaxx + gxy*betayx + gxz*betazx) & 
       / (TWO * (lapse + ONE) )

  Axy = ( betax*gxyx + betay*gxyy + betaz*gxyz &
       + gxx*betaxy + gxy*betayy + gxz*betazy &
       + gxy*betaxx + gyy*betayx + gyz*betazx) & 
       / (TWO * (lapse + ONE) )

  Axz = ( betax*gxzx + betay*gxzy + betaz*gxzz &
       + gxx*betaxz + gxy*betayz + gxz*betazz &
       + gxz*betaxx + gyz*betayx + gzz*betazx) & 
       / (TWO * (lapse + ONE) )
  Ayy = ( betax*gyyx + betay*gyyy + betaz*gyyz &
       + gxy*betaxy + gyy*betayy + gyz*betazy &
       + gxy*betaxy + gyy*betayy + gyz*betazy) & 
       / (TWO * (lapse + ONE) )
  Ayz = ( betax*gyzx + betay*gyzy + betaz*gyzz &
       + gxy*betaxz + gyy*betayz + gyz*betazz &
       + gxz*betaxy + gyz*betayy + gzz*betazy) &
       / (TWO * (lapse + ONE) )
  Azz = ( betax*gzzx + betay*gzzy + betaz*gzzz &
       + gxz*betaxz + gyz*betayz + gzz*betazz &
       + gxz*betaxz + gyz*betayz + gzz*betazz) & 
       / (TWO * (lapse + ONE) )

  temp =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz &
       - gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz )/ temp
  gupxy = - ( gxy * gzz - gyz * gxz )/ temp
  gupxz =   ( gxy * gyz - gyy * gxz )/ temp
  gupyy =   ( gxx * gzz - gxz * gxz )/ temp
  gupyz = - ( gxx * gyz - gxy * gxz )/ temp
  gupzz =   ( gxx * gyy - gxy * gxy )/ temp  
  trK = (gupxx*Axx + gupyy*Ayy + gupzz*Azz + & 
       TWO*(gupxy*Axy + gupxz*Axz + gupyz*Ayz))

  Axx = Axx - F1o3 * gxx * trK
  Axy = Axy - F1o3 * gxy * trK
  Axz = Axz - F1o3 * gxz * trK
  Ayy = Ayy - F1o3 * gyy * trK
  Ayz = Ayz - F1o3 * gyz * trK
  Azz = Azz - F1o3 * gzz * trK

! In axisymmetry, the A_ij values off the y/=0 plane will
! need to be corrected. 

!  if (Symmetry == AXISYM) then
!     call axibc_scalar(ex,X,Y,Z,trK)
!     call axibc_tensor(ex,X,Y,Z,Axx,Axy,Axz,Ayy,Ayz,Azz)
!  endif
end subroutine kset_c_v2
