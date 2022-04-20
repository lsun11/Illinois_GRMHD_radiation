#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine particle_tracer(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  integer, dimension(3)                       :: ext,global_ext
  real*8                                      :: dT,dTo2,dTo6,dX,dY,dZ,xmax
  integer                                     :: ii,jj
  ext = cctk_lsh
  global_ext = cctk_gsh
  
  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  if (CCTK_TIME .ne. 0) then
     call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
          ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
          Xglobal,Yglobal,Zglobal,pos,1.d0/u0_p,slope(:,1))
     call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
          ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
          Xglobal,Yglobal,Zglobal,pos,vx_p,slope(:,2))
     call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
          ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
          Xglobal,Yglobal,Zglobal,pos,vy_p,slope(:,3)) 
     call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
          ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
          Xglobal,Yglobal,Zglobal,pos,vz_p,slope(:,4)) 
     
     dTo2=dT*0.5d0
     dTo6=dT/6.0d0
     coordt=coord+dTo2*slope
     
     do ii=1,3
        do jj=1,narr 
           pos(jj,ii) = coordt(jj,ii+1)
        end do
     end do
     
 

  call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
       ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
       Xglobal,Yglobal,Zglobal,pos,0.5d0*(1.d0/u0_p+1.d0/u0),slopet(:,1))
  call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
       ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
       Xglobal,Yglobal,Zglobal,pos,0.5d0*(vx_p+vx),slopet(:,2))
  call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
       ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
       Xglobal,Yglobal,Zglobal,pos,0.5d0*(vy_p+vy),slopet(:,3)) 
  call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
       ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
       Xglobal,Yglobal,Zglobal,pos,0.5d0*(vz_p+vz),slopet(:,4)) 
  
  coordt=coord+dTo2*slopet


  do ii=1,3
     do jj=1,narr 
        pos(jj,ii) = coordt(jj,ii+1)
     end do
  end do
 
  
  call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
       ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
       Xglobal,Yglobal,Zglobal,pos,0.5d0*(1.d0/u0_p+1.d0/u0),slopem(:,1))
  call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
       ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
       Xglobal,Yglobal,Zglobal,pos,0.5d0*(vx_p+vx),slopem(:,2))
  call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
       ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
       Xglobal,Yglobal,Zglobal,pos,0.5d0*(vy_p+vy),slopem(:,3)) 
  call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
       ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
       Xglobal,Yglobal,Zglobal,pos,0.5d0*(vz_p+vz),slopem(:,4)) 
  
  coordt=coord+dT*slopem
  slopem=slopet+slopem

  do ii=1,3
     do jj=1,narr 
        pos(jj,ii) = coordt(jj,ii+1)
     end do
  end do
 
  call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
       ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
       Xglobal,Yglobal,Zglobal,pos,1.d0/u0,slopet(:,1))
  call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
       ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
       Xglobal,Yglobal,Zglobal,pos,vx,slopet(:,2))
  call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
       ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
       Xglobal,Yglobal,Zglobal,pos,vy,slopet(:,3)) 
  call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
       ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
       Xglobal,Yglobal,Zglobal,pos,vz,slopet(:,4)) 



  coord=coord+dTo6*(slope+slopet+2.0d0*slopem)
  

endif

vx_p=vx
vy_p=vy
vz_p=vz
u0_p=u0


end subroutine particle_tracer
