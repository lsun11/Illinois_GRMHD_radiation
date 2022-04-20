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
  real*8, allocatable,dimension(:)            :: u0_p_int,vx_p_int,vy_p_int,vz_p_int
  real*8, allocatable,dimension(:)            :: u0_int,vx_int,vy_int,vz_int
  integer                                     :: vindex
  ext = cctk_lsh
  global_ext = cctk_gsh

  dT = CCTK_DELTA_TIME
  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  allocate(u0_p_int(narr))
  allocate(vx_p_int(narr))
  allocate(vy_p_int(narr))  
  allocate(vz_p_int(narr))  

  allocate(u0_int(narr)) 
  allocate(vx_int(narr))
  allocate(vy_int(narr))
  allocate(vz_int(narr)) 
  
  
  if (CCTK_TIME .ne. 0) then
     call CCTK_VarIndex(vindex,"OS_collapse::u0_p")  
     call interp_driver_carp(cctkGH,narr,pos,vindex,u0_p_int)
     call CCTK_VarIndex(vindex,"OS_collapse::vx_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vx_p_int)
     call CCTK_VarIndex(vindex,"OS_collapse::vy_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vy_p_int)
     call CCTK_VarIndex(vindex,"OS_collapse::vz_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vz_p_int)

     
     slope(:,1)=1.d0/u0_p_int
     slope(:,2)=vx_p_int
     slope(:,3)=vy_p_int
     slope(:,4)=vz_p_int

     ! call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !      ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !      Xglobal,Yglobal,Zglobal,pos,1.d0/u0_p,slope(:,1))
     ! call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !      ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !      Xglobal,Yglobal,Zglobal,pos,vx_p,slope(:,2))
     ! call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !      ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !      Xglobal,Yglobal,Zglobal,pos,vy_p,slope(:,3)) 
     ! call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !      ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !      Xglobal,Yglobal,Zglobal,pos,vz_p,slope(:,4)) 

     dTo2=dT*0.5d0
     dTo6=dT/6.0d0
     coordt=coord+dTo2*slope
     write(*,*) "dTo2: ",dTo2
    
     
     do ii=1,3
        do jj=1,narr 
           pos(jj,ii) = coordt(jj,ii+1)
        end do
     end do


     call CCTK_VarIndex(vindex,"OS_collapse::u0_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,u0_p_int)
     call CCTK_VarIndex(vindex,"OS_collapse::vx_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vx_p_int)
     call CCTK_VarIndex(vindex,"OS_collapse::vy_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vy_p_int)
     call CCTK_VarIndex(vindex,"OS_collapse::vz_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vz_p_int)


     call CCTK_VarIndex(vindex,"mhd_evolve::u0")
     call interp_driver_carp(cctkGH,narr,pos,vindex,u0_int)
     call CCTK_VarIndex(vindex,"mhd_evolve::vx")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vx_int)
     call CCTK_VarIndex(vindex,"mhd_evolve::vy")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vy_int)
     call CCTK_VarIndex(vindex,"mhd_evolve::vz")
     
     slopet(:,1)=0.5d0*(1.d0/u0_p_int+1.d0/u0_int)
     slopet(:,2)=0.5d0*(vx_p_int+vx_int)
     slopet(:,3)=0.5d0*(vy_p_int+vy_int)
     slopet(:,4)=0.5d0*(vz_p_int+vz_int)
         
     ! call interp_driver_carp(cctkGH,narr,pos,vindex,vz_p_int)    
     ! call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !      ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !      Xglobal,Yglobal,Zglobal,pos,0.5d0*(1.d0/u0_p+1.d0/u0),slopet(:,1))
     ! call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !      ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !      Xglobal,Yglobal,Zglobal,pos,0.5d0*(vx_p+vx),slopet(:,2))
     ! call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !      ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !      Xglobal,Yglobal,Zglobal,pos,0.5d0*(vy_p+vy),slopet(:,3)) 
     ! call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !      ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !      Xglobal,Yglobal,Zglobal,pos,0.5d0*(vz_p+vz),slopet(:,4)) 

     coordt=coord+dTo2*slopet

     
     do ii=1,3
        do jj=1,narr 
           pos(jj,ii) = coordt(jj,ii+1)
        end do
     end do
     
     call CCTK_VarIndex(vindex,"OS_collapse::u0_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,u0_p_int)
     call CCTK_VarIndex(vindex,"OS_collapse::vx_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vx_p_int)
     call CCTK_VarIndex(vindex,"OS_collapse::vy_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vy_p_int)
     call CCTK_VarIndex(vindex,"OS_collapse::vz_p")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vz_p_int)

     call CCTK_VarIndex(vindex,"mhd_evolve::u0")
     call interp_driver_carp(cctkGH,narr,pos,vindex,u0_int)
     call CCTK_VarIndex(vindex,"mhd_evolve::vx")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vx_int)
     call CCTK_VarIndex(vindex,"mhd_evolve::vy")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vy_int)
     call CCTK_VarIndex(vindex,"mhd_evolve::vz")

     slopem(:,1)=0.5d0*(1.d0/u0_p_int+1.d0/u0_int)
     slopem(:,2)=0.5d0*(vx_p_int+vx_int)
     slopem(:,3)=0.5d0*(vy_p_int+vy_int)
     slopem(:,4)=0.5d0*(vz_p_int+vz_int)


     !call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !     ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !     Xglobal,Yglobal,Zglobal,pos,0.5d0*(1.d0/u0_p+1.d0/u0),slopem(:,1))
     !call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !     ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !     Xglobal,Yglobal,Zglobal,pos,0.5d0*(vx_p+vx),slopem(:,2))
     !call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !     ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !     Xglobal,Yglobal,Zglobal,pos,0.5d0*(vy_p+vy),slopem(:,3)) 
     !call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !     ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !     Xglobal,Yglobal,Zglobal,pos,0.5d0*(vz_p+vz),slopem(:,4)) 

     coordt=coord+dT*slopem
     slopem=slopet+slopem

     do ii=1,3
        do jj=1,narr 
           pos(jj,ii) = coordt(jj,ii+1)
        end do
     end do

     call CCTK_VarIndex(vindex,"mhd_evolve::u0")
     call interp_driver_carp(cctkGH,narr,pos,vindex,u0_int)
     call CCTK_VarIndex(vindex,"mhd_evolve::vx")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vx_int)
     call CCTK_VarIndex(vindex,"mhd_evolve::vy")
     call interp_driver_carp(cctkGH,narr,pos,vindex,vy_int)
     call CCTK_VarIndex(vindex,"mhd_evolve::vz")

     slopet(:,1)=1.d0/u0_int
     slopet(:,2)=vx_int
     slopet(:,3)=vy_int
     slopet(:,4)=vz_int
     !call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !  ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !  Xglobal,Yglobal,Zglobal,pos,1.d0/u0,slopet(:,1))
     !call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !     ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !     Xglobal,Yglobal,Zglobal,pos,vx,slopet(:,2))
     !call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !     ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !     Xglobal,Yglobal,Zglobal,pos,vy,slopet(:,3)) 
     !call interpolate_pointset(cctkGH,cctk_nghostzones,narr,2,&
     !     ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,&
     !     Xglobal,Yglobal,Zglobal,pos,vz,slopet(:,4)) 



     coord=coord+dTo6*(slope+slopet+2.0d0*slopem)


  endif

  vx_p=vx
  vy_p=vy
  vz_p=vz
  u0_p=u0


end subroutine particle_tracer
