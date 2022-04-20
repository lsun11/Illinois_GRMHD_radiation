#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
! Predictor
!
!-----------------------------------------------------------------------------
subroutine Center_Predictor(CCTK_ARGUMENTS)
  
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8 				:: dT, dX,dY,dZ,val
  integer, dimension(3)			:: ext,global_ext
  real*8, dimension(1,3) :: pointcoords

  dT = CCTK_DELTA_TIME
  dX=CCTK_DELTA_SPACE(1)
  dY=CCTK_DELTA_SPACE(2)
  dZ=CCTK_DELTA_SPACE(3)
  ext = cctk_lsh
  global_ext = cctk_gsh

  if(bigp.ne.0.) then
     
     xbhold=xbh1
     ybhold=ybh1
     zbhold=zbh1
     
     pointcoords(1,1)=xbh1
     pointcoords(1,2)=ybh1
     pointcoords(1,3)=zbh1

     if(Symmetry.ne.4) then
        call interpolate_pointset(cctkGH,cctk_nghostzones,1,2,ext, &
             global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d, &
             Xglobal,Yglobal,Zglobal,pointcoords,shiftx_p,val)
        xbh1=xbh1-dT*val
        call interpolate_pointset(cctkGH,cctk_nghostzones,1,2,ext, &
             global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d, &
             Xglobal,Yglobal,Zglobal,pointcoords,shifty_p,val)
        ybh1=ybh1-dT*val
     endif

     if (Symmetry.ne.1) then
        call interpolate_pointset(cctkGH,cctk_nghostzones,1,2,ext, &
             global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d, &
             Xglobal,Yglobal,Zglobal,pointcoords,shiftz_p,val)
        zbh1=zbh1-dT*val     
     endif
     
  endif

end subroutine Center_Predictor
  
!===============================================================================
! Next, the Corrector step...
!===============================================================================
subroutine Center_Corrector(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8 				:: dT, dX,dY,dZ,val1,val2
  integer, dimension(3)			:: ext,global_ext
  real*8, dimension(1,3) :: pointcoords,pointcoords2

  dT = CCTK_DELTA_TIME
  dX=CCTK_DELTA_SPACE(1)
  dY=CCTK_DELTA_SPACE(2)
  dZ=CCTK_DELTA_SPACE(3)
  ext = cctk_lsh
  global_ext = cctk_gsh
  global_ext = cctk_gsh

  if(bigP.ne.0.) then
     
     if (iter_count == 1) then
        pointcoords(1,1)=xbhold
        pointcoords(1,2)=ybhold
        pointcoords(1,3)=zbhold
        pointcoords2(1,1)=xbh1
        pointcoords2(1,2)=ybh1
        pointcoords2(1,3)=zbh1
        if(Symmetry.ne.4) then
           call interpolate_pointset(cctkGH,cctk_nghostzones,1,2,ext, &
                global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d, &
                Xglobal,Yglobal,Zglobal,pointcoords,shiftx_p,val1)
           call interpolate_pointset(cctkGH,cctk_nghostzones,1,2,ext, &
                global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d, &
                Xglobal,Yglobal,Zglobal,pointcoords2,shiftx,val2)
           xbh1=xbhold-dT*0.5*(val1+val2)
           call interpolate_pointset(cctkGH,cctk_nghostzones,1,2,ext, &
                global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d, &
                Xglobal,Yglobal,Zglobal,pointcoords,shifty_p,val1)
           call interpolate_pointset(cctkGH,cctk_nghostzones,1,2,ext, &
                global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d, &
                Xglobal,Yglobal,Zglobal,pointcoords2,shifty,val2)
           ybh1=ybhold-dT*0.5*(val1+val2)
        endif
        if (Symmetry.ne.1) then
           call interpolate_pointset(cctkGH,cctk_nghostzones,1,2,ext, &
                global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d, &
                Xglobal,Yglobal,Zglobal,pointcoords,shiftz_p,val1)
           call interpolate_pointset(cctkGH,cctk_nghostzones,1,2,ext, &
                global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d, &
                Xglobal,Yglobal,Zglobal,pointcoords2,shiftz,val2)
           zbh1=zbhold-dT*0.5*(val1+val2)
        endif
!        write(6,*)'Correct_center1:',xbh1,ybh1,zbh1
        
     else ! iter_count=2 
        pointcoords(1,1)=xbhold
        pointcoords(1,2)=ybhold
        pointcoords(1,3)=zbhold
        pointcoords2(1,1)=xbh1
        pointcoords2(1,2)=ybh1
        pointcoords2(1,3)=zbh1
        if(Symmetry.ne.4) then
           call interpolate_pointset(cctkGH,cctk_nghostzones,1,2,ext, &
                global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d, &
                Xglobal,Yglobal,Zglobal,pointcoords,shiftx_p,val1)
           call interpolate_pointset(cctkGH,cctk_nghostzones,1,2,ext, &
                global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d, &
                Xglobal,Yglobal,Zglobal,pointcoords2,shiftx_t,val2)
           xbh1=xbhold-dT*0.5*(val1+val2)
           call interpolate_pointset(cctkGH,cctk_nghostzones,1,2,ext, &
                global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d, &
                Xglobal,Yglobal,Zglobal,pointcoords,shifty_p,val1)
           call interpolate_pointset(cctkGH,cctk_nghostzones,1,2,ext, &
                global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d, &
                Xglobal,Yglobal,Zglobal,pointcoords2,shifty_t,val2)
           ybh1=ybhold-dT*0.5*(val1+val2)
        endif
        if (Symmetry.ne.1) then
           call interpolate_pointset(cctkGH,cctk_nghostzones,1,2,ext, &
                global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d, &
                Xglobal,Yglobal,Zglobal,pointcoords,shiftz_p,val1)
           call interpolate_pointset(cctkGH,cctk_nghostzones,1,2,ext, &
                global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d, &
                Xglobal,Yglobal,Zglobal,pointcoords2,shiftz_t,val2)
           zbh1=zbhold-dT*0.5*(val1+val2)
        endif
!        write(6,*)'Correct_center2:',xbh1,ybh1,zbh1
        
     endif

  endif

 end subroutine Center_Corrector
