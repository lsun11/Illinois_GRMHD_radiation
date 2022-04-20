!-------------------------------------------------------
!    :: Driver routine for MHD timestepping, v2.0 ::
! (i.e., computing RHS's of all conservative variables)
!-------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine driver_constrained_transport(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  real*8                :: dX,dY,dZ,b2bt
  integer               :: index,ierr,handle,dummy
  CCTK_REAL             :: reduction_value
  integer               :: AXISYM,i,j,k
  parameter(AXISYM = 4)

  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  if(em_evolve_enable==1) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !          \||                              ||/          !
     !           ::    CONSTRAINED TRANSPORT,    ::           !
     !           ::COMPUTE \partial_t \tilde{B}^i::           !
     !          /||                              ||\          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if(constrained_transport_scheme == 1) then
        ! Constrained transport: Compute \tilde{f}^{ij}
        !           call flux_ct(ext,X,Y,Z,fxy,fxz,fyx,fyz,fzx,fzy, &
        !                ftxy,ftxz,ftyx,ftyz,ftzx,ftzy, &
        !                Sym_Bz,Symmetry)  
        call flux_ct_cpp(cctkGH,ext,X,Y,Z,fxy,fxz,fyx,fyz,fzx,fzy, &
             ftxy,ftxz,ftyx,ftyz,ftzx,ftzy, &
             Sym_Bz,Symmetry)  
     else if(constrained_transport_scheme == 2) then
        write(*,*) "SORRY, FLUX_CD SCHEME HAS NOT YET BEEN PORTED TO CARPET VERSION OF ABE!"
        stop
        !        ftyz = 0.D0
        !        ftzx = 0.D0
        !        ftzy = 0.D0
        !        call flux_new_cd(ext,X,Y,Z,ftxy,ftxz,ftyx,fxy,fxz,fyx,fyz,fzx,fzy,Sym_Bz,Symmetry)
        !        ! CAN'T HAVE SYNC'S HERE:
        !        !!call CCTK_SyncGroup(dummy,cctkGH,'em_evolve_ct::em_misc')
        !        if(excision_enable == 1) then
        !           call vector_excision_bc(ext,X,Y,Z, &
        !                ftxy,ftxz,ftyx,Symmetry,excision_zone_gf);
        !           call vector_excision_bc(ext,X,Y,Z, &
        !                ftyz,ftzx,ftzy,Symmetry,excision_zone_gf);
        !        end if
        !
        !        call advect_flux_new_cd(ext,dX,dY,dZ,X,Y,Z,Bxtilde_or_Ax_rhs,Bytilde_or_Ay_rhs,Bztilde_or_Az_rhs,ftxy,ftxz,ftyx,fxy,fxz,fyx,fyz,fzx,fzy,Sym_Bz,Symmetry)
        !
     else if(constrained_transport_scheme == 0) then
        ftxy = fxy
        ftxz = fxz
        ftyx = fyx
        ftyz = fyz
        ftzx = fzx
        ftzy = fzy
     end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !          \||                              ||/          !
     !           ::             END              ::           !
     !           ::    CONSTRAINED TRANSPORT,    ::           !
     !           ::COMPUTE \partial_t \tilde{B}^i::           !
     !          /||                              ||\          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end if

end subroutine driver_constrained_transport
