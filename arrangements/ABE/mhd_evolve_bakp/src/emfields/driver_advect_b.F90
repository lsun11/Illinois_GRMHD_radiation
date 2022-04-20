!-------------------------------------------------------
!    :: Driver routine for MHD timestepping, v2.0 ::
! (i.e., computing RHS's of all conservative variables)
!-------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine driver_advect_b(CCTK_ARGUMENTS)

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
        ! CAN'T HAVE SYNC'S HERE:
        !!call CCTK_SyncGroup(dummy,cctkGH,'em_evolve_ct::em_misc')
        !!if(excision_enable == 1) then
        !!   call vector_excision_bc(ext,X,Y,Z, &
        !!        ftxy,ftxz,ftyx,Symmetry,excision_zone_gf);
        !!   call vector_excision_bc(ext,X,Y,Z, &
        !!        ftyz,ftzx,ftzy,Symmetry,excision_zone_gf);
        !!end if

        ! Compute \partial_t B^i
        !call advect_b_ct(ext,X,Y,Z,Bxtilde_or_Ax_rhs,Bytilde_or_Ay_rhs,Bztilde_or_Az_rhs, &
        call advect_b_ct_cpp(cctkGH,ext,X,Y,Z,Bxtilde_or_Ax_rhs,Bytilde_or_Ay_rhs,Bztilde_or_Az_rhs, &
             ftxy,ftxz,ftyx,ftyz,ftzx,ftzy, &
             Symmetry)

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
        ! CAN'T HAVE SYNC'S HERE:
        !!call CCTK_SyncGroup(dummy,cctkGH,'em_evolve_ct::em_misc')
        !!if(excision_enable == 1) then
        !!   call vector_excision_bc(ext,X,Y,Z, &
        !!        ftxy,ftxz,ftyx,Symmetry,excision_zone_gf);
        !!   call vector_excision_bc(ext,X,Y,Z, &
        !!        ftyz,ftzx,ftzy,Symmetry,excision_zone_gf);
        !!end if

        ! Compute \partial_t B^i
        call advect_b_ct_cpp(cctkGH,ext,X,Y,Z,Bxtilde_or_Ax_rhs,Bytilde_or_Ay_rhs,Bztilde_or_Az_rhs, &
             ftxy,ftxz,ftyx,ftyz,ftzx,ftzy, &
             Symmetry)
        !        call advect_b_ct(ext,X,Y,Z,Bxtilde_or_Ax_rhs,Bytilde_or_Ay_rhs,Bztilde_or_Az_rhs, &
        !             ftxy,ftxz,ftyx,ftyz,ftzx,ftzy, &
        !             Symmetry)
     end if

     ! Hyperbolic divergence cleaning: adding an extra MHD mode that satisfies
     !  the telegraph equation, propagating & dissipating away the divergence
     !  of the magnetic fields.
     if(hyperbolic_divergence_cleaning_enable==1) then
        call hyperbolic_divergence_cleaning(CCTK_ARGUMENTS)
     else
        Blagrangemultiplier_p_p = 0.D0
        Blagrangemultiplier_p = 0.D0
        Blagrangemultiplier = 0.D0
        Blagrangemultiplier_rhs = 0.D0
     end if

     ! FIXME: NEED THIS?!
     !  call CCTK_SyncGroup(dummy,cctkGH,'em_evolve_ct::em_rhs')
     !  call CartSymGN(dummy,cctkGH,'em_evolve_ct::em_rhs')

     !!if(excision_enable == 1) then
     !!   call vector_excision_bc(ext,X,Y,Z, &
     !!        Bxtilde_or_Ax_rhs,Bytilde_or_Ay_rhs,Bztilde_or_Az_rhs, &
     !!        Symmetry,excision_zone_gf);
     !!end if

     if(excision_enable == 1) then
        call remove_interior(ext,X,Y,Z,Bxtilde_or_Ax_rhs,excision_zone_gf,Symmetry)
        call remove_interior(ext,X,Y,Z,Bytilde_or_Ay_rhs,excision_zone_gf,Symmetry)
        call remove_interior(ext,X,Y,Z,Bztilde_or_Az_rhs,excision_zone_gf,Symmetry)
        if (hyperbolic_divergence_cleaning_enable==1) &
          call remove_interior(ext,X,Y,Z,Blagrangemultiplier_rhs,excision_zone_gf,Symmetry)
     end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !          \||                              ||/          !
     !           ::             END              ::           !
     !           ::    CONSTRAINED TRANSPORT,    ::           !
     !           ::COMPUTE \partial_t \tilde{B}^i::           !
     !          /||                              ||\          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !First set the fij's and ftildeij's to zero in the ghostzones.  If we don't, then
     !  when Bitilde_rhs is set, we get undefined values.
     if (constrained_transport_scheme .ne. 3) then 
        !$omp parallel do
        do k=1,cctk_lsh(3)
           do j=1,cctk_lsh(2)
              do i=1,cctk_lsh(1)
                 if( &
                      i.le.cctk_nghostzones(1) .or. &
                      j.le.cctk_nghostzones(2) .or. &
                      k.le.cctk_nghostzones(3) .or. &
                      i.gt.cctk_lsh(1)-cctk_nghostzones(1) .or. &
                      j.gt.cctk_lsh(2)-cctk_nghostzones(2) .or. &
                      k.gt.cctk_lsh(3)-cctk_nghostzones(3)) then
                    Bxtilde_or_Ax_rhs(i,j,k) = 0.D0
                    Bytilde_or_Ay_rhs(i,j,k) = 0.D0
                    Bztilde_or_Az_rhs(i,j,k) = 0.D0
                 end if
              end do
           end do
        end do
        !$omp end parallel do
     end if


  end if

end subroutine driver_advect_b
