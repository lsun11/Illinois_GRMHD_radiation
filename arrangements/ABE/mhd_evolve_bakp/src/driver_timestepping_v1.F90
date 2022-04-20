!-------------------------------------------------------
!    :: Driver routine for MHD timestepping, v1.0 ::
! (i.e., computing RHS's of all conservative variables)
!-------------------------------------------------------

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
subroutine mhd_timestepping_v1(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: ext
  real*8                :: dX,dY,dZ
  integer               :: index,ierr,handle,dummy
  CCTK_REAL             :: reduction_value
  integer               :: AXISYM
  parameter(AXISYM = 4)

  ext = cctk_lsh

  dX = CCTK_DELTA_SPACE(1)
  dY = CCTK_DELTA_SPACE(2)
  dZ = CCTK_DELTA_SPACE(3)

  !The following lines are only needed if we have emfields turned off, since emfields does exactly the same thing at iter_count==1!
  if(iter_count==1.and.em_evolve_enable==0) then
     ! Note that h_p == [h from previous timeSTEP, not previous MoL step].  h_p MUST be set for primitive solver.
     h_p = h

     ! In the MHD code, rho_b and P are NOT updated on the boundaries in Axisymmetry.  
     !   As a result, the maximum value of these quantities may occur at the boundary and not change!
     !   Here we correct for this problem by zeroing out the boundary values of rho_b and P.
     if(Symmetry==AXISYM) then
        if(Z(1,1,1).lt.0.D0) then
           rho_b(:,:,1) = rho_b(:,:,2)
           P(:,:,1) = P(:,:,2)
        end if
        if(X(1,1,1).lt.0.D0) then
           rho_b(1,:,:) = rho_b(2,:,:)
           P(1,:,:) = P(2,:,:)
        end if
     end if

     ! Also, rho_b_max, P_max, and rhos_max must also be set for primitive solver!
     call CCTK_VarIndex(index,"mhd_evolve::rho_b")
     call CCTK_ReductionHandle(handle,"maximum")
     if (handle .gt. 0) then
        call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,reduction_value,1,index)
        if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH).eq.0) then
           !        print *,"1: Maximum value of rho_b is ",reduction_value
        end if
     else
        call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
     end if
     rho_b_max = reduction_value  
     call CCTK_VarIndex(index,"mhd_evolve::P")
     call CCTK_ReductionHandle(handle,"maximum")
     if (handle .gt. 0) then
        call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,reduction_value,1,index)
        if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH).eq.0) then
           !        print *,"1: Maximum value of P is ",reduction_value
        end if
     else
        call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
     end if
     P_max = reduction_value  
     call CCTK_VarIndex(index,"mhd_evolve::rho_star")
     call CCTK_ReductionHandle(handle,"maximum")
     if (handle .gt. 0) then
        call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,reduction_value,1,index)
        if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH).eq.0) then
           !        print *,"1: Maximum value of rho_star is ",reduction_value
        end if
     else
        call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
     end if
     rhos_max = reduction_value  
  end if

  ! initialize right hand sides to zero <-- Do we really need to zero these?
  rho_star_rhs = 0.D0
  tau_rhs = 0.D0
  mhd_st_x_rhs = 0.D0
  mhd_st_y_rhs = 0.D0
  mhd_st_z_rhs = 0.D0


  m = 1
  call reconstruction_v1(CCTK_PASS_FTOF)
  call mhd_flux_and_rhs1D_v1(CCTK_PASS_FTOF)

  if (Symmetry .ne. AXISYM) then
     m = 2
     call reconstruction_v1(CCTK_PASS_FTOF)
     call mhd_flux_and_rhs1D_v1(CCTK_PASS_FTOF)
  end if

  m = 3
  call reconstruction_v1(CCTK_PASS_FTOF)
  call mhd_flux_and_rhs1D_v1(CCTK_PASS_FTOF)

  if(enable_HARM_energyvariable==0) then
     if(1==0) then
     call mhd_tau_curvature_cpp(cctkGH,cctk_lsh,cctk_nghostzones,Symmetry, &
          tau_rhs,rho_star,h,P, &
          sbt,sbx,sby,sbz, &
          u0,vx,vy,vz, &
          lapm1,shiftx,shifty,shiftz, &
          phi,trK, &
          gxx,gxy,gxz,gyy,gyz,gzz, &
          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
          Axx,Axy,Axz,Ayy,Ayz,Azz, &
          dX,dY,dZ)
  end if

        call mhd_tau_curvature(ext,tau_rhs,rho_star,h,P, &
             sbt,sbx,sby,sbz, &
             u0,vx,vy,vz,lapm1, &
             shiftx,shifty,shiftz, phi, &
             gxx,gxy,gxz,gyy,gyz,gzz, &
             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
             trK,Axx,Axy,Axz,Ayy,Ayz,Azz, &
             dX,dY,dZ)
  end if


  if(excision_enable==1) then
     call remove_interior(ext,X,Y,Z,rho_star_rhs,excision_zone_gf,Symmetry)
     call remove_interior(ext,X,Y,Z,tau_rhs,excision_zone_gf,Symmetry)
     call remove_interior(ext,X,Y,Z,mhd_st_x_rhs,excision_zone_gf,Symmetry)
     call remove_interior(ext,X,Y,Z,mhd_st_y_rhs,excision_zone_gf,Symmetry)
     call remove_interior(ext,X,Y,Z,mhd_st_z_rhs,excision_zone_gf,Symmetry)
  end if

end subroutine mhd_timestepping_v1
