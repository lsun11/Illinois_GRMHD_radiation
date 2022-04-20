#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
! Compute shock error integrands
!-----------------------------------------------------------------------------
subroutine shock_error_integrand(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: num_rows
  parameter(num_rows = 50)
  real*8, dimension(num_rows)              :: R_out,rho_out,P_out,vx_out

  real*8                                   :: old_surf_radius,SymmFactor,gammashock,xmin,xmax,fudgefactor
  real*8                                   :: rhoshock_out,Pshock_out,vxshock_out
  integer :: header_flag,handle,index,ierr,ii,j,k

  Shock_Which_Int = 1000

  if(MOD(cctk_iteration,out_every)==0) then
     !     call CCTK_VarIndex(index,"grid::X")
     !     call CCTK_ReductionHandle(handle,"minimum")
     !     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xmin,1,index)
     !     call CCTK_ReductionHandle(handle,"maximum")
     !     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,xmax,1,index)

     gammashock = 1.D0 + 1.D0/npolyshock

     xmin = -1
     xmax = 1

     write(*,*) "HIIIII>>>",cctk_time,xmin,xmax,gammashock,rho_l,rho_r
     call shocktests_rieman_analytic(cctk_time,xmin,xmax,gammashock,rho_l,rho_r,num_rows,R_out,rho_out,P_out,vx_out)

     rho_shock_analytic = 0.D0
     rho_error = 0.D0
     rho_error_reflection = 0.D0

     P_shock_analytic = 0.D0
     P_error = 0.D0
     P_error_reflection = 0.D0

     vx_shock_analytic = 0.D0
     vx_error = 0.D0
     vx_error_reflection = 0.D0

     ! We must multiply by cctk_levfac's here because the Carpet weight function assumes 3D integrals.
     !    You see, when a CCTK_Reduce() sum is called, it actually sums the integrand*weight function.
     !    The weight function is set to 1/(cctk_levfac(1)*cctk_levfac(2)*cctk_levfac(3)) in non-
     !    ghostzone regions of the grid, where e.g., cctk_levfac(1)=1 on the coarsest level and
     !    it increases from there by factors of 2 as the grids get finer.
     !    Thus for a 1D integral along the x direction, we need to multiply by cctk_levfac(2)*cctk_levfac(3).

     ! In addition, we must multiply by 4.D0 on all levels because we set up the grid as follows:
     !
     !CoordBase::boundary_size_z_lower = 3
     !CoordBase::boundary_shiftout_z_lower = 1
     !CoordBase::boundary_size_y_lower = 3
     !CoordBase::boundary_shiftout_y_lower = 1
     !
     !    Any time there is a boundary_shiftout_z_lower, Carpet assumes that this a symmetry boundary, and will 
     !    multiply the weight function by a factor of 0.5 for each "symmetry".  Here, Carpet assumes that there
     !    are 2 symmetries here: one in y and one in z.  Thus the weight function on the coarsest level will be 
     !    set to 0.25, and to 0.25*0.25 on the finer level (due to the additional cctk_levfac factors of 2.  Note
     !    that these symmetry factors work well in the case of, e.g., equatorial symmetry, where the z=0 plane is
     !    included on the grid, since we multiply the final result by symmfactor=2, and we don't want to double-count
     !    the z=0 points.
     fudgefactor = 4.D0*cctk_levfac(2)*cctk_levfac(3)

     j=cctk_nghostzones(2)+1
     k=cctk_nghostzones(3)+1
     do ii=1,cctk_lsh(1)
        call get_analytic_interp(X(ii,1,1),num_rows,R_out,rho_out,rhoshock_out)
        rho_shock_analytic(ii,j,k) = rhoshock_out
        rho_error(ii,j,k) = abs(rho_b(ii,j,k) - rhoshock_out)*fudgefactor

        call get_analytic_interp(X(ii,1,1),num_rows,R_out,P_out,Pshock_out)
        P_shock_analytic(ii,j,k) = Pshock_out
        P_error(ii,j,k) = abs(P(ii,j,k) - Pshock_out)*fudgefactor

        call get_analytic_interp(X(ii,1,1),num_rows,R_out,vx_out,vxshock_out)
        vx_shock_analytic(ii,j,k) = vxshock_out
        vx_error(ii,j,k) = abs(vx(ii,j,k) - vxshock_out)*fudgefactor

        ! following if statement singles out reflection error!
        if(X(ii,1,1).lt.0.2D0 .and. X(ii,1,1) .gt. -0.2D0) then
           call get_analytic_interp(X(ii,1,1),num_rows,R_out,rho_out,rhoshock_out)
           rho_error_reflection(ii,j,k) = abs(rho_b(ii,j,k) - rhoshock_out)*fudgefactor

           call get_analytic_interp(X(ii,1,1),num_rows,R_out,P_out,Pshock_out)
           P_error_reflection(ii,j,k) = abs(P(ii,j,k) - Pshock_out)*fudgefactor

           call get_analytic_interp(X(ii,1,1),num_rows,R_out,vx_out,vxshock_out)
           vx_error_reflection(ii,j,k) = abs(vx(ii,j,k) - vxshock_out)*fudgefactor
        end if

     end do
  end if
end subroutine shock_error_integrand

subroutine get_analytic_interp(x,num_rows,R_out,rho_out,rho_shock_analytic)
  implicit none

  real*8 :: x,rho_shock_analytic
  integer :: num_rows
  real*8,dimension(num_rows) :: R_out,rho_out

  real*8 :: y1,y2
  real*8 :: x1,x2
  real*8 :: a,b
  integer i,j

  j=0

  do i=2,num_rows
     if(j.eq.0) then
        if(x.ge.R_out(i-1) .and. x.le.R_out(i)) then
           j=i
        end if
     end if
  end do

  if(j.eq.0) then
     write(*,*) "COULD NOT FIND VALUE FOR X INTERPOLATION!, x=",x
     stop
  end if

  x1 = R_out(j-1)
  x2 = R_out(j)

  y1 = rho_out(j-1)
  y2 = rho_out(j)

  a = (y2-y1) / (x2-x1)
  b = y1 - a*x1

  rho_shock_analytic = a*x + b

end subroutine get_analytic_interp
