      subroutine IDScalarWaveMoL_Errors (cctk_dim,cctk_gsh,cctk_lsh,cctk
     &_lbnd,cctk_ubnd,cctk_lssh,cctk_from,cctk_to,cctk_bbox,cctk_delta_t
     &ime,cctk_time,cctk_delta_space,cctk_origin_space,cctk_levfac,cctk_
     &levoff,cctk_levoffdenom,cctk_timefac,cctk_convlevel,cctk_convfac,c
     &ctk_nghostzones,cctk_iteration,cctkGH, X0scalarevolveerror,X1scala
     &revolveerror,X2scalarevolveerror,phierror,psierror, X0coordinates,
     &X0scalarevolve,X1coordinates,X1scalarevolve,X2coordinates,X2scalar
     &evolve,coarse_dx,coarse_dy,coarse_dz,phi,phi_p,phi_p_p,psi,psi_p,p
     &si_p_p,r,x,y,z)
      implicit none
      INTEGER cctk_dim
      INTEGER cctk_gsh(cctk_dim),cctk_lsh(cctk_dim)
      INTEGER cctk_lbnd(cctk_dim),cctk_ubnd(cctk_dim)
      INTEGER cctk_lssh(3*cctk_dim)
      INTEGER cctk_from(cctk_dim),cctk_to(cctk_dim)
      INTEGER cctk_bbox(2*cctk_dim)
      REAL*8 cctk_delta_time, cctk_time
      REAL*8 cctk_delta_space(cctk_dim)
      REAL*8 cctk_origin_space(cctk_dim)
      INTEGER cctk_levfac(cctk_dim)
      INTEGER cctk_levoff(cctk_dim)
      INTEGER cctk_levoffdenom(cctk_dim)
      INTEGER cctk_timefac
      INTEGER cctk_convlevel
      INTEGER cctk_convfac
      INTEGER cctk_nghostzones(cctk_dim)
      INTEGER cctk_iteration
      integer*8 cctkGH
      INTEGER X0scalarevolveerror
      INTEGER X1scalarevolveerror
      INTEGER X2scalarevolveerror
      REAL*8 phierror(X0scalarevolveerror,X1scalarevolveerror,X2scalarev
     &olveerror)
      REAL*8 psierror(X0scalarevolveerror,X1scalarevolveerror,X2scalarev
     &olveerror)
      INTEGER X0coordinates
      INTEGER X0scalarevolve
      INTEGER X1coordinates
      INTEGER X1scalarevolve
      INTEGER X2coordinates
      INTEGER X2scalarevolve
      REAL*8 coarse_dx
      REAL*8 coarse_dy
      REAL*8 coarse_dz
      REAL*8 phi(X0scalarevolve,X1scalarevolve,X2scalarevolve)
      REAL*8 phi_p(X0scalarevolve,X1scalarevolve,X2scalarevolve)
      REAL*8 phi_p_p(X0scalarevolve,X1scalarevolve,X2scalarevolve)
      REAL*8 psi(X0scalarevolve,X1scalarevolve,X2scalarevolve)
      REAL*8 psi_p(X0scalarevolve,X1scalarevolve,X2scalarevolve)
      REAL*8 psi_p_p(X0scalarevolve,X1scalarevolve,X2scalarevolve)
      REAL*8 r(X0coordinates,X1coordinates,X2coordinates)
      REAL*8 x(X0coordinates,X1coordinates,X2coordinates)
      REAL*8 y(X0coordinates,X1coordinates,X2coordinates)
      REAL*8 z(X0coordinates,X1coordinates,X2coordinates)
      
      integer      CCTK_Equals, CCTK_MyProc, CCTK_nProcs, CCTK_IsThornAc
     &tive
      external     CCTK_Equals, CCTK_MyProc, CCTK_nProcs, CCTK_IsThornAc
     &tive
      integer*8 CCTK_PointerTo, CCTK_NullPointer
      external     CCTK_PointerTo, CCTK_NullPointer
      
      REAL*8  excision_radius
      REAL*8  run_time
      INTEGER*4 Symmetry
      INTEGER*4 bssn_enable
      INTEGER*4 cowling_enable
      INTEGER*4 excision_enable
      INTEGER*4 fisheye_enable
      INTEGER*4 iter_count
      INTEGER*4 number_of_mol_ministeps
      INTEGER*4 rot_metric
      INTEGER*4 trA_detg_enforce
      COMMON /cctk_params_global/excision_radius,run_time,Symmetry,bssn_
     &enable,cowling_enable,excision_enable,fisheye_enable,iter_count,nu
     &mber_of_mol_ministeps,rot_metric,trA_detg_enforce
      REAL*8  amplitude
      REAL*8  origin(3)
      REAL*8  phase_offset(3)
      REAL*8  pulse_direction(3)
      REAL*8  pulse_offset(3)
      REAL*8  radius
      REAL*8  sigma
      REAL*8  time_offset
      REAL*8  wave_number(3)
      integer*8  initial_data
      COMMON /IDScalarWaveMoLpriv/amplitude,origin,phase_offset,pulse_di
     &rection,pulse_offset,radius,sigma,time_offset,wave_number,initial_
     &data
      
      REAL*8 tmp
      integer i, j, k
      do k=1,cctk_lsh(3)
         do j=1,cctk_lsh(2)
            do i=1,cctk_lsh(1)
               phierror(i,j,k) = phi(i,j,k)
               psierror(i,j,k) = psi(i,j,k)
            end do
         end do
      end do
      call IDScalarWaveMol_InitialData (cctk_dim,cctk_gsh,cctk_lsh,cctk_
     &lbnd,cctk_ubnd,cctk_lssh,cctk_from,cctk_to,cctk_bbox,cctk_delta_ti
     &me,cctk_time,cctk_delta_space,cctk_origin_space,cctk_levfac,cctk_l
     &evoff,cctk_levoffdenom,cctk_timefac,cctk_convlevel,cctk_convfac,cc
     &tk_nghostzones,cctk_iteration,cctkGH, X0scalarevolveerror,X1scalar
     &evolveerror,X2scalarevolveerror,phierror,psierror, X0coordinates,X
     &0scalarevolve,X1coordinates,X1scalarevolve,X2coordinates,X2scalare
     &volve,coarse_dx,coarse_dy,coarse_dz,phi,phi_p,phi_p_p,psi,psi_p,ps
     &i_p_p,r,x,y,z)
      do k=1,cctk_lsh(3)
         do j=1,cctk_lsh(2)
            do i=1,cctk_lsh(1)
               tmp = phi(i,j,k)
               phi(i,j,k) = phierror(i,j,k)
               phierror(i,j,k) = phi(i,j,k) - tmp
               tmp = psi(i,j,k)
               psi(i,j,k) = psierror(i,j,k)
               psierror(i,j,k) = psi(i,j,k) - tmp
            end do
         end do
      end do
      end
