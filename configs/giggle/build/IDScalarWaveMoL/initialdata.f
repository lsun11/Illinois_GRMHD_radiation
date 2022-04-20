      subroutine IDScalarWaveMoL_InitialData (cctk_dim,cctk_gsh,cctk_lsh
     &,cctk_lbnd,cctk_ubnd,cctk_lssh,cctk_from,cctk_to,cctk_bbox,cctk_de
     &lta_time,cctk_time,cctk_delta_space,cctk_origin_space,cctk_levfac,
     &cctk_levoff,cctk_levoffdenom,cctk_timefac,cctk_convlevel,cctk_conv
     &fac,cctk_nghostzones,cctk_iteration,cctkGH, X0scalarevolveerror,X1
     &scalarevolveerror,X2scalarevolveerror,phierror,psierror, X0coordin
     &ates,X0scalarevolve,X1coordinates,X1scalarevolve,X2coordinates,X2s
     &calarevolve,coarse_dx,coarse_dy,coarse_dz,phi,phi_p,phi_p_p,psi,ps
     &i_p,psi_p_p,r,x,y,z)
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
      
      REAL*8 one, two
      parameter (one = 1)
      parameter (two = 2)
      REAL*8 pi
      parameter (pi = 3.141592653589793238462643383279502884197169399375
     &105820974944592307816406286208998628034825342117068d0)
      REAL*8 omega
      REAL*8 tmp
      integer i, j, k
      if ((CCTK_Equals(initial_data, "plane wave") .ne. 0)) then
         omega = sqrt(wave_number(1)**2 + wave_number(2)**2
     $        + wave_number(3)**2)
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  tmp =  wave_number(1)*(x(i,j,k)-phase_offset(1))
     $                 + wave_number(2)*(y(i,j,k)-phase_offset(2))
     $                 + wave_number(3)*(z(i,j,k)-phase_offset(3))
     $                 + omega*(cctk_time-time_offset)
                  phi(i,j,k) = amplitude * cos (2*pi * tmp)
                  psi(i,j,k) = - amplitude * 2*pi * omega * sin (2*pi * 
     &tmp)
               end do
            end do
         end do
      else if ((CCTK_Equals(initial_data, "Gaussian pulse") .ne. 0)) the
     &n
         omega = sqrt(pulse_direction(1)**2 + pulse_direction(2)**2
     $        + pulse_direction(3)**2)
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  tmp =  pulse_direction(1)*(x(i,j,k)-pulse_offset(1))
     $                 + pulse_direction(2)*(y(i,j,k)-pulse_offset(2))
     $                 + pulse_direction(3)*(z(i,j,k)-pulse_offset(3))
     $                 + omega*(cctk_time-time_offset)
                  phi(i,j,k) = amplitude * exp (-0.5d0 * tmp**2)
                  psi(i,j,k) = - tmp * omega * phi(i,j,k)
               end do
            end do
         end do
      else if ((CCTK_Equals(initial_data, "Gaussian") .ne. 0)) then
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  tmp = sqrt(  (x(i,j,k)-origin(1))**2
     $                       + (y(i,j,k)-origin(2))**2
     $                       + (z(i,j,k)-origin(3))**2)
                  if (tmp .gt. 1.0d-6 * sigma) then
                     phi(i,j,k) = amplitude/2
     $                    * (  (cctk_time + tmp)/tmp * exp (-0.5d0 * ((c
     &ctk_time + tmp - radius) / sigma)**2)
     $                       - (cctk_time - tmp)/tmp * exp (-0.5d0 * ((c
     &ctk_time - tmp - radius) / sigma)**2))
                     psi(i,j,k) = amplitude/2
     $                    * (+ (1 + (cctk_time + tmp) * (cctk_time + tmp
     & - radius) / sigma**2) / tmp * exp (-0.5d0 * ((cctk_time + tmp - r
     &adius) / sigma)**2)
     $                       - (1 + (cctk_time - tmp) * (cctk_time - tmp
     & - radius) / sigma**2) / tmp * exp (-0.5d0 * ((cctk_time - tmp - r
     &adius) / sigma)**2))
                  else
                     phi(i,j,k) = amplitude
     $                    * (  (radius * cctk_time - cctk_time**2 + sigm
     &a**2) / sigma**2
     $                       + (radius**3 * cctk_time - cctk_time**4 + 6
     & * cctk_time**2 * sigma**2 - 3 * sigma**4 - 3 * radius**2 * (cctk_
     &time - sigma)**2 + 3 * radius * (cctk_time**3 - 3 * cctk_time * si
     &gma**2)) / (6 * sigma**6) * tmp**2)
     $                    * exp (-0.5d0 * (radius - cctk_time)**2 / sigm
     &a**2)
                     psi(i,j,k) = amplitude
     $                    * (  (radius**2 * cctk_time - 2 * radius * cct
     &k_time**2 + cctk_time**3 + 2 * radius * sigma**2 - 3 * cctk_time *
     & sigma**2) / sigma**4
     $                       + (radius**4 * cctk_time + cctk_time**5 - 1
     &0 * cctk_time**3 * sigma**2 + 15 * cctk_time * sigma**4 - 4 * radi
     &us**4 * (cctk_time**2 - sigma**2) + 6 * radius**2 * (cctk_time**3 
     &- 3 * cctk_time * sigma**2) - 4 * radius * (cctk_time**4 - 6 * cct
     &k_time**2 * sigma**2 + 3 * sigma**4)) / (6 * sigma**8) * tmp**2)
     $                    * exp (-0.5d0 * (radius - cctk_time)**2 / sigm
     &a**2)
                  end if
               end do
            end do
         end do
      else if ((CCTK_Equals(initial_data, "level index") .ne. 0)) then
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  phi(i,j,k) = log (one * cctk_levfac(1)) / log (two)
                  psi(i,j,k) = log (one * cctk_levfac(1)) / log (two)
               end do
            end do
         end do
      end if
      end
