#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
! Particle tracer for BHNS thorn.
! This routine has been optimized, and should not slow down
! your run with 1000 particles.
!
! Read the code & comments carefully to figure out the timestepping.
!-----------------------------------------------------------------------------
  
subroutine bhns_diagnostics_trace_particles(CCTK_ARGUMENTS)
  
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
 
  ! Local variables needed for file output:
  character                                  :: filename*50
  real*8, dimension(N_particles_to_trace,5) :: export_data
  character, dimension(50)                   :: data_headers*20
  integer                                    :: num_cols,header_flag

  ! Local variables needed for RK4 timestepping:
  real*8                                   :: dT
  real*8                                   :: k1_x,k1_y,k1_z
  real*8                                   :: k2_x,k2_y,k2_z
  real*8                                   :: k3_x,k3_y,k3_z
  real*8                                   :: k4_x,k4_y,k4_z

  ! Local variables needed for interpolation:
  real*8,dimension(N_particles_to_trace,3) :: pointcoords
  integer                                  :: vindex

  ! Loop dummy variable
  integer                                  :: which_particle

  !------------------------------------------------------------!
  ! MRI wavelength calculator
  ! Loop dummy variables
  integer                                  :: i,j,which_gridpoint
  real*8,dimension(MRI_wavelength_calculator_Nxy*MRI_wavelength_calculator_Nxy,3) :: MRI_wavelength_grid
  real*8,dimension(MRI_wavelength_calculator_Nxy*MRI_wavelength_calculator_Nxy,8) :: MRI_export_data
  !------------------------------------------------------------!
  ! Max B poloidal, B toroidal diagnostic, as well as rho_star weighted b^2/P diagnostic.
  real*8, dimension(50) :: export_data_singleline
  real*8                :: numerator,denominator
  integer :: handle,ierr
  !------------------------------------------------------------!

  header_flag=0

  if(cctk_iteration.ge.bhns_particle_tracer_start) then

     if(Symmetry.ne.0 .and. Symmetry.ne.1) then
        write(*,*) "Sorry, particle tracers only work when Symmetry==0 or 1."
        stop
     end if

     !-----------------------------------------------------------------------------
     ! INITIAL DATA

     if(cctk_iteration.eq.bhns_particle_tracer_start) then
        ! Don't change the default random seed!
        !CALL RANDOM_SEED                    ! Initialize random number generator

        do which_particle=1,N_particles_to_trace
           ! Seed particles inside sphere centered on coords of NS, with radius r = bhns_tracer_r

           ! First choose a guess that is very far from the star, so that we go inside the do while loop below.
           tracer_x(which_particle) = bhns_tracer_r + 1000.D0
           tracer_y(which_particle) = 0.D0
           tracer_z(which_particle) = 0.D0

           ! Next we throw particles into a cube, and if they fall within the sphere of the NS, we keep them:
           do while(sqrt(tracer_x(which_particle)**2 + tracer_y(which_particle)**2 + tracer_z(which_particle)**2) > bhns_tracer_r)

              CALL RANDOM_NUMBER(tracer_x(which_particle))
              CALL RANDOM_NUMBER(tracer_y(which_particle))
              CALL RANDOM_NUMBER(tracer_z(which_particle))

              ! tracer_xyz(which_particle) is now between 0 and 1. We want it between -1 and 1:
              tracer_x(which_particle) = 2.D0*(tracer_x(which_particle) - 0.5D0)
              tracer_y(which_particle) = 2.D0*(tracer_y(which_particle) - 0.5D0)
              tracer_z(which_particle) = 2.D0*(tracer_z(which_particle) - 0.5D0)

              ! Multiply tracer_xyz(which_particle) by bhns_tracer_r, which will give us a point within a cube of 
              !    half sidelength bhns_tracer_r, centered at the origin.
              tracer_x(which_particle) = tracer_x(which_particle)*bhns_tracer_r
              tracer_y(which_particle) = tracer_y(which_particle)*bhns_tracer_r
              tracer_z(which_particle) = tracer_z(which_particle)*bhns_tracer_r

           end do

           ! Shift the center of the ball of particles to the NS center of mass
           tracer_x(which_particle) = tracer_x(which_particle) + position_x(2)
           tracer_y(which_particle) = tracer_y(which_particle) + position_y(2)
           tracer_z(which_particle) = tracer_z(which_particle) + position_z(2)

           ! If we have a particle at z<0 with equatorial symmetry turned on, set it to z>0:
           if(Symmetry.eq.1 .and. tracer_z(which_particle).lt.0.D0) then
              tracer_z(which_particle) = abs(tracer_z(which_particle))
           end if

        end do

        ! We are now at t = t_n, so we go ahead and perform the first RK4 substep.

        !Warning: Don't adjust the value of bhns_rk4_particle_tracer_step after
        ! evaluating the first RK4 step below. Setting it to one forces us to dump
        ! the initial particle distribution and go to the next RK4 step
        ! properly.
        bhns_rk4_particle_tracer_step=1
        
        do which_particle=1,N_particles_to_trace
           pointcoords(which_particle,1) = tracer_x(which_particle)
           pointcoords(which_particle,2) = tracer_y(which_particle)
           pointcoords(which_particle,3) = tracer_z(which_particle)
        end do

!!$           x_rhs(:,bhns_rk4_particle_tracer_step) = 10.D0*sin(CCTK_TIME)
!!$           y_rhs(:,bhns_rk4_particle_tracer_step) = 10.D0*sin(CCTK_TIME)
!!$           z_rhs(:,bhns_rk4_particle_tracer_step) = 10.D0*sin(CCTK_TIME)
!!$
        call CCTK_VarIndex(vindex,"mhd_evolve::vx")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,x_rhs(:,bhns_rk4_particle_tracer_step))
        call CCTK_VarIndex(vindex,"mhd_evolve::vy")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,y_rhs(:,bhns_rk4_particle_tracer_step))
        call CCTK_VarIndex(vindex,"mhd_evolve::vz")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,z_rhs(:,bhns_rk4_particle_tracer_step))


     end if ! End of initial data routine.
     !-----------------------------------------------------------------------------

     !-----------------------------------------------------------------------------
     ! RK4 TIMESTEPPING

     if(MOD(cctk_iteration,particle_tracer_substep_every)==0) then
        ! RK4 Timestepping procedure to advance each particle to the next timestep.
        ! For notes on the RK4 algorithm used here, visit:
        !   http://en.wikipedia.org/w/index.php?title=Runge%E2%80%93Kutta_methods&oldid=471711997
        !
        ! We are solving y' = f(t,y), where
        !           y' = dx^i/dt , i=1,2,3
        !           f(t,y) = v^i , i=1,2,3
        ! Thus, the basic ODE we're solving is:
        !           dx^i/dt = v^i.
        !
        ! To evaluate the RHS of the above ODE, we need to interpolate v^i at 
        ! values x^i at times t, t+dt/2, and t+dt -- all as specified by RK4.
        ! 
        ! Note that these interpolations are expensive, so dt is some multiple
        ! of the time between _full_ simulation timesteps, dT.
        ! We should NOT need dt<dT to get reasonable particle positions!

        dT = CCTK_DELTA_TIME*particle_tracer_substep_every*2

        write(*,*) "****************"
        write(*,*) " Particle Tracer substep # ",bhns_rk4_particle_tracer_step,"dt=",dT
        write(*,*) " Particle Tracer thinks the time is =",CCTK_TIME

        write(*,*) "howdy!",CCTK_DELTA_TIME,cctk_delta_time,cctk_timefac

        if(bhns_rk4_particle_tracer_step==2) then
           ! Within this if() statement, we perform two rhs evaluations at t = t_n + 1/2 dT,
           ! performing the 2nd and 3rd RK4 steps.

           ! We are now at t = t_n + 1/2 dT

           do which_particle=1,N_particles_to_trace
              ! x_rhs(...,1) == f(t_n,y_n) == k1_x / dT
              k1_x = dT*x_rhs(which_particle,1)
              k1_y = dT*y_rhs(which_particle,1)
              k1_z = dT*z_rhs(which_particle,1)

              pointcoords(which_particle,1) = tracer_x(which_particle) + 0.5D0 * k1_x
              pointcoords(which_particle,2) = tracer_y(which_particle) + 0.5D0 * k1_y
              pointcoords(which_particle,3) = tracer_z(which_particle) + 0.5D0 * k1_z
           end do

!!$           x_rhs(:,bhns_rk4_particle_tracer_step) = 10.D0*sin(CCTK_TIME)
!!$           y_rhs(:,bhns_rk4_particle_tracer_step) = 10.D0*sin(CCTK_TIME)
!!$           z_rhs(:,bhns_rk4_particle_tracer_step) = 10.D0*sin(CCTK_TIME)
!!$
           call CCTK_VarIndex(vindex,"mhd_evolve::vx")
           call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,x_rhs(:,bhns_rk4_particle_tracer_step))
           call CCTK_VarIndex(vindex,"mhd_evolve::vy")
           call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,y_rhs(:,bhns_rk4_particle_tracer_step))
           call CCTK_VarIndex(vindex,"mhd_evolve::vz")
           call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,z_rhs(:,bhns_rk4_particle_tracer_step))

           bhns_rk4_particle_tracer_step=3

           write(*,*) " Particle Tracer substep # ",bhns_rk4_particle_tracer_step

           do which_particle=1,N_particles_to_trace
              ! x_rhs(...,2) == f(t_n + 1/2 dT , y_n + 1/2 k_1) == k2_x / dT
              k2_x = dT*x_rhs(which_particle,2)
              k2_y = dT*y_rhs(which_particle,2)
              k2_z = dT*z_rhs(which_particle,2)

              pointcoords(which_particle,1) = tracer_x(which_particle) + 0.5D0 * k2_x
              pointcoords(which_particle,2) = tracer_y(which_particle) + 0.5D0 * k2_y
              pointcoords(which_particle,3) = tracer_z(which_particle) + 0.5D0 * k2_z
           end do

!!$           x_rhs(:,bhns_rk4_particle_tracer_step) = 10.D0*sin(CCTK_TIME)
!!$           y_rhs(:,bhns_rk4_particle_tracer_step) = 10.D0*sin(CCTK_TIME)
!!$           z_rhs(:,bhns_rk4_particle_tracer_step) = 10.D0*sin(CCTK_TIME)
!!$
           call CCTK_VarIndex(vindex,"mhd_evolve::vx")
           call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,x_rhs(:,bhns_rk4_particle_tracer_step))
           call CCTK_VarIndex(vindex,"mhd_evolve::vy")
           call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,y_rhs(:,bhns_rk4_particle_tracer_step))
           call CCTK_VarIndex(vindex,"mhd_evolve::vz")
           call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,z_rhs(:,bhns_rk4_particle_tracer_step))

           bhns_rk4_particle_tracer_step=4

        else if(bhns_rk4_particle_tracer_step==4) then
           ! Within this if() statement, we perform two rhs evaluations at t = t_n + dT, 
           ! performing the 4th RK4 step for this RK4 timestep, and the 1st RK4 step for 
           ! the next RK4 timestep

           ! We are now at t = t_n + dT
           write(*,*) " Particle Tracer substep # ",bhns_rk4_particle_tracer_step

           do which_particle=1,N_particles_to_trace
              ! x_rhs(...,1) == f(t_n,y_n)
              k1_x = dT*x_rhs(which_particle,1)
              k1_y = dT*y_rhs(which_particle,1)
              k1_z = dT*z_rhs(which_particle,1)

              ! x_rhs(...,2) == f(t_n + 1/2 dT , y_n + 1/2 k_1) == k2_x / dT
              k2_x = dT*x_rhs(which_particle,2)
              k2_y = dT*y_rhs(which_particle,2)
              k2_z = dT*z_rhs(which_particle,2)

              ! x_rhs(...,3) == f(t_n + 1/2 dT , y_n + 1/2 k_2) == k3_x / dT
              k3_x = dT*x_rhs(which_particle,3)
              k3_y = dT*y_rhs(which_particle,3)
              k3_z = dT*z_rhs(which_particle,3)

              pointcoords(which_particle,1) = tracer_x(which_particle) + k3_x
              pointcoords(which_particle,2) = tracer_y(which_particle) + k3_y
              pointcoords(which_particle,3) = tracer_z(which_particle) + k3_z
           end do

!!$           x_rhs(:,bhns_rk4_particle_tracer_step) = 10.D0*sin(CCTK_TIME)
!!$           y_rhs(:,bhns_rk4_particle_tracer_step) = 10.D0*sin(CCTK_TIME)
!!$           z_rhs(:,bhns_rk4_particle_tracer_step) = 10.D0*sin(CCTK_TIME)
!!$
           call CCTK_VarIndex(vindex,"mhd_evolve::vx")
           call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,x_rhs(:,bhns_rk4_particle_tracer_step))
           call CCTK_VarIndex(vindex,"mhd_evolve::vy")
           call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,y_rhs(:,bhns_rk4_particle_tracer_step))
           call CCTK_VarIndex(vindex,"mhd_evolve::vz")
           call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,z_rhs(:,bhns_rk4_particle_tracer_step))

           ! Now that we have xyz_rhs at all 4 RK4 substeps, 
           !  we next construct k1, k2, k3, and k4, then 
           !  update tracer_xyz
           do which_particle=1,N_particles_to_trace
              ! x_rhs(...,1) == f(t_n,y_n)
              k1_x = dT*x_rhs(which_particle,1)
              k1_y = dT*y_rhs(which_particle,1)
              k1_z = dT*z_rhs(which_particle,1)

              ! x_rhs(...,2) == f(t_n + 1/2 dT , y_n + 1/2 k_1) == k2_x / dT
              k2_x = dT*x_rhs(which_particle,2)
              k2_y = dT*y_rhs(which_particle,2)
              k2_z = dT*z_rhs(which_particle,2)

              ! x_rhs(...,3) == f(t_n + 1/2 dT , y_n + 1/2 k_2) == k3_x / dT
              k3_x = dT*x_rhs(which_particle,3)
              k3_y = dT*y_rhs(which_particle,3)
              k3_z = dT*z_rhs(which_particle,3)

              ! x_rhs(...,4) == f(t_n + dT , y_n + k_3) == k4_x / dT
              k4_x = dT*x_rhs(which_particle,4)
              k4_y = dT*y_rhs(which_particle,4)
              k4_z = dT*z_rhs(which_particle,4)

              tracer_x(which_particle) = tracer_x(which_particle) + 1.D0/6.D0 * (k1_x + 2.D0*k2_x + 2.D0*k3_x + k4_x)
              tracer_y(which_particle) = tracer_y(which_particle) + 1.D0/6.D0 * (k1_y + 2.D0*k2_y + 2.D0*k3_y + k4_y)
              tracer_z(which_particle) = tracer_z(which_particle) + 1.D0/6.D0 * (k1_z + 2.D0*k2_z + 2.D0*k3_z + k4_z)

              ! Truncation error might conceivably push particle to z<0 with equatorial symmetry turned on, so correct the problem:
              if(Symmetry.eq.1 .and. tracer_z(which_particle).lt.0.D0) then
                 tracer_z(which_particle) = abs(tracer_z(which_particle))
              end if

           end do

           bhns_rk4_particle_tracer_step = 1

           write(*,*) " Particle Tracer substep # ",bhns_rk4_particle_tracer_step
           
           ! Next do the first RK4 substep!
           do which_particle=1,N_particles_to_trace
              pointcoords(which_particle,1) = tracer_x(which_particle)
              pointcoords(which_particle,2) = tracer_y(which_particle)
              pointcoords(which_particle,3) = tracer_z(which_particle)
           end do

!!$           x_rhs(:,bhns_rk4_particle_tracer_step) = 10.D0*sin(CCTK_TIME)
!!$           y_rhs(:,bhns_rk4_particle_tracer_step) = 10.D0*sin(CCTK_TIME)
!!$           z_rhs(:,bhns_rk4_particle_tracer_step) = 10.D0*sin(CCTK_TIME)
!!$
           call CCTK_VarIndex(vindex,"mhd_evolve::vx")
           call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,x_rhs(:,bhns_rk4_particle_tracer_step))
           call CCTK_VarIndex(vindex,"mhd_evolve::vy")
           call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,y_rhs(:,bhns_rk4_particle_tracer_step))
           call CCTK_VarIndex(vindex,"mhd_evolve::vz")
           call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,z_rhs(:,bhns_rk4_particle_tracer_step))
           
        end if
     end if
     !--------------------------------------------------------------------------------!
     !  FILE OUTPUT
     if((MOD(cctk_iteration,particle_tracer_substep_every)==0 .and. bhns_rk4_particle_tracer_step==1) &
          .or. cctk_iteration.eq.bhns_particle_tracer_start) then
        ! If we are at bhns_rk4_particle_tracer_step==1, it is time to dump the output files.
        ! Yes, this seems weird, but keep in mind that we must update the RHS twice every
        !  time this routine is called (except when setting up initial data): 
        !  The second & third RK4 substeps occur at t + 1/2 delta t
        !  The fourth RK4 substep occurs at t + delta t, and the first RK4 substep
        !  for the next timestep uses the same velocity data as the previous fourth RK4 
        !  substep (note that each time the rhs's are in general evaluated at different
        !  spatial points).

        ! These files contain tracer_xyz() data at the latest timestep.
        
        write(*,*) "Particle tracer: Dumping data files..."
        filename = 'bhns-particles.mon'

        ! First interpolate ( -1 - u_0 ) at all the up-to-date particle locations:
        do which_particle=1,N_particles_to_trace
           pointcoords(which_particle,1) = tracer_x(which_particle)
           pointcoords(which_particle,2) = tracer_y(which_particle)
           pointcoords(which_particle,3) = tracer_z(which_particle)
        end do

        call CCTK_VarIndex(vindex,"mhd_evolve::temp1")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,m1_minus_u_lower_zero)

        ! Then output everything to file.
        do which_particle=1,N_particles_to_trace
!           write(filename,41)which_particle
!41         FORMAT('bhns-particle.',I4.4)

           num_cols = 5
           export_data(which_particle,:) = 0.D0

           data_headers(1) = '# Time'
           export_data(which_particle,1) = CCTK_TIME

           data_headers(2) = 'x position'
           export_data(which_particle,2) = tracer_x(which_particle)

           data_headers(3) = 'y position'
           export_data(which_particle,3) = tracer_y(which_particle)

           data_headers(4) = 'z position'
           export_data(which_particle,4) = tracer_z(which_particle)

           data_headers(5) = '-1 - u_0'
           export_data(which_particle,5) = m1_minus_u_lower_zero(which_particle)
        end do
        if(cctk_iteration.eq.bhns_particle_tracer_start) then
           header_flag=1
        end if
                                       !output_multiline_data_to_file(filename,num_cols,col_names,num_rows,output_data,header_flag)
        if(CCTK_MyProc(CCTKGH)==0) call output_multiline_data_to_file(filename,num_cols,data_headers,N_particles_to_trace,export_data,header_flag)
        
        write(*,*) "Particle tracer: Finished dumping data files at time",CCTK_TIME
        write(*,*) "Particle tracer: Previous dump was at time",CCTK_TIME-dT

        ! We have already updated the xyz_rhs at the first RK4 substep (see above),
        ! so the next substep will be 2.
        bhns_rk4_particle_tracer_step = 2
     end if
     write(*,*) "****************"
  end if
  
  ! Next output MRI wavelength grid.
  if(cctk_iteration.ge.ITERATION_TO_output_MRI_wavelength .and. MOD(cctk_iteration,out_every)==0) then
     filename = 'bhns-MRI_wavelengths.mon'

     which_gridpoint=1
     do i=1,MRI_wavelength_calculator_Nxy
        do j=1,MRI_wavelength_calculator_Nxy
           MRI_wavelength_grid(which_gridpoint,1) = bh_posn_x(1) + (i - MRI_wavelength_calculator_Nxy*0.5D0)*MRI_wavelength_calculator_dxy
           MRI_wavelength_grid(which_gridpoint,2) = bh_posn_y(1) + (j - MRI_wavelength_calculator_Nxy*0.5D0)*MRI_wavelength_calculator_dxy
           MRI_wavelength_grid(which_gridpoint,3) = bh_posn_z(1)

           MRI_export_data(which_gridpoint,1) = CCTK_TIME
           MRI_export_data(which_gridpoint,2) = MRI_wavelength_grid(which_gridpoint,1)
           MRI_export_data(which_gridpoint,3) = MRI_wavelength_grid(which_gridpoint,2)
           MRI_export_data(which_gridpoint,4) = MRI_wavelength_grid(which_gridpoint,3)

           which_gridpoint=which_gridpoint+1
        end do
     end do
     which_gridpoint=which_gridpoint-1 ! fix off-by-one error

     call CCTK_VarIndex(vindex,"mhd_evolve::temp7")
     call interp_driver_carp(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,5))
     call CCTK_VarIndex(vindex,"mhd_evolve::temp8")
     call interp_driver_carp(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,6))
     call CCTK_VarIndex(vindex,"mhd_evolve::MONOPOLE")
     call interp_driver_carp(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,7))
     call CCTK_VarIndex(vindex,"mhd_evolve::rho_b")
     call interp_driver_carp(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,8))
 
     num_cols = 8
     
     data_headers(1) = '# Time'
     data_headers(2) = 'x position'
     data_headers(3) = 'y position'
     data_headers(4) = 'z position'
     data_headers(5) = 'Omega'
     data_headers(6) = 'MRI wavelength'
     data_headers(7) = 'MRI wavelength/dX'
     data_headers(8) = 'rho_b'

!     if(cctk_iteration.eq.ITERATION_TO_output_MRI_wavelength) then
     header_flag=1
 !    end if

     if(CCTK_MyProc(CCTKGH)==0) call output_multiline_data_to_file(filename,num_cols,data_headers,which_gridpoint,MRI_export_data,header_flag)
  end if

  if(1==0 .and. MOD(cctk_iteration,out_every)==0) then
     filename = 'bhns-Pandrhob.mon'

     which_gridpoint=1
     do i=1,MRI_wavelength_calculator_Nxy
        do j=1,MRI_wavelength_calculator_Nxy
           MRI_wavelength_grid(which_gridpoint,1) = bh_posn_x(1) + (i - MRI_wavelength_calculator_Nxy*0.5D0)*MRI_wavelength_calculator_dxy
           MRI_wavelength_grid(which_gridpoint,2) = bh_posn_y(1) + (j - MRI_wavelength_calculator_Nxy*0.5D0)*MRI_wavelength_calculator_dxy
           MRI_wavelength_grid(which_gridpoint,3) = bh_posn_z(1)

           MRI_export_data(which_gridpoint,1) = CCTK_TIME
           MRI_export_data(which_gridpoint,2) = MRI_wavelength_grid(which_gridpoint,1)
           MRI_export_data(which_gridpoint,3) = MRI_wavelength_grid(which_gridpoint,2)
           MRI_export_data(which_gridpoint,4) = MRI_wavelength_grid(which_gridpoint,3)

           which_gridpoint=which_gridpoint+1
        end do
     end do
     which_gridpoint=which_gridpoint-1 ! fix off-by-one error

     call CCTK_VarIndex(vindex,"mhd_evolve::rho_b")
     call interp_driver_carp(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,5))
     call CCTK_VarIndex(vindex,"mhd_evolve::P")
     call interp_driver_carp(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,6))
 
     num_cols = 6
     
     data_headers(1) = '# Time'
     data_headers(2) = 'x position'
     data_headers(3) = 'y position'
     data_headers(4) = 'z position'
     data_headers(5) = 'rho_b'
     data_headers(6) = 'P'

     if(cctk_iteration.eq.ITERATION_TO_output_MRI_wavelength) then
        header_flag=1
     end if

     if(CCTK_MyProc(CCTKGH)==0) call output_multiline_data_to_file(filename,num_cols,data_headers,which_gridpoint,MRI_export_data,header_flag)
  end if

  if(MOD(cctk_iteration,out_every)==0 .and. em_evolve_enable==1) then
!     if(cctk_iteration.eq.0) then
        header_flag=1
!     end if

     filename = 'bhns-max_tor_pol.mon'
     num_cols = 8
     export_data_singleline = 0.D0
     data_headers = 'unknowndude'

     data_headers(1) = '# Time'
     export_data_singleline(1) = CCTK_TIME

     data_headers(2) = 'B_tor_rhostar'
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_VarIndex(vindex,"mhd_evolve::temp3")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,numerator  ,1,vindex)

     export_data_singleline(2) = numerator

     data_headers(3) = 'B_pol_rhostar'
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_VarIndex(vindex,"mhd_evolve::temp4")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,numerator  ,1,vindex)

     export_data_singleline(3) = numerator

     data_headers(4) = 'rhostar'
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_VarIndex(vindex,"mhd_evolve::temp9")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,denominator,1,vindex)
     export_data_singleline(4) = denominator

     data_headers(5) = 'b^2/P_3'
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_VarIndex(vindex,"mhd_evolve::temp5")
     call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, export_data_singleline(5), 1, vindex)
     export_data_singleline(5) = export_data_singleline(4)*cctk_delta_space(1)*cctk_delta_space(2)*cctk_delta_space(3)
 
     data_headers(6) = 'rho_star_3'
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_VarIndex(vindex,"mhd_evolve::temp6")
     call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, export_data_singleline(6), 1, vindex) 
     export_data_singleline(6) = export_data_singleline(5)*cctk_delta_space(1)*cctk_delta_space(2)*cctk_delta_space(3)

     data_headers(7) = 'B_tor_max'
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex(vindex,"mhd_evolve::temp10")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data_singleline(7),1,vindex)

     data_headers(8) = 'B_pol_max'
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex(vindex,"mhd_evolve::temp11")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data_singleline(8),1,vindex)


     
!!$     data_headers(6) = 'b^2/P_4'
!!$     call CCTK_ReductionHandle(handle,"sum")
!!$     call CCTK_VarIndex(vindex,"mhd_evolve::temp9")
!!$     call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, export_data_singleline(6), 1, vindex)
!!$     export_data_singleline(6) = export_data_singleline(6)*cctk_delta_space(1)*cctk_delta_space(2)*cctk_delta_space(3)
!!$ 
!!$     data_headers(7) = 'rho_star_4'
!!$     call CCTK_ReductionHandle(handle,"sum")
!!$     call CCTK_VarIndex(vindex,"mhd_evolve::temp11")
!!$     call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, export_data_singleline(7), 1, vindex) 
!!$     export_data_singleline(7) = export_data_singleline(7)*cctk_delta_space(1)*cctk_delta_space(2)*cctk_delta_space(3)
!!$
!!$     data_headers(8) = 'b^2/P_2'
!!$     call CCTK_ReductionHandle(handle,"sum")
!!$     call CCTK_VarIndex(vindex,"mhd_evolve::temp10")
!!$     call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, export_data_singleline(8), 1, vindex)
!!$     export_data_singleline(8) = export_data_singleline(8)*cctk_delta_space(1)*cctk_delta_space(2)*cctk_delta_space(3)
!!$ 
!!$     data_headers(9) = 'rho_star_2'
!!$     call CCTK_ReductionHandle(handle,"sum")
!!$     call CCTK_VarIndex(vindex,"mhd_evolve::temp2")
!!$     call CCTK_Reduce(ierr,cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, export_data_singleline(9), 1, vindex) 
!!$     export_data_singleline(9) = export_data_singleline(9)*cctk_delta_space(1)*cctk_delta_space(2)*cctk_delta_space(3)

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data_singleline,header_flag)
 
  end if

end subroutine bhns_diagnostics_trace_particles
