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
  real*8, dimension(N_particles_to_trace,6) :: export_data
  character, dimension(50)                   :: data_headers*20
  integer                                    :: num_cols,header_flag

  ! Local variables needed for RK4 timestepping:
  real*8                                   :: dT,phiangle
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
  integer                                  :: i,j,k,which_gridpoint
  real*8,dimension(MRI_wavelength_calculator_Nxy*MRI_wavelength_calculator_Nxy,3) :: MRI_wavelength_grid
  real*8,dimension(MRI_wavelength_calculator_Nxy*MRI_wavelength_calculator_Nxy,16) :: MRI_export_data
  !------------------------------------------------------------!
  ! Max B poloidal, B toroidal diagnostic, as well as rho_star weighted b^2/P diagnostic.
  real*8, dimension(50) :: export_data_singleline
  real*8                :: numerator,denominator, x0, y0, z0
  integer :: handle,ierr
  !------------------------------------------------------------!
  integer :: NO_SYMM, EQUATORIAL
  parameter(NO_SYMM = 0, EQUATORIAL = 1)

  ! define bad particles array
!  real,dimension(N_particles_to_trace) :: bad_particles


  header_flag=0

  select case (particle_center)
     case(0)  ! Choose center in the par file
        x0=bhns_tracer_x0   
        y0=bhns_tracer_y0 
        z0=bhns_tracer_z0 
     case(1)  ! Center on box 1
        x0=position_x(1)
        y0=position_y(1) 
        z0=position_z(1)
     case(2)  ! Center on box 2
        x0=position_x(2)
        y0=position_y(2)
        z0=position_z(2)
     case(3)  ! Center on BH
        x0=bh_posn_x(1)
        y0=bh_posn_y(1)
        z0=bh_posn_z(1)
     case(4)  ! Center on NS
        x0=CoMx_VolInt/CoM_VolInt_denominator
        y0=CoMy_VolInt/CoM_VolInt_denominator
        if(Symmetry==EQUATORIAL) then
           z0=0.d0
        else
           z0=CoMz_VolInt/CoM_VolInt_denominator
        end if
  end select


  if(cctk_iteration.ge.bhns_particle_tracer_start) then

     if(Symmetry.ne.0 .and. Symmetry.ne.1) then
        write(*,*) "Sorry, particle tracers only work when Symmetry==0 or 1."
        stop
     end if
     
     ! set bad particles to zero
!     bad_particles  = 0.d0

     !-----------------------------------------------------------------------------
     ! INITIAL DATA

     if(cctk_iteration.eq.bhns_particle_tracer_start) then  !start ID
        write(*,*) "In particle tracer initial data."
        write(*,*) "Particle distribution center x,y,z=",x0,y0,z0
        write(*,*) "Particle distribution Geometry =", initial_particle_geometry
        write(*,*) "Particle distribution Large radius =", bhns_tracer_r
        write(*,*) "Particle distribution Small radius =", bhns_tracer_rin
        if (initial_particle_geometry.eq.1) then
           write(*,*) "particle_cylinder_cone_zmax =", particle_cylinder_cone_zmax
           write(*,*) "particle_cylinder_cone_zmin =", particle_cylinder_cone_zmin
        end if


        ! Don't change the default random seed!
        !CALL RANDOM_SEED                    ! Initialize random number generator

        do which_particle=1,N_particles_to_trace
           ! Seed particles inside sphere centered on coords of NS, with radius r = bhns_tracer_r

           ! First choose a guess that is very far from the star, so that we go inside the do while loop below.
           tracer_x(which_particle) = bhns_tracer_r + 10000.D0
           tracer_y(which_particle) = 0.D0
           tracer_z(which_particle) = 0.D0

           ! Next we throw particles into a cube, and if they fall within the specified geometry, we keep them:
           if (initial_particle_geometry.eq.0) then ! spherical geometry
              do while(sqrt((tracer_x(which_particle)-x0)**2 + (tracer_y(which_particle)-y0)**2 + (tracer_z(which_particle)-z0)**2) > bhns_tracer_r)

                 CALL RANDOM_NUMBER(tracer_x(which_particle))
                 CALL RANDOM_NUMBER(tracer_y(which_particle))
                 CALL RANDOM_NUMBER(tracer_z(which_particle))
                 
                 ! tracer_xyz(which_particle) is now between 0 and 1. We want it between -bhns_tracer_r + x0 to bhns_tracer_r + x0:
                 tracer_x(which_particle) = 2.D0*bhns_tracer_r*tracer_x(which_particle) + x0 - bhns_tracer_r
                 tracer_y(which_particle) = 2.D0*bhns_tracer_r*tracer_y(which_particle) + y0 - bhns_tracer_r
                 tracer_z(which_particle) = 2.D0*bhns_tracer_r*tracer_z(which_particle) + z0 - bhns_tracer_r
                 
              end do
           else if (initial_particle_geometry.eq.1) then ! cylindrical geometry
              do while(sqrt((tracer_x(which_particle)-x0)**2 + (tracer_y(which_particle)-y0)**2) > bhns_tracer_r .or. sqrt((tracer_x(which_particle)-x0)**2 + (tracer_y(which_particle)-y0)**2) < bhns_tracer_rin)
                 
                 call random_number(tracer_x(which_particle))
                 CALL RANDOM_NUMBER(tracer_y(which_particle))
                 CALL RANDOM_NUMBER(tracer_z(which_particle))
                 
                 ! tracer_xyz(which_particle) is now between 0 and 1. We want it between -bhns_tracer_r + x0 to bhns_tracer_r + x0:
                 tracer_x(which_particle) = 2.D0*bhns_tracer_r*tracer_x(which_particle) + x0 - bhns_tracer_r
                 tracer_y(which_particle) = 2.D0*bhns_tracer_r*tracer_y(which_particle) + y0 - bhns_tracer_r
                 
                 ! Need to change z to be within zmin and zmax
                 tracer_z(which_particle) = (particle_cylinder_cone_zmax-particle_cylinder_cone_zmin)*tracer_z(which_particle) + particle_cylinder_cone_zmin
              end do
           else 
              write(*,*) "Stopping: initial particle geometry not supported: change bhns::initial_particle_geometry parameter"
           end if
           

           ! Shift the center of the ball of particles to the NS center of mass
!           tracer_x(which_particle) = tracer_x(which_particle) + position_x(2)
!           tracer_y(which_particle) = tracer_y(which_particle) + position_y(2)
!           tracer_z(which_particle) = tracer_z(which_particle) + position_z(2)

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

           if(pointcoords(which_particle,1).gt.tracer_x_max) pointcoords(which_particle,1)=tracer_x_max
           if(pointcoords(which_particle,2).gt.tracer_y_max) pointcoords(which_particle,2)=tracer_y_max
           if(pointcoords(which_particle,3).gt.tracer_z_max) pointcoords(which_particle,3)=tracer_z_max

           !If we have a particle at z<0 with equatorial symmetry turned on, set it to z>0:
           if(Symmetry.eq.1 .and. tracer_z(which_particle).lt.0.D0) then
              tracer_z(which_particle) = abs(tracer_z(which_particle))
           end if

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

        write(*,*) "Done with particle tracer initial data"
     end if ! End of initial data routine.
     !-----------------------------------------------------------------------------

     !-----------------------------------------------------------------------------
     ! RK4 TIMESTEPPING

     if(mod(cctk_iteration,particle_tracer_substep_every)==0) then
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
              
              if((pointcoords(which_particle,1).gt.tracer_x_max).or.(pointcoords(which_particle,1).lt.tracer_x_min).or.(isnan(pointcoords(which_particle,1))).or. &
                 (pointcoords(which_particle,2).gt.tracer_y_max).or.(pointcoords(which_particle,2).lt.tracer_y_min).or.(isnan(pointcoords(which_particle,2))).or. &
                 (pointcoords(which_particle,3).gt.tracer_z_max).or.(pointcoords(which_particle,3).lt.tracer_z_min).or.(isnan(pointcoords(which_particle,3)))) then
                 pointcoords(which_particle,1)  = 0.d0
                 pointcoords(which_particle,2)  = 0.d0
                 pointcoords(which_particle,3)  = 0.d0
              end if
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

              if((pointcoords(which_particle,1).gt.tracer_x_max).or.(pointcoords(which_particle,1).lt.tracer_x_min).or.(isnan(pointcoords(which_particle,1))).or. &
                 (pointcoords(which_particle,2).gt.tracer_y_max).or.(pointcoords(which_particle,2).lt.tracer_y_min).or.(isnan(pointcoords(which_particle,2))).or. &
                 (pointcoords(which_particle,3).gt.tracer_z_max).or.(pointcoords(which_particle,3).lt.tracer_z_min).or.(isnan(pointcoords(which_particle,3)))) then
                 pointcoords(which_particle,1)  = 0.d0
                 pointcoords(which_particle,2)  = 0.d0
                 pointcoords(which_particle,3)  = 0.d0
              end if

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
              
              pointcoords(which_particle,1) = tracer_x(which_particle) + 0.5D0 * k2_x
              pointcoords(which_particle,2) = tracer_y(which_particle) + 0.5D0 * k2_y
              pointcoords(which_particle,3) = tracer_z(which_particle) + 0.5D0 * k2_z

              if((pointcoords(which_particle,1).gt.tracer_x_max).or.(pointcoords(which_particle,1).lt.tracer_x_min).or.(isnan(pointcoords(which_particle,1))).or. &
                 (pointcoords(which_particle,2).gt.tracer_y_max).or.(pointcoords(which_particle,2).lt.tracer_y_min).or.(isnan(pointcoords(which_particle,2))).or. &
                 (pointcoords(which_particle,3).gt.tracer_z_max).or.(pointcoords(which_particle,3).lt.tracer_z_min).or.(isnan(pointcoords(which_particle,3)))) then
                 pointcoords(which_particle,1)  = 0.d0
                 pointcoords(which_particle,2)  = 0.d0
                 pointcoords(which_particle,3)  = 0.d0
              end if

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
              if((Symmetry.eq.1 .and. tracer_z(which_particle).lt.0.D0)) then
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

              pointcoords(which_particle,1) = tracer_x(which_particle) + 0.5D0 * k2_x
              pointcoords(which_particle,2) = tracer_y(which_particle) + 0.5D0 * k2_y
              pointcoords(which_particle,3) = tracer_z(which_particle) + 0.5D0 * k2_z

              if((pointcoords(which_particle,1).gt.tracer_x_max).or.(pointcoords(which_particle,1).lt.tracer_x_min).or.(isnan(pointcoords(which_particle,1))).or. &
                 (pointcoords(which_particle,2).gt.tracer_y_max).or.(pointcoords(which_particle,2).lt.tracer_y_min).or.(isnan(pointcoords(which_particle,2))).or. &
                 (pointcoords(which_particle,3).gt.tracer_z_max).or.(pointcoords(which_particle,3).lt.tracer_z_min).or.(isnan(pointcoords(which_particle,3)))) then
                 pointcoords(which_particle,1)  = 0.d0
                 pointcoords(which_particle,2)  = 0.d0
                 pointcoords(which_particle,3)  = 0.d0
              end if

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
     if((mod(cctk_iteration,particle_tracer_substep_every)==0 .and. bhns_rk4_particle_tracer_step==1) &
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

        !  check the position of the particles

        do which_particle=1,N_particles_to_trace
           if((isnan(tracer_x(which_particle))).or.(isnan(tracer_y(which_particle))).or.(isnan(tracer_z(which_particle)))) then
              tracer_x(which_particle)=tracer_x_max 
              tracer_y(which_particle)=tracer_y_max 
              tracer_z(which_particle)=tracer_z_max 
              
!              bad_particles(which_particle)  = 1.d0
           end if
        end do


        ! Frist interpolate the density
        call CCTK_VarIndex(vindex,"mhd_evolve::rho_b")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,local_density)

        ! First interpolate ( -1 - u_0 ) at all the up-to-date particle locations:
        do which_particle=1,N_particles_to_trace
           pointcoords(which_particle,1) = tracer_x(which_particle)
           pointcoords(which_particle,2) = tracer_y(which_particle)
           pointcoords(which_particle,3) = tracer_z(which_particle)

!           if(pointcoords(which_particle,1).gt.tracer_x_max) pointcoords(which_particle,1)=tracer_x_max
!           if(pointcoords(which_particle,2).gt.tracer_y_max) pointcoords(which_particle,2)=tracer_y_max
!           if(pointcoords(which_particle,3).gt.tracer_z_max) pointcoords(which_particle,3)=tracer_z_max

           if((pointcoords(which_particle,1).gt.tracer_x_max).or.(pointcoords(which_particle,1).lt.tracer_x_min).or.(isnan(pointcoords(which_particle,1))).or. &
                (pointcoords(which_particle,2).gt.tracer_y_max).or.(pointcoords(which_particle,2).lt.tracer_y_min).or.(isnan(pointcoords(which_particle,2))).or. &
                (pointcoords(which_particle,3).gt.tracer_z_max).or.(pointcoords(which_particle,3).lt.tracer_z_min).or.(isnan(pointcoords(which_particle,3)))) then
              pointcoords(which_particle,1)  = 0.d0
              pointcoords(which_particle,2)  = 0.d0
              pointcoords(which_particle,3)  = 0.d0
           end if
        end do

        call CCTK_VarIndex(vindex,"mhd_evolve::temp1")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,m1_minus_u_lower_zero)


        ! reset the position of bad particles to a given value       

!        do which_particle=1,N_particles_to_trace
!           if(bad_particles(which_particle)==1.d0) then
!              tracer_x(which_particle) = 1.0d10
!              tracer_y(which_particle) = 1.0d10
!              tracer_z(which_particle) = 1.0d10
!           end if
!        end do

        ! Then output everything to file.
        do which_particle=1,N_particles_to_trace
           !           write(filename,41)which_particle
           !41         FORMAT('bhns-particle.',I4.4)

           num_cols = 6
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

           data_headers(6) = 'rest-mass density'
           export_data(which_particle,6) = local_density(which_particle)
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

           data_headers(1) = '# Time'
           MRI_export_data(which_gridpoint,1) = CCTK_TIME
           data_headers(2) = 'x position'
           MRI_export_data(which_gridpoint,2) = MRI_wavelength_grid(which_gridpoint,1)
           data_headers(3) = 'y position'
           MRI_export_data(which_gridpoint,3) = MRI_wavelength_grid(which_gridpoint,2)
           data_headers(4) = 'z position'
           MRI_export_data(which_gridpoint,4) = MRI_wavelength_grid(which_gridpoint,3)

           which_gridpoint=which_gridpoint+1
        end do
     end do
     which_gridpoint=which_gridpoint-1 ! fix off-by-one error

     ! temp2 = b^2/(2*P)*rho_star, for rho>1e3 rho_atm
     ! temp3 = B_toroidal*rho_star, for rho>1e3 rho_atm
     ! temp4 = B_poloidal*rho_star, for rho>1e3 rho_atm
     ! temp5 = dX, for rho>1e3 rho_atm
     ! temp6 = K, for rho>1e3 rho_atm
     ! temp7 = Omega, for rho>1e3 rho_atm
     ! temp8 = lambda_MRI, for rho>1e3 rho_atm
     ! temp9 = rho_star, for rho>1e3 rho_atm
     ! temp10 = B_toroidal
     ! temp11 = B_poloidal
     ! MONOPOLE = rho_star
     ! P_thermal = b^2

!here

     data_headers(5) = 'dX'
     call CCTK_VarIndex(vindex,"mhd_evolve::temp5")
     call bhns_interp_driver_carp_arborder(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,5),1)
     data_headers(6) = 'K'
     call CCTK_VarIndex(vindex,"mhd_evolve::temp6")
     call bhns_interp_driver_carp_arborder(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,6),1)
     data_headers(7) = 'Omega'
     call CCTK_VarIndex(vindex,"mhd_evolve::temp7")
     call bhns_interp_driver_carp_arborder(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,7),1)
     data_headers(8) = 'lambda_MRI'
     call CCTK_VarIndex(vindex,"mhd_evolve::temp8")
     call bhns_interp_driver_carp_arborder(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,8),1)
     data_headers(9) = 'rho_b'
     call CCTK_VarIndex(vindex,"mhd_evolve::rho_b")
     call bhns_interp_driver_carp_arborder(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,9),1)
     data_headers(10) = 'B_toroidal'
     call CCTK_VarIndex(vindex,"mhd_evolve::temp10")
     call bhns_interp_driver_carp_arborder(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,10),1)
     data_headers(11) = 'B_toroidal'
     call CCTK_VarIndex(vindex,"mhd_evolve::temp11")
     call bhns_interp_driver_carp_arborder(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,11),1)
     data_headers(12) = 'P'
     call CCTK_VarIndex(vindex,"mhd_evolve::P")
     call bhns_interp_driver_carp_arborder(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,12),1)
     data_headers(13) = 'b^2'
     call CCTK_VarIndex(vindex,"mhd_evolve::P_thermal")
     call bhns_interp_driver_carp_arborder(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,13),1)
     data_headers(14) = 'Ax'
     call CCTK_VarIndex(vindex,"mhd_evolve::Ax")
     call bhns_interp_driver_carp_arborder(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,14),1)
     data_headers(15) = 'Ay'
     call CCTK_VarIndex(vindex,"mhd_evolve::Ay")
     call bhns_interp_driver_carp_arborder(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,15),1)
     data_headers(16) = 'Az'
     call CCTK_VarIndex(vindex,"mhd_evolve::Az")
     call bhns_interp_driver_carp_arborder(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,16),1)

     num_cols = 16

     !     if(cctk_iteration.eq.ITERATION_TO_output_MRI_wavelength) then
     header_flag=1
     !    end if

     if(CCTK_MyProc(CCTKGH)==0) call output_multiline_data_to_file(filename,num_cols,data_headers,which_gridpoint,MRI_export_data,header_flag)

     !************************************************************
          filename = 'bhns-b2_8pirho.mon'

     which_gridpoint=1
     do i=1,MRI_wavelength_calculator_Nxy
        do k=1,MRI_wavelength_calculator_Nxy
           MRI_wavelength_grid(which_gridpoint,1) = bh_posn_x(1) + (i - MRI_wavelength_calculator_Nxy*0.5D0)*MRI_wavelength_calculator_dxy
           MRI_wavelength_grid(which_gridpoint,2) = bh_posn_y(1) 
           
           if(bhns_domain .eq. 0) then
              MRI_wavelength_grid(which_gridpoint,3) = bh_posn_z(1) + k*MRI_wavelength_calculator_dxy           
           else if(bhns_domain .eq. 1) then
              MRI_wavelength_grid(which_gridpoint,3) = bh_posn_z(1) + (k - MRI_wavelength_calculator_Nxy*0.5D0)*MRI_wavelength_calculator_dxy
           else
              print *,"Check your domain. bhns-b2_8pirho.mon assumes either full or bitant domains "
              stop
           end if

           data_headers(1) = '# Time'
           MRI_export_data(which_gridpoint,1) = CCTK_TIME
           data_headers(2) = 'x position'
           MRI_export_data(which_gridpoint,2) = MRI_wavelength_grid(which_gridpoint,1)
           data_headers(3) = 'y position'
           MRI_export_data(which_gridpoint,3) = MRI_wavelength_grid(which_gridpoint,2)
           data_headers(4) = 'z position'
           MRI_export_data(which_gridpoint,4) = MRI_wavelength_grid(which_gridpoint,3)

           which_gridpoint=which_gridpoint+1
        end do
     end do
     which_gridpoint=which_gridpoint-1 ! fix off-by-one error

     ! temp2 = b^2/(2*P)*rho_star, for rho>1e3 rho_atm
     ! temp3 = B_toroidal*rho_star, for rho>1e3 rho_atm
     ! temp4 = B_poloidal*rho_star, for rho>1e3 rho_atm
     ! temp5 = dX, for rho>1e3 rho_atm
     ! temp6 = K, for rho>1e3 rho_atm
     ! temp7 = Omega, for rho>1e3 rho_atm
     ! temp8 = lambda_MRI, for rho>1e3 rho_atm
     ! temp9 = rho_star, for rho>1e3 rho_atm
     ! temp10 = B_toroidal
     ! temp11 = B_poloidal
     ! MONOPOLE = rho_star
     ! P_thermal = b^2

     data_headers(5) = 'b^2'
     call CCTK_VarIndex(vindex,"mhd_evolve::P_thermal")
     call bhns_interp_driver_carp_arborder(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,5),1)
     data_headers(6) = 'rho_b'
     call CCTK_VarIndex(vindex,"mhd_evolve::rho_b")
     call bhns_interp_driver_carp_arborder(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,6),1)
     data_headers(7) = 'vx'
     call CCTK_VarIndex(vindex,"mhd_evolve::vx")
     call bhns_interp_driver_carp_arborder(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,7),1)
     data_headers(8) = 'vy'
     call CCTK_VarIndex(vindex,"mhd_evolve::vy")
     call bhns_interp_driver_carp_arborder(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,8),1)
     data_headers(9) = 'vz'
     call CCTK_VarIndex(vindex,"mhd_evolve::vz")
     call bhns_interp_driver_carp_arborder(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,9),1)
     data_headers(10) = 'P'
     call CCTK_VarIndex(vindex,"mhd_evolve::P")
     call bhns_interp_driver_carp_arborder(cctkGH,which_gridpoint,MRI_wavelength_grid,vindex,MRI_export_data(:,10),1)

     num_cols = 10
     header_flag=1
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
     num_cols = 7
     export_data_singleline = 0.D0
     data_headers = 'unknowndude'

     data_headers(1) = '# Time'
     export_data_singleline(1) = CCTK_TIME

     ! temp2 = b^2/(2*P)*rho_star, for rho>1e3 rho_atm
     ! temp3 = B_toroidal*rho_star, for rho>1e3 rho_atm
     ! temp4 = B_poloidal*rho_star, for rho>1e3 rho_atm
     ! temp5 = dX, for rho>1e3 rho_atm
     ! temp6 = K, for rho>1e3 rho_atm
     ! temp7 = Omega, for rho>1e3 rho_atm
     ! temp8 = lambda_MRI, for rho>1e3 rho_atm
     ! temp9 = rho_star, for rho>1e3 rho_atm
     ! temp10 = B_toroidal
     ! temp11 = B_poloidal
     ! MONOPOLE = rho_star
     ! P_thermal = lambda_MRI*rho_star, for rho>1e3 rho_atm

     data_headers(2) = 'B_tor_rhostarN_3'
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_VarIndex(vindex,"mhd_evolve::temp3")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data_singleline(2),1,vindex)

     data_headers(3) = 'B_pol_rhostarN_3'
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_VarIndex(vindex,"mhd_evolve::temp4")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data_singleline(3),1,vindex)

     data_headers(4) = 'B^2/(2P)rhostarN_3'
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_VarIndex(vindex,"mhd_evolve::temp2")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data_singleline(4),1,vindex)

     data_headers(5) = 'rhostar_D_3'
     call CCTK_ReductionHandle(handle,"sum")
     call CCTK_VarIndex(vindex,"mhd_evolve::temp9")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data_singleline(5),1,vindex)

     data_headers(6) = 'B_tor_max'
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex(vindex,"mhd_evolve::temp10")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data_singleline(6),1,vindex)

     data_headers(7) = 'B_pol_max'
     call CCTK_ReductionHandle(handle,"maximum")
     call CCTK_VarIndex(vindex,"mhd_evolve::temp11")
     call CCTK_Reduce (ierr,cctkGH, -1,handle, 1, CCTK_VARIABLE_REAL,export_data_singleline(7),1,vindex)

     if(CCTK_MyProc(CCTKGH)==0) call output_data_to_file(filename,num_cols,data_headers,export_data_singleline,header_flag)

  end if

end subroutine bhns_diagnostics_trace_particles
