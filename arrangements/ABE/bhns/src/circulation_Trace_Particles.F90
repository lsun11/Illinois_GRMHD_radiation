#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
!-----------------------------------------------------------------------------
! Particle tracer to compute the circulation
! This routine has been optimized, and should not slow down
! your run 
!
! Read the code & comments carefully to figure out the timestepping.
!-----------------------------------------------------------------------------
  
subroutine circulation_trace_particles(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  ! Local variables needed for file output:
  character                                  :: filename*50
  real*8, dimension(N_particles_to_trace,6)  :: export_data
  character, dimension(50)                   :: data_headers*20
  integer                                    :: num_cols,header_flag

  ! Local variables needed for RK4 timestepping:
  real*8                                   :: dT
  real*8                                   :: k1_x,k1_y,k1_z,k1_taul
  real*8                                   :: k2_x,k2_y,k2_z,k2_taul
  real*8                                   :: k3_x,k3_y,k3_z,k3_taul
  real*8                                   :: k4_x,k4_y,k4_z,k4_taul

  integer, dimension(3)                    :: ext

  ! Local variables needed for interpolation:
  real*8,dimension(N_particles_to_trace,4) :: pointcoords
  integer                                  :: vindex

  ! Loop dummy variable
  integer                                  :: which_particle
  real*8                                   :: x0, y0,phiangle,output_integral
  integer :: handle,ierr
  !------------------------------------------------------------!

  header_flag=0



  if(cctk_iteration.ge.bhns_particle_tracer_start) then 
  
  ! -----------------------
  ! ---  INITIAL DATA   ---
  ! -----------------------

  if(cctk_iteration.eq.0) then     

  ! position at t =0
  x0 = initial_ns_coord_x
  y0 = initial_ns_coord_y

  print *,"x0,y0=",x0,y0

  ! define angle
  phiangle  =  2.0d0*acos(-1.0d0)/(N_particles_to_trace-1)


     do which_particle=1,N_particles_to_trace
        ! define particles in a ring of radius at proper time = tau.
        ! notice that I'm using the z coordinate as the proper time,
        ! which means each point has coordinates (x,y,tau)
        tracer_x(which_particle)      =  bhns_R_NS*cos(which_particle*phiangle) + x0
        tracer_y(which_particle)      =  bhns_R_NS*sin(which_particle*phiangle) + y0
        tracer_z(which_particle)      =  0.0d0
        tracer_taul(which_particle)   =  0.0d0
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
        pointcoords(which_particle,4) = tracer_z(which_particle)
     end do
     
     ! compute 3-velocities 
     call CCTK_VarIndex(vindex,"mhd_evolve::vx")
     call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,x_rhs(:,bhns_rk4_particle_tracer_step))
     call CCTK_VarIndex(vindex,"mhd_evolve::vy")
     call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,y_rhs(:,bhns_rk4_particle_tracer_step))
     call CCTK_VarIndex(vindex,"mhd_evolve::vz")
     call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,z_rhs(:,bhns_rk4_particle_tracer_step))

     ! compute 1/u^t
     call CCTK_VarIndex(vindex,"mhd_evolve::temp18")
     call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,taul_rhs(:,bhns_rk4_particle_tracer_step))

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
        write(*,*) " Particle Tracer circle substep # ",bhns_rk4_particle_tracer_step,"dt=",dT
        write(*,*) " Particle Tracer circle thinks the time is =",CCTK_TIME

        write(*,*) "howdy!",CCTK_DELTA_TIME,cctk_delta_time,cctk_timefac

     if(bhns_rk4_particle_tracer_step==2) then
        ! Within this if() statement, we perform two rhs evaluations at t = t_n + 1/2 dT,
        ! performing the 2nd and 3rd RK4 steps.
        
        ! We are now at t = t_n + 1/2 dT
        
        do which_particle=1,N_particles_to_trace
           ! x_rhs(...,1) == f(t_n,y_n) == k1_x / dT
           k1_x     = dT*   x_rhs(which_particle,1)
           k1_y     = dT*   y_rhs(which_particle,1)
           k1_z     = dT*   z_rhs(which_particle,1)
           k1_taul  = dT*taul_rhs(which_particle,1)

           pointcoords(which_particle,1) =    tracer_x(which_particle) + 0.5D0 * k1_x
           pointcoords(which_particle,2) =    tracer_y(which_particle) + 0.5D0 * k1_y
           pointcoords(which_particle,3) =    tracer_z(which_particle) + 0.5D0 * k1_z
           pointcoords(which_particle,4) = tracer_taul(which_particle) + 0.5D0 * k1_taul
        end do
        
        ! recompute velocities           
        call CCTK_VarIndex(vindex,"mhd_evolve::vx")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,x_rhs(:,bhns_rk4_particle_tracer_step))
        call CCTK_VarIndex(vindex,"mhd_evolve::vy")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,y_rhs(:,bhns_rk4_particle_tracer_step))
        call CCTK_VarIndex(vindex,"mhd_evolve::vz")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,z_rhs(:,bhns_rk4_particle_tracer_step))

        call CCTK_VarIndex(vindex,"mhd_evolve::temp18")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,taul_rhs(:,bhns_rk4_particle_tracer_step))
        
        bhns_rk4_particle_tracer_step=3
        
        write(*,*) " Particle Tracer substep # ",bhns_rk4_particle_tracer_step
        
        do which_particle=1,N_particles_to_trace
           ! x_rhs(...,2) == f(t_n + 1/2 dT , y_n + 1/2 k_1) == k2_x / dT
           k2_x    = dT*   x_rhs(which_particle,2)
           k2_y    = dT*   y_rhs(which_particle,2)
           k2_z    = dT*   z_rhs(which_particle,2)
           k2_taul = dT*taul_rhs(which_particle,2)
           
           pointcoords(which_particle,1) = tracer_x(which_particle)    + 0.5D0 * k2_x
           pointcoords(which_particle,2) = tracer_y(which_particle)    + 0.5D0 * k2_y
           pointcoords(which_particle,3) = tracer_z(which_particle)    + 0.5D0 * k2_z
           pointcoords(which_particle,4) = tracer_taul(which_particle) + 0.5D0 * k2_taul
        end do
        
        
        ! recompute velocities           
        call CCTK_VarIndex(vindex,"mhd_evolve::vx")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,x_rhs(:,bhns_rk4_particle_tracer_step))
        call CCTK_VarIndex(vindex,"mhd_evolve::vy")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,y_rhs(:,bhns_rk4_particle_tracer_step))
        call CCTK_VarIndex(vindex,"mhd_evolve::vz")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,z_rhs(:,bhns_rk4_particle_tracer_step))
        
        ! recompute 1/u^t
        call CCTK_VarIndex(vindex,"mhd_evolve::temp18")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,taul_rhs(:,bhns_rk4_particle_tracer_step))

        bhns_rk4_particle_tracer_step=4
        
     else if(bhns_rk4_particle_tracer_step==4) then
        ! Within this if() statement, we perform two rhs evaluations at t = t_n + dT, 
        ! performing the 4th RK4 step for this RK4 timestep, and the 1st RK4 step for 
        ! the next RK4 timestep
        
        ! We are now at t = t_n + dT
        write(*,*) " Particle Tracer substep # ",bhns_rk4_particle_tracer_step
        
        do which_particle=1,N_particles_to_trace
           ! x_rhs(...,1) == f(t_n,y_n)
           k1_x     = dT*   x_rhs(which_particle,1)
           k1_y     = dT*   y_rhs(which_particle,1)
           k1_z     = dT*   z_rhs(which_particle,1)
           k1_taul  = dT*taul_rhs(which_particle,1)
           
           ! x_rhs(...,2) == f(t_n + 1/2 dT , y_n + 1/2 k_1) == k2_x / dT
           k2_x      = dT*   x_rhs(which_particle,2)
           k2_y      = dT*   y_rhs(which_particle,2)
           k2_z      = dT*   z_rhs(which_particle,2)
           k2_taul   = dT*taul_rhs(which_particle,2)
           
           ! x_rhs(...,3) == f(t_n + 1/2 dT , y_n + 1/2 k_2) == k3_x / dT
           k3_x     = dT*   x_rhs(which_particle,3)
           k3_y     = dT*   y_rhs(which_particle,3)
           k3_z     = dT*   z_rhs(which_particle,3)
           k3_taul  = dT*taul_rhs(which_particle,3)
           
           pointcoords(which_particle,1) = tracer_x(which_particle)    + k3_x
           pointcoords(which_particle,2) = tracer_y(which_particle)    + k3_y
           pointcoords(which_particle,3) = tracer_z(which_particle)    + k3_z
           pointcoords(which_particle,4) = tracer_taul(which_particle) + k3_taul
           
           pointcoords(which_particle,1) = tracer_x(which_particle)    + 0.5D0 * k2_x
           pointcoords(which_particle,2) = tracer_y(which_particle)    + 0.5D0 * k2_y
           pointcoords(which_particle,3) = tracer_z(which_particle)    + 0.5D0 * k2_z
           pointcoords(which_particle,4) = tracer_taul(which_particle) + 0.5D0 * k2_taul
        end do
        
        
        ! recompute velocities           
        call CCTK_VarIndex(vindex,"mhd_evolve::vx")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,x_rhs(:,bhns_rk4_particle_tracer_step))
        call CCTK_VarIndex(vindex,"mhd_evolve::vy")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,y_rhs(:,bhns_rk4_particle_tracer_step))
        call CCTK_VarIndex(vindex,"mhd_evolve::vz")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,z_rhs(:,bhns_rk4_particle_tracer_step))
        
        ! recompute 1/u^t
        call CCTK_VarIndex(vindex,"mhd_evolve::temp18")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,taul_rhs(:,bhns_rk4_particle_tracer_step))


        ! Now that we have xyz_rhs at all 4 RK4 substeps, 
        !  we next construct k1, k2, k3, and k4, then 
        !  update tracer_xyz
        do which_particle=1,N_particles_to_trace
           ! x_rhs(...,1) == f(t_n,y_n)
           k1_x      = dT*   x_rhs(which_particle,1)
           k1_y      = dT*   y_rhs(which_particle,1)
           k1_z      = dT*   z_rhs(which_particle,1)
           k1_taul   = dT*taul_rhs(which_particle,1)
           
           ! x_rhs(...,2) == f(t_n + 1/2 dT , y_n + 1/2 k_1) == k2_x / dT
           k2_x      = dT*   x_rhs(which_particle,2)
           k2_y      = dT*   y_rhs(which_particle,2)
           k2_z      = dT*   z_rhs(which_particle,2)
           k2_taul   = dT*taul_rhs(which_particle,2)
           
           ! x_rhs(...,3) == f(t_n + 1/2 dT , y_n + 1/2 k_2) == k3_x / dT
           k3_x      = dT*   x_rhs(which_particle,3)
           k3_y      = dT*   y_rhs(which_particle,3)
           k3_z      = dT*   z_rhs(which_particle,3)
           k3_taul   = dT*taul_rhs(which_particle,3)
           
           ! x_rhs(...,4) == f(t_n + dT , y_n + k_3) == k4_x / dT
           k4_x      = dT*   x_rhs(which_particle,4)
           k4_y      = dT*   y_rhs(which_particle,4)
           k4_z      = dT*   z_rhs(which_particle,4)
           k4_taul   = dT*taul_rhs(which_particle,4)
           
           tracer_x(which_particle)    = tracer_x(which_particle)    + 1.D0/6.D0 * (k1_x    + 2.D0*k2_x + 2.D0*k3_x + k4_x)
           tracer_y(which_particle)    = tracer_y(which_particle)    + 1.D0/6.D0 * (k1_y    + 2.D0*k2_y + 2.D0*k3_y + k4_y)
           tracer_z(which_particle)    = tracer_z(which_particle)    + 1.D0/6.D0 * (k1_z    + 2.D0*k2_z + 2.D0*k3_z + k4_z)
           tracer_taul(which_particle) = tracer_taul(which_particle) + 1.D0/6.D0 * (k1_taul + 2.D0*k2_taul + 2.D0*k3_taul + k4_taul)
        end do
        
        
        bhns_rk4_particle_tracer_step = 1
        
        write(*,*) " Particle Tracer substep # ",bhns_rk4_particle_tracer_step
        
        ! Next do the first RK4 substep!
        do which_particle=1,N_particles_to_trace
           pointcoords(which_particle,1) =     tracer_x(which_particle)
           pointcoords(which_particle,2) =     tracer_y(which_particle)
           pointcoords(which_particle,3) =     tracer_z(which_particle)
           pointcoords(which_particle,4) =  tracer_taul(which_particle)
           
           
           pointcoords(which_particle,1) =    tracer_x(which_particle) + 0.5D0 * k2_x
           pointcoords(which_particle,2) =    tracer_y(which_particle) + 0.5D0 * k2_y
           pointcoords(which_particle,3) =    tracer_z(which_particle) + 0.5D0 * k2_z
           pointcoords(which_particle,4) = tracer_taul(which_particle) + 0.5D0 * k2_taul
        end do
        
        
        call CCTK_VarIndex(vindex,"mhd_evolve::vx")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,x_rhs(:,bhns_rk4_particle_tracer_step))
        call CCTK_VarIndex(vindex,"mhd_evolve::vy")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,y_rhs(:,bhns_rk4_particle_tracer_step))
        call CCTK_VarIndex(vindex,"mhd_evolve::vz")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,z_rhs(:,bhns_rk4_particle_tracer_step))

        call CCTK_VarIndex(vindex,"mhd_evolve::temp18")
        call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,taul_rhs(:,bhns_rk4_particle_tracer_step))
     end if
  end if
  !--------------------------------------------------------------------------------!
  !  FILE OUTPUT
  if((mod(cctk_iteration,particle_tracer_substep_every)==0 .and. bhns_rk4_particle_tracer_step==1) &
        .or. cctk_iteration.eq.0) then
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
     filename = 'bhns_circulation_ring.mon'
     
     ! First update the position of the particles 

     do which_particle=1,N_particles_to_trace
        pointcoords(which_particle,1) =    tracer_x(which_particle)
        pointcoords(which_particle,2) =    tracer_y(which_particle)
        pointcoords(which_particle,3) =    tracer_z(which_particle)
        pointcoords(which_particle,4) = tracer_taul(which_particle)
     end do
     

     ! Compute  u_x, u_y, and h
     call CCTK_VarIndex(vindex,"mhd_evolve::temp_16")
     call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,u_x_local)
     
     call CCTK_VarIndex(vindex,"mhd_evolve::temp_17")
     call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,u_y_local)
     
     
     call CCTK_VarIndex(vindex,"mhd_evolve::h")
     call interp_driver_carp(cctkGH,N_particles_to_trace,pointcoords,vindex,m1_minus_u_lower_zero)
 

     
     ! Then output everything to file.
     do which_particle=1,N_particles_to_trace
        
        num_cols = 8
        export_data(which_particle,:) = 0.D0
        
        data_headers(1) = '# Time'
        export_data(which_particle,1) = CCTK_TIME

        data_headers(2) = 'tau'
        export_data(which_particle,2) = tracer_taul(which_particle)
        
        data_headers(3) = 'x position'
        export_data(which_particle,3) = tracer_x(which_particle)
        
        data_headers(4) = 'y position'
        export_data(which_particle,4) = tracer_y(which_particle)
        
        data_headers(5) = 'u_x'
        export_data(which_particle,5) = u_x_local(which_particle)
        
        data_headers(6) = 'u_y'
        export_data(which_particle,6) = u_y_local(which_particle)
        
        data_headers(7) = 'h'
        export_data(which_particle,7) =  m1_minus_u_lower_zero(which_particle)

        data_headers(8) = 'circulation'
        export_data(which_particle,8) =  output_integral
        
     end do
     
     if(cctk_iteration.eq.0) then
        header_flag=1
     end if
     
     if(CCTK_MyProc(CCTKGH)==0) call output_multiline_data_to_file(filename,num_cols,data_headers,N_particles_to_trace,export_data,header_flag)
     
     write(*,*) "Particle tracer: Finished dumping data files at time",CCTK_TIME
     write(*,*) "Particle tracer: Previous dump was at time",CCTK_TIME-dT
     
     bhns_rk4_particle_tracer_step = 2
  end if
  write(*,*) "****************"
end if
end subroutine circulation_trace_particles
