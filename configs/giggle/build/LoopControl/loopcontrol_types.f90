module loopcontrol_types
  ! Note: this type must correspond to the corresponding C type
  ! declared in loopcontrol.h
  type lc_statmap_t
     sequence
     integer*8 next
     ! Name
     integer*8 name
     integer*8 statset_list
  end type lc_statmap_t
  ! Note: this type must correspond to the corresponding C type
  ! declared in loopcontrol.h
  type lc_control_t
     sequence
     integer*8 statmap
     integer*8 statset
     integer*8 stattime
     ! Copy of arguments (useful for debugging)
     integer imin, jmin, kmin
     integer imax, jmax, kmax
     integer ilsh, jlsh, klsh
     ! Control settings for thread parallelism (useful for debugging)
     integer iiimin, jjjmin, kkkmin
     integer iiimax, jjjmax, kkkmax
     integer iiistep, jjjstep, kkkstep
     ! Control settings for current thread (useful for debugging)
     integer thread_num
     integer iii, jjj, kkk
     ! Control settings for tiling loop
     integer iimin, jjmin, kkmin
     integer iimax, jjmax, kkmax
     integer iistep, jjstep, kkstep
     ! Timing statistics
     double precision time_setup_begin, time_calc_begin
  end type lc_control_t
end module loopcontrol_types
