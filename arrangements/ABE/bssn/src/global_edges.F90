!---------------------------------------------------------------------------------
! Figure out whether or not the local processor is on the edge of the global grid
!   Used with code that updates outer boundaries (e.g., update_boundary)
!   This routine is used so that outer boundaries are updated only when local
!      processor grid touches the outer boundary
!---------------------------------------------------------------------------------

subroutine global_edges(dx,dy,dz, &
     Xglobal_min,Yglobal_min,Zglobal_min,Xlocal_min,Ylocal_min,Zlocal_min, &
     Xglobal_max,Yglobal_max,Zglobal_max,Xlocal_max,Ylocal_max,Zlocal_max, &
     have_bdry_min,have_bdry_max,Symmetry)

  real*8  :: Xglobal_min,Yglobal_min,Zglobal_min,Xlocal_min,Ylocal_min,Zlocal_min
  real*8  :: Xglobal_max,Yglobal_max,Zglobal_max,Xlocal_max,Ylocal_max,Zlocal_max
  real*8  :: dx,dy,dz
  integer, dimension(3) :: have_bdry_min,have_bdry_max
  integer :: Symmetry
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)

  if(abs(Xglobal_max-Xlocal_max).lt.0.1*dx)then
     have_bdry_max(1) = 1
  else
     have_bdry_max(1) = 0
  end if

  if((abs(Yglobal_max-Ylocal_max).lt.0.1*dy).and. &
       Symmetry.ne.AXISYM)then
     have_bdry_max(2) = 1
  else
     have_bdry_max(2) = 0
  end if

  if(abs(Zglobal_max-Zlocal_max).lt.0.1*dz)then
     have_bdry_max(3) = 1
  else
     have_bdry_max(3) = 0
  end if

  if((abs(Xglobal_min-Xlocal_min).lt.0.1*dx).and. &
       (Symmetry==EQUATORIAL.or.Symmetry==NO_SYMM.or.Symmetry==PI_SYMM))then
     have_bdry_min(1) = 1
  else
     have_bdry_min(1) = 0
  end if
  
  if((abs(Yglobal_min-Ylocal_min).lt.0.1*dy).and. &
       (Symmetry==EQUATORIAL.or.Symmetry==NO_SYMM))then
     have_bdry_min(2) = 1
  else
     have_bdry_min(2) = 0
  end if

  if((abs(Zglobal_min-Zlocal_min).lt.0.1*dz).and.Symmetry==NO_SYMM)then
     have_bdry_min(3) = 1
  else
     have_bdry_min(3) = 0
  end if

  return
end subroutine global_edges
