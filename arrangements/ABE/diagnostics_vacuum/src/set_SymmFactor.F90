!-------------------------------------------------------------------------------------------
! Set the factor we should multiply our integrals by.  E.g., SymmFactor=2 for equatorial symmetry.
!--------------------------------------------------------------------------------------------
subroutine set_SymmFactor(Symmetry,SymmFactor)
  implicit none
!
! I/O parameters:
!
  integer                                     :: Symmetry
  real*8                                      :: SymmFactor

! Other variables:
!
  integer :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  integer                            :: i,j,k
  real*8                             :: ZERO, ONE, SIX
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter( ZERO = 0.D0 )

  SymmFactor = 1.D0

  if(Symmetry == AXISYM) then
     SymmFactor = 1.D0
     return

  else if(Symmetry == OCTANT) then
     SymmFactor = 8.D0
     return
     
  else if(Symmetry == EQUATORIAL) then
     SymmFactor = 2.D0
     return
  end if
end subroutine set_SymmFactor

