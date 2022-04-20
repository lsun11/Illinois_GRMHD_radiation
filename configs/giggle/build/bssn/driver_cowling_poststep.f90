!-----------------------------------------------------------------------------
! Cowling approximation: set rhss == 0
!-----------------------------------------------------------------------------
!
subroutine cowling_poststep()
  implicit none
  ! Do nothing here.
  ! Just want Carpet to do prolongation for BSSN_refbd
end subroutine cowling_poststep
