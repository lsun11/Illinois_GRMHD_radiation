!---------------------------------------------------------------------
!
! radius
!
!---------------------------------------------------------------------

  real*8 function radius(x,y,z,Symmetry)
    implicit none
    real*8,  intent(in) :: x,y,z
    integer, intent(in) :: Symmetry
    integer, parameter  :: AXISYM = 4

    if (Symmetry == AXISYM) then
       radius = sqrt(x*x + z*z)
    else
       radius = sqrt(x*x + y*y + z*z)
    end if
  end function radius
