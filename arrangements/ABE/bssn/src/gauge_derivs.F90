!-----------------------------------------------------
! time derivatives of lapse and shift
! (for when we restart with hyperbolic gauges)
!-----------------------------------------------------
subroutine gauge_derivs(ex,dT,lapse_old,lapse_new,dtlapse, &
     betax_old,betax_new,dtbetax,betay_old,betay_new,dtbetay, &
     betaz_old,betaz_new,dtbetaz)
  implicit none
  integer, dimension(3)                :: ex
  real*8  :: dT
  real*8, dimension(ex(1),ex(2),ex(3)) :: lapse_old,lapse_new,dtlapse
  real*8, dimension(ex(1),ex(2),ex(3)) :: betax_old,betax_new,dtbetax
  real*8, dimension(ex(1),ex(2),ex(3)) :: betay_old,betay_new,dtbetay
  real*8, dimension(ex(1),ex(2),ex(3)) :: betaz_old,betaz_new,dtbetaz
  dtlapse = (lapse_new - lapse_old)/dT
  dtbetax = (betax_new - betax_old)/dT
  dtbetay = (betay_new - betay_old)/dT
  dtbetaz = (betaz_new - betaz_old)/dT
end subroutine gauge_derivs
