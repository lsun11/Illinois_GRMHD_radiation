SUBROUTINE bremsstrahlung(theta,F_ei,F_ee)
  IMPLICIT NONE
  real*8 :: theta
  real*8 :: F_ei,F_ee
  real*8 :: pi
  pi=3.14159265358979323846
  ! (Narayan and Yi, 1995)
  if (theta < 1) then
     F_ei = 4.d0 * sqrt(2.d0*theta/pi**3) &
          *(1.d0+1.781d0*(theta)**(1.34))
  else
     F_ei =  (9.d0*theta)/(2.d0*pi)*(log(1.123*theta+0.48)+1.5)
  endif
  if (theta < 1) then
     F_ee = 20.d0/(9.d0*sqrt(pi))*(44.d0-3.d0*pi**2)*&
          (theta)**(1.5d0) &
          *(1.d0+1.1d0*theta+(theta)**2-1.25d0*(theta)**2.5d0)
  else
     F_ee = 24.d0*theta*(log(1.123*theta)+1.28d0)
  endif
END SUBROUTINE bremsstrahlung
