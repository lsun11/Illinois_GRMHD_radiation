FUNCTION funcv_bisect_step(xlo,xhi,fval,outputdiff,lo_or_hi) 
  IMPLICIT NONE 
  REAL*8, INTENT(INOUT) :: xlo,xhi,fval
  REAL*8 :: outputdiff
  REAL*8 :: ftmplo,ftmphi
  REAL*8 :: M,R
  INTEGER :: funcv_bisect_step,lo_or_hi
  M = 1.D0

  if(lo_or_hi == 1) then
     if(xlo .le. 3.D0*M*0.5) xlo = 3.D0*M*0.5+1.D-5
     R = xlo
     
     ftmplo = fval - (2.D0*R + M + sqrt(4.D0*R*R + 4.D0*M*R + 3.D0*M*M))*0.25D0 * &
          (((4.D0+3.D0*sqrt(2.D0))*(2.D0*R-3.D0*M)) / &
          (8.D0*R+6.D0*M+3.D0*sqrt(8.D0*R*R+8.D0*M*R+6.D0*M*M)))**(1.D0/sqrt(2.D0))
     ftmphi = outputdiff
  else if(lo_or_hi == 0) then
     R = xhi
     ftmphi = fval - (2.D0*R + M + sqrt(4.D0*R*R + 4.D0*M*R + 3.D0*M*M))*0.25D0 * &
          (((4.D0+3.D0*sqrt(2.D0))*(2.D0*R-3.D0*M)) / &
          (8.D0*R+6.D0*M+3.D0*sqrt(8.D0*R*R+8.D0*M*R+6.D0*M*M)))**(1.D0/sqrt(2.D0))
     ftmplo = outputdiff
  else
     if(xlo .le. 3.D0*M*0.5) xlo = 3.D0*M*0.5+1.D-5
     R = xlo
     
     ftmplo = fval - (2.D0*R + M + sqrt(4.D0*R*R + 4.D0*M*R + 3.D0*M*M))*0.25D0 * &
          (((4.D0+3.D0*sqrt(2.D0))*(2.D0*R-3.D0*M)) / &
          (8.D0*R+6.D0*M+3.D0*sqrt(8.D0*R*R+8.D0*M*R+6.D0*M*M)))**(1.D0/sqrt(2.D0))

     R = xhi
     ftmphi = fval - (2.D0*R + M + sqrt(4.D0*R*R + 4.D0*M*R + 3.D0*M*M))*0.25D0 * &
          (((4.D0+3.D0*sqrt(2.D0))*(2.D0*R-3.D0*M)) / &
          (8.D0*R+6.D0*M+3.D0*sqrt(8.D0*R*R+8.D0*M*R+6.D0*M*M)))**(1.D0/sqrt(2.D0))
  end if
     
  if(abs(ftmplo) .le. abs(ftmphi)) then
     outputdiff = abs(ftmplo)
     lo_or_hi = 0
     funcv_bisect_step = 0
  else
     outputdiff = abs(ftmphi)
     lo_or_hi = 1
     funcv_bisect_step = 1
  end if

END FUNCTION funcv_bisect_step
