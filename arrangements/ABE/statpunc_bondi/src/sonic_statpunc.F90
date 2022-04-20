SUBROUTINE sonic_statpunc(M,a2_inf,gamma_th,a2_s,Rs_s,u2_s)
  IMPLICIT NONE
  real*8 :: M,a2_inf,gamma_th,a2_s,Rs_s,u2_s,a2_s_prev
  a2_s = a2_inf
  a2_s_prev=a2_s
  
  if (abs(gamma_th - 5.d0/3.d0) .gt. 1.d-4) then
     do
        a2_s = (2.d0*(gamma_th-1.d0)*a2_inf+&
             (1.d0-6.d0*(gamma_th-1.d0))*a2_s*a2_s+&
             3.d0*a2_s*a2_s*a2_s-a2_inf*a2_inf)/(2.d0*(gamma_th-1.d0)-&
             3.d0*(gamma_th-1.d0)**2)
        if (abs((a2_s-a2_s_prev)/a2_s_prev) .lt. 1.d-12) exit
        a2_s_prev = a2_s
     end do
  else
     do 
        a2_s = sqrt(a2_s*a2_s*a2_s+4.d0/9.d0*a2_inf-1.d0/3.d0*a2_inf*a2_inf)
        if (abs((a2_s-a2_s_prev)/a2_s_prev) .lt. 1.d-14) exit
        a2_s_prev = a2_s
     end do
  endif
 
  !inward 4-vel squared at sonic radius
  u2_s = a2_s/(1.d0+3.d0*a2_s)

  !sonic radius in schwarzchild coords
  Rs_s = (1.d0+3.d0*a2_s)/(2.d0*a2_s)*M 
 
END SUBROUTINE sonic_statpunc
