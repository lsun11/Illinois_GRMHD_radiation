FUNCTION get_a2(rho0_derivs_s,Rs_s,a2_s,M,Rs,Mdot,K_poly,Gamma,G29_rhs)
  USE nr
  USE nrtype
  IMPLICIT NONE
  
  INTERFACE
     FUNCTION a2_function(a2,M,Rs,Mdot,K_poly,Gamma,G29_rhs)
       USE nrtype
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: a2,M,Rs,Mdot,K_poly,Gamma,G29_rhs
       REAL(DP) :: a2_function
     END FUNCTION a2_function
  END INTERFACE
  real(DP), INTENT(IN) :: Rs_s,a2_s,M,Rs,Mdot,K_poly,Gamma,G29_rhs,rho0_derivs_s
  real(DP)             :: x1,x2,xacc,a2_low,a2_high,tiny,delta
  real(DP)             :: get_a2
  REAL(DP), DIMENSION(:), POINTER :: xb1,xb2
  INTEGER(I4B) :: num,nb,found,i
  
  !set xacc
  xacc = 1.d-12
  tiny = 1.d-2
  !initialize 
  found=0

  !there are problems with finding a solution close to the sonic point, which is why I extrapolate using the derivative.  This could probably be improved
  if (Rs .gt. Rs_s+tiny) then
   num=100
     do while (found .eq. 0)
        num=num*10
        delta=a2_s/num
        
        do i=1,num
           x1=a2_s-i*delta
           x2=a2_s-(i-1)*delta
           if (a2_function(x1,M,Rs,Mdot,K_poly,Gamma,G29_rhs)*a2_function(x2,M,Rs,Mdot,K_poly,Gamma,G29_rhs) .lt. 0) then
              found=1
              exit
           endif
        end do
     end do
     get_a2 = rtbis(a2_function,x1,x2,xacc,M,Rs,Mdot,K_poly,Gamma,G29_rhs)
  endif
  
  if (Rs .lt. Rs_s-tiny) then
     num=100
     do while (found .eq. 0)
        num=num*10
        delta=abs(0.99*1.d0/3.d0-a2_s)/num

        do i=1,num
           x1=a2_s+(i-1)*delta
           x2=a2_s+(i)*delta
           if (a2_function(x1,M,Rs,Mdot,K_poly,Gamma,G29_rhs)*a2_function(x2,M,Rs,Mdot,K_poly,Gamma,G29_rhs) .lt. 0) then
              found=1
              exit
           endif
        end do
     end do
     get_a2 = rtbis(a2_function,x1,x2,xacc,M,Rs,Mdot,K_poly,Gamma,G29_rhs)
  endif
  
  if (Rs .gt. Rs_s-tiny .and. Rs .lt. Rs_s+tiny) then
     get_a2=a2_s+rho0_derivs_s*(Rs-Rs_s)
  endif
  
END FUNCTION get_a2

FUNCTION a2_function(a2,M,Rs,Mdot,K_poly,Gamma,G29_rhs)
  USE nrtype
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: a2,M,Rs,Mdot,K_poly,Gamma,G29_rhs
  REAL(DP) :: a2_function

  a2_function = (1.d0 - 2.d0*M/Rs+(Mdot/(4*PI_D*(1.d0/K_poly*a2/(Gamma-Gamma/(Gamma-1)*a2))**(1.d0/(Gamma-1.d0))*Rs*Rs))**2)*(1.d0+a2/(Gamma-1.d0-a2))**2-G29_rhs
END FUNCTION a2_function


