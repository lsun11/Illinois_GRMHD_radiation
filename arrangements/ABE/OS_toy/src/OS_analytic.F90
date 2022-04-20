SUBROUTINE OS_analytic(tau,rho,Q,R_edge,M)
use nr
use nrtype
IMPLICIT NONE 
interface 
   FUNCTION os_func(eta,tau,R_edge,M)
     USE nrtype
     IMPLICIT NONE 
     REAL(DP), INTENT(IN) :: eta,tau,R_edge,M
     REAL(DP) :: os_func
   END FUNCTION os_func
end interface
real(DP) :: x1,x2,xacc
real(DP) :: eta,Q,rho,tau,R_edge,M
x1 = -1.d-6
x2 = pi
xacc = 1.d-6
write(*,*) "Begin OS_analytic. tau_center_OS=", tau
write(*,*) "about to calculate eta"
eta = abs(rtbis(os_func,x1,x2,xacc,tau,R_edge,M))
write(*,*) "just calculated eta, eta=", eta
Q=0.5d0*(1.d0+cos(eta))
rho=M/(4.d0/3.d0*pi*R_edge**3)*Q**(-3)

write(*,*) "done OS_analytic, rho, Q=", rho,Q
END SUBROUTINE OS_analytic

FUNCTION os_func(eta,tau,R_edge,M)
USE nrtype
IMPLICIT NONE
REAL(DP), INTENT(IN) :: eta,tau,R_edge,M
REAL(DP) :: os_func
os_func = sqrt(R_edge**3*M*M/8.d0)*(eta + sin(eta)) - tau
END FUNCTION os_func   
