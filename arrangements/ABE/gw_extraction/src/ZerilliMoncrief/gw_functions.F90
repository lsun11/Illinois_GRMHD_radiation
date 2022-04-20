module gw_functions
contains
  FUNCTION plgndr(l,m,x)
    !USE nrtype; 
    USE nrutil_gw_extraction, ONLY : arth ,assert
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: l,m
    REAL*8, INTENT(IN) :: x
    REAL*8 :: plgndr
    INTEGER :: i,ll
    REAL*8 :: fact,pll,pmm,pmmp1,somx2
    !call assert(m >= 0, m <= l, abs(x) <= 1.D0, 'plgndr args')
    if (m .lt. 0) then 
       write(*,*) 'm = ',m
       write(*,*) 'Error! m<0 in plgndr.'
       stop
    end if
    if (abs(m) .gt. l) then 
       write(*,*) 'l, m = ',l,m
       write(*,*) 'Error! |m| > l in plgndr.'
       stop
    end if
    if (abs(x) .gt. 1.d0) then 
       write(*,*) 'x = ',x
       write(*,*) 'Error! |x| > 1 in plgndr.'
       stop
    end if
    pmm=1.D0
    if (m > 0) then
       somx2=sqrt((1.D0-x)*(1.D0+x))
!following lines added by Zach
       fact = 1.D0
       do i=1,m
          pmm = pmm*(-fact*somx2)
          fact = fact + 2.D0
       end do
       
!following lines removed by Zach
!       pmm=product(arth(1.D0,2.0D0,m))*somx2**m
!       if (mod(m,2) == 1) pmm=-pmm
    end if
    if (l == m) then
       plgndr=pmm
    else
       pmmp1=x*(2*m+1)*pmm
       if (l == m+1) then
          plgndr=pmmp1
       else
          do ll=m+2,l
             pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
             pmm=pmmp1
             pmmp1=pll
          end do
          plgndr=pll
       end if
    end if
  END FUNCTION plgndr
  
  FUNCTION factorial(n)
    integer i
    real*8 factorial
    
    factorial = 1.D0
    do i=2,n
       factorial = factorial * i
    end do
    
  END FUNCTION factorial

end module gw_functions
