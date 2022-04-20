FUNCTION bessk2(x)
  IMPLICIT NONE
  INTERFACE
     FUNCTION bessk0(x)
       IMPLICIT NONE
       INTERFACE
          FUNCTION bessi0(x)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: x
            REAL*8 :: bessi0
            REAL*8 :: ax
            REAL*8 :: y
          END FUNCTION bessi0
       END INTERFACE
       REAL*8, INTENT(IN) :: x
       REAL*8 :: bessk0
       REAL*8 :: y
     END FUNCTION bessk0
  END INTERFACE
  INTERFACE
     FUNCTION bessk1(x)
       IMPLICIT NONE
       INTERFACE
          FUNCTION bessi1(x)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: x
            REAL*8 :: bessi1
            REAL*8 :: ax
          END FUNCTION bessi1
       END INTERFACE
       REAL*8, INTENT(IN) :: x
       REAL*8 :: bessk1
       REAL*8 :: y
     END FUNCTION bessk1
  END INTERFACE
  REAL*8, INTENT(IN) :: x
  REAL*8 :: bessk2
  REAL*8 :: bk,bkm,bkp,tox
  tox=2.d0/x
  bkm=bessk0(x)
  bk=bessk1(x)

  bkp=bkm+tox*bk
  bkm=bk
  bk=bkp

  bessk2=bk
END FUNCTION bessk2


FUNCTION bessk0(x)
  IMPLICIT NONE
  INTERFACE
     FUNCTION bessi0(x)
       IMPLICIT NONE
       REAL*8, INTENT(IN) :: x
       REAL*8 :: bessi0
       REAL*8 :: ax
       REAL*8 :: y
     END FUNCTION bessi0
  END INTERFACE
  REAL*8, INTENT(IN) :: x
  REAL*8 :: bessk0
  REAL*8 :: y
  if (x <= 2.0) then
     y=x*x/4.d0
     bessk0=(-log(x/2.d0)*bessi0(x))+(-0.57721566d0+y*(0.42278420d0 &
                        +y*(0.23069756d0+y*(0.3488590d-1+y*(0.262698d-2 &
                        +y*(0.10750d-3+y*0.74d-5))))))
  else
     y=(2.d0/x)
     bessk0=(exp(-x)/sqrt(x))*(1.25331414d0+y*(-0.7832358d-1 &
                        +y*(0.2189568d-1+y*(-0.1062446d-1+y*(0.587872d-2 &
                        +y*(-0.251540d-2+y*0.53208d-3))))))
  end if
END FUNCTION bessk0


FUNCTION bessk1(x)
  IMPLICIT NONE
  INTERFACE
     FUNCTION bessi1(x)
       IMPLICIT NONE
       REAL*8, INTENT(IN) :: x
       REAL*8 :: bessi1
       REAL*8 :: ax
       REAL*8 :: y
     END FUNCTION bessi1
  END INTERFACE
  REAL*8, INTENT(IN) :: x
  REAL*8 :: bessk1
  REAL*8 :: y
  if (x <= 2.0) then
     y=x*x/4.d0
     bessk1=(log(x/2.d0)*bessi1(x))+(1.d0/x)*(1.d0+y*(0.15443144d0 &
                        +y*(-0.67278579d0+y*(-0.18156897d0+y*(-0.1919402d-1 &
                        +y*(-0.110404d-2+y*(-0.4686d-4)))))))
  else
     y=2.d0/x
     bessk1=(exp(-x)/sqrt(x))*(1.25331414d0+y*(0.23498619d0 &
                        +y*(-0.3655620d-1+y*(0.1504268d-1+y*(-0.780353d-2 &
                        +y*(0.325614d-2+y*(-0.68245d-3)))))))
  end if
END FUNCTION bessk1


FUNCTION bessi0(x)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: x
  REAL*8 :: bessi0
  REAL*8 :: ax
  REAL*8 :: y
  ax=abs(x)
  if (ax < 3.75) then
     y=x/3.75
     y=y*y
     bessi0= 1.0+y*(3.5156229d0+y*(3.0899424d0+y*(1.2067492d0 &
          +y*(0.2659732d0+y*(0.360768d-1+y*0.45813d-2)))))
 
  else
     y=3.75d0/ax
     bessi0=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592d-1 &
                        +y*(0.225319d-2+y*(-0.157565d-2+y*(0.916281d-2 &
                        +y*(-0.2057706d-1+y*(0.2635537d-1+y*(-0.1647633d-1 &
                        +y*0.392377d-2))))))))
  end if
END FUNCTION bessi0

FUNCTION bessi1(x)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: x
  REAL*8 :: bessi1
  REAL*8 :: ax
  REAL*8 :: y
  ax=abs(x)
  if (ax < 3.75) then
     y=x/3.75
     y=y*y;
     bessi1=ax*(0.5+y*(0.87890594d0+y*(0.51498869d0+y*(0.15084934d0 &
          +y*(0.2658733d-1+y*(0.301532d-2+y*0.32411d-3))))))
  else
     y=3.75/ax
     bessi1=0.2282967d-1+y*(-0.2895312d-1+y*(0.1787654d-1-y*0.420059d-2))
     bessi1=0.39894228d0+y*(-0.3988024d-1+y*(-0.362018d-2 &
          +y*(0.163801d-2+y*(-0.1031555d-1+y*bessi1))))
     bessi1 = bessi1*(exp(ax)/sqrt(ax))
  end if
  if (x < 0.0) bessi1=-bessi1
END FUNCTION bessi1
