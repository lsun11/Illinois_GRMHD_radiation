MODULE nrutil_gw_extraction
        IMPLICIT NONE
        INTEGER, PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
        INTERFACE assert
                MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
        END INTERFACE
        INTERFACE arth
                MODULE PROCEDURE arth_d, arth_i
        END INTERFACE
CONTAINS
!BL
        SUBROUTINE assert1(n1,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1
        if (.not. n1) then
                write (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert1'
        end if
        END SUBROUTINE assert1
!BL
        SUBROUTINE assert2(n1,n2,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1,n2
        if (.not. (n1 .and. n2)) then
                write (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert2'
        end if
        END SUBROUTINE assert2
!BL
        SUBROUTINE assert3(n1,n2,n3,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1,n2,n3
        if (.not. (n1 .and. n2 .and. n3)) then
                write (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert3'
        end if
        END SUBROUTINE assert3
!BL
        SUBROUTINE assert4(n1,n2,n3,n4,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1,n2,n3,n4
        if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
                write (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert4'
        end if
        END SUBROUTINE assert4
!BL
        SUBROUTINE assert_v(n,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, DIMENSION(:), INTENT(IN) :: n
        if (.not. all(n)) then
                write (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert_v'
        end if
        END SUBROUTINE assert_v
!BL
        FUNCTION assert_eqn(nn,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER, DIMENSION(:), INTENT(IN) :: nn
        INTEGER :: assert_eqn
        if (all(nn(2:) == nn(1))) then
                assert_eqn=nn(1)
        else
                write (*,*) 'nrerror: an assert_eq failed with this tag:', &
                        string
                STOP 'program terminated by assert_eqn'
        end if
        END FUNCTION assert_eqn
!BL
        FUNCTION arth_d(first,increment,n)
        REAL*8, INTENT(IN) :: first,increment
        INTEGER, INTENT(IN) :: n
        REAL*8, DIMENSION(n) :: arth_d
        INTEGER :: k,k2
        REAL*8 :: temp
        if (n > 0) arth_d(1)=first
        if (n <= NPAR_ARTH) then
                do k=2,n
                        arth_d(k)=arth_d(k-1)+increment
                end do
        else
                do k=2,NPAR2_ARTH
                        arth_d(k)=arth_d(k-1)+increment
                end do
                temp=increment*NPAR2_ARTH
                k=NPAR2_ARTH
                do
                        if (k >= n) exit
                        k2=k+k
                        arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
                        temp=temp+temp
                        k=k2
                end do
        end if
        END FUNCTION arth_d
!BL
        FUNCTION arth_i(first,increment,n)
        INTEGER, INTENT(IN) :: first,increment,n
        INTEGER, DIMENSION(n) :: arth_i
        INTEGER :: k,k2,temp
        if (n > 0) arth_i(1)=first
        if (n <= NPAR_ARTH) then
                do k=2,n
                        arth_i(k)=arth_i(k-1)+increment
                end do
        else
                do k=2,NPAR2_ARTH
                        arth_i(k)=arth_i(k-1)+increment
                end do
                temp=increment*NPAR2_ARTH
                k=NPAR2_ARTH
                do
                        if (k >= n) exit
                        k2=k+k
                        arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
                        temp=temp+temp
                        k=k2
                end do
        end if
        END FUNCTION arth_i
!BL
END MODULE nrutil_gw_extraction
