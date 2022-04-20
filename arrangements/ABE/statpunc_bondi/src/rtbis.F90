	FUNCTION rtbis(func,x1,x2,xacc,M,r,Mdot,K,Gamma,G29_rhs)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: x1,x2,xacc
	REAL(DP) :: rtbis
        REAL(DP) :: M,r,Mdot,K,Gamma,G29_rhs
        INTERFACE
           FUNCTION func(a2,M,r,Mdot,K,Gamma,G29_rhs)
             USE nrtype
             IMPLICIT NONE
             REAL(DP), INTENT(IN) :: a2,M,r,Mdot,K,Gamma,G29_rhs
             REAL(DP) :: func
           END FUNCTION func
        END INTERFACE
        INTEGER(I4B), PARAMETER :: MAXIT=40
	INTEGER(I4B) :: j
	REAL(DP) :: dx,f,fmid,xmid
	fmid=func(x2,M,r,Mdot,K,Gamma,G29_rhs)
	f=func(x1,M,r,Mdot,K,Gamma,G29_rhs)
	if (f*fmid >= 0.0) call nrerror('rtbis: root must be bracketed')
	if (f < 0.0) then
		rtbis=x1
		dx=x2-x1
	else
		rtbis=x2
		dx=x1-x2
	end if
	do j=1,MAXIT
		dx=dx*0.5_sp
		xmid=rtbis+dx
		fmid=func(xmid,M,r,Mdot,K,Gamma,G29_rhs)
		if (fmid <= 0.0) rtbis=xmid
		if (abs(dx) < xacc .or. fmid == 0.0) RETURN
	end do
	call nrerror('rtbis: too many bisections')
	END FUNCTION rtbis
