c lapack_wrapper.f -- wrapper routines for LAPACK [sd]gecon()
c $Header: /numrelcvs/AEIThorns/AHFinderDirect/src/elliptic/lapack_wrapper.F77,v 1.5 2003/03/24 14:50:12 jthorn Exp $
c
c These subroutines are wrappers around the LAPACK [sd]gecon() subroutines.
c These subroutines take only integer/real/double precision arguments,
c avoiding problems with C/C++ --> Fortran passing of the character string
c arguments used by [sd]gecon().
c
c Arguments:
c norm_int = (in) 0 ==> infinity-norm
c                 1 ==> 1-norm
c
c
c Note that the Compaq f90 compiler complains about empty (or all-comment)
c files, so we still have to define an empty subroutine even if it is never
c called.  :( Oh well, memory is cheap...
c
        subroutine sgecon_wrapper(norm_int,
     $                                  N, A, LDA, anorm, rcond,
     $                                  WORK, IWORK, info)
        integer norm_int
        integer N, LDA
        real A(LDA,N)
        real anorm, rcond
        real WORK(*)
        integer iwork(*)
        integer info
        return
        end
        subroutine dgecon_wrapper(norm_int,
     $                                  N, A, LDA, anorm, rcond,
     $                                  WORK, IWORK, info)
        integer norm_int
        integer N, LDA
        double precision A(LDA,N)
        double precision anorm, rcond
        double precision WORK(*)
        integer iwork(*)
        integer info
        return
        end
