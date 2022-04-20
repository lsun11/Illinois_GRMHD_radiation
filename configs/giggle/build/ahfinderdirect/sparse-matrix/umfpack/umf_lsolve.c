/* ========================================================================== */
/* === UMF_lsolve =========================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*  solves Lx = b, where L is the lower triangular factor of a matrix */
/*  B is overwritten with the solution X. */
/*  Returns the floating point operation count */

#include "umf_internal.h"

GLOBAL double UMF_lsolve
(
    NumericType *Numeric,
    Entry X [ ],		/* b on input, solution x on output */
    Int Pattern [ ]		/* a work array of size n */
)
{
    Int k, deg, *ip, j, row, *Lpos, *Lilen, *Lip, llen, lp, newLchain,
	pos, npiv ;
    Entry *xp, xk ;

    /* ---------------------------------------------------------------------- */

    if (Numeric->n_row != Numeric->n_col) return (0.) ;
    npiv = Numeric->npiv ;
    Lpos = Numeric->Lpos ;
    Lilen = Numeric->Lilen ;
    Lip = Numeric->Lip ;
    deg = 0 ;

#ifndef NDEBUG
    DEBUG4 (("Lsolve start:\n")) ;
    for (j = 0 ; j < Numeric->n_row ; j++)
    {
	DEBUG4 (("Lsolve start "ID": ", j)) ;
	EDEBUG4 (X [j]) ;
	DEBUG4 (("\n")) ;
    }
#endif

    for (k = 0 ; k < npiv ; k++)
    {

	/* ------------------------------------------------------------------ */
	/* make column of L in Pattern [0..deg-1] */
	/* ------------------------------------------------------------------ */

	lp = Lip [k] ;
	newLchain = (lp < 0) ;
	if (newLchain)
	{
	    lp = -lp ;
	    deg = 0 ;
	    DEBUG4 (("start of chain for column of L\n")) ;
	}

	/* remove pivot row */
	pos = Lpos [k] ;
	if (pos != EMPTY)
	{
	    DEBUG4 (("  k "ID" removing row "ID" at position "ID"\n",
	    k, Pattern [pos], pos)) ;
	    ASSERT (!newLchain) ;
	    ASSERT (deg > 0) ;
	    ASSERT (pos >= 0 && pos < deg) ;
	    ASSERT (Pattern [pos] == k) ;
	    Pattern [pos] = Pattern [--deg] ;
	}

	/* concatenate the pattern */
	ip = (Int *) (Numeric->Memory + lp) ;
	llen = Lilen [k] ;
	for (j = 0 ; j < llen ; j++)
	{
	    row = *ip++ ;
	    DEBUG4 (("  row "ID"  k "ID"\n", row, k)) ;
	    ASSERT (row > k) ;
	    Pattern [deg++] = row ;
	}

	/* ------------------------------------------------------------------ */
	/* use column k of L */
	/* ------------------------------------------------------------------ */

	xk = X [k] ;
	if (IS_NONZERO (xk))
	{
	    xp = (Entry *) (Numeric->Memory + lp + UNITS (Int, llen)) ;
	    for (j = 0 ; j < deg ; j++)
	    {
		DEBUG4 (("  row "ID"  k "ID" value", Pattern [j], k)) ;
		EDEBUG4 (*xp) ;
		DEBUG4 (("\n")) ;
		/* X [Pattern [j]] -= xk * (*xp) ; */
		MULT_SUB (X [Pattern [j]], xk, *xp) ;
		xp++ ;
	    }
	}
    }

#ifndef NDEBUG
    for (j = 0 ; j < Numeric->n_row ; j++)
    {
	DEBUG4 (("Lsolve done "ID": ", j)) ;
	EDEBUG4 (X [j]) ;
	DEBUG4 (("\n")) ;
    }
    DEBUG4 (("Lsolve done.\n")) ;
#endif

    return (MULTSUB_FLOPS * ((double) Numeric->lnz)) ;
}
