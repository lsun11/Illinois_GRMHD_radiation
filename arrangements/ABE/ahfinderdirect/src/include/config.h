/* config.h -- compile-time configuration for AHFinderDirect */
/* $Header: /numrelcvs/AEIThorns/AHFinderDirect/src/include/config.h,v 1.16 2003/06/04 14:04:10 jthorn Exp $ */

#ifndef AHFINDERDIRECT__CONFIG_H
#define AHFINDERDIRECT__CONFIG_H

/*
 * This header file is #included by C++ and C files, and also by
 * Fortran 77 files which are run through a C-like preprocessor.  Alas,
 * due to a Cactus bug on some systems, this means that all preprocessor
 * commands must start in column 1. :( :( :(
 */

/*
 * prerequisites:
 *	"cctk.h" or "fake_cctk.h"
 */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/*
 * definitions for C/C++ only
 */

#ifdef CCODE
typedef CCTK_REAL fp;

/* Fortran 'integer' type */
typedef CCTK_INT integer;
#endif	/* CCODE */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/*
 * definitions for C++, C, and Fortran
 */

/* CCTK_VWarn() "warning level" for fatal errors (terminate execution) */
/* note low-level software uses error_exit() instead of CCTK_VWarn(), */
/* with exit codes defined in "stdc.h" */
#define FATAL_ERROR	(-1)

/******************************************************************************/

/*
 * ***** definitions for integer and floating point datatypes *****
 */

/*
 * FIXME: this assumes fp == C 'double'
 * CCTK_REAL_PRECISION_{4,8,16} are helpful, but not quite enough)
 */
#undef  FP_IS_FLOAT
#define FP_IS_DOUBLE

#if   defined(FP_IS_FLOAT)
#define FP_SCANF_FORMAT "%f"
#define FP_EPSILON	FLT_EPSILON	/* from <float.h> */
#elif defined(FP_IS_DOUBLE)
#define FP_SCANF_FORMAT "%lf"
#define FP_EPSILON	DBL_EPSILON	/* from <float.h> */
#else
#error "don't know fp datatype!"
#endif

/*
 * FIXME: this should be set somehow at configure time
 */
#define FORTRAN_INTEGER_IS_INT
#undef  FORTRAN_INTEGER_IS_LONG

/******************************************************************************/

/*
 * The angular finite differencing in our multipatch system can be
 * either 2nd order or 4th order, configurable here.
 */
#ifdef FINITE_DIFF_ORDER
/* someone (probably a Makefile) has already defined this */
/* ==> leave it alone */
#else
#define FINITE_DIFF_ORDER	4
#endif

/******************************************************************************/

/*
 * What types of Jacobian matrix storage and linear solvers do we want to
 * compile in support for?  Note that each of these requires linking with
 * the corresponding linear-solver library; see "../elliptic/Jacobian.hh"
 * for details on the storage formats, or "../make.configuration.defn"
 * (if it exists) for details on the libraries.
 */

/* store as (Fortran) dense matrix, solve with LAPACK */
#undef  HAVE_DENSE_JACOBIAN__LAPACK

/* store as row-oriented sparse matrix, solve with ILUCG */
#define HAVE_ROW_SPARSE_JACOBIAN__ILUCG

/* store as row-oriented sparse matrix, solve with UMFPACK */
#define HAVE_ROW_SPARSE_JACOBIAN__UMFPACK

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/*
 * The following #ifdefs and #defines shouldn't need changing; they
 * decode the previous set of Jacobian storage format and linear solver
 * #defines to determine which storage formats we are using.
 */
#ifdef HAVE_DENSE_JACOBIAN__LAPACK
#define HAVE_DENSE_JACOBIAN
#endif

#ifdef HAVE_ROW_SPARSE_JACOBIAN__ILUCG
#define HAVE_ROW_SPARSE_JACOBIAN
#endif

#ifdef HAVE_ROW_SPARSE_JACOBIAN__UMFPACK
#define HAVE_ROW_SPARSE_JACOBIAN
#endif

#endif	/* AHFINDERDIRECT__CONFIG_H */
