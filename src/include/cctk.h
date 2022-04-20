/*@@
   @header    cctk.h
   @date      Tue Jan 26 17:29:34 1999
   @author    Tom Goodale
   @desc
              Main include file for the CCTK.
              All thorns should include this.
   @enddesc
   @version   $Header$
 @@*/

#ifndef _CCTK_H_
#define _CCTK_H_ 1

/* Grab the main configuration info. */
#include "cctk_Config.h"

/* Define which thorn the file is in */
#include "definethisthorn.h"

/* Include the constants */
#include "cctk_Constants.h"


/* Define some stuff */

#ifdef FCODE

#include "cctk_Faces.h"
#include "cctk_Interp.h"
#include "cctk_Loop.h"
#include "cctk_WarnLevel.h"

#define CCTK_PRINTSEPARATOR\
  print '("--------------------------------------------------------------------------------")'

#define _CCTK_FARGUMENTS \
cctk_dim,cctk_gsh,cctk_lsh,cctk_lbnd,cctk_ubnd,cctk_lssh,cctk_from,cctk_to,\
cctk_bbox,cctk_delta_time,cctk_time,cctk_delta_space,cctk_origin_space,\
cctk_levfac,cctk_levoff,cctk_levoffdenom,cctk_timefac,cctk_convlevel,\
cctk_convfac,cctk_nghostzones,cctk_iteration,cctkGH

#define _DECLARE_CCTK_ARGUMENTS _DECLARE_CCTK_FARGUMENTS
#define _DECLARE_CCTK_FARGUMENTS INTEGER cctk_dim&&\
                           INTEGER cctk_gsh(cctk_dim),cctk_lsh(cctk_dim)&&\
                           INTEGER cctk_lbnd(cctk_dim),cctk_ubnd(cctk_dim)&&\
                           INTEGER cctk_lssh(CCTK_NSTAGGER*cctk_dim)&&\
                           INTEGER cctk_from(cctk_dim),cctk_to(cctk_dim)&&\
                           INTEGER cctk_bbox(2*cctk_dim)&&\
                           CCTK_REAL cctk_delta_time, cctk_time&&\
                           CCTK_REAL cctk_delta_space(cctk_dim)&&\
                           CCTK_REAL cctk_origin_space(cctk_dim)&&\
                           INTEGER cctk_levfac(cctk_dim)&&\
                           INTEGER cctk_levoff(cctk_dim)&&\
                           INTEGER cctk_levoffdenom(cctk_dim)&&\
                           INTEGER cctk_timefac&&\
                           INTEGER cctk_convlevel&&\
                           INTEGER cctk_convfac&&\
                           INTEGER cctk_nghostzones(cctk_dim)&&\
                           INTEGER cctk_iteration&&\
                           CCTK_POINTER cctkGH&&

#define CCTK_WARN(a,b) CCTK_Warn(a,__LINE__,__FORTRANFILE__,CCTK_THORNSTRING,b)

#define CCTK_CoordRegisterSystem(a,b,c) CCTKi_CoordRegisterSystem(a,b,CCTK_THORNSTRING,c)

/* John Shalf says that the operator .ne. needs to be enclosed by
   spaces, because ANSI C preprocessors otherwise interpret the
   character sequence ".0" as preprocessor token */
#define CCTK_EQUALS(a,b) (CCTK_Equals(a,b) .ne. 0)

#define CCTK_PASS_FTOF CCTK_FARGUMENTS

#define CCTK_ORIGIN_SPACE(x) (cctk_origin_space(x)+cctk_delta_space(x)/cctk_levfac(x)*cctk_levoff(x)/cctk_levoffdenom(x))
#define CCTK_DELTA_SPACE(x) (cctk_delta_space(x)/cctk_levfac(x))
#define CCTK_DELTA_TIME (cctk_delta_time/cctk_timefac)
#define CCTK_LSSH(stag,dim) cctk_lssh((stag)+CCTK_NSTAGGER*(dim))
#define CCTK_LSSH_IDX(stag,dim) ((stag)+CCTK_NSTAGGER*(dim))

#ifdef F90CODE

#define _DECLARE_CCTK_FUNCTIONS                   \
  external     CCTK_PointerTo                   &&\
  CCTK_POINTER CCTK_PointerTo                   &&\
  interface                                     &&\
     integer function CCTK_Equals (arg1, arg2)  &&\
       implicit none                            &&\
       CCTK_POINTER_TO_CONST arg1               &&\
       character(*) arg2                        &&\
     end function CCTK_Equals                   &&\
     integer function CCTK_MyProc (cctkGH)      &&\
       implicit none                            &&\
       CCTK_POINTER_TO_CONST cctkGH             &&\
     end function CCTK_MyProc                   &&\
     integer function CCTK_nProcs (cctkGH)      &&\
       implicit none                            &&\
       CCTK_POINTER_TO_CONST cctkGH             &&\
     end function CCTK_nProcs                   &&\
     integer function CCTK_IsThornActive (name) &&\
       implicit none                            &&\
       character(*) name                        &&\
     end function CCTK_IsThornActive            &&\
     CCTK_POINTER function CCTK_NullPointer ()  &&\
       implicit none                            &&\
     end function CCTK_NullPointer              &&\
  end interface                                 &&

#else /* ! F90CODE */

#define _DECLARE_CCTK_FUNCTIONS \
  integer      CCTK_Equals, CCTK_MyProc, CCTK_nProcs, CCTK_IsThornActive &&\
  external     CCTK_Equals, CCTK_MyProc, CCTK_nProcs, CCTK_IsThornActive &&\
  CCTK_POINTER CCTK_PointerTo, CCTK_NullPointer &&\
  external     CCTK_PointerTo, CCTK_NullPointer &&

#endif /* ! F90CODE */

#endif /*FCODE*/

#ifdef CCODE

/* get the definition of ptrdiff_t */
#include <stddef.h>

#include "cGH.h"

#include "cctk_ActiveThorns.h"
#include "cctk_Banner.h"
#include "cctk_Cache.h"
#include "cctk_Coord.h"
#include "cctk_Comm.h"
#include "cctk_CommandLine.h"
#include "cctk_Complex.h"
#include "cctk_DebugDefines.h"
#include "cctk_Faces.h"
#include "cctk_File.h"
#include "cctk_Flesh.h"
#include "cctk_FortranString.h"
#include "cctk_Functions.h"
#include "cctk_GHExtensions.h"
#include "cctk_Groups.h"
#include "cctk_GroupsOnGH.h"
#include "cctk_Interp.h"
#include "cctk_IO.h"
#include "cctk_IOMethods.h"
#include "cctk_Loop.h"
#include "cctk_Main.h"
#include "cctk_Malloc.h"
#include "cctk_Misc.h"
#include "cctk_Parameter.h"
#include "cctk_Reduction.h"
#include "cctk_Stagger.h"
#include "cctk_Sync.h"
#include "cctk_Timers.h"
#include "cctk_WarnLevel.h"


/*
 * routines/macros to compute the linear index
 * of a grid funtion element from its i/j/k dimensions
 *
 * These are defined as inline functions when the language supports this,
 * otherwise they are defined as macros.
 * For CCTK_DEBUG these are external C routines defined in DebugDefines.c.
 */

#ifdef __cplusplus
#define HAVE_INLINE
#else
#ifndef inline
#define HAVE_INLINE
#endif
#endif

#ifdef HAVE_INLINE

/* The "inline" keyword is supported */

static inline int CCTK_GFINDEX1D (const cGH *GH, int i);
static inline int CCTK_GFINDEX1D (const cGH *GH, int i)
{
  GH = GH;
  return (i);
}

static inline int CCTK_GFINDEX2D (const cGH *GH, int i, int j);
static inline int CCTK_GFINDEX2D (const cGH *GH, int i, int j)
{
  return (i + GH->cctk_lsh[0]*j);
}

static inline int CCTK_GFINDEX3D (const cGH *GH, int i, int j, int k);
static inline int CCTK_GFINDEX3D (const cGH *GH, int i, int j, int k)
{
  return (i + GH->cctk_lsh[0]*(j + GH->cctk_lsh[1]*k));
}

static inline int CCTK_GFINDEX4D (const cGH *GH, int i, int j, int k, int l);
static inline int CCTK_GFINDEX4D (const cGH *GH, int i, int j, int k, int l)
{
  return (i + GH->cctk_lsh[0]*(j + GH->cctk_lsh[1]*(k + GH->cctk_lsh[2] * l)));
}

#else /* ! defined(__cplusplus) && defined(inline) */

#ifdef CCTK_DEBUG
/* The "inline" keyword is not supported, and we want to debug */

#define CCTK_GFINDEX1D CCTK_GFIndex1D
#define CCTK_GFINDEX2D CCTK_GFIndex2D
#define CCTK_GFINDEX3D CCTK_GFIndex3D
#define CCTK_GFINDEX4D CCTK_GFIndex4D

#else /* ! defined(CCTK_DEBUG) */
/* The "inline" keyword is not supported, and we want to optimise */

#define CCTK_GFINDEX1D(GH, i)          (i)
#define CCTK_GFINDEX2D(GH, i, j)       ((i) + (GH)->cctk_lsh[0] * (j))
#define CCTK_GFINDEX3D(GH, i, j, k)    ((i) + (GH)->cctk_lsh[0] *             \
                                        ((j) + (GH)->cctk_lsh[1] * (k)))
#define CCTK_GFINDEX4D(GH, i, j, k, l) ((i) + (GH)->cctk_lsh[0] *             \
                                        ((j) + (GH)->cctk_lsh[1] *            \
                                         ((k) + (GH)->cctk_lsh[2] * (l))))

#endif /* ! defined(CCTK_DEBUG) */

#endif /* ! defined(__cplusplus) && defined(inline) */


#define CCTK_PRINTSEPARATOR \
  printf("--------------------------------------------------------------------------------\n");

#define _DECLARE_CCTK_ARGUMENTS _DECLARE_CCTK_CARGUMENTS
#define _DECLARE_CCTK_CARGUMENTS \
            ptrdiff_t   cctki_dummy_int;\
            int         cctk_dim = cctkGH->cctk_dim;\
            int        *cctk_gsh = cctkGH->cctk_gsh;\
            int        *cctk_lsh = cctkGH->cctk_lsh;\
            int        *cctk_lbnd = cctkGH->cctk_lbnd;\
            int        *cctk_ubnd = cctkGH->cctk_ubnd;\
            int        *cctk_lssh = cctkGH->cctk_lssh;\
            int        *cctk_from = cctkGH->cctk_from;\
            int        *cctk_to = cctkGH->cctk_to;\
            int        *cctk_bbox = cctkGH->cctk_bbox;\
            CCTK_REAL   cctk_delta_time = cctkGH->cctk_delta_time;\
            CCTK_REAL   cctk_time = cctkGH->cctk_time;\
            CCTK_REAL  *cctk_delta_space = cctkGH->cctk_delta_space;\
            CCTK_REAL  *cctk_origin_space = cctkGH->cctk_origin_space;\
            int        *cctk_levfac = cctkGH->cctk_levfac;\
            int        *cctk_levoff = cctkGH->cctk_levoff;\
            int        *cctk_levoffdenom = cctkGH->cctk_levoffdenom;\
            int         cctk_timefac = cctkGH->cctk_timefac;\
            int         cctk_convlevel = cctkGH->cctk_convlevel;\
            int         cctk_convfac = cctkGH->cctk_convfac;\
            int        *cctk_nghostzones = cctkGH->cctk_nghostzones;\
            int         cctk_iteration = cctkGH->cctk_iteration;\
            enum { \
		    cctk_dummyi_dummy_int = sizeof( cctki_dummy_int ), \
		    cctk_dummy_dim = sizeof( cctk_dim ), \
		    cctk_dummy_gsh = sizeof( cctk_gsh ), \
		    cctk_dummy_lsh = sizeof( cctk_lsh ), \
		    cctk_dummy_lbnd = sizeof( cctk_lbnd ), \
		    cctk_dummy_ubnd = sizeof( cctk_ubnd ), \
		    cctk_dummy_lssh = sizeof( cctk_lssh ), \
		    cctk_dummy_from = sizeof( cctk_from ), \
		    cctk_dummy_to = sizeof( cctk_to ), \
		    cctk_dummy_bbox = sizeof( cctk_bbox ), \
		    cctk_dummy_delta_time = sizeof( cctk_delta_time ), \
		    cctk_dummy_time = sizeof( cctk_time ), \
		    cctk_dummy_delta_space = sizeof( cctk_delta_space ), \
		    cctk_dummy_origin_space = sizeof( cctk_origin_space ), \
		    cctk_dummy_levoff = sizeof( cctk_levoff ), \
		    cctk_dummy_levfac = sizeof( cctk_levfac ), \
		    cctk_dummy_levoffdenom = sizeof( cctk_levoffdenom ), \
		    cctk_dummy_timefac = sizeof( cctk_timefac ), \
		    cctk_dummy_convlevel = sizeof( cctk_convlevel ), \
		    cctk_dummy_convfac = sizeof( cctk_convfac ), \
		    cctk_dummy_nghostzones = sizeof( cctk_nghostzones ), \
		    cctk_dummy_iteration = sizeof( cctk_iteration ) \
	    };


#define _INITIALISE_CCTK_C2F
#define _DECLARE_CCTK_C2F
#define _PASS_CCTK_C2F(xGH) &((xGH)->cctk_dim),\
                            (xGH)->cctk_gsh,(xGH)->cctk_lsh,\
                            (xGH)->cctk_lbnd,(xGH)->cctk_ubnd,\
                            (xGH)->cctk_lssh,\
                            (xGH)->cctk_from,(xGH)->cctk_to,\
                            (xGH)->cctk_bbox,\
                            &((xGH)->cctk_delta_time),\
                            &((xGH)->cctk_time), (xGH)->cctk_delta_space,\
                            (xGH)->cctk_origin_space,\
                            (xGH)->cctk_levfac,\
                            (xGH)->cctk_levoff,\
                            (xGH)->cctk_levoffdenom,\
                            &((xGH)->cctk_timefac),\
                            &((xGH)->cctk_convlevel),\
                            &((xGH)->cctk_convfac),\
                            (xGH)->cctk_nghostzones,\
                            &((xGH)->cctk_iteration),\
                            &(xGH)
#define _CCTK_C2F_PROTO     int *,\
                            int *,\
                            int *,int *, int *, int *, int *,int *,int *,\
                            CCTK_REAL *, CCTK_REAL *, CCTK_REAL *,\
                            CCTK_REAL *,\
                            int *,\
                            int *,\
                            int *,\
                            int *,\
                            int *,\
                            int *,\
                            int *,\
                            int *,\
                            cGH **

#define CCTK_EQUALS(a,b) (CCTK_Equals((a),(b)))

#define CCTK_PASS_CTOC cctkGH

#define CCTK_ORIGIN_SPACE(x) (cctk_origin_space[x]+cctk_delta_space[x]/cctk_levfac[x]*cctk_levoff[x]/cctk_levoffdenom[x])
#define CCTK_DELTA_SPACE(x) (cctk_delta_space[x]/cctk_levfac[x])
#define CCTK_DELTA_TIME (cctk_delta_time/cctk_timefac)
#define CCTK_LSSH(stag,dim) cctk_lssh[(stag)+CCTK_NSTAGGER*(dim)]
#define CCTK_LSSH_IDX(stag,dim) ((stag)+CCTK_NSTAGGER*(dim))

#define CCTK_WARN(a,b) CCTK_Warn(a,__LINE__,__FILE__,CCTK_THORNSTRING,b)

#define CCTK_MALLOC(s) CCTKi_Malloc(s,__LINE__,__FILE__)
#define CCTK_FREE(p) CCTKi_Free(p)

#endif /*CCODE*/

#define CCTK_INFO(a) CCTK_Info(CCTK_THORNSTRING,(a))
#define CCTK_PARAMWARN(a) CCTK_ParamWarn(CCTK_THORNSTRING,(a))

/*
#define CCTK_MAKESTRING(x) CCTK_REALSTRING(x)
#define CCTK_REALSTRING(x) #x
#define CCTK_WARN(a,b) CCTK_Warn(a,CCTK_MAKESTRING(CCTK_THORN),b,__LINE__,__FILE__)
*/


#endif
