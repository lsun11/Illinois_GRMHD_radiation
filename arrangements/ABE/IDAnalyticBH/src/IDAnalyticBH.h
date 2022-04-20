 /*@@
   @file      IDAnalyticBH.h
   @date      23.Mar.2003
   @author    Jonathan Thornburg <jthorn@aei.mpg.de>
   @desc      Overall header file for this thorn; defines prototypes
              externally visible routines.
   @enddesc 
   @version $Header: /cactus/CactusEinstein/IDAnalyticBH/src/IDAnalyticBH.h,v 1.3 2004/05/17 12:40:46 goodale Exp $
   @history
   @hdate 6.May.2003
   @hauthor Jonathan Thornburg
   @hdesc   add lots more prototypes
   @endhistory
 @@*/

/******************************************************************************/

#define SQR(a) ((a)*(a))

/******************************************************************************/

/*
 * prototypes for routines called directly by the Cactus scheduler
 */

void IDAnalyticBH_ParamChecker(CCTK_ARGUMENTS);         /* ParamChecker.c */
void Schwarzschild(CCTK_ARGUMENTS);                     /* Schwarzschild.c */
void KerrID(CCTK_ARGUMENTS);                            /* Kerr.c */
void BrillLindquist(CCTK_ARGUMENTS);                    /* BrillLindquist.c */
void Misner_standard(CCTK_ARGUMENTS);                   /* Misner_standard.c */
void Misner_multiple(CCTK_ARGUMENTS);                   /* Misner_multiple.c */

/******************************************************************************/

/*
 * prototypes for other routines visible outside their source files
 */

/* Misner_points.c */
void Misner_init(int n, CCTK_REAL mu, int terms);
CCTK_REAL MisnerEvalPsi(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);

/* misc.c */
void IDAnalyticBH_zero_CCTK_REAL_array(int N, CCTK_REAL A[]);
