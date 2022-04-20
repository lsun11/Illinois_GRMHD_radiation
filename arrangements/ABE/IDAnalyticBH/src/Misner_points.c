 /*@@
   @file      Misner_points.c
   @date
   @author    Steve Brandt
   @desc
              This calculates the conformal factor for nbh black holes,
              with naked mass m0 = 2 csch(mu) each, and placed on a circle in the
              xy plane around the origin, of radius coth(mu).
              One of them sits on the positive x axis, the others are evenly spaced.
              Naked mass here corresponds to the term m0 / (2 |r - r0|) in the expansion.
   @history
   @hdate 6.May.2003
   @hauthor Jonathan Thornburg
   @hdesc code cleanup, add prototypes
   @endhistory
   @enddesc
 @@*/

#include <stdlib.h>
#include <stdio.h>

#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "IDAnalyticBH.h"

static const char *rcsid = "$Header: /cactus/CactusEinstein/IDAnalyticBH/src/Misner_points.c,v 1.13 2004/08/12 22:10:39 yye00 Exp $";

CCTK_FILEVERSION(CactusEinstein_IDAnalyticBH_Misner_points_c)

/******************************************************************************/

/*
 * data structures local to this file
 */

/* Basic data about a brill-lindquist black hole term. */
struct bhole {
  CCTK_REAL x,y;
  CCTK_REAL mass;

  /* i gives either the number of the seed
     black hole we are starting with, or
     the number of the seed black hole that
     was used to isometrize this term. */
  int i;

  /* isos points to an array of nbhole (or is it nbhole-1??)
     structures, malloc()-ed in  fill_iso() .
     At present this array is never freed, and in fact the pointers may
     be overwritten with other malloc()s in fill_iso() recursive calls! :( */
  struct bhole *isos;
};

/******************************************************************************/

/*
 * static data
 */

static int nbholes;

#define MAXBHOLES 10
/* The seed black holes. */
struct bhole bholes[MAXBHOLES];

/******************************************************************************/

/*
 * prototypes for functions local to this file
 */
static CCTK_REAL csch(CCTK_REAL theta);
static CCTK_REAL coth(CCTK_REAL theta);
static void iso(struct bhole *a1, struct bhole *a2, struct bhole *a3);
static void fill_iso(struct bhole *b, int n);
static CCTK_REAL eval_bh_psi(const struct bhole *b,
                             CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);

/******************************************************************************/
/***** functions visible outside this file ************************************/
/******************************************************************************/

 /*@@
   @routine    Misner_init
   @date
   @author     Steve Brandt
   @desc
      Initialises the black holes then makes the isometry black holes
   @enddesc
   @calls
   @history
   @endhistory
@@*/
void Misner_init(int n, CCTK_REAL mu, int terms)
{

  int i;
  CCTK_REAL pi,ang;

  if (! ((nbholes=n) < MAXBHOLES))
    CCTK_VWarn (0, __LINE__, __FILE__,"IDAnalyticBH", "Number of black holes exceed maximum number of black holes");

  pi = 4.0*atan(1.);

  ang = 2.*pi/(n);

  for(i=0;i<n;i++)
  {
    bholes[i].x = coth(mu)*cos(ang*i);
    bholes[i].y = coth(mu)*sin(ang*i);
    bholes[i].mass = csch(mu);
    bholes[i].i = i;
    bholes[i].isos = 0;
  }

  CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   Misner_init(): warning: about to call  fill_iso() ; at present\n"
"                  this routine seems to leak quasi-infinite amounts of\n"
"                  memory at the rate of O(100 megabytes/second) :( :( :(\n");

  for(i=0;i<n;i++)
    fill_iso(&bholes[i],terms);

}

/******************************************************************************/

 /*@@
   @routine    MisnerEvalPsi
   @date
   @author     Steve Brandt
   @desc
      Evaluate psi at a point
   @enddesc
   @calls      eval_bh_psi

   @history
   @hdate      6.May.2003
   @hauthor    Jonathan Thornburg <jthorn@aei.mpg.de>
   @hdesc      Return result directly rather than via a pointer argument.
   @endhistory
 @@*/
CCTK_REAL MisnerEvalPsi(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z)
{
  int i;
  CCTK_REAL sum = 1.0;

  for(i=0;i<nbholes;i++)
  {
    sum += eval_bh_psi(&bholes[i],x,y,z);
  }

  return sum;
}

/******************************************************************************/
/***** functions local to this file *******************************************/
/******************************************************************************/

static CCTK_REAL csch(CCTK_REAL theta) {
  return 1.0/sinh(theta);
}

/******************************************************************************/

static CCTK_REAL coth(CCTK_REAL theta) {
  return cosh(theta)/sinh(theta);
}

/******************************************************************************/

 /*@@
   @routine    iso
   @date
   @author     Steve Brandt
   @desc
               Isometrize black hole a1 through hole a2
   @enddesc
   @calls      iso
   @history

   @endhistory

@@*/

static void iso(struct bhole *a1, struct bhole *a2, struct bhole *a3)
{
  CCTK_REAL rad,radtwo;
  radtwo=(
      (a1->x - a2->x)*(a1->x - a2->x)+
      (a1->y - a2->y)*(a1->y - a2->y)
    );
  rad=sqrt(radtwo);
  a3->mass = a1->mass*a2->mass/rad;
  a3->x = a2->x+(a2->mass*a2->mass)*(a1->x - a2->x)/radtwo;
  a3->y = a2->y+(a2->mass*a2->mass)*(a1->y - a2->y)/radtwo;
}

/******************************************************************************/

 /*@@
   @routine    fill_iso
   @date
   @author     Steve Brandt
   @desc
               Fills in the iso structure of a given black hole.
               Applies recursively to the number of terms desired.
   @enddesc
   @calls      fill_iso
   @history

   @endhistory

@@*/

static void fill_iso(struct bhole *b, int n)
{
  int i,j;
  if(n==0)
  {
    b->isos = 0;
    return;
  }
  b->isos = (struct bhole *)malloc(sizeof(struct bhole)*(nbholes-1));
  if (! (b->isos != 0))
    CCTK_VWarn (0, __LINE__, __FILE__,"IDAnalyticBH", "error in function fill_iso");
  for(j=0, i=0;i<nbholes;i++)
  {
    if(i != b->i) {
      iso(b,&bholes[i],&b->isos[j]);
      b->isos[j].i = i;
      fill_iso(&b->isos[j],n-1);
      j++;
    }
  }
}

/******************************************************************************/

 /*@@
   @routine    eval_bh_psi
   @date
   @author     Steve Brandt
   @desc

   @enddesc
   @calls      eval_bh_psi
   @history

   @endhistory

@@*/
static CCTK_REAL eval_bh_psi(const struct bhole *b,
                             CCTK_REAL x, CCTK_REAL y, CCTK_REAL z)
{
  int i;
  CCTK_REAL res;
  res = 0.0;
  if(b->isos != 0)
  {
    for(i=0;i<nbholes-1;i++)
    {
      res += eval_bh_psi(&b->isos[i],x,y,z);
    }
  }
  res += b->mass/sqrt(
      (x - b->x)*(x - b->x)+
      (y - b->y)*(y - b->y)+
      z*z
    );
  return res;
}
