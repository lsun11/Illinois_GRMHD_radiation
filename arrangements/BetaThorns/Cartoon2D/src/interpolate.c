/* interpolate.c -- interpolation functions */
/* $Header: /cactusdevcvs/BetaThorns/Cartoon2D/src/interpolate.c,v 1.12 2005/11/10 14:05:35 hawke Exp $ */

/*
 * interpolate_local - ... (local)
 * interpolate_local_order1 - ... (order 1) (linear) (2-point)
 * interpolate_local_order2 - ... (order 2) (3-point)
 * interpolate_local_order3 - ... (order 3) (4-point)
 * interpolate_local_order4 - ... (order 4) (5-point)
 * interpolate_local_order5 - ... (order 5) (6-point)
 */

/* To make porting easy, only the interpolate_local* functions 
 * have been retained in this file. JT's original code (kept in 
 * the archive directory) also has functions for derivatives and
 * interpolation on non-uniform grids. ---Sai.
 */

#include <math.h>
#include <stdio.h>

#include "cctk.h"

CCTK_REAL interpolate_local(int order, CCTK_REAL x0, CCTK_REAL dx, 
                            CCTK_REAL y[],
                            CCTK_REAL x,
                            CCTK_INT excised_points[]);

CCTK_REAL interpolate_eno(CCTK_INT order, CCTK_REAL x0, CCTK_REAL dx,
                          CCTK_REAL y[],
                          CCTK_REAL x,
                            CCTK_INT excised_points[]);

/* prototypes for private functions defined in this file */
static CCTK_REAL interpolate_local_order1(CCTK_REAL x0, CCTK_REAL dx,
                                          CCTK_REAL y[],
                                          CCTK_REAL x);
static CCTK_REAL interpolate_local_order2(CCTK_REAL x0, CCTK_REAL dx,
                                          CCTK_REAL y[],
                                          CCTK_REAL x);
static CCTK_REAL interpolate_local_order3(CCTK_REAL x0, CCTK_REAL dx,
                                          CCTK_REAL y[],
                                          CCTK_REAL x);
static CCTK_REAL interpolate_local_order4(CCTK_REAL x0, CCTK_REAL dx,
                                          CCTK_REAL y[],
                                          CCTK_REAL x);
static CCTK_REAL interpolate_local_order5(CCTK_REAL x0, CCTK_REAL dx,
                                          CCTK_REAL y[],
                                          CCTK_REAL x);

/*****************************************************************************/

/*
 * This function does local Lagrange polynomial interpolation on a
 * uniform grid.  That is, given  order+1  uniformly spaced data points
 *  (x0 + i*dx, y[i]) , i = 0...order, it interpolates a degree-(order)
 * (Lagrange) polynomial through them and evaluates this polynomial at
 * a specified  x  coordinate.  In general the error in this computed
 * function value is  O(h^(order+1)) .
 *
 * Except for possible end effects (see  end_action  (below)), the local
 * interpolation is internally centered to improve its conditioning.  Even
 * so, the interpolation is highly ill-conditioned for small  dx  and/or
 * large  order  and/or  x  outside the domain of the data points.
 *
 * The interpolation formulas were (are) all derived via Maple, see
 * "./interpolate.in" and "./interpolate.out".
 *
 * Arguments:
 * order = (in) The order of the local interpolation, i.e. the degree
 *		of the interpolated polynomial.
 * x0 = (in) The x coordinate corresponding to y[0].
 * dx = (in) The x spacing between the data points.
 * y[0...order] = (in) The y data array.
 * x = (in) The x coordinate for the interpolation.
 */
CCTK_REAL interpolate_local(int order,
                            CCTK_REAL x0, CCTK_REAL dx,
                            CCTK_REAL y[],
                            CCTK_REAL x,
                            CCTK_INT excised_points[])
{

  /* 
     We assume that the excised region has one of the following 4 forms:
     
     o o o o o         (No excision)
     x o o o o         (left excision)
     o o o o x         (right excision)
     x o o o x         (left and right excision)

     Anything else (e.g., x o x o x) will go horribly wrong.
  */

  int start, end, available_order;
  
  start = 0;
  end   = order;
  
  while ((start < order + 1)&&(excised_points[start]))
  {
    ++start;
  }
  while ((end > -1)&&(excised_points[end]))
  {
    --end;
  }
  
  if (end < start)
  {
    /* The whole block is excised */
    return y[0];
  }

  /* The possible interpolation order that can be used */

  available_order = end - start;

  if (available_order == 0)
  {
    /* There is only one non-excised point */
    return y[start];
  }

  if (dx * start + x0 > x)
  {
    /* We would have to extrapolate to the left */
    return y[start];
  }
  
  if (dx * end + x0 < x)
  {
    /* We would have to extrapolate to the right */
    return y[end];
  }

  switch (available_order) 
  {
  case 1: return interpolate_local_order1(x0 + start * dx, dx, &y[start], x);
  case 2: return interpolate_local_order2(x0 + start * dx, dx, &y[start], x);
  case 3: return interpolate_local_order3(x0 + start * dx, dx, &y[start], x);
  case 4: return interpolate_local_order4(x0 + start * dx, dx, &y[start], x);
  case 5: return interpolate_local_order5(x0 + start * dx, dx, &y[start], x);
  default: CCTK_WARN(0, "Interpolation order not supported");
  }
}

/*****************************************************************************/

/*
 * This function interpolates a 1st-order (linear) Lagrange polynomial
 * in a 2-point [0...1] data array centered about the [0.5] grid position.
 */
static CCTK_REAL interpolate_local_order1(CCTK_REAL x0, CCTK_REAL dx,
                                          CCTK_REAL y[],
                                          CCTK_REAL x)
{
CCTK_REAL c0 = (+1.0/2.0)*y[0] + (+1.0/2.0)*y[1];
CCTK_REAL c1 =          - y[0] +          + y[1];

CCTK_REAL xc = x0 + 0.5*dx;
CCTK_REAL xr = (x - xc) / dx;

return(  c0 + xr*c1  );
}

/*****************************************************************************/

/*
 * This function interpolates a 2nd-order (quadratic) Lagrange polynomial
 * in a 3-point [0...2] data array centered about the [1] grid position.
 */
static CCTK_REAL interpolate_local_order2(CCTK_REAL x0, CCTK_REAL dx,
                                          CCTK_REAL y[],
                                          CCTK_REAL x)
{
CCTK_REAL c0 = y[1];
CCTK_REAL c1 = (-1.0/2.0)*y[0]        + (+1.0/2.0)*y[2];
CCTK_REAL c2 = (+1.0/2.0)*y[0] - y[1] + (+1.0/2.0)*y[2];

CCTK_REAL xc = x0 + 1.0*dx;
CCTK_REAL xr = (x - xc) / dx;

return(  c0 + xr*(c1 + xr*c2)  );
}

/*****************************************************************************/

/*
 * This function interpolates a 3rd-order (cubic) Lagrange polynomial
 * in a 4-point [0...3] data array centered about the [1.5] grid position.
 */
static CCTK_REAL interpolate_local_order3(CCTK_REAL x0, CCTK_REAL dx,
                                          CCTK_REAL y[],
                                          CCTK_REAL x)
{
CCTK_REAL c0 =   (-1.0/16.0)*y[0] + (+9.0/16.0)*y[1]
	    + (+9.0/16.0)*y[2] + (-1.0/16.0)*y[3];
CCTK_REAL c1 =   (+1.0/24.0)*y[0] + (-9.0/ 8.0)*y[1]
	    + (+9.0/ 8.0)*y[2] + (-1.0/24.0)*y[3];
CCTK_REAL c2 =   (+1.0/ 4.0)*y[0] + (-1.0/ 4.0)*y[1]
	    + (-1.0/ 4.0)*y[2] + (+1.0/ 4.0)*y[3];
CCTK_REAL c3 =   (-1.0/ 6.0)*y[0] + (+1.0/ 2.0)*y[1]
	    + (-1.0/ 2.0)*y[2] + ( 1.0/ 6.0)*y[3];

CCTK_REAL xc = x0 + 1.5*dx;
CCTK_REAL xr = (x - xc) / dx;

return(  c0 + xr*(c1 + xr*(c2 + xr*c3))  );
}

/*****************************************************************************/

/*
 * This function interpolates a 4th-order (quartic) Lagrange polynomial
 * in a 5-point [0...4] data array centered about the [2] grid position.
 */
static CCTK_REAL interpolate_local_order4(CCTK_REAL x0, CCTK_REAL dx,
                                          CCTK_REAL y[],
                                          CCTK_REAL x)
{
CCTK_REAL c0 = y[2];
CCTK_REAL c1 =   (+1.0/12.0)*y[0] + (-2.0/ 3.0)*y[1]
	    + (+2.0/ 3.0)*y[3] + (-1.0/12.0)*y[4];
CCTK_REAL c2 = + (-1.0/24.0)*y[0] + (+2.0/ 3.0)*y[1] + (-5.0/ 4.0)*y[2]
	    + (+2.0/ 3.0)*y[3] + (-1.0/24.0)*y[4];
CCTK_REAL c3 = + (-1.0/12.0)*y[0] + (+1.0/ 6.0)*y[1]
	    + (-1.0/ 6.0)*y[3] + (+1.0/12.0)*y[4];
CCTK_REAL c4 = + (+1.0/24.0)*y[0] + (-1.0/ 6.0)*y[1] + (+1.0/ 4.0)*y[2]
	    + (-1.0/ 6.0)*y[3] + (+1.0/24.0)*y[4];

CCTK_REAL xc = x0 + 2.0*dx;
CCTK_REAL xr = (x - xc) / dx;

return(  c0 + xr*(c1 + xr*(c2 + xr*(c3 + xr*c4)))  );
}

/*****************************************************************************/

/*
 * This function interpolates a 5th-order (quintic) Lagrange polynomial
 * in a 6-point [0...5] data array centered about the [2.5] grid position.
 */
static CCTK_REAL interpolate_local_order5(CCTK_REAL x0, CCTK_REAL dx,
                                          CCTK_REAL y[],
                                          CCTK_REAL x)
{
CCTK_REAL c0 =   (+ 3.0/256.0)*y[0] + (-25.0/256.0)*y[1]
	    + (+75.0/128.0)*y[2] + (+75.0/128.0)*y[3]
	    + (-25.0/256.0)*y[4] + (+ 3.0/256.0)*y[5];
CCTK_REAL c1 =   (- 3.0/640.0)*y[0] + (+25.0/384.0)*y[1]
	    + (-75.0/ 64.0)*y[2] + (+75.0/ 64.0)*y[3]
	    + (-25.0/384.0)*y[4] + (+ 3.0/640.0)*y[5];
CCTK_REAL c2 =   (- 5.0/ 96.0)*y[0] + (+13.0/ 32.0)*y[1]
	    + (-17.0/ 48.0)*y[2] + (-17.0/ 48.0)*y[3]
	    + (+13.0/ 32.0)*y[4] + (- 5.0/ 96.0)*y[5];
CCTK_REAL c3 =   (+ 1.0/ 48.0)*y[0] + (-13.0/ 48.0)*y[1]
	    + (+17.0/ 24.0)*y[2] + (-17.0/ 24.0)*y[3]
	    + (+13.0/ 48.0)*y[4] + (- 1.0/ 48.0)*y[5];
CCTK_REAL c4 =   (+ 1.0/ 48.0)*y[0] + (- 1.0/ 16.0)*y[1]
	    + (+ 1.0/ 24.0)*y[2] + (+ 1.0/ 24.0)*y[3]
	    + (- 1.0/ 16.0)*y[4] + (+ 1.0/ 48.0)*y[5];
CCTK_REAL c5 =   (- 1.0/120.0)*y[0] + (+ 1.0/ 24.0)*y[1]
	    + (- 1.0/ 12.0)*y[2] + (+ 1.0/ 12.0)*y[3]
	    + (- 1.0/ 24.0)*y[4] + (+ 1.0/120.0)*y[5];

CCTK_REAL xc = x0 + 2.5*dx;
CCTK_REAL xr = (x - xc) / dx;

return(  c0 + xr*(c1 + xr*(c2 + xr*(c3 + xr*(c4 + xr*c5))))  );
}


/*****************************************************************************/

CCTK_REAL interpolate_eno(CCTK_INT order, 
                          CCTK_REAL x0, CCTK_REAL dx,
                          CCTK_REAL y[],
                          CCTK_REAL x,
                          CCTK_INT excised_points[])
{

/*   CCTK_REAL diff[4]; */
/*   CCTK_INT seed = 0; */
/*   CCTK_INT j,r;   */

/*   CCTK_REAL c0; */
/*   CCTK_REAL c1; */
/*   CCTK_REAL c2; */

/*   CCTK_REAL xc; */
/*   CCTK_REAL xr; */

  CCTK_REAL result;
  CCTK_REAL undiv_diff [13][6];
  CCTK_INT diff, i, ii, r;
  CCTK_REAL fp_i, fp_ii, x_rel;

  /* Find seed index */
/*   while (x > x0 + dx * ((CCTK_REAL)seed+0.5)) seed++; */

/*   if (seed!=2) { */
    /* Not enough stencil, only perform linear interpolation */   
/*     seed = 0; */
/*     while (x > x0 + dx * ((CCTK_REAL)(seed+1))) seed++; */
/*     if (seed==4) seed=3; */
/*     return y[seed] + (y[seed+1]-y[seed])/dx * (x-x0-(CCTK_REAL)seed*dx); */
/*   } */

  /*  Calculate the undivided differences */
/*   for (j=0; j<=3; j++) */
/*     diff[j] = y[j+1] - y[j]; */

  /*  Find the stencil */
/*   if ( fabs(diff[1]) < fabs(diff[2]) ) { */
/*     if ( fabs(diff[1]-diff[0]) < fabs(diff[2]-diff[1]) ) */
/*       r = 0; */
/*     else */
/*       r = 1; */
/*   } else { */
/*     if ( fabs(diff[2]-diff[1]) < fabs(diff[3]-diff[2]) ) */
/*       r = 1; */
/*     else */
/*       r = 2; */
/*   } */

  /* Interpolate second order */
/*   c0 = y[r+1]; */
/*   c1 = (-1.0/2.0)*y[r] + (+1.0/2.0)*y[r+2]; */
/*   c2 = (+1.0/2.0)*y[r] - y[r+1] + (+1.0/2.0)*y[r+2]; */

/*   xc = x0 + dx * (CCTK_REAL)(r+1); */
/*   xr = (x - xc) / dx; */
/*   result = (  c0 + xr*(c1 + xr*c2)  ); */

  for (i = 1; i < 2 * order + 2; ++i)
  {
    undiv_diff[i][0] = y[i-1];
  }
  for (i = 0; i < 6; ++i)
  {
    undiv_diff[0][i] = 1e10;
    undiv_diff[2 * order + 2][i] = 1e10;
  }

  for (diff = 1; diff < order + 1; ++diff)
  {
    for (i = 1; i < 2 * order + 2; ++i)
    {
      undiv_diff[i][diff] = undiv_diff[i][diff-1] - undiv_diff[i-1][diff-1];
    }
    undiv_diff[0][diff] = -10 * undiv_diff[1][diff];
    undiv_diff[2 * order + 2][diff]= -10 * undiv_diff[2 * order + 1][diff];
  }
  
  fp_i = (x - x0) / dx;
  fp_ii = floor(fp_i);
  x_rel = fp_i - fp_ii;
  
  ii = (CCTK_INT)(fp_ii) + 1;
  
  if (ii < 1)
  {
    ii = 1;
    x_rel = fp_i - (CCTK_REAL)(ii) + 1.0;
  }
  else if (ii > 2 * order + 1)
  {
    ii = 2 * order + 1;
    x_rel = fp_i - (CCTK_REAL)(ii) + 1.0;
  }
  
  r = 0;
  for (diff = 1; diff < order + 1; ++diff)
  {
    if (fabs(undiv_diff[ii-r][diff]) < fabs(undiv_diff[ii-r+1][diff]) )
    {
      ++r;
    }
    
    if (ii - r < diff + 1)
    {
      --r;
    }
    else if (ii - r + diff > 2 * order + 1)
    {
      ++r;
    }
  }

  result = undiv_diff[ii-r][0];

  switch (order)
  {
    case 1:
      result += (x_rel + r - 0.0) / 1.0 * undiv_diff[ii-r+1][1];
      break;
    case 2:
      result += (x_rel + r - 0.0) / 1.0 * (undiv_diff[ii-r+1][1] +
                (x_rel + r - 1.0) / 2.0 * (undiv_diff[ii-r+2][2]));
      break;
    case 3:
      result += (x_rel + r - 0.0) / 1.0 * (undiv_diff[ii-r+1][1] +
                (x_rel + r - 1.0) / 2.0 * (undiv_diff[ii-r+2][2] +
                (x_rel + r - 2.0) / 3.0 * (undiv_diff[ii-r+3][3])));
      break;
    case 4:
      result += (x_rel + r - 0.0) / 1.0 * (undiv_diff[ii-r+1][1] +
                (x_rel + r - 1.0) / 2.0 * (undiv_diff[ii-r+2][2] +
                (x_rel + r - 2.0) / 3.0 * (undiv_diff[ii-r+3][3] +
                (x_rel + r - 3.0) / 4.0 * (undiv_diff[ii-r+4][4]))));
      break;
    case 5:
      result += (x_rel + r - 0.0) / 1.0 * (undiv_diff[ii-r+1][1] +
                (x_rel + r - 1.0) / 2.0 * (undiv_diff[ii-r+2][2] +
                (x_rel + r - 2.0) / 3.0 * (undiv_diff[ii-r+3][3] +
                (x_rel + r - 3.0) / 4.0 * (undiv_diff[ii-r+4][4] +
                (x_rel + r - 4.0) / 5.0 * (undiv_diff[ii-r+5][5])))));
      break;
  }
  
/* #define CARTOON_ENO_DBG 1 */
#ifdef CARTOON_ENO_DBG

  printf("Cartoon ENO interpolation.\n"
         "Input: order %d x0 %g dx %g x %g.\n",
         order, x0, dx, x);
  printf("Data is\n");
  for (i = 0; i < 2 * order + 1; ++i)
  {
    printf(" %g", y[i]);
  }
  printf("\nEntries used were\n");
  printf("0: %g   1: %g   2: %g\n",
         undiv_diff[ii-r][0], undiv_diff[ii-r+1][1], undiv_diff[ii-r+2][2]);
  printf("Parameters are fp_i %g fp_ii %g x_rel %g ii %d r %d\n",
         fp_i, fp_ii, x_rel, ii, r);
  printf("Result is %g\n", result);
  
#endif

  return result;  

}
