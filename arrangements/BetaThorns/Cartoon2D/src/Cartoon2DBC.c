/*@@
   @file      Cartoon2DBC.c
   @date      Thu Nov  4 13:35:00 MET 1999
   @author    Sai Iyer
   @desc 
              Apply Cartoon2D boundary conditions
              An implementation of Steve Brandt's idea for doing
              axisymmetry with 3d stencils.
              Cactus 4 version of Bernd Bruegmann's code.             
   @enddesc
 @@*/

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_FortranString.h"

#include "Cartoon2D_tensors.h"
#include "SpaceMask.h"

static const char *rcsid = "$Id: Cartoon2DBC.c,v 1.26 2008/01/16 15:05:14 schnetter Exp $";

CCTK_FILEVERSION(BetaThorns_Cartoon2D_Cartoon2DBC_c);

int BndCartoon2DGN(const cGH *GH, int tensortype, const char *group);
int BndCartoon2DVN(const cGH *GH, int tensortype, const char *var);
int BndCartoon2DVI(const cGH *GH, int tensortype, int prolongtype, int vi);
static CCTK_REAL Cartoon2DInterp(const cGH *GH, CCTK_REAL *f, int i, int di,
                             int ijk, CCTK_REAL x);
static CCTK_REAL Cartoon2DInterp_ENO(const cGH *cctkGH, 
                                    CCTK_REAL *f, int i, int di, 
                                    int ijk, CCTK_REAL x);
CCTK_REAL interpolate_local(int order, CCTK_REAL x0, CCTK_REAL dx, 
                            CCTK_REAL y[], CCTK_REAL x, 
                            CCTK_INT excised_points[]);
CCTK_REAL interpolate_eno(CCTK_INT order, CCTK_REAL x0, CCTK_REAL dx,
                          CCTK_REAL y[], CCTK_REAL x,
                          CCTK_INT excised_points[]);
void CCTK_FCALL CCTK_FNAME(BndCartoon2DVI)(int *retval, const cGH **GH, 
                int *tensortype, int *vi);
void CCTK_FCALL CCTK_FNAME(BndCartoon2DVN)(int *retval, const cGH **GH,
                int *tensortype, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(BndCartoon2DGN)
                          (int *retval, const cGH **GH, int *tensortype, ONE_FORTSTRING_ARG);


void Cartoon2D_InitExcisionVars(CCTK_ARGUMENTS);
void SetupExcisionVars(CCTK_ARGUMENTS);

/* set boundaries of a grid tensor assuming axisymmetry 
   - handles lower boundary in x
   - does not touch other boundaries
   - coordinates and rotation coefficients are independent of z and should
     be precomputed
   - this is also true for the constants in interpolator, but this may 
     be too messy
   - minimizes conceptual warpage, not computationally optimized
*/
/* uses rotation matrix and its inverse as linear transformation on
   arbitrary tensor indices -- I consider this a good compromise
   between doing index loops versus using explicit formulas in
   cos(phi) etc with messy signs 
*/

/* BndCartoon2DVI is called only on the first variable of a group.  It
   then is automatically applied to the remaining variables in the
   group, using the tensortype argument to determine how many there
   are and how to treat them. */

int BndCartoon2DVI(const cGH *cctkGH, int tensortype, int prolongtype, int vi)

{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  char * groupname;

  int i, j, k, di, dj, ijk;
  int jj, n, s;
  int lnx, lnz, ny;
  int gi, var0;
  int n_vars;
  int timelevel;
  CCTK_REAL x, y, rho;
  CCTK_REAL lx0, dx0, dy0;
  CCTK_REAL rxx, rxy, ryx, ryy;
  CCTK_REAL sxx, sxy, syx, syy;
  CCTK_REAL f[100], fx, fy, fz, fxx, fxy, fxz, fyy, fyz, fzz;
  CCTK_REAL *t, *tx, *ty, *tz, *txx, *txy, *txz, *tyy, *tyz, *tzz;

  gi = CCTK_GroupIndexFromVarI (vi);
  assert (gi >= 0);
  n_vars = CCTK_NumVarsInGroupI (gi);
  assert (n_vars > 0);
  var0 = CCTK_FirstVarIndexI (gi);
  assert (var0 >= 0);

  if (n_vars > 100)
  {
    groupname = CCTK_GroupName (gi);
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Group \"%s\" has more than 100 variables -- cannot apply Cartoon boundary condition",
                groupname);
    free (groupname);
    return -1;
  }

  switch (tensortype) {
  case TENSORTYPE_SCALAR:
    break;
  case TENSORTYPE_U:
    if (n_vars != 3)
    {
      groupname = CCTK_GroupName (gi);
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Group \"%s\" is declared with the tensor type U, but does not contain 3 variables -- cannot apply Cartoon boundary condition",
                  groupname);
      free (groupname);
      return -1;
    }
    break;
  case TENSORTYPE_DDSYM:
    if (n_vars != 6)
    {
      groupname = CCTK_GroupName (gi);
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Group \"%s\" is declared with the tensor type DDSYM, but does not contain 6 variables -- cannot apply Cartoon boundary condition",
                  groupname);
      free (groupname);
      return -1;
    }
    break;
  default:
    CCTK_WARN(0,"Tensor type not recognized by Cartoon2D.");
  }
    
  s = cctkGH->cctk_nghostzones[0];
  lnx = cctkGH->cctk_lsh[0];
  lnz = cctkGH->cctk_lsh[2];
  ny  = cctkGH->cctk_gsh[1];

  dx0 = CCTK_DELTA_SPACE(0);
  dy0 = CCTK_DELTA_SPACE(1);
  lx0 = CCTK_ORIGIN_SPACE(0) + dx0 * cctk_lbnd[0];

  timelevel = 0; 

  /* Cactus loop in C:  grid points with y = 0 */
  /* compare thorn_BAM_Elliptic */
  /* strides used in stencils, should be provided by cactus */
  di = 1;
  dj = lnx;

  /* y = 0 */
  assert (ny % 2 == 1);
  j = ny/2;

  /* make sure that the input data is synced */
  CCTK_SyncGroupI(cctkGH, gi);

  /* z-direction: include lower and upper boundary */
  for (k = 0; k < lnz; k++)

  /* y-direction: as many zombies as the evolution stencil needs */
#if 0
  for (jj = 1, dj = jj*lnx; jj <= ny/2; jj++, dj = jj*lnx) 
#else
  for (jj = 0, dj = jj*lnx; jj <= ny/2; jj++, dj = jj*lnx) 
#endif

  /* x-direction: zombie for x < 0, including upper boundary for extrapol */
#if 0
  for (i = -s; i < lnx; i++) {
#else
  for (i = 0; i < lnx; i++)
    if (! (i > s && jj == 0)) {
#endif

    /* index into linear memory */
    ijk = CCTK_GFINDEX3D(cctkGH,i,j,k);
   
    /* what a neat way to hack in the zombie for x < 0, y = 0 */
    if (i < 0) {i += s; ijk += s; dj = 0;}

    /* compute coordinates (could also use Cactus grid function) */
    x = lx0 + dx0 * i;
#if 0
    y = (dj) ? dy0 * jj : 0;
#else
    y = dy0 * jj;
#endif
    rho = sqrt(x*x + y*y);

    /* compute rotation matrix
       note that this also works for x <= 0 (at lower boundary in x)
       note that rho is nonzero by definition if cactus checks dy ... 
    */
    rxx = x/rho;
    ryx = y/rho;
    rxy = -ryx;
    ryy = rxx;

    /* inverse rotation matrix, assuming detr = 1 */
    sxx = ryy;
    syx = -ryx;
    sxy = -rxy;
    syy = rxx;

    /* interpolate grid functions */
    switch (prolongtype)
    {
      case PROLONG_NONE:
        /* no-op */
        return (0);
        break;
      case PROLONG_LAGRANGE:
        for (n = 0; n < n_vars; n++)
          f[n] = Cartoon2DInterp(cctkGH, cctkGH->data[var0+n][timelevel], 
                                 i, di, ijk, rho); 
        break;
      case PROLONG_ENO:
        for (n = 0; n < n_vars; n++)
          f[n] = Cartoon2DInterp_ENO(cctkGH, cctkGH->data[var0+n][timelevel], i, 
                                     di, ijk, rho); 
        break;
      default:
        CCTK_WARN(0, "Prolongation type not recognized by Cartoon2D.");
    }

    /* rotate grid tensor by matrix multiplication */
    if (tensortype == TENSORTYPE_SCALAR) {
      /* Scalar groups can have arbitrary many components; these
         "components" are then all scalars.  */
      for (n = 0; n < n_vars; n++)
      {
        t = cctkGH->data[var0 + n][timelevel];
        t[ijk+dj] = f[n];
        t[ijk-dj] = f[n];
      }
    }
    else if (tensortype == TENSORTYPE_U) {
      tx = cctkGH->data[vi][timelevel];
      ty = cctkGH->data[vi+1][timelevel];
      tz = cctkGH->data[vi+2][timelevel];
      fx = f[0];
      fy = f[1];
      fz = f[2];

      tx[ijk+dj] = rxx * fx + rxy * fy;
      ty[ijk+dj] = ryx * fx + ryy * fy;
      tz[ijk+dj] = fz;

      tx[ijk-dj] = sxx * fx + sxy * fy;
      ty[ijk-dj] = syx * fx + syy * fy;
      tz[ijk-dj] = fz;
    }
    else if (tensortype == TENSORTYPE_DDSYM) {
      txx = cctkGH->data[vi][timelevel];
      txy = cctkGH->data[vi+1][timelevel];
      txz = cctkGH->data[vi+2][timelevel];
      tyy = cctkGH->data[vi+3][timelevel];
      tyz = cctkGH->data[vi+4][timelevel];
      tzz = cctkGH->data[vi+5][timelevel];
      fxx = f[0];
      fxy = f[1];
      fxz = f[2];
      fyy = f[3];
      fyz = f[4];
      fzz = f[5];
      
      txx[ijk+dj] = fxx*sxx*sxx + 2*fxy*sxx*syx + fyy*syx*syx;
      tyy[ijk+dj] = fxx*sxy*sxy + 2*fxy*sxy*syy + fyy*syy*syy;
      txy[ijk+dj] = fxx*sxx*sxy + fxy*(sxy*syx+sxx*syy) + fyy*syx*syy;
      txz[ijk+dj] = fxz*sxx + fyz*syx;
      tyz[ijk+dj] = fxz*sxy + fyz*syy;
      tzz[ijk+dj] = fzz;

      txx[ijk-dj] = fxx*rxx*rxx + 2*fxy*rxx*ryx + fyy*ryx*ryx;
      tyy[ijk-dj] = fxx*rxy*rxy + 2*fxy*rxy*ryy + fyy*ryy*ryy;
      txy[ijk-dj] = fxx*rxx*rxy + fxy*(rxy*ryx+rxx*ryy) + fyy*ryx*ryy;
      txz[ijk-dj] = fxz*rxx + fyz*ryx;
      tyz[ijk-dj] = fxz*rxy + fyz*ryy;
      tzz[ijk-dj] = fzz;
    } 
    else {
      return(-1);
    }

#if 0
    /* what a neat way to hack out the zombies for x < 0, y = 0 */
    if (dj == 0) {i -= s; ijk -= s; dj = jj*lnx;}
#endif
  }

  /* syncs needed after interpolation (actually only for x direction) */
  CCTK_SyncGroupI(cctkGH, gi);

  return(0);
}

void CCTK_FCALL CCTK_FNAME(BndCartoon2DVI)
                          (int *retval, const cGH **GH, int *tensortype, int *vi)
{
  *retval = BndCartoon2DVI(*GH, *tensortype, PROLONG_LAGRANGE, *vi);
}

int BndCartoon2DVN(const cGH *GH, int tensortype, const char *inpvarname)
{
  DECLARE_CCTK_PARAMETERS

  int vi;

  if (verbose) printf("Cartoon2D called for %s\n", inpvarname);
  vi = CCTK_VarIndex(inpvarname);

  return(BndCartoon2DVI(GH, tensortype, PROLONG_LAGRANGE, vi));
}

void CCTK_FCALL CCTK_FNAME(BndCartoon2DVN)
                          (int *retval, const cGH **GH, int *tensortype, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(impvarname)
  *retval = BndCartoon2DVN(*GH, *tensortype, impvarname);
  free(impvarname);
}

int BndCartoon2DGN(const cGH *GH, int tensortype, const char *group)
{
  DECLARE_CCTK_PARAMETERS

  if (verbose && GH->cctk_iteration==1) CCTK_VInfo(CCTK_THORNSTRING,"Cartoon2D called for %s (message appears only once)", group);
  
  return(BndCartoon2DVI(GH, tensortype, PROLONG_LAGRANGE, 
                        CCTK_FirstVarIndex(group)));
}

void CCTK_FCALL CCTK_FNAME(BndCartoon2DGN)
                          (int *retval, const cGH **GH, int *tensortype, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(group)
  *retval = BndCartoon2DGN(*GH, *tensortype, group);
  free(group);
}

/* interpolation on x-axis */
static CCTK_REAL Cartoon2DInterp(const cGH *cctkGH, CCTK_REAL *f, 
                                 int i, int di, int ijk, CCTK_REAL x)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL x0;
  CCTK_REAL lx0, dx0;
  CCTK_REAL y[6], r;
  int lnx;
  int n, offset;

  CCTK_REAL *old_mask;
  CCTK_INT *new_mask;

  CCTK_INT excised_points[6];

  if (*excision_active < 0)
  {
    SetupExcisionVars((cGH*)cctkGH);
  }

  if (old_excision)
  {
    old_mask = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 0, *old_mask_vi);
  }
  else
  {
    old_mask = NULL;
  }
  
  if (new_excision)
  {
    new_mask = (CCTK_INT*)CCTK_VarDataPtrI(cctkGH, 0, *new_mask_vi);
  }
  else
  {
    new_mask = NULL;
  }

  /* 
     If this boundary point is excised then we return 
     without doing anything 
  */

  if (*excision_active)
  {
    if ((*excision_active == 1)||(*excision_active == 3))
    {
      if (old_mask[ijk] < 0.75)
      {
        return f[ijk];
      }
    }
    if ((*excision_active == 2)||(*excision_active == 3))
    {
      if (SpaceMask_CheckStateBits(new_mask, (ijk),
                                   *new_excision_field,
                                   *new_excision_descriptor))
      {
        return f[ijk];
      }
    }
  }  

  lnx = cctkGH->cctk_lsh[0];
    
  dx0 = CCTK_DELTA_SPACE(0);
  lx0 = CCTK_ORIGIN_SPACE(0) + dx0 * cctk_lbnd[0];

  /* find i such that x(i) < x <= x(i+1)
     for rotation on entry always x > x(i), but sometimes also x > x(i+1) */
  while (x > lx0 + dx0 * (i+1)) {i++; ijk++;}

  /* first attempt to interface to JT's interpolator */

  /* offset of stencil, note that rotation leads to x close to x0 
       for large x */
  offset = order/2;
  
  /* shift stencil at boundaries */
  /* note that for simplicity we don't distinguish between true boundaries
     and decomposition boundaries: the sync fixes that! */
  if (i-offset < 0) offset = i;
  if (i-offset+order >= lnx) offset = i+order-lnx+1;
  
  /* fill in info */
  /* fills in old data near axis, but as long as stencil at axis is 
     centered, no interpolation happens anyway */
  x0 = lx0 + dx0 * (i-offset);
  for (n = 0; n <= order; n++) {
    y[n] = f[ijk-offset+n];
    excised_points[n] = 0;
    if (*excision_active)
    {
      if ((*excision_active == 1)||(*excision_active == 3))
      {
        excised_points[n] = (old_mask[ijk-offset+n] > 0.75) ? 0 : 1;
      }
      if ((*excision_active == 2)||(*excision_active == 3))
      {
        excised_points[n] |= 
          (SpaceMask_CheckStateBits(new_mask, (ijk-offset+n),
                                    *new_excision_field,
                                    *new_excision_descriptor));
      }
    }
  }
  
  /* call interpolator */
  r = interpolate_local(order, x0, dx0, y, x, excised_points);
  
  return r;
}

static CCTK_REAL Cartoon2DInterp_ENO(const cGH *cctkGH, 
                                    CCTK_REAL *f, int i, int di, 
                                    int ijk, CCTK_REAL x)
{
  
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;
  
  CCTK_REAL x0;
  CCTK_REAL lx0, dx0;
  CCTK_REAL y[11], r;
  int lnx;
  int n, offset;
  int tmp_order = 2;

  CCTK_REAL *old_mask;
  CCTK_INT *new_mask;

  CCTK_INT excised_points[11];

  if (*excision_active < 0)
  {
    SetupExcisionVars((cGH*)cctkGH);
  }

  if (old_excision)
  {
    old_mask = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 0, *old_mask_vi);
  }
  else
  {
    old_mask = NULL;
  }
  
  if (new_excision)
  {
    new_mask = (CCTK_INT*)CCTK_VarDataPtrI(cctkGH, 0, *new_mask_vi);
  }
  else
  {
    new_mask = NULL;
  }

  /* 
     If this boundary point is excised then we return 
     without doing anything 
  */

  if (*excision_active)
  {
    if ((*excision_active == 1)||(*excision_active == 3))
    {
      if (old_mask[ijk] < 0.75)
      {
        return f[ijk];
      }
    }
    if ((*excision_active == 2)||(*excision_active == 3))
    {
      if (SpaceMask_CheckStateBits(new_mask, (ijk),
                                   *new_excision_field,
                                   *new_excision_descriptor))
      {
        return f[ijk];
      }
    }
  }
  
  lnx = cctkGH->cctk_lsh[0];
  
  dx0 = CCTK_DELTA_SPACE(0);
  lx0 = CCTK_ORIGIN_SPACE(0) + dx0 * cctk_lbnd[0];
  
  /* find i such that x(i) < x <= x(i+1)
     for rotation on entry always x > x(i), but sometimes also x > x(i+1) */
  while (x > lx0 + dx0 * (i+1)) {i++; ijk++;}
  
  /* first attempt to interface to JT's interpolator */
  
  /* offset of stencil, note that rotation leads to x close to x0 
     for large x */
  offset = eno_order;
/*   offset = tmp_order; */
  
  /* shift stencil at boundaries */
  /* note that for simplicity we don't distinguish between true boundaries
     and decomposition boundaries: the sync fixes that! */
  if (i-offset < 0) offset = i;
/*   if (i-offset+4 >= lnx) offset = i+5-lnx; */
/*   if (i-offset+2*tmp_order >= lnx) offset = i+2*tmp_order+1-lnx; */
  if (i-offset+2*eno_order >= lnx) offset = i+2*eno_order+1-lnx;
  
  /* fill in info */
  /* fills in old data near axis, but as long as stencil at axis is 
     centered, no interpolation happens anyway */  
  x0 = lx0 + dx0 * (i-offset);
/*   for (n = 0; n < 2 * tmp_order + 1; ++n)  */
  for (n = 0; n < 2 * eno_order + 1; ++n) 
  {
    y[n] = f[ijk-offset+n];
    excised_points[n] = 0;
    if (*excision_active)
    {
      if ((*excision_active == 1)||(*excision_active == 3))
      {
        excised_points[n] = (old_mask[ijk-offset+n] > 0.75) ? 1 : 0;
      }
      if ((*excision_active == 2)||(*excision_active == 3))
      {
        excised_points[n] |= 
          (SpaceMask_CheckStateBits(new_mask, (ijk-offset+n),
                                    *new_excision_field,
                                    *new_excision_descriptor));
      }
    }
  }
  
  /* call interpolator */
/*   r = interpolate_eno(tmp_order, x0, dx0, y, x); */
  r = interpolate_eno(eno_order, x0, dx0, y, x, excised_points);
  
  return r;
}

void Cartoon2D_InitExcisionVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *excision_active = -1;
  *old_mask_vi = -1;
  *new_mask_vi = -1;
  *new_excision_field = -1;
  *new_excision_descriptor = -1;

}

void SetupExcisionVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (old_excision)
  {
    *old_mask_vi = CCTK_VarIndex(old_style_excision_var);
  }
  else
  {
    *old_mask_vi = -1;
  }
  if (new_excision)
  {
    *new_mask_vi = CCTK_VarIndex(new_style_excision_var);
  }
  else
  {
    *new_mask_vi = -1;
  }

  if (*old_mask_vi > -1)
  {
    if (*new_mask_vi > -1)
    {
      *excision_active = 3;
      *new_excision_field = SpaceMask_GetTypeBits(new_mask_field_name);
      *new_excision_descriptor = SpaceMask_GetStateBits(new_mask_field_name,
                                                        new_mask_excised_name);
    }
    else
    {
      *excision_active = 1;
    }
  }
  else
  {
    if (*new_mask_vi > -1)
    {
      *excision_active = 2;
      *new_excision_field = SpaceMask_GetTypeBits(new_mask_field_name);
      *new_excision_descriptor = SpaceMask_GetStateBits(new_mask_field_name,
                                                        new_mask_excised_name);
    }
    else
    {
      *excision_active = 0;
    }
  }
  
}
