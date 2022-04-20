/*@@
   @file      SetGrid.c
   @date      April 2002
   @author    Denis Pollney
   @desc
              Reset the grid sizes before they are set by CartGrid3D, so
              that they have sensible cartoon-compatible values.
   @enddesc
   @version   $Id: SetGrid.c,v 1.7 2008/01/16 15:04:27 schnetter Exp $
 @@*/

#include <stdio.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameter.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header: /cactusdevcvs/BetaThorns/Cartoon2D/src/SetGrid.c,v 1.7 2008/01/16 15:04:27 schnetter Exp $";

CCTK_FILEVERSION(BetaThorns_Cartoon2D_SetGrid_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

int Cartoon2D_SetGrid(void);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    Cartoon_SetGrid
   @date       April 2002
   @author     Denis Pollney
   @desc
               Resets the grid dimensions in a cartoon-compatible way if
               the user has specified grid::type="byspacing".

               Generally, "byspacing" would put the z-axis in the middle
               of the x range. However, for cartoon, it is required to
               be at the edge of the range with only some zombie-zones
               extending across. This routine ensures that this is the
               case, by fixing the x range and re-specifying nx so that
               the dx value which was specified by the user is unchanged.
               The y range is reset to be a width of 1 plus space for the
               zombie points.

               In the process, a couple of other parameters such as
               grid::bitant_plane and grid::avoid_originy are also checked
               to ensure that they are cartoon-compatible.

               If the grid type is specified to be "byrange", then this
               routine does nothing and it is up to the user to choose
               appropriate ranges for what they want to do.

               Note that this routine currently needs to be scheduled at
               CCTK_RECOVER_PARAMETERS, because this is the only place where
               it is still possible to modify non-steerable parameters.
   @enddesc
   @@*/
int Cartoon2D_SetGrid(void)
{
  DECLARE_CCTK_PARAMETERS

  const CCTK_INT *cctk_int_ptr;
  const CCTK_REAL *cctk_real_ptr;
  int nx, ny, nz, ghost_size_y, avoid_y0, cartoon_ny;
  double dx, dy, dz;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  char p_val[80];
  const char *domain, *plane;


  CCTK_WARN (0, "Erik Schnetter, 2006-05-11: This routine does not work.  The schedule bin RECOVER_PARAMETERS is special; not all routines in this bin are executed.");

  /*
   * Determine the y ghostzone size.
   */
  cctk_int_ptr = (const CCTK_INT *) CCTK_ParameterGet("ghost_size_y", "pugh", NULL);
  if (cctk_int_ptr == NULL)
    CCTK_WARN(0, "pugh::ghost_size_y must be set explicitly");

  ghost_size_y = *cctk_int_ptr;
  cartoon_ny = 2 * ghost_size_y + 1;

  /*
   * Get the x,y,z grid sizes.
   */
  cctk_int_ptr = (const CCTK_INT *) CCTK_ParameterGet("global_nsize", "pugh", NULL);
  if ((cctk_int_ptr == NULL) || (*cctk_int_ptr == -1))
  {
    nx = *(const CCTK_INT *) CCTK_ParameterGet("global_nx", "pugh", NULL);
    ny = *(const CCTK_INT *) CCTK_ParameterGet("global_ny", "pugh", NULL);
    nz = *(const CCTK_INT *) CCTK_ParameterGet("global_nz", "pugh", NULL);
  }
  else
  {
    nx = ny = nz = *cctk_int_ptr;
  }

  /*
   * Reset the y grid size to be one layer thick, plus twice the ghostzone
   * size.
   */
  if (ny != cartoon_ny)
  {
    CCTK_VInfo(CCTK_THORNSTRING, "Resetting pugh::global_ny to %d", cartoon_ny);
    sprintf(p_val, "%d", cartoon_ny);
    CCTK_ParameterSet("global_ny", "pugh",  p_val);
  }


  /*
   * Reset the grid sizes if the user has specified a "byspacing" grid.
   */
  if (strcmp(CCTK_ParameterValString("type", "cartgrid3d"), "byrange"))
  {
    /*
     * We need to set the grid spacings explicitly, so from now on
     * consider the grid to be 'byrange'.
     */
    CCTK_ParameterSet("type", "cartgrid3d", "byrange");

    /*
     * Get the x,y,z grid spacing.
     */
    cctk_real_ptr = (const CCTK_REAL *) CCTK_ParameterGet("dxyz", "cartgrid3d", NULL);
    if (cctk_real_ptr == NULL)
    {
      dx = *(const CCTK_REAL *) CCTK_ParameterGet("dx", "cartgrid3d", NULL);
      dy = *(const CCTK_REAL *) CCTK_ParameterGet("dy", "cartgrid3d", NULL);
      dz = *(const CCTK_REAL *) CCTK_ParameterGet("dz", "cartgrid3d", NULL);
    }
    else
    {
      dx = dy = dz = *cctk_real_ptr;
    }

    /*
     * Set the x grid dimensions.
     */
    nx += ghost_size_y + 1;
    CCTK_VInfo(CCTK_THORNSTRING, "Adding pugh::ghost_size_y+1 to pugh::global_nx");
    sprintf(p_val, "%d", nx);
    CCTK_ParameterSet("global_nx", "pugh",  p_val);

    xmin = -(ghost_size_y-0.5) * dx;
    xmax = xmin + (nx-1) * dx;

    CCTK_VInfo(CCTK_THORNSTRING, "Setting x-range to [%f, %f]", xmin, xmax);

    sprintf(p_val, "%f", xmin);
    if (CCTK_ParameterSet("xmin", "cartgrid3d",  p_val) < 0)
    {
       CCTK_WARN(0,"Error setting parameter");
    }
    sprintf(p_val, "%f", xmax);
    if (CCTK_ParameterSet("xmax", "cartgrid3d",  p_val) < 0)
    {
       CCTK_WARN(0,"Error setting parameter");
    }

    /*
     * Set the y grid dimensions.
     */
    avoid_y0 = *(const CCTK_INT *) CCTK_ParameterGet("avoid_originy", "cartgrid3d", NULL);
    if (avoid_y0)
      CCTK_WARN(0, "Cartoon2D requires grid::avoid_originy=\"no\"");

    ymax = ghost_size_y * dy;
    ymin = -ymax;

    CCTK_VInfo(CCTK_THORNSTRING, "Setting y-range to [%f, %f]", ymin, ymax);

    sprintf(p_val, "%f", ymin);
    if (CCTK_ParameterSet("ymin", "cartgrid3d",  p_val) < 0)
    {
       CCTK_WARN(0,"Error setting parameter");
    }
    sprintf(p_val, "%f", ymax);
    if (CCTK_ParameterSet("ymax", "cartgrid3d",  p_val) < 0)
    {
       CCTK_WARN(0,"Error setting parameter");
    }

    /*
     * Set the z grid dimensions for bitant mode.
     */

    zmin = zmax = 0;
    domain = *((const char * const *) CCTK_ParameterGet("domain", "cartgrid3d", NULL));
    if (!strcmp(domain, "bitant"))
    {
      plane = *(const char * const *) CCTK_ParameterGet("bitant_plane", "cartgrid3d", NULL);
      if (strcmp(plane, "xy"))
        CCTK_WARN(0, "Cartoon2D requires grid::bitant_plane=\"xy\"");

      cctk_int_ptr = (const CCTK_INT *) CCTK_ParameterGet("ghost_size", "pugh", NULL);
      if ((cctk_int_ptr == NULL) || (*cctk_int_ptr < 0))
        cctk_int_ptr = (const CCTK_INT *) CCTK_ParameterGet("ghost_size_z", "pugh", NULL);

      nz += 1;
      CCTK_VInfo(CCTK_THORNSTRING, "Increasing PUGH::global_nz by 1");
      sprintf(p_val, "%d", nz);
      CCTK_ParameterSet("global_nz", "pugh",  p_val);

      zmax = (nz-0.5) * dz;
      zmin = -zmax;

      if ((cctk_int_ptr != NULL) && (*cctk_int_ptr > 0))
      {
        nz += (*cctk_int_ptr);
        CCTK_VInfo(CCTK_THORNSTRING, "Increasing PUGH::global_nz by "
                                     "PUGH::ghost_size_z");
        sprintf(p_val, "%d", nz);
        CCTK_ParameterSet("global_nz", "pugh",  p_val);
      }
    }
    else if (!strcmp(domain, "full"))
    {
      nz += 2;
      CCTK_VInfo(CCTK_THORNSTRING, "Increasing PUGH::global_nz by 2");
      sprintf(p_val, "%d", nz);
      CCTK_ParameterSet("global_nz", "pugh",  p_val);

      zmax = (nz-1) * dz/2;
      zmin = -zmax;
    }
    else
    {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Cartoon2D is unable to handle grid::domain=\"%s\"",
                   domain);
    }

    CCTK_VInfo(CCTK_THORNSTRING, "Setting z-range to [%f,%f]", zmin, zmax);

    sprintf(p_val, "%f", zmin);
    if (CCTK_ParameterSet("zmin", "cartgrid3d",  p_val) < 0)
    {
       CCTK_WARN(0,"Error setting parameter");
    }
    sprintf(p_val, "%f", zmax);
    if (CCTK_ParameterSet("zmax", "cartgrid3d",  p_val) < 0)
    {
       CCTK_WARN(0,"Error setting parameter");
    }
  }

  return 0;
}
